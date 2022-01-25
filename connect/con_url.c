/*  $Id: con_url.c,v 6.6 1999/04/09 22:27:24 vakatov Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Author:  Denis Vakatov
*
* File Description:
*  Implement CONNECTOR to hit an URL.
*
*  NOTE:  See in "connectr.h" for the detailed specification of the underlying
*         connector("CONNECTOR", "SConnectorTag") methods and structures.
*
* --------------------------------------------------------------------------
* $Log: con_url.c,v $
* Revision 6.6  1999/04/09 22:27:24  vakatov
* Added flag "URLC_URL_ENCODE_ARGS";  thus do not URL-encode CGI args by
* default
*
* Revision 6.5  1999/04/09 21:36:00  vakatov
* Split former CGI "args" into "path?args":
* - added "path" arg. to Ncbi_ConnectURL(), have it URL-encoded automagically;
* - added "path" field to SNetConnInfo
* - split ..._ENGINE_URL defs and config. names to ..._ENGINE_PATH / ARG
*
* Revision 6.4  1999/04/08 18:20:53  vakatov
* Added the on-the-fly URL-encoding/decoding functionality
* (see flags "URLC_URL_*" for URL_CreateConnector())
*
* Revision 6.3  1999/04/02 20:41:52  vakatov
* Got rid of occasional comment inside comment(in CVS Log)
*
* Revision 6.2  1999/04/01 21:52:11  vakatov
* Fixed for the change in spec:  "n_written/n_read" args in
* FConnectorWrite/Read to be non-NULL and "*n_written / *n_read" := 0
*
* Revision 6.1  1999/03/26 19:50:45  vakatov
* tiny type casts
*
* Revision 6.0  1999/03/25 23:04:58  vakatov
* Initial revision
*
* ==========================================================================
*/

#include <ncbi.h>
#include <ncbisock.h>
#include <ncbibuf.h>
#include <con_url.h>


/***********************************************************************
 *  INTERNAL -- Auxiliary types and static functions
 ***********************************************************************/

/* All internal data necessary to perform the (re)connect and i/o
 */
typedef struct {
  SNetConnInfo*   info;        /* connection configuration */
  Nlm_Char*       user_header; /* user-defined part of the HTTP header */
  SOCK            sock;        /* socket;  NULL if not in the "READ" mode */
  Nlm_Uint4       n_connect;   /* # of "WRITE"/"READ" auto-reconnect cycles */
  Nlm_Boolean     sure_flush;  /* send at least a header on CLOSE/re-CONNECT */
  Nlm_Boolean     strip_header;/* strip HTTP header from the response */
  Nlm_Boolean     url_input;   /* URL-decode all input data but HTTP header */
  Nlm_Boolean     url_output;  /* URL-encode all output data but HTTP header */
  Nlm_Boolean     encode_args; /* URL-encode "info->args" */
  Nlm_Boolean     first_read;  /* the 1st attempt to read after connect? */
  BUF             buf;         /* storage to accumulate output data */
  const STimeout* c_timeout;   /* NULL(infinite) or points to "cc_timeout" */
  STimeout        cc_timeout;  /* storage for a (finite) connect timeout */
  const STimeout* w_timeout;   /* NULL(infinite) or points to "ww_timeout" */
  STimeout        ww_timeout;  /* storage for a (finite) write timeout */
} SUrlConnector;


/* Reset the accumulated output data and close socket
 */
static void s_Close(SUrlConnector* uuu)
{
  ASSERT( uuu->sock );
  BUF_Read(uuu->buf, 0, UINT4_MAX);
  SOCK_Close(uuu->sock);
  uuu->sock = 0;
}


/* Flush the connector's accumulated output data.
 * NOTE:  on both success and any error, all internally stored accumulated
 *        output data will be discarded(BUF_Size(uuu->buf) --> zero)!
 */
static EConnStatus s_FlushData(SUrlConnector* uuu)
{
  ASSERT( uuu->sock );
  if ( !BUF_Size(uuu->buf) )
    return eCONN_Success;

  SOCK_SetTimeout(uuu->sock, eSOCK_OnWrite, uuu->w_timeout, 0, 0);
  do {
    char          buf[4096];
    Nlm_Uint4     n_written;
    Nlm_Uint4     size = BUF_Peek(uuu->buf, buf, sizeof(buf));
    ESOCK_ErrCode err_code = SOCK_Write(uuu->sock, buf, size, &n_written);

    /* on any error, just discard all data and fail */
    if (err_code != eSOCK_ESuccess) {
      BUF_Read(uuu->buf, 0, UINT4_MAX);
      return ESOCK2ECONN(err_code);
    }
    
    /* on success, discard the succesfully written data and continue */
    BUF_Read(uuu->buf, 0, n_written);
  } while ( BUF_Size(uuu->buf) );

  return eCONN_Success;
}


/* Connect to the "port:host" specified in "uuu->info", then
 * compose and form relevant HTTP header and
 * flush the accumulated output data("uuu->buf") after the HTTP header.
 * On error, all accumulated output data will be lost, and the connector
 * socket will be NULL.
 */
static EConnStatus s_ConnectAndSend(SUrlConnector* uuu)
{
  EConnStatus status;
  ASSERT( !uuu->sock );

  /* connect & send HTTP header */
  uuu->sock = Ncbi_ConnectURL
    (uuu->info->host, uuu->info->port, uuu->info->path, uuu->info->args,
     (size_t)BUF_Size(uuu->buf), uuu->c_timeout, uuu->w_timeout,
     uuu->user_header, uuu->encode_args);
  if ( !uuu->sock ) {
    BUF_Read(uuu->buf, 0, UINT4_MAX);
    return eCONN_Unknown;
  }

  /* flush the accumulated output data */
  if ((status = s_FlushData(uuu)) != eCONN_Success) {
    if (status == eCONN_Timeout)
      status = eCONN_Closed;  /* fake it to avoid somebody to re-try */
    s_Close(uuu);
    return status;
  }

  uuu->first_read = TRUE;
  return eCONN_Success;
}


/* Send the accumulated output data(if any) to server, then close socket
 */
static void s_FlushAndClose(SUrlConnector* uuu)
{
  if ( uuu->sock ) {
    ASSERT( !BUF_Size(uuu->buf) );
    SOCK_Close(uuu->sock);
    uuu->sock = 0;
  }
  else if (uuu->sure_flush  ||  BUF_Size(uuu->buf)) {
    if (uuu->n_connect != UINT4_MAX  &&
        s_ConnectAndSend(uuu) == eCONN_Success)
      s_Close(uuu);
  }
}


/***********************************************************************
 *  INTERNAL -- "s_VT_*" functions for the "virt.table" of connector methods
 ***********************************************************************/

static const Nlm_Char* s_VT_GetType
(void* connector)
{
  return "HTTP_URL";
}


static EConnStatus s_VT_Connect
(void*           connector,
 const STimeout* timeout)
{
  SUrlConnector* uuu = (SUrlConnector*)connector;

  /* flush the accumulated output data and close current socket */
  s_FlushAndClose(uuu);

  /* reset the auto-reconnect counter */
  uuu->n_connect = (uuu->n_connect != UINT4_MAX-1) ? 0 : UINT4_MAX-1;

  /* NOTE: the real connect will be performed on the first "READ"(or "CLOSE"),
   * or on the "WAIT" on read -- see in "s_ConnectAndSend()";
   * here we just store the timeout value -- for the future connect */
  if ( timeout ) {
    uuu->cc_timeout = *timeout;
    uuu->c_timeout  = &uuu->cc_timeout;
  } else {
    uuu->c_timeout = 0;
  }

  return eCONN_Success;
}


static EConnStatus s_VT_Wait
(void*           connector,
 EConnDirection  direction,
 const STimeout* timeout)
{
  SUrlConnector* uuu = (SUrlConnector*)connector;
  if (uuu->n_connect == UINT4_MAX)
    return eCONN_Closed; /* no more i/o permitted */


  if (direction == eCONN_Read) {
    if ( !uuu->sock ) {
      EConnStatus status = s_ConnectAndSend(uuu);
      if (status != eCONN_Success)
        return status;
    }
    return ESOCK2ECONN( SOCK_Select(uuu->sock, eSOCK_OnRead, timeout) );
  } else if (direction == eCONN_Write) {
    return (uuu->sock  &&  uuu->n_connect == UINT4_MAX-1) ?
      eCONN_Closed : eCONN_Success;
  }

  return eCONN_InvalidArg;
}


#ifdef IMPLEMENTED__CONN_WaitAsync
static EConnStatus s_VT_WaitAsync
(void*                   connector,
 FConnectorAsyncHandler  func,
 SConnectorAsyncHandler* data)
{
  return eCONN_NotSupported;
}
#endif


static EConnStatus s_VT_Write
(void*           connector,
 const void*     buf,
 Nlm_Uint4       size,
 Nlm_Uint4*      n_written,
 const STimeout* timeout)
{
  SUrlConnector* uuu = (SUrlConnector*)connector;
  if (uuu->n_connect == UINT4_MAX)
    return eCONN_Closed; /* no more connects permitted */

  /* if trying to "WRITE" after "READ"... */
  if ( uuu->sock ) {
    uuu->n_connect++;
    if (uuu->n_connect == UINT4_MAX)
      return eCONN_Closed; /* no more connects permitted */
    /* close the socket, and thus switch to the "WRITE" mode */
    SOCK_Close(uuu->sock);
    uuu->sock = 0;
  }

  /* accumulate all output in the memory buffer */
  if ( uuu->url_output ) {
    /* with URL-encoding */
    Nlm_Uint4 dst_size = 3 * size;
    void* dst = Nlm_MemNew(dst_size);
    Nlm_Uint4 dst_written;
    URL_Encode(buf, size, n_written, dst, dst_size, &dst_written);
    ASSERT( *n_written == size );
    if ( !BUF_Write(&uuu->buf, dst, dst_written) ) {
      Nlm_MemFree(dst);
      return eCONN_Unknown;
    }
    Nlm_MemFree(dst);
  }
  else {
    /* "as is" (without URL-encoding) */
    if ( !BUF_Write(&uuu->buf, buf, size) )
      return eCONN_Unknown;
    *n_written = size;
  }

  /* store the write timeout */
  if ( timeout ) {
    uuu->ww_timeout = *timeout;
    uuu->w_timeout  = &uuu->ww_timeout;
  } else {
    uuu->w_timeout = 0;
  }

  return eCONN_Success;
}


static EConnStatus s_VT_Flush
(void*           connector,
 const STimeout* timeout)
{
  SUrlConnector* uuu = (SUrlConnector*)connector;

  /* The real flush will be performed on the first "READ"(or "CLOSE"),
   * or on the "WAIT".
   * We just store the write timeout here, that's all...
   */
  if ( timeout ) {
    uuu->ww_timeout = *timeout;
    uuu->w_timeout  = &uuu->ww_timeout;
  } else {
    uuu->w_timeout = 0;
  }

  return eCONN_Success;
}


static EConnStatus s_VT_Read
(void*           connector,
 void*           buf,
 Nlm_Uint4       size,
 Nlm_Uint4*      n_read,
 const STimeout* timeout)
{
  SUrlConnector* uuu = (SUrlConnector*)connector;
  if (uuu->n_connect == UINT4_MAX)
    return eCONN_Closed; /* no more connects permitted */

  if ( !uuu->sock ) {
    /* not in the "READ" mode yet... so "CONNECT" and "WRITE" first */
    EConnStatus status = s_ConnectAndSend(uuu);
    if (status != eCONN_Success)
      return status;
    ASSERT( uuu->sock );
  }

  /* if it is a "fake" read -- just to connect and flush the output data */
  if ( !size )
    return eCONN_Success;

  /* set timeout */
  SOCK_SetTimeout(uuu->sock, eSOCK_OnRead, timeout, 0, 0);

  /* skip HTTP header */
  if ( uuu->first_read ) {
    uuu->first_read = FALSE;
    if ( uuu->strip_header ) {
      EConnStatus status = SOCK_StripToPattern(uuu->sock, "\r\n\r\n", 4, 0, 0);
      if (status != eCONN_Success)
        return status;
    }
  }

  /* just read, no URL-decoding */
  if ( !uuu->url_input )
    return ESOCK2ECONN( SOCK_Read(uuu->sock, buf, size, n_read) );

  /* read and URL-decode */
  {{
    EConnStatus status;
    Nlm_Uint4   n_peeked, n_decoded;
    Nlm_Uint4   peek_size = 3 * size;
    void*       peek_buf = Nlm_MemNew(peek_size);

    /* peek the data */
    status = ESOCK2ECONN(SOCK_Peek(uuu->sock, peek_buf, peek_size, &n_peeked));
    if (status != eCONN_Success) {
      Nlm_MemFree(peek_buf);
      return status;
    }

    /* decode, then discard the successfully decoded data from the input */
    if ( URL_Decode(peek_buf, n_peeked, &n_decoded, buf, size, n_read) ) {
      Nlm_Uint4 x_read;
      SOCK_Read(uuu->sock, peek_buf, n_decoded, &x_read);
      ASSERT( x_read == n_decoded );
      status = eCONN_Success;
    }
    else {
      status = eCONN_Unknown;
    }

    Nlm_MemFree(peek_buf);
    return status;
  }}
}


static EConnStatus s_VT_Close
(CONNECTOR       connector,
 const STimeout* timeout)
{
  SUrlConnector* uuu = (SUrlConnector*)connector->handle;

  s_FlushAndClose(uuu);
  NetConnInfo_Destroy(&uuu->info);
  Nlm_MemFree(uuu->user_header);
  BUF_Destroy(uuu->buf);
  Nlm_MemFree(uuu);
  Nlm_MemFree(connector);
  return eCONN_Success;
}



/***********************************************************************
 *  EXTERNAL -- the connector's "constructor"
 ***********************************************************************/


NLM_EXTERN CONNECTOR URL_CreateConnector
(const SNetConnInfo* info,
 const Nlm_Char*     user_header,
 URLC_Flags          flags)
{
  CONNECTOR      ccc = (SConnector   *)Nlm_MemNew(sizeof(SConnector   ));
  SUrlConnector* uuu = (SUrlConnector*)Nlm_MemNew(sizeof(SUrlConnector));

  /* initialize internal data structures */
  uuu->info = info ? NetConnInfo_Clone(info) : NetConnInfo_Create(0, 0);
  NetConnInfo_AdjustForCernProxy(uuu->info);
  uuu->user_header = user_header ? Nlm_StringSave(user_header) : 0;
  uuu->sock = 0;
  uuu->n_connect    = (flags & URLC_AUTO_RECONNECT) ? 0 : (UINT4_MAX - 1);
  uuu->sure_flush   = (Nlm_Boolean)(flags & URLC_SURE_FLUSH);
  uuu->url_input    = (Nlm_Boolean)(flags & URLC_URL_DECODE_INP);
  uuu->url_output   = (Nlm_Boolean)(flags & URLC_URL_ENCODE_OUT);
  uuu->encode_args  = (Nlm_Boolean)(flags & URLC_URL_ENCODE_ARGS);
  uuu->strip_header = (Nlm_Boolean)
     (uuu->url_input  ||  !(flags & URLC_HTTP_HEADER));
    
  uuu->c_timeout    = 0;

  /* initialize handle */
  ccc->handle = uuu;

  /* initialize virtual table */ 
  ccc->vtable.get_type    = s_VT_GetType;
  ccc->vtable.connect     = s_VT_Connect;
  ccc->vtable.wait        = s_VT_Wait;
#ifdef IMPLEMENTED__CONN_WaitAsync
  ccc->vtable.wait_async  = s_VT_WaitAsync;
#endif
  ccc->vtable.write       = s_VT_Write;
  ccc->vtable.flush       = s_VT_Flush;
  ccc->vtable.read        = s_VT_Read;
  ccc->vtable.close       = s_VT_Close;

  /* done */
  return ccc;
}


#ifdef TEST_MODULE__CON_URL

#include <conntest.h>

extern
#ifdef __cplusplus
"C"
#endif
Int2 Main(void)
{
  SNetConnInfo* info = NetConnInfo_Create(0, 0);
  const Char*   user_header = 0;
  Uint4         flags;
  STimeout      timeout;
  CONNECTOR     connector;
  FILE*         log_file;

  if ( info->host )
    MemFree(info->host);
  info->host = StringSave("ray.nlm.nih.gov");

  info->port = 6224;

  if ( info->path )
    MemFree(info->path);
  info->path = StringSave("/cgi-bin/tools/vakatov/con_url.cgi");

  if ( info->args )
    MemFree(info->args);
  info->args = StringSave("arg1 arg2 arg3");

  timeout.sec  = 5;
  timeout.usec = 123456;
  log_file = FileOpen("con_url.out", "wb");

  flags = URLC_HTTP_HEADER | URLC_URL_CODEC | URLC_URL_ENCODE_ARGS;
  connector = URL_CreateConnector(info, user_header, flags);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_SINGLE_BOUNCE_PRINT);

  MemFree(info->args);
  info->args = StringSave("arg1+arg2+arg3");

  flags = 0;
  connector = URL_CreateConnector(info, user_header, flags);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_SINGLE_BOUNCE_CHECK);

  flags = URLC_AUTO_RECONNECT;
  connector = URL_CreateConnector(info, user_header, flags);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_ALL);

  flags = URLC_AUTO_RECONNECT | URLC_URL_CODEC;
  connector = URL_CreateConnector(info, user_header, flags);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_ALL);

  FileClose(log_file);
  NetConnInfo_Destroy(&info);
  return 0;
}

#endif  /* TEST_MODULE__CON_URL */
