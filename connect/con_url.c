/*  $Id: con_url.c,v 6.13 1999/12/09 21:07:17 vakatov Exp $
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
* Revision 6.13  1999/12/09 21:07:17  vakatov
* #ifdef TEST__HitURL::Main():  use VERIFY rather then ASSERT
*
* Revision 6.12  1999/08/04 21:06:44  vakatov
* s_VT_Read() -- check for the returned HTTP reply status
*
* Revision 6.11  1999/07/26 18:03:25  vakatov
* Added new test (see #ifdef TEST__HitURL)
* Redesigned the old test (see #ifdef TEST_MODULE__CON_URL) to use
* transient app.parameters rather than the "info" structure directly
*
* Revision 6.10  1999/07/16 22:21:56  vakatov
* + URL_CreateConnectorEx() to switch the peer URL and other parameters
* "on-the-fly", on every reconnect
* Retry(up to "info->conn_try" times) when cannot establish the connection.
* Printout the "info" content in the "info->debug_printout" mode.
*
* Revision 6.9  1999/07/12 16:39:07  vakatov
* s_VT_Read() -- added debug printout of the HTTP header
*
* Revision 6.8  1999/07/09 22:55:39  vakatov
* s_VT_Read() -- using "SOCK_Eof()", catch EOF on PEEK and so escape an
* infinite loop
*
* Revision 6.7  1999/06/28 16:28:01  vakatov
* SNetConnInfo:: separated the HTTP and the CERL-like(non-transparent)
* firewall proxies;  renamed config./env. parameters accordingly:
*   SRV_HTTP_PROXY_HOST/PORT replaced the former SRV_PROXY_HOST/PORT
*   SRV_PROXY_HOST now specifies the host name of non-transparent CERN proxy
*   SRV_PROXY_PORT is obsolete
*   SRV_PROXY_TRANSPARENT is obsolete
* Also:  NetConnInfo_AdjustForCernProxy --> ...HttpProxy
*
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
  FAdjustInfo     adjust_info;
  void*           adjust_data;
  FAdjustCleanup  adjust_cleanup;
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


/* Adjust the "uuu->info" with "uuu->adjust_info()";
 * connect to the "port:host" specified in "uuu->info", then
 * compose and form relevant HTTP header and
 * flush the accumulated output data("uuu->buf") after the HTTP header.
 * On error, all accumulated output data will be lost, and the connector
 * socket will be NULL.
 */
static EConnStatus s_ConnectAndSend(SUrlConnector* uuu)
{
  EConnStatus status;
  Nlm_Uint4 i;
  Nlm_Uint4 conn_try = uuu->info->conn_try ? uuu->info->conn_try : 1;
  ASSERT( !uuu->sock );

  /* the re-try loop... */
  for (i = 0;  i < conn_try  &&  !uuu->sock;  i++) {
    /* adjust the connection info */
    if ( uuu->adjust_info ) {
      uuu->adjust_info(uuu->info, uuu->adjust_data, i);
      if ( uuu->info->debug_printout )
        NetConnInfo_Print(uuu->info, stderr);
    }

    /* connect & send HTTP header */
    uuu->sock = Ncbi_ConnectURL
      (uuu->info->host, uuu->info->port, uuu->info->path, uuu->info->args,
       (size_t)BUF_Size(uuu->buf), uuu->c_timeout, uuu->w_timeout,
       uuu->user_header, uuu->encode_args);
  }

  /* ? error */
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

  /* first read::  skip HTTP header (and check reply status) */
  if ( uuu->first_read ) {
    uuu->first_read = FALSE;
    if ( uuu->strip_header ) {
      Nlm_Boolean server_error = FALSE;
      EConnStatus status;
      BUF         buf = 0;

      /* check status (assume the reply status is in the first line) */
      status = SOCK_StripToPattern(uuu->sock, "\r\n", 2, &buf, 0);
      if (status == eCONN_Success) {
          char str[64];
          int http_v1, http_v2, http_status = 0;
          Uint4 n_peek = BUF_Peek(buf, str, sizeof(str)-1);
          ASSERT( 2 <= n_peek  &&  n_peek < sizeof(str) );
          str[n_peek] = '\0';
          if (sscanf(str, " HTTP/%d.%d %d ", &http_v1, &http_v2, &http_status)
              != 3  ||  http_status < 200  ||  299 < http_status)
            server_error = TRUE;
      }

      /* skip HTTP header */
      if (status == eCONN_Success) {
        if ( uuu->info->debug_printout ) {
          char  data[256];
          Uint4 n_read;

          /* skip & printout the HTTP header */
          status = SOCK_StripToPattern(uuu->sock, "\r\n\r\n", 4, &buf, 0);
          fprintf(stderr, "\
\n\n----- [BEGIN] CON_URL HTTP Header(%ld bytes follow after \\n) -----\n",
                  (long)BUF_Size(buf));
          while ((n_read = BUF_Read(buf, data, sizeof(data))) > 0)
            fwrite(data, 1, (size_t)n_read, stderr);
          fprintf(stderr, "\
\n----- [END] CON_URL HTTP Header -----\n\n");
          fflush(stderr);

          /* skip & printout the content, if server error is detected */
          if ( server_error ) {
            fprintf(stderr, "\
\n\n----- [BEGIN] Detected a server error -----\n");
            for (;;) {
              ESOCK_ErrCode err_code =
                SOCK_Read(uuu->sock, (void*)data, (Uint4)sizeof(data),
                          &n_read);
              if (err_code != eSOCK_ESuccess)
                break;
              FileWrite((const void*)data, 1, (size_t)n_read, stderr);
            }
            fprintf(stderr, "\
\n----- [END] Detected a server error -----\n\n");
            fflush(stderr);
          }
        } else if ( !server_error ) {
          /* skip HTTP header */
          status = SOCK_StripToPattern(uuu->sock, "\r\n\r\n", 4, 0, 0);
        }
      }

      buf = BUF_Destroy(buf);
      if ( server_error )
        return eCONN_Unknown;
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
        if ( n_decoded ) {
            Nlm_Uint4 x_read;
            SOCK_Read(uuu->sock, peek_buf, n_decoded, &x_read);
            ASSERT( x_read == n_decoded );
            status = eCONN_Success;
        } else if ( SOCK_Eof(uuu->sock) ) {
            /* we are at EOF, and the remaining data cannot be decoded */
            status = eCONN_Unknown;
        } 
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
  if ( uuu->adjust_cleanup )
    uuu->adjust_cleanup(uuu->adjust_data);
  BUF_Destroy(uuu->buf);
  Nlm_MemFree(uuu);
  Nlm_MemFree(connector);
  return eCONN_Success;
}



/***********************************************************************
 *  EXTERNAL -- the connector's "constructor"
 ***********************************************************************/


NLM_EXTERN CONNECTOR URL_CreateConnectorEx
(const SNetConnInfo* info,
 const Nlm_Char*     user_header,
 URLC_Flags          flags,
 FAdjustInfo         adjust_info,
 void*               adjust_data,
 FAdjustCleanup      adjust_cleanup)
{
  CONNECTOR      ccc = (SConnector   *)Nlm_MemNew(sizeof(SConnector   ));
  SUrlConnector* uuu = (SUrlConnector*)Nlm_MemNew(sizeof(SUrlConnector));

  /* initialize internal data structures */
  uuu->info = info ? NetConnInfo_Clone(info) : NetConnInfo_Create(0, 0);
  NetConnInfo_AdjustForHttpProxy(uuu->info);
  if ( uuu->info->debug_printout )
    NetConnInfo_Print(uuu->info, stderr);

  uuu->user_header = user_header ? Nlm_StringSave(user_header) : 0;
  uuu->adjust_info    = adjust_info;
  uuu->adjust_data    = adjust_data;
  uuu->adjust_cleanup = adjust_cleanup;
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


NLM_EXTERN CONNECTOR URL_CreateConnector
(const SNetConnInfo* info,
 const Nlm_Char*     user_header,
 URLC_Flags          flags)
{
  return URL_CreateConnectorEx(info, user_header, flags, 0, 0, 0);
}



/***********************************************************************
 *  TEST 1
 ***********************************************************************/

#ifdef TEST_MODULE__CON_URL

#include <conntest.h>

#define TEST_ENGINE_HOST     "ray.nlm.nih.gov"
#define TEST_ENGINE_PORT     "6224"
#define TEST_ENGINE_PATH     "/cgi-bin/tools/vakatov/con_url.cgi"
#define TEST_ENGINE_ARGS     "arg1+arg2+arg3"

#define TEST_LOGFILE         "con_url.out"
#define TEST_DEBUG_PRINTOUT  "yes"


extern
#ifdef __cplusplus
"C"
#endif
Int2 Main(void)
{
  const Char*   user_header = 0;
  Uint4         flags;
  STimeout      timeout;
  CONNECTOR     connector;
  FILE*         log_file;

  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_HOST, TEST_ENGINE_HOST);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_PORT, TEST_ENGINE_PORT);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_PATH, TEST_ENGINE_PATH);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_ARGS, TEST_ENGINE_ARGS);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_DEBUG_PRINTOUT, TEST_DEBUG_PRINTOUT);

  timeout.sec  = 5;
  timeout.usec = 123456;
  log_file = FileOpen(TEST_LOGFILE, "wb");

  flags = URLC_HTTP_HEADER | URLC_URL_CODEC | URLC_URL_ENCODE_ARGS;
  connector = URL_CreateConnector(0, user_header, flags);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_SINGLE_BOUNCE_PRINT);

  flags = 0;
  connector = URL_CreateConnector(0, user_header, flags);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_SINGLE_BOUNCE_CHECK);

  flags = URLC_AUTO_RECONNECT;
  connector = URL_CreateConnector(0, user_header, flags);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_ALL);

  flags = URLC_AUTO_RECONNECT | URLC_URL_CODEC;
  connector = URL_CreateConnector(0, user_header, flags);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_ALL);

  FileClose(log_file);
  return 0;
}

#endif  /* TEST_MODULE__CON_URL */



/***********************************************************************
 *  TEST 2
 ***********************************************************************/


#ifdef TEST__HitURL

#ifdef __cplusplus
extern "C"
#endif
Int2 Main(void)
{
  /* Prepare to connect:  parse and check cmd.-line args, etc. */
  Int4   argc = GetArgc();
  Char **argv = GetArgv();

  const Char* host         = (argc > 1) ? argv[1] : "";
  const Char* port_str     = (argc > 2) ? argv[2] : "0";
  Uint2       port         = (Uint2) (atoi(port_str));
  const Char* path         = (argc > 3) ? argv[3] : "";
  const Char* args         = (argc > 4) ? argv[4] : "";
  const Char* inp_file     = (argc > 5) ? argv[5] : "";
  const Char* user_header  = (argc > 6) ? argv[6] : "";

  CONN conn;
  EConnStatus status;
  char buffer[10000];

  ErrSetLogfile("stderr", 0);
  ErrSetFatalLevel(SEV_FATAL);
  ErrSetMessageLevel(SEV_MIN);
  ErrSetOptFlags(EO_SHOW_FILELINE | EO_SHOW_ERRTEXT | EO_SHOW_MSGTEXT);

  fprintf(stderr, "Running...\n"
          "  Executable:      '%s'\n"
          "  URL host:        '%s'\n"
          "  URL port:         %hu\n"
          "  URL path:        '%s'\n"
          "  URL args:        '%s'\n"
          "  Input data file: '%s'\n"
          "  User header:     '%s'\n"
          " Reply(if any) from the hit URL goes to the standard output.\n\n",
          argv[0],
          host, (unsigned short)port, path, args, inp_file, user_header);

  if (argc < 4) {
    fprintf(stderr,
            "Usage:   %s host port path args inp_file [user_header]\n"
            "Example: %s ............\n",
            argv[0], argv[0]);
    ErrPostEx(SEV_ERROR, 0, 0, "Two few arguments.");
    return 1;
  }

  if (FileLengthEx(inp_file) < 0) {
    ErrPostEx(SEV_ERROR, 0, 0, "Non-existent file '%s'", inp_file);
    return 2;
  }

  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_HOST, host);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_PORT, port_str);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_PATH, path);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_ARGS, args);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_DEBUG_PRINTOUT, "yes");

  /* Connect */
  {{
    CONNECTOR connector = URL_CreateConnector(0, user_header, 0);
    ASSERT( connector );
    VERIFY( CONN_Create(connector, &conn) == eCONN_Success );
  }}

  {{ /* Pump data from the input file to URL */
    FILE* fp = FileOpen(inp_file, "rb");
    if ( !fp ) {
      ErrPostEx(SEV_ERROR, 0, 0, "Cannot open file '%s' for read", inp_file);
      return 4;
    }

    for (;;) {
      Uint4 n_written;
      size_t n_read = FileRead((void *)buffer, 1, sizeof(buffer), fp);
      if (n_read <= 0)
        break; /* EOF */


      status = CONN_Write(conn, buffer, (Uint4)n_read, &n_written);
      if (status != eCONN_Success) {
        ErrPostEx(SEV_ERROR, 0, 0, "Error writing to te URL(%s)",
                  CONN_StatusString(status));
        return 6;
      }
    }

    FileClose(fp);
  }}

  /* Read reply from connection, write it to standard output */
  {{
    Uint4 n_read;
    for (;;) {
      status = CONN_Read(conn, (void*)buffer, (Uint4)sizeof(buffer),
                         &n_read, eCR_Read);
      if (status != eCONN_Success)
        break;

      FileWrite((const void*)buffer, 1, (size_t)n_read, stdout);
    }

    if (status != eCONN_Closed) {
      ErrPostEx(SEV_WARNING, 0, 0,
                "Error reading from URL(%s)", CONN_StatusString(status));
    }
    fprintf(stdout, "\n");
  }}

  /* Success:  close the connection and exit */
  CONN_Close(conn);
  return 0;
}

#endif /* TEST__HitURL */
