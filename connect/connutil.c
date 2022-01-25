/*  $Id: connutil.c,v 6.15 2000/06/29 17:32:21 vakatov Exp $
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
*   Auxiliary API to get "standard" connection-related configuration
*   from the program config.-file(s) and environment
*
* --------------------------------------------------------------------------
* $Log: connutil.c,v $
* Revision 6.15  2000/06/29 17:32:21  vakatov
* Added more MIME sub-types
* Allow non-NCBI ("standard") MIME types
* Added MIME_*ContentTypeEx(), and tests for them
*
* Revision 6.14  1999/09/13 15:54:31  vakatov
* Added URL_DecodeEx() -- for "relaxed" URL decoding: let the user to
* allow some of the symbols prohibited by the standard
*
* Revision 6.13  1999/07/26 17:58:35  vakatov
* Use "\r\n" rather than just "\n" as the HTTP header line terminator.
* Made the hard-coded names and defautl values of config.parameters
* (#define CFG_CONN_, DEF_CONN_) be public.
* +eMIME_Fasta.
*
* Revision 6.12  1999/07/23 13:16:08  vakatov
* MIME_ParseContentType() - dont crash if passed NULL string
*
* Revision 6.11  1999/07/23 00:37:20  vakatov
* Final version of NCBI MIME functions
*
* Revision 6.10  1999/07/21 21:32:32  vakatov
* Get rid of the "nph-" prefixes for DISPD/NCBID -- use new DISPD/NCBID
*
* Revision 6.9  1999/07/16 21:19:10  vakatov
* + NetConnInfo_Print()
* + NCBI MIME-type handling routines and typedefs for a future use (MIME_*)
* NetConnInfo_Create():  typo fixed ("info->ncbid_*")
* Ncbi_ConnectURL():     use standard HTTP header terminator("\r\n\r\n")
*
* Revision 6.8  1999/07/09 22:53:37  vakatov
* Added more conf. parameters:  SRV_NCBID_PORT, SRV_NCBID_PATH, SRV_LB_DISABLE
*
* Revision 6.7  1999/06/28 16:28:02  vakatov
* SNetConnInfo:: separated the HTTP and the CERL-like(non-transparent)
* firewall proxies;  renamed config./env. parameters accordingly:
*   SRV_HTTP_PROXY_HOST/PORT replaced the former SRV_PROXY_HOST/PORT
*   SRV_PROXY_HOST now specifies the host name of non-transparent CERN proxy
*   SRV_PROXY_PORT is obsolete
*   SRV_PROXY_TRANSPARENT is obsolete
* Also:  NetConnInfo_AdjustForCernProxy --> ...HttpProxy
*
* Revision 6.6  1999/04/09 22:27:01  vakatov
* Ncbi_ConnectURL():  added "encode_args"
*
* Revision 6.5  1999/04/09 21:36:01  vakatov
* Split former CGI "args" into "path?args":
* - added "path" arg. to Ncbi_ConnectURL(), have it URL-encoded automagically;
* - added "path" field to SNetConnInfo
* - split ..._ENGINE_URL defs and config. names to ..._ENGINE_PATH / ARG
*
* Revision 6.4  1999/04/08 18:15:47  vakatov
* URL_Decode():  a more accurate handling of incomplete "%.." sequences
*
* Revision 6.3  1999/04/07 20:43:52  vakatov
* Added URL_Encode() and URL_Decode()
*
* Revision 6.2  1999/04/01 21:38:20  vakatov
* DEF_ENGINE_URL:  added leading '/'
* ==========================================================================
*/

#include <ncbi.h>
#include <connutil.h>


/***********************************************************************
 *  EXTERNAL
 ***********************************************************************/


NLM_EXTERN SNetConnInfo* NetConnInfo_Create
(const Char* conf_file,
 const Char* conf_section)
{
  SNetConnInfo* info = (SNetConnInfo*)MemNew(sizeof(SNetConnInfo));

  /* fallbacks for the conf. file & section names */
  if (!conf_file  ||  !*conf_file)
    conf_file = DEF_CONN_CONF_FILE;
  if (!conf_section  ||  !*conf_section)
    conf_section = DEF_CONN_CONF_SECTION;

  {{ /* client host */
    info->client_host = (Char*)MemNew(MAX_CONN_HOST_LEN);
    if ( !GetHostName(info->client_host, MAX_CONN_HOST_LEN) ) {
      ErrPostEx(SEV_WARNING, 0, 0,
                "[NetConnInfo_Create]  Cannot get local host name");
      info->client_host[0] = '\0';
    }
  }}

  {{ /* alternate server host name */
    info->host = (Char*)MemNew(MAX_CONN_HOST_LEN);
    GetEnvParam(conf_file, conf_section, CFG_CONN_ENGINE_HOST,
                info->host, MAX_CONN_HOST_LEN, DEF_CONN_ENGINE_HOST);
  }}

  {{ /* alternate service port */
    Char str[32];
    int      val;
    GetEnvParam(conf_file, conf_section, CFG_CONN_ENGINE_PORT,
                str, sizeof(str), "");
    val = atoi(str);
    info->port = (Uint2)(val > 0 ? val : DEF_CONN_ENGINE_PORT);
  }}

  {{ /* alternate service path */
    info->path = (Char*)MemNew(MAX_CONN_PATH_LEN);
    GetEnvParam(conf_file, conf_section, CFG_CONN_ENGINE_PATH,
                info->path, MAX_CONN_PATH_LEN, DEF_CONN_ENGINE_PATH);
  }}

  {{ /* alternate args */
    info->args = (Char*)MemNew(MAX_CONN_ARGS_LEN);
    GetEnvParam(conf_file, conf_section, CFG_CONN_ENGINE_ARGS,
                info->args, MAX_CONN_ARGS_LEN, DEF_CONN_ENGINE_ARGS);
  }}

  {{ /* alternate connection timeout */
    Char   str[32];
    double val;
    GetEnvParam(conf_file, conf_section, CFG_CONN_TIMEOUT,
                str, sizeof(str), "");
    val = atof(str);
    if (val <= 0)
      val = DEF_CONN_TIMEOUT;
    info->timeout.sec  = (Uint4)val;
    info->timeout.usec = (Uint4)((val - info->timeout.sec) * 1000000);
  }}

  {{ /* alternate the max. number of attempts to establish a connection */
    Char str[32];
    int  val;
    GetEnvParam(conf_file, conf_section, CFG_CONN_TRY,
                str, sizeof(str), "");
    val = atoi(str);
    info->conn_try = (Uint4)((val > 0) ? val : DEF_CONN_TRY);
  }}

  {{ /* HTTP proxy server? */
    Char http_proxy_host[MAX_CONN_HOST_LEN];
    GetEnvParam(conf_file, conf_section, CFG_CONN_HTTP_PROXY_HOST,
                http_proxy_host, sizeof(http_proxy_host),
                DEF_CONN_HTTP_PROXY_HOST);

    if ( *http_proxy_host ) {
      /* yes, use the specified HTTP proxy server */
      Char  str[32];
      int   val;
      GetEnvParam(conf_file, conf_section, CFG_CONN_HTTP_PROXY_PORT,
                  str, sizeof(str), "");
      val = atoi(str);

      info->http_proxy_port = (Uint2)(val>0 ? val : DEF_CONN_HTTP_PROXY_PORT);
      info->http_proxy_host = StringSave(http_proxy_host);
    }
  }}

  /* non-transparent CERN-like firewall proxy server? */
  {{
    Char proxy_host[MAX_CONN_HOST_LEN];
    GetEnvParam(conf_file, conf_section, CFG_CONN_PROXY_HOST,
                proxy_host, sizeof(proxy_host), DEF_CONN_PROXY_HOST);
    if ( *proxy_host )
      info->proxy_host = StringSave(proxy_host);
  }}

  {{ /* alternate the debug printout feature */
    Char str[32];
    GetEnvParam(conf_file, conf_section, CFG_CONN_DEBUG_PRINTOUT,
                str, sizeof(str), DEF_CONN_DEBUG_PRINTOUT);
    info->debug_printout = (Boolean)
      (*str  &&
       (StringICmp(str, "1"   ) == 0  ||
        StringICmp(str, "true") == 0  ||
        StringICmp(str, "yes" ) == 0));
  }}

  {{ /* alternate the firewall mode, if not set already */
    if ( !info->firewall ) {
      Char str[32];
      GetEnvParam(conf_file, conf_section, CFG_CONN_FIREWALL,
                  str, sizeof(str), DEF_CONN_FIREWALL);
      info->firewall = (Boolean)
        (*str  &&
         (StringICmp(str, "1"   ) == 0  ||
          StringICmp(str, "true") == 0  ||
          StringICmp(str, "yes" ) == 0));
    }
  }}

  {{ /* alternate NCBID port */
    Char str[32];
    int      val;
    GetEnvParam(conf_file, conf_section, CFG_CONN_NCBID_PORT,
                str, sizeof(str), "");
    val = atoi(str);
    info->ncbid_port = (Uint2)(val > 0 ? val : DEF_CONN_NCBID_PORT);
  }}

  {{ /* alternate NCBID path */
    info->ncbid_path = (Char*)MemNew(MAX_CONN_PATH_LEN);
    GetEnvParam(conf_file, conf_section, CFG_CONN_NCBID_PATH,
                info->ncbid_path, MAX_CONN_PATH_LEN, DEF_CONN_NCBID_PATH);
  }}

  {{ /* if to prohibit the use of local load balancer */
    Char str[32];
    GetEnvParam(conf_file, conf_section, CFG_CONN_LB_DISABLE,
                str, sizeof(str), DEF_CONN_LB_DISABLE);
    info->lb_disable = (Boolean)
        (*str  &&
         (StringICmp(str, "1"   ) == 0  ||
          StringICmp(str, "true") == 0  ||
          StringICmp(str, "yes" ) == 0));
  }}

  return info;
}  /* CONN_ComposeConnectInfo */


NLM_EXTERN Boolean NetConnInfo_AdjustForHttpProxy
(SNetConnInfo* info)
{
  Char* x_path;
  if (info->http_proxy_adjusted  ||  !info->http_proxy_host)
    return FALSE;

  x_path = (Char*)MemNew(16 + StrLen(info->host) + StrLen(info->path));
  sprintf(x_path, "http://%s:%hu/%s",
          info->host, (unsigned short)info->port, info->path);
  MemFree(info->host);
  MemFree(info->path);

  info->host = StringSave(info->http_proxy_host);
  info->port = info->http_proxy_port;
  info->path = x_path;
  info->http_proxy_adjusted = TRUE;
  return TRUE;
}


NLM_EXTERN SNetConnInfo* NetConnInfo_Clone
(const SNetConnInfo* info)
{
  SNetConnInfo* x_info;

  if ( !info )
    return 0;

  x_info = (SNetConnInfo*) MemNew(sizeof(SNetConnInfo));

  if ( info->client_host )
    x_info->client_host        = StringSave(info->client_host);
  if ( info->host )
    x_info->host               = StringSave(info->host);
  x_info->port                 = info->port;
  if ( info->path )
    x_info->path               = StringSave(info->path);
  if ( info->args )
    x_info->args               = StringSave(info->args);
  x_info->timeout              = info->timeout;
  x_info->conn_try             = info->conn_try;
  if ( info->http_proxy_host )
    x_info->http_proxy_host    = StringSave(info->http_proxy_host);
  x_info->http_proxy_port      = info->http_proxy_port;
  if ( info->proxy_host )
    x_info->proxy_host         = StringSave(info->proxy_host);
  x_info->debug_printout       = info->debug_printout;
  x_info->firewall             = info->firewall;
  x_info->ncbid_port           = info->ncbid_port;
  if ( info->ncbid_path )
    x_info->ncbid_path         = StringSave(info->ncbid_path);
  x_info->lb_disable           = info->lb_disable;
  x_info->http_proxy_adjusted  = info->http_proxy_adjusted;

  return x_info;
}


static void s_PrintString(FILE* fp, const char* name, const char* str) {
  if ( str )
    fprintf(fp, "%-16s: \"%s\"\n", name, str);
  else
    fprintf(fp, "%-16s: <NULL>\n", name);
}
static void s_PrintULong(FILE* fp, const char* name, unsigned long lll) {
  fprintf(fp, "%-16s: %lu\n", name, lll);
}
static void s_PrintBool(FILE* fp, const char* name, Nlm_Boolean bbb) {
  fprintf(fp, "%-16s: %s\n", name, bbb ? "TRUE" : "FALSE");
}

NLM_EXTERN void NetConnInfo_Print
(const SNetConnInfo* info,
 FILE*               fp)
{
  if ( !fp )
    return;

  fprintf(fp, "\n\n----- [BEGIN] NetConnInfo_Print -----\n");

  if ( info ) {
    s_PrintString(fp, "client_host",     info->client_host);
    s_PrintString(fp, "host",            info->host);
    s_PrintULong (fp, "port",            info->port);
    s_PrintString(fp, "path",            info->path);
    s_PrintString(fp, "args",            info->args);
    s_PrintULong (fp, "timeout(sec)",    info->timeout.sec);
    s_PrintULong (fp, "timeout(usec)",   info->timeout.usec);
    s_PrintULong (fp, "conn_try",        info->conn_try);
    s_PrintString(fp, "http_proxy_host", info->http_proxy_host);
    s_PrintULong (fp, "http_proxy_port", info->http_proxy_port);
    s_PrintString(fp, "proxy_host",      info->proxy_host);
    s_PrintBool  (fp, "debug_printout",  info->debug_printout);
    s_PrintBool  (fp, "firewall",        info->firewall);
    s_PrintULong (fp, "ncbid_port",      info->ncbid_port);
    s_PrintString(fp, "ncbid_path",      info->ncbid_path);
    s_PrintBool  (fp, "lb_disable",      info->lb_disable);
    s_PrintBool  (fp, "proxy_adjusted",  info->http_proxy_adjusted);
  } else {
    fprintf(fp, "<NULL>\n");
  }

  fprintf(fp, "----- [END] NetConnInfo_Print -----\n\n");
}


NLM_EXTERN void NetConnInfo_Destroy
(SNetConnInfo** info)
{
  if (!info  ||  !*info)
    return;

  MemFree((*info)->client_host);
  MemFree((*info)->host);
  MemFree((*info)->path);
  MemFree((*info)->args);
  MemFree((*info)->http_proxy_host);
  MemFree((*info)->proxy_host);
  MemFree((*info)->ncbid_path);
  MemFree(*info);
  *info = 0;
}


NLM_EXTERN EConnStatus ESOCK2ECONN
(ESOCK_ErrCode err_code)
{
  switch ( err_code ) {
  case eSOCK_ESuccess:
    return eCONN_Success;
  case eSOCK_ETimeout:
    return eCONN_Timeout;
  case eSOCK_EClosed:
    return eCONN_Closed;
  case eSOCK_EUnknown:
    return eCONN_Unknown;
  }
  ASSERT(0);
  return eCONN_Unknown;
}


NLM_EXTERN SOCK Ncbi_ConnectURL
(const Char*     host,
 Uint2           port,
 const Char*     path,
 const Char*     args,
 size_t          content_length,
 const STimeout* c_timeout,
 const STimeout* rw_timeout,
 const Char*     user_header,
 Boolean         encode_args
 )
{
  static const Char X_POST_1[] = "POST ";
  static const Char X_POST_Q[] = "?";
  static const Char X_POST_E[] = " HTTP/1.0\r\n";

  SOCK  sock;
  Char  buffer[128];
  Char* x_args = 0;

  /* check the args */
  if (!host  ||  !*host  ||  !port  ||  !path  ||  !*path  ||
      (user_header  &&  *user_header  &&
       user_header[StrLen(user_header)-1] != '\n')) {
    ErrPostEx(SEV_ERROR, 0, 0,
              "[Ncbi_ConnectURL]  Bad arguments");
    return 0;
  }

  /* connect to HTTPD */
  if (SOCK_Create(host, port, c_timeout, &sock) != eSOCK_ESuccess) {
    ErrPostEx(SEV_ERROR, 0, 0,
              "[Ncbi_ConnectURL]  Cannot connect to host \"%s\", port %d;",
              host, (int)port);
    return 0;
  }

  /* setup i/o timeout for the connection */
  if (SOCK_SetTimeout(sock, eSOCK_OnReadWrite, rw_timeout, 0, 0)
      != eSOCK_ESuccess) {
    ErrPostEx(SEV_ERROR, 0, 0,
              "[Ncbi_ConnectURL]  Cannot setup timeout for the connection"
              " handshake with host \"%s\", port %d",
              host, (int)port);
    SOCK_Close(sock);
    return 0;
  }

  /* URL-encode "args", if any specified */
  if (args  &&  *args) {
    if ( encode_args ) {
      Uint4 src_size = StrLen(args);
      Uint4 dst_size = 3 * src_size;
      Uint4 src_read, dst_written;
      x_args = (Char*)MemNew(dst_size + 1);
      URL_Encode(args, src_size, &src_read, x_args, dst_size, &dst_written);
      x_args[dst_written] = '\0';
      ASSERT( src_read == src_size );
    } else {
      x_args = StringSave(args);
    }
  }

  /* compose and send HTTP header */
  if (/*  POST <path>?<args> HTTP/1.0\r\n */
      SOCK_Write(sock, (const void*)X_POST_1,  StrLen(X_POST_1 ), 0)
      != eSOCK_ESuccess  ||
      SOCK_Write(sock, (const void*)path, StrLen(path), 0)
      != eSOCK_ESuccess  ||
      (x_args  &&
       (SOCK_Write(sock, (const void*)X_POST_Q, StrLen(X_POST_Q), 0)
        != eSOCK_ESuccess  ||
        SOCK_Write(sock, (const void*)x_args, StrLen(x_args), 0)
        != eSOCK_ESuccess
        )
       )  ||
      SOCK_Write(sock, (const void*)X_POST_E,  StrLen(X_POST_E ), 0)
      != eSOCK_ESuccess  ||

      /*  <user_header> */
      (user_header  &&
       SOCK_Write(sock, (const void*)user_header, StrLen(user_header), 0)
       != eSOCK_ESuccess)  ||

      /*  Content-Length: <content_length>\r\n\r\n */
      sprintf(buffer, "Content-Length: %lu\r\n\r\n",
              (unsigned long)content_length) <= 0  ||
      SOCK_Write(sock, (const void*)buffer, StrLen(buffer), 0)
      != eSOCK_ESuccess)
    {
      /* error */
      ErrPostEx(SEV_ERROR, 0, 0,
                "[Ncbi_ConnectURL]  Error sending HTTP header");
      MemFree(x_args);
      SOCK_Close(sock);
      return 0;
    }

  /* success */
  MemFree(x_args);
  return sock;
}


/* Code for the "*_StripToPattern()" functions
 */
typedef  EConnStatus (*FPeekRead)
     (void*   source,
      void*   buffer,
      Uint4   size,
      Uint4*  n_read,
      Boolean do_peek
      );

static EConnStatus s_StripToPattern
(void*       source,
 FPeekRead   PeekRead,
 const void* pattern,
 Uint4       pattern_size,
 BUF*        buf,
 Uint4*      n_discarded)
{
  EConnStatus status;
  char*       buffer;
  Uint4       buffer_size;
  Uint4       n_read = 0;

  /* check args */
  if ( n_discarded )
    *n_discarded = 0;
  if (!source  ||  !pattern  ||  !pattern_size)
    return eCONN_InvalidArg;

  /* allocate a temporary read buffer */
  buffer_size = 2 * pattern_size;
  if (buffer_size < 4096)
    buffer_size = 4096;
  buffer = (char*)MemNew((size_t)buffer_size);

  /* peek/read;  search for the pattern;  maybe, store the discarded data */
  for (;;) {
    /* peek */
    Uint4 n_peeked, n_stored, x_discarded;
    ASSERT( n_read < pattern_size );
    status = PeekRead(source, buffer + n_read, buffer_size - n_read, &n_peeked,
                      TRUE);
    if ( !n_peeked ) {
      ASSERT( status != eCONN_Success );
      break; /* error */
    }

    n_stored = n_read + n_peeked;

    if (n_stored >= pattern_size) {
      /* search for the pattern */
      Uint4 n_check = n_stored - pattern_size + 1;
      const char* b;
      for (b = buffer;  n_check;  b++, n_check--) {
        if (*b != *((char*)pattern))
          continue;
        if (MemCmp(b, pattern, (size_t)pattern_size) == 0)
          break; /* found */
      }
      /* pattern found */
      if ( n_check ) {
        Uint4 x_read =  b - buffer + pattern_size;
        ASSERT( MemCmp(b, pattern, pattern_size) == 0 );
        status = PeekRead(source, buffer + n_read, x_read - n_read,
                          &x_discarded, FALSE);
        ASSERT( status == eCONN_Success );
        ASSERT( x_discarded == x_read - n_read );
        if ( buf )
          BUF_Write(buf, buffer + n_read, x_read - n_read);
        if ( n_discarded )
          *n_discarded += x_read - n_read;
        break; /* success */
      }
    }

    /* pattern not found yet */
    status = PeekRead(source, buffer + n_read, n_peeked, &x_discarded, FALSE);
    ASSERT( status == eCONN_Success );
    ASSERT( x_discarded == n_peeked );
    if ( buf )
      BUF_Write(buf, buffer + n_read, n_peeked);
    if ( n_discarded )
      *n_discarded += n_peeked;
    n_read = n_stored;

    if (n_read > pattern_size) {
      Uint4 n_cut = n_read - pattern_size + 1;
      n_read = pattern_size - 1;
      MemMove(buffer, buffer + n_cut, (size_t)n_read);
    }
  }

  /* cleanup & exit */
  MemFree(buffer);
  return status;
}

static EConnStatus s_CONN_Read
(void*   source,
 void*   buffer,
 Uint4   size,
 Uint4*  n_read,
 Boolean do_peek)
{
  return CONN_Read((CONN)source, buffer, size, n_read,
                   do_peek ? eCR_Peek : eCR_Read);
}

NLM_EXTERN EConnStatus CONN_StripToPattern
(CONN        conn,
 const void* pattern,
 Uint4       pattern_size,
 BUF*        buf,
 Uint4*      n_discarded)
{
  return s_StripToPattern(conn, s_CONN_Read, pattern, pattern_size, buf,
                          n_discarded);
}


static EConnStatus s_SOCK_Read
(void*   source,
 void*   buffer,
 Uint4   size,
 Uint4*  n_read,
 Boolean do_peek)
{
  return ESOCK2ECONN(do_peek ?
                     SOCK_Peek((SOCK)source, buffer, size, n_read) :
                     SOCK_Read((SOCK)source, buffer, size, n_read));
}

NLM_EXTERN EConnStatus SOCK_StripToPattern
(SOCK        sock,
 const void* pattern,
 Uint4       pattern_size,
 BUF*        buf,
 Uint4*      n_discarded)
{
  return s_StripToPattern(sock, s_SOCK_Read, pattern, pattern_size, buf,
                          n_discarded);
}


/* Return integer (0..15) corresponding to the "ch" as a hex digit
 * Return -1 on error
 */
static int s_HexChar(char ch)
{
  if ('0' <= ch  &&  ch <= '9')
    return ch - '0';
  if ('a' <= ch  &&  ch <= 'f')
    return 10 + (ch - 'a');
  if ('A' <= ch  &&  ch <= 'F')
    return 10 + (ch - 'A');
  return -1;
}

/* The URL-encoding table
 */
static const char s_Encode[256][4] = {
  "%00", "%01", "%02", "%03", "%04", "%05", "%06", "%07",
  "%08", "%09", "%0A", "%0B", "%0C", "%0D", "%0E", "%0F",
  "%10", "%11", "%12", "%13", "%14", "%15", "%16", "%17",
  "%18", "%19", "%1A", "%1B", "%1C", "%1D", "%1E", "%1F",
  "+",   "!",   "%22", "%23", "$",   "%25", "%26", "'",
  "(",   ")",   "*",   "%2B", ",",   "-",   ".",   "%2F",
  "0",   "1",   "2",   "3",   "4",   "5",   "6",   "7",
  "8",   "9",   "%3A", "%3B", "%3C", "%3D", "%3E", "%3F",
  "%40", "A",   "B",   "C",   "D",   "E",   "F",   "G",
  "H",   "I",   "J",   "K",   "L",   "M",   "N",   "O",
  "P",   "Q",   "R",   "S",   "T",   "U",   "V",   "W",
  "X",   "Y",   "Z",   "%5B", "%5C", "%5D", "%5E", "_",
  "%60", "a",   "b",   "c",   "d",   "e",   "f",   "g",
  "h",   "i",   "j",   "k",   "l",   "m",   "n",   "o",
  "p",   "q",   "r",   "s",   "t",   "u",   "v",   "w",
  "x",   "y",   "z",   "%7B", "%7C", "%7D", "%7E", "%7F",
  "%80", "%81", "%82", "%83", "%84", "%85", "%86", "%87",
  "%88", "%89", "%8A", "%8B", "%8C", "%8D", "%8E", "%8F",
  "%90", "%91", "%92", "%93", "%94", "%95", "%96", "%97",
  "%98", "%99", "%9A", "%9B", "%9C", "%9D", "%9E", "%9F",
  "%A0", "%A1", "%A2", "%A3", "%A4", "%A5", "%A6", "%A7",
  "%A8", "%A9", "%AA", "%AB", "%AC", "%AD", "%AE", "%AF",
  "%B0", "%B1", "%B2", "%B3", "%B4", "%B5", "%B6", "%B7",
  "%B8", "%B9", "%BA", "%BB", "%BC", "%BD", "%BE", "%BF",
  "%C0", "%C1", "%C2", "%C3", "%C4", "%C5", "%C6", "%C7",
  "%C8", "%C9", "%CA", "%CB", "%CC", "%CD", "%CE", "%CF",
  "%D0", "%D1", "%D2", "%D3", "%D4", "%D5", "%D6", "%D7",
  "%D8", "%D9", "%DA", "%DB", "%DC", "%DD", "%DE", "%DF",
  "%E0", "%E1", "%E2", "%E3", "%E4", "%E5", "%E6", "%E7",
  "%E8", "%E9", "%EA", "%EB", "%EC", "%ED", "%EE", "%EF",
  "%F0", "%F1", "%F2", "%F3", "%F4", "%F5", "%F6", "%F7",
  "%F8", "%F9", "%FA", "%FB", "%FC", "%FD", "%FE", "%FF"
};
#define VALID_URL_SYMBOL(ch)  (s_Encode[(unsigned char)ch][0] != '%')


NLM_EXTERN Boolean URL_DecodeEx
(const void* src_buf,
 Uint4   src_size,
 Uint4*  src_read,
 void*   dst_buf,
 Uint4   dst_size,
 Uint4*  dst_written,
 Char*   allow_symbols)
{
  unsigned char *src = (unsigned char*)src_buf;
  unsigned char *dst = (unsigned char*)dst_buf;

  *src_read    = 0;
  *dst_written = 0;
  if (!src_size  ||  !dst_size)
    return TRUE;

  for ( ;  *src_read != src_size  &&  *dst_written != dst_size;
        (*src_read)++, (*dst_written)++, src++, dst++) {
    switch ( *src ) {
    case '%': {
      int i1, i2;
      if (*src_read + 2 > src_size)
        return TRUE;
      if ((i1 = s_HexChar(*(++src))) == -1)
        return (Boolean)(*dst_written ? TRUE : FALSE);
      if (*src_read + 3 > src_size)
        return TRUE;
      if ((i2 = s_HexChar(*(++src))) == -1)
        return (Boolean)(*dst_written ? TRUE : FALSE);

      *dst = (unsigned char)((i1 << 4) + i2);
      *src_read += 2;
      break;
    }

    case '+': {
      *dst = ' ';
      break;
    }

    default:
      if (VALID_URL_SYMBOL(*src)  ||
          (allow_symbols  &&  strchr(allow_symbols, *src)))
        *dst = *src;
      else
        return (Boolean)(*dst_written ? TRUE : FALSE);
    }
  }

  ASSERT( src == (unsigned char*)src_buf + *src_read    );
  ASSERT( dst == (unsigned char*)dst_buf + *dst_written );
  return TRUE;
}


NLM_EXTERN Boolean URL_Decode
(const void* src_buf,
 Uint4   src_size,
 Uint4*  src_read,
 void*   dst_buf,
 Uint4   dst_size,
 Uint4*  dst_written)
{
  return URL_DecodeEx
    (src_buf, src_size, src_read, dst_buf, dst_size, dst_written, 0);
}


NLM_EXTERN void URL_Encode
(const void* src_buf,
 Uint4   src_size,
 Uint4*  src_read,
 void*   dst_buf,
 Uint4   dst_size,
 Uint4*  dst_written)
{
  unsigned char *src = (unsigned char*)src_buf;
  unsigned char *dst = (unsigned char*)dst_buf;

  *src_read    = 0;
  *dst_written = 0;
  if (!src_size  ||  !dst_size)
    return;

  for ( ;  *src_read != src_size  &&  *dst_written != dst_size;
        (*src_read)++, (*dst_written)++, src++, dst++) {
    const char* subst = s_Encode[*src];
    if (*subst != '%') {
      *dst = *subst;
    } else if (*dst_written < dst_size - 2) {
      *dst = '%';
      *(++dst) = *(++subst);
      *(++dst) = *(++subst);
      *dst_written += 2;
    }
    else {
      return;
    }
  }
  ASSERT( src == (unsigned char*)src_buf + *src_read    );
  ASSERT( dst == (unsigned char*)dst_buf + *dst_written );
}



/****************************************************************************
 * NCBI-specific MIME content type and sub-types
 */

static const char* s_MIME_Type[eMIME_T_Unknown+1] = {
  "x-ncbi-data",
  "text",
  "application",
  "unknown"
};

static const char* s_MIME_SubType[eMIME_Unknown+1] = {
  "x-asn-text",
  "x-asn-binary",
  "x-fasta",
  "x-www-form",
  "html",
  "x-unknown"
};

static const char* s_MIME_Encoding[eENCOD_None+1] = {
  "urlencoded",
  ""
};


NLM_EXTERN char* MIME_ComposeContentTypeEx
(EMIME_Type     type,
 EMIME_SubType  subtype,
 EMIME_Encoding encoding,
 char*          buf,
 size_t         buflen)
{
  static const char s_ContentType[] = "Content-Type: ";
  const char*       x_Type          = s_MIME_Type    [(int) type];
  const char*       x_SubType       = s_MIME_SubType [(int) subtype];
  const char*       x_Encoding      = s_MIME_Encoding[(int) encoding];
  char              x_buf[MAX_CONTENT_TYPE_LEN];

  if ( *x_Encoding ) {
    ASSERT(sizeof(s_ContentType) + strlen(x_Type) + strlen(x_SubType) + strlen(x_Encoding) + 4 < MAX_CONTENT_TYPE_LEN);
    sprintf(x_buf, "%s%s/%s-%s\r\n",
            s_ContentType, x_Type, x_SubType, x_Encoding);
  } else {
      ASSERT(sizeof(s_ContentType) + strlen(x_Type) + strlen(x_SubType) + 3 < MAX_CONTENT_TYPE_LEN);
    sprintf(x_buf, "%s%s/%s\r\n", s_ContentType, x_Type, x_SubType);
  }

  ASSERT( strlen(x_buf) < sizeof(x_buf) );
  ASSERT( strlen(x_buf) < buflen );
  StringNCpy_0(buf, x_buf, buflen);
  return buf;
}


NLM_EXTERN char* MIME_ComposeContentType
(EMIME_SubType  subtype,
 EMIME_Encoding encoding,
 char*          buf,
 size_t         buflen)
{
    return MIME_ComposeContentTypeEx(eMIME_T_NcbiData,
                                     subtype, encoding, buf, buflen);
}


NLM_EXTERN Boolean MIME_ParseContentTypeEx
(const char*     str,
 EMIME_Type*     type,
 EMIME_SubType*  subtype,
 EMIME_Encoding* encoding)
{
  char*   x_buf;
  char*   x_type;
  char*   x_subtype;
  int     i;

  if ( type )
    *type = eMIME_T_Unknown;
  if ( subtype )
    *subtype = eMIME_Unknown;
  if ( encoding )
    *encoding = eENCOD_None;

  if (!str  ||  !*str)
    return FALSE;

  {{
    size_t  x_size = strlen(str) + 1;
    x_buf     = (char*) malloc(2 * x_size);
    x_type    = x_buf  + x_size;
  }}

  strcpy(x_buf, str);
  StrLower(x_buf);
  if ((sscanf(x_buf, " content-type: %s ", x_type) != 1  &&
       sscanf(x_buf, " %s ", x_type) != 1)  ||
      (x_subtype = strchr(x_type, '/')) == 0) {
    free(x_buf);
    return FALSE;
  }
  *x_subtype++ = '\0';

  if ( type ) {
    for (i = 0;  i < (int) eMIME_T_Unknown;  i++) {
      if ( !strncmp(x_type, s_MIME_Type[i], strlen(s_MIME_Type[i])) )
        *type = (EMIME_Type) i;
    }
  }

  if ( subtype ) {
    for (i = 0;  i < (int) eMIME_Unknown;  i++) {
      if ( !strncmp(x_subtype, s_MIME_SubType[i], strlen(s_MIME_SubType[i])) )
        *subtype = (EMIME_SubType) i;
    }
  }

  if ( encoding ) {
    for (i = 0;  i < (int)eENCOD_None;  i++) {
      if (strstr(x_subtype, s_MIME_Encoding[i]) != 0)
        *encoding = (EMIME_Encoding)i;
    }
  }
  
  free(x_buf);
  return TRUE;
}


NLM_EXTERN Nlm_Boolean MIME_ParseContentType
(const char*     str,
 EMIME_SubType*  subtype,
 EMIME_Encoding* encoding)
{
    EMIME_Type type;
    if ( !MIME_ParseContentTypeEx(str, &type, subtype, encoding) )
        return FALSE;

    if (type != eMIME_T_NcbiData) {
        if ( subtype )
            *subtype  = eMIME_Unknown;
        if ( encoding )
            *encoding = eENCOD_None;
        return FALSE;
    }

    return TRUE;
}


/***********************************************************************
 *  TEST 1:  Ncbi_ConnectURL()
 ***********************************************************************/


#ifdef TEST__Ncbi_ConnectURL

#ifdef __cplusplus
extern "C"
#endif
Int2 Main(void)
{
  /* Prepare to connect:  parse and check cmd.-line args, etc. */
  Int4   argc = GetArgc();
  Char **argv = GetArgv();

  const Char* host         = (argc > 1) ? argv[1] : "";
  Uint2       port         = (argc > 2) ? (Uint2)(atoi(argv[2])) : 0;
  const Char* path         = (argc > 3) ? argv[3] : "";
  const Char* args         = (argc > 4) ? argv[4] : "";
  const Char* inp_file     = (argc > 5) ? argv[5] : "";
  const Char* user_header  = (argc > 6) ? argv[6] : "";

  size_t   content_length;
  STimeout timeout;

  SOCK sock;
  ESOCK_ErrCode err_code;
  char buffer[10000];

  ErrSetLogfile("stderr", 0);
  ErrSetFatalLevel(SEV_ERROR);
  ErrSetMessageLevel(SEV_MIN);
  ErrSetOptFlags(EO_SHOW_FILELINE | EO_SHOW_ERRTEXT | EO_SHOW_MSGTEXT);

  fprintf(stderr, "Running...\n"
          "  Executable:      '%s'\n"
          "  URL host:        '%s'\n"
          "  URL port:        %hu\n"
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

  {{
    Int4 file_len = FileLengthEx(inp_file);
    if (file_len < 0) {
      ErrPostEx(SEV_ERROR, 0, 0, "Non-existent file '%s'", inp_file);
      return 2;
    }
    content_length = (size_t)file_len;
  }}

  timeout.sec  = 10;
  timeout.usec = 0;

  /* Connect */
  sock = Ncbi_ConnectURL(host, port, path, args, content_length,
                         &timeout, &timeout, user_header, TRUE);
  if ( !sock )
    return 3;

  {{ /* Pump data from the input file to socket */
    FILE* fp = FileOpen(inp_file, "rb");
    if ( !fp ) {
      ErrPostEx(SEV_ERROR, 0, 0, "Cannot open file '%s' for read", inp_file);
      return 4;
    }

    for (;;) {
      Uint4 n_written;
      size_t n_read = FileRead((void *)buffer, 1, sizeof(buffer), fp);
      if (n_read <= 0) {
        if ( content_length ) {
          ErrPostEx(SEV_ERROR, 0, 0,
                    "Cannot read last %ld bytes from file '%s'",
                    (long)content_length, inp_file);
          return 5;
        }
        break;
      }

      content_length -= n_read;
      err_code = SOCK_Write(sock, buffer, (Uint4)n_read, &n_written);
      if (err_code != eSOCK_ESuccess) {
        ErrPostEx(SEV_ERROR, 0, 0, "Error writing to socket(%s)",
                  SOCK_ErrCodeStr(err_code));
        return 6;
      }
    }

    FileClose(fp);
  }}

  /* Read reply from socket, write it to standard output */
  {{
    Uint4 n_read;
    for (;;) {
      err_code = SOCK_Read(sock, (void*)buffer, (Uint4)sizeof(buffer),
                           &n_read);
      if (err_code != eSOCK_ESuccess)
        break;

      FileWrite((const void*)buffer, 1, (size_t)n_read, stdout);
    }

    if (err_code != eSOCK_EClosed) {
      ErrPostEx(SEV_WARNING, 0, 0,
                "Error occured after reading %ld bytes from socket(%s)",
                (long)content_length, SOCK_ErrCodeStr(err_code));
    }
    fprintf(stdout, "\n");
  }}

  /* Success:  close the socket and exit */
  SOCK_Close(sock);
  return 0;
}

#endif /* TEST__Ncbi_ConnectURL */



/***********************************************************************
 *  TEST 2:  URL_Encode(), URL_Decode()
 ***********************************************************************/

#ifdef TEST__URL_Encode

#ifdef __cplusplus
extern "C"
#endif
Int2 Main(void)
{
  typedef struct {
    char*       src_buf;
    Uint4   src_size;
    Uint4   src_read;
    char*       dst_buf;
    Uint4   dst_size;
    Uint4   dst_written;
    Boolean ok;
  } STestArg;

  static const STestArg s_TestEncode[] = {
    { "",           0, 0,  "",    0, 0, TRUE },
    { "abc",        3, 3,  "abc", 3, 3, TRUE },
    { "_ _%_;_\n_", 7, 0,  "",    0, 0, TRUE },
    { "_ _%_;_\n_", 0, 0,  "",    0, 0, TRUE },
    { "_ _%_;_\n_", 0, 0,  "",    7, 0, TRUE },
    { "_ _%_;_\n_:_\\_\"_", 15, 15,
      "_+_%25_%3B_%0A_%3A_%5C_%22_", 27, 27, TRUE },
    { "_ _%_;_\n_:_\\_\"_", 15, 13,
      "_+_%25_%3B_%0A_%3A_%5C_%22_", 25, 23, TRUE },
    { "_%_", 3, 1,  "_%25_", 2, 1, TRUE },
    { "_ _%_;_\n_", 7, 7,
      "_+_%25_%3B_%0A", 100, 11, TRUE }
  };

  static const STestArg s_TestDecode[] = {
    { "",    0, 0,   "", 0, 0,  TRUE },
    { "%25", 1, 0,   "", 0, 0,  TRUE },
    { "%25", 2, 0,   "", 0, 0,  TRUE },
    { "%25", 3, 3,  "%", 1, 1,  TRUE },
    { "%25", 3, 0,  "%", 0, 0,  TRUE },
    { "%%%", 2, 0,   "", 1, 0,  FALSE },
    { "%%%", 3, 0,   "", 1, 0,  FALSE },
    { "%xy", 3, 0,   "", 1, 0,  FALSE },
    { "\n",  1, 0,   "", 1, 0,  FALSE },
    { "a\t", 2, 1,  "a", 1, 1,  TRUE },
    { "#\n", 1, 0,   "", 0, 0,  TRUE },
    { "%a-", 3, 0,   "", 1, 0,  FALSE },
    { "%a-", 3, 0,   "", 0, 0,  TRUE },
    { "_+_%25_%3B_%0A_%3A_%5C_%22_", 27, 27,
      "_ _%_;_\n_:_\\_\"_", 15, 15, TRUE },
    { "_+_%25_%3B_%0A_%3A_%5C_%22_", 25, 23,
      "_ _%_;_\n_:_\\_\"_", 13, 13, TRUE },
    { "_+_%25_%3B_%0A_%3A_%5C_%22_", 27, 23,
      "_ _%_;_\n_:_\\_\"_", 13, 13, TRUE }
  };

  static const STestArg s_TestDecodeEx[] = {
    { "",    0, 0,    "", 0, 0,  TRUE },
    { "%25", 3, 0,   "%", 0, 0,  TRUE },
    { "%%%", 2, 0,    "", 1, 0,  FALSE },
    { "%xy", 3, 0,    "", 1, 0,  FALSE },
    { "\n",  1, 0,    "", 1, 0,  FALSE },
    { ">>a", 3, 3, ">>a", 3, 3,  TRUE },
    { ">b[", 3, 3, ">b[", 4, 3,  TRUE },
    { ">b]", 3, 2, ">b",  3, 2,  TRUE },
    { "[b]", 3, 2, "[b",  3, 2,  TRUE },
    { "<b>", 3, 0,   "",  3, 0,  FALSE },
    { "<e>", 3, 0,   "",  5, 0,  FALSE }
  };

  size_t i;
  Uint4 src_read, dst_written;
  char dst[1024];

  for (i = 0;  i < DIM(s_TestEncode);  i++) {
    const STestArg* arg = &s_TestEncode[i];
    URL_Encode(arg->src_buf, arg->src_size, &src_read,
               dst, arg->dst_size, &dst_written);
    ASSERT( src_read == arg->src_read );
    ASSERT( dst_written == arg->dst_written );
    ASSERT( !dst_written  ||  !memcmp(dst, arg->dst_buf, dst_written) );
  }

  for (i = 0;  i < DIM(s_TestDecode);  i++) {
    const STestArg* arg = &s_TestDecode[i];
    Boolean ok = URL_Decode(arg->src_buf, arg->src_size, &src_read,
                            dst, arg->dst_size, &dst_written);
    ASSERT( ok == arg->ok );
    ASSERT( src_read == arg->src_read );
    ASSERT( dst_written == arg->dst_written );
    ASSERT( !dst_written  ||  !memcmp(dst, arg->dst_buf, dst_written) );
  }

  for (i = 0;  i < DIM(s_TestDecodeEx);  i++) {
    const STestArg* arg = &s_TestDecodeEx[i];
    Boolean ok = URL_DecodeEx(arg->src_buf, arg->src_size, &src_read,
                              dst, arg->dst_size, &dst_written, "[>");
    ASSERT( ok == arg->ok );
    ASSERT( src_read == arg->src_read );
    ASSERT( dst_written == arg->dst_written );
    ASSERT( !dst_written  ||  !memcmp(dst, arg->dst_buf, dst_written) );
  }

  return 0;
}

#endif /* TEST__URL_Encode */



/***********************************************************************
 *  TEST 3:  Miscellaneous
 ***********************************************************************/

#ifdef TEST_MODULE__CONNUTIL

static Boolean s_CheckMIME
(const char* str,
 EMIME_Type type, EMIME_SubType subtype, EMIME_Encoding encoding)
{
    EMIME_Type     x_type;
    EMIME_SubType  x_subtype;
    EMIME_Encoding x_encoding;

    if (type == eMIME_T_NcbiData) {
        if (!MIME_ParseContentType(str, &x_subtype, 0)  ||
            x_subtype  != subtype) {
            return FALSE;
        }
        if (!MIME_ParseContentType(str, 0, &x_encoding)  ||
            x_encoding != encoding) {
            return FALSE;
        }
        if (!MIME_ParseContentType(str, &x_subtype, &x_encoding)  ||
            x_subtype != subtype  ||  x_encoding != encoding) {
            return FALSE;
        }
    }

    if (!MIME_ParseContentTypeEx(str, &x_type, &x_subtype, 0)  ||
        x_type != type  ||  x_subtype != subtype) {
        return FALSE;
    }
    if (!MIME_ParseContentTypeEx(str, &x_type, 0, &x_encoding)  ||
        x_type != type  ||  x_encoding != encoding) {
        return FALSE;
    }
    if (!MIME_ParseContentTypeEx(str, &x_type, &x_subtype, &x_encoding)  ||
        x_type != type  ||  x_subtype != subtype  ||  x_encoding != encoding) {
        return FALSE;
    }
    str = strchr(str, ':');
    if ( str ) {
        str++;
        return s_CheckMIME(str, type, subtype, encoding);
    }

    return TRUE;
}


#ifdef __cplusplus
extern "C"
#endif
Int2 Main(void)
{
  int i,j,k;

  /* MIME API */
  EMIME_Type     type;
  EMIME_SubType  subtype;
  EMIME_Encoding encoding;
  char str[MAX_CONTENT_TYPE_LEN];
  *str = '\0';
  for (k = 0, type = (EMIME_Type) k;
       k <= (int) eMIME_T_Unknown;  k++, type = (EMIME_Type) k) {
    for (i = 0, subtype = (EMIME_SubType) i;
         i <= (int)eMIME_Unknown;  i++, subtype = (EMIME_SubType) i) {
      for (j = 0, encoding = (EMIME_Encoding) j; 
           j <= (int)eENCOD_None;  j++, encoding = (EMIME_Encoding) j) {
        ASSERT( !s_CheckMIME(str, type, subtype, encoding) );
        MIME_ComposeContentTypeEx(type, subtype, encoding, str, sizeof(str));
        ASSERT( s_CheckMIME(str, type, subtype, encoding) );
      }
    }
  }

  ASSERT( s_CheckMIME("content-type:  x-ncbi-data/x-asn-binary ",
                      eMIME_T_NcbiData, eMIME_AsnBinary, eENCOD_None) );
  ASSERT( s_CheckMIME("content-type:  application/x-www-form-urlencoded ",
                      eMIME_T_Application, eMIME_WwwForm, eENCOD_Url) );
  ASSERT( s_CheckMIME("content-TYPE: \t x-ncbi-data/x-asn-text-urlencoded\r",
                      eMIME_T_NcbiData, eMIME_AsnText, eENCOD_Url) );
  ASSERT( s_CheckMIME("x-ncbi-data/x-eeee",
                      eMIME_T_NcbiData, eMIME_Unknown, eENCOD_None) );

  ASSERT( !s_CheckMIME("content-TYPE : x-ncbi-data/x-unknown\r",
                       eMIME_T_NcbiData, eMIME_Unknown, eENCOD_None) );
  ASSERT( s_CheckMIME("text/html",
                       eMIME_T_Text, eMIME_Html, eENCOD_None) );
  ASSERT( !s_CheckMIME("", eMIME_T_NcbiData, eMIME_Unknown, eENCOD_None) );
  ASSERT( !s_CheckMIME(0, eMIME_T_NcbiData, eMIME_Unknown, eENCOD_None) );

  return 0;
}

#endif /* TEST_MODULE__CONNUTIL */
