/*  $Id: con_stateless.c,v 6.3 1999/08/02 22:20:52 vakatov Exp $
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
*   Implement a CONNECTOR to perform a stateless connection to NCBI services.
*
*  NOTE: See in "connectr.h" for the detailed specification of the underlying
*        connector("CONNECTOR", "SConnectorTag") methods and structures.
*
* --------------------------------------------------------------------------
* $Log: con_stateless.c,v $
* Revision 6.3  1999/08/02 22:20:52  vakatov
* In the HTTP header, use "\r\n" instead of just "\n"
*
* Revision 6.2  1999/07/26 18:07:05  vakatov
* Redesigned the test (see #ifdef TEST_MODULE__CON_STATELESS) to use
* public standard name of conf.parameters(as #def'd in "connutil.h" now)
*
* Revision 6.1  1999/07/21 21:43:26  vakatov
* Initial revision
*
* ==========================================================================
*/

#include <ncbi.h>
#include <con_url.h>
#include <con_stateless.h>


/***********************************************************************
 *  STATIC -- for the LB mode only
 ***********************************************************************/

#if defined(LB_DIRECT)

#  include <lbapi.h>

typedef struct {
  char*     service_name;
  unsigned* skip_ip;  /* addresses to ignore */
  size_t    size;     /* of "skip_ip" */
  size_t    n_skip;   /* # of valid entries in "skip_ip" */
} SUserdata;

static SUserdata* s_UserdataCreate(const char* service_name)
{
  SUserdata* data = (SUserdata*) MemNew(sizeof(SUserdata));
  data->service_name = StringSave(service_name);
  data->size   = 0;
  data->n_skip = 0;
  return data;
}

static void s_UserdataDestroy(SUserdata* data)
{
  MemFree(data->skip_ip);
  MemFree(data->service_name);
  MemFree(data);
}


static void s_AdjustInfo(SNetConnInfo* info, void* data, Uint4 conn_try)
{
  SUserdata* x_data = (SUserdata*) data;
  unsigned   ip_addr;
  Char       str_addr[16];

  if ( info->lb_disable )
    return;

  if (conn_try == 0)
    x_data->n_skip = 0;

  ip_addr = (Uint4)LBGetIPAddress(x_data->service_name,
                                  0, x_data->skip_ip, x_data->n_skip);
  if (ip_addr == 0) {
    if (x_data->n_skip == 0)
      return;
    x_data->n_skip = 0;
    ip_addr = LBGetIPAddress(x_data->service_name, 0, 0, 0);
  }
  if (ip_addr == 0 )
    return;

  if (x_data->size == x_data->n_skip) {
    x_data->size += 32;
    x_data->skip_ip = (unsigned*)
      Realloc(x_data->skip_ip, x_data->size * sizeof(unsigned));
  }
  x_data->skip_ip[++x_data->n_skip] = ip_addr;

  if ( !Uint4toInaddr((Uint4)ip_addr, str_addr, sizeof(str_addr)) )
    return;

  MemFree(info->host);  info->host = StringSave(str_addr);
  info->port = info->ncbid_port;
  MemFree(info->path);  info->path = StringSave(info->ncbid_path);
  info->http_proxy_adjusted = FALSE;
  NetConnInfo_AdjustForHttpProxy(info);
}

static void s_AdjustCleanup(void* data)
{
  s_UserdataDestroy((SUserdata*) data);
}

#endif /* LB_DIRECT */



/***********************************************************************
 *  EXTERNAL -- the connector's "constructor"
 ***********************************************************************/


NLM_EXTERN CONNECTOR STATELESS_CreateConnector
(const char* service_name)
{
  /* get default connection parameters */
  SNetConnInfo* info = NetConnInfo_Create(0, 0);

  /* set connection flags */
  URLC_Flags flags = URLC_AUTO_RECONNECT | URLC_SURE_FLUSH;

  /* compose the user's part of HTTP header */
  static const char fmt_user_header[] = "\
User-Agent: NCBI_Client from %s\r\n\
Content-Type: x-ncbi-data/x-unknown\r\n";
#define MAX_HOST_LEN 64
  char client_hostname[MAX_HOST_LEN];
  char user_header[sizeof(fmt_user_header) + sizeof(client_hostname) - 3];
  if ( !GetHostName(client_hostname, sizeof(client_hostname)) )
       StringNCpy_0(client_hostname, "UNKNOWN", sizeof(client_hostname));
  sprintf(user_header, fmt_user_header, client_hostname);
  ASSERT( strlen(user_header) < sizeof(user_header) );


  /* compose the URL cmd-line args */
  {{
    static const char fmt_args[] = "service=%s&address=%s&platform=%s";
    const char* client_platform = Nlm_PlatformName();
    size_t args_size = sizeof(fmt_args) - 6 +
      strlen(service_name) + strlen(client_hostname) + strlen(client_platform);

    MemFree(info->args);
    info->args = (char*) MemNew(args_size);
    sprintf(info->args, fmt_args,
            service_name, client_hostname, client_platform);
    ASSERT( strlen(info->args) == args_size - 1 );
  }}

  /* create the connector, then cleanup & return */
  {{
    CONNECTOR connector = URL_CreateConnectorEx
      (info, user_header, flags,
#if defined(LB_DIRECT)
       s_AdjustInfo, s_UserdataCreate(service_name), s_AdjustCleanup
#else
       0, 0, 0
#endif
       );

    NetConnInfo_Destroy(&info);
    return connector;
  }}
}



/***********************************************************************
 *  TEST
 ***********************************************************************/

#ifdef TEST_MODULE__CON_STATELESS

#define TEST_ALTERNATE_URLS   0
#define USE_HTTP_PROXY        0

#define DISPD_SCRIPT          "dispd"
#define NCBID_SCRIPT          "ncbid"

#define TEST_CONF_FILE        "ncbi"
#define TEST_CONF_SECTION     "NET_SERV"

#define TEST_HTTP_PROXY_HOST  "www.ncbi.nlm.nih.gov"
#define TEST_ENGINE_HOST      "www.ncbi.nlm.nih.gov"
#define TEST_ENGINE_PORT      "80"

#define TEST_ENGINE_PATH \
  "/Service/" DISPD_SCRIPT ".cgi"
#define TEST_NCBID_PATH \
  "/Service/" NCBID_SCRIPT ".cgi"

#define TEST_LOGFILE          "con_stateless.log"
#define TEST_SERVICE_NAME     "IO_BOUNCE"
#define TEST_DEBUG_PRINTOUT   "yes"


#include <conntest.h>

extern
#ifdef __cplusplus
"C"
#endif
Int2 Main(void)
{
  STimeout  timeout;
  CONNECTOR connector;
  FILE*     log_file;

#if (TEST_ALTERNATE_URLS > 0)
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_HOST, TEST_ENGINE_HOST);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_PORT, TEST_ENGINE_PORT);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_ENGINE_PATH, TEST_ENGINE_PATH);
  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           DEF_CONN_ENGINE_ARGS, TEST_ENGINE_ARGS);

  Nlm_TransientSetAppParam(TEST_CONF_FILE, TEST_CONF_SECTION,
                           CFG_CONN_NCBID_PATH, TEST_NCBID_PATH);
#endif

#if (USE_HTTP_PROXY > 0)
  Nlm_TransientSetAppParam(TEST_CONF_FILE, TEST_CONF_SECTION,
                           CFG_CONN_HTTP_PROXY_HOST, TEST_HTTP_PROXY_HOST);
#endif

  Nlm_TransientSetAppParam(DEF_CONN_CONF_FILE, DEF_CONN_CONF_SECTION,
                           CFG_CONN_DEBUG_PRINTOUT, TEST_DEBUG_PRINTOUT);

  timeout.sec  = 5;
  timeout.usec = 123456;
  log_file = FileOpen(TEST_LOGFILE, "wb");

  connector = STATELESS_CreateConnector(TEST_SERVICE_NAME);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_SINGLE_BOUNCE_PRINT);

  connector = STATELESS_CreateConnector(TEST_SERVICE_NAME);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_SINGLE_BOUNCE_CHECK);

  connector = STATELESS_CreateConnector(TEST_SERVICE_NAME);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_ALL);

  connector = STATELESS_CreateConnector(TEST_SERVICE_NAME);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_ALL);

  FileClose(log_file);
  return 0;
}

#endif  /* TEST_MODULE__CON_STATELESS */
