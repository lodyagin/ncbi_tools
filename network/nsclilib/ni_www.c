/*  $RCSfile: ni_www.c,v $  $Revision: 4.19 $  $Date: 1999/02/18 18:46:48 $
* ==========================================================================
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
* ==========================================================================
*
* Authors:  Denis Vakatov, Sergei Shavirin, Eugene Yaschenko
*
* File Description:
*   Basic client code for the new(HTTPD-based) NCBI dispatchers
*   (regular and connectionless protocols)
*
* --------------------------------------------------------------------------
* $Log: ni_www.c,v $
* Revision 4.19  1999/02/18 18:46:48  shavirin
* Added few new functions, those do not interfere with existing API.
*
* Revision 4.18  1998/12/15 17:28:52  vakatov
* Added config parameter "SRV_PROXY_TRANSPARENT" -- to indicate the use
* of transparent/non-transparent CERN-like proxy
*
* Revision 4.17  1998/12/08 17:03:58  vakatov
* Look for SRV_*** connection parameters in #SRV_SECTION(="NET_SERV")
* rather than in "configSection"(service-specific)
*
* Revision 4.16  1998/09/08 17:59:08  vakatov
* Added WWW/Firewall network interface
*
* Revision 4.15  1998/08/05 20:19:37  vakatov
* s_AsnRead():   Always make sure that the network conection is
* established by the moment of data read;
* + minor logic fix due to the fix in "ncbicli.c" R5.0
*
* Revision 4.14  1998/05/20 20:22:05  vakatov
* Added two new config. options SRV_PROXY_HOST and SRV_PROXY_PORT in order
* to allow one to connect via a proxy HTTP server
*
* Revision 4.13  1998/05/08 15:26:59  vakatov
* [LB_DIRECT]  now can skip IP addresses(for the dispatcher hosts) which
* were already offered by the load-balancing daemon but failed by some reason
*
* Revision 4.12  1998/05/01 19:16:56  vakatov
* Connection establishing defaults: 2 attempts with 30 sec timeout(was 3 & 20)
*
* Revision 4.11  1998/04/30 23:02:25  vakatov
* s_ServiceDisconnect():  check "sinfo->nic" for NULL
*
* Revision 4.10  1998/04/09 15:36:36  vakatov
* [LB_DIRECT]  Introduced SRV_NO_LB_DIRECT env/conf variable to let one
* disable the use of local "nsdaemon" load-balance info
*
* Revision 4.9  1998/04/08 23:04:21  vakatov
* Count references to the Dispatcher struct to avoid its premature freeing
*
* Revision 4.8  1998/04/03 22:46:43  vakatov
* [#LB_DIRECT]  Allowed shortcut around micasa/sucasa(for the hosts
* running "lbdaemon")
*
* Revision 4.7  1998/04/02 20:21:10  vakatov
* Added possibility to printout service reply headers
*
* Revision 4.6  1998/03/30 17:50:20  vakatov
* Ingrafted to the main NCBI CVS tree
*
* ==========================================================================
*/

#include <ncbi.h>
#include <ncbinet.h>
#include <ncbicli.h>
#include <lbapi.h>


/*********************************
 * INTERNALS
 */

/* Hard-coded constants, environment parameter names & defaults
 */

#define SRV_SECTION         "NET_SERV"

#define ENV_PROXY_HOST      "SRV_PROXY_HOST"
#define DEF_PROXY_HOST      ""

#define ENV_PROXY_PORT      "SRV_PROXY_PORT"
#define DEF_PROXY_PORT      80

#define ENV_PROXY_TRANSPARENT  "SRV_PROXY_TRANSPARENT"
#define DEF_PROXY_TRANSPARENT  ""

#define ENV_ENGINE_HOST     "SRV_ENGINE_HOST"
#define DEF_ENGINE_HOST     "www.ncbi.nlm.nih.gov"

#define ENV_ENGINE_PORT     "SRV_ENGINE_PORT"
#define DEF_ENGINE_PORT     80

#define ENV_ENGINE_URL      "SRV_ENGINE_URL"
#define DEF_ENGINE_URL      "/Service/nph-dispd.cgi"

#define ENV_CONN_TIMEOUT    "SRV_CONN_TIMEOUT"
#define DEF_CONN_TIMEOUT    30

#define ENV_CONN_TRY        "SRV_CONN_TRY"
#define DEF_CONN_TRY        3

#define ENV_DEBUG_PRINTOUT  "SRV_DEBUG_PRINTOUT"
#define DEF_DEBUG_PRINTOUT  ""

#ifdef LB_DIRECT
#define DEF_ENGINE_LB_URL   "/Service/nph-ncbid.cgi"
#define ENV_NO_LB_DIRECT    "SRV_NO_LB_DIRECT"
#endif /* LB_DIRECT */


/* Local typedefs
 */

typedef struct {
  Char       client_host[64];  /* effective client hostname */
  Char       service[256];     /* requested service */
  Char       disp_host[64];    /* dispatcher:  host */
  Char       disp_path[256];   /* dispatcher:  path(to CGI script) */
  Uint2      disp_port;        /* dispatcher:  service port */
  NIC        nic;              /* handle to the connection */ 
  STimeout   timeout;          /* i/o(connection and read/write) timeout */
  Uint4      conn_try;         /* max. number of attempts to establish conn. */
  ENIC_Agent agent;
  Uint4      flags;            /* to be passed to NIC_GetService() */
  Boolean    no_lb_direct;     /* prohibit the use of local "nsdaemon" */
  AsnIoBSPtr aibsp;
} ServiceInfo;


/* Static functions
 */

#ifdef LB_DIRECT
static unsigned s_AlternateDispatcher(ServiceInfo *sinfo,
                                      unsigned *skip_ip, size_t n_skip)
{
  Uint4 ip_addr;
  Char  str_addr[16];

  if ( sinfo->no_lb_direct )
    return 0;

  ip_addr = (Uint4)LBGetIPAddress(sinfo->service, 0, skip_ip, n_skip);
  if (!ip_addr  ||  !Uint4toInaddr(ip_addr, str_addr, sizeof(str_addr)))
    return 0;

  StringNCpy_0(sinfo->disp_host, str_addr, sizeof(sinfo->disp_host));
  sinfo->disp_port = DEF_ENGINE_PORT;
  StringNCpy_0(sinfo->disp_path, DEF_ENGINE_LB_URL, sizeof(DEF_ENGINE_LB_URL));
  return ip_addr;
}
#endif

static Boolean s_Connect(ServiceInfo *sinfo)
{
  Uint4 conn_try;
#ifdef LB_DIRECT
  unsigned skip_ip[MAX_NUM_HOSTS];
  MemSet(skip_ip, '\0', sizeof(skip_ip));
#endif

  ASSERT( !sinfo->nic );
  for (conn_try = 0;  !sinfo->nic  &&  conn_try < sinfo->conn_try;
       conn_try++) {
#ifdef LB_DIRECT
    skip_ip[conn_try] = s_AlternateDispatcher(sinfo, skip_ip,
                                              conn_try%MAX_NUM_HOSTS);
#endif
    sinfo->nic = NIC_GetService
      (sinfo->service,
       sinfo->disp_host, sinfo->disp_port, sinfo->disp_path,
       &sinfo->timeout, sinfo->agent, sinfo->client_host,
       sinfo->aibsp->bsp, sinfo->flags);
  }

  return (Boolean)(sinfo->nic != 0);
}

Boolean NET_Connect(VoidPtr handle)
{
    Uint4 conn_try;
    ServiceInfo *sinfo;
        
#ifdef LB_DIRECT
    unsigned skip_ip[MAX_NUM_HOSTS];
    MemSet(skip_ip, '\0', sizeof(skip_ip));
#endif

    sinfo = (ServiceInfo *) handle;
    
    if(sinfo->nic != NULL) {
        NIC_CloseService(sinfo->nic);
        sinfo->nic = NULL;
    }

    for (conn_try = 0;  !sinfo->nic  &&  conn_try < sinfo->conn_try;
         conn_try++) {
#ifdef LB_DIRECT
        skip_ip[conn_try] = s_AlternateDispatcher(sinfo, skip_ip,
                                                  conn_try%MAX_NUM_HOSTS);
#endif
        sinfo->nic = NIC_GetService
            (sinfo->service,
             sinfo->disp_host, sinfo->disp_port, sinfo->disp_path,
             &sinfo->timeout, sinfo->agent, sinfo->client_host,
             NULL, sinfo->flags);
    }
    
    return (Boolean)(sinfo->nic != 0);
}


static Int2 LIBCALLBACK s_AsnRead(Pointer p, CharPtr buff, Uint2 len)
{
  ServiceInfo *sinfo = (ServiceInfo *)p;
    
  if (BSLen(sinfo->aibsp->bsp) > 0) {
    /* connect, if needed;  flush the content of output BSP, if any */
    if ( sinfo->nic ) {
      if (sinfo->agent == eNIC_WWWDirect) {
        NIC_CloseService(sinfo->nic);
        sinfo->nic = 0;
        if ( !s_Connect(sinfo) ) /* it also flushes the content of out. BSP */
          return 0;
      } else { /* flush the content of output BSP */
        Nlm_BSUnitPtr bsup;
        SOCK sock = NIC_GetSOCK(sinfo->nic);
        ASSERT(sinfo->agent == eNIC_WWWClient);
        for (bsup = sinfo->aibsp->bsp->chain;  bsup  &&  bsup->len;
             bsup = bsup->next) {
          const void *ptr = Nlm_HandLock(bsup->str);
          ESOCK_ErrCode err_code = SOCK_Write(sock, ptr, (Uint4)bsup->len, 0);
          Nlm_HandUnlock(bsup->str);
          if (err_code != eSOCK_ESuccess)
            return 0;
        }
      }
    } else {
      if ( !s_Connect(sinfo) ) /* it also flushes the content of output BSP */
        return 0;
    }

    /* cleanup the output BSP content */
    BSSeek(sinfo->aibsp->bsp, 0, SEEK_SET);
    BSDelete(sinfo->aibsp->bsp, BSLen(sinfo->aibsp->bsp));
  }
  else { /* no data to (pre-)write;  just connect if necessary */
    if (!sinfo->nic  &&  !s_Connect(sinfo))
      return 0;
  }


  {{ /* read */
    Uint4 n_read;
    SOCK_Read(NIC_GetSOCK(sinfo->nic), buff, (Uint4)len, &n_read);
    return (Int2)n_read;
  }}
}


static void LIBCALLBACK s_AsnErrorFunc(Int2 type, CharPtr message)
{
  ErrPostEx(SEV_ERROR, 77, type, message);
}


/* The interface implementaion functions
 */

static NI_DispatcherPtr s_GenericInit
(ENIInterface interface,
 CharPtr configFile, CharPtr configSection, Boolean showMonitor,
 CharPtr lastDispatcher, Int2 lastDispLen)
{
  NI_DispatcherPtr disp = (NI_DispatcherPtr)MemNew(sizeof(NI_Dispatcher));
  disp->interface = interface;

  if ( lastDispatcher )
    StringNCpy_0(lastDispatcher, "WWW Based Dispatcher", lastDispLen);

  disp->motd = StringSave("This system uses HTTPD servers on NCBI to "
                          "initialize network connection");
  disp->adminInfo = StringSave("");
  disp->referenceCount = 1;
  return disp;
}


static NI_DispatcherPtr s_GenericInitWWW
(CharPtr configFile, CharPtr configSection, Boolean showMonitor,
 CharPtr lastDispatcher, Int2 lastDispLen)
{
  return s_GenericInit(eNII_WWW, configFile, configSection, showMonitor,
                       lastDispatcher, lastDispLen);
}


static NI_DispatcherPtr s_GenericInitWWWFirewall
(CharPtr configFile, CharPtr configSection, Boolean showMonitor,
 CharPtr lastDispatcher, Int2 lastDispLen)
{
  return s_GenericInit(eNII_WWWFirewall, configFile, configSection,
                       showMonitor, lastDispatcher, lastDispLen);
}


static NI_DispatcherPtr s_GenericInitWWWDirect
(CharPtr configFile, CharPtr configSection, Boolean showMonitor,
 CharPtr lastDispatcher, Int2 lastDispLen)
{
  return s_GenericInit(eNII_WWWDirect, configFile, configSection, showMonitor,
                       lastDispatcher, lastDispLen);
}


static NI_DispatcherPtr s_SetDispatcherWWW
(NI_DispatcherPtr disp, CharPtr host, CharPtr svc, int timeout,
 Int4 dispserialnum, ValNodePtr encryption, Boolean useOutServ)
{
  return s_GenericInitWWW(0, 0, 0, 0, 0);
}

static NI_DispatcherPtr s_SetDispatcherWWWFirewall
(NI_DispatcherPtr disp, CharPtr host, CharPtr svc, int timeout,
 Int4 dispserialnum, ValNodePtr encryption, Boolean useOutServ)
{
  return s_GenericInitWWWFirewall(0, 0, 0, 0, 0);
}

static NI_DispatcherPtr s_SetDispatcherWWWDirect
(NI_DispatcherPtr disp, CharPtr host, CharPtr svc, int timeout,
 Int4 dispserialnum, ValNodePtr encryption, Boolean useOutServ)
{
  return s_GenericInitWWWDirect(0, 0, 0, 0, 0);
}


static NI_HandPtr s_GenericGetService
(NI_DispatcherPtr disp, CharPtr configFile, CharPtr configSection,
 CharPtr defService, Boolean hasResource)
{
  NI_HandPtr   result = (NI_HandPtr)MemNew(sizeof(NI_Handle));
  ServiceInfo *sinfo  = (ServiceInfo *)MemNew(sizeof(ServiceInfo));

  {{ /* alternate service name */
    static const Char ENV_PREFIX[] = "NI_SERVICE_NAME_";
    CharPtr envName = (CharPtr)MemNew(sizeof(ENV_PREFIX) +
                                      StringLen(configSection));
    StringCpy(envName, ENV_PREFIX);
    StringCat(envName, configSection);
    NI_GetEnvParamEx(configFile, configSection, envName, "SERVICE_NAME",
                     sinfo->service, sizeof(sinfo->service), defService);
    MemFree(envName);
  }}

  /* alternate the dispatcher's host name */
  NI_GetEnvParam(configFile, SRV_SECTION, ENV_ENGINE_HOST,
                 sinfo->disp_host, sizeof(sinfo->disp_host), DEF_ENGINE_HOST);

  {{ /* alternate the dispatcher's port */
    Char str[32];
    int  val;
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_ENGINE_PORT,
                   str, sizeof(str), "");
    val = atoi(str);
    sinfo->disp_port = (Uint2)(val > 0 ? val : DEF_ENGINE_PORT);
  }}

  /* alternate the dispatcher's CGI path */
  NI_GetEnvParam(configFile, SRV_SECTION, ENV_ENGINE_URL,
                 sinfo->disp_path, sizeof(sinfo->disp_path), DEF_ENGINE_URL);

  {{ /* CERN-like firewall proxy server? */
    Char proxy_host[sizeof(sinfo->disp_host)];
    /* get the PROXY server name, if specified */
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_PROXY_HOST,
                   proxy_host, sizeof(proxy_host), DEF_PROXY_HOST);

    if ( *proxy_host ) { /* use PROXY server */
      Char str[32 + sizeof(sinfo->disp_host) + sizeof(sinfo->disp_path)];
      int  val;
      Uint2 proxy_port;
      NI_GetEnvParam(configFile, SRV_SECTION, ENV_PROXY_PORT,
                     str, sizeof(str), "");
      val = atoi(str);
      proxy_port = (Uint2)(val > 0 ? val : DEF_PROXY_PORT);
      sprintf(str, "http://%s:%u%s",
              sinfo->disp_host, (unsigned)sinfo->disp_port, sinfo->disp_path);
      StringNCpy_0(sinfo->disp_host, proxy_host, sizeof(sinfo->disp_host));
      sinfo->disp_port = proxy_port;
      StringNCpy_0(sinfo->disp_path, str, sizeof(sinfo->disp_path));

      /* non-transparent proxy? */
      NI_GetEnvParam(configFile, SRV_SECTION, ENV_PROXY_TRANSPARENT,
                     str, sizeof(str), DEF_PROXY_TRANSPARENT);
      if (!*str  ||  (StringICmp(str, "1"   )  &&
                      StringICmp(str, "true")  &&
                      StringICmp(str, "yes" )))
        sinfo->flags |= NIC_CERN_PROXY;

      /* auto-set to FIREWALL mode */
      sinfo->flags |= NIC_FIREWALL;
    }
  }}

  {{ /* alternate the connection timeout */
    Char   str[32];
    double val;
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_CONN_TIMEOUT,
                   str, sizeof(str), "");
    val = atof(str);
    if (val <= 0)
      val = DEF_CONN_TIMEOUT;
    sinfo->timeout.sec  = (Uint4)val;
    sinfo->timeout.usec = (Uint4)((val - sinfo->timeout.sec) * 1000000);
  }}

  {{ /* alternate the max. number of attemts to establish a connection */
    Char str[32];
    int  val;
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_CONN_TRY,
                   str, sizeof(str), "");
    val = atoi(str);
    sinfo->conn_try = (Uint4)((val > 0) ? val : DEF_CONN_TRY);
  }}

  {{ /* alternate the connection flags */
    Char  str[32];
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_DEBUG_PRINTOUT,
                   str, sizeof(str), DEF_DEBUG_PRINTOUT);
    if (*str  &&  (!StringICmp(str, "1"   )  ||
                   !StringICmp(str, "true")  ||
                   !StringICmp(str, "yes" )))
      sinfo->flags |= NIC_DEBUG_PRINTOUT;
  }}

#ifdef LB_DIRECT
  {{ /* if to prohibit the use of local "nsdaemon" info */
    Char  str[32];
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_NO_LB_DIRECT,
                   str, sizeof(str), "");
    sinfo->no_lb_direct = (Boolean)
      (*str  &&  StringICmp(str, "0")  &&
       StringICmp(str, "false")  &&  StringICmp(str, "no"));
  }}
#endif

  /* open ASN i/o, etc. */
  result->extra_proc_info = sinfo;
  sinfo->aibsp = AsnIoBSOpen("wb", BSNew(1024));
  result->waip = sinfo->aibsp->aip;
  result->raip = AsnIoNew((ASNIO_BIN | ASNIO_IN), (FILE *)0,
                          (void *)sinfo, s_AsnRead, (IoFuncType)0);

  AsnIoSetErrorMsg(result->raip, s_AsnErrorFunc);
  AsnIoSetErrorMsg(result->waip, s_AsnErrorFunc);
  switch ( disp->interface ) {
  case eNII_WWWFirewall:
    sinfo->flags |= NIC_FIREWALL;
  case eNII_WWW:
    sinfo->agent = eNIC_WWWClient;
    break;
  case eNII_WWWDirect:
    sinfo->agent = eNIC_WWWDirect;
    break;
  default:
    ASSERT ( 0 );
    return 0;
  }
  result->hostname = StringSave("");
  result->disp = disp;
  disp->referenceCount++;

  return result;
}

/*  ---- Here are set of temporary functions to be deleted */

VoidPtr NET_GenericGetService (CharPtr defService, 
                               CharPtr configSection, 
                               Int4 interface)
{
    ServiceInfo *sinfo  = (ServiceInfo *)MemNew(sizeof(ServiceInfo));
    
    {{ /* alternate service name */
        static const Char ENV_PREFIX[] = "NI_SERVICE_NAME_";
        CharPtr envName = (CharPtr)MemNew(sizeof(ENV_PREFIX) +
                                          StringLen(configSection));
        StringCpy(envName, ENV_PREFIX);
        StringCat(envName, configSection);
        NI_GetEnvParamEx(NULL, configSection, envName, "SERVICE_NAME",
                         sinfo->service, sizeof(sinfo->service), defService);
        MemFree(envName);
    }}
    
    /* alternate the dispatcher's host name */
    NI_GetEnvParam(NULL, SRV_SECTION, ENV_ENGINE_HOST,
                   sinfo->disp_host, sizeof(sinfo->disp_host), DEF_ENGINE_HOST);
    {{ /* alternate the dispatcher's port */
        Char str[32];
        int  val;
        NI_GetEnvParam(NULL, SRV_SECTION, ENV_ENGINE_PORT,
                       str, sizeof(str), "");
        val = atoi(str);
        sinfo->disp_port = (Uint2)(val > 0 ? val : DEF_ENGINE_PORT);
    }}
    
    /* alternate the dispatcher's CGI path */
    NI_GetEnvParam(NULL, SRV_SECTION, ENV_ENGINE_URL,
                   sinfo->disp_path, sizeof(sinfo->disp_path), DEF_ENGINE_URL);
    
    {{ /* CERN-like firewall proxy server? */
        Char proxy_host[sizeof(sinfo->disp_host)];
        /* get the PROXY server name, if specified */
        NI_GetEnvParam(NULL, SRV_SECTION, ENV_PROXY_HOST,
                       proxy_host, sizeof(proxy_host), DEF_PROXY_HOST);
        
        if ( *proxy_host ) { /* use PROXY server */
            Char str[32 + sizeof(sinfo->disp_host) + sizeof(sinfo->disp_path)];
            int  val;
            Uint2 proxy_port;
            NI_GetEnvParam(NULL, SRV_SECTION, ENV_PROXY_PORT,
                           str, sizeof(str), "");
            val = atoi(str);
            proxy_port = (Uint2)(val > 0 ? val : DEF_PROXY_PORT);
            sprintf(str, "http://%s:%u%s",
                    sinfo->disp_host, (unsigned)sinfo->disp_port, sinfo->disp_path);
            StringNCpy_0(sinfo->disp_host, proxy_host, sizeof(sinfo->disp_host));
            sinfo->disp_port = proxy_port;
            StringNCpy_0(sinfo->disp_path, str, sizeof(sinfo->disp_path));
            
            /* non-transparent proxy? */
            NI_GetEnvParam(NULL, SRV_SECTION, ENV_PROXY_TRANSPARENT,
                           str, sizeof(str), DEF_PROXY_TRANSPARENT);
            if (!*str  ||  (StringICmp(str, "1"   )  &&
                            StringICmp(str, "true")  &&
                            StringICmp(str, "yes" )))
                sinfo->flags |= NIC_CERN_PROXY;
            
            /* auto-set to FIREWALL mode */
            sinfo->flags |= NIC_FIREWALL;
        }
    }}
    
    {{ /* alternate the connection timeout */
        Char   str[32];
        double val;
        NI_GetEnvParam(NULL, SRV_SECTION, ENV_CONN_TIMEOUT,
                       str, sizeof(str), "");
        val = atof(str);
        if (val <= 0)
            val = DEF_CONN_TIMEOUT;
        sinfo->timeout.sec  = (Uint4)val;
        sinfo->timeout.usec = (Uint4)((val - sinfo->timeout.sec) * 1000000);
    }}
    
    {{ /* alternate the max. number of attemts to establish a connection */
        Char str[32];
        int  val;
        NI_GetEnvParam(NULL, SRV_SECTION, ENV_CONN_TRY,
                       str, sizeof(str), "");
        val = atoi(str);
        sinfo->conn_try = (Uint4)((val > 0) ? val : DEF_CONN_TRY);
    }}
    
    {{ /* alternate the connection flags */
        Char  str[32];
        NI_GetEnvParam(NULL, SRV_SECTION, ENV_DEBUG_PRINTOUT,
                       str, sizeof(str), DEF_DEBUG_PRINTOUT);
        if (*str  &&  (!StringICmp(str, "1"   )  ||
                       !StringICmp(str, "true")  ||
                       !StringICmp(str, "yes" )))
            sinfo->flags |= NIC_DEBUG_PRINTOUT;
    }}
    
#ifdef LB_DIRECT 
    {{ /* if to prohibit the use of local "nsdaemon" info */
        Char  str[32];
        NI_GetEnvParam(NULL, SRV_SECTION, ENV_NO_LB_DIRECT,
                       str, sizeof(str), "");
        sinfo->no_lb_direct = (Boolean)
            (*str  &&  StringICmp(str, "0")  &&
             StringICmp(str, "false")  &&  StringICmp(str, "no"));
    }}
#endif
    
    switch (interface) {
    case eNII_WWWFirewall:
        sinfo->flags |= NIC_FIREWALL;
    case eNII_WWW:
        sinfo->agent = eNIC_WWWClient;
        break;
    case eNII_WWWDirect:
        sinfo->agent = eNIC_WWWDirect;
        break;
    default:
        ASSERT ( 0 );
        return 0;
    }
    
    return (VoidPtr) sinfo;
}

static Int2 s_EndServices(NI_DispatcherPtr disp)
{
    ASSERT ( disp->referenceCount > 0 );
    if (--disp->referenceCount == 0) {
        MemFree(disp->adminInfo);
        MemFree(disp->motd);
        MemFree(disp);
    }
    return 0;
}

SOCK NET_GetSOCK(VoidPtr handle)
{
    ServiceInfo *sinfo = (ServiceInfo *) handle;

    if(sinfo != NULL && sinfo->nic != NULL)
        return NIC_GetSOCK(sinfo->nic);

    return NULL;
}

Int2 NET_ServiceDisconnect(VoidPtr handle)
{
    ServiceInfo *sinfo = (ServiceInfo *) handle;
    if (sinfo->nic)
        NIC_CloseService(sinfo->nic);
    MemFree(sinfo);
    return 0;
}

static Int2 s_ServiceDisconnect(NI_HandPtr mhp)
{
  ServiceInfo *sinfo = (ServiceInfo *)mhp->extra_proc_info;
  ByteStorePtr bsp   = sinfo->aibsp->bsp;

  s_EndServices(mhp->disp);
  if ( sinfo->nic )
    NIC_CloseService(sinfo->nic);
  AsnIoBSClose(sinfo->aibsp);
  BSFree(bsp);
  MemFree(sinfo);
  AsnIoClose(mhp->raip);
  MemFree(mhp->hostname);
  MemFree(mhp);
  return 0;
} 


/* Exported table of interface functions
 */
static const NIInterface s_NII_WWW = {
  s_GenericInitWWW,
  s_SetDispatcherWWW,
  s_GenericGetService,
  s_ServiceDisconnect,
  s_EndServices
};
const NIInterface *g_NII_WWW = &s_NII_WWW;

static const NIInterface s_NII_WWWFirewall = {
  s_GenericInitWWWFirewall,
  s_SetDispatcherWWWFirewall,
  s_GenericGetService,
  s_ServiceDisconnect,
  s_EndServices
};
const NIInterface *g_NII_WWWFirewall = &s_NII_WWWFirewall;

static const NIInterface s_NII_WWWDirect = {
  s_GenericInitWWWDirect,
  s_SetDispatcherWWWDirect,
  s_GenericGetService,
  s_ServiceDisconnect,
  s_EndServices
};
const NIInterface *g_NII_WWWDirect = &s_NII_WWWDirect;


/* EOF */
