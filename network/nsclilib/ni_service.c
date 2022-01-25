/*  $RCSfile: ni_service.c,v $  $Revision: 6.6 $  $Date: 2002/04/23 17:57:54 $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   NCBI Named Service Client (based on SERVICE connector)
 *
 * --------------------------------------------------------------------------
 * $Log: ni_service.c,v $
 * Revision 6.6  2002/04/23 17:57:54  lavr
 * Recognize "INFINITE" as a timeout from registry/environment
 *
 * Revision 6.5  2002/04/16 21:58:06  lavr
 * Few fixes after test compilation and run
 *
 * Revision 6.4  2002/04/16 21:33:24  lavr
 * Add compatibility for service parameters taken as for WWW disp
 *
 * Revision 6.3  2002/03/23 04:21:04  lavr
 * Typo corrected
 *
 * Revision 6.2  2002/03/22 22:22:45  lavr
 * Try to do the best in setting up proper timeouts
 *
 * Revision 6.1  2001/02/21 22:09:15  lavr
 * Initial revision
 *
 * ==========================================================================
 */

#include <ncbi.h>
#include <ncbinet.h>
#include <connect/ncbi_connection.h>
#include <connect/ncbi_service_connector.h>


/*********************************
 * INTERNALS
 */

/* Hard-coded constants, environment parameter names & defaults
 * NOTE:: These are taken from ni_www.c for backward compatibility.
 *        Their use will eventually get deprecated in the future...
 */

#define SRV_SECTION         "NET_SERV"

#define ENV_ENGINE_HOST     "SRV_ENGINE_HOST"
#define DEF_ENGINE_HOST     DEF_CONN_HOST

#define ENV_ENGINE_PORT     "SRV_ENGINE_PORT"
#define DEF_ENGINE_PORT     DEF_CONN_PORT

#define ENV_ENGINE_URL      "SRV_ENGINE_URL"
#define DEF_ENGINE_URL      DEF_CONN_PATH

#define ENV_TIMEOUT         "SRV_CONN_TIMEOUT"
#define DEF_TIMEOUT         DEF_CONN_TIMEOUT

#define ENV_CONN_TRY        "SRV_CONN_TRY"
#define DEF_CONN_TRY        DEF_CONN_MAX_TRY

#define ENV_HTTP_PROXY_HOST "SRV_HTTP_PROXY_HOST"
#define DEF_HTTP_PROXY_HOST DEF_CONN_HTTP_PROXY_HOST

#define ENV_HTTP_PROXY_PORT "SRV_HTTP_PROXY_PORT"
#define DEF_HTTP_PROXY_PORT DEF_CONN_HTTP_PROXY_PORT

#define ENV_PROXY_HOST      "SRV_PROXY_HOST"
#define DEF_PROXY_HOST      DEF_CONN_PROXY_HOST

#define ENV_DEBUG_PRINTOUT  "SRV_DEBUG_PRINTOUT"
#define DEF_DEBUG_PRINTOUT  DEF_CONN_DEBUG_PRINTOUT

#define ENV_NO_LB_DIRECT    "SRV_NO_LB_DIRECT"
#define DEF_NO_LB_DIRECT    DEF_CONN_LB_DISABLE


/* Static functions
 */

static Int2 LIBCALLBACK s_AsnRead(Pointer p, CharPtr buff, Uint2 len)
{
    size_t n_read = 0;
    CONN_Read((CONN) p, buff, len, &n_read, eIO_Plain);
    return (Int2) n_read;
}


static Int2 LIBCALLBACK s_AsnWrite(Pointer p, CharPtr buff, Uint2 len)
{
    size_t n_written = 0;
    CONN_Write((CONN) p, buff, len, &n_written);
    return (Int2) n_written;
}


static void LIBCALLBACK s_AsnErrorFunc(Int2 type, CharPtr message)
{
    ErrPostEx(SEV_ERROR, 88, type, message);
}


/* The interface implementaion functions
 */

static NI_DispatcherPtr s_GenericInit
(CharPtr configFile, CharPtr configSection, Boolean showMonitor,
 CharPtr lastDispatcher, Int2 lastDispLen)
{
    NI_DispatcherPtr disp = (NI_DispatcherPtr)MemNew(sizeof(NI_Dispatcher));
    disp->interface = eNII_Service;
    
    if ( lastDispatcher )
        StringNCpy_0(lastDispatcher, "NCBI Named Service", lastDispLen);
    
    disp->motd = StringSave("Load-balancing service mapping facility");
    disp->adminInfo = StringSave("Anton Lavrentiev (lavr@ncbi.nlm.nih.gov)");
    disp->referenceCount = 1;
    return disp;
}


static NI_DispatcherPtr s_SetDispatcher
(NI_DispatcherPtr disp, CharPtr host, CharPtr svc, int timeout,
 Int4 dispserialnum, ValNodePtr encryption, Boolean useOutServ)
{
    return s_GenericInit(0, 0, 0, 0, 0);
}


static NI_HandPtr s_GenericGetService
(NI_DispatcherPtr disp, CharPtr configFile, CharPtr configSection,
 CharPtr defService, Boolean hasResource)
{
    SConnNetInfo* net_info;
    Char          str[64];
    NI_HandPtr    result;
    double        valf;
    CONN          conn;
    int           val;
    CONNECTOR     c;

    {{ /* alternate service name */
        static const Char ENV_PREFIX[] = "NI_SERVICE_NAME_";
        CharPtr envName = (CharPtr)MemNew(sizeof(ENV_PREFIX) +
                                          StringLen(configSection));
        StringCpy(envName, ENV_PREFIX);
        StringCat(envName, configSection);
        NI_GetEnvParamEx(configFile, configSection, envName, "SERVICE_NAME",
                         str, sizeof(str), defService);
        MemFree(envName);
    }}
    if (!(net_info = ConnNetInfo_Create(str))) {
        ErrPostEx(SEV_ERROR, 0, 1, "[Service NI Client] "
                  " Cannot set parameters for service \"%s\"", str);
        return 0;
    }

    /* Now override default parameters with proprietary parameters
     * of older WWW service dispatcher -- should go obsolete soon... */

    /* alternate dispatcher's host name & port */
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_ENGINE_HOST,
                   net_info->host, sizeof(net_info->host), DEF_ENGINE_HOST);

    NI_GetEnvParam(configFile, SRV_SECTION, ENV_ENGINE_PORT,
                   str, sizeof(str), "");
    val = atoi(str);
    net_info->port = val > 0 ? val : DEF_ENGINE_PORT;

    /* alternate the dispatcher's CGI path */
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_ENGINE_URL,
                   net_info->path, sizeof(net_info->path), DEF_ENGINE_URL);

    /* alternate HTTP proxy host & port */
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_HTTP_PROXY_HOST,
                   net_info->http_proxy_host,
                   sizeof(net_info->http_proxy_host), DEF_HTTP_PROXY_HOST);
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_HTTP_PROXY_PORT,
                   str, sizeof(str), "");
    val = atoi(str);
    net_info->http_proxy_port = val > 0 ? val : DEF_HTTP_PROXY_PORT;

    /* alternate non-transparent CERN-like firewall proxy server */
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_PROXY_HOST,
                   net_info->proxy_host, sizeof(net_info->proxy_host),
                   DEF_PROXY_HOST);

    /* alternate connection timeout */
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_TIMEOUT,
                   str, sizeof(str), "");
    if (strlen(str) < 3  ||  StringNICmp(str, "infinite", strlen(str)) != 0) {
        valf = atof(str);
        if (valf <= 0)
            valf = DEF_TIMEOUT;
        if (!net_info->timeout)
            net_info->timeout   = &net_info->tmo;
        net_info->timeout->sec  =
            (unsigned int) valf;
        net_info->timeout->usec =
            (unsigned int) ((valf - net_info->timeout->sec) * 1000000);
    } else
        net_info->timeout = 0;

    /* alternate max. number of attemts to establish connection */
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_CONN_TRY,
                   str, sizeof(str), "");
    val = atoi(str);
    net_info->max_try = val > 0 ? val : DEF_CONN_TRY;

    NI_GetEnvParam(configFile, SRV_SECTION, ENV_DEBUG_PRINTOUT,
                   str, sizeof(str), DEF_DEBUG_PRINTOUT);
    if (*str  &&  (StringICmp(str, "1"   ) == 0 ||
                   StringICmp(str, "true") == 0 ||
                   StringICmp(str, "yes" ) == 0))
        net_info->debug_printout = 1/*true*/;

    /* whether to prohibit the use of local LBSMD */
    NI_GetEnvParam(configFile, SRV_SECTION, ENV_NO_LB_DIRECT,
                   str, sizeof(str), DEF_NO_LB_DIRECT);
    if (*str  &&  (StringICmp(str, "0"   )  != 0 &&
                   StringICmp(str, "false") != 0 &&
                   StringICmp(str, "no" )   != 0))
        net_info->lb_disable = 1/*true*/;

    /* establish connection to the server */
    if (!(c = SERVICE_CreateConnectorEx(defService, fSERV_Any, net_info, 0)) ||
        CONN_Create(c, &conn) != eIO_Success) {
        ErrPostEx(SEV_ERROR, 0, 1, "[Service NI Client] "
                  " Service \"%s\" unusable", net_info->service);
        ConnNetInfo_Destroy(net_info);
        return 0;
    }
    CONN_SetTimeout(conn, eIO_Open,      net_info ? net_info->timeout : 0);
    CONN_SetTimeout(conn, eIO_ReadWrite, net_info ? net_info->timeout : 0);
    CONN_SetTimeout(conn, eIO_Close,     net_info ? net_info->timeout : 0);

    /* open ASN i/o, etc. */
    result = (NI_HandPtr) MemNew(sizeof(NI_Handle));
    result->extra_proc_info = conn;
    result->raip = AsnIoNew(ASNIO_BIN | ASNIO_IN,  (FILE*) 0,
                            (void*) conn, s_AsnRead,  (IoFuncType) 0);
    result->waip = AsnIoNew(ASNIO_BIN | ASNIO_OUT, (FILE*) 0,
                            (void*) conn, (IoFuncType) 0, s_AsnWrite);
    AsnIoSetErrorMsg(result->raip, s_AsnErrorFunc);
    AsnIoSetErrorMsg(result->waip, s_AsnErrorFunc);
    result->hostname = StringSave(net_info ? net_info->client_host : "");
    result->disp = disp;
    disp->referenceCount++;
    ConnNetInfo_Destroy(net_info);

    return result;
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


static Int2 s_ServiceDisconnect(NI_HandPtr mhp)
{
    s_EndServices(mhp->disp);
    CONN_Close((CONN)mhp->extra_proc_info);
    AsnIoClose(mhp->raip);
    AsnIoClose(mhp->waip);
    MemFree(mhp->hostname);
    MemFree(mhp);
    return 0;
} 


/* Exported table of interface functions
 */

static const NIInterface s_NII_Service = {
    s_GenericInit,
    s_SetDispatcher,
    s_GenericGetService,
    s_ServiceDisconnect,
    s_EndServices
};
const NIInterface *g_NII_Service = &s_NII_Service;

/* EOF */
