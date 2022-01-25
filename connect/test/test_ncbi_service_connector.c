/*  $Id: test_ncbi_service_connector.c,v 6.22 2002/10/28 15:47:12 lavr Exp $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   Standard test for the service connector
 *
 */

#include "../ncbi_ansi_ext.h"
#include "../ncbi_priv.h"
#include <connect/ncbi_service_connector.h>
#include <stdlib.h>
/* This header must go last */
#include "test_assert.h"


int main(int argc, const char* argv[])
{
    static char obuf[128] = "UUUUUZZZZZZUUUUUUZUZUZZUZUZUZUZUZ\n";
    const char* service = argc > 1 ? argv[1] : "bounce";
    const char* host = argc > 2 ? argv[2] : "www.ncbi.nlm.nih.gov";
    CONNECTOR connector;
    SConnNetInfo *info;
    EIO_Status status;
    STimeout  timeout;
    char ibuf[1024];
    CONN conn;
    size_t n;

    CORE_SetLOGFILE(stderr, 0/*false*/);

    info = ConnNetInfo_Create(service);
    strcpy(info->host, host);
    if (argc > 3) {
        strncpy0(obuf, argv[3], sizeof(obuf) - 2);
        obuf[n = strlen(obuf)] = '\n';
        obuf[++n]              = 0;
    }

    connector = SERVICE_CreateConnectorEx(service, fSERV_Any, info, 0);
    ConnNetInfo_Destroy(info);

    if (!connector)
        CORE_LOG(eLOG_Fatal, "Failed to create service connector");

    if (CONN_Create(connector, &conn) != eIO_Success)
        CORE_LOG(eLOG_Fatal, "Failed to create connection");

    timeout.sec  = 5;
    timeout.usec = 123456;

    CONN_SetTimeout(conn, eIO_ReadWrite, &timeout);

    if (CONN_Write(conn, obuf, strlen(obuf), &n) != eIO_Success ||
        n != strlen(obuf)) {
        CONN_Close(conn);
        CORE_LOG(eLOG_Fatal, "Error writing to connection");
    }

    if (CONN_Wait(conn, eIO_Read, &timeout) != eIO_Success) {
        CONN_Close(conn);
        CORE_LOG(eLOG_Fatal, "Error waiting for reading");
    }

    status = CONN_Read(conn, ibuf, n, &n, eIO_ReadPersist);
    if (status != eIO_Success) {
        if (!n)
            CONN_Close(conn);
        CORE_LOG(n ? eLOG_Error : eLOG_Fatal, "Error reading from connection");
    }

    CORE_LOGF(eLOG_Note,
              ("%d bytes read from service (%s):\n%.*s",
               (int) n, CONN_GetType(conn), (int) n, ibuf));
    CONN_Close(conn);

#if 0
    CORE_LOG(eLOG_Note, "Trying ID1 service");

    info = ConnNetInfo_Create(service);
    connector = SERVICE_CreateConnectorEx("ID1", fSERV_Any, info);
    ConnNetInfo_Destroy(info);

    if (!connector)
        CORE_LOG(eLOG_Fatal, "Service ID1 not available");

    if (CONN_Create(connector, &conn) != eIO_Success)
        CORE_LOG(eLOG_Fatal, "Failed to create connection");

    if (CONN_Write(conn, "\xA4\x80\x02\x01\x02\x00\x00", 7, &n) !=
        eIO_Success || n != 7) {
        CONN_Close(conn);
        CORE_LOG(eLOG_Fatal, "Error writing to service ID1");
    }

    if (CONN_Read(conn, ibuf, sizeof(ibuf), &n, eIO_ReadPlain) != eIO_Success){
        CONN_Close(conn);
        CORE_LOG(eLOG_Fatal, "Error reading from service ID1");
    }

    CORE_LOGF(eLOG_Note, ("%d bytes read from service ID1", n));
    CONN_Close(conn);
#endif

    return 0/*okay*/;
}


/*
 * --------------------------------------------------------------------------
 * $Log: test_ncbi_service_connector.c,v $
 * Revision 6.22  2002/10/28 15:47:12  lavr
 * Use "ncbi_ansi_ext.h" privately and use strncpy0()
 *
 * Revision 6.21  2002/09/24 15:10:09  lavr
 * Fix test not to dereference NULL pointer resulting from failed connection
 *
 * Revision 6.20  2002/08/07 16:38:08  lavr
 * EIO_ReadMethod enums changed accordingly; log moved to end
 *
 * Revision 6.19  2002/03/22 19:47:41  lavr
 * Test_assert.h made last among the include files
 *
 * Revision 6.18  2002/03/21 22:02:16  lavr
 * Change default server from "ray" into "www.ncbi.nlm.nih.gov"
 *
 * Revision 6.17  2002/02/05 21:45:55  lavr
 * Included header files rearranged
 *
 * Revision 6.16  2002/01/16 21:23:15  vakatov
 * Utilize header "test_assert.h" to switch on ASSERTs in the Release mode too
 *
 * Revision 6.15  2001/09/24 20:36:22  lavr
 * Adjusted parameters in SERVICE_CreateConnectorEx()
 *
 * Revision 6.14  2001/06/11 22:17:28  lavr
 * Wait-for-reading timeout made finite
 *
 * Revision 6.13  2001/06/07 17:53:47  lavr
 * Persistent reading from test connection
 *
 * Revision 6.12  2001/06/01 16:19:10  lavr
 * Added (ifdef'ed out) an internal test connection to service ID1
 *
 * Revision 6.11  2001/05/11 16:05:41  lavr
 * Change log message corrected
 *
 * Revision 6.10  2001/05/11 15:38:01  lavr
 * Print connector type along with read data
 *
 * Revision 6.9  2001/04/24 21:42:43  lavr
 * Brushed code to use CORE_LOG facility only.
 *
 * Revision 6.8  2001/01/25 17:13:22  lavr
 * Added: close/free everything on program exit: useful to check memory leaks
 *
 * Revision 6.7  2001/01/23 23:22:34  lavr
 * debug_printout was given an enum value (instead of "boolean" 1)
 *
 * Revision 6.6  2001/01/09 15:35:20  lavr
 * Removed header <unistd.h>, unknown on WinNT
 *
 * Revision 6.5  2001/01/08 23:13:19  lavr
 * C/C++ -> C only
 *
 * Revision 6.4  2001/01/08 22:42:42  lavr
 * Further development of the test-suite
 *
 * Revision 6.3  2001/01/03 22:40:24  lavr
 * Minor adjustment
 *
 * Revision 6.2  2000/12/29 18:25:27  lavr
 * More tests added (still not yet complete).
 *
 * Revision 6.1  2000/10/20 17:31:07  lavr
 * Initial revision
 *
 * ==========================================================================
 */
