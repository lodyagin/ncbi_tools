/*  $Id: ncbi_service_lbsmd.c,v 6.26 2001/07/05 17:00:23 lavr Exp $
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
 *   Low-level API to resolve NCBI service name to the server meta-address
 *   with the use of NCBI Load-Balancing Service Mapper (LBSMD).
 *   UNIX only!
 *
 * --------------------------------------------------------------------------
 * $Log: ncbi_service_lbsmd.c,v $
 * Revision 6.26  2001/07/05 17:00:23  lavr
 * Severity change eLOG_Warning->eLOG_Error in log for failed heap lock
 *
 * Revision 6.25  2001/07/03 20:50:21  lavr
 * RAND_MAX included in the interval search.
 * Log error message if unable to lock the heap.
 *
 * Revision 6.24  2001/06/25 15:36:53  lavr
 * s_GetNextInfo now takes one additional argument for host environment
 *
 * Revision 6.23  2001/06/19 19:12:01  lavr
 * Type change: size_t -> TNCBI_Size; time_t -> TNCBI_Time
 *
 * Revision 6.22  2001/05/24 21:27:37  lavr
 * Skip pre-expired servers (with zero expiration time)
 *
 * Revision 6.21  2001/05/03 16:58:16  lavr
 * FIX: Percent is taken of local bonus coef instead of the value itself
 *
 * Revision 6.20  2001/05/03 16:35:58  lavr
 * Local bonus coefficient modified: meaning of negative value changed
 *
 * Revision 6.19  2001/04/26 20:20:01  lavr
 * Better way of choosing local server with a tiny (e.g. penalized) status
 *
 * Revision 6.18  2001/04/24 21:35:03  lavr
 * Penalty interface added; treatment of bonus coefficient for local servers
 *
 * Revision 6.17  2001/03/29 21:14:35  lavr
 * 'penalty' changed to 'fine'
 *
 * Revision 6.16  2001/03/27 23:36:13  lavr
 * Uptime and penalty members added to service struct in shared memory
 *
 * Revision 6.15  2001/03/07 22:23:52  lavr
 * BUGFIX: Open returns 0 if service not found in local LBSM table
 *
 * Revision 6.14  2001/03/07 20:45:37  lavr
 * s_CalculateStatus moved to ncbi_lbsm.c and made global LBSM_CalculateStatus
 *
 * Revision 6.13  2001/03/05 23:11:22  lavr
 * Special treatment of negative rate (for statically configured servers)
 *
 * Revision 6.12  2001/03/02 20:10:21  lavr
 * Typos fixed
 *
 * Revision 6.11  2001/01/08 22:41:39  lavr
 * No important mods here; only some cosmetic changes
 *
 * Revision 6.10  2000/12/29 18:07:32  lavr
 * Modified to increase rate for locally running services (BONUS);
 * 'expiration period' replaced with 'expiration time' returned by
 * SERV_GetNextInfo.
 *
 * Revision 6.9  2000/12/06 22:07:36  lavr
 * LBSM_LOCAL_SVC_BONUS introduced to increase rate of a service, which
 * runs locally, when we calculate balancing information
 *
 * Revision 6.8  2000/10/20 17:21:40  lavr
 * VTable updated (empty 'Update' added)
 * Uninited variable use fixed
 *
 * Revision 6.7  2000/10/05 21:39:40  lavr
 * Marked UNIX-only
 *
 * Revision 6.6  2000/06/20 19:01:35  lavr
 * Code cleaned to removed some warnings from gcc compiler
 *
 * Revision 6.5  2000/05/31 23:12:23  lavr
 * First try to assemble things together to get working service mapper
 *
 * Revision 6.4  2000/05/23 19:02:49  lavr
 * Server-info now includes rate; verbal representation changed
 *
 * Revision 6.3  2000/05/22 16:53:12  lavr
 * Rename service_info -> server_info everywhere (including
 * file names) as the latter name is more relevant
 *
 * Revision 6.2  2000/05/12 21:43:00  lavr
 * Cleaned up for the C++ compilation, etc.
 *
 * Revision 6.1  2000/05/12 18:43:16  lavr
 * First working revision
 *
 * ==========================================================================
 */

#include "ncbi_lbsm.h"
#include "ncbi_lbsm_ipc.h"
#include "ncbi_priv.h"
#include "ncbi_servicep_lbsmd.h"
#include <connect/ncbi_ansi_ext.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>


/* Default rate increase if svc runs locally */
#define SERV_LBSMD_LOCAL_SVC_BONUS 1.1


#ifdef __cplusplus
extern "C" {
#endif

    static SSERV_Info* s_GetNextInfo(SERV_ITER iter, char** env);
    static int/*bool*/ s_Penalize(SERV_ITER iter, double penalty);
    
    static const SSERV_VTable s_op = {
        s_GetNextInfo, 0, s_Penalize, 0, "LBSMD"
    };

#ifdef __cplusplus
} /* extern "C" */
#endif


static SSERV_Info* s_GetNextInfo(SERV_ITER iter, char** env)
{
    double total = 0.0, point = -1.0, access = 0.0, p = 0.0, status;
    struct SSERV_List {
        double               status;
        const SLBSM_Host*    host;
        const SLBSM_Service* svc;
    }* list = 0;
    size_t i, n_list, n_max_list;
    int/*bool*/ stateless_only;
    SSERV_Info* info = 0;
    TSERV_Type type;
    HEAP heap;

    if (!LBSM_Shmem_Lock()) {
        CORE_LOG_SYS_ERRNO(eLOG_Error, "Cannot lock LBSMD heap");
        return 0;
    }

    stateless_only = (iter->type & fSERV_StatelessOnly) != 0;
    type = iter->type & ~fSERV_StatelessOnly;
    n_list = n_max_list = 0;
    if ((heap = LBSM_Shmem_Attach()) != 0) {
        const SLBSM_Service *svc = 0;
        const SLBSM_Host *host;
        LBSM_Setup(heap);

        while ((svc = LBSM_LookupService(iter->service, svc)) != 0) {
            if (!svc->info.rate/*not working*/ || !svc->info.time/*expired*/)
                continue;

            if (type && !(type & (int)svc->info.type))
                continue; /* type doesn't match */

            if (stateless_only && svc->info.sful) {
                /* Skip stateful-only non-CGI (NCBID and standalone) svc */
                if (!(svc->info.type & fSERV_Http))
                    continue;
            }

            for (i = 0; i < iter->n_skip; i++)
                if (SERV_EqualInfo(&svc->info, iter->skip[i]))
                    break;
            if (i < iter->n_skip)
                continue; /* excluded */

            if (!(host = LBSM_LookupHost(svc->info.host)) &&
                svc->info.rate > 0.0) {
                CORE_LOG(eLOG_Error, "No host entry for dynamic service");
                continue; /* no host information for non-static service */
            }

            status = LBSM_CalculateStatus(&host->load, svc->info.rate,
                                          svc->info.flag, svc->fine);
            if (status <= 0.0)
                continue;

            /* This server should be taken into consideration */
            if (n_list == n_max_list) {
                size_t n = n_max_list + 10;

                struct SSERV_List* temp;
                if (!(temp = (struct SSERV_List*)
                      (list ? realloc(list, sizeof(*temp)*n) :
                       malloc(sizeof(*temp)*n))))
                    break;

                list = temp;
                n_max_list = n;
            }
            if (svc->info.host == iter->preferred_host) {
                if (svc->info.coef > 0.0)
                    status *= svc->info.coef;
                else {
                    status *= SERV_LBSMD_LOCAL_SVC_BONUS;
                    if (svc->info.coef < 0.0 && access < status) {
                        access = status;
                        point  = total; /* Latch this local server */
                        p      = -svc->info.coef;
                    }
                }
            }
            total              += status;
            list[n_list].status = total;
            list[n_list].host   = host;
            list[n_list].svc    = svc;
            n_list++;
        }

        if (list) {
            size_t info_size;

            /* We will take pre-chosen local server only if its status is
               not less than p% of the average rest status; otherwise, we
               ignore the server, and apply the general procedure by seeding
               a random point. */
            if (point < 0.0 || access*(n_list - 1) < p*0.01*(total - access))
                point = (total * rand()) / (double) RAND_MAX;
            for (i = 0; i < n_list; i++) {
                if (point <= list[i].status)
                    break;
            }
            assert(i < n_list);

            info_size = SERV_SizeOfInfo(&list[i].svc->info);
            if ((info = (SSERV_Info*) malloc(info_size)) != 0) {
                memcpy(info, &list[i].svc->info, info_size);
                info->rate = list[i].status - (i ? list[i - 1].status : 0.0);
                info->time += (TNCBI_Time) time(0);  /* Expiration time now */
            }
            if (env) {
                if ((host = list[i].host) != 0 && host->env)
                    *env = strdup((const char*) host + host->env);
                else
                    *env = 0;
            }

            free(list);
        }

        LBSM_Shmem_Detach(heap);
    }

    LBSM_Shmem_Unlock();

    return info;
}


static int/*bool*/ s_Penalize(SERV_ITER iter, double penalty)
{
    return LBSM_SubmitPenalty(iter->service, iter->last->type,
                              penalty, iter->last->host);
}


/***********************************************************************
 *  EXTERNAL
 ***********************************************************************/

const SSERV_VTable *SERV_LBSMD_Open(SERV_ITER iter)
{
    SSERV_Info* info;

    /* Daemon is running if LBSM_LBSMD returns 1: mutex exists but read
       operation failed with errno == EAGAIN (the mutex is busy) */
    if (LBSM_LBSMD(0) <= 0 || errno != EAGAIN)
        return 0;
    srand((int)time(0) + (int)getpid());
    if (!(info = s_GetNextInfo(iter, 0)))
        return 0;
    free(info);
    return &s_op;
}
