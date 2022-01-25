/*  $Id: ncbi_service_lbsmd.c,v 6.17 2001/03/29 21:14:35 lavr Exp $
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
#include "ncbi_servicep_lbsmd.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>


#ifdef __cplusplus
extern "C" {
#endif

    static SSERV_Info* s_GetNextInfo(SERV_ITER iter);

    static const SSERV_VTable s_op = { s_GetNextInfo, 0, 0 };

#ifdef __cplusplus
} /* extern "C" */
#endif


static SSERV_Info* s_GetNextInfo(SERV_ITER iter)
{
    double total = 0.0, point, status;
    struct SSERV_List {
        double               status;
        const SLBSM_Service* svc;
    }* list = 0;
    size_t i, n_list, n_max_list;
    int/*bool*/ stateless_only;
    SSERV_Info* info = 0;
    TSERV_Type type;
    HEAP heap;

    if (!LBSM_Shmem_Lock())
        return 0;

    stateless_only = (iter->type & fSERV_StatelessOnly) != 0;
    type = iter->type & ~fSERV_StatelessOnly;
    n_list = n_max_list = 0;
    if ((heap = LBSM_Shmem_Attach()) != 0) {
        const SLBSM_Service *svc = 0;
        LBSM_Setup(heap);

        while ((svc = LBSM_Search(iter->service, svc)) != 0) {
            if (!svc->info.rate)
                continue;

            if (type && !(type & (int)svc->info.type))
                continue;

            if (stateless_only && svc->info.sful) {
                /* Skip stateful-only non-CGI (NCBID and standalone) svc */
                if (!(svc->info.type & fSERV_Http))
                    continue;
            }

            for (i = 0; i < iter->n_skip; i++)
                if (SERV_EqualInfo(&svc->info, iter->skip[i]))
                    break;
            if (i < iter->n_skip)
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
            list[n_list].svc    = svc;
            status              = LBSM_CalculateStatus(&svc->load,
                                                       svc->info.rate,
                                                       svc->info.flag,
                                                       svc->fine);
            if (svc->info.host == iter->preferred_host)
                status *= SERV_LBSMD_LOCAL_SVC_BONUS;
            total              += status;
            list[n_list].status = total;
            n_list++;
        }

        if (list) {
            size_t info_size;

            point = (total * rand()) / (double)RAND_MAX;
            for (i = 0; i < n_list; i++) {
                if (point < list[i].status)
                    break;
            }
            assert(i < n_list);
            
            info_size = SERV_SizeOfInfo(&list[i].svc->info);
            if ((info = (SSERV_Info*) malloc(info_size)) != 0) {
                memcpy(info, &list[i].svc->info, info_size);
                info->rate = list[i].status - (i ? list[i - 1].status : 0.0);
                info->time += time(0);  /* Set 'expiration time' */
            }
            
            free(list);
        }

        LBSM_Shmem_Detach(heap);
    }

    LBSM_Shmem_Unlock();

    return info;
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
    if (!(info = s_GetNextInfo(iter)))
        return 0;
    free(info);
    return &s_op;
}
