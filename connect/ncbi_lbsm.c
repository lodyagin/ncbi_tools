/*  $Id: ncbi_lbsm.c,v 6.15 2001/07/05 16:53:00 lavr Exp $
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
 *    LBSM client-server data exchange API
 *
 * --------------------------------------------------------------------------
 * $Log: ncbi_lbsm.c,v $
 * Revision 6.15  2001/07/05 16:53:00  lavr
 * Common exit point in LBSM_SubmitPenalty after malloc
 *
 * Revision 6.14  2001/07/03 20:44:33  lavr
 * Removed redundand NULL-pointer checks.
 * Added new consistency check in LBSM_ServiceLookup().
 *
 * Revision 6.13  2001/06/25 15:34:57  lavr
 * Renamed and redesigned functions:
 *   LBSM_Search -> LBSM_LookupService
 *   LBSM_RevokeAll removed
 *   LBSM_LookupHost added, new
 *   LBSM_Publish -> LBSM_PublishService
 *   LBSM_PublishHost added, new
 *
 * Revision 6.12  2001/05/24 21:28:36  lavr
 * LBSM_SERV_RATE_MIN, LBSM_SERV_RATE_MAX -> LBSM_DEFAULT_SERV_RATE
 *
 * Revision 6.11  2001/05/11 15:33:01  lavr
 * BUGFIX: Restore signal handler after penalty submission; made C++-compliant
 *
 * Revision 6.10  2001/04/24 21:40:00  lavr
 * Missing #include <stdio.h> put in.
 *
 * Revision 6.9  2001/04/05 23:11:36  lavr
 * Added function: LBSM_SubmitPenalty
 *
 * Revision 6.8  2001/03/29 21:14:18  lavr
 * 'penalty' changed to 'fine'
 *
 * Revision 6.7  2001/03/27 23:36:13  lavr
 * Uptime and penalty members added to service struct in shared memory
 *
 * Revision 6.6  2001/03/07 20:45:00  lavr
 * Added: LBSM_CalculateStatus (moved from ncbi_service_lbsmd.c and made global)
 *
 * Revision 6.5  2001/03/05 23:09:49  lavr
 * More comments to LBSM helper functions
 *
 * Revision 6.4  2001/03/02 20:08:55  lavr
 * Typo fixed
 *
 * Revision 6.3  2000/05/17 16:13:49  lavr
 * NCBI_* (for ANSI ext) undone - now "ncbi_ansi_ext.h" does good prototyping
 *
 * Revision 6.2  2000/05/15 19:06:08  lavr
 * Use home-made ANSI extensions (NCBI_***)
 *
 * Revision 6.1  2000/05/12 19:14:03  lavr
 * First working revision
 *
 * ==========================================================================
 */

#include <connect/ncbi_ansi_ext.h>
#include <connect/ncbi_socket.h>
#include "ncbi_lbsm.h"
#include "ncbi_priv.h"
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>


static HEAP s_Heap = 0;


void LBSM_Setup(HEAP heap)
{
    s_Heap = heap;
}


const SLBSM_Service* LBSM_LookupService(const char* name,
                                        const SLBSM_Service* s)
{
    const SLBSM_Entry* e = &s->entry;
    if (e && e->type != eLBSM_Service) {
        CORE_LOG(eLOG_Error, "Invalid preceding entry in service lookup");
        return 0;
    }
    while ((e = (SLBSM_Entry*) HEAP_Walk(s_Heap, &e->head)) != 0) {
        if ((short) e->head.flag && e->type == eLBSM_Service) {
            s = (const SLBSM_Service*) e;
            if (strcasecmp((char*) s + s->name, name) == 0)
                return s;
        }
    }
    return 0;
}


const SLBSM_Host* LBSM_LookupHost(unsigned addr)
{
    const SLBSM_Entry* e = 0;
    while ((e = (SLBSM_Entry*) HEAP_Walk(s_Heap, &e->head)) != 0) {
        if ((short) e->head.flag && e->type == eLBSM_Host) {
            const SLBSM_Host* h = (const SLBSM_Host*) e;
            if (h->addr == addr)
                return h;
        }
    }
    return 0;
}


int/*bool*/ LBSM_PublishService(const SLBSM_Service* s)
{
    const SLBSM_Service* f = 0;
    SLBSM_Service* b;
    const char* name;
    size_t size;

    if (!s || s->entry.type != eLBSM_Service)
        return 0;

    name = (const char*) s + s->name;
    while ((f = LBSM_LookupService(name, f)) != 0)
        if (SERV_EqualInfo(&s->info, &f->info))
            break;
    if (f)
        HEAP_Free(s_Heap, (SHEAP_Block*) &f->entry.head);

    size = sizeof(*s) - sizeof(s->entry.head) - sizeof(s->info) +
        SERV_SizeOfInfo(&s->info) + strlen(name)+1;
    if (!(b = (SLBSM_Service*) HEAP_Alloc(s_Heap, size)))
        return 0;

    /* Copy service to the heap, preserving block header */
    memcpy((char*) b + sizeof(b->entry.head),
           (char*) s + sizeof(s->entry.head), size);
    if (!b->entry.ttl)
        b->entry.ttl = LBSM_DEFAULT_TTL;

    return 1;
}


int/*bool*/ LBSM_PublishHost(const SLBSM_Host* h)
{
    const SLBSM_Host* f = 0;
    const char* env;
    SLBSM_Host* b;
    size_t size;

    if (!h || h->entry.type != eLBSM_Host)
        return 0;

    if ((f = LBSM_LookupHost(h->addr)) != 0)
        HEAP_Free(s_Heap, (SHEAP_Block*) &f->entry.head);

    env = h->env ? (const char*) h + h->env : 0;
    size = sizeof(*h) - sizeof(h->entry.head) + (env ? strlen(env)+1 : 0);
    if (!(b = (SLBSM_Host*) HEAP_Alloc(s_Heap, size)))
        return 0;

    /* Copy host to the heap, preserving block header */
    memcpy((char*) b + sizeof(b->entry.head),
           (char*) h + sizeof(h->entry.head), size);
    if (!b->entry.ttl)
        b->entry.ttl = LBSM_DEFAULT_TTL;

    return 1;
}


void LBSM_AgeAll(void)
{        
    SLBSM_Entry* e = 0;
    while ((e = (SLBSM_Entry*) HEAP_Walk(s_Heap, &e->head)) != 0) {
        if ((short)e->head.flag && !--(e->ttl))
            HEAP_Free(s_Heap, &e->head);
    }
}


static double s_CalculateStatusBLAST(const SLBSM_Load* load)
{
    return 1.0/LBSM_DEFAULT_SERV_RATE +
        LBSM_DEFAULT_SERV_RATE*(load->nrproc - load->avgBLAST)*
        (load->queuesBLAST[0] ?
         exp(-8.0*load->queuesBLAST[2]/load->queuesBLAST[0]) : 1.0);
}


double LBSM_CalculateStatus(const SLBSM_Load* load, double rate,
                            ESERV_Flags flags, double fine)
{
    double status;

    if (rate >= 0.0) {
        if (!(flags & fSERV_Blast)) {
            double cpuload = load->avg;

            /* Penalize machine if loadavg > nr.of CPUs */
            if (load->avg > load->nrproc && load->nrproc)
                cpuload *= exp(load->avg/load->nrproc - 1.0);
            status = 1.0/LBSM_DEFAULT_SERV_RATE +
                LBSM_DEFAULT_SERV_RATE*load->nrproc/(cpuload*load->nrtask+1.0);
        } else
            status = s_CalculateStatusBLAST(load);
    } else
        status = -LBSM_DEFAULT_SERV_RATE;
    status *= rate/LBSM_DEFAULT_SERV_RATE;
    status *= 1.0 - (fine < 0.0 ? 0.0 :
                     (fine > 100.0 ? 100.0 : fine))/100.0;
    return status;
}


#ifdef __cplusplus
extern "C" {
    static void (*oldhandler)(int);
}
#else
static void (*oldhandler)(int);
#endif /* __cplusplus */


int/*bool*/ LBSM_SubmitPenalty(const char* sname, ESERV_Type type,
                               double penalty, unsigned int host)
{
    const char* tname = SERV_TypeStr(type);
    struct sockaddr_un un;
    char hname[32];
    int s, retval;
    char* msg;

    if (!sname || !tname)
        return 0;
    if (host && SOCK_ntoa(host, hname, sizeof(hname)) != 0)
        return 0;
    if (!(msg = (char*) malloc(strlen(sname) + 1 + strlen(tname) + 1 +
                               4/*penalty*/ + (host ? strlen(hname) : 0) +
                               2/*\n\0*/)))
        return 0;

    if (penalty < 0.0)
        penalty = 0.0;
    if (penalty > 100.0)
        penalty = 100.0;

    retval = 0; /* assume worst */
    if (sprintf(msg, "%s %s %#.0f%s%s\n", sname, tname, penalty,
                host ? " " : "", host ? hname : "") > 0) {
        if ((oldhandler = signal(SIGPIPE, SIG_IGN)) != SIG_ERR) {
            un.sun_family = PF_UNIX;
            strcpy(un.sun_path, LBSM_DEFAULT_FEEDFILE);
            if ((s = socket(PF_UNIX, SOCK_STREAM, 0)) >= 0) {
                if (connect(s, (struct sockaddr*) &un, sizeof(un)) == 0 &&
                    write(s, msg, strlen(msg)) > 0)
                    retval = 1/*success*/;
                close(s);
            }
            signal(SIGPIPE, oldhandler);
        }
    }
    free(msg);
    return retval;
}
