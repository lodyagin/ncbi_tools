/*  $Id: ncbi_lbsm.c,v 6.8 2001/03/29 21:14:18 lavr Exp $
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
#include "ncbi_lbsm.h"
#include <math.h>
#include <string.h>


static HEAP s_Heap = 0;


void LBSM_Setup(HEAP heap)
{
    s_Heap = heap;
}


const SLBSM_Service* LBSM_Search(const char* name, const SLBSM_Service* s)
{
    while ((s = (SLBSM_Service*) HEAP_Walk(s_Heap, s ? &s->head : 0)) != 0) {
        if ((short)s->head.flag && strcasecmp((char*) s + s->name, name) == 0)
            return s;
    }
    return 0;
}


int/*bool*/ LBSM_Publish(const SLBSM_Service* s)
{
    SLBSM_Service* b;
    const SLBSM_Service* f = 0;
    size_t size = sizeof(*s) - sizeof(s->head) - sizeof(s->info) +
        SERV_SizeOfInfo(&s->info) + strlen((char*) s + s->name)+1;

    while ((f = LBSM_Search((char*) s + s->name, f)) != 0) {
        if (SERV_EqualInfo(&s->info, &f->info))
            break;
    }
    if (f)
        HEAP_Free(s_Heap, (SHEAP_Block*) &f->head);

    if (!(b = (SLBSM_Service*) HEAP_Alloc(s_Heap, size)))
        return 0;

    /* Copy service to the heap, preserving block header */
    memcpy((char*) b + sizeof(b->head), (char*) s + sizeof(s->head), size);
    if (!b->ttl)
        b->ttl = LBSM_DEFAULT_TTL;

    return 1;
}


void LBSM_AgeAll(void)
{        
    SLBSM_Service* s = 0;
    while ((s = (SLBSM_Service*) HEAP_Walk(s_Heap, s ? &s->head : 0)) != 0) {
        if ((short)s->head.flag && !--(s->ttl))
            HEAP_Free(s_Heap, &s->head);
    }
}


void LBSM_RevokeAll(unsigned int host)
{
    SLBSM_Service* s = 0;
    while ((s = (SLBSM_Service*) HEAP_Walk(s_Heap, s ? &s->head : 0)) != 0) {
        if ((short)s->head.flag && s->info.host == host && s->info.rate >= 0.0)
            HEAP_Free(s_Heap, &s->head);
    }
}


static double s_CalculateStatusBLAST(const SLBSM_Load* load)
{
    return 1.0/LBSM_SERV_RATE_MAX +
        LBSM_SERV_RATE_MAX*(load->nrproc - load->avgBLAST)*
        (load->queuesBLAST[2] > 0 ?
         exp(-4.0*load->queuesBLAST[2]/load->queuesBLAST[0]) : 1.0);
}


double LBSM_CalculateStatus(const SLBSM_Load* load, double rate,
                            ESERV_Flags flags, double fine)
{
    double status;

    if (rate >= 0.0) {
        if (!(flags & fSERV_Blast)) {
            double cpuload = load->avg;

            /* Penalize machine if loadavg > nr.of CPUs */
            if (load->avg > load->nrproc)
                cpuload *= exp(load->avg/load->nrproc - 1.0);
            status = 1.0/LBSM_SERV_RATE_MAX +
                LBSM_SERV_RATE_MAX*load->nrproc/(cpuload*load->nrtask + 1.0);
        } else
            status = s_CalculateStatusBLAST(load);
    } else
        status = -LBSM_SERV_RATE_MAX;
    status *= rate/LBSM_SERV_RATE_MAX;
    status *= 1.0 - (fine < 0.0 ? 0.0 :
                     (fine > 100.0 ? 100.0 : fine))/100.0;
    return status;
}
