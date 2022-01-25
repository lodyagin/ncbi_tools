#ifndef NCBI_LBSM__H
#define NCBI_LBSM__H

/*  $Id: ncbi_lbsm.h,v 6.12 2001/03/29 21:13:56 lavr Exp $
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
 * $Log: ncbi_lbsm.h,v $
 * Revision 6.12  2001/03/29 21:13:56  lavr
 * Default feedback file name added; 'penalty' changed to 'fine'
 *
 * Revision 6.11  2001/03/27 23:36:13  lavr
 * Uptime and penalty members added to service struct in shared memory
 *
 * Revision 6.10  2001/03/07 20:44:27  lavr
 * Added: LBSM_CalculateStatus
 *
 * Revision 6.9  2001/03/06 23:54:57  lavr
 * Minor fix: floating point added to floating point values
 *
 * Revision 6.8  2001/03/05 23:09:58  lavr
 * More comments to LBSM helper functions
 *
 * Revision 6.7  2000/12/29 17:58:05  lavr
 * Pretty printed; BONUS constant added for services running locally.
 *
 * Revision 6.6  2000/12/06 22:06:53  lavr
 * LBSM_LOCAL_SVC_BONUS introduced to increase rate of a service, which
 * runs locally, when we calculate balancing information
 *
 * Revision 6.5  2000/05/23 21:05:29  lavr
 * Memory leaks fixed (appeared after recent server-info structure rearrangement)
 *
 * Revision 6.4  2000/05/23 19:02:48  lavr
 * Server-info now includes rate; verbal representation changed
 *
 * Revision 6.3  2000/05/22 16:53:10  lavr
 * Rename service_info -> server_info everywhere (including
 * file names) as the latter name is more relevant
 *
 * Revision 6.2  2000/05/18 14:11:19  lavr
 * Yet another cosmetic update
 *
 * Revision 6.1  2000/05/12 19:11:59  lavr
 * First working revision
 *
 * ==========================================================================
 */

#include <connect/ncbi_heapmgr.h>
#include <connect/ncbi_server_info.h>

#ifdef __cplusplus
extern "C" {
#endif


#if defined(_DEBUG) && !defined(NDEBUG)
/* #define LBSM_DEBUG 1 */
#endif

#define LBSM_DEFAULT_CFGFILE  "lbsmd.cfg"
#define LBSM_DEFAULT_LOGFILE  "lbsmd.log"
#define LBSM_DEFAULT_RUNFILE  "lbsmd.pid"
#define LBSM_DEFAULT_FEEDFILE "/tmp/lbsm"
#define LBSM_DEFAULT_PORT     0x7321    /* Port number, host byte order      */
#define LBSM_DEFAULT_TIME     30        /* Exp.time for svc. to be used      */
#define LBSM_DEFAULT_TTL      5         /* Svc unpublished after TTL expires */


/* Load parameters of the machine
 */
typedef  struct {
    double       avg;
    unsigned int nrtask;
    unsigned int nrproc;
    double       avgBLAST;
    unsigned int nrtaskBLAST[5];
    unsigned int queuesBLAST[3];
} SLBSM_Load;


/* Uptime parameters of the machine
 */
typedef struct {
    time_t       sys_uptime;      /* Time the system booted up  */
    time_t       start_time;      /* Time the daemon started up */
} SLBSM_Uptime;


/* Service rate limits (for "SLBSM_Service.rate" below)
 */
#define LBSM_SERV_RATE_MIN  0.     /* Not available, e.g. down */
#define LBSM_SERV_RATE_MAX  1000.  /* Perfectly available      */


/* Full service info structure as kept in the heap
 */
typedef struct {
    SHEAP_Block    head;   /* heap manager stuff             */
    size_t         name;   /* Name of this service (offset)  */
    unsigned short ttl;    /* time-to-live (rundown counter) */
    SLBSM_Load     load;   /* load information               */
    SLBSM_Uptime   uptime; /* uptime information             */
    double         fine;   /* feedback on unreachability, %  */
    SSERV_Info     info;   /* server descriptor              */
} SLBSM_Service;


/* Setup the LBSM heap to work with
 */
void LBSM_Setup(HEAP heap);


/* Get next service info from the LBSM heap.
 */
const SLBSM_Service* LBSM_Search
(const char*          service_name,  /* name of service we are looking for */
 const SLBSM_Service* prev_service   /* previously found info (or NULL)    */
 );


/* Put (allocate + copy) new service info to the LBSM heap,
 * and assign this service a default ttl (if original ttl == 0).
 * If same service exists already, replace the information.
 */
int/*bool*/ LBSM_Publish(const SLBSM_Service* service);


/* Traverse through all service info blocks in the LBSM heap,
 * "age" them (decrement "SLBSM_Service.ttl" field),
 * and remove expired service infos (those with "ttl" == 0).
 */
void LBSM_AgeAll(void);


/* Revoke all services, running on 'host' and having non-negative
 * rating coefficients, that is all dynamic, aka lb, (vs. static) services.
 */
void LBSM_RevokeAll(unsigned int host);


/* Calculate status of the service based on machine load and rating.
 */
double LBSM_CalculateStatus
(const SLBSM_Load* load,   /* [in]  machine load parameters                */
 double            rate,   /* [in]  service rate                           */
 ESERV_Flags       flags,  /* [in]  status calculation flags (Reg | Blast) */
 double            fine    /* [in]  feedback penalty from application(s)   */
);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* NCBI_LBSM__H */
