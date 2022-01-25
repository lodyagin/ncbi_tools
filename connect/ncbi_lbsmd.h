#ifndef CONNECT___NCBI_LBSMD__H
#define CONNECT___NCBI_LBSMD__H

/*  $Id: ncbi_lbsmd.h,v 6.16 2007/04/20 01:55:30 kazimird Exp $
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
 *
 */

#include "ncbi_servicep.h"
#include <connect/ncbi_heapmgr.h>


#ifdef __cplusplus
extern "C" {
#endif


const SSERV_VTable* SERV_LBSMD_Open(SERV_ITER    iter,
                                    SSERV_Info** info,
                                    HOST_INFO*   host_info,
                                    int/*bool*/  no_dispd);


/* Get configuration file name. Returned '\0'-terminated string
 * is to be free()'d by a caller when no longer needed.
 * Return NULL if no configuration file name is available.
 * LBSMD_FastHeapAccess() was set to "eOff" and there is a cached copy
 * of LBSM heap kept in-core, it will be released by this call.
 */
extern NCBI_XCONNECT_EXPORT const char* LBSMD_GetConfig(void);


/* Get (perhaps cached) copy of LBSM heap, which is guaranteed to be
 * current for given time "time".  If "time" passed as 0, current time
 * of the day is assumed.  Return NULL if the copy operation cannot
 * be performed (due to various reasons, including the original
 * LBSM shmem to be obsolete).
 * Returned heap (if non-NULL) has a serial number reflecting which
 * shmem segment has been used to get the snapshot.  The serial number
 * is negated for newer heap structure, which has dedicated version
 * entry format.  Older heap structure uses SLBSM_OldEntry instead,
 * and has TTLs for entries instead of expiration times.  The returned
 * copy must be passed to (MT-locked by the caller) HEAP_Destroy() when
 * no longer needed.
 * The copy may be cached in-core, the only way to release it is to
 * call LBSMD_GetConfig() provided that LBSM_FastHeapAccess() has
 * been set to "eOff" (which is the default setting).
 */
extern NCBI_XCONNECT_EXPORT HEAP LBSMD_GetHeapCopy(TNCBI_Time time);


extern NCBI_XCONNECT_EXPORT ESwitch LBSMD_FastHeapAccess(ESwitch onoff);


/* Host info getters */
typedef const void* LBSM_HINFO;

int LBSM_HINFO_CpuCount(LBSM_HINFO hinfo);


int LBSM_HINFO_TaskCount(LBSM_HINFO hinfo);


int/*bool*/ LBSM_HINFO_LoadAverage(LBSM_HINFO hinfo, double lavg[2]);


int/*bool*/ LBSM_HINFO_Status(LBSM_HINFO hinfo, double status[2]);


#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* CONNECT___NCBI_LBSMD__H */
