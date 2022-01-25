/*  $Id: ncbi_lbsm_ipc.c,v 6.17 2001/03/26 18:38:43 lavr Exp $
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
 *    Implementation of the LBSM client-server data exchange API
 *    with the use of SYSV IPC
 *
 * --------------------------------------------------------------------------
 * $Log: ncbi_lbsm_ipc.c,v $
 * Revision 6.17  2001/03/26 18:38:43  lavr
 * Suppression of compiler warning about uninited parameter use (dummy)
 *
 * Revision 6.16  2001/03/19 23:07:51  lavr
 * Log typo corrected :-(
 *
 * Revision 6.15  2001/03/19 23:02:38  lavr
 * LBSM_UnLBSMD now accepts one argument and returns a value
 *
 * Revision 6.14  2001/01/12 23:51:39  lavr
 * Message logging modified for use LOG facility only
 *
 * Revision 6.13  2000/12/06 22:08:49  lavr
 * Added a kludge to avoid double definition of 'union semun' on FreeBSD
 *
 * Revision 6.12  2000/10/20 17:09:30  lavr
 * Minor cosmetic update
 *
 * Revision 6.11  2000/10/06 18:28:04  lavr
 * Conditional definition of 'union semun' (undefined everywhere but IRIX)
 *
 * Revision 6.6  2000/06/20 19:01:35  lavr
 * Code cleaned to removed some warnings from gcc compiler
 *
 * Revision 6.5  2000/06/05 20:44:51  lavr
 * 'Debug' macros modified
 *
 * Revision 6.4  2000/05/17 16:14:32  lavr
 * A tiny little cosmetic change
 *
 * Revision 6.3  2000/05/16 16:36:53  lavr
 * Debug defines modified to reduce gcc warnings
 *
 * Revision 6.2  2000/05/12 21:42:58  lavr
 * Cleaned up for the C++ compilation, etc.
 *
 * Revision 6.1  2000/05/12 19:31:59  lavr
 * First working revision
 *
 * ==========================================================================
 */

#include "ncbi_lbsm_ipc.h"
#include "ncbi_priv.h"
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifndef HAVE_SEMUN
/* This sequence of defines causes 'union semun' be undefined on IRIX */
#ifdef _XOPEN_SOURCE
#  undef _XOPEN_SOURCE
#endif
#define _XOPEN_SOURCE 1
#endif
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>


static int    s_Shmid     = -1;
static int    s_Muxid     = -1;
static char*  s_Shmem     =  0;
static size_t s_ShmemSize =  0;


int LBSM_LBSMD(int/*bool*/ check_n_lock)
{
    struct sembuf lock[2];
    int id = semget(LBSM_MUTEX_KEY, check_n_lock ? 2 : 0,
                    check_n_lock ? (IPC_CREAT | LBSM_SHARED_PROT) : 0);
    if (id < 0)
        return -1;

    s_Muxid = id;
    /* Check & lock daemon presence: done atomically */
    lock[0].sem_num = 0;
    lock[0].sem_op  = 0;
    lock[0].sem_flg = IPC_NOWAIT;
    lock[1].sem_num = 0;
    lock[1].sem_op  = 1;
    lock[1].sem_flg = SEM_UNDO;
    if (semop(id, lock, check_n_lock ? sizeof(lock)/sizeof(lock[0]) : 1) < 0)
        return 1;
    return 0;
}


/* Daemon use: undaemon != 0; Client use: undaemon == 0 => return LBSMD PID */
pid_t LBSM_UnLBSMD(int/*bool*/ undaemon)
{
    static union semun dummy;
    pid_t pid = 0;
    
    if (s_Muxid >= 0) {
        if (!undaemon) {
            if (LBSM_Shmem_Lock()) {
                int shmid = shmget(LBSM_SHMEM_KEY, 0, 0);
                struct shmid_ds shmem_ds;
                
                if (shmid > 0 && shmctl(shmid, IPC_STAT, &shmem_ds) == 0)
                    pid = shmem_ds.shm_cpid;
                LBSM_Shmem_Unlock();
            }
        } else
            semctl(s_Muxid, 0, IPC_RMID, dummy);
    }
    s_Muxid = -1;
    return pid;
}


int/*bool*/ LBSM_Shmem_Lock(void)
{
    struct sembuf lock[2];
    int i;

    for (i = 0; i < 5; i++) {
        lock[0].sem_num = 1;
        lock[0].sem_op  = 0;
        lock[0].sem_flg = 0;
        lock[1].sem_num = 1;
        lock[1].sem_op  = 1;
        lock[1].sem_flg = SEM_UNDO;
        
        if (semop(s_Muxid, lock, sizeof(lock)/sizeof(lock[0])) >= 0)
            return 1;

        if (errno != EINTR)
            break;
    }
    return 0;
}


int/*bool*/ LBSM_Shmem_Unlock(void)
{
    struct sembuf unlock[1];

    unlock[0].sem_num = 1;
    unlock[0].sem_op  = -1;
    unlock[0].sem_flg = IPC_NOWAIT | SEM_UNDO;
    if (semop(s_Muxid, unlock, sizeof(unlock)/sizeof(unlock[0])) < 0)
        return 0;
    return 1;
}


#ifdef __cplusplus
extern "C" {
    static char *s_ExpandShmem(char *, size_t);
}
#endif

static char *s_ExpandShmem(char *mem, size_t newsize)
{
    char *hugebuf = 0;

    if (mem && newsize) {
        /* Reallocation: make a copy of current content */
        assert(mem == s_Shmem);
        if (!(hugebuf = (char *)malloc(s_ShmemSize)))
            CORE_LOG_SYS_ERRNO(eLOG_Error,
                               "Expand: Unable to allocate buffer");
        memcpy(hugebuf, mem, s_ShmemSize);
    }
 
    /* Delete current shared memory region */
    if (mem && shmdt(mem) < 0)
        CORE_LOG_SYS_ERRNO(eLOG_Error, "Expand: Cannot detach shared memory");
    s_Shmem = 0;
    if (s_Shmid >= 0) {
        if (shmctl(s_Shmid, IPC_RMID, 0) < 0)
            CORE_LOG_SYS_ERRNO(eLOG_Warning,
                               "Expand: Unable to delete shared memory");
        s_Shmid = -1;
    }

    if (newsize) {
        /* New allocation or reallocation */
        if ((s_Shmid = shmget(LBSM_SHMEM_KEY, newsize,
                              IPC_CREAT | LBSM_SHARED_PROT)) < 0 ||
            !(s_Shmem = (char *)shmat(s_Shmid, 0, 0)))
            CORE_LOG_SYS_ERRNO(eLOG_Error,
                               "Expand: Unable to create shared memory");
        if (mem) {
            /* Reallocation: restore original content */
            memcpy(s_Shmem, hugebuf, s_ShmemSize);
            free(hugebuf);
        }
    }
    s_ShmemSize = newsize;

    return s_Shmem;
}


HEAP LBSM_Shmem_Create(void)
{
    HEAP heap;

    if (s_Shmem)
        CORE_LOG(eLOG_Warning, "Shared memory already exists");

    if (!LBSM_Shmem_Lock())
        return 0;

    s_ShmemSize = 0;
    heap = HEAP_Create(0, 0, getpagesize(), s_ExpandShmem);
    
    LBSM_Shmem_Unlock();
    
    return heap;
}


void LBSM_Shmem_Destroy(HEAP heap)
{
    HEAP_Destroy(heap);
    /* 'expand' has done all local housekeeping */
}


HEAP LBSM_Shmem_Attach(void)
{
    HEAP heap;

    if (s_Shmem)
        CORE_LOG(eLOG_Warning, "Shared memory already attached");

    if ((s_Shmid = shmget(LBSM_SHMEM_KEY, 0, 0)) < 0 ||
        !(s_Shmem = (char *)shmat(s_Shmid, 0, SHM_RDONLY))) {
        heap = 0;
        CORE_LOG_SYS_ERRNO(eLOG_Warning, "Unable to attach shared memory");
    } else {
        heap = HEAP_Attach(s_Shmem);
    }

    return heap;
}


void LBSM_Shmem_Detach(HEAP heap)
{
    HEAP_Detach(heap);
    if (s_Shmem) {
        shmdt(s_Shmem);
        /* Attached heap didn't have smart 'expand', thus the statements */
        s_Shmem = 0;
    }
    s_Shmid = -1;
}
