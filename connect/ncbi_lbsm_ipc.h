#ifndef NCBI_LBSM_IPC__H
#define NCBI_LBSM_IPC__H

/*  $Id: ncbi_lbsm_ipc.h,v 6.15 2001/07/05 16:58:12 lavr Exp $
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
 * Author:  Anton Lavrentiev, Denis Vakatov
 *
 * File Description:
 *    Implementation of the LBSM client-server data exchange API
 *    with the use of SYSV IPC (shared memory and semaphores)
 *
 * --------------------------------------------------------------------------
 * $Log: ncbi_lbsm_ipc.h,v $
 * Revision 6.15  2001/07/05 16:58:12  lavr
 * Comments updated to show calls for use by daemon only
 *
 * Revision 6.14  2001/07/03 20:46:45  lavr
 * Added function: LBSM_Shmem_Retry() to delete hanging locks.
 * Modified function: LBSM_Shmem_Create() - timeout argument added
 *
 * Revision 6.13  2001/06/04 20:55:03  lavr
 * Widened the trick to get rid of 'union semun' definition on IRIX
 * (the former trick has a side effect of not defining hton*() and ntoh*())
 *
 * Revision 6.12  2001/05/14 14:01:37  lavr
 * Fixed: Linux compilation fails due to 'union semun' undefined
 *
 * Revision 6.11  2001/05/11 15:29:39  lavr
 * Added additional condition whether to define union semun with GNU Library
 *
 * Revision 6.10  2001/03/19 23:07:51  lavr
 * Log typo corrected :-(
 *
 * Revision 6.9  2001/03/19 23:02:27  lavr
 * LBSM_UnLBSMD now accepts one argument and returns a value
 *
 * Revision 6.8  2000/12/06 22:08:50  lavr
 * Added a kludge to avoid double definition of 'union semun' on FreeBSD
 *
 * Revision 6.7  2000/10/20 17:11:16  lavr
 * Access rights to shared memory changed: 'other' now may write
 * (to eliminate problem with UID/GID mismatch for LBSMD and CGIs)
 *
 * Revision 6.6  2000/10/06 17:38:23  lavr
 * Dirty trick undone - didn't come out anyway
 *
 * Revision 6.5  2000/10/06 17:35:25  lavr
 * Dirty trick for elimination of declaration of 'union semun',
 * which is used in a call to 'semctl'
 *
 * Revision 6.4  2000/10/05 21:29:40  lavr
 * ncbiconf.h removed
 *
 * Revision 6.3  2000/05/18 14:11:35  lavr
 * Yet another cosmetic update
 *
 * Revision 6.2  2000/05/15 19:06:09  lavr
 * Tiny little cosmetic change
 *
 * Revision 6.1  2000/05/12 19:19:16  lavr
 * First working revision
 *
 * ==========================================================================
 */

#include <connect/ncbi_heapmgr.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef HAVE_SEMUN
/* This sequence of defines causes 'union semun' be undefined on IRIX */
#  ifdef _XOPEN_SOURCE
#    define _XOPEN_SOURCE_SAVE _XOPEN_SOURCE
#    undef  _XOPEN_SOURCE
#  endif
#  define _XOPEN_SOURCE 1
#endif
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#ifndef HAVE_SEMUN
#  undef _XOPEN_SOURCE
#  ifdef _XOPEN_SOURCE_SAVE
#    define _XOPEN_SOURCE _XOPEN_SOURCE_SAVE
#  endif
#endif


#ifdef __cplusplus
extern "C" {
#endif


/* Keys to access the LBSM shared memory and semaphores
 */
#define LBSM_SHMEM_KEY   0x130DFB2
#define LBSM_MUTEX_KEY   0x12CC359

#define LBSM_SHARED_PROT ( S_IRUSR | S_IWUSR | \
                           S_IRGRP | S_IWGRP | \
                           S_IROTH | S_IWOTH )

#if !defined(HAVE_SEMUN) && !defined(__FreeBSD__) && \
    (!defined(__GNU_LIBRARY__) || defined(_SEM_SEMUN_UNDEFINED))
union semun {
    int              val;
    struct semid_ds* buf;
    unsigned short*  array;
};
#endif


/* Check (and lock, if "check_n_lock" is TRUE) an instance of LBSMD.
 * Return value:
 * -1 means the mutex could not be acquired;
 *  0 means the mutex was vacant, and operation went through;
 *  1 means the lock operation failed (mutex was already locked, or whatever).
 * In cases of non-zero return code, "errno" must be analyzed to
 * figure out the problem (if any). Locking is reserved for use by daemon
 * only and undone automatically on program exit.
 *
 * This must be the first call prior to any shared resources use, as it
 * sets up an internal semaphore ID for all locking/unlocking code.
 */
int LBSM_LBSMD(int/*bool*/ check_n_lock);


/* Remove LBSMD-specific internal semaphore ID (LBSM_LBSMD undo).
 * This call is mostly for use by daemon (called with argument != 0),
 * but if called with its argument zero then could also be used by
 * clients to obtain PID of the daemon running (LBSM_LBSMD yet must be
 * called before). Return value (if differs from 0) is most likely the PID
 * of running LBSM daemon (actually, the PID of LBSM shared memory creator).
 */
pid_t LBSM_UnLBSMD(int/*bool*/ undaemon);


/* Create shared memory based LBSM heap to work with.
 * Designed for use solely by LBSM daemon.
 */
HEAP LBSM_Shmem_Create(int timeout);


/* Destroy the shared memory based LBSM heap (created by LBSM_Shmem_Create)
 */
void LBSM_Shmem_Destroy(HEAP heap);


/* Attach to a shared memory based LBSM heap
 */
HEAP LBSM_Shmem_Attach(void);


/* Detach the shared memory based LBSM heap (attached by LBSM_Shmem_Attach)
 */
void LBSM_Shmem_Detach(HEAP heap);


/* Return a snapshot of LBSMD shared memory (must be later freed by caller).
 * Shared memory gets firstly attached and lastly detached by this call.
 */
HEAP LBSM_Shmem_Copy(void);


/* Lock the shared memory for exclusive use
 */
int/*bool*/ LBSM_Shmem_Lock(void);


/* Modification of LBSM_Shmem_Lock, which will lock the shared memory even
 * if existing lock is pending at least timeout seconds by removing the lock.
 * This call is intended for exlcusive use by LBSM daemon.
 */
int/*bool*/ LBSM_Shmem_Retry(int timeout);


/* Unlock the shared memory
 */
int/*bool*/ LBSM_Shmem_Unlock(void);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* NCBI_LBSM_IPC__H */
