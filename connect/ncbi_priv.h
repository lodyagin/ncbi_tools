#ifndef CONNECT___NCBI_PRIV__H
#define CONNECT___NCBI_PRIV__H

/*  $Id: ncbi_priv.h,v 6.20 2007/06/25 14:25:36 kazimird Exp $
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
 * Author:  Denis Vakatov
 *
 * File Description:
 *    Private aux. code for the "ncbi_*.[ch]"
 *
 *********************************
 * Random generator seeding support
 *    private global:  g_NCBI_ConnectRandomSeed
 *    macro:           NCBI_CONNECT_SRAND_ADDEND
 * Critical section (basic multi-thread synchronization):
 *    private global:  g_CORE_MT_Lock
 *    macros:          CORE_LOCK_WRITE, CORE_LOCK_READ, CORE_UNLOCK
 * Tracing and logging:
 *    private global:  g_CORE_Log
 *    macros:          CORE_LOG[F](), CORE_DATA[F](), CORE_LOG[F]_ERRNO[_EX]()
 *
 */

#include "ncbi_assert.h"
#include <connect/ncbi_util.h>


#ifdef __cplusplus
extern "C" {
#endif


/******************************************************************************
 *  Random generator seeding support
 */

extern NCBI_XCONNECT_EXPORT int g_NCBI_ConnectRandomSeed;
extern NCBI_XCONNECT_EXPORT int g_NCBI_ConnectSrandAddend();
#define NCBI_CONNECT_SRAND_ADDEND g_NCBI_ConnectSrandAddend()


/******************************************************************************
 *  Multi-Thread SAFETY
 */

/* Always use the following macros and functions to access "g_CORE_MT_Lock",
 * dont access/change it directly!
 */
extern NCBI_XCONNECT_EXPORT MT_LOCK g_CORE_MT_Lock;

#define CORE_LOCK_WRITE  verify( MT_LOCK_Do(g_CORE_MT_Lock, eMT_Lock    ) )
#define CORE_LOCK_READ   verify( MT_LOCK_Do(g_CORE_MT_Lock, eMT_LockRead) )
#define CORE_UNLOCK      verify( MT_LOCK_Do(g_CORE_MT_Lock, eMT_Unlock  ) )



/******************************************************************************
 *  ERROR HANDLING and LOGGING
 */

/* Always use the following macros and functions to access "g_CORE_Log",
 * dont access/change it directly!
 */
extern NCBI_XCONNECT_EXPORT LOG g_CORE_Log;

extern NCBI_XCONNECT_EXPORT const char* g_CORE_Sprintf(const char* fmt, ...)
#ifdef __GNUC__
         __attribute__((format(printf, 1, 2)))
#endif
;

#define CORE_LOG(level, message)  do { \
    if (g_CORE_Log  ||  (level) == eLOG_Fatal) { \
        CORE_LOCK_READ; \
        LOG_WRITE(g_CORE_Log, level, message); \
        CORE_UNLOCK; \
    } \
} while (0)

#define CORE_LOGF(level, fmt_args)  do { \
    if (g_CORE_Log  ||  (level) == eLOG_Fatal) { \
        CORE_LOCK_READ; \
        LOG_WRITE(g_CORE_Log, level, g_CORE_Sprintf fmt_args); \
        CORE_UNLOCK; \
    } \
} while (0)

#ifdef _DEBUG
#  define CORE_TRACE(message)    CORE_LOG(eLOG_Trace, message)
#  define CORE_TRACEF(fmt_args)  CORE_LOGF(eLOG_Trace, fmt_args)
#  define CORE_DEBUG_ARG(arg)    arg
#else
#  define CORE_TRACE(message)    ((void) 0)
#  define CORE_TRACEF(fmt_args)  ((void) 0)
#  define CORE_DEBUG_ARG(arg)    /*arg*/
#endif /*_DEBUG*/

#define CORE_DATA(data, size, message)  do { \
    if ( g_CORE_Log ) { \
        CORE_LOCK_READ; \
        LOG_DATA(g_CORE_Log, data, size, message); \
        CORE_UNLOCK; \
    } \
} while (0)

#define CORE_DATAF(data, size, fmt_args)  do { \
    if ( g_CORE_Log ) { \
        CORE_LOCK_READ; \
        LOG_DATA(g_CORE_Log, data, size, g_CORE_Sprintf fmt_args); \
        CORE_UNLOCK; \
    } \
} while (0)

#define CORE_LOG_ERRNO(level, x_errno, message)  do { \
    if (g_CORE_Log  ||  (level) == eLOG_Fatal) { \
        CORE_LOCK_READ; \
        LOG_WRITE_ERRNO_EX(g_CORE_Log, level, message, x_errno, 0); \
        CORE_UNLOCK; \
    } \
} while (0)

#define CORE_LOGF_ERRNO(level, x_errno, fmt_args)  do { \
    if (g_CORE_Log  ||  (level) == eLOG_Fatal) { \
        CORE_LOCK_READ; \
        LOG_WRITE_ERRNO_EX(g_CORE_Log, level, g_CORE_Sprintf fmt_args, \
                           x_errno, 0); \
        CORE_UNLOCK; \
    } \
} while (0)

#define CORE_LOG_ERRNO_EX(level, x_errno, x_descr, message)  do { \
    if (g_CORE_Log  ||  (level) == eLOG_Fatal) { \
        CORE_LOCK_READ; \
        LOG_WRITE_ERRNO_EX(g_CORE_Log, level, message, x_errno, x_descr); \
        CORE_UNLOCK; \
    } \
} while (0)

#define CORE_LOGF_ERRNO_EX(level, x_errno, x_descr, fmt_args)  do { \
    if (g_CORE_Log  ||  (level) == eLOG_Fatal) { \
        CORE_LOCK_READ; \
        LOG_WRITE_ERRNO_EX(g_CORE_Log, level, g_CORE_Sprintf fmt_args, \
                           x_errno, x_descr); \
        CORE_UNLOCK; \
    } \
} while (0)



/******************************************************************************
 *  REGISTRY
 */

/* Always use the following macros and functions to access "g_CORE_Registry",
 * dont access/change it directly!
 */
extern NCBI_XCONNECT_EXPORT REG g_CORE_Registry;

#define CORE_REG_GET(section, name, value, value_size, def_value) \
    g_CORE_RegistryGET(section, name, value, value_size, def_value)

#define CORE_REG_SET(section, name, value, storage)  do { \
    CORE_LOCK_READ; \
    REG_Set(g_CORE_Registry, section, name, value, storage); \
    CORE_UNLOCK; \
} while (0)


/* (private, to be used exclusively by the above macro CORE_REG_GET) */
extern NCBI_XCONNECT_EXPORT char* g_CORE_RegistryGET
(const char* section,
 const char* name,
 char*       value,
 size_t      value_size,
 const char* def_value);


#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* CONNECT___NCBI_PRIV__H */
