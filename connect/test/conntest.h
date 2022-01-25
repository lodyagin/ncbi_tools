#ifndef CONNTEST__H
#define CONNTEST__H

/*  $Id: conntest.h,v 6.0 1999/03/25 23:05:00 vakatov Exp $
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
*   Test suite for the NCBI connector(see in "connectn.[ch]", "connectr.h")
*
* --------------------------------------------------------------------------
* $Log: conntest.h,v $
* Revision 6.0  1999/03/25 23:05:00  vakatov
* Initial revision
*
* ==========================================================================
*/

#include <connectr.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif


/* Create connection based on the passed "connector" and then:
 * 1) write some data to the connection and expect them the same data come
 *    back from the connection;
 * 2) after reading back all the "bounced" data, read any extra data
 *    coming from the connection and write the extra data to the
 *    "extra_data_file" output file(if non-NULL).
 * 3) close the connection
 */

#define TESTCONN_ALL 0x0
#define TESTCONN_SINGLE_BOUNCE_PRINT 0x1
#define TESTCONN_MULTI_BOUNCE_PRINT  0x2
#define TESTCONN_SINGLE_BOUNCE_CHECK 0x4

NLM_EXTERN void Ncbi_TestConnector
(CONNECTOR       connector,  /* [in] connector handle */
 const STimeout* timeout,    /* [in] timeout for all i/o */
 FILE*           log_file,   /* [in] log file */
 Nlm_Uint4       flags       /* [in] tests to run */
 );


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* CONNTEST__H */
