#ifndef CON_SOCK__H
#define CON_SOCK__H

/*  $Id: con_sock.h,v 6.0 1999/04/01 22:09:32 vakatov Exp $
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
*  Implement CONNECTOR for a network socket(based on the NCBI "SOCK").
*
*  See in "connectr.h" for the detailed specification of the underlying
*  connector("CONNECTOR", "SConnectorTag") methods and structures.
*
* --------------------------------------------------------------------------
* $Log: con_sock.h,v $
* Revision 6.0  1999/04/01 22:09:32  vakatov
* Initial revision
*
* ==========================================================================
*/

#include <connectr.h>
#include <connutil.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif


/* This is equivalent to SOCK_CreateConnectorEx(host, port, conn_try, 0,0).
 */
NLM_EXTERN CONNECTOR SOCK_CreateConnector
(const Nlm_Char* host,     /* server:  host */
 Nlm_Uint2       port,     /* server:  service port */
 Nlm_Uint4       conn_try  /* max.number of attempts to establish conn */
 );


/* Create new CONNECTOR structure to handle connection to a socket.
 * Make up to "conn_try" attempts to connect to the "host:port" before
 * giving up.
 * On successful connect, send the first "init_size" bytes from buffer
 * "init_data"(can be NULL -- then send nothing) to the connection.
 * NOTE:  the connector makes(and then uses) its own copy of the "init_data".
 * Return NULL on error.
 */
NLM_EXTERN CONNECTOR SOCK_CreateConnectorEx
(const Nlm_Char* host,      /* server:  host */
 Nlm_Uint2       port,      /* server:  service port */
 Nlm_Uint4       conn_try,  /* max.number of attempts to establish conn */
 const void*     init_data, /* data to send to server on connect */
 Nlm_Uint4       init_size  /* size of the "init_data" buffer */
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

#endif /* CON_SOCK__H */

