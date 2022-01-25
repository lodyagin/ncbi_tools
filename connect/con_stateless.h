#ifndef CON_STATELESS__H
#define CON_STATELESS__H

/*  $Id: con_stateless.h,v 6.1 1999/07/21 21:43:27 vakatov Exp $
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
*   Implement a CONNECTOR to perform a stateless connection to NCBI services.
*
*  NOTE: See in "connectr.h" for the detailed specification of the underlying
*        connector("CONNECTOR", "SConnectorTag") methods and structures.
*
* --------------------------------------------------------------------------
* $Log: con_stateless.h,v $
* Revision 6.1  1999/07/21 21:43:27  vakatov
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


/* Create new CONNECTOR structure to connect to an NCBI service using
 * NCBI HTTPD-based dispatcher.
 *
 * The HTTP request is sent to <host>:<port>, and its HTTP header
 * has the following fixed format:
 *
 *    POST <path>?service=<service_name>&
 *                address=<client_hostname>&
 *                platform=<client_platform> HTTP/1.0\n
 *    User-Agent: NCBIDirClientHTTP from <client_hostname>\n
 *    Content-Type: x-ncbi-data/x-unknown\n
 *    Content-Length: <accumulated_data_length>\r\n\r\n
 *
 * where the following is generated automagically:
 *    <client_hostname> := name of the local host,     see Nlm_GetHostName()
 *    <client_platform> := name of the local platform, see Nlm_PlatformName()
 *    <accumulated_data_length> := calculated length of the accumulated data
 *
 * NOTE:  This implementation is based on the URLC_*** API(see in
 * "con_url.[ch]"), with the following flags always set:
 *     URLC_AUTO_RECONNECT, URLC_SURE_FLUSH;
 * and the URLC_HTTP_HEADER is unset.
 */

NLM_EXTERN CONNECTOR STATELESS_CreateConnector
(const char* service_name  /* can contain alpha-numeric and '_' only */
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

#endif /* CON_STATELESS__H */

