#ifndef NCBICLI__H
#define NCBICLI__H

/*  $RCSfile: ncbicli.h,v $  $Revision: 4.4 $  $Date: 1998/12/15 17:24:58 $
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
*   API to open(using HTTPD-based dispatcher) and handle socket-based
*   connection to an NCBI service
*
* --------------------------------------------------------------------------
* $Log: ncbicli.h,v $
* Revision 4.4  1998/12/15 17:24:58  vakatov
* +NIC_CERN_PROXY flag(to indicate the use of CERN-like non-transparent proxy)
*
* Revision 4.3  1998/09/08 17:59:05  vakatov
* Added WWW/Firewall network interface
*
* Revision 4.2  1998/04/02 20:21:10  vakatov
* Added possibility to printout service reply headers
*
* Revision 4.1  1998/03/30 17:50:09  vakatov
* Ingrafted to the main NCBI CVS tree
*
* ==========================================================================
*/

#include <ncbisock.h>

#define NIC_ERROR   888

#define NIC                Nlm_NIC
#define ENIC_Agent         Nlm_ENIC_Agent
#define NIC_GetService     Nlm_NIC_GetService
#define NIC_GetSOCK        Nlm_NIC_GetSOCK
#define NIC_CloseService   Nlm_NIC_CloseService


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct Nlm_NICtag;              /* internal storage  */
typedef struct Nlm_NICtag *NIC; /* handle */


/* Open read/write connection to the service specified by "service_name"
 * using dispatcher located at URL "http://disp_host:disp_port/disp_path".
 * If both "timeout_sec" and "timeout_usec" are equal to zero then wait
 * ad infinitum; otherwise, fail if the specified timeout expires.
 * Return handle to the opened connection; NULL on error.
 */

typedef enum {
  eNIC_WWWClient = 0,
  eNIC_WWWDirect,

  eNIC_Default /* use default agent */
} ENIC_Agent;

/* the "flags" arg for NIC_GetService()
 */
/* printout reply from dispatcher */
#define NIC_DEBUG_PRINTOUT 0x1
/* use an NCBI firewall daemon to make a connection */
#define NIC_FIREWALL       0x2
/* non-transparent CERN-like proxy */
#define NIC_CERN_PROXY     0x4


NLM_EXTERN NIC NIC_GetService
(const Nlm_Char  *service_name,
 const Nlm_Char  *disp_host,
 Nlm_Uint2        disp_port,
 const Nlm_Char  *disp_path,
 const STimeout  *timeout,
 ENIC_Agent       client_agent,
 const Char      *client_host,
 const ByteStore *service_query,
 Uint4            flags
 );

/* Get SOCK structure associated with the connection
 */ 
NLM_EXTERN SOCK NIC_GetSOCK
(NIC nic
 );

/* Close the connection, destroy relevant internal data
 */
NLM_EXTERN Nlm_Boolean NIC_CloseService
(NIC nic
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

#endif /* NCBI_CLI__H */