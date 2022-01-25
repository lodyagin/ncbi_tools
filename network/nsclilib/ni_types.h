/*      
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
* File Name:    ni_types.h
*
* Author:       Beatty, Gish
*
* Version Creation Date:        1/1/92
*
* $Revision: 6.8 $
*
* File Description: 
*
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: ni_types.h,v $
* Revision 6.8  2012/02/19 03:45:25  lavr
* Cleanup of obsolete features
*
* Revision 6.7  2004/11/19 14:11:03  lavr
* Reinstate OBSOLETED eNII_ constants which still may be in use in some DEAD code in the toolkit
*
* Revision 6.5  2002/08/08 01:52:28  lavr
* Default dispatcher set to SERVICE
*
* Revision 6.4  2001/02/21 22:09:27  lavr
* SERVICE connector included
*
* Revision 6.3  1998/09/08 17:59:07  vakatov
* Added WWW/Firewall network interface
*
* Revision 6.2  1998/05/05 22:45:39  vakatov
* Added "eNII_Debug" network interface
*
* Revision 6.1  1998/03/30 17:56:06  vakatov
* Added ENIInterface enumerator definition and added "interface" field
* to the NI_Dispatcher struct
*
* Revision 5.3  1997/01/28 21:24:33  epstein
* move NodePtr definition to ncbimisc.h
*
* Revision 5.2  1996/06/28  17:14:39  epstein
* add job-penalty
*
* Revision 5.1  1996/06/27  17:18:17  epstein
* add load-balancing
*
* Revision 4.2  1996/04/29  15:29:10  epstein
* add disp to NI_HandPtr so that service-handle can encapsulate greater context
*
* Revision 4.1  1995/11/27  20:59:31  epstein
* add client support for direct-connection services
* 
* Revision 1.19  1995/05/24  12:09:04  epstein
* add support for tracking of how many times a client IP has used a service within a time interval
*/

#ifndef _NI_TYPES_
#define _NI_TYPES_

#include <ncbi.h>
#include <asn.h>

#define NI_Handle       MHandle    /* for API use */
#define NI_HandPtr      MHandPtr   /* for API use */


typedef struct MHandle {
        CharPtr         hostname;       /* name of peer machine */
        AsnIoPtr        raip;           /* ASNtool IO read pointer */
        AsnIoPtr        waip;           /* ASNtool IO write pointer */
        VoidPtr         extra_proc_info; /* extra processing info, used externally */
        struct NI_Dispatcher PNTR disp;
} MHandle, PNTR MHandPtr;


/* The available connection interfaces
 */
typedef enum {
  /* Refer to "s_NII" in "ni_lib_.c" when changing the enumerator ordering
   * or adding new interfaces */
  eNII_Dispatcher = 0,  /* old-fashioned NCBI dispatched-based connection  | OBSOLETE */
  eNII_WWW,             /* WWW-based connection                            | OBSOLETE */
  eNII_WWWFirewall,     /* eNII_WWW + pass through the NCBI firewall daemon| OBSOLETE */
  eNII_WWWDirect,       /* WWW-based stateless connection                  | OBSOLETE */
  eNII_Service,         /* SERVICE-based connection                         */
  eNII_Debug,           /* direct client-server connection                  */

  /* FEATURE: add new interfaces *above* this point(i.e. above eNII_Default),
   * so that eNII_Default be equal to the number of avail. interfaces */
  eNII_Default      /* let program try environment and config files */

/* NII_DEFAULT will be used if user did not explicitly specify the interface
 * and if application failed to find it in environment and config files */
#define NII_DEFAULT eNII_Service
} ENIInterface;


typedef struct NI_Dispatcher {
        ENIInterface     interface;
        Int2             referenceCount;     /* # of services connected via this dispatcher */
        CharPtr          adminInfo;          /* info. regarding administrator */
        CharPtr          motd;               /* message of the day for user */
} NI_Dispatcher, PNTR NI_DispatcherPtr;


#endif
