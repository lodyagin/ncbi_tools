/*   tax0.h
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
* File Name:  tax0.h
*
* Author:  Vladimir Soussov
*
* Version Creation Date:   01/22/96
*
* $Revision: 6.0 $
*
* File Description: 
*       API for Taxonomy Archive service
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*
*
* RCS Modification History:
* $Log: tax0.h,v $
* Revision 6.0  1997/08/25 18:41:23  madden
* Revision changed to 6.0
*
* Revision 5.1  1996/07/03 18:48:42  epstein
* eliminate double-quotes in favor of angle-brackets
*
 * Revision 5.0  1996/05/28  14:15:13  ostell
 * Set to revision 5.0
 *
 * Revision 1.1  1996/03/06  17:05:19  soussov
 * Initial revision
 *
*/

#include <ncbi.h>
#include <asn.h>
#include <objfeat.h>
#include <objtaxc0.h>

Boolean Tax0Init PROTO((void));
Boolean Tax0Fini PROTO((void));

TaxonIdListPtr Tax0GetTaxId PROTO((TaxonNamePtr tnp));
TaxonIdListPtr Tax0GetChildren PROTO((Int4 id_tax));
TaxonIdListPtr Tax0GetParents PROTO((Int4 id_tax));
OrgRefPtr Tax0GetRef PROTO((Int4 id_tax));
TaxCompleteListPtr Tax0GetComplete PROTO((TaxonIdNamePtr tinp));

/* # of retries to get a server */
#define TAXARCH_SERV_RETRIES 2

#define MAXIDLIST 1   /* was 50 */

