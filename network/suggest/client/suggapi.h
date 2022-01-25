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
* File Name: suggapi.h
*
* Author:  Yuri Sadykov
*
* Version Creation Date:   08/14/95
*
* $Revision: 6.0 $
*
* File Description: 
*	Header file for API for Suggest service
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
*/

#ifndef __suggapi_h__
#define __suggapi_h__

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */

#ifndef _NCBI_Seq_
#include <objseq.h>
#endif	/* _NCBI_Seq_ */

#ifndef _suggen_
#include "suggen.h"
#endif	/* _suggen_ */

Boolean SuggestInit PROTO((void));
Boolean SuggestFini PROTO((void));
SeqAnnotPtr SuggestFindIntervals PROTO((SuggestIntervalsPtr pIntervals));

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */

#endif	/* __suggapi_h__ */