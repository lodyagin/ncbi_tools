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
* File Name: findrepl.h
*
* Author:  Yuri Sadykov
*
* Version Creation Date:   10/17/95
*
* $Id: findrepl.h,v 6.0 1997/08/25 18:05:40 madden Exp $
* $Revision: 6.0 $
*
* File Description:
*	Interface to find/replace.
*
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
* -------------------------
* $Log: findrepl.h,v $
* Revision 6.0  1997/08/25 18:05:40  madden
* Revision changed to 6.0
*
* Revision 5.3  1997/06/19 18:37:43  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.2  1997/03/17 23:44:39  kans
* added whole_word parameter to FindInEntity and FindInEntityX
*
 * Revision 5.1  1996/09/06  20:20:41  kans
 * keeps going even if ObjMgrTypeFind returns NULL (e.g., on OBJ_BIOSEQ_SEG),
 * and adds a case_counts parameter for case sensitive/insensitive searches.
 *
 * Revision 5.0  1996/05/28  13:23:23  ostell
 * Set to revision 5.0
 *
 * Revision 1.3  1996/01/03  23:06:32  ostell
 * support for longer replaces, controlled updating
 *
 * Revision 1.1  1995/12/31  18:13:14  kans
 * Initial revision
 *
* Revision 1.1.1.1  1995/10/19 18:42:10  sad
* Initial version
*
*/

#ifndef __FINDREPL_H__
#define __FINDREPL_H__

#include <ncbi.h>
#include <gather.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */


typedef Uint2 (LIBCALLBACK *FindStructFunc) PROTO ((Pointer findstruct));

#define FS_FLAG_CONTINUE    1	/* keep looking */
#define FS_FLAG_REPLACE     2   /* replace it? */

#define FS_CONTINUE	    1
#define FS_STOP		    0
#define FS_REPLACE_CONTINUE 3
#define FS_REPLACE_STOP     2

typedef struct _FindStruct {
    GatherContextPtr	gcp;
    Pointer		udp;	    /* Data supplied by user */
    CharPtr		findstr;    /* String to find */
	CharPtr     replstr;    /* String to replace with */
} FindStruct, PNTR FindStructPtr;

/*** defines for send_update parameter below ***/

#define UPDATE_NEVER 0		/* never send ObjMgrMessage (OM_MSG_UPDATE... */
#define UPDATE_EACH 1		/* send it on each replace */
#define UPDATE_ONCE 2		/* send once for whole entityID, if any replacements occur */

NLM_EXTERN Int4 LIBCALL FindInEntity(Uint2 EntityID, CharPtr findstr, CharPtr replstr, Boolean select_item,
						  Int2 send_update, Boolean case_counts, Boolean whole_word, Boolean replace_all);

NLM_EXTERN Int4 LIBCALL FindInEntityX(Uint2 EntityID, CharPtr findstr, CharPtr replstr, FindStructFunc userfunc,
	      Pointer userdata, Boolean select_item, Int2 send_update, Boolean case_counts, Boolean whole_word);

#ifdef __cplusplus
extern "C" }
#endif	/* __cplusplus */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif	/* __FINDREPL_H__ */