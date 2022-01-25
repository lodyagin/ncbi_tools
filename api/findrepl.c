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
* File Name: findrepl.c
*
* Author:  Yuri Sadykov
*
* Version Creation Date:   10/17/95
*
* $Id: findrepl.c,v 6.0 1997/08/25 18:05:38 madden Exp $
* $Revision: 6.0 $
*
* File Description:
*	The implementation of find/replace
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
* $Log: findrepl.c,v $
* Revision 6.0  1997/08/25 18:05:38  madden
* Revision changed to 6.0
*
* Revision 5.3  1997/06/19 18:37:41  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.2  1997/03/17 23:44:39  kans
* added whole_word parameter to FindInEntity and FindInEntityX, and protected
* against multiple ObjMgrAlsoSelects on a single itemID
*
 * Revision 5.1  1996/09/06  20:20:41  kans
 * keeps going even if ObjMgrTypeFind returns NULL (e.g., on OBJ_BIOSEQ_SEG),
 * and adds a case_counts parameter for case sensitive/insensitive searches.
 *
 * Revision 5.0  1996/05/28  13:23:23  ostell
 * Set to revision 5.0
 *
 * Revision 1.7  1996/02/28  04:53:06  ostell
 * fix to prevernt recursion on substring replaces
 *
 * Revision 1.6  1996/02/26  20:24:05  kans
 * replace needs MemCopy instead of StringMove (JO), and set dirty flag
 *
 * Revision 1.5  1996/01/03  23:06:32  ostell
 * support for longer replaces, controlled updating
 *
 * Revision 1.3  1996/01/02  18:40:07  ostell
 * simplified code.
 *
 * Revision 1.2  1996/01/01  00:05:14  kans
 * replaced StringStr with StringISearch to ignore case
 *
 * Revision 1.1  1995/12/31  18:13:14  kans
 * Initial revision
 *
* Revision 1.1.1.1  1995/10/19 18:42:10  sad
* Initial version
*
*/

#include <findrepl.h>
#include <objsub.h>

typedef struct _GatherFindStruct {
    CharPtr	    pchFindStr;
	CharPtr     pchReplStr;
	Int4        flen, replen;
    FindStructFunc  UserFunc;
    Pointer	    pUserData;
    Int4	    iFoundCount;
    Uint2	    iStop;
	AsnIoPtr	aip;          /* AsnIoNull, used for find, and replace if replace < find length */
	AsnExpOptPtr aeop;
	ObjMgrPtr   omp;
	Boolean	    bSelect,	  /* if TRUE, send an ObjMgrSelect message */
		        replace_all,  /* if TRUE, replace all instances found */
				do_replace,   /* signal used when replace string longer than find */
				did_replace,  /* TRUE if any replaces were done */
				substring,    /* TRUE if findstr is a substring of replstr */
				case_counts,  /* if TRUE use StringSearch, else StringISearch */
				whole_word;   /* if TRUE must not have flanking printable characters */
	Boolean		notYetFound;
	Int2 send_update;  /* if 1, send ObjMgrUpdate message on each item, 2 only on input entityID at end */
} GatherFindStruct, PNTR GatherFindStructPtr;

typedef struct _AsnWriteFindStruct {
    GatherFindStructPtr	pGFS;
    GatherContextPtr	pGC;
} AsnWriteFindStruct, PNTR AsnWriteFindStructPtr;

static Boolean
GatherFindCallBack(GatherContextPtr pGContext);

static void
LIBCALLBACK AsnWriteFindCallBack(AsnExpOptStructPtr pAEOS);

static void
LIBCALLBACK AsnWriteReplaceCallBack(AsnExpOptStructPtr pAEOS);

static Uint2 LIBCALLBACK
DefaultUserCallBack(Pointer findstruct);

static Pointer NEAR FindReplMemCopy (Pointer from, AsnReadFunc readfunc, AsnWriteFunc writefunc,
									 GatherFindStructPtr gfsp);

static CharPtr SearchForString (CharPtr str, CharPtr sub, Boolean case_counts, Boolean whole_word)

{
  CharPtr  ptr = NULL;
  CharPtr  tmp;

  if (case_counts) {
    ptr = StringSearch (str, sub);
  } else {
    ptr = StringISearch (str, sub);
  }
  if (ptr == NULL) return NULL;
  if (whole_word) {
    if (ptr > str) {
      tmp = ptr - 1;
      if (! IS_WHITESP (*tmp)) return NULL;
    }
    tmp = ptr + StringLen (sub);
    if (*tmp != '\0' && (! IS_WHITESP (*tmp))) return NULL;
  }
  return ptr;
}

NLM_EXTERN Int4 LIBCALL FindInEntity(Uint2 EntityID, CharPtr findstr, CharPtr replstr, Boolean select_item,
						  Int2 send_update, Boolean case_counts, Boolean whole_word, Boolean replace_all)
{
    GatherFindStruct	GFStruct;

    if (EntityID == 0 || findstr == NULL)
	return -1;

    GFStruct.pchFindStr = findstr;
	GFStruct.flen = StringLen(findstr);
	GFStruct.pchReplStr = replstr;
	GFStruct.replen = StringLen(replstr);
    GFStruct.UserFunc = NULL;
    GFStruct.pUserData = &GFStruct;
    GFStruct.iFoundCount = 0;
    GFStruct.iStop = FS_CONTINUE;
    GFStruct.aip = AsnIoNullOpen();
	GFStruct.aeop = AsnExpOptNew(GFStruct.aip, NULL, NULL, AsnWriteFindCallBack);
	GFStruct.omp = ObjMgrGet();
	GFStruct.bSelect = select_item;
	GFStruct.replace_all = replace_all;
	GFStruct.send_update = send_update;
	GFStruct.do_replace = FALSE;
	GFStruct.did_replace = FALSE;
	GFStruct.substring = FALSE;
	GFStruct.case_counts = case_counts;
	GFStruct.whole_word = whole_word;
	if (GFStruct.flen <= GFStruct.replen)  /* possible infinite loop */
	{
		if ((SearchForString(replstr, findstr, case_counts, whole_word)) != NULL)
			GFStruct.substring = TRUE;
	}

    GatherEntity(EntityID, &GFStruct, GatherFindCallBack, NULL);

	AsnIoClose(GFStruct.aip);

	if ((GFStruct.did_replace) && (GFStruct.send_update == UPDATE_ONCE))
		ObjMgrSendMsg(OM_MSG_UPDATE, EntityID, 0,0); 

    return GFStruct.iFoundCount;
}

NLM_EXTERN Int4 LIBCALL FindInEntityX(Uint2 EntityID, CharPtr findstr, CharPtr replstr, FindStructFunc userfunc,
	      Pointer userdata, Boolean select_item, Int2 send_update, Boolean case_counts, Boolean whole_word)
{
    GatherFindStruct	GFStruct;

    if (EntityID == 0 || findstr == NULL || userfunc == NULL)
	return -1;

    GFStruct.pchFindStr = findstr;
	GFStruct.pchReplStr = replstr;
    GFStruct.UserFunc = userfunc;
    GFStruct.pUserData = userdata;
    GFStruct.iFoundCount = 0;
    GFStruct.iStop = FS_CONTINUE;
    GFStruct.aip = AsnIoNullOpen();
	GFStruct.aeop = AsnExpOptNew(GFStruct.aip, NULL, NULL, AsnWriteFindCallBack);
	GFStruct.omp = ObjMgrGet();
	GFStruct.bSelect = select_item;
	GFStruct.replace_all = FALSE;
	GFStruct.send_update = send_update;
	GFStruct.do_replace = FALSE;
	GFStruct.did_replace = FALSE;
	GFStruct.substring = FALSE;
	GFStruct.case_counts = case_counts;
	GFStruct.whole_word = whole_word;
	if (GFStruct.flen <= GFStruct.replen)  /* possible infinite loop */
	{
		if ((SearchForString(replstr, findstr, case_counts, whole_word)) != NULL)
			GFStruct.substring = TRUE;
	}

    GatherEntity(EntityID, &GFStruct, GatherFindCallBack, NULL);
	AsnIoClose(GFStruct.aip);

	if ((GFStruct.did_replace) && (GFStruct.send_update == UPDATE_ONCE))
		ObjMgrSendMsg(OM_MSG_UPDATE, EntityID, 0,0); 

    return GFStruct.iFoundCount;
}

static Boolean NEAR FindReplFunc ( ObjMgrTypePtr pOMType, Pointer oldptr, GatherFindStructPtr gfsp)
{
	Boolean done = FALSE, replaced = FALSE;
	Pointer newptr;

	gfsp->notYetFound = TRUE;
	while (! done)
	{
		(pOMType->asnwrite)(oldptr, gfsp->aip, NULL);   /* look for strings  */
		                                                             /* replace if new shorter than old */
		if (gfsp->do_replace)  /* have to replace by copy */
		{
			 newptr = FindReplMemCopy (oldptr,(pOMType->asnread), (pOMType->asnwrite), gfsp);
			 if (newptr != NULL)
			 {
				 GatherOverWrite(oldptr, newptr, pOMType->datatype);
				 (*(pOMType->freefunc))(newptr); /* free it */
				 replaced = TRUE;
				 gfsp->did_replace = TRUE;
			 }
			 gfsp->do_replace = FALSE;
			 if (gfsp->iStop == FS_STOP)
				 done = TRUE;
			 if (gfsp->substring)
				 done = TRUE;   /* avoid recursion */
		}
		else
			done = TRUE;
	}

	return replaced;
}

static Boolean
GatherFindCallBack(GatherContextPtr pGContext)
{
    ObjMgrTypePtr	pOMType;
	GatherFindStructPtr gfsp;
	Boolean replaced = FALSE;
    /*
    AsnWriteFindStruct	pAWFS;
    */

	gfsp = (GatherFindStructPtr)(pGContext->userdata);
	if (gfsp->iStop == FS_STOP)
		return FALSE;
	gfsp->aeop->user_data = pGContext;

    pOMType = ObjMgrTypeFind(gfsp->omp, pGContext->thistype, NULL, NULL);
	if (pOMType == NULL)
		return TRUE;

    switch (pGContext->thistype) {

		case OBJ_BIOSEQ: {
			BioseqPtr   pBioseq;
			Pointer	    pSeq_ext;
			SeqHistPtr  pSeqHist;
			ValNodePtr  pDescr;
			SeqAnnotPtr pAnnot;

			pBioseq = (BioseqPtr)pGContext->thisitem;
	
			pSeq_ext = pBioseq->seq_ext;
			pBioseq->seq_ext = NULL;
			pSeqHist = pBioseq->hist;
			pBioseq->hist = NULL;
			pDescr = pBioseq->descr;
			pBioseq->descr = NULL;
			pAnnot = pBioseq->annot;
			pBioseq->annot = NULL;

			replaced = FindReplFunc (pOMType, pGContext->thisitem, gfsp);
	
			pBioseq->seq_ext = pSeq_ext;
			pBioseq->hist = pSeqHist;
			pBioseq->descr = pDescr;
			pBioseq->annot = pAnnot;
	
			break;
		    }

		case OBJ_BIOSEQSET: {
			BioseqSetPtr	pBS;
			ValNodePtr	pDescr;
			SeqAnnotPtr	pSeqAnnot;
			SeqEntryPtr	pSeqEntry;

			pBS = (BioseqSetPtr)pGContext->thisitem;

			pDescr = pBS->descr;
			pBS->descr = NULL;
			pSeqAnnot = pBS->annot;
			pBS->annot = NULL;
			pSeqEntry = pBS->seq_set;
			pBS->seq_set = NULL;

			replaced = FindReplFunc (pOMType, pGContext->thisitem, gfsp);

			pBS->descr = pDescr ;
			pBS->annot = pSeqAnnot;
			pBS->seq_set = pSeqEntry;

			break;
		   }

		case OBJ_SEQANNOT: {
			SeqAnnotPtr	pSeqAnnot;
			Pointer data;
	
			pSeqAnnot = (SeqAnnotPtr)pGContext->thisitem;
			data = pSeqAnnot->data;
			pSeqAnnot->data = NULL;

			replaced = FindReplFunc (pOMType, pGContext->thisitem, gfsp);
	
			pSeqAnnot->data = data;
			break;
			}

	    case OBJ_SUBMIT_BLOCK:
		case OBJ_SEQHIST:
		case OBJ_SEQID:
		case OBJ_SEQLOC:
		case OBJ_SEQCODE:
		case OBJ_PUB:
		case OBJ_SEQFEAT:
		case OBJ_SEQALIGN:
		case OBJ_SEQGRAPH:
		case OBJ_SEQDESC:
		case OBJ_BIOSEQ_MAPFEAT:
		case OBJ_BIOSEQ_SEG:
		case OBJ_GENETIC_CODE:
			replaced = FindReplFunc (pOMType, pGContext->thisitem, gfsp);
			break;

		case OBJ_SEQSUB_CONTACT:
		case OBJ_SEQSUB_CIT:
		case OBJ_ANNOTDESC:
		case OBJ_MEDLINE_ENTRY:
		case OBJ_SEQFEAT_CIT:
			if (pGContext->indent == 0)  /* these are normally included in other things */
				replaced = FindReplFunc (pOMType, pGContext->thisitem, gfsp);
			break;

		case OBJ_SEQENTRY:    /* Bioseq or BioseqSet will get this */
		default:
			break;
	}

	if (replaced)
		ObjMgrSetDirtyFlag(pGContext->entityID, TRUE);
	if ((gfsp->send_update == UPDATE_EACH) && replaced)  /* only true with a copy type replace */
		ObjMgrSendMsg(OM_MSG_UPDATE, pGContext->entityID, pGContext->itemID, pGContext->thistype); 
	if (gfsp->bSelect && replaced)		/* only true with a copy type replace */
		ObjMgrAlsoSelect(pGContext->entityID, pGContext->itemID, pGContext->thistype, 0, NULL);

	if (gfsp->iStop == FS_STOP)
		return FALSE;
	else
		return TRUE;
}

static void
LIBCALLBACK AsnWriteFindCallBack(AsnExpOptStructPtr pAEOS)
{
    GatherContextPtr    pGC;
    GatherFindStructPtr pGFS;
	CharPtr foundit, tmp;
	Boolean replaceit = FALSE, foundone = FALSE, done, replaced = FALSE;

    pGC = (GatherContextPtr)pAEOS->data;
	pGFS = (GatherFindStructPtr)pGC->userdata;
    if ((pGFS->iStop == FS_STOP) || (pGFS->do_replace))
		return;


    if (ISA_STRINGTYPE(AsnFindBaseIsa(pAEOS->atp))) {
	CharPtr		    pchSource;
	CharPtr		    pchFind, pchReplace, ptr, lastptr;
	Int2			iStop;

	pchSource = (CharPtr)pAEOS->dvp->ptrvalue;
	pchFind = pGFS->pchFindStr;
	pchReplace = pGFS->pchReplStr;
	done = FALSE;
	ptr = pchSource;
	lastptr = ptr;
	while (((foundit = SearchForString(ptr, pchFind, pGFS->case_counts, pGFS->whole_word)) != NULL) && (! done)) {
	    FindStruct	FS;

		foundone = TRUE;
 		if (pGFS->replace_all)
			replaceit = TRUE;
	    (pGFS->iFoundCount)++;

		if (pGFS->UserFunc != NULL)
		{
		    FS.gcp = pGC;
		    FS.udp = ((GatherFindStructPtr)pGC->userdata)->pUserData;
		    FS.findstr = pchFind;
			FS.replstr = pchReplace;
	/*	    FS.select_item = pGFS->bSelect; */

		    iStop = (pGFS->UserFunc)(&FS);
			pGFS->iStop = (iStop & FS_FLAG_CONTINUE);   /* only set the continue bit */
			if (iStop & FS_FLAG_REPLACE)
				replaceit = TRUE;
		}

		if (replaceit)
		{
			if (pGFS->replen <= pGFS->flen)   /* replace in_situ */
			{
				tmp = MemCopy(foundit, pchReplace, pGFS->replen);
				tmp += pGFS->replen;
				foundit += pGFS->flen;
				ptr = foundit;
				tmp = StringMove(tmp, foundit);
				replaced = TRUE;
				pGFS->did_replace = TRUE;
			}
			else
			{
				pGFS->do_replace = TRUE;   /* signal to replace by copy */
				return;
			}
			
		}
		else
			ptr = foundit + pGFS->flen;

		if (pGFS->iStop == FS_STOP)
			done = TRUE;

	}
    }

	if (pGFS->bSelect && foundone && pGFS->notYetFound) {
		ObjMgrAlsoSelect(pGC->entityID, pGC->itemID, pGC->thistype, 0, NULL);
		pGFS->notYetFound = FALSE;
	}
	if ((pGFS->send_update == UPDATE_EACH) && replaced)
		ObjMgrSendMsg(OM_MSG_UPDATE, pGC->entityID, pGC->itemID, pGC->thistype); 
}

static void
LIBCALLBACK AsnWriteReplaceCallBack(AsnExpOptStructPtr pAEOS)
{
    GatherFindStructPtr pGFS;
	CharPtr foundit, tmp, tmp2, newstring;
	Int4 diff;

	pGFS = (GatherFindStructPtr)pAEOS->data;
    if (! pGFS->do_replace)
		return;

    if (ISA_STRINGTYPE(AsnFindBaseIsa(pAEOS->atp))) {
	CharPtr		    pchSource;
	CharPtr		    pchFind, pchReplace;

	pchSource = (CharPtr)pAEOS->dvp->ptrvalue;
	pchFind = pGFS->pchFindStr;
	pchReplace = pGFS->pchReplStr;
	if ((foundit = SearchForString(pchSource, pchFind, pGFS->case_counts, pGFS->whole_word)) == NULL)
		return;

	diff = pGFS->replen - pGFS->flen;
	
	newstring = MemNew(StringLen(pchSource) + diff + 1);
	tmp = pchSource;
	tmp2 = newstring;
	while (tmp != foundit)
	{
		*tmp2 = *tmp;
		tmp++; tmp2++;
	}

	tmp2 = MemCopy(tmp2, pchReplace, pGFS->replen);
	tmp2 += pGFS->replen;
	tmp += pGFS->flen;
	tmp2 = StringMove(tmp2, tmp);

	pGFS->do_replace = FALSE;
	MemFree(pAEOS->dvp->ptrvalue);
	pAEOS->dvp->ptrvalue = newstring;
	}

	return;
}


static Uint2 LIBCALLBACK
DefaultUserCallBack(Pointer findstruct)
{
    FindStructPtr   pFS;
	GatherFindStructPtr gfsp;
	GatherContextPtr gcp;

    if (findstruct == NULL)
	return FS_CONTINUE;

    pFS = (FindStructPtr)findstruct;
	gfsp = (GatherFindStructPtr)(pFS->udp);
	gcp = pFS->gcp;

	if (gfsp->replace_all)
	{

	}
    if (*(Boolean*)pFS->udp) {
	GatherContextPtr pGC;

	pGC= pFS->gcp;
	ObjMgrAlsoSelect(pGC->entityID, pGC->itemID, pGC->thistype, 0, NULL);
    }

    return FS_CONTINUE;
}

/*****************************************************************
*
*   This function is a lift of AsnIoMemCopy, modified to support
*      string substitution
*
******************************************************************/
static Pointer NEAR FindReplMemCopy (Pointer from, AsnReadFunc readfunc, AsnWriteFunc writefunc, GatherFindStructPtr gfsp)
{
	Pointer res;
	AsnIoBSPtr aibp;
	ByteStorePtr bsp;
	AsnExpOptPtr aeop;
	
	if ((from == NULL) || (readfunc == NULL) || (writefunc == NULL))
		return NULL;

	bsp = BSNew(5000);
	aibp = AsnIoBSOpen("wb", bsp);
	if (aibp == NULL)
		return NULL;

    if (! (*writefunc)(from, aibp->aip, NULL))
	{
		AsnIoBSClose(aibp);
		BSFree(bsp);
		return NULL;
	}

	AsnIoBSClose(aibp);

	aibp = AsnIoBSOpen("rb", bsp);
	aeop = AsnExpOptNew(aibp->aip, NULL, (Pointer)gfsp, AsnWriteReplaceCallBack); 
	res = (*readfunc)(aibp->aip, NULL);

	AsnIoBSClose(aibp);
	BSFree(bsp);

	return res;
}
