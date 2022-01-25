/*   $Id: samutil.c,v 1.10 2000/01/24 20:54:34 vakatov Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  $Id: samutil.c,v 1.10 2000/01/24 20:54:34 vakatov Exp $
*
* Author:  Lewis Geer
*
* Version Creation Date:   8/12/99
*
* $Revision: 1.10 $
*
* File Description: Utility functions for AlignIds and SeqAlignLocs
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: samutil.c,v $
* Revision 1.10  2000/01/24 20:54:34  vakatov
* SAM_ViewString::  made #define to fix for the DLL build on PC
*
* Revision 1.9  2000/01/24 16:11:13  lewisg
* speed up seqid comparison in color manager, make fast windows version of SetColor()
*
* Revision 1.8  1999/12/11 01:30:34  lewisg
* fix bugs with sharing colors between ddv and cn3d
*
* Revision 1.7  1999/12/03 23:17:24  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 1.6  1999/11/24 21:24:30  vakatov
* Fixed for the C++ and/or MSVC DLL compilation
*
* Revision 1.5  1999/10/05 23:18:15  lewisg
* add ddv and udv to cn3d with memory management
*
* Revision 1.4  1999/09/27 17:53:08  kans
* seqalign entityID/itemID/itemtype now in GatherIndex substructure
*
* Revision 1.3  1999/09/27 17:49:12  lewisg
* fix denseseg constructor, bug in valnode loops, add SAM_ValNodeByPosition
*
* Revision 1.2  1999/09/21 19:33:53  lewisg
* fix broken declarations
*
* Revision 1.1  1999/09/21 18:09:14  lewisg
* binary search added to color manager, various bug fixes, etc.
*
* Revision 1.4  1999/09/03 23:27:32  lewisg
* minor speedups by avoiding casts
*
* Revision 1.3  1999/09/03 14:01:40  lewisg
* use faster seqid compare SAM_CompareID
*
* Revision 1.2  1999/09/01 23:02:59  lewisg
* binary search in color functions
*
* Revision 1.1  1999/08/13 22:08:16  lewisg
* color manager updated to use alignment coords
*
*
* ==========================================================================
*/

#include <samutil.h>
#include <sequtil.h>


/*****************************************************************************

Function: SAM_ExtractSips()

Purpose: Return a ValNode list containing a SeqId for each Bioseq
         contained in a SeqEntry.
  
Parameters: sep, pointer to the SeqEntry to explore

Returns: ValNode list of pointers to SeqId's.  Do NOT deallocate, these are
         not duplicates!

*****************************************************************************/

typedef struct _SAM_ExtractStruc {
    ValNode *sips;
} SAM_ExtractStruc;

static void SAM_ExtractSipsCallback(SeqEntryPtr sep, Pointer mydata,
                            Int4 index, Int2 indent)
{
    Bioseq *bsp;
    SAM_ExtractStruc *pExtract;

    pExtract = (SAM_ExtractStruc *) mydata;   
    if(sep == NULL || pExtract == NULL) return;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if(bsp == NULL) return;
    ValNodeAddPointer(&pExtract->sips, 0, bsp->id);
}

ValNode * SAM_ExtractSips(SeqEntry *sep)
{
    SAM_ExtractStruc Extract;
    
    Extract.sips = NULL;
    BioseqExplore(sep, (void *)&Extract, SAM_ExtractSipsCallback);
    
    return Extract.sips;
}


/*****************************************************************************

Function: SAM_MakeViewerFree()

Purpose: Set an object to use OM_OPT_FREE_IF_NO_VIEW flag.
  
Parameters: data, pointer to the object to be flagged.

Returns: 1 if OK, 0 otherwise.

*****************************************************************************/

NLM_EXTERN Int4 SAM_MakeViewerFree (void *data)
{
    ObjMgrData *omdp;
    ObjMgr *omp;
    Uint2 entityID;
    Uint2 options;

    omp = ObjMgrReadLock();
    if(omp == NULL) ThrowError;
    omdp = ObjMgrFindByData(omp, data);
    if(omdp == NULL) ThrowError;
    if(!ObjMgrUnlock()) ThrowError;
    entityID = omdp->EntityID;
    options = ObjMgrGetOptions(entityID);
    options |= OM_OPT_FREE_IF_NO_VIEW;
    ObjMgrSetOptions(options, entityID);

    return 1;

error:
    ErrPostEx(SEV_ERROR, 0, 0, "Error");
    return 0;
}

/*****************************************************************************

Function: SAM_MakeTemp()

Purpose: Make an object temporary loaded.
  
Parameters: data, pointer to the object to be flagged.

Returns: 1 if OK, 0 otherwise.

*****************************************************************************/

NLM_EXTERN Int4 SAM_MakeTemp (void *data)
{
    ObjMgr *omp;
 
    omp = ObjMgrWriteLock();
    if(omp == NULL) ThrowError;
    if(!ObjMgrSetTempLoad (omp, data)) ThrowError;
    if(!ObjMgrUnlock()) ThrowError;

    return 1;

error:
    ErrPostEx(SEV_ERROR, 0, 0, "Error");
    return 0;
}

/*****************************************************************************

Function: SAM_ValNodePut()

Purpose: Put a ValNode in a list of ValNodes by cardinal order.
  
Parameters: ppvnHead, pointer to head of ValNode list
            Num, the number of the ValNode to insert before.  begins with 0.
            pvnInsert, the valnode to insert

Returns: 1 if OK, 0 otherwise.

Notes: Hangs on loops that don't include the head.

*****************************************************************************/

NLM_EXTERN Int4 SAM_ValNodePut
(ValNode **ppvnHead, ValNode *pvnInsert, Int4 Num)
{
    ValNode *pvn, *pvnPrevious = NULL;
    Boolean First = TRUE;
    Int4 i;
    
    if(ppvnHead == NULL || pvnInsert == NULL) return 0;
    
    for(pvn = *ppvnHead, i = 0; pvn != NULL && i != Num;
            pvn = pvn->next, i++) {
        if(!First && pvn == *ppvnHead) return 0;  /* loop */
        First = FALSE;
        pvnPrevious = pvn;
    }

    if(pvn == NULL && i != Num) return 0;  /* nobody home */
    if(pvn == NULL && i == Num) {  /* at the end */
        pvnPrevious->next = pvnInsert;
        pvnInsert->next = NULL;
        return 1;
    }
    if(i == 0) { /* at the beginning */
        pvnInsert->next = *ppvnHead;
        *ppvnHead = pvnInsert;
        return 1;
    }
    /* somewhere in the middle */
    pvnPrevious->next = pvnInsert;
    pvnInsert->next = pvn;
    return 1;
}

/*****************************************************************************

Function: SAM_ValNodeExtract()

Purpose: Takes a ValNode out of a list of ValNodes by cardinal order.
  
Parameters: ppvnHead, pointer to head of ValNode list
            Num, the number of the ValNode to extract (the 2nd, 3rd, etc.)

Returns: The extracted ValNode. NULL otherwise.

Notes: Hangs on loops that don't include the head.

*****************************************************************************/

NLM_EXTERN ValNode * SAM_ValNodeExtract(ValNode **ppvnHead, Int4 Num)
{
    ValNode *pvn, *pvnPrevious = NULL;
    Boolean First = TRUE;
    Int4 i;
    
    if(ppvnHead == NULL) return NULL;
    
    for(pvn = *ppvnHead, i = 0; pvn != NULL; pvn = pvn->next, i++) {
        if(!First && pvn == *ppvnHead) return NULL;  /* loop */
        First = FALSE;
        if (i == Num) break;
        pvnPrevious = pvn;
    }
    if(pvn == NULL) return NULL;  /* nobody home */
    
    if(pvnPrevious != NULL) pvnPrevious->next = pvn->next;
    else *ppvnHead = pvn->next;
    
    pvn->next = NULL;
    return pvn;
}

/*****************************************************************************

Function: SAM_ValNodeByPosition()

Purpose: Return pointer to ValNode in a list of ValNodes by cardinal order.
  
Parameters: ppvnHead, pointer to head of ValNode list
            Num, the number of the ValNode to get (the 2nd, 3rd, etc.).
                 starts at 0.

Returns: The ValNode at position NUM. NULL otherwise.

Notes: Hangs on loops that don't include the head.

*****************************************************************************/

NLM_EXTERN ValNode * SAM_ValNodeByPosition(ValNode **ppvnHead, Int4 Num)
{
    ValNode *pvn;
    Boolean First = TRUE;
    Int4 i;
    
    if(ppvnHead == NULL) return NULL;
    
    for(pvn = *ppvnHead, i = 0; pvn != NULL; pvn = pvn->next, i++) {
        if(!First && pvn == *ppvnHead) return NULL;  /* loop */
        First = FALSE;
        if (i == Num) break;
    }
    return pvn;
}


/*****************************************************************************

Function: SAM_SeqAlignExtract()

Purpose: Takes a SeqAlign out of a list of SeqAligns.  The extracted SeqAlign
         is the first one pointer matches the passed SeqAlign pointer.
  
Parameters: psalpHead, pointer to head of SeqAlign list
            salpCheck, points to the SeqAlign to be extracted

Returns: The extracted SeqAlign. NULL otherwise

Notes: Hangs on loops that don't include the head.

*****************************************************************************/

NLM_EXTERN SeqAlign * SAM_SeqAlignExtract
(SeqAlign **psalpHead, SeqAlign *salpCheck)
{
    SeqAlign *salp, *salpPrevious = NULL;
    Boolean fFirst = TRUE;
    
    if(psalpHead == NULL) return NULL;
    
    for(salp = *psalpHead; salp != NULL; salp = salp->next) {
        if(!fFirst && salp == *psalpHead) return NULL;  /* loop */
        fFirst = FALSE;
        if (salp == salpCheck) break;
        salpPrevious = salp;
    }
    if(salp == NULL) return NULL;  /* nobody home */
    
    if(salpPrevious != NULL) salpPrevious->next = salp->next;
    else *psalpHead = salp->next;
    
    salp->next = NULL;
    return salp;
}

/*****************************************************************************
*
*   Adds a SeqAlign newnode to the end of a SeqAlign chain started by head.
*   
*   If the head is NULL, makes the newnode the head.
*   Returns the head of the SeqAlign chain, otherwise returns NULL on error.
*
*****************************************************************************/

NLM_EXTERN SeqAlign * SAM_Add2SeqAlign(SeqAlign ** head, SeqAlign *newnode)
{   
    SeqAlign *salp;
    
    if (head == NULL)
        return NULL;
    salp = *head;
    if (salp != NULL )   {
        while (salp->next != NULL) salp = salp->next;
        salp->next = newnode;
    }
    else
        *head = newnode;
    return *head;
}

/*****************************************************************************
*
*   Frees a list of SeqId's.  Returns the remaining SeqId * if fails,
*   Otherwise NULL. 
*
*****************************************************************************/

NLM_EXTERN SeqId * SAM_FreeSeqIdSet(SeqId *sip)
{
    for(;sip != NULL; sip = sip->next) {
        if(SeqIdFree(sip) != NULL) return sip;
    }
    return NULL;
}

/*****************************************************************************
*
*   Retrieves SeqIds from a set of seqlocs, duplicates the first SeqId, then
*   appends it to a list of SeqId's, which is returned.  
*
*   Returns NULL if the seqlocs are of any type other that SEQLOC_INT or the
*   SeqLoc doesn't contain a SeqId.
*
*****************************************************************************/

NLM_EXTERN SeqId * SAM_SeqIdFromSeqLoc(SeqLoc *slp, Int4 * NumSeqs)
{
    SeqInt *pSeqInt;
    SeqId *sip, *sipHead;

    sip = NULL;
    sipHead = NULL;
    *NumSeqs = 0;
    for(; slp != NULL; slp = slp->next) {

        if(slp->choice != SEQLOC_INT) goto error;

        pSeqInt = (SeqInt *)slp->data.ptrvalue;
        if (pSeqInt == NULL || pSeqInt->id == NULL) goto error;
        
        sip = SeqIdDup(pSeqInt->id);
        if(sip == NULL) goto error;
        sip->next = NULL;
        
        ValNodeLink(&sipHead, sip);
        (*NumSeqs)++;
    }    
    return sipHead;

error:
    ErrPostEx(SEV_ERROR, 0, 0, 
        "SAM_SeqIdFromSeqLoc: Error");
    SAM_FreeSeqIdSet(sipHead);
    return NULL;
}


/*****************************************************************************

Function: SAM_NewDenseSeg()

Purpose: Constructs a new DenseSeg, all arrays initialized to sizes based on
          NumSeqs and NumSegs

Parameters: NumSeqs, the number of sequences
            NumSegs, the number of Segments
            Strands, if TRUE, initialize the strands array.

Returns: The new DenseSeg. NULL otherwise.

*****************************************************************************/

NLM_EXTERN DenseSeg *SAM_NewDenseSeg
(Int4 NumSeqs, Int4 NumSegs, Boolean Strands)
{
    DenseSeg *pDenseSeg = NULL;

    pDenseSeg = DenseSegNew();
    if(pDenseSeg == NULL) goto error;

    pDenseSeg->starts = MemNew(sizeof(Int4)*NumSeqs*NumSegs);
    pDenseSeg->lens = MemNew(sizeof(Int4)*NumSegs);
    if(pDenseSeg->starts == NULL || pDenseSeg->lens == NULL) goto error;
    if(Strands) {
        pDenseSeg->strands = MemNew(sizeof(Uint1)*NumSeqs);
        if(pDenseSeg->strands == NULL) goto error;
    }
    else pDenseSeg->strands = NULL;

    pDenseSeg->dim = NumSeqs;
    pDenseSeg->numseg = NumSegs;
    pDenseSeg->scores = NULL;
    pDenseSeg->ids = NULL;
    return pDenseSeg;

error:
    ErrPostEx(SEV_ERROR, 0, 0, 
        "SAM_NewDenseSeg: Error");
    DenseSegFree(pDenseSeg);
    return NULL;

}

/*****************************************************************************
*
*   Constructs a new SeqAlign of type type and segment type segtype.  Add segs
*   to the seg pointer and sets the dimension of the SeqAlign to dim.
*
*   Returns NULL on error
*
*****************************************************************************/

NLM_EXTERN SeqAlign *SAM_NewSeqAlign
(Uint1 type, Uint1 segtype, Pointer segs, Int2 dim)
{
    SeqAlign *pSeqAlign;

    pSeqAlign = MemNew(sizeof(SeqAlign));

    pSeqAlign->type = type;

    pSeqAlign->segtype = segtype;
    pSeqAlign->dim = dim;
    pSeqAlign->segs = segs;

    pSeqAlign->score = NULL;
    pSeqAlign->next = NULL;
	pSeqAlign->bounds = NULL;
    pSeqAlign->master = NULL;
    pSeqAlign->saip = NULL;
    pSeqAlign->idx.entityID = 0;
    pSeqAlign->idx.itemID = 0;
    pSeqAlign->idx.itemtype = 0;

    return pSeqAlign;
}

/*****************************************************************************
*
*   Checks to see how a postion is inside or is in front of or in back of a
*   range.
*   If it is inside, return SAM_TOTALLAP
*   If it doesn't, return SAM_NOLAP & SAM_NOLAPFRONT if in front
*   If it doesn't and is in back return SAM_NOLAP & SAM_NOLAPBACK
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_InRange(Int4 Position, Int4 From, Int4 To)
{
    return SAM_RangeOverlap(Position, Position, From, To);
}

/*****************************************************************************
*
*   Checks to see if Range1 overlaps Range2.
*   If it does completely, return SAM_TOTALLAP
*   If it doesn't, return SAM_NOLAP & SAM_NOLAPFRONT if in front
*   If it doesn't and is in back return SAM_NOLAP & SAM_NOLAPBACK
*   If pRange1 overlaps the front of pRange2, return SAM_FRONTLAP
*   If pRange1 overlaps the rear of pRange2, return SAM_BACKLAP
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_RangeOverlap(Int4 From1, Int4 To1, Int4 From2, Int4 To2)
{
    if(To1 < From2) return SAM_NOLAPFRONT;
    if(From1 > To2) return SAM_NOLAPBACK;
    if(From1 >= From2) {
        if(To1 <= To2) return SAM_TOTALLAP;
        else return SAM_BACKLAP;
    }
    return SAM_FRONTLAP;
}

/*****************************************************************************
*
*   Lexically compare two SeqId's.  Checks ALL sips on both chains.
*   Returns -1 if sip1 <  sip2
*            0 if sip1 == sip2
*            1 if sip1 >  sip2
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_LexicalComp(SeqId *sip1, SeqId *sip2)
{
    Char Id1[SAM_SIPBUF], Id2[SAM_SIPBUF];
    SeqId *sipThis1, *sipThis2;
    Int4 Compare;

    for(sipThis1 = sip1, sipThis2 = sip2;;
        sipThis1 = sipThis1->next, sipThis2 = sipThis2->next) {

        MakeReversedSeqIdString (sipThis1, Id1, (size_t) SAM_SIPBUF);
        MakeReversedSeqIdString (sipThis2, Id2, (size_t) SAM_SIPBUF);
        
        Compare = StrCmp(Id1, Id2);
        
        if(Compare == 0) {
            if(sipThis1->next == NULL && sipThis2->next == NULL) return 0;
            else if(sipThis1->next == NULL) return -1;
            else if(sipThis2->next == NULL) return 1;
            continue;
        }
        return Compare;
    }
}


/*****************************************************************************
*
*   Orders two SeqId's for binary searches, etc.  DOES check the full
*   ValNode lists.  The ordering is arbitrary but consistent.
*   Returns -1 if sip1 <  sip2
*            0 if sip1 == sip2
*            1 if sip1 >  sip2
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_OrderSeqIDChain (SeqId *sip1, SeqId *sip2)
{
    Int4 retval = 1;

    for(;sip1 != NULL && sip2 != NULL; sip1 = sip1->next, sip2 = sip2->next) {        
        retval = SAM_OrderSeqID(sip1, sip2);
        if(retval != 0) return retval;
    }
    if(sip1 == NULL && sip2 == NULL) return 0;
    else if(sip2 == NULL) return -1;
    else return 1;
}

/*****************************************************************************
*
*   Orders two SeqId's for binary searches, etc. Does NOT check the full
*   ValNode lists.  The ordering is arbitrary but consistent.
*   Returns -1 if sip1 <  sip2
*            0 if sip1 == sip2
*            1 if sip1 >  sip2
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_OrderSeqID(SeqId *sip1, SeqId *sip2)
{
    Char Buf1[SAM_SIPBUF], Buf2[SAM_SIPBUF];
    Int4 retval = 1;
    
    if(sip1->choice == sip2->choice ) goto check;
    
    if((sip1->choice == SEQID_GENBANK ||
        sip1->choice == SEQID_EMBL || sip1->choice == SEQID_DDBJ) && 
        (sip2->choice == SEQID_GENBANK || sip2->choice == SEQID_EMBL ||
        sip2->choice == SEQID_DDBJ)) goto check;
    goto nocheck;
    
check:
    switch (sip1->choice) {
    case SEQID_GI:
    case SEQID_GIBBSQ:
    case SEQID_GIBBMT:
        retval = sip1->data.intvalue - sip2->data.intvalue;
        if (retval == 0) break;
        if (retval > 0) return 1;
        return -1;
        
    case SEQID_LOCAL:
    case SEQID_GIIM:
    case SEQID_GENERAL:
    case SEQID_PDB:
    case SEQID_PATENT:
    case SEQID_PRF:
    case SEQID_DDBJ:
    case SEQID_OTHER:
    case SEQID_EMBL:
    case SEQID_GENBANK:
    case SEQID_PIR:
    case SEQID_SWISSPROT:
        SeqIdWrite (sip1, Buf1, PRINTID_FASTA_SHORT, SAM_SIPBUF);
        SeqIdWrite (sip2, Buf2, PRINTID_FASTA_SHORT, SAM_SIPBUF);
        retval = StrCmp(Buf1, Buf2);
        if(retval != 0) return retval;
        break;
        
    default:
        retval = 1;
        break;
    }
    return retval;
    
nocheck:
    if (sip1->choice > sip2->choice) return 1;
    return -1;
}


/*****************************************************************************
*
*   Compare two SeqId to make sure all valnodes compare exactly
*
*****************************************************************************/

NLM_EXTERN Boolean SAM_SeqIdCompareAll(SeqId *sip1, SeqId *sip2)
{
    SeqId *sip;
    Boolean retval = TRUE;

    if(sip1 == NULL || sip2 == NULL) return FALSE;
    if(ValNodeLen(sip1) != ValNodeLen(sip2)) return FALSE;

    for(sip = sip1; sip != NULL; sip = sip->next)
        if(!SeqIdIn(sip, sip2)) retval = FALSE;

    return retval;
}


/*****************************************************************************
*
*   Compare two AlignId to make sure all valnodes compare exactly
*
*****************************************************************************/

NLM_EXTERN Boolean SAM_AlignIdCompare(AlignId *saip1, AlignId *saip2)
{
    SeqId *saip;
    Boolean retval = TRUE;

    if(saip1 == NULL || saip2 == NULL) return FALSE;
    if(ValNodeLen(saip1) != ValNodeLen(saip2)) return FALSE;

    for(saip = saip1; saip != NULL; saip = saip->next)
        if(!SAM_AlignIdIn(saip, saip2)) retval = FALSE;

    return retval;
}

/*****************************************************************************
*
*     Looks for single AlignId, "a" in chain of AlignIds, "b"
*
*****************************************************************************/

NLM_EXTERN Boolean SAM_AlignIdIn (AlignId *a, AlignId *b)
{
	AlignId *now;
	Uint1 retval;

	if (a == NULL)
	    return FALSE;

	for (now =b; now != NULL; now = now -> next)
	{
        retval = SAM_AlignIdComp(a, now);
		switch (retval)
		{
			case SIC_YES:
				return TRUE;
			case SIC_NO:
				return FALSE;
		}
    }
    return FALSE;
}


/*****************************************************************************
*
*   	Compares a to b and returns
*
*   SIC_DIFF   = different types, could not be compared
*   SIC_NO     = types could be compared, and ids are different
*   SIC_YES    = types could be compared, and ids are the same
*
*****************************************************************************/

NLM_EXTERN Uint1 SAM_AlignIdComp (AlignId *a, AlignId *b)
{
    Uint1 choice;

    if ((a == NULL) || (b == NULL))
        return SIC_DIFF;

	choice = a->choice;
	if (choice != b->choice) return SIC_DIFF;
    switch (choice)
    {
        case AlignId_id:   
            if (ObjectIdMatch((ObjectIdPtr)a->data.ptrvalue, (ObjectIdPtr)b->data.ptrvalue))
				return SIC_YES;
			else
				return SIC_NO;
        case AlignId_gi:  /* gi */
            if (a->data.intvalue == b->data.intvalue)
				return SIC_YES;
			else
				return SIC_NO;
        case AlignId_itemid:  /* "permanent" itemid */
            if (a->data.intvalue == b->data.intvalue)
				return SIC_YES;
			else
				return SIC_NO;
		default:
			ErrPostEx(SEV_ERROR, 0,0, "AlignIdComp: unsupported type [%d]",
				(int)choice);
			return SIC_DIFF;
     }
}

/*******************************************************
*
*   duplicate a list of AlignId *
*
*******************************************************/

NLM_EXTERN AlignId * SAM_AlignIdDupList (AlignId *id_list)
{
  SeqId *sip=NULL;
  SeqId *sid;

  for (sid = id_list; sid != NULL; sid = sid->next) {
         ValNodeLink(&sip, SAM_AlignIdDup(sid));  
  }
  return sip;
}

/*******************************************************
*
*   Duplicates one AlignId
*
*******************************************************/

NLM_EXTERN AlignId * SAM_AlignIdDup (AlignId *oldid)
{
	AlignId *newid = NULL;

    if (oldid == NULL)
        return oldid;

	newid = ValNodeNew(NULL);
	if (newid == NULL) return newid;
	MemCopy(newid, oldid, sizeof(ValNode));
	newid->next = NULL;    /* not in chain */
    switch (oldid->choice)
    {
        case AlignId_id:
			newid->data.ptrvalue = ObjectIdDup((ObjectIdPtr)oldid->data.ptrvalue);
			break;
        case AlignId_itemid:
        case AlignId_gi:
            break;
     }
	return newid;
}
