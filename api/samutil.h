/*   $Id: samutil.h,v 1.9 2000/01/24 20:54:35 vakatov Exp $
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
* File Name:  $Id: samutil.h,v 1.9 2000/01/24 20:54:35 vakatov Exp $
*
* Author:  Lewis Geer
*
* Version Creation Date:   8/12/99
*
* $Revision: 1.9 $
*
* File Description: Utility functions for AlignIds and SeqAlignLocs
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: samutil.h,v $
* Revision 1.9  2000/01/24 20:54:35  vakatov
* SAM_ViewString::  made #define to fix for the DLL build on PC
*
* Revision 1.8  1999/12/11 01:30:34  lewisg
* fix bugs with sharing colors between ddv and cn3d
*
* Revision 1.7  1999/12/03 23:17:23  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 1.6  1999/11/24 21:24:30  vakatov
* Fixed for the C++ and/or MSVC DLL compilation
*
* Revision 1.5  1999/11/15 18:30:07  lewisg
* get rid of extra redraws when selecting
*
* Revision 1.4  1999/10/05 23:18:16  lewisg
* add ddv and udv to cn3d with memory management
*
* Revision 1.3  1999/09/27 17:49:12  lewisg
* fix denseseg constructor, bug in valnode loops, add SAM_ValNodeByPosition
*
* Revision 1.2  1999/09/21 19:38:59  lewisg
* update broken declarations
*
* Revision 1.1  1999/09/21 18:09:14  lewisg
* binary search added to color manager, various bug fixes, etc.
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


#ifndef SAMUTIL_H
#define SAMUTIL_H

#include <ncbistd.h>
#include <objall.h>
#include <sequtil.h>
#include <objgen.h>
#include <objalignloc.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif


#define SAM_SIPBUF 80

#define SAM_NOLAP 0x1
#define SAM_NOLAPFRONT 0x3
#define SAM_NOLAPBACK 0x5
#define SAM_TOTALLAP 0x8
#define SAM_FRONTLAP 0x10
#define SAM_BACKLAP 0x20

#define ThrowError if(!Nlm_ErrSetContext(THIS_MODULE,THIS_FILE,__LINE__,DBFLAG,0,0,0)) goto error

typedef SeqEntryPtr (*SAM_SeqFetchProc) (Int4 uid, Int2 retcode);

/* same as Vibrant RecT */
typedef  struct  _SAM_RecT {
  Nlm_Int2  left;
  Nlm_Int2  top;
  Nlm_Int2  right;
  Nlm_Int2  bottom;
} SAM_RecT, PNTR SAM_RecTPtr;

/* defines so that views can be identified */
#define SAMVIEWSELF 0
#define SAMVIEWCN3D 1
#define SAMVIEWDDV 2
#define SAMVIEWUDV 3
#define SAMVIEWSEQUIN 4
#define SAMVIEWNENTREZ 5

#define SAM_ViewString "Viewer Global"

typedef  struct  _SAM_ViewerGlobal {
  SAM_RecT Rect;  /* where the master viewer wants the called viewer to go */
  SAM_SeqFetchProc  FetchProc;  /* the fetch function the master has set */
  Int4 MasterViewer;  /* the name of the master viewer */
} SAM_ViewGlobal, PNTR SAM_ViewGlobalPtr;


/*****************************************************************************

Function: SAM_ExtractSips()

Purpose: Return a ValNode list containing a SeqId for each Bioseq
         contained in a SeqEntry.
  
Parameters: sep, pointer to the SeqEntry to explore

Returns: ValNode list of pointers to SeqId's.  Do NOT deallocate, these are
         not duplicates!

*****************************************************************************/

NLM_EXTERN ValNode * SAM_ExtractSips(SeqEntry *sep);


/*****************************************************************************

Function: SAM_MakeViewerFree()

Purpose: Set an object to use OM_OPT_FREE_IF_NO_VIEW flag.
  
Parameters: data, pointer to the object to be flagged.

Returns: 1 if OK, 0 otherwise.

*****************************************************************************/

NLM_EXTERN Int4 SAM_MakeViewerFree (void *data);


/*****************************************************************************

Function: SAM_MakeTemp()

Purpose: Make an object temporary loaded.
  
Parameters: data, pointer to the object to be flagged.

Returns: 1 if OK, 0 otherwise.

*****************************************************************************/

NLM_EXTERN Int4 SAM_MakeTemp (void *data);


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
(ValNode **ppvnHead, ValNode *pvnInsert, Int4 Num);


/*****************************************************************************

Function: SAM_ValNodeExtract()

Purpose: Takes a ValNode out of a list of ValNodes by cardinal order.
  
Parameters: ppvnHead, pointer to head of ValNode list
            Num, the number of the ValNode to extract (the 2nd, 3rd, etc.)

Returns: The extracted ValNode. NULL otherwise.

Notes: Hangs on loops that don't include the head.

*****************************************************************************/

NLM_EXTERN ValNode * SAM_ValNodeExtract(ValNode **ppvnHead, Int4 Num);

/*****************************************************************************

Function: SAM_ValNodeByPosition()

Purpose: Return pointer to ValNode in a list of ValNodes by cardinal order.
  
Parameters: ppvnHead, pointer to head of ValNode list
            Num, the number of the ValNode to get (the 2nd, 3rd, etc.).
                 starts at 0.

Returns: The ValNode at position NUM. NULL otherwise.

Notes: Hangs on loops that don't include the head.

*****************************************************************************/

NLM_EXTERN ValNode * SAM_ValNodeByPosition(ValNode **ppvnHead, Int4 Num);


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
(SeqAlign **psalpHead, SeqAlign *salpCheck);

/*****************************************************************************
*
*   Adds a SeqAlign newnode to the end of a SeqAlign chain started by head.
*   
*   If the head is NULL, makes the newnode the head.
*   Returns the head of the SeqAlign chain, otherwise returns NULL on error.
*
*****************************************************************************/

NLM_EXTERN SeqAlign * SAM_Add2SeqAlign(SeqAlign ** head, SeqAlign *newnode);

/*****************************************************************************
*
*   Frees a list of SeqId's.  Returns the remaining SeqId * if fails,
*   Otherwise NULL. 
*
*****************************************************************************/

NLM_EXTERN SeqId * SAM_FreeSeqIdSet(SeqId *sip);

/*****************************************************************************
*
*   Retrieves SeqIds from a set of seqlocs, duplicates the first SeqId, then
*   appends it to a list of SeqId's, which is returned.  
*
*   Returns NULL if the seqlocs are of any type other that SEQLOC_INT or the
*   SeqLoc doesn't contain a SeqId.
*
*****************************************************************************/

NLM_EXTERN SeqId * SAM_SeqIdFromSeqLoc(SeqLoc *slp, Int4 *NumSeqs);


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
(Int4 NumSeqs, Int4 NumSegs, Boolean Strands);


/*****************************************************************************
*
*   Constructs a new SeqAlign of type type and segment type segtype.  Add segs
*   to the seg pointer and sets the dimension of the SeqAlign to dim.
*
*   Returns NULL on error
*
*****************************************************************************/

NLM_EXTERN SeqAlign *SAM_NewSeqAlign
(Uint1 type, Uint1 segtype, Pointer segs, Int2 dim);


/*****************************************************************************
*
*   Checks to see how a postion is inside or is in front of or in back of a
*   range.
*   If it is inside, return SAM_TOTALLAP
*   If it doesn't, return SAM_NOLAP & SAM_NOLAPFRONT if in front
*   If it doesn't and is in back return SAM_NOLAP & SAM_NOLAPBACK
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_InRange(Int4 Position, Int4 From, Int4 To);

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

NLM_EXTERN Int4 SAM_RangeOverlap(Int4 From1, Int4 To1, Int4 From2, Int4 To2);

/*****************************************************************************
*
*   Lexically compare two SeqId's.  Checks ALL sips on both chains.
*   Returns -1 if sip1 <  sip2
*            0 if sip1 == sip2
*            1 if sip1 >  sip2
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_LexicalComp(SeqId *sip1, SeqId *sip2);

/*****************************************************************************
*
*   Orders two SeqId's for binary searches, etc.  DOES check the full
*   ValNode lists.  The ordering is arbitrary but consistent.
*   Returns -1 if sip1 <  sip2
*            0 if sip1 == sip2
*            1 if sip1 >  sip2
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_OrderSeqIDChain (SeqId *sip1, SeqId *sip2);

/*****************************************************************************
*
*   Orders two SeqId's for binary searches, etc. Does NOT check the full
*   ValNode lists.  The ordering is arbitrary but consistent.
*   Returns -1 if sip1 <  sip2
*            0 if sip1 == sip2
*            1 if sip1 >  sip2
*
*****************************************************************************/

NLM_EXTERN Int4 SAM_OrderSeqID(SeqId *sip1, SeqId *sip2);

/*****************************************************************************
*
*   Orders two SeqId's for binary searches, etc.  If any id's in chains are
*   equivalent, returns 0.
*   Returns -1 if sip1 <  sip2
*            0 if sip1 == sip2
*            1 if sip1 >  sip2
*
*****************************************************************************/

NLM_EXTERN Boolean SAM_SeqIdCompareAll(SeqId *sip1, SeqId *sip2);

/*****************************************************************************
*
*   Compare two AlignId to make sure all valnodes compare exactly
*
*****************************************************************************/

NLM_EXTERN Boolean SAM_AlignIdCompare(AlignId *saip1, AlignId *saip2);


/*****************************************************************************
*
*     Looks for single AlignId, "a" in chain of AlignIds, "b"
*
*****************************************************************************/

NLM_EXTERN Boolean SAM_AlignIdIn (AlignId *a, AlignId *b);


/*****************************************************************************
*
*   	Compares a to b and returns
*
*   SIC_DIFF   = different types, could not be compared
*   SIC_NO     = types could be compared, and ids are different
*   SIC_YES    = types could be compared, and ids are the same
*
*****************************************************************************/

NLM_EXTERN Uint1 SAM_AlignIdComp (AlignId *a, AlignId *b);


/*******************************************************
*
*   duplicate a list of AlignId *
*
*******************************************************/

NLM_EXTERN AlignId * SAM_AlignIdDupList (AlignId *id_list);


/*******************************************************
*
*   Duplicates one AlignId
*
*******************************************************/

NLM_EXTERN AlignId * SAM_AlignIdDup (AlignId *oldid);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* SAMUTIL_H */
