/* ===========================================================================
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
*  any work or product based on this material.
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
* File Name:  alignmgr.h
*
* Author:  Sarah Wheelan
*
* Version Creation Date:   7/99
*
* $Revision: 6.14 $
*
* File Description: SeqAlign indexing and messaging functions 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: alignmgr.h,v $
* Revision 6.14  1999/09/17 16:55:59  wheelan
* added AlnMgrPropagateSeqIdsBySapList
*
* Revision 6.13  1999/09/13 14:34:13  wheelan
* added more functions to support accessing the alignment by row number
*
* Revision 6.12  1999/09/06 16:37:36  wheelan
* added AlnMgrGetNextLengthBit and associated function
*
* Revision 6.11  1999/09/06 15:51:31  wheelan
* added row management functions, new structure for storing row information
*
* Revision 6.10  1999/09/01 20:12:18  wheelan
* added new merge function and the typedef for the structure it uses
*
* Revision 6.9  1999/09/01 14:38:46  wheelan
* added AlnMgrGetStrand
*
* Revision 6.8  1999/08/30 19:28:52  wheelan
* Added AlnMgrGetSapForSip for master-slave alignments
*
* Revision 6.7  1999/08/26 20:35:54  wheelan
* added parent indexing and pairwise-to-multiple functions
*
* Revision 6.6  1999/08/19 17:24:46  wheelan
* changed AMAlignIndex structure, added more api functions
*
* Revision 6.5  1999/08/12 20:56:57  vakatov
* [WIN32] Added missed LIBCALLBACK
*
* Revision 6.4  1999/08/12 12:41:56  wheelan
* added comments, and functions to index the parent
*
* Revision 6.3  1999/08/06 13:43:48  wheelan
* added several api functions; changed all names to AlnMgr..
*
* Revision 6.2  1999/07/30 14:08:41  wheelan
* added api functions to access indexes
*
* Revision 6.1  1999/07/29 12:56:45  wheelan
* initial checkin
*

* ==========================================================================
*/

#ifndef _ALIGNMGR_
#define _ALIGNMGR_

#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>
#include <seqmgr.h>
#include <salutil.h>
#include <sequtil.h>
#include <salpedit.h>

/* for SAIndex indextype field */
#define INDEX_SEGS 1
#define INDEX_PARENT 2

/* return values for AlnMgrCheckAlignForParent */
#define AM_CHILD 1
#define AM_PARENT 2

/* return values for AlnMgrCheckOverlapping */
#define CHECK_ERROR -2
#define NO_OVERLAP -1


/***************************************************************************
*
*  Each child seqalign (segtype SAS_DENSEG) has an SAIndex structure, 
*  (in the sap->saip field -- the SAIndex is saip->indextype INDEX_SEGMENT
*  and may be accessed by casting the sap->saip pointer to an SAIndexPtr)
*  which records the alignment coordinates (each child must be a 
*  continuous alignment, so each will have obvious coordinates).  
*  The SAIndex structure has an array of SASeqDat structures, one 
*  for each "row" of the alignment, and each SASeqDat structure 
*  contains an array of Uint2s representing the segment numbers that 
*  the sequence in that row participates in (not gapped).
*
***************************************************************************/

typedef struct saseqdat {
   Uint2Ptr  sect;
   Uint2Ptr  unsect;
   Uint2     numsect;
   Uint2     numunsect;
} SASeqDat, PNTR SASeqDatPtr;

typedef struct saindex{
   Uint1 indextype;
   SeqAlignIndexFreeFunc freefunc;
   Uint4Ptr     aligncoords;
   Int4 master;
   SASeqDatPtr PNTR  ssdp;
} SAIndex, PNTR SAIndexPtr;

SASeqDatPtr SASeqDatNew(void);
SAIndexPtr SAIndexNew(void);
Boolean SAIndexFree(VoidPtr index);

/***************************************************************************
*
*  Each parent seqalign (segtype SAS_DISC) has an AMAlignIndex structure
*  (in the sap->saip field -- the AMAlignIndex is saip->indextype 
*  INDEX_PARENT and may be accessed by casting the sap->saip pointer to
*  an AMAlignIndexPtr).  The AMAlignIndex has alignment coordinates across
*  all children, if possible, or across a subset of the children.  The
*  saps field is an array of SeqAlignPtrs, each of which appears in the 
*  order that it occurs in the combined alignment -- this means that a given
*  pointer may occur more than once.  The AMAlignDatPtr PNTR
*  points to an array of AMAlignDat structure, one per bioseq, each of
*  which has a sorted list of saps that the bioseq participates in.
*
***************************************************************************/
typedef struct row_source {
   SeqIdPtr    id;
   Uint4Ptr    which_saps;
   Uint4Ptr    num_in_sap;
   Uint4       numsaps;
} RowSource, PNTR RowSourcePtr;

RowSourcePtr RowSourceNew(void);
RowSourcePtr RowSourceFree(RowSourcePtr rsp);

typedef struct amaligndat {
   SeqAlignPtr PNTR saps;
   Int4        numsaps;
   Uint2Ptr    segments;
   Uint2       numseg;
} AMAlignDat, PNTR AMAlignDatPtr;

typedef struct amalignindex {
   Uint1 indextype;
   SeqAlignIndexFreeFunc freefunc;
   Uint4Ptr aligncoords;
   SeqAlignPtr PNTR  saps;
   Uint4 numseg;
   Int4Ptr lens;
   Int4Ptr starts;
   Int4 alnsaps; /* the number of child seqaligns contained in the multiple */
   Int4 numsaps; /* the total number of child seqaligns */
   SeqIdPtr ids;
   Int4 numbsqs;
   RowSourcePtr PNTR rowsource;
   Uint4 numrows;
   Int4 master; /*tells which row is the master row*/
   AMAlignDatPtr PNTR amadp;
} AMAlignIndex, PNTR AMAlignIndexPtr;

AMAlignIndexPtr AMAlignIndexNew(void);
Boolean AMAlignIndexFree(VoidPtr index);
AMAlignDatPtr AMAlignDatNew(void);
AMAlignDatPtr AMAlignDatFree(AMAlignDatPtr amadp);


/***************************************************************************
*
*  The AlnMsg structure is used by the function AlnMgrGetNextAlnBit.  The
*  calling function sets the fields which_master (NULL means alignment
*  coordinates), from_m, to_m (in alignment or master coordinates) and
*  either which_bsq or row_num (row_num avoids the repeated SeqId problem).  If
*  from_m and to_m are not set, the entire alignment is returned.  The function
*  fills in the AlnMsg structure and returns TRUE if there are more
*  segments of that bioseq in the region requested, FALSE if not. Be sure
*  not to set any other field -- many of these are used by the function
*  to indicate what has already been returned.  Also, make sure to call
*  AlnMsgNew instead of MemNew to allocate a new structure, as a couple of
*  fields must be initialized for the function to work correctly.
*
**************************************************************************/

typedef struct messagestruct {
/* set by function asking for alignment info */
   SeqIdPtr  which_master;
   Int4      from_m;
   Int4      to_m;
   SeqIdPtr  which_bsq;
   Int4      row_num;
/* ONLY set by AlnMgrGetNextAlnBit*/
   Int4      from_b;
   Int4      to_b;
   Uint1     gap;
   Uint1     strand;
   Int4      prev;
   Int4      real_from;
   Int4      prev_sap;
   Int4      len_left;
   Boolean   send_space;
} AlnMsg, PNTR AlnMsgPtr;

AlnMsgPtr AlnMsgNew(void);

typedef struct am_tinyinfo {
   Int4  start;
   Uint4  numgap;
   Uint2  which;
} AMTinyInfo, PNTR AMTinyInfoPtr;


SeqAlignIndexPtr SeqAlignIndexNew(void);


Boolean AlnMgrIndexSeqAlign(SeqAlignPtr sap);

/***************************************************************************
*
*  AlnMgrAnythingToSeg takes any SeqAlign and does an in-place transformation
*  to the parent-child structure.  Each dense-seg, dense-diag and std-seg 
*  is put into its own seqalign, and the child seqaligns are linked 
*  together in no particular order and put in the sap->segs field of the
*  new parent (which takes over the pointer passed in).  The parent
*  has segtype SAS_DISC, and each child has segtype SAS_DENSEG or SAS_STD.
*  Each child, then, is a continuous, nonoverlapping alignment and therefore
*  may be indexed.
*
***************************************************************************/

Boolean AlnMgrAnythingToSeg (SeqAlignPtr sap);

/***********************************************************************
*
*  AlnMgrIndexLinkedSegs and AlnMgrIndexParentSA create and fill in the
*  SASeqIndex and AMAlignIndex structures on the children and the parent,
*  respectively.  IndexLinkedSegs is called on the sap->segs field of
*  the parent, so that the pointer of the first child in the list
*  gets passed in.  AlnMgrIndexParentSA is called on the parent, and 
*  the children must already be indexed (the function does check) in order
*  for it to work.  AlnMgrIndexParentSA calls AlnMgrPropagateUpSeqIdPtrs
*  to create a list of all SeqIdPtrs present in all the children (each
*  is only listed once, in the order that its AMAlignDat structure occurs
*  in).
*
***********************************************************************/
Boolean AlnMgrIndexLinkedSegs (SeqAlignPtr sap);
Boolean AlnMgrIndexParentSA(SeqAlignPtr sap);
SeqIdPtr AlnMgrPropagateUpSeqIdPtrs(SeqAlignPtr sap, Int4Ptr num);
SeqIdPtr  AlnMgrPropagateSeqIdsBySapList(AMAlignIndexPtr amaip);


/**********************************************************************
*
*  am_print_seqalign_indexes is a useful function for debugging; 
*  given a parent or child seqalign it prints a text summary of the
*  indexes for that seqalign.
*
***********************************************************************/
void am_print_seqalign_indexes(SeqAlignPtr sap);


/**********************************************************************
*
*  AlnMgrCheckAlignForParent is called by almost every other function;
*  it checks whether the parents and children are correctly labeled
*  and indexed, and it returns AM_CHILD to indicate that the sap 
*  passed in is a child, AM_PARENT for a parent, or -1 for an error.
*  If the indexes aren't present, it will try to create them.
*
**********************************************************************/
Int4 AlnMgrCheckAlignForParent(SeqAlignPtr sap);


/***********************************************************************
*
*  AlnMgrSortSeqAligns is a variant of the ValNodeSort function, and
*  calls very similar heapsort functions.  It can take a comparison
*  function that needs userdata, so more specific sorts are possible
*  without defining special structures for every type of sort.
*
***********************************************************************/
SeqAlignPtr PNTR AlnMgrSortSeqAligns (SeqAlignPtr sap, int (LIBCALLBACK *compare)(VoidPtr, VoidPtr, VoidPtr), VoidPtr userdata, Int4Ptr numsap);
static void heapsort_with_userdata (VoidPtr b, size_t nel, size_t width, int (LIBCALLBACK *compar)PROTO((VoidPtr, VoidPtr, VoidPtr)), VoidPtr userdata);
static void heapify_with_userdata(CharPtr base0, CharPtr base, CharPtr lim, CharPtr last, size_t width, int(LIBCALLBACK *compar)PROTO((VoidPtr, VoidPtr, VoidPtr)), VoidPtr userdata);

/**********************************************************************
*
*  AlnMgrCompareIncreasingBySeqIdPtr takes a SeqIdPtr as userdata,
*  and sorts the alignments in increasing order according to the 
*  region of the bioseq indicated that is contained in the alignment.
*  If the bioseq is not in the alignment, the alignment will be put
*  first, so all alignments in which the given bioseq does not 
*  participate occur at the beginning of the list, making it easy to
*  check for them and remove them.
*
**********************************************************************/
int LIBCALLBACK AlnMgrCompareIncreasingBySeqIdPtr (VoidPtr base, VoidPtr large_son, VoidPtr userdata);

/*********************************************************************
*
*  AlnMgrFindFirst is crucial to the AlnMgrMakeFakeMultiple function;
*  it uses the sorted order of the seqaligns in each AMAlignDat
*  structure to guide a heapsort of all the seqaligns.
*
*********************************************************************/
int LIBCALLBACK AlnMgrFindFirst(VoidPtr base, VoidPtr large_son, VoidPtr userdata);
int LIBCALLBACK AlnMgrCompareTips(VoidPtr base, VoidPtr large_son);

/************************************************************************
*
*  AlnMgrGetNextLengthBit should be called iteratively on an alignment
*  to return the lengths of all the aligned and unaligned pieces in
*  the alignment.  Don't change the value in r, just pass in a pointer 
*  to an allocated Int4 set to 0 initially.  length will be positive
*  for aligned regions and negative for unaligned regions (it is
*  set to the negative of the largest unaligned sequence in the region).
*  AlnMgrGetNextLengthBit calls AlnMgrGetMaxUnalignedLength to get
*  the lengths of the unaligned regions.  This function makes a lot of
*  assumptions about the relationship between sap1 and sap2; call it
*  with much hesitation.
*
************************************************************************/
Boolean AlnMgrGetNextLengthBit(SeqAlignPtr sap, Int4Ptr length, Int4Ptr r);
Int4 AlnMgrGetMaxUnalignedLength(SeqAlignPtr sap1, SeqAlignPtr sap2);


/*************************************************************************
*
*  AlnMgrGetNextAlnBit takes an AlnMsgPtr, with (at the minimum) the
*  which_bsq field filled in to indicate which bioseq should be returned.
*  The function returns the segments of the bioseq which span the region
*  of the alignment indicated, and can return them according to either
*  alignment coordinates (if which_master is NULL) or a master coordinate
*  system (need to fill in the SeqIdPtr of the master).  The function
*  returns TRUE if there are more segments of the bioseq to retrieve,
*  and FALSE if not.  It uses the two binary search functions to quickly
*  retrieve the required data from the indexes.
*
*************************************************************************/
Boolean AlnMgrGetNextAlnBit(SeqAlignPtr sap, AlnMsgPtr amp);
Uint4 binary_search_on_uint4_list(Uint4Ptr list, Uint4 pos, Uint4 listlen);
Int4 binary_search_on_uint2_list(Uint2Ptr list, Uint2 ele, Uint2 listlen);
Int4 binary_search_by_chunk(Int4Ptr list, Int4 ele, Int4 listlen, Int4 chunksize, Int4 offset);
Int4 binary_search_segment_array(SASeqDatPtr ssdp, Int4 pos, Int4 numseq, Int4 offset, DenseSegPtr dsp);



/************************************************************************
*
*  These are several utility functions which get needed data from the
*  indexes.  "N" refers to row number.
*
************************************************************************/
Int4 AlnMgrGetAlnLength(SeqAlignPtr sap, Boolean fill_in);
Int4 AlnMgrGetNumSeqs(SeqAlignPtr sap);
SeqIdPtr AlnMgrGetNthSeqIdPtr(SeqAlignPtr sap, Int4 n);
void AlnMgrGetNthSeqRangeInSA(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop);
Uint1 AlnMgrGetStrand(SeqAlignPtr sap, SeqIdPtr sip);
Uint1 AlnMgrGetNthStrand(SeqAlignPtr sap, Int4 n);
Int4 AlnMgrGetNForSip(SeqAlignPtr sap, SeqIdPtr sip);
Int4 AlnMgrGetSapForSip(AMAlignIndexPtr amaip, SeqIdPtr sip, Int4 which);
Int4 AlnMgrMapToBsqCoords(SeqAlignPtr sap, Uint4 pos, SeqIdPtr sip, SeqIdPtr master); 
Int4 AlnMgrMapRowCoords(SeqAlignPtr sap, Uint4 pos, Int4 row, SeqIdPtr master);


/***********************************************************************
*
*  AlnMgrMakeFakeMultiple calls AlnMgrCheckOverlapping to decide whether
*  an alignment is linear.  Then, if possible, it calls AlnMgrMakeAlignCoords
*  to create alignment coordinates across all children contained in the
*  parent.
*
***********************************************************************/
Boolean AlnMgrMakeFakeMultiple(SeqAlignPtr sap);
void AlnMgrMakeAlignCoords(SeqAlignPtr sap);
Boolean AlnMgrMergeIntoMSMultByMaster(AMAlignIndexPtr amaip, Int4Ptr lens, Uint4Ptr numseg);
Boolean AlnMgrMergeSegments(Int4Ptr lens, SeqAlignPtr sap, SeqIdPtr master, Int4Ptr where, Int4 which);
Boolean AlnMgrFillInStarts(SeqAlignPtr PNTR saparray, Int4Ptr starts, Int4 numseg, Int4Ptr lens, Int4 numsaps, Uint4Ptr aligncoords);
Int4 AlnMgrGetStartFromMaster(SeqAlignPtr sap, Int4 start);
Uint4 AlnMgrGetMasterGapStartForSeg(SeqAlignPtr sap, Int4 which_gap, Uint4Ptr aligncoord);
Boolean AlnMgrReconcileGaps(Int4Ptr lens, Uint4Ptr aligncoords, Int4 num);
Boolean AlnMgrMakeMultSegments(AMAlignIndexPtr amaip);

Int4 AlnMgrCheckOverlapping(SeqAlignPtr sap);

/******************************************************************************
*
*  AlnMgrGetMaxSegments simply adds up the number of segments for each
*  SeqAlign in a linked list, to get the maximum number of segments
*  for the merge of the list (for memory allocation in AlnMgrMakeFakeMultiple).
*
*******************************************************************************/
Int4 AlnMgrGetMaxSegments(SeqAlignPtr sap);

/*******************************************************************************
*
*  Row Management functions:
*
*******************************************************************************/
Int4 AlnMgrGetNumRows(SeqAlignPtr sap);
Int4 AlnMgrGetMaxRowsForParentPartial(SeqAlignPtr sap);
Boolean AlnMgrGetRowsForPartial(SeqAlignPtr sap);
Boolean AlnMgrGetRowsForMasterSlave(SeqAlignPtr sap);


/*******************************************************************************
*
*  AlnMgrFindMaster returns the SeqIdPtr of the first bioseq present in all
*  child alignments, unless the sap->master field is set in the child alignments,
*  in which case that SeqIdPtr is returned (if it's the same in all children).
*
*******************************************************************************/
SeqIdPtr AlnMgrFindMaster(SeqAlignPtr sap);

/*******************************************************************************
*
*  AlnMgrCheckRealMaster makes sure that the master seqid given appears
*  once and only once in each seqalign in the set if a parent is given,
*  or once and only one in the seqalign if a child is given.
*
*******************************************************************************/
Boolean AlnMgrCheckRealMaster(SeqAlignPtr sap, SeqIdPtr master);

void AlnMgrSetMaster(SeqAlignPtr sap, SeqIdPtr master);
void AlnMgrMakeMasterPlus(SeqAlignPtr sap);

SeqAlignPtr AlnMgrGetSubAlign(SeqAlignPtr sap, SeqIdPtr which_master, Int4 from, Int4 to);

#ifdef __cplusplus
}
#endif

#endif
