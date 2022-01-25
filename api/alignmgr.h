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
* $Revision: 6.33 $
*
* File Description: SeqAlign indexing and messaging functions 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: alignmgr.h,v $
* Revision 6.33  2000/01/12 17:43:20  wheelan
* added AlnMgrGetNumSegments, AlnMgrDeleteRow
*
* Revision 6.32  1999/11/30 14:36:27  wheelan
* added AlnMgrMakeMultipleByScore
*
* Revision 6.31  1999/11/26 15:42:21  vakatov
* Fixed for the C++ and/or MSVC DLL compilation
*
* Revision 6.30  1999/11/18 19:30:24  wheelan
* added AlnMgrDeleteChildByPointer
*
* Revision 6.29  1999/10/25 18:16:44  wheelan
* Added AlnMgrGetUniqueSeqs
*
* Revision 6.28  1999/10/19 19:25:57  wheelan
* added is_aligned return to AlnMgrGetNextNthSeqRange; removed static functions; changed structures for AlnMgrMakeSegmentedMasterSlave
*
* Revision 6.27  1999/10/15 21:51:03  durand
* add AlnMgrIsSAPDiscAli()
*
* Revision 6.26  1999/10/15 15:41:20  wheelan
* fixed typo
*
* Revision 6.25  1999/10/15 15:15:05  wheelan
* added defines for AlnMgrGetNthRowTail
*
* Revision 6.24  1999/10/15 13:48:32  wheelan
* added AlnMgrGetNthRowTail
*
* Revision 6.23  1999/10/14 16:10:31  kans
* new includes and prototypes added
*
* Revision 6.22  1999/10/13 19:28:35  wheelan
* added speedup for segmented master-slave creation
*
* Revision 6.21  1999/10/07 13:37:33  wheelan
* added AlnMgrIndexSingleSeqAlign
*
* Revision 6.20  1999/10/06 19:34:49  wheelan
* added several viewer and editor management functions (AlnMgrCopy . . . )
*
* Revision 6.19  1999/10/05 15:15:18  wheelan
* added AlnMgrGetNthUnalignedForNthRow
*
* Revision 6.18  1999/10/04 14:58:17  wheelan
* added AlnMgrMapBioseqToSeqAlign
*
* Revision 6.17  1999/09/23 16:03:28  wheelan
* Added structures and functions to support segmented master-slave alignments
*
* Revision 6.16  1999/09/22 13:20:04  wheelan
* Added AlnMgrGetNextNthSeqRange and structures for segmented master-slave handling
*
* Revision 6.15  1999/09/21 19:14:15  wheelan
* added functions to make segmented master-slave, added fields to AlnMsg and AMAlignIndex structures
*
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

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>
#include <seqmgr.h>
#include <salutil.h>
#include <sequtil.h>
#include <salpedit.h>
#include <samutil.h>

/* for SAIndex indextype field */
#define INDEX_SEGS 1
#define INDEX_PARENT 2

/* return values for AlnMgrCheckAlignForParent */
#define AM_CHILD 1
#define AM_PARENT 2

/* return values for AlnMgrCheckOverlapping */
#define CHECK_ERROR -2
#define NO_OVERLAP -1

/* values for amaip->mstype */
#define AM_MASTERSLAVE 1
#define AM_SEGMENTED_MASTERSLAVE 2

/* values for AlnMgrGetNthRowTail */
#define LEFT_TAIL 1
#define RIGHT_TAIL 0

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

NLM_EXTERN SASeqDatPtr SASeqDatNew(void);
NLM_EXTERN SAIndexPtr SAIndexNew(void);
NLM_EXTERN Boolean SAIndexFree(VoidPtr index);

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
   Uint1       strand;
   struct row_source PNTR next;
} RowSource, PNTR RowSourcePtr;

NLM_EXTERN RowSourcePtr RowSourceNew(void);
NLM_EXTERN RowSourcePtr RowSourceFree(RowSourcePtr rsp);

typedef struct amaligndat {
   SeqIdPtr    sip;
   SeqAlignPtr PNTR saps;
   Int4        numsaps;
   Uint2Ptr    segments;
   Uint2       numseg;
} AMAlignDat, PNTR AMAlignDatPtr;

typedef struct amalignindex {
   Uint1 indextype;
   SeqAlignIndexFreeFunc freefunc;
   Int2 mstype;
   Uint4Ptr aligncoords;
   SeqAlignPtr PNTR  saps;
   Uint4 numseg;
   Int4Ptr lens;
   Int4Ptr ulens;
   Int4Ptr starts;
   Int4 alnsaps; /* the number of child seqaligns contained in the multiple */
   Int4 numsaps; /* the total number of child seqaligns */
   SeqIdPtr ids;
   Int4 numbsqs;
   RowSourcePtr PNTR rowsource;
   Uint4 numrows;
   Int4 master; /*tells which row is the master row*/
   AMAlignDatPtr PNTR amadp;
   SeqAlignPtr parent;
} AMAlignIndex, PNTR AMAlignIndexPtr;

NLM_EXTERN AMAlignIndexPtr AMAlignIndexNew(void);
NLM_EXTERN Boolean AMAlignIndexFree(VoidPtr index);
NLM_EXTERN AMAlignDatPtr AMAlignDatNew(void);
NLM_EXTERN AMAlignDatPtr AMAlignDatFree(AMAlignDatPtr amadp);


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
   Int4      row_num;  /*  1-based numbering  */
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
   Uint2     place;
} AlnMsg, PNTR AlnMsgPtr;

NLM_EXTERN AlnMsgPtr AlnMsgNew(void);
NLM_EXTERN AlnMsgPtr AlnMsgReNew(AlnMsgPtr amp);

/* used by AlnMgrMakeMultSegments */
typedef struct am_tinyinfo {
   Int4   start;
   Int4   stop;
   Uint4  numgap;
   Uint2  which;
   Int4   numsap;
} AMTinyInfo, PNTR AMTinyInfoPtr;

typedef struct am_aligninfo {
   SeqAlignPtr  align;
   Int4 align_len;
} AMAlignInfo, PNTR AMAlignInfoPtr;


/* used by AlnMgrMakeSegmentedMasterSlave */
typedef struct am_msms {
   Int4         start;
   Int4         stop;
   Int4         sstart;
   Int4         sstop;
   Uint1        strand;
   SeqIdPtr     sip;
   SeqAlignPtr  sap;
   Int4         nsap;
   Int4         n;
   Int4         count;
   Int4         masternum;
   struct am_msms PNTR next;
} AMmsms, PNTR AMmsmsPtr;

typedef struct am_tinymsinfo {
   Int4       start;
   AMmsmsPtr  ams;
} AMTinyMSInfo, PNTR AMTinyMSInfoPtr;

typedef struct am_siplist {
   SeqIdPtr  sip;
   Int4      first_row;
   struct am_siplist PNTR next;
} AMsiplist, PNTR AMsiplistPtr;

NLM_EXTERN SeqAlignIndexPtr SeqAlignIndexNew(void);

/********************************************************************************
*
*  AlnMgrIndexSingleSeqAlign indexes (in place) only the first seqalign or
*  seqalign set in the chain that is passed in.  It will extensively
*  rearrange the first seqalign given.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrIndexSingleSeqAlign(SeqAlignPtr sap);

NLM_EXTERN Boolean AlnMgrIndexSingleChildSeqAlign(SeqAlignPtr sap);

/********************************************************************************
*
*  AlnMgrReIndexSeqAlign frees the parent indexes, indexes any child
*  seqaligns that are not indexed (it assumes that any indexed child 
*  seqaligns are correctly indexed), and reindexes the set.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrReIndexSeqAlign(SeqAlignPtr sap);

/********************************************************************************
*
*  AlnMgrIndexSeqAlign indexes (in place) the ENTIRE chain of seqaligns
*  and seqalign sets passed in, and extensively rearranges the seqalign.   
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrIndexSeqAlign(SeqAlignPtr sap);


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

NLM_EXTERN Boolean AlnMgrAnythingToSeg (SeqAlignPtr sap);

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
NLM_EXTERN Boolean AlnMgrIndexLinkedSegs (SeqAlignPtr sap);
NLM_EXTERN Boolean AlnMgrIndexParentSA(SeqAlignPtr sap);
NLM_EXTERN SeqIdPtr AlnMgrPropagateUpSeqIdPtrs(SeqAlignPtr sap, Int4Ptr num);
NLM_EXTERN SeqIdPtr  AlnMgrPropagateSeqIdsBySapList(AMAlignIndexPtr amaip);
NLM_EXTERN SeqIdPtr  AlnMgrPropagateSeqIdsByRow(AMAlignIndexPtr amaip);

NLM_EXTERN Boolean AlnMgrMakeMultipleByScore(SeqAlignPtr sap);


/**********************************************************************
*
*  am_print_seqalign_indexes is a useful function for debugging; 
*  given a parent or child seqalign it prints a text summary of the
*  indexes for that seqalign.
*
***********************************************************************/
NLM_EXTERN void am_print_seqalign_indexes(SeqAlignPtr sap);


/**********************************************************************
*
*  AlnMgrCheckAlignForParent is called by almost every other function;
*  it checks whether the parents and children are correctly labeled
*  and indexed, and it returns AM_CHILD to indicate that the sap 
*  passed in is a child, AM_PARENT for a parent, or -1 for an error.
*  If the indexes aren't present, it will try to create them.
*
**********************************************************************/
NLM_EXTERN Int4 AlnMgrCheckAlignForParent(SeqAlignPtr sap);


/***********************************************************************
*
*  AlnMgrSortSeqAligns is a variant of the ValNodeSort function, and
*  calls very similar heapsort functions.  It can take a comparison
*  function that needs userdata, so more specific sorts are possible
*  without defining special structures for every type of sort.
*
***********************************************************************/
NLM_EXTERN SeqAlignPtr PNTR AlnMgrSortSeqAligns (SeqAlignPtr sap, int (LIBCALLBACK *compare)(VoidPtr, VoidPtr, VoidPtr), VoidPtr userdata, Int4Ptr numsap);

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
NLM_EXTERN int LIBCALLBACK AlnMgrCompareIncreasingBySeqIdPtr (VoidPtr base, VoidPtr large_son, VoidPtr userdata);

/*********************************************************************
*
*  AlnMgrFindFirst is crucial to the AlnMgrMakeFakeMultiple function;
*  it uses the sorted order of the seqaligns in each AMAlignDat
*  structure to guide a heapsort of all the seqaligns.
*
*********************************************************************/
NLM_EXTERN int LIBCALLBACK AlnMgrFindFirst(VoidPtr base, VoidPtr large_son, VoidPtr userdata);
NLM_EXTERN int LIBCALLBACK AlnMgrCompareTips(VoidPtr base, VoidPtr large_son);


/************************************************************************
*
*  AlnMgrGetNextLengthBit should be called iteratively on an alignment
*  to return the lengths of all the aligned and unaligned pieces in
*  the alignment.  Don't change the value in r, just pass in a pointer
*  to an allocated Int4 set to 0 initially.  The lengths of the unaligned
*  pieces are precomputed using AlnMgrGetMaxUnalignedLength; if no
*  precomputed values are found, this function is used to compute the
*  lengths on the fly.
*
************************************************************************/
NLM_EXTERN Boolean AlnMgrGetNextLengthBit(SeqAlignPtr sap, Int4Ptr length, Int4Ptr r);
NLM_EXTERN Int4 AlnMgrGetMaxUnalignedLength(SeqAlignPtr sap1, SeqAlignPtr sap2);


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
NLM_EXTERN Boolean AlnMgrGetNextAlnBit(SeqAlignPtr sap, AlnMsgPtr amp);

/********************************************************************************
*
*  binary search functions
*
********************************************************************************/
NLM_EXTERN Uint4 binary_search_on_uint4_list(Uint4Ptr list, Uint4 pos, Uint4 listlen);
NLM_EXTERN Int4 binary_search_on_uint2_list(Uint2Ptr list, Uint2 ele, Uint2 listlen);
NLM_EXTERN Int4 binary_search_by_chunk(Int4Ptr list, Int4 ele, Int4 listlen, Int4 chunksize, Int4 offset);
NLM_EXTERN Int4 binary_search_segment_array(SASeqDatPtr ssdp, Int4 pos, Int4 numseq, Int4 offset, DenseSegPtr dsp);


/************************************************************************
*
*  These are several utility functions which get needed data from the
*  indexes.  "N" refers to row number (using a 1-based numbering system).
*
************************************************************************/
NLM_EXTERN Int4 AlnMgrGetAlnLength(SeqAlignPtr sap, Boolean fill_in);
NLM_EXTERN Int4 AlnMgrGetNumSeqs(SeqAlignPtr sap);
NLM_EXTERN SeqIdPtr AlnMgrGetUniqueSeqs(SeqAlignPtr sap, Int4Ptr n);
NLM_EXTERN SeqIdPtr AlnMgrGetNthSeqIdPtr(SeqAlignPtr sap, Int4 n);
NLM_EXTERN void AlnMgrGetNthSeqRangeInSA(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop);
NLM_EXTERN Int4 AlnMgrGetNumSegments(SeqAlignPtr sap);

/********************************************************************************
*
*  AlnMgrGetNextNthSeqRange is called recursively to return the lengths of
*  all aligned and all internal unaligned regions of any row in a seqalign.
*  If there is an error, or if the function is called past the last block,
*  the function returns FALSE.  Set where to point to an allocated integer
*  equal to 0 on the first call and don't change it during the loop.  Set
*  the boolean unaligned to FALSE to get only the aligned regions, and TRUE to
*  get the aligned regions plus all internal unaligned regions.  For unaligned
*  regions, *is_aligned will be FALSE.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrGetNextNthSeqRange(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop, Int4Ptr where, BoolPtr is_aligned, Boolean unaligned);


/********************************************************************************
*
*  AlnMgrGetNthRowTail retrieves the blocks of sequence on either end of the
*  alignment, by row.  which_tail is LEFT_TAIL to retrieve the ends which come
*  before alignment coordinate 0, and RIGHT_TAIL to retrieve the other ends.
*  The function returns TRUE if successful, FALSE for an error.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrGetNthRowTail(SeqAlignPtr sap, Int4 n, Uint1 which_tail, Int4Ptr start, Int4Ptr stop, Uint1Ptr strand);
NLM_EXTERN Boolean AlnMgrGetNthUnalignedForNthRow(SeqAlignPtr sap, Int4 unaligned, Int4 row, Int4Ptr start, Int4Ptr stop);
NLM_EXTERN Uint1 AlnMgrGetStrand(SeqAlignPtr sap, SeqIdPtr sip);
NLM_EXTERN Uint1 AlnMgrGetNthStrand(SeqAlignPtr sap, Int4 n);
NLM_EXTERN Int4 AlnMgrGetNForSip(SeqAlignPtr sap, SeqIdPtr sip);
NLM_EXTERN Int4 AlnMgrGetNForSap(AMAlignIndexPtr amaip, SeqAlignPtr sap);

/********************************************************************************
*
*  AlnMgrGetAllNForSip is called in a while loop to return all the rows that a
*  seqid appears in in a given seqalign.  Use n = 0 to start, and then on
*  return, if the return is TRUE, n will be the row number of the next row
*  that the seqid appears in.  If the return is FALSE, either there was an
*  error or there are no (more) rows containing that seqid.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrGetAllNForSip(SeqAlignPtr sap, SeqIdPtr sip, Int4Ptr n);

NLM_EXTERN Int4 AlnMgrGetSapForSip(AMAlignIndexPtr amaip, SeqIdPtr sip, Int4 which);
NLM_EXTERN Int4 AlnMgrMapToBsqCoords(SeqAlignPtr sap, Uint4 pos, SeqIdPtr sip, SeqIdPtr master); 

/********************************************************************************
*
*  AlnMgrMapRowCoords maps a position in a given row to the bioseq coordinate
*  of that row.  If master is NULL, the alignment is taken to be flattened;
*  otherwise it is an alignment according to that master (this will change the
*  correspondence between row coordinates and bioseq coordinates).  The return
*  value will be either a positive bioseq coordinate, or -1 if the bioseq is
*  gapped at that row position.
*
********************************************************************************/
NLM_EXTERN Int4 AlnMgrMapRowCoords(SeqAlignPtr sap, Uint4 pos, Int4 row, SeqIdPtr master);

/********************************************************************************
*
*  AlnMgrMapBioseqToSeqAlign takes a position in bioseq coordinates in a
*  row and maps it to seqalign coordinates, using the given master as
*  the alignment master (if master is NULL the alignment is flat).  A
*  return value of -1 indicates an error; a return value of -2 indicates
*  that the given bioseq coordinates are not contained in the alignment
*  specified.
*
********************************************************************************/
NLM_EXTERN Int4 AlnMgrMapBioseqToSeqAlign(SeqAlignPtr sap, Uint4 pos, Int4 row_num, SeqIdPtr master);


/***********************************************************************
*
*  AlnMgrMakeFakeMultiple calls AlnMgrCheckOverlapping to decide whether
*  an alignment is linear.  Then, if possible, it calls AlnMgrMakeAlignCoords
*  to create alignment coordinates across all children contained in the
*  parent.
*
***********************************************************************/
NLM_EXTERN Boolean AlnMgrMakeFakeMultiple(SeqAlignPtr sap);
NLM_EXTERN void AlnMgrMakeAlignCoords(SeqAlignPtr sap);
NLM_EXTERN Boolean AlnMgrMergeIntoMSMultByMaster(AMAlignIndexPtr amaip, Int4Ptr lens, Uint4Ptr numseg);
NLM_EXTERN Boolean AlnMgrMergeSegments(Int4Ptr lens, SeqAlignPtr sap, SeqIdPtr master, Int4Ptr where, Int4 which);
NLM_EXTERN Boolean AlnMgrFillInStarts(SeqAlignPtr PNTR saparray, Int4Ptr starts, Int4 numseg, Int4Ptr lens, Int4 numsaps, Uint4Ptr aligncoords);
NLM_EXTERN Int4 AlnMgrGetStartFromMaster(SeqAlignPtr sap, Int4 start);
NLM_EXTERN Uint4 AlnMgrGetMasterGapStartForSeg(SeqAlignPtr sap, Int4 which_gap, Uint4Ptr aligncoord);
NLM_EXTERN Boolean AlnMgrReconcileGaps(Int4Ptr lens, Uint4Ptr aligncoords, Int4 num);
NLM_EXTERN Boolean AlnMgrMakeMultSegments(AMAlignIndexPtr amaip);

NLM_EXTERN Int4 AlnMgrCheckOverlapping(SeqAlignPtr sap);

/******************************************************************************
*
*  AlnMgrGetMaxSegments simply adds up the number of segments for each
*  SeqAlign in a linked list, to get the maximum number of segments
*  for the merge of the list (for memory allocation in AlnMgrMakeFakeMultiple).
*
*******************************************************************************/
NLM_EXTERN Int4 AlnMgrGetMaxSegments(SeqAlignPtr sap);

/*******************************************************************************
*
*  Row Management functions:
*
*******************************************************************************/
NLM_EXTERN Int4 AlnMgrGetNumRows(SeqAlignPtr sap);
NLM_EXTERN Int4 AlnMgrGetMaxRowsForParentPartial(SeqAlignPtr sap);
NLM_EXTERN Boolean AlnMgrGetRowsForPartial(SeqAlignPtr sap);
NLM_EXTERN Boolean AlnMgrGetRowsForMasterSlave(SeqAlignPtr sap);


/*******************************************************************************
*
*  AlnMgrFindMaster returns the SeqIdPtr of the first bioseq present in all
*  child alignments, unless the sap->master field is set in the child alignments,
*  in which case that SeqIdPtr is returned (if it's the same in all children).
*
*******************************************************************************/
NLM_EXTERN SeqIdPtr AlnMgrFindMaster(SeqAlignPtr sap);

/*******************************************************************************
*
*  AlnMgrCheckRealMaster makes sure that the master seqid given appears
*  once and only once in each seqalign in the set if a parent is given,
*  or once and only one in the seqalign if a child is given.
*
*******************************************************************************/
NLM_EXTERN Boolean AlnMgrCheckRealMaster(SeqAlignPtr sap, SeqIdPtr master);

NLM_EXTERN Boolean AlnMgrMakeSegmentedMasterSlave(SeqAlignPtr sap);
NLM_EXTERN int LIBCALLBACK AlnMgrCompareAMS(VoidPtr base, VoidPtr large_son);
NLM_EXTERN int LIBCALLBACK AlnMgrCompareMasterAMS(VoidPtr base, VoidPtr large_son);

NLM_EXTERN Boolean AlnMgrForceMasterSlave(SeqAlignPtr sap);

NLM_EXTERN void AlnMgrSetMaster(SeqAlignPtr sap, SeqIdPtr master);
NLM_EXTERN void AlnMgrMakeMasterPlus(SeqAlignPtr sap);

NLM_EXTERN SeqAlignPtr AlnMgrGetSubAlign(SeqAlignPtr sap, SeqIdPtr which_master, Int4 from, Int4 to);

/********************************************************************************
*
*   viewer and editor management functions  
* 
********************************************************************************/ 

NLM_EXTERN SeqAlignPtr AlnMgrCopyIndexedParentSeqAlign(SeqAlignPtr sap);
NLM_EXTERN RowSourcePtr AlnMgrCopyRowSource(RowSourcePtr rsp);
NLM_EXTERN AMAlignDatPtr AlnMgrCopyamadp(AMAlignDatPtr amadp, SeqAlignPtr sap_tmp, SeqAlignPtr seg_head);
NLM_EXTERN SeqAlignIndexPtr AlnMgrCopyIndexesForChildSeqAlign(SeqAlignPtr sap);
NLM_EXTERN SASeqDatPtr  AlnMgrCopySASeqDat(SASeqDatPtr ssdp);
NLM_EXTERN SeqAlignPtr AlnMgrCopyAndIndexSingleAlignment(SeqAlignPtr sap);
NLM_EXTERN Boolean AlnMgrCopyIndexedParentIntoSap(SeqAlignPtr sap, SeqAlignPtr target);
NLM_EXTERN Boolean AlnMgrDeleteChildByPointer(SeqAlignPtr parent, SeqAlignPtr child);
NLM_EXTERN void AlnMgrDeleteRow(SeqAlignPtr sap, Int4 row);

/********************************************************************************
*
*   SeqAlign types  
* 
********************************************************************************/ 
NLM_EXTERN  Boolean AlnMgrIsSAPDiscAli(SeqAlignPtr sap);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif