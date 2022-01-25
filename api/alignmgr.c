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
* File Name:  alignmgr.c
*
* Author:  Sarah Wheelan
*
* Version Creation Date:   7/99
*
* $Revision: 6.56 $
*
* File Description: SeqAlign indexing and messaging functions
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: alignmgr.c,v $
* Revision 6.56  2000/01/19 15:45:09  wheelan
* many, many bug fixes in AlnMgrGetSubAlign and AlnMgrGetNextAlnBit
*
* Revision 6.55  2000/01/14 18:50:36  wheelan
* fixed bug in AlnMgrGetSubAlign
*
* Revision 6.54  2000/01/12 17:43:19  wheelan
* added AlnMgrGetNumSegments, AlnMgrDeleteRow
*
* Revision 6.53  1999/12/02 20:31:59  lewisg
* put seqentries into bioseqset and fix calling convention in alignmgr.c
*
* Revision 6.52  1999/11/30 14:36:39  wheelan
* added AlnMgrMakeMultipleByScore; bug fixes
*
* Revision 6.51  1999/11/26 15:42:19  vakatov
* Fixed for the C++ and/or MSVC DLL compilation
*
* Revision 6.50  1999/11/24 11:29:52  wheelan
* added missing return values
*
* Revision 6.49  1999/11/18 19:30:33  wheelan
* added AlnMgrDeleteChildByPointer, bug fixes
*
* Revision 6.48  1999/11/03 12:47:05  wheelan
* added code to correctly handle internal gaps in segmented master-slave alignments
*
* Revision 6.47  1999/11/02 12:38:38  wheelan
* bug fixes when only one child
*
* Revision 6.46  1999/10/25 18:17:23  wheelan
* Added AlnMgrGetUniqueSeqs, fixed merge function to handle single child seqalign correctly
*
* Revision 6.45  1999/10/19 19:27:03  wheelan
* added static defines; changed behavior of AlnMgrGetNextNthSeqRange; rewrote AlnMgrMakeSegmentedMasterSlave to handle more cases
*
* Revision 6.44  1999/10/15 21:51:02  durand
* add AlnMgrIsSAPDiscAli()
*
* Revision 6.43  1999/10/15 18:19:05  wheelan
* added rudimentary ability to default to master-slave type if possible
*
* Revision 6.42  1999/10/15 13:48:47  wheelan
* added AlnMgrGetNthRowTail, extended capability of AlnMgrGetNthStrand
*
* Revision 6.41  1999/10/14 16:10:30  kans
* new includes and prototypes added
*
* Revision 6.40  1999/10/13 19:29:03  wheelan
* added speedup for segmented master-slave creation
*
* Revision 6.39  1999/10/07 13:37:16  wheelan
* added AlnMgrIndexSingleSeqAlign, which only indexes the first seqalign in a list; also added automatic computation of max length of unaligned regions for time savings
*
* Revision 6.38  1999/10/06 19:35:09  wheelan
* added several viewer and editor management functions; fixed many bugs in AlnMgrGetNextAlnBit
*
* Revision 6.37  1999/10/05 15:15:31  wheelan
* added AlnMgrGetNthUnalignedForNthRow
*
* Revision 6.36  1999/10/05 14:02:31  wheelan
* bug fixes in AlnMgrGetNextAlnBit
*
* Revision 6.35  1999/10/04 14:58:08  wheelan
* bug fixes; added AlnMgrMapBioseqToSeqAlign
*
* Revision 6.34  1999/09/24 15:04:55  lewisg
* AlnMgrGetNextAlnBit: amp->to_m changed when calling child
*
* Revision 6.33  1999/09/24 14:29:58  wheelan
* changed behavior of AlnMgrGetNextLengthBit to mimic other GetNext functions, completed functionality of AlnMgrGetSubAlign, bug fixes
*
* Revision 6.32  1999/09/23 16:03:32  wheelan
* Added structures and functions to support segmented master-slave alignments
*
* Revision 6.31  1999/09/22 13:19:15  wheelan
* made AlnMsg row_num field 1-based, added AlnMgrGetNextNthSeqRange, started adding functions to handle a segmented master-slave alignment
*
* Revision 6.30  1999/09/21 19:15:28  wheelan
* changed AlnMgrGetNextAlnBit to return FALSE if called once more past the end; various bug fixes; implemented part of AlnMgrGetSubAlign
*
* Revision 6.29  1999/09/20 12:12:58  wheelan
* added safety checks in case input seqalign has no strand or score information
*
* Revision 6.28  1999/09/20 11:58:52  wheelan
* modified AlnMgrGetNthSeqRange to use new row information structures
*
* Revision 6.27  1999/09/17 16:55:33  wheelan
* bug fixes, added AlnMgrPropagateSeqIdsBySapList to correctly associate seqids with rows
*
* Revision 6.26  1999/09/14 15:48:50  kans
* AlnMgrMapRowCoords returns -1 on failure at end of function
*
* Revision 6.25  1999/09/13 19:57:10  sicotte
* Make AlnMgrMapBsqCoord work for continous alignments
*
* Revision 6.24  1999/09/13 19:43:09  sicotte
* bug fixes
*
* Revision 6.23  1999/09/13 14:33:24  wheelan
* added support for row numbers in AlnMgrGetNextAlnBit
*
* Revision 6.22  1999/09/08 13:36:16  wheelan
* fixed bugs found by Patrick Durand
*
* Revision 6.21  1999/09/08 11:55:35  sicotte
* fix bug that was missing end segments
*
* Revision 6.20  1999/09/08 11:49:13  wheelan
* added capability to return length of unaligned regions
*
* Revision 6.19  1999/09/07 12:11:17  wheelan
* fixed bugs pointed out by Hugues
*
* Revision 6.18  1999/09/06 16:37:44  wheelan
* added AlnMgrGetNextLengthBit and associated function
*
* Revision 6.17  1999/09/06 15:55:55  wheelan
* IndexSeqAlign now makes the fake multiple if possible
*
* Revision 6.16  1999/09/06 15:52:25  wheelan
* added row management functions, made most functions minus-strand compliant, added smarter test for master-slave vs partial
*
* Revision 6.15  1999/09/01 20:11:56  wheelan
* added new merge function and the typedef for the structure it uses
*
* Revision 6.14  1999/09/01 14:40:06  wheelan
* added AlnMgrGetStrand, fixed bugs in GetNextAlnBit, added more cases to AlnMgrIndexSeqAlign
*
* Revision 6.13  1999/08/30 19:28:06  wheelan
* modified AlnMgrGetNextAlnBit to handle master-slave alignments
*
* Revision 6.12  1999/08/26 20:35:21  wheelan
* added parent indexing and pairwise-to-multiple functions
*
* Revision 6.11  1999/08/20 11:23:53  wheelan
* fixed AlnMgrGetNthSeqRange for minus strands
*
* Revision 6.10  1999/08/19 19:30:26  wheelan
* made case for SAT_PARTIAL in AlnMgrGetNextAlnBit
*
* Revision 6.9  1999/08/19 17:24:50  wheelan
* changed AMAlignIndex structure, added more api functions
*
* Revision 6.8  1999/08/12 20:56:56  vakatov
* [WIN32] Added missed LIBCALLBACK
*
* Revision 6.7  1999/08/12 12:41:53  wheelan
* added comments, and functions to index the parent
*
* Revision 6.6  1999/08/06 18:31:19  wheelan
* fixed compiler error
*
* Revision 6.5  1999/08/06 16:38:43  kans
* fixed Mac compiler complaints
*
* Revision 6.4  1999/08/06 13:44:14  wheelan
* added several functions; changed all function names to AlnMgr..
*
* Revision 6.3  1999/07/30 14:17:52  wheelan
* fixes to keep Mac compiler happy
*
* Revision 6.2  1999/07/30 14:08:37  wheelan
* added api functions to access indexes
*
* Revision 6.1  1999/07/29 12:56:25  wheelan
* initial checkin
*

* ==========================================================================
*/



#include <alignmgr.h>
#include <samutil.h>

/***************************************************************************
*
*  static functions
*
***************************************************************************/
static void heapsort_with_userdata (VoidPtr b, size_t nel, size_t width, int (LIBCALLBACK *compar)PROTO((VoidPtr, VoidPtr, VoidPtr)), VoidPtr userdata);
static void heapify_with_userdata(CharPtr base0, CharPtr base, CharPtr lim, CharPtr last, size_t width, int(LIBCALLBACK *compar)PROTO((VoidPtr, VoidPtr, VoidPtr)), VoidPtr userdata);
static Boolean am_get_nth_range_for_partial(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop, Int4Ptr where, BoolPtr is_aligned, Boolean unaligned);
static AMmsmsPtr am_sort_ammsms(AMmsmsPtr ams_head, Int4 n);
static AMmsmsPtr am_sort_masterams(AMmsmsPtr ams_head, Int4 n);
static Int4 am_get_first_rsp_for_sip(SeqIdPtr sip, AMsiplistPtr siplist);
static int LIBCALLBACK AMCompareAlignInfoProc(VoidPtr ptr1, VoidPtr ptr2);
static Int4 AlnMgrMapSegmentCoords(SeqAlignPtr sap, Uint4 pos, Int4 row, SeqIdPtr master, Int4Ptr len);



/*******************************************************************
*
*  all the memory allocation/deallocation functions
*
*******************************************************************/

NLM_EXTERN SeqAlignIndexPtr SeqAlignIndexNew(void)
{
   return (SeqAlignIndexPtr)(MemNew(sizeof(SeqAlignIndex)));
}

static Boolean LIBCALLBACK SAIndexFreeFunc(VoidPtr index)
{
    return SAIndexFree(index);
}

NLM_EXTERN SAIndexPtr SAIndexNew(void)
{
   SAIndexPtr  saip;

   saip = (SAIndexPtr)MemNew(sizeof(SAIndex));
   saip->master = -1;
   saip->freefunc = (SeqAlignIndexFreeFunc)(SAIndexFreeFunc);
   return saip;
}

NLM_EXTERN Boolean SAIndexFree(VoidPtr index)
{
   Boolean     retval;
   SAIndexPtr  saip;

   retval = FALSE;
   if (!index)
      return retval;
   saip = (SAIndexPtr)index;
   if (saip->indextype != INDEX_SEGS)
      return retval;
   MemFree(saip->aligncoords);
   MemFree(saip->ssdp);
   retval = TRUE;
   return retval;
}

NLM_EXTERN SASeqDatPtr SASeqDatNew(void)
{
   return (SASeqDatPtr)(MemNew(sizeof(SASeqDat)));
}

NLM_EXTERN RowSourcePtr RowSourceNew(void)
{
   return (RowSourcePtr)(MemNew(sizeof(RowSource)));
}

NLM_EXTERN RowSourcePtr RowSourceFree(RowSourcePtr rsp)
{
   if (rsp == NULL)
      return NULL;
   rsp->id = SeqIdSetFree(rsp->id);
   MemFree(rsp->which_saps);
   MemFree(rsp->num_in_sap);
   return NULL;
}

static Boolean LIBCALLBACK AMAlignIndexFreeFunc (VoidPtr data)
{
    return AMAlignIndexFree(data);
}


NLM_EXTERN AMAlignIndexPtr AMAlignIndexNew(void)
{
   AMAlignIndexPtr  amaip;

   amaip = (AMAlignIndexPtr)MemNew(sizeof(AMAlignIndex));
   amaip->freefunc = (SeqAlignIndexFreeFunc)(AMAlignIndexFreeFunc);
   amaip->master = -2;
   return amaip;
}

NLM_EXTERN Boolean AMAlignIndexFree(VoidPtr index)
{
   AMAlignIndexPtr  amaip;         
   Int4             i;
   Boolean          retval;

   retval = FALSE;
   amaip = (AMAlignIndexPtr)(index);
   if (!amaip)
      return retval;
   if (amaip->indextype != INDEX_PARENT)
      return retval;
   amaip->ids = SeqIdSetFree(amaip->ids);
   for (i=0; i<(amaip->numbsqs); i++)
   {
      amaip->amadp[i] = AMAlignDatFree(amaip->amadp[i]);
   }
   MemFree(amaip->amadp);
   MemFree(amaip->aligncoords);
   MemFree(amaip->lens);
   MemFree(amaip->ulens);
   MemFree(amaip->starts);
   for (i=0; i<(amaip->numrows); i++)
   {
      amaip->rowsource[i] = RowSourceFree(amaip->rowsource[i]);
   }
   MemFree(amaip->rowsource);
   retval = TRUE;
   return retval;
}

NLM_EXTERN AMAlignDatPtr AMAlignDatNew(void)
{
   return (AMAlignDatPtr)(MemNew(sizeof(AMAlignDat)));
}

NLM_EXTERN AMAlignDatPtr AMAlignDatFree(AMAlignDatPtr amadp)
{
   if (amadp == NULL)
      return NULL;
   SeqIdFree(amadp->sip);
   MemFree(amadp->saps);
   MemFree(amadp->segments);
   return NULL;
}

NLM_EXTERN AlnMsgPtr AlnMsgNew(void)
{
   AlnMsgPtr  amp;

   amp = (AlnMsgPtr)MemNew(sizeof(AlnMsg));
   amp->send_space = FALSE;
   amp->row_num = -1;
   amp->prev = -2;
   amp->prev_sap = -2;
   amp->place = 0;
   return amp;
}

NLM_EXTERN AlnMsgPtr AlnMsgReNew(AlnMsgPtr amp)
{
   amp = (AlnMsgPtr)MemNew(sizeof(AlnMsg));
   amp->send_space = FALSE;
   amp->row_num = -1;
   amp->prev = -2;
   amp->prev_sap = -2;
   amp->place = 0;
   return amp;
}

/********************************************************************************
*
*  AlnMgrIndexSingleSeqAlign indexes (in place) only the first seqalign or
*  seqalign set in the chain that is passed in.  It will extensively
*  rearrange the first seqalign given.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrIndexSingleSeqAlign(SeqAlignPtr sap)
{
   SeqAlignPtr  sap_next;

   if (sap == NULL)
      return TRUE;
   sap_next = NULL;
   if (sap->next)
      sap_next = sap->next;
   sap->next = NULL;
   AlnMgrIndexSeqAlign(sap);
   sap->next = sap_next;
   if (sap->saip)
      return TRUE;
   else
      return FALSE;
}

NLM_EXTERN Boolean AlnMgrIndexSingleChildSeqAlign(SeqAlignPtr sap)
{
   SeqAlignPtr sap_next;

   if (sap->segtype == SAS_DISC)
      return FALSE;
   sap_next = NULL;
   if (sap->next)
      sap_next = sap->next;
   sap->next = NULL;
   if (sap->saip != NULL)
   {
      if (sap->saip->indextype == INDEX_SEGS)
         SAIndexFree(sap->saip);
   }
   if (sap->segtype == SAS_DENSEG)
      AlnMgrIndexLinkedSegs(sap);
   else if (sap->segtype == SAS_DENDIAG)
      AlnMgrIndexSingleSeqAlign(sap);
   sap->next = sap_next;
   if (sap->saip)
      return TRUE;
   else
      return FALSE;
}

/********************************************************************************
*
*  AlnMgrReIndexSeqAlign frees the parent indexes, indexes any child
*  seqaligns that are not indexed (it assumes that any indexed child 
*  seqaligns are correctly indexed), and reindexes the set.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrReIndexSeqAlign(SeqAlignPtr sap)
{
   SeqAlignPtr      sap_tmp;
   SeqAlignPtr      tmp_next;

   if (sap == NULL)
      return FALSE;
   if (sap->segtype != SAS_DISC) /*  we don't know what we're dealing with */
      return FALSE;
   if (!AMAlignIndexFree((Pointer)sap->saip))
      return FALSE;
   sap->saip = NULL;
   sap_tmp = (SeqAlignPtr)sap->segs;
   while (sap_tmp)
   {
      if (sap_tmp->saip == NULL)
      {
         tmp_next = sap_tmp->next;
         sap_tmp->next = NULL;
         if (!AlnMgrIndexLinkedSegs(sap_tmp))
            return FALSE;
         sap_tmp->next = tmp_next;
      }
      sap_tmp = sap_tmp->next;
   }
   if (!AlnMgrIndexParentSA(sap))
      return FALSE;
   if (!AlnMgrMakeFakeMultiple(sap))
      return FALSE;
   return TRUE;
}

/********************************************************************************
*
*  AlnMgrIndexSeqAlign indexes (in place) the ENTIRE chain of seqaligns
*  and seqalign sets passed in, and extensively rearranges the seqalign.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrIndexSeqAlign(SeqAlignPtr sap)
{
   Boolean      retval;
   SeqAlignPtr  salp;
   SeqAlignPtr  salp_head;
   SeqAlignPtr  salp_prev;
   SeqAlignPtr  sap_tmp;

   retval = FALSE;
   if (!sap)
      return retval;
   if (sap->saip != NULL)
   {
      return TRUE;
   }
   if (sap->segtype == SAS_DISC)
   {
      if (sap->next == NULL)
      {
         if (!AlnMgrIndexLinkedSegs((SeqAlignPtr)(sap->segs)))
            return retval;
         if (!AlnMgrIndexParentSA(sap))
            return retval;
         if (!AlnMgrMakeFakeMultiple(sap))
            return retval;
      } else
      {
         sap_tmp = sap;
         salp_head = salp_prev = NULL;
         while (sap_tmp)
         {
            if (sap_tmp->segtype == SAS_DISC)
            {
               salp = (SeqAlignPtr)sap_tmp->segs;
               if (salp_head)
               {
                  salp_prev->next = salp;
                  while (salp->next)
                  {
                     salp = salp->next;
                  }
                  salp_prev = salp;
               } else
               {
                  salp_head = salp;
                  while (salp->next)
                  {
                     salp = salp->next;
                  }
                  salp_prev = salp;
               }
               sap_tmp = sap_tmp->next;
            } else if (sap_tmp->segtype == SAS_DENDIAG || sap_tmp->segtype == SAS_DENSEG)
            {
               salp = sap_tmp->next;
               sap_tmp->next = NULL;
               if (salp_head)
               {
                  salp_prev->next = sap_tmp;
                  salp_prev = sap_tmp;
               } else
               {
                  salp_head = salp_prev = sap_tmp;
               }
               sap_tmp = salp;
            }
         }
         sap->next = NULL;
         sap->segs = (Pointer)salp_head;
         if (!AlnMgrIndexLinkedSegs((SeqAlignPtr)(sap->segs)))
            return retval;
         if (!AlnMgrIndexParentSA(sap))
            return retval;
         if (!AlnMgrMakeFakeMultiple(sap))
            return retval;
      }
   } else
   {
      if (!AlnMgrAnythingToSeg(sap))
         return retval;
      if (!AlnMgrIndexLinkedSegs((SeqAlignPtr)(sap->segs)))
         return retval;
      if (!AlnMgrIndexParentSA(sap))
         return retval;
      if (!AlnMgrMakeFakeMultiple(sap))
         return retval;
   }
   retval = TRUE;
   return retval;
}


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
NLM_EXTERN Boolean AlnMgrAnythingToSeg (SeqAlignPtr sap)
{
   DenseDiagPtr  ddp;
   DenseDiagPtr  ddp_prev;
   DenseSegPtr   dsp;
   Int4          i;
   Boolean       retval;
   SeqAlignPtr   salp;
   SeqAlignPtr   salp_tmp;
   SeqAlignPtr   sap_head;
   SeqAlignPtr   sap_new;
   SeqAlignPtr   sap_prev;
   StdSegPtr     ssp;
   StdSegPtr     ssp_next;

   retval = FALSE;
   if (!sap)
      return retval;
   sap_new = SeqAlignNew();
   sap_new->type = SAT_GLOBAL;
   sap_new->segtype = sap->segtype;
   sap_new->dim = sap->dim;
   sap_new->segs = sap->segs;
   sap_new->master = sap->master;
   sap_new->bounds = sap->bounds;
   sap_new->next = sap->next;
   sap_new->score = sap->score;
   sap->next = NULL;
   sap->segtype = SAS_DISC;
   sap->type = 0;
   sap->dim = 0;
   sap->master = NULL;
   sap->bounds = NULL;
   sap->score = NULL;
   salp = sap_new;
   sap_head = sap_prev = NULL;
   while (salp)
   {
      if (salp->segtype < 1)
      {
         return retval;
      } else if (salp->segtype == SAS_DENDIAG)
      {
         ddp = (DenseDiagPtr)salp->segs;
         while (ddp)
         {
            sap_new = SeqAlignNew();
            sap_new->type = SAT_GLOBAL;
            sap_new->segtype = SAS_DENSEG;
            sap_new->dim = ddp->dim;
            dsp = DenseSegNew();
            dsp->dim = sap_new->dim;
            dsp->numseg = 1;
            dsp->starts = ddp->starts;
            ddp->starts = NULL;
            dsp->lens = (Int4Ptr)MemNew(sizeof(Int4));
            dsp->lens[0] = ddp->len;
            ddp->len = 0;
            dsp->scores = ddp->scores;
            ddp->scores = NULL;
            dsp->strands = ddp->strands;
            ddp->strands = NULL;
            if (dsp->strands == NULL)
            {
               dsp->strands = (Uint1Ptr)MemNew(dsp->dim * sizeof(Uint1));
               for (i=0; i<dsp->dim; i++)
               {
                  dsp->strands[i] = Seq_strand_plus;
               }
            }
            dsp->ids = SeqIdDupList(ddp->id);
            sap_new->segs = (Pointer)dsp;
            if (dsp->scores)
               sap_new->score = ScoreDup(dsp->scores);
            if (!sap_head)
            {
               sap_head = sap_prev = sap_new;
            } else
            {
               sap_prev->next = sap_new;
               sap_prev = sap_new;
            }
            ddp_prev = ddp;
            ddp = ddp->next;
            DenseDiagFree(ddp_prev);
         }
         salp_tmp = salp->next;
         sap_prev->next = salp_tmp;
         salp = salp_tmp;
         retval = TRUE;
      } else if (salp->segtype == SAS_DENSEG)
      {
         if (!sap_head)
            sap_head = sap_prev = salp;
         else
         {
            sap_prev->next = salp;
            sap_prev = salp;
         }
         dsp = (DenseSegPtr)salp->segs;
         if (dsp->strands == NULL)
         {
            dsp->strands = (Uint1Ptr)MemNew((dsp->dim)*(dsp->numseg)* sizeof(Uint1));
            for (i=0; i<(dsp->dim)*(dsp->numseg); i++)
            {
               dsp->strands[i] = Seq_strand_plus;
            }
         }
         salp = salp->next;
         retval = TRUE;
      } else if (salp->segtype == SAS_STD)
      {
         sap_prev = sap_head = NULL;
         ssp = (StdSegPtr)salp->segs;
         while (ssp)
         {
            sap_new = SeqAlignNew();
            if (sap_head)
            {
               sap_prev->next = sap_new;
               sap_prev = sap_new;
            } else
            {
               sap_head = sap_prev = sap_new;
            }
            sap_new->segtype = SAS_STD;
            sap_new->type = SAT_GLOBAL;
            sap_new->segs = (Pointer)ssp;
            ssp_next = ssp->next;
            ssp->next = NULL;
            ssp = ssp_next;
         }
         salp_tmp = salp->next;
         salp = (Pointer)sap_head;
         sap_prev->next = salp_tmp;
         salp = salp_tmp;
      }
   }
   sap->segs = (Pointer)sap_head;
   return retval;
}


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
NLM_EXTERN Boolean AlnMgrIndexLinkedSegs (SeqAlignPtr sap)
{ /* all the Uint2's may have to be changed to Uint4's */
   Int4         currseq;
   DenseSegPtr  dsp;
   Uint2        i;
   Uint4        qlen;
   Boolean      retval;
   SAIndexPtr   saip;
   SASeqDatPtr  ssdp;

   retval = FALSE;
   while (sap)
   {
      if (sap->segtype == SAS_DENSEG)
      {
         dsp = (DenseSegPtr)sap->segs;
         saip = SAIndexNew();
         saip->aligncoords = (Uint4Ptr)MemNew((dsp->numseg+1)*sizeof(Uint4));
         qlen = 0;
         saip->ssdp = (SASeqDatPtr PNTR)MemNew((dsp->dim+1)*sizeof(SASeqDatPtr));
         for (i = 0; i<(dsp->dim); i++)
         {
            ssdp = SASeqDatNew();
            saip->ssdp[i] = ssdp;
         }
         for (i = 0; i<(dsp->numseg); i++)
         {
            saip->aligncoords[i] = qlen;
            qlen += dsp->lens[i];
            for (currseq = 0; currseq<(dsp->dim); currseq++)
            {
               if ((dsp->starts[dsp->dim*i+currseq]) != -1)
               {
                  saip->ssdp[currseq]->numsect++;
               }
            }
         }
         for (currseq = 0; currseq<(dsp->dim); currseq++)
         {
            saip->ssdp[currseq]->sect = (Uint2Ptr)MemNew((saip->ssdp[currseq]->numsect)*sizeof(Uint2));
            saip->ssdp[currseq]->unsect = (Uint2Ptr)MemNew((dsp->numseg - saip->ssdp[currseq]->numsect)*sizeof(Uint2));
            saip->ssdp[currseq]->numsect = 0;
         }
         for (i=0; i<(dsp->numseg); i++)
         {
            for (currseq=0; currseq<(dsp->dim); currseq++)
            {
               if ((dsp->starts[dsp->dim*i+currseq]) != -1)
               {
                  saip->ssdp[currseq]->sect[saip->ssdp[currseq]->numsect] = i;
                  saip->ssdp[currseq]->numsect++;
               } else
               {
                  saip->ssdp[currseq]->unsect[saip->ssdp[currseq]->numunsect]=i;
                  saip->ssdp[currseq]->numunsect++;
               }
            }
         }
         saip->indextype = INDEX_SEGS;
         sap->saip = (SeqAlignIndexPtr)saip;
      }
      sap = sap->next;
      retval = TRUE;
   }
   return retval;
}

NLM_EXTERN Boolean AlnMgrIndexParentSA(SeqAlignPtr sap)
{
   AMAlignDatPtr    amadp;
   AMAlignIndexPtr  amaip;
   Int4             count;
   Boolean          done;
   Int4             i;
   Int4             notfound;
   Int4             numsap;
   Boolean          retval;
   SeqAlignPtr      salp;
   SeqIdPtr         sip;

   retval = FALSE;
   if (!sap)
      return retval;
   if (sap->segtype != SAS_DISC)
      return retval;
   if (((SeqAlignPtr)(sap->segs))->saip == NULL)
   {
      if (!AlnMgrIndexLinkedSegs((SeqAlignPtr)(sap->segs)))
         return retval;
   }
   amaip = (AMAlignIndexPtr)sap->saip;
   if (amaip)
      sap->saip = (Pointer)AMAlignIndexFree(amaip);
   sap->saip = NULL;
   amaip = AMAlignIndexNew();
   count = 0;
   amaip->indextype = INDEX_PARENT;
   amaip->ids = AlnMgrPropagateUpSeqIdPtrs(sap, &count);
   sip = amaip->ids;
   amaip->numbsqs = count;
   amaip->amadp = (AMAlignDatPtr PNTR)MemNew((count+1)*sizeof(AMAlignDatPtr));
   for (count = 0; count < amaip->numbsqs; count++)
   {
      amadp = AMAlignDatNew();
      amaip->amadp[count] = amadp;
      numsap = 0;
      amadp->saps = AlnMgrSortSeqAligns((SeqAlignPtr)sap->segs, AlnMgrCompareIncreasingBySeqIdPtr, sip, &numsap);
      done = FALSE;
      notfound = 0;
      for (i=0; i<numsap && !done; i++)
      {
         if (AlnMgrGetNForSip(amadp->saps[i], sip) < 0)
         {
            notfound++;
         } else
         {
            done = TRUE;
         }
      }
      amadp->numsaps = numsap - notfound;
      for (i=0; i<(numsap - notfound); i++)
      {
         amadp->saps[i] = amadp->saps[i+notfound];
      }
      for (i=(numsap - notfound); i<numsap; i++)
      {
         amadp->saps[i] = NULL;
      }
      amadp->sip = SeqIdDup(sip);
      sip = sip->next;
   }
   i = 0;
   salp = (SeqAlignPtr)sap->segs;
   while (salp)
   {
      i++;
      salp = salp->next;
   }
   amaip->numsaps = i;
   amaip->parent = sap;
   sap->saip = (Pointer)amaip;
   retval = TRUE;
   return retval;
}

NLM_EXTERN SeqIdPtr AlnMgrPropagateUpSeqIdPtrs(SeqAlignPtr sap, Int4Ptr num)
{
   Int4             count;
   DenseSegPtr      dsp;
   Boolean          found;
   SeqAlignPtr      salp;
   SeqIdPtr         sip_head;
   SeqIdPtr         sip_list;
   SeqIdPtr         sip_tmp;
   SeqIdPtr         sip_tmp2;

   if (!sap)
      return NULL;
   if (sap->segtype == SAS_DISC)
      salp = (SeqAlignPtr)(sap->segs);
   else
      salp = sap;
   count = 0;
   sip_list = sip_head = NULL;
   while (salp)
   {
      dsp = (DenseSegPtr)salp->segs;
      sip_tmp = dsp->ids;
      if (!sip_list)
      {
         sip_head = sip_list = SeqIdDup(sip_tmp);
         sip_tmp = sip_tmp->next;
         count++;
      }
      while (sip_tmp)
      {
         sip_tmp2 = sip_head;
         found = FALSE;
         while (sip_tmp2 && !found)
         {
            if (SeqIdComp(sip_tmp, sip_tmp2) == SIC_YES)
               found = TRUE;
            sip_tmp2 = sip_tmp2->next;
         }
         if (!found)
         {
            sip_list->next = SeqIdDup(sip_tmp);
            sip_list = sip_list->next;
            sip_list->next = NULL;
            count++;
         }
         sip_tmp = sip_tmp->next;
      }
      salp = salp->next;
   }
   if (num)
      *num = count;
   return sip_head;
}

NLM_EXTERN SeqIdPtr AlnMgrPropagateSeqIdsBySapList(AMAlignIndexPtr amaip)
{
   DenseSegPtr  dsp;
   Int4         i;
   Int4         j;
   SAIndexPtr   saip;
   SeqAlignPtr  salp;
   SeqIdPtr     sip;
   SeqIdPtr     sip_head;
   SeqIdPtr     sip_tmp;
   SeqIdPtr     sip_tmp2;

   if (amaip == NULL)
      return NULL;
   if (amaip->saps == NULL)
      return NULL;
   sip_head = NULL;
   for (i=0; i<(amaip->alnsaps); i++)
   {
      j=1;
      salp = amaip->saps[i];
      saip = (SAIndexPtr)salp->saip;
      dsp = (DenseSegPtr)(salp->segs);
      sip_tmp = dsp->ids;
      while (j<saip->master)
      {
         sip_tmp = sip_tmp->next;
         j++;
      }
      if (sip_head == NULL)
         sip_head = sip = SeqIdDup(sip_tmp);
      sip_tmp = dsp->ids;
      j=0;
      while (sip_tmp)
      {
         j++;
         if (j!=saip->master)
         {
            sip_tmp2 = SeqIdDup(sip_tmp);
            sip->next = sip_tmp2;
            sip = sip->next;
         }
         sip_tmp = sip_tmp->next; 
      }
   }
   return sip_head;
}

NLM_EXTERN SeqIdPtr  AlnMgrPropagateSeqIdsByRow(AMAlignIndexPtr amaip)
{
   Int4      i;
   SeqIdPtr  sip;
   SeqIdPtr  sip_head;
   SeqIdPtr  sip_tmp;

   if (amaip->rowsource == NULL)
      return NULL;
   sip_head = sip = SeqIdDup(amaip->rowsource[0]->id);
   for (i=1; i<amaip->numrows; i++)
   {
      sip_tmp = SeqIdDup(amaip->rowsource[i]->id);
      sip->next = sip_tmp;
      sip = sip->next;
   }
   return sip_head;
}

NLM_EXTERN Boolean AlnMgrMakeMultipleByScore(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   FloatHi          bit_score;
   Boolean          conflict;
   FloatHi          evalue;
   Int4             i;
   Int4             j;
   Int4             n;
   Int4             number;
   SAIndexPtr       saip1;
   SAIndexPtr       saip2;
   SeqAlignPtr      salp;
   AMAlignInfoPtr   salp_list;
   SeqAlignPtr      PNTR saparray;
   Int4             score;
   SeqIdPtr         sip;
   Int4             start1;
   Int4             start2;
   Int4             startm1;
   Int4             startm2;
   Int4             stop1;
   Int4             stop2;
   Int4             stopm1;
   Int4             stopm2;
   Uint1            strand;
   Uint1            strand_curr;
   AMTinyInfoPtr    PNTR tiparray;

   if (sap == NULL)
      return FALSE;
   i = AlnMgrCheckAlignForParent(sap);
   if (i != AM_PARENT)
      return FALSE;
   amaip = (AMAlignIndexPtr)sap->saip;
   if (amaip == NULL)
      return FALSE;
   if (amaip->numbsqs > 2)
      return FALSE;
   if (sap->master == NULL)
      return FALSE;
   salp = (SeqAlignPtr)sap->segs;
   n = amaip->numsaps;
   salp_list = Calloc(n, sizeof (AMAlignInfo));
   for (salp, i=0; i<n; i++, salp=salp->next)
   {
      salp_list[i].align=salp;
      GetScoreAndEvalue(salp, &score, &bit_score, &evalue, &number);
      salp_list[i].align_len = score;
   }
   HeapSort (salp_list, n, sizeof (AMAlignInfo), AMCompareAlignInfoProc);
   saip1 = (SAIndexPtr)salp_list[0].align->saip;
   if (saip1 == NULL)
      return FALSE;
   strand = AlnMgrGetNthStrand(salp_list[0].align, 3-saip1->master);
   if (strand != Seq_strand_minus)
      strand = Seq_strand_plus;
   amaip->alnsaps = 0;
   for (i=0; i<n; i++)
   {
      if ((saip1 = (SAIndexPtr)salp_list[i].align->saip) == NULL)
         return FALSE;
      AlnMgrGetNthSeqRangeInSA(salp_list[i].align, saip1->master, &startm1, &stopm1);
      AlnMgrGetNthSeqRangeInSA(salp_list[i].align, 3-saip1->master, &start1, &stop1);
      strand_curr = AlnMgrGetNthStrand(salp_list[i].align, 3-saip1->master);
      if (strand_curr != Seq_strand_minus)
         strand_curr = Seq_strand_plus;
      if (strand_curr != strand)
         conflict = TRUE;
      else
         conflict = FALSE;
      for (j=0; j<i && !conflict; j++)
      {
         if ((saip2 = (SAIndexPtr)(amaip->saps[j]->saip)) == NULL)
            return FALSE;
         AlnMgrGetNthSeqRangeInSA(amaip->saps[j], saip2->master, &startm2, &stopm2);
         AlnMgrGetNthSeqRangeInSA(amaip->saps[j], 3-saip2->master, &start2, &stop2);
         if (startm1 < startm2)
         {
            if (stopm1 >= startm2)
               conflict = TRUE;
            if (strand == Seq_strand_minus)
            {
               if (start1 <= stop2)
                  conflict = TRUE;
            } else
            {
               if (stop1 >= start2)
                  conflict = TRUE;
            }
         } else if (startm1 > startm2)
         {
            if (startm1 <= stopm2)
               conflict = TRUE;
            if (strand == Seq_strand_minus)
            {
               if (stop1 >= start2)
                  conflict = TRUE;
            } else
            {
               if (stop2 >= start1)
                  conflict = TRUE;
            }
         } else if (startm1 == startm2)
            conflict = TRUE;
      }
      if (!conflict)
      {
         amaip->saps[amaip->alnsaps] = salp_list[i].align;
         amaip->alnsaps++;
      }
   }
   tiparray = (AMTinyInfoPtr PNTR)MemNew((amaip->alnsaps)*sizeof(AMTinyInfoPtr));
   for (i=0; i<amaip->alnsaps; i++)
   {
      saip1 = (SAIndexPtr)salp_list[i].align->saip;
      AlnMgrGetNthSeqRangeInSA(amaip->saps[i], saip1->master, &start1, &stop1);
      tiparray[i] = (AMTinyInfoPtr)MemNew(sizeof(AMTinyInfo));
      tiparray[i]->start = start1;
      tiparray[i]->stop = stop1;
      tiparray[i]->numgap = saip1->master;
      tiparray[i]->numsap = i;
   }
   HeapSort((Pointer)tiparray, (size_t)(amaip->alnsaps), sizeof(AMTinyInfoPtr), AlnMgrCompareTips);
   saparray = (SeqAlignPtr PNTR)(MemNew((amaip->alnsaps)*sizeof(SeqAlignPtr)));
   for (i=0; i<amaip->alnsaps; i++)
   {
      saparray[i] = amaip->saps[i];
   }
   for (i=0; i<amaip->alnsaps; i++)
   {
      amaip->saps[i] = saparray[tiparray[i]->numsap];
      tiparray[i]->numsap = i;
   }
   MemFree(saparray);
   amaip->numseg = amaip->alnsaps;
   amaip->aligncoords = (Uint4Ptr)MemNew((amaip->alnsaps)*sizeof(Uint4));
   amaip->lens = (Int4Ptr)MemNew((amaip->alnsaps)*sizeof(Int4));
   amaip->rowsource = (RowSourcePtr PNTR)MemNew(2*sizeof(RowSourcePtr));
   amaip->rowsource[0] = (RowSourcePtr)MemNew(sizeof(RowSource));
   amaip->rowsource[0]->id = SeqIdDup(sap->master);
   amaip->rowsource[0]->which_saps = (Uint4Ptr)MemNew((amaip->alnsaps+1)*sizeof(Uint4));
   amaip->rowsource[0]->num_in_sap = (Uint4Ptr)MemNew((amaip->alnsaps+1)*sizeof(Uint4));
   amaip->rowsource[1] = (RowSourcePtr)MemNew(sizeof(RowSource));
   sip = SeqIdDup(AlnMgrGetNthSeqIdPtr(amaip->saps[tiparray[0]->numsap], 3-tiparray[0]->numgap));
   amaip->rowsource[1]->id = sip;
   amaip->rowsource[1]->which_saps = (Uint4Ptr)MemNew((amaip->alnsaps+1)*sizeof(Uint4));
   amaip->rowsource[1]->num_in_sap = (Uint4Ptr)MemNew((amaip->alnsaps+1)*sizeof(Uint4));
   for (i=0; i<amaip->alnsaps; i++)
   {
      amaip->rowsource[0]->which_saps[i] = amaip->rowsource[1]->which_saps[i] = tiparray[i]->numsap + 1;
      amaip->rowsource[0]->num_in_sap[i] = tiparray[i]->numgap;
      amaip->rowsource[1]->num_in_sap[i] = 3-tiparray[i]->numgap;
      amaip->lens[i] = AlnMgrGetAlnLength(amaip->saps[tiparray[i]->numsap], FALSE);;
      if (i>0)
         amaip->aligncoords[i] = amaip->aligncoords[i-1] + amaip->lens[i-1];
      else
         amaip->aligncoords[i] = 0;
   }
   amaip->rowsource[0]->numsaps = amaip->rowsource[1]->numsaps = amaip->alnsaps;
   amaip->master = 1;
   amaip->numrows = 2;
   for (i=0; i<amaip->alnsaps; i++)
   {
      MemFree(tiparray[i]);
   }
   MemFree(tiparray);
   MemFree(salp_list);
   sap->type = SAT_MASTERSLAVE;
   amaip->mstype = AM_SEGMENTED_MASTERSLAVE;
   return TRUE;
}

static int LIBCALLBACK AMCompareAlignInfoProc(VoidPtr ptr1, VoidPtr ptr2)
{
   AMAlignInfoPtr aip_1;
   AMAlignInfoPtr  aip_2;
   if (ptr1 != NULL && ptr2 != NULL)
   {
      aip_1 = (AMAlignInfoPtr) ptr1;
      aip_2 = (AMAlignInfoPtr) ptr2;
      if(aip_1->align_len > aip_2->align_len)
         return -1;
      else if(aip_1->align_len < aip_2->align_len)
         return 1;
      else
         return 0;
  }
  return 0;
}

NLM_EXTERN void am_print_seqalign_indexes(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             i;
   Int4             j;
   SAIndexPtr       saip;

   if (!sap)
      return;
   if (!sap->saip)
      return;
   while (sap)
   {
      if (sap->segtype == SAS_DENSEG && sap->saip)
      {
         dsp = (DenseSegPtr)sap->segs;
         if (sap->saip->indextype == INDEX_SEGS)
            saip = (SAIndexPtr)(sap->saip);
         printf("\naligncoords: ");
         for (i=0; i<(dsp->numseg); i++)
         {
            printf("%d ", saip->aligncoords[i]);
         }
         fflush(stdout);
         for (i=0; i<(dsp->dim); i++)
         {
            printf("\n");
            printf("Sequence %d:", i);
            for (j=0; j<(saip->ssdp[i]->numsect); j++)
            {
               printf("%d ", saip->ssdp[i]->sect[j]);
            }
            fflush(stdout);
         }
      } else if (sap->segtype == SAS_DISC && sap->saip)
      {
         if (sap->saip->indextype == INDEX_PARENT)
            amaip = (AMAlignIndexPtr)(sap->saip);
         if (sap->type == SAT_PARTIAL)
            printf("SAT_PARTIAL\n");
         else if (sap->type == SAT_MASTERSLAVE)
            printf("SAT_MASTERSLAVE\n");
         printf("Parent info:\n");
         printf("numbsqs = %d\n", amaip->numbsqs);
         printf("numsaps = %d\n", amaip->numsaps);
         printf("alnsaps = %d\n", amaip->alnsaps);
         printf("numseg = %d\n", amaip->numseg);
         fflush(stdout);
         for (i=0; i<amaip->numbsqs; i++)
         {
            printf("Sequence %d:", i);
            printf(" %d saps\n", amaip->amadp[i]->numsaps);
            fflush(stdout);
         }
         printf("Starts: ");
         if (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_SEGMENTED_MASTERSLAVE)
         {
            printf("Segmented\n");
         } else 
         {
            for (i=0; i<(amaip->numseg*amaip->numsaps); i++)
            {
               printf("%d ", amaip->starts[i]);
            }
         }
         fflush(stdout);
         printf("\nTotal Length: %d \n", AlnMgrGetAlnLength(sap, TRUE));
         printf("Alignment Length: %d\n", AlnMgrGetAlnLength(sap, FALSE));
         if (amaip->lens)
         {
            printf("lens: ");
            for (i=0; i<amaip->numseg; i++)
            {
               printf("%i ", amaip->lens[i]);
            }
            printf("\n");
            fflush(stdout);
            printf("aligncoords: ");
            for (i=0; i<amaip->numseg; i++)
            {
               printf("%i ", amaip->aligncoords[i]);
            }
            printf("\n");
            fflush(stdout);
         }
         if (amaip->saps)
         {
            for (i=0; i<amaip->numbsqs; i++)
            {
               printf("Segments: ");
               for (j=0; j<(amaip->amadp[i]->numseg); j++)
               {
                  printf("%d ", amaip->amadp[i]->segments[j]);
               }
               printf("\n");
               fflush(stdout);
            }
         }
         if (amaip->rowsource)
         {
            printf("Rowsource arrays:\n");
            for (i=0; i<(amaip->numrows); i++)
            {
               printf("row %d  ", (i+1));
               for (j=0; j<(amaip->rowsource[i]->numsaps); j++)
               {
                  printf("%d: %d   ", amaip->rowsource[i]->which_saps[j], amaip->rowsource[i]->num_in_sap[j]);
               }
               printf("\n");
            }
         }
         am_print_seqalign_indexes((SeqAlignPtr)sap->segs);
      }
      sap = sap->next;
   }
   return;
}

/*CHECK*/
NLM_EXTERN Int4 AlnMgrCheckAlignForParent(SeqAlignPtr sap)
{
   if (sap->segtype == SAS_DISC)
   {
      if (!sap->saip)
      {
         if (!AlnMgrIndexParentSA(sap))
            return -1;
         else
            return AM_PARENT;
      } else if (sap->saip->indextype == INDEX_PARENT)
      {
         return AM_PARENT;
      } else
      {
         return -1;
      }
   } else if (sap->segtype == SAS_DENSEG)
   {
      if (!sap->saip)
      {
         if (sap->segs == NULL)
            return -1;
         AlnMgrAnythingToSeg(sap);
         if (!AlnMgrIndexLinkedSegs((SeqAlignPtr)sap->segs))
            return -1;
         return AM_PARENT;
      } else if (sap->saip->indextype == INDEX_SEGS)
      {
         return AM_CHILD;
      } else
      {
         return -1;
      }
   }
   return -1;
}


/***********************************************************************
*
*  AlnMgrSortSeqAligns is a variant of the ValNodeSort function, and
*  calls very similar heapsort functions.  It can take a comparison
*  function that needs userdata, so more specific sorts are possible
*  without defining special structures for every type of sort.
*
***********************************************************************/
NLM_EXTERN SeqAlignPtr PNTR AlnMgrSortSeqAligns (SeqAlignPtr sap, int (LIBCALLBACK *compar)(VoidPtr, VoidPtr, VoidPtr), VoidPtr userdata, Int4Ptr numsap)
{
   SeqAlignPtr  PNTR head;
   Int4         i;
   Int4         num;
   SeqAlignPtr  tmp;

   if (!sap)
      return NULL;
   tmp = sap;
   num = 0;
   while (tmp)
   {
      num++;
      tmp = tmp->next;
   }
   head = MemNew(((size_t) num + 1)*sizeof(SeqAlignPtr));
   tmp = sap;
   for (i = 0; i<num; i++)
   {
      head[i]=tmp;
      tmp = tmp->next;
      if (!tmp)
         break;
   }
   heapsort_with_userdata(head, (size_t)num, sizeof(SeqAlignPtr), compar, userdata);
   if (numsap)
      *numsap = num;
   return head;
}

static void heapsort_with_userdata (VoidPtr b, size_t nel, size_t width, int (LIBCALLBACK *compar)PROTO((VoidPtr, VoidPtr, VoidPtr)), VoidPtr userdata)
{
   register CharPtr base = (CharPtr)b;
   register size_t i;
   register char ch;
   register CharPtr base0=(CharPtr)base, lim, basef;

   if (nel<2)
      return;
   lim = &base[((nel-2)/2)*width];
   basef = &base[(nel-1)*width];
   i = nel/2;
   for (base = &base0[(i-1)*width]; i>0; base=base-width)
   {
      heapify_with_userdata(base0, base, lim, basef, width, compar, userdata);
      i--;
   }
   for (base=&base0[(nel-1)*width]; base>base0; base -= width)
   {
      for (i = 0; i<width; i++)
      {
         ch = base0[i];
         base0[i] = base[i];
         base[i] = ch;
      }
      lim = base0 + ((base-base0)/2 - width);
      if (base> (base0+width))
         heapify_with_userdata(base0, base0, lim, base-width, width, compar, userdata);
   }
   return;
}

static void heapify_with_userdata(CharPtr base0, CharPtr base, CharPtr lim, CharPtr last, size_t width, int(LIBCALLBACK *compar)PROTO((VoidPtr, VoidPtr, VoidPtr)), VoidPtr userdata)
{
   register size_t i;
   register char ch;
   register CharPtr left_son, large_son;

   left_son = base0 + 2*(base-base0) + width;
   while (base<=lim)
   {
      if (left_son == last)
      {
         large_son = left_son;
      } else
      {
         if((*compar)(left_son, left_son+width, userdata) >= 0)
            large_son = left_son;
         else
            large_son = left_son + width;
      }
      if ((*compar)(base, large_son, userdata) < 0)
      {
         for (i = 0; i<width; i++)
         {
            ch = base[i];
            base[i] = large_son[i];
            large_son[i] = ch;
         }
         base = large_son;
         left_son = base0 + 2*(base-base0) + width;
      } else
      {
         break;
      }
   }
   return;
}

/*************************************************************************
*
*  sorting comparison functions
*
*************************************************************************/
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
NLM_EXTERN int LIBCALLBACK AlnMgrCompareIncreasingBySeqIdPtr (VoidPtr base, VoidPtr large_son, VoidPtr userdata)
{
   Boolean      done;
   DenseSegPtr  dsp1;
   DenseSegPtr  dsp2;
   Int4         n1;
   Int4         n2;
   SeqAlignPtr  sap1;
   SeqAlignPtr  sap2;
   SeqIdPtr     sip;
   SeqIdPtr     sip_tmp;
   Int4         start1;
   Int4         start2;
   Int4         stop1;
   Int4         stop2;
   Uint2        strand1;
   Uint2        strand2;

   sap1 = *((SeqAlignPtr PNTR) base);
   sap2 = *((SeqAlignPtr PNTR) large_son);
   sip = (SeqIdPtr)userdata;
   if (!sap1||!sap2||!sip)
      return 0;
   dsp1 = (DenseSegPtr)sap1->segs;
   dsp2 = (DenseSegPtr)sap2->segs;
   if (!dsp1||!dsp2)
      return 0;
   n1 = n2 = 0;
   done = FALSE;
   sip_tmp = dsp1->ids;
   while (sip_tmp && !done)
   {
      n1 += 1;
      if (SeqIdComp(sip_tmp, sip) == SIC_YES)
         done = TRUE;
      sip_tmp = sip_tmp->next;
   }
   if (!done)
      return -1;
   done = FALSE;
   sip_tmp = dsp2->ids;
   while (sip_tmp && !done)
   {
      n2+=1;
      if (SeqIdComp(sip_tmp, sip) == SIC_YES)
         done = TRUE;
      sip_tmp = sip_tmp->next;
   }
   if (!done)
      return 1;
   AlnMgrGetNthSeqRangeInSA(sap1, n1, &start1, &stop1);
   AlnMgrGetNthSeqRangeInSA(sap2, n2, &start2, &stop2);
   strand1 = AlnMgrGetNthStrand(sap1, n1);
   strand2 = AlnMgrGetNthStrand(sap2, n2);
   if ((strand1 == strand2) && strand1 != Seq_strand_minus)
   {
      if (start1 < start2)
         return -1;
      else if (start2 < start1)
         return 1;
      else if (start1 == start2)
      {
        if (stop1 < stop2)
           return -1;
        else if (stop2 < stop1)
           return 1;
        else
           return 0;
      }
   } else if ((strand1 == strand2) && strand1 ==  Seq_strand_minus)
   {
      if (start1 > start2)
         return -1;
      else if (start2 > start1)
         return 1;
      else if (start1 == start2)
      {
         if (stop1 < stop2)
            return -1;
         else if (stop2 < stop1)
            return 1;
         else
            return 0;
      }
   }
   else
      return 0;
   return 0;
}

/*********************************************************************
*
*  AlnMgrFindFirst is crucial to the AlnMgrMakeFakeMultiple function;
*  it uses the sorted order of the seqaligns in each AMAlignDat
*  structure to guide a heapsort of all the seqaligns.
*
*********************************************************************/
NLM_EXTERN int LIBCALLBACK AlnMgrFindFirst(VoidPtr base, VoidPtr large_son, VoidPtr userdata)
{
   AMAlignDatPtr    amadp;
   AMAlignIndexPtr  amaip;
   Int4             i;
   SeqAlignPtr      sap1;
   SeqAlignPtr      sap2;
   Int4             x;
   Int4             y;
   Int4             z;

   amaip = (AMAlignIndexPtr)userdata;
   if (amaip == NULL || base == NULL || large_son == NULL)
      return 0;
   sap1 = *((SeqAlignPtr PNTR) base);
   sap2 = *((SeqAlignPtr PNTR) large_son);
   if (base == large_son)
      return 0;
   x = y = -1;
   z = amaip->numbsqs;
   while (z)
   {
      amadp = amaip->amadp[(amaip->numbsqs - z)];
      for (i=0; i<(amadp->numsaps); i++)
      {
         if (amadp->saps[i] == sap1)
            x=i;
         else if (amadp->saps[i] == sap2)
            y=i;
      }
      if (x!=-1 && y!=-1)
      {
         if (x < y)
            return -1;
         else if (y < x)
            return 1;
      }
      z--;
   }
   return 0;
}

NLM_EXTERN int LIBCALLBACK AlnMgrCompareTips(VoidPtr base, VoidPtr large_son)
{
   AMTinyInfoPtr  tip1;
   AMTinyInfoPtr  tip2;

   tip1 = *((AMTinyInfoPtr PNTR) base);
   tip2 = *((AMTinyInfoPtr PNTR) large_son);
   if (tip1 == NULL || tip2 == NULL)
      return 0;
   if (tip1->start < tip2->start)
      return -1;
   else if (tip1->start > tip2->start)
      return 1;
   else
      return 0;
}


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
NLM_EXTERN Boolean AlnMgrGetNextLengthBit(SeqAlignPtr sap, Int4Ptr length, Int4Ptr r)
{
   AMAlignIndexPtr  amaip;
   Int4             i;

   if (sap == NULL || length == NULL || r == NULL)
      return FALSE;
   i = AlnMgrCheckAlignForParent(sap);
   if (i == AM_CHILD)
   {
      if (*r == 1)
         return FALSE;
      *length = AlnMgrGetAlnLength(sap, FALSE);
      *r = 1;
      return TRUE;
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      if (sap->type == SAT_PARTIAL)
      {
         if (*r < 0)
         {
            if ((-*r) >= amaip->numsaps)
               return FALSE;
            if (amaip->ulens == NULL)
               *length = -(AlnMgrGetMaxUnalignedLength(amaip->saps[(-*r)-1], amaip->saps[(-*r)]));
            else
               *length = -(amaip->ulens[(-*r)-1]);
            *r = -(*r);
            return TRUE;
         } else
         {
            if (*r >= amaip->numsaps)
               return FALSE;
            *length = AlnMgrGetAlnLength(amaip->saps[*r], FALSE);
            *r = -((*r)+1);
            return TRUE;
         }
      } else if (sap->type == SAT_MASTERSLAVE)
      {
         if (amaip->mstype == AM_MASTERSLAVE)
         {
            if (*r == 1)
               return FALSE;
            *length = amaip->aligncoords[amaip->numseg-1] + amaip->lens[amaip->numseg-1];
            *r = 1;
            return TRUE;
         } else if (amaip->mstype == AM_SEGMENTED_MASTERSLAVE)
         {
            if (*r < 0)
            {
               if ((-*r) >= amaip->numsaps)
                  return FALSE;
               if (amaip->ulens == NULL)
                  *length = -(AlnMgrGetMaxUnalignedLength(amaip->saps[(-*r)-1], amaip->saps[(-*r)]));
               else
                  *length = -(amaip->ulens[(-*r)-1]);
               *r = -(*r);
               return TRUE;
            } else
            {
               if (*r >= amaip->numsaps)
                  return FALSE;
               *length = AlnMgrGetAlnLength(amaip->saps[*r], FALSE);
               *r = -((*r)+1);
               return TRUE;
            }
         }
      }
   }
   return FALSE;
}

NLM_EXTERN Int4 AlnMgrGetMaxUnalignedLength(SeqAlignPtr sap1, SeqAlignPtr sap2)
{
   DenseSegPtr  dsp;
   Int4         max;
   Int4         n1;
   Int4         n2;
   SeqIdPtr     sip;
   Int4         start1;
   Int4         start2;
   Int4         stop1;
   Int4         stop2;

   if (sap1 == NULL || sap2 == NULL)
      return 0;
   dsp = (DenseSegPtr)sap1->segs;
   sip = dsp->ids;
   max = 0;
   while (sip)
   {
      n1 = AlnMgrGetNForSip(sap1, sip);
      AlnMgrGetNthSeqRangeInSA(sap1, n1, &start1, &stop1);
      n2 = AlnMgrGetNForSip(sap2, sip);
      if (n2 >= 0)
      {
         AlnMgrGetNthSeqRangeInSA(sap2, n2, &start2, &stop2);
         if (start2 > stop1)
         {
            if (start2 - stop1 - 1 > max)
               max = start2 - stop1 - 1;
         } else
         {
            if (start1 - stop2 - 1 > max)
               max = start1 - stop2 - 1;
         }
      }
      sip = sip->next;
   }
   return max;
}


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
*  retrieve the required data from the indexes.  (NEXT)
*
*************************************************************************/
NLM_EXTERN Boolean AlnMgrGetNextAlnBit (SeqAlignPtr sap, AlnMsgPtr amp)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             endoffset;
   Boolean          found;
   Int4             i;
   Int4             len;
   Boolean          more;
   Int4             offset;
   Boolean          retval;
   Int4             rf_tmp;
   SAIndexPtr       saip;
   SASeqDatPtr      ssdp;
   Int4             start_b;
   Uint4            start_m;
   Uint4            start_tmp;
   Uint4            stop_m;
   Uint4            stop_tmp;

   retval = FALSE;
   if (!sap)
      return retval;
   if (!amp)
      return retval;
   i = AlnMgrCheckAlignForParent(sap);
   if (i == AM_CHILD)
   {
      saip = (SAIndexPtr)sap->saip;
      dsp = (DenseSegPtr)sap->segs;
      if (!dsp)
         return retval;
      if (!amp->which_master)
      {
         if (amp->place == 1)
            return retval;
         len = AlnMgrGetAlnLength(sap, FALSE);
         if (amp->to_m > len-1)
            return retval;
         if (amp->to_m < 0)
            amp->to_m = len-1;
         if (amp->row_num == -1)
         {
            if (!amp->which_bsq)
               return retval;
            amp->row_num = AlnMgrGetNForSip(sap, amp->which_bsq) - 1;
            if (amp->row_num == -1)
               return retval;
         }
         if (amp->prev != -2)
         {
            start_m = amp->prev;
         } else
         {
            start_m = binary_search_on_uint4_list(saip->aligncoords, amp->from_m, dsp->numseg);
            amp->real_from = amp->from_m;
         }
         stop_m = binary_search_on_uint4_list(saip->aligncoords, amp->to_m, dsp->numseg);
         ssdp = saip->ssdp[amp->row_num-1];
         offset = amp->real_from - saip->aligncoords[start_m];
         start_b = binary_search_on_uint2_list(ssdp->sect, start_m, ssdp->numsect);
         if (dsp->strands[start_b*(dsp->dim)+amp->row_num-1] == Seq_strand_minus)
            amp->strand = Seq_strand_minus;
         else
            amp->strand = Seq_strand_plus;
         if ((stop_m - start_m) > 0)
         {
            retval = TRUE;
            if (start_b >= 0)
            {
               if (amp->strand != Seq_strand_minus)
               {
                  amp->from_b = dsp->starts[start_b*(dsp->dim)+amp->row_num-1] + offset;
                  amp->to_b = dsp->starts[start_b*(dsp->dim)+amp->row_num-1] + dsp->lens[start_b] - 1;
               } else
               {
                  amp->from_b = dsp->starts[start_b*(dsp->dim)+amp->row_num-1];
                  amp->to_b = amp->from_b + dsp->lens[start_b] - 1 - offset;
               }
               amp->gap = 0;
            } else
            {
               amp->from_b = amp->real_from;
               amp->to_b = saip->aligncoords[start_m + 1] - 1;
               amp->gap = 1; 
            }
            amp->real_from = saip->aligncoords[start_m + 1];
            amp->prev = start_m + 1;
         } else
         {
            amp->place = 1;
            endoffset = amp->to_m - saip->aligncoords[start_m];
            if (start_b >= 0)
            {
               if (amp->strand != Seq_strand_minus)
               {
                  amp->from_b = dsp->starts[start_b*(dsp->dim)+amp->row_num-1] + offset;
                  amp->to_b = dsp->starts[start_b*(dsp->dim)+amp->row_num-1] + endoffset;
               } else
               {
                  amp->from_b = dsp->starts[start_b*(dsp->dim)+amp->row_num-1] + dsp->lens[start_b]  - endoffset - 1;
                  amp->to_b = amp->from_b + amp->to_m - amp->real_from;
               }
               amp->gap = 0;
            } else
            {
               amp->from_b = amp->real_from;
               amp->to_b = amp->to_m;
               amp->gap = 1;
            }
            amp->real_from = 0;
            amp->prev = -2;
            retval = TRUE;
         }
      }
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (!amaip->saps)
         return retval;
      if (amp->place == 1)
         return retval;
      if (!amp->which_bsq && amp->row_num==-1)
         return retval;
      if (sap->type == SAT_PARTIAL && amp->which_master == NULL)
      {
         len = AlnMgrGetAlnLength(sap, FALSE);
         if (amp->to_m < 0)
            amp->to_m = len-1;
         if (amp->to_m > len-1)
            return FALSE;
         if (amp->prev_sap != -2)
         {
            start_m = amp->prev_sap;
            amp->len_left = amp->to_m - amp->real_from + 1;
         } else
         {
            start_m = binary_search_on_uint4_list(amaip->aligncoords, amp->from_m, amaip->alnsaps);
            amp->real_from = amp->from_m;
            amp->prev_sap = start_m;
            amp->len_left = amp->to_m - amp->from_m + 1;
         }
         stop_m = binary_search_on_uint4_list(amaip->aligncoords, amp->to_m, amaip->alnsaps);
         offset = amp->real_from - amaip->aligncoords[start_m];
         if (amp->len_left > (amaip->lens[start_m]-offset))
         {
            endoffset = amaip->lens[start_m] - offset;
         } else
         {
            endoffset = amp->len_left;
         }
         stop_tmp = amp->to_m;
         start_tmp = amp->from_m;
         if ((stop_m - start_m) == 0)
         {
            amp->from_m = offset + amaip->starts[start_m];
            amp->to_m = amp->from_m + endoffset - 1;
            amp->prev = -2;
            rf_tmp = amp->real_from;
            AlnMgrGetNextAlnBit((amaip->saps[start_m]), amp);
            amp->len_left = amp->len_left - (amp->to_b - amp->from_b + 1);
            amp->to_m = stop_tmp;
            amp->from_m = start_tmp;
            if (amp->len_left == 0)
            {
               amp->real_from = amp->to_m + 1;
               amp->prev_sap = -2;
               amp->place = 1;
            } else
            {
               amp->real_from = rf_tmp + (amp->to_b - amp->from_b + 1);
               amp->place = 0;
            }
            retval = TRUE;
         } else
         {
            retval = TRUE;
            amp->from_m = offset + amaip->starts[start_m];
            amp->to_m = amp->from_m + endoffset - 1;
            more = AlnMgrGetNextAlnBit((amaip->saps[start_m]), amp);
            amp->len_left = amp->len_left - (amp->to_m - amp->from_m + 1);
            amp->to_m = stop_tmp;
            amp->real_from = amp->to_m - amp->len_left + 1;
            amp->from_m = start_tmp;
            if (amp->place == 1)
            {
               amp->prev_sap += 1;
               amp->send_space = TRUE;
               if (amp->len_left > 0)
                  amp->place = 0;
            }
            if (amp->len_left == 0 || amp->real_from >= amp->to_m)
            {
               amp->place = 1;
               retval = FALSE;
               amp->prev_sap = -2;
            }
         }
      } else if (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_MASTERSLAVE && amp->which_master == NULL)
      {
         if (amp->place == 1)
            return retval;
         len = AlnMgrGetAlnLength(sap, FALSE);
         if (amp->to_m > len-1)
            return retval;
         if (amp->to_m < 0)
            amp->to_m = len-1;
         if (amp->row_num == -1) 
         { 
            if(!amp->which_bsq)
                return retval;
            else 
            {
                amp->row_num = AlnMgrGetNForSip(sap,amp->which_bsq);
                if(amp->row_num == -1)
                    return retval;
            }
         }
         if (amp->row_num == amaip->master)
         {
            amp->strand = Seq_strand_plus;
            if (amp->prev != -2)
            {
               amp->prev += 1;
               start_m = amp->prev;
            } else
            {
               start_m = binary_search_on_uint4_list(amaip->aligncoords, amp->from_m, amaip->numseg);
               amp->real_from = amp->from_m;
               amp->prev = start_m;
            }
            offset = amp->real_from - amaip->aligncoords[start_m];
            endoffset = amaip->lens[start_m] - offset - (amp->to_m - amp->real_from + 1);
            if (endoffset < 0 && (start_m+1) < amaip->numseg)
               retval = TRUE;
            else
            {
               retval = TRUE;
               amp->place = 1;
               amp->row_num = -1;
               amp->prev = -2;
            }
            i=0;
            found = FALSE;
            while (!found && i < amaip->numsaps)
            {
               if (amaip->starts[i+(amaip->numsaps)*start_m] >= 0)
                  found = TRUE;
               else
                  i++;
            }
            amp->from_b = AlnMgrMapToBsqCoords(amaip->saps[i], amaip->starts[i+(amaip->numsaps)*start_m]+offset, NULL, NULL);
            if (amp->from_b >= 0)
            {
               if (endoffset >=0)
                  amp->to_b = amp->from_b + amaip->lens[start_m] - 1 - offset  - endoffset;
               else
                  amp->to_b = amp->from_b + amaip->lens[start_m] - offset - 1;
               amp->gap = 0;
            } else
            { 
               amp->from_b = amp->real_from;
               amp->gap = 1;
               if (endoffset >= 0)
                  amp->to_b = amp->from_b + amaip->lens[start_m] - 1 - offset  - endoffset;
               else
                  amp->to_b = amp->from_b + amaip->lens[start_m] - offset - 1;
            }
            amp->real_from += amp->to_b - amp->from_b + 1;
         } else
         {
            if (amp->prev != -2)
            {
               amp->prev += 1;
               start_m = amp->prev;
            } else
            {
               start_m = binary_search_on_uint4_list(amaip->aligncoords, amp->from_m, amaip->numseg);
               amp->real_from = amp->from_m;
               amp->prev = start_m;
            }
            if (amp->prev_sap == -2)
               amp->prev_sap=amaip->rowsource[amp->row_num-1]->which_saps[0];
            i = amp->prev_sap-1;
            amp->strand = AlnMgrGetNthStrand(amaip->saps[i], amaip->rowsource[amp->row_num-1]->num_in_sap[0]);
            offset = amp->real_from - amaip->aligncoords[start_m];
            endoffset = amaip->lens[start_m] - offset - (amp->to_m - amp->real_from + 1);
            if (endoffset <= 0 && (start_m + 1) < amaip->numseg)
               retval = TRUE;
            else
            {
               retval = TRUE;
               amp->place = 1;
               amp->prev = amp->prev_sap = -2;
            }
            if (amaip->starts[i+(amaip->numsaps)*start_m] < 0)
               amp->from_b = -1;
            else
               amp->from_b = AlnMgrMapRowCoords(amaip->saps[i], amaip->starts[i+(amaip->numsaps)*start_m]+offset, amaip->rowsource[amp->row_num-1]->num_in_sap[0], NULL);
            if (amp->from_b >= 0)
            {
               if (amp->strand != Seq_strand_minus)
               {
                  if (endoffset >=0)
                     amp->to_b = amp->from_b + amaip->lens[start_m] - offset - 1 -
endoffset;
                  else
                     amp->to_b = amp->from_b + amaip->lens[start_m] - offset -1;
               } else
               {
                  amp->to_b = amp->from_b;
                  if (endoffset >= 0)
                     amp->from_b = amp->to_b + amaip->lens[start_m] - 1 + endoffset
;
                  else
                     amp->from_b = amp->to_b + amaip->lens[start_m] - 1;
               }
               amp->gap = 0;
            } else
            {
               amp->from_b = amp->real_from;
               amp->gap = 1;
               if (endoffset >= 0)
                  amp->to_b = amp->from_b + amaip->lens[start_m] - offset - 1 - endoffset;
               else
                  amp->to_b = amp->from_b + amaip->lens[start_m] - offset - 1;
            }
            amp->real_from += amp->to_b - amp->from_b + 1;
            if (amp->real_from > amp->to_m)
            {
               retval = TRUE;
               amp->place = 1;
               amp->row_num = -1;
               amp->prev = -2;
            }
         }
      } else if (sap->type == SAT_MASTERSLAVE && amp->which_master)
      {
      } else if (sap->type == SAT_DIAGS && amp->which_master)
      {
      } else if (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_SEGMENTED_MASTERSLAVE && amp->which_master == NULL)
      {
         if (amp->place == 1)
            return retval;
         len = AlnMgrGetAlnLength(sap, FALSE);
         if (amp->to_m > len-1)
            return retval;
         if (amp->to_m < 0)
            amp->to_m = len-1;
         if (amp->row_num == -1)
         {
            if(!amp->which_bsq)
                return retval;
            else
            {
                amp->row_num = AlnMgrGetNForSip(sap,amp->which_bsq);
                if(amp->row_num == -1)
                    return retval;
            }
         }
         if (amp->prev == -2)
         {
            start_m = binary_search_on_uint4_list(amaip->aligncoords, amp->from_m, amaip->numseg);
            amp->real_from = amp->from_m;
            amp->prev = start_m;
         } else
            start_m = amp->prev;
         offset = amp->real_from - amaip->aligncoords[start_m];
         if (offset < 0)
            offset = 0;
         if (amaip->rowsource[amp->row_num-1]->which_saps[start_m] == 0)
         {
            amp->from_b = amaip->aligncoords[start_m]+offset;
            amp->gap = 2;
            amp->strand = Seq_strand_unknown;
         } else
         {
            len = 0;
            amp->strand = AlnMgrGetNthStrand(amaip->saps[amaip->rowsource[amp->row_num-1]->which_saps[start_m]-1], amaip->rowsource[amp->row_num-1]->num_in_sap[start_m]);
            amp->from_b = AlnMgrMapSegmentCoords(amaip->saps[amaip->rowsource[amp->row_num-1]->which_saps[start_m]-1], offset, amaip->rowsource[amp->row_num-1]->num_in_sap[start_m], NULL, &len);
            if (amp->from_b == -1)
            {
               amp->from_b = amaip->aligncoords[start_m]+offset;
               amp->gap = 1;
            } else
               amp->gap = 0;
         }
         endoffset = amp->to_m - (amaip->aligncoords[start_m] + len + offset -1);
         if (endoffset <= 0)
            amp->place = 1;
         else if (len >= amaip->lens[start_m] - offset)
         {
            amp->prev++;
            amp->send_space = 1;
         } else
            amp->send_space = 0;
         if (endoffset <= 0)
            amp->to_b = amp->from_b + len + endoffset -1;
         else
         {
            amp->to_b = amp->from_b + len - 1;
            amp->real_from = amp->real_from + amp->to_b - amp->from_b + 1;
         }
         if (amp->strand == Seq_strand_minus && amp->gap == 0)
         {
            offset = amp->to_b - amp->from_b;
            amp->to_b = amp->from_b;
            amp->from_b = amp->to_b - offset;
         }
         retval = TRUE;
      }
   }
   return retval;
}

NLM_EXTERN Uint4 binary_search_on_uint4_list(Uint4Ptr list, Uint4 pos, Uint4 listlen)
{
   Uint4  L;
   Uint4  mid;
   Uint4  R;

   if (list == NULL || listlen == 0)
      return -1;
   L = 0;
   R = listlen - 1;
   while (L < R)
   {
      mid = (L+R)/2;
      if (list[mid + 1] <= pos)
      {
         L = mid + 1;
      } else
      {
         R = mid;
      }
   }
   return R;
}

NLM_EXTERN Int4 binary_search_on_uint2_list(Uint2Ptr list, Uint2 ele, Uint2 listlen)
{
   Uint2  L;
   Uint2  mid;
   Uint2  R;

   if (list == NULL || listlen == 0)
      return -1;
   L = 0;
   R = listlen - 1;
   while (L < R)
   {
      mid = (L+R)/2;
      if (ele <= list[mid])
      {
         R = mid;
      } else
      {
         L = mid+1;
      }
   }
   if (ele == list[R])
      return list[R];
   else
      return -1;
}

NLM_EXTERN Int4 binary_search_by_chunk(Int4Ptr list, Int4 ele, Int4 listlen, Int4 chunksize, Int4 offset)
{
   Int4  L;
   Int4  mid;
   Int4  R;

   if (list == NULL || listlen == 0)
      return -1;
   L = 0;
   R = (listlen/chunksize) - 1;
   while (L < R)
   {
      mid = (L+R)/2;
      if (ele <= list[(mid)*chunksize + offset] && list[(mid)*chunksize + offset] >= 0)
      {
         R = mid;
      } else
      {
         L = mid + 1;
      }
   }
   return R;
}

NLM_EXTERN Int4 binary_search_segment_array(SASeqDatPtr ssdp, Int4 pos, Int4 numseq, Int4 offset, DenseSegPtr dsp)
{
   Int4  L;
   Int4  mid;
   Int4  R;

   if (ssdp == NULL || numseq == 0)
      return -1;
   L = 0;
   R = ssdp->numsect - 1;
   while (L < R)
   {
      mid = (L+R)/2;
      if (pos <= (dsp->starts[(ssdp->sect[mid])*numseq + offset]))
      {
         R = mid;
      } else
      {
         L = mid+1;
      }
   }
   return (ssdp->sect[R]);
}

/************************************************************************
*
*  These are several utility functions which get needed data from the
*  indexes.  "N" refers to row number.
*
************************************************************************/
NLM_EXTERN Int4 AlnMgrGetAlnLength(SeqAlignPtr sap, Boolean fill_in)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             i;
   Int4             length;
   SAIndexPtr       saip;

   if (!sap)
      return 0;
   i = AlnMgrCheckAlignForParent(sap);
   if (i<0)
   {
      return 0;
   } else if (i == AM_CHILD)
   {
      dsp = (DenseSegPtr)sap->segs;
      if (!dsp)
         return 0;
      saip = (SAIndexPtr)sap->saip;
      return ((saip->aligncoords[dsp->numseg-1])+dsp->lens[dsp->numseg-1]);
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      if (!amaip)
         return 0;
      if (!amaip->saps)
      {
         if (!AlnMgrMakeFakeMultiple(sap))
            return 0;
      }
      if (fill_in)
      {
         if (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_MASTERSLAVE)
            return (amaip->lens[(amaip->numseg)-1] + amaip->aligncoords[amaip->numseg-1]);
         else if (sap->type == SAT_PARTIAL)
         {
            length = 0;
            for (i=0; i<(amaip->numsaps-1); i++)
            {
               length += AlnMgrGetMaxUnalignedLength(amaip->saps[i], amaip->saps[i+1]);
            }
            return (length + amaip->lens[(amaip->numseg)-1] + amaip->aligncoords[amaip->numseg-1]);
         }
      } else
      {
         return (amaip->lens[(amaip->numseg)-1] + amaip->aligncoords[amaip->numseg-1]);
      }
   }
   return 0;
}

NLM_EXTERN Int4 AlnMgrGetNumSeqs(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             i;

   if (!sap)
      return 0;
   i = AlnMgrCheckAlignForParent(sap);
   if (i<0)
      return 0;
   if (i == AM_CHILD)
   {
      dsp = (DenseSegPtr)sap->segs;
      if (!dsp)
         return 0;
      return (dsp->dim);
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (!amaip)
         return 0;
      return (amaip->numbsqs);
   }
   return 0;
}

NLM_EXTERN SeqIdPtr AlnMgrGetUniqueSeqs(SeqAlignPtr sap, Int4Ptr n)
{
   AMAlignIndexPtr  amaip;
   Int4             c;
   DenseSegPtr      dsp;
   Boolean          found;
   Int4             i;
   Int4             m;
   SeqIdPtr         sip;
   SeqIdPtr         sip_head;
   SeqIdPtr         sip_prev;
   SeqIdPtr         sip_tmp;

   if (sap == NULL)
      return 0;
   i = AlnMgrCheckAlignForParent(sap);
   if (i<0)
      return 0;
   sip_head = sip_prev = NULL;
   if (i == AM_CHILD)
   {
      dsp = (DenseSegPtr)sap->segs;
      if (dsp == NULL)
         return 0;
      sip = dsp->ids;
      m = 0;
      while (sip)
      {
         sip_tmp = sip_head;
         found = FALSE;
         while (!found && sip_tmp != NULL)
         {
            if (SAM_OrderSeqID(sip, sip_tmp) == 0)
               found = TRUE;
            sip_tmp = sip_tmp->next;
         }
         if (!found)
         {
            m++;
            if (sip_head == NULL)
            {
               sip_head = sip_prev = SeqIdDup(sip);
            } else
            {
               sip_prev->next = SeqIdDup(sip);
               sip_prev = sip_prev->next;
            }
         }
         sip = sip->next;
      }
      if (n)
         *n = m;
      return sip_head;
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (amaip == NULL)
         return 0;
      m = 0;
      if (amaip->alnsaps == 1)
      {
         return (AlnMgrGetUniqueSeqs((SeqAlignPtr)sap->segs, n));
      }
      for (c=0; c<amaip->numrows; c++)
      {
         sip = amaip->rowsource[c]->id;
         sip_tmp = sip_head;
         found = FALSE;
         while (!found && sip_tmp != NULL)
         {
            if (SAM_OrderSeqID(sip, sip_tmp) == 0)
               found = TRUE;
            sip_tmp = sip_tmp->next;
         }
         if (!found)
         {
            m++;
            if (sip_head == NULL)
            {
               sip_head = sip_prev = SeqIdDup(sip);
            } else
            {
               sip_prev->next = SeqIdDup(sip);
               sip_prev = sip_prev->next;
            }
         }
      }
      if (n)
         *n = m;
      return sip_head;
   }
   return NULL;
}

NLM_EXTERN SeqIdPtr AlnMgrGetNthSeqIdPtr(SeqAlignPtr sap, Int4 n)
{
   AMAlignIndexPtr  amaip;
   Int4             count;
   DenseSegPtr      dsp;
   Int4             i;
   SeqIdPtr         sip;

   if (!sap)
      return NULL;
   i = AlnMgrCheckAlignForParent(sap);
   if (i<0)
      return NULL;
   else if (i == AM_CHILD)
   {
      dsp = (DenseSegPtr)sap->segs;
      if (!dsp)
         return NULL;
      sip = dsp->ids;
      count = 0;
      while (sip)
      {
         count++;
         if (count == n)
            return (SeqIdDup(sip));
         sip = sip->next;
      }
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (n <= amaip->numrows)
      {
         return (SeqIdDup(amaip->rowsource[n-1]->id));
      } else
         return NULL;
   }
   return NULL;
}

/* (RANGE) */
NLM_EXTERN void AlnMgrGetNthSeqRangeInSA(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop)
{
   AMAlignDatPtr    amadp;
   AMAlignIndexPtr  amaip;
   Uint2            beg;
   Int4             bsq;
   DenseSegPtr      dsp;
   Uint2            end;
   Int4             i;
   Int4             j;
   RowSourcePtr     rsp;
   SAIndexPtr       saip;
   SeqIdPtr         sip;
   Uint2            strand;
   Int4             tmp_beg;
   Int4             tmp_end;
   Int4             tmp_start;
   Int4             tmp_stop;

   if (!sap)
      return;
   i = AlnMgrCheckAlignForParent(sap);
   if (i < 0)
   {
      return;
   }
   else if (i == AM_CHILD)
   {
      if (n<1)
         return;
      saip = (SAIndexPtr)(sap->saip);
      bsq = n-1;
      dsp = (DenseSegPtr)sap->segs;
      if (n > dsp->dim)
         return;
      if (!dsp)
         return;
      strand = dsp->strands[bsq];
      if (strand != Seq_strand_minus)
      {
         beg = saip->ssdp[bsq]->sect[0];
         if (start)
            *start = dsp->starts[beg*(dsp->dim)+bsq];
         end = saip->ssdp[bsq]->sect[(saip->ssdp[bsq]->numsect)-1];
         if (stop)
            *stop = (dsp->starts[end*(dsp->dim)+bsq] + dsp->lens[end] - 1);
         return;
      } else
      {
         beg = saip->ssdp[bsq]->sect[(saip->ssdp[bsq]->numsect)-1];
         if (start)
            *start = dsp->starts[beg*(dsp->dim)+bsq];
         end = saip->ssdp[bsq]->sect[0];
         if (stop)
            *stop = (dsp->starts[end*(dsp->dim)+bsq] + dsp->lens[end] - 1);
         return;
      }
   } else if (i == AM_PARENT)
   {
      if (n<1)
         return;
      bsq = n-1;
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (amaip->rowsource == NULL)
      {
         amadp = amaip->amadp[bsq];
         sip = amaip->ids;
         for (j = 0; j<bsq; j++)
         {
            sip = sip->next;
            if (sip == NULL)
               return;
         }
         for (j = 0; j<(amadp->numsaps); j++)
         {
            tmp_start = tmp_stop = 0;
            AlnMgrGetNthSeqRangeInSA(amadp->saps[j], AlnMgrGetNForSip(amadp->saps[j], sip), &tmp_start, &tmp_stop);
            if (j == 0)
            {
               tmp_beg = tmp_start;
               tmp_end = tmp_stop;
            } else
            {
               if (tmp_start < tmp_beg)
                  tmp_beg = tmp_start;
               if (tmp_stop > tmp_end)
                  tmp_end = tmp_stop;
            }
         }
         if (start)
            *start = tmp_beg;
         if (stop)
            *stop = tmp_end;
         return;
      } else
      {
         sip = amaip->ids;
         if (n > amaip->numrows)
            return;
         rsp = (RowSourcePtr)amaip->rowsource[n-1];
         for (j=0; j<(rsp->numsaps); j++)
         {
            tmp_start = tmp_stop = 0;
            AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[j]-1], rsp->num_in_sap[j], &tmp_start, &tmp_stop);
            if (j==0)
            {
               tmp_beg = tmp_start;
               tmp_end = tmp_stop;
            } else
            {
               if (tmp_start < tmp_beg)
                  tmp_beg = tmp_start;
               if (tmp_stop > tmp_end)
                  tmp_end = tmp_stop;
            }
         }
         if (start)
            *start = tmp_beg;
         if (stop)
            *stop = tmp_end;
         return;
      }
   }
   return;
}

NLM_EXTERN Int4 AlnMgrGetNumSegments(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             i;

   if (sap == NULL)
      return -1;
   i = AlnMgrCheckAlignForParent(sap);
   if (i == AM_CHILD)
   {
      dsp = (DenseSegPtr)(sap->segs);
      return (dsp->numseg);
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      return (amaip->numseg);
   } else
      return -1;
}

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
NLM_EXTERN Boolean AlnMgrGetNextNthSeqRange(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop, Int4Ptr where, BoolPtr is_aligned, Boolean unaligned)
{
   if (sap == NULL || n <= 0)
      return FALSE;
   if (sap->saip == NULL)
      return FALSE;
   if (sap->saip->indextype == INDEX_PARENT && sap->type == SAT_PARTIAL)
   {
      return (am_get_nth_range_for_partial(sap, n, start, stop, where, is_aligned, unaligned));
   } else
   {
      if (*where == 0)
      {
         AlnMgrGetNthSeqRangeInSA(sap, n, start, stop);
         *where = 1;
         return TRUE;
      } else
         return FALSE;
   }
}

static Boolean am_get_nth_range_for_partial(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop, Int4Ptr where, BoolPtr is_aligned, Boolean unaligned)
{
   AMAlignIndexPtr  amaip;
   RowSourcePtr     rsp;
   Uint2            strand;
   Int4             tmp_start;
   Int4             tmp_stop;
   Int4             tmp_where;

   amaip = (AMAlignIndexPtr)sap->saip;
   rsp = amaip->rowsource[n-1];
   tmp_where = *where;
   if (tmp_where >= 0)
   {
      if (tmp_where >= rsp->numsaps)
         return FALSE;
      if (is_aligned)
         *is_aligned = TRUE;
      AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[tmp_where]-1], rsp->num_in_sap[tmp_where], start, stop);
      if (unaligned && (sap->type == SAT_PARTIAL || (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_SEGMENTED_MASTERSLAVE)))
         tmp_where = -(tmp_where+1);
      else
         tmp_where += 1;
   } else if (tmp_where < 0 && unaligned == TRUE)
   {
      if (-tmp_where >= rsp->numsaps)
         return FALSE;
      if (is_aligned)
         *is_aligned = FALSE;
      strand = AlnMgrGetNthStrand(amaip->saps[rsp->which_saps[(-tmp_where)]-1], n);
      tmp_start = tmp_stop = 0;
      if (start)
      { 
         if (strand == Seq_strand_minus)
         {
            AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[(-tmp_where)]-1], rsp->num_in_sap[(-tmp_where)], &tmp_start, NULL);
            *start = tmp_start + 1;
         } else
         {
             AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[(-tmp_where)-1]-1], rsp->num_in_sap[(-tmp_where)-1], NULL, &tmp_start);
             *start = tmp_start + 1;
         }
      }
      if (stop)
      {
         if (strand == Seq_strand_minus)
         {
            AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[(-tmp_where)-1]-1], rsp->num_in_sap[(-tmp_where)-1], NULL, &tmp_stop);
            *stop = tmp_stop - 1;
         } else
         {
            AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[-tmp_where]-1], rsp->num_in_sap[-tmp_where], &tmp_stop, NULL);
            *stop = tmp_stop - 1;
         }
      }
      if (tmp_start + 1 > tmp_stop - 1)
      {
         if (start)
            *start = -1;
         if (stop)
            *stop = -1;
      }
      tmp_where = -tmp_where;
   }
   *where = tmp_where;
   return TRUE;
}

/********************************************************************************
*
*  AlnMgrGetNthRowTail retrieves the blocks of sequence on either end of the
*  alignment, by row.  which_tail is LEFT_TAIL to retrieve the ends which come
*  before alignment coordinate 0, and RIGHT_TAIL to retrieve the other ends.
*  The function returns TRUE if successful, FALSE for an error.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrGetNthRowTail(SeqAlignPtr sap, Int4 n, Uint1 which_tail, Int4Ptr start, Int4Ptr stop, Uint1Ptr strand)
{
   BioseqPtr  bsp;
   SeqIdPtr   sip;
   Int4       tmp_start;
   Int4       tmp_stop;
   Uint1      tmp_strand;

   if (sap == NULL || n < 1)
      return FALSE;
   tmp_start = tmp_stop = -1;
   AlnMgrGetNthSeqRangeInSA(sap, n, &tmp_start, &tmp_stop);
   if (tmp_start == -1 || tmp_stop == -1)
      return FALSE;
   tmp_strand = AlnMgrGetNthStrand(sap, n);
   if (which_tail == LEFT_TAIL)
   {
      if (tmp_strand == Seq_strand_minus)
      {
         sip = AlnMgrGetNthSeqIdPtr(sap, n);
         bsp = BioseqLockById(sip);
         if (tmp_stop == bsp->length-1 || stop == NULL)
         {
            if (start)
               *start = -1;
            if (stop)
               *stop = -1;
         } else
         {
            if (bsp == NULL)
               return FALSE;
            if (start)
               *start = tmp_stop-1;
            if (stop)
               *stop = bsp->length-1;
         }
         BioseqUnlockById(sip);
         if (strand)
            *strand = tmp_strand;
      } else
      {
         if (tmp_start >= 1)
         {
            if (start)
               *start = 0;
            if (stop)
               *stop = tmp_start - 1;
         } else
         {
            if (start)
               *start = -1;
            if (stop)
               *stop = -1;
         }
         if (strand)
            *strand = tmp_strand;
      }
   } else if (which_tail == RIGHT_TAIL)
   {
      if (tmp_strand == Seq_strand_minus)
      {
         if (tmp_start >= 1)
         {
            if (start)
               *start = 0;
            if (stop)
               *stop = tmp_start - 1;
         } else
         {
            if (start)
               *start = -1;
            if (stop)
               *stop = -1;
         }
         if (strand)
            *strand = tmp_strand;
      } else
      {
         sip = AlnMgrGetNthSeqIdPtr(sap, n);
         bsp = BioseqLockById(sip);
         if (bsp == NULL)
            return FALSE;
         if (bsp->length-1 == tmp_stop)
         {
            if (start)
               *start = -1;
            if (stop)
               *stop = -1;
         } else
         {
            if (start)
               *start = tmp_stop + 1;
            if (stop)
               *stop = bsp->length-1;
         }
         if (strand)
            *strand = tmp_strand;
         BioseqUnlockById(sip);
      }
   }
   return TRUE;
}


NLM_EXTERN Boolean AlnMgrGetNthUnalignedForNthRow(SeqAlignPtr sap, Int4 unaligned, Int4 row, Int4Ptr start, Int4Ptr stop)
{
   AMAlignIndexPtr  amaip;
   Int4             i;
   RowSourcePtr     rsp;
   Uint2            strand;
   Int4             tmp_start;
   Int4             tmp_stop;

   i = AlnMgrCheckAlignForParent(sap);
   if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (sap->type == SAT_PARTIAL || (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_SEGMENTED_MASTERSLAVE))
      {
         if (unaligned > amaip->alnsaps - 1)
            return FALSE;
         tmp_start = tmp_stop = 0;
         rsp = amaip->rowsource[row-1];
         strand = AlnMgrGetNthStrand(amaip->saps[rsp->which_saps[unaligned]-1], row);
         if (start)
         {
            if (strand == Seq_strand_minus)
            {
               AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[unaligned]-1], rsp->num_in_sap[unaligned], NULL, &tmp_start);
               *start = tmp_start + 1;
            } else
            {
               AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[unaligned-1]-1], rsp->num_in_sap[unaligned-1], NULL, &tmp_start);
               *start = tmp_start + 1;
            }
         }
         if (stop)
         {
            if (strand == Seq_strand_minus)
            {
               AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[unaligned-1]-1], rsp->num_in_sap[unaligned-1], &tmp_stop, NULL);
               *stop = tmp_stop - 1;
            } else
            {
               AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[unaligned]-1], rsp->num_in_sap[unaligned], &tmp_stop, NULL);
               *stop = tmp_stop - 1;
            }
         }
         if (tmp_start + 1 > tmp_stop - 1)
         {
            if (start)
               *start = -1;
            if (stop)
               *stop = -1;
         }
         return TRUE;
      } else
         return FALSE;
   } else
      return FALSE;
}


NLM_EXTERN Uint1 AlnMgrGetStrand(SeqAlignPtr sap, SeqIdPtr sip)
{
   Int4  i;

   i = AlnMgrGetNForSip(sap, sip);
   return (AlnMgrGetNthStrand(sap, i));
}

NLM_EXTERN Uint1 AlnMgrGetNthStrand(SeqAlignPtr sap, Int4 n)
{
   AMAlignIndexPtr  amaip;
   Int4             c;
   DenseSegPtr      dsp;
   Int4             m;
   SeqAlignPtr      salp;

   if (!sap || n < 1)
      return 0;
   if (sap->segtype != SAS_DENSEG)
   {
      if (sap->saip == NULL)
         return 0;
      amaip = (AMAlignIndexPtr)sap->saip;
      if (n > amaip->numrows)
         return 0;
      c = 0;
      while (amaip->rowsource[n-1]->which_saps[c] == 0)
      {
         c++;
         if (c >= amaip->alnsaps)
            return (Seq_strand_unknown);
      }
      salp = amaip->saps[amaip->rowsource[n-1]->which_saps[c]-1];
      dsp = (DenseSegPtr)salp->segs;
      m = amaip->rowsource[n-1]->num_in_sap[0];
      if (m > dsp->dim)
         return 0;
      return (dsp->strands[m-1]);
   } else
   {
      dsp = (DenseSegPtr)sap->segs;
      if (!dsp)
         return 0;
      if (n==0)
         return 0;
      if (dsp->strands)
         return (dsp->strands[n-1]);
      else
         return (Seq_strand_plus);
   }
}

NLM_EXTERN Int4 AlnMgrGetNForSip(SeqAlignPtr sap, SeqIdPtr sip)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             i;
   Int4             n;
   SeqIdPtr         sip_tmp;

   i = AlnMgrCheckAlignForParent(sap);
   if (i<0)
      return -1;
   if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      sip_tmp = amaip->ids;
      n = 0;
      while (sip_tmp)
      {
         n++;
         if (SeqIdComp(sip_tmp, sip) == SIC_YES)
            return n;
         sip_tmp = sip_tmp->next;
      }
   } else if (i == AM_CHILD)
   {
      dsp = (DenseSegPtr)(sap->segs);
      sip_tmp = dsp->ids;
      n = 0;
      while (sip_tmp)
      {
         n++;
         if (SeqIdComp(sip_tmp, sip) == SIC_YES)
            return n;
         sip_tmp = sip_tmp->next;
      }
   }
   return -1;
}

NLM_EXTERN Int4 AlnMgrGetNForSap(AMAlignIndexPtr amaip, SeqAlignPtr sap)
{
   Int4  i;

   if (sap == NULL || amaip == NULL)
      return -1;
   if (sap->saip->indextype != INDEX_SEGS)
      return -1;
   i = 0;
   while (i<amaip->alnsaps)
   {
      if (amaip->saps[i] == sap)
         return (i+1);
      i++;
   }
   return -1;
}


/********************************************************************************
*
*  AlnMgrGetAllNForSip is called in a while loop to return all the rows that a
*  seqid appears in in a given seqalign.  Use n = 0 to start, and then on
*  return, if the return is TRUE, n will be the row number of the next row
*  that the seqid appears in.  If the return is FALSE, either there was an
*  error or there are no (more) rows containing that seqid.
*
********************************************************************************/
NLM_EXTERN Boolean AlnMgrGetAllNForSip(SeqAlignPtr sap, SeqIdPtr sip, Int4Ptr n)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             i;
   SeqIdPtr         sip_tmp;

   if (sap == NULL || sip == NULL || n == NULL)
      return FALSE;
   if (sap->saip == NULL)
      return FALSE;
   if (sap->saip->indextype == INDEX_SEGS)
   {
      i = 1;
      dsp = (DenseSegPtr)sap->segs;
      sip_tmp = dsp->ids;
      while (i <= *n)
      {
         sip_tmp = sip_tmp->next;
         i++;
      }
      while (sip_tmp)
      {
         if (SeqIdComp(sip_tmp, sip) == SIC_YES)
         {
            *n = i;
            return TRUE;
         }
         i++;
         sip_tmp = sip_tmp->next;
      }
   } else if (sap->saip->indextype == INDEX_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      i = *n + 1;
      while (i <= amaip->numrows)
      {
         if (SeqIdComp(amaip->rowsource[i-1]->id, sip) == SIC_YES)
         {
            *n = i;
            return TRUE;
         }
         i++;
      }
   }
   return FALSE;
}

NLM_EXTERN Int4 AlnMgrGetSapForSip(AMAlignIndexPtr amaip, SeqIdPtr sip, Int4 which)
{
   Int4  i;
   Int4  j;
   Int4  n;

   i = 0;
   for (n=0; n<(amaip->numsaps); n++)
   {
      j = AlnMgrGetNForSip(amaip->saps[n], sip);
      if (j != -1)
      {
         if (i==which)
            return n;
         else
            i++;
      }
   }
   return -1;
}

/********************************************************************************
*
*  AlnMgrMapToBsqCoords returns the bioseq coordinate for an alignment
*  position.  If master is NULL, the alignment position is taken to be from
*  a flattened alignment; otherwise, the function returns the corresponding
*  position in the given master.
*
********************************************************************************/

NLM_EXTERN Int4 AlnMgrMapToBsqCoords(SeqAlignPtr sap, Uint4 pos, SeqIdPtr sip, SeqIdPtr master)
{
   DenseSegPtr  dsp;
   Int4         n;
   Int4         offset;
   SAIndexPtr   saip;
   Int4         start;

   if (!sap)
      return -1;
   if (sap->segtype == SAS_DENSEG)
   {
      saip = (SAIndexPtr)(sap->saip);
      dsp = (DenseSegPtr)(sap->segs);
      if (sip == NULL)
         n = saip->master;
      else
         n = AlnMgrGetNForSip(sap, sip);
      if (!master)
      {
         start = binary_search_on_uint4_list(saip->aligncoords, pos, dsp->numseg);
         offset = pos - saip->aligncoords[start];
         if (dsp->starts[(dsp->dim*start) + n - 1] == -1)
            return -1;
         else
            if (dsp->strands[(dsp->dim*start) + n - 1] != Seq_strand_minus)
               return (dsp->starts[(dsp->dim*start) + n - 1] + offset);
            else
               return (dsp->starts[(dsp->dim*start) + n - 1] + dsp->lens[start] - 1 - offset);
      } else
      {
      }
   } else if (sap->segtype == SAS_DISC)
   {
       SeqAlignPtr salp;
       salp = (SeqAlignPtr)sap->segs;
       if(salp->next==NULL)
           return AlnMgrMapToBsqCoords(salp, pos, sip, master);
   }
   return -1;
}

static Int4 AlnMgrMapSegmentCoords(SeqAlignPtr sap, Uint4 pos, Int4 row, SeqIdPtr master, Int4Ptr len)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             offset;
   SAIndexPtr       saip;
   Int4             start;

   if (sap == NULL || row < 0 || len == NULL)
      return -1;
   if (sap->saip == NULL)
      return -1;
   if (sap->saip->indextype == INDEX_SEGS)
   {
      saip = (SAIndexPtr)sap->saip;
      dsp = (DenseSegPtr)sap->segs;
      if (master == NULL)
      {
         start = binary_search_on_uint4_list(saip->aligncoords, pos, dsp->numseg);
         offset = pos - saip->aligncoords[start];
         *len = dsp->lens[start]-offset;
         if (dsp->starts[(dsp->dim*start) + row - 1] == -1)
            return -1;
         else
            if (dsp->strands[(dsp->dim*start) + row - 1] != Seq_strand_minus)
               return (dsp->starts[(dsp->dim*start) + row - 1] + offset);
            else
               return (dsp->starts[(dsp->dim*start) + row - 1] + dsp->lens[start] - 1 - offset);
      } else
      {
      }
   }
   return -1;
}


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
NLM_EXTERN Int4 AlnMgrMapRowCoords(SeqAlignPtr sap, Uint4 pos, Int4 row, SeqIdPtr master)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             offset;
   SAIndexPtr       saip;
   Int4             start;

   if (sap == NULL || row < 0)
      return -1;
   if (sap->saip == NULL)
      return -1;
   if (sap->saip->indextype == INDEX_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      if (row > amaip->numrows)
         return -1;
      if (master == NULL)
      {
      } else
      {
      }
   } else if (sap->saip->indextype == INDEX_SEGS)
   {
      saip = (SAIndexPtr)sap->saip;
      dsp = (DenseSegPtr)sap->segs;
      if (master == NULL)
      {
         start = binary_search_on_uint4_list(saip->aligncoords, pos, dsp->numseg);
         offset = pos - saip->aligncoords[start];
         if (dsp->starts[(dsp->dim*start) + row - 1] == -1)
            return -1;
         else
            if (dsp->strands[(dsp->dim*start) + row - 1] != Seq_strand_minus)
               return (dsp->starts[(dsp->dim*start) + row - 1] + offset);
            else
               return (dsp->starts[(dsp->dim*start) + row - 1] + dsp->lens[start] - 1 - offset);
      } else
      {
      }
   }
   return -1;
}


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
NLM_EXTERN Int4 AlnMgrMapBioseqToSeqAlign(SeqAlignPtr sap, Uint4 pos, Int4 row_num, SeqIdPtr master)
{
   AMAlignIndexPtr  amaip;
   Boolean          done;
   DenseSegPtr      dsp;
   Int4             i;
   Int4             j;
   Int4             k;
   Uint2            L;
   Int4             m;
   Uint2            mid;
   Uint1            n;
   Int4             offset;
   Uint2            R;
   SAIndexPtr       saip;
   SASeqDatPtr      ssdp;
   Int4             start;
   Int4             stop;

   if (sap == NULL || row_num < 0)
      return -1;
   AlnMgrGetNthSeqRangeInSA(sap, row_num, &start, &stop);
   if (pos < start || pos > stop)
      return -2;
   i = AlnMgrCheckAlignForParent(sap);
   if (i == AM_CHILD)
   {
      if (master == NULL)
      {
         saip = (SAIndexPtr)sap->saip;
         ssdp = saip->ssdp[row_num-1];
         if (ssdp == NULL)
            return -1;
         dsp = (DenseSegPtr)sap->segs;
         L = 0;
         R = ssdp->numsect - 1;
         while (L < R)
         {
            mid = (L + R)/2;
            if (dsp->starts[ssdp->sect[mid+1]+row_num-1] <= pos)
            {
               L = mid + 1;
            } else
               R = mid;
         }
         n = AlnMgrGetNthStrand(sap, row_num);
         offset = pos - dsp->starts[ssdp->sect[R]+row_num-1];
         if (offset > dsp->lens[ssdp->sect[R]] || offset < 0)
            return -2;
         if (n != Seq_strand_minus)
            return (saip->aligncoords[R] + offset);
         else
            return (saip->aligncoords[R] + dsp->lens[R] - 1 - offset);
      } else
      {
      }
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      if (row_num > amaip->numrows)
         return -1;
      j = k = 0;
      m = -1;
      done = FALSE;
      while (!done && j < amaip->numseg)
      {
         k = AlnMgrMapRowCoords(sap, amaip->aligncoords[j], row_num, master);
         if (k == -1)
            j++;
         else if (k > pos)
            done = TRUE;
         else if (k <= pos)
         {
            m = j;
            offset = pos - k;
         }
      }
      if (m == -1)
         return -2;
      n = AlnMgrGetNthStrand(sap, row_num);
      if (n != Seq_strand_minus)
      {
         return (amaip->aligncoords[m] + offset);
      } else
      {
         return (amaip->aligncoords[m] + amaip->lens[m] - 1 - offset);
      }
   } else
      return -1;
   return -1;
}


/***********************************************************************
*
*  AlnMgrMakeFakeMultiple calls AlnMgrCheckOverlapping to decide whether
*  an alignment is linear.  Then, if possible, it calls AlnMgrMakeAlignCoords
*  to create alignment coordinates across all children contained in the
*  parent.  (MULT)
*
***********************************************************************/
NLM_EXTERN Boolean AlnMgrMakeFakeMultiple(SeqAlignPtr sap)
{
   AMAlignDatPtr    amadp;
   AMAlignIndexPtr  amaip;
   Int4             i;
   Int4             j;
   Boolean          ms;
   Int4             n;
   Boolean          retval;

   retval = FALSE;
   if (!sap)
      return retval;
   i = AlnMgrCheckAlignForParent(sap);
   if (i<0)
   {
      return retval;
   }
   if (i==AM_PARENT)
   {
      n = AlnMgrCheckOverlapping(sap);
      if (n == NO_OVERLAP)
      {
         sap->type = SAT_PARTIAL;
         amaip = (AMAlignIndexPtr)sap->saip;
         if (amaip->saps)
            MemFree(amaip->saps);
         amaip->saps = AlnMgrSortSeqAligns((SeqAlignPtr)(sap->segs), AlnMgrFindFirst, amaip, &amaip->numsaps);
         amaip->alnsaps = amaip->numsaps;
         amaip->starts = (Int4Ptr)MemNew((amaip->alnsaps)*(amaip->alnsaps)*sizeof(Int4));
         amaip->lens = (Int4Ptr)MemNew((amaip->alnsaps)*sizeof(Int4));
         amaip->ulens = (Int4Ptr)MemNew((amaip->alnsaps)*sizeof(Int4));
         amaip->numseg = amaip->alnsaps;
         for (j=0; j<(amaip->alnsaps); j++)
         {
            amaip->lens[j] = AlnMgrGetAlnLength(amaip->saps[j], FALSE);
            amaip->starts[j] = 0;
         }
         AlnMgrMakeAlignCoords(sap);
         if (!AlnMgrGetRowsForPartial(sap))
            return retval;
         for (j=0; j<(amaip->alnsaps-1); j++)
         {
            amaip->ulens[j] = AlnMgrGetMaxUnalignedLength(amaip->saps[j], amaip->saps[j+1]);
         }
         retval = TRUE;
      } else /*should add function to check for pairwise multiple vs. diags*/
      {
         amaip = (AMAlignIndexPtr)sap->saip;
         if (amaip->saps)
            MemFree(amaip->saps);
         sap->master = AlnMgrFindMaster(sap);
         amaip->alnsaps = amaip->numsaps;
         ms = FALSE;
         ms = AlnMgrCheckRealMaster(sap, sap->master);
         if (sap->master && ms == TRUE)
         {
            retval = TRUE;
            AlnMgrSetMaster(sap, sap->master);
            AlnMgrMakeMasterPlus(sap);
            n = AlnMgrGetNForSip(sap, sap->master);
            sap->type = SAT_MASTERSLAVE;
            amaip->master = 1;
            amaip->numseg = AlnMgrGetMaxSegments((SeqAlignPtr)(sap->segs));
            amaip->alnsaps = amaip->numsaps;
            amaip->lens = (Int4Ptr)MemNew((amaip->numseg)*sizeof(Int4));
            amadp = amaip->amadp[n-1];
            amaip->saps = amadp->saps;
            if (amaip->numsaps < amaip->numbsqs)
            {
               amaip->ids = SeqIdSetFree(amaip->ids);
               amaip->ids = AlnMgrPropagateSeqIdsBySapList(amaip);
               if (!AlnMgrMergeIntoMSMultByMaster(amaip, amaip->lens, &amaip->numseg))
                  retval = FALSE;
               amaip->starts = (Int4Ptr)MemNew((amaip->numseg)*(amaip->numsaps)*sizeof(Int4));
               amaip->aligncoords = (Uint4Ptr)MemNew((amaip->numseg)*sizeof(Uint4));
               if (!AlnMgrFillInStarts(amadp->saps, amaip->starts, amaip->numseg, amaip->lens, amaip->numsaps, amaip->aligncoords))
                  retval = FALSE;
               if (amaip->numseg > 1)
                  amaip->numseg -= 1;
               if (!AlnMgrMakeMultSegments(amaip))
                  retval = FALSE;
               if (!AlnMgrGetRowsForMasterSlave(sap))
                  retval = FALSE;
            } else
               retval = FALSE;
         }
         if (retval == FALSE && sap->master != NULL)
         {
            if (AlnMgrMakeSegmentedMasterSlave(sap))
            {
               sap->type = SAT_MASTERSLAVE;
               amaip->ids = SeqIdSetFree(amaip->ids);
               amaip->ids = AlnMgrPropagateSeqIdsByRow(amaip);
               retval = TRUE;
            } else
            {
               if (AlnMgrForceMasterSlave(sap))
               {
                  amaip->ids = SeqIdSetFree(amaip->ids);
                  amaip->ids = AlnMgrPropagateSeqIdsByRow(amaip);
                  amaip->mstype = AM_MASTERSLAVE;
                  retval = TRUE;
               }
            }
         } else
         {
            amaip->mstype = AM_MASTERSLAVE;
         }
      }
   }
   return retval;
}


NLM_EXTERN void AlnMgrMakeAlignCoords(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   Int4             i;
   Int4             j;

   i = AlnMgrCheckAlignForParent(sap);
   if (i < 0 || i==AM_CHILD)
      return;
   amaip = (AMAlignIndexPtr)(sap->saip);
   if (!amaip->saps)
      return;
   amaip->aligncoords = (Uint4Ptr)MemNew((amaip->alnsaps)*sizeof(Uint4));
   amaip->aligncoords[0] = 0;
   for (j=0; j<((amaip->alnsaps)-1); j++)
   {
      amaip->aligncoords[j+1] = AlnMgrGetAlnLength(amaip->saps[j], FALSE) + amaip->aligncoords[j];
   }
   return;
}

NLM_EXTERN Boolean AlnMgrMergeIntoMSMultByMaster(AMAlignIndexPtr amaip, Int4Ptr lens, Uint4Ptr numseg)
{
   Uint4          count;
   DenseSegPtr    dsp;
   Int4           gap;
   Int4           i;
   Int4           j;
   Int4           n;
   Boolean        retval;
   SAIndexPtr     saip;
   AMTinyInfoPtr  tip;
   AMTinyInfoPtr  PNTR tiparray;

   retval = FALSE;
   if (numseg == NULL)
      return retval;
   tiparray = (AMTinyInfoPtr PNTR)MemNew((*numseg+1)*sizeof(AMTinyInfoPtr));
   j = 0;
   count = 0;
   for (i=0; i<(amaip->alnsaps); i++)
   {
      dsp = (DenseSegPtr)(amaip->saps[i]->segs);
      saip = (SAIndexPtr)amaip->saps[i]->saip;
      gap = 0;
      for (n=0; n<(dsp->numseg); n++)
      {
         if (dsp->starts[n*(dsp->dim)+saip->master-1] != -1)
         {
            tip = (AMTinyInfoPtr)MemNew(sizeof(AMTinyInfo));
            tip->start = dsp->starts[n*(dsp->dim)+saip->master-1];
            tip->which = i+1;
            tip->numgap = gap;
            tiparray[j] = tip;
            j++;
            count++;
            gap = 0;
         } else
         {
            gap++;
         }
      }
      tip = (AMTinyInfoPtr)MemNew(sizeof(AMTinyInfo));
      AlnMgrGetNthSeqRangeInSA(amaip->saps[i], saip->master, NULL, &tip->start);
      tip->start += 1;
      tip->which = i+1;
      tip->numgap = gap;
      tiparray[j] = tip;
      j++;
      count++;
   }
   *numseg = count;
   HeapSort((Pointer)tiparray, (size_t)(*numseg), sizeof(AMTinyInfoPtr), AlnMgrCompareTips);
   *numseg = count-1;
   count = 0;
   for (i=0; i<=(*numseg); i++)
   {
      if (count!=0 && (tiparray[i]->start == lens[count-1])) 
         count--;
      for (j=1; j<=(tiparray[i]->numgap); j++)
      {
         lens[count] = -(tiparray[i]->which);
         count++;
      }
      lens[count] = tiparray[i]->start;
      count++;  
   }
   for (i=0; i<=(*numseg); i++)
   {
      MemFree(tiparray[i]);
   }
   MemFree(tiparray);
   *numseg = count;
   return TRUE;
}

NLM_EXTERN Boolean AlnMgrMergeSegments(Int4Ptr lens, SeqAlignPtr sap, SeqIdPtr master, Int4Ptr where, Int4 which)
{
   DenseSegPtr  dsp;
   Boolean      found;
   Int4         i;
   Int4         j;
   Int4         n;
   Int4         num;
   Int4         r;
   Boolean      retval;
   Int4         s;
   SAIndexPtr   saip;
   Int4Ptr      tmp;
   Int4         z;

   retval = FALSE;
   if (!sap || !master || !lens)
      return retval;
   if (!where)
      return retval;
   n = AlnMgrGetNForSip(sap, master);
   if (n<0)
      return retval;
   if (sap->segtype == SAS_DENSEG)
   {
      dsp = (DenseSegPtr)(sap->segs);
      if (!dsp)
         return retval;
      saip = (SAIndexPtr)(sap->saip);
      if (!saip)
         return retval;
   } else
   {
      return retval;
   }
   if (*where == 0)
   {
      for(j=0; j<(dsp->numseg); j++)
      {
         if (dsp->starts[(j*(dsp->dim)) + n - 1] < 0)
         {
            s = -(which);
         } else
         {
            s = dsp->starts[(j*(dsp->dim)) + n - 1];
         }
         lens[*where] = s;
         *where = *where + 1;
      }
      AlnMgrGetNthSeqRangeInSA(sap, saip->master, NULL, &lens[dsp->numseg]);
      lens[dsp->numseg] += 1;
      *where = *where + 1;
   } else
   {
      tmp = (Int4Ptr)MemNew((dsp->numseg+1)*sizeof(Int4));
      for(j=0; j<(dsp->numseg); j++)
      {
         if (dsp->starts[(j*(dsp->dim)) + n - 1] < 0)
         {
            s = -(which);
         } else
         {
            s = dsp->starts[(j*(dsp->dim)) + n - 1];
         }
         tmp[j] = s;
      }
      AlnMgrGetNthSeqRangeInSA(sap, saip->master, NULL, &tmp[dsp->numseg]);
      tmp[dsp->numseg] += 1;
      s = 0;
      for (j=0; j<=(dsp->numseg); j++)
      {
         num = 0;
         while (tmp[j] < 0 && num<(dsp->numseg))
         {
            num++;
            j++;
         }
         num++;
         found = FALSE;
         for (i=s; !found && i<*where; i++)
         {
            r = 0;
            if (lens[i] < 0)
            {
            } else if (tmp[j] < lens[i])
            {
               if (i>0)
               {
                  while (((i-r-1)>=0) && (lens[i-r-1] < 0))
                  {
                     r++;
                  }
               }
               s = i;
               for (z = *where-1; z >= i-r; z--)
               {
                  lens[z+num] = lens[z];
               }
               for (z = num; z > 0; z--)
               {
                  lens[i-r] = tmp[j-z+1];
                  i++;
               }
               found = TRUE;
               *where = *where + num;
            } else if (tmp[j] == lens[i])
            {
               s = i;
               for (z = *where-1; z >= i; z--)
               {
                  lens[z+num-1] = lens[z];
               }
               for (z = num-1; z > 0; z--)
               {
                  lens[i] = tmp[j-z];
                  i++;
               }
               found = TRUE;
               *where = *where + num - 1;
            }
         }
         if (!found)
         {
            s = *where;
            for (z = *where+num-1; z >= *where; z--)
            {
               lens[z+num] = lens[z];
            }
            for (z = num-1; z >= 0; z--)
            {
               lens[i] = tmp[j-z];
               i++;
            }
            found = TRUE;
            *where = *where + num;
         }
      }
      MemFree(tmp);
   }
   retval = TRUE;
   return retval;
}


NLM_EXTERN Boolean AlnMgrFillInStarts(SeqAlignPtr PNTR saparray, Int4Ptr starts, Int4 numseg, Int4Ptr lens, Int4 numsaps, Uint4Ptr aligncoords)
{
   Int4Ptr  alnlen;
   Boolean  done;
   Int4     gap_pos;
   Int4     i;
   Int4     j;
   Int4     length;
   Boolean  retval;

   retval = FALSE;
   for (i=0; i<numsaps; i++)
   {
      gap_pos = 0;
      for (j=0; j<numseg; j++)
      {
         if(lens[j] >= 0)
         {
            starts[(numsaps*j)+i] = AlnMgrGetStartFromMaster(saparray[i], lens[j]);
         } else
         {
            if (lens[j] == -(i+1))
            {
               starts[(numsaps*j)+i] = AlnMgrGetMasterGapStartForSeg(saparray[i], gap_pos, &aligncoords[j]);
               gap_pos += 1;
            } else
            {
               starts[(numsaps*j)+i] = -1;
            }
         }
      }
   }
   if (!AlnMgrReconcileGaps(lens, aligncoords, numseg))
      return retval;
   alnlen = (Int4Ptr)MemNew(numsaps*sizeof(Int4));
   for (i=0; i<numsaps; i++)
   {
      alnlen[i] = AlnMgrGetAlnLength(saparray[i], FALSE);
   }
   for (i=0; i<numsaps; i++)
   {
      length = 0;
      done = FALSE;
      for (j=0; j<numseg; j++)
      {
         if (starts[(numsaps*j)+i] == -2)
         {
            if (length > 0)
            {
               if (lens[j]+length-1 < alnlen[i])
               {
                  starts[(numsaps*j)+i] = length;
                  length += lens[j];
               } else
               {
                  done = TRUE;
               }
            }
         } else if (starts[(numsaps*j)+i] == -1)
         {
            if (length == 0)
               starts[(numsaps*j)+i] = -2;
            else if (done)
               starts[(numsaps*j)+i] = -2;
         } else
         {
            length = starts[(numsaps*j)+i] + lens[j];
         }
      }
   }
   j = 0;
   numseg -= 1;
   done = FALSE;
   if (numseg != 0)
      done = FALSE;
   else
      done = TRUE;
   for (i=(numsaps*(numseg-1)); (!done && i<(numsaps*numseg)); i++)
   {
      if (starts[i] >= 0)
      {
         done = TRUE;
         lens[numseg-1] = alnlen[j]-starts[i];
      }
      else
         j++;
   }
   fflush(stdout);
   MemFree(alnlen);
   retval = TRUE;
   return retval;
}

NLM_EXTERN Int4 AlnMgrGetStartFromMaster(SeqAlignPtr sap, Int4 pos)
{
   DenseSegPtr  dsp;
   SAIndexPtr   saip;
   Int4         start;

   saip = (SAIndexPtr)(sap->saip);
   dsp = (DenseSegPtr)(sap->segs);
   start = binary_search_segment_array(saip->ssdp[saip->master-1], pos, dsp->dim, saip->master - 1, (DenseSegPtr)sap->segs);
   if (dsp->starts[(start*dsp->dim)+saip->master-1] != pos)
   {
      return -2;
   } else
   {
      return (saip->aligncoords[start]);
   }
}

NLM_EXTERN Uint4 AlnMgrGetMasterGapStartForSeg(SeqAlignPtr sap, Int4 which_gap, Uint4Ptr aligncoord)
{
   DenseSegPtr  dsp;
   SAIndexPtr   saip;

   saip = (SAIndexPtr)(sap->saip);
   dsp = (DenseSegPtr)(sap->segs);
   if (aligncoord)
      *aligncoord = dsp->lens[saip->ssdp[saip->master-1]->unsect[which_gap]];
   return saip->aligncoords[saip->ssdp[saip->master-1]->unsect[which_gap]];
}


NLM_EXTERN Boolean AlnMgrReconcileGaps(Int4Ptr lens, Uint4Ptr aligncoords, Int4 num)
{
   Int4  i;
   Int4  j;
   Int4  r;

   for (i=0; i<num; i++)
   {
      if (lens[i] < 0)
      {
         r = 1;
         while (lens[i+r] < 0)
         {
            r++;
         }
         lens[i] = lens[i+r];
         for (j=i+1; j<num; j++)
         {
            if (lens[j] >= 0)
               lens[j] = lens[j] + aligncoords[i];
         }
      }
   }
   for (i=0; i<num; i++)
   {
      aligncoords[i] = lens[i];
   }
   for (i=0; i<num-1; i++)
   {
      lens[i] = lens[i+1] - lens[i];
   }
   return TRUE;
}

NLM_EXTERN Boolean AlnMgrMakeMultSegments(AMAlignIndexPtr amaip)
{
   Int4      i;
   Int4      j;
   Uint2     n;
   Boolean   retval;
   Uint2Ptr  segments;
   Uint2Ptr  tmp;

   retval = FALSE;
   tmp = (Uint2Ptr)MemNew((amaip->numseg)*sizeof(Uint2));
   for (i=0; i<amaip->numsaps; i++)
   {
      n = 0;
      for (j=0; j<amaip->numseg; j++)
      {
         if (amaip->starts[((amaip->numsaps)*j)+i] >= 0)
         {
            tmp[n] = j;
            n++;
         }
      }
      segments = (Uint2Ptr)MemNew(n*sizeof(Uint2));
      for (j=0; j<n; j++)
      {
         segments[j] = tmp[j];
      }
      if (!amaip->amadp[i])
         return retval;
      amaip->amadp[i]->segments = segments;
      amaip->amadp[i]->numseg = n;
      amaip->amadp[i]->numseg = n;
   }
   MemFree(tmp);
   retval = TRUE;
   return retval;
}

NLM_EXTERN Int4 AlnMgrCheckOverlapping(SeqAlignPtr sap)
{
   AMAlignDatPtr    amadp;
   AMAlignIndexPtr  amaip;
   Int4             end;
   Int4             c;
   Int4             i;
   Int4             j;
   Int4             n;
   Int4             prevstrand;
   SeqIdPtr         sip;
   Int4             start;
   Int4             stop;
   Uint2            strand;

   i = AlnMgrCheckAlignForParent(sap);
   if (i<0)
      return CHECK_ERROR;
   else if (i==AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      if (amaip->numsaps == 1)
         return 1;
      sip = amaip->ids;
      for (j=0; j<(amaip->numbsqs); j++)
      {
         end = -1;
         amadp = amaip->amadp[j];
         prevstrand = -1;
         for (c=0; c<(amadp->numsaps); c++)
         {
            n = AlnMgrGetNForSip(amadp->saps[c], sip);
            strand = AlnMgrGetNthStrand(amadp->saps[c], n);
            if (prevstrand != -1)
            {
               if (strand != prevstrand)
                  return j;
            } else
               prevstrand = strand;
            AlnMgrGetNthSeqRangeInSA(amadp->saps[c], n, &start, &stop);
            if (strand != Seq_strand_minus)
            {
               if (start <= end)
                  return j;
               else
                  end = stop;
            } else
            {
               if (end != -1 && stop >= end)
                  return j;
               else
                  end = start;
            }
         }
         sip = sip->next;
      }
   } else if (i==AM_CHILD)
   {
      return NO_OVERLAP;
   }
   return NO_OVERLAP;
}

/*****************************************************************************
*
*  AlnMgrGetMaxSegments simply adds up the number of segments for each
*  SeqAlign in a linked list, to get the maximum number of segments 
*  for the merge of the list (for memory allocation in AlnMgrMakeFakeMultiple).
*
******************************************************************************/

NLM_EXTERN Int4 AlnMgrGetMaxSegments(SeqAlignPtr sap)
{
   DenseSegPtr  dsp;
   Int4         ernie; /* the running total, also a happy hamster */

   ernie = 0;
   while (sap)
   {
      if (sap->segtype == SAS_DENSEG)
      {
         dsp = (DenseSegPtr)(sap->segs);
         ernie += dsp->numseg;
      } else if (sap->segtype == SAS_STD)
      {
         ernie += 1;
      } else
         return 0;
      sap = sap->next;
      ernie += 1;
   }
   return ernie;
}

/*******************************************************************************
*
*  Row Management functions:
*
*******************************************************************************/
NLM_EXTERN Int4 AlnMgrGetNumRows(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;

   if (sap == NULL || sap->saip == NULL)
      return -1;
   if (sap->saip->indextype == INDEX_SEGS)
   {
      dsp = (DenseSegPtr)sap->segs;
      if (dsp == NULL)
         return -1;
      return (dsp->dim);
   } else if (sap->saip->indextype == INDEX_PARENT)
   {
      if ((amaip = (AMAlignIndexPtr)sap->saip) == NULL)
         return -1;
      if (amaip->numrows)
         return (amaip->numrows);
   }
   return 0;
}

NLM_EXTERN Int4 AlnMgrGetMaxRowsForParentPartial(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   Int4             i;
   Int4             j;
   Int4             max;

   if (sap == NULL || sap->saip == NULL)
      return -1;
   max = -1;
   if (sap->saip->indextype == INDEX_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      for (i=0; i<(amaip->alnsaps); i++)
      {
         j = AlnMgrGetNumRows(amaip->saps[i]);
         if (j==-1)
            return -1;
         if (j>max)
            max = j;
      }
   }
   return max;
}

NLM_EXTERN Boolean AlnMgrGetRowsForPartial(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   Int4             curr;
   DenseSegPtr      dsp;
   Boolean          found;
   Int4             i;
   Int4             j;
   Int4             k;
   Boolean          retval;
   RowSourcePtr     PNTR rowsource;
   RowSourcePtr     rsp;
   SeqAlignPtr      salp;
   SeqIdPtr         sip;

   retval = FALSE;
   if (sap == NULL || sap->saip == NULL)
      return retval;
   if (sap->saip->indextype != INDEX_PARENT)
      return retval;
   if (sap->type != SAT_PARTIAL)
      return retval;
   amaip = (AMAlignIndexPtr)sap->saip;
   i = AlnMgrGetMaxRowsForParentPartial(sap);
   if (i < 0)
      return retval;
   else
      amaip->numrows = i;
   rowsource = (RowSourcePtr PNTR)MemNew((amaip->numrows)*sizeof(RowSourcePtr));
   curr = -1;
   for (i=0; i<(amaip->alnsaps); i++)
   {
      salp = amaip->saps[i];
      dsp = (DenseSegPtr)salp->segs;
      sip = dsp->ids;
      for (j=0; j<(dsp->dim); j++)
      {
         found = FALSE;
         k = 0;
         while (!found && k <= curr)
         {
            if (SeqIdComp(sip, rowsource[k]->id) == SIC_YES)
               found = TRUE;
            else
               k++;
         }
         if (!found)
         {
            curr++;
            rsp = RowSourceNew();
            rsp->which_saps = (Uint4Ptr)MemNew((amaip->alnsaps)*sizeof(Uint4));
            rsp->num_in_sap = (Uint4Ptr)MemNew((amaip->alnsaps)*sizeof(Uint4));
            rowsource[curr] = rsp;
            rsp->id = SeqIdDup(sip);
            rsp->which_saps[rsp->numsaps] = i+1;
            rsp->num_in_sap[rsp->numsaps] = AlnMgrGetNForSip(salp, sip);
            (rsp->numsaps)++;
         } else
         {
            rowsource[k]->which_saps[rowsource[k]->numsaps] = i+1;
            rowsource[k]->num_in_sap[rowsource[k]->numsaps] = AlnMgrGetNForSip(salp, sip);
            (rowsource[k]->numsaps)++;
         }
         sip = sip->next;
      }
   }
   amaip->numrows = curr+1;
   amaip->rowsource = rowsource;
   amaip->master = -2;
   return TRUE;
}

NLM_EXTERN Boolean AlnMgrGetRowsForMasterSlave(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             i;
   Int4             j;
   Int4             k;
   Boolean          retval;
   RowSourcePtr     PNTR rowsource;
   RowSourcePtr     rsp;
   SAIndexPtr       saip;
   SeqAlignPtr      salp;
   SeqIdPtr         sip;

   retval = FALSE;
   if (sap == NULL || sap->saip == NULL)
      return retval;
   if (sap->saip->indextype != INDEX_PARENT)
      return retval;
   if (sap->type != SAT_MASTERSLAVE)
      return retval;
   amaip = (AMAlignIndexPtr)sap->saip;
   i = 1;
   salp = (SeqAlignPtr)sap->segs;
   while (salp)
   {
      j = AlnMgrGetNumRows(salp);
      if (j < 0)
         return retval;
      else
         i += (j-1); /*don't count the master over and over*/
      salp = salp->next;
   }
   rowsource = (RowSourcePtr PNTR)MemNew((i+1)*sizeof(RowSourcePtr));
   rsp = RowSourceNew();
   rsp->id = SeqIdDup(sap->master);
   rsp->which_saps = (Uint4Ptr)MemNew((amaip->alnsaps)*sizeof(Uint4));
   rsp->num_in_sap = (Uint4Ptr)MemNew((amaip->alnsaps)*sizeof(Uint4));
   rsp->numsaps = amaip->alnsaps;
   rowsource[0] = rsp;
   amaip->numrows = 1;
   for (j=0; j<(amaip->alnsaps); j++)
   {
      salp = amaip->saps[j];
      saip = (SAIndexPtr)salp->saip;
      dsp = (DenseSegPtr)(salp->segs);
      sip = dsp->ids;
      k=1;
      while (sip)
      {
         if (k != saip->master)
         {
            rsp = RowSourceNew();
            rsp->id = SeqIdDup(sip);
            rsp->which_saps = (Uint4Ptr)MemNew(sizeof(Uint4));
            rsp->num_in_sap = (Uint4Ptr)MemNew(sizeof(Uint4));
            rsp->numsaps = 1;
            rsp->which_saps[0] = j+1;
            rsp->num_in_sap[0] = k;
            rowsource[amaip->numrows] = rsp;
            amaip->numrows++;
         } else
         {
            rowsource[0]->which_saps[j] = j+1;
            rowsource[0]->num_in_sap[j] = k;
            amaip->master = 1;
         }
         k++;
         sip = sip->next;
      }
   }
   amaip->rowsource = rowsource;
   return TRUE;
}


/*******************************************************************************
*
*  AlnMgrFindMaster returns the (duplicated) SeqIdPtr of the first bioseq 
*  that is present in every child alignment, unless the sap->master field 
*  is set in the child alignments, in which case that SeqIdPtr is returned 
*  (if it's the same in all children).
*
*******************************************************************************/

NLM_EXTERN SeqIdPtr AlnMgrFindMaster(SeqAlignPtr sap)
{
   AMAlignDatPtr    amadp;
   AMAlignIndexPtr  amaip;
   Boolean          done;
   Int4             i;
   SeqAlignPtr      salp;
   SeqIdPtr         sip;

   i = AlnMgrCheckAlignForParent(sap);
   if (i<0)
      return NULL;
   else if (i==AM_CHILD)
   {
      return SeqIdDup(sap->master);
   } else if (i==AM_PARENT)
   {
      salp = (SeqAlignPtr)(sap->segs);
      sip = NULL;
      done = FALSE;
      while (salp && !done)
      {
         if (salp->master)
         {
            if (sip)
            {
               if (SeqIdComp(sip, salp->master) != SIC_YES)
               {
                  done = TRUE;
                  sip = SeqIdFree(sip);
               }
            } else
            {
               sip = SeqIdDup(salp->master);
            }
         } else
         {
            sip = SeqIdFree(sip);
            done = TRUE;
         }
         salp = salp->next;
      }
      if (sip)
         return sip;
      amaip = (AMAlignIndexPtr)(sap->saip);
      sip = amaip->ids;
      for (i=0; i<(amaip->numbsqs); i++)
      {
         amadp = amaip->amadp[i];
         if (!amadp || !sip)
            return NULL;
         else 
         {
            if (amadp->numsaps == amaip->numsaps)
               return (SeqIdDup(sip));
         }
         sip = sip->next;
      }
      return NULL;
   }
   return NULL;
}


/*******************************************************************************
*
*  AlnMgrCheckRealMaster makes sure that the master seqid given appears
*  once and only once in each seqalign in the set if a parent is given,
*  or once and only one in the seqalign if a child is given.
*
*******************************************************************************/
NLM_EXTERN Boolean AlnMgrCheckRealMaster(SeqAlignPtr sap, SeqIdPtr master)
{
   DenseSegPtr  dsp;
   Int4         i;
   Boolean      retval;
   SeqAlignPtr  salp;
   SeqIdPtr     sip;

   retval = FALSE;
   if (!sap || !master)
      return retval;
   if (sap->segtype == SAS_DISC)
   {
      salp = (SeqAlignPtr)sap->segs;
      while (salp)
      {
         dsp = (DenseSegPtr)salp->segs;
         sip = dsp->ids;
         i = 0;
         while (sip)
         {
            if (SeqIdComp(sip, master) == SIC_YES)
            {
               i++;
               if (i > 1)
                  return retval;
            }
            sip = sip->next;
         }
         salp = salp->next;
      }
   } else if (sap->segtype == SAS_DENSEG)
   {
      dsp = (DenseSegPtr)sap->segs;
      sip = dsp->ids;
      i = 0;
      while (sip)
      {
         if (SeqIdComp(sip, master) == SIC_YES)
         {
            i++;
            if (i > 1)
               return retval;
         }
         sip = sip->next;
      }
   }
   return TRUE;
}

NLM_EXTERN Boolean AlnMgrMakeSegmentedMasterSlave(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   AMmsmsPtr        ams;
   AMmsmsPtr        PNTR amsarray;
   AMmsmsPtr        ams_head;
   AMmsmsPtr        ams_master;
   AMmsmsPtr        ams_mtmp;
   AMmsmsPtr        ams_tmp;
   Int4             c;
   Boolean          done;
   DenseSegPtr      dsp;
   Boolean          found;
   Int4             i;
   Int4             j;
   Int4             n;
   Int4             max;
   Boolean          ok;
   RowSourcePtr     rsp;
   Int4             rspnum;
   SAIndexPtr       saip;
   SeqAlignPtr      salp;
   SeqIdPtr         sip;
   AMsiplistPtr     siplist;
   AMsiplistPtr     siplist_new;
   AMsiplistPtr     siplist_tmp;
   Int4             sstart;
   Int4             sstop;
   Int4             start;
   Int4             stop;
   Int4Ptr          tmparray;

   if (sap == NULL)
      return FALSE;
   amaip = (AMAlignIndexPtr)sap->saip;
   if (amaip == NULL)
      return FALSE;
   if (amaip->master < 0)
      return FALSE;
   ams_head = NULL;
   n = 0;
   for (i=0; i<(amaip->numsaps); i++)
   {
      salp = amaip->saps[i];
      saip = (SAIndexPtr)(salp->saip);
      if (saip->master < 0)
         return FALSE;
      AlnMgrGetNthSeqRangeInSA(salp, saip->master, &start, &stop);
      dsp = (DenseSegPtr)salp->segs;
      sip = dsp->ids;
      j = 1;
      while (sip != NULL)
      {
         if (j != saip->master)
         {
            n++;
            ams = (AMmsmsPtr)MemNew(sizeof(AMmsms));
            ams->start = start; 
            ams->stop = stop;
            ams->sap = salp;
            ams->nsap = i+1;
            ams->sip = sip;
            ams->n = j;
            AlnMgrGetNthSeqRangeInSA(salp, j, &sstart, &sstop);
            ams->sstart = sstart;
            ams->sstop = sstop;
            ams->strand = AlnMgrGetNthStrand(salp, j);
            if (ams_head == NULL)
            {
               ams_head = ams_tmp = ams;
            } else
            {
               ams_tmp->next = ams;
               ams_tmp = ams;
            }
         }
         sip = sip->next;
         j++;
      }
   }
   ams_head = am_sort_ammsms(ams_head, n);
   ams_master = NULL;
   ams = ams_head;
   n = 0;
   while (ams)
   {
      if (ams_master)
      {
         ams_mtmp = ams_master;
         found = FALSE;
         while (!found && ams_mtmp)
         {
            if (ams->start == ams_mtmp->start && ams->stop == ams_mtmp->stop)
            {
               found = TRUE;
               ams->masternum = ams_mtmp->masternum;
               ams_mtmp->count++;
            }
            else
               ams_mtmp = ams_mtmp->next;
         }
         if (!found)
         {
            n++;
            ams_tmp = (AMmsmsPtr)MemNew(sizeof(AMmsms));
            ams_tmp->start = ams->start;
            ams_tmp->stop = ams->stop;
            ams_tmp->sap = ams->sap;
            ams_tmp->nsap = ams->nsap;
            ams_tmp->sip = sap->master;
            ams_tmp->count = 1;
            ams_tmp->masternum = ams->masternum = n;
            saip = (SAIndexPtr)(ams->sap->saip);
            ams_tmp->n = saip->master;
            ams_tmp->next = ams_master;
            ams_master = ams_tmp;
         }
      } else
      {
         n++;
         ams_tmp = (AMmsmsPtr)MemNew(sizeof(AMmsms));
         ams_tmp->start = ams->start;
         ams_tmp->stop = ams->stop;
         ams_tmp->sap = ams->sap;
         ams_tmp->nsap = ams->nsap;
         ams_tmp->sip = sap->master;
         ams_tmp->count = 1;
         ams_tmp->masternum = ams->masternum = n;
         saip = (SAIndexPtr)(ams->sap->saip);
         ams_tmp->n = saip->master;
         ams_master = ams_tmp;
      }
      ams = ams->next;
   }
   ams_master = am_sort_masterams(ams_master, n);
   max = c = 0;
   ams = ams_master;
   ams_tmp = NULL;
   amsarray = (AMmsmsPtr PNTR)MemNew((n+1)*sizeof(AMmsmsPtr));
   while (ams)
   {
      amsarray[c] = ams;
      if (ams_tmp)
      {
         if (ams->start <= ams_tmp->stop)
         {
            MemFree(amsarray);
            return FALSE; /* add code here to compress all lines??? */
         }
      }
      max += ams->count;
      c++;
      ams_tmp = ams;
      ams = ams->next;
   }
   amaip->mstype = AM_SEGMENTED_MASTERSLAVE;
   amaip->rowsource = (RowSourcePtr PNTR)MemNew((max+1)*sizeof(RowSourcePtr));
   if (amaip->aligncoords)
      MemFree(amaip->aligncoords);
   amaip->aligncoords = (Uint4Ptr)MemNew((c+1)*sizeof(Uint4));
   amaip->lens = (Int4Ptr)MemNew((c+1)*sizeof(Int4));
   amaip->numseg = c;
   tmparray = (Int4Ptr)MemNew((c+1)*sizeof(Int4));
   ams = ams_master;
   for (j=0; ams && j < c; j++)
   {
      amaip->lens[j] = AlnMgrGetAlnLength(ams->sap, FALSE);
      amaip->aligncoords[j+1] = amaip->aligncoords[j] + amaip->lens[j];
      tmparray[ams->masternum] = j;
      ams = ams->next;
   }
   rsp = RowSourceNew();
   rsp->id = SeqIdDup(ams_master->sip);
   rsp->which_saps = (Uint4Ptr)MemNew(c*sizeof(Uint4));
   rsp->num_in_sap = (Uint4Ptr)MemNew(c*sizeof(Uint4));
   rsp->numsaps = 0;
   ams = ams_master;
   while (ams)
   {
      rsp->which_saps[rsp->numsaps] = ams->nsap;
      rsp->num_in_sap[rsp->numsaps] = ams->n;
      rsp->numsaps++;
      ams = ams->next;
   }
   amaip->rowsource[0] = rsp;
   amaip->numrows = 1;
   siplist = (AMsiplistPtr)MemNew(sizeof(AMsiplist));
   siplist->sip = rsp->id;
   siplist->first_row = 0;
   siplist_tmp = siplist;
   ams = ams_head;
   rsp = RowSourceNew();
   rsp->which_saps = (Uint4Ptr)MemNew(c*sizeof(Uint4));
   rsp->num_in_sap = (Uint4Ptr)MemNew(c*sizeof(Uint4));
   amaip->rowsource[amaip->numrows] = rsp;
   amaip->numrows++;
   while (ams && amaip->numrows <= max)
   {
      if (rsp->id == NULL) /* new rsp */
      {
         rsp->id = SeqIdDup(ams->sip);
         rsp->strand = ams->strand;
         rsp->which_saps[tmparray[ams->masternum]] = ams->nsap;
         rsp->num_in_sap[tmparray[ams->masternum]] = ams->n;
         rsp->numsaps = c;
         rspnum = am_get_first_rsp_for_sip(ams->sip, siplist);
         if (rspnum == -1) /* need to add to seqid list */
         {
            siplist_new = (AMsiplistPtr)MemNew(sizeof(AMsiplist));
            siplist_new->sip = ams->sip;
            siplist_new->first_row = amaip->numrows-1;
            siplist_tmp->next = siplist_new;
            siplist_tmp = siplist_new;
         }
      } else /* some fields already filled -- check for conflicts or new row */
      {
         n = SeqIdComp(rsp->id, ams->sip);
         if (n == SIC_YES && ams->strand == rsp->strand) /* could be same row -- check for conflicts */
         {
            ok = FALSE;
            if (rsp->which_saps[tmparray[ams->masternum]] == 0) /* put in same row */
            {
               done = FALSE;
               i = 0;
               while (!done && i<c)
               {
                  if (rsp->which_saps[i] != 0)
                     done = TRUE;
                  else
                     i++;
               }
               if (done)
               {
                  AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[i] - 1], rsp->num_in_sap[i], &start, &stop);
                  if (ams->strand == Seq_strand_minus)
                  {
                     if (tmparray[ams->masternum] < i)
                     {
                        if (stop >= ams->sstart)
                           ok = FALSE;
                        else
                           ok = TRUE;
                     } else
                     {
                        if (start <= ams->sstop)
                           ok = FALSE;
                        else
                           ok = TRUE;
                     }
                  } else
                  {
                     if (tmparray[ams->masternum] < i)
                     {
                        if (start <= ams->sstop)
                           ok = FALSE;
                        else
                           ok = TRUE; 
                     } else
                     {
                        if (stop >= ams->sstart)
                           ok = FALSE;
                        else
                           ok = TRUE;
                     }
                  }
               }
            }
            if (ok)
            {
               rsp->which_saps[tmparray[ams->masternum]] = ams->nsap;
               rsp->num_in_sap[tmparray[ams->masternum]] = ams->n;
               rsp->numsaps=c;
            } else
            {
               rspnum = am_get_first_rsp_for_sip(ams->sip, siplist);
               if (rspnum == -1) /* make a new row */
               {
                  rsp = RowSourceNew();
                  rsp->strand = ams->strand;
                  rsp->which_saps = (Uint4Ptr)MemNew(c*sizeof(Uint4));
                  rsp->num_in_sap = (Uint4Ptr)MemNew(c*sizeof(Uint4));
                  amaip->rowsource[amaip->numrows] = rsp;
                  amaip->numrows++;
                  rsp->id = SeqIdDup(ams->sip);
                  rsp->which_saps[tmparray[ams->masternum]] = ams->nsap;
                  rsp->num_in_sap[tmparray[ams->masternum]] = ams->n;
                  rsp->numsaps = c;
                  siplist_new = (AMsiplistPtr)MemNew(sizeof(AMsiplist));
                  siplist_new->sip = ams->sip;
                  siplist_new->first_row = amaip->numrows-1;
                  siplist_tmp->next = siplist_new;
                  siplist_tmp = siplist_new;
               } else
               {
                  done = FALSE;
                  while (rspnum < amaip->numrows && !done && SAM_OrderSeqID(ams->sip, amaip->rowsource[rspnum]->id) == 0)
                  {
                     rsp = amaip->rowsource[rspnum];
                     if (rsp->which_saps[tmparray[ams->masternum]] == 0) /* fits here */
                     {
                        done = TRUE;
                        found = FALSE;
                        i = 0;
                        while (!found && i<c)
                        {
                           if (rsp->which_saps[i] != 0)
                              found = TRUE;
                           else
                              i++;
                        }
                        if (found)
                        {
                           AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[i] - 1], rsp->num_in_sap[i], &start, &stop);
                           if (ams->strand == Seq_strand_minus)
                           {
                              if (tmparray[ams->masternum] < i)
                              {
                                 if (stop >= ams->sstart)
                                    ok = FALSE;
                                 else
                                    ok = TRUE;
                              } else
                              {
                                 if (start <= ams->sstop)
                                    ok = FALSE;
                                 else
                                    ok = TRUE;
                              }
                           } else
                           {
                              if (tmparray[ams->masternum] < i)
                              {
                                 if (start <= ams->sstop)
                                    ok = FALSE;
                                 else
                                    ok = TRUE;
                              } else
                              {
                                 if (stop >= ams->sstart)
                                    ok = FALSE;
                                 else
                                    ok = TRUE;
                              }
                           }
                        }
                        if (ok && found)
                        {
                           rsp->which_saps[tmparray[ams->masternum]] = ams->nsap;
                           rsp->num_in_sap[tmparray[ams->masternum]] = ams->n;
                           rsp->numsaps = c;
                        } else
                        {
                           rsp = RowSourceNew();
                           rsp->which_saps = (Uint4Ptr)MemNew(c*sizeof(Uint4));
                           rsp->num_in_sap = (Uint4Ptr)MemNew(c*sizeof(Uint4));
                           rsp->strand = ams->strand;
                           amaip->rowsource[amaip->numrows] = rsp;
                           amaip->numrows++;
                           rsp->id = SeqIdDup(ams->sip);
                           rsp->which_saps[tmparray[ams->masternum]] = ams->nsap;
                           rsp->num_in_sap[tmparray[ams->masternum]] = ams->n;
                           rsp->numsaps = c;
                        }
                     }
                     rspnum++;
                  }
                  if (!done) /* didn't fit */
                  {
                     rsp = RowSourceNew();
                     rsp->which_saps = (Uint4Ptr)MemNew(c*sizeof(Uint4));
                     rsp->num_in_sap = (Uint4Ptr)MemNew(c*sizeof(Uint4));
                     amaip->rowsource[amaip->numrows] = rsp;
                     amaip->numrows++;
                     rsp->id = SeqIdDup(ams->sip);
                     rsp->which_saps[tmparray[ams->masternum]] = ams->nsap;
                     rsp->num_in_sap[tmparray[ams->masternum]] = ams->n;
                     rsp->numsaps=c;
                  }
               }
            }
         } else  /* make a new row */
         {
            rsp = RowSourceNew();
            rsp->which_saps = (Uint4Ptr)MemNew(c*sizeof(Uint4));
            rsp->num_in_sap = (Uint4Ptr)MemNew(c*sizeof(Uint4));
            amaip->rowsource[amaip->numrows] = rsp;
            amaip->numrows++;
            rsp->id = SeqIdDup(ams->sip);
            rsp->strand = ams->strand;
            rsp->which_saps[tmparray[ams->masternum]] = ams->nsap;
            rsp->num_in_sap[tmparray[ams->masternum]] = ams->n;
            rsp->numsaps=c;
            siplist_new = (AMsiplistPtr)MemNew(sizeof(AMsiplist));
            siplist_new->sip = ams->sip;
            siplist_new->first_row = amaip->numrows-1;
            siplist_tmp->next = siplist_new;
            siplist_tmp = siplist_new;
         }
      }
      ams = ams->next;
   }
   siplist_tmp = siplist;
   while (siplist_tmp)
   {
      siplist_new = siplist_tmp->next;
      siplist_tmp->sip = NULL;
      siplist_tmp->next = NULL;
      MemFree(siplist_tmp);
      siplist_tmp = siplist_new;
   }
   ams = ams_master;
   while (ams)
   {
      ams_tmp = ams->next;
      ams->sap = NULL;
      ams->sip = NULL;
      ams->next = NULL;
      MemFree(ams);
      ams = ams_tmp;
   }
   ams = ams_head;
   while (ams)
   {
      ams_tmp = ams->next;
      ams->sap = NULL;
      ams->sip = NULL;
      ams->next = NULL;
      MemFree(ams);
      ams = ams_tmp;
   }
   MemFree(amsarray);
   MemFree(tmparray);
   amaip->starts = (Int4Ptr)MemNew(amaip->numseg*sizeof(Int4));
   return TRUE;
}

static Int4 am_get_first_rsp_for_sip(SeqIdPtr sip, AMsiplistPtr siplist)
{
   AMsiplistPtr  siplist_tmp;

   if (sip == NULL || siplist == NULL)
      return -1;
   siplist_tmp = siplist;
   while (siplist_tmp)
   {
      if (SeqIdComp(sip, siplist_tmp->sip) == SIC_YES)
      {
         return (siplist_tmp->first_row);
      }
      siplist_tmp = siplist_tmp->next;
   }
   return -1;
}

static AMmsmsPtr am_sort_ammsms(AMmsmsPtr ams_head, Int4 n)
{
   AMmsmsPtr        ams;
   AMmsmsPtr        ams_tmp;
   AMmsmsPtr PNTR   ams_array;
   Int4             i;

   if (ams_head == NULL || n == 0)
      return NULL;
   if (n == 1)
      return ams_head;
   ams_array = (AMmsmsPtr PNTR)MemNew((n+1)*sizeof(AMmsmsPtr));
   ams = ams_head;
   for (i=0; ams!=NULL && i<n; i++)
   {
      ams_array[i] = ams;
      ams = ams->next;
   }
   HeapSort((Pointer)ams_array, (size_t)(n), sizeof(AMmsmsPtr), AlnMgrCompareAMS);
   ams_tmp = NULL;
   for (i=0; i<n; i++)
   {
      if (ams_tmp != NULL)
      {
         ams->next = ams_array[i];
         ams = ams->next;
         ams->next = NULL;
      } else
      {
         ams_tmp = ams = ams_array[i];
         ams_tmp->next = NULL;
      }
   }
   return ams_tmp;
}

NLM_EXTERN int LIBCALLBACK AlnMgrCompareAMS(VoidPtr base, VoidPtr large_son)
{
   AMmsmsPtr  ams1;
   AMmsmsPtr  ams2;

   ams1 = *((AMmsmsPtr PNTR) base);
   ams2 = *((AMmsmsPtr PNTR) large_son);
   if (ams1 == NULL || ams2 == NULL)
      return 0;
   return (SAM_OrderSeqID(ams1->sip, ams2->sip));
}

static AMmsmsPtr am_sort_masterams(AMmsmsPtr ams_head, Int4 n)
{
   AMmsmsPtr        ams;
   AMmsmsPtr        ams_tmp;
   AMmsmsPtr PNTR   ams_array;
   Int4             i;

   if (ams_head == NULL || n == 0)
      return NULL;
   if (n == 1)
      return ams_head;
   ams_array = (AMmsmsPtr PNTR)MemNew((n+1)*sizeof(AMmsmsPtr));
   ams = ams_head;
   for (i=0; ams!=NULL && i<n; i++)
   {
      ams_array[i] = ams;
      ams = ams->next;
   }
   HeapSort((Pointer)ams_array, (size_t)(n), sizeof(AMmsmsPtr), AlnMgrCompareMasterAMS);
   ams_tmp = NULL;
   for (i=0; i<n; i++)
   {
      if (ams_tmp != NULL)
      {
         ams->next = ams_array[i];
         ams = ams->next;
         ams->next = NULL;
      } else
      {
         ams_tmp = ams = ams_array[i];
         ams_tmp->next = NULL;
      }
   }
   return ams_tmp;
}

NLM_EXTERN int LIBCALLBACK AlnMgrCompareMasterAMS(VoidPtr base, VoidPtr large_son)
{
   AMmsmsPtr  ams1;
   AMmsmsPtr  ams2;

   ams1 = *((AMmsmsPtr PNTR) base);
   ams2 = *((AMmsmsPtr PNTR) large_son);
   if (ams1 == NULL || ams2 == NULL)
      return 0;
   if (ams1->start < ams2->start)
      return -1;
   else if (ams1->start > ams2->start)
      return 1;
   else if (ams1->stop < ams2->stop)
      return -1;
   else
      return 0;
}


NLM_EXTERN void AlnMgrSetMaster(SeqAlignPtr sap, SeqIdPtr master)
{
   SAIndexPtr   saip;
   SeqAlignPtr  salp;

   if (sap->segtype != SAS_DISC || !master)
      return;
   sap->master = SeqIdDup(master);
   salp = (SeqAlignPtr)(sap->segs);
   while (salp)
   {
      if (!salp->saip)
         return;
      saip = (SAIndexPtr)(salp->saip);
      saip->master = AlnMgrGetNForSip(salp, master);
      salp = salp->next;
   }
   return;
}

NLM_EXTERN void AlnMgrMakeMasterPlus(SeqAlignPtr sap)
{
   DenseSegPtr      dsp;
   Int4             i;
   Int4             master;
   SAIndexPtr       saip;
   SeqAlignPtr      sap_tmp;

   i = AlnMgrCheckAlignForParent(sap);
   if (i==AM_CHILD)
   {
      saip = (SAIndexPtr)(sap->saip);
      if (saip->master < 0)
         return;
      else
         master = saip->master;
      dsp = (DenseSegPtr)(sap->segs);
      if (dsp->strands[saip->master-1] == Seq_strand_minus)
      {
         sap_tmp = sap;
         sap = sap->next;
         sap_tmp->next = NULL;
         sap_tmp = SeqAlignListReverseStrand(sap_tmp);
         if (!AlnMgrIndexSingleChildSeqAlign(sap_tmp))
            return;
         saip = (SAIndexPtr)(sap_tmp->saip);
         saip->master = master;
         sap_tmp->next = sap;
         sap = sap_tmp;
      }
   } else if (i==AM_PARENT)
   {
      sap_tmp = (SeqAlignPtr)(sap->segs);
      while (sap_tmp)
      {
         AlnMgrMakeMasterPlus(sap_tmp);
         sap_tmp = sap_tmp->next;
      }
   }
   return;
}

NLM_EXTERN Boolean AlnMgrForceMasterSlave(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   AMAlignDatPtr    amadp;
   Int4             n;

   if (sap == NULL || sap->master == NULL || sap->saip == NULL)
      return FALSE;
   amaip = (AMAlignIndexPtr)sap->saip;
   n = AlnMgrGetNForSip(sap, sap->master);
   if (n < 1)
      return FALSE;
   amadp = amaip->amadp[n-1];
   if (!AlnMgrMergeIntoMSMultByMaster(amaip, amaip->lens, &amaip->numseg))
      return FALSE;
   amaip->starts = (Int4Ptr)MemNew((amaip->numseg)*(amaip->numsaps)*sizeof(Int4));
   amaip->aligncoords = (Uint4Ptr)MemNew((amaip->numseg)*sizeof(Uint4));
   if (!AlnMgrFillInStarts(amadp->saps, amaip->starts, amaip->numseg, amaip->lens, amaip->numsaps, amaip->aligncoords))
      return FALSE;
   if (amaip->numseg > 1)
      amaip->numseg -= 1;
   if (!AlnMgrGetRowsForMasterSlave(sap))
      return FALSE;
   return TRUE;
}

NLM_EXTERN SeqAlignPtr AlnMgrGetSubAlign(SeqAlignPtr sap, SeqIdPtr which_master, Int4 from, Int4 to)
{
   AMAlignIndexPtr  amaip;
   AlnMsgPtr        amp;
   DenseSegPtr      dsp;
   Int4             i;
   Int4             j;
   Boolean          more;
   Uint4            n;
   Int4             numaln;
   SeqAlignPtr      salp;
   SeqAlignPtr      salp_head;
   SeqAlignPtr      salp_prev;
   SeqAlignPtr      sap_parent;
   SeqIdPtr         sip;
   SeqIdPtr         sip_curr;
   SeqIdPtr         sip_prev;
   Uint4            start;
   Uint4            stop;
   Int4Ptr          trackarray;

   i = AlnMgrCheckAlignForParent(sap);
   if (i == AM_CHILD)
   {
      salp = SeqAlignDup(sap);
      return salp;
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (amaip == NULL)
         return NULL;
      if (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_MASTERSLAVE)
      {
         salp = SeqAlignNew();
         salp->type = SAT_MASTERSLAVE;
         salp->segtype = SAS_DENSEG;
         salp->dim = amaip->numrows;
         dsp = DenseSegNew();
         dsp->dim = amaip->numrows;
         dsp->numseg = amaip->numseg;
         dsp->starts = (Int4Ptr)MemNew((amaip->numseg+1)*(amaip->numrows)*sizeof(Int4));
         dsp->lens = (Int4Ptr)MemNew((amaip->numseg+1)* sizeof(Int4));
         dsp->strands = (Uint1Ptr)MemNew((amaip->numseg+1)*(amaip->numrows)*sizeof(Uint1));
         sip_curr = amaip->ids;
         dsp->ids = SeqIdDupList(sip_curr);
         amp = AlnMsgNew();
         for (j=0; j<(amaip->numrows); j++)
         {
            if (j == amaip->master - 1)
               salp->master = SeqIdDup(sip_curr);
            sip_curr = sip_curr->next;
            amp->which_master = which_master;
            amp->from_m = from;
            amp->to_m = to;
            amp->row_num = j + 1;
            more = TRUE;
            n = 0;
            while (more)
            {
               more = AlnMgrGetNextAlnBit(sap, amp);
               if (amp->gap == 0)
               {
                  dsp->starts[n*(dsp->dim) + j] = amp->from_b;
               } else
               {
                  dsp->starts[n*(dsp->dim) + j] = -1;
               }
               if (j == 0)
                  dsp->lens[n] = amp->to_b - amp->from_b + 1;
               dsp->strands[n*(dsp->dim) + j] = amp->strand;
               n++;
            }
            amp = AlnMsgReNew(amp);
         }
         salp->segs = (Pointer)dsp;
         return salp;
      } else if (sap->type == SAT_PARTIAL || (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_SEGMENTED_MASTERSLAVE))
      {
         amp = AlnMsgNew();
         amp->which_master = which_master;
         amp->from_m = from;
         amp->to_m = to;
         amp->row_num = 1;
         trackarray = (Int4Ptr)MemNew((amaip->numseg+1)*sizeof(Int4));
         numaln = 0;
         while (more = AlnMgrGetNextAlnBit(sap, amp))
         {
            if (amp->send_space)
            {
               numaln++;
               amp->send_space = FALSE;
            } else
               trackarray[numaln]++;
         }
         numaln++;
         salp_head = NULL;
         sip_curr = NULL;
         for (j=0; j<amaip->numrows; j++)
         {
            sip = AlnMgrGetNthSeqIdPtr(sap, j+1);
            if (sip_curr != NULL)
            {
               sip_prev->next = sip;
               sip_prev = sip;
            } else
               sip_curr = sip_prev = sip;
         }
         for (j=0; j<numaln; j++)
         {
            salp = SeqAlignNew();
            if (salp_head != NULL)
            {
               salp_prev->next = salp;
               salp_prev = salp;
            } else
               salp_prev = salp_head = salp;
            salp->type = SAT_PARTIAL;
            salp->segtype = SAS_DENSEG;
            salp->dim = amaip->numrows;
            dsp = DenseSegNew();
            dsp->dim = amaip->numrows;
            dsp->numseg = trackarray[j]+1;
            dsp->starts = (Int4Ptr)MemNew((dsp->dim)*(trackarray[j]+1)*sizeof(Int4));
            dsp->lens = (Int4Ptr)MemNew((dsp->dim)*(trackarray[j]+1)*sizeof(Int4));
            dsp->strands = (Uint1Ptr)MemNew((dsp->dim)*(trackarray[j]+1)*sizeof(Uint1));
            dsp->ids = SeqIdDupList(sip_curr);
            salp->segs = (Pointer)dsp;
         }
         amp = AlnMsgReNew(amp);
         for (j=0; j<(amaip->numrows); j++)
         {
            salp = salp_head;
            dsp = (Pointer)(salp->segs);
            if (j == amaip->master - 1)
               salp->master = SeqIdDup(sip_curr);
            sip_curr = sip_curr->next;
            amp->which_master = which_master;
            amp->from_m = from;
            amp->to_m = to;
            amp->row_num = j + 1;
            more = TRUE;
            n = 0;
            while (more)
            {
               more = AlnMgrGetNextAlnBit(sap, amp);
               if (amp->gap == 0)
               {
                  dsp->starts[n*(dsp->dim) + j] = amp->from_b;
               } else
               {
                  dsp->starts[n*(dsp->dim) + j] = -1;
               }
               if (j == 0)
                  dsp->lens[n] = amp->to_b - amp->from_b + 1;
               dsp->strands[n*(dsp->dim) + j] = amp->strand;
               n++;
               if (amp->send_space == TRUE)
               {
                  salp = salp->next;
                  dsp = (DenseSegPtr)(salp->segs);
                  amp->send_space = FALSE;
                  n=0;
               }
            }
            amp = AlnMsgReNew(amp);
         }
         /*if (salp_head->next != NULL) 
         {
            sap_parent = SeqAlignNew();
            sap_parent->type = SAT_MASTERSLAVE;
            sap_parent->segtype = SAS_DISC;
            sap_parent->dim = amaip->numrows;
            sap_parent->segs = (Pointer)salp_head;
            MemFree(trackarray);
            MemFree(amp);
            return sap_parent;
         } else
         {*/
            MemFree(trackarray);
            MemFree(amp);
            return salp_head;
        /* }*/
      } else if (sap->type == SAT_DIAGS)
      {
         salp = SeqAlignDup(sap);
         return salp;
      }
   }
   return NULL;
}


/********************************************************************************
*
*   viewer and editor management functions
*
********************************************************************************/

NLM_EXTERN SeqAlignPtr AlnMgrCopyIndexedParentSeqAlign(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   AMAlignIndexPtr  amaip_new;
   Boolean          found;
   Int4             i;
   Int4Ptr          orderarray;
   Int4             r;
   SeqAlignPtr      sap_new;
   SeqAlignPtr      sap_tmp;
   SeqAlignPtr      seg_head;
   SeqAlignPtr      seg_new;
   SeqAlignPtr      seg_prev;
   SeqAlignPtr      seg_tmp;

   if (sap->saip == NULL)
      return NULL;
   if (sap->saip->indextype != INDEX_PARENT)
      return NULL;
   amaip = (AMAlignIndexPtr)sap->saip;
   amaip_new = AMAlignIndexNew();
   sap_new = SeqAlignDup(sap);
   sap_new->saip = (SeqAlignIndexPtr)amaip_new;
   amaip_new->indextype = amaip->indextype;
   amaip_new->freefunc = amaip->freefunc;
   amaip_new->mstype = amaip->mstype;
   amaip_new->aligncoords = (Uint4Ptr)MemNew((amaip->numseg+1)*sizeof(Uint4));
   amaip_new->numseg = amaip->numseg;
   amaip_new->lens = (Int4Ptr)MemNew((amaip->numseg+1)*sizeof(Int4));
   for (i=0; i<amaip->numseg; i++)
   {
      amaip_new->aligncoords[i] = amaip->aligncoords[i];
      amaip_new->lens[i] = amaip->lens[i];
   }
   amaip_new->starts = (Int4Ptr)MemNew((amaip->numseg)*(amaip->numsaps+1)*sizeof(Int4));
   for(i=0; i<(amaip->numseg)*(amaip->numsaps); i++)
   {
      amaip_new->starts[i] = amaip->starts[i];
   }
   amaip_new->alnsaps = amaip->alnsaps;
   amaip_new->numsaps = amaip->numsaps;
   amaip_new->ids = SeqIdDupList(amaip->ids);
   amaip_new->numbsqs = amaip->numbsqs;
   amaip_new->rowsource = (RowSourcePtr PNTR)MemNew((amaip->numrows+1)*sizeof(RowSourcePtr));
   for (i=0; i<amaip->numrows; i++)
   {
      amaip_new->rowsource[i] = AlnMgrCopyRowSource(amaip->rowsource[i]);
   }
   amaip_new->numrows = amaip->numrows;
   amaip_new->master = amaip->master;
   seg_head = NULL;
   sap_tmp = (SeqAlignPtr)sap->segs;
   while (sap_tmp != NULL)
   {
      seg_new = SeqAlignDup(sap_tmp);
      if (seg_head != NULL)
      {
         seg_prev->next = seg_new;
         seg_prev = seg_new;
      } else
         seg_head = seg_prev = seg_new;
      sap_tmp = sap_tmp->next;
   }
   sap_new->segs = (Pointer)seg_head;
   i = 0;
   orderarray = (Int4Ptr)MemNew((amaip->numsaps)*sizeof(Int4));
   seg_new = seg_head;
   sap_tmp = (SeqAlignPtr)sap->segs;
   while (sap_tmp != NULL && seg_new != NULL)
   {
      seg_new->saip = AlnMgrCopyIndexesForChildSeqAlign(sap_tmp);
      found = FALSE;
      r = 0;
      while (!found && r < amaip->numsaps)
      {
         if (sap_tmp == amaip->saps[r])
         {
            orderarray[i] = r;
            found = TRUE;
         }
         r++;
      }
      i++;
      seg_new = seg_new->next;
      sap_tmp = sap_tmp->next;
   }
   amaip_new->saps = (SeqAlignPtr PNTR)MemNew((amaip->numsaps+1)*sizeof(SeqAlignPtr));
   seg_tmp = (SeqAlignPtr)sap_new->segs;
   i = 0;
   while (seg_tmp)
   {
      amaip_new->saps[orderarray[i]] = seg_tmp;
      i++;
      seg_tmp = seg_tmp->next;
   }
   sap_tmp = (SeqAlignPtr)sap->segs;
   amaip_new->amadp = (AMAlignDatPtr PNTR)MemNew((amaip->numbsqs+1)*sizeof(AMAlignDatPtr));
   seg_head = (SeqAlignPtr)sap_new->segs;
   for (i=0; i<amaip->numbsqs; i++)
   {
      amaip_new->amadp[i] = AlnMgrCopyamadp(amaip->amadp[i], sap_tmp, seg_head);
   }
   MemFree(orderarray);
   return sap_new;
}

NLM_EXTERN RowSourcePtr AlnMgrCopyRowSource(RowSourcePtr rsp)
{
   Int4          i;
   RowSourcePtr  rsp_new;

   rsp_new = RowSourceNew();
   rsp_new->id = SeqIdDup(rsp->id);
   rsp_new->which_saps = (Uint4Ptr)MemNew((rsp->numsaps+1)*sizeof(Uint4));
   rsp_new->num_in_sap = (Uint4Ptr)MemNew((rsp->numsaps+1)*sizeof(Uint4));
   for (i=0; i<rsp->numsaps; i++)
   {
      rsp_new->which_saps[i] = rsp->which_saps[i];
      rsp_new->num_in_sap[i] = rsp->num_in_sap[i];
   }
   rsp_new->numsaps = rsp->numsaps;
   return rsp_new;
}

NLM_EXTERN AMAlignDatPtr AlnMgrCopyamadp(AMAlignDatPtr amadp, SeqAlignPtr sap_tmp, SeqAlignPtr seg_head)
{
   AMAlignDatPtr  amadp_new;
   Boolean        found;
   Int4           i;
   Int4           j;
   Int4Ptr        orderarray;
   SeqAlignPtr    sap_old;
   SeqAlignPtr    sap_new;

   if (sap_tmp == NULL || amadp == NULL || seg_head == NULL)
      return NULL;
   amadp_new = AMAlignDatNew();
   amadp_new->sip = SeqIdDup(amadp->sip);
   amadp_new->numsaps = amadp->numsaps;
   amadp_new->saps = (SeqAlignPtr PNTR)MemNew((amadp->numsaps+1)*sizeof(SeqAlignPtr));
   orderarray = (Int4Ptr)MemNew((amadp->numsaps+1)*sizeof(Int4));
   sap_old = sap_tmp;
   j = 0;
   while (sap_old)
   {
      i=0;
      found = FALSE;
      while (!found && i<amadp->numsaps)
      {
         if (sap_old == amadp->saps[i])
         {
            orderarray[i] = j;
            found = TRUE;
         }
         i++;
      }
      sap_old = sap_old->next;
      j++;
   }
   for (i=0; i<amadp->numsaps; i++)
   {
      sap_new = seg_head;
      j=0;
      while (j<orderarray[i])
      {
         sap_new = sap_new->next;
         j++;
      }
      amadp_new->saps[i] = sap_new;
   }
   amadp_new->segments = (Uint2Ptr)MemNew((amadp->numseg+1)*sizeof(Uint2));
   for (i=0; i<amadp->numseg; i++)
   {
      amadp_new->segments[i] = amadp->segments[i];
   }
   amadp_new->numseg = amadp->numseg;
   MemFree(orderarray);
   return amadp_new;
}

NLM_EXTERN SeqAlignIndexPtr AlnMgrCopyIndexesForChildSeqAlign(SeqAlignPtr sap)
{
   DenseSegPtr  dsp;
   Int4         i;
   SAIndexPtr   saip;
   SAIndexPtr   saip_new;

   if (sap == NULL || sap->saip == NULL)
      return NULL;
   dsp = (DenseSegPtr)sap->segs;
   saip = (SAIndexPtr)sap->saip;
   saip_new = SAIndexNew();
   saip_new->indextype = saip->indextype;
   saip_new->freefunc = saip->freefunc;
   saip_new->master = saip->master;
   saip_new->aligncoords = (Uint4Ptr)MemNew((dsp->numseg + 1)*sizeof(Uint4));
   for (i=0; i<dsp->numseg; i++)
   {
      saip_new->aligncoords[i] = saip->aligncoords[i];
   }
   saip_new->ssdp = (SASeqDatPtr PNTR)MemNew((dsp->dim+1)*sizeof(SASeqDatPtr));
   for (i=0; i<dsp->dim; i++)
   {
      saip_new->ssdp[i] = AlnMgrCopySASeqDat(saip->ssdp[i]);
   }
   return (SeqAlignIndexPtr)saip_new;
}

NLM_EXTERN SASeqDatPtr  AlnMgrCopySASeqDat(SASeqDatPtr ssdp)
{
   Int4         i;
   SASeqDatPtr  ssdp_new;

   if (ssdp == NULL)
      return NULL;
   ssdp_new = SASeqDatNew();
   ssdp_new->numsect = ssdp->numsect;
   ssdp_new->numunsect = ssdp->numunsect;
   ssdp_new->sect = (Uint2Ptr)MemNew((ssdp->numsect+1)*sizeof(Uint2));
   for (i=0; i<ssdp->numsect; i++)
   {
      ssdp_new->sect[i] = ssdp->sect[i];
   }
   ssdp_new->unsect = (Uint2Ptr)MemNew((ssdp->numunsect+1)*sizeof(Uint2));
   for (i=0; i<ssdp->numunsect; i++)
   {
      ssdp_new->unsect[i] = ssdp->unsect[i];
   }
   return ssdp_new;
}

NLM_EXTERN SeqAlignPtr AlnMgrCopyAndIndexSingleAlignment(SeqAlignPtr sap)
{
   SeqAlignPtr  sap_new;

   if (sap == NULL)
      return NULL;
   sap_new = SeqAlignDup(sap);
   sap_new->type = SAT_MASTERSLAVE;
   sap_new->saip = AlnMgrCopyIndexesForChildSeqAlign(sap);
   return sap_new;
}

NLM_EXTERN Boolean AlnMgrCopyIndexedParentIntoSap(SeqAlignPtr sap, SeqAlignPtr target)
{
   AMAlignIndexPtr  amaip;
   AMAlignIndexPtr  amaip_new;
   DenseSegPtr      dsp_tmp;
   Boolean          found;
   Int4             i;
   Int4Ptr          orderarray;
   Int4             r;
   SeqAlignPtr      sap_tmp;
   SeqAlignPtr      seg_head;
   SeqAlignPtr      seg_new;
   SeqAlignPtr      seg_prev;
   SeqAlignPtr      seg_tmp;

   if (sap->saip == NULL || target == NULL)
      return FALSE;
   if (sap->saip->indextype != INDEX_PARENT)
      return FALSE;
   AMAlignIndexFree((Pointer)target->saip);
   target->saip = NULL;
   amaip = (AMAlignIndexPtr)sap->saip;
   amaip_new = AMAlignIndexNew();
   target->type = sap->type;
   target->segtype = sap->segtype;
   target->dim = sap->dim;
   target->score = ScoreSetFree(target->score);
   target->score = ScoreDup(sap->score);
   target->master = SeqIdFree(target->master);
   target->master = SeqIdDup(sap->master);
   dsp_tmp = (DenseSegPtr)target->segs;
   target->segs = DenseSegFree(dsp_tmp);
   target->segs = NULL;
   target->saip = (SeqAlignIndexPtr)amaip_new;
   amaip_new->indextype = amaip->indextype;
   amaip_new->freefunc = amaip->freefunc;
   amaip_new->mstype = amaip->mstype;
   amaip_new->aligncoords = (Uint4Ptr)MemNew((amaip->numseg+1)*sizeof(Uint4));
   amaip_new->numseg = amaip->numseg;
   amaip_new->lens = (Int4Ptr)MemNew((amaip->numseg+1)*sizeof(Int4));
   for (i=0; i<amaip->numseg; i++)
   {
      amaip_new->aligncoords[i] = amaip->aligncoords[i];
      amaip_new->lens[i] = amaip->lens[i];
   }
   amaip_new->starts = (Int4Ptr)MemNew((amaip->numseg)*(amaip->numsaps+1)*sizeof(Int4));
   for(i=0; i<(amaip->numseg)*(amaip->numsaps); i++)
   {
      amaip_new->starts[i] = amaip->starts[i];
   }
   amaip_new->alnsaps = amaip->alnsaps;
   amaip_new->numsaps = amaip->numsaps;
   amaip_new->ids = SeqIdDupList(amaip->ids);
   amaip_new->numbsqs = amaip->numbsqs;
   amaip_new->rowsource = (RowSourcePtr PNTR)MemNew((amaip->numrows+1)*sizeof(RowSourcePtr));
   for (i=0; i<amaip->numrows; i++)
   {
      amaip_new->rowsource[i] = AlnMgrCopyRowSource(amaip->rowsource[i]);
   }
   amaip_new->numrows = amaip->numrows;
   amaip_new->master = amaip->master;
   sap_tmp = (SeqAlignPtr)sap->segs;
   seg_head = NULL;
   while (sap_tmp)
   {
      seg_new = SeqAlignDup(sap_tmp);
      sap_tmp = sap_tmp->next;
      if (seg_head != NULL)
      {
         seg_prev->next = seg_new;
         seg_prev = seg_new;
      } else
         seg_head = seg_prev = seg_new;
   }
   sap_tmp = (SeqAlignPtr)sap->segs;
   i = 0;
   orderarray = (Int4Ptr)MemNew((amaip->numsaps)*sizeof(Int4));
   target->segs = (SeqAlignPtr)seg_head;
   seg_new = seg_head;
   seg_head = NULL;
   while (sap_tmp && seg_new)
   {
      seg_new->saip = AlnMgrCopyIndexesForChildSeqAlign(sap_tmp);
      found = FALSE;
      r = 0;
      while (!found && r < amaip->numsaps)
      {
         if (sap_tmp == amaip->saps[r])
         {
            orderarray[i] = r;
            found = TRUE;
         }
         r++;
      }
      i++;
      seg_new = seg_new->next;
      sap_tmp = sap_tmp->next;
   }
   amaip_new->saps = (SeqAlignPtr PNTR)MemNew((amaip->numsaps+1)*sizeof(SeqAlignPtr));
   seg_tmp = (SeqAlignPtr)target->segs;
   i = 0;
   while (seg_tmp)
   {
      amaip_new->saps[orderarray[i]] = seg_tmp;
      i++;
      seg_tmp = seg_tmp->next;
   }
   sap_tmp = (SeqAlignPtr)sap->segs;
   amaip_new->amadp = (AMAlignDatPtr PNTR)MemNew((amaip->numbsqs+1)*sizeof(AMAlignDatPtr));
   seg_head = (SeqAlignPtr)target->segs;
   for (i=0; i<amaip->numbsqs; i++)
   {
      amaip_new->amadp[i] = AlnMgrCopyamadp(amaip->amadp[i], sap_tmp, seg_head);
   }
   MemFree(orderarray);
   return TRUE;
}

NLM_EXTERN Boolean AlnMgrDeleteChildByPointer(SeqAlignPtr parent, SeqAlignPtr child)
{
   Boolean      found;
   Int4         i;
   SeqAlignPtr  salp;
   SeqAlignPtr  salp_head;
   SeqAlignPtr  salp_prev;

   if (parent == NULL || child == NULL)
      return FALSE;
   i = AlnMgrCheckAlignForParent(parent);
   if (i != INDEX_PARENT)
      return FALSE;
   salp_head = salp_prev = NULL;
   salp = (SeqAlignPtr)(parent->segs);
   found = FALSE;
   while (salp && !found)
   {
      if (salp == child)
         found = TRUE;
      else
      {
         if (salp_head)
            salp_prev = salp;
         else
            salp_head = salp_prev = salp;
         salp = salp->next;
      }
   }
   if (!found)
      return FALSE;
   if (salp_head != NULL)
   {
      salp_prev->next = salp->next;
      salp->next = NULL;
      SeqAlignFree(salp);
   } else
   {
      salp_head = salp->next;
      salp->next = NULL;
      SeqAlignFree(salp);
   }
   parent->segs = (Pointer)salp_head;
   return (AlnMgrReIndexSeqAlign(parent));
}

NLM_EXTERN void AlnMgrDeleteRow(SeqAlignPtr sap, Int4 row)
{
   AMAlignIndexPtr  amaip;
   Int4             i;
   Int4             j;
   RowSourcePtr     rsp;

   if (sap == NULL)
      return;
   i = AlnMgrCheckAlignForParent(sap);
   if (i == INDEX_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      if (amaip == NULL)
         return;
      if (row > amaip->numrows)
         return;
      rsp = amaip->rowsource[row-1];
      amaip->rowsource[row-1] = NULL;
      for (j=row-1; j<amaip->numrows-1; j++)
      {
         amaip->rowsource[j] = amaip->rowsource[j+1];
      }
      amaip->numrows--;
   } else
      return;
   return;
}

/*******************************************************************************

  Function : AlnMgrIsSAPDiscAli()
  
  Purpose : check if a SeqAlign is discontinuous
  
  Parameters : SeqAlignPtr
  
  Return value : TRUE if discontinous, FALSE otherwise

*******************************************************************************/
NLM_EXTERN  Boolean AlnMgrIsSAPDiscAli(SeqAlignPtr sap)
{
AMAlignIndexPtr  amaip;
Boolean          bRet=FALSE;

	if (!sap || !sap->saip) return(bRet);

	if (sap->saip->indextype == INDEX_PARENT){
		amaip = (AMAlignIndexPtr)sap->saip;
		if (sap->type == SAT_PARTIAL || (sap->type == SAT_MASTERSLAVE && 
			amaip->mstype == AM_SEGMENTED_MASTERSLAVE)){
			bRet=TRUE;
		}
	}
	
	return(bRet);
}

