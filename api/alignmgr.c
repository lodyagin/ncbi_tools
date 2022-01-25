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
* $Revision: 6.27 $
*
* File Description: SeqAlign indexing and messaging functions
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: alignmgr.c,v $
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

/*******************************************************************
*
*  all the memory allocation/deallocation functions
*
*******************************************************************/

SeqAlignIndexPtr SeqAlignIndexNew(void)
{
   return (SeqAlignIndexPtr)(MemNew(sizeof(SeqAlignIndex)));
}

SAIndexPtr SAIndexNew(void)
{
   SAIndexPtr  saip;

   saip = (SAIndexPtr)MemNew(sizeof(SAIndex));
   saip->master = -1;
   saip->freefunc = (SeqAlignIndexFreeFunc)(SAIndexFree);
   return saip;
}

Boolean SAIndexFree(VoidPtr index)
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

SASeqDatPtr SASeqDatNew(void)
{
   return (SASeqDatPtr)(MemNew(sizeof(SASeqDat)));
}

RowSourcePtr RowSourceNew(void)
{
   return (RowSourcePtr)(MemNew(sizeof(RowSource)));
}

RowSourcePtr RowSourceFree(RowSourcePtr rsp)
{
   if (rsp == NULL)
      return NULL;
   rsp->id = SeqIdSetFree(rsp->id);
   MemFree(rsp->which_saps);
   MemFree(rsp->num_in_sap);
   return NULL;
}

AMAlignIndexPtr AMAlignIndexNew(void)
{
   AMAlignIndexPtr  amaip;

   amaip = (AMAlignIndexPtr)MemNew(sizeof(AMAlignIndex));
   amaip->freefunc = (SeqAlignIndexFreeFunc)(AMAlignIndexFree);
   amaip->master = -2;
   return amaip;
}

Boolean AMAlignIndexFree(VoidPtr index)
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
   MemFree(amaip->starts);
   for (i=0; i<(amaip->numrows); i++)
   {
      amaip->rowsource[i] = RowSourceFree(amaip->rowsource[i]);
   }
   MemFree(amaip->rowsource);
   retval = TRUE;
   return retval;
}

AMAlignDatPtr AMAlignDatNew(void)
{
   return (AMAlignDatPtr)(MemNew(sizeof(AMAlignDat)));
}

AMAlignDatPtr AMAlignDatFree(AMAlignDatPtr amadp)
{
   MemFree(amadp->saps);
   MemFree(amadp->segments);
   return NULL;
}

AlnMsgPtr AlnMsgNew(void)
{
   AlnMsgPtr  amp;

   amp = (AlnMsgPtr)MemNew(sizeof(AlnMsg));
   amp->send_space = FALSE;
   amp->row_num = -1;
   amp->prev = -2;
   amp->prev_sap = -2;
   return amp;
}

Boolean AlnMgrIndexSeqAlign(SeqAlignPtr sap)
{
   Boolean      retval;
   SeqAlignPtr  salp;
   SeqAlignPtr  salp_head;
   SeqAlignPtr  salp_prev;
   SeqAlignPtr  sap_tmp;

   retval = FALSE;
   if (!sap)
      return retval;
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
Boolean AlnMgrAnythingToSeg (SeqAlignPtr sap)
{
   DenseDiagPtr  ddp;
   DenseDiagPtr  ddp_prev;
   DenseSegPtr   dsp;
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
   sap_new->type = sap->type;
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
            dsp->ids = SeqIdDupList(ddp->id);
            sap_new->segs = (Pointer)dsp;
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
Boolean AlnMgrIndexLinkedSegs (SeqAlignPtr sap)
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

Boolean AlnMgrIndexParentSA(SeqAlignPtr sap)
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
   sap->saip = (Pointer)amaip;
   retval = TRUE;
   return retval;
}

SeqIdPtr AlnMgrPropagateUpSeqIdPtrs(SeqAlignPtr sap, Int4Ptr num)
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
            if (SeqIdComp(sip_tmp, sip_tmp2))
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

SeqIdPtr AlnMgrPropagateSeqIdsBySapList(AMAlignIndexPtr amaip)
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

void am_print_seqalign_indexes(SeqAlignPtr sap)
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
         for (i=0; i<(amaip->numseg*amaip->numsaps); i++)
         {
            printf("%d ", amaip->starts[i]);
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

Int4 AlnMgrCheckAlignForParent(SeqAlignPtr sap)
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
SeqAlignPtr PNTR AlnMgrSortSeqAligns (SeqAlignPtr sap, int (LIBCALLBACK *compar)(VoidPtr, VoidPtr, VoidPtr), VoidPtr userdata, Int4Ptr numsap)
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
int LIBCALLBACK AlnMgrCompareIncreasingBySeqIdPtr (VoidPtr base, VoidPtr large_son, VoidPtr userdata)
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
      if (SeqIdComp(sip_tmp, sip))
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
      if (SeqIdComp(sip_tmp, sip))
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
int LIBCALLBACK AlnMgrFindFirst(VoidPtr base, VoidPtr large_son, VoidPtr userdata)
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

int LIBCALLBACK AlnMgrCompareTips(VoidPtr base, VoidPtr large_son)
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
*  to an allocated Int4 set to 0 initially.
*  AlnMgrGetNextLengthBit calls AlnMgrGetMaxUnalignedLength to get
*  the lengths of the unaligned regions.  This function makes a lot of
*  assumptions about the relationship between sap1 and sap2; call it
*  with much hesitation.
*
************************************************************************/
Boolean AlnMgrGetNextLengthBit(SeqAlignPtr sap, Int4Ptr length, Int4Ptr r)
{
   AMAlignIndexPtr  amaip;
   Int4             i;
   Boolean          retval;

   retval = FALSE;
   if (sap == NULL || length == NULL || r == NULL)
      return retval;
   i = AlnMgrCheckAlignForParent(sap);
   if (i == AM_CHILD)
   {
      *length = AlnMgrGetAlnLength(sap, FALSE);
      return FALSE;
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      if (sap->type == SAT_PARTIAL)
      {
         if (*r < 0)
         {
            *length = -(AlnMgrGetMaxUnalignedLength(amaip->saps[*r-1], amaip->saps[*r]));
            *r = -(*r);
         } else
         {
            *length = AlnMgrGetAlnLength(amaip->saps[*r], FALSE);
            if (*r < amaip->numsaps-1)
            {
               *r = -((*r)+1);
               return TRUE;
            } else
               return FALSE;
         }
      } else if (sap->type == SAT_MASTERSLAVE)
      {
         *length = amaip->aligncoords[amaip->numseg-1] + amaip->lens[amaip->numseg-1];
         return retval;
      }
   }
   return retval;
}

Int4 AlnMgrGetMaxUnalignedLength(SeqAlignPtr sap1, SeqAlignPtr sap2)
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
Boolean AlnMgrGetNextAlnBit (SeqAlignPtr sap, AlnMsgPtr amp)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             endoffset;
   Boolean          found;
   Int4             i;
   Int4             len;
   Boolean          more;
   Uint4            offset;
   Boolean          retval;
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
         len = AlnMgrGetAlnLength(sap, FALSE);
         if (amp->to_m == 0 || amp->to_m > len)
            amp->to_m = len;
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
         }
         else
         {
            start_m = binary_search_on_uint4_list(saip->aligncoords, amp->from_m, dsp->numseg);
            amp->real_from = amp->from_m;
         }
         stop_m = binary_search_on_uint4_list(saip->aligncoords, amp->to_m, dsp->numseg);
         ssdp = saip->ssdp[amp->row_num];
         offset = amp->real_from - saip->aligncoords[start_m];
         start_b = binary_search_on_uint2_list(ssdp->sect, start_m, ssdp->numsect);
         if (dsp->strands[start_b*(dsp->dim)+amp->row_num] == Seq_strand_minus)
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
                  amp->from_b = dsp->starts[start_b*(dsp->dim)+amp->row_num] + offset;
                  amp->to_b = dsp->starts[start_b*(dsp->dim)+amp->row_num] + dsp->lens[start_b] - 1;
               } else
               {
                  amp->from_b = dsp->starts[start_b*(dsp->dim)+amp->row_num];
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
            endoffset = amp->to_m - saip->aligncoords[start_m];
            if (start_b >= 0)
            {
               if (amp->strand != Seq_strand_minus)
               {
                  amp->from_b = dsp->starts[start_b*(dsp->dim)+amp->row_num] + offset;
                  amp->to_b = dsp->starts[start_b*(dsp->dim)+amp->row_num] + endoffset -1;
               } else
               {
                  amp->from_b = dsp->starts[start_b*(dsp->dim)+amp->row_num] + dsp->lens[start_b]  - endoffset;
                  amp->to_b = amp->from_b + endoffset - offset -1;
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
         }
      }
   } else if (i == AM_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (!amaip->saps)
         return retval;
      if (!amp->which_bsq && amp->row_num==-1)
         return retval;
      if (sap->type == SAT_PARTIAL && amp->which_master == NULL)
      {
         len = AlnMgrGetAlnLength(sap, FALSE);
         if (amp->to_m == 0 ||amp->to_m > len)
            amp->to_m = len;
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
         if (amp->len_left > amaip->lens[start_m])
         {
            endoffset = amaip->lens[start_m];
         } else
         {
            endoffset = amp->len_left;
         }
         stop_tmp = amp->to_m;
         start_tmp = amp->from_m;
         if ((stop_m - start_m) == 0)
         {
            amp->from_m = offset + amaip->starts[start_m];
            amp->to_m = amp->from_m + endoffset;
            AlnMgrGetNextAlnBit((amaip->saps[start_m]), amp);
            amp->len_left = amp->len_left - (amp->to_m - amp->from_m + 1)+1;
            amp->real_from = amp->to_m + 1;
            amp->to_m = stop_tmp;
            amp->from_m = start_tmp;
            amp->prev_sap = -2;
         } else
         {
            retval = TRUE;
            amp->from_m = offset + amaip->starts[start_m];
            amp->to_m = amp->from_m + endoffset;
            more = AlnMgrGetNextAlnBit((amaip->saps[start_m]), amp);
            amp->len_left = amp->len_left - (amp->to_m - amp->from_m + 1)+1;
            amp->to_m = stop_tmp;
            amp->real_from = amp->to_m - amp->len_left + 1;
            amp->from_m = start_tmp;
            if (more == FALSE)
            {
               amp->prev_sap += 1;
               amp->send_space = TRUE;
            }
            if (amp->len_left == 0)
            {
               retval = FALSE;
               amp->prev_sap = -2;
            }
         }
      } else if (sap->type == SAT_MASTERSLAVE && amp->which_master == NULL)
      {
         len = AlnMgrGetAlnLength(sap, FALSE);
         if (amp->to_m < 0 ||amp->to_m > len)
            amp->to_m = len;
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
            amp->strand = AlnMgrGetStrand(amaip->saps[i], amp->which_bsq);
            offset = amp->real_from - amaip->aligncoords[start_m];
            endoffset = amaip->lens[start_m] - offset - (amp->to_m - amp->real_from + 1);
            if (endoffset <= 0 && (start_m + 1) < amaip->numseg)
               retval = TRUE;
            else
            {
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
               retval = FALSE;
               amp->row_num = -1;
               amp->prev = -2;
            }
         }
      } else if (sap->type == SAT_MASTERSLAVE && amp->which_master)
      {
      } else if (sap->type == SAT_DIAGS && amp->which_master)
      {
      }
   }
   return retval;
}

Uint4 binary_search_on_uint4_list(Uint4Ptr list, Uint4 pos, Uint4 listlen)
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

Int4 binary_search_on_uint2_list(Uint2Ptr list, Uint2 ele, Uint2 listlen)
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

Int4 binary_search_by_chunk(Int4Ptr list, Int4 ele, Int4 listlen, Int4 chunksize, Int4 offset)
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

Int4 binary_search_segment_array(SASeqDatPtr ssdp, Int4 pos, Int4 numseq, Int4 offset, DenseSegPtr dsp)
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
Int4 AlnMgrGetAlnLength(SeqAlignPtr sap, Boolean fill_in)
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
         if (sap->type == SAT_MASTERSLAVE)
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

Int4 AlnMgrGetNumSeqs(SeqAlignPtr sap)
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

SeqIdPtr AlnMgrGetNthSeqIdPtr(SeqAlignPtr sap, Int4 n)
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
      sip = amaip->ids;
      count = 0;
      while (sip)
      {
         count++;
         if (count == n)
            return (SeqIdDup(sip));
         sip = sip->next;
      }
   }
   return NULL;
}

void AlnMgrGetNthSeqRangeInSA(SeqAlignPtr sap, Int4 n, Int4Ptr start, Int4Ptr stop)
{
   AMAlignDatPtr    amadp;
   AMAlignIndexPtr  amaip;
   Uint2            beg;
   Int4             bsq;
   DenseSegPtr      dsp;
   Uint2            end;
   Int4             i;
   Int4             j;
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
      amadp = amaip->amadp[bsq];
      sip = amaip->ids;
      for (j = 0; j<bsq; j++)
      {
         sip = sip->next;
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
   }
   return;
}

Uint1 AlnMgrGetStrand(SeqAlignPtr sap, SeqIdPtr sip)
{
   Int4  i;

   i = AlnMgrGetNForSip(sap, sip);
   return (AlnMgrGetNthStrand(sap, i));
}

Uint1 AlnMgrGetNthStrand(SeqAlignPtr sap, Int4 n)
{
   DenseSegPtr  dsp;

   if (!sap)
      return 0;
   if (sap->segtype != SAS_DENSEG)
      return 0;
   dsp = (DenseSegPtr)sap->segs;
   if (!dsp)
      return 0;
   if (n==0)
      return 0;
   return (dsp->strands[n-1]);
}

Int4 AlnMgrGetNForSip(SeqAlignPtr sap, SeqIdPtr sip)
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
         if (SeqIdComp(sip_tmp, sip))
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
         if (SeqIdComp(sip_tmp, sip))
            return n;
         sip_tmp = sip_tmp->next;
      }
   }
   return -1;
}

Int4 AlnMgrGetSapForSip(AMAlignIndexPtr amaip, SeqIdPtr sip, Int4 which)
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

Int4 AlnMgrMapToBsqCoords(SeqAlignPtr sap, Uint4 pos, SeqIdPtr sip, SeqIdPtr master)
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


Int4 AlnMgrMapRowCoords(SeqAlignPtr sap, Uint4 pos, Int4 row, SeqIdPtr master)
{
   AMAlignIndexPtr  amaip;
   DenseSegPtr      dsp;
   Int4             offset;
   SAIndexPtr       saip;
   Int4             start;

   if (sap == NULL || pos < 0 || row < 0)
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


/***********************************************************************
*
*  AlnMgrMakeFakeMultiple calls AlnMgrCheckOverlapping to decide whether
*  an alignment is linear.  Then, if possible, it calls AlnMgrMakeAlignCoords
*  to create alignment coordinates across all children contained in the
*  parent.  (MULT)
*
***********************************************************************/
Boolean AlnMgrMakeFakeMultiple(SeqAlignPtr sap)
{
   AMAlignDatPtr    amadp;
   AMAlignIndexPtr  amaip;
   Int4             i;
   Int4             j;
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
         amaip->numseg = amaip->alnsaps;
         for (j=0; j<(amaip->alnsaps); j++)
         {
            amaip->lens[j] = AlnMgrGetAlnLength(amaip->saps[j], FALSE);
            amaip->starts[j] = 0;
         }
         AlnMgrMakeAlignCoords(sap);
         if (!AlnMgrGetRowsForPartial(sap))
            return retval;
         retval = TRUE;
      } else /*should add function to check for pairwise multiple vs. diags*/
      {
         amaip = (AMAlignIndexPtr)sap->saip;
         if (amaip->saps)
            MemFree(amaip->saps);
         sap->master = AlnMgrFindMaster(sap);
         if (sap->master && AlnMgrCheckRealMaster(sap, sap->master))
         {
            AlnMgrSetMaster(sap, sap->master);
            AlnMgrMakeMasterPlus(sap);
            n = AlnMgrGetNForSip(sap, sap->master);
            sap->type = SAT_MASTERSLAVE;
            amaip->numseg = AlnMgrGetMaxSegments((SeqAlignPtr)(sap->segs));
            amaip->alnsaps = amaip->numsaps;
            amaip->lens = (Int4Ptr)MemNew((amaip->numseg)*sizeof(Int4));
            amadp = amaip->amadp[n-1];
            amaip->saps = amadp->saps;
            if (!AlnMgrMergeIntoMSMultByMaster(amaip, amaip->lens, &amaip->numseg))
               return retval;
            amaip->ids = SeqIdSetFree(amaip->ids);
            amaip->ids = AlnMgrPropagateSeqIdsBySapList(amaip);
            amaip->starts = (Int4Ptr)MemNew((amaip->numseg)*(amaip->numsaps)*sizeof(Int4));
            amaip->aligncoords = (Uint4Ptr)MemNew((amaip->numseg)*sizeof(Uint4));
            if (!AlnMgrFillInStarts(amadp->saps, amaip->starts, amaip->numseg, amaip->lens, amaip->numsaps, amaip->aligncoords))
               return retval;
            if (amaip->numseg > 1)
               amaip->numseg -= 1;
            if (!AlnMgrMakeMultSegments(amaip))
               return retval;
            if (!AlnMgrGetRowsForMasterSlave(sap))
               return retval;
            retval = TRUE;
         } else
         {
         }
      }
   }
   return retval;
}


void AlnMgrMakeAlignCoords(SeqAlignPtr sap)
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

Boolean AlnMgrMergeIntoMSMultByMaster(AMAlignIndexPtr amaip, Int4Ptr lens, Uint4Ptr numseg)
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
   for (i=0; i<=(count-1); i++)
   {
      printf("%d %d %d\n", tiparray[i]->start, tiparray[i]->numgap, tiparray[i]->which);
   }
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

Boolean AlnMgrMergeSegments(Int4Ptr lens, SeqAlignPtr sap, SeqIdPtr master, Int4Ptr where, Int4 which)
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


Boolean AlnMgrFillInStarts(SeqAlignPtr PNTR saparray, Int4Ptr starts, Int4 numseg, Int4Ptr lens, Int4 numsaps, Uint4Ptr aligncoords)
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

Int4 AlnMgrGetStartFromMaster(SeqAlignPtr sap, Int4 pos)
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

Uint4 AlnMgrGetMasterGapStartForSeg(SeqAlignPtr sap, Int4 which_gap, Uint4Ptr aligncoord)
{
   DenseSegPtr  dsp;
   SAIndexPtr   saip;

   saip = (SAIndexPtr)(sap->saip);
   dsp = (DenseSegPtr)(sap->segs);
   if (aligncoord)
      *aligncoord = dsp->lens[saip->ssdp[saip->master-1]->unsect[which_gap]];
   return saip->aligncoords[saip->ssdp[saip->master-1]->unsect[which_gap]];
}


Boolean AlnMgrReconcileGaps(Int4Ptr lens, Uint4Ptr aligncoords, Int4 num)
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

Boolean AlnMgrMakeMultSegments(AMAlignIndexPtr amaip)
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

Int4 AlnMgrCheckOverlapping(SeqAlignPtr sap)
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

Int4 AlnMgrGetMaxSegments(SeqAlignPtr sap)
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
Int4 AlnMgrGetNumRows(SeqAlignPtr sap)
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

Int4 AlnMgrGetMaxRowsForParentPartial(SeqAlignPtr sap)
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

Boolean AlnMgrGetRowsForPartial(SeqAlignPtr sap)
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
            if (SeqIdComp(sip, rowsource[k]->id))
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
            rsp->id = sip;
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

Boolean AlnMgrGetRowsForMasterSlave(SeqAlignPtr sap)
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
   rsp->id = sap->master;
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
            rsp->id = sip;
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

SeqIdPtr AlnMgrFindMaster(SeqAlignPtr sap)
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
               if (!SeqIdComp(sip, salp->master))
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
Boolean AlnMgrCheckRealMaster(SeqAlignPtr sap, SeqIdPtr master)
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
            if (SeqIdComp(sip, master))
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
         if (SeqIdComp(sip, master))
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

void AlnMgrSetMaster(SeqAlignPtr sap, SeqIdPtr master)
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

void AlnMgrMakeMasterPlus(SeqAlignPtr sap)
{
   DenseSegPtr      dsp;
   Int4             i;
   SAIndexPtr       saip;
   SeqAlignPtr      sap_tmp;

   i = AlnMgrCheckAlignForParent(sap);
   if (i==AM_CHILD)
   {
      saip = (SAIndexPtr)(sap->saip);
      if (saip->master < 0)
         return;
      dsp = (DenseSegPtr)(sap->segs);
      if (dsp->strands[saip->master] == Seq_strand_minus)
      {
         sap_tmp = sap;
         sap_tmp->next = NULL;
         sap_tmp = SeqAlignListReverseStrand(sap_tmp);
         if (!AlnMgrIndexLinkedSegs(sap_tmp))
            return;
         sap_tmp->next = sap->next;
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
