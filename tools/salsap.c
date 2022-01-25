/*   salsap.c
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
* File Name:  salsap.c
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.61 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/
#include <salsap.h>
#include <salutil.h>
#include <sqnutils.h>
#include <edutil.h>
#include <objmgr.h>
#include <seqmgr.h>

extern Boolean is_dim1seqalign (SeqAlignPtr salp)
{
  if (salp->dim == 1) { 
     return TRUE;
  }
  return FALSE;
}  

extern Boolean is_dim2seqalign (SeqAlignPtr salp)
{
  SeqAlignPtr nextsalp;
  DenseSegPtr dsp;
  SeqIdPtr    sip;

  if (salp->dim == 2 && salp->next != NULL) {
     if (salp->segtype == 2) {
        dsp = (DenseSegPtr) salp->segs;
        sip = dsp->ids;
        nextsalp = salp->next;
        while (nextsalp) {
           if (nextsalp->dim != 2)
              break;
           if (nextsalp->segtype != 2) 
              break;
           dsp = (DenseSegPtr) nextsalp->segs;
           if (! SeqIdForSameBioseq (sip, dsp->ids))
              break;
           nextsalp = nextsalp->next;
        }
        if (nextsalp == NULL) 
           return TRUE;
     }
  }
  return FALSE;
}

extern Boolean is_fasta_seqalign (SeqAlignPtr salp)
{
  DenseSegPtr dsp;
  Int4Ptr     startp;
  Boolean     gap;
  Int4        k;
  Int2        j;
  Boolean     one_gap = FALSE;

  if (salp->dim > 2) {
     if (salp->segtype == 2) {
        dsp = (DenseSegPtr) salp->segs;
        for (j=0; j<dsp->dim; j++) {
           for (k=0; k<dsp->numseg; k++) {
              startp=dsp->starts;
              if (startp[dsp->dim*k + j] < 0) {
                 gap = TRUE;
                 one_gap = TRUE;
              }
              else if (gap)
                 return FALSE;
           }
        }
     }
  }
  else 
     return FALSE;
  if (!one_gap)
     return FALSE;
  return TRUE;
}

extern SeqAnnotPtr SeqAnnotForSeqAlign (SeqAlignPtr salp)
{
  SeqAnnotPtr sap;

  if (salp==NULL)
     return NULL;
  sap = SeqAnnotNew ();
  sap->type = 2;
  sap->data = (Pointer) salp;
/*
  seqannot->desc = ValNodeCopyStr(NULL, Annot_descr_name, "blastp");
*/
  return sap;
}

extern SeqAlignPtr is_salp_in_sap (SeqAnnotPtr sap, Uint1 choice)
{
  SeqAlignPtr      salp = NULL;

  if (sap != NULL) {
     for (; sap!= NULL; sap=sap->next) {
        if (sap->type == choice) {
           salp = (SeqAlignPtr) sap->data;
           return salp;
        }
     }   
  }
  return NULL;
}


extern SeqAlignPtr seqalign_list_free (SeqAlignPtr salp)
{
  SeqAlignPtr tmp,
              next;
  tmp = salp;
  while (tmp!=NULL) {
     next = tmp->next;
     tmp->next = NULL;
     SeqAlignFree (tmp);
     tmp = next;
  }
  return NULL;
}


static SeqIdPtr SeqAlignIDList (SeqAlignPtr salp)
{
  DenseSegPtr  dsp;
  DenseDiagPtr ddp;
  CompSegPtr   csp;
  SeqIdPtr     sip = NULL;

  if (salp!=NULL) {
     if (salp->segtype == 1) {
        ddp = (DenseDiagPtr) salp->segs;
        if(ddp!=NULL) 
           sip = ddp->id;
     }
     else if (salp->segtype == 2) {
        dsp = (DenseSegPtr) salp->segs;
        if (dsp!=NULL) 
           sip = dsp->ids;
     } 
     else if (salp->segtype == 4) {
        csp = (CompSegPtr) salp->segs;
        if (csp!=NULL)
           sip = csp->ids;
     }
  }
  return sip;
}

extern SeqIdPtr SeqAlignId (SeqAlignPtr salp, Int2 index)
{
  SeqIdPtr    sip = NULL;
  Int2        j;

  if (salp!=NULL) {
     sip = SeqAlignIDList (salp);
     if (sip != NULL)
     {
        for (j=0; j<index; j++)
           sip = sip->next;
     }
  }
  return sip;
}

extern Boolean FindSeqIdinSeqAlign (SeqAlignPtr salphead, SeqIdPtr sip)
{
  SeqAlignPtr salp;
  SeqIdPtr    tmpsip;

  for (salp = salphead; salp!= NULL; salp=salp->next)
  {
     tmpsip = SeqAlignIDList (salp);
     if ((position_inIdlist(sip, tmpsip)) > 0)
     {
        return TRUE;
     }
  }
  return FALSE;
}


extern SeqAlignPtr SeqAlignIdReplace (SeqAlignPtr salp, Int2 index, SeqIdPtr newsip)
{
  DenseSegPtr  dsp;
  DenseDiagPtr ddp;
  SeqIdPtr     sip = NULL,
               presip;
  Int2         j;

  if (newsip==NULL)
     return salp;
  newsip = SeqIdDup(newsip);
  if (salp!=NULL) {
     if (salp->segtype == 2) {
        dsp = (DenseSegPtr) salp->segs;
        if (dsp!=NULL) {
           sip = dsp->ids;
           if (index==0)
              dsp->ids = newsip;
           else
           {
              presip = dsp->ids;
              sip=sip->next;
              for (j=0; j<index-1; j++, sip=sip->next)
                 presip = presip->next;
              presip->next = newsip;
           }
           newsip->next = sip->next;
           sip->next = NULL;
           SeqIdFree (sip);
        }
     } 
     else if (salp->segtype == 1) {
        ddp = (DenseDiagPtr) (salp->segs);
        if(ddp!=NULL) {
           sip = ddp->id;
           if (index==0)
              ddp->id = newsip;
           else
           {
              presip = ddp->id;
              sip=sip->next;
              for (j=0; j<index-1; j++, sip=sip->next)
                 presip = presip->next;
              presip->next = newsip;
           }
           newsip->next = sip->next;
           sip->next = NULL;
           SeqIdFree (sip);
        }
     }
  }
  return salp;
}

extern Boolean SeqAlignSeqLocComp (SeqAlignPtr salphead, ValNodePtr vnp)
{
  SeqAlignPtr salp;
  ValNodePtr tmp;
  SeqLocPtr  slp;
  SeqIdPtr   sip,
             tmpsip;

  for (tmp=vnp; tmp!=NULL; tmp=tmp->next)
  {
     slp= (SeqLocPtr)tmp->data.ptrvalue;
     if (slp)
     {
        sip=SeqLocId(slp);
        for (salp = salphead; salp!= NULL; salp=salp->next)
        {
           tmpsip = SeqAlignIDList (salp);
           if ((position_inIdlist(sip, tmpsip)) > 0)
           {
              return TRUE;
           }
        }
     }
  } 
  return FALSE;
}


extern Int4 SeqAlignLength (SeqAlignPtr salp)
{
  DenseSegPtr  dsp;
  DenseDiagPtr ddp;
  CompSegPtr   csp;
  Int4Ptr      lenp;
  Int4         lens = 0;
  Int2         j;

  if (salp!=NULL) {
     if (salp->segtype == 2) {
        dsp = (DenseSegPtr) salp->segs;
        if (dsp) { 
           lenp = (Int4Ptr) dsp->lens;
           for (j=0; j<dsp->numseg; j++, lenp++)
              lens += *lenp;
        }
     } 
     else if (salp->segtype == 1) {
       ddp = (DenseDiagPtr)salp->segs;
       if(ddp!=NULL) 
          lens=(Int4)ddp->len;
     }
     else if (salp->segtype == 4) {
        csp = (CompSegPtr)salp->segs;      
        if (csp) {
          lenp = (Int4Ptr) csp->lens;
          for (j=0; j<csp->numseg; j++, lenp++)
             lens += *lenp;
        }
     }
  }
  return lens;
}

extern Uint1 SeqAlignStrand (SeqAlignPtr salp, Int2 index)
{
  DenseDiagPtr ddp;
  DenseSegPtr  dsp;
  Uint1        strand = Seq_strand_unknown;

  if (salp == NULL)
     return 0;
  if (salp->segtype == 1)
  {
     ddp = (DenseDiagPtr) salp->segs;
     if (ddp != NULL) {
        if ((Boolean)(ddp->strands != NULL))
           strand = *(ddp->strands);
     }
  }
  else if (salp->segtype == 2)
  {
     dsp = (DenseSegPtr) salp->segs;
     if (dsp!=NULL)
     {
        if ((Boolean)(dsp->strands != NULL)) {
           strand = (dsp->strands)[index];
           while (strand!=Seq_strand_plus && strand!=Seq_strand_minus) {
              index+=dsp->dim;
              if (index >= dsp->dim*dsp->numseg)
                 break;
              strand = (dsp->strands)[index];
           }
        }
     }
  }
  return strand;
}

extern Int4 SeqAlignStart (SeqAlignPtr salp, Int2 index)
{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  CompSegPtr    csp;
  Int4Ptr       startp;
  Int4          j, 
                val = (Int4)-1;
  Uint1         strand;
  
  if (salp == NULL)
     return -1;
  if (salp->segtype == 1) 
  {
     ddp = (DenseDiagPtr) salp->segs;
     if (ddp != NULL) {
        startp = ddp->starts;
        if (index > 0)
           startp += index;
        val = *startp;
     }
  } 
  else if (salp->segtype == 2) 
  {
     strand = SeqAlignStrand (salp, index);
     dsp = (DenseSegPtr) salp->segs;
     if (dsp!=NULL) 
     {
        startp = dsp->starts;
        if (index > 0)
           startp += index;
        for (j = 0; j < dsp->numseg; j++, startp += dsp->dim) 
           if (*startp > -1)
              break;
        if (j < dsp->numseg) {
           if (strand == Seq_strand_minus)
              val = *startp + dsp->lens[j] - 1;
           else
              val = *startp;              
        }
     }
  }
  else if (salp->segtype == 4)
  {
     csp = (CompSegPtr) salp->segs;
     val = csp->from[index];
  }
  return val;
}

extern Int4 SeqAlignStop (SeqAlignPtr salp, Int2 index)
{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  CompSegPtr    csp;
  Int4Ptr       startp;
  BoolPtr       startb;
  Int4          j, 
                val = (Int4)-1;
  Uint1         strand;
  
  if (salp == NULL)
     return -1;
  if (salp->segtype == 1) 
  {
     ddp = (DenseDiagPtr) salp->segs;
     if (ddp != NULL) {
        startp = ddp->starts;
        if (index > 0)
           startp += index;
        val = *startp + ddp->len;
     }
  } 
  else if (salp->segtype == 2) 
  {
     dsp = (DenseSegPtr) salp->segs;
     if (dsp!=NULL) 
     {
        if ((Boolean)(dsp->strands != NULL))
           strand = dsp->strands[index];
        else
           strand = Seq_strand_plus;
        startp = dsp->starts + ((dsp->dim * dsp->numseg) - dsp->dim);
        if (index > 0)
           startp += index;
        for (j = dsp->numseg-1; j >= 0; j--, startp-=dsp->dim) 
           if (*startp > -1)
              break;
        if (j >= 0) {
           if (strand == Seq_strand_minus)
              val = *startp;
           else
              val = *startp + dsp->lens[j] - 1;              
        }
     }
  }
  else if (salp->segtype == 4)
  {
     csp = (CompSegPtr) salp->segs;
     val = csp->from[index];
     startb=csp->starts;
     for (j=0; j<csp->numseg; j++)
     {
        if (*startb)
           val+=csp->lens[j]; 
        startb+=csp->dim;
     } 
  }
  return val;
}

extern Uint1 SeqAlignMolType (SeqAlignPtr salp)
{
  DenseSegPtr  dsp;
  DenseDiagPtr ddp;
  StdSegPtr    ssp;
  BioseqPtr    bsp;
  SeqIdPtr     sip = NULL;
  Int2         k;
  Int2         dim;
  Uint1        moltype = 0;
  Boolean      molb;
 
  if (salp==NULL)
     return FALSE;
  if (salp->segtype == 1) {
     ddp = (DenseDiagPtr) salp->segs;
     sip = ddp->id;
     dim = ddp->dim;
  }
  else if (salp->segtype == 2) {
     dsp = (DenseSegPtr) salp->segs;
     sip = dsp->ids;
     dim = dsp->dim;
  }
  else if (salp->segtype == 3) {
     ssp = (StdSegPtr) salp->segs;
     sip = ssp->ids;
     dim = ssp->dim; 
  }
  if (sip!=NULL) {
   for (k = 0; k < dim && sip!=NULL; k++, sip = sip->next)
   {
     bsp = BioseqLockById (sip);
     if (bsp!=NULL)
     {   
        break;
     }
   }
   if (bsp!=NULL) {
      moltype = bsp->mol;
      molb = (Boolean) (ISA_aa(moltype));
      BioseqUnlock (bsp);
      for (; k < dim && sip!=NULL; k++, sip = sip->next)
      {
         bsp = BioseqLockById (sip);
         if (bsp!=NULL)
         {   
            if ((Boolean) (ISA_aa(bsp->mol)) != molb) {
               BioseqUnlock (bsp);
               moltype = 0;
               ErrPostEx (SEV_ERROR, 0, 0, "All Source Bioseq's are not Nucleotide");
               break;
            }
            BioseqUnlock (bsp);
         }
      }
   }
  }
  return moltype;
}

extern Int4 SeqAlignBestScore (SeqAlignPtr salp)
{
  ScorePtr    score;
  ObjectIdPtr oip;
  Int4        val = 0;

  if (salp!=NULL) {
     if (salp->segtype == 2) {
        for (score = (ScorePtr) salp->score; score!=NULL; score=score->next)
        {
           oip = (ObjectIdPtr) score->id;
           if (oip != NULL) {
              if (StringStr(oip->str, "score")!=NULL) {
                 if (score->choice == 1) {
                    val = (Int4)score->value.intvalue; 
                 }
              }
           }
        }
     }
  }
  return val;
}

extern SeqAlignPtr SeqAlignBestHit (SeqAlignPtr salp, Int4 length, Int4 threshold)
{
  SeqAlignPtr  tmp;
  Int4         len,
               gap;
  float        val;

  for (tmp=salp; tmp!=NULL; tmp=tmp->next)
  {
     len = SeqAlignLength(salp);
     if (100*((float)len/(float)length) >= threshold)
        return tmp;
  }
  return NULL;
}

extern Int4 SeqAlignGapCount (SeqAlignPtr salp)
{
  DenseSegPtr   dsp;
  Int4          sum = 0,
                j;
  Int2          k;

  if (salp == NULL)
     return 0;
  if (salp->segtype != 2)
     return 0;
  if (salp->segtype == 2)
  {
     dsp = (DenseSegPtr) salp->segs;
     if (dsp!=NULL)
     {
        for (j=0; j<dsp->numseg; j++)
           for (k=0; k<dsp->dim; k++)
           {
              if (dsp->starts[(dsp->dim*j)+k] < 0)
                 sum+= dsp->lens[j];
           }
     }
  }
  return sum;
}


extern void SeqAlignStartUpdate (SeqAlignPtr salp, SeqIdPtr target_sip, Int4 offset, Uint1 strand)
{
  SeqAlignPtr salptmp;
  DenseSegPtr dsp;
  SeqIdPtr    pre, sip, next;
  Int4Ptr     lenp,
              startp;
  Uint1Ptr    strandp;
  Int4        len_sum,
              start;
  Int2        index, k, j;

  if (salp==NULL || offset<=0)
     return;
  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next)
  {
     if (salptmp->segtype == 2)
     {
        dsp = (DenseSegPtr) salptmp->segs;
        pre = NULL;
        index=0;
        sip=dsp->ids;
        while (sip) 
        {
           next=sip->next;
           if (SeqIdForSameBioseq(target_sip, sip))
           {
              if (strand == Seq_strand_minus)
              {
                 strandp=dsp->strands;
                 strandp+=index;
                 for (j=0; j < dsp->numseg && strandp!=NULL; j++)
                 {
                    if (*strandp == Seq_strand_minus)
                       *strandp = Seq_strand_plus;
                    else if (*strandp == Seq_strand_plus)
                       *strandp = Seq_strand_minus;
                    strandp+=dsp->dim;
                 }
                 lenp=dsp->lens;
                 startp=dsp->starts;
                 start=startp[index];
                 len_sum=0;
                 j=dsp->dim*dsp->numseg-dsp->dim;
                 k=dsp->numseg-1;
                 for (; j>=0; j-=dsp->dim) {
                    startp[j+index]=start+len_sum;
                    len_sum+=lenp[k];
                    k--;
                 }
              }
              for (j=0; j<dsp->numseg; j++) {
                 if (dsp->starts[dsp->dim*j+index] != -1)
                    dsp->starts[dsp->dim*j+index] += offset;
              }
           }
           pre=sip;
           sip=next;
           index++;
        }
     }
  }
}

/**********************************************************************
*
*       SeqAlignLink(head, new)
*       link the new align to the end of head align. return the
*       start of the linked chain
*
**********************************************************************/

extern SeqAlignPtr SeqAlignLink(SeqAlignPtr head, SeqAlignPtr a_new)
{
        SeqAlignPtr curr;

        curr = head;
        if(curr!= NULL){
                while(curr->next != NULL)
                        curr= curr->next;
                curr->next = a_new;
                return head;
        }
        else
                return a_new;

}

/***********************************************************************
***    
***    SeqLocListFromSeqAlign
***      read SeqAlignPtr
***      return list of ValNodePtr-SeqLocPtr
***
************************************************************************/
extern ValNodePtr SeqLocListFromSeqAlign (SeqAlignPtr salp)
{
  SeqAlignPtr   salptmp;
  ValNodePtr    vnp = NULL;
  DenseSegPtr   dsp;
  DenseDiagPtr  ddp;
  SeqLocPtr     slp;
  SeqIdPtr      sip;
  Int4          start, 
                stop;
  Uint1         strand; 
  Int2          offset;
  
  if (salp==NULL) 
     return NULL;
  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next)
  {
     sip = SeqAlignId (salptmp, 0);
     if (sip != NULL)
     {
        offset = 0;
        for (; sip != NULL; sip = sip->next)
        {
           start = SeqAlignStart (salptmp, offset);
           stop = SeqAlignStop (salptmp, offset);
           strand = SeqAlignStrand (salptmp, offset);
           slp = SeqLocIntNew (start, stop, strand, sip);
           if (slp!=NULL)
              ValNodeAddPointer (&vnp, 0, (Pointer) slp);
           offset++;
        }
     }
  }
  return vnp;
}



extern SeqLocPtr SeqLocFromSeqAlign (SeqAlignPtr salp, SeqIdPtr sip)
{
  SeqLocPtr     slp = NULL;
  SeqIdPtr      sip1 = NULL, 
                siptmp;
  Int4          start, 
                stop;
  Uint1         strand; 
  Int2          offset;
 
  if (salp != NULL)
  {
     offset = 0;
     sip1=SeqAlignId(salp, 0);
     if (sip  == NULL)
        siptmp = sip1;
     else 
     {
        sip = SeqIdFindBest (sip, 0);
        for (siptmp=sip1; siptmp!=NULL; siptmp=siptmp->next, offset++) {
           if (SeqIdForSameBioseq (sip,siptmp))
              break;
        }
     }
     if (siptmp!=NULL)
     {
        start = SeqAlignStart (salp, offset);
        stop = SeqAlignStop (salp, offset);
        strand = SeqAlignStrand (salp, offset);
        slp = SeqLocIntNew (start, stop, strand, siptmp);
     }
     
  }
  return slp; 
}

extern SeqLocPtr SeqLocMixFromSeqAlign (SeqAlignPtr salp, SeqIdPtr sip)
{
  DenseSegPtr dsp;
  SeqLocPtr   headslp = NULL,
              slp = NULL,
              tmp;
  SeqIdPtr    siptmp;
  Int4Ptr     lens, start,
              start2;
  Int4        salp_start, salp_stop;
  Int2        offset, j;

  if (salp == NULL)
     return NULL;
  if (salp->segtype != 2)
     return NULL;
  dsp =  (DenseSegPtr) salp->segs;
  if (dsp == NULL)
     return NULL;
  offset = 0;
  if (sip  == NULL)
     sip = dsp->ids;
  else {
     sip = SeqIdFindBest (sip, 0);
     for (siptmp=dsp->ids; siptmp!=NULL; siptmp=siptmp->next, offset++) {
        if (SeqIdForSameBioseq (sip,siptmp))
           break;
     }
     if (siptmp == NULL)
        return NULL;
  }
  start=dsp->starts + offset;
  if (offset==0)
     start2 = start+1;
  else
     start2 = dsp->starts;
  lens=dsp->lens;
  for (j=0; j<dsp->numseg; j++) { 
     if (*start > -1) 
        break;
     start +=dsp->dim;
     start2+=dsp->dim;
     lens++;
  }
  for (; j<dsp->numseg; j++, lens++) {
     if (*start > -1 && *start2 > -1) {
        salp_start = *start;
        salp_stop = salp_start + *lens -1;
        slp = SeqLocIntNew (salp_start, salp_stop, Seq_strand_plus, sip);
        if (headslp == NULL) {
           headslp = slp;
        }
        else if (headslp->choice == SEQLOC_MIX) {
           tmp = (ValNodePtr)(headslp->data.ptrvalue);
           while (tmp->next != NULL)
              tmp = tmp->next;
           tmp->next = slp;
        }
        else {
           tmp = ValNodeNew(NULL);
           tmp->choice = SEQLOC_MIX;
           tmp->data.ptrvalue = (Pointer)headslp;
           headslp = tmp;
           tmp = (ValNodePtr)(headslp->data.ptrvalue);
           tmp->next = slp;
        }
     }
     start +=dsp->dim;
     start2+=dsp->dim;
  }
  return headslp; 
}

static Boolean is_id_in (SeqIdPtr sip, ValNodePtr vnp)
{
  ValNodePtr vnptmp;
  SeqIdPtr   siptmp;

  for (vnptmp=vnp; vnptmp!=NULL; vnptmp=vnptmp->next)
  {
     siptmp = SeqLocId((SeqLocPtr)vnptmp->data.ptrvalue);
     if (siptmp!=NULL) {
        if (SeqIdForSameBioseq(sip, siptmp))
           return TRUE;
     }
  }
  return FALSE;
}

extern ValNodePtr WholeSeqLocListFromSeqAlign (SeqAlignPtr salp)
{
  SeqAlignPtr   salptmp;
  BioseqPtr     bsp;
  ValNodePtr    vnp = NULL;
  DenseSegPtr   dsp;
  DenseDiagPtr  ddp;
  SeqIdPtr      sip,
                siptmp;
  SeqLocPtr     slp;

  if (salp==NULL) 
     return NULL;
  for (salptmp = salp; salptmp!=NULL; salptmp=salptmp->next)
  {
     sip = NULL;
     if (salptmp->segtype == 1)
     {
        ddp = (DenseDiagPtr) salptmp->segs;
        sip = ddp->id;        
     }
     else if (salptmp->segtype == 2)
     {
        dsp = (DenseSegPtr) salptmp->segs;
        sip = dsp->ids;
     }
     if (sip != NULL) 
     {
        for (; sip != NULL; sip = sip->next)
        {
           siptmp = SeqIdDup (sip); 
           if (vnp==NULL || (vnp!=NULL && !is_id_in (siptmp, vnp)))
           {
              bsp = BioseqLockById (siptmp);
              if (bsp!= NULL) 
              {
                 slp = SeqLocIntNew (0, bsp->length-1, Seq_strand_plus, siptmp);
                 if (slp != NULL)
                 {
                    ValNodeAddPointer (&vnp, 0, (Pointer) slp);
                 }
                 BioseqUnlock (bsp);
              }
           }
        }
     }
  }
  return vnp;
}


static SeqIdPtr MySeqIdSetFree (SeqIdPtr sip)
{
  SeqIdPtr next;
 
    while (sip != NULL)
    {
        next = sip->next;
        sip->next = NULL;
        if (sip!=NULL) {
            if (sip->choice > 0 && sip->choice < 30) SeqIdFree(sip);
        }
        else break;
        sip = next;
    }
    return NULL;
}

extern SeqAlignPtr CompSeqAlignFree (SeqAlignPtr salp)
{
  CompSegPtr csp;
 
  csp = (CompSegPtr) salp->segs;
  if (csp->ids != NULL) MySeqIdSetFree (csp->ids);
  csp->ids=NULL;
  if (csp->from != NULL) MemFree (csp->from);
  csp->from=NULL;
  if (csp->lens != NULL) MemFree (csp->lens);
  csp->lens=NULL;
  if (csp->starts!= NULL) MemFree (csp->starts);
  csp->starts=NULL;
  MemFree ( csp);
  salp->segs = NULL;
  SeqAlignFree (salp);
  return NULL;
}

extern SeqAnnotPtr CompSeqAnnotFree (SeqAnnotPtr sap)
{
  SeqAlignPtr salp;

  if ( sap == NULL ) return NULL; 
  if ( sap->type != 2 ) {
      ErrPostEx(SEV_WARNING, 0, 0, "fail in CompSeqAnnotFree [1] %d", (int)sap->type); 
      return NULL; 
  }
  if ( ( salp = (SeqAlignPtr) sap->data ) == NULL ) {
      ErrPostEx(SEV_WARNING, 0, 0, "fail in CompSeqAnnotFree [2]"); 
      return NULL; 
  }
  if ( salp->segtype == 4 ) 
      CompSeqAlignFree (salp);
  else if ( salp->segtype == 2 )
      SeqAnnotFree (sap);
  else {
      ErrPostEx(SEV_WARNING, 0, 0, "fail in CompSeqAnnotFree [3] %d", salp->segtype); 
      return NULL; 
  }
  sap->data = NULL;
  return NULL;
}

extern SeqAnnotPtr SeqAlignBoolSegCpy (SeqAnnotPtr sap, Int4 from, Int4 to)
{
  SeqAnnotPtr sapnew =NULL;
  SeqAlignPtr salp, salpnew =NULL, salpre =NULL;
  CompSegPtr  csp =NULL;
  CompSegPtr  newcsp =NULL;
  BoolPtr     cspstart =NULL;
  Int4Ptr     csplens =NULL, newcsplens =NULL;
  Int4Ptr     cspfrom =NULL, newcspfrom =NULL;
  Int2        j, k, nbseq, numseg, newnumseg;
  Int4        sumlens, offset;

  if ( sap == NULL ) return NULL; 
  if ( sap->type != 2 ) {
      ErrPostEx(SEV_WARNING, 0, 0, "fail in SeqAlignBoolSegCpy [1]"); 
      return NULL; 
  }
  if ( ( salp = (SeqAlignPtr) sap->data ) == NULL ) {
      ErrPostEx(SEV_WARNING, 0, 0, "fail in SeqAlignBoolSegCpy [1-2]"); 
      return NULL; 
  }
  if ( salp->segtype != 4 ) {
      ErrPostEx(SEV_WARNING, 0, 0, "fail in SeqAlignBoolSegCpy [2]"); 
      return NULL; 
  }
  sapnew = SeqAnnotNew ();
  if (sapnew == NULL) {
      return NULL; 
  }
  sapnew->type = 2;

  while ( salp != NULL )
  {
      salpnew = SeqAlignNew ();
      if ( salpnew == NULL ) {
          SeqAnnotFree (sapnew);
          return NULL; 
      }
      salpnew->type = 3;
      salpnew->segtype = 4;
      salpnew->dim = salp->dim;
      csp = (CompSegPtr) salp->segs;
      newcsp = (CompSegPtr) MemNew ( sizeof (CompSeg) );
      if ( newcsp == NULL ) {
          ErrPostEx(SEV_WARNING, 0, 0, "fail in SeqAlignBoolSegCpy [4]"); 
          SeqAnnotFree (sapnew);
          return NULL; 
      }
      salpnew->segs = (Pointer) newcsp;
      nbseq = newcsp->dim = csp->dim;
      newcsp->ids = SeqIdDupList (csp->ids);
      numseg = newcsp->numseg = csp->numseg;

          /* copy the first position (from) + sum of length */
      newcsp->from = (Int4Ptr) MemNew ((size_t) ((nbseq +2) * sizeof (Int4)));
      if ( newcsp->from == NULL ) {
          SeqAnnotFree (sapnew);
          return NULL; 
      }
      for (j = 0; j < nbseq; j++) {
          offset = -1;
          sumlens = 0;
          csplens = csp->lens;
          cspstart = csp->starts;
          cspstart += j;
          for (k = 0; k < numseg ; k++, csplens++, cspstart +=nbseq) {
               if ( from >= sumlens && from < sumlens + *csplens ) {
                    if ( (Boolean)(*cspstart) ) offset =(sumlens +*csplens -from);
                    break;
               }
               if ( (Boolean)(*cspstart) ) sumlens += (*csplens);
          }
          if ( offset < 0 ) 
               for (; k < numseg && offset < 0; k++) {
                    if ( (Boolean)(*cspstart) ) offset = 0;
               }
          newcsp->from [j] = csp->from [j] +offset;
      }
          /* count number of new segments between from and to */
      sumlens =0;
      newnumseg =0;
      cspstart = csp->starts;
      for (k=0; k<numseg; k++, sumlens += (Int4)(*csplens), csplens++, cspstart += (Int4)(nbseq))
      {
          if ( from >= sumlens && from < sumlens + *csplens ) break;
      }
      for (; k<numseg; k++, sumlens += *csplens, csplens++, cspstart +=nbseq)
      {
          newnumseg++;
          if ( to >= sumlens && to < sumlens + *csplens ) break;
      }
          /* copy segment lengths within from and to */
      newcsp->lens =(Int4Ptr) MemNew ((size_t) ((newnumseg+2) * sizeof (Int4)));
      if ( newcsp->lens == NULL ) {
          SeqAnnotFree (sapnew);
          return NULL; 
      }
      sumlens =0;
      csplens = csp->lens;
      cspstart = csp->starts;
      newcsplens = newcsp->lens;
      for (k=0; k<numseg; k++, sumlens += *csplens, csplens++, cspstart +=nbseq)
      {
          if ( from >= sumlens && from < sumlens + *csplens ) break;
      }
      *newcsplens = sumlens +*csplens -from;
      newcsplens++;
      csplens++;
      k++;
      for (; k<numseg; k++, sumlens += *csplens, csplens++, cspstart +=nbseq)
      {
          *newcsplens = *csplens;
          if ( to >= sumlens && to < sumlens + *csplens ) break;
          newcsplens++;
      }
      *newcsplens = sumlens -to;
      newcsp->starts =(BoolPtr) MemNew((size_t)((nbseq *newnumseg +2) * sizeof(Boolean)));
      if ( newcsp->starts == NULL ) {
          SeqAnnotFree (sapnew);
          return NULL; 
      }
      for (j = 0; j < nbseq; j++) 
      {
          sumlens = 0;
          csplens = csp->lens;
          cspfrom = csp->from;
          cspfrom += j;
          newcspfrom = newcsp->from;
          newcspfrom += j;
          for(k=0;k<numseg; k++,sumlens +=*csplens, csplens++, cspstart +=nbseq)
          {
               if ( from >= sumlens && from < sumlens + *csplens ) break;
          }
          for (; k<numseg; k++, sumlens +=*csplens, csplens++, cspstart +=nbseq)
          {
               *newcsplens = *csplens;
               if ( to >= sumlens && to < sumlens + *csplens ) break;
          }
      }
/*
      csp->strands
      csp->scores
*/
      if ( salpre == NULL ) sapnew->data = (Pointer) salpnew;
      else salpre->next = salpnew;
      salpre = salpnew;
      salp = salp->next;
  }
  return sapnew; 
}

extern SeqAlignPtr SeqAlignDenseSegToBoolSeg (SeqAlignPtr salp)
{
  SeqAlignPtr salpnew =NULL;
  DenseSegPtr dsp =NULL;
  Int4Ptr     dspstarts =NULL, dspstarts2 = NULL;
  Int4Ptr     dsplens =NULL;
  Uint1Ptr    strandp;
  Uint1       strand;
  CompSegPtr  csp =NULL;
  Int4        j; 
  Int2        nbseq, numseg;

  salpnew = SeqAlignNew ();
  if ( salpnew == NULL ) {
          return NULL; 
  }
  salpnew->type = 3;
  salpnew->segtype = 4;
  salpnew->dim = salp->dim;

  salpnew->score = salp->score;
  salp->score = NULL;

  dsp = (DenseSegPtr) salp->segs;
  csp = (CompSegPtr) MemNew ( sizeof (CompSeg) );
  if ( csp == NULL ) {
      return NULL; 
  }
  salpnew->segs = (Pointer) csp;
  csp->dim = dsp->dim;
  nbseq = dsp->dim;

  csp->ids = SeqIdDupBestList (dsp->ids);
  csp->numseg = dsp->numseg;
  numseg = dsp->numseg;

  if (dsp->strands != NULL) {
     strandp = dsp->strands;
     csp->strands=(Uint1Ptr)MemNew((size_t)((nbseq*numseg+4)*sizeof (Uint1)));
     for (j=0; j<nbseq*numseg; j++, strandp++)
        csp->strands[j] = *strandp;
  }
  csp->lens = (Int4Ptr) MemNew ((size_t) ((numseg+2) * sizeof (Int4)));
  if ( csp->lens == NULL ) {
     return NULL; 
  }
  for (j = 0, dsplens = dsp->lens; j < numseg ; j++, dsplens++ ) 
     csp->lens [j] = *dsplens;
  csp->from = (Int4Ptr) MemNew ((size_t) ((nbseq +2) * sizeof (Int4)));
  if ( csp->from == NULL ) {
     return NULL; 
  }
  for (j = 0; j < nbseq +2; j++) 
     csp->from [j] = 0;
  dspstarts=dsp->starts;
  dsplens = dsp->lens;
  strandp = dsp->strands;
  if (strandp!=NULL)
     strand = *strandp;
  else 
     strand = Seq_strand_unknown;
  for (j = 0; j < nbseq ; j++, dspstarts++)
  {
     if (*dspstarts < 0) {
        for(dspstarts2=dspstarts;dspstarts2!=NULL;dspstarts2+=nbseq) {
           if (*dspstarts2 > -1) 
              break;
           dsplens++;
        }
        if (dspstarts2!=NULL) { 
           if (strand==Seq_strand_minus)
              csp->from [j] = *dspstarts2 + *dsplens;
           else
              csp->from [j] = *dspstarts2;
        }
     }
     else {
        if (strand==Seq_strand_minus) 
           csp->from [j] = *dspstarts + *dsplens;
        else
           csp->from [j] = *dspstarts;
     }
     if (strandp!=NULL) {
        strandp++;
        strand = *strandp;
     } else strand = Seq_strand_unknown;
  }
  csp->starts=(BoolPtr) MemNew((size_t)((nbseq *numseg +2) * sizeof(Boolean)));
  if ( csp->starts == NULL ) {
     return NULL; 
  }
  for (j =0, dspstarts =dsp->starts; j <nbseq*numseg ; j++, dspstarts++ ) {
          csp->starts [j] = (Boolean) (*dspstarts >= 0);
  }
  return salpnew;
}

extern SeqAnnotPtr SeqAnnotDenseSegToBoolSeg (SeqAnnotPtr sap)
{
  SeqAnnotPtr sapnew =NULL;
  SeqAlignPtr salp=NULL, salpnew=NULL, salpre=NULL;

  if ( sap == NULL ) return NULL; 
  if ( sap->type != 2 ) {
      return NULL; 
  }
  if ( ( salp = (SeqAlignPtr) sap->data ) == NULL ) {
      return NULL; 
  }
  if ( salp->segtype == 4 ) return sap;
  if ( salp->segtype != 2 ) {
      return NULL; 
  }
  sapnew = SeqAnnotNew ();
  if (sapnew == NULL) {
      return NULL; 
  }
  sapnew->type = 2;
  while ( salp != NULL )
  {
      salpnew = SeqAlignDenseSegToBoolSeg (salp);
      if ( salpre == NULL ) sapnew->data = (Pointer) salpnew;
      else salpre->next = salpnew;
      salpre = salpnew;
      salp = salp->next;
  }
  return sapnew; 
}

extern SeqAlignPtr SeqAlignBoolSegToDenseSeg (SeqAlignPtr salp)
{
  SeqAlignPtr salpnew =NULL;
  DenseSegPtr dsp =NULL;
  CompSegPtr  csp =NULL;
  BoolPtr     cspstarts =NULL;
  Int4Ptr     cspfrom =NULL;
  Int4Ptr     csplens =NULL;
  Uint1Ptr    strandp;
  Uint1       strand;
  Int4Ptr     dspstarts =NULL;
  Int4        sumlens, k;
  Int2        j, nbseq, numseg;

  salpnew = SeqAlignNew ();
  if (salpnew != NULL) {
     salpnew->type = 3;
     salpnew->segtype = 2;
     salpnew->dim = salp->dim;

     salpnew->score = salp->score;
     salp->score = NULL;

     csp = (CompSegPtr) salp->segs;
     dsp = DenseSegNew ();
     if ( dsp == NULL ) {
        return NULL; 
     }
     salpnew->segs = (Pointer) dsp;
     nbseq = csp->dim;
     dsp->dim = csp->dim;
     dsp->ids = SeqIdDupList (csp->ids);
     dsp->numseg = csp->numseg;
     numseg = csp->numseg;

     dsp->lens  =(Int4Ptr)MemNew((size_t)((numseg + 2) * sizeof (Int4))); 
     for (j = 0, csplens = csp->lens; j < numseg ; j++, csplens++ ) 
        dsp->lens [j] = *csplens;
     if (csp->strands != NULL) {
        strandp=csp->strands;
        dsp->strands =(Uint1Ptr)MemNew((size_t)((nbseq*numseg+4)*sizeof (Uint1)));
        for (j=0; j<nbseq*numseg; j++, strandp++)
           dsp->strands[j] = *strandp;
     }
     dsp->starts=(Int4Ptr)MemNew((size_t)((nbseq *numseg +4) * sizeof (Int4)));
     for (k = 0; k < nbseq *numseg +4; k++) 
        dsp->starts[k] = -1;
     cspstarts = csp->starts;
     dspstarts = dsp->starts;
     strandp = csp->strands;
     if (strandp!=NULL)
        strand = *strandp;
     else 
        strand = Seq_strand_unknown;
     for (j = 0, cspfrom = csp->from; j < nbseq ; j++, cspfrom++) 
     {
        csplens = csp->lens;
        sumlens = 0;
        for (k = 0; k < numseg; k++, csplens++) {
           if ( (Boolean)(cspstarts [nbseq *k +j]) ) {
              if (strand == Seq_strand_minus) {
                 sumlens += *csplens;
                 dspstarts [nbseq *k +j] = *cspfrom - sumlens;
              }
              else {
                 dspstarts [nbseq *k +j] = *cspfrom + sumlens;
                 sumlens += *csplens;
              }
           }
        }
        if (strandp!=NULL) {
           strandp++;
           strand = *strandp;
        } 
        else 
           strand = Seq_strand_unknown;
     }
  }
  return salpnew;
}

extern SeqAnnotPtr SeqAnnotBoolSegToDenseSeg (SeqAnnotPtr sap)
{
  SeqAnnotPtr sapnew =NULL;
  SeqAlignPtr salp, salpnew =NULL, salpre =NULL;

  if ( sap == NULL ) return NULL; 
  if ( sap->type != 2 ) {
      return NULL; 
  }
  if ( ( salp = (SeqAlignPtr) sap->data ) == NULL ) {
      return NULL; 
  }
  if ( salp->segtype == 2 ) return sap;
  if ( salp->segtype != 4 ) {
      return NULL; 
  }
  sapnew = SeqAnnotNew ();
  if (sapnew == NULL) return NULL;
  sapnew->type = 2;

  while ( salp != NULL )
  {
      salpnew = SeqAlignBoolSegToDenseSeg (salp);
      if ( salpre == NULL ) sapnew->data = (Pointer) salpnew;
      else salpre->next = salpnew;
      salpre = salpnew;
      salp = salp->next;
  }
  return sapnew; 
}


extern void CompSeqAlignPrint (SeqAlignPtr salp)
{
  CompSegPtr  csp =NULL;
  SeqIdPtr    sip;
  BoolPtr     startp;
  Int4        j, k;
  FILE        *fout;
  Char    strLog[128];

  csp = (CompSegPtr) salp->segs;
  if (csp!=NULL) {
     fout = FileOpen ("LogFile", "w");
     if (fout!=NULL) {
        fprintf (fout, "\n");
        for (j=0, sip=csp->ids; sip!=NULL; sip=sip->next, j++)
        {
            SeqIdWrite (sip, strLog, PRINTID_FASTA_LONG, 120);
            fprintf (fout, "%d %s \n", (int)(j+1), strLog);
        } 
        fprintf (fout, "\n");
        for (j=0; j<csp->dim; j++) 
           fprintf (fout, " %d ", (int)csp->from[j]); 
        fprintf (fout, "\n");
        startp = csp->starts;
        for (j = 0; j < csp->numseg; j++) {
            fprintf (fout, "%3d lg %6ld ", (int)j, (long)csp->lens[j]);
            for (k = 0; k < csp->dim; k++) { 
               fprintf (fout, " %d", (int)*startp);
               startp++;
            }
/*
            for ( k = 0; k < csp->dim; k++) { 
               if (csp->strands!=NULL)
                 fprintf (fout, " %d", (int)csp->strands[(Int4)csp->dim*k+j]);
            }
*/
            fprintf (fout, "\n"); 
        } 
        fprintf (fout, "\n"); 
        FileClose (fout);
     }
  }
}

extern SeqAlignPtr build_seqalign_fromstart (Int2 dim, Int2 numseg, SeqIdPtr sip, Int4Ptr starts, Int4Ptr lens)
{
  SeqAlignPtr salp;
  DenseSegPtr dsp;
  SeqIdPtr    next;

  salp = SeqAlignNew ();
  if (salp != NULL) {
     salp->type = 3;
     salp->segtype = 2;
     salp->dim = dim;
     dsp = DenseSegNew ();
     if (dsp != NULL) {
        salp->segs = (Pointer) dsp;
        dsp->dim = dim;
        dsp->ids = sip;
        dsp->numseg = numseg;
        dsp->starts = starts;
        dsp->lens = lens;
        return salp;
     }
  }
  MemFree (starts);
  while (sip != NULL) {
     next = sip->next;
     sip->next = NULL;
     SeqIdFree (sip);
     sip = next;
  }
  return NULL;
}

static Int4 getfirstpos_inseqalign (SeqAlignPtr salp, Int4 position, SeqIdPtr sip)
{
  CompSegPtr  dsp;
  BoolPtr     dspstart;
  Int4Ptr     dsplens;
  Int4        sumlens = 0;
  Int4        sumstart = 0;
  Int2        numseg = 0;
  Int2        index;
  Boolean     seen = FALSE;
  Boolean     seq_end;

  dsp = (CompSegPtr) salp->segs;
  if (dsp == NULL) {
         return -1;
  }
  index = position_inIdlist (sip, dsp->ids);
  if (index == 0) {
         return -1;
  }
  index -= 1;
  numseg = 0;
  for (dsplens =dsp->lens; dsplens !=NULL && numseg <dsp->numseg; dsplens++) {
         sumlens += *dsplens;
         numseg++; 
  }
  seq_end = (Boolean) (position >= sumlens);
  if (seq_end) return -1;
  sumlens = 0;
  numseg = 0;
  dspstart = dsp->starts + index;
  dsplens = dsp->lens;
  if (dspstart == NULL || dsplens == NULL ) {
         return -1;
  }
  while ( !seen && numseg < dsp->numseg ) {
         numseg++;
         if (position >=sumlens && position <sumlens +*dsplens ) {
                seen = TRUE;
         }
         else if ( numseg == dsp->numseg ) 
         {
                if ( salp->next == NULL ) break; 
                else 
                { 
                       sumstart += sumlens + *dsplens;
                       salp = salp->next;
                       dsp = (CompSegPtr) salp->segs;
                       dspstart = dsp->starts + index;
                       dsplens = dsp->lens;
                       numseg = 0;
                }
         }
         else 
         {
                sumlens += *dsplens;
                dspstart += dsp->dim; 
                dsplens++;
         }
  }
  if ( !seen ) return -1;
  if ((Boolean) *dspstart) return position;
  if (numseg == dsp->numseg) return -1;
  while (numseg < dsp->numseg) {
         numseg++;
         sumlens += *dsplens;
         dspstart += dsp->dim; 
         if ((Boolean) *dspstart) break;
         dsplens++;
  }
  return sumlens;
}

static Int4 getlastpos_inseqalign (SeqAlignPtr salp, Int4 position, SeqIdPtr sip)
{
  CompSegPtr  dsp;
  BoolPtr     dspstart;
  Int4Ptr     dsplens;
  Int4        sumlens = 0;
  Int4        sumstart = 0;
  Int4        lastpos;
  Int2        numseg = 0;
  Int2        index;
  Boolean     seen = FALSE;
  Boolean     seq_end;

  dsp = (CompSegPtr) salp->segs;
  if (dsp == NULL) {
         return -1;
  }
  index = position_inIdlist (sip, dsp->ids);
  if (index == 0) {
         return -1;
  }
  index -= 1;
  numseg = 0;
  for (dsplens =dsp->lens; dsplens !=NULL && numseg < dsp->numseg; dsplens++) {
         sumlens += *dsplens;
         numseg++;
  }
  seq_end = (Boolean) (position >= sumlens);
  if (position < 0 || seq_end) 
     return sumlens-1;
  sumlens = 0;
  numseg = 0;
  dspstart = dsp->starts + index;
  dsplens = dsp->lens;
  if (dspstart == NULL || dsplens == NULL ) {
         return -1;
  }
  lastpos = 0;
  while ( !seen && numseg < dsp->numseg ) {
         numseg++;
         if ((Boolean) *dspstart) {
                lastpos = sumlens + *dsplens - 1;
         }
         if (position >=sumlens && position <sumlens +*dsplens ) {
                seen = TRUE;
         }
         else if ( numseg == dsp->numseg ) 
         {
                if ( salp->next == NULL ) break; 
                else 
                { 
                       sumstart += sumlens + *dsplens;
                       salp = salp->next;
                       dsp = (CompSegPtr) salp->segs;
                       dspstart = dsp->starts + index;
                       dsplens = dsp->lens;
                       numseg = 0;
                }
         }
         else 
         {
                sumlens += *dsplens;
                dspstart += dsp->dim; 
                dsplens++;
         }
  }
  if ( !seen ) return -1;
  if ((Boolean) *dspstart) return position;
  return lastpos;
}

/*********************************************************
***
***  SeqAlignPtr procedures
***
**********************************************************/
extern SeqAlignPtr SeqAlignDup (SeqAlignPtr salp)
{
  SeqAnnotPtr sap, 
              sap2;
  SeqAlignPtr salp2 = NULL,
              next; 

  next = salp->next;
  salp->next = NULL;
  sap = SeqAnnotForSeqAlign (salp);
  sap2 = (SeqAnnotPtr) AsnIoMemCopy ((Pointer) sap, (AsnReadFunc) SeqAnnotAsnRead, (AsnWriteFunc) SeqAnnotAsnWrite);
  if (sap2!=NULL) {
     salp2 = sap2->data;
     sap2->data = NULL;
     SeqAnnotFree (sap2);
  }
  if (next != NULL)
     salp->next = next;
  return salp2;
}

extern SeqAlignPtr SeqAlignDupRegion (SeqAlignPtr salp, Int2 to_numseg, Int4 subseg, Boolean first_part)
{
  SeqAlignPtr newsalp = NULL;
  DenseSegPtr dsp = NULL, newdsp = NULL;
  Int4Ptr     dspstart = NULL;
  Int4Ptr     dsplen = NULL;
  Int4Ptr     newstart = NULL;
  Int4Ptr     newlens = NULL;
  Int2        newseg;
  Int2        n;
  Int2        k;
  Int4        j;

  dsp = (DenseSegPtr) salp->segs;
  newsalp = SeqAlignNew ();
  if ( newsalp == NULL ) return NULL;
  newsalp->type = salp->type;
  newsalp->segtype = salp->segtype;
  newsalp->dim = salp->dim;
  newdsp = DenseSegNew ();
  if ( newdsp == NULL ) return NULL;
  newsalp->segs = (Pointer) newdsp;
  n = newdsp->dim = dsp->dim;
  newdsp->ids = SeqIdDupList (dsp->ids);
  if ( to_numseg > dsp->numseg ) to_numseg = 0;
  if ( to_numseg == 0 )
         newseg = dsp->numseg;
  else if ( first_part )
         newseg = to_numseg;
  else if ( !first_part )
         newseg = dsp->numseg - to_numseg;
  newdsp->numseg = newseg;
  newdsp->starts = (Int4Ptr) MemNew ((size_t) ((n*newseg+4) * sizeof (Int4)));
  if ( newdsp->starts == NULL ) return NULL;
  for (j = 0; j < n * newseg + 4; j++) newdsp->starts [j] = -1;
  newdsp->lens  = (Int4Ptr) MemNew ((size_t) ((newseg + 2) * sizeof (Int4))); 
  if ( newdsp->lens == NULL ) return NULL;
  for (j = 0; j < newseg + 2; j++)     newdsp->lens [j] = 0;
  newstart = newdsp->starts;
  newlens  = newdsp->lens;
  dspstart = dsp->starts;
  dsplen   = dsp->lens;
  if ( first_part )
  {
         for ( k = 0; k < newseg - 1; k++ ) 
         {
                for (j = 0; j < n; j++) newstart [j] = dspstart[j];
                *newlens = *dsplen;
                newstart += n;
                newlens++;
                dspstart += n;
                dsplen++;
         }
         for (j = 0; j < n; j++) newstart [j] = dspstart[j];
         *newlens = subseg;
  }
  else {
         for ( k = 0; k < to_numseg-1; k++ ) 
         {
                dspstart += n;
                dsplen++;
         } 
         for (j = 0; j < n; j++) 
                       if ( dspstart[j] > -1 )
                              newstart [j] = dspstart[j] + subseg;
                       else 
                              newstart [j] = dspstart[j];
         *newlens = *dsplen - subseg;
         newstart += n;
         newlens++;
         dspstart += n;
         dsplen++;
         k++;
         for ( ; k < to_numseg + newseg; k++ ) 
         {
                for (j = 0; j < n; j++) newstart [j] = dspstart[j];
                *newlens = *dsplen;
                newstart += n;
                newlens++;
                dspstart += n;
                dsplen++;
         }
  }
  return newsalp;
}

extern SeqAlignPtr SeqAlignAdd (SeqAlignPtr *salp_head, SeqAlignPtr salp)
{
  SeqAlignPtr salp_tmp;

  if (salp==NULL)
     return *salp_head;
  salp_tmp = *salp_head;
  if (salp_tmp != NULL ) 
  {
         while ( salp_tmp->next != NULL ) 
                salp_tmp = salp_tmp->next;
         salp_tmp->next = salp;
  }
  else *salp_head = salp;
  return *salp_head;
}

extern SeqAlignPtr SeqAlignDupAdd (SeqAlignPtr *salp_head, SeqAlignPtr salp, Int2 to_numseg, Int4 subseg, Boolean first_part)
{
  SeqAlignPtr salp_tmp, salp_copy;

  salp_copy = SeqAlignDupRegion (salp, to_numseg, subseg, first_part);
  if ( salp_copy == NULL )
         return *salp_head;
  if ( (salp_tmp = *salp_head) != NULL ) 
  {
         while ( salp_tmp->next != NULL ) 
                salp_tmp = salp_tmp->next;
         salp_tmp->next = salp_copy;
  }
  else *salp_head = salp_copy;
  return *salp_head;
}

extern SeqAlignPtr SeqAlignEndExtend (SeqAlignPtr sap, Int4 start1, Int4 start2, Int4 stop1, Int4 stop2, Int4 x1, Int4 y1, Int4 x2, Int4 y2, Uint1 strand1, Uint1 strand2)
{
   DenseDiagPtr  ddp;
   DenseSegPtr   dsp;
   Int4Ptr       n_starts, n_lens;
   Uint1Ptr      n_strands;
   Int4          index1,
                 minlen = 0;
   Int2          n;
   Int2          offset=0;
   Boolean       is_strands = FALSE;

   if (sap==NULL)
      return NULL;
   if (start1==x1 && start2==y1 && stop1==x2 && stop2==y2)
      return sap;
   if (sap->segtype == 1) {
      ddp = (DenseDiagPtr) sap->segs;
   } 
   else if (sap->segtype == 2)
    {
      dsp = (DenseSegPtr) sap->segs;
      is_strands = (Boolean) (dsp->strands!=NULL);
      n = dsp->numseg;
      n_starts = (Int4Ptr) MemNew((2*(n+2)+4)*sizeof(Int4));
      n_lens = (Int4Ptr) MemNew((n+4)*sizeof(Int4));
      if (is_strands)
         n_strands = MemNew((2*(n+2)+4)*sizeof(Int1));
      if (x1 > start1 && y1 > start2) { 
         minlen = MIN ((x1-start1), (y1-start2)); 
         n = 0;
      }
      if ((x1 > start1 || y1 > start2)  && ((x1-start1) != (y1-start2))) { 
         if (x1 > start1 && (x1-start1) > minlen) {
            n_starts[0] = 0;
            n_starts[1] = -1;
            n_lens[0] = (x1 - start1) - minlen;
            if (is_strands) {
               n_strands[0] = strand1;
               n_strands[1] = strand2; 
            }
         } 
         else if (y1 > start2 && (y1-start2) > minlen) {
            n_starts[0] = -1;
            n_starts[1] = 0;
            n_lens[0] = (y1 - start2) - minlen;
            if (is_strands) {
               n_strands[0] = strand1;
               n_strands[1] = strand2; 
            }
         }
         offset += 1;
      }
      if (minlen > 0) {
         n += offset;
         n_starts[2*n] = x1 - minlen;
         n_starts[2*n+1] = y1 - minlen;
         n_lens[n] = minlen;
         if (is_strands) {
            n_strands[2*n] = strand1;
            n_strands[2*n+1] = strand2;
         }
         offset += 1;
      }
      n = dsp->numseg;
      for (index1=0; index1<n; index1++) {
         n_starts [2*(index1+offset)] = dsp->starts [2*index1];
         if (is_strands)
            n_strands[2*(index1+offset)] = dsp->strands[2*index1];
         n_starts [2*(index1+offset)+1]=dsp->starts [2*index1+1];
         if (is_strands)
            n_strands[2*(index1+offset)+1]=dsp->strands[2*index1+1];
         n_lens[index1+offset] = dsp->lens[index1];
      }
      if (x2 < stop1 &&  y2 < stop2) {
         n += offset;
         minlen = MIN ((stop1-x2), (stop2-y2));
         n_starts[2*n] = x2 + 1;
         n_starts[2*n+1] = y2 +1;
         n_lens[n] = minlen;
         if (is_strands) {
            n_strands[2*n] = strand1;
            n_strands[2*n+1] = strand2;
         }
         x2 += minlen;
         y2 += minlen;
         offset += 1;
      }
      if (x2 < stop1 || y2 < stop2) { 
         n += offset;
         if (x2 < stop1) {
            n_starts[2*n] = x2 +1;
            n_starts[2*n+1] = -1;
            n_lens[n] = stop1 - x2;
            if (is_strands) {
               n_strands[2*n] = strand1;
               n_strands[2*n+1] = strand2;
            }
         } 
         else if (y2 < stop2) {
            n_starts[2*n] = -1;
            n_starts[2*n+1] = y2 +1;
            n_lens[n] = stop2 - y2;
            if (is_strands) {
               n_strands[2*n] = strand1;
               n_strands[2*n+1] = strand2;
            }
         }
         offset += 1;
      }
      dsp->numseg = n+1;
      MemFree(dsp->starts);
      if (is_strands)
         MemFree(dsp->strands);
      MemFree(dsp->lens);
      dsp->starts = n_starts;
      dsp->lens = n_lens;
      if (is_strands)
         dsp->strands = n_strands;
   }
   return sap;
}


extern SeqAlignPtr SeqAlignTrunc (SeqAlignPtr salp, Int4 from, Int4 to)
{
  DenseSegPtr dsp;
  Int4Ptr     int4p;
  Int2        j;

  if (salp == NULL)
     return NULL;
  if (salp->segtype == 2) {
     dsp = (DenseSegPtr) salp->segs; 
     if (from != 0) {
        int4p = dsp->starts;
        for (j=0; j<dsp->dim; j++, int4p++)
           if (*int4p > -1)
              *int4p += from; 
        int4p = dsp->lens;
        *int4p -= from;
     }
     if (to != 0) {
        int4p = dsp->lens;
        for (j=0; j<(dsp->numseg-1); j++)
           int4p++;
        *int4p += to;
     }
  }
  return salp;
}

extern SeqAlignPtr SeqAlignMap (SeqAlignPtr salp)
{
  SeqAlignPtr tmp;
  DenseSegPtr dsp;
  Int4Ptr     startp;
  Int4Ptr     lenp;
  Int2        j;
  
  if (salp == NULL)
     return NULL;
  for (tmp=salp; tmp!=NULL; tmp=tmp->next) {
     if (tmp->segtype == 2) {
        dsp = (DenseSegPtr) tmp->segs;
        for (j=0, startp = dsp->starts; j<(dsp->numseg-1); j++)
           startp += dsp->dim;
        while (*startp < 0) {
           startp -= dsp->dim;
           dsp->numseg--;
        }
        startp = dsp->starts;
        lenp = dsp->lens;
        while (*startp < 0) {
           startp += dsp->dim;
           lenp++;
           dsp->numseg--;
        }
        dsp->starts = startp;
        dsp->lens = lenp;
     }     
  }
  return salp;
}

extern SeqAlignPtr SeqAlignMerge (SeqAlignPtr salp1, SeqAlignPtr salp2, Boolean return_salp)
{
  SeqAlignPtr salp_toreturn;
  DenseSegPtr dsp1= NULL, 
              dsp2 = NULL;
  Int4Ptr     dspstarts= NULL;
  Int4Ptr     dsptmp= NULL;
  Int4Ptr     dsplens= NULL;
  Uint1Ptr    dspstrands;
  Uint1Ptr    dsptmp1;
  Int2        j, k, n, newseg;
  Uint1       st1, st2;
  Boolean     merge_segment = FALSE;

  if (salp1==NULL) {
     if (salp2 == NULL)
        return NULL;
     else
        return salp2;
  }
  else if (salp2 == NULL)
     return salp1;
  
  if (return_salp)
     salp_toreturn = salp2;
  else 
     salp_toreturn = salp1;

  if ( salp1->segtype != 2 || salp2->segtype != 2 ) {
     return salp_toreturn; 
  }
  if (return_salp) {
     dsp1 = (DenseSegPtr) salp1->segs;
     dsp2 = (DenseSegPtr) salp2->segs;
  }
  else {
     dsp1 = (DenseSegPtr) salp2->segs;
     dsp2 = (DenseSegPtr) salp1->segs;
  }
  if ( dsp1==NULL || dsp2==NULL || dsp1->dim != dsp2->dim) {
     return salp_toreturn; 
  }
  n = dsp1->dim;
  newseg = dsp1->numseg + dsp2->numseg; 
  dspstarts = (Int4Ptr) MemNew ((size_t) ((n*newseg+4) * sizeof (Int4)));
  if ( dspstarts == NULL ) {
     return salp_toreturn; 
  }
  st1 = SeqAlignStrand (salp1, 0);
  st2 = SeqAlignStrand (salp1, 1);
  dsptmp = dsp1->starts;
  for (j = 0; j < n * dsp1->numseg; j++, dsptmp++) {
     dspstarts [j] = *dsptmp;
  }
  dsptmp = dsp2->starts;
  if (n==2) {
     if (dspstarts [j-2]> -1 && dsptmp[0] > -1) {
        if (dspstarts [j-1] > -1 && dsptmp[1] > -1) 
           merge_segment = TRUE;
     }    
  }
  if (merge_segment) {
     if (st1==Seq_strand_minus)
        dspstarts [j-2] = dsptmp[0];
     if (st2==Seq_strand_minus)
        dspstarts [j-1] = dsptmp[1];
     newseg--;
     k=n;
     dsptmp += n;
  }
  else 
     k = 0;
  for (; k < n * dsp2->numseg; k++, j++, dsptmp++) {
     dspstarts [j] = *dsptmp;
  }
  dsplens  = (Int4Ptr) MemNew ((size_t) ((newseg + 2) * sizeof (Int4))); 
  if ( dsplens == NULL ) {
     return salp_toreturn; 
  }
  dsptmp = dsp1->lens;
  for (j = 0; j < dsp1->numseg; j++, dsptmp++) 
     dsplens [j] = *dsptmp;
  dsptmp = dsp2->lens;
  if (merge_segment) {
     dsplens [j-1] += *dsptmp;
     k=1;
     dsptmp ++;
  } 
  else
     k = 0; 
  for (; k < dsp2->numseg; k++, j++, dsptmp++) 
     dsplens [j] = *dsptmp;
  dspstrands = (Uint1Ptr) MemNew ((size_t) ((n*newseg+4) * sizeof (Uint1)));
  if ( dspstrands == NULL ) {
     return salp_toreturn;
  }
  dsptmp1=dsp1->strands;
  for (j = 0; j < n * dsp1->numseg; j++, dsptmp1++) {
     dspstrands [j] = *dsptmp1;
  }
  dsptmp1 = dsp2->strands; 
  if (merge_segment) {
     k=n; 
     dsptmp1 += n;
  } 
  else 
     k = 0; 
  for (; k < n * dsp2->numseg; k++, j++, dsptmp1++) {
     dspstrands [j] = *dsptmp1;
  }
  dsp1 = (DenseSegPtr) salp1->segs;
  dsp1->numseg = newseg;
  MemFree(dsp1->starts);
  dsp1->starts = dspstarts;
  MemFree(dsp1->lens);
  dsp1->lens = dsplens;
  MemFree(dsp1->strands);
  dsp1->strands = dspstrands; 
  return salp1;
}

extern SeqAnnotPtr SeqAnnotMerge (SeqAnnotPtr sap1, SeqAnnotPtr sap2, Boolean return_salp)
{
  SeqAlignPtr salp1=NULL, salp2=NULL;

  if ( sap1 == NULL ) return sap2;
  if ( sap2 == NULL ) return sap1;
  if ( sap1->type != 2 ||  sap2->type != 2 ) {
         return NULL; 
  }
  salp1 = (SeqAlignPtr) sap1->data;
  salp2 = (SeqAlignPtr) sap2->data;
  if (return_salp) 
     salp1 = SeqAlignMerge (salp1, salp2, TRUE);
  else 
     salp1 = SeqAlignMerge (salp1, salp2, FALSE);
  return sap1;
}

static SeqAlignPtr SeqAlignAddSeg (SeqAlignPtr salp, Int4 pos1, Int4 pos2, Int4 len)
{
  DenseSegPtr dsp;
  Int4Ptr     startp;
  Int4Ptr     lenp;
  Uint1Ptr    dspstrands,
              dsptmp1;
  Int4Ptr     dsptmp;
  Int2        j;

  dsp = (DenseSegPtr) salp->segs;

  startp = (Int4Ptr) MemNew ((size_t) ((dsp->dim*dsp->numseg+4) * sizeof (Int4)));
  if ( startp == NULL ) {
     return NULL;
  }
  dsptmp = dsp->starts;
  for (j = 0; j < dsp->dim*dsp->numseg; j++, dsptmp++) {
     startp [j] = *dsptmp;
  }
  startp[j] = pos1;
  startp[j+1] = pos2;
  MemFree(dsp->starts);
  dsp->starts = startp;

  lenp  = (Int4Ptr) MemNew ((size_t) ((dsp->numseg + 2) * sizeof (Int4)));
  if ( lenp == NULL ) {
     return NULL;
  }
  dsptmp = dsp->lens;
  for (j = 0; j < dsp->numseg; j++, dsptmp++)
     lenp [j] = *dsptmp;
  lenp [j] = len;
  MemFree(dsp->lens);
  dsp->lens = lenp;

  dspstrands = (Uint1Ptr) MemNew ((size_t) ((dsp->dim*dsp->numseg+4) * sizeof (Uint1)));
  if ( dspstrands == NULL ) {
     return NULL;
  }
  dsptmp1=dsp->strands;
  for (j = 0; j < dsp->dim * dsp->numseg; j++, dsptmp1++) {
     dspstrands [j] = *dsptmp1;
  }
  dspstrands [j] = Seq_strand_unknown;
  dspstrands [j+1] = Seq_strand_unknown;
  MemFree(dsp->strands);
  dsp->strands = dspstrands;

  dsp->numseg += 1;
  return salp; 
}

static SeqAlignPtr shrinksap5 (SeqAlignPtr salp, Int4 offset)
{
  DenseSegPtr dsp;
  Int4Ptr     dsptmp;
  Int2        j;

  dsp = (DenseSegPtr) salp->segs;
  dsptmp = dsp->starts;
  for (j = 0; j < dsp->dim; j++)
     dsptmp[j] += offset; 
  dsptmp = dsp->lens;
  *dsptmp -= offset;
  return salp;
}

static SeqAlignPtr shrinksap3 (SeqAlignPtr salp, Int4 offset)
{
  DenseSegPtr dsp;
  Int4Ptr     dsptmp;

  dsp = (DenseSegPtr) salp->segs;
  dsptmp = dsp->lens;
  *dsptmp -= offset;
  return salp;
}

extern SeqAlignPtr SeqAlignExtend (SeqAlignPtr salp1, SeqAlignPtr salp2)
{
  SeqAlignPtr salptmp;
  SeqLocPtr   slp1, slp2,
              slp1b, slp2b;
  ValNodePtr  vnp1, vnp2,
              slpt1, slpt2;
  Int4        gaplen;
  Int4        gaplen1, stop, stopb;
  Int2        j;
  Uint1       choice;
  Boolean     goOn=TRUE;

  slpt1 = SeqLocListFromSeqAlign (salp1);
  while (goOn) {
     goOn = FALSE;
     for (salptmp = salp2; salptmp!=NULL; salptmp=salptmp->next) {
        slpt2 = SeqLocListFromSeqAlign (salptmp);   
        choice = 0;
        for (vnp1=slpt1, vnp2=slpt2, j=0; vnp1!=NULL && vnp2!=NULL; j++) {
           slp1=(SeqLocPtr)vnp1->data.ptrvalue;   
           slp2=(SeqLocPtr)vnp2->data.ptrvalue;   
           gaplen1 = SeqLocStart(slp2) - SeqLocStop(slp1);
           if (gaplen1 ==0) {
              if (j==0) choice =1;
              else choice =2;
              break;
           }
           gaplen1 = SeqLocStart(slp1) - SeqLocStop(slp2);
           if (gaplen1 ==0 || gaplen1==1 || gaplen1==-1 || gaplen1==-2) {
              if (j==0) choice = 3;
              else choice =4;
              break;
           }
           vnp1=vnp1->next;
           vnp2=vnp2->next;
        } 
        if (choice > 0) {
           if (choice==1 || choice==3) {
              vnp1=vnp1->next;
              vnp2=vnp2->next;
           }
           else {
              vnp1=slpt1;
              vnp2=slpt2;
           }
           slp1b=(SeqLocPtr)vnp1->data.ptrvalue;
           slp2b=(SeqLocPtr)vnp2->data.ptrvalue;
           if (choice ==1 || choice == 2) {
              gaplen = SeqLocStart(slp2b) - SeqLocStop(slp1b) -1;   
              if (gaplen >= 1) {
                 stop = SeqLocStop(slp1);
                 stopb = SeqLocStop(slp1b);
                 if (gaplen1==0) {
                    shrinksap5 (salptmp, 1);
                    gaplen-=1;
                 }
                 if (gaplen > 1) {
                    if (choice ==1) {
                       SeqAlignAddSeg(salp1, -1, stopb+1, gaplen);
                    }
                    else {
                       SeqAlignAddSeg(salp1,stop+1, -1, gaplen);
                    }
                 }      
                 salp1 = SeqAlignMerge (salp1, salptmp, TRUE);
                 goOn=TRUE;
              } 
           } else {
              gaplen = SeqLocStart(slp1b) - SeqLocStop(slp2b) -1;
              if (gaplen >= 1) { 
                 stop = SeqLocStop(slp2);
                 stopb = SeqLocStop(slp2b);
                 if (gaplen1==0 || gaplen1==-1 || gaplen1==-2) {
                    gaplen1=ABS(gaplen1)+1;
                    shrinksap3 (salptmp, gaplen1);
                    stop-=gaplen1;
                    stopb-=gaplen1;
                    gaplen+=gaplen1;
                 }
                 if (gaplen > 1) {      
                    if (choice ==3) {
                       SeqAlignAddSeg(salptmp,-1,stop+1, gaplen); 
                    }
                    else {
                       SeqAlignAddSeg(salptmp, stopb+1,-1, gaplen); 
                    }
                 }
                 salp1 = SeqAlignMerge (salp1, salptmp, FALSE);
                 goOn=TRUE;
              }
           }   
           if (goOn) {
              ValNodeFreeType (&slpt1, TypeSeqLoc);
              slpt1 = SeqLocListFromSeqAlign (salp1); 
              break;
           }
        } 
        ValNodeFreeType (&slpt2, TypeSeqLoc);
     }   
  }
  ValNodeFreeType (&slp1, TypeSeqLoc);
  return salp1;
}

static Uint1 check_salplength (DenseSegPtr dsp)
{
  BioseqPtr   bsp;
  SeqIdPtr    sip;
  Int4Ptr     lenp;
  Int4Ptr     startp;
  Uint1Ptr    strandp;
  Uint1       strand;
  Int4        bsplength,
              sumlens,
              from,
              startpre, lenpre;
  Int2        k, j;
  Uint1       modif = 0;

  sip = dsp->ids;
  strandp = dsp->strands;
  for (k = 0; k < dsp->dim && sip!=NULL; k++, sip = sip->next)
  {
     bsp = BioseqLockById (sip);
     if (bsp!=NULL) 
     {
        bsplength = bsp->length;
        BioseqUnlock(bsp); 
        startp = dsp->starts;
        startp += k;
        lenp = dsp->lens;
        sumlens = 0;
        from = -1;
        if (strandp!=NULL)
           strand = *strandp;
        else 
           strand = Seq_strand_unknown;
        if (strand == Seq_strand_minus) {
           for (j = dsp->numseg-1; j >= 0; j--) {
              if (startp[dsp->dim * j] >= 0) {
                 if (from < 0) from = startp[dsp->dim * j];
                 sumlens += lenp[j];    
                 if (from + sumlens > bsplength) {
                    if (j==0)
                       return 2; 
                    startp[j] = -1;
                    modif = 1;     
                 }
              }   
           }
        }
        else {
           startpre = lenpre = -1;
           for (j = 0; j < dsp->numseg; j++) {
              if (*startp >= 0) {
                 if (from < 0) 
                    from = *startp;
                 sumlens += *lenp;
                 if (from + sumlens > bsplength || (startpre>=0 && *startp != startpre+lenpre)) {
                    if (j==0)
                       return 2;
                    *startp = -1;
                    modif = 1;
                 }
                 startpre = *startp;
                 lenpre = *lenp;
              }
              startp += dsp->dim;
              lenp++;
           }
        }
     }
     if (strandp!=NULL) strandp++;
  }
  return modif;
}

static Uint1 check_salpcomplength(CompSegPtr csp)
{
  BioseqPtr   bsp;
  SeqIdPtr    sip;
  Int4Ptr     lenp,
              fromp;
  BoolPtr     startp;
  Uint1Ptr    strandp;
  Uint1       strand;
  Int4        bsplength,
              sumlens;
  Int2        k, j;
  Uint1       modif = 0;

  strandp = csp->strands;
  sip = csp->ids;
  for (k = 0; k < csp->dim && sip!=NULL; k++, sip = sip->next)
  {
     bsp = BioseqLockById (sip);
     if (bsp!=NULL)
     {
        bsplength = bsp->length;
        BioseqUnlock(bsp);

        startp = csp->starts;
        startp += k;
        fromp = csp->from;
        fromp += k;
        sumlens = 0;
        lenp = csp->lens;  
        if (strandp!=NULL)
           strand = *strandp;
        else 
           strand = Seq_strand_unknown;
        if (strand==Seq_strand_minus) {
           for (j = 0; j < csp->numseg; j++) {
              if (*startp) {
                 sumlens += *lenp;
                 if (*fromp - sumlens > bsplength) {
                    if (j==csp->numseg-1) 
                       return 2;
                    *startp = FALSE;
                    modif = 1;
                 }
              }
              startp += csp->dim;
              lenp++;
           }
        }
        else {
           for (j = 0; j < csp->numseg; j++) {
              if (*startp) {
                 sumlens += *lenp;
                 if (*fromp + sumlens > bsplength) {
                    if (j==0) 
                       return 2;
                    *startp = FALSE;
                    modif = 1;

                 }
              }
              startp += csp->dim;
              lenp++;
           }
        }
        if (strandp!=NULL) strandp++;
     }
  }
  return modif;
}

extern SeqAlignPtr check_salp_forlength (SeqAlignPtr salp)
{
/***********************************************
  SeqAlignPtr salptmp, 
              presalp=NULL,
              salp1, out;
  Uint1       modif = 0;

  salp1 = salptmp = salp;
  while (salptmp!=NULL) {
     if (salptmp->segtype == 4) 
     {   
        modif = check_salpcomplength((CompSegPtr) salptmp->segs);
     }
     else if (salptmp->segtype == 2) 
     {  
        modif = check_salplength((DenseSegPtr) salptmp->segs);
     }
     if (modif==2) {
        out = salptmp;
        if (presalp!=NULL) {
           presalp->next = salptmp->next;
        } 
        else salp1 = salptmp->next; 
        salptmp = salptmp->next; 
        out->next = NULL;
        if (out->segtype == 4) 
           CompSeqAlignFree (out); 
        else if (out->segtype == 2)
           SeqAlignFree (out);
     }
     else {
        presalp = salptmp;
        salptmp=salptmp->next;
     }
  }
*******************************/
  return salp;
}

extern SeqAlignPtr check_salp_forstrand (SeqAlignPtr salp)
{
  SeqAlignPtr salptmp;
  DenseSegPtr dsp;
  CompSegPtr  csp;
  Int4Ptr     lenp;
  Int4Ptr     startp;
  Uint1Ptr    strandp = NULL;
  Int4        numseg;
  Int2        dim;
  Int4        j, k, n, tmp;
  Boolean     retourne = FALSE;

  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next) 
  {
     if (salptmp->segtype == 4) 
     {   
        csp=(CompSegPtr) salptmp->segs;
        strandp = csp->strands;
        numseg = csp->numseg;
        dim = csp->dim;
        lenp = csp->lens;
     }
     else if (salptmp->segtype == 2) 
     {  
        dsp = (DenseSegPtr) salptmp->segs;
        strandp = dsp->strands;
        numseg = dsp->numseg;
        dim = dsp->dim;
        lenp = dsp->lens;
        startp = dsp->starts;
     }
     if (strandp!=NULL) {
        if (*strandp != Seq_strand_plus && *strandp!=Seq_strand_minus)
           for (j=0; j<numseg; j++, strandp+=dim) {
              if (*strandp == Seq_strand_plus || *strandp==Seq_strand_minus)
                 break;
           }
        if (strandp!=NULL) 
           retourne = (Boolean) (*strandp == Seq_strand_minus);
     }
     if (retourne) {
        for (j=0; j < numseg*dim && strandp!=NULL; j++, strandp++) 
        {
           if (*strandp == Seq_strand_minus) 
              *strandp = Seq_strand_plus;
           else if (*strandp == Seq_strand_plus)
              *strandp = Seq_strand_minus;
        }
        for (j=0, k=numseg-1; j<numseg/2; j++, k--) {
           tmp=lenp[j];
           lenp[j]=lenp[k];
           lenp[k]=tmp;
        }
        for (j=0, k=(dim*numseg-dim); j<(dim*numseg-1)/2; j+=dim, k-=dim) {
           for (n=0; n<dim; n++) {
              tmp=startp[j+n];
              startp[j+n]=startp[k+n];
              startp[k+n]=tmp;
           }
        }
     }
  }
  return salp;
}

/*************************************************
***  
***      LocalAlignToSeqAnnotDimn
***              
*************************************************/
extern SeqAnnotPtr LocalAlignToSeqAnnotDimn (ValNodePtr seqvnp, SeqIdPtr seqsip, ValNodePtr fromp, Int2 nbseq, Int4 lens, ValNodePtr strands, Boolean trunc_emptyends)
{
  SeqAnnotPtr  sap;
  SeqAlignPtr  salp;
  DenseSegPtr  dsp;
  ValNodePtr   vnp;
  CharPtr      seqstr;
  Boolean PNTR filter;
  Boolean PNTR begun;
  Int4Ptr      from;
  Int4Ptr      lgseq;
  Int4Ptr      start;
  Uint1Ptr     strandp;
  Int4         lenstmp, the_lens;
  Int4         j;
  Int2         seg, numseg;
  Int2         k;
  Boolean      nogap;

  for (k =0, vnp =seqvnp; k < nbseq && vnp !=NULL; k++, vnp =vnp->next)
  {
         seqstr = (CharPtr) vnp->data.ptrvalue;
         lenstmp = StringLen (seqstr);
         if (k==0) 
            the_lens = lenstmp;
         else if (lenstmp != the_lens) {
            ErrPostEx (SEV_ERROR, 0, 0, "Sequence alignment of different lengths"); 
            return NULL;
         }
  }
  if (lens > 0 && lens != the_lens) {
/*
     WriteLog ("Length problem in the sequence alignment %ld %ld\n", lens, the_lens); 
     if (lens > the_lens)
*/
  }
  lens = the_lens;

         /*****************************
         ** count number of segments **
         *****************************/
  filter= MemNew ((size_t) ((nbseq + 1) * sizeof (Boolean)));
  j = 0;
  for (k =0, vnp =seqvnp; k < nbseq && vnp !=NULL; k++, vnp =vnp->next) 
  {
         seqstr = (CharPtr) vnp->data.ptrvalue;
         filter [k] = (Boolean)( seqstr [j] != '-' );
  }
  numseg = 1;
  while ( j < lens ) 
  {
         seg = 0;
         for (k =0, vnp =seqvnp; k <nbseq && vnp !=NULL; k++, vnp =vnp->next) 
         {
                seqstr = (CharPtr) vnp->data.ptrvalue;
                nogap = (Boolean) ( seqstr [j] != '-' );
                if ( filter [k] != nogap ) {
                   seg++;
                   filter [k] = nogap;
                }
         } 
         if ( seg > 0 ) ++numseg;
         j++;
  }
         /********************************************
         ** allocate SeqAnnot + SeqAlign + DenseSeg  **
         *********************************************/
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;
  sap->type = 2;
  salp = SeqAlignNew ();
  if (salp == NULL) return NULL;
  salp->type = 3;
  salp->segtype = 2;
  salp->dim = nbseq;
  sap->data = (Pointer) salp;
  dsp = DenseSegNew ();
  salp->segs = (Pointer) dsp;
  dsp->dim = nbseq;
  dsp->ids = seqsip;
  dsp->numseg = numseg;
  dsp->starts = (Int4Ptr)MemNew((size_t)((nbseq *numseg + 4) * sizeof (Int4)));
  for (j = 0; j < nbseq *numseg + 4; j++) 
     dsp->starts[j] = -1;
  dsp->lens   = (Int4Ptr) MemNew ((size_t) ((numseg + 2) * sizeof (Int4))); 
  for (k = 0; k < numseg + 2; k++) 
     dsp->lens[k] = 0;
  dsp->strands = (Uint1Ptr)MemNew((size_t) ((numseg*nbseq+4)*sizeof (Uint1)));
  strandp = dsp->strands;
  for (j = 0; j < numseg*nbseq ; j++, strandp++) 
     *(strandp) = Seq_strand_unknown; 
  if (strands!=NULL) {
     strandp = dsp->strands;
     for (k=0; k<numseg; k++) {
        for (j = 0, vnp=strands; j < nbseq && vnp!=NULL; j++, vnp=vnp->next) {
           *strandp = (Uint1)vnp->data.intvalue;
           strandp++;
        }
     }
  }
  j = 0;
  for (k =0, vnp =seqvnp; k < nbseq && vnp !=NULL; k++, vnp =vnp->next)
  {
         seqstr = (CharPtr) vnp->data.ptrvalue;
         filter [k] = (Boolean)( seqstr [j] != '-' );
  }
  lenstmp= 0;
  numseg = 0;
  while ( j < lens )
  {
         seg = 0;
         for (k =0, vnp =seqvnp; k <nbseq && vnp !=NULL; k++, vnp =vnp->next)
         {
            seqstr = (CharPtr) vnp->data.ptrvalue;
            nogap = (Boolean) ( seqstr [j] != '-' );
            if ( filter [k] != nogap ) {
               seg++;
               filter [k] = nogap;
            }
         }
         if ( seg > 0 ) {
            dsp->lens[numseg] = lenstmp;
            ++numseg;
            lenstmp = 0;
         }
         lenstmp++;
         j++;
  }
  if (lenstmp > 0)
     dsp->lens[numseg] = lenstmp;
         /******************************
         ***  store the segments      **
         ******************************/
  lgseq = MemNew ((size_t) ((nbseq + 1) * sizeof (Int4)));
  from = MemNew ((size_t) ((nbseq + 1) * sizeof (Int4)));
  begun = MemNew ((size_t) ((nbseq + 1) * sizeof (Boolean)));
  if ( lgseq == NULL || from == NULL || begun == NULL )
     return NULL;
  for (k = 0; k < nbseq; k++) 
     lgseq[k] = 0;
  if (fromp == NULL)
     for (k = 0; k < nbseq; k++) from [k] = 0;
  else {
     for (k=0, vnp=fromp; k<nbseq && vnp!=NULL; vnp=vnp->next, k++)
     {
        from [k] = (Int4)vnp->data.intvalue; 
     }
  }
  for (k = 0; k < nbseq; k++) 
     begun[k] = FALSE;
  start = dsp->starts;
  strandp = dsp->strands;
  j = 0;
  for (k =0, vnp =seqvnp; k <nbseq && vnp !=NULL; k++, vnp =vnp->next) 
  {
     seqstr = (CharPtr) vnp->data.ptrvalue;
     filter [k] = (Boolean)( seqstr [j] != '-' );
     if ( filter [k] ) {
        start [k] = lgseq [k] + from [k];
        begun [k] = TRUE;
     }
  }
  j++;
  numseg = 0;
  while ( j < lens ) 
  {
         seg = 0;
         for (k =0, vnp =seqvnp; k <nbseq && vnp !=NULL; k++, vnp =vnp->next) 
         {
                seqstr = (CharPtr) vnp->data.ptrvalue;
                nogap = ( seqstr [j] != '-' );
                if  ( nogap && begun [k] ) 
                   ++lgseq [k]; 
                if ( filter [k] != nogap ) {
                   seg++;
                   filter[k] = nogap;
                   begun [k] = TRUE;
                }
         }   
         if ( seg > 0 ) {
                start += nbseq;
                for (k = 0; k < nbseq; k++) {
                   if ( filter[k] ) {
                      if (strandp[k]==Seq_strand_minus) {
                             start[k] = from[k] - lgseq[k];
                      } else {
                             start[k] = from[k] + lgseq[k];
                      }
                   }
                }
                ++numseg;
         }
         j++;
  }
  MemFree (filter);
  MemFree (lgseq);
  MemFree (from);
  MemFree (begun);
  if (trunc_emptyends && salp!=NULL) {
     nogap = TRUE;
     dsp = (DenseSegPtr) salp->segs;
     while (nogap && dsp->numseg>0) {
        start = dsp->starts;
        start += (dsp->dim * (dsp->numseg-1) );
        for (j=0; j < dsp->dim; j++, start++) {
           if (*start > -1) 
              break;
        }   
        if (j == dsp->dim)
           dsp->numseg --;
        else nogap = FALSE;
     }
  }
  for (k=0, vnp=seqvnp; k <nbseq && vnp!=NULL; k++, vnp=vnp->next) {
     if (dsp->strands[k] == Seq_strand_minus) {
        start = dsp->starts; 
        start += k;
        lgseq = dsp->lens;
        for (j=0; j<dsp->numseg; j++, lgseq++) {
           if (*start > -1) {
              *start -= *lgseq;
           }
           start += nbseq;
        }
     }
  }
  if (salp == NULL || dsp->numseg == 0) {
     sap = SeqAnnotFree (sap);
     return NULL;
  }
  salp = check_salp_forlength (salp);
  if (salp == NULL) {  
     sap = SeqAnnotFree (sap);
     return NULL; 
  } 
  return sap;
}

/*************************************************
***  Alignment
***      LocalAlignToSeqAnnotCompDimn
***              
*************************************************/
extern SeqAnnotPtr LocalAlignToSeqAnnotCompDimn (ValNodePtr seqvnp, SeqIdPtr seqsip, Int2 nbseq, Int4 lens)
{
  SeqAnnotPtr  sap;
  SeqAlignPtr  salp;
  CompSegPtr   csp;
  ValNodePtr   vnp;
  CharPtr      seqstr;
  BoolPtr      start;
  Boolean PNTR filter;
  Boolean PNTR begun;
  Int4Ptr      lgseq;
  Int2         seg, numseg;
  Int4         j;
  Int2         k;
  Boolean      nogap;

         /*****************************
         ** count number of segments **
         *****************************/
  filter= MemNew ((size_t) ((nbseq + 1) * sizeof (Boolean)));
  j = 0;
  nogap = FALSE;
  while ( ! nogap && j < lens-1 ) 
  {
         for (k =0, vnp =seqvnp; k < nbseq && vnp !=NULL; k++, vnp =vnp->next) 
         {
                seqstr = (CharPtr) vnp->data.ptrvalue;
                nogap = ( seqstr [j] != '-' );
                if  ( nogap ) break;
         }
         if  ( nogap ) break;
         j++;
  }
  for (k =0, vnp =seqvnp; k < nbseq && vnp !=NULL; k++, vnp =vnp->next) 
  {
         seqstr = (CharPtr) vnp->data.ptrvalue;
         filter [k] = ( seqstr [j] != '-' );
  }
  numseg = 1;
  while ( j < lens ) 
  {
         seg = 0;
         for (k =0, vnp =seqvnp; k <nbseq && vnp !=NULL; k++, vnp =vnp->next) 
         {
                seqstr = (CharPtr) vnp->data.ptrvalue;
                nogap = ( seqstr [j] != '-' );
                if ( filter [k] != nogap ) {
                       seg++;
                       filter [k] = nogap;
                }
         } 
         if ( seg > 0 ) ++numseg;
         j++;
  }
         /********************************************
         ** allocate SeqAnnot + SeqAlign + CompSeg  **
         *********************************************/
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;
  sap->type = 2;
  salp = SeqAlignNew ();
  if (salp == NULL) return NULL;
  salp->type = 3;
  salp->segtype = 4;
  salp->dim = nbseq;
  sap->data = (Pointer) salp;
  csp = (CompSegPtr) MemNew ( sizeof (CompSeg) );
  salp->segs = (Pointer) csp;
  csp->dim = nbseq;
  csp->ids = ValNodeSeqIdListDup (seqsip);
  csp->numseg = numseg;
  csp->from = (Int4Ptr) MemNew ((size_t) ((nbseq + 2) * sizeof (Int4)));
  for (j = 0; j < nbseq + 2; j++) csp->from [j] = 0;
  csp->starts =(BoolPtr)MemNew((size_t)((nbseq*numseg+ 4) * sizeof (Boolean)));
  for (j = 0; j < nbseq *numseg + 4; j++) csp->starts[j] = FALSE;
  csp->lens   = (Int4Ptr) MemNew ((size_t) ((numseg + 2) * sizeof (Int4))); 
  for (j = 0; j < numseg + 2; j++) csp->lens[j] = 0;

         /******************************
         ***  store the segments      **
         ******************************/
  lgseq = MemNew ((size_t) ((nbseq + 1) * sizeof (Int4)));
  begun = MemNew ((size_t) ((nbseq + 1) * sizeof (Boolean)));
  if ( lgseq == NULL || begun == NULL ) return NULL;
  for (k = 0; k < nbseq; k++) lgseq[k] = 0;
  for (k = 0; k < nbseq; k++) begun[k] = FALSE;
  start = csp->starts;
  j = 0;
  nogap = FALSE;
  while ( !nogap && j < lens-1 ) 
  {
         for (k =0, vnp =seqvnp; k <nbseq && vnp !=NULL; k++, vnp =vnp->next) 
         {
                seqstr = (CharPtr) vnp->data.ptrvalue;
                nogap = ( seqstr [j] != '-' );
                if  ( nogap ) break;
         }
         if  ( nogap ) break;
         j++;
  }
  for (k =0, vnp =seqvnp; k <nbseq && vnp !=NULL; k++, vnp =vnp->next) 
  {
         seqstr = (CharPtr) vnp->data.ptrvalue;
         filter[k] = ( seqstr [j] != '-' );
         if (filter [k] ) 
                begun [k] = start [k] = TRUE;
  }
  j++;
  numseg = 0;
  csp->lens[numseg]++;
  while ( j < lens ) 
  {
         seg = 0;
         for (k =0, vnp =seqvnp; k <nbseq && vnp !=NULL; k++, vnp =vnp->next) 
         {
                seqstr = (CharPtr) vnp->data.ptrvalue;
                nogap = ( seqstr [j] != '-' );
                if  ( nogap && begun [k] ) ++lgseq [k]; 
                if ( filter [k] != nogap ) {
                       seg++;
                       filter [k] = nogap;
                       begun [k] = TRUE;
                }
         }   
         if ( seg > 0 ) {
                start += nbseq;
                for (k = 0; k < nbseq; k++) {
                       start [k] = filter [k];
                }
                ++numseg;
         }
         csp->lens[numseg]++;
         j++;
  }
  MemFree (filter);
  MemFree (lgseq);
  MemFree (begun);
  return sap;
}

/*******************************************
***
***  SeqAlignIDDelete 
***
********************************************/
extern SeqAlignPtr SeqAlignIDDelete (SeqAlignPtr salphead, SeqIdPtr sip)
{
  SeqAlignPtr salp,
              presalp = NULL,
              nextsalp;
  SeqIdPtr    presip, nextsip, tmpsip;
  DenseSegPtr dsp;
  Int4Ptr     startp;
  Int4        j, k;
  Int2        index = 0;

  salp = salphead; 
  while (salp!= NULL) 
  {
     nextsalp = salp->next; 
     tmpsip = SeqAlignIDList (salp);
     if ((position_inIdlist(sip, tmpsip)) > 0) 
     {
        if (salp->segtype == 1 || salp->dim == 2) 
        {
           if (salphead==NULL)
              salphead = nextsalp;
           if (presalp!=NULL)
              presalp->next = nextsalp;
           salp->next = NULL;
           salp = SeqAlignFree (salp); 
           if (presalp!=NULL)
              salp = presalp;
        } 
        else if (salp->segtype == 2) 
        {
           dsp = (DenseSegPtr) salp->segs;
           presip = NULL;
           while (tmpsip!=NULL) 
           {
            nextsip = tmpsip->next;
            if (SeqIdForSameBioseq (tmpsip, sip)) 
            {
              if (presip==NULL) {
                 dsp->ids = nextsip;
              } else {
                 presip->next = nextsip;
              }
              tmpsip->next = NULL;
              SeqIdFree (tmpsip);
              startp = dsp->starts;
              k=0;
              for (j=0; j<dsp->dim*dsp->numseg; j++) {
                 if (((j-index) % (dsp->dim))!=0) {
                    startp[k] = startp[j];
                    k++;
                 }
              } 
              dsp->dim--;
              salp->dim--;
            }
            else 
              presip = tmpsip;
            tmpsip = nextsip;
            index++;
           } 
        }
     }
     presalp = salp;
     salp = nextsalp;
  }
  return salphead; 
}

extern Boolean SeqAlignIDCache (SeqAlignPtr salphead, SeqIdPtr sip)
{
  SeqAlignPtr salp;
  SeqIdPtr    tmpsip;
  Boolean     ok = FALSE;

  for (salp = salphead; salp!= NULL; salp=salp->next)
  {
     tmpsip = SeqAlignIDList (salp);
     if ((position_inIdlist(sip, tmpsip)) > 0) 
     {
        if (salp->segtype == 1 || salp->dim == 2) 
        {
           salp->type = 0;
           ok = TRUE;
        } 
        else if (salp->segtype == 2)
        {
           SeqAlignIDDelete (salp, sip); 
           ok = TRUE;
        }
     }
  }
  return ok; 
}

extern SeqAlignPtr SeqAlignIDUncache (SeqAlignPtr salphead, SeqIdPtr sip)
{
  SeqAlignPtr salp;
  SeqIdPtr    tmpsip;

  for (salp = salphead; salp!= NULL; salp=salp->next)
  {
     tmpsip = SeqAlignIDList (salp);
     if ((position_inIdlist(sip, tmpsip)) > 0) 
     {
        if (salp->segtype == 1 || salp->dim == 2) 
        {
           salp->type = 3;
        } 
     }
  }
  return salphead; 
}

extern SeqAlignPtr SeqAlignIDUncacheAll (SeqAlignPtr salphead)
{
  SeqAlignPtr salp;
  SeqIdPtr    tmpsip;

  for (salp = salphead; salp!= NULL; salp=salp->next)
  {
     if (salp->type < 1)
        if (salp->segtype == 1 || salp->dim == 2) 
           salp->type = 3;
  }
  return salphead; 
}

typedef struct ccid1 {
  SeqIdPtr       sip;
} CcId1, PNTR CcId1Ptr;

static void FindDeleteSeqAlignCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqAlignPtr        salp;
  CcId1Ptr           cip;

  cip = (CcId1Ptr)mydata;
  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           salp=is_salp_in_sap(bsp->annot, 2);
           SeqAlignIDDelete (salp, cip->sip);
        }
     }   
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           salp=is_salp_in_sap(bssp->annot, 2);
           SeqAlignIDDelete (salp, cip->sip);
        }
     }   
  }
}

extern void DelAlignItem (SeqEntryPtr sep, SeqIdPtr sip)
{
  CcId1            ci;

  ci.sip = SeqIdDup (sip);
  SeqEntryExplore (sep, (Pointer)&ci, FindDeleteSeqAlignCallback);
  if (ci.sip != NULL)
     SeqIdFree (ci.sip);
}

/**************************************************
***
***************************************************/
extern SeqAlignPtr SeqAlignDeleteByLoc (SeqLocPtr slp, SeqAlignPtr salp)
{
  SeqIdPtr    sip;
  DenseSegPtr dsp;
  Int4Ptr     dspstart;
  Int4Ptr     newstartp, newstart;
  Int4Ptr     dsplens;
  Int4Ptr     newlensp, newlens;
  Uint1Ptr    strandp;
  Uint1Ptr    newstrdp, newstrd;
  Int4        from;
  Int4        sumlens = 0;
  Int4        seqlens = 0;
  Int4        lensplus;
  Int4        position;
  Int2        newseg;
  Int2        j, tmp;
  Int2        numseg;
  Int2        inter_salp = 0;
  Int2        index;
  Boolean     seen = FALSE;
  Boolean     delete_before = FALSE;
  Int2    intersalpwidth=0;
  Int2    is_end=0;

  if (salp == NULL)
     return NULL;
  sip = SeqLocId(slp);
  dsp = (DenseSegPtr) salp->segs;
  if (dsp == NULL) {
         return salp;
  }
  index = position_inIdlist (sip, dsp->ids);
  if (index == 0) {
         return salp;
  }
  index -= 1;
  from = SeqAlignStart(salp, index);
/**/     
  delete_before = (Boolean) (SeqLocStart(slp) <= from);
  if (delete_before) 
     position = SeqLocStop(slp) +1;
  else
     position = SeqLocStart(slp);
/**/
  dspstart = dsp->starts + index;
  dsplens = dsp->lens;
  if (dspstart == NULL || dsplens == NULL ) {
         return salp;
  }
  numseg = 1;
  while ( !seen && numseg <= dsp->numseg ) 
  {
     if (position>=from+seqlens && position<from+seqlens+*dsplens)
     {
        if (*dspstart > -1)
           lensplus = (Int4)ABS(position - from - seqlens);
        seen = TRUE;
     }
     else if (*dspstart > -1 && position<=from+seqlens+*dsplens 
     && is_end==APPEND_RESIDUE) {
                lensplus = (Int4)ABS(position - from - seqlens);
                seen = TRUE;
     }
     else if ( numseg == dsp->numseg )
     {
                break;
    }
     else if (numseg < dsp->numseg)
     {
                sumlens += *dsplens;
                if (*dspstart > -1) 
                   seqlens += *dsplens;
                dspstart += dsp->dim;
                dsplens++;
     }
     if ( !seen ) numseg++;
  }       
  if ( !seen ) {
    return salp;
  }
  if (position != from+seqlens) 
  {
     newseg = dsp->numseg+1;
     newstart =(Int4Ptr)MemNew((size_t)((dsp->dim*newseg+4)*sizeof (Int4)));
     newlens  =(Int4Ptr)MemNew((size_t) ((newseg + 2) * sizeof (Int4)));
     if (dsp->strands!=NULL) 
        newstrd  =(Uint1Ptr)MemNew((size_t)((dsp->dim*newseg+4)*sizeof(Uint1)));
     if (newstart!=NULL && newlens!=NULL) 
     {
        if (dsp->strands!=NULL) 
        {
           strandp = dsp->strands;
           newstrdp= newstrd;
           for (tmp=0; tmp<dsp->dim*dsp->numseg; newstrdp++, tmp++) {
              *newstrdp = *strandp;
              strandp++;
           }
           strandp = dsp->strands;
           for (tmp=0; tmp<dsp->dim; newstrdp++, tmp++) {
              *newstrdp = *strandp;
              strandp++;
           }
        }
        dspstart = dsp->starts;
        dsplens = dsp->lens;
        newstartp = newstart;
        newlensp = newlens;
/*------*/
        if (delete_before) 
        {
           for(tmp=0; tmp<numseg-1; tmp++) {
              for (j = 0; j < dsp->dim; j++) {
                 if (j == index)
                    *newstartp = -1;
                 else
                    *newstartp = *dspstart;
                 newstartp++; 
                 dspstart++;
              }
              *newlensp = *dsplens;
              newlensp++; dsplens++;
           }        
        }
        else {
           for(tmp=0; tmp<numseg -1; tmp++) {
              for (j = 0; j < dsp->dim; j++) {
                 *newstartp = *dspstart;
                 newstartp++; 
                 dspstart++;
              }
              *newlensp = *dsplens;
              newlensp++; dsplens++;
           }
        }
/*------*/
        *newlensp = lensplus+intersalpwidth*inter_salp;
        for (j = 0; j < dsp->dim; j++) {
           if (delete_before)
           {
              if (j == index) {
                 *newstartp = -1;
                 *(newstartp + dsp->dim) = *dspstart + *newlensp;
              }
              else {
                 *newstartp = *dspstart;
                 if (*dspstart < 0)
                    *(newstartp + dsp->dim) = -1;
                 else
                    *(newstartp + dsp->dim) = *dspstart + *newlensp;
              }
           }
           else {           
              *newstartp = *dspstart;
              if (j == index || *dspstart < 0)
                 *(newstartp + dsp->dim) = -1;
              else 
                 *(newstartp + dsp->dim) = *dspstart + *newlensp;
           }
           newstartp++; 
           dspstart++;
        }
        newstartp += dsp->dim;
        newlensp++; 
        tmp++;
        *newlensp = *dsplens-(lensplus+intersalpwidth*inter_salp);
        newlensp++; 
        dsplens++;
/*------*/
        if (delete_before)
        {
           for(; tmp < dsp->numseg; tmp++) {
              for (j = 0; j < dsp->dim; j++) {
                 *newstartp = *dspstart;
                 newstartp++; 
                 dspstart++;
              }
              *newlensp = *dsplens;
              newlensp++; dsplens++;
           }
        }
        else 
        {
           for(; tmp < dsp->numseg; tmp++) {
              for (j = 0; j < dsp->dim; j++) {
                 if (j == index)
                    *newstartp = -1; 
                 else
                    *newstartp = *dspstart;
                 newstartp++; 
                 dspstart++;
              }
              *newlensp = *dsplens;
              newlensp++; dsplens++;
           }
        }
/*------*/
        dsp->numseg = newseg;
        dspstart = dsp->starts;
        dsp->starts = newstart;
        MemFree (dspstart);
        dsplens = dsp->lens;
        dsp->lens = newlens;
        MemFree (dsplens);
        strandp = dsp->strands;
        dsp->strands = newstrd;
        MemFree (strandp);
     }
     if (delete_before)
     {
        sumlens = from;
        dspstart = dsp->starts + index;
        dsplens = dsp->lens;
        for(tmp=0; tmp< dsp->numseg; tmp++) {
           if (*dspstart > -1) {
              *dspstart = sumlens;
              sumlens+=*dsplens;
           }
           dspstart += dsp->dim;
           dsplens++; 
        }
     }
  }
  else if (position == from+seqlens && numseg <= dsp->numseg) {
     if (delete_before) 
     {
        dspstart = dsp->starts + index;
        for(tmp=0; tmp< numseg-1; tmp++) {
           *dspstart = -1;
           dspstart += dsp->dim;
        }     
        /*------------------*/ 
        sumlens = from;
        dspstart = dsp->starts + index;
        dsplens = dsp->lens;
        for(tmp=0; tmp< dsp->numseg; tmp++) {
           if (*dspstart > -1) {
              *dspstart = sumlens;
              sumlens+=*dsplens;
           }
           dspstart += dsp->dim;
           dsplens++; 
        }
     }
     else {
        while (numseg <= dsp->numseg) {
           *dspstart = -1;
           numseg++;
           dspstart += dsp->dim;
        }
     }
  }
  return salp;
}

/*******************************************
***
***   DeleteRegion
***  !@!!!!!!!!!!!!!!!!!!!!!!! > CompSegPtr
********************************************/
extern SeqAlignPtr DeleteRegion (SeqIntPtr sip, SeqAlignPtr salp)
{
  CompSegPtr  dsp;
  SeqAlignPtr salp1 =  NULL;
  BoolPtr     dspstart = NULL;
  Int4Ptr     dsplens = NULL;
  Int4        delete_from;
  Int4        delete_to;
  Int2        numseg_before = 0;
  Int2        subseglens = 0;
  Int4        sumseglens = 0;
  Boolean     seen = FALSE;

  if ( sip == NULL ) return salp;
  delete_from = sip->from;
  delete_to = sip->to;

        /*****************************************
         *** copy salp(s) until delete_from
         *****************************************/
  if ( (dsp = (CompSegPtr) salp->segs ) == NULL) {
         return NULL;
  }
  dspstart = dsp->starts;
  dsplens = dsp->lens;
  while ( !seen )
  {
         if ( !(seen = locate_in_seqalign (delete_from, dsp->dim, dsp->numseg, &dspstart, &dsplens, &numseg_before, &subseglens, &sumseglens)) ) 
         {
                salp1 = SeqAlignDupAdd (&salp1, salp, 0, 0, 0);
                if ( salp->next == NULL ) break;
                else { 
                      salp = salp->next;
                      dsp = (CompSegPtr) salp->segs;
                      dspstart = dsp->starts;
                      dsplens = dsp->lens;
                }
         }
  }
  if ( !seen ) {
    return NULL;
  }
  salp1 = SeqAlignDupAdd (&salp1, salp, numseg_before, subseglens, TRUE);
        /*****************************************
         *** delete salp until delete_to
         *****************************************/
  seen = FALSE;
  while ( !seen )
  {
         if ( !(seen = locate_in_seqalign (delete_to, dsp->dim, dsp->numseg, &dspstart, &dsplens, &numseg_before, &subseglens, &sumseglens)) ) 
         {
                if ( salp->next == NULL ) break;
                else { 
                      salp = salp->next;
                      dsp = (CompSegPtr) salp->segs;
                      dspstart = dsp->starts;
                      dsplens = dsp->lens;
                }
         }
  }
  if ( !seen ) {
         return NULL;
  }
        /*****************************************
         *** copy salp from delete_to to the end
         *****************************************/
  salp1 = SeqAlignDupAdd (&salp1, salp, numseg_before, subseglens, FALSE);
  if ( ( salp = salp->next) != NULL )
     while ( salp != NULL )
     {
        salp1 = SeqAlignDupAdd (&salp1, salp, 0, 0, 0);
        salp = salp->next;
     }
  return salp1;
}

/*********************************************************
***
***  DenseDiagPtr procedures
***
**********************************************************/

static SeqAlignPtr DenseDiagToDenseSegFunc (SeqAlignPtr salp, Boolean add_ends)
{
  BioseqPtr    bsp;
  SeqAlignPtr  newsalp;
  DenseSegPtr  dsp;
  DenseDiagPtr ddp;
  DenseDiagPtr nextddp;
  ValNodePtr   head;
  ValNodePtr   vnp;
  SeqIdPtr     sip1, sip2;
  Int4Ptr      ddp_starts,
               nextddp_starts,
               dspstart;
  Int4Ptr      dsplen;
  Int4         minstart, laststart;
  Int4         ddp_st1, ddp_st2;
  Int4         nextddp_st1, nextddp_st2;
  Int4         ddlen;
  Int4         interlen1, interlen2;
  Int4         intermax, intermin;
  Int4         j;
  Int4         start1 = 0, start2 = 0,
               stop1 = 0, stop2 = 0;
  Int4         slpstart1 = 0, slpstart2 = 0,
               slpstop1 = 0, slpstop2 = 0;
  Uint1        strand1, strand2;
  Int2         numseg;
  Boolean      delete;

  newsalp = SeqAlignNew ();
  newsalp->type = 3;
  newsalp->segtype = 2;
  newsalp->dim = 2;
  dsp = DenseSegNew ();
  newsalp->segs = (Pointer) dsp;
  dsp->dim = 2;

  ddp =(DenseDiagPtr)salp->segs;
  dsp->ids = SeqIdDupList (ddp->id);

  numseg = 0;
  for (ddp =(DenseDiagPtr)salp->segs; ddp != NULL; ddp = ddp->next) 
     numseg++;
  numseg = numseg *3 -2;
 
  head = NULL; 
  laststart = -10;
  for (j=0; j<numseg; j++) {
     minstart = 99999; 
     nextddp = NULL;
     for (ddp = (DenseDiagPtr)salp->segs; ddp!=NULL; ddp=ddp->next) {
        if (laststart < *(ddp->starts) && minstart > *(ddp->starts)) {
           minstart = *(ddp->starts);
           nextddp = ddp;
        } 
     }
     if (nextddp!=NULL) {
        ValNodeAddPointer(&head, 0, nextddp);
        laststart = *(nextddp->starts);
     }
  }
  if (head==NULL)
     return NULL;
  for (vnp=head; vnp!=NULL; vnp=vnp->next) {
     ddp = (DenseDiagPtr) vnp->data.ptrvalue;
     ddp->next = NULL;
  }
  salp->segs = (Pointer) head->data.ptrvalue;
  head->data.ptrvalue = NULL;
  vnp = head->next;
  nextddp = (DenseDiagPtr) salp->segs;
  for (; vnp!=NULL; vnp=vnp->next) {
     ddp = (DenseDiagPtr) vnp->data.ptrvalue;
     nextddp->next = ddp;
     nextddp = nextddp->next;
     vnp->data.ptrvalue = NULL;
  }
  head = ValNodeFree (head);
  ddp = (DenseDiagPtr) salp->segs;
  ddlen = ddp->len;
  ddp_starts = ddp->starts;
  ddp_st1 = *ddp_starts;
  ddp_starts++;
  ddp_st2 = *ddp_starts;
  ddp_starts++;
  while (ddp!=NULL) {
     delete=FALSE;
     if (ddp->next != NULL)
     {   
        nextddp = ddp->next;
        nextddp_starts = nextddp->starts;
        nextddp_st1 = *nextddp_starts;
        nextddp_starts++;
        nextddp_st2 = *nextddp_starts;
        interlen1 = nextddp_st1 - ddp_st1 - ddlen;
        interlen2 = nextddp_st2 - ddp_st2 - ddlen;
        if (interlen1 < 0 || interlen2 < 0) {
           if (interlen1 < 0 && interlen2 < 0) {
              ddp->next = nextddp->next;
              nextddp->next = NULL;
              DenseDiagFree(nextddp);
              delete=TRUE;
           }
           else if (interlen1 < 0) {
              if (ABS(interlen1) >= ddlen) {
                ddp->next = nextddp->next;
                nextddp->next = NULL;
                DenseDiagFree(nextddp);
                delete=TRUE;
              }  
           }
           else if (interlen2 < 0) {
              if (ABS(interlen2) >= ddlen) {
                ddp->next = nextddp->next;
                nextddp->next = NULL;
                DenseDiagFree(nextddp);
                delete=TRUE;
              }  
           }
        }
     } 
     if (!delete) {
        ddp = ddp->next;
        if (ddp != NULL)
        {
         ddlen = ddp->len;
         ddp_starts = ddp->starts;
         ddp_st1 = *ddp_starts;
         ddp_starts++;
         ddp_st2 = *ddp_starts;
         ddp_starts++;
        }
     }
  }
  dsp->starts = (Int4Ptr) MemNew ((size_t) ((2*numseg + 4) * sizeof (Int4)));
  dsp->lens  = (Int4Ptr) MemNew ((size_t) ((numseg + 2) * sizeof (Int4))); 
  for (j = 0; j < 2*numseg + 4; j++) dsp->starts [j] = -1;
  for (j = 0; j < numseg + 2; j++)   dsp->lens [j] = 0;
  dspstart = dsp->starts;
  dsplen = dsp->lens;

  ddp =(DenseDiagPtr)salp->segs;
  ddlen = ddp->len;
  ddp_starts = ddp->starts;
  ddp_st1 = *ddp_starts;
  ddp_starts++;
  ddp_st2 = *ddp_starts;
  ddp_starts++;
  numseg = 0;
  while (ddp != NULL) 
  {
     numseg++;
     *dspstart = ddp_st1;
     dspstart++;
     *dspstart = ddp_st2;
     dspstart++;
     if (ddp->next != NULL) 
     {
        nextddp = ddp->next;
        nextddp_starts = nextddp->starts;
        nextddp_st1 = *nextddp_starts;
        nextddp_starts++;
        nextddp_st2 = *nextddp_starts;
        interlen1 = nextddp_st1 - ddp_st1 - ddlen;
        interlen2 = nextddp_st2 - ddp_st2 - ddlen;
        if (interlen1 < 0 || interlen2 < 0) {
           if (interlen1 < 0 && interlen2 < 0) {
              return NULL;
           }
           if (interlen1 < 0) {
              if (ABS(interlen1) < ddlen) ddlen += interlen1;
              else {
                 return NULL;
              }
           }
           if (interlen2 < 0) {
              if (ABS(interlen2) < ddlen) ddlen += interlen2;
              else {
                 return NULL;
              }
           }
           interlen1 = nextddp_st1 - ddp_st1 - ddlen;
           interlen2 = nextddp_st2 - ddp_st2 - ddlen;
        }
        *dsplen = ddlen;
        dsplen++;
        if (interlen1 <= 0)
        {
           numseg++;
           *dspstart = -1;
           dspstart++;
           *dspstart = ddp_st2 + ddlen;
           dspstart++;
           *dsplen = interlen2;
           dsplen++;
        }
        else if (interlen2 <= 0)
        {
           numseg++;
           *dspstart = ddp_st1 + ddlen;
           dspstart++;
           *dspstart = -1;
           dspstart++;
           *dsplen = interlen1;
           dsplen++;
        }
        else if (interlen1 == interlen2)
        {
           numseg++;
           *dspstart = ddp_st1 + ddlen;
           dspstart++;
           *dspstart = ddp_st2 + ddlen;
           dspstart++;
           *dsplen = interlen1;
           dsplen++;
        }
        else 
        {
           if (interlen1 > interlen2) {
              intermax = interlen1;
              intermin = interlen2;
           }
           else  {
              intermax = interlen2;
              intermin = interlen1;
           }
           numseg++;
           *dspstart = ddp_st1 + ddlen;
           dspstart++;
           *dspstart = ddp_st2 + ddlen;
           dspstart++;
           *dsplen = intermin;
           dsplen++;
           numseg++;
           if (interlen1 > interlen2) {
              *dspstart = ddp_st1 + ddlen + intermin;
              dspstart++;
              *dspstart = -1; 
           }
           else {
              *dspstart = -1; 
              dspstart++;
              *dspstart = ddp_st2 + ddlen + intermin;
           }
           dspstart++;
           *dsplen = intermax - intermin;
           dsplen++;
        }
     }
     else {
        *dsplen = ddlen;
        dsplen++;
     }
     ddp = ddp->next;
     if (ddp != NULL) 
     {
        ddlen = ddp->len;
        ddp_starts = ddp->starts;
        ddp_st1 = *ddp_starts;
        ddp_starts++;
        ddp_st2 = *ddp_starts;
        ddp_starts++;
     }
  }
  dsp->numseg = numseg;
  if (add_ends && newsalp!=NULL) 
  {
     strand1 = SeqAlignStrand (newsalp, 0);
     strand2 = SeqAlignStrand (newsalp, 1);
     start1 = SeqAlignStart (newsalp, 0);
     start2 = SeqAlignStart (newsalp, 1);
     sip1 = SeqAlignId (newsalp, 0);
     sip2 = SeqAlignId (newsalp, 1);
     slpstart1= 0;
     bsp = BioseqLockById(sip1);
     if (bsp!=NULL) {
        slpstop1 = bsp->length-1;
        BioseqUnlock (bsp);
     }
     else
        slpstop1 = stop1;
     slpstart2 = 0; 
     bsp = BioseqLockById(sip2);
     if (bsp!=NULL) {
        slpstop2 = bsp->length-1;
        BioseqUnlock (bsp);
     }  
     else
        slpstop2 = stop2;
     newsalp = SeqAlignEndExtend (newsalp, slpstart1, slpstart2, -1, -1, start1, start2, -1, -1, strand1, strand2); 
     stop1 = SeqAlignStop (newsalp, 0);
     stop2 = SeqAlignStop (newsalp, 1); 
     newsalp = SeqAlignEndExtend (newsalp, -1, -1, slpstop1, slpstop2, -1, -1, stop1, stop2, strand1, strand2);
  } 
  return newsalp;
}

extern SeqAlignPtr DenseDiagToDenseSeg (SeqAlignPtr salp, Boolean add_ends)
{
  SeqAlignPtr cur = NULL, 
              newsalp = NULL, 
              new1 = NULL, 
              new = NULL;

  for (cur=salp; cur!= NULL; cur = cur->next) {
     newsalp = DenseDiagToDenseSegFunc(cur, add_ends);
     if (newsalp != NULL) {
        if (new1 == NULL) { 
           new1 = newsalp;
           new = new1;
        } else {
           new->next = newsalp; 
           new = new->next;
        }
     }
  }
  return new1;
}

extern SeqAlignPtr  DenseSegToDenseDiag (SeqAlignPtr salp)
{
  SeqAlignPtr     salptmp,
                  salp2 = NULL,
                  salp2tmp;
  DenseSegPtr     dsp;
  DenseDiagPtr    ddp,
                  ddphead = NULL;
  SeqIdPtr        sip;
  Int4Ptr         startp,
                  lenp;
  Int2            j;

  if (salp!=NULL)
  {
     for (salptmp = salp; salptmp!=NULL; salptmp = salptmp->next)
     {
        if (salptmp->segtype == 2)
        {
           ddphead = NULL;
           dsp = (DenseSegPtr) salptmp->segs;
           if (dsp!=NULL)
           {
              sip = dsp->ids;
              lenp = (Int4Ptr) dsp->lens;
              startp = dsp->starts;
              for (j=0; j<dsp->numseg; j++, lenp++)
              {
                 if (*startp > -1 && *(startp+1) > -1)
                 {
                    ddp = DenseDiagCreate (2, sip, startp, *lenp, NULL, NULL);
                    DenseDiagAdd (&ddphead, ddp);
                 }
                 startp += dsp->dim;
              }
           }
           if (ddphead != NULL)
           {
              salp2tmp = SeqAlignNew ();
              if (salp2tmp != NULL) {
                 salp2tmp->type = 3;
                 salp2tmp->segtype = 1;
                 salp2tmp->dim = dsp->dim;
                 salp2tmp->segs = (Pointer) ddphead;
                 salp2 = SeqAlignAdd (&salp2, salp2tmp);
              }
           }
        }
        else if (salptmp->segtype == 1)
        {
           salp2tmp = SeqAlignDup (salptmp);
           salp2 = SeqAlignAdd (&salp2, salp2tmp);
        }
     }
  }
  return salp2;
}

extern DenseDiagPtr DenseDiagCreate (Int4 dim, SeqIdPtr id, Int4Ptr starts, Int4 len, Uint1Ptr strands, ScorePtr scores)

{
  DenseDiagPtr ddp_copy;
  Int4         j;

  ddp_copy = DenseDiagNew();
  ddp_copy->dim = dim;
  if (id != NULL) {
     ddp_copy->id = SeqIdDupList (id);
  }
  ddp_copy->starts = (Int4Ptr)MemNew((size_t)((dim+1)*sizeof(Int4)));
  for (j = 0; j < dim; j++, starts++) {
         ddp_copy->starts [j] = *starts;   
  }
  ddp_copy->len = len;
  if ( strands != NULL ) 
     ddp_copy->strands = strands;
  if ( scores != NULL ) 
     ddp_copy->scores = scores;
  return ddp_copy;
}

extern DenseDiagPtr DenseDiagDup (DenseDiagPtr ddp)

{
  DenseDiagPtr ddp_copy;
  SeqIdPtr     sip;
  ScorePtr     sp;
  Int4         dim;
  Int4         j;

  ddp_copy = DenseDiagNew();
  ddp_copy->dim = dim = ddp->dim;
  ddp_copy->id = NULL;
  for (sip = ddp->id, j = 0; sip != NULL && j < dim; sip = sip->next, j++) {
         AddSeqId (&ddp_copy->id, sip);
  }
  ddp_copy->starts = MemNew((size_t)((dim+1)*sizeof(Int4)));
  for ( j = 0; j < dim; j++) 
         ddp_copy->starts [j] = ddp->starts [j];
  ddp_copy->len = ddp->len;
  ddp_copy->strands = MemNew((size_t)((dim+1)*sizeof(Uint1)));
  for ( j = 0; j < dim; j++) 
         ddp_copy->strands [j] = ddp->strands [j];
  ddp_copy->scores = NULL;
  for (sp = ddp->scores; sp != NULL; sp = sp->next ) {
         ScoreDupAdd ( &ddp_copy->scores, sp);
  }
  ddp_copy->next = NULL;
  return ddp_copy;
}

extern DenseDiagPtr DenseDiagAdd (DenseDiagPtr *ddp_head, DenseDiagPtr ddp)
{
  DenseDiagPtr ddp_tmp;

  if (ddp!=NULL)
  {
     if (*ddp_head != NULL)
     {
        ddp_tmp = *ddp_head;
        while (ddp_tmp->next != NULL) 
           ddp_tmp = ddp_tmp->next; 
        ddp_tmp->next = ddp;
     } 
     else 
        *ddp_head = ddp;
  }
  return *ddp_head;
}

extern DenseDiagPtr DenseDiagInsert (DenseDiagPtr ddp_before, DenseDiagPtr ddp)

{
  DenseDiagPtr ddp_tmp;

  if ( (ddp_tmp = ddp_before->next) == NULL) {
         ddp_before->next = ddp;
         return ddp_before;
  }
  ddp_before->next = ddp;
  ddp->next = ddp_tmp;
  return ddp_before;
}

extern DenseDiagPtr DenseDiagPrecede (DenseDiagPtr ddp_after, DenseDiagPtr *ddp)

{
  DenseDiagPtr ddp_tmp;

  ddp_tmp = *ddp;
  ddp_tmp->next = ddp_after;
  return *ddp;
}

extern DenseDiagPtr DenseDiagSortAdd (DenseDiagPtr *ddp_head, DenseDiagPtr ddp)

{
  DenseDiagPtr ddp_tmp;

  if ( (ddp_tmp = *ddp_head) != NULL ) {
         if ( *(ddp->starts) < *(ddp_tmp->starts) ) 
                *ddp_head = DenseDiagPrecede (ddp_tmp, &ddp);

         else if ( *(ddp->starts) > *(ddp_tmp->starts) && ddp_tmp->next == NULL)
                ddp_tmp->next = ddp;

         else {
                while ( ddp_tmp->next != NULL ) {
                       if ( *(ddp->starts) < *(ddp_tmp->next->starts) ) break;
                       ddp_tmp = ddp_tmp->next; 
                }
                ddp_tmp = DenseDiagInsert (ddp_tmp, ddp);
         }
  } 
  else *ddp_head = ddp;
  return *ddp_head;
}

extern void DenseDiagPrint (ValNodePtr ddp)
{
  ValNodePtr   vnp;
  DenseDiagPtr ddp_cur;
  SeqIdPtr     ddp_sip;
  Int4Ptr      ddp_starts;
  Char    strLog[128];

  for ( vnp = ddp; vnp != NULL; vnp = vnp->next)
  {
        if ( vnp->choice != 0 ) break;
        if ( (ddp_cur = (DenseDiagPtr) vnp->data.ptrvalue) != NULL) 
        {
              ddp_starts = ddp_cur->starts;
              for (ddp_sip = ddp_cur->id; ddp_sip != NULL; ddp_sip = ddp_sip->next, ddp_starts++) 
              {
                    SeqIdWrite (ddp_sip, strLog, PRINTID_FASTA_LONG, 120);
              }
        }
  }
}

/*****************************************************************
***
***    IS_sipindensediag
***
*****************************************************************/
extern Boolean IS_seqidindensediag (SeqIdPtr sip, ValNodePtr ddia_list, SeqAlignPtr salp, Int2 index, Int4 from, Int4 to, DenseDiagPtr *block, Int2 intersalpwidth)
{
  ValNodePtr   vnp;
  DenseDiagPtr ddp_cur;
  DenseDiagPtr ddp_new;
  Int4Ptr      ddp_starts;
  SeqIdPtr     ddp_sip;
  Int4         start;
  Boolean      same = FALSE;
  Boolean      is_in = FALSE;

  *block = NULL;
  for ( vnp = ddia_list; vnp != NULL; vnp = vnp->next)
  {
        if ( vnp->choice != 0 ) break;
        if ( (ddp_cur = (DenseDiagPtr) vnp->data.ptrvalue) != NULL) 
        {
              ddp_starts = ddp_cur->starts;
              for (ddp_sip = ddp_cur->id; ddp_sip != NULL; ddp_sip = ddp_sip->next) 
              {
                    same = SeqIdForSameBioseq (sip, ddp_sip);
                    if ( same ) 
                    {
                          start = SeqCoordToAlignCoord (*(ddp_starts), sip, salp, intersalpwidth, 0);
                          if (start > to || start + ddp_cur->len <= from ) {
                                continue;
                          }
                          else {
                                is_in = TRUE;
                                ddp_new = DenseDiagCreate (1, sip, &start, ddp_cur->len, NULL, NULL);
                                DenseDiagSortAdd (block, ddp_new);
                          }
                    }
                    ddp_starts++;
                    /*  
                    {
                    Char  str [128];
                    Char  str2[128];
                    SeqIdPtr (sip, strLog, PRINTID_FASTA_LONG);
                    SeqIdPtr (ddp_sip, strLog, PRINTID_FASTA_LONG);
                    }
                    */
              }
        }
  }
  return is_in;
}


extern DenseDiagPtr GetDenDiag (SeqAlignPtr salp, Int2 index, Int2 *index_entry)
{
  DenseDiagPtr ddp;
  Int2         j = 0;
  Int2         k = 0;

  while ( salp != NULL ) {
        if ( salp->segtype == 1 ) 
        {
              for (ddp =(DenseDiagPtr)salp->segs; ddp != NULL; ddp = ddp->next) 
              {
                    j++;
                    /*
                    {
                    Char  str [128];
                    Char  str2 [128];
                    SeqIdPtr (ddp->id, strLog, PRINTID_FASTA_LONG);
                    WriteLog ("GetDenDiag %s \n", strLog);
                    SeqIdPtr ((ddp->id->next, strLog, PRINTID_FASTA_LONG);
                    WriteLog ("and GetDenDiag %d %d %s \n", j,index, strLog);
                    }
                    */
                    if ( index == j ) break;
                    if ( ddp->next == NULL ) break;
              }
              if ( index == j ) break;
              if ( salp->next == NULL ) break;
              salp = salp->next;
              k++; 
        }
  }
  if ( index != j ) return NULL;
  *index_entry = k;
  return ddp;
}

extern SeqAlignPtr SeqAlignDiagAdd (SeqAlignPtr headp, Int4 pos, Int4 len)
{
  DenseDiagPtr ddp,
               ddpp;
  SeqAlignPtr  salp;

  ddp = DenseDiagCreate (1, NULL, &pos, len, NULL, NULL);
  if (headp == NULL) {
     salp = SeqAlignNew();
     salp->type = 2;
     salp->dim = 1;
     salp->segtype = 3;
     salp->segs = (Pointer)ddp;
  }
  else {
     salp = headp;
     ddpp = (DenseDiagPtr)salp->segs;
     for (; ddpp->next!=NULL; ddpp=ddpp->next)
        continue;
     ddpp->next = ddp;
  }
  return salp;
}


static SeqAlignPtr AddDenseDiagToSeqAlign (SeqAlignPtr salp, DenseDiagPtr dendia, Int2 numseg)

{
  SeqAlignPtr newsalp;
  DenseSegPtr newdsp;
  DenseSegPtr dsp;
  Int4Ptr     gap;
  Int4Ptr     dspstart;
  Int4Ptr     dsplen;
  Int4Ptr     denstart;
  Int4Ptr     newstart;
  Int4Ptr     newlens;
  Int2        n;
  Int2        gapi, newseg;
  Int4        min, minlens,
              max, maxlens;
  Int4        j;
  Int2        k;

  dsp = (DenseSegPtr) salp->segs;
  n       = dsp->dim;
  dspstart= dsp->starts;
  dsplen  = dsp->lens;
  denstart= dendia->starts;

                /* locate the segment before the dendia */
  if ( numseg > 1) {
    for ( k = 1; k < numseg; k++ ) {
      dspstart += n;
      dsplen++;
    }
  }
                /* array of gap lenghts */

  gap = MemNew ((size_t) ((n+1) * sizeof(Int4)));
  max = 0;
  for ( j = 0; j < n; j++) {
    if  ( denstart [j] >= 0 && dspstart [j] >= 0 )
      if ( denstart [j] - dspstart [j] > max ) 
        max = denstart [j] - dspstart [j];
  }
  maxlens = max;
  min = maxlens;
  for ( j = 0; j < n; j++) {
    if  ( denstart [j] >= 0 && dspstart [j] >= 0 )
      if ( denstart [j] - dspstart [j] < min ) 
        min = denstart [j] - dspstart [j];
  }
  minlens = min;
  gapi = 0;
  gap [gapi] = maxlens;
  for (k = 1; k < n; k++ ) {
    max = minlens;
    for ( j = 0; j < n; j++) {
      if  ( denstart [j] >= 0 && dspstart [j] >= 0 ) {
        if ( denstart [j] - dspstart [j] >= max 
          && denstart [j] - dspstart [j] < gap [gapi]) 
          max = denstart [j] - dspstart [j];
      }
    }
    gapi++;
    gap [gapi] = max;
    if ( max == minlens ) break;
  }
  for (k = 0; k <= gapi; k++ ) {
    gap[k] = maxlens - gap [k];
  }
  newseg = dsp->numseg + 2 * gapi + 2;
  newsalp = SeqAlignNew ();
  newsalp->type = salp->type;
  newsalp->segtype = salp->segtype;
  newsalp->dim = salp->dim;
  newdsp = DenseSegNew ();
  newsalp->segs = (Pointer) newdsp;
  newdsp->dim = n;
  newdsp->numseg = newseg;
  newdsp->starts = (Int4Ptr) MemNew ((size_t) ((n*newseg + 4) * sizeof (Int4)));
  newdsp->lens  = (Int4Ptr) MemNew ((size_t) ((newseg + 2) * sizeof (Int4))); 
  for (j = 0; j < n*newseg + 4; j++) newdsp->starts [j] = -1;
  for (j = 0; j < newseg + 2; j++)   newdsp->lens [j] = 0;

  newstart = newdsp->starts;
  newlens  = newdsp->lens;
  dspstart = dsp->starts;
  dsplen   = dsp->lens;
  for ( k = 1; k < numseg; k++ ) {
      for (j = 0; j < n; j++) newstart [j] = dspstart[j];
      *newlens = *dsplen;
      newstart += n;
      newlens++;
      dspstart += n;
      dsplen++;
  }
                     /* block dsp before dendia */

  for (j = 0; j < n; j++) newstart [j] = dspstart[j];
  *newlens = minlens;
  newstart += n;
  newlens++;
                     /* blocks with gaps before dendia */

  for (k = 0; k < gapi; k++ ) {
    for (j = 0; j < n; j++) {
      if (  ( denstart [j] - dspstart [j] >= maxlens - gap[k] ) 
        || ( denstart [j] == SALSA_ND && dspstart [j] >= 0 ) )
        newstart [j] = dspstart [j] + minlens + gap[k];
    }
    *newlens =  gap[k+1];
    newstart += n;
    newlens++;

  }
                     /* block dendia */

  for (j = 0; j < n; j++) 
    if ( denstart[j] >= 0 )
        newstart [j] = denstart[j];
    else if ( dspstart[j] >= 0 )
        newstart [j] = dspstart [j] + maxlens;
  *newlens = dendia->len;
  dspstart += n;
  dsplen++;
  denstart = dspstart;
  dspstart = newstart;
  newstart+= n;
  newlens++;
  
                /* array of gap lenghts */
  max = 0;
  for ( j = 0; j < n; j++) {
    if  ( denstart [j] >= 0 && dspstart [j] >= 0 )
      if ( denstart [j] - dspstart [j] > max ) 
        max = denstart [j] - dspstart [j];
  }
  maxlens = max;
  min = maxlens;
  for ( j = 0; j < n; j++) {
    if  ( denstart [j] >= 0 && dspstart [j] >= 0 )
      if ( denstart [j] - dspstart [j] < min ) 
        min = denstart [j] - dspstart [j];
  }
  minlens = min;
  gapi = 0;
  gap [gapi] = maxlens;
  for (k = 1; k < n; k++ ) {
    max = minlens;
    for ( j = 0; j < n; j++) {
      if  ( denstart [j] >= 0 && dspstart [j] >= 0 ) {
        if ( denstart [j] - dspstart [j] >= max 
          && denstart [j] - dspstart [j] < gap [gapi]) 
          max = denstart [j] - dspstart [j];
      }
    }
    gapi++;
    gap [gapi] = max;
    if ( max == minlens ) break;
  }
  for (k = 0; k <= gapi; k++ ) {
    gap[k] = maxlens - gap [k];
  }
                     /* blocks after dendia with gaps */

  for (k = 0; k < gapi; k++ ) {
    for (j = 0; j < n; j++) {
      if (  ( denstart [j] - dspstart [j] >= maxlens - gap[k] ) 
        || ( denstart [j] == SALSA_ND && dspstart [j] >= 0 ) )
        newstart [j] = dspstart [j] + minlens + gap[k];
    }
    *newlens =  gap[k+1];
    newstart += n;
    newlens++;

  }
                     /* blocks dsp after dendia */

  for ( k = numseg+1; k <= dsp->numseg; k++ ) {
    for (j = 0; j < n; j++) newstart [j] = denstart[j];
    *newlens = *dsplen;
    if ( k < dsp->numseg ) {
      newstart += n;
      newlens++;
      denstart += n;
      dsplen++;
    }    
  }
  return newsalp;
}

static SeqAlignPtr DiagDel (SeqAlignPtr blocp, DenseDiagPtr delp)
{
  DenseDiagPtr ddp = NULL, ddpp = NULL;
 
  if (blocp!=NULL) {
     ddp = blocp->segs;
     if (ddp == delp) {
        if (ddp->next == NULL) {
           blocp->segs = NULL;
           DenseDiagFree (ddp);
           blocp = SeqAlignFree (blocp);
        } else {
           blocp->segs = ddp->next;
           ddp->next = NULL;
           DenseDiagFree (ddp);
        }
     } 
     else {
        ddpp = ddp;
        ddp = ddp->next;
        for (; ddp!=NULL; ddp=ddp->next) {
           if (ddp == delp) {
              ddpp->next = ddp->next;
              ddp->next = NULL;
              DenseDiagFree (ddp);
           }
           if (ddpp->next != NULL)
              ddpp = ddpp->next;
        }
     }
  }
  return blocp;
}

static Int2 Check_compatibility (DenseSegPtr dsp, DenseDiagPtr dendia, Int2 *numseg_before)
{
  Int4Ptr     dspstart;
  Int4Ptr     dspstart_before;
  Int4Ptr     dsplen;
  Int4Ptr     dsplen_before;
  Int4Ptr     denstart;
  Boolean     curpos, pos1;
  Int2        curdiff, diff1;
  Boolean     compat = TRUE;
  Boolean     first = TRUE;
  Int2        n, numseg;
  Int4        j;

  dspstart_before= dspstart = dsp->starts;
  dsplen_before  = dsplen   = dsp->lens;
  denstart = dendia->starts;
  n = dsp->dim;

  /* *********** check twist */      
  for ( numseg = 0; numseg < dsp->numseg; numseg++ ) {
    j = 0;
    first = compat = TRUE;
    while ( j < n ) {
      if ( denstart[j] >= 0 && dspstart[j] >= 0) {
        curpos = (Boolean) (denstart[j] >= dspstart[j]);
        if ( first ) {
          first = FALSE;
          pos1 = curpos;
        } 
        else {
          compat = ( curpos == pos1 );
          if ( !compat ) break;
        }
      }
      j++;
    }
    if ( !compat ) break;
    if ( curpos ) {
      dspstart_before = dspstart; dsplen_before = dsplen;
    }
    else break;
    dspstart += n;
    dsplen++;
  }
  *numseg_before = numseg;
  dspstart = dspstart_before;
  dsplen = dsplen_before;
  if ( !compat ) return 1;

  /* *********** check prolongation */
  first = TRUE;
  for ( j = 0; j < n; j++ ) {
         if ( denstart[j] >= 0 && dspstart[j] >= 0) {
                curdiff = (Int4)ABS(denstart[j] - dspstart[j]);
                if ( first ) {
                       first = FALSE;
                       diff1 = curdiff;
                } 
                else {
                       compat = (curdiff == diff1);        
                       if ( !compat ) break;
                }
         }
  }
  if ( compat ) return 2;
  if ( !compat ) return 3;

  return 0; 
}

extern SeqAlignPtr DenseDiagAlign (SeqAlignPtr salp, DenseDiagPtr dendia)

{
  DenseDiagPtr ddia_dimn;
  DenseSegPtr  dsp;
  SeqIdPtr     sip, sip2;
  Int2         dim;
  Int2         compatibility_mess;
  Int2         numseg;
  Boolean      same = FALSE;
  Int4Ptr      start;
  Uint1Ptr     strand;
  Int4         j;


  if ( salp->type != 3 ) {
    return salp;
  }
                         /*  validate dendiag */
  dsp = (DenseSegPtr) salp->segs;
  ddia_dimn = DenseDiagNew ();
  dim = salp->dim;
  ddia_dimn->dim = dim;
  ddia_dimn->id = dsp->ids;
  ddia_dimn->len = dendia->len;
  ddia_dimn->starts = MemNew((size_t)((dim+1)*sizeof(Int4)));
  ddia_dimn->strands = MemNew((size_t)((dim+1)*sizeof(Uint1)));
  for ( j = 0; j < dim; j++) {
         ddia_dimn->starts [j] = SALSA_ND;
         ddia_dimn->strands [j] = 0;
  }
  start = dendia->starts;
  strand = dendia->strands;
  for ( sip = dendia->id; sip != NULL; sip = sip->next)
  {
         for(sip2=ddia_dimn->id, j=0; sip2!=NULL && j<dim; sip2=sip2->next, j++) 
         {
                same = SeqIdForSameBioseq (sip, sip2);
                if ( same )
                {
                       ddia_dimn->starts [j] = *start;
                       ddia_dimn->strands[j] = *strand;
                       break;
                 }
         }
         start++;
         strand++;
  }
  compatibility_mess = Check_compatibility (dsp, ddia_dimn, &numseg);
  switch ( compatibility_mess ) {
  
    case 1:
         Message (MSG_OK, "twist"); 
         break;

    case 2:
         Message (MSG_OK, "target included: no necessary change"); 
         break;

    case 3:
         salp = AddDenseDiagToSeqAlign (salp, ddia_dimn, numseg);
         break;

    default:
         break;
  }
  
  return salp;
}

/***********************************************************************
***    
***    DenDiagToSeqLoc
***      read SeqAnnotPtr-densediag
***      n: number of sip
***      return list of ValNodePtr-SeqLocPtr
***
************************************************************************/
extern ValNodePtr DenDiagToSeqLoc (SeqAnnotPtr sap, ValNodePtr adpslp, Int2 blastscore_threshold, Int2 *n)
{
  SeqAlignPtr   salp;
  SeqAlignPtr   pre_salp = NULL;
  DenseDiagPtr  ddp;
  SeqIdPtr      sip;
  Int4Ptr       start; 
  Uint1Ptr      strand;
  ValNodePtr    vnp1 = NULL, vnptmp = NULL;
  SeqLocPtr     sltmp;
  SeqIntPtr     sitmp;
  Boolean       same;
  Int2          nseqloc;
  Int2          j = 0;
  Char    strLog[128];
/*
  ScorePtr      ddscore;
  Char          str[128];
  int           sc;
  SeqAlignPtr   salptmp;
*/

  nseqloc = *n;
  vnp1 = adpslp;
  salp = (SeqAlignPtr) sap->data;
  while ( salp != NULL ) {
             ddp = (DenseDiagPtr) salp->segs;
             /**************************************
             ***  score ? blastscore_threshold
             ***************************************/
/*
             sc = -1;
             while ( ( ddscore = ddp->scores ) != NULL ) 
             {
                   if ( ddscore->id == NULL ) break;
                   StringCpy (str, ddscore->id->str);
                   if ( StringNCmp (str, "Score", 5) == 0) 
                   {
                         sc = (int) ddscore->value.intvalue;
                         break;
                   }
                   ddscore = ddscore->next;
             }
*/
             /**************************************
             ***  no ddscore->id->str  --> read value
             ***************************************/
/*
             if  ( sc == -1 && ddp->scores != NULL) {
                   sc = (int) ddp->scores->value.intvalue;
             }
*/
             /**************************************
             ***  delete dense diag < blastscore_threshold
             ***************************************/
/*
             if ( sc < blastscore_threshold ) {
                   if ( pre_salp != NULL ) pre_salp->next = NULL;
                   while ( salp != NULL ) {
                          salptmp = salp;
                          salp = salp->next;
                          SeqAlignFree (salptmp);
                   }
                   salp = pre_salp;
                   break;
             } 
*/
             /**************************************
             ***  list of ids in densediag
             ***************************************/
             sip   = ddp->id;
             start = ddp->starts;
             strand= ddp->strands;
             for ( j = 0; j < ddp->dim; j++ ) 
             {
                    if ( nseqloc > 0 ) 
                    {
                          for (vnptmp=vnp1; vnptmp!=NULL; vnptmp=vnptmp->next) 
                          {
                              sltmp = (SeqLocPtr) vnptmp->data.ptrvalue;
                              sitmp = (SeqIntPtr) sltmp->data.ptrvalue;
                              same = SeqIdForSameBioseq(sip, SeqLocId(sltmp));
                              if ( same )
                              {
                                 if (*start < SeqLocStart (sltmp))
                                         sitmp->from = *start;
                                 if (*start + ddp->len - 1 > SeqLocStop (sltmp))
                                         sitmp->to = *start + ddp->len - 1;
                                 /* TO DO: comparer strands??? ? */
                                 if ( same ) break;
                              }
                          }
                          if (!same) {
                              sltmp = (SeqLocPtr) ValNodeNew (NULL);
                              sltmp->choice = SEQLOC_INT;
                              sitmp = SeqIntNew ();
                              sitmp->from = *start;
                              sitmp->to = *start + ddp->len - 1;
                              sitmp->strand = *strand;
                              sitmp->id = SeqIdDup (sip);
                              sltmp->data.ptrvalue = (Pointer) sitmp; 
                              ValNodeAddPointer (&vnp1, 0, (Pointer) sltmp); 
                              nseqloc++;

SeqIdWrite (SeqLocId (sltmp), strLog, PRINTID_FASTA_LONG, 120);
Message (MSG_OK, "si-new %d %s %ld %ld\n", nseqloc, strLog, (long) SeqLocStart (sltmp), (long) SeqLocStop (sltmp));
                          } 
                    }
                    else {
                          sltmp = (SeqLocPtr) ValNodeNew (NULL);
                          sltmp->choice = SEQLOC_INT;
                          sitmp = SeqIntNew ();
                          sitmp->from = *start;
                          sitmp->to = *start + ddp->len - 1;
                          sitmp->strand = *strand;
                          sitmp->id = SeqIdDup (sip);
                          sltmp->data.ptrvalue = (Pointer) sitmp; 
                          ValNodeAddPointer (&vnp1, 0, (Pointer) sltmp); 
                          nseqloc++;
                    }
                    sip = sip->next; 
                    start++; 
                    strand++;
             }
             pre_salp = salp;
             if ( salp == NULL ) break;
             salp = salp->next;
  }
  *n = nseqloc;
  return vnp1;
}

/***********************************************************************
***
***********************************************************************/
extern SeqAlignPtr SeqLocToFastaSeqAlign (ValNodePtr vnp)
{
  SeqAlignPtr  salp;
  ValNodePtr   vnptmp;
  DenseSegPtr  dsp;
  Int4Ptr      lengthsort;
  Int4         maxlen,
               len, min, pre_min;
  Int4         j;
  Int2         nseq,
               k, numseg;

  nseq = 0;
  for (vnptmp=vnp; vnptmp!=NULL;vnptmp=vnptmp->next)
     if (vnptmp->data.ptrvalue != NULL) nseq++;
  if (nseq == 0)
     return NULL;
  salp = SeqAlignNew ();
  salp->type = 3;
  salp->segtype = 2;
  salp->dim = nseq;
  dsp = DenseSegNew ();
  salp->segs = (Pointer) dsp;
  dsp->dim = nseq;
  dsp->ids = SeqIdListfromSeqLoc (vnp);

         /****************************
         ** count nb of segments
         ****************************/
  maxlen = MaxLengthSeqLoc (vnp);
  lengthsort = MemNew((size_t) ((nseq+1)*sizeof(Int4)));
  pre_min = 0;
  numseg = 1;
  lengthsort [numseg] = 0;
  for ( j = 0; j < nseq; j++ ) 
  {
         vnptmp = vnp;
         min = maxlen;
         for ( k = 0; k < nseq; k++ ) 
         {
                 len = SeqLocLen ((SeqLocPtr) vnptmp->data.ptrvalue);
                 if ( len < min && len > pre_min ) min = len;
                 vnptmp = vnptmp->next;
                 if ( vnptmp == NULL ) break;
         }
         if ( min > pre_min ) 
         {
                  lengthsort [numseg] = min;
                  pre_min = min;
                  numseg++;
         }
  }
         /****************************
         ** copy starts, lens
         ****************************/
  dsp->starts = (Int4Ptr) MemNew ((size_t) ((nseq*numseg + 4) * sizeof (Int4)));
  dsp->lens   = (Int4Ptr) MemNew ((size_t) ((numseg + 2) * sizeof (Int4))); 
  for (j = 0; j < nseq*numseg + 4; j++) dsp->starts[j] = -1;
  for (j = 0; j < numseg + 2; j++) dsp->lens[j] = 0;
  vnptmp = vnp;
  for ( j = 0; j < nseq; j++ ) 
  {
         dsp->starts[j] = SeqLocStart ((SeqLocPtr) vnptmp->data.ptrvalue);
         vnptmp = vnptmp->next;
         if ( vnptmp == NULL ) break;
  }
  for ( k = 1; k < numseg; k++ ) 
  {
         vnptmp = vnp;
         for ( j = 0; j < nseq; j++ ) 
         {
              if (lengthsort[k] < SeqLocLen((SeqLocPtr) vnptmp->data.ptrvalue)) 
                        dsp->starts [j+k*nseq]= dsp->starts[j]+ lengthsort[k];
              vnptmp = vnptmp->next;
              if ( vnptmp == NULL ) break;
         }
         dsp->lens [k-1] = lengthsort [k] - lengthsort [k-1];
  }
  numseg--;
  dsp->numseg = numseg;
  dsp->strands= (Uint1Ptr) MemNew ((size_t) ((numseg*nseq+4)*sizeof (Uint1)));
  for (j = 0; j < numseg*nseq + 4; j++) 
     dsp->strands[j] = Seq_strand_plus;
  MemFree (lengthsort);
  return salp;
}

/*******************************************************
*** SeqEntryToSeqAlignFunc
***    aligns the bioseqs that are present in a SeqEntry (sep) 
***    returns a SeqAlign
*** SeqEntryToSeqAlign
***    calls SeqEntryToSeqAlignFunc
***    returns a SeqAnnot
******************************************************/

extern SeqAlignPtr SeqEntryToSeqAlignFunc (SeqEntryPtr sep, SeqLocPtr master, Uint1 bsp_mol, Int2 method)
{
  ValNodePtr         vnp = NULL,
                     vnp2 = NULL;
  SeqAlignPtr        salp = NULL;
  Int2               nb;

  vnp = SeqEntryToSeqLoc (sep, &nb, bsp_mol);
  if (vnp != NULL) {
     if (master != NULL) 
     {
        ValNodeAddPointer (&vnp2, 0, master);
        vnp2->next = vnp;
        vnp = vnp2; 
     }
     salp = SeqLocListToSeqAlign (vnp, method, NULL);
  } 
  return salp;
}

extern SeqAnnotPtr SeqEntryToSeqAlign (SeqEntryPtr sep, Uint1 bsp_mol)
{
  SeqAlignPtr salp;

  salp = SeqEntryToSeqAlignFunc (sep, NULL, bsp_mol, PRG_ANYALIGN);
  return SeqAnnotForSeqAlign (salp);
}

static Boolean sap_replace (SeqAnnotPtr sap, SeqAlignPtr salp, Uint1 choice)
{
  if (sap != NULL) {
     for (; sap!= NULL; sap=sap->next) {
        if (sap->type == choice) {
           SeqAlignFree ((SeqAlignPtr)sap->data);
           sap->data = (Pointer)salp;
           return TRUE;
        }
     }   
  }
  return FALSE;
}

extern void ReplaceSeqAlignInSeqEntry (Uint2 entityID, Uint2 itemID, SeqAlignPtr salp)
{
  SeqEntryPtr      sep,
                   sep1 = NULL;
  SeqEntryPtr      sept = NULL;
  BioseqSetPtr     bssp = NULL;
  BioseqPtr        bsp = NULL;
  SeqAnnotPtr      sap = NULL;

  sep = GetBestTopParentForItemID (entityID, itemID, OBJ_BIOSEQ);
  if (sep != NULL) {
     if (IS_Bioseq(sep)) {
        entityID = ObjMgrGetEntityIDForChoice (sep);
        sep1 = GetTopSeqEntryForEntityID (entityID);
        bsp = (BioseqPtr) sep->data.ptrvalue;
     }   
     else if(IS_Bioseq_set(sep)) {
        sep1 = sep;
     }   
     if (sep1 != NULL) {
        bssp = NULL; bsp = NULL;
        if (IS_Bioseq(sep1)) {
           bsp = (BioseqPtr) sep1->data.ptrvalue;
           sap_replace(bsp->annot, salp, 2);
        }
        else if(IS_Bioseq_set(sep1)) {
           bssp = (BioseqSetPtr)sep1->data.ptrvalue;
           while (bssp!=NULL && bssp->_class == 7) {
              sept = bssp->seq_set;
              bssp = NULL; bsp = NULL;
              if (IS_Bioseq(sept))  {
                 bsp = (BioseqPtr) sept->data.ptrvalue;
                 break;
              }
              else if (IS_Bioseq_set(sept))
                 bssp = (BioseqSetPtr) sept->data.ptrvalue;
           }
           if (bssp!=NULL) {
              sap = bssp->annot;
              if((sap==NULL || salp==NULL) && IS_Bioseq(sep)) {
                 bsp = (BioseqPtr) sep->data.ptrvalue;
                 sap_replace(bsp->annot, salp, 2);
              }
              else 
                 sap_replace(sap, salp, 2);
              if (sap==NULL && IS_Bioseq_set(sep)) {
                 bssp = (BioseqSetPtr) sep->data.ptrvalue;
                 for (sept = bssp->seq_set; sept!=NULL; sept=sept->next) {
                    if (IS_Bioseq(sept)) {
                       bsp = (BioseqPtr) sept->data.ptrvalue;
                       sap_replace(bsp->annot, salp, 2);
                    }
                 }  
              }
           }
           else if (bsp!=NULL) {
              sap_replace(bsp->annot, salp, 2);
           }
        }
     }   
  }     
  return;
}

/*********************************************************/
static Pointer sap_empty (SeqAnnotPtr sap, Uint1 type, Pointer PNTR ptr)
{
  SeqAlignPtr      salp = NULL;

  if (sap != NULL) {
     for (; sap!= NULL; sap=sap->next) {
        if (sap->type == type) {
           salp = (SeqAlignPtr) sap->data;
           *ptr = (Pointer) sap;
           break;
        }
     }   
  }
  return salp;
}

typedef struct ccid2 {
  Uint1      choice;
  SeqIdPtr   sip;
  Pointer    sap;
} CcId2, PNTR CcId2Ptr;


static void FindSeqAlignCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqAlignPtr        salp, 
                     salptmp;
  DenseSegPtr        dsp;
  CcId2Ptr           cip;
  Boolean            found;
  Pointer            this_sap;

  cip = (CcId2Ptr)mydata;
  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           salp=sap_empty(bsp->annot, 2, &this_sap);
           if (salp!=NULL) {
              found=FALSE;
              salptmp=salp;
              if (cip->sip!=NULL) {
                 while (!found && salptmp!=NULL) 
                 {
                    dsp = salptmp->segs;
                    found = (Boolean)(position_inIdlist(cip->sip, dsp->ids)>0);
                    salptmp=salptmp->next;
                 }
              }
              if (found || cip->sip==NULL) {
                 if (cip->sap==NULL) {
                    if (cip->choice==OBJ_SEQALIGN)
                       cip->sap = (Pointer)salp;
                    else if (cip->choice==OBJ_SEQANNOT)
                       cip->sap = (Pointer) this_sap;
                 }
              }
           }
        }
     }
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           salp=sap_empty(bssp->annot, 2, &this_sap);
           if (salp!=NULL) {
              found=FALSE;
              salptmp=salp;
              if (cip->sip!=NULL) {
                 while (!found && salptmp!=NULL) 
                 {
                    dsp = salptmp->segs;
                    found = (Boolean)(position_inIdlist(cip->sip, dsp->ids)>0);
                    salptmp=salptmp->next;
                 }
              }
              if (found || cip->sip==NULL) {
                 if (cip->sap==NULL) {
                    if (cip->choice==OBJ_SEQALIGN)
                       cip->sap = (Pointer)salp;
                    else if (cip->choice==OBJ_SEQANNOT)
                       cip->sap = (Pointer) this_sap;
                 }
              }
           }
        }
     }
  }
}

extern Pointer FindSeqAlignInSeqEntry (SeqEntryPtr sep, Uint1 choice)
{
  SeqEntryPtr      sep_head;
  BioseqPtr        bsp;
  Uint2            entityID;
  CcId2            ci;

  if (sep==NULL)
     return NULL;
  ci.sap = NULL;
  ci.sip = NULL;
  if (choice != OBJ_SEQALIGN && choice != OBJ_SEQANNOT)
     return NULL;
  ci.choice = choice;
  if (IS_Bioseq(sep)) {
     bsp = (BioseqPtr) sep->data.ptrvalue;
     if (bsp!=NULL)
        ci.sip = SeqIdDup (bsp->id);
  }
  entityID = ObjMgrGetEntityIDForChoice (sep);
  sep_head = GetTopSeqEntryForEntityID (entityID);
  SeqEntryExplore (sep_head, (Pointer)&ci, FindSeqAlignCallback);
  if (ci.sip != NULL)
     SeqIdFree (ci.sip);
  return ci.sap;
}

/***********************************************
***  ReadBuffer from spp+sap  
***    in : spp, sap, from + to in seq coordinates
***    out: length of buffer + buffer
************************************************/
extern Int4 readbuff_fromseqalign (SeqPortPtr spp, SeqAlignPtr salp,  Int2 index, CharPtr buffer, Int4 from, Int4 to, Int4 offset, Boolean strand)
{
  BioseqPtr   bsp;
  DenseSegPtr dsp;
  SeqIdPtr    sip;
  Int4Ptr     dspstart;
  Int4Ptr     dsplens;
  Int4        sumlens, sumstart;
  Int4        bufflen, buffstart;
  Int4        seglenstobuffer = 0;
  Int4        j;
  Int4        maxlen = 0;
  Int2        numseg;
  Boolean     seen = FALSE;
  Boolean     nogap;
  Boolean     ok = TRUE;

  if (spp == NULL) {
    return 0;
  }
  if (buffer == NULL) {
    return 0;
  }
                   /**********************************
                    ***  locate segment including 'from'
                    ***********************************/
  if ( (dsp = (DenseSegPtr) salp->segs ) == NULL) {
         return 0;
  }  
  if (strand == Seq_strand_minus) {
     sip = dsp->ids;
     for (j = 0; sip!=NULL && j < index; j++) 
        sip=sip->next;
     if (sip!=NULL) {
        bsp = BioseqLockById (sip);
        if (bsp!=NULL) {
           maxlen = bsp->length;
           BioseqUnlock(bsp);
        }
     }     
     if (maxlen==0)
        return 0;
  }
  dsplens = dsp->lens;
  dspstart = dsp->starts + index;
  seen = LocateInSeqAlignDenSeg (from, dsp->dim, dsp->numseg, &dspstart, &dsplens, &numseg, &seglenstobuffer);
  if (!seen) {
    ErrPostEx (SEV_ERROR, 0, 0, "fail in readbuff_fromseqalign_sap [4] %ld %ld %ld %ld %ld %ld", from, dsp->dim, dsp->numseg, *dspstart, *dsplens, seglenstobuffer);
    return 0;
  }
                    /***********************************
                    ***  read segments until 'to'
                    ***********************************/
  bufflen = MIN((Int4)(*dsplens-seglenstobuffer),(Int4)(to - from));
  if (strand == Seq_strand_minus) {
     buffstart = *dspstart + seglenstobuffer;
     if (buffstart>-1)
        buffstart = maxlen - buffstart - bufflen;
  } else {
     buffstart = *dspstart + seglenstobuffer;
  }
  nogap = (*dspstart >= 0);
  if (offset < 0 && !nogap)
     offset = 0;
  sumlens = bufflen;
  sumstart = 0;
/**
        WriteLog ("/ %d %d %d %d %d %d %d\n", (int)seglenstobuffer,
        (int)buffstart, (int)bufflen, (int) *dsplens, (int)*dspstart, to, from);
**/
  while ( ok ) 
  {
    if ( nogap ) 
    {
       if (strand == Seq_strand_minus) {
        offset = ReadBufferFromSep (spp, buffer, (Int4)buffstart,
                                   (Int4)(buffstart+bufflen), offset);
       } else {
        offset = ReadBufferFromSep (spp, buffer, (Int4)buffstart,
                                   (Int4)(buffstart+bufflen), offset);
       }
    } 
    else 
    {
        for (j = 0; j < bufflen; j++, offset++) 
                buffer[offset] = '-';
        buffer[offset] = '\0';
    }
    if ( ! (ok = (numseg < dsp->numseg)) ) 
    {
        if ( salp->next == NULL ) break; 
    }
    else {
        numseg++;
        dspstart += dsp->dim; 
        dsplens++;
    }
    nogap = (*dspstart >= 0);
    bufflen = MIN ((Int4) (*dsplens), (Int4) (to - from));
    if ( nogap ) {
       if (strand == Seq_strand_minus) {
          buffstart = *dspstart + sumstart;
          buffstart = maxlen - buffstart - bufflen; 
       } else {
          buffstart = *dspstart + sumstart;
       }
    }
    sumlens += bufflen;
/**
    if (index != 0) 
        WriteLog ("- %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld \n", 
        (long)offset, (long)seglenstobuffer,
        (long)sumstart, (long)*dspstart, (long)buffstart, 
        (long)*dsplens, (long)from, (long)to, (long)bufflen, (long)sumlens);
**/
  }
  return ((Int4)offset);
}

/*******************************************************
*** DElete dsp->ids
********************************************************/
extern SeqAlignPtr aaSeqAlign_to_dnaSeqAlign (SeqAlignPtr salp, ValNodePtr vnp, ValNodePtr framep)
{
  DenseSegPtr      dsp;
  ValNodePtr       vnptmp,
                   fvnp;
  SeqIdPtr         sip,
                   siptmp;
  Int4Ptr          startp;
  Int4Ptr          lenp;
  Int4             from, sumlens;
  Int2             j, k;
  Uint1            frame;

  if (salp == NULL) 
     return NULL;
  dsp = (DenseSegPtr) salp->segs;
  if (dsp == NULL) 
     return NULL;
  
  j = 0;
  for (vnptmp = vnp; vnptmp != NULL; vnptmp = vnptmp->next)
     if (vnptmp->data.ptrvalue != NULL) j++;
     else break;
  if (j != dsp->dim)
     return NULL;
  sip = dsp->ids;
  while (sip!=NULL) {
     siptmp = sip->next;
     sip = SeqIdFree (sip);
     sip = siptmp;
  }
  sip = SeqIdDup (SeqLocId ((SeqLocPtr) vnp->data.ptrvalue));
  vnptmp = vnp->next;
  dsp->ids = sip;
  for (; vnptmp != NULL; vnptmp = vnptmp->next)
  {
     siptmp = SeqIdDup (SeqLocId ((SeqLocPtr) vnptmp->data.ptrvalue));
     sip->next = siptmp;
     sip = sip->next;
  }
  lenp = dsp->lens;
  for (j = 0; j < dsp->numseg; j++, lenp++) {   
     *lenp = *lenp * 3;
  }
  frame=0;
  if (framep!=NULL) {
     fvnp = framep;
  }
  for (k = 0; k < dsp->dim; k++) 
  {  
     lenp = dsp->lens;
     startp = dsp->starts;
     startp += k;
     while (*startp < 0) startp += dsp->dim;
     from = *startp;
     if (fvnp!=NULL)
        frame=(Uint1) fvnp->data.intvalue;
     if (frame == 2)
        from +=1;
     else if (frame == 3)
        from+=2;
     startp = dsp->starts;
     startp += k;
     sumlens = 0;
     for (j = 0; j < dsp->numseg; j++) {
        if (*startp > -1) {
           *startp = from + sumlens;
           sumlens += *lenp;
        }
        startp += dsp->dim;
        lenp++;
     }
     if (framep!=NULL) 
        fvnp=fvnp->next;
  }
  return salp;
}

extern SeqAnnotPtr aaSeqAnnot_to_dnaSeqAnnot (SeqAnnotPtr sap, ValNodePtr vnp, ValNodePtr framep)
{
  if (sap!=NULL && vnp!=NULL)
     sap->data = (Pointer)aaSeqAlign_to_dnaSeqAlign((SeqAlignPtr) sap->data, vnp, framep);
  return sap;
}


extern SeqAlignPtr SortSeqAlign (SeqAlignPtr PNTR salp)
{
  SeqAlignPtr      salp1, salptmp,
                   newsalp = NULL, 
                   presalp = NULL, 
                   minsalp = NULL, 
                   preminsalp = NULL;
  DenseSegPtr      dsp,
                   dsptmp;
  SeqIdPtr         sip, siptmp;
  Int4Ptr          start;
  Int4             minstart;
  BioseqPtr        bsp;

  salp1=*salp;
  while (salp1!=NULL)
  {
     minstart = 999999;
     presalp=NULL;
     minsalp = preminsalp = NULL;
     salptmp=salp1;
     while (salptmp!=NULL) {
        dsp=(DenseSegPtr) salptmp->segs;
        if (dsp!=NULL) {
           start = dsp->starts;
           if (*start < minstart) {
              minstart = *start;
              minsalp = salptmp;
              preminsalp = presalp;
           }
        }
        presalp=salptmp;
        salptmp=salptmp->next;
     } 
     if (minsalp != NULL) {
        dsp=(DenseSegPtr) minsalp->segs;
        sip=dsp->ids->next;
        bsp=BioseqLockById(sip);
        if (bsp==NULL) {
           minsalp=NULL;
        }
        else BioseqUnlock (bsp); 
     }
     if (minsalp != NULL && newsalp != NULL) {
        dsp=(DenseSegPtr) minsalp->segs;
        sip=dsp->ids->next;
        salptmp=newsalp;
        while (salptmp!=NULL) {
           dsptmp=(DenseSegPtr) salptmp->segs;
           siptmp=dsptmp->ids->next;
           if (SeqIdForSameBioseq(sip, siptmp)) {
              break;
           }
           salptmp=salptmp->next;
        }
        if (salptmp!=NULL)
           minsalp=NULL;
     }
     if (minsalp != NULL) {
        if (preminsalp==NULL) {
           salp1 = salp1->next;
           minsalp->next = NULL;
        }
        else {
           preminsalp->next = minsalp->next;
           minsalp->next = NULL;
        }
        SeqAlignAdd (&newsalp, minsalp);
     }
     else break;
  }
  *salp = newsalp;
  return newsalp;
}

extern SeqAlignPtr SortSeqAlignFromList (SeqAlignPtr salp, Int2Ptr sortlst)
{
  SeqAlignPtr newsalp;
  CompSegPtr  csp, 
              newcsp;
  SeqIdPtr    sip;
  Int4Ptr     lenp;
  Uint1Ptr    stp, st2p;
  Int2        j, k, m;

  csp = (CompSegPtr) salp->segs;
  newsalp = SeqAlignNew ();
  if ( newsalp == NULL ) return NULL;
  newsalp->type = salp->type;
  newsalp->segtype = salp->segtype;
  newsalp->dim = salp->dim;
  newcsp = (CompSegPtr) MemNew ( sizeof (CompSeg) );
  newsalp->segs = newcsp;
  newcsp->dim = csp->dim;
  newcsp->numseg = csp->numseg;
  newcsp->ids = NULL;
  for (j=0; j < csp->dim; j++)
  {
     for (k=0, sip = csp->ids; k < csp->dim && sip!=NULL; k++) {
        if (sortlst[k] == (j+1)) {
           newcsp->ids = AddSeqId (&(newcsp->ids), sip);
           break;
        }
        sip = sip->next;
     }
  }
  newcsp->from = (Int4Ptr) MemNew ((size_t) ((csp->dim + 2) * sizeof (Int4)));
  for (j=0; j < csp->dim; j++)
  {
     for (k=0, lenp = csp->from; k < csp->dim; k++, lenp++)
        if (sortlst[k] == (j+1)) {
           newcsp->from [j] = *lenp;
           break;
        }
  }
  newcsp->starts =(BoolPtr)MemNew((size_t)((csp->dim*csp->numseg+ 4) * sizeof (Boolean)));
  for (j=0; j < csp->dim; j++)
  {
     for (k=0, stp = csp->starts; k < csp->dim; k++, stp ++)
        if (sortlst[k] == (j+1)) {
           for (m=0, st2p = stp; m < csp->numseg; m++, st2p+=csp->dim) {
              newcsp->starts [m*csp->dim + j] = *st2p;
           }
           break;
        }
  }
  newcsp->lens=(Int4Ptr) MemNew((size_t)((csp->numseg + 2) * sizeof (Int4)));
  for (j = 0; j < csp->numseg + 2; j++) 
     newcsp->lens[j] = csp->lens[j];
  
  return newsalp;
}
