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
* File Name:  actutils.c
*
* Author:  Sarah Wheelan
*
* Version Creation Date:   2/00
*
* $Revision: 6.22 $
*
* File Description: utility functions for alignments
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: actutils.c,v $
* Revision 6.22  2001/05/22 12:01:07  wheelan
* changes to avoid overflow on alpha platform
*
* Revision 6.21  2001/03/26 16:45:57  wheelan
* fixed uninitialized variables
*
* Revision 6.20  2001/01/09 23:18:55  lewisg
* fix memory leaks
*
* Revision 6.19  2000/10/23 18:43:30  wheelan
* minor bug fix
*
* Revision 6.18  2000/09/07 04:53:40  sicotte
* fix alignment calls, bad matrix calls, and misc alignments problems for sequence update
*
* Revision 6.16  2000/08/28 16:19:11  sicotte
* moved AlnMgrSeqAlignMergeTwoPairwiseEx AlnMgrSeqAlignMergeTwoPairwise AlnMgrSeqAlignMergePairwiseSet to actutils.c from alignmgr.c
*
* Revision 6.15  2000/07/21 21:41:04  sicotte
* fix bug in AlnMgrForcePairwiseContinuousEx, that was inserting two copies of all the seqaligns when trying to realign both ends
*
* Revision 6.14  2000/05/05 11:53:10  wheelan
* bug fixes in ACT_MakeProfileFromSA
*
* Revision 6.13  2000/05/04 16:45:19  wheelan
* changes to profile builder to accomodate IBMed BLAST results
*
* Revision 6.12  2000/05/03 19:59:42  wheelan
* added fix for NULL alignments
*
* Revision 6.11  2000/05/02 12:00:21  wheelan
* fixed memory leaks
*
* Revision 6.10  2000/05/01 19:54:26  wheelan
* fixed memory leak
*
* Revision 6.9  2000/04/22 15:54:57  wheelan
* bug fixes in profile maker
*
* Revision 6.8  2000/04/18 13:57:14  wheelan
* added AlnMgrForcePairwiseContinuousEx
*
* Revision 6.7  2000/04/11 14:50:28  wheelan
* bug fix in AlnMgrForceContinuous
*
* Revision 6.6  2000/03/16 12:53:45  wheelan
* bug fix in AlnMgrForceContinuous
*
* Revision 6.5  2000/03/14 11:25:47  wheelan
* added ACT_ProfileFree functions
*
* Revision 6.4  2000/03/10 16:57:37  wheelan
* fixed AlnMgrForceContinuous
*
* Revision 6.3  2000/03/07 17:55:02  wheelan
* bug fixes in AlnMgrForcePairwiseContinuous
*
* Revision 6.2  2000/03/02 21:11:06  lewisg
* use bandalign for import sequence, make standalone ddv use viewmgr, make dialogs modal, send color update
*
* Revision 6.1  2000/02/11 17:31:44  kans
* initial checkin of functions depending upon blast/bandalign (SW)
*
* ==========================================================================
*/

#include <actutils.h>
#include <viewmgr.h>

static void StateTableSearch (TextFsaPtr tbl, CharPtr txt, Int2Ptr state, Int4 pos, ACT_sitelistPtr PNTR asp_prev, ACT_sitelistPtr PNTR asp_head);
static Boolean am_isa_gap(Int4 start, Int4 prevstop, Uint1 strand);
static void am_fix_strand(SeqAlignPtr sap, Uint1 strand1, Uint1 strand2);


NLM_EXTERN ACT_CGInfoPtr ACT_FindCpG(BioseqPtr bsp)
{
   ACT_CGInfoPtr    acg;
   ACT_CGInfoPtr    acg_head;
   ACT_CGInfoPtr    acg_prev;
   ACT_sitelistPtr  asp;
   ACT_sitelistPtr  asp_head;
   ACT_sitelistPtr  asp_prev;
   ACT_sitelistPtr  asp_tmp;
   Uint1Ptr         buf;
   CharPtr          c;
   Int2             j;
   Int4             max_len;
   Int4             offset;
   Int4             pos;
   Uint1            prev;
   Uint1            r;
   Uint1            res1;
   Uint1            res2;
   Uint1            residue;
   Int4             start;
   Int2             state;
   Int4             state_r;
   Int2             state_test;
   TextFsaPtr       tbl;
   Int4             x;
   FloatHi          y;

   if (bsp == NULL)
      return NULL;
   if (bsp->mol == Seq_mol_aa)
   {
      Message(SEV_WARNING, "Must use nucleotide sequence\n");
      return NULL;
   }
   if (bsp->length < MAX_LEN)
      max_len = bsp->length+1;
   else
      max_len = MAX_LEN;
   buf = (Uint1Ptr)MemNew((max_len)*sizeof(Uint1));
   state = 0;
   prev = 0;
   pos = 0;
   state_r = 0;
   asp_prev = asp_head = NULL;
   tbl = TextFsaNew();
   if (tbl == NULL)
      return NULL;
   TextFsaAdd(tbl, "CCCGGG");
   TextFsaAdd(tbl, "CCGCGG");
   state_test = 0;
   acg = (ACT_CGInfoPtr)MemNew(sizeof(ACT_CGInfo));
   acg_head = NULL;
   offset = 0;
   while ((residue = ACT_GetResidue(pos, buf, &offset, bsp)) != 0)
   {
      if (residue == 65)
         c = "A";
      else if (residue == 67)
         c = "C";
      else if (residue == 71)
         c = "G";
      else if (residue == 84)
         c = "T";
      else
         c = "N";
      StateTableSearch(tbl, c, &state_test, pos, &asp_prev, &asp_head);
      acg->length++;
      acg->to = pos;
      pos++;
      if (residue == 65)
         acg->a++;
      else if (residue == 67)
         acg->c++;
      else if (residue == 71)
      {
         if (prev == 67)
            acg->cg++;
         acg->g++;
      } else if (residue == 84)
         acg->t++;
      else
         acg->n++;
      prev = residue;
      if (acg->length >= 200)
      {
         if (state == 0)
         {
            if (100*(acg->cg)*(acg->length) > 60*(acg->c)*(acg->g) && 10*(acg->c+acg->g) > 6*(acg->length))
            {
               state = 1;
            } else
            {
               res1 = ACT_GetResidue(acg->from, buf, &offset, bsp);
               if (res1 == 67)
               {
                  res2 = ACT_GetResidue(acg->from+1, buf, &offset, bsp);
                  if (res2 == 71)
                     acg->cg--;
               }
               if (res1 == 65)
                  acg->a--;
               else if (res1 == 67)
                  acg->c--;
               else if (res1 == 71)
                  acg->g--;
               else if (res1 == 84)
                  acg->t--;
               else
                  acg->n--;
               acg->from++;
               acg->length--;
            }
         } else if (state == 1)
         {
            if (100*(acg->cg)*(acg->length) <= 60*(acg->c)*(acg->g) || 10*(acg->c+acg->g) < 6*(acg->length))
            {
               state = 0;
               if (acg_head)
               {
                  acg_prev->next = acg;
                  acg_prev = acg;
               } else
                  acg_head = acg_prev = acg;
               j=0;
               if (acg->from - 2000 < 0)
                  start = 0;
               else
                  start = acg->from - 2000;
               r = 1;
               acg->sequence = (CharPtr)MemNew(20000*sizeof(Char));
               for (x=start; x<(acg->to+2000) && r > 0; x++, j++)
               {
                  r = ACT_GetResidue(x, buf, &offset, bsp);
                  if (r == 65)
                     acg->sequence[j] = 'A';
                  else if (r == 67)
                     acg->sequence[j] = 'C';
                  else if (r == 71)
                     acg->sequence[j] = 'G';
                  else if (r == 84)
                     acg->sequence[j] = 'T';
                  else
                     acg->sequence[j] = 'N';
               }
               y = (FloatHi)((acg->cg)*(acg->length))/(FloatHi)((acg->c)*(acg->g));
               printf("Coordinates: %d to %d CpG: %d to %d conf: %f\n%s\n", start, x, acg->from, acg->to, y, acg->sequence);
               MemFree(acg->sequence);
               acg = (ACT_CGInfoPtr)MemNew(sizeof(ACT_CGInfo));
               acg->from = pos;
               acg->length = 1;
               if (residue == 65)
                  acg->a++;
               else if (residue == 67)
                  acg->c++;
               else if (residue == 71)
                  acg->g++;
               else if (residue == 84)
                  acg->t++;
               else
                  acg->n++;
            } else if (acg->length > 1000)
            {
               res1 = ACT_GetResidue(acg->from, buf, &offset, bsp);
               if (res1 == 65 || res1 == 84)
               {
                  acg->from++;
                  acg->length--;
                  if (res1 == 65)
                     acg->a--;
                  else if (res1 == 84)
                     acg->t--;
               }
            }
         }
      }
   }
   /* check for restriction sites in potential islands found */
   acg = acg_head;
   return acg_head;
   acg_prev = NULL; /* Statement not reached... */
   while (acg)
   {
      asp_prev = NULL;
      asp = asp_head;
      while (asp)
      {
         if (asp->start >= acg->from && asp->start < acg->to - 9)
         {
            if (asp_prev)
            {
               asp_prev->next = asp->next;
               asp_tmp = asp->next;
            } else
               asp_head = asp_tmp = asp->next;
            asp->next = acg->asp;
            acg->asp = asp;
            asp = asp_tmp;
         } else
         {
            asp_prev = asp;
            asp = asp->next;
         }
      }
      /*if (acg->asp == NULL)
      {
         if (acg_prev != NULL)
         {
            acg_prev->next = acg->next;
            acg->next = NULL;
            MemFree(acg);
            acg = acg_prev->next;
         } else
         {
            acg_head = acg->next;
            acg->next = NULL;
            MemFree(acg);
            acg = acg_head;
         }
      } else
      {
         if (acg->asp->next == NULL)
         {
            if (acg_prev != NULL)
            {
               acg_prev->next = acg->next;
               acg->next = NULL;
               MemFree(acg);
               acg = acg_prev->next;
            } else
            {
               acg_head = acg->next;
               acg->next = NULL;
               MemFree(acg);
               acg = acg_head;
            }
         } else
         {
            acg_prev = acg;
            acg = acg->next;
         }
      }*/
   }
   MemFree(buf);
   return acg_head;
}

NLM_EXTERN Uint1 ACT_GetResidue(Int4 pos, Uint1Ptr buf, Int4Ptr offset, BioseqPtr bsp)
{
   Int4        bufsize;
   SeqPortPtr  spp;

   if (offset == NULL)
      return 0;
   if (pos > bsp->length - 1)
      return 0;
   if (buf[0] == 0 || pos > (*offset + MAX_LEN - 1) || pos < *offset)
   {
      if (pos > *offset + MAX_LEN - 1)
         *offset = *offset + MAX_LEN;
      else if (pos < *offset)
         *offset = pos;
      else
         *offset = 0;
      if (bsp->length < MAX_LEN)
      {
         bufsize = bsp->length;
         spp = SeqPortNew(bsp, 0, -1, 0, Seq_code_iupacna);
      } else if (bsp->length < *offset + MAX_LEN-1)
      {
         bufsize = bsp->length - *offset + 1;
         spp = SeqPortNew(bsp, *offset, bsp->length-1, 0, Seq_code_iupacna);
      } else
      {
         bufsize = MAX_LEN;
         spp = SeqPortNew(bsp, *offset, *offset+MAX_LEN-1, 0, Seq_code_iupacna);
      }
      if (spp == NULL)
      {
         Message(SEV_WARNING, "Couldn't create SeqPort\n");
         return 0;
      }
      SeqPortRead(spp, buf, bufsize);
   }
   return (buf[pos-(*offset)]);
}

static void StateTableSearch (TextFsaPtr tbl, CharPtr txt, Int2Ptr state, Int4 pos, ACT_sitelistPtr PNTR asp_prev, ACT_sitelistPtr PNTR asp_head)
{
   ACT_sitelistPtr  asp;
   Char             ch;
   ValNodePtr       matches;
   CharPtr          ptr;
   ValNodePtr       vnp;

   if (tbl == NULL || txt == NULL) return;

   for (ptr = txt, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
      *state = TextFsaNext (tbl, *state, ch, &matches);
      for (vnp = matches; vnp != NULL; vnp = vnp->next) {
         asp = (ACT_sitelistPtr)MemNew(sizeof(ACT_sitelist));
         if (asp_prev)
         {
            if (*asp_prev)
            {
               (*asp_prev)->next = asp;
               *asp_prev = asp;
            } else
               *asp_prev = *asp_head = asp;
            asp->name = (CharPtr) vnp->data.ptrvalue;
            asp->start = pos;
         }
      }
   }
}

NLM_EXTERN ACTProfilePtr ACT_ProfileNew(Boolean nuc)
{
   ACTProfilePtr  app;
   FloatHiPtr     PNTR freq;

   app = (ACTProfilePtr)MemNew(sizeof(ACTProfile));
   if (nuc)
   {
      freq = (FloatHiPtr PNTR)MemNew(ACT_NUCLEN*sizeof(FloatHiPtr));
      app->freq = freq;
      app->nuc = TRUE;
   } else
   {
      freq = (FloatHiPtr PNTR)MemNew(ACT_PROTLEN*sizeof(FloatHiPtr));
      app->freq = freq;
      app->nuc = FALSE;
   }
   return app;
}

/***************************************************************************
*
*  ACT_ProfileFree frees a single profile; ACT_ProfileSetFree frees an
*  entire linked list of profiles.
*
***************************************************************************/
NLM_EXTERN ACTProfilePtr ACT_ProfileFree(ACTProfilePtr app)
{
   Int4  i;
   Int4  j;

   if (app == NULL)
      return NULL;
   if (app->nuc)
      j = ACT_NUCLEN;
   else
      j = ACT_PROTLEN;
   for (i=0; i<j; i++)
   {
      MemFree(app->freq[i]);
   }
   MemFree(app->freq);
   app->next = NULL;
   MemFree(app);
   return NULL;
}

NLM_EXTERN ACTProfilePtr ACT_ProfileSetFree(ACTProfilePtr app)
{
   ACTProfilePtr  app_next;

   while (app != NULL)
   {
      app_next = app->next;
      app->next = NULL;
      ACT_ProfileFree(app);
      app = app_next;
   }
   return NULL;
}

NLM_EXTERN void ACT_BuildProfile(SeqLocPtr slp, ACTProfilePtr PNTR app, Int4Ptr count, Int4 length)
{
   Int4        i;
   Int4        len;
   Uint1       res;
   SeqPortPtr  spp;

   if (app == NULL)
      return;
   if (slp == NULL)
   {
      *count = *count+length;
      if ((*app)->len <= *count)
      {
         *count = 0;
         *app = (*app)->next;
      }
      return;
   }
   len = SeqLocLen(slp);
   if (len <= 0)
      return;
   if ((*app)->len == 0)
   {
      (*app)->len = len;
      if ((*app)->nuc)
      {
         for (i=0; i<ACT_NUCLEN; i++)
         {
            (*app)->freq[i] = (FloatHiPtr)MemNew((*app)->len*sizeof(FloatHi));
         }
      } else
      {
         for (i=0; i<ACT_PROTLEN; i++)
         {
            (*app)->freq[i] = (FloatHiPtr)MemNew((*app)->len*sizeof(FloatHi));
         }
      }
   } else
   {
      if (len > (*app)->len) /* seqloc is longer than the */
         return;          /* existing profile -- don't add it     */
   }
   if ((*app)->nuc)
      spp = SeqPortNewByLoc(slp, Seq_code_ncbi4na);
   else
      spp = SeqPortNewByLoc(slp, Seq_code_ncbistdaa);
   if (spp == NULL)
      return;
   if (*count == 0)
     (*app)->numseq++;
   i=0;
   if ((*app)->nuc == FALSE)
   {
      while ((res = SeqPortGetResidue(spp)) != SEQPORT_EOF && i+*count<((*app)->len))
      {
         (*app)->freq[res][i+*count]++;
         i++;
      }
   } else
   {
      while ((res = SeqPortGetResidue(spp)) != SEQPORT_EOF && i+*count<((*app)->len))
      {
         if (res == 1)
         {
            (*app)->freq[0][i+*count]++;
         } else if (res == 2)
         {
            (*app)->freq[1][i+*count]++;
         } else if (res == 4)
         {
            (*app)->freq[2][i+*count]++;
         } else if (res == 8)
         {
            (*app)->freq[3][i+*count]++;
         } else
         {
            (*app)->freq[4][i+*count]++;
         }
         i++;
      }
   }
   SeqPortFree(spp);
   if (len+*count == (*app)->len)
   {
      *app = (*app)->next;
      *count = 0;
   } else
      *count = *count + len;
   return;
}

NLM_EXTERN FloatHi ACT_ScoreProfile(BioseqPtr bsp, Int4 pos, Uint1 strand, ACTProfilePtr app)
{
   Int4        i;
   Uint1       res;
   FloatHi     retval;
   SeqPortPtr  spp;

   retval = 0;
   if (bsp == NULL || app == NULL || pos < 0)
      return retval;
   if (pos + app->len-1 >= bsp->length)
      return retval;
   if (ISA_na(bsp->mol))
   {
      spp = SeqPortNew(bsp, pos, pos+app->len-1, strand, Seq_code_ncbi4na);
      if (spp == NULL)
         return retval;
      i = 0;
      while ((res = SeqPortGetResidue(spp)) != SEQPORT_EOF && i<app->len)
      {
         if (res == 1)
         {
            retval += app->freq[0][i];
         } else if (res == 2)
         {
            retval += app->freq[1][i];
         } else if (res == 4)
         {
            retval += app->freq[2][i];
         } else if (res == 8)
         {
            retval += app->freq[3][i];
         } else
         {
            retval += app->freq[4][i];
         }
         i++;
      }
      retval = retval / app->len;
      return retval;
   } else
   {
      spp = SeqPortNew(bsp, pos, pos+app->len-1, strand, Seq_code_ncbistdaa);
      if (spp == NULL)
         return retval;
      i = 0;
      while ((res = SeqPortGetResidue(spp)) != SEQPORT_EOF && i<app->len)
      {
         retval += app->freq[res][i];
         i++;
      }
      retval = retval / app->len;
      return retval;
   }
}

NLM_EXTERN void ACT_EstimateConfidence(ACTProfilePtr app)
{
   FloatHi  conf;
   Int4     i;
   Int4     j;
   Int4     max;
   Int4     numres;

   if (app == NULL)
      return;
   if (app->nuc)
      numres = ACT_NUCLEN;
   else
      numres = ACT_PROTLEN;
   conf = 1;
   while (app)
   {
      for (i=0; i<app->len; i++)
      {
         max = 0;
         for (j=0; j<numres; j++)
         {
            if (app->freq[j][i] > max)
               max = app->freq[j][i];
         }
         if (max > 0)
            conf = conf*max;
         if (conf > INT4_MAX)
            conf = INT4_MAX;
      }
      app->confidence = conf;
      app = app->next;
   }
   return;
}

NLM_EXTERN ACTProfilePtr ACT_SortProfilesByConfidence(ACTProfilePtr app)
{
   ACTProfilePtr  app_head;
   ACTProfilePtr  PNTR array;
   Int4           count;
   Int4           i;

   if (app == NULL)
      return NULL;
   app_head = app;
   count = 0;
   while (app != NULL)
   {
      count++;
      app = app->next;
   }
   array = (ACTProfilePtr PNTR)MemNew(count*sizeof(ACTProfilePtr));
   app = app_head;
   count = 0;
   while (app != NULL)
   {
      array[count] = app;
      count++;
      app = app->next;
   }
   HeapSort((Pointer)array, (size_t)count, sizeof(ACTProfilePtr), ACT_CompareProfileConfidence);
   app_head = app = array[0];
   for (i=1; i<count; i++)
   {
      app->next = array[i];
      app = app->next;
   }
   return app_head;
}

NLM_EXTERN int LIBCALLBACK ACT_CompareProfileConfidence(VoidPtr base, VoidPtr large_son)
{
   ACTProfilePtr  app1;
   ACTProfilePtr  app2;

   app1 = *((ACTProfilePtr PNTR) base);
   app2 = *((ACTProfilePtr PNTR) large_son);
   if (app1 == NULL || app2 == NULL)
      return 0;
   if (app1->confidence > app2->confidence)
      return -1;
   else if (app1->confidence < app2->confidence)
      return 1;
   else
      return 0;
}

NLM_EXTERN ACTProfilePtr ACT_MakeProfileFromSA(SeqAlignPtr sap)
{
   AMAlignIndexPtr  amaip;
   AlnMsgPtr        amp;
   ACTProfilePtr    app;
   ACTProfilePtr    app_head;
   ACTProfilePtr    app_prev;
   BioseqPtr        bsp;
   Int4             count;
   Int4             i;
   Int4             j;
   Boolean          more;
   Boolean          nuc;
   Int4             numrows;
   SeqIdPtr         sip;
   SeqLocPtr        slp;

   if (sap == NULL)
      return NULL;
   if (sap->saip == NULL)
      return NULL;
   if (sap->saip->indextype == INDEX_PARENT)
   {
      amaip = (AMAlignIndexPtr)(sap->saip);
      if (amaip->mstype == AM_NEATINDEX || amaip->mstype == AM_LITE || amaip->mstype == AM_NULL)
         return NULL;
   }
   sip = AlnMgrGetNthSeqIdPtr(sap, 1);
   bsp = BioseqLockById(sip);
   if (bsp == NULL)
      return NULL;
   if (ISA_na(bsp->mol))
      nuc = TRUE;
   else
      nuc = FALSE;
   BioseqUnlockById(sip);
   sip = SeqIdFree(sip);
   amp = AlnMsgNew();
   amp->to_m = -1;
   amp->row_num = 1;
   app_head = NULL;
   if (sap->saip->indextype == INDEX_PARENT)
   {
      for (i=0; i<amaip->numseg; i++)
      {
         app = ACT_ProfileNew(nuc);
         app->len = amaip->lens[i];
         if (nuc)
         {
            for (j=0; j<ACT_NUCLEN; j++)
            {
               app->freq[j] = (FloatHiPtr)MemNew(app->len*sizeof(FloatHi));
            }
         } else
         {
            for (j=0; j<ACT_PROTLEN; j++)
            {
               app->freq[j] = (FloatHiPtr)MemNew(app->len*sizeof(FloatHi));
            }
         }
         if (app_head != NULL)
         {
            app_prev->next = app;
            app_prev = app;
         } else
            app_head = app_prev = app;
      }
   } else
   {
      while ((Boolean) (more = AlnMgrGetNextAlnBit(sap, amp)))
      {
         app = ACT_ProfileNew(nuc);
         app->len = amp->to_b - amp->from_b + 1;;
         if (nuc)
         {
            for (j=0; j<ACT_NUCLEN; j++)
            {
               app->freq[j] = (FloatHiPtr)MemNew(app->len*sizeof(FloatHi));
            }
         } else
         {
            for (j=0; j<ACT_PROTLEN; j++)
            {
               app->freq[j] = (FloatHiPtr)MemNew(app->len*sizeof(FloatHi));
            }
         }
         if (app_head != NULL)
         {
            app_prev->next = app;
            app_prev = app;
         } else
            app_head = app_prev = app;
      }
   }
   numrows = AlnMgrGetNumRows(sap);
   for (i=1; i<=numrows; i++)
   {
      amp = AlnMsgReNew(amp);
      amp->to_m = -1;
      amp->row_num = i;
      app = app_head;

      sip = AlnMgrGetNthSeqIdPtr(sap, i);
      bsp = BioseqLockById(sip);
      count = 0;
      while ((Boolean) (more = AlnMgrGetNextAlnBit(sap, amp)) && app != NULL)
      {
         if (amp->gap == 0 && bsp != NULL)
         {
            slp = SeqLocIntNew(amp->from_b, amp->to_b, amp->strand, sip);
            ACT_BuildProfile(slp, &app, &count, 0);
            SeqLocFree(slp);
         } else if (amp->gap != 0)
            ACT_BuildProfile(NULL, &app, &count, (amp->to_b - amp->from_b + 1));
      }
      BioseqUnlockById(sip);
      sip = SeqIdFree(sip);
   }
   ACT_EstimateConfidence(app_head);
   AlnMsgFree(amp);
   return app_head;
}

NLM_EXTERN Boolean ACT_AddBioseqToSAByProfile(SeqAlignPtr sap, BioseqPtr bsp)
{
   ACT_PlaceBoundsPtr  abp;
   ACTProfilePtr       app;
   ACTProfilePtr       app_head;
   ACT_PositionPtr     apos;
   ACT_PositionPtr     aposminus;
   ACT_TopScorePtr     PNTR ats;
   ACT_TopScorePtr     PNTR atsminus;
   ACT_TopScorePtr     ats_tmp;
   Int4                c;
   Int4                i;
   Int4                j;
   Int4                len;
   Int4                nprof;
   FloatHiPtr          scorearray;
   FloatHiPtr          scorearrayminus;

   if (sap == NULL || bsp == NULL)
      return FALSE;
   app = ACT_MakeProfileFromSA(sap);
   if (app == NULL)
      return FALSE;
   app_head = app;
   nprof = 0;
   len = 0;
   while (app)
   {
      nprof += 1;
      len += app->len;
      app = app->next;
   }
   app = app_head;
   if (len > bsp->length)
   {
      ErrPostEx(SEV_WARNING, 0, 0, "Sequence too short to align");
      return FALSE;
   }
   scorearray = (FloatHiPtr)MemNew((bsp->length)*sizeof(FloatHi));
   scorearrayminus = NULL;
   ats = (ACT_TopScorePtr PNTR)MemNew(nprof*sizeof(ACT_TopScorePtr));
   c = 0;
   while (app)
   {
      for (i=0; i<bsp->length; i++)
      {
         scorearray[i] = ACT_ScoreProfile(bsp, i, Seq_strand_plus, app);
      }
      ats[c] = ACT_FindPeakScores(scorearray, bsp->length);
      c++;
      app = app->next;
   }
   ats = ACT_SortAndTruncate(ats);
   apos = (ACT_PositionPtr)MemNew(sizeof(ACT_Position));
   abp = (ACT_PlaceBoundsPtr)MemNew(sizeof(ACT_PlaceBounds));
   abp->ats = ats;
   abp->app = app_head;
   abp->len = bsp->length;
   abp->apos = apos;
   abp->nprof = nprof;
   abp->currprof = nprof-1;
   abp->currpos = (Int4Ptr)MemNew(nprof*sizeof(Int4));
   abp->boundarray = (Int4Ptr)MemNew(nprof*sizeof(Int4));
   abp->numats = (Int4Ptr)MemNew(nprof*sizeof(Int4));
   abp->currats = (ACT_TopScorePtr PNTR)MemNew(nprof*sizeof(ACT_TopScorePtr));
   for (i=0; i<nprof; i++)
   {
      abp->boundarray[i] = abp->ats[i]->pos;
      ats_tmp = abp->ats[i];
      abp->currats[i] = ats_tmp;
      j=0;
      while (ats_tmp != NULL)
      {
         ats_tmp = ats_tmp->next;
         j++;
      }
      abp->numats[i] = j;
      if (i>0)
      {
         ats_tmp = abp->currats[i];
         while (abp->boundarray[i] <= abp->boundarray[i-1])
         {
            ats_tmp = ats_tmp->next;
            if (ats_tmp == NULL)
               return FALSE;
            abp->boundarray[i] = ats_tmp->pos;
            abp->currats[i] = ats_tmp;
            abp->currpos[i]++;
         }
      }
   }
   apos->posarray = (Int4Ptr)MemNew(nprof*sizeof(Int4));
   apos = ACT_PlaceByScore(abp);
   if (ISA_na(bsp->mol)) /* if it's a nucleotide, try the minus strand too */
   {
   }
   for (i=0; i<nprof; i++)
   {
      MemFree(ats[i]);
   }
   MemFree(ats);
   MemFree(scorearray);
   MemFree(apos->posarray);
   MemFree(apos);
   MemFree(abp->boundarray);
   MemFree(abp->currpos);
   MemFree(abp->numats);
   MemFree(abp);
   return FALSE;
}

NLM_EXTERN ACT_TopScorePtr PNTR ACT_SortAndTruncate(ACT_TopScorePtr PNTR ats)
{
   return NULL;
}

NLM_EXTERN ACT_PositionPtr ACT_PlaceByScore(ACT_PlaceBoundsPtr abp)
{
   Int4             i;
   FloatHi          score;

   score = ACT_CalcScore(abp);
   if (score > abp->apos->score)
   {
      abp->apos->score = score;
      for (i=0; i<abp->nprof; i++)
      {
         abp->apos->posarray[i] = abp->boundarray[i];
      }
   }
   while (abp->currpos[abp->currprof] < abp->numats[abp->currprof] - 1)
   {
      abp->currpos[abp->currprof]++;
      abp->currats[abp->currprof] = abp->currats[abp->currprof]->next;
      abp->boundarray[abp->currprof] = abp->currats[abp->currprof]->pos;
      score = ACT_CalcScore(abp);
      if (score > abp->apos->score)
      {
         abp->apos->score = score;
         for (i=0; i<abp->nprof; i++)
         {
            abp->apos->posarray[i] = abp->boundarray[i];
         }
      }
   }
   while(abp->currpos[abp->currprof] >= abp->numats[abp->currprof]-1 && abp->currprof >= 0)
   {
      abp->currprof--;
   }
   if (abp->currprof < 0)
      return (abp->apos);
   for (i=abp->currprof+1; i<abp->nprof; i++)
   {
      abp->currpos[i] = 0;
      abp->currats[i] = abp->ats[i];
      abp->boundarray[i] = abp->currats[i]->pos;
      while (abp->boundarray[i] <= abp->boundarray[i-1])
      {
         if (abp->currpos[abp->currprof] >= abp->numats[abp->currprof]-1)
            return (abp->apos);
         abp->currpos[i]+=1;
         abp->currats[i] = abp->currats[i]->next;
         abp->boundarray[i] = abp->currats[i]->pos;
      }
   }
   abp->currpos[abp->currprof]++;
   abp->currats[abp->currprof] = abp->currats[abp->currprof]->next;
   abp->boundarray[abp->currprof] = abp->currats[abp->currprof]->pos;
   abp->currprof = abp->nprof-1;
   abp->apos = ACT_PlaceByScore(abp);
   return (abp->apos);
}

NLM_EXTERN FloatHi ACT_CalcScore(ACT_PlaceBoundsPtr abp)
{
   ACTProfilePtr    app;
   ACT_TopScorePtr  ats;
   Int4             i;
   Int4             j;
   FloatHi          score;

   app = abp->app;
   score = 0;
   for (i=1; i<abp->nprof; i++)
   {
      if (app == NULL)
         return 0;
      if (abp->boundarray[i-1] + app->len >= abp->boundarray[i])
         return 0;
      j=0;
      ats = abp->ats[i-1];
      while (j<abp->currpos[i-1])
      {
         if (ats == NULL)
            return 0;
         ats = ats->next;
      }
      score += app->confidence*ats->score;
      app = app->next;
   }
   return score;
}

NLM_EXTERN ACT_TopScorePtr ACT_FindPeakScores(FloatHiPtr scorearray, Int4 len)
{
   ACT_TopScorePtr  ats;
   ACT_TopScorePtr  ats_head;
   ACT_TopScorePtr  ats_new;
   ACT_TopScorePtr  ats_newhead;
   ACT_TopScorePtr  ats_newprev;
   ACT_TopScorePtr  ats_prev;
   FloatHi          diff;
   FloatHi          diff_prev;
   Int4             i;

   if (scorearray == NULL)
      return NULL;
   diff = diff_prev = 0;
   diff_prev = scorearray[1] - scorearray[0];
   ats_head = NULL;
   for (i=1; i<len-1; i++)
   {
      diff = scorearray[i+1]-scorearray[i];
      if (diff < 0 && diff_prev >= 0) /* peak */
      {
         ats = (ACT_TopScorePtr)MemNew(sizeof(ACT_TopScore));
         ats->score = scorearray[i-1];
         ats->pos = i-1;
         if (ats_head != NULL)
         {
            ats_prev->next = ats;
            ats_prev = ats;
         } else
            ats_head = ats_prev = ats;
      }
      diff_prev = diff;
   }
   ats = ats_prev = ats_head;
   ats = ats->next;
   ats_newhead = NULL;
   diff_prev = 0;
   while (ats)
   {
      diff = ats->score - ats_prev->score;
      if (diff < 0 && diff_prev >= 0)
      {
         ats_new = (ACT_TopScorePtr)MemNew(sizeof(ACT_TopScore));
         ats_new->score = ats_prev->score;
         ats_new->pos = ats_prev->pos;
         if (ats_newhead != NULL)
         {
            ats_newprev->next = ats_new;
            ats_newprev = ats_new;
         } else
            ats_newhead = ats_newprev = ats_new;
      }
      diff_prev = diff;
      ats_prev = ats;
      ats = ats->next;
   }
   ats_prev = ats_head;
   while (ats_prev)
   {
      ats = ats_prev->next;
      ats_prev->next = NULL;
      MemFree(ats_prev);
      ats_prev = ats;
   }
   return ats_newhead;
}

static FloatHi act_get_eval(Int4 exp)
{
  FloatHi eval;
  Int4 i;

  eval = 1;
  for (i=1; i<=exp; i++)
  {
     eval = eval/10;
  }
  return eval;
}

NLM_EXTERN SeqAlignPtr ACT_GlobalAlignSimple(BioseqPtr bsp1, BioseqPtr bsp2,
                                             Boolean Default)
{
   Char *program;
   BLAST_OptionsBlkPtr options;
   SeqAlign *sap, *sap_final;

   if (bsp1 == NULL || bsp2 == NULL)
      return NULL;

   if (ISA_aa(bsp1->mol))
   {
      if (ISA_aa(bsp2->mol))
         program = StringSave("blastp");
      else
         return NULL;
   } else if (ISA_na(bsp1->mol))
   {
      if (ISA_na(bsp2->mol))
         program = StringSave("blastn");
      else
         return NULL;
   }
   else return NULL;

   options = BLASTOptionNew(program, TRUE);
   if(!Default) {
       options->gapped_calculation = TRUE;
       options->expect_value = act_get_eval(60);
       options->wordsize = 20;
       options->filter_string = StringSave("m L");
   }
   sap = BlastTwoSequences(bsp1, bsp2, program, options);
   MemFree(program);
   BLASTOptionDelete(options);
   if (sap == NULL)
   {
      return NULL;
   }
   if(!AlnMgrIndexSeqAlign(sap)) goto error;
   if (!AlnMgrMakeMultipleByScore(sap)) goto error;
   AlnMgrDeleteHidden(sap, FALSE);
   sap_final = AlnMgrForcePairwiseContinuous(sap);
error:
   SeqAlignSetFree(sap);
   return sap_final;
}


NLM_EXTERN SeqAlignPtr ACT_GlobalAlignTwoSeq(BioseqPtr bsp1, BioseqPtr bsp2)
{
   SeqAlignPtr          sap_final;

   sap_final = ACT_GlobalAlignSimple(bsp1, bsp2, FALSE);
   if(sap_final == NULL) return NULL;
   sap_final->idx.entityID = ObjMgrRegister(OBJ_SEQALIGN, sap_final);
   AssignIDsInEntity(sap_final->idx.entityID, OBJ_SEQALIGN, (Pointer)sap_final);
   return sap_final;
}

NLM_EXTERN SeqAlignPtr AlnMgrForcePairwiseContinuous(SeqAlignPtr sap)
{
   return (AlnMgrForcePairwiseContinuousEx(sap, 0, -1, 0, -1));
}

NLM_EXTERN SeqAlignPtr AlnMgrForcePairwiseContinuousEx(SeqAlignPtr sap, Int4 start_1, Int4 stop_1, Int4 start_2, Int4 stop_2)
{
   AMAlignIndexPtr      amaip;
   AlnMsgPtr            amp;
   BioseqPtr            bsp1;
   BioseqPtr            bsp2;
   Int4                 currstart2;
   DenseSegPtr          dsp;
   DenseSegPtr          dsp_old;
   Int4                 end1;
   Int4                 end2;
   GlobalBandStructPtr  gbsp;
   Int4                 i;
   Boolean              is_prot;
   Int4                 j;
   Boolean              more;
   Int4                 numseg;
   Int4                 prevstop1;
   Int4                 prevstop2;
   RowSourcePtr         rsp;
   SeqAlignPtr          sap_new;
   SeqAlignPtr          sap_tmp;
   SeqIdPtr             sip1;
   SeqIdPtr             sip2;
   SeqLocPtr            slp1;
   SeqLocPtr            slp2;
   Int4                 start1;
   Int4                 start2;
   Int4                 stop1;
   Int4                 stop2;
   Uint1                strand1;
   Uint1                strand2;

   if (sap == NULL)
      return NULL;
   if (sap->saip == NULL)
   {
      if (!AlnMgrIndexSeqAlign(sap))
         return NULL;
   }
   if (sap->saip->indextype != INDEX_PARENT)
      return NULL;
   i = AlnMgrGetNumRows(sap);
   if (i != 2)
      return NULL;
   sip1 = AlnMgrGetNthSeqIdPtr(sap, 1);
   sip2 = AlnMgrGetNthSeqIdPtr(sap, 2);
   bsp1 = BioseqLockById(sip1);
   if (bsp1 == NULL)
      return NULL;
   bsp2 = BioseqLockById(sip2);
   if (bsp2 == NULL)
   {
      BioseqUnlock(bsp1);
      return NULL;
   }
   if (ISA_na(bsp1->mol))
      is_prot = FALSE;
   else
      is_prot = TRUE;
   amaip = (AMAlignIndexPtr)(sap->saip);
   rsp = amaip->rowsource[0];
   prevstop1 = start_1 - 1;
   prevstop2 = start_2 - 1;
   strand1 = AlnMgrGetNthStrand(sap, 1);
   strand2 = AlnMgrGetNthStrand(sap, 2);
   if (strand2 == Seq_strand_minus)
   {
      if (stop_2 == -1)
         prevstop2 = bsp2->length;
      else
         prevstop2 = stop_2+1;
   }
   for (i=0; i<rsp->numsaps; i++)
   {
      sap_new = NULL;
      AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[i]-1], rsp->num_in_sap[i], &start1, &stop1);
      AlnMgrGetNthSeqRangeInSA(amaip->saps[rsp->which_saps[i]-1], 3-rsp->num_in_sap[i], &start2, &stop2);
      if (strand2 == Seq_strand_minus)
         currstart2 = stop2;
      else
         currstart2 = start2;
      if (am_isa_gap(start1, prevstop1, strand1))
      {
         if (am_isa_gap(currstart2, prevstop2, strand2))
         {
             /* HS: gaps in both, try to realign ends */
            slp1 = SeqLocIntNew(prevstop1+1, start1-1, strand1, sip1);
            if (strand2 != Seq_strand_minus)
               slp2 = SeqLocIntNew(prevstop2+1, start2-1, strand2, sip2);
            else
               slp2 = SeqLocIntNew(currstart2+1, prevstop2-1, strand2, sip2);
            gbsp = CreatBandStruct(slp1, slp2, NULL, is_prot, G_BAND_QUADRATIC);
            sap_new = GlobalBandByLoc(gbsp, slp1, slp2, is_prot, G_BAND_QUADRATIC);
            am_fix_strand(sap_new, strand1, strand2);
            gbsp = GlobalBandStructDelete(gbsp);
            SeqLocFree(slp1);
            SeqLocFree(slp2);
            if(sap_new) {
                sap_tmp = (SeqAlignPtr)(sap->segs);
                while (sap_tmp->next != NULL)
                    {
                        sap_tmp = sap_tmp->next;
                    }
                sap_tmp->next = sap_new;
                amaip->numsaps++;/* HS */
            }
         } else  /* extend first alignment ; 
                   HS: extend length of gap on 2nd sequence to cover first sequence */
         { 
            sap_tmp = amaip->saps[0];
            SAIndexFree(sap_tmp->saip);
            if (sap_tmp->score != NULL)
            {
               ScoreSetFree(sap_tmp->score);
               sap_tmp->score = NULL;
            }
            sap_tmp->saip = NULL;
            dsp_old = (DenseSegPtr)(sap_tmp->segs);
            dsp = DenseSegNew();
            dsp->starts = (Int4Ptr)MemNew((dsp_old->numseg+1)*2*sizeof(Int4));
            dsp->lens = (Int4Ptr)MemNew((dsp_old->numseg+1)*sizeof(Int4));
            dsp->starts[1] = -1;
            dsp->starts[0] = 0;
            dsp->lens[0] = start1 - 1;
            dsp->strands = (Uint1Ptr)MemNew((dsp_old->numseg+1)*2*sizeof(Uint1));
            dsp->dim = 2;
            dsp->numseg = dsp_old->numseg+1;
            dsp->ids = dsp_old->ids;
            dsp_old->ids = NULL;
            dsp->strands[0] = strand1;
            dsp->strands[1] = strand2;
            for (j=0; j<(dsp_old->numseg); j++)
            {
               dsp->starts[2*(j+1)] = dsp_old->starts[2*j];
               dsp->starts[2*(j+1)+1] = dsp_old->starts[2*j+1];
               dsp->lens[j+1] = dsp_old->lens[j];
               dsp->strands[2*(j+1)] = strand1;
               dsp->strands[2*j] = strand2;
            }
            sap_tmp->segs = (Pointer)(dsp);
            DenseSegFree(dsp_old);
         }
      } else if (am_isa_gap(currstart2, prevstop2, strand2)) /*extend first alignment*/
          { /* HS extend gap on 1st sequence to cover second sequence */
         sap_tmp = amaip->saps[0];
         SAIndexFree(sap_tmp->saip);
         if (sap_tmp->score != NULL)
         {
            ScoreSetFree(sap_tmp->score);
            sap_tmp->score = NULL;
         }
         sap_tmp->saip = NULL;
         dsp_old = (DenseSegPtr)(sap_tmp->segs);
         dsp = DenseSegNew();
         dsp->starts = (Int4Ptr)MemNew((dsp_old->numseg+1)*2*sizeof(Int4));
         dsp->lens = (Int4Ptr)MemNew((dsp_old->numseg+1)*sizeof(Int4));
         dsp->starts[0] = -1;
         if (strand2 != Seq_strand_minus)
         {
            dsp->lens[0] = start2 - prevstop2 - 1;
            dsp->starts[1] = prevstop2+1;
         } else
         {
            dsp->starts[1] = currstart2+1;
            dsp->lens[0] = prevstop2-currstart2 - 1;
         }
         dsp->strands = (Uint1Ptr)MemNew((dsp_old->numseg+1)*2*sizeof(Uint1));
         dsp->dim = 2;
         dsp->numseg = dsp_old->numseg+1;
         dsp->ids = dsp_old->ids;
         dsp_old->ids = NULL;
         dsp->strands[0] = strand1;
         dsp->strands[1] = strand2;
         for (j=0; j<(dsp_old->numseg); j++)
         {
            dsp->starts[2*(j+1)] = dsp_old->starts[2*j];
            dsp->starts[2*(j+1)+1] = dsp_old->starts[2*j+1];
            dsp->lens[j+1] = dsp_old->lens[j];
            dsp->strands[2*(j+1)] = strand1;
            dsp->strands[2*j] = strand2;
         }
         sap_tmp->segs = (Pointer)(dsp);
         DenseSegFree(dsp_old);
      }
      prevstop1 = stop1;
      if (strand2 != Seq_strand_minus)
         prevstop2 = stop2;
      else
         prevstop2 = start2;
   }
   if (strand2 == Seq_strand_minus)
   {
      if (start_1 == 0)
         end2 = -1;
      else
         end2 = start_1+1;
   } else
   {
      if (stop_2 == -1)
         end2 = bsp2->length;
      else
         end2 = stop_2+1;
   }
   if (stop_1 == -1)
      end1 = bsp1->length;
   else
      end1 = stop_1+1;
   if (am_isa_gap(end1, prevstop1, strand1))
   {
      if (am_isa_gap(end2, prevstop2, strand2))
      {
          /* HS gaps in both, realign ends */
         gbsp = GlobalBandStructCreate(G_BAND_QUADRATIC);
         slp1 = SeqLocIntNew(prevstop1+1, bsp1->length-1, strand1, sip1);
         if (strand2 != Seq_strand_minus)
            slp2 = SeqLocIntNew(prevstop2+1, end2-1, strand2, sip2);
         else
            slp2 = SeqLocIntNew(end2+1, prevstop2-1, strand2, sip2);
         gbsp->seqloc1 = slp1;
         gbsp->seqloc2 = slp2;
         sap_new = GlobalBandByLoc(gbsp, slp1, slp2, is_prot, G_BAND_QUADRATIC);
         am_fix_strand(sap_new, strand1, strand2);
         gbsp = GlobalBandStructDelete(gbsp);
         SeqLocFree(slp1);
         SeqLocFree(slp2);
         if(sap_new) {
             sap_tmp = (SeqAlignPtr)(sap->segs);
             while (sap_tmp->next != NULL)
                 {
                     sap_tmp = sap_tmp->next;
                 }
             sap_tmp->next = sap_new;
             amaip->numsaps++; /* HS */
         }
      } else  /* extend last alignment */
      {
         sap_tmp = amaip->saps[amaip->alnsaps-1];
         SAIndexFree(sap_tmp->saip);
         if (sap_tmp->score != NULL)
         {
            ScoreSetFree(sap_tmp->score);
            sap_tmp->score = NULL;
         }
         sap_tmp->saip = NULL;
         dsp_old = (DenseSegPtr)(sap_tmp->segs);
         dsp = DenseSegNew();
         dsp->numseg = dsp_old->numseg+1;
         dsp->starts = (Int4Ptr)MemNew((dsp->numseg+1)*2*sizeof(Int4));
         dsp->lens = (Int4Ptr)MemNew((dsp->numseg+1)*sizeof(Int4));
         dsp->strands = (Uint1Ptr)MemNew((dsp->numseg+1)*2*sizeof(Uint1));
         for (j=0; j<dsp_old->numseg; j++)
         {
            dsp->starts[2*j] = dsp_old->starts[2*j];
            dsp->starts[2*j+1] = dsp_old->starts[2*j+1];
            dsp->lens[j] = dsp_old->lens[j];
            dsp->strands[2*j] = strand1;
            dsp->strands[2*j+1] = strand2;
         }
         dsp->starts[2*j+1] = -1;
         dsp->starts[2*j] = prevstop1+1;
         dsp->lens[j] = bsp1->length - prevstop1 - 1;
         dsp->strands[2*j] = strand1;
         dsp->strands[2*j+1] = strand2;
         dsp->dim = 2;
         dsp->ids = dsp_old->ids;
         dsp_old->ids = NULL;
         sap_tmp->segs = (Pointer)(dsp);
         DenseSegFree(dsp_old);
      }
   } else if (am_isa_gap(end2, prevstop2, strand2)) /* extend last alignment */
   {
      sap_tmp = amaip->saps[amaip->alnsaps-1];
      SAIndexFree(sap_tmp->saip);
      if (sap_tmp->score != NULL)
      {
         ScoreSetFree(sap_tmp->score);
         sap_tmp->score = NULL;
      }
      sap_tmp->saip = NULL;
      dsp_old = (DenseSegPtr)(sap_tmp->segs);
      dsp = DenseSegNew();
      dsp->numseg = dsp_old->numseg+1;
      dsp->starts = (Int4Ptr)MemNew((dsp->numseg+1)*2*sizeof(Int4));
      dsp->lens = (Int4Ptr)MemNew((dsp->numseg+1)*sizeof(Int4));
      dsp->strands = (Uint1Ptr)MemNew((dsp->numseg+1)*2*sizeof(Uint1));
      for (j=0; j<dsp_old->numseg; j++)
      {
         dsp->starts[2*j] = dsp_old->starts[2*j];
         dsp->starts[2*j+1] = dsp_old->starts[2*j+1];
         dsp->lens[j] = dsp_old->lens[j];
         dsp->strands[2*j] = strand1;
         dsp->strands[2*j+1] = strand2;
      }
      dsp->starts[2*j] = -1;
      if (strand2 != Seq_strand_minus)
      {
         dsp->starts[2*j+1] = prevstop2 + 1;
         dsp->lens[j] = bsp2->length - prevstop2 - 1;
      } else
      {
         dsp->starts[2*j+1] = end2+1;
         dsp->lens[j] = prevstop2 - end2 - 1;
      }
      dsp->strands[2*j] = strand1;
      dsp->strands[2*j+1] = strand2;
      dsp->dim = 2;
      dsp->ids = dsp_old->ids;
      dsp_old->ids = NULL;
      sap_tmp->segs = (Pointer)(dsp);
      DenseSegFree(dsp_old);
   }
   BioseqUnlock(bsp1);
   BioseqUnlock(bsp2);
   AlnMgrReIndexSeqAlign(sap);
   amaip = (AMAlignIndexPtr)(sap->saip);
   if (amaip == NULL)
      return NULL;
   amp = AlnMsgNew();
   amp->row_num = 1;
   amp->from_m = 0;
   amp->to_m = -1;
   numseg = 0;
   while ((Boolean) (more = AlnMgrGetNextAlnBit(sap, amp)))
   {
      numseg++;
   }
   sap_new = SeqAlignNew();
   dsp = DenseSegNew();
   dsp->dim = 2;
   dsp->numseg = numseg;
   dsp->starts = (Int4Ptr)MemNew((dsp->numseg*2)*sizeof(Int4));
   dsp->lens = (Int4Ptr)MemNew((dsp->numseg)*sizeof(Int4));
   dsp->strands = (Uint1Ptr)MemNew((dsp->numseg*2)*sizeof(Int4));
   for (i=0; i<2; i++)
   {
      amp = AlnMsgReNew(amp);
      amp->row_num = i+1;
      amp->from_m = 0;
      amp->to_m = -1;
      numseg = 0;
      while ((Boolean) (more = AlnMgrGetNextAlnBit(sap, amp)))
      {
         if (amp->gap == 0)
         {
            dsp->starts[numseg*2+i] = amp->from_b;
            dsp->lens[numseg] = amp->to_b - amp->from_b + 1;
            if (i == 0)
            {
               dsp->strands[numseg*2] = strand1;
            } else
               dsp->strands[numseg*2+i] = strand2;
         } else
         {
            dsp->starts[numseg*2+i] = -1;
            dsp->lens[numseg] = amp->to_b - amp->from_b + 1;
            if (i==0)
               dsp->strands[numseg*2] = strand1;
            else
               dsp->strands[numseg*2+i] = strand2;
         }
         numseg++;
         if (numseg > dsp->numseg)
            return NULL;
      }
   }
   dsp->ids = SeqIdDup(sip1);
   dsp->ids->next = SeqIdDup(sip2);
   sap_new->type = SAT_GLOBAL;
   sap_new->segtype = SAS_DENSEG;
   sap_new->dim = 2;
   sap_new->segs = (Pointer)dsp;
   amp = AlnMsgFree(amp);
   return sap_new;
}

static Boolean am_isa_gap(Int4 start, Int4 prevstop, Uint1 strand)
{
   if (strand != Seq_strand_minus)
   {
      if (start > prevstop+1)
         return TRUE;
      else
         return FALSE;
   } else
   {
      if (prevstop > start+1)
         return TRUE;
      else
         return FALSE;
   }
}

static void am_fix_strand(SeqAlignPtr sap, Uint1 strand1, Uint1 strand2)
{
   DenseSegPtr  dsp;
   Int4         i;

   if (sap == NULL || strand1 == 0 || strand2 == 0)
      return;
   if (sap->segtype != SAS_DENSEG)
      return;
   dsp = (DenseSegPtr)(sap->segs);
   if (dsp->dim != 2)
      return;
   for (i=0; i<dsp->numseg; i++)
   {
      dsp->strands[i*2] = strand1;
      dsp->strands[(i*2) + 1] = strand2;
   }
   return;
}

/*
  Algorithm to merge SeqAligns resulting from two blast hits.
  Untested results If try to merge very non-diagonal SeqAligns.
  
  
 */

NLM_EXTERN SeqAlignPtr AlnMgrSeqAlignMergeTwoPairwiseEx(SeqAlignPtr sap_global,SeqAlignPtr salp1,SeqAlignPtr salp2,Int4 which_master, Uint1 strand_master, Uint1 strand_subject, Int4 startm1,Int4 stopm1, Int4 start1, Int4 stop1,Int4 startm2,Int4 stopm2, Int4 start2, Int4 stop2) {
    /* Overlaps/slight-gaps between two segments can be divided in ten categories
       (five per strand) depending on the relative pentants of the two
       end-points of the alignment.
       .. for ++ alignments (The diagonal is the first Hit,
                   the "end" of the other diagonal falls into one of the
                   5 regions (pentant).
       \    |
        \II |
         \  |   III            region I & II yield overlaps in both master&slave
    I     \ |                  regions III & V will require a gap in one sequence
           \|                  region IV will require gaps in both through a 
 -----------------------                     needleman-wunsch alignment.
            | 
      V     |    IV
            |
            |

for +- alignments
            |    /
            | II/ 
    III     |  /
            | /   I
            |/
 ------------------------
            |
    IV      |    V
            |
            |
                       For regions II and V, could use a simpler "sliding" method
                       to try a simple long gap.. to bridge the two diagonals.
                       XXX For now, simply use Needleman-Wunsch.

    Diagrammatic conventions have the query vertical, and the subject horizontal.

    */
    Int4 gapOpen = 6;
    Int4 gapExtend = 1;
    Int4 reward=1,penalty=-3;
    Int4 score;
    Nlm_FloatHi Kappa=0.0,Lambda=0.0;
    FloatHi H;
    Int4 ** matrix;
    BLAST_MatrixPtr blmat;
    SeqLocPtr slp1,slp2;
    BioseqPtr bsp1,bsp2;
    AlnMsgPtr amp1,amp2;
    Int4 start2_in1,PostGap;
    SAIndexPtr saip;
    Int4 top_c,bottom_c,left_c,right_c;
    Int4 start2_proj_m1,start1_proj_m2;
    Int4 stop2_proj_m1,stop1_proj_m2;
    SeqIdPtr sip1,sip2;
    char * matrix_name;
    SeqAlignPtr salp_new,salp_merging;
    AMAlignIndexPtr amaip;
    Int4 numrows;
    Uint1 mol1,mol2;
    
    amaip = (AMAlignIndexPtr)(sap_global->saip);
    numrows = amaip->numrows;
    if(numrows<2 || numrows>3) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: Invalid call with SeqAlign with too %s row (%d instead of 2 or 3 )\n", numrows<2 ? "little" : "many",numrows);
        return NULL;
    }
    /* Find Alignment coordinates of start2 .. and see if in seqalign */

    sip1 = SeqAlignId(salp1,0);
    sip2 = SeqAlignId(salp1,1);
    if(!sip1 || !sip2) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: SeqAlign with NULL SeqIds\n");
        return NULL;
    }
    bsp1 = BioseqLockById(sip1);
    if(!bsp1) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: SeqAlign with NULL SeqIds\n");
        return NULL;
    }
    bsp2 = BioseqLockById(sip2);
    if(!bsp2) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: SeqAlign with NULL SeqIds\n");
        return NULL;
    }
    mol1 = bsp1->mol;
    mol2 = bsp2->mol;
    if((bsp1->mol == Seq_mol_aa && bsp2->mol != Seq_mol_aa) || (bsp1->mol != Seq_mol_aa && bsp2->mol == Seq_mol_aa)) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: Invalid call with protein-dna alignment \n");
        return NULL;
    }
    BioseqUnlockById(sip1);
    BioseqUnlockById(sip2);

    if(strand_subject != Seq_strand_minus && strand_master != Seq_strand_plus) {
        /*
          Top Left Corner 
        */
        top_c = startm2;
        start2_proj_m1 = AlnMgrMapBioseqToBioseq(sap_global,start2,numrows,1,(Boolean) TRUE,&PostGap);
        if(start2_proj_m1!=-2) {
            if(start2_proj_m1==-1) {
                if(PostGap!=-1)
                    start2_proj_m1 = PostGap;
            }
            if(start2_proj_m1!=-1 && start2_proj_m1<top_c)
                top_c = start2_proj_m1;
        }
        if(stopm1<top_c) {
            top_c = stopm1;
        }
        left_c = AlnMgrMapBioseqToBioseq(sap_global,top_c,1,2,(Boolean) TRUE,&PostGap);

        if(left_c!=-2) {
            if(left_c==-1) {
                if(PostGap!=-1)
                    left_c = PostGap;
            }
        }
        if(left_c<0) {
            ErrPostEx(SEV_WARNING,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: Cannot project coordinate back in SeqAlign : bug in input SeqAlign or AlnMgrMapBioseqToBioseq or  AlnMgrMapBioseqToSeqAlign\n");
            return NULL;
        }
        
        /*
          Bottom Right Corner 
        */
        bottom_c = stopm1;
        stop1_proj_m2 = AlnMgrMapBioseqToBioseq(sap_global,stop1,2,1,(Boolean) TRUE,&PostGap);
        if(stop1_proj_m2!=-2) {
            if(stop1_proj_m2==-1) {
            if(PostGap!=-1)
                stop1_proj_m2 = PostGap;
            }
            if(stop1_proj_m2!=-1 && stop1_proj_m2>bottom_c)
                bottom_c = stop1_proj_m2;
        }
        if(startm2>bottom_c)
            bottom_c = startm2;
        
        right_c = AlnMgrMapBioseqToBioseq(sap_global,bottom_c,1,numrows,(Boolean) TRUE,&PostGap);
        if(right_c!=-2) {
            if(right_c==-1) {
                if(PostGap!=-1)
                    right_c = PostGap;
            }
        }
        if(right_c<0) {
            ErrPostEx(SEV_WARNING,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: Cannot project coordinate back in SeqAlign : bug in input SeqAlign or AlnMgrMapBioseqToBioseq or  AlnMgrMapBioseqToSeqAlign\n");
            return NULL;
        }
        slp1 = SeqLocIntNew(top_c,bottom_c,Seq_strand_plus, SeqAlignId(salp1,0));
        slp2 = SeqLocIntNew(left_c,right_c,Seq_strand_plus, SeqAlignId(salp1,1));
    } else if (strand_subject == Seq_strand_minus && strand_master != Seq_strand_minus) {
        /*
          Top Right Corner 
        */
        top_c = startm2;
        stop2_proj_m1 = AlnMgrMapBioseqToBioseq(sap_global,stop2,numrows,1,(Boolean) TRUE,&PostGap);
        if(stop2_proj_m1!=-2) {
            if(stop2_proj_m1==-1) {
                if(PostGap!=-1)
                    stop2_proj_m1 = PostGap;
            }
            if(stop2_proj_m1!=-1 && stop2_proj_m1<top_c)
                top_c = stop2_proj_m1;
        }
        if(stopm1<top_c)
            top_c = stopm1;
        
        right_c = AlnMgrMapBioseqToBioseq(sap_global,top_c,1,2,(Boolean) TRUE,&PostGap);
        if(right_c!=-2) {
            if(right_c==-1) {
                if(PostGap!=-1)
                    right_c = PostGap;
            }
        }
        if(right_c<0) {
            ErrPostEx(SEV_WARNING,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: Cannot project coordinate back in SeqAlign : bug in input SeqAlign or AlnMgrMapBioseqToBioseq or  AlnMgrMapBioseqToSeqAlign\n");
            return NULL;
        }
        
        /*
          Bottom Left Corner 
        */
        bottom_c = stopm1;
        start1_proj_m2 = AlnMgrMapBioseqToBioseq(sap_global,start1,2,1,(Boolean) TRUE,&PostGap);
        if(start1_proj_m2!=-2) {
            if(start1_proj_m2==-1) {
            if(PostGap!=-1)
                start1_proj_m2 = PostGap;
            }
            if(start1_proj_m2!=-1 && start1_proj_m2>bottom_c)
                bottom_c = start1_proj_m2;
        }
        if(startm2>bottom_c)
            bottom_c = startm2;
        
        left_c = AlnMgrMapBioseqToBioseq(sap_global,bottom_c,1,numrows,(Boolean) TRUE,&PostGap);
        if(left_c!=-2) {
            if(left_c==-1) {
                if(PostGap!=-1)
                    left_c = PostGap;
            }
        }
        if(left_c<0) {
            ErrPostEx(SEV_WARNING,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: Cannot project coordinate back in SeqAlign : bug in input SeqAlign or AlnMgrMapBioseqToBioseq or  AlnMgrMapBioseqToSeqAlign\n");
            return NULL;
        }
        slp1 = SeqLocIntNew(top_c,bottom_c,Seq_strand_plus, SeqAlignId(salp1,0));
        slp2 = SeqLocIntNew(left_c,right_c,Seq_strand_minus, SeqAlignId(salp1,1));
        
    } else {
        ErrPostEx(SEV_WARNING,0,0,"%c %c alignment or unsupported\n",strand_master != Seq_strand_minus ? '+' : '-',strand_subject != Seq_strand_minus ? '+' : '-');
        return NULL;
    }
    if(mol1 == Seq_mol_aa) {
        BlastKarlinGetDefaultMatrixValues("BLOSUM62", &gapOpen, &gapExtend, &Lambda, &Kappa, &H);
        blmat = BLAST_MatrixFetch("BLOSUM62");
        matrix = (Int4Ptr PNTR) blmat->matrix;
        blmat->name = MemFree(blmat->name);
        if(blmat->posFreqs != NULL) {
            Int4 index;
            for (index = 0; index < blmat->rows; index++) {
                MemFree(blmat->posFreqs[index]);
            }
            MemFree(blmat->posFreqs);
        }
        MemFree(blmat);
    } else {
        Int4 index;
        BLASTMatrixStructurePtr matrix_struct;
        matrix_struct = (BLASTMatrixStructurePtr) MemNew(sizeof(BLASTMatrixStructure));
        for (index=0; index<BLAST_MATRIX_SIZE-1; index++)
            {
                matrix_struct->matrix[index] = matrix_struct->long_matrix + index*BLAST_MATRIX_SIZE;
            }
        matrix = (Int4Ptr PNTR ) matrix_struct->matrix;
        BlastScoreBlkMatCreateEx((BLAST_ScorePtr PNTR)matrix,penalty,reward);
        
    }

    salp_merging = NeedlemanWunschQuadraticByLoc(slp1,slp2, matrix, gapOpen, gapExtend,&score, Kappa, Lambda);
    MemFree(matrix);
    AlnMgrReIndexSeqAlign(salp_merging);

    /*
      NOW Merge all 3 SeqAligns..
      Since we insisted that we know the coordinates on BOTH seqaligns..
      all we have to do is Make a new SeqAlign.. and
      use GetNextAlnBit 
     */
    salp_new = AlnMgrMerge3OverlappingSeqAligns(salp1,salp_merging,salp2,top_c,bottom_c);
    SeqAlignFree(salp_merging);
    return salp_new;
}





/*
  Take a Single local SeqAlign and align ends of alignment, or add gaps,
  to make it into a global alignment.
 */

NLM_EXTERN SeqAlignPtr AlnMgrSeqAlignLocalToGlobal(SeqAlignPtr sap) {
    Int4 gapOpen = 5;
    Int4 gapExtend = 1;
    Int4 reward=1,penalty=-3;
    Int4 score;
    Nlm_FloatHi Kappa=0.0,Lambda=0.0;
    FloatHi H;
    BLAST_Score ** matrix=NULL;
    BLAST_MatrixPtr blmat;
    SeqLocPtr slps1,slps2,slpe1,slpe2;
    Int4 len1,len2;
    BioseqPtr bsp1,bsp2;
    SAIndexPtr saip;
    SeqIdPtr sip1,sip2;
    char * matrix_name;
    SeqAlignPtr salp_new,salp_merging,salp_begin,salp_end,sap0;
    Int4 start1,start2,stop1,stop2;
    Uint1 mol1,mol2;
    Uint1 qstrand,sstrand;

    /* Find Alignment coordinates of start2 .. and see if in seqalign */
    if(sap->saip)
        sap0 = ((AMAlignIndexPtr)(sap->saip))->saps[0];
    else
        sap0 = sap;

    sip1 = SeqIdDup(SeqAlignId(sap0,0));
    sip2 = SeqIdDup(SeqAlignId(sap0,1));
    if(!sip1 || !sip2) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: SeqAlign with NULL SeqIds\n");
        return NULL;
    }
    bsp1 = BioseqLockById(sip1);
    if(!bsp1) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: SeqAlign with NULL SeqIds\n");
        return NULL;
    }
    bsp2 = BioseqLockById(sip2);
    if(!bsp2) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignMergeTwoPairwiseEx: SeqAlign with NULL SeqIds\n");
        return NULL;
    }
    mol1 = bsp1->mol;
    mol2 = bsp2->mol;
    len1 = bsp1->length;
    len2 = bsp2->length;
    if((bsp1->mol == Seq_mol_aa && bsp2->mol != Seq_mol_aa) || (bsp1->mol != Seq_mol_aa && bsp2->mol == Seq_mol_aa)) {
        ErrPostEx(SEV_ERROR,0,0,"AlnMgrSeqAlignLocalToGlobal: Invalid call with protein-dna alignment \n");
        return NULL;
    }
    BioseqUnlockById(sip1);
    BioseqUnlockById(sip2);
    SeqIdFree(sip1);
    SeqIdFree(sip2);
    qstrand = SeqAlignStrand(sap0,0);
    sstrand = SeqAlignStrand(sap0,1);
    slps1=NULL;
    slps2=NULL;
    slpe1=NULL;
    slpe2=NULL;
    start1 = SeqAlignStart(sap0,0);
    stop1 = SeqAlignStop(sap0,0);
    start2 = SeqAlignStart(sap0,1);
    stop2 = SeqAlignStop(sap0,1);

    if(qstrand!=Seq_strand_minus) {
        if(start1>0) 
            slps1 = SeqLocIntNew(0,start1-1,Seq_strand_plus, SeqAlignId(sap0,0));
        if(stop1<len1-1) 
            slpe1 = SeqLocIntNew(stop1+1,len1-1,Seq_strand_plus, SeqAlignId(sap0,0));
    } else {
        if(stop1<len1-1) 
            slps1 = SeqLocIntNew(stop1+1,len1-1,Seq_strand_minus, SeqAlignId(sap0,0));
        if(start1>0) 
            slpe1 = SeqLocIntNew(0,start1-1,Seq_strand_minus, SeqAlignId(sap0,0));
    }
    if(sstrand!=Seq_strand_minus) {
        if(start2>0) 
            slps2 = SeqLocIntNew(0,start2-1,Seq_strand_plus, SeqAlignId(sap0,1));
        if(stop2<len2-1) 
            slpe2 = SeqLocIntNew(stop2+1,len2-1,Seq_strand_plus, SeqAlignId(sap0,1));
    } else {
        if(stop2<len2-1) 
            slps2 = SeqLocIntNew(stop2+1,len2-1,Seq_strand_minus, SeqAlignId(sap0,1));
        if(start2>0) 
            slpe2 = SeqLocIntNew(0,start2-1,Seq_strand_minus, SeqAlignId(sap0,1));
    }

    /*
       If either end is not an end-gap 
    */
    if(slps1&&slps2 || slpe1&&slpe2) {
        if(mol1 == Seq_mol_aa) {
            BlastKarlinGetDefaultMatrixValues("BLOSUM62", &gapOpen, &gapExtend, &Lambda, &Kappa, &H);
            blmat = BLAST_MatrixFetch("BLOSUM62");
            matrix = blmat->matrix;
            blmat->name = MemFree(blmat->name);
            if(blmat->posFreqs != NULL) {
                Int4 index;
                for (index = 0; index < blmat->rows; index++) {
                    MemFree(blmat->posFreqs[index]);
                }
                MemFree(blmat->posFreqs);
            }
            MemFree(blmat);
        } else {
            Int4 index;
            BLASTMatrixStructurePtr matrix_struct;
            matrix_struct = (BLASTMatrixStructurePtr) MemNew(sizeof(BLASTMatrixStructure));
            for (index=0; index<BLAST_MATRIX_SIZE-1; index++)
                {
                    matrix_struct->matrix[index] = matrix_struct->long_matrix + index*BLAST_MATRIX_SIZE;
                }
            matrix = matrix_struct->matrix;
            BlastScoreBlkMatCreateEx(matrix,penalty,reward);
        }
    }

    /* 
       Align 1st end
    */
    if(slps1 && slps2) {
        salp_begin = NeedlemanWunschQuadraticByLoc(slps1,slps2, matrix, gapOpen, gapExtend,&score, Kappa, Lambda);
    if(salp_begin)
        AlnMgrIndexSeqAlign(salp_begin);

    } else if (slps1) {
        /* align left end of query with gaps 
          Q 123456789QQQQQQQQQQQQ
          S  ------------12345678
         */
        DenseSegPtr dsp;
        salp_begin = SeqAlignNew();
        salp_begin->segtype = SAS_DENSEG;
        salp_begin->type = SAT_PARTIAL;

        dsp = DenseSegNew();
        salp_begin->segs = dsp;
        dsp->dim=2;
        dsp->numseg=1;
        dsp->starts = (Int4Ptr)MemNew(2*sizeof(Int4));
        dsp->lens = (Int4Ptr)MemNew(2*sizeof(Int4));
        dsp->strands = (Uint1Ptr)MemNew(2*sizeof(Uint1));
        dsp->ids = SeqIdDupList(SeqAlignId(sap0,0));
        dsp->starts[1]=-1;
        /* works for +- strand */
        dsp->starts[0]=SeqLocStart(slps1); 
        dsp->lens[0]=SeqLocStop(slps1)+1-dsp->starts[0];

        dsp->strands[0]=qstrand;
        dsp->strands[1]=sstrand;
    } else if (slps2) {
        /*
          align left end of subject with gaps
          Q  -----------123456789
          S  123456789SSSSSSSSSSS

         */
        DenseSegPtr dsp;
        salp_begin = SeqAlignNew();
        salp_begin->segtype = SAS_DENSEG;
        salp_begin->type = SAT_PARTIAL;

        dsp = DenseSegNew();
        salp_begin->segs = dsp;
        dsp->dim=2;
        dsp->numseg=1;
        dsp->starts = (Int4Ptr)MemNew(2*sizeof(Int4));
        dsp->lens = (Int4Ptr)MemNew(2*sizeof(Int4));
        dsp->strands = (Uint1Ptr)MemNew(2*sizeof(Uint1));
        dsp->ids = SeqIdDupList(SeqAlignId(sap0,0));
        dsp->starts[0]=-1;
        dsp->starts[1]=SeqLocStart(slps2);
        dsp->lens[0]=SeqLocStop(slps2)+1-dsp->starts[1];
        dsp->strands[0]=qstrand;
        dsp->strands[1]=sstrand;
    } else {
        /* Left-end completely aligned. */
        salp_begin = NULL;
    }

    if(slpe1 && slpe2) {
        salp_end = NeedlemanWunschQuadraticByLoc(slpe1,slpe2, matrix, gapOpen, gapExtend,&score, Kappa, Lambda);
    if(salp_end)
        AlnMgrIndexSeqAlign(salp_end);
    } else if (slpe1) {
        DenseSegPtr dsp;
        salp_end = SeqAlignNew();
        salp_end->segtype = SAS_DENSEG;
        salp_end->type = SAT_PARTIAL;

        dsp = DenseSegNew();
        salp_end->segs = dsp;
        dsp->dim=2;
        dsp->numseg=1;
        dsp->starts = (Int4Ptr)MemNew(2*sizeof(Int4));
        dsp->lens = (Int4Ptr)MemNew(2*sizeof(Int4));
        dsp->strands = (Uint1Ptr)MemNew(2*sizeof(Uint1));
        dsp->ids = SeqIdDupList(SeqAlignId(sap0,0));
        dsp->starts[1]=-1;
        dsp->starts[0]=SeqLocStart(slpe1);
        dsp->lens[0]=SeqLocStop(slpe1)+1-dsp->starts[0];
        dsp->strands[0]=qstrand;
        dsp->strands[1]=sstrand;
    } else if (slpe2) {
        DenseSegPtr dsp;
        salp_end = SeqAlignNew();
        salp_end->segtype = SAS_DENSEG;
        salp_end->type = SAT_PARTIAL;

        dsp = DenseSegNew();
        salp_end->segs = dsp;
        dsp->dim=2;
        dsp->numseg=1;
        dsp->starts = (Int4Ptr)MemNew(2*sizeof(Int4));
        dsp->lens = (Int4Ptr)MemNew(2*sizeof(Int4));
        dsp->strands = (Uint1Ptr)MemNew(2*sizeof(Uint1));
        dsp->ids = SeqIdDupList(SeqAlignId(sap0,0));
        dsp->starts[1]=SeqLocStart(slpe2);
        dsp->starts[0]=-1;
        dsp->lens[0]=SeqLocStop(slpe2)+1-dsp->starts[1];
        dsp->strands[0]=qstrand;
        dsp->strands[1]=sstrand;
    } else {
        salp_end = NULL;
    }

    if(sap) {
        AlnMgrIndexSeqAlign(sap);
        sap0 = ((AMAlignIndexPtr)(sap->saip))->saps[0];
    }

    if(salp_begin || salp_end) {
        Int4 cut_start,cut_stop;
        cut_start = SeqAlignStart(sap0,0)-1;
        if(cut_start<0 && salp_begin!=NULL) {
            cut_start=0;
        }
        cut_stop = SeqAlignStop((SeqAlignPtr)sap0,0)+1;
        if(cut_stop>len1-1 && salp_end!=NULL)
            cut_stop = len1-1;
        salp_new = AlnMgrMerge3OverlappingSeqAligns(salp_begin,sap,salp_end,cut_start,cut_stop);
        if(salp_begin)
            SeqAlignFree(salp_begin);
        if(salp_end)
            SeqAlignFree(salp_end);
        if(matrix)
            MemFree(matrix);
    } else {
        salp_new =(SeqAlignPtr) AsnIoMemCopy (sap0, (AsnReadFunc) SeqAlignAsnRead,(AsnWriteFunc) SeqAlignAsnWrite);
      /*  AlnMgrIndexSeqAlign(salp_new);*/
    }
    if(slps1)
        SeqLocFree(slps1);
    if(slps2)
        SeqLocFree(slps2);
    if(slpe1)
        SeqLocFree(slpe1);
    if(slpe2)
        SeqLocFree(slpe2);

    if(salp_new)
        salp_new->type = SAT_GLOBAL;
    return salp_new;
}

NLM_EXTERN SeqAlignPtr AlnMgrSeqAlignMergeTwoPairwise(SeqAlignPtr sap_global,SeqAlignPtr salp1,SeqAlignPtr salp2,Int4 which_master) {
    Uint1 strand_subject, strand_master;
    Int4 start1,start2,stop1,stop2,startm1,startm2,stopm1,stopm2;
    AlnMgrGetNthSeqRangeInSA(salp1, 3-which_master, &start1, &stop1);
    AlnMgrGetNthSeqRangeInSA(salp2, 3-which_master, &start2, &stop2);
    AlnMgrGetNthSeqRangeInSA(salp1, which_master, &startm1, &stopm1);
    AlnMgrGetNthSeqRangeInSA(salp2, which_master, &startm2, &stopm2);
    
    strand_subject = AlnMgrGetNthStrand(salp1, 3-which_master);
    strand_master = AlnMgrGetNthStrand(salp1, which_master);
    if (strand_subject != Seq_strand_minus)
        strand_subject = Seq_strand_plus;

    return AlnMgrSeqAlignMergeTwoPairwiseEx(sap_global,salp1,salp2,which_master,strand_master, strand_subject,startm1, stopm1, start1, stop1, startm2, stopm2, start2, stop2);

}

/*
  Take a list of Pairwise SeqAligns, output from AlnMgrRemoveInconsistentFromPairwiseSet
  and Merges them if possible (if they don't have gaps of AM_MPS_MAXGAP(~1000 nt) or more )
  return the number of SeqAligns it merged the results in, or -1 if failed.
*/

#define AM_MPS_MAX_GAP 1000 
NLM_EXTERN Int4 AlnMgrSeqAlignMergePairwiseSet(SeqAlignPtr PNTR sap_ptr ) {
    AMAlignIndexPtr amaip;
    SAIndexPtr saip1;
    SeqAlignPtr sap,salp;
    AMTinyInfoPtr PNTR tiparray;
    Int4 i;
    FloatHi bit_score,evalue,score;
    Int4 numsegs;
    SeqAlignPtr salp_merged, PNTR saparray,salp_head,salp_last;
    Int4 startm1,stopm1,startm2,stopm2;
    Int4 start1,stop1,start2,stop2;
    Boolean DeleteAln,MergeAlns;
    Uint1 strand,strand_master;

   if (sap_ptr == NULL || *sap_ptr==NULL)
       return 0;
   sap= *sap_ptr;
   if((sap->saip != NULL && sap->saip->indextype != INDEX_PARENT))
      return -1;
   if (sap->saip == NULL)
   {
      if (!AlnMgrIndexLite(sap))
         return -1;
   }
   amaip = (AMAlignIndexPtr)(sap->saip);
   if(!amaip)
       return -1;
   if (amaip->numbsqs > 2)
      return -1;
   if(amaip->numsaps>=2) {
       salp = salp_head = salp_last = (SeqAlignPtr)(sap->segs);
       while(salp_last->next!=NULL) {
           salp_last = salp_last->next;
       }
       
       tiparray = (AMTinyInfoPtr PNTR)MemNew((amaip->alnsaps)*sizeof(AMTinyInfoPtr));
       for (i=0; i<amaip->alnsaps; i++)
           {
               saip1 = (SAIndexPtr)amaip->saps[i]->saip;
               AlnMgrGetNthSeqRangeInSA(amaip->saps[i], saip1->master, &start1, &stop1);
               tiparray[i] = (AMTinyInfoPtr)MemNew(sizeof(AMTinyInfo));
               tiparray[i]->start = start1;
               tiparray[i]->stop = stop1;
               tiparray[i]->numsap = i;
               tiparray[i]->which = 0; /* Important */
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
       
       
       if ((saip1 = (SAIndexPtr)amaip->saps[0]->saip)==NULL) 
           return 0;
       strand = AlnMgrGetNthStrand(amaip->saps[tiparray[0]->numsap], 3-saip1->master);
       strand_master = AlnMgrGetNthStrand(amaip->saps[tiparray[0]->numsap], saip1->master);
       if (strand != Seq_strand_minus)
           strand = Seq_strand_plus;
       numsegs=0;
       startm1 = tiparray[0]->start;
       stopm1 = tiparray[0]->stop;
       AlnMgrGetNthSeqRangeInSA(amaip->saps[0], 3-saip1->master, &start1, &stop1);
       i=1;
       while(i<amaip->numsaps) {
           /*
             if two current alignments overlap or touch.. package
             seqaligns in one dsp
           */
           startm2 = tiparray[i]->start;
           stopm2 = tiparray[i]->stop;
           AlnMgrGetNthSeqRangeInSA(amaip->saps[i], 3-saip1->master, &start2, &stop2);
           DeleteAln = FALSE;
           if (startm1 != startm2) {
               if((strand!=Seq_strand_minus && start2<=start1) || (strand == Seq_strand_minus && stop2<=stop1)) {
                   ErrPostEx(SEV_WARNING,0,0,"BUG: Input SeqAlign to AlnMgrSeqAlignMergePairwiseSet violated input conditions\n");
               } else {
                   Int4 vert_gap,horiz_gap,gap_area;
                   if(strand!=Seq_strand_minus) {
                       vert_gap=startm2-stopm1;
                       horiz_gap=start2-stop1;
                   } else {
                       vert_gap=startm2-stopm1;
                       horiz_gap=start1-stop2;
                   }
                   gap_area = vert_gap*horiz_gap;
                   if(gap_area<0) 
                       gap_area = - gap_area;
                   if(gap_area<AM_MPS_MAX_GAP*AM_MPS_MAX_GAP) {
                       MergeAlns = TRUE;
                   } else {
                       if(((vert_gap>=-20 && vert_gap<=20) || (horiz_gap >=-20 && horiz_gap<=20))
                          && gap_area <= 10*AM_MPS_MAX_GAP*AM_MPS_MAX_GAP) {
                           MergeAlns = TRUE;
                       } else
                           MergeAlns = FALSE;
                   }
                   if(MergeAlns) {
                       salp_merged = AlnMgrSeqAlignMergeTwoPairwiseEx(sap,amaip->saps[i-1],amaip->saps[i],saip1->master, strand_master, (Uint1) Seq_strand_plus, startm1, stopm1, start1, stop1,startm2, stopm2, start2,  stop2);
                       amaip->saps[i-1]=salp_merged;
                       if(salp_head)
                           salp_last->next = salp_merged;
                       else
                           salp_head = salp_merged;
                       salp_last = salp_merged;
                       DeleteAln=TRUE;
                   } else { /* Too much or too little Overlap to merge */
                       ErrPostEx(SEV_WARNING,0,0,"AlnMgrSeqAlignMergePairwiseSet: Input alignments overlap too much/little to merge(Indicative of repeats/large deletions):\nMaster Overlap = %d , Subject Overlap = %d\n",vert_gap,horiz_gap);
                       DeleteAln = FALSE;
                   }
               }
           } else {
               DeleteAln = TRUE;
               /* Overlapping starts!
                  Since Guaranteed that the 
                  first one is longer because of the sort Fn. 
                  --> Delete the shorter one (second alignment )
               */
           }
           if(DeleteAln) {
               amaip->numsaps--;
               amaip->saps[i]=NULL;
               if(i<amaip->numsaps) {
                   Int4 j;
                   for(j=i+1;j<amaip->numsaps;j++) {
                       MemMove(tiparray+j-1,tiparray+j,sizeof(tiparray[0]));
                       amaip->saps[j-1]=amaip->saps[j];
                   }
               }
               /* Deletes SeqAligns not in amaip->saps[0..amaip->numsaps-1] array */
               AlnMgrDeleteHidden(sap,FALSE); /* Delete unused SeqAligns and Re-Index */
               amaip = (AMAlignIndexPtr)(sap->saip);
               saip1 = (SAIndexPtr)amaip->saps[0]->saip;
               /* Try again re-merging current (don't increment i) */
           } else {
               i++;
               startm1 = startm2;
               stopm1 = stopm2;
               start1 = start2;
               stop1 = stop2;
           }
       }
   }
   if(amaip->numsaps==1) {
       SeqAlignPtr sap_tmp;
       sap_tmp = SeqAlignDup((SeqAlignPtr)sap->segs);
       sap->segs=NULL;
       SeqAlignFree(sap);
       sap = sap_tmp;
   }
   *sap_ptr = sap;
   AlnMgrReIndexSeqAlign(sap);
   return amaip->numsaps;
}
