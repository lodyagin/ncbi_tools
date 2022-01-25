/*   salfiles.c
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
* File Name:  salfiles.c
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.49 $
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
#include <salfiles.h>
#include <salstruc.h>
#include <salutil.h>
#include <salsap.h>
#include <salpanel.h>
#include <salparam.h>
#include <biosrc.h>
#include <cdrgn.h>
#include <fstyle.h>
#include <subutil.h>
#include <satutil.h>
#include <tofasta.h>
#include <objacces.h>
#include <accutils.h>
#include <seqmgr.h>

#define MAXSTR          512
#define OBJ_VIRT        254


static void FindBioseqCB3 (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  ValNodePtr  vnp;

  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     if (IS_Bioseq(sep)) {
        vnp = (ValNodePtr)mydata;
        if (vnp->data.ptrvalue==NULL)
        {
           vnp->data.ptrvalue=(BioseqPtr) sep->data.ptrvalue;
        }
     }
  }
}
 
static SelEdStructPtr is_sip_inseqinfo (SeqIdPtr sip, SelEdStructPtr seq_info)
{
  SelEdStructPtr tmp;
  SeqLocPtr      slp;

  for (tmp=seq_info; tmp!=NULL; tmp=tmp->next)
  {
     slp=(SeqLocPtr)tmp->region;
     if (SeqIdForSameBioseq (sip, SeqLocId(slp)))
     {
        return tmp;
     }
  }
  return NULL;
}

extern ValNodePtr CCReadAnythingLoop (CharPtr filename, SelEdStructPtr seq_info)
{
  Char         name [PATH_MAX];
  Pointer      dataptr;
  FILE         *fp;
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqLocPtr    slp;
  ValNodePtr   head = NULL,
               vnp,
               slphead = NULL;
  Uint2        datatype; 
  Uint2        entityID;

  ValNode      vn;

SelEdStructPtr tmp;
SeqIdPtr new_sip;

  if (filename == NULL) 
  {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }
     filename = name;
  }
  fp = FileOpen (filename, "r");
  if (fp != NULL) {
    while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE)) != NULL) {
      ValNodeAddPointer (&head, datatype, dataptr);
    }
    FileClose (fp);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      datatype = vnp->choice;
      dataptr = vnp->data.ptrvalue;
      entityID = ObjMgrRegister (datatype, dataptr);
      if (datatype == OBJ_BIOSEQ)
      {
         bsp=(BioseqPtr)vnp->data.ptrvalue;
         slp = SeqLocIntNew (0, bsp->length-1, Seq_strand_plus, SeqIdFindBest (bsp->id, 0));

if (seq_info) {
 if ((tmp=is_sip_inseqinfo(bsp->id, seq_info)) != NULL) {
   sep=SeqEntryNew ();
   sep->choice=1;
   sep->data.ptrvalue = (Pointer)bsp; 
   SeqEntryReplaceSeqID (sep, SeqLocId(slp));
   ValNodeFree(slp);
   slp = SeqLocIntNew(0, bsp->length-1, Seq_strand_plus, bsp->id);
   sep->data.ptrvalue = NULL;
   SeqEntryFree(sep);
 }
}
         ValNodeAddPointer (&slphead, 0, (Pointer) slp);
      }
      else if (datatype == OBJ_SEQENTRY)
      {
         sep=(SeqEntryPtr)vnp->data.ptrvalue;
         vn.data.ptrvalue=NULL;
         SeqEntryExplore (sep, &vn, FindBioseqCB3);
         if (vn.data.ptrvalue!=NULL) {
            bsp=(BioseqPtr)vn.data.ptrvalue;
            slp = SeqLocIntNew (0, bsp->length-1, Seq_strand_plus, SeqIdFindBest (bsp->id, 0));

if (seq_info) {
 if ((tmp=is_sip_inseqinfo(bsp->id, seq_info)) != NULL) {
   SeqEntryReplaceSeqID (sep, SeqLocId(slp));
   ValNodeFree(slp);
   slp = SeqLocIntNew(0, bsp->length-1, Seq_strand_plus, bsp->id);
 }
}
            ValNodeAddPointer (&slphead, 0, (Pointer) slp);
         }
      }
    }
  }
  return slphead;
}

/*************************************************
***   Sequence: 
***       FastaRead
***
*************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR last_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol     /* Returns special symbol if no SeqEntry */
 );


static SeqEntryPtr NewFastaRead (FILE *fp, Boolean is_na, Boolean parseSeqId, Int2 *seqnumber, Int2 *segnumber, SeqIdPtr PNTR siplst, Int4 *lengthmax)
{
  BioseqPtr     bsp;
  SeqEntryPtr   sep = NULL;
  SeqEntryPtr   lastsep = NULL;
  SeqEntryPtr   nextsep = NULL;
  SeqEntryPtr   last = NULL;
  CharPtr       errormsg = NULL;
  ValNodePtr    head = NULL;
  Char          lastchar;
  ObjectIdPtr   oid = NULL;
  SeqIdPtr      sip;
  SeqIdPtr      siphead = NULL, 
                siptmp;
  Char          str [32];
  ValNodePtr    vnp;
  Int4          count;
  Int4          lensmax = 0;
  Int2          nseq = 0;
  Int2          segcount = 0, 
                segtotal = 0;
  Boolean       insegset;
  Boolean       isLocalUnknownID;
  
  SeqEntryPtr sepnuc;
  BioseqPtr   segbsp;
  SeqIdPtr    segsip=NULL, lastsegsip=NULL;

  count = 0;
  last = sep;
  lastsep = NULL;
  insegset = FALSE;
  nextsep = FastaToSeqEntryInternal ((void *)fp, 2, NULL, is_na, &errormsg, parseSeqId, &lastchar);
  while (nextsep != NULL || (lastchar != EOF && lastchar != NULLB && lastchar != 255)) {
          if (nextsep != NULL) {
            count++;
            if (IS_Bioseq (nextsep) && nextsep->data.ptrvalue != NULL) {
              bsp = (BioseqPtr) nextsep->data.ptrvalue;
              if (bsp->length > lensmax)
                 lensmax = bsp->length;
              isLocalUnknownID = FALSE;
              sip = bsp->id;
              if (sip != NULL && sip->choice == SEQID_LOCAL) {
                oid = (ObjectIdPtr) sip->data.ptrvalue;
                if (oid != NULL && oid->str != NULL) {
                  isLocalUnknownID = (Boolean) (StringICmp (oid->str, "unknown") == 0);
                }
              }
              if ((! parseSeqId) || isLocalUnknownID) {
                oid = ObjectIdNew ();
                if (oid != NULL) {
                  if (is_na) {
                    sprintf (str, "nuc %ld", (long) count);
                  } else {
                    sprintf (str, "prot %ld", (long) count);
                  }
                  oid->str = StringSave (str);
                  sip = ValNodeNew (NULL);
                  if (sip != NULL) {
                    sip->choice = SEQID_LOCAL;
                    sip->data.ptrvalue = (Pointer) oid;
                    bsp->id = SeqIdFree (bsp->id);
                    bsp->id = sip;
                    SeqMgrReplaceInBioseqIndex (bsp);
                  } else {
                    ObjectIdFree (oid);
                  }
                }
              }
            }
            SeqEntryPack (nextsep);
            if (sep != NULL) {     
              if (insegset) {
                if (lastsep != NULL) {
                  AddSeqEntryToSeqEntry (lastsep, nextsep, TRUE);
                  segcount ++;
                  if (segcount > segtotal)
                     segtotal = segcount;
                  sepnuc = FindNucSeqEntry (lastsep);
                  if (IS_Bioseq(sepnuc)) {
                     segbsp=(BioseqPtr)sepnuc->data.ptrvalue;
                     segsip=segbsp->id;
                     if (segsip != NULL) {
                       if (lastsegsip==NULL || !SeqIdMatch(segsip,lastsegsip)) 
                       {
                          siptmp = SeqIdDup (segsip);
                          siphead = AddSeqId (&siphead, siptmp);
                          lastsegsip = segsip;
                       }
                     }
                  }
                } 
                else {
                  lastsep = nextsep;
                  last->next = nextsep;
                  last = nextsep;
                  segcount=1;
                  if (segcount > segtotal)
                    segtotal = segcount;
                  nseq++;
                }
              } 
              else {
                last->next = nextsep;
                last = nextsep;
                segcount=1;
                if (segcount > segtotal)
                  segtotal = segcount;
                nseq++;
                if (lastsegsip==NULL || !SeqIdMatch(sip,lastsegsip)) 
                {
                  siptmp = SeqIdDup (sip);
                  siphead = AddSeqId (&siphead, siptmp);
                  lastsegsip = sip;
                }
              }
            } 
            else {
              if (insegset && lastsep == NULL) {
                lastsep = nextsep;
                sep = nextsep;
                last = sep;
                segcount=1;
                if (segcount > segtotal)
                  segtotal = segcount;
                nseq++;
              } 
              else {
                sep = nextsep;
                last = sep;
                segcount=1;
                if (segcount > segtotal)
                  segtotal = segcount;
                nseq++;
                if (lastsegsip==NULL || !SeqIdMatch(sip,lastsegsip)) 
                {
                  siptmp = SeqIdDup (sip);
                  siphead = AddSeqId (&siphead, siptmp);
                  lastsegsip = sip;
                }
              }
            }
            vnp = ValNodeNew (head);
            if (head == NULL) {
              head = vnp;
            }
            if (vnp != NULL) {
              vnp->data.ptrvalue = errormsg;
              errormsg = NULL;
            }
          } else if (lastchar == '[') {
            insegset = TRUE;
            lastsep = NULL;
          } else if (lastchar == ']') {
            insegset = FALSE;
          }
          nextsep = FastaToSeqEntryInternal ((void *)fp, 2, NULL, is_na, &errormsg, parseSeqId, &lastchar);
  }
  if (segnumber !=NULL) 
     *segnumber = segtotal;
  if (sip!=NULL)
     *siplst = siphead;
  else 
     SeqIdFree (siphead);
  if (lengthmax != NULL)
     *lengthmax = lensmax;
  if(seqnumber != NULL)
     *seqnumber = nseq;
  return sep;
}

static ValNodePtr ReadAlignmentToStrings (CharPtr path, Int4 length, Int2 segnumber)
{
  Char         name[PATH_MAX];
  FILE         *fp;
  ValNodePtr   vnpal, tmp, vnp;
  Char         str[255]; 
  CharPtr      strp,
               seqstr;
  Int4         strlens, 
               lmax=0,
               lgseq=0;
  Int2         inseg = 0;
  Boolean       insegb = FALSE;
  Boolean      goOn,
               startp;
  
  
  Int2         j = 0;   
  Int2         nseq = 0;
  
  if (path == NULL) {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }
     path = name;
  }
  vnpal = ValNodeNew (NULL);
  tmp = vnpal;
  for (j=1; j<segnumber; j++) {
     vnp = ValNodeNew (NULL);
     tmp->next = vnp;
     tmp = tmp->next;
  }
  if ( (fp = FileOpen (path, "r")) == NULL)  {
     ValNodeFree (vnpal);
     return NULL;
  }
  vnp=NULL;
  lmax = length + length/2;
  goOn = (Boolean)(fgets(str, sizeof(str), fp)!=NULL);
  if (goOn) {
     strp = str;
     while (*strp == ' ' && *strp!='\0' && *strp!='\n')
        strp++;
     if (*strp!='\0' && *strp!='\n')
        strlens = StringLen (strp);
     else 
        goOn=FALSE;
  }
  while (goOn) 
  {
     if (strlens > 0) {
        if (*strp == '>') {
           if (!insegb) {
              vnp = vnpal;
           }
           else {
              inseg++;
              if (inseg==1)
                 vnp = vnpal;
              else
                 vnp = vnp->next;
           }
           startp = FALSE;
        }
        else if (StringStr(strp, "[")!= NULL) {
           if (vnp!=NULL) {
           }
           insegb = TRUE;
           inseg = 0;
           startp= FALSE;
        } 
        else if (StringStr(strp, "]")!= NULL) {
           insegb = FALSE;
           inseg = 0;
           startp= FALSE;
        } 
        else {
           if (!startp) {
              seqstr=(CharPtr)MemNew((size_t)((lmax + 1) * sizeof(Char)));
              for (strlens=0; strlens<lmax; strlens++) 
                 seqstr[strlens] = ' ';
              if (vnp->data.ptrvalue==NULL) {
                 tmp = NULL;
                 ValNodeAddPointer (&tmp, 0, (Pointer)seqstr);
                 vnp->data.ptrvalue = (Pointer) tmp;
              } else {
                 tmp = (ValNodePtr)vnp->data.ptrvalue;
                 ValNodeAddPointer (&tmp, 0, (Pointer)seqstr);
              }
              lgseq = 0;
              startp = TRUE;
           }              
           for (j=0; j<strlens; j++)
           {
              if (strp[j]=='\n' || strp[j]=='\0' || strp[j]=='\r' )
                 break;
              strp[j] = TO_UPPER(strp[j]);
              if (StringChr("ABCDEFGHIKLMNPQRSTUVWXYZ-*", strp[j]) != NULL) {
                 seqstr [lgseq] = strp[j];
                 lgseq++;
              }
           }
           seqstr [lgseq] = '\0';
        }
     }
     goOn = (Boolean)(fgets(str, sizeof(str), fp)!=NULL);
     if (goOn) {
        strp = str;
        while (*strp == ' ' && *strp!='\0' && *strp!='\n')
           strp++;
        if (*strp!='\0' && *strp!='\n')
           strlens = StringLen (strp);
        else 
           goOn=FALSE;
     }
  }
  fclose (fp);
  return vnpal;
}

static ValNodePtr get_lens_fromseqalign (SeqAlignPtr salp)
{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  ValNodePtr    fromp = NULL;
  Int4Ptr       startp;
  Int4          j,
                val = (Int4)-1;
  Int2          index;
  Uint1         strand;
 
  if (salp == NULL)
     return NULL;
  if (salp->segtype == 1)
  {
     ddp = (DenseDiagPtr) salp->segs;
     if (ddp != NULL) {
      for (index=0; index<ddp->dim; index++) {
        startp = ddp->starts;
        if (index > 0)
           startp += index;
        val = *startp + ddp->len;
        ValNodeAddInt (&fromp, 1, (Int4)(val+1));
      }
     }   
  }  
  else if (salp->segtype == 2)
  {
     dsp = (DenseSegPtr) salp->segs;
     if (dsp!=NULL)
     {   
      for (index=0; index<dsp->dim; index++) 
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
           ValNodeAddInt (&fromp, 1, (Int4)(val+1));
        }
        else 
           ValNodeAddInt (&fromp, 1, (Int4)(-1));
      }
     }
  }     
  return fromp;
}

static SeqAnnotPtr LocalAlign1ToSeqAnnotDimn (ValNodePtr vnpal, SeqIdPtr seqsip, ValNodePtr fromp, Int2 nbseq, Int4 lens, ValNodePtr strands, Boolean trunc_emptyends)
{
  SeqAnnotPtr  sap1=NULL;
  ValNodePtr   tmp;

  if (vnpal!=NULL && vnpal->data.ptrvalue != NULL) {
     tmp = (ValNodePtr) vnpal->data.ptrvalue;
     sap1 = LocalAlignToSeqAnnotDimn (tmp, seqsip, fromp, nbseq, lens, NULL, FALSE);
  }
  return sap1;
}

static SeqAnnotPtr LocalAlignsToSeqAnnotDimn (ValNodePtr vnpal, SeqIdPtr seqsip, ValNodePtr fromp, Int2 nbseq, Int2 nbseg, Int4 lens, ValNodePtr strands, Boolean trunc_emptyends)
{
  SeqAnnotPtr  sap1 = NULL, 
               sap = NULL;
  SeqAlignPtr  salphead = NULL,
               salptmp;
  ValNodePtr   vnp, 
               tmp;
  SeqIdPtr     siplst;

  vnp = vnpal; 
  salphead = NULL;
  while (salphead == NULL && vnp != NULL) 
  {
     siplst = SeqIdDupList (seqsip);
     tmp = (ValNodePtr) vnp->data.ptrvalue;
     sap1 = LocalAlignToSeqAnnotDimn (tmp, siplst, fromp, nbseq, lens, NULL, FALSE);
     if (sap1!=NULL && sap1->data!=NULL)
        salphead = (SeqAlignPtr) sap1->data;
     vnp = vnp->next;
  }
  if (fromp!=NULL)
     ValNodeFree (fromp);
  salptmp = salphead;
  while (vnp!=NULL) 
  {
     fromp = get_lens_fromseqalign (salptmp);
     siplst = SeqIdDupList (seqsip);
     tmp = (ValNodePtr) vnp->data.ptrvalue;
     sap = LocalAlignToSeqAnnotDimn (tmp, siplst, fromp, nbseq, lens, NULL, FALSE);
     if (sap!=NULL && sap->data!=NULL) {
        salptmp->next = (SeqAlignPtr)sap->data;
        salptmp = salptmp->next;
     }
     vnp = vnp->next;
     if (fromp!=NULL)
        ValNodeFree (fromp);
  }
  SeqIdFree (seqsip);
  return sap1;
}  

/*************************************************
***   Sequence:
***       FastaRead
***
*************************************************/
extern SeqEntryPtr FastaRead (CharPtr path, Uint2 mol_type)
{
  Char         name[PATH_MAX];
  SeqEntryPtr  sep_list = NULL, sep = NULL, pre_sep = NULL;
  BioseqSetPtr bssp;
  BioseqPtr    bsp;
  FILE         *fpin;
  Int2         j = 0;

  if (path == NULL) {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }
     path = name;
  }
  if ( (fpin = FileOpen (path, "r")) == NULL)  {
     return NULL;
  }
  while ((sep = FastaToSeqEntry (fpin, (Boolean)ISA_na (mol_type) ) ) != NULL)
  {
     if (j == 0) sep_list = sep;
     else  pre_sep->next = sep;
     pre_sep = sep;
     j++;
  }
  FileClose(fpin);
  if ( j == 0 )  {
     return NULL;
  }
  else if (j == 1) {
     sep_list->choice = 1;
     return sep_list;
  }
  bssp = BioseqSetNew ();
  if ( bssp == NULL )
     return NULL;
  bssp->_class = 14;
  bssp->seq_set = sep_list;
  for (sep = sep_list; sep!=NULL; sep=sep->next) {
     bsp = (BioseqPtr) sep->data.ptrvalue;
     ObjMgrConnect (OBJ_BIOSEQ, (Pointer) bsp, OBJ_BIOSEQSET, (Pointer) bssp);
  }
  sep = SeqEntryNew ();
  if ( sep  == NULL )
     return NULL;
  sep->choice = 2;
  sep->data.ptrvalue = (Pointer) bssp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bssp, sep);
  return sep;
}

static SeqEntryPtr FastaReadAdvanced (CharPtr path, Uint2 mol_type, Int2 *seqnumber, Int2 *segnumber, SeqIdPtr PNTR sip, Int4 *lengthmax)
{
  Char         name[PATH_MAX];
  SeqEntryPtr  sep = NULL, 
               sep1 = NULL;
  BioseqSetPtr bssp;
  FILE         *fpin;
   
  if (path == NULL) {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }
     path = name;
  }
  if ( (fpin = FileOpen (path, "r")) == NULL)  {
     return NULL;
  }
  if (segnumber != NULL)
     sep = NewFastaRead (fpin, (Boolean)ISA_na (mol_type), TRUE, seqnumber, segnumber, sip, lengthmax);
  else 
     sep = NewFastaRead (fpin, (Boolean)ISA_na (mol_type), TRUE, NULL, NULL, NULL, NULL);
  FileClose (fpin);
  if (sep != NULL) {
     bssp = BioseqSetNew ();
     bssp->_class = 14;
     bssp->seq_set = sep;
     sep1 = SeqEntryNew ();
     sep1->choice = 2;
     sep1->data.ptrvalue = bssp;
     SeqMgrLinkSeqEntry (sep1, 0, NULL);
  }
  return sep1;
}

/*****************************************************
*** GapFastaRead
***    first read the sequences as FASTA
***    2d read the sequence text with the gaps (-)
***       the max length allocated for the char array
***       that is the max length of the sequences plus
***       a 1/2 of gaps.
***
******************************************************/
extern SeqEntryPtr GapFastaRead (CharPtr path, Uint2 mol_type)
{
  Char         name[PATH_MAX];
  SeqAnnotPtr  sap = NULL;
  SeqEntryPtr  sephead = NULL;
  SeqEntryPtr  sep = NULL;
  BioseqSetPtr bssp;
  BioseqPtr    bsp;
  ValNodePtr   vnp;
  SeqIdPtr     sip1 = NULL,
               sipnew = NULL,
               siptmp = NULL;
  Int4         lmax;
  Int2         nseq = 0,
               seqnumber = 0,
               segnumber;

  if (path == NULL) {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }
     path = name;
  }
  sephead = FastaReadAdvanced (path, mol_type, &seqnumber, &segnumber, &sip1, &lmax);
  if (sephead == NULL) {
     return NULL;
  }
  nseq=0;
  for (siptmp=sip1; siptmp!=NULL; siptmp=siptmp->next) {
     nseq++;
  }
  if (nseq != seqnumber) {
     ErrPostEx (SEV_ERROR, 0, 0, "The submission contains an error");
     return NULL;
  }
  if (matching_seqid (sip1)) {
     ErrPostEx (SEV_ERROR, 0, 0, "The submission contains identical sequence IDs");
     return NULL;
  } 
  vnp = ReadAlignmentToStrings (path, lmax, segnumber);
  if (segnumber > 1)
     sap=LocalAlignsToSeqAnnotDimn(vnp,sip1,NULL,seqnumber,segnumber, 0, NULL, FALSE);
  else  
     sap=LocalAlign1ToSeqAnnotDimn (vnp, sip1, NULL, seqnumber, 0, NULL, FALSE);
  if (sap!=NULL && sap->data !=NULL) {
     if (IS_Bioseq(sephead)) {
        bsp=(BioseqPtr)sephead->data.ptrvalue;
        bsp->annot = sap;
     }
     else if (IS_Bioseq_set(sephead)) {
        bssp = (BioseqSetPtr)sephead->data.ptrvalue;
        bssp->annot = sap;
     }
  }
  return sephead;
}

static Boolean seq_char (Char car, Char missingchar, Char gapchar)
{
  if (car == 'A') return TRUE;
  if (car == 'T') return TRUE;
  if (car == 'G') return TRUE;
  if (car == 'C') return TRUE;
  if (car == 'U') return TRUE;
  if (car == 'a') return TRUE;
  if (car == 't') return TRUE;
  if (car == 'g') return TRUE;
  if (car == 'c') return TRUE;
  if (car == 'u') return TRUE;
  if (car == missingchar) return TRUE;
  if (car == gapchar) return TRUE;
  if (car == '*') return TRUE;
  return FALSE;
}

static Boolean has_extrachar (CharPtr str, Char missingchar, Char gapchar)
{
  Int2     j;
  Boolean  ret = FALSE;
 
  if (str==NULL)
     return TRUE;
  if (*str=='\0' || *str=='\n')
     return TRUE;
  for (j=0; j<StrLen(str); j++) {
     if (str[j]!='\n' && str[j]!='\0' && str[j]!='\r' && str[j]!=' ') {
        if (!isdigit(str[j])) {
           if ( (StringChr ("ABCDGHKMNRSTUVWY", str[j])) == NULL &&
             (StringChr ("abcdghkmnrstuvwy", str[j])) == NULL &&
             str[j]!=gapchar && str[j] != missingchar)  {
              ret = TRUE;
              break;
           }
        }
     }
  }
  return ret;
}

static Char nexustoseq (Char car, Char missingchar, Char gapchar)
{
  if (car == ':')
     return ('-');
  if (car == '.')
     return ('-');
  if (car == missingchar)
     return ('N');
  if (car == gapchar)
     return ('-');
  if (isalpha (car))
     return car;
  return ('\0');
}

static Boolean ConvertPaupToFastaGap (CharPtr path, CharPtr tmpfile)
{
  FILE       *fp, *fpout;
  CharPtr    tmp;
  Char       str [MAXSTR];
  Char       str2 [MAXSTR];
  Char       gapchar = '-';
  Char       missingchar = '?';
  Char       car;
  Int4       strlens;
  Int4       lg_seq = 0;
  Int2       n_seq = 0;
  Int2       n_tmp = 0;
  Int4       j, j1, 
             k, 
             k1=0;
  Boolean    goOn, first_line;
 
  if ( (fp = FileOpen (path, "r")) == NULL) {
     return FALSE;
  }
  goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
  while (goOn) {
     if (! stringhasnotext (str)) {
        if (StringLen (str) > 0 && str [0] != '>')
           break;
     }   
     goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
  }
  if (!goOn) {
    FileClose(fp); 
    return FALSE;
  }
  while (goOn) {
        tmp = StringStr(str, "INTERLEAVE");
        if (tmp == NULL)
           tmp = StringStr(str, "interleave");
        if (tmp != NULL) {
           n_seq = 0;
           lg_seq= 0;
           ErrPostEx (SEV_ERROR, 0, 0, "This is a NEXUS interleave format"); 
           break;
        }
        tmp = StringStr(str, "GAP=");
        if (tmp == NULL)
           tmp = StringStr(str, "gap=");
        if (tmp != NULL) {
           while (*tmp!='\0' && *tmp!='\n' && *tmp!='=')
              tmp++;
           if (*tmp!='\0' && *tmp!='\n') 
              tmp++;
           while (*tmp!='\0' && *tmp!='\n' && *tmp==' ') 
              tmp++;
           if (*tmp!='\0' && *tmp!='\n')
              gapchar = *tmp;
        }
        tmp = StringStr(str, "MISSING=");
        if (tmp == NULL)
           tmp = StringStr(str, "missing=");
        if (tmp != NULL) {
           while (*tmp!='\0' && *tmp!='\n' && *tmp!='=')
              tmp++;
           if (*tmp!='\0' && *tmp!='\n')
              tmp++;
           while (*tmp!='\0' && *tmp!='\n' && *tmp==' ') 
              tmp++;
           if (*tmp!='\0' && *tmp!='\n')
              missingchar = *tmp;
        }
        if (n_seq == 0) {
           tmp = StringStr(str, "NTAX");
           if (tmp == NULL)
              tmp = StringStr(str, "ntax");
           if (tmp != NULL) {
              while (*tmp!='\0' && *tmp!='\n' && !isdigit (*tmp))
                 tmp++;
              if (*tmp!='\0' && *tmp!='\n')
                 n_seq = (Int2) atoi(tmp);
           }          
        }
        if (lg_seq == 0) {
           tmp = StringStr(str, "NCHAR");
           if (tmp == NULL)
              tmp = StringStr(str, "nchar");
           if (tmp != NULL) {
              while (*tmp!='\0' && !isdigit (*tmp))
                 tmp++;
              if (*tmp!='\0')
                 lg_seq = (Int4) atol(tmp);
           }
        }
        tmp = StringStr(str, "MATRIX");
        if (tmp == NULL)
           tmp = StringStr(str, "matrix");
        if (tmp!=NULL) {
           break;
        }   
        goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
  }   
  if (n_seq == 0 || lg_seq == 0) {
     FileClose(fp); 
     return FALSE;
  }
  while (goOn) {
     tmp = StringStr(str, "MATRIX");
     if (tmp == NULL)
        tmp = StringStr(str, "matrix");
     if (tmp != NULL)
        break;
     goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
  }
  if (!goOn) {
     FileClose(fp); 
     return FALSE;
  }
  if ( (fpout = FileOpen (tmpfile, "w")) == NULL) {
     FileClose(fp); 
     return FALSE;
  }
  goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
  first_line = TRUE;
  n_tmp = 0;
  k=0;
  while (goOn) {
     strlens = StringLen (str);
     if (strlens > 0) {
        if (str[0] == ';')
           break;
        if (has_extrachar (str, missingchar, gapchar)) {
           first_line = TRUE;
        }
        j1=j=0;
        while (j < strlens) 
        { 
           if (str[j]=='\0' || str[j] == '\n' || str[j] == '\r' ) {
              str[j1]='\0';
              break;
           }
           if (!first_line) {
              car = nexustoseq (str[j], missingchar, gapchar);
              if (car != '\0') {
                 str[j1] = car;
                 j1++;
              }
              j++;
           }
           else if (first_line) {
              if (str[j]!=' ') {
                 str[j1] = str[j];
                 j1++;
                 j++;
              }
              else {
                 while (str[j] == ' ')
                    j++;
                 k1=0;
                 while (str[j]!='\0' && str[j]!='\n' && j < strlens) {
                    car = nexustoseq (str[j], missingchar, gapchar);
                    if (car != '\0') {
                       str2[k1] = car;
                       k1++;
                    }
                    str[j] = ' ';
                    j++;
                 }
                 if (k1>0)
                    str2[k1] = '\0';       
              }
           }
        }
        strlens = StringLen (str);
        if (strlens > 0 && !stringhasnocharplus (str)) {
           if (!first_line && has_extrachar (str, missingchar, gapchar)) {
              first_line = TRUE; 
           }
           if (first_line) {
              if (strlens > 1) {
                 fprintf(fpout, ">%s\n", str);
                 first_line = FALSE; 
                 k=0;
                 n_tmp++;
                 if (k1 > 0) {
                    fprintf(fpout, "%s\n", str2);
                    k1 = 0;   
                    k += StringLen (str2);
                 }
              }
           }
           else {
              fprintf(fpout, "%s\n", str);
              k += strlens; 
              if (k >= lg_seq)  {
                 if (n_tmp == n_seq)
                    break;
                 first_line = TRUE;
              }
           }
        }
     }      
     k1=0;
     goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
  }    
  FileClose(fp);
  fprintf(fpout, "\n");
  FileClose(fpout);
  return TRUE;
}

/*************************************************
***   Sequence: 
***       IdRead
***
*************************************************/
extern ValNodePtr IdRead (CharPtr path)
{
  Char         name[PATH_MAX];
  FILE        *fp;
  ValNodePtr   sqloc_list = NULL;
  SeqLocPtr    slp = NULL,
               slpnew = NULL;
  BioseqPtr    bsp;
  SeqIdPtr     sip = NULL;
  SeqIntPtr    sit = NULL;
  Char         str [256];
  Int4         pos, 
               from ,to;
  Boolean      goOn;
  Boolean      stop;
  CharPtr      ptr;
  CharPtr      chptr;
  Char         ch;

  Uint2        choice = SEQID_GI;

  if (path == NULL) {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }
     path = name;
  }
  if ( (fp = FileOpen (path, "r")) == NULL)  {
      return NULL;
  }
  if (fp != NULL) {
    str [0] = '\0';
    pos = ftell(fp);
    goOn = (Boolean)(fgets (str, sizeof (str), fp) != NULL);
    if (!goOn) {
      FileClose (fp);
      return NULL; 
    }
    if (StringLen(str) == 0) {
      FileClose (fp);
      return NULL; 
    }
    while (goOn) {
      ptr = str;
      ch = *ptr;
      while (ch != '\n' && ch != '\r' && ch != '\0') {
        ptr++;
        ch = *ptr;
      }
      *ptr = '\0';
      if (str [0] == '&') {
        goOn = FALSE;
      } else if (str [0] != '?') {
              if (str [0] == '>')           
                ptr = str + 1;
              else 
                ptr = str;
              while (*ptr == ' ' || *ptr == '?') {
                ptr++;
              }
              stop = FALSE;
              if (*ptr == '"') {
                ptr++;
                chptr = StringChr (ptr, '"');
              } else {
                chptr = StringChr (ptr, ' ');
              }
              if (chptr != NULL) {
                *chptr = '\0';
                chptr++;
                if (choice != 0) ptr = check_seqid (choice, ptr);
                sip = MakeSeqID (ptr);
              } else if (*ptr != '\0') {
                if (choice != 0) ptr = check_seqid (choice, ptr);
                sip = MakeSeqID (ptr);
                stop = TRUE;
              } else {
                sip = MakeSeqID ("lcl|unknown");
                stop = TRUE;
              }
              from = to = -1;
              if ( !stop ) {
               from = to = -1;
               if ( *chptr != '\0' ) {
                  while ( *chptr != '\0' && *chptr == ' ' ) chptr++;
               }
               if ( *chptr != '\0' ) {
                  ptr = chptr;
                  while ( *chptr != ' ' && *chptr != '\0' ) chptr++;
                  *chptr = '\0';
                  chptr++;
                  from = (Int4) atoi (ptr);
               }
               if ( from >= 0 ) {
                  while ( *chptr != '\0' && *chptr == ' ' ) chptr++;
               }
               if ( *chptr != '\0' ) {
                  ptr = chptr;
                  while ( *chptr != ' ' && *chptr != '\0' ) chptr++;
                  *chptr = '\0';
                  chptr++;
                  to = (Int4) atoi (ptr);
               }
              }
              if (from >= 0 && to >= 0) {
                 slp = SeqLocIntNew (from, to, Seq_strand_plus, sip);
                 ValNodeAddPointer (&sqloc_list, 0, (Pointer) slp);
              } else {
                 bsp = BioseqLockById (sip);
                 if (bsp!=NULL) {
                    slp = SeqLocIntNew (0, bsp->length-1, Seq_strand_plus, sip);
                    BioseqUnlock (bsp);
                    ValNodeAddPointer (&sqloc_list, 0, (Pointer) slp);
                 }
              }
      }
      str [0] = '\0';
      pos = ftell(fp);
      goOn = (Boolean) (goOn && (fgets (str, sizeof (str), fp) != NULL));
    }
  }
  FileClose(fp);
  return sqloc_list;
}

/*******************************************************
*** AsnReadForSalsa
***   copied from Jonathan's code
***   without the following lines:
***
          rsult = SeqEntryNew ();
          if (rsult != NULL) {
            rsult->choice = sep->choice;
            rsult->data.ptrvalue = sep->data.ptrvalue;
            sep->data.ptrvalue = NULL;
            if (datatype == OBJ_SEQSUB) {
              SeqSubmitFree ((SeqSubmitPtr) dataptr);
            } else {
              SeqEntryFree (sep);
            }
            if (!ObjMgrRegister (OBJ_SEQENTRY, (Pointer) rsult))
               rsult = SeqEntryFree (rsult);
          }  
***********************************************************/
extern SeqEntryPtr AsnReadForSalsa (CharPtr path)
{
  Char         name[PATH_MAX];
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr;
  Uint2         datatype;
  Uint2         entityID;
  SeqEntryPtr   rsult;
  SeqEntryPtr   sep;

  rsult = NULL;
  if (path == NULL) {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }
     path = name;
  }
  if (path != NULL && path [0] != '\0') {
    dataptr = ObjMgrGenericAsnTextFileRead (path, &datatype, &entityID);
    if (dataptr != NULL && entityID > 0) {
      if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
          datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {
        sep = GetTopSeqEntryForEntityID (entityID);
        if (sep == NULL) {
          sep = SeqEntryNew ();
          if (sep != NULL) {
            if (datatype == OBJ_BIOSEQ) {
              bsp = (BioseqPtr) dataptr;
              sep->choice = 1;
              sep->data.ptrvalue = bsp;
              SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
            } else if (datatype == OBJ_BIOSEQSET) {
              bssp = (BioseqSetPtr) dataptr;
              sep->choice = 2;
              sep->data.ptrvalue = bssp;
              SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
            } else {
              sep = SeqEntryFree (sep);
            }
          }  
          sep = GetTopSeqEntryForEntityID (entityID);
        }
        if (sep != NULL) {
           rsult = sep;
        }
      }
    }
  }
  return rsult;
}

/*************************************************
***  Sequence
***     BioseqSetFileWrite
***
*************************************************/
extern Int2 BioseqSetFileWrite (BioseqSetPtr bssp)
{
  AsnIoPtr     aip;
  AsnTypePtr   atp;
  AsnModulePtr amp;
  Char         path [PATH_MAX]; 

  if ( bssp == NULL ) {
    return 0;
  }
  if (! GetInputFileName (path, PATH_MAX,"","TEXT"))  {
    return 0;
  }
  amp = AsnAllModPtr ();
  atp = AsnTypeFind (amp,"Bioseq-set");
  if ((aip = AsnIoOpen (path,"w")) == NULL) {
    return 0;
  }
  if ( ! BioseqSetAsnWrite ( bssp, aip, atp ) ) {
  }
  aip = AsnIoClose (aip);
  return 1;
}

/*************************************************
***  Sequence
***     BioseqFileWrite
***
*************************************************/
extern Int2 BioseqFileWrite (BioseqPtr bsp)
{
  AsnIoPtr     aip;
  AsnTypePtr   atp;
  AsnModulePtr amp;
  Char         path [PATH_MAX]; 

  if ( bsp == NULL ) {
    return 0;
  }
  if (! GetInputFileName (path, PATH_MAX,"","TEXT"))  {
    return 0;
  }
  amp = AsnAllModPtr ();
  atp = AsnTypeFind (amp,"Bioseq");
  if ((aip = AsnIoOpen (path,"w")) == NULL) {
    return 0;
  }
  if ( ! BioseqAsnWrite ( bsp, aip, atp ) ) {
  }
  aip = AsnIoClose (aip);
  return 1;
}

/*************************************************
***  Sequence
***     seqentry_read
***     seqentry_write
***
*************************************************/
extern SeqEntryPtr seqentry_read (CharPtr path)
{
  Char         name[PATH_MAX];
  AsnIoPtr     aip;
  AsnTypePtr   atp;
  AsnModulePtr amp;
  SeqEntryPtr  sep;

  if (path == NULL )
  {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }
     path = name;
  }
  amp = AsnAllModPtr ();
  atp = AsnTypeFind (amp,"SeqEntry");
  if ((aip = AsnIoOpen (path,"r")) == NULL) {
    return NULL;
  }
  sep = SeqEntryAsnRead ( aip, atp );
  aip = AsnIoClose (aip);
  return sep;
}

extern Boolean seqentry_write (SeqEntryPtr sep, CharPtr path)
{
  Char         name[PATH_MAX];
  AsnIoPtr     aip;
  AsnTypePtr   atp;
  AsnModulePtr amp;

  if ( sep == NULL ) {
    return 0;
  }
  if (path == NULL )
  {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return 0;
     }
     path = name;
  }
  amp = AsnAllModPtr ();
  atp = AsnTypeFind (amp,"SeqEntry");
  if ((aip = AsnIoOpen (path,"w")) == NULL) {
    return 0;
  }
  if ( ! SeqEntryAsnWrite ( sep, aip, atp ) ) {
  }
  aip = AsnIoClose (aip);
  return 1;
}

/*************************************************
***  Sequence
***     WriteSeqToFasta
***
*************************************************/
static void selalignnode_tofasta (CharPtr path, ValNodePtr anp_list, Int2 Width_Page, Uint2 mol_type, Boolean writeseq)
{
  AlignNodePtr     anp;
  SelStructPtr     ssp;
  SeqLocPtr        slp;
  SeqPortPtr       spp;
  Char             buffer[128];
  CharPtr          str, str2;
  Int4             j;
  FILE             *fout;
  Char    strLog[128];

  if ( (fout = FileOpen (path, "w")) == NULL) {
    return;
  }
  ssp = ObjMgrGetSelected();  
  for (; ssp != NULL; ssp = ssp->next) 
     if ( checkssp_for_editor (ssp) ) 
     {
         anp = (AlignNodePtr) AlignNodeFind (anp_list, ssp->entityID, ssp->itemID, ssp->itemtype);
         if ( anp != NULL ) {
                slp = CollectSeqLocFromAlignNode (anp);
                SeqIdWrite ( SeqLocId (slp), strLog, PRINTID_FASTA_LONG, 120);
                str = strLog;
                str2 = StringStr(str,"lcl|");
                if ( str2 != NULL ) str+=4;
                fprintf (fout, "> %s   length %ld  from %ld to %ld \n", strLog,
                            (long) (SeqLocStop (slp) - SeqLocStart (slp)), 
                            (long) SeqLocStart (slp), (long) SeqLocStop (slp) );
                if (writeseq)
                {
                   if ( mol_type == Seq_mol_aa )
                       spp = SeqPortNewByLoc (slp, Seq_code_ncbieaa);
                   else
                       spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
                   j = 0;
                   while ( j < (SeqLocStop (slp) - SeqLocStart (slp)) ) 
                   {
                       j +=ReadBufferFromSep (spp, buffer, j, j +Width_Page, 0);
                       fprintf(fout, "%s\n", buffer); 
                   }
                   SeqPortFree (spp);
                }
         }
     }
  FileClose (fout);
}

/*************************************************
***  Sequence
***     BioseqSetToFasta 
***
*************************************************/
extern void EditBioseqToFasta (BioseqPtr bsp, FILE *fout, Boolean is_na, Int4 from, Int4 to)
{
  SeqLocPtr        slp;
  SeqPortPtr       spp;
  Char             buffer[128];
  Char             str [128];
  Int4             txt_out;
  Int4             Width_Page = 60;
  Int4             j;

  if (bsp == NULL) 
     return;
  if (fout == NULL)
     return; 
  SeqIdWrite (SeqIdFindBest(bsp->id, 0), str, PRINTID_FASTA_LONG, 120);
  if (from < 0) 
     from = 0;
  if (to < 0) 
     to =  bsp->length-1;
  fprintf (fout, ">%s   (%ld - %ld)\n", str, (long)(from+1), (long)(to+1));
  slp = SeqLocIntNew (from, to, Seq_strand_plus, SeqIdFindBest(bsp->id, 0));
  if ( bsp->mol == Seq_mol_aa )
     spp = SeqPortNewByLoc (slp, Seq_code_ncbieaa);
  else
     spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
  j = 0;
  while ( j < SeqLocStop (slp) - SeqLocStart (slp) +1)
  {
     txt_out = ReadBufferFromSep (spp, buffer, j, j +Width_Page, 0);
     if (txt_out == 0) break;
     j += txt_out;
     fprintf(fout, "%s\n", buffer);
  }
  SeqPortFree (spp);
  return;
}

static Boolean BioseqSetToFasta (BioseqSetPtr bssp, Boolean firstout)
{
  SeqEntryPtr      sep;
  BioseqPtr        bsp;
  FILE             *fout;
  Int2             count = 0;

  if ( (fout = FileOpen ("ffile", "w")) != NULL) {
     sep = bssp->seq_set;
     while (sep != NULL)
     {
         count++;
         if (count == 1 && !firstout) {}
         else {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            EditBioseqToFasta (bsp, fout, (Boolean)ISA_na(bsp->mol), -1, -1);
         }
         sep = sep->next;
     }
     FileClose (fout);
     return TRUE;
  }
  return FALSE;
}

/*************************************************
***  SeqAnnot
***
*************************************************/
extern SeqAnnotPtr seqannot_read (CharPtr path)
{
  AsnIoPtr     aip;
  AsnTypePtr   atp;
  AsnModulePtr amp;
  SeqAnnotPtr  sap_head = NULL, pre_sap = NULL, sap;
  Boolean      first =  TRUE;

  amp = AsnAllModPtr ();
  atp = AsnTypeFind (amp,"Seq-annot");
  if ((aip = AsnIoOpen(path,"r")) == NULL) {
         return NULL;
  }
  while ((atp = AsnReadId (aip, amp, atp)) != NULL) 
  {
         sap = SeqAnnotNew ();
         if ((sap = SeqAnnotAsnRead (aip, atp)) == NULL) 
         {
                return NULL;
         }
         if ( first ) {
                sap_head = sap;
                first = FALSE;
         } else pre_sap->next = sap;
         pre_sap = sap;
  }
  aip = AsnIoClose (aip);
  return sap_head;
}

/*************************************************
***  SeqAnnot
***
*************************************************/
extern Int2 seqannot_write (SeqAnnotPtr sap, CharPtr path)
{
  Char         name[PATH_MAX];
  AsnIoPtr     aip;
  AsnTypePtr   atp;
  AsnModulePtr amp;
  SeqAnnotPtr  saptmp;
  SeqAlignPtr  salp;

  if ( sap == NULL ) {
         return 0;
  }
  if ( ( salp = (SeqAlignPtr) sap->data ) == NULL ) {
         return 0;
  }
  if ( salp->segtype == 4 ) {
         saptmp = SeqAnnotBoolSegToDenseSeg (sap); 
  }
  else   saptmp = sap;
  if (path == NULL) {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return 0;
     }
     path = name;
  }
  amp = AsnAllModPtr ();
  atp = AsnTypeFind (amp,"Seq-annot");
  if ((aip = AsnIoOpen (path, "w")) == NULL) {
         return 0;
  }
  while ( sap != NULL ) {
         if ( ! SeqAnnotAsnWrite ( sap, aip, atp ) ) {
                break;
         }
         sap = sap->next;
  }
  AsnIoReset(aip);
  aip = AsnIoClose (aip);
  if ( salp->segtype == 4 ) CompSeqAnnotFree (saptmp);
  return 1;
}

extern void seqalign_write (SeqAlignPtr salp, CharPtr path)
{
  SeqAnnotPtr  sap;

  if (salp!=NULL) {
     sap = SeqAnnotNew (); 
     if (sap != NULL) {
        sap->type = 2; 
        sap->data = (Pointer) salp; 
        seqannot_write (sap, path); 
        sap->data = NULL; 
        sap = SeqAnnotFree (sap); 
     }
  }
}

/*************************************************
***  SeqAnnot
***     SeqAnnotFileRead
***
*************************************************/
extern SeqAnnotPtr SeqAnnotFileRead (void)
{
  Char         path [PATH_MAX];  /* path+filename */ 
  SeqAnnotPtr  sap_head = NULL;
  SeqAlignPtr  salp;
  DenseSegPtr  dsp;
  DenseDiagPtr ddp;
   
  if (!GetInputFileName (path, PATH_MAX, "", "TEXT")) 
     return NULL;
  sap_head = seqannot_read (path);

  if ( sap_head->type == 2 ) 
  {
         salp = (SeqAlignPtr) sap_head->data;
         if ( salp == NULL ) return NULL;
         if ( salp->dim < 2 || salp->type == 0 )
                return NULL;
         if ( salp->segtype == 1 ) {
                ddp = (DenseDiagPtr) salp->segs;
                if ( ddp->dim < 2 ){
                       SeqAnnotFree (sap_head);
                       return NULL;
                }
         } 
         else if ( salp->segtype == 2 ) {
                dsp = (DenseSegPtr) salp->segs;
                if ( dsp->dim < 2 ){
                       SeqAnnotFree (sap_head);
                       return NULL;
                }
         }
         else return NULL;
  }
  else {
         return NULL; 
  }
  return (sap_head);
}

/*************************************************
***
***  Alignment
***              
*************************************************/
static Boolean seq_line (CharPtr str)
{
  Int4 lens;
  Int4 val1, val2, j;

  if (str != NULL) {
     lens = StringLen(str);
     val1 = 0;
     val2 = 0;
     for (j = lens; j > 0; j--) {
        str[j] = TO_UPPER(str[j]);
        if (str[j] >= 'A' && str[j] <= 'Z') {
           val1++;
           if (str[j]=='A' || str[j]=='C' || str[j]=='T' || str[j]=='G')
              val2++;
        }
     }
     if (val2 > (2*val1/3))
        return TRUE;
  }
  return FALSE;
}

static ValNodePtr ReadLocalAlign (CharPtr path, Int2 align_format, Int2 n_seq)
{
  FILE       *fp;
  ValNodePtr seqvnp = NULL, vnp;
  CharPtr    tmp,
             tmp1;
  Char       str [MAXSTR];
  Int4 PNTR  lgseq;
  Int4       lmax;
  Int4       strlens;
  Int2       i_seq, j;	
  Int2       leftmargin;
  Int4       lg_seq = 0;
  int        val1;
  long       val2;
  Boolean    found_seq;
  Boolean    goOn;
  Boolean    first;
  Char       gapchar = '-';
  Char       missingchar = '?';

  if ( (fp = FileOpen (path, "r")) == NULL) {
         return NULL;
  }
  goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
  while (goOn) {
     if (! stringhasnotext (str)) {
        if (StringLen (str) > 0 && str [0] != '>') 
           break;
     }
     goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
  }
  if (align_format == SALSAA_GCG){
     n_seq = 1;
     goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
     while (goOn) {
        n_seq++;
        goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
     } 
     FileClose(fp);
     fp = FileOpen (path, "r");
     goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
     while (goOn) {
        if (! stringhasnotext (str)) {
           if (StringLen (str) > 0 && str [0] != '>')
              break;
        }   
        goOn=(Boolean)(fgets (str, sizeof (str), fp) != NULL);
     }
     leftmargin = SALSAA_GCG;
  }
  else if (align_format == SALSA_NEXUS) {
     found_seq = FALSE;
     lg_seq = 0;
     n_seq = 0;
     while (goOn) {
        if (n_seq == 0) {
           tmp = StringStr(str, "NTAX");
           if (tmp == NULL)  
              tmp = StringStr(str, "ntax");
           if (tmp != NULL) {   
              while (tmp!='\0' && !isdigit (*tmp)) 
                 tmp++;
              if (tmp!='\0') 
                 n_seq = (Int2) atoi(tmp); 
           }
        }
        if (lg_seq == 0) {
           tmp = StringStr(str, "NCHAR");
           if (tmp == NULL)  
              tmp = StringStr(str, "nchar");
           if (tmp != NULL) {   
              while (tmp!='\0' && !isdigit (*tmp)) 
                 tmp++;
              if (tmp!='\0') 
                 lg_seq = (Int4) atol(tmp); 
           }
        }
        tmp = StringStr(str, "GAP=");
        if (tmp == NULL)
           tmp = StringStr(str, "gap=");
        if (tmp != NULL) {
           while (*tmp!='\0' && *tmp!='\n' && *tmp!='=')
              tmp++;
           if (*tmp!='\0' && *tmp!='\n') 
              tmp++;
           while (*tmp!='\0' && *tmp!='\n' && *tmp==' ') 
              tmp++;
           if (*tmp!='\0' && *tmp!='\n')
              gapchar = *tmp;
        }
        tmp = StringStr(str, "MISSING=");
        if (tmp == NULL)
           tmp = StringStr(str, "missing=");
        if (tmp != NULL) {
           while (*tmp!='\0' && *tmp!='\n' && *tmp!='=')
              tmp++;
           if (*tmp!='\0' && *tmp!='\n')
              tmp++;
           while (*tmp!='\0' && *tmp!='\n' && *tmp==' ') 
              tmp++;
           if (*tmp!='\0' && *tmp!='\n')
              missingchar = *tmp;
        }
        if (n_seq>0 && lg_seq>0 && seq_line (str)) {
           if (seq_char(str[0], missingchar, gapchar) 
           && seq_char(str[1], missingchar, gapchar) 
           && seq_char(str[2], missingchar, gapchar)) {
              leftmargin = 0;
              found_seq = TRUE;
              break;
           }
           for (leftmargin = 0; leftmargin<MAXSTR-1; leftmargin++) {
              if (str[leftmargin] == ' ' 
              && seq_char(str[leftmargin+1], missingchar, gapchar)) {
                 found_seq = TRUE;
                 break;
              }
           }
           break;
        }
        goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
     }
     if (!found_seq) 
        n_seq = 0;
  }
  else if (align_format == SALSA_PHYLIP) {
     if (sscanf (str, "%d %ld", &val1, &val2) == 2) {
        n_seq = (Int2) val1;
        lg_seq = (Int4) val2;
     }
     goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
     leftmargin = SALSA_PHYLIP;
  }
  else if (align_format == SALSA_CLUSTALV) {
     if (n_seq == 0) {
        FileClose(fp);
        return NULL;
     }
     for ( j =0; j < 4; j++) 
         goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
     leftmargin = SALSA_CLUSTALV;
  }
  if (n_seq == 0) {
     FileClose(fp);
     return NULL;
  }
  if (lg_seq == 0) 
     lmax = LENGTHMAX;
  else 
     lmax = lg_seq;
  
  for (j = 0; j < n_seq; j++) {
     tmp = (CharPtr) MemNew((size_t) ((lmax + 1) * sizeof(Char)));
     for (strlens = 0; strlens < lmax; strlens++) 
         tmp [strlens] = ' ';
     tmp [lmax] = '\0';
     ValNodeAddPointer (&seqvnp, 0, (Pointer)tmp);
  }
  lgseq = (Int4Ptr) MemNew((size_t) ((n_seq + 1) * sizeof(Int4)));
  for (j = 0; j < n_seq; j++) lgseq [j] = 0; 
  
  tmp1 = (CharPtr) seqvnp->data.ptrvalue;
  i_seq = 0;
  vnp = seqvnp;
  first = TRUE;
  while (goOn)
  {
     strlens = StringLen (str);
     if (strlens > 0) {
           if (str[0] == ';') 
              break;
           if (! stringhasnocharplus (str) && str[0]!='>') 
           {
                tmp = (CharPtr) vnp->data.ptrvalue;
                for (j = leftmargin; j < strlens && lgseq [i_seq] <= lmax; j++) 
                { 
                   if (str[j] == '\n' || str[j] == '\r' ) break;
                   str[j] = TO_UPPER (str[j]);  
                   if (str[j] == gapchar)
                      str[j] = '-';
                   else if (str[j] == ':')
                      str[j] = '-';
                   else if (str[j] == '.') 
                   { 
                      if (align_format == SALSA_PHYLIP && i_seq != 0)
                      {
/**
                         if (tmp1[])
                            str [j]= tmp1[lgseq[i_seq]];
                         else
**/
                            str [j] = '-'; 
                      }
                      else
                         str[j] = '-';
                   }
                   else if (str[j] == '?')
                      str[j] = 'N';
                   if ((str[j] >= 'A' && str[j] <= 'Z') || str[j]=='*' || str[j] == '-') { 
                      tmp [lgseq[i_seq]] = str[j]; 
                      ++lgseq [i_seq];
                   }
                }
                ++i_seq;
                if (i_seq == n_seq) {
                   i_seq = 0;
                   vnp = seqvnp;
                   if (align_format == SALSA_PHYLIP && first) {
                      leftmargin=0; 
                      first = FALSE;
                   }
                } 
                else vnp = vnp->next;
            }
     }
     goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
  }
  FileClose(fp);
  for (lmax = 0, j = 0; j < n_seq; j++) 
     if (lgseq[j] > lmax) 
        lmax = lgseq[j];
  for (vnp = seqvnp, j = 0 && vnp != NULL; j < n_seq; j++, vnp = vnp->next) 
  {
     tmp = (CharPtr) vnp->data.ptrvalue;
     tmp [lmax] = '\0';
  }
  if (lg_seq == 0 ) 
     lg_seq = lmax;
  else if (lmax < lg_seq) 
  {
     if (lg_seq < LENGTHMAX )
        Message(MSG_OK, "Length in file %d != alignment length %d", (int) lg_seq, (int) lmax);
     lg_seq = lmax;
  }
  return seqvnp;  
}

/*************************************************
***  Alignment
***              
*************************************************/
static SeqIdPtr ReadLocalName (CharPtr path, Int2 nbseq, Int2 format)
{
  FILE       *fp;
  SeqIdPtr   sip1 = NULL,
             sipnew = NULL, siptmp;
  Char       str [MAXSTR];
  Int2       leftmargin;
  Int2       j;
  int        i_seq = 0;	
  Boolean    found_seq;
  Boolean    goOn;
  Char       gapchar = '-';
  Char       missingchar = '?';

  if ( (fp = FileOpen (path, "r")) == NULL) {
    return NULL;
  }
  goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
  while (goOn) {
     if (! stringhasnotext (str)) {
        if (StringLen (str) > 0 && str [0] != '>') 
           break;
     }
     goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
  }
  if (format == SALSA_NEXUS) {
     found_seq = FALSE;
     while (goOn) {
        if (seq_line (str)) {
           for (j = 0; j<MAXSTR-1; j++) {
              if (str[j] == ' ' 
              && seq_char(str[j+1], missingchar, gapchar)) {
                 found_seq = TRUE;
                 break;
              }
           }
           break;
        }
        goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
     }
     if (!found_seq) 
        return NULL;
     leftmargin = j;
  }
  else if (format == SALSA_PHYLIP)  {
     goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
     leftmargin = SALSA_PHYLIP;
  }
  else if (format == SALSA_CLUSTALV) {
     for ( j =0; j < 4; j++) 
         fgets (str, MAXSTR, fp);
     leftmargin = SALSA_CLUSTALV;
  }
  else 
     leftmargin = format;
  while (goOn && i_seq < nbseq ) 
  {
         if ( StringLen (str) > 0 ) 
         {                        
                str [leftmargin] = '\0';
                for (j=leftmargin-1; j>0 && str[j] == ' '; j--) 
                   str[j] = '\0';
                sipnew = MakeSeqID (str);
                if (sip1 == NULL)
                   sip1 = sipnew;
                else
                   siptmp->next = sipnew;
                siptmp = sipnew;
         }
         goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
         ++i_seq;
  }
  FileClose(fp);
  if (matching_seqid (sip1))
     ErrPostEx (SEV_ERROR, 0, 0, "MATCHING SEQIDs");
  return sip1;
}

/*************************************************
***  Alignment
***      ReadAlignmentFunc
***              
*************************************************/
extern void ReadAlignView (BioseqPtr target_bsp, Uint2 import_format)
{
  ValNodePtr  vnp=NULL;
  SeqLocPtr   slp=NULL;
  SeqLocPtr   source_slp;
  SeqIdPtr    sip;
  SeqAlignPtr salp;
  BioseqPtr   source_bsp = NULL;
  SeqEntryPtr source_sep = NULL;

  if (target_bsp==NULL)
     return;
  switch (import_format) {
     case SALSA_FASTA:
        source_sep = FastaRead (NULL, Seq_mol_na);
        break;
     case SALSA_ASN1:
        source_sep = AsnReadForSalsa (NULL);
        break;
     default:
        break;
  }
  if (source_sep!=NULL) {
     if (IS_Bioseq_set(source_sep))
        source_sep = FindNucSeqEntry(source_sep);
     if (IS_Bioseq(source_sep))
     {
        sip = SeqIdFindBest(target_bsp->id, 0);
        slp = SeqLocIntNew (0, target_bsp->length - 1, Seq_strand_plus, sip);
        ValNodeAddPointer(&vnp, 0, (Pointer)slp);
        source_bsp = (BioseqPtr)source_sep->data.ptrvalue;
        source_slp = SeqLocIntNew (0, source_bsp->length-1, Seq_strand_plus, source_bsp->id);
        ValNodeAddPointer(&vnp, 0, (Pointer)source_slp);
        salp = SeqLocListToSeqAlign (vnp, PRGALIGNDEFAULT, NULL);
        if (salp != NULL) {
           LaunchAlignViewer (salp);
        }
     }
  }
  return;
}


extern SeqEntryPtr ReadAlignmentFunc (CharPtr path, Uint1 mol_type, Uint1 format, Int2 n_seq, Boolean save_seqentry, Boolean save_sap, SeqIdPtr sqloc_list)
{
  SeqEntryPtr  sep_list, 
               sep, 
               pre_sep = NULL;
  BioseqSetPtr bssp;
  BioseqPtr    bsp;
  ValNodePtr   seqvnp , vnp;
  SeqIdPtr     seqsip , sip; 
  SeqAnnotPtr  sap;
  Int4         lens;
  Int2         k;

         /************************************
         **  read sequences into Charptr seqvnp 
         **  read names into Charptr seqsip 
         *************************************/
  seqvnp = ReadLocalAlign (path, format, n_seq);
  if (seqvnp == NULL)
  {
         ValNodeFree (seqvnp);
         return NULL;
  }
  k = 0;
  for (vnp=seqvnp; vnp!=NULL; vnp=vnp->next) k++;
  if (n_seq == 0)
     n_seq = k;
  else 
     if (k != n_seq) {
         ValNodeFree (seqvnp);
         return NULL;
     }

  if (sqloc_list != NULL)
         seqsip = sqloc_list;
  else 
         seqsip = ReadLocalName (path, n_seq, format); 

         /************************************
         **  sequences + names into SeqEntry    
         *************************************/
  if ( save_seqentry )
  {
         bssp = BioseqSetNew ();
         if (bssp == NULL) 
         {
                ValNodeFree (seqvnp);
                return NULL;
         }
         bssp->_class = 14;     
         vnp = seqvnp;
         sip = seqsip;
         for (k = 0; k < n_seq && vnp!=NULL && sip!=NULL; k++, vnp=vnp->next) 
         {
                lens = (Int4) StringLen ((CharPtr) vnp->data.ptrvalue);
                sep = StringToSeqEntry ((CharPtr) vnp->data.ptrvalue, sip, lens, mol_type);
                if (sep != NULL) {
                       if (pre_sep == NULL) bssp->seq_set = sep;
                       else pre_sep->next = sep;
                       pre_sep = sep;
                }
                sip = sip->next;
         }
         sep_list = SeqEntryNew ();
         if ( sep_list  == NULL ) {
                ValNodeFree (seqvnp);
                return NULL;
         }
         sep_list->choice = 2;
         sep_list->data.ptrvalue = (Pointer) bssp; 
         SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bssp, sep_list);

         for (sep = bssp->seq_set; sep!=NULL; sep=sep->next) {
                bsp = (BioseqPtr) sep->data.ptrvalue;
                ObjMgrConnect (OBJ_BIOSEQ, (Pointer) bsp, OBJ_BIOSEQSET, (Pointer) bssp);
         }
  }

         /*********************************
         **  alignment annot into sap 
         *********************************/
  if ( save_sap )
  {
         sap = LocalAlignToSeqAnnotDimn (seqvnp, seqsip, NULL, n_seq, 0, NULL, FALSE);
         if ( sap==NULL ) {
                ValNodeFree (seqvnp);
                return NULL;
         }
         bssp->annot = sap;
  }
  ValNodeFree (seqvnp);
  return sep_list;
}


extern SeqEntryPtr ReadLocalAlignment (Uint1 format, CharPtr path)
{
  SeqEntryPtr sep = NULL;
  Char        name [PATH_MAX];
  Char        tmpfile [PATH_MAX];
 
  if (path == NULL)
  {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }   
     path = name;
  }
  switch (format) 
  {
     case SALSA_FASTGAP:
            sep = GapFastaRead (path, Seq_mol_na);
            break;
     case SALSA_PAUP:
            TmpNam (tmpfile);
            if (ConvertPaupToFastaGap (path, tmpfile) ) 
            {
               sep = GapFastaRead (tmpfile, Seq_mol_na);
               FileRemove (tmpfile);
            }
            break;
     case SALSA_PHYLIP:
     case SALSA_NEXUS:
            sep = ReadAlignmentFunc (path, Seq_mol_na, format, 0, TRUE, TRUE, NULL);
            break;
     case SALSA_FASTA:
            sep = FastaReadAdvanced (path, Seq_mol_aa, NULL, NULL, NULL, NULL); 
            break;
     case SALSAA_FASTGAP:
            sep = GapFastaRead (path, Seq_mol_aa);
            break;
     case SALSAA_PHYLIP:
     case SALSAA_NEXUS:
            format = 11;
     case SALSAA_GCG:
            sep = ReadAlignmentFunc (path, Seq_mol_aa, format, 0, TRUE, TRUE, NULL);
            break;
     default:
            break;
  }
  return sep;
}

extern SeqEntryPtr ReadAnyAlignment (Boolean is_prot, CharPtr path)
{
  FILE        *fp;
  SeqEntryPtr sep = NULL;
  Char        name [PATH_MAX];
  CharPtr     tmp;
  Char        str [MAXSTR];
  Boolean     goOn;
  int         val1;
  long        val2;
 
  if (path == NULL)
  {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return NULL;
     }   
     path = name;
  }
  if ( (fp = FileOpen (path, "r")) == NULL) {
    return NULL;
  }
  goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
  while (goOn) {
     if (! stringhasnotext (str)) {
        if (StringLen (str) > 0)
           break;
     }
     goOn = (Boolean)(fgets (str, MAXSTR, fp) != NULL);
  }
  FileClose (fp);     

  if (str[0] == '>')
  {
     if (is_prot)
        sep = ReadLocalAlignment (SALSAA_FASTGAP, path);
     else 
        sep = ReadLocalAlignment (SALSA_FASTGAP, path);
     return sep;
  }
  tmp = StringStr(str, "NEXUS");
  if (tmp == NULL)
     tmp = StringStr(str, "nexus");
  if (tmp)
  {
     if (is_prot)
        sep = ReadLocalAlignment (SALSAA_NEXUS, path);
     else
        sep = ReadLocalAlignment (SALSA_NEXUS, path);
     return sep;
  }
  if (sscanf (str, "%d %ld", &val1, &val2) == 2) {
     if (val1 > 0 && val2 > 0)
     {
        if (is_prot)
           sep = ReadLocalAlignment (SALSAA_PHYLIP, path);
        else 
           sep = ReadLocalAlignment (SALSA_PHYLIP, path);
     } 
     return sep;
  }
  return sep;
}

extern void showtextalign_fromalign (SeqAlignPtr salp, CharPtr path, FILE *fp)
{
  Char        name [PATH_MAX];
  SeqAnnotPtr sap;
  Int4        line = 80;     
  Uint4       option = TXALIGN_MASTER;
  Boolean     do_close = TRUE;

  if (salp == NULL)
     return;
  if (path == NULL && fp == NULL) 
  {
     if (! GetInputFileName (name, PATH_MAX,"","TEXT"))  {
        return;
     }
     path = name;
  }
  if (path != NULL && fp == NULL) {
     fp = FileOpen (path, "w");
  }
  else
     do_close = FALSE;
  if (fp != NULL) {
     sap = SeqAnnotForSeqAlign (salp);  
     ShowTextAlignFromAnnot (sap, line, fp, NULL, NULL, option, NULL, NULL, NULL);
     if (do_close) {
        FileClose(fp);
     }
     sap->data = NULL;
     SeqAnnotFree (sap);
  } 
}

/***********************************************************
***
*** Import functions:
***    from any file (ASN.1, FASTA, gi/acc#), 
***       calls ReadAsnFastaOrFlatFile 
***    Download from Entrez
***       copy of FetchFromNet from sequin2.c (JK)
***
************************************************************/

extern SeqAlignPtr ImportFromFile (EditAlignDataPtr adp)
{
  SeqAlignPtr      salp = NULL,
                   salp_original=NULL;
  SeqAnnotPtr      sap = NULL;
  SeqEntryPtr      sep = NULL;
  ValNodePtr       importslp = NULL,
                   sqloc = NULL;
  SeqLocPtr        slp;
  Boolean          new_seqalign=FALSE,
                   ok=FALSE,
                   replace_salp=FALSE;

  if (adp==NULL)
     return NULL;

  importslp = CCReadAnythingLoop (NULL, adp->seq_info);
  if (importslp != NULL) 
  {
     if (adp->sap_original != NULL) {
        salp_original = (SeqAlignPtr)(adp->sap_original->data);
        if (salp_original->dim==2 || is_dim2seqalign (salp_original))
           salp_original=salp_original;
        else {
           if (salp_original->dim == 1)
              replace_salp = TRUE;
           salp_original=NULL;
        }
     }
     if (salp_original)
        ok=SeqAlignSeqLocComp (salp_original, importslp); 
     if (!ok)
     {
        slp = (SeqLocPtr) adp->master.region;
        if (slp!=NULL) 
        {
           ValNodeAddPointer (&sqloc, 0, (Pointer)slp);
           sqloc->next = importslp;
           salp = SeqLocListToSeqAlign (sqloc, (Int2)adp->align_format, NULL);
           if (salp!=NULL) 
           {
              if (salp_original != NULL) {
                 salp = SeqAlignLink (salp_original, salp);
                 new_seqalign = TRUE;
              }
           }
           else {
              if (adp->align_format==PRG_BLAST || adp->align_format==PRGALIGNDEFAULT)
                 Message (MSG_OK, "Blast detected no sequence similarity and could not construct an alignment");
              else
                 Message (MSG_OK, "No significant similarity detected. No alignment produced");
           }
        }
        if (!new_seqalign && !replace_salp)
           salp = SeqAlignFree (salp);
     }
     else
        Message(MSG_OK, "Can not import a sequence already in the editor"); 
  }
  return salp;
}

/*------------------------------------------------------------*/
typedef struct fetchform {
  FORM_MESSAGE_BLOCK
  GrouP           accntype;
  TexT            accession;
  ButtoN          accept;
  EditAlignDataPtr adp;
  WindoW          editor_window;
} FetchForm, PNTR FetchFormPtr;

static void FetchFormMessage (ForM f, Int2 mssg)
{
  FetchFormPtr  ffp;

  ffp = (FetchFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    switch (mssg) {
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (ffp->appmessage != NULL) {
          ffp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void FetchTextProc (TexT t)
{
  FetchFormPtr  ffp;

  ffp = (FetchFormPtr) GetObjectExtra (t);
  if (ffp == NULL) return;
  if (TextHasNoText (t)) {
    SafeDisable (ffp->accept);
  } else {
    SafeEnable (ffp->accept);
  }
}

static Int4 CCAccessionToGi (CharPtr string, Int2 seqtype)
{
   CharPtr str;
   LinkSetPtr lsp;
   Int4 gi;

   str = MemNew (StringLen (string) + 10);
   sprintf (str, "\"%s\" [ACCN]", string);
   lsp = EntrezTLEvalString (str, seqtype, -1, NULL, NULL);
   MemFree (str);
   if (lsp == NULL) return 0;
   if (lsp->num <= 0) {
       LinkSetFree (lsp);
       return 0;
   }
   gi = lsp->uids [0];
   LinkSetFree (lsp);
   return gi;
}

static SeqAlignPtr align_this (SeqEntryPtr sep, SeqLocPtr master, SeqAnnotPtr sap, WindoW editor_window, EditAlignDataPtr adp)
{
  SeqAlignPtr  salp = NULL, 
               salp_original = NULL;
  ValNodePtr   vnp=NULL,
               vnp2=NULL;
  Int2         n;
  Boolean      ok = FALSE,
               new_seqalign = FALSE,
               replace_salp = FALSE;

  if (sep==NULL)
     return NULL;
  if (!IS_Bioseq(sep))
     return NULL;

  vnp = SeqEntryToSeqLoc (sep, &n, adp->mol_type); 
  if (vnp)
  {
     if (sap)
     {
        salp_original = (SeqAlignPtr)(sap->data);
        if (salp_original->dim==2 || is_dim2seqalign (salp_original))
           salp_original=salp_original;
        else {
           if (salp_original->dim == 1)
              replace_salp = TRUE;
           salp_original=NULL;
        }
     } 
     if (salp_original)
        ok=SeqAlignSeqLocComp (salp_original, vnp);
     if (!ok)
     { 
        if (master!=NULL)
        {
           ValNodeAddPointer(&vnp2, 0, (Pointer)master);
           vnp2->next=vnp;
           vnp=vnp2;
        }
        salp = SeqLocListToSeqAlign (vnp, adp->align_format, NULL);
        if (salp!=NULL) 
        {
           if (salp_original) 
           {
              salp = SeqAlignLink (salp_original, salp);
              new_seqalign = TRUE;
           }
        }
        else {
           if (adp->align_format==PRG_BLAST || adp->align_format==PRGALIGNDEFAULT)
              Message (MSG_OK, "Blast detected no sequence similarity and could not construct an alignment");
           else
              Message (MSG_OK, "No significant similarity detected. No alignment produced");
        }
        if (new_seqalign || replace_salp) {
           repopulate_panel (editor_window, adp, salp);
        }
        else
           salp = SeqAlignFree (salp);
     }
     else
        Message(MSG_OK, "Can not import a sequence already in the editor"); 
  }
  return salp;
}

static SeqEntryPtr SeqEntryNewForBioseq (BioseqPtr bsp)
{
  SeqEntryPtr new_sep;
  BioseqPtr   new_bsp;
  SeqLocPtr   slp;
  SeqIdPtr    new_sip;

  slp=SeqLocIntNew(0, bsp->length-1, Seq_strand_plus, bsp->id);
  new_sip = MakeNewProteinSeqId (slp, NULL);
  new_bsp=BioseqCopy(new_sip, bsp->id, 0, bsp->length-1, Seq_strand_plus, TRUE);              
  if (new_bsp) {
     new_sep=SeqEntryNew();
     new_sep->choice = 1;
     new_sep->data.ptrvalue=(Pointer)new_bsp;
     SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) new_bsp, new_sep);
  }
  ValNodeFree(slp);
  return new_sep;
} 

static void CCDownloadProc (ButtoN b)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  CharPtr       dbname;
  Uint2         entityID;
  FetchFormPtr  ffp;
  Boolean       idTypes [NUM_SEQID];
  SeqEntryPtr   sep=NULL;
  ValNodePtr    sip2=NULL;
  SeqIdPtr      sip;
  Char          str [32];
  Int4          uid;
  ForM          w;
  Boolean       is_prot,
                is_newbsp=FALSE;

  ffp = (FetchFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  w = ffp->form;
  Hide (w);
  Update ();
  GetTitle (ffp->accession, str, sizeof (str));
  if (StringHasNoText (str)) {
    Message (MSG_OK, "Please enter an accession number or gi");
    Show (w);
    Select (w);
    Select (ffp->accession);
    return;
  }
  is_prot = (Boolean) (ISA_aa(ffp->adp->mol_type));
  WatchCursor ();
  sep = NULL;
  uid = 0;
/**/
  EntrezInit("Salsa", TRUE, NULL);
/**/
  if (GetValue (ffp->accntype) == 1) {
    if (is_prot)
       uid = CCAccessionToGi (str, TYP_AA);
    else 
       uid = CCAccessionToGi (str, TYP_NT);
  } else {
    if (! StrToLong (str, &uid)) {
     uid = 0;
    }
  }
  if (uid > 0) {
    ValNodeAddInt (&sip2, SEQID_GI, (Int4) uid);
    bsp=BioseqLockById(sip2);
    ValNodeFree(sip2);
    if (bsp) {
       sep = SeqMgrGetSeqEntryForData ((Pointer)bsp);
    }
    else {
       sep = EntrezSeqEntryGet (uid, -2);
       is_newbsp=(Boolean) (sep!=NULL);
    }
/**/
    EntrezFini (); 
/**/
    if (sep == NULL) {
      Message (MSG_OK, "Unable to find this record in the database.");
    }
    else { 
       Remove (w);
       if (IS_Bioseq(sep))
       {
          if (!is_newbsp) 
          {
             sep = SeqEntryNewForBioseq ((BioseqPtr)sep->data.ptrvalue);
          } 
          if (sep!=NULL) {
             bsp=(BioseqPtr)sep->data.ptrvalue;
             if (bsp) 
             {
                entityID = ObjMgrRegister (OBJ_BIOSEQ, (Pointer)bsp);
                if (entityID>0)
                   align_this(sep, (SeqLocPtr)ffp->adp->master.region, ffp->adp->sap_original, ffp->editor_window, ffp->adp);
                else 
                   SeqEntryFree (sep); 
             }
             else 
                SeqEntryFree (sep); 
          }
       }
       else {
          ErrPostEx (SEV_ERROR, 0, 0, "Unable to process object type.");
          SeqEntryFree (sep); 
       } 
    } 
  }
  else {
    EntrezFini (); 
    Message (MSG_OK, "Unable to find this record in the database.");
  }
  ArrowCursor ();
  Show (w);
  Select (w);
  Select (ffp->accession);
  ArrowCursor ();
  return;
}

static void CCCommonFetchFromNet (BtnActnProc actn, BtnActnProc cancel, EditAlignDataPtr adp, WindoW editor_window)

{
  GrouP              c;
  FetchFormPtr       ffp;
  GrouP              g;
  WindoW             w;

  Update ();
  w = NULL;
  ffp = MemNew (sizeof (FetchForm));
  if (ffp != NULL) {
    w = FixedWindow (-50, -33, -10, -10, "Download From Entrez", NULL);
    SetObjectExtra (w, ffp, StdCleanupFormProc);
    ffp->form = (ForM) w;
    ffp->formmessage = FetchFormMessage;
    ffp->adp= adp;
    ffp->editor_window = editor_window;
    SetGroupSpacing (w, 10, 10);

    g = HiddenGroup (w, -3, 0, NULL);
    StaticPrompt (g, "Type", 0, stdLineHeight, programFont, 'l');
    ffp->accntype = HiddenGroup (g, 4, 0, NULL);
    RadioButton (ffp->accntype, "Accession");
    RadioButton (ffp->accntype, "GI");
    SetValue (ffp->accntype, 1);
    ffp->accession = DialogText (g, "", 6, FetchTextProc);
    SetObjectExtra (ffp->accession, ffp, NULL);

    c = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    ffp->accept = DefaultButton (c, "Retrieve", actn);
    SetObjectExtra (ffp->accept, ffp, NULL);
    Disable (ffp->accept);
    PushButton (c, "Cancel", cancel);

    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
    Show (w);
    Select (w);
    Update ();
  }
}
extern void CCFetchFromNet (EditAlignDataPtr adp, WindoW editor_window)
{
  CCCommonFetchFromNet (CCDownloadProc, StdCancelButtonProc, adp, editor_window);
}

/*******************
*** FEATURES
********************/

static void get_client_rect (PaneL p, RectPtr prc)
{
  ObjectRect (p, prc);
  InsetRect (prc, HRZ_BORDER_WIDTH, VER_BORDER_WIDTH);
}

/***************************************************************
***  switch_featOrder
***
*****************************************************************/
static void switch_featOrder (EditAlignDataPtr adp, Uint1 choice)
{
  Int2  oldstyle;
  Int2  j;
  Int4  groupNum;

  if (choice > 0) {
     oldstyle = GetMuskCurrentSt ();
     SetMuskCurrentSt (GetMuskStyleName (adp->styleNum));
     for(j =0; j<FEATDEF_ANY; j++) 
     {   
        adp->featOrder[j] = (Uint1)GetMuskCParam(j, MSM_FORDER, MSM_NUM); 
        groupNum = (Uint1)GetMuskCParam(j, MSM_FGROUP, MSM_NUM); 
        adp->groupOrder[j] = (Uint1)GetMuskCParam(MSM_GROUPS, (Int2)groupNum, MSM_NUM); 
     } 
     SetMuskCurrentSt (GetMuskStyleName (oldstyle));
  } 
  else
     for(j=0; j<FEATDEF_ANY; ++j) adp->featOrder[j] = choice;
}


/*********************************************
***   sesp_to_pept
***
***
*********************************************/
typedef struct ccid {
  SeqIdPtr    sip;
  SeqEntryPtr sep;
  BioseqPtr   bsp;
} CcId, PNTR CcIdPtr;
 
static void FindSeqEntryForSeqIdCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  CcIdPtr            cip;
 
  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     cip = (CcIdPtr)mydata;
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL && ISA_na (bsp->mol)) {
           sip = SeqIdFindBest(bsp->id, 0);
           if (SeqIdForSameBioseq(cip->sip, sip))
              cip->sep = sep;
              cip->bsp = bsp;
        }
     }   
  }
}
 
static Int2 CC_SeqEntryToGeneticCode (Uint2 entityID, SeqIdPtr sip)
{
  SeqEntryPtr sep_head,
              sep;
  CcId        ci;
  Int2        genCode = 0;

  sep_head  = GetTopSeqEntryForEntityID (entityID);
  ci.sip = SeqIdDup (sip);
  ci.sep = NULL;
  ci.bsp = NULL;
  SeqEntryExplore(sep_head,(Pointer)&ci, FindSeqEntryForSeqIdCallback);
  sep = ci.sep;
  SeqIdFree (ci.sip);
  if (sep!=NULL)
     genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
  return genCode;
}


extern Boolean sesp_to_pept (SelEdStructPtr cds, SeqAlignPtr salp, ValNodePtr sqlocs, Boolean partial)
{
  SelEdStructPtr cdsp;
  SelEdStructPtr cds1;
  SeqLocPtr      slp;
  SeqIntPtr      sit;
  ValNodePtr     pept;
  ByteStorePtr   bsp;
  CharPtr        pep = NULL;
  CharPtr        pepPtr = NULL;
  CharPtr        str = NULL, 
                 strPtr = NULL;
  CharPtr        buffer = NULL, 
                 bufferPtr = NULL;
  Int4           k, 
                 strlens, slplens;
  Int4           sumlens;
  Int2           genCode;
  Int2           cb;
  Uint1          codonbase;
  Uint1          codonstart;
  Uint1          strand;

  if ( cds == NULL || salp == NULL )
     return FALSE;
  if (cds->regiontype == 0 || cds->region == NULL)
     return FALSE;
  cds1 = cds;
  while (cds1->prev != NULL) {
     cds1 = cds1->prev;
  }
  slp = sesp_to_slp (cds1, salp, sqlocs, partial);
  slplens = SeqLocLen (slp);
  if ( slplens < 3 ) {
     SeqLocFree (slp);
     return FALSE;
  }
  if (SeqLocStart(slp) > 0)
     codonstart = 1;
  else 
     codonstart = cds1->codonstart;
  strand = SeqLocStrand (slp);

  genCode = CC_SeqEntryToGeneticCode (cds1->entityID, SeqLocId(slp));
  if (genCode == 0)
     genCode = Seq_code_ncbieaa;

  codonbase = 0;
  sit = (SeqIntPtr) slp->data.ptrvalue;
  if (strand == Seq_strand_minus && codonstart > 1) {
     cb = (Int2)(slplens % (Int4) 3);
     if (cb == 1) {
          sit->to --;
     }
  }
  else if (strand == Seq_strand_minus) {
     cb = (Int2)(slplens % (Int4) 3);
     if (cb == 1 && sit->from >0) {
          sit->from --;
     } else if (cb == 2) {
          sit->from ++;
     }
     if (cb == 0) codonbase = 0;
     else if (cb == 1) codonbase = 1;
     else if (cb == 2) codonbase = 2;
  }  
  slplens = SeqLocLen (slp);
  bsp = cds_to_pept (slp, codonstart, genCode, TRUE);
  str = (CharPtr) BSMerge (bsp, NULL);
  BSFree (bsp);
  pep = MemNew ((size_t) ((slplens + 5) *sizeof(Char)));
  pep = emptystring (pep, (Int4)(slplens + 5));
  pep [slplens + 3] = '\0';
  pepPtr = pep;
  *pepPtr = ' ';
  pepPtr += codonbase +1;
  strlens = 3*StringLen(str);
  if (slplens < strlens) {
     strlens=(Int4)(slplens/(Int4)3);
     str [strlens] ='\0';
  }
  if (strand == Seq_strand_minus)
     reverse_string (str);
  strlens = StringLen(str);
  strPtr = str;
  for (k = 0; k < strlens; k++, pepPtr += 3, strPtr++) {
          *pepPtr = *strPtr; 
  }
  MemFree (str);
/*
  strlens = SeqLocLen (slp) + 5;
  buffer = MemNew ((size_t) (strlens *sizeof(Char)));
  buffer = emptystring (buffer, strlens);
  buffer [strlens -1] = '\0';
  bufferPtr = buffer;
  *bufferPtr = ' ';
  sip = SeqLocId (slp);
  for (cdsp= cds1; cdsp != NULL; cdsp = cdsp->next)
  {
     slp = (SeqLocPtr) cdsp->region; 
     buffer = ReadBufferFromSap (pep, buffer, salp, sip, SeqLocStart(slp), SeqLocStop(slp));
  }
  MemFree (pep);
*/
  buffer = pep;
  SeqLocFree (slp);
  if (cds1->data != NULL) {
     pept = cds1->data;
     cds1->data = NULL;
     pept->data.ptrvalue = MemFree (pept->data.ptrvalue);
     pept = ValNodeFree (pept);
  }
  sumlens = 0;
  for (cdsp= cds1; cdsp != NULL; cdsp = cdsp->next)
  {
     pept = ValNodeNew (NULL);
     pept->choice = 0;
     pept->data.ptrvalue = (Pointer) buffer;
     cdsp->data = pept;
     cdsp->offset = sumlens;
     sumlens += SeqLocLen ((SeqLocPtr) cdsp->region);
     pept = NULL;
  }
  return TRUE;
}


/*******************************************************************
***   TranslateProc, TranslateButton
***
***   CdRgnToProtProc
***
********************************************************************/
extern void CdRgnToProtProc (PaneL pnl,  EditAlignDataPtr adp)
{
  WindoW           temport;
  SelStructPtr     ssp;
  SelEdStructPtr   cds;
  ValNodePtr       feathead = NULL, 
                   vnp = NULL;
  Uint2            itemsubtype;
  Boolean          seq_select = FALSE;

  ssp = ObjMgrGetSelected(); 
  for (; ssp != NULL; ssp = ssp->next)
  {  
     if ( checkssp_for_editor (ssp) && ssp->itemtype == OBJ_VIRT) { 
        feathead = adp->feat;
        itemsubtype = FEATDEF_CDS;
     }
     else if ( checkssp_for_editor (ssp) && ssp->itemtype == OBJ_SEQFEAT) {
        feathead = adp->seqfeat;
        itemsubtype = SEQFEAT_CDREGION;
     }
     else feathead = NULL;
     if (feathead != NULL)
     {
        for (vnp = feathead; vnp != NULL; vnp = vnp->next)
        {
           if (vnp->choice == itemsubtype) {
              cds = (SelEdStructPtr) vnp->data.ptrvalue;
              if (cds->entityID == ssp->entityID && cds->itemID == ssp->itemID) 
              {
                 if (sesp_to_pept(cds, (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, TRUE))
                    seq_select = TRUE;
                 break;
              }
           }
        }
     }
  }
  if (!seq_select) return;
  data_collect_arrange (adp, TRUE);
  SeqEdSetCorrectBarMax (pnl, adp->nlines, adp->voffset);
  temport = SavePort(ParentWindow(pnl));
  Select (pnl);
  inval_panel (pnl, -1, -1);
  RestorePort (temport);
  return; 
}

/*******************************************************************
***
***   UntranslateFunc
***
********************************************************************/
extern void UntranslateFunc (PaneL pnl,  EditAlignDataPtr adp)
{
  WindoW           temport;
  SelStructPtr     ssp;
  SelEdStructPtr   cds;
  ValNodePtr       pept;
  ValNodePtr       feathead = NULL, vnp = NULL;
  Uint2            itemsubtype;
  Boolean          seq_select = FALSE;

  ssp = ObjMgrGetSelected(); 
  for (; ssp != NULL; ssp = ssp->next)
  {  
     if ( checkssp_for_editor (ssp) && ssp->itemtype == OBJ_VIRT) { 
        feathead = adp->feat;
        itemsubtype = FEATDEF_CDS;
     }
     else if ( checkssp_for_editor (ssp) && ssp->itemtype == OBJ_SEQFEAT) {
        feathead = adp->seqfeat;
        itemsubtype = SEQFEAT_CDREGION;
     }
     else feathead = NULL;
     if (feathead != NULL) {
        for (vnp = feathead; vnp != NULL; vnp = vnp->next) {
           if (vnp->choice == itemsubtype) 
         {
              cds = (SelEdStructPtr) vnp->data.ptrvalue;
              if (cds->entityID == ssp->entityID && cds->itemID == ssp->itemID) 
              {
                 if (cds->data != NULL) {
                    pept = cds->data;
                    cds->data = NULL;
                    pept->data.ptrvalue = MemFree (pept->data.ptrvalue);
                    ValNodeFree (pept);
                    for (; cds != NULL; cds = cds->next) {
                       cds->data = NULL;
                    }
                    seq_select = TRUE;
                 }
                 break;
              }
           }
        }
     }
  }
  if (!seq_select) return;
  data_collect_arrange (adp, TRUE);
  SeqEdSetCorrectBarMax (pnl, adp->nlines, adp->voffset);
  temport = SavePort(ParentWindow(pnl));
  Select (pnl);
  inval_panel (pnl, -1, -1);
  RestorePort (temport);
  return; 
}



extern Boolean ShowFeatureFunc (EditAlignDataPtr adp)
{
  AlignNodePtr     anp;
  ValNodePtr       vnp;
  SeqLocPtr        slp;
  Boolean          seq_select = FALSE;

  switch_featOrder (adp, 1);
  for (vnp = adp->anp_list; vnp != NULL; vnp = vnp->next) {
     anp = (AlignNodePtr) vnp->data.ptrvalue;
     if ( anp != NULL ) {
        if ( anp->segs->cnp == NULL ) {
           slp = CollectSeqLocFromAlignNode (anp);
           CollectFeatureForAlign (slp, anp, adp->featOrder, adp->groupOrder);
           adp->seqfeat=CollectFeatureForEditor (slp, adp->seqfeat, anp->seq_entityID, anp->bsp_itemID, adp->featOrder, FALSE);
           seq_select = TRUE;
        } 
     }
  }
  if (!seq_select) {
     switch_featOrder (adp, 0);
     return FALSE;
  }
  OrderFeatProc (adp->anp_list);
  if (adp->seqfeat != NULL)
     checkselectsequinfeature_for_editor (adp->seqfeat);
  return TRUE;
}

/***********************************************************
***
***  HideFeatureProc
***
*** loop on Bioseq to delete the features in those selected only.
***
***********************************************************/
extern Boolean HideFeatureFunc (EditAlignDataPtr adp)
{
  AlignNodePtr     anp;
  SelStructPtr     ssp;
  ValNodePtr       vnp;
  AlignSegPtr      asp, aspnext;
  Boolean          seq_select = FALSE;

  switch_featOrder (adp, 0);
  if (adp->input_format == OBJ_BIOSEQ) 
  {
     if ( checkOMss_for_itemtype (OBJ_BIOSEQ) == 0 ) 
          ssp = &(adp->master);
     else ssp = ObjMgrGetSelected();  
     for (; ssp != NULL; ssp = ssp->next)  {
         if ( checkssp_for_editor (ssp) && ssp->itemtype == OBJ_BIOSEQ )  {
            adp->seqfeat =SeqfeatlistFree_fromID (adp->seqfeat, ssp->entityID);
            anp = (AlignNodePtr) AlignNodeFind (adp->anp_list, ssp->entityID, ssp->itemID, ssp->itemtype);
            if ( anp != NULL ) {
               asp = anp->segs;
               while(asp !=NULL)
               {
                     aspnext = asp->next;
                     if(asp->cnp != NULL)
                           FreeFeatureList(asp->cnp);
                     asp->cnp = NULL;
                     if(asp->mismatch)
                           ValNodeFree(asp->mismatch);
                     asp->mismatch = NULL;
                     asp = aspnext;
               }
               seq_select = TRUE;
            }
         }
     }
  }
  else if (adp->input_format == OBJ_SEQALIGN) 
  {
     for (vnp = adp->anp_list; vnp != NULL; vnp = vnp->next) {
            anp = (AlignNodePtr) vnp->data.ptrvalue;
            if ( anp != NULL ) {
               asp = anp->segs;
               while(asp !=NULL)
               {
                     aspnext = asp->next;
                     if(asp->cnp != NULL)
                           FreeFeatureList(asp->cnp);
                     asp->cnp = NULL;
                     if(asp->mismatch)
                           ValNodeFree(asp->mismatch);
                     asp->mismatch = NULL;
                     asp = aspnext;
               }
               seq_select = TRUE;
            }
     }
     if (seq_select) adp->seqfeat =SeqfeatlistFree (adp->seqfeat);
  }
  if (!seq_select) {
     switch_featOrder (adp, 1);
     return FALSE;
  }
  return TRUE;
}

/***********************************************************
***
***  ResetFeatureProc
***
***********************************************************/
extern Boolean ResetFeatureFunc (EditAlignDataPtr adp)
{
/*
  AlignNodePtr     anp;
  ValNodePtr       vnp;
  SeqLocPtr        slp;
  AlignSegPtr      asp, aspnext;
  Boolean          seq_select = FALSE;

  switch_featOrder (adp, 1);
  ssp = Gettranslation (adp->seqfeat);
  adp->seqfeat = SeqfeatlistFree (adp->seqfeat);
  for (vnp = adp->anp_list; vnp != NULL; vnp = vnp->next) {
     anp = (AlignNodePtr) vnp->data.ptrvalue;
     if ( anp != NULL ) {
        asp = anp->segs;
        while(asp !=NULL)
        {
           aspnext = asp->next;
           if(asp->cnp != NULL)
                           FreeFeatureList(asp->cnp);
           asp->cnp = NULL;
           if(asp->mismatch)
                           ValNodeFree(asp->mismatch);
           asp->mismatch = NULL;
           asp = aspnext;
           slp = CollectSeqLocFromAlignNode (anp);
           CollectFeatureForAlign (slp, anp, adp->featOrder, adp->groupOrder);
           adp->seqfeat=CollectFeatureForEditor (slp, adp->seqfeat, anp->seq_entityID, anp->bsp_itemID, adp->featOrder, FALSE);
        } 
     }
  }
  if (adp->seqfeat !=NULL) {
        for (vnp = adp->seqfeat; vnp != NULL; vnp = vnp->next)
        {
           if (vnp->choice == SEQFEAT_CDREGION) {
              cds = (SelEdStructPtr) vnp->data.ptrvalue;
              if (cds->entityID == ssp->entityID && cds->itemID == ssp->itemID)
              {
                 if (sesp_to_pept(cds, (SeqAlignPtr) adp->sap_align->data, adp->
sqloc_list, TRUE))
              }
           }
        }
  }
*/
  return TRUE;
}

/**********************************************************************
***   PropagateFeatureProc
***      build features taking one selected feature as template
***
*** 
LIST starts with 15 
***********************************************************************/

#define first_GBFeat  15
#define number_GBFeat 58
static CharPtr GBFeat[number_GBFeat] = {
"allele", "attenuator", "C_region", "CAAT_signal", "CDS", 
"conflict", "D-loop", "D_segment", "enhancer",  "exon",  
"GC_signal", "gene",  "intron",  "J_segment",  "LTR",  
"mat_peptide", "misc_binding",   "misc_difference",  
"misc_feature", "misc_recomb",  "misc_RNA",   "misc_signal",   
"misc_structure", "modified_base",  "mutation", "N_region", 
"old_sequence", "polyA_signal",  "polyA_site", "precursor_RNA",   
"prim_transcript", "primer_bind",  "promoter",   "protein_bind",  "RBS",  
"repeat_region", "repeat_unit",  "rep_origin",  "S_region",  "satellite",  
"sig_peptide", "source",  "stem_loop",  "STS",   "TATA_signal", 
"terminator", "transit_peptide",  "unsure",   "V_region",   "V_segment",   
"variation", "virion",   "3'clip",   "3'UTR",   "5'clip",  "5'UTR",  
"-10_signal", "-35_signal"};

static Boolean FindSqFeatItem (GatherContextPtr gcp)
{
  SeqFeatPtr PNTR sfpp;
 
  sfpp = (SeqFeatPtr PNTR) gcp->userdata;
  if (sfpp != NULL && gcp->thistype == OBJ_SEQFEAT) {
    *sfpp = (SeqFeatPtr) gcp->thisitem;
  }
  return TRUE;
}


static void PropagateFeatureProc (ButtoN b)
{
  WindoW           wdialog;
  PaneL            pnl;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  SeqAnnotPtr      sap;
  SeqAlignPtr      salp;
  SelEdStructPtr   feat;
  SelEdStructPtr   sesp;
  SeqIdPtr         featsip;
  SeqLocPtr        featslp;
  SeqLocPtr        new_slp;
  ValNodePtr       vnp,
                   vnpf,
                   vnpfeat = NULL;
  ValNodePtr       vnpsfp = NULL;
  AlignNodePtr     anp;
  SeqFeatPtr       source_sfp;
  SeqFeatPtr       source_dup;
  Uint2            eID, iID,
                   subtype;
  Int2             j, jmax,
                   k, kmax;
  Uint1            frame;
  Boolean          val;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update ();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (b));
  if ( ( pnl = GetPanelFromWindow (dbdp->w) ) != NULL ) 
  {
   if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) { 
    if (adp->seqnumber > 0)  
    {
     sap = SeqAnnotBoolSegToDenseSeg (adp->sap_align);
     salp = (SeqAlignPtr) sap->data;
     vnpfeat = NULL;
     vnpf = adp->seqfeat;
     kmax = CountItems (dbdp->lst3);
     jmax = CountItems (dbdp->lst2);
     for (k=1; k<=kmax; k++)
     {
        val = GetItemStatus (dbdp->lst3, k);
        if (val) {
           vnp = adp->anp_list; 
           for (j=1; j<=jmax; j++) 
           {
              val = GetItemStatus (dbdp->lst2, j);
              if (val) {
                    feat  = (SelEdStructPtr) vnpf->data.ptrvalue;
                    featslp = (SeqLocPtr) feat->region;
                    featsip = SeqLocId (featslp);
                    subtype = vnpf->choice;
                    GatherItem (feat->entityID, feat->itemID, OBJ_SEQFEAT, (Pointer) (&source_sfp), FindSqFeatItem); 
                    if (source_sfp != NULL) {
                       anp = (AlignNodePtr) vnp->data.ptrvalue;
                       new_slp = CopySeqLocFromSeqAlign (source_sfp, anp->sip, featsip, salp, adp->gap_choice, &frame);
                       if (new_slp != NULL) {
                          eID = anp->seq_entityID;
                          iID = anp->bsp_itemID;
                          if (is_newfeat (adp->seqfeat, eID, new_slp) )
                          {
                             source_dup = (SeqFeatPtr) AsnIoMemCopy((Pointer) source_sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
                             sesp = new_seledstruct_fromseqloc (eID, iID, subtype, iID, new_slp, feat->label, NULL, 0, frame);
                             if (sesp != NULL) {
                              ValNodeAddPointer(&vnpfeat, 0, (Pointer) sesp);
                              ValNodeAddPointer(&vnpsfp,0, (Pointer)source_dup);
                             }
                          }
                       }
                    }
              }
              vnp = vnp->next;
           }
        }
        vnpf = vnpf->next;
     }
     val = ApplyNewSeqFeat (vnpfeat, vnpsfp, adp->stoptransl);
     if (val) {    
        sap = SeqAnnotFree (sap);
        ObjMgrSendMsg (OM_MSG_UPDATE, adp->master.entityID, adp->master.itemID, OBJ_BIOSEQ); 
     }
    }
   }
  }
  Remove (wdialog);
  return; 
}


static void getchoicegaps (GrouP c)
{
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Int2             j;

  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (c));
  if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) {
     j = GetValue (c);
     if (j == 2) 
        adp->gap_choice = IGNORE_GAP_CHOICE;
     else
        adp->gap_choice = DEFAULT_GAP_CHOICE;
  }
}

static void select_lst_sseq (LisT lst)
{
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  ValNodePtr       vnp,
                   vnpf;
  SelEdStructPtr   feat;
  SeqIdPtr         sip;
  SeqLocPtr        slp,
                   slpseq;
  Int2             j, k, jmax, kmax;
  Boolean          val;

  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (lst));
  if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) 
  { 
     vnp = adp->sqloc_list; 
     jmax = CountItems (lst);
     for (j=1; j<=jmax; j++) 
     {
        val = GetItemStatus (lst, j);
        vnpf = adp->seqfeat;
        kmax = CountItems (dbdp->lst3);
        for (k=1; k<=kmax; k++)
        {
           slpseq = (SeqLocPtr) vnp->data.ptrvalue;
           feat  = (SelEdStructPtr) vnpf->data.ptrvalue;
           slp = (SeqLocPtr) feat->region;
           sip = SeqLocId (slp);
           if (SeqIdForSameBioseq (SeqLocId(slpseq), sip)) {
                 SetItemStatus (dbdp->lst3, k, val);
           }
           else { 
              if (val) {
                 SetItemStatus (dbdp->lst3, k, FALSE);
              }
           }
           vnpf = vnpf->next;
        }
        vnp = vnp->next;
     }
  }
}

static void select_lst_feat (LisT lst)
{
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  ValNodePtr       vnp,
                   vnpf;
  SelEdStructPtr   feat;
  SeqLocPtr        slp,
                   slpseq;
  SeqIdPtr         sip;
  Int2             j, k, jmax, kmax;
  Boolean          val;

  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (lst));
  if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) { 
     vnp = adp->sqloc_list; 
     jmax = CountItems (dbdp->lst1);
     for (j=1; j<=jmax; j++) 
     {
        vnpf = adp->seqfeat;
        kmax = CountItems (lst);
        val = FALSE;
        for (k=1; k<=kmax; k++)
        {
           slpseq = (SeqLocPtr) vnp->data.ptrvalue;
           feat  = (SelEdStructPtr) vnpf->data.ptrvalue;
           slp = (SeqLocPtr) feat->region;
           sip = SeqLocId (slp);
           if (SeqIdForSameBioseq (SeqLocId(slpseq), sip)) {
              val = GetItemStatus (lst, k);
           }
           if (val) break;
           vnpf = vnpf->next;
        }
        SetItemStatus (dbdp->lst1, j, val);
        vnp = vnp->next;
     }
  }
}

static void ExtTranslButton (ButtoN bn)
{
  WindoW           wdialog;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  adp->stoptransl = !(GetStatus (bn));
  return;
}

static void selectall (LisT lst)
{
  Int2 j, max;
  max = CountItems (lst);
  for (j=1; j<=max; j++) {
     SetItemStatus (lst, j, TRUE);
  }
}

static void selectall1 (ButtoN b)
{
  DialogBoxDataPtr dbdp;
  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (b));
  selectall (dbdp->lst1);
}

static void selectall2 (ButtoN b)
{
  DialogBoxDataPtr dbdp;
  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (b));
  selectall (dbdp->lst2);
}

static void selectall3 (ButtoN b)
{
  DialogBoxDataPtr dbdp;
  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (b));
  selectall (dbdp->lst3);
}

static ValNodePtr ShowAllFeatureFunc (EditAlignDataPtr adp)
{
  AlignNodePtr     anp;
  ValNodePtr       vnp,
                   allseqfeat = NULL;
  SeqLocPtr        slp;

  for (vnp = adp->anp_list; vnp != NULL; vnp = vnp->next) {
     anp = (AlignNodePtr) vnp->data.ptrvalue;
     if ( anp != NULL ) {
           slp = CollectSeqLocFromAlignNode (anp);
           CollectFeatureForAlign (slp, anp, adp->featOrder, adp->groupOrder);
           allseqfeat = CollectFeatureForEditor (slp, allseqfeat, anp->seq_entityID, anp->bsp_itemID, adp->featOrder, TRUE);
     }
  }
  return allseqfeat;
}

static ValNodePtr AddAllFeatureFunc (EditAlignDataPtr adp, ValNodePtr allseqfeat)
{
  AlignNodePtr     anp;
  ValNodePtr       vnp;
  SeqLocPtr        slp;

  for (vnp = adp->anp_list; vnp != NULL; vnp = vnp->next) {
     anp = (AlignNodePtr) vnp->data.ptrvalue;
     if ( anp != NULL ) {
           slp = CollectSeqLocFromAlignNode (anp);
           CollectFeatureForAlign (slp, anp, adp->featOrder, adp->groupOrder);
           allseqfeat = CollectFeatureForEditor (slp, allseqfeat, anp->seq_entityID, anp->bsp_itemID, adp->featOrder, TRUE);
     }
  }
  return allseqfeat;
}

extern void PropagateFeatDialog (IteM i)
{
  WindoW           w, wdialog;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  GrouP            g1, g2, g2b, g3, g4, g5,
                   d, c, h;
  LisT             lst_seq,
                   lst_sseq,
                   lst_feat; 
  Char             str [128];
  Char             str2 [24];
  CharPtr          tmp;
  CharPtr          strp;
  ValNodePtr       vnp;
  SelEdStructPtr   feat;
  SeqLocPtr        slp,
                   slpseq;
  Int4             start, stop;
  Uint1            strand;
  ButtoN           b;
  PrompT           p;

  ValNodePtr  allseqfeat = NULL;

  w = ParentWindow (i);
  if ( ( adp = GetAlignEditData (w) ) == NULL ) 
     return;
  adp->seqfeat = AddAllFeatureFunc (adp, adp->seqfeat);

  allseqfeat = adp->seqfeat;
  wdialog = FixedWindow (-50, -33, -10, -10, "Feature Propagation", StdCloseWindowProc);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;

  g1 = HiddenGroup (wdialog, 0, -5, NULL);

  g2 = HiddenGroup (g1, 2, 0, NULL);

  g4 = HiddenGroup (g2, 0, -3, NULL);
  p = StaticPrompt (g4, "Select source sequences", 0, dialogTextHeight,systemFont, 'c');
  lst_sseq = ExtendedList (g4, 12, 5, select_lst_sseq);
  for (vnp = adp->sqloc_list; vnp!=NULL; vnp = vnp->next)
  {
     str [0] = '\0';
     slpseq = (SeqLocPtr) vnp->data.ptrvalue;
     SeqIdWrite (SeqLocId(slpseq), str, PRINTID_REPORT, sizeof (str));
     ListItem (lst_sseq, str);
  }
  dbdp->lst1 = lst_sseq;
  b = PushButton (g4, "Select all", selectall1);
  AlignObjects (ALIGN_CENTER, (HANDLE) lst_sseq, (HANDLE) b, (HANDLE) p, NULL);

  g5 = HiddenGroup (g2, 0, -3, NULL);
  p = StaticPrompt (g5, "Select target sequences", 0, dialogTextHeight,systemFont, 'c');
  lst_seq = ExtendedList (g5, 12, 5, NULL);
  for (vnp = adp->sqloc_list; vnp!=NULL; vnp = vnp->next)
  {
     str [0] = '\0';
     slpseq = (SeqLocPtr) vnp->data.ptrvalue;
     SeqIdWrite (SeqLocId(slpseq), str, PRINTID_REPORT, sizeof (str));
     ListItem (lst_seq, str);
  }
  dbdp->lst2 = lst_seq;
  b = PushButton (g5, "Select all", selectall2);
  AlignObjects (ALIGN_CENTER, (HANDLE) lst_seq, (HANDLE) b, (HANDLE) p, NULL);

  g2b = HiddenGroup (g1, 1, 0, NULL);
  g3 = HiddenGroup (g2b, -1, 0, NULL);
  p = StaticPrompt (g3, "Select source Features", 0, dialogTextHeight,systemFont, 'c');
  lst_feat = ExtendedList (g3, 31, 5, select_lst_feat);
  for (vnp = allseqfeat; vnp!=NULL; vnp=vnp->next)
  {
     str [0] = '\0';
     tmp = str;
     feat  = (SelEdStructPtr) vnp->data.ptrvalue;
     slp = (SeqLocPtr) feat->region;
     start = SeqLocStart (slp) +1;
     stop = SeqLocStop (slp) +1;
     strand = SeqLocStrand (slp);
     if (vnp->choice == FEATDEF_GENE)
        tmp = StringMove (tmp, "GENE: ");
     else  if (vnp->choice == FEATDEF_mRNA)
        tmp = StringMove (tmp, "mRNA: ");
     else  if (vnp->choice == FEATDEF_CDS)
        tmp = StringMove (tmp, "CDS: ");
     else if (vnp->choice>=first_GBFeat && vnp->choice<number_GBFeat) {
        strp = GBFeat[(vnp->choice-first_GBFeat)];
        tmp = StringMove (tmp, strp); 
        tmp = StringMove (tmp, ": ");
     }
     if (feat->label != NULL)
        tmp = StringMove (tmp, feat->label);
     if (strand == Seq_strand_minus) {
        sprintf (str2, " (%ld..%ld) minus strand", (long)start, (long)stop);
     } else
        sprintf (str2, " (%ld..%ld)", (long)start, (long)stop);
     tmp = StringMove (tmp, str2);
     ListItem (lst_feat, str);
  }
  dbdp->lst3 = lst_feat;
  b = PushButton (g3, "Select all", selectall3);
  AlignObjects (ALIGN_CENTER, (HANDLE) lst_feat, (HANDLE) b, (HANDLE) p, NULL);

  c = HiddenGroup (g1, 2, 0, getchoicegaps);
  RadioButton (c, "split at gaps");
  RadioButton (c, "extend over gaps");
  SetValue (c, (Int2)(adp->gap_choice + 1));

  d = HiddenGroup (g1, 1, 0, NULL);
  CheckBox (d, "extend translation after internal stop codon", ExtTranslButton);

  h = HiddenGroup (g1, 2, 0, NULL);
  PushButton (h, "Propagate", PropagateFeatureProc);
  PushButton (h, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g2, (HANDLE) g3, (HANDLE) c, (HANDLE) d, (HANDLE) h, NULL);

  RealizeWindow (wdialog);
  Show (wdialog);
  return;
}

/******************************************************************/
static SelEdStructPtr split_feat (SelEdStructPtr feat, Int4 pos, Int4 changevalue)
{
  SeqLocPtr      slpfeat;
  SeqIntPtr      sitfeat;
  SelEdStructPtr new, next;
  Int4           from , to;

  slpfeat = (SeqLocPtr) feat->region;
  next = feat->next;
  if (changevalue >= 0)
     from = (Int4)(pos + changevalue);
  else {
     from = pos;
  }
  to = (Int4)(SeqLocStop(slpfeat) + changevalue);
  sitfeat = (SeqIntPtr) slpfeat->data.ptrvalue;
  sitfeat->to = pos -1;
  new = new_seledstruct (feat->entityID, feat->itemID, feat->itemtype, feat->bsp_itemID, from, to, SeqLocId (slpfeat), SeqLocStrand (slpfeat), FALSE, feat->label, feat->data, feat->offset + SeqLocLen(slpfeat), 1);
  feat->next = new;
  new->next = next;
  new->prev = feat;
  return new;
}

extern ValNodePtr update_featpept (EditAlignDataPtr adp, ValNodePtr feathead, RecT *rp, SelStructPtr ssp, Int4 changevalue, Uint2 itemsubtype)
{
  ValNodePtr     vnpfeat,
                 vnpfeatnext;
  SelEdStructPtr feat,
                 feat1, next;
  SeqLocPtr      slpfeat;
  SeqLocPtr      slpssp;
  SeqIntPtr      sitfeat;
  Int4           width;
  Int4           lg;
  Boolean        overlap, precede, succeed, deletefeat;

  if (ssp == NULL) return feathead;
  if (ssp->regiontype == 0 || ssp->region == NULL) return feathead;
  slpssp = (SeqLocPtr) ssp->region;
  if (SeqLocStart(slpssp) == SeqLocStop(slpssp)) return feathead;  
  width = adp->visibleWidth;
  if (adp->columnpcell > 0) 
         width += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;
  vnpfeat = feathead;
  while (vnpfeat != NULL)
  {
    vnpfeatnext = vnpfeat->next;   
    feat1 = (SelEdStructPtr) vnpfeat->data.ptrvalue;
    if( feat1 != NULL && (vnpfeat->choice ==itemsubtype || itemsubtype == 255)) 
    {
      feat = feat1;
      while (feat != NULL) 
      {
         next = feat->next;
         slpfeat = (SeqLocPtr) feat->region;
         overlap = overlapp_ssp (slpssp, slpfeat); 
         precede = precede_ssp (slpssp, slpfeat);
         succeed = succeed_ssp (slpssp, slpfeat);
         deletefeat = FALSE;
         if (overlap || precede || succeed ) 
         {
            sitfeat = (SeqIntPtr) slpfeat->data.ptrvalue;
            if (precede) 
            {
               sitfeat->from = sitfeat->from + changevalue ;
               sitfeat->to = sitfeat->to + changevalue ;
            }
            else if (succeed) 
            {
               if (changevalue < 0 ) sitfeat->to =sitfeat->to +changevalue ;
            }
            else if (overlap) 
            {
               if (changevalue < 0)
               {
                  if ( include_ssp (slpssp, slpfeat) ) {
                     deletefeat = TRUE;
                     feathead = del_ssp_fromid (feathead, itemsubtype, feat);
                     if (rp != NULL)
                        inval_rect (rp->left, rp->top, rp->right, rp->bottom);
                  }
                  else if ( include_ssp (slpfeat, slpssp) ) {
                     if (!adp->spliteditmode)
                          sitfeat->to = sitfeat->to + changevalue;
                     else
                          feat=split_feat(feat,SeqLocStart(slpssp), changevalue);
                  }
                  else if ((lg = overlapp_startssp (slpssp, slpfeat)) > 0) {
                     if (changevalue < 0) {
                          sitfeat->from = sitfeat->from - (abs(changevalue)-lg);
                          sitfeat->to = sitfeat->to + changevalue;
                     }
                     else 
                        ErrPostEx (SEV_ERROR, 0, 0, "Cut what ?");
                  }
                  else if ((lg = overlapp_startssp (slpfeat, slpssp)) > 0) {
                     if (changevalue < 0) {
                          sitfeat->to = sitfeat->to - lg ;
                     } 
                     else sitfeat->to = sitfeat->to + lg ;
                  }
               }
               else {
                  if (!adp->spliteditmode)
                     sitfeat->to = sitfeat->to + changevalue ;
                  else {
                     feat = split_feat (feat, SeqLocStart(slpssp), changevalue);
                  }
               }
            }
            if (!deletefeat && rp != NULL)
            {
               inval_selstruct(adp, feat->entityID, feat->itemID, feat->itemtype, itemsubtype, rp, adp->margin.left,(Int2)(width *adp->charw));
               inval_selstruct(adp, feat->entityID, feat->itemID, feat->itemtype, itemsubtype, rp, adp->margin.left, (Int2)(width *adp->charw));
            }
         }
         feat = next;
      }
      if (feat1 != NULL)
         if (feat1->data != NULL)
            sesp_to_pept (feat1, (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, TRUE);
    }
    vnpfeat = vnpfeatnext;
  }
  return feathead;
}

extern void ShowFeatureProc (PaneL pnl, Boolean invalidate) 
{
  WindoW             temport;
  SeqEditViewFormPtr wdp;
  EditAlignDataPtr   adp;
  Boolean            ok;

  wdp = (SeqEditViewFormPtr)GetObjectExtra (ParentWindow(pnl));
  if ( wdp == NULL ) return;
  if ( ( adp = GetAlignDataPanel (pnl) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  adp->showfeat = (!adp->showfeat);
  ResetClip ();
  WatchCursor ();
  if ( adp->showfeat ) {
        ok = (Boolean) ShowFeatureFunc (adp);
        if (ok) {
           Enable (wdp->hidefeatitem);
           Disable(wdp->showfeatitem);
           SetTitle (wdp->showfeatbt, "Hide Feat.");
        }
  }
  else  {
        ok = (Boolean) HideFeatureFunc (adp);
        if (ok) {
           Disable (wdp->hidefeatitem);
           Enable(wdp->showfeatitem);
           SetTitle (wdp->showfeatbt, "Show Feat.");
        }
  }
  if (ok) {
     data_collect_arrange (adp, TRUE);
     SeqEdSetCorrectBarMax (pnl, adp->nlines, adp->voffset);
     if (invalidate) {
        temport = SavePort(ParentWindow(pnl));
        Select (pnl);
        inval_panel (pnl, -1, -1);
        RestorePort (temport);
     }
  }
  ArrowCursor ();
}

/******************************************************
***
***   LaunchCDSEditor on a Bioseq (input_itemID)
***
*******************************************************/
static void LaunchCDSEditor (Uint2 input_entityID, Uint2 input_itemID, SeqLocPtr slp, Uint1 codonstart)
{
  WindoW         w;
  SeqEntryPtr    top_sep;
  FeatureFormPtr cfp;

  if (slp != NULL && input_entityID != 0) 
  {
     top_sep = GetTopSeqEntryForEntityID (input_entityID);
     input_entityID = SeqMgrGetEntityIDForSeqEntry (top_sep);
   
     w = (WindoW) CreateCdRgnForm (-50, -33, "Coding Region", NULL, top_sep, CdRgnFeatFormActnProc);
     cfp = (FeatureFormPtr) GetObjectExtra (w);
     if (cfp != NULL) {
        cfp->input_entityID = input_entityID;
        cfp->input_itemID = input_itemID;
        cfp->input_itemtype = OBJ_BIOSEQ;
        cfp->this_itemtype = OBJ_SEQFEAT;
        cfp->this_subtype = FEATDEF_CDS;
        PointerToForm (cfp->form, NULL);
        SendMessageToForm (cfp->form, VIB_MSG_INIT);
        PointerToDialog (cfp->location, (Pointer) slp);
        CdRgnTranslateWithFrame (cfp->form, 1);
     }
     Show (w);
     Select (w);
  }
  return;
}
/***************************************************************
***  slpfeatreplacefunc
***
*****************************************************************/
static Boolean slpfeatreplacefunc(GatherContextPtr gcp)
{
  SeqFeatPtr sfp;
  SeqLocPtr  slp;
  Boolean    p3, p5; 

  if(gcp->thistype != OBJ_SEQFEAT) 
     return FALSE;
  sfp = (SeqFeatPtr)(gcp->thisitem);
  slp = (SeqLocPtr) gcp->userdata;
  CheckSeqLocForPartial (sfp->location, &p5, &p3);
  SetSeqLocPartial (slp, p5, p3);
  sfp->location = SeqLocFree (sfp->location);  
  sfp->location = slp;
  return TRUE;
}

/******************************************************************
***
*** SaveFeatProc 
***     look at the selected items
***     if new feature, attaches it (AttachDataForProc)
***     other, replaces it (GatherItem)
***     sends a message to ObjMgr (ObjMgrSendMsg (OM_MSG_UPDATE..))
***     write the seqenrty in the temporary file
***
*** SaveFeatureProc, SaveFeatureButton : call SaveFeatProc
***     sends a message to ObjMgr (ObjMgrSendMsg (OM_MSG_UPDATE..))
***
*******************************************************************/
extern void SaveFeatProc (PaneL pnl)
{
  EditAlignDataPtr   adp;
  SeqEntryPtr        sep;
  SelStructPtr       ssp = NULL;
  SelEdStructPtr     feat;
  ValNodePtr         vnp,
                     next;
  SeqLocPtr          slp;
  RecT               rp;
  Int2               width;
  Uint2              bsp_eID, bsp_iID;
  Int2               handled;

  if ( ( adp = GetAlignDataPanel (pnl) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  if ( checkOMss_for_itemtype (OBJ_VIRT) == 0 
    && checkOMss_for_itemtype (OBJ_SEQFEAT) == 0 ) return;
  sep = GetBestTopParentForItemID (adp->master.entityID, adp->master.itemID, adp->master.itemtype);
  if (sep == NULL)
     return;
  ssp = ObjMgrGetSelected(); 
  for (; ssp != NULL; ssp = ssp->next)
  {  
     if ( checkssp_for_editor (ssp) && ssp->itemtype == OBJ_VIRT) 
     {
        vnp = adp->feat; 
        while  (vnp != NULL) {
           next = vnp->next;
           if (vnp->choice == SEQFEAT_CDREGION) {
              feat = (SelEdStructPtr) vnp->data.ptrvalue;
              if ( is_samess_ses (ssp, feat) ) 
              {
                 if (SeqLocStart((SeqLocPtr)feat->region)==0 || feat->codonstart == 1) 
                 {
                    adp->curfeat = feat;
                    slp = sesp_to_slp (feat, (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE);
                    bsp_eID = SeqMgrGetEntityIDForSeqEntry (sep);
                    bsp_iID = feat->bsp_itemID; 
                    LaunchCDSEditor (bsp_eID, bsp_iID, slp, feat->codonstart);
                 }
                 else 
                    ErrPostEx (SEV_ERROR, 0, 0, "Codon start must be 1");
                 break;
              }
           }
        }
        vnp = next;
     }
     else if ( checkssp_for_editor (ssp) && ssp->itemtype == OBJ_SEQFEAT) 
     {
        vnp = adp->seqfeat; 
        while  (vnp != NULL) 
        {
           next = vnp->next;
           if (vnp->choice == SEQFEAT_CDREGION
            || vnp->choice == SEQFEAT_GENE || vnp->choice == SEQFEAT_RNA) 
           {
              feat = (SelEdStructPtr) vnp->data.ptrvalue;
              if ( is_samess_ses (ssp, feat) ) 
              {
                 slp = sesp_to_slp (feat,(SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE);
                 if (slp != NULL) 
                 {
                    bsp_eID = feat->entityID;
                    bsp_iID = feat->itemID;
                    GatherItem (bsp_eID, bsp_iID, OBJ_SEQFEAT, (Pointer)(slp), slpfeatreplacefunc);
                    ObjMgrSendMsg (OM_MSG_UPDATE, adp->master.entityID,  adp->master.itemID, adp->master.itemtype);
                    HideFeatureFunc (adp);
                    adp->showfeat = FALSE;
                    ShowFeatureProc(pnl, FALSE);
                    get_client_rect (pnl, &rp);
                    width = adp->visibleWidth;
                    if (adp->columnpcell > 0)
                      width +=(Int2)adp->visibleWidth/(Int2) adp->columnpcell;
                    inval_all (adp, &rp, (Uint2)255, OBJ_VIRT, OBJ_SEQFEAT, width);
 
                    WatchCursor ();
                    Update ();
                    handled = GatherProcLaunch (OMPROC_EDIT, FALSE, bsp_eID, bsp_iID, OBJ_SEQFEAT, 0, 0, OBJ_SEQFEAT, 0);
                    ArrowCursor ();
                    Update ();
                 }
                 break;
              }
           }
           vnp = next;
        }
        seqentry_write (sep, adp->tmpfile);
     }
  }
  adp->dirty = TRUE;
  if (!adp->showfeat) {
     ShowFeatureProc(pnl, TRUE);
  }
  return;
}

/******************************************************************
***
*** SaveAllFeatProc 
***     looks at the features items
***     attaches the new features (AttachDataForProc)
***     replaces the old features (GatherItem)
***
OLD VERSION:

  OMProcControl      ompc;
  EditAlignDataPtr   adp;
  SelEdStructPtr     feat;
  ValNodePtr         vnp;
  SeqLocPtr          slp;

  if ( ( adp = GetAlignDataPanel (pnl) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;

  MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
  ompc.input_entityID = adp->master.entityID;
  ompc.input_itemID = adp->master.itemID;
  ompc.input_itemtype = adp->master.itemtype;
  ompc.output_itemtype = OBJ_SEQFEAT;

  for  (vnp=adp->seqfeat; vnp != NULL; vnp = vnp->next) 
  {
     if (vnp->choice == SEQFEAT_CDREGION
      || vnp->choice == SEQFEAT_GENE || vnp->choice == SEQFEAT_RNA) {
        feat = (SelEdStructPtr) vnp->data.ptrvalue;
        slp = sesp_to_slp (feat, (SeqAlignPtr) adp->sap_align->data, FALSE);
        if (slp != NULL) {
           GatherItem (feat->entityID, feat->itemID, OBJ_SEQFEAT,
                      (Pointer)(slp), slpfeatreplacefunc);
       }
     }
  }
*******************************************************************/
extern void SaveAllFeatProc (PaneL pnl)
{
  return;
}


static void MakeFeatFunc (EditAlignDataPtr adp, SelStructPtr ssp, Uint2 itemsubtype, Uint1 strand)
{
  SelEdStructPtr   feat = NULL;
  SelStructPtr     ssptmp;
  SeqLocPtr        slp;
  SeqIntPtr        sit;
  SeqLocPtr        slpfeat;
  Int4             from, to;
  Uint2            itemID;

  slp = (SeqLocPtr) ssp->region;
  ssptmp = is_selectedbyID (ssp->entityID, 255, OBJ_VIRT);
  if (ssptmp == NULL) {
         adp->nfeat++;
         itemID = adp->nfeat;
  } else {
         itemID = ssptmp->itemID;
  }
  from = SeqLocStart (slp);
  to = SeqLocStop (slp);
  if (to == APPEND_RESIDUE) {
     slpfeat = ValNodeNew (NULL);
     slpfeat->choice = SEQLOC_WHOLE;
     slpfeat->data.ptrvalue = (Pointer) SeqIdDup (SeqLocId(slp));
     to = SeqLocLen(slpfeat) -1;
     SeqLocFree (slpfeat);
  }
  if (from >= 0 && to > 0) {
  sit = (SeqIntPtr) slp->data.ptrvalue;
  sit->strand = strand;
  feat = ss_to_ses (ssp); 
  slpfeat = (SeqLocPtr)feat->region;
  setposition_toses (feat, from, to);

  /*  feat->entityID !!!!!!!!!!!!!!!!!!!!!*/
  feat->bsp_itemID = feat->itemID;
  feat->itemID = itemID;
  feat->itemtype = OBJ_VIRT;
  feat->codonstart = 1;
  feat->offset = 0;
  feat->dirty = TRUE;
  feat->next = NULL;
  feat->prev = NULL;
  adp->feat = AddFeatFunc (feat, &(adp->feat), itemsubtype);
  }
}

extern void MakeFeatProc (PaneL pnl, Uint2 itemsubtype, Uint1 strand)
{
  WindoW           temport;
  EditAlignDataPtr adp;
  SelStructPtr     ssp;

  if ( ( adp = GetAlignDataPanel (pnl) ) != NULL ) {
     if (adp->seqnumber > 0 || ISA_na(adp->mol_type)) {
        ssp = ObjMgrGetSelected();  
        for (; ssp != NULL; ssp = ssp->next) 
        {
           if (checkssp_for_editor (ssp) && ssp->itemtype == OBJ_BIOSEQ ) {
              MakeFeatFunc (adp, ssp, itemsubtype, strand);
           }
        }
        data_collect_arrange (adp, TRUE);
        SeqEdSetCorrectBarMax (pnl, adp->nlines, adp->voffset);
        adp->dirty = TRUE;
        temport = SavePort(ParentWindow(pnl));
        Select (pnl);
        inval_panel (pnl, -1, -1);
        RestorePort (temport);
     }
  }
  return; 
}

extern void TranslateAllBioseq (PaneL pnl,  EditAlignDataPtr adp)
{
  WindoW           temport;
  SelEdStructPtr   cds;
  AlignNodePtr     anp;
  SelStructPtr     ssp;
  SeqLocPtr        slp;
  ValNodePtr       vnp = NULL;


  if (adp->seqnumber > 0 || ISA_na(adp->mol_type)) {
     for (vnp = adp->anp_list; vnp != NULL; vnp = vnp->next) {
        if ( (anp = (AlignNodePtr) vnp->data.ptrvalue) != NULL)
        {
           slp = CollectSeqLocFromAlignNode(anp);
           if (slp!=NULL) {
              ssp = SelStructNew (anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ, SeqLocStart(slp), SeqLocStop(slp), SeqLocId(slp), SeqLocStrand(slp), FALSE);
              if (ssp!=NULL)
                 MakeFeatFunc (adp, ssp, SEQFEAT_CDREGION, Seq_strand_plus);
           }
        }
     }
     vnp = adp->feat;
     if (vnp != NULL)
     {   
        for (; vnp != NULL; vnp = vnp->next)
        {
           if (vnp->choice == FEATDEF_CDS) {
              cds = (SelEdStructPtr) vnp->data.ptrvalue;
              sesp_to_pept(cds, (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, TRUE); 
           }   
        }
        data_collect_arrange (adp, TRUE);
        SeqEdSetCorrectBarMax (pnl, adp->nlines, adp->voffset);
        temport = SavePort(ParentWindow(pnl));
        Select (pnl);
        inval_panel (pnl, -1, -1);
        RestorePort (temport);
     }   
  }
  return;
}

