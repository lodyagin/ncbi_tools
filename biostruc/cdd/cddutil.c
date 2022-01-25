/* $Id: cddutil.c,v 1.13 2000/10/10 18:55:06 shoemake Exp $
*===========================================================================
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
* File Name:  cddutil.c
*
* Author:  Aron Marchler-Bauer
*
* Initial Version Creation Date: 10/18/1999
*
* $Revision: 1.13 $
*
* File Description: CDD utility routines
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddutil.c,v $
* Revision 1.13  2000/10/10 18:55:06  shoemake
* fixed memory error in CddAssignProfileRange
*
* Revision 1.12  2000/09/08 21:43:51  hurwitz
* adding PSSM calculation to DDE
*
* Revision 1.11  2000/09/01 21:59:04  hurwitz
* re-order columns from PSSM of CDs to column order expected in threading
*
* Revision 1.10  2000/08/30 21:33:55  hurwitz
* added new and free functions for Seq_Mtf and Qry_Seq
*
* Revision 1.9  2000/08/30 16:03:57  bauer
* fixed GCC compiler warning fo CddGetPairId
*
* Revision 1.8  2000/08/14 21:52:04  hurwitz
* added CddCalcPSSM
*
* Revision 1.7  2000/08/14 19:36:26  hurwitz
* got CddCposComp working and tested
*
* Revision 1.6  2000/08/11 19:54:00  hurwitz
* restored CddDenDiagCposComputation and CddCposComputation to original, added CddCposComp which combines the 2
*
* Revision 1.5  2000/08/10 18:18:59  kans
* commented out direct calls to ID1 services
*
* Revision 1.4  2000/08/10 16:50:06  kans
* use StringSave instead of unavailable strdup
*
* Revision 1.3  2000/08/09 21:29:08  hurwitz
* adding cddutil.c to VC++ build
*
* Revision 1.2  2000/07/19 19:32:55  bauer
* added modification logging
*
*
* ==========================================================================
*/


#include <stdio.h>
#include <ncbi.h>
#include <accentr.h>
#include <lsqfetch.h>
/* #include <netentr.h> */
/* #include <www.h> */
/* #include <sys/resource.h> */
#include <accutils.h>
#include <mmdbapi.h>
#include <mmdbapi1.h>
/* #include <asnmime.h> */
#include <objseq.h>
#include <objmime.h>
#include <strimprt.h>
#include <cdd.h>
/* #include <bjcdd.h> */
#include <cddutil.h>
#include <edutil.h>
#include <posit.h>
#include <cddposutil.h>
#include <blastpri.h>
#include <accid1.h>
#include <thrddecl.h>

static void CddCposCompPart1(SeqAlignPtr listOfSeqAligns, BioseqPtr bsp,
                             compactSearchItems* compactSearch, ValNodePtr* LetterHead,
                             posSearchItems* posSearch);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Cdd asn1 reader and writer wrappers                                       */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Boolean LIBCALL CddWriteToFile(CddPtr pcdd, CharPtr cFile, Boolean bBin)
{
  AsnIoPtr aip = NULL;
  Boolean  bWriteOK;
  
  if (!pcdd) return(FALSE);
  if (bBin) {
     aip = AsnIoOpen(cFile,"wb");
  } else {
     aip = AsnIoOpen(cFile,"w");
  }
  bWriteOK = CddAsnWrite(pcdd,aip,NULL);
  AsnIoClose(aip);
  return(bWriteOK);
}


CddPtr LIBCALL CddReadFromFile(CharPtr cFile, Boolean bBin)
{
  AsnIoPtr aip = NULL;
  CddPtr   pcdd;
  
  if (bBin) {
     aip = AsnIoOpen(cFile,"rb");
     pcdd = CddAsnRead(aip,NULL);
  } else {
     aip = AsnIoOpen(cFile,"r");
     pcdd = CddAsnRead(aip,NULL);
  }

  AsnIoClose(aip);
  return(pcdd);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Cdd-tree asn1 reader and writer wrappers                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Boolean LIBCALL CddTreeWriteToFile(CddTreePtr pcddt, CharPtr cFile, Boolean bBin)
{
  AsnIoPtr aip = NULL;
  Boolean  bWriteOK;
  
  if (!pcddt) return(FALSE);
  if (bBin) {
     aip = AsnIoOpen(cFile,"wb");
  } else {
     aip = AsnIoOpen(cFile,"w");
  }
  bWriteOK = CddTreeAsnWrite(pcddt,aip,NULL);
  AsnIoClose(aip);
  return(bWriteOK);
}


CddTreePtr LIBCALL CddTreeReadFromFile(CharPtr cFile, Boolean bBin)
{
  AsnIoPtr     aip = NULL;
  CddTreePtr   pcddt;
  
  if (bBin) {
     aip = AsnIoOpen(cFile,"rb");
     pcddt = (CddTreePtr) CddTreeAsnRead(aip,NULL);
  } else {
     aip = AsnIoOpen(cFile,"r");
     pcddt = (CddTreePtr) CddTreeAsnRead(aip,NULL);
  }

  AsnIoClose(aip);
  return(pcddt);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Functions to add to a Cdd data structure                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LIBCALL CddAssignDescr(CddPtr pcdd, Pointer pThis, Int4 iWhat, Int4 iIval)
{
  ValNodePtr  vnp;

  vnp = ValNodeNew(NULL);
  vnp->choice = iWhat;
  switch (iWhat) {
    case CddDescr_othername:
      vnp->data.ptrvalue = (CharPtr) pThis;
      break;
    case CddDescr_category:
      vnp->data.ptrvalue = (CharPtr) pThis;
      break;
    case CddDescr_comment:
      vnp->data.ptrvalue = (CharPtr) pThis;
      break;
    case CddDescr_reference:
      vnp->data.ptrvalue = (ValNodePtr) pThis;
      break;
    case CddDescr_create_date:
      vnp->data.ptrvalue = (DatePtr) pThis;
      break;
    case CddDescr_tax_source:
      vnp->data.ptrvalue = (OrgRefPtr) pThis;
      break;
    case CddDescr_source:
      vnp->data.ptrvalue = (CharPtr) pThis;
      break;
    case CddDescr_status:
      vnp->data.intvalue = (Int4) iIval;
      break;
    default:
      vnp->data.ptrvalue = pThis;
      break;
  }
  ValNodeLink(&(pcdd->description),vnp);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Query the status of a CD                                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int4 LIBCALL CddGetStatus(CddPtr pcdd)
{
  CddDescrPtr    pCddesc;
  
  pCddesc = pcdd->description;
  while (pCddesc) {
    if (pCddesc->choice == CddDescr_status) {
      return(pCddesc->data.intvalue);
    }
    pCddesc = pCddesc->next;
  }
  return 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* report Errors in processing and exit immediately                          */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LIBCALL CddHtmlError(CharPtr cErrTxt) 
{
  printf("Content-type: text/html\n\n");
  printf("<h2>Error:</h2>\n");
  printf("<h3>%s</h3>\n",cErrTxt);
  exit(1);
}

void LIBCALL CddSevError(CharPtr cErrTxt) 
{
  printf(" Error: %s\n",cErrTxt);
  exit(1);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Find bioseqptr in Cdd internal sequence set                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
BioseqPtr LIBCALL CddFindSip(SeqIdPtr sip, SeqEntryPtr sep)
{
  BioseqPtr    bsp = NULL;
  BioseqSetPtr bssp;
  SeqEntryPtr  sepThis;

  sepThis = sep;
  
  while (sepThis) {
    if (sepThis->choice == 2) {
      bssp = sepThis->data.ptrvalue;
      bsp = CddFindSip(sip, bssp->seq_set);
    } else if (sepThis->choice == 1) {
      bsp = (BioseqPtr) sepThis->data.ptrvalue;
      if (CddSameSip(bsp->id, sip)) return(bsp);
    }
    sepThis = sepThis->next;
  }
  return(NULL);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* extract from a SeqEntryPtr a "minimum bioseq" with a bit of description   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
BioseqPtr LIBCALL CddExtractBioseq(SeqEntryPtr sep, SeqIdPtr sip)
{
  BioseqPtr    bsp, bspNew, bspTemp;
  BioseqSetPtr bssp;
  SeqIdPtr     sipNew;
  ValNodePtr   vnpThis, vnp = NULL;
  SeqAnnotPtr  sanp;

  bsp = BioseqFindInSeqEntry(sip,sep);
  sipNew = SeqIdDup(sip);
  bspNew = (BioseqPtr) BioseqCopy(sipNew,bsp->id,0,bsp->length-1,0,FALSE);
  if (sep->choice == 2) {
    bssp = sep->data.ptrvalue;
    vnp = bssp->descr;
    bssp->descr = NULL;
  } else if (sep->choice == 1) {
    bspTemp = sep->data.ptrvalue;
    vnp = bspTemp->descr;
    bspTemp->descr = NULL;
  }
  sanp = bsp->annot;
  vnpThis = bsp->descr;
  bsp->descr = NULL;
  ValNodeLink(&(vnp),vnpThis);
  bspNew->descr = vnp;
  if (sanp) {
    bspNew->annot = sanp;
    bsp->annot = NULL; 
  }
  return(bspNew);
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Cdd specific alignment format converters                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* make a DenseSeg Alignment-Set turning each diag into a separate alignment */
/*---------------------------------------------------------------------------*/
SeqAlignPtr LIBCALL CddMSLDenDiagToMSLDenSeg(SeqAlignPtr salp)
{
  SeqAlignPtr         salpnew, salphead, salptail;
  DenseDiagPtr        ddp;
  DenseSegPtr         dsp;
  Int4                lastm, lasts;

  if (!salp) return(NULL);
  salphead = NULL;
  salptail = NULL;
  while (salp) {
    if (salp->segtype) {
      if (salp->segtype != SAS_DENDIAG) {
        CddSevError("CddMSLDenDiagToMSLDenSeg: Wrong segtype specified!");
      }
    }
    if (salp->dim) {
      if (salp->dim != 2) {
        CddSevError("CddMSLDenDiagToMSLDenSeg: Expect alignments of dimension 2!");
      }
    }
    ddp = salp->segs;
    lastm = -1;
    lasts = -1;
    while (ddp) {
      salpnew = SeqAlignNew();
      if (!salphead) salphead = salpnew;
      salpnew->type = SAT_PARTIAL;
      salpnew->segtype = SAS_DENSEG;
      salpnew->dim = 2;
      salpnew->score = ddp->scores;
      dsp = DenseSegNew();
      salpnew->segs = dsp;
      dsp->dim = 2;
      dsp->numseg = 1;
      dsp->ids = ddp->id;
      dsp->starts = ddp->starts;
      dsp->lens = &(ddp->len);
      dsp->strands = ddp->strands;
      dsp->scores = ddp->scores;
      if (salptail) salptail->next = salpnew;
      salptail = salpnew;
      ddp = ddp->next;
    }
    salp = salp->next;
  }
  return(salphead);
}

/*---------------------------------------------------------------------------*/
/* convert a dendiag pairwise alignment set into a multiple dendiag alignment*/
/* if possible - i.e. check that the number of segments and their extents on */
/* the common master is the same for all pairwise alignments                 */
/*---------------------------------------------------------------------------*/
SeqAlignPtr LIBCALL CddMSLDenDiagToMULDenDiag(SeqAlignPtr salp)
{
  SeqAlignPtr     salpHead;
  SeqAlignPtr     salpNew;
  Int4            iNumSegs, iCnt, iDim;
  Int4            *iSegStart;
  Int4            *iSegLen;
  DenseDiagPtr    ddp, ddpNew = NULL, ddpTail = NULL;
  Boolean         bIsOk = TRUE;
  SeqIdPtr        sipNew, sipTail;

  if (!salp) return NULL;
  if (salp->dim && salp->dim != 2) {
    printf(" CddMSLDenDiagToMulDenDiag Warning: Can't convert alignment of dim!=2\n");
    return(salp);
  }
  iNumSegs = 0;
  ddp = (DenseDiagPtr) salp->segs;
  while (ddp) {
    iNumSegs++;
    ddp = ddp->next;
  }
  iSegStart = calloc(iNumSegs,sizeof(Int4));
  iSegLen = calloc(iNumSegs,sizeof(Int4));
  iCnt = 0;
  ddp = (DenseDiagPtr) salp->segs;
  while (ddp) {
    iSegStart[iCnt] = ddp->starts[0];
    iSegLen[iCnt] = ddp->len;    
    iCnt++;
    ddp = ddp->next;
  }
  iDim = 1;
  salpHead = salp;
  while (salpHead) {
    iDim++;
    if (iDim > 2) {
      ddp = (DenseDiagPtr) salpHead->segs;
      iCnt = 0;
      while (ddp) {
        if (iCnt > iNumSegs || ddp->starts[0] != iSegStart[iCnt] ||
            ddp->len != iSegLen[iCnt]) {
          bIsOk = FALSE;
          break;
        }
        iCnt++;  
        ddp = ddp->next;
      }
      if (iCnt != iNumSegs) {
        bIsOk = FALSE;
        break;
      }
    }
    salpHead = salpHead->next;
  }
  
  if (!bIsOk) {
    printf(" CddMSLDenDiagToMulDenDiag Warning: Can't convert alignment!\n");
    return(salp);
  }
  salpNew = SeqAlignNew();
  salpNew->type = salp->type;
  salpNew->segtype = salp->segtype;
  salpNew->dim = iDim;

  sipNew = NULL;
  salpHead = salp;
  while (salpHead) {
    ddp = (DenseDiagPtr)salpHead->segs;
    if (!sipNew) {
      sipNew = SeqIdDup(ddp->id);
      sipNew->next = SeqIdDup(ddp->id->next);
      sipTail = sipNew->next;
    } else {
      sipTail->next = SeqIdDup(ddp->id->next);
      sipTail = sipTail->next;
    }
    salpHead = salpHead->next;
  }

  salpHead = salp;
  iCnt = 0;
  while (salpHead) {
    ddp = (DenseDiagPtr) salpHead->segs;
    while (ddp) {
      if (!ddpNew) {
        ddpNew = (DenseDiagPtr) DenseDiagNew();
        ddpNew->dim = iDim;
        ddpNew->id = sipNew;
        ddpNew->starts = calloc(iDim,sizeof(Int4));
        ddpNew->starts[0]=ddp->starts[0];
        ddpNew->starts[1]=ddp->starts[1];
        ddpNew->len=ddp->len;
        if (!ddpTail) {
          salpNew->segs = ddpNew;
          ddpTail = ddpNew;
        } else {
          ddpTail->next = ddpNew;
          ddpTail = ddpTail->next;
        }
      } else {
        ddpNew->starts[iCnt+1]=ddp->starts[1];
      }
      if (ddpTail->next) {
        ddpTail = ddpTail->next;
        ddpNew = ddpTail;
      } else {
        ddpNew = NULL;
      }
      ddp = ddp->next;
      if (!ddp) {
        ddpNew = (DenseDiagPtr) salpNew->segs;
        ddpTail = ddpNew;
      }
    }
    iCnt++;
    salpHead = salpHead->next;
  }
  free(iSegLen);
  free(iSegStart);
  return(salpNew);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Modified after PGPReadBlastOptions (from blastpgp.c)                      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
PGPBlastOptionsPtr LIBCALL CddGetBlastOptions(BioseqPtr bsp)
{
  PGPBlastOptionsPtr  bop;
  BLAST_OptionsBlkPtr options;

  bop = MemNew (sizeof(PGPBlastOptions));
  bop->blast_database = StringSave("nr");
  bop->number_of_descriptions = 500;
  bop->number_of_alignments = 500;
  bop->believe_query = FALSE;
  bop->query_bsp = bsp;
  options = BLASTOptionNew("blastp",TRUE);
  bop->options = options;
  BLASTOptionSetGapParams(options,"BLOSUM62",0,0);
  bop->align_options = 0;
  bop->print_options = 0;
  options->required_start = 0;
  options->required_end = -1;
  options->window_size = 40;
  options->threshold_second = 11;
  options->dropoff_2nd_pass = (int)7.0;
  options->expect_value = 10.0;
  options->hitlist_size = 500;
  options->two_pass_method = FALSE;
  options->multiple_hits_only = TRUE;
  options->gap_open = 11;
  options->gap_extend = 1;
  options->decline_align = INT2_MAX;
  options->gap_x_dropoff = 15;
  options->gap_x_dropoff_final =25;
  options->gap_trigger = 22.0;
  options->filter_string = StringSave("F");
  options->number_of_cpus = 1;
  options->isPatternSearch = FALSE;
  options->ethresh = 0.000;
  options->pseudoCountConst = 10;
  options->maxNumPasses = 1;
  options->hsp_range_max = 0;
  options->block_width = 33;
  options = BLASTOptionValidate(options,"blastp");
  bop->fake_bsp = BioseqNew();
  bop->fake_bsp->descr = bop->query_bsp->descr;
  bop->fake_bsp->repr = bop->query_bsp->repr;
  bop->fake_bsp->mol = bop->query_bsp->mol;
  bop->fake_bsp->length = bop->query_bsp->length;
  bop->fake_bsp->seq_data_type = bop->query_bsp->seq_data_type;
  bop->fake_bsp->seq_data = bop->query_bsp->seq_data;
  bop->fake_bsp->id = SeqIdDup(bop->query_bsp->id);
  bop->query_bsp->id = SeqIdSetFree(bop->query_bsp->id);
  return(bop);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Generation of position-specific scoring matrices for database scanning    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LIBCALL CddDenDiagCposComputation(SeqAlignPtr listOfSeqAligns, BioseqPtr bsp,
                                       BioseqPtr bspF, CddPtr pcdd)
{
    Int4                numalign, numseq;    /*number of alignments and seqs */
    posSearchItems      *posSearch;          /*holds position-specific info  */
    compactSearchItems  *compactSearch = NULL;
    BLAST_ScoreBlkPtr   sbp;
    SeqLocPtr           private_slp = NULL;
    SeqPortPtr          spp = NULL;
    Uint1Ptr            query_seq, query_seq_start;
    Uint1               residue;
    Int4                index, a, KarlinReturn, array_size, iNew;
    Nlm_FloatHiPtr      lambda, K, H;
    Int4Ptr             open, extend;
    BLAST_ResFreqPtr    stdrfp;
    Int2                iStatus;
    ValNodePtr          error_return = NULL;
    ValNodePtr          ColumnHead;
    ValNodePtr          newRow, RowHead;
    ValNodePtr          newLetter, LetterHead;
    SeqCodeTablePtr     sctp;
    Char                ckptFileName[PATH_MAX];
    Char                cseqFileName[PATH_MAX];
    FILE                *fp;
    
    numalign = CddCountDenDiagSeqAligns(listOfSeqAligns, &numseq);
    posSearch = (posSearchItems *) MemNew(1*sizeof(posSearchItems));
    sbp = BLAST_ScoreBlkNew(Seq_code_ncbistdaa,1);
    sbp->read_in_matrix = TRUE;
    sbp->protein_alphabet = TRUE;
    sbp->posMatrix = NULL;
    sbp->number_of_contexts = 1;
    iStatus = BlastScoreBlkMatFill(sbp,"BLOSUM62");
    compactSearch = compactSearchNew(compactSearch);
    compactSearch->qlength = bsp->length;
    compactSearch->alphabetSize = sbp->alphabet_size;
    compactSearch->matrix = sbp->matrix;
    compactSearch->gapped_calculation = TRUE;
    compactSearch->pseudoCountConst = 10;
    compactSearch->ethresh = 0.001;

/*---------------------------------------------------------------------------*/
/* get the query sequence                                                    */
/*---------------------------------------------------------------------------*/
    private_slp = SeqLocIntNew(0, bspF->length-1 , Seq_strand_plus, SeqIdFindBest(bspF->id, SEQID_GI));
		spp = SeqPortNewByLoc(private_slp, Seq_code_ncbistdaa);
    SeqPortSet_do_virtual(spp, TRUE);
		query_seq_start = (Uint1Ptr) MemNew(((bspF->length)+2)*sizeof(Char));
		query_seq_start[0] = NULLB;
		query_seq = query_seq_start+1;
    index=0;
    while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
      if (IS_residue(residue)) {
				query_seq[index] = residue;
				index++;
			}
		}
		query_seq[index] = NULLB;
		spp = SeqPortFree(spp);
		private_slp = SeqLocFree(private_slp);
    compactSearch->query = query_seq;

    BlastScoreBlkFill(sbp,(CharPtr) compactSearch->query,compactSearch->qlength,0);

    sbp->kbp_gap_std[0] = BlastKarlinBlkCreate();
    KarlinReturn = BlastKarlinBlkGappedCalc(sbp->kbp_gap_std[0],11,1,sbp->name,&error_return);
    if (1 == KarlinReturn) {
      BlastErrorPrint(error_return);
    }
    sbp->kbp_gap_psi[0] = BlastKarlinBlkCreate();
    KarlinReturn = BlastKarlinBlkGappedCalc(sbp->kbp_gap_psi[0],11,1,sbp->name,&error_return);
    if (1 == KarlinReturn) {
      BlastErrorPrint(error_return);
    }

    array_size = BlastKarlinGetMatrixValues(sbp->name,&open,&extend,&lambda,&K,&H,NULL);
    if (array_size > 0) {
      for (index = 0; index < array_size; index++) {
        if (open[index] == INT2_MAX && extend[index] == INT2_MAX) {
          sbp->kbp_ideal = BlastKarlinBlkCreate();
          sbp->kbp_ideal->Lambda = lambda[index];
          sbp->kbp_ideal->K = K[index];
          sbp->kbp_ideal->H = H[index];
        }
      }
      MemFree(open);
      MemFree(extend);
      MemFree(lambda);
      MemFree(K);
      MemFree(H);
    }
    if (sbp->kbp_ideal == NULL) {
      sbp->kbp_ideal = BlastKarlinBlkStandardCalcEx(sbp);
    }
    compactSearch->lambda       = sbp->kbp_gap_std[0]->Lambda;
    compactSearch->kbp_std      = sbp->kbp_std;
    compactSearch->kbp_psi      = sbp->kbp_psi;
    compactSearch->kbp_gap_std  = sbp->kbp_gap_std;
    compactSearch->kbp_gap_psi  = sbp->kbp_gap_psi;
    compactSearch->lambda_ideal = sbp->kbp_ideal->Lambda;
    compactSearch->K_ideal      = sbp->kbp_ideal->K;

    stdrfp = BlastResFreqNew(sbp);
    BlastResFreqStdComp(sbp,stdrfp); 
    compactSearch->standardProb = MemNew(compactSearch->alphabetSize * sizeof(Nlm_FloatHi));
    for(a = 0; a < compactSearch->alphabetSize; a++)
      compactSearch->standardProb[a] = stdrfp->prob[a];
    stdrfp = BlastResFreqDestruct(stdrfp);
    posSearch->posInformation = NULL;
/*---------------------------------------------------------------------------*/
/* numseq is replaced with numalign (last argument) - each alignment is a    */
/* set of diags and represents an independent sequence/domain                */
/*---------------------------------------------------------------------------*/
    CddposAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numalign);
    CddfindDenseDiagThreshSequences(posSearch, listOfSeqAligns, numalign, numalign);
    sbp->kbp     = sbp->kbp_psi;
    sbp->kbp_gap = sbp->kbp_gap_psi;
    CddposInitializeInformation(posSearch,sbp,compactSearch);
    CddposDenseDiagDemographics(posSearch, compactSearch, listOfSeqAligns);
    CddposPurgeMatches(posSearch, compactSearch);
    CddposComputeExtents(posSearch, compactSearch);
    CddposComputeSequenceWeights(posSearch, compactSearch);
    CddposCheckWeights(posSearch, compactSearch);
    posSearch->posFreqs = (Nlm_FloatHi **) CddposComputePseudoFreqs(posSearch, compactSearch, TRUE);
    CddposFreqsToMatrix(posSearch,compactSearch);
    CddposScaling(posSearch, compactSearch);

    sctp = SeqCodeTableFind(Seq_code_ncbistdaa);
    LetterHead = NULL;
    for (a=0;a<compactSearch->alphabetSize;a++) {
      newLetter = ValNodeNew(NULL);
      newLetter->data.ptrvalue = MemNew(2*sizeof(Char));
      Nlm_StrNCpy(newLetter->data.ptrvalue,&(sctp->letters[a]),1);
      ValNodeLink(&(LetterHead),newLetter);
    }

    pcdd->posfreq = (MatrixPtr) MatrixNew();
    pcdd->posfreq->ncolumns = compactSearch->qlength;
    pcdd->posfreq->nrows = compactSearch->alphabetSize;
    pcdd->posfreq->scale_factor = 10000;
    pcdd->posfreq->row_labels = LetterHead;    
    ColumnHead = NULL;
    for (index = 0; index<compactSearch->qlength;index++) {
      RowHead = NULL;
      for (a=0;a<compactSearch->alphabetSize;a++) {
         newRow = ValNodeNew(NULL);
         iNew = (Int4) (0.5 + (Nlm_FloatHi) pcdd->posfreq->scale_factor * posSearch->posFreqs[index][a]);
         newRow->data.intvalue = iNew;
           ValNodeLink(&(ColumnHead),newRow);
      }
    }
    pcdd->posfreq->columns = ColumnHead;

    pcdd->scoremat = (MatrixPtr) MatrixNew();
    pcdd->scoremat->ncolumns = compactSearch->qlength;
    pcdd->scoremat->nrows = compactSearch->alphabetSize;
    pcdd->scoremat->scale_factor = 1;
    pcdd->scoremat->row_labels = LetterHead;    
    ColumnHead = NULL;
    for (index = 0; index<compactSearch->qlength;index++) {
      RowHead = NULL;
      for (a=0;a<compactSearch->alphabetSize;a++) {
         newRow = ValNodeNew(NULL);
         iNew = (Int4) posSearch->posMatrix[index][a];
         newRow->data.intvalue = iNew;
           ValNodeLink(&(ColumnHead),newRow);
      }
    }
    pcdd->scoremat->columns = ColumnHead;

/*---------------------------------------------------------------------------*/
/* Construct name for checkpoint file                                        */
/*---------------------------------------------------------------------------*/
    strcpy(ckptFileName,pcdd->name);
    strcat(ckptFileName,CKPTEXT);

    CddposTakeCheckpoint(posSearch, compactSearch, ckptFileName, &error_return);
    strcpy(cseqFileName,pcdd->name);
    strcat(cseqFileName,CSEQEXT);
  
    fp = FileOpen(cseqFileName, "w");
    if (NULL == fp) {
      BlastConstructErrorMessage("CddDenDiagCposComputation", "Could not open fasta file", 1, &error_return);
    }
    BioseqToFasta(bsp,fp,FALSE);
    FileClose(fp);
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Generation of position-specific scoring matrices for database scanning    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LIBCALL CddCposComputation(SeqAlignPtr listOfSeqAligns, BioseqPtr bsp,
                                CddPtr pcdd)
{
    Int4                numalign, numseq;    /*number of alignments and seqs */
    posSearchItems      *posSearch;          /*holds position-specific info  */
    compactSearchItems  *compactSearch = NULL;
    BLAST_ScoreBlkPtr   sbp;
    SeqLocPtr           private_slp = NULL;
    SeqPortPtr          spp = NULL;
    Uint1Ptr            query_seq, query_seq_start;
    Uint1               residue;
    Int4                index, a, KarlinReturn, array_size, iNew;
    Nlm_FloatHiPtr      lambda, K, H;
    Int4Ptr             open, extend;
    BLAST_ResFreqPtr    stdrfp;
    Int2                iStatus;
    ValNodePtr          error_return = NULL;
    ValNodePtr          ColumnHead;
    ValNodePtr          newRow, RowHead;
    ValNodePtr          newLetter, LetterHead;
    SeqCodeTablePtr     sctp;
    Char                ckptFileName[PATH_MAX];
    Char                cseqFileName[PATH_MAX];
    FILE                *fp;

    numalign = CddCountSeqAligns(listOfSeqAligns, &numseq);
    posSearch = (posSearchItems *) MemNew(1*sizeof(posSearchItems));
    sbp = BLAST_ScoreBlkNew(Seq_code_ncbistdaa,1);
    sbp->read_in_matrix = TRUE;
    sbp->protein_alphabet = TRUE;
    sbp->posMatrix = NULL;
    sbp->number_of_contexts = 1;
    iStatus = BlastScoreBlkMatFill(sbp,"BLOSUM62");
    compactSearch = compactSearchNew(compactSearch);
    compactSearch->qlength = bsp->length;
    compactSearch->alphabetSize = sbp->alphabet_size;
    compactSearch->matrix = sbp->matrix;
    compactSearch->gapped_calculation = TRUE;
    compactSearch->pseudoCountConst = 10;
    compactSearch->ethresh = 0.001;

/*---------------------------------------------------------------------------*/
/* get the query sequence                                                    */
/*---------------------------------------------------------------------------*/
    private_slp = SeqLocIntNew(0, bsp->length-1 , Seq_strand_plus, SeqIdFindBest(bsp->id, SEQID_GI));
		spp = SeqPortNewByLoc(private_slp, Seq_code_ncbistdaa);
    SeqPortSet_do_virtual(spp, TRUE);
		query_seq_start = (Uint1Ptr) MemNew(((bsp->length)+2)*sizeof(Char));
		query_seq_start[0] = NULLB;
		query_seq = query_seq_start+1;
    index=0;
    while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
      if (IS_residue(residue)) {
				query_seq[index] = residue;
				index++;
			}
		}
		query_seq[index] = NULLB;
		spp = SeqPortFree(spp);
		private_slp = SeqLocFree(private_slp);
    compactSearch->query = query_seq;

    BlastScoreBlkFill(sbp,(CharPtr) compactSearch->query,compactSearch->qlength,0);

    sbp->kbp_gap_std[0] = BlastKarlinBlkCreate();
    KarlinReturn = BlastKarlinBlkGappedCalc(sbp->kbp_gap_std[0],11,1,sbp->name,&error_return);
    if (1 == KarlinReturn) {
      BlastErrorPrint(error_return);
    }
    sbp->kbp_gap_psi[0] = BlastKarlinBlkCreate();
    KarlinReturn = BlastKarlinBlkGappedCalc(sbp->kbp_gap_psi[0],11,1,sbp->name,&error_return);
    if (1 == KarlinReturn) {
      BlastErrorPrint(error_return);
    }

    array_size = BlastKarlinGetMatrixValues(sbp->name,&open,&extend,&lambda,&K,&H,NULL);
    if (array_size > 0) {
      for (index = 0; index < array_size; index++) {
        if (open[index] == INT2_MAX && extend[index] == INT2_MAX) {
          sbp->kbp_ideal = BlastKarlinBlkCreate();
          sbp->kbp_ideal->Lambda = lambda[index];
          sbp->kbp_ideal->K = K[index];
          sbp->kbp_ideal->H = H[index];
        }
      }
      MemFree(open);
      MemFree(extend);
      MemFree(lambda);
      MemFree(K);
      MemFree(H);
    }
    if (sbp->kbp_ideal == NULL) {
      sbp->kbp_ideal = BlastKarlinBlkStandardCalcEx(sbp);
    }
    compactSearch->lambda       = sbp->kbp_gap_std[0]->Lambda;
    compactSearch->kbp_std      = sbp->kbp_std;
    compactSearch->kbp_psi      = sbp->kbp_psi;
    compactSearch->kbp_gap_std  = sbp->kbp_gap_std;
    compactSearch->kbp_gap_psi  = sbp->kbp_gap_psi;
    compactSearch->lambda_ideal = sbp->kbp_ideal->Lambda;
    compactSearch->K_ideal      = sbp->kbp_ideal->K;

    stdrfp = BlastResFreqNew(sbp);
    BlastResFreqStdComp(sbp,stdrfp); 
    compactSearch->standardProb = MemNew(compactSearch->alphabetSize * sizeof(Nlm_FloatHi));
    for(a = 0; a < compactSearch->alphabetSize; a++)
      compactSearch->standardProb[a] = stdrfp->prob[a];
    stdrfp = BlastResFreqDestruct(stdrfp);
    posSearch->posInformation = NULL;
/*---------------------------------------------------------------------------*/
/* numseq is replaced with numalign (last argument) - each alignment is a    */
/* single segment and represents a separate sequence                         */
/*---------------------------------------------------------------------------*/
    CddposAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numalign);
    CddfindThreshSequences(posSearch, listOfSeqAligns, numalign, numalign);
    sbp->kbp     = sbp->kbp_psi;
    sbp->kbp_gap = sbp->kbp_gap_psi;
    CddposInitializeInformation(posSearch,sbp,compactSearch);
    CddposDemographics(posSearch, compactSearch, listOfSeqAligns);
    CddposPurgeMatches(posSearch, compactSearch);
    CddposComputeExtents(posSearch, compactSearch);
    CddposComputeSequenceWeights(posSearch, compactSearch);
    CddposCheckWeights(posSearch, compactSearch);
    posSearch->posFreqs = (Nlm_FloatHi **) CddposComputePseudoFreqs(posSearch, compactSearch, TRUE);
    CddposFreqsToMatrix(posSearch,compactSearch);
    CddposScaling(posSearch, compactSearch);

    sctp = SeqCodeTableFind(Seq_code_ncbistdaa);
    LetterHead = NULL;
    for (a=0;a<compactSearch->alphabetSize;a++) {
      newLetter = ValNodeNew(NULL);
      newLetter->data.ptrvalue = MemNew(2*sizeof(Char));
      Nlm_StrNCpy(newLetter->data.ptrvalue,&(sctp->letters[a]),1);
      ValNodeLink(&(LetterHead),newLetter);
    }

    pcdd->posfreq = (MatrixPtr) MatrixNew();
    pcdd->posfreq->ncolumns = compactSearch->qlength;
    pcdd->posfreq->nrows = compactSearch->alphabetSize;
    pcdd->posfreq->scale_factor = 10000;
    pcdd->posfreq->row_labels = LetterHead;    
    ColumnHead = NULL;
    for (index = 0; index<compactSearch->qlength;index++) {
      RowHead = NULL;
      for (a=0;a<compactSearch->alphabetSize;a++) {
         newRow = ValNodeNew(NULL);
         iNew = (Int4) (0.5 + (Nlm_FloatHi) pcdd->posfreq->scale_factor * posSearch->posFreqs[index][a]);
         newRow->data.intvalue = iNew;
           ValNodeLink(&(ColumnHead),newRow);
      }
    }
    pcdd->posfreq->columns = ColumnHead;

    pcdd->scoremat = (MatrixPtr) MatrixNew();
    pcdd->scoremat->ncolumns = compactSearch->qlength;
    pcdd->scoremat->nrows = compactSearch->alphabetSize;
    pcdd->scoremat->scale_factor = 1;
    pcdd->scoremat->row_labels = LetterHead;    
    ColumnHead = NULL;
    for (index = 0; index<compactSearch->qlength;index++) {
      RowHead = NULL;
      for (a=0;a<compactSearch->alphabetSize;a++) {
         newRow = ValNodeNew(NULL);
         iNew = (Int4) posSearch->posMatrix[index][a];
         newRow->data.intvalue = iNew;
           ValNodeLink(&(ColumnHead),newRow);
      }
    }
    pcdd->scoremat->columns = ColumnHead;

/*---------------------------------------------------------------------------*/
/* Construct name for checkpoint file                                        */
/*---------------------------------------------------------------------------*/
    strcpy(ckptFileName,pcdd->name);
    strcat(ckptFileName,CKPTEXT);

    CddposTakeCheckpoint(posSearch, compactSearch, ckptFileName, &error_return);
    strcpy(cseqFileName,pcdd->name);
    strcat(cseqFileName,CSEQEXT);
  
    fp = FileOpen(cseqFileName, "w");
    if (NULL == fp) {
      BlastConstructErrorMessage("CddCposComputation", "Could not open fasta file", 1, &error_return);
    }
    BioseqToFasta(bsp,fp,FALSE);
    FileClose(fp);
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Count the number of seqAligns in a list (returned) and count the number   */
/* of distinct target sequences represented (passed back in numSequences);   */
/* Important assumption: All SeqAligns with  the same target sequence are    */
/* consecutive in the list                                                   */
/* ONLY works for lists of Seqaligns containing DenseSegs!!!                 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int4 LIBCALL CddCountSeqAligns(SeqAlignPtr listOfSeqAligns, Int4 * numSequences)
{
  SeqAlignPtr curSeqAlign, prevSeqAlign;
  Int4        seqAlignCounter;
  DenseSegPtr curSegs;
  SeqIdPtr    curId, prevId;

  seqAlignCounter = 0;
  *numSequences = 0;
  curSeqAlign = listOfSeqAligns;
  prevSeqAlign = NULL;
  while (NULL != curSeqAlign) {
    curSegs = (DenseSegPtr) curSeqAlign->segs;
    if(curSegs->ids == NULL) break;
    curId = curSegs->ids->next; 
    seqAlignCounter++;
    if ((NULL == prevSeqAlign) || (!(SeqIdMatch(curId, prevId))))
      (*numSequences)++;
    prevSeqAlign = curSeqAlign;
    prevId = curId;
    curSeqAlign = curSeqAlign->next;
  }
  return(seqAlignCounter);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Count the number of seqAligns in a list (returned) and count the number   */
/* of distinct target sequences represented (passed back in numSequences);   */
/* Important assumption: All SeqAligns with  the same target sequence are    */
/* consecutive in the list                                                   */
/* ONLY works for lists of Seqaligns containing DenseDiags!!!                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int4 LIBCALL CddCountDenDiagSeqAligns(SeqAlignPtr listOfSeqAligns, Int4 * numSequences)
{
  SeqAlignPtr  curSeqAlign;
  Int4         seqAlignCounter;
  DenseDiagPtr curSegs;

  seqAlignCounter = 0;
  *numSequences = 0;
  curSeqAlign = listOfSeqAligns;
  while (NULL != curSeqAlign) {
    curSegs = (DenseDiagPtr) curSeqAlign->segs;
    if(curSegs->id == NULL) break;
    seqAlignCounter++;
    (*numSequences)++;
    curSeqAlign = curSeqAlign->next;
  }
  return(seqAlignCounter);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* assign the range of the master sequence involved in alignments            */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LIBCALL CddAssignProfileRange(CddPtr pcdd, SeqIdPtr sip)
{
  SeqAlignPtr  salp;
  DenseDiagPtr ddp;
  SeqIntPtr    sintp;
  Int4         iMax = -10000000;
  Int4         iMin =  10000000;

  sintp = (SeqIntPtr) SeqIntNew();

  salp = (SeqAlignPtr) pcdd->seqannot->data;

  while (salp) {
    ddp = (DenseDiagPtr) salp->segs;
    while (ddp) {
      if (ddp->starts[0] < iMin) iMin = ddp->starts[0];
      if ((ddp->starts[0] + ddp->len - 1) > iMax) 
        iMax = ddp->starts[0] + ddp->len - 1;
      ddp = ddp->next;
    }
    salp = salp->next;
  }
  sintp->from = iMin;
  sintp->to = iMax;
  sintp->id = (SeqIdPtr) SeqIdDup(sip);
  if (iMax >= iMin) {
    pcdd->profile_range = (struct seqint PNTR) sintp;
  } else {
    printf(" CddAssignProfileRange: Warning: Can not assign alignment range!\n");
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* return a pointer to a specific bioseq from a seq-entry set                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
BioseqPtr LIBCALL CddFindBioseqInSeqEntry(SeqEntryPtr sep, SeqIdPtr sip)
{
  SeqEntryPtr  sepThis;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;

  sepThis = sep;

  while (sepThis) {
    if (IS_Bioseq(sepThis)) {
      bsp = sep->data.ptrvalue;
      if (CddSameSip(bsp->id, sip)) return(bsp);
    } else if (IS_Bioseq_set(sepThis)) {
      bssp = sep->data.ptrvalue;
      bsp = CddFindBioseqInSeqEntry(bssp->seq_set,sip);
      if (CddSameSip(bsp->id, sip)) return(bsp);
    }
    sepThis = sepThis->next;
  }
  return(NULL);
}

/*---------------------------------------------------------------------------*/
/* rips out and returns a PDBSeqId from a SeqId                              */
/*---------------------------------------------------------------------------*/
PDBSeqIdPtr LIBCALL CddGetPdbSeqId(SeqIdPtr sip)
                    /* may need to be modified according to how bioseq id is */
{
  SeqIdPtr    seq_id = NULL;
  PDBSeqIdPtr pdb_seq_id = NULL;

  seq_id = sip;
  while(seq_id != NULL){
    if(seq_id->choice == 15) {
      pdb_seq_id = seq_id->data.ptrvalue;
      break;
    }
    seq_id = seq_id->next;
  }
  return(pdb_seq_id);
}


/*---------------------------------------------------------------------------*/
/* test whether two SeqIds are the same, only considers PDB-Ids and integers */
/* as well as local sequence-id's                                            */
/*---------------------------------------------------------------------------*/
Boolean LIBCALL CddSameSip(SeqIdPtr sip1, SeqIdPtr sip2)
{
  PDBSeqIdPtr psip1, psip2;
  SeqIdPtr    tsip1, tsip2;
  ObjectIdPtr oidp1, oidp2;

  if (sip1 == sip2) return(TRUE);
  tsip1 = sip1;
  while (tsip1) {
    tsip2 = sip2;
    while (tsip2) {
      if (tsip1->choice == tsip2->choice) {
        if (tsip1->choice == 15) {
          psip1 = (PDBSeqIdPtr) CddGetPdbSeqId(tsip1);
          psip2 = (PDBSeqIdPtr) CddGetPdbSeqId(tsip2);
          if (strcmp(psip1->mol,psip2->mol)==0 && psip1->chain==psip2->chain)
            return(TRUE);
        } else if (tsip1->choice == SEQID_LOCAL) {
          oidp1 = tsip1->data.ptrvalue;
          oidp2 = tsip2->data.ptrvalue;
          if (oidp1 && oidp1->id == oidp2->id) return(TRUE);
          if (strcmp(oidp1->str,oidp2->str)==0) return(TRUE);
        }else if (tsip1->data.intvalue == tsip2->data.intvalue) return(TRUE);
      }
      tsip2 = tsip2->next;
    }
    tsip1 = tsip1->next;
  }
  return(FALSE);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* reindex master of a MSL-DenSeg Alignment according to a n-terminal offset */
/* must allocate new arrays for "starts", these may just be pointers to DDG  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LIBCALL CddReindexMSLDenSegMaster(SeqAlignPtr salp, Int4 offset)
{
  SeqAlignPtr salpThis;
  DenseSegPtr dsp;
  Int4        *newstarts, i;

  if (salp->type == SAT_PARTIAL && salp->segtype == SAS_DENSEG) {
    salpThis = salp;
    while (salpThis) {
      dsp = salpThis->segs;
      if (dsp) {
        newstarts = (Int4 *)MemNew(dsp->dim * sizeof(Int4));
        for (i=0;i<dsp->dim;i++) newstarts[i] = dsp->starts[i];
        newstarts[0] -= offset;
        dsp->starts = newstarts;
      }
      salpThis = salpThis->next;
    }
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* reindex master of a MSL-DenDiag Alignment according to a n-terminal offset*/
/* must allocate new arrays for "starts", these may just be pointers to DDG  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LIBCALL CddReindexMSLDenDiagMaster(SeqAlignPtr salp, Int4 offset)
{
  SeqAlignPtr  salpThis;
  DenseDiagPtr ddp;

  if (salp->type == SAT_PARTIAL && salp->segtype == SAS_DENDIAG) {
    salpThis = salp;
    while (salpThis) {
      ddp = salpThis->segs;
      while (ddp) {
        ddp->starts[0] -= offset;
        ddp = ddp->next;
      }
      salpThis = salpThis->next;
    }
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* makes a copy of one master-slave dense-diag alignment - does NOT copy the */
/* identifiers. Should only be used for a reindexing the master for CposComp.*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SeqAlignPtr LIBCALL CddCopyMSLDenDiag(SeqAlignPtr salp)
{
  SeqAlignPtr  salpNew = NULL, salpHead = NULL, salpTail = NULL;
  DenseDiagPtr ddp, ddpNew, ddpTail;
  Int4         i;
  
  if (salp->type == SAT_PARTIAL && salp->segtype == SAS_DENDIAG) {
    while (salp) {
      salpNew = SeqAlignNew();
      salpNew->type = salp->type;
      salpNew->segtype = salp->segtype;
      salpNew->dim = salp->dim;
      salpNew->segs = NULL;
      ddp = salp->segs;
      while (ddp) {
        ddpNew = DenseDiagNew();
        ddpNew->dim = ddp->dim;
        ddpNew->id = ddp->id;
        ddpNew->starts = MemNew(ddpNew->dim * sizeof(Int4));
        for (i =0;i<ddpNew->dim;i++) ddpNew->starts[i] = ddp->starts[i];
        ddpNew->len = ddp->len;
        if (!salpNew->segs) {
          salpNew->segs = ddpNew;
          ddpTail = ddpNew;
          ddpTail->next = NULL;
        } else {
          ddpTail->next = ddpNew;
          ddpTail = ddpNew;
          ddpTail->next = NULL;
        }
        ddp = ddp->next;
      }
      if (!salpHead) {
        salpHead = salpNew;
        salpTail = salpHead;
        salpTail->next = NULL;
      } else {
        salpTail->next = salpNew;
        salpTail = salpNew;
        salpTail->next = NULL;
      }
      salp = salp->next;
    }
  }
  return salpHead;
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* allocate and free a new ExplicitAlignment object                          */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
CddExpAlignPtr CddExpAlignNew()
{
  CddExpAlignPtr pCDea;
  
  pCDea = MemNew(sizeof(CddExpAlign));
  if (!pCDea) return NULL;
  pCDea->length = 0;
  pCDea->ids = NULL;
  pCDea->adata = NULL;
  return(pCDea);
}

CddExpAlignPtr CddExpAlignFree(CddExpAlignPtr pCDea)
{
  if (!pCDea) return NULL;
  if (pCDea->adata) free(pCDea->adata);
  free(pCDea);
  return(pCDea);
}

void           CddExpAlignAlloc(CddExpAlignPtr pCDea, Int4 iLength)
{
  Int4    i;
  
  if (!pCDea || !iLength) return;
  pCDea->length = iLength;
  if (pCDea->adata) free(pCDea->adata);
  pCDea->adata = MemNew(sizeof(Int4)*pCDea->length);  
  for (i=0;i<pCDea->length;i++) pCDea->adata[i] = -1;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* reindex a m-s1, m-s2 alignment pair to a s1-s2 alignment pair             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
CddExpAlignPtr CddReindexExpAlign(CddExpAlignPtr pCDea1, Int4 newlength,
                                  CddExpAlignPtr pCDea2, Int4 iOuter, Int4 iInner)
{
  CddExpAlignPtr     pCDea = NULL;
  Int4               i, iP1, iP2;
  
  if (!pCDea1 || !pCDea2) return pCDea;
  pCDea = CddExpAlignNew();
  CddExpAlignAlloc(pCDea,newlength);
  for (i=0;i<pCDea1->length;i++) {
    iP1 = pCDea1->adata[i]; iP2 = pCDea2->adata[i];
    if (iP1 > -1 && iP2 > -1){
      if (iP1 >= newlength) {
        printf("Warning: index error between sequences %d and %d\n",
	        iInner, iOuter);
        CddSevError("Error while reindexing explicit alignments!");
      } else  pCDea->adata[iP1] = iP2;
    }
  }
  return pCDea;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* returns a pairwise % id for two sequences in a CD (ordered by alignmt)    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
FloatHi CddGetPairId(TrianglePtr pTri, Int4 idx1, Int4 idx2)
{
  Int4 i,  index;
  float    findex;
  ScorePtr psc;
  
  if (idx1 == idx2) return 100.0;
  if (idx1 > idx2) { i = idx1; idx1 = idx2; idx2 = i; }
  findex=(float)idx1*((float)pTri->nelements-((float)idx1+1.0)/2.0-1.0)+(float)idx2-1.0;
  index = (Int4) findex;
  if (index < 0) CddSevError("Impossible index in CddGetPairId!");
  psc = (ScorePtr) pTri->scores; i=0;
  while (psc && (i<index)) {
    i++; psc = psc->next;
  }
  if (psc) return psc->value.realvalue;
  else CddSevError("Ran out of scores in CddGetPairId!");
  return((FloatHi) 0.0);
}

static Boolean HitYet(Int4Ptr retlist, Int4 index, Int4 i)
{
  Int4 j;
  for (j=0;j<index;j++) {
    if (retlist[j]==i) {
      return TRUE; 
    }
  }
  return FALSE;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* return a list of indices for the n most dissimilar sequences in a CD      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int4Ptr CddMostDiverse(TrianglePtr pTri, Int4 length)
{
  Int4Ptr     retlist;
  Int4        index, winner, i, j, ncomp;
  FloatHi     sumcomp, mincomp;
  FloatHi     **iMat;
  ScorePtr    psc;
  
  if (!pTri) return NULL;
  if (length >= pTri->nelements) return NULL;
  iMat = (FloatHi **) calloc(pTri->nelements,sizeof (FloatHi *));
  for (i=0;i<pTri->nelements;i++) {
    iMat[i] = (FloatHi *) calloc(pTri->nelements,sizeof(FloatHi));
    iMat[i][i] = 100;
  }
  psc = (ScorePtr) pTri->scores;
  for (i=0;i<pTri->nelements;i++) {
    for (j=i+1;j<pTri->nelements;j++) {
      iMat[i][j] = iMat[j][i] = psc->value.realvalue;
      psc = psc->next;
    }
  }
  
  retlist = MemNew(length * sizeof(Int4));
  retlist[0] = 0;
  index = 0; for (index = 1; index < length; index++) {
    mincomp = 100.0, winner = -1;
    for (i=1;i<pTri->nelements;i++) {
      sumcomp = 0.0; ncomp = 0;
      if (!HitYet(retlist,index,i)) {
        for (j=0;j<pTri->nelements;j++) {
	  if (HitYet(retlist,index,j)) {
/*	    sumcomp += CddGetPairId(pTri, i, j); */
            sumcomp += iMat[i][j];
	    ncomp ++;
	  }
	}
	if (ncomp) sumcomp /= (FloatHi) ncomp;
	if (sumcomp < mincomp) {
	  mincomp = sumcomp; winner = i;
	}
      }
    }
    retlist[index] = winner;
  }
  for (i=0;i<pTri->nelements;i++) free(iMat[i]);
  free(iMat);
  return (retlist);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* return a list of indices for the n most similar sequences for the query   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int4Ptr CddMostSimilarToQuery(ScorePtr psc, Int4 length)
{
  Int4Ptr    retlist;
  Int4       winner, index, i, nelements;
  FloatHi    sumcomp, mincomp;
  ScorePtr   thispsc;
  
  if (!psc) return NULL;
  nelements = 0; thispsc = psc;
  while (thispsc) { nelements++; thispsc = thispsc->next; }
  if (length >= nelements) return NULL;
  retlist = MemNew(length * sizeof(Int4));
  retlist[0] = 0;
  for (index=1;index<length;index++) {
    mincomp = 0.0; winner = -1;
    i = 0;
    thispsc = psc->next;
    while (thispsc) {
      i++;    
      if (!HitYet(retlist,index,i)) {
        sumcomp = thispsc->value.realvalue;
        if (sumcomp > mincomp) {
	  mincomp = sumcomp;
	  winner = i;
	}
      }
      thispsc = thispsc->next;
    }
    retlist[index] = winner;
  }  
  return(retlist);
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* retrieve a sequence by trying to locate it in the CDD bioseq-set or other */
/* wise through the sequence/object-manager or ID1                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
BioseqPtr CddRetrieveBioseqById(SeqIdPtr sip, SeqEntryPtr sep)
{
  BioseqPtr   bsp = NULL;
  Int4        uid = 0;
/*  SeqEntryPtr sepThis; */
  Int2        retcode = 1;

  bsp = (BioseqPtr) CddFindSip(sip, sep);
  if (!bsp) {
    bsp = BioseqLockById(sip);
  }
  /*
  if (!bsp) {
    uid = ID1FindSeqId(sip);
    if (uid > 0) {
      sepThis = (SeqEntryPtr) ID1SeqEntryGet(uid,retcode);
      if (sepThis) {
        bsp = CddExtractBioseq(sepThis, sip);
      }
    }
  }
  */
  return bsp;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Calculate pairwise percentages of identical residues for sequences in a CD*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
TrianglePtr CddCalculateTriangle(CddPtr pcdd)
{
   BioseqPtr      bsp1 = NULL, bsp2, bsp3;
   BioseqSetPtr   bssp = NULL;
   CddExpAlignPtr pCDea, pCDea1, pCDea2;
   DenseDiagPtr   ddp;
   SeqAlignPtr    salpOuter, salpInner, salpHead;
   SeqIdPtr       sip1, sip2, sip3;
   SeqPortPtr     spp1, spp2, spp3;
   ScorePtr       pscHead=NULL, psc, pscTail;
   TrianglePtr    pTri;
   Uint1Ptr       buffer1, buffer2, buffer3;
   Int4           i, alen, nid, iOuter, iInner;
   

   if (!pcdd->seqannot) return NULL;
   salpHead = pcdd->seqannot->data;
   pTri = TriangleNew();
   pTri->nelements = 1;

/*---------------------------------------------------------------------------*/
/* Loop around Outer pointer to get master-slave alignments                  */
/*---------------------------------------------------------------------------*/
   salpOuter = salpHead; bsp1 = NULL; sip1 = NULL;
   while (salpOuter) {
     pTri->nelements++;
     ddp = salpOuter->segs;
     if (!sip1) sip1 = SeqIdDup(ddp->id);
     sip2 = SeqIdDup(ddp->id->next);
     if (!bssp) bssp = pcdd->sequences->data.ptrvalue;
     if (!bsp1) {
       bsp1 = (BioseqPtr) CddRetrieveBioseqById(sip1, bssp->seq_set);
       buffer1 = MemNew((bsp1->length)*sizeof(Uint1));
       spp1 = SeqPortNew(bsp1, 0, bsp1->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
       for (i=0; i<bsp1->length;i++) buffer1[i] = SeqPortGetResidue(spp1);
       spp1 = SeqPortFree(spp1);
     }
     bsp2 = NULL;
     bsp2 = (BioseqPtr) CddRetrieveBioseqById(sip2, bssp->seq_set);
     buffer2 = MemNew((bsp2->length)*sizeof(Uint1));
     spp2 = SeqPortNew(bsp2, 0, bsp2->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
     for (i=0; i<bsp2->length;i++) buffer2[i] = SeqPortGetResidue(spp2);
     spp2 = SeqPortFree(spp2);
     pCDea = CddExpAlignNew();
     CddExpAlignAlloc(pCDea,bsp1->length);
     while (ddp) {
       for (i=0;i<ddp->len;i++) {
         pCDea->adata[ddp->starts[0]+i]=ddp->starts[1]+i;
       }
       ddp = ddp->next;
     }
     alen = 0; nid = 0;
     for (i=0;i<bsp1->length;i++) {
       if (pCDea->adata[i] > -1) {
         alen++; if (buffer1[i] == buffer2[pCDea->adata[i]]) nid++;
       }
     }
     psc = ScoreNew(); psc->id=ObjectIdNew(); psc->id->str=StringSave("%id");
     psc->choice = 2;
     if (alen) psc->value.realvalue = 100.0 * (FloatHi) nid / (FloatHi) alen;
     else psc->value.realvalue = 0.0;
     psc->next = NULL;
     if (!pscHead) {
       pscHead = psc; pscTail = psc;
     } else {
       pscTail->next = psc; pscTail = psc; psc = NULL;
     }
     pCDea = CddExpAlignFree(pCDea);
     MemFree(buffer2); SeqIdFree(sip2);
     salpOuter = salpOuter->next;
   }
   
/*---------------------------------------------------------------------------*/
/* Loop around Outer and Inner Pointer to get slave-slave similarities       */
/*---------------------------------------------------------------------------*/
   salpOuter = salpHead; iOuter = 1;
   while (salpOuter->next) {
     iOuter++;
     ddp = salpOuter->segs;
     sip2 = SeqIdDup(ddp->id->next); bsp2 = NULL;
     bsp2 = (BioseqPtr) CddRetrieveBioseqById(sip2, bssp->seq_set);
     buffer2 = MemNew((bsp2->length)*sizeof(Uint1));
     spp2 = SeqPortNew(bsp2, 0, bsp2->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
     for (i=0; i<bsp2->length;i++) buffer2[i] = SeqPortGetResidue(spp2);
     spp2 = SeqPortFree(spp2);
     pCDea1 = CddExpAlignNew();
     CddExpAlignAlloc(pCDea1,bsp1->length);
     while (ddp) {
       for (i=0;i<ddp->len;i++) pCDea1->adata[ddp->starts[0]+i]=ddp->starts[1]+i;
       ddp = ddp->next;
     }
     salpInner = salpOuter->next;
     iInner = iOuter;
     while (salpInner) {
       iInner++;
       ddp = salpInner->segs;
       sip3 = SeqIdDup(ddp->id->next); bsp3 = NULL;
       bsp3 = (BioseqPtr) CddRetrieveBioseqById(sip3,bssp->seq_set);
       buffer3  =MemNew((bsp3->length)*sizeof(Uint1));
       spp3 = SeqPortNew(bsp3, 0, bsp3->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
       for (i=0; i<bsp3->length;i++) buffer3[i] = SeqPortGetResidue(spp3);
       spp3 = SeqPortFree(spp3);
       pCDea2 = CddExpAlignNew();
       CddExpAlignAlloc(pCDea2,bsp1->length);
       while (ddp) {
         for (i=0;i<ddp->len;i++) pCDea2->adata[ddp->starts[0]+i]=ddp->starts[1]+i;
         ddp = ddp->next;
       }
       pCDea = CddReindexExpAlign(pCDea1, bsp2->length, pCDea2, iOuter, iInner);
       alen = 0; nid = 0;
       for (i=0;i<bsp2->length;i++) {
         if (pCDea->adata[i] > -1) {
           alen++; if (buffer2[i] == buffer3[pCDea->adata[i]]) nid++;
         }
       }
       psc = ScoreNew(); psc->id=ObjectIdNew(); psc->id->str=StringSave("%id");
       psc->choice = 2;
       if (alen) psc->value.realvalue = 100.0 * (FloatHi) nid / (FloatHi) alen;
       else psc->value.realvalue = 0.0;
       psc->next = NULL;
       pscTail->next = psc; pscTail = psc; psc = NULL;
       pCDea2 = CddExpAlignFree(pCDea2);
       pCDea = CddExpAlignFree(pCDea);
       MemFree(buffer3); SeqIdFree(sip3); 
       salpInner = salpInner->next;
     }
     pCDea1 = CddExpAlignFree(pCDea1);
     MemFree(buffer2); SeqIdFree(sip2);
     salpOuter = salpOuter->next;
   }
   if (sip1) SeqIdFree(sip1); if (bsp1) bsp1 = NULL;   
   if (buffer1) MemFree(buffer1);
   pTri->scores = (struct score PNTR) pscHead;
   return(pTri);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Calculate pairwise similarities for the query added to a cd alignment     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
ScorePtr CddCalculateQuerySim(CddPtr pcdd, SeqAlignPtr salp)
{
   BioseqPtr      bsp1, bsp2, bsp3;
   BioseqSetPtr   bssp = NULL;
   CddExpAlignPtr pCDea, pCDea1, pCDea2;
   DenseDiagPtr   ddp;
   SeqAlignPtr    salpOuter, salpInner, salpHead;
   SeqIdPtr       sip1, sip2, sip3;
   SeqPortPtr     spp1, spp2, spp3;
   ScorePtr       pscHead=NULL, psc, pscTail;
   Uint1Ptr       buffer1, buffer2, buffer3;
   Int4           i, alen, nid;
   

   if (!salp || !pcdd) return NULL;
   salpHead = salp;
/*---------------------------------------------------------------------------*/
/* master-query alignment                                                    */
/*---------------------------------------------------------------------------*/
   salpOuter = salpHead; bsp1 = NULL; sip1 = NULL;
   pCDea = CddExpAlignNew();
   ddp = salpOuter->segs;
   sip1 = SeqIdDup(ddp->id);
   sip2 = SeqIdDup(ddp->id->next);
   bssp = pcdd->sequences->data.ptrvalue;
   bsp1 = (BioseqPtr) CddRetrieveBioseqById(sip1, bssp->seq_set);
   buffer1 = MemNew((bsp1->length)*sizeof(Uint1));
   spp1 = SeqPortNew(bsp1, 0, bsp1->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
   for (i=0; i<bsp1->length;i++) buffer1[i] = SeqPortGetResidue(spp1);
   spp1 = SeqPortFree(spp1);
   bsp2 = (BioseqPtr) CddRetrieveBioseqById(sip2, bssp->seq_set);
   buffer2 = MemNew((bsp2->length)*sizeof(Uint1));
   spp2 = SeqPortNew(bsp2, 0, bsp2->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
   for (i=0; i<bsp2->length;i++) buffer2[i] = SeqPortGetResidue(spp2);
   spp2 = SeqPortFree(spp2);
   CddExpAlignAlloc(pCDea,bsp1->length);
   while (ddp) {
     for (i=0;i<ddp->len;i++) {
       pCDea->adata[ddp->starts[0]+i]=ddp->starts[1]+i;
     }
     ddp = ddp->next;
   }
   alen = 0; nid = 0;
   for (i=0;i<bsp1->length;i++) {
     if (pCDea->adata[i] > -1) {
       alen++; if (buffer1[i] == buffer2[pCDea->adata[i]]) nid++;
     }
   }
   pCDea = CddExpAlignFree(pCDea);
   psc = ScoreNew(); psc->id=ObjectIdNew(); psc->id->str=StringSave("%id");
   psc->choice = 2;
   if (alen) psc->value.realvalue = 100.0 * (FloatHi) nid / (FloatHi) alen;
   else psc->value.realvalue = 0.0;
   psc->next = NULL;
   pscHead = psc; pscTail = psc; psc = NULL;
   MemFree(buffer2); SeqIdFree(sip2);
/*---------------------------------------------------------------------------*/
/* Loop around Inner Pointer to get query-slave similarities                 */
/*---------------------------------------------------------------------------*/
   salpOuter = salpHead;
   pCDea1 = CddExpAlignNew();
   ddp = salpOuter->segs;
   sip2 = SeqIdDup(ddp->id->next);
   bsp2 = (BioseqPtr) CddRetrieveBioseqById(sip2, bssp->seq_set);
   buffer2 = MemNew((bsp2->length)*sizeof(Uint1));
   spp2 = SeqPortNew(bsp2, 0, bsp2->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
   for (i=0; i<bsp2->length;i++) buffer2[i] = SeqPortGetResidue(spp2);
   spp2 = SeqPortFree(spp2);
   CddExpAlignAlloc(pCDea1,bsp1->length);
   while (ddp) {
     for (i=0;i<ddp->len;i++) pCDea1->adata[ddp->starts[0]+i]=ddp->starts[1]+i;
     ddp = ddp->next;
   }
   salpInner = salpOuter->next;
   while (salpInner) {
     pCDea2 = CddExpAlignNew();
     ddp = salpInner->segs;
     sip3 = SeqIdDup(ddp->id->next);
     bsp3 = (BioseqPtr) CddRetrieveBioseqById(sip3,bssp->seq_set);
     buffer3  =MemNew((bsp3->length)*sizeof(Uint1));
     spp3 = SeqPortNew(bsp3, 0, bsp3->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
     for (i=0; i<bsp3->length;i++) buffer3[i] = SeqPortGetResidue(spp3);
     spp3 = SeqPortFree(spp3);
     CddExpAlignAlloc(pCDea2,bsp1->length);
     while (ddp) {
       for (i=0;i<ddp->len;i++) pCDea2->adata[ddp->starts[0]+i]=ddp->starts[1]+i;
       ddp = ddp->next;
     }
     pCDea = CddReindexExpAlign(pCDea1, bsp2->length, pCDea2,0,0);
     alen = 0; nid = 0;
     for (i=0;i<bsp2->length;i++) {
       if (pCDea->adata[i] > -1) {
         alen++; if (buffer2[i] == buffer3[pCDea->adata[i]]) nid++;
       }
     }
     pCDea = CddExpAlignFree(pCDea);
     psc = ScoreNew(); psc->id=ObjectIdNew(); psc->id->str=StringSave("%id");
     psc->choice = 2;
     if (alen) psc->value.realvalue = 100.0 * (FloatHi) nid / (FloatHi) alen;
     else psc->value.realvalue = 0.0;
     psc->next = NULL;
     pscTail->next = psc; pscTail = psc; psc = NULL;
     pCDea2 = CddExpAlignFree(pCDea2);
     MemFree(buffer3); SeqIdFree(sip3); 
     salpInner = salpInner->next;
   }
   pCDea1 = CddExpAlignFree(pCDea1);
   MemFree(buffer2); SeqIdFree(sip2);
   if (sip1) SeqIdFree(sip1); if (bsp1) bsp1 = NULL;   
   if (buffer1) MemFree(buffer1);
   return(pscHead);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Calculate a weighted 50/50 consensus sequence and make it new master      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SeqAlignPtr LIBCALL CddConsensus(SeqAlignPtr salp,
                                 SeqEntryPtr sep,
                                 BioseqPtr   bsp)
{
  Int4Ptr  trunc_on_virtual;
  
  trunc_on_virtual = MemNew(sizeof(Int4)*bsp->length);

  free (trunc_on_virtual);

  return( (SeqAlignPtr) NULL);
}



static void CddCposCompPart1(SeqAlignPtr listOfSeqAligns, BioseqPtr bsp,
                             compactSearchItems* compactSearch, ValNodePtr* LetterHead,
                             posSearchItems* posSearch) {
/*---------------------------------------------------------------------------*/
/* first part of CddCposComp()                                               */
/* basically returns compactSearch, LetterHead, and posSearch                */
/*---------------------------------------------------------------------------*/
    Int4                numalign, numseq;    /*number of alignments and seqs */
    BLAST_ScoreBlkPtr   sbp;
    SeqLocPtr           private_slp = NULL;
    SeqPortPtr          spp = NULL;
    Uint1Ptr            query_seq, query_seq_start;
    Uint1               residue;
    Int4                index, a, KarlinReturn, array_size;
    Nlm_FloatHiPtr      lambda, K, H;
    Int4Ptr             open, extend;
    BLAST_ResFreqPtr    stdrfp;
    Int2                iStatus;
    ValNodePtr          error_return = NULL;
    ValNodePtr          newLetter;
    SeqCodeTablePtr     sctp;
    BioseqPtr           bspFake, temp_bsp;

    /* this is used for the DenseDiag calculations */
    bspFake  = (BioseqPtr) BioseqCopy(NULL,bsp->id,0,bsp->length-1,0,FALSE);

    if (listOfSeqAligns->segtype == SAS_DENDIAG) {
      numalign = CddCountDenDiagSeqAligns(listOfSeqAligns, &numseq);
    }
    else {
      numalign = CddCountSeqAligns(listOfSeqAligns, &numseq);
    }
    sbp = BLAST_ScoreBlkNew(Seq_code_ncbistdaa,1);
    sbp->read_in_matrix = TRUE;
    sbp->protein_alphabet = TRUE;
    sbp->posMatrix = NULL;
    sbp->number_of_contexts = 1;
    iStatus = BlastScoreBlkMatFill(sbp,"BLOSUM62");
    ASSERT(iStatus == 0);
    compactSearch->qlength = bsp->length;
    compactSearch->alphabetSize = sbp->alphabet_size;
    compactSearch->matrix = sbp->matrix;
    compactSearch->gapped_calculation = TRUE;
    compactSearch->pseudoCountConst = 10;
    compactSearch->ethresh = 0.001;

/*---------------------------------------------------------------------------*/
/* get the query sequence                                                    */
/*---------------------------------------------------------------------------*/
    temp_bsp = bsp;
    if (listOfSeqAligns->segtype == SAS_DENDIAG) {
      bsp = bspFake;  /* to make call compatible with CddDenDiagCposComputation */
    }
    private_slp = SeqLocIntNew(0, bsp->length-1 , Seq_strand_plus, SeqIdFindBest(bsp->id, SEQID_GI));
		spp = SeqPortNewByLoc(private_slp, Seq_code_ncbistdaa);
    SeqPortSet_do_virtual(spp, TRUE);
    query_seq_start = (Uint1Ptr) MemNew(((bsp->length)+2)*sizeof(Char));
    bsp = temp_bsp;
		query_seq_start[0] = NULLB;
    query_seq = query_seq_start+1;
    index=0;
    while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
      if (IS_residue(residue)) {
				query_seq[index] = residue;
				index++;
			}
		}
		query_seq[index] = NULLB;
		spp = SeqPortFree(spp);
		private_slp = SeqLocFree(private_slp);
    compactSearch->query = query_seq;

    BlastScoreBlkFill(sbp,(CharPtr) compactSearch->query,compactSearch->qlength,0);

    sbp->kbp_gap_std[0] = BlastKarlinBlkCreate();
    KarlinReturn = BlastKarlinBlkGappedCalc(sbp->kbp_gap_std[0],11,1,sbp->name,&error_return);
    if (1 == KarlinReturn) {
      BlastErrorPrint(error_return);
    }
    sbp->kbp_gap_psi[0] = BlastKarlinBlkCreate();
    KarlinReturn = BlastKarlinBlkGappedCalc(sbp->kbp_gap_psi[0],11,1,sbp->name,&error_return);
    if (1 == KarlinReturn) {
      BlastErrorPrint(error_return);
    }

    array_size = BlastKarlinGetMatrixValues(sbp->name,&open,&extend,&lambda,&K,&H,NULL);
    if (array_size > 0) {
      for (index = 0; index < array_size; index++) {
        if (open[index] == INT2_MAX && extend[index] == INT2_MAX) {
          sbp->kbp_ideal = BlastKarlinBlkCreate();
          sbp->kbp_ideal->Lambda = lambda[index];
          sbp->kbp_ideal->K = K[index];
          sbp->kbp_ideal->H = H[index];
        }
      }
      MemFree(open);
      MemFree(extend);
      MemFree(lambda);
      MemFree(K);
      MemFree(H);
    }
    if (sbp->kbp_ideal == NULL) {
      sbp->kbp_ideal = BlastKarlinBlkStandardCalcEx(sbp);
    }
    compactSearch->lambda       = sbp->kbp_gap_std[0]->Lambda;
    compactSearch->kbp_std      = sbp->kbp_std;
    compactSearch->kbp_psi      = sbp->kbp_psi;
    compactSearch->kbp_gap_std  = sbp->kbp_gap_std;
    compactSearch->kbp_gap_psi  = sbp->kbp_gap_psi;
    compactSearch->lambda_ideal = sbp->kbp_ideal->Lambda;
    compactSearch->K_ideal      = sbp->kbp_ideal->K;

    stdrfp = BlastResFreqNew(sbp);
    BlastResFreqStdComp(sbp,stdrfp); 
    compactSearch->standardProb = MemNew(compactSearch->alphabetSize * sizeof(Nlm_FloatHi));
    for(a = 0; a < compactSearch->alphabetSize; a++)
      compactSearch->standardProb[a] = stdrfp->prob[a];
    stdrfp = BlastResFreqDestruct(stdrfp);
    posSearch->posInformation = NULL;
/*---------------------------------------------------------------------------*/
/* numseq is replaced with numalign (last argument) - each alignment is a    */
/* single segment and represents a separate sequence                         */
/*---------------------------------------------------------------------------*/
    CddposAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numalign);
    if (listOfSeqAligns->segtype == SAS_DENDIAG) {
      CddfindDenseDiagThreshSequences(posSearch, listOfSeqAligns, numalign, numalign);
    }
    else {
      CddfindThreshSequences(posSearch, listOfSeqAligns, numalign, numalign);
    }
    sbp->kbp     = sbp->kbp_psi;
    sbp->kbp_gap = sbp->kbp_gap_psi;
    CddposInitializeInformation(posSearch,sbp,compactSearch);
    if (listOfSeqAligns->segtype == SAS_DENDIAG) {
      CddposDenseDiagDemographics(posSearch, compactSearch, listOfSeqAligns);
    }
    else {
      CddposDemographics(posSearch, compactSearch, listOfSeqAligns);
    }
    CddposPurgeMatches(posSearch, compactSearch);
    CddposComputeExtents(posSearch, compactSearch);
    CddposComputeSequenceWeights(posSearch, compactSearch);
    CddposCheckWeights(posSearch, compactSearch);
    posSearch->posFreqs = (Nlm_FloatHi **) CddposComputePseudoFreqs(posSearch, compactSearch, TRUE);
    CddposFreqsToMatrix(posSearch,compactSearch);
    CddposScaling(posSearch, compactSearch);

    sctp = SeqCodeTableFind(Seq_code_ncbistdaa);
    for (a=0;a<compactSearch->alphabetSize;a++) {
      newLetter = ValNodeNew(NULL);
      newLetter->data.ptrvalue = MemNew(2*sizeof(Char));
      Nlm_StrNCpy(newLetter->data.ptrvalue,&(sctp->letters[a]),1);
      ValNodeLink(LetterHead,newLetter);
    }
}


Int4 LIBCALL CddGetNewIndexForThreading(char InChar, char* Output) {
/*---------------------------------------------------------------------------*/
/*  get the index in Output where InChar is found.                           */
/*  return -1 if char is not found in the Output array.                      */
/*---------------------------------------------------------------------------*/
    Int4  i, Index=0;

    /* look through Output for InChar */
    for (i=0; i<StrLen(Output); i++) {
      if (Output[i] == InChar) {
        return(i);
      }
    }
    return(-1);
}


Seq_Mtf* LIBCALL CddCalcPSSM(SeqAlignPtr listOfSeqAligns, BioseqPtr bsp) {
/*---------------------------------------------------------------------------*/
/* Generation of position-specific scoring matrices.                         */
/* Put results in a Seq_Mtf structure in preparation for calling atd().      */
/* allocate space for 2-d arrays pssm->ww and pssm->freq here.               */
/* pssm->ww[i][j] is pssm, where i is sequence index, j is residue index.    */
/* pssm->freqs[i][j] is table of residue frequencies.                        */
/*---------------------------------------------------------------------------*/
    posSearchItems      posSearch;           /*holds position-specific info  */
    compactSearchItems  compactSearch;
    ValNodePtr          LetterHead = NULL;
    Int4                i, j, jj, scale_factor;
    char*               Input;     /* order of residue-types from CD's */
    char*               Output;    /* order of residue-types needed for threading */
    Boolean*            Coverage;  /* for making sure all columns are filled */
    Int4                OutLen;
    Seq_Mtf*            pssm;

    /* allocate arrays for strings that determine input and output order */
    Input = (char*) MemNew(sizeof(InputOrder));
    Output = (char*) MemNew(sizeof(OutputOrder));
    StrCpy(Input, InputOrder);
    StrCpy(Output, OutputOrder);
    
    OutLen = StrLen(Output);
    CddCposCompPart1(listOfSeqAligns, bsp, &compactSearch, &LetterHead, &posSearch);

    /* construct new pssm */
    pssm = NewSeqMtf(compactSearch.qlength, OutLen);
    Coverage = (Boolean*) MemNew(sizeof(Boolean) * OutLen);  /* set to 0 by MemNew */
    
    /* copy results to pssm */
    scale_factor = 10000;
    for (j=0; j<compactSearch.alphabetSize;j++) {
      jj = CddGetNewIndexForThreading(Input[j], Output);
      if (jj != -1) {
        for (i=0; i<compactSearch.qlength;i++) {
          pssm->ww[i][jj] = (Int4) posSearch.posMatrix[i][j];
          pssm->freqs[i][jj] = (Int4) (0.5 + (Nlm_FloatHi) scale_factor * posSearch.posFreqs[i][j]);
          Coverage[jj] = TRUE;
        }
      }
    }

    /* make sure all columns are filled */
    for (i=0; i<OutLen; i++) {
      if (Coverage[i] == FALSE) {
        ASSERT(FALSE);  /* should never get here */
      }
    }
    return(pssm);
}


void LIBCALL CddCposComp(SeqAlignPtr listOfSeqAligns, BioseqPtr bsp, CddPtr pcdd) {
/*---------------------------------------------------------------------------*/
/* Generation of position-specific scoring matrices for database scanning    */
/*---------------------------------------------------------------------------*/
    posSearchItems      posSearch;           /*holds position-specific info  */
    compactSearchItems  compactSearch;
    Int4                index, a, iNew;
    ValNodePtr          error_return = NULL;
    ValNodePtr          ColumnHead;
    ValNodePtr          newRow, RowHead;
    ValNodePtr          LetterHead = NULL;
    Char                ckptFileName[PATH_MAX];
    Char                cseqFileName[PATH_MAX];
    FILE                *fp;

    CddCposCompPart1(listOfSeqAligns, bsp, &compactSearch, &LetterHead, &posSearch);

    /* put results in proper format... */

    pcdd->posfreq = (MatrixPtr) MatrixNew();
    pcdd->posfreq->ncolumns = compactSearch.qlength;
    pcdd->posfreq->nrows = compactSearch.alphabetSize;
    pcdd->posfreq->scale_factor = 10000;
    pcdd->posfreq->row_labels = LetterHead;    
    ColumnHead = NULL;
    for (index = 0; index<compactSearch.qlength;index++) {
      RowHead = NULL;
      for (a=0;a<compactSearch.alphabetSize;a++) {
         newRow = ValNodeNew(NULL);
         iNew = (Int4) (0.5 + (Nlm_FloatHi) pcdd->posfreq->scale_factor * posSearch.posFreqs[index][a]);
         newRow->data.intvalue = iNew;
           ValNodeLink(&(ColumnHead),newRow);
      }
    }
    pcdd->posfreq->columns = ColumnHead;

    pcdd->scoremat = (MatrixPtr) MatrixNew();
    pcdd->scoremat->ncolumns = compactSearch.qlength;
    pcdd->scoremat->nrows = compactSearch.alphabetSize;
    pcdd->scoremat->scale_factor = 1;
    pcdd->scoremat->row_labels = LetterHead;    
    ColumnHead = NULL;
    for (index = 0; index<compactSearch.qlength;index++) {
      RowHead = NULL;
      for (a=0;a<compactSearch.alphabetSize;a++) {
         newRow = ValNodeNew(NULL);
         iNew = (Int4) posSearch.posMatrix[index][a];
         newRow->data.intvalue = iNew;
           ValNodeLink(&(ColumnHead),newRow);
      }
    }
    pcdd->scoremat->columns = ColumnHead;

/*---------------------------------------------------------------------------*/
/* Construct name for checkpoint file                                        */
/*---------------------------------------------------------------------------*/
    strcpy(ckptFileName,pcdd->name);
    strcat(ckptFileName,CKPTEXT);

    CddposTakeCheckpoint(&posSearch, &compactSearch, ckptFileName, &error_return);
    strcpy(cseqFileName,pcdd->name);
    strcat(cseqFileName,CSEQEXT);
  
    fp = FileOpen(cseqFileName, "w");
    if (NULL == fp) {
      BlastConstructErrorMessage("CddCposComp", "Could not open fasta file", 1, &error_return);
    }
    BioseqToFasta(bsp,fp,FALSE);
    FileClose(fp);
}
