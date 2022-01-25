/* $Id: cddutil.c,v 1.38 2001/03/07 20:27:37 bauer Exp $
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
* $Revision: 1.38 $
*
* File Description: CDD utility routines
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddutil.c,v $
* Revision 1.38  2001/03/07 20:27:37  bauer
* fixed bug in CddBioseqCopy
*
* Revision 1.37  2001/03/07 16:30:33  bauer
* fixed alignment reindexing bug
*
* Revision 1.36  2001/03/05 20:25:09  bauer
* fixed gap-penalty selection when changing scoring matrix
*
* Revision 1.35  2001/02/15 15:41:41  shoemake
* fixed small memory leak
*
* Revision 1.34  2001/02/14 17:11:03  bauer
* relaced calls to BioseqCopy with CddBioseqCopy
*
* Revision 1.33  2001/02/13 20:55:10  bauer
* fixed bug when comparing local ids
*
* Revision 1.32  2001/02/07 11:36:51  thiessen
* fix minor memory leak
*
* Revision 1.31  2001/02/06 22:55:35  bauer
* Scoring Matrix now a function parameter in CddDenDiagCposComp2
*
* Revision 1.30  2001/02/06 20:56:10  hurwitz
* get rid of problematic BioSeqNew
*
* Revision 1.29  2001/02/05 22:58:46  bauer
* added alignment reindexing function
*
* Revision 1.28  2001/02/01 17:15:15  bauer
* removed SeqIdDup in CddReadBlastOptions
*
* Revision 1.27  2001/01/26 15:06:40  lewisg
* use entrez2 to retrieve structures
*
* Revision 1.26  2001/01/24 03:08:09  bauer
* fixed memory leaks
*
* Revision 1.25  2001/01/22 21:25:06  hurwitz
* bioseq id shouldn't be freed
*
* Revision 1.24  2001/01/19 23:34:55  bauer
* fixed memory leaks
*
* Revision 1.23  2001/01/19 21:43:44  bauer
* removed dependency from NR for threading PSSM calculations
*
* Revision 1.22  2001/01/17 21:32:02  bauer
* changes to PSSM calculation
*
* Revision 1.21  2001/01/12 01:33:32  bauer
* removed CddGetBlastOptions
*
* Revision 1.20  2001/01/12 00:17:01  bauer
* added routines for information content calculation
*
* Revision 1.19  2001/01/11 23:48:24  bauer
* added check for consensus-Cd
*
* Revision 1.18  2000/12/29 00:52:51  hurwitz
* cleaning memory leaks
*
* Revision 1.17  2000/12/20 18:56:40  hurwitz
* new random num gen, more debug printing
*
* Revision 1.16  2000/12/01 19:36:57  hurwitz
* added scale factor to PSSM calcs
*
* Revision 1.15  2000/11/22 22:34:49  hurwitz
* changed pssm scale factor to match contact potential scale factor
*
* Revision 1.14  2000/11/14 22:08:44  hurwitz
* added weighting for pssm calculation
*
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
/*#include <accentr.h>*/
#include <lsqfetch.h>
/* #include <netentr.h> */
/* #include <www.h> */
/* #include <sys/resource.h> */
/*#include <accutils.h>*/
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
/* Check if alignment master is a consensus Sequence                         */
/*---------------------------------------------------------------------------*/
Boolean LIBCALL SeqAlignHasConsensus(SeqAlignPtr salp)
{
  DenseDiagPtr      ddp;
  SeqIdPtr          sip;
  ObjectIdPtr       oidp;

  if (!salp) return(FALSE);
  ASSERT(salp->segtype == SAS_DENDIAG);
  ddp = salp->segs; if (!ddp) return (FALSE);
  sip = ddp->id;    if (!sip) return (FALSE);
  if (sip->choice != SEQID_LOCAL) return (FALSE);
  oidp = sip->data.ptrvalue; if (!oidp) return(FALSE);
  if (StringCmp(oidp->str,"consensus")==0) return (TRUE);
  return(FALSE);
}

/*---------------------------------------------------------------------------*/
/* Check if Cdd has a consensus Sequence                                     */
/*---------------------------------------------------------------------------*/
Boolean LIBCALL CddHasConsensus(CddPtr pcdd)
{
  SeqIntPtr   sintp;
  SeqIdPtr    sip;
  ObjectIdPtr oidp;
  
  if (!pcdd) return (FALSE);
  sintp = (SeqIntPtr) pcdd->profile_range;
  if (!sintp) return (FALSE);
  sip = (SeqIdPtr) sintp->id;
  if (!sip) return (FALSE);
  if (sip->choice != SEQID_LOCAL) return (FALSE);
  oidp = sip->data.ptrvalue; if (!oidp) return(FALSE);
  if (StringCmp(oidp->str,"consensus")==0) return (TRUE);
  return(FALSE);
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
/* version of BioseqCopy that doesn't use BioseqFind and expects the old bsp */
/* adapted by Ben                                                            */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
BioseqPtr LIBCALL CddBioseqCopy (SeqIdPtr newid, BioseqPtr oldbsp, Int4 from,
                                 Int4 to, Uint1 strand, Boolean do_feat)
{
  BioseqPtr   newbsp=NULL, tmpbsp;
  Uint1       seqtype;
  ValNodePtr  tmp;
  ObjectIdPtr oid;
  Int4        len, i;
  Int2        residue;
  ValNode     fake;
  SeqLocPtr   the_segs, head, curr;
  Boolean     handled = FALSE, split;
  SeqFeatPtr  sfp, newsfp, lastsfp;
  DeltaSeqPtr dsp;
  SeqEntryPtr oldscope;


  if (from < 0) return FALSE;
  if (oldbsp == NULL) return NULL;
  len = to - from + 1;
  if (len <= 0) return NULL;
  newbsp = BioseqNew();
  if (newid != NULL) newbsp->id = SeqIdDup(newid);
  else {
    tmp = ValNodeNew(NULL);
    tmp->choice = SEQID_LOCAL;
    oid = ObjectIdNew();
    tmp->data.ptrvalue = (Pointer)oid;
    oid->str = StringSave("Clipboard");
    tmpbsp = BioseqFind(tmp);   /* old clipboard present? */
    if (tmpbsp == NULL) {
      oldscope = SeqEntrySetScope (NULL);
      if (oldscope != NULL) {
        tmpbsp = BioseqFind(tmp);
        SeqEntrySetScope (oldscope);
      }
    }
    if (tmpbsp != NULL) BioseqFree(tmpbsp);
    newbsp->id = tmp;
  }

  newbsp->repr = oldbsp->repr;
  newbsp->mol = oldbsp->mol;
  newbsp->length = len;
  newbsp->seq_ext_type = oldbsp->seq_ext_type;

  if (newbsp->repr == Seq_repr_virtual)
    handled = TRUE;               /* no more to do */

  if ((newbsp->repr == Seq_repr_raw) ||
    (newbsp->repr == Seq_repr_const)) {
      AsnTypePtr  atp;
      AsnIoPtr  aip;
    if (ISA_aa(newbsp->mol)) {
      seqtype = Seq_code_ncbieaa;
    } else {
      seqtype = Seq_code_iupacna;
    }
/*---------------------------------------------------------------------------*/
/* need to check whether seqtypes agree, otherwise the ByteStore is invalid! */
/*---------------------------------------------------------------------------*/
    if (seqtype != oldbsp->seq_data_type) {
      seqtype = oldbsp->seq_data_type;
    }
    newbsp->seq_data_type = seqtype;
    newbsp->seq_data = BSDup(oldbsp->seq_data);
    handled = TRUE;
  }

  if ((newbsp->repr == Seq_repr_seg) ||
    (newbsp->repr == Seq_repr_ref) ||
    (newbsp->repr == Seq_repr_delta)) {
    if (newbsp->repr == Seq_repr_seg) { /* segmented */
      fake.choice = SEQLOC_MIX;   /* make SEQUENCE OF Seq-loc, into one */
      fake.data.ptrvalue = oldbsp->seq_ext;
      fake.next = NULL;
      the_segs = (SeqLocPtr)&fake;
      head = (SeqLocPtr)SeqLocCopyPart (the_segs, from, to, strand, FALSE, NULL, NULL);
    } else if (newbsp->repr == Seq_repr_ref) { /* reference: is a Seq-loc */
      head = (SeqLocPtr)SeqLocCopyPart ((SeqLocPtr)(oldbsp->seq_ext), from, to,
                               strand, TRUE, NULL, NULL);
    } else if (newbsp->repr == Seq_repr_delta) {
      dsp = (DeltaSeqPtr)(oldbsp->seq_ext);  /* real data is here */
      the_segs = (SeqLocPtr)DeltaSeqsToSeqLocs(dsp);
      head = (SeqLocPtr)SeqLocCopyPart (the_segs, from, to, strand, FALSE, NULL, NULL);
      SeqLocFree (the_segs);
    }
    newbsp->seq_ext = (Pointer)head;
    handled = TRUE;
  }

  if (newbsp->repr == Seq_repr_map) {
    lastsfp = NULL;
    for (sfp = (SeqFeatPtr)(oldbsp->seq_ext); sfp != NULL; sfp = sfp->next) {
      split = FALSE;
      curr = (SeqLocPtr)SeqLocCopyRegion(newbsp->id, sfp->location, oldbsp, from, to, strand, &split);
      if (curr != NULL) {  /* got one */
        newsfp = (SeqFeatPtr)AsnIoMemCopy((Pointer)sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
        SeqLocFree(newsfp->location);
        newsfp->location = curr;
        if (split)
          newsfp->partial = TRUE;
        if (lastsfp == NULL)  /* first one */
          newbsp->seq_ext = (Pointer)newsfp;
        else
          lastsfp->next = newsfp;
        lastsfp = newsfp;
      }
    }
    handled = TRUE;
  }
  if (! handled) goto erret;
  if (do_feat)
    SeqFeatsCopy (newbsp, oldbsp, from, to, strand);

  return newbsp;

erret:
  BioseqFree(newbsp);
  return NULL;
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
  bspNew = (BioseqPtr) CddBioseqCopy(sipNew,bsp,0,bsp->length-1,0,FALSE);
  sipNew = SeqIdFree(sipNew);
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
/*---------------------------------------------------------------------------*/
/* Create a Dense-Diag seqalign from a Dense-Seg SeqAlign                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SeqAlignPtr LIBCALL CddMSLDenSegToMSLDenDiag(SeqAlignPtr salp)
{
  SeqAlignPtr         salpnew, salphead, salptail;
  DenseDiagPtr        ddp, ddphead, ddptail;
  DenseSegPtr         dsp;
  Int4                i, s1, s2;


  if (!salp) return(NULL);
  salphead = NULL;
  salptail = NULL;
  while (salp) {
    if (salp->segtype) {
      if (salp->segtype != SAS_DENSEG) {
        CddSevError("CddMSLDenSegToMSLDenDiag: Wrong segtype specified!");
      }
    }
    if (salp->dim) {
      if (salp->dim != 2) {
        CddSevError("CddMSLDenSegToMSLDenDiag: Expect alignments of dimension 2!");
      }
    }
    dsp = salp->segs;
    ddphead = NULL;
    ddptail = NULL;
    for (i=0;i<dsp->numseg;i++) {
      s1 = dsp->starts[i*2];
      s2 = dsp->starts[i*2+1]; 
      if (s1 >=0 && s2>= 0) {
        ddp = DenseDiagNew();
	ddp->starts = MemNew(2*sizeof(Int4));
	ddp->starts[0] = s1;
	ddp->starts[1] = s2;
	ddp->len=dsp->lens[i];
        ddp->id = dsp->ids;
	ddp->dim = 2;
        if (!ddphead) {
	  ddphead = ddp;
	  ddptail = ddp;
	  ddptail->next = NULL;
	} else {
	  ddptail->next = ddp;
	  ddptail = ddp;
	  ddptail->next = NULL;
	}
      }
    }
    salpnew = SeqAlignNew();
    salpnew->type = SAT_PARTIAL;
    salpnew->segtype = SAS_DENDIAG;
    salpnew->dim = 2;
    salpnew->segs = ddphead;
    if (!salphead) {
      salphead = salpnew;
      salptail = salpnew;
      salptail->next = NULL;
    } else {
      salptail->next = salpnew;
      salptail = salpnew;
      salptail->next = NULL;
    }
    salp = salp->next;
  }
  return(salphead);
}

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
/* Calculate per column information content for a SeqAlign                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Nlm_FloatHiPtr LIBCALL SeqAlignInform(SeqAlignPtr salp, BioseqPtr bsp_master,
                                      Boolean bHasConsensus)
{
  BLAST_ScoreBlkPtr  sbp;
  BLAST_ResFreqPtr   stdrfp;
  Nlm_FloatHiPtr     standardProb;
  Nlm_FloatHiPtr     Informativeness;
  Int2               iStatus;
  Int4               i, c, a;
  Int4               qlength;
  Int4               alphabetSize = 26;
  MatrixPtr          posfreq;
  ValNodePtr         vnp;
  Int4               nColumns;
  Nlm_FloatHi        posEpsilon = 0.0001;
  Nlm_FloatHi        QoverPest;
  Int4               *ntypes;
  Nlm_FloatHi        **typefreq, sumfreq;
  DenseDiagPtr       ddp;
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  SeqPortPtr         spp;
  Uint1Ptr           buffer;
  
  if (!salp) return(NULL);
  if (!bsp_master) return(NULL);
  qlength = bsp_master->length;
  Informativeness = MemNew(sizeof(Nlm_FloatHi) * qlength);
  sbp = BLAST_ScoreBlkNew(Seq_code_ncbistdaa,1);
  sbp->read_in_matrix = TRUE;
  sbp->protein_alphabet = TRUE;
  sbp->posMatrix = NULL;
  sbp->number_of_contexts = 1;
  iStatus = BlastScoreBlkMatFill(sbp,"BLOSUM62");
  stdrfp = BlastResFreqNew(sbp);
  BlastResFreqStdComp(sbp,stdrfp); 
  standardProb = MemNew(alphabetSize * sizeof(Nlm_FloatHi));
  for(a = 0; a < alphabetSize; a++) standardProb[a] = stdrfp->prob[a];
  stdrfp = BlastResFreqDestruct(stdrfp);
  ntypes = (Int4 *) MemNew(qlength*sizeof(Int4));
  for (i=0;i<qlength;i++) ntypes[i] = 0;
  typefreq = (Nlm_FloatHi**) MemNew(qlength * sizeof(Nlm_FloatHi *));
  for (i=0;i<qlength;i++) {
    typefreq[i] = (Nlm_FloatHi *)MemNew(alphabetSize * sizeof(Nlm_FloatHi));
    for (a=0;a<alphabetSize;a++) typefreq[i][a] = 0.0;
  }
  if (!bHasConsensus) {   /* count residues in the master/representative too */
    ddp = salp->segs;
    sip = ddp->id;
    bsp = bsp_master;
    buffer = MemNew((bsp->length)*sizeof(Uint1));
    spp = SeqPortNew(bsp, 0, bsp->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
    for (i=0; i<bsp->length;i++) buffer[i] = SeqPortGetResidue(spp);
    spp = SeqPortFree(spp);
    while (ddp) {
      for (c=ddp->starts[0];c<ddp->starts[0]+ddp->len;c++) {
        if (buffer[c] >= 0 && buffer[c] < alphabetSize)
          typefreq[c][buffer[c]] += 1.0;
      }
      ddp = ddp->next;
    }
    MemFree(buffer);
  }  
  while (salp) {
    ddp = salp->segs;
    sip = ddp->id->next;
    bsp = CddRetrieveBioseqById(sip,NULL);
    buffer = MemNew((bsp->length)*sizeof(Uint1));
    spp = SeqPortNew(bsp, 0, bsp->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
    for (i=0; i<bsp->length;i++) buffer[i] = SeqPortGetResidue(spp);
    spp = SeqPortFree(spp);
    while (ddp) {
      for (c=ddp->starts[1];c<ddp->starts[1]+ddp->len;c++) {
        i = ddp->starts[0]+c-ddp->starts[1];
        if (buffer[c] >= 0 && buffer[c] < alphabetSize)
	  typefreq[i][buffer[c]] += 1.0;
      }
      ddp = ddp->next;
    }
    MemFree(buffer);
    salp = salp->next;
  }
  for (i=0;i<qlength;i++) {
    sumfreq = 0.0;
    for (a=0;a<alphabetSize;a++) {
      if (typefreq[i][a] > 0.0) {
        ntypes[i]++;
	sumfreq += typefreq[i][a];
      }
    }
    if (sumfreq > 0.0) for (a=0;a<alphabetSize;a++) {
      typefreq[i][a] /= sumfreq; 
    }
  }

  for (c=0;c<qlength;c++) {
    Informativeness[c] = 0.0;
    for (a=0;a<alphabetSize;a++) {
      if (standardProb[a] > posEpsilon) {
        QoverPest = typefreq[c][a] / standardProb[a];
        if (QoverPest > posEpsilon) {
	  Informativeness[c] += typefreq[c][a] * log(QoverPest) / NCBIMATH_LN2;        
        }
      }
    }
  }
  for (i=0;i<qlength;i++) MemFree(typefreq[i]);
  MemFree(typefreq);
  MemFree(ntypes);
  MemFree(standardProb);
  sbp = (BLAST_ScoreBlkPtr) BLAST_ScoreBlkDestruct(sbp);
  return(Informativeness);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Calculate per column information content for a CDD from the seqalign      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Nlm_FloatHiPtr LIBCALL CddAlignInform(CddPtr pcdd, Nlm_FloatHi * Niobs)
{
  Int4               offset;
  Boolean            bHasConsensus = CddHasConsensus(pcdd);
  BLAST_ScoreBlkPtr  sbp;
  BLAST_ResFreqPtr   stdrfp;
  Nlm_FloatHiPtr     standardProb;
  Nlm_FloatHiPtr     Informativeness;
  Int2               iStatus;
  Int4               i, c, a;
  Int4               qlength;
  Int4               alphabetSize = 26;
  MatrixPtr          posfreq;
  ValNodePtr         vnp;
  Int4               nColumns;
  Nlm_FloatHi        posEpsilon = 0.0001;
  Nlm_FloatHi        QoverPest;
  Int4               *ntypes;
  Nlm_FloatHi        **typefreq, sumfreq;
  SeqAlignPtr        salp;
  DenseDiagPtr       ddp;
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  BioseqSetPtr       bssp;
  SeqPortPtr         spp;
  Uint1Ptr           buffer;
  
  bssp = pcdd->sequences->data.ptrvalue;
  offset = pcdd->profile_range->from;
  if (!pcdd) return(NULL);
  if (!pcdd->trunc_master) return(NULL);
  qlength = pcdd->trunc_master->length;
  Informativeness = MemNew(sizeof(Nlm_FloatHi) * qlength);
  sbp = BLAST_ScoreBlkNew(Seq_code_ncbistdaa,1);
  sbp->read_in_matrix = TRUE;
  sbp->protein_alphabet = TRUE;
  sbp->posMatrix = NULL;
  sbp->number_of_contexts = 1;
  iStatus = BlastScoreBlkMatFill(sbp,"BLOSUM62");
  stdrfp = BlastResFreqNew(sbp);
  BlastResFreqStdComp(sbp,stdrfp); 
  standardProb = MemNew(alphabetSize * sizeof(Nlm_FloatHi));
  for(a = 0; a < alphabetSize; a++) standardProb[a] = stdrfp->prob[a];
  stdrfp = BlastResFreqDestruct(stdrfp);
  ntypes = (Int4 *) MemNew(qlength*sizeof(Int4));
  for (i=0;i<qlength;i++) ntypes[i] = 0;
  typefreq = (Nlm_FloatHi**) MemNew(qlength * sizeof(Nlm_FloatHi *));
  for (i=0;i<qlength;i++) {
    typefreq[i] = (Nlm_FloatHi *)MemNew(alphabetSize * sizeof(Nlm_FloatHi));
    for (a=0;a<alphabetSize;a++) typefreq[i][a] = 0.0;
  }
  salp = pcdd->seqannot->data;
  if (!bHasConsensus) {   /* count residues in the master/representative too */
    ddp = salp->segs;
    sip = ddp->id;
    bsp = CddRetrieveBioseqById(sip,bssp->seq_set);
    buffer = MemNew((bsp->length)*sizeof(Uint1));
    spp = SeqPortNew(bsp, 0, bsp->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
    for (i=0; i<bsp->length;i++) buffer[i] = SeqPortGetResidue(spp);
    spp = SeqPortFree(spp);
    while (ddp) {
      for (c=ddp->starts[0];c<ddp->starts[0]+ddp->len;c++) {
        typefreq[c-offset][buffer[c]] += 1.0;
      }
      ddp = ddp->next;
    }
    MemFree(buffer);
  }  
  while (salp) {
    ddp = salp->segs;
    sip = ddp->id->next;
    bsp = CddRetrieveBioseqById(sip,bssp->seq_set);
    buffer = MemNew((bsp->length)*sizeof(Uint1));
    spp = SeqPortNew(bsp, 0, bsp->length-1, Seq_strand_unknown, Seq_code_ncbistdaa);
    for (i=0; i<bsp->length;i++) buffer[i] = SeqPortGetResidue(spp);
    spp = SeqPortFree(spp);
    while (ddp) {
      for (c=ddp->starts[1];c<ddp->starts[1]+ddp->len;c++) {
        i = ddp->starts[0]-offset+c-ddp->starts[1];
        typefreq[i][buffer[c]] += 1.0;
      }
      ddp = ddp->next;
    }
    MemFree(buffer);
    salp = salp->next;
  }
  for (i=0;i<qlength;i++) {
    sumfreq = 0.0;
    for (a=0;a<alphabetSize;a++) {
      if (typefreq[i][a] > 0.0) {
        ntypes[i]++;
	sumfreq += typefreq[i][a];
      }
    }
    if (sumfreq > 0.0) for (a=0;a<alphabetSize;a++) {
      typefreq[i][a] /= sumfreq; 
    }
  }

  for (c=0;c<qlength;c++) {
    Informativeness[c] = 0.0;
    for (a=0;a<alphabetSize;a++) {
      if (standardProb[a] > posEpsilon) {
        QoverPest = typefreq[c][a] / standardProb[a];
        if (QoverPest > posEpsilon) {
	  Informativeness[c] += typefreq[c][a] * log(QoverPest) / NCBIMATH_LN2;        
        }
      }
    }
  }

  *Niobs = 0.0;
  for (i=0;i<qlength;i++) *Niobs += (Nlm_FloatHi) ntypes[i];
  *Niobs /= (Nlm_FloatHi) qlength;

  for (i=0;i<qlength;i++) MemFree(typefreq[i]);
  MemFree(typefreq);
  MemFree(ntypes);
  MemFree(standardProb);
  return(Informativeness);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Calculate per column information content for a CDD pssm                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Nlm_FloatHiPtr LIBCALL CddPssmInform(CddPtr pcdd)
{
  BLAST_ScoreBlkPtr  sbp;
  BLAST_ResFreqPtr   stdrfp;
  Nlm_FloatHiPtr     standardProb;
  Nlm_FloatHiPtr     Informativeness;
  Nlm_FloatHi        thisposFreq;
  Int2               iStatus;
  Int4               i, c, a;
  Int4               qlength;
  Int4               alphabetSize = 26;
  MatrixPtr          posfreq;
  ValNodePtr         vnp;
  Int4               nColumns;
  Nlm_FloatHi        posEpsilon = 0.0001;
  Nlm_FloatHi        QoverPest;

  if (!pcdd) return(NULL);
  if (!pcdd->trunc_master) return(NULL);
  qlength = pcdd->trunc_master->length;
  Informativeness = MemNew(sizeof(Nlm_FloatHi) * qlength);
  sbp = BLAST_ScoreBlkNew(Seq_code_ncbistdaa,1);
  sbp->read_in_matrix = TRUE;
  sbp->protein_alphabet = TRUE;
  sbp->posMatrix = NULL;
  sbp->number_of_contexts = 1;
  iStatus = BlastScoreBlkMatFill(sbp,"BLOSUM62");
  stdrfp = BlastResFreqNew(sbp);
  BlastResFreqStdComp(sbp,stdrfp); 
  standardProb = MemNew(alphabetSize * sizeof(Nlm_FloatHi));
  for(a = 0; a < alphabetSize; a++) standardProb[a] = stdrfp->prob[a];
  stdrfp = BlastResFreqDestruct(stdrfp);
  posfreq = pcdd->posfreq;
  if (posfreq->nrows != alphabetSize || posfreq->ncolumns != qlength) {
    CddSevError("posfreq matrix size mismatch in CddPssmInform");  
  }
  vnp = posfreq->columns;
  for (c=0;c<qlength;c++) {
    Informativeness[c] = 0.0;
    for (a=0;a<alphabetSize;a++) {
      if (standardProb[a] > posEpsilon) {
        thisposFreq = (Nlm_FloatHi) vnp->data.intvalue / (Nlm_FloatHi) posfreq->scale_factor;
        QoverPest = thisposFreq / standardProb[a];
        if (QoverPest > posEpsilon) {
	  Informativeness[c] += thisposFreq * log(QoverPest) / NCBIMATH_LN2;        
        }
      }
      vnp = vnp->next;
    }
  }
  MemFree(standardProb);
  return(Informativeness);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Calculate per column information content for a frequency table            */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Nlm_FloatHiPtr LIBCALL CddPosFreqInform(Nlm_FloatHi **posFreq, Int4 ncol, Int4 nrow)
{
  BLAST_ScoreBlkPtr  sbp;
  BLAST_ResFreqPtr   stdrfp;
  Nlm_FloatHiPtr     standardProb;
  Nlm_FloatHiPtr     Informativeness;
  Nlm_FloatHi        thisposFreq;
  Int2               iStatus;
  Int4               i, c, a;
  Int4               qlength;
  Int4               alphabetSize = 26;
  Int4               nColumns;
  Nlm_FloatHi        posEpsilon = 0.0001;
  Nlm_FloatHi        QoverPest;

  if (!posFreq) return(NULL);
  qlength = ncol;
  if (nrow != alphabetSize) return(NULL);
  Informativeness = MemNew(sizeof(Nlm_FloatHi) * qlength);
  sbp = BLAST_ScoreBlkNew(Seq_code_ncbistdaa,1);
  sbp->read_in_matrix = TRUE;
  sbp->protein_alphabet = TRUE;
  sbp->posMatrix = NULL;
  sbp->number_of_contexts = 1;
  iStatus = BlastScoreBlkMatFill(sbp,"BLOSUM62");
  stdrfp = BlastResFreqNew(sbp);
  BlastResFreqStdComp(sbp,stdrfp); 
  standardProb = MemNew(alphabetSize * sizeof(Nlm_FloatHi));
  for(a = 0; a < alphabetSize; a++) standardProb[a] = stdrfp->prob[a];
  stdrfp = BlastResFreqDestruct(stdrfp);

  for (c=0;c<qlength;c++) {
    Informativeness[c] = 0.0;
    for (a=0;a<alphabetSize;a++) {
      if (standardProb[a] > posEpsilon) {
        thisposFreq = posFreq[c][a];
        QoverPest = thisposFreq / standardProb[a];
        if (QoverPest > posEpsilon) {
	  Informativeness[c] += thisposFreq * log(QoverPest) / NCBIMATH_LN2;        
        }
      }
    }
  }
  MemFree(standardProb);
  return(Informativeness);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Generation of position-specific scoring matrices for database scanning    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void LIBCALL CddDenDiagCposComputation(SeqAlignPtr listOfSeqAligns, BioseqPtr bsp,
                                       BioseqPtr bspF, CddPtr pcdd, Int4 pseudoCnt)
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
    Boolean             bHasConsensus;
    
    bHasConsensus = CddHasConsensus(pcdd);
    
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
    compactSearch->pseudoCountConst = pseudoCnt;
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
    posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numalign);
    CddfindDenseDiagThreshSequences(posSearch, listOfSeqAligns, numalign, numalign);
    sbp->kbp     = sbp->kbp_psi;
    sbp->kbp_gap = sbp->kbp_gap_psi;
    CddposInitializeInformation(posSearch,sbp,compactSearch);
    CddposDenseDiagDemographics(posSearch, compactSearch, listOfSeqAligns);
    posPurgeMatches(posSearch, compactSearch);
/*---------------------------------------------------------------------------*/
/* remove the consensus from the CposComputation if present in the CDD       */
/*---------------------------------------------------------------------------*/
    if (bHasConsensus) posCancel(posSearch,compactSearch,0,0,0,compactSearch->qlength);
    posComputeExtents(posSearch, compactSearch);
    posComputeSequenceWeights(posSearch, compactSearch, 1.0);
    posCheckWeights(posSearch, compactSearch);
    posSearch->posFreqs = (Nlm_FloatHi **) posComputePseudoFreqs(posSearch, compactSearch, TRUE);
    CddposFreqsToMatrix(posSearch,compactSearch);
    posScaling(posSearch, compactSearch);

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
    posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numalign);
    CddfindThreshSequences(posSearch, listOfSeqAligns, numalign, numalign);
    sbp->kbp     = sbp->kbp_psi;
    sbp->kbp_gap = sbp->kbp_gap_psi;
    CddposInitializeInformation(posSearch,sbp,compactSearch);
    CddposDemographics(posSearch, compactSearch, listOfSeqAligns);
    posPurgeMatches(posSearch, compactSearch);
    posComputeExtents(posSearch, compactSearch);
    posComputeSequenceWeights(posSearch, compactSearch, 1.0);
    posCheckWeights(posSearch, compactSearch);
    posSearch->posFreqs = (Nlm_FloatHi **) posComputePseudoFreqs(posSearch, compactSearch, TRUE);
    CddposFreqsToMatrix(posSearch,compactSearch);
    posScaling(posSearch, compactSearch);

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
    if(seq_id->choice == SEQID_PDB) {
      pdb_seq_id = seq_id->data.ptrvalue;
      break;
    }
    seq_id = seq_id->next;
  }
  return(pdb_seq_id);
}


/*---------------------------------------------------------------------------*/
/* test whether two SeqIds are the same, only considers PDB-Ids and GIs      */
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
        if (tsip1->choice == SEQID_PDB) {
          psip1 = (PDBSeqIdPtr) CddGetPdbSeqId(tsip1);
          psip2 = (PDBSeqIdPtr) CddGetPdbSeqId(tsip2);
          if (strcmp(psip1->mol,psip2->mol)==0 && psip1->chain==psip2->chain)
            return(TRUE);
        } else if (tsip1->choice == SEQID_LOCAL) {
          oidp1 = tsip1->data.ptrvalue;
          oidp2 = tsip2->data.ptrvalue;
          if (NULL != oidp1 && strcmp(oidp1->str,oidp2->str)==0) return(TRUE);
        } else if (tsip1->choice == SEQID_GI) {
	  if (tsip1->data.intvalue == tsip2->data.intvalue) return(TRUE);
	}
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
  pCDea->bIdAlloc = FALSE;
  return(pCDea);
}

CddExpAlignPtr CddExpAlignFree(CddExpAlignPtr pCDea)
{
  SeqIdPtr sip, sipNext;

  if (!pCDea) return NULL;
  if (NULL != pCDea->adata) free(pCDea->adata);
  if (pCDea->bIdAlloc) {
    sip = pCDea->ids;
    while (sip) {
      sipNext = sip->next;
      SeqIdFree(sip);
      sip = sipNext;
    }
  }
  free(pCDea);
  return(pCDea);
}

void           CddExpAlignAlloc(CddExpAlignPtr pCDea, Int4 iLength)
{
  Int4    i;
  
  if (!pCDea || !iLength) return;
  pCDea->length = iLength;
  if (NULL != pCDea->adata) free(pCDea->adata);
  pCDea->adata = MemNew(sizeof(Int4)*pCDea->length);  
  for (i=0;i<pCDea->length;i++) pCDea->adata[i] = -1;
}

/*---------------------------------------------------------------------------*/
/* Convert a SeqAlign (pairwise, DenseDiag) into an ExplicitAlignmentObject  */
/*---------------------------------------------------------------------------*/
CddExpAlignPtr SeqAlignToCddExpAlign(SeqAlignPtr salp, SeqEntryPtr sep)
{
  CddExpAlignPtr   pCDea;
  BioseqPtr        bsp1;
  DenseDiagPtr     ddp;
  Int4             i;
  SeqIdPtr         sip;
  
  pCDea = CddExpAlignNew();
  ddp = salp->segs;
  pCDea->ids = ddp->id;
  sip = SeqIdDup(pCDea->ids);
  bsp1 = CddRetrieveBioseqById(sip,sep);
  SeqIdFree(sip);
  CddExpAlignAlloc(pCDea,bsp1->length);
  while (ddp) {
    for (i=0;i<ddp->len;i++) {
      if (ddp->starts[0]+i >= bsp1->length) {
        printf("Warning: Indexing error!\n");
      }
      pCDea->adata[ddp->starts[0]+i]=ddp->starts[1]+i;
    }
    ddp = ddp->next;
  }
  return(pCDea);
}

/*---------------------------------------------------------------------------*/
/* Invert an explicit alignment object (flip master-slave)                   */
/*---------------------------------------------------------------------------*/
CddExpAlignPtr InvertCddExpAlign(CddExpAlignPtr pCDea, SeqEntryPtr sep)
{
  BioseqPtr      bsp2;
  CddExpAlignPtr pCDeaNew;
  Int4           i;
 
  if (!pCDea->ids || !pCDea->ids->next) CddSevError("No SeqId in InvertCddExpAlign");
  pCDeaNew = CddExpAlignNew();
  pCDeaNew->ids = SeqIdDup(pCDea->ids->next);
  pCDeaNew->ids->next = SeqIdDup(pCDea->ids);
  pCDeaNew->bIdAlloc = TRUE;
  bsp2 = CddRetrieveBioseqById(pCDea->ids->next, sep);
  CddExpAlignAlloc(pCDeaNew, bsp2->length);
  for (i=0;i<pCDea->length;i++) {
    if (pCDea->adata[i] != -1) {
      pCDeaNew->adata[pCDea->adata[i]] = i;
    }
  }
  return(pCDeaNew);
}

/*---------------------------------------------------------------------------*/
/* Convert an Explicit Alignment Object (with proper IDs) to a SeqAlign      */
/*---------------------------------------------------------------------------*/
SeqAlignPtr    CddExpAlignToSeqAlign(CddExpAlignPtr pCDea)
{
  SeqAlignPtr   salp;
  DenseDiagPtr  ddp, ddplast = NULL;
  Int4          i, last;
  
  if (!pCDea->ids || !pCDea->ids->next)
   CddSevError("Missing Sequence ID in CddExpAlignToSeqAlign");
  salp = SeqAlignNew();
  salp->type = SAT_PARTIAL;
  salp->segtype = SAS_DENDIAG;
  salp->dim = 2;
  salp->segs = NULL;
  
  i= -1;
  while (++i < pCDea->length) {
    if (pCDea->adata[i] != -1) {
      ddp = DenseDiagNew();
      ddp->id = SeqIdDup(pCDea->ids);
      ddp->id->next = SeqIdDup(pCDea->ids->next);
      ddp->dim = 2;
      ddp->starts = MemNew(2*sizeof(Int4));
      ddp->starts[0] = i;
      ddp->starts[1] = pCDea->adata[i];
      ddp->len = 1;
      while (i<(pCDea->length-1) &&  pCDea->adata[i+1]==(pCDea->adata[i]+1)) {
	ddp->len++;
	i++;
      }
      if (ddplast) {
        ddplast->next = ddp;
      } else {
        salp->segs = ddp;
      }
      ddplast = ddp;
    }
  }
  return(salp);
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
  SeqIdPtr           sip1, sip2;
  
  if (!pCDea1 || !pCDea2) return pCDea;
  if (pCDea1->ids && pCDea1->ids->next) {
    sip1 = SeqIdDup(pCDea1->ids->next);  
  } else sip1 = NULL;
  if (pCDea2->ids && pCDea2->ids->next) {
    sip2 = SeqIdDup(pCDea2->ids->next);  
  } else sip2 = NULL;
  
  pCDea = CddExpAlignNew();
  CddExpAlignAlloc(pCDea,newlength);
  for (i=0;i<pCDea1->length && i<pCDea2->length;i++) {
    iP1 = pCDea1->adata[i]; iP2 = pCDea2->adata[i];
    if (iP1 > -1 && iP2 > -1){
      if (iP1 >= newlength) {
        printf("Warning: index error between sequences %d and %d\n",
	        iInner, iOuter);
        CddSevError("Error while reindexing explicit alignments!");
      } else  pCDea->adata[iP1] = iP2;
    }
  }
  if (sip1 && sip2) {
    pCDea->ids = sip1;
    pCDea->ids->next = sip2;
    pCDea->bIdAlloc = TRUE;
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
				 BioseqPtr   bsp,
				 SeqIntPtr   range,
				 BioseqPtr   *bspCons,
				 SeqAlignPtr *salpCons)
{
  BioseqPtr           bspThis;
  CddAlignWeightPtr   pCAW;
  CddExpAlignPtr      pCDea;
  CddExtAlignCellPtr  *pCEACell;
  CddExtAlignColPtr   pCEACol;
  CharPtr             cbuffer, outptr;
  DenseDiagPtr        ddp;
  Int4Ptr             trunc_on_virtual;
  Int4Ptr             trunc_aligned;
  Int4Ptr             virtual_on_trunc;
  ObjectIdPtr         oidp;
  SeqAlignPtr         salpThis, salpNew, salpTail;
  SeqEntryPtr         sepNew;
  SeqIdPtr            sipThis;
  SeqMapTablePtr      smtp;
  SeqPortPtr          spp;
  Uint1Ptr            buffer;
  FloatHi             SumWeight;
  Int4                maxextlen, nalign, ConsensusLen;
  Int4                index, midpt, ContWeight = 0;
  Int4                offset, i, j , k, ipt, seqpos;
  Int4                m_last, m_frst, s_insert, m_end;
  Int4                s_last, s_frst, maxpos = 0;
  Uint1               aatype;
  
  if (!salp) return (NULL);
  if (salp->type != SAT_PARTIAL || salp->segtype != SAS_DENDIAG) return(NULL);
  if (!range) return (NULL);
  offset = range->from;
/*---------------------------------------------------------------------------*/
/* trunc_on_virtual holds, for every residue of the trunc_master, its posi-  */
/* tion on the maximal expanded new master sequence                          */
/* bsp is supposed to be the "truncated master", i.e. a bioseq covering the  */
/* interval from the first to the last aligned residue on the master sequence*/
/*---------------------------------------------------------------------------*/
  trunc_on_virtual = MemNew(sizeof(Int4)*bsp->length);
  for (i=0;i<bsp->length;i++) trunc_on_virtual[i] = i;
  maxextlen = bsp->length-1;
  trunc_aligned = MemNew(sizeof(Int4)*bsp->length);
  for (i=0;i<bsp->length;i++) trunc_aligned[i] = 0;

/*---------------------------------------------------------------------------*/
/* expand values on trunc_on_virtual to consider insertions relativ to master*/
/*---------------------------------------------------------------------------*/
  nalign = 1;
  salpThis = salp;
  while (salpThis) {
    nalign++;
    ddp = salpThis->segs;
    if (ddp) for (j=ddp->starts[0]-offset;j<ddp->starts[0]+ddp->len-offset;j++) trunc_aligned[j]=1;
    while (ddp && ddp->next) {
      for (j=ddp->next->starts[0]-offset;j<ddp->next->starts[0]+ddp->next->len-offset;j++) trunc_aligned[j]=1;
      m_last = ddp->starts[0]+ddp->len - 1 - offset;
      m_frst = ddp->next->starts[0] - offset;
      m_end  = ddp->next->starts[0] + ddp->next->len - 1 - offset;
      ASSERT(m_frst >= 0 && m_frst < bsp->length);
      ASSERT(m_last >= 0 && m_last < bsp->length);
      ASSERT(m_end >= 0 && m_end < bsp->length);
/*---------------------------------------------------------------------------*/
/* shift block from the c- to n-terminus if necessary                        */
/*---------------------------------------------------------------------------*/
/*    for (j = m_end; j > m_frst; j--) {
        if (trunc_on_virtual[j-1] < trunc_on_virtual[j]-1) {
	  trunc_on_virtual[j-1] = trunc_on_virtual[j] - 1;
	}
      }
*/
      s_insert = ddp->next->starts[1] - ddp->starts[1] - ddp->len;
      if (s_insert > (m_frst-m_last- 1)) {
        if (s_insert > (trunc_on_virtual[m_frst]-trunc_on_virtual[m_last]-1)) {
/*---------------------------------------------------------------------------*/
/* need to insert additional space, find insertion point first               */
/*---------------------------------------------------------------------------*/
	  ipt = m_frst;
	  for (j=m_frst;j>m_last;j--) {
	    if (trunc_on_virtual[j-1] < trunc_on_virtual[j]-1) {
	      ipt = j; break;
	    }
	  }
/*---------------------------------------------------------------------------*/
/* do the shift starting at the insertion point and moving to the end        */
/*---------------------------------------------------------------------------*/
          s_insert -= trunc_on_virtual[m_frst]-trunc_on_virtual[m_last]-1;
	  ASSERT(s_insert > 0);
	  for (j=ipt;j<bsp->length;j++) trunc_on_virtual[j]+=s_insert;
	}
      }
      ddp = ddp->next;
    }
    salpThis = salpThis->next;
  }
/*---------------------------------------------------------------------------*/
/* take care of unaligned regions on the master - distribute evenly on both  */
/* sides of the gap. If uneven number, one more residue shifted to N-terminus*/
/*---------------------------------------------------------------------------*/
  i=-1;
  while (++i < bsp->length-1) {
    if (trunc_aligned[i]==0) {
      j=i; while(trunc_aligned[j]==0 && j < (bsp->length)) j++;
      j--;
      if (j > i && j < (bsp->length)) {
        midpt = i + (j-i)/2;
	for (k=j;k>midpt;k--) {
	  if (trunc_on_virtual[k]<trunc_on_virtual[k+1]-1) {
	    trunc_on_virtual[k] = trunc_on_virtual[k+1]-1;
	  }
	}
        i = j;  
      }
    }
  }
/*  for (j=0;j<bsp->length;j++) printf(" %d",trunc_on_virtual[j]);
  printf("\n"); */
  maxextlen = trunc_on_virtual[bsp->length-1]+1;
  for (i=1;i<bsp->length;i++) {
    if (trunc_on_virtual[i] <= trunc_on_virtual[i-1]) {
      CddSevError("Error determining alignment coordinate system");
    }
  }
/*---------------------------------------------------------------------------*/
/* allocate space for storing the extended alignment plus weights and columns*/
/*---------------------------------------------------------------------------*/
  pCAW = MemNew(nalign * sizeof(CddAlignWeight));
  pCEACol = MemNew(maxextlen * sizeof(CddExtAlignCol));
  pCEACell = MemNew(nalign * sizeof(CddExtAlignCellPtr));
  for (j=0;j<nalign;j++) {
    pCEACell[j] = MemNew(maxextlen * sizeof(CddExtAlignCell));
    for (k=0;k<maxextlen;k++) {
      pCEACell[j][k].seqpos = -1;
      pCEACell[j][k].aatype = 26;
    } 
  }
/*---------------------------------------------------------------------------*/
/* start to fill in the extended alignment data, row 0 is taken by old       */
/* master which is added as fully extended sequence (i.e. not truncated)     */
/*---------------------------------------------------------------------------*/
  i = 0;
  salpThis = salp;
  ddp = salpThis->segs;
  sipThis = ddp->id;
  bspThis = CddRetrieveBioseqById(sipThis,sep);
  pCAW[i].bsp = bspThis;
  spp = SeqPortNew(bspThis, 0, LAST_RESIDUE, Seq_strand_unknown,
                   Seq_code_ncbistdaa);
  buffer = MemNew((bspThis->length)*sizeof(Uint1));
  for (index=0; index<bspThis->length; index++)
    buffer[index] = SeqPortGetResidue(spp);
  spp = SeqPortFree(spp);
  for (j=0;j<bsp->length;j++) {
    k = trunc_on_virtual[j];
    seqpos = j + offset;
    ASSERT (k>=0 && k<maxextlen);
    pCEACell[i][k].seqpos = seqpos;
    pCEACell[i][k].aatype = buffer[seqpos];
    if (pCEACell[i][k].aatype < 0 || pCEACell[i][k].aatype > 26) {
      pCEACell[i][k].aatype = 26;
    }
  }
  MemFree(buffer);
/*---------------------------------------------------------------------------*/
/* continue to fill in alignment data for the remaining sequences            */
/*---------------------------------------------------------------------------*/
  salpThis = salp;
  while (salpThis) {
    i++;
    ddp = salpThis->segs;
    sipThis = ddp->id->next;
    bspThis = CddRetrieveBioseqById(sipThis,sep);
    pCAW[i].bsp = bspThis;
    spp = SeqPortNew(bspThis, 0, LAST_RESIDUE, Seq_strand_unknown,
                     Seq_code_ncbistdaa);
    buffer = MemNew((bspThis->length)*sizeof(Uint1));
    for (index=0; index<bspThis->length; index++)
      buffer[index] = SeqPortGetResidue(spp);
    spp = SeqPortFree(spp);
    while (ddp) {
      for (j = ddp->starts[0]; j<ddp->starts[0]+ddp->len;j++) {
        k = trunc_on_virtual[j - offset];
        ASSERT (k>=0 && k<maxextlen);
	pCEACell[i][k].seqpos = ddp->starts[1] + j - ddp->starts[0];
	pCEACell[i][k].aatype = buffer[pCEACell[i][k].seqpos];
        if (pCEACell[i][k].aatype < 0 || pCEACell[i][k].aatype > 26) {
          pCEACell[i][k].aatype = 26;
        }
      }
      if (ddp->next) {
        m_last = ddp->starts[0]+ddp->len-1-offset;
	m_frst = ddp->next->starts[0]-offset;
	s_last = ddp->starts[1]+ddp->len-1;
	s_frst = ddp->next->starts[1];
	if (s_frst - s_last > 1) {   /* have unaligned residues */
          midpt = s_last + (s_frst - s_last)/2;
	  for (j = s_last+1;j<=midpt;j++) {
	    k = trunc_on_virtual[m_last] + j - s_last;
            ASSERT (k>=0 && k<maxextlen);
	    pCEACell[i][k].seqpos = j;
	    pCEACell[i][k].aatype = buffer[j];
            if (pCEACell[i][k].aatype < 0 || pCEACell[i][k].aatype > 26) {
              pCEACell[i][k].aatype = 26;
            }
	  }
          for (j = s_frst-1; j>midpt; j--) {
	    k = trunc_on_virtual[m_frst] - (s_frst - j);
            ASSERT (k>=0 && k<maxextlen);
	    pCEACell[i][k].seqpos = j;
	    pCEACell[i][k].aatype = buffer[j];
            if (pCEACell[i][k].aatype < 0 || pCEACell[i][k].aatype > 26) {
              pCEACell[i][k].aatype = 26;
            }
	  }	
	}
      }
      ddp = ddp->next;
    }
    MemFree(buffer);
    salpThis = salpThis->next;
  }
  
/*---------------------------------------------------------------------------*/
/* now do statistics on the columns in this matrix to determine weights      */
/*---------------------------------------------------------------------------*/
  for (k=0;k<maxextlen;k++) {
    pCEACol[k].conpos = -1;
    pCEACol[k].occpos = 0;
    pCEACol[k].ntypes = 0;
    pCEACol[k].aatype = 0;
    pCEACol[k].w_occpos = 0.0;
    pCEACol[k].typecount = (Int4 *) MemNew(26*sizeof(Int4));
    pCEACol[k].wtypefreq = (FloatHi *) MemNew(26*sizeof(FloatHi));
    for (j=0;j<26;j++) {
      pCEACol[k].typecount[j] = 0;
      pCEACol[k].wtypefreq[j] = 0.0;
    }
    for (i=0;i<nalign;i++) {
      if (pCEACell[i][k].seqpos != -1 && pCEACell[i][k].aatype < 26) {
        pCEACol[k].typecount[(Int4)pCEACell[i][k].aatype]++;
	pCEACol[k].occpos++;
	if (pCEACol[k].occpos > maxpos) maxpos = pCEACol[k].occpos;
      }
    }
    for (j=0;j<26;j++) if (pCEACol[k].typecount[j]) pCEACol[k].ntypes++;
  }
/*---------------------------------------------------------------------------*/
/* use blocks with the highest position counts to determine weight           */
/*---------------------------------------------------------------------------*/
  for (i=0;i<nalign;i++) {
    pCAW[i].weight = 0.0;
    pCAW[i].nposaligned = 0;
  }
  for (k=0;k<maxextlen;k++) {
    if (pCEACol[k].occpos >= maxpos) {
      ContWeight++;
      for (i=0;i<nalign;i++) {
        if (pCEACell[i][k].seqpos != -1 && pCEACell[i][k].aatype < 26) {
          pCAW[i].nposaligned++;
	  if (pCEACol[k].typecount[pCEACell[i][k].aatype] > 0) {
  	    pCAW[i].weight += 1.0 / (pCEACol[k].ntypes * 
	                             pCEACol[k].typecount[pCEACell[i][k].aatype]);
	  }	
        }
      }    
    }
  }
/*---------------------------------------------------------------------------*/
/* if sequence not covered by weight calculation, give it the expected weight*/
/*---------------------------------------------------------------------------*/
  for (i=0;i<nalign;i++) {
    if (pCAW[i].weight <= 0.0) {
      pCAW[i].weight = (FloatHi) ContWeight * 1.0 / (FloatHi) nalign;
    }
  }
  SumWeight = 0.0;
  for(i=0;i<nalign;i++) {
    SumWeight += pCAW[i].weight;
  }
  for(i=0;i<nalign;i++) pCAW[i].weight /= SumWeight;
/*---------------------------------------------------------------------------*/
/* now determine columns with at least 50% representation and the consensus  */
/*---------------------------------------------------------------------------*/
  ConsensusLen = 0;
  for (k=0;k<maxextlen;k++) {
    for (i=0;i<nalign;i++) {
      if (pCEACell[i][k].seqpos != -1 && pCEACell[i][k].aatype < 26) {
        pCEACol[k].w_occpos += pCAW[i].weight;
	pCEACol[k].wtypefreq[pCEACell[i][k].aatype]+=pCAW[i].weight;
      }
    }
    if (pCEACol[k].w_occpos >= 0.5) {
      pCEACol[k].conpos = ConsensusLen;
      SumWeight = 0.0;
      ConsensusLen++;
      for (j=0;j<26;j++) {
        if (pCEACol[k].wtypefreq[j] > SumWeight) {
	  SumWeight = pCEACol[k].wtypefreq[j];
	  pCEACol[k].aatype = (Uint1) j;
	}
      }
    }
  }  
  buffer = (Uint1 *) MemNew((ConsensusLen+1) * sizeof(Uint1));
  cbuffer = (Char *) MemNew((ConsensusLen+1) * sizeof(Char));
  for (k=0;k<maxextlen;k++) {
    if (pCEACol[k].conpos != -1) {
      buffer[pCEACol[k].conpos] = pCEACol[k].aatype;
    }
  }
/*---------------------------------------------------------------------------*/
/* translate sequence into ascii/character and make new bioseq               */
/*---------------------------------------------------------------------------*/
  smtp = SeqMapTableFind(Seq_code_ncbieaa, Seq_code_ncbistdaa);
  for (i=0;i<ConsensusLen;i++) {
    cbuffer[i] = (Char) SeqMapTableConvert(smtp, buffer[i]);
  }
  sepNew = FastaToSeqBuffEx(cbuffer,&outptr,FALSE,NULL,FALSE);
  *bspCons = (BioseqPtr) sepNew->data.ptrvalue;
  sipThis = (*bspCons)->id;
  oidp = ObjectIdNew();
  oidp->str = MemNew(10*sizeof(Char));
  strcpy(oidp->str,"consensus");
  oidp->str[9]='\0';
  MemFree(sipThis->data.ptrvalue);
  sipThis->data.ptrvalue = oidp;
  sipThis->choice = SEQID_LOCAL;
  sipThis->next = NULL;
  
/*---------------------------------------------------------------------------*/
/* create the pairwise alignment old-master(not truncated) new-master        */
/*---------------------------------------------------------------------------*/
  bspThis = pCAW[0].bsp;
  pCDea = CddExpAlignNew();
  CddExpAlignAlloc(pCDea,bspThis->length);
  pCDea->ids = SeqIdDup(bspThis->id);
  pCDea->ids->next = SeqIdDup((*bspCons)->id);
  pCDea->bIdAlloc = TRUE;
  for (k=0;k<maxextlen;k++) {
    if (pCEACol[k].conpos != -1) {
      pCDea->adata[pCEACell[0][k].seqpos]=pCEACol[k].conpos;
    }
  }
  *salpCons = CddExpAlignToSeqAlign(pCDea);
  pCDea = CddExpAlignFree(pCDea);

/*---------------------------------------------------------------------------*/
/* create a multiple alignment salpNew where the consensus is the new master */
/*---------------------------------------------------------------------------*/
  salpNew = NULL;
  for (i=0;i<nalign;i++) {
    pCDea = CddExpAlignNew();
    CddExpAlignAlloc(pCDea,ConsensusLen);
    pCDea->ids = SeqIdDup((*bspCons)->id);
    pCDea->ids->next = SeqIdDup(pCAW[i].bsp->id);
    pCDea->bIdAlloc = TRUE;
    for (k=0;k<maxextlen;k++) {
      if (pCEACol[k].conpos != -1) {
        pCDea->adata[pCEACol[k].conpos]=pCEACell[i][k].seqpos;
      }
    }
    salpThis = CddExpAlignToSeqAlign(pCDea);
    pCDea = CddExpAlignFree(pCDea);
    if (!salpNew) {
      salpNew = salpThis;
    } else {
      salpTail->next = salpThis;
    }
    salpTail = salpThis;
  }

/*---------------------------------------------------------------------------*/
/* deallocate memory for this function                                       */
/*---------------------------------------------------------------------------*/
  MemFree(cbuffer);
  MemFree(buffer);
  for (j=0;j<nalign;j++) MemFree(pCEACell[j]);
  MemFree(pCEACell);
  for (j=0;j<maxextlen;j++) {
    MemFree(pCEACol[j].typecount);
    MemFree(pCEACol[j].wtypefreq);    
  }
  MemFree(pCEACol);
  MemFree(pCAW);
  MemFree(trunc_aligned);
  MemFree(trunc_on_virtual);
  return(salpNew);
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
    bspFake  = (BioseqPtr) CddBioseqCopy(NULL,bsp,0,bsp->length-1,0,FALSE);

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
    posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numalign);
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
    posPurgeMatches(posSearch, compactSearch);
    posComputeExtents(posSearch, compactSearch);
    posComputeSequenceWeights(posSearch, compactSearch, 1.0);
    posCheckWeights(posSearch, compactSearch);
    posSearch->posFreqs = (Nlm_FloatHi **) posComputePseudoFreqs(posSearch, compactSearch, TRUE);
    CddposFreqsToMatrix(posSearch,compactSearch);
    posScaling(posSearch, compactSearch);

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


Seq_Mtf* LIBCALL CddCalcPSSM(SeqAlignPtr listOfSeqAligns, BioseqPtr bsp,
                             double Weight, double ScaleFactor) {
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
    Int4                i, j, jj;
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
    for (j=0; j<compactSearch.alphabetSize;j++) {
      jj = CddGetNewIndexForThreading(Input[j], Output);
      if (jj != -1) {
        for (i=0; i<compactSearch.qlength;i++) {
          pssm->ww[i][jj] = ThrdRound(posSearch.posMatrix[i][j]*ScaleFactor*Weight);
          pssm->freqs[i][jj] = ThrdRound(posSearch.posFreqs[i][j]*ScaleFactor);
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

    /* clean up */
    MemFree(Input);
    MemFree(Output);
    MemFree(Coverage);
    CddposFreeMemory(&posSearch);

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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* modified after PGPReadBlastOptions from blastpgp.c, as of 12/27/00        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static PGPBlastOptionsPtr CddReadBlastOptions(BioseqPtr bsp, Int4 iPseudo, CharPtr matrix_name)
{
    PGPBlastOptionsPtr  bop;
    BLAST_OptionsBlkPtr options;
    SeqEntryPtr         sep;
    Boolean             is_dna;
    ObjectIdPtr         obidp;

    bop = MemNew(sizeof(PGPBlastOptions));
    bop->blast_database   = StringSave("nr");
    bop->number_of_descriptions = 500;
    bop->number_of_alignments = 500;
    bop->believe_query = FALSE;
    options = BLASTOptionNew("blastp",TRUE);                /* gapped Blast! */
    bop->options = options;
    bop->query_bsp = bsp; 
    BLASTOptionSetGapParams(options, matrix_name, 0, 0);
    bop->is_xml_output = FALSE;
    bop->align_options = 0;
    bop->print_options = 0;
    options->required_start = 0;
    options->required_end = -1;
    options->window_size = 40;
    options->dropoff_2nd_pass  = 7.0;
    options->expect_value  = 10.0;
    options->hitlist_size = 500;
    options->two_pass_method  = FALSE;
    options->multiple_hits_only  = TRUE;
    if (strcmp(matrix_name,"BLOSUM62") == 0) {
      options->gap_open = 11;
      options->gap_extend = 1;
    }
    options->decline_align = INT2_MAX;
    options->gap_x_dropoff = 15;
    options->gap_x_dropoff_final = 25;
    options->gap_trigger = 22.0;
    options->filter_string = StringSave("F");
    options->number_of_cpus = 1;
    options->isPatternSearch = FALSE;
    options->ethresh = 0.0;                  /* zeroed out, will not be used */
    options->pseudoCountConst = iPseudo;
    options->maxNumPasses = 1;
    options->hsp_range_max  = 0;
    options->block_width  = 20; /* Default value - previously '-L' parameter */
    options->tweak_parameters = TRUE;
    options->smith_waterman = FALSE;
    if (bop->options->tweak_parameters) {
      /*allows for extra matches in first pass of screening,
        hitlist_size */
      bop->options->original_expect_value = bop->options->expect_value;
      bop->options->hitlist_size *= 2; 
    }

    options = BLASTOptionValidate(options, "blastp");
    if (options == NULL) return NULL;
    bop->fake_bsp = bop->query_bsp;
#if 0
      BioseqNew();
    bop->fake_bsp->descr = bop->query_bsp->descr;
    bop->fake_bsp->repr = bop->query_bsp->repr;
    bop->fake_bsp->mol = bop->query_bsp->mol;
    bop->fake_bsp->length = bop->query_bsp->length;
    bop->fake_bsp->seq_data_type = bop->query_bsp->seq_data_type;
    bop->fake_bsp->seq_data = bop->query_bsp->seq_data;
    bop->fake_bsp->id = bop->query_bsp->id;
#endif
    /* bop->fake_bsp->seq_data = BSDup(bop->query_bsp->seq_data); */
    /* bop->fake_bsp->id = SeqIdDup(bop->query_bsp->id); */
    /* bug fix -- Lewis, Dave and Aron */
    /* bop->query_bsp->id = SeqIdSetFree(bop->query_bsp->id); */
    /* BLASTUpdateSeqIdInSeqInt(options->query_lcase_mask, bop->fake_bsp->id); */
    return bop;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* updated version of the PSSM calculation routine, modified from blastpgp.c */
/* as of 12/27/2000                                                          */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Seq_Mtf * LIBCALL CddDenDiagCposComp2(BioseqPtr bspFake, Int4 iPseudo,
                                      SeqAlignPtr salp, CddPtr pcdd,
				      BioseqPtr bspOut, double Weight,
				      double ScaleFactor, CharPtr matrix_name)
{
  PGPBlastOptionsPtr      bop;
  BlastSearchBlkPtr       search;
  posSearchItems         *posSearch = NULL;
  compactSearchItems     *compactSearch = NULL;
  ValNodePtr              ColumnHead, LetterHead, newRow, RowHead, newLetter;
  SeqCodeTablePtr         sctp;
  Int4                    numSeqs, numASeqs, alignLength, index, a, iNew, i;
  Int4                    scale_factor, j, OutLen, jj;
  Char                    takeCkptFileName[PATH_MAX];
  Char                    cseqFileName[PATH_MAX];
  FILE                   *fp;
  Boolean                 bHasConsensus;
  Nlm_FloatHi            *AInf, SumAInf = 0.0;
  char*                   Input;    /* order of residue-types from CD's      */
  char*                   Output;   /*order of res-types needed for threading*/
  Boolean*                Coverage; /* for making sure all columns are filled*/
  Seq_Mtf*                pssm;

  if (pcdd) {
    bHasConsensus = CddHasConsensus(pcdd);
  } else {                                        /* determine from seqalign */
    bHasConsensus = SeqAlignHasConsensus(salp);
  }
  if (iPseudo <= 0) {                             /* need to determine first */
    AInf = SeqAlignInform(salp, bspFake, bHasConsensus); 
    for (i=0;i<bspFake->length;i++) SumAInf += AInf[i];
    if (SumAInf > 84) iPseudo = 10;                    /* purely empirical   */
    else if (SumAInf > 55) iPseudo = 7;
    else if (SumAInf > 43) iPseudo = 5;
    else if (SumAInf > 41.5) iPseudo = 4;
    else if (SumAInf > 40) iPseudo = 3;
    else if (SumAInf > 39) iPseudo = 2;
    else iPseudo = 1;
    MemFree(AInf);
    if (pcdd) printf("%s AInf:%6.1f PseudoCt: %d\n",pcdd->name,SumAInf, iPseudo);
  }
  if (NULL == matrix_name) {
    static CharPtr defaultMatrix = "BLOSUM62";
    matrix_name = defaultMatrix;
  }
  bop = (PGPBlastOptionsPtr) CddReadBlastOptions(bspFake, iPseudo, matrix_name);
  if (pcdd) {          /* set up search with nr database size for CD dumping */
    search = BLASTSetUpSearchWithReadDb(bop->fake_bsp, "blastp", 
                                        bop->query_bsp->length, 
                                        bop->blast_database,
                                        bop->options, NULL);
  } else {          /* or avoid using nr when calculating PSSM for threading */
    search = BLASTSetUpSearch(bop->fake_bsp, "blastp",
                              bop->query_bsp->length,
			      0, NULL, bop->options, NULL);
  
  }

  posSearch = NULL;
  compactSearch = NULL;
  posSearch = (posSearchItems *)MemNew(sizeof(posSearchItems));

  compactSearch = compactSearchNew(compactSearch);
  copySearchItems(compactSearch, search, bop->options->matrix);
  posInitializeInformation(posSearch,search);

/*---------------------------------------------------------------------------*/
/* section modelled after BposComputation, posit.c as of 12/28/2000          */
/*---------------------------------------------------------------------------*/
  search->posConverged = FALSE;
  numSeqs = CddCountDenDiagSeqAligns(salp, &numASeqs) + 1; 
  alignLength = bspFake->length;
  posAllocateMemory(posSearch, compactSearch->alphabetSize, 
                       compactSearch->qlength, numSeqs);
  CddposProcessAlignment(posSearch,compactSearch,salp,numSeqs,alignLength,bspFake);
  posSearch->posNumSequences = numSeqs;
  posPurgeMatches(posSearch, compactSearch);
  if (bHasConsensus) posCancel(posSearch,compactSearch,0,0,0,compactSearch->qlength);
  posComputeExtents(posSearch, compactSearch);
  posComputeSequenceWeights(posSearch, compactSearch, 1.0);
  posCheckWeights(posSearch, compactSearch);
/*  printf("Calling posComputePseudoFreqs with %s\n",compactSearch->standardMatrixName);
  fflush(stdout); */
  posSearch->posFreqs = posComputePseudoFreqs(posSearch,compactSearch,TRUE);
  if (NULL == search->sbp->posFreqs)
    search->sbp->posFreqs =  allocatePosFreqs(compactSearch->qlength, compactSearch->alphabetSize);
  copyPosFreqs(posSearch->posFreqs,search->sbp->posFreqs, compactSearch->qlength, compactSearch->alphabetSize);

/*---------------------------------------------------------------------------*/
/* Construct name for checkpoint file and write out (if in a CD context)     */
/*---------------------------------------------------------------------------*/
  if (pcdd) {
    strcpy(takeCkptFileName,pcdd->name);
    strcat(takeCkptFileName,CKPTEXT);
    if (NULL != takeCkptFileName)
      posTakeCheckpoint(posSearch, compactSearch, takeCkptFileName, NULL);
  }

  posFreqsToMatrix(posSearch,compactSearch, NULL, 1);
  posScaling(posSearch, compactSearch);

  if (pcdd) {
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

  } else  { /* if pcdd is empty, fill in the Seq_Mtf data structure */
    Input = (char*) MemNew(sizeof(InputOrder));
    Output = (char*) MemNew(sizeof(OutputOrder));
    StrCpy(Input, InputOrder);
    StrCpy(Output, OutputOrder);
    OutLen = StrLen(Output);
    pssm = NewSeqMtf(compactSearch->qlength, OutLen);
    Coverage = (Boolean*) MemNew(sizeof(Boolean) * OutLen);  /* set to 0 by MemNew */
    for (j=0; j<compactSearch->alphabetSize;j++) {
      jj = CddGetNewIndexForThreading(Input[j], Output);
      if (jj != -1) {
        for (i=0; i<compactSearch->qlength;i++) {
          pssm->ww[i][jj] = ThrdRound(posSearch->posMatrix[i][j]*ScaleFactor*Weight);
          pssm->freqs[i][jj] = ThrdRound(posSearch->posFreqs[i][j]*ScaleFactor);
          Coverage[jj] = TRUE;
        }
      }
    }
    for (i=0; i<OutLen; i++) {
      if (Coverage[i] == FALSE) {
        ASSERT(FALSE);  /* should never get here */
      }
    }
    MemFree(Input);
    MemFree(Output);
    MemFree(Coverage);
  }

  if (bspOut && pcdd) {
    strcpy(cseqFileName,pcdd->name);
    strcat(cseqFileName,CSEQEXT);
    fp = FileOpen(cseqFileName, "w");
    if (NULL == fp) {
      BlastConstructErrorMessage("CddDenDiagCposComp2", "Could not open fasta file", 1, NULL);
    }
    BioseqToFasta(bspOut,fp,FALSE);
    FileClose(fp);
  }

/*---------------------------------------------------------------------------*/
/* clean up                                                                  */
/*---------------------------------------------------------------------------*/
  posSearch->NumSequences = posSearch->posNumSequences;
  posSearch->QuerySize = alignLength;
  CddposFreeMemory(posSearch);
  MemFree(posSearch);
  CddposFreeMemory2(compactSearch);
  MemFree(compactSearch);
  search = BlastSearchBlkDestruct(search);
  bop->options = BLASTOptionDelete(bop->options);
  MemFree(bop->blast_database);
  MemFree(bop);
  if (NULL == pcdd) {        /* return pssm if used from within the threader */
    return(pssm);
  }
  return NULL;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Free a Cdd data structure - works on a single CD only (not linked lists)  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
CddPtr LIBCALL CddFreeCarefully(CddPtr pcdd)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  ScorePtr     nxtscore, thisscore;
  SeqAnnotPtr  sap, sapThis;
  SeqEntryPtr  sep, sepThis;
  ValNodePtr   vnp, vnpThis;
  
  if(pcdd == NULL) return NULL;
  if (NULL != pcdd->name) pcdd->name = MemFree(pcdd->name);
  if (NULL != pcdd->id)   pcdd->id   = (CddIdSetPtr) CddIdSetFree(pcdd->id);

  sap = pcdd->seqannot;
  while (sap) {
    sapThis = sap->next;
    sap = SeqAnnotFree(sap);
    sap = sapThis;
  }

  vnp = pcdd->description;
  while (vnp) {
    vnpThis = vnp->next;
    switch (vnp->choice) {
      case CddDescr_othername:
        vnp->data.ptrvalue = MemFree(vnp->data.ptrvalue);
	break;
      case CddDescr_category:
        vnp->data.ptrvalue = MemFree(vnp->data.ptrvalue);
	break;
      case CddDescr_comment:
        vnp->data.ptrvalue = MemFree(vnp->data.ptrvalue);
	break;
      case CddDescr_reference:
        vnp->data.ptrvalue = PubFree(vnp->data.ptrvalue);
	break;
      case CddDescr_create_date:
        vnp->data.ptrvalue = DateFree(vnp->data.ptrvalue);
	break;
      case CddDescr_tax_source:
        vnp->data.ptrvalue = OrgRefFree(vnp->data.ptrvalue);
	break;
      case CddDescr_source:
        vnp->data.ptrvalue = MemFree(vnp->data.ptrvalue);
	break;
      default:
	break;
    }
    vnp = MemFree(vnp);
    vnp = vnpThis;
  }

  if (NULL != pcdd->features) {
    pcdd->features = BiostrucAnnotSetFree(pcdd->features);
  }

  if (NULL != pcdd->sequences) {
    bssp = pcdd->sequences->data.ptrvalue;
    if (NULL != bssp) {
      sep = bssp->seq_set;
      while (sep) {
        sepThis = sep->next;
	if (sep->choice == 1) {
          bsp = sep->data.ptrvalue;
	  if (bsp->id->choice != SEQID_LOCAL)BioseqUnlock(bsp);
	  BioseqFree(bsp);
        }
        MemFree(sep);
	sep = sepThis;
      }
      bssp = MemFree(bssp);
    }
    pcdd->sequences = ValNodeFree(pcdd->sequences);
  }

  if(NULL != pcdd->profile_range) {
    pcdd->profile_range = SeqIntFree(pcdd->profile_range);
  }
  if (NULL != pcdd->trunc_master) {
    pcdd->trunc_master = BioseqFree(pcdd->trunc_master);
  }
  if (NULL != pcdd->posfreq) {
    vnp = pcdd->posfreq->row_labels; if (NULL != vnp) ValNodeFree(vnp);
    vnp = pcdd->posfreq->columns;    if (NULL != vnp) ValNodeFree(vnp);
    pcdd->posfreq = MemFree(pcdd->posfreq);
  }
  if (NULL != pcdd->scoremat) {
    vnp = pcdd->scoremat->columns;    if (NULL != vnp) ValNodeFree(vnp);
    pcdd->scoremat = MemFree(pcdd->scoremat);
  }
  if (NULL != pcdd->distance) {
    nxtscore = pcdd->distance->scores;
    while (nxtscore) {
      if (NULL != nxtscore->id) ObjectIdFree(nxtscore->id);
      thisscore = nxtscore;
      nxtscore = nxtscore->next;
      thisscore = MemFree(thisscore);
    }
    pcdd->distance = MemFree(pcdd->distance);
  }
  return MemFree(pcdd);
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* reindex a SeqAlign (a linked list of pairwise Den-Diag Master-Slave Seq-  */
/* aligns) to a new Master, pointed out by it's SeqId. Will pick the first   */
/* instance of that SeqId in the alignment, and return a SeqAlign of the     */
/* same format as the input was                                              */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SeqAlignPtr LIBCALL CddReindexSeqAlign(SeqAlignPtr salp, SeqIdPtr sipMaster,
                                       BioseqSetPtr bssp)
{
  SeqAlignPtr    salpThis, salpNew, salpHead, salpTail;
  DenseDiagPtr   ddp;
  CddExpAlignPtr pCDea, pCDea1, pCDeaR;
  Int4           iNotThis = -1, i = -1, mlength;
  BioseqPtr      bsp;

  bsp = CddRetrieveBioseqById(sipMaster, bssp->seq_set);
  if (!bsp) return NULL;
  mlength = bsp->length;
  salpThis = salp; while(salpThis) {
    iNotThis++;
    ddp = salpThis->segs;
    if (CddSameSip(sipMaster,ddp->id->next)) break;
    else salpThis = salpThis->next; 
  }
  pCDea1 = (CddExpAlignPtr) SeqAlignToCddExpAlign(salpThis,bssp->seq_set);
  pCDeaR = (CddExpAlignPtr) InvertCddExpAlign(pCDea1, bssp->seq_set);
  salpNew = (SeqAlignPtr) CddExpAlignToSeqAlign(pCDeaR);
  pCDeaR = CddExpAlignFree(pCDeaR);
  salpTail = salpNew;
  salpHead = salp;
  while (salpHead) {
    i++; 
    if (i!= iNotThis) {    
      pCDea = (CddExpAlignPtr) SeqAlignToCddExpAlign(salpHead, bssp->seq_set);
      pCDeaR = (CddExpAlignPtr) CddReindexExpAlign(pCDea1,mlength,pCDea,0,i);
      salpThis = (SeqAlignPtr) CddExpAlignToSeqAlign(pCDeaR);
      pCDeaR = CddExpAlignFree(pCDeaR);
      pCDea = CddExpAlignFree(pCDea);
      salpTail->next = salpThis;
      salpTail = salpThis;
    }
    salpHead = salpHead->next;
  }
  
  pCDea1 = CddExpAlignFree(pCDea1);
  return (salpNew);
}
