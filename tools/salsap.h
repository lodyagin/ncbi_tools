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
* File Name:  salsap.h
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.35 $
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

#ifndef _SALSAP_
#define _SALSAP_

#include <salsa.h>
#include <seqport.h>

extern Boolean is_dim1seqalign (SeqAlignPtr salp);
extern Boolean is_dim2seqalign (SeqAlignPtr salp);

extern Boolean is_fasta_seqalign (SeqAlignPtr salp);

extern SeqAnnotPtr SeqAnnotForSeqAlign (SeqAlignPtr salp);
extern SeqAlignPtr is_salp_in_sap (SeqAnnotPtr sap, Uint1 choice);

extern SeqAlignPtr seqalign_list_free (SeqAlignPtr salp);

/**********************************************************************
*
*   Functions for SeqAlign (1 or linked list)
*
***********************************************************************/
extern SeqIdPtr    SeqAlignId (SeqAlignPtr salp, Int2 index);
extern Boolean     FindSeqIdinSeqAlign (SeqAlignPtr salphead, SeqIdPtr sip);
extern SeqAlignPtr SeqAlignIdReplace (SeqAlignPtr salp, Int2 index, SeqIdPtr newsip);
extern Boolean     SeqAlignSeqLocComp (SeqAlignPtr salphead, ValNodePtr vnp);
extern Int4        SeqAlignLength (SeqAlignPtr salp);
extern Uint1       SeqAlignStrand (SeqAlignPtr salp, Int2 index);
extern Int4        SeqAlignStart (SeqAlignPtr salp, Int2 index);
extern Int4        SeqAlignStop (SeqAlignPtr salp, Int2 index);
extern Uint1       SeqAlignMolType (SeqAlignPtr salp);
extern Int4        SeqAlignBestScore (SeqAlignPtr salp);
extern SeqAlignPtr SeqAlignBestHit (SeqAlignPtr salp, Int4 length, Int4 threshold);

extern Int4        SeqAlignGapCount (SeqAlignPtr salp);

/**********************************************************************
*   Update SeqAlign
***********************************************************************/
extern void        SeqAlignStartUpdate (SeqAlignPtr salp, SeqIdPtr target_sip, Int4 offset, Uint1 strand);

/**********************************************************************
*
*       SeqAlignLink(head, new)
*       link the new align to the end of head align. return the
*       start of the linked chain
*
**********************************************************************/
extern SeqAlignPtr SeqAlignLink PROTO((SeqAlignPtr head, SeqAlignPtr a_new));

extern ValNodePtr  SeqLocListFromSeqAlign (SeqAlignPtr salp);
extern SeqLocPtr   SeqLocFromSeqAlign (SeqAlignPtr salp, SeqIdPtr sip);
extern SeqLocPtr   SeqLocMixFromSeqAlign (SeqAlignPtr salp, SeqIdPtr sip);
extern ValNodePtr  WholeSeqLocListFromSeqAlign (SeqAlignPtr salp);

extern SeqAnnotPtr CompSeqAnnotFree (SeqAnnotPtr sap);
extern SeqAlignPtr CompSeqAlignFree (SeqAlignPtr salp);
extern SeqAnnotPtr SeqAlignBoolSegCpy (SeqAnnotPtr sap, Int4 from, Int4 to);
extern SeqAlignPtr SeqAlignDenseSegToBoolSeg (SeqAlignPtr salp);
extern SeqAnnotPtr SeqAnnotDenseSegToBoolSeg (SeqAnnotPtr sap);
extern SeqAlignPtr SeqAlignBoolSegToDenseSeg (SeqAlignPtr salp);
extern SeqAnnotPtr SeqAnnotBoolSegToDenseSeg (SeqAnnotPtr sap);
extern void CompSeqAlignPrint (SeqAlignPtr salp);

extern SeqAlignPtr build_seqalign_fromstart (Int2 dim, Int2 numseg, SeqIdPtr sip, Int4Ptr starts, Int4Ptr lens);

/*********************************************************
***
***  SeqAlignDup
***
**********************************************************/
extern SeqAlignPtr SeqAlignDup (SeqAlignPtr salp);
extern SeqAlignPtr SeqAlignDupRegion (SeqAlignPtr salp, Int2 to_numseg, Int4 subseg, Boolean first_part);

/*********************************************************
***
***  SeqAlignAdd
***
**********************************************************/
extern SeqAlignPtr SeqAlignAdd (SeqAlignPtr *salp_head, SeqAlignPtr salp);
extern SeqAlignPtr SeqAlignDupAdd (SeqAlignPtr *salp_head, SeqAlignPtr salp, Int2 to_numseg, Int4 subseg, Boolean first_part);
extern SeqAlignPtr SeqAlignEndExtend (SeqAlignPtr sap, Int4 start1, Int4 start2, Int4 stop1, Int4 stop2, Int4 x1, Int4 y1, Int4 x2, Int4 y2, Uint1 strand1, Uint1 strand2);

/*********************************************************
***
***  SeqAlignTrunc = truncates the extremitites of seqalign salp
***  SeqAlignMap = delete the segments at extremities when 1rst sequence has gaps
***
**********************************************************/
extern SeqAlignPtr SeqAlignTrunc (SeqAlignPtr salp, Int4 from, Int4 to);
extern SeqAlignPtr SeqAlignMap (SeqAlignPtr salp);

/*********************************************************
***
***  SeqAnnotMerge
***      return_salp =  TRUE  if seqalign1 precedes seqalign2
***                     FALSE if otherwise
**********************************************************/
extern SeqAlignPtr SeqAlignMerge (SeqAlignPtr salp1, SeqAlignPtr salp2, Boolean return_salp);
extern SeqAnnotPtr SeqAnnotMerge (SeqAnnotPtr sap1, SeqAnnotPtr sap2, Boolean return_salp);

extern SeqAlignPtr SeqAlignExtend (SeqAlignPtr salp1, SeqAlignPtr salp2);

extern SeqAlignPtr check_salp_forlength (SeqAlignPtr salp);
extern SeqAlignPtr check_salp_forstrand (SeqAlignPtr salp);
extern SeqAnnotPtr LocalAlignToSeqAnnotDimn (ValNodePtr seqvnp, SeqIdPtr seqsip, ValNodePtr fromp, Int2 nbseq, Int4 lens, ValNodePtr strands, Boolean trunc_emptyends);
extern SeqAnnotPtr LocalAlignToSeqAnnotCompDimn (ValNodePtr seqvnp, SeqIdPtr seqsip, Int2 nbseq, Int4 lens);

/*******************************************
***
***   Delete
***
********************************************/
extern Boolean     SeqAlignIDCache (SeqAlignPtr salp, SeqIdPtr sip);
extern SeqAlignPtr SeqAlignIDUncache (SeqAlignPtr salphead, SeqIdPtr sip);
extern SeqAlignPtr SeqAlignIDUncacheAll (SeqAlignPtr salphead);
extern SeqAlignPtr SeqAlignIDDelete (SeqAlignPtr salp, SeqIdPtr sip);
extern void        DelAlignItem (SeqEntryPtr sep, SeqIdPtr sip);

/*******************************************
***
***   DeleteRegion
***
********************************************/
extern SeqAlignPtr SeqAlignDeleteByLoc (SeqLocPtr slp, SeqAlignPtr salp);
extern SeqAlignPtr DeleteRegion (SeqIntPtr vnp, SeqAlignPtr salp);

/*********************************************************
***
***  DenseDiagPtr procedures
***
**********************************************************/
extern SeqAlignPtr  DenseDiagToDenseSeg (SeqAlignPtr salp, Boolean add_ends);
extern SeqAlignPtr  DenseSegToDenseDiag (SeqAlignPtr salp);

extern DenseDiagPtr DenseDiagCreate (Int4 dim, SeqIdPtr id, Int4Ptr start, Int4 lens, Uint1Ptr strands, ScorePtr scores);
extern DenseDiagPtr DenseDiagDup (DenseDiagPtr ddp);
extern DenseDiagPtr DenseDiagAdd (DenseDiagPtr *ddp_head, DenseDiagPtr ddp);
extern DenseDiagPtr DenseDiagInsert (DenseDiagPtr ddp_before, DenseDiagPtr ddp);
extern DenseDiagPtr DenseDiagPrecede (DenseDiagPtr ddp_after, DenseDiagPtr *ddp);
extern DenseDiagPtr DenseDiagSortAdd (DenseDiagPtr *ddp_head, DenseDiagPtr ddp);
extern void DenseDiagPrint (ValNodePtr ddp);
extern Boolean IS_seqidindensediag (SeqIdPtr sip, ValNodePtr ddia_list, SeqAlignPtr salp, Int2 index, Int4 from, Int4 to, DenseDiagPtr *block, Int2 intersalpwidth);
extern DenseDiagPtr GetDenDiag (SeqAlignPtr salp, Int2 index, Int2 *index_entry);
extern SeqAlignPtr SeqAlignDiagAdd (SeqAlignPtr headp, Int4 pos, Int4 len);
extern SeqAlignPtr DenseDiagAlign (SeqAlignPtr salp, DenseDiagPtr dendia);

/***********************************************************************
***    
***    DenDiagToSeqLoc
***      read SeqAnnotPtr-densediag
***      n: number of sip
***      return list of ValNodePtr-SeqLocPtr
***
************************************************************************/
extern ValNodePtr DenDiagToSeqLoc (SeqAnnotPtr sap, ValNodePtr adpslp, Int2 blastscore_threshold, Int2 *n);

/***********************************************************************
***
***   SeqLocToFastaSeqAlign
***      creates a SeqAlign where all gaps are at the sequence ends
***      FASTA-style
***
***********************************************************************/
extern SeqAlignPtr SeqLocToFastaSeqAlign (ValNodePtr vnp);

/***********************************************************************
***
***   SeqEntryToSeqAlignFunc
***      creates a SeqAlign from all bioseqs in SeqEntry sep
***      if master NOT NULL, all Bioseqs are aligned to it
***
***   SeqEntryToSeqAlign
***      calls SeqEntryToSeqAlignFunc, return SeqAnnotPtr
***
*************************************************************************/
extern SeqAlignPtr SeqEntryToSeqAlignFunc (SeqEntryPtr sep, SeqLocPtr master, Uint1 bsp_mol, Int2 method);
extern SeqAnnotPtr SeqEntryToSeqAlign (SeqEntryPtr sep, Uint1 bsp_mol);

extern void        ReplaceSeqAlignInSeqEntry (Uint2 entityID, Uint2 itemID, SeqAlignPtr salp);
extern Pointer     FindSeqAlignInSeqEntry (SeqEntryPtr sep, Uint1 choice);

extern Int4        readbuff_fromseqalign (SeqPortPtr spp, SeqAlignPtr salp,  Int2 index, CharPtr buffer, Int4 from, Int4 to, Int4 offset, Boolean strand);

extern SeqAlignPtr aaSeqAlign_to_dnaSeqAlign (SeqAlignPtr sap, ValNodePtr vnp, ValNodePtr framep);
extern SeqAnnotPtr aaSeqAnnot_to_dnaSeqAnnot (SeqAnnotPtr sap, ValNodePtr vnp, ValNodePtr framep);


extern SeqAlignPtr SortSeqAlign (SeqAlignPtr PNTR salp);
extern SeqAlignPtr SortSeqAlignFromList (SeqAlignPtr salp, Int2Ptr sortlst);


#endif
