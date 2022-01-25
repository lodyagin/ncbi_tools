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
* $Revision: 6.39 $
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


	
extern SeqAlignPtr build_seqalign_fromstart (Int2 dim, Int2 numseg, 
	SeqIdPtr sip, Int4Ptr starts, Int4Ptr lens);

extern SeqAlignPtr SeqLocToFastaSeqAlign (ValNodePtr vnp);

extern SeqAnnotPtr LocalAlignToSeqAnnotDimn (ValNodePtr seqvnp, SeqIdPtr seqsip, 
	ValNodePtr fromp, Int2 nbseq, Int4 lens, ValNodePtr strands, 
	Boolean trunc_emptyends);

extern DenseDiagPtr DenseDiagCreate (Int4 dim, SeqIdPtr id, Int4Ptr starts, 
	Int4 len, Uint1Ptr strands, ScorePtr scores);

extern SeqAlignPtr SeqEntryAlignToMasterFunc (SeqEntryPtr sep, SeqLocPtr master, 
	Uint1 bsp_mol, Int2 method);
	
extern SeqAnnotPtr SeqEntryToSeqAlign (SeqEntryPtr sep, Uint1 bsp_mol);

extern Pointer FindSeqAlignInSeqEntry (SeqEntryPtr sep, Uint1 choice);

extern SeqAlignPtr is_salp_in_sap (SeqAnnotPtr sap, Uint1 choice);

extern Boolean is_dim1seqalign (SeqAlignPtr salp);
	
extern Boolean is_dim2seqalign (SeqAlignPtr salp);

extern Int4 SeqAlignLength (SeqAlignPtr salp);

extern Int4 SeqAlignLengthForId (SeqAlignPtr salp, SeqIdPtr sip);

extern Uint1 SeqAlignMolType (SeqAlignPtr salp);

extern SeqIdPtr SeqAlignIDList (SeqAlignPtr salp);

extern SeqIdPtr SeqAlignId (SeqAlignPtr salp, Int2 index);

extern Boolean FindSeqIdinSeqAlign (SeqAlignPtr salphead, SeqIdPtr sip);

extern Boolean SeqAlignSeqLocComp (SeqAlignPtr salphead, ValNodePtr vnp);

extern Uint1 SeqAlignStrand (SeqAlignPtr salp, Int2 index);

extern Int4 SeqAlignStart (SeqAlignPtr salp, Int2 index);

extern Boolean get_pos_from_salp (SeqAlignPtr salp, Int4 pos, 
	Int4 PNTR offset, Int4Ptr PNTR startp, Int4Ptr PNTR lenp, 
	Int4 PNTR numseg);

extern Int4 SeqAlignStop (SeqAlignPtr salp, Int2 index);

/**SeqLoc**/
extern ValNodePtr SeqLocListFromSeqAlign (SeqAlignPtr salp);

extern SeqLocPtr SeqLocFromSeqAlign (SeqAlignPtr salp, SeqIdPtr sip);

extern SeqLocPtr SeqLocMixFromSeqAlign (SeqAlignPtr salp, SeqIdPtr sip);

extern ValNodePtr SeqLocListOfBioseqsFromSeqAlign (SeqAlignPtr salp);

/**Alignment Score**/
extern Int4 SeqAlignBestScore (SeqAlignPtr salp);

extern SeqAlignPtr SeqAlignBestHit (SeqAlignPtr salp, Int4 length, Int4 threshold);
	
extern Int4 SeqAlignGapCount (SeqAlignPtr salp);

/**Sequence**/
extern Int4 readbuff_fromseqalign (SeqPortPtr spp, SeqAlignPtr salp, Int2 index, 
	CharPtr buffer, Int4 from, Int4 to, Int4 offset, Boolean strand);

/**misc**/
extern Boolean is_fasta_seqalign (SeqAlignPtr salp);

/**
*Functions taking a SeqAlign and returning a pointer to the "same" SeqAlign
**/

/**Link, Merge**/
extern SeqAlignPtr SeqAlignLink(SeqAlignPtr head, SeqAlignPtr a_new);

extern SeqAlignPtr SeqAlignMerge (SeqAlignPtr salp1, SeqAlignPtr salp2, 
	Boolean return_salp);
	
extern SeqAnnotPtr SeqAnnotMerge (SeqAnnotPtr sap1, SeqAnnotPtr sap2, 
	Boolean return_salp);

extern DenseDiagPtr DenseDiagLink (DenseDiagPtr *ddp_head, DenseDiagPtr ddp);

extern DenseDiagPtr DenseDiagInsert (DenseDiagPtr ddp_before, DenseDiagPtr ddp);

extern DenseDiagPtr DenseDiagPrecede (DenseDiagPtr ddp_after, DenseDiagPtr *ddp);

extern DenseDiagPtr DenseDiagLinkSort (DenseDiagPtr *ddp_head, DenseDiagPtr ddp);

/**Sort**/
extern SeqAlignPtr SortSeqAlign (SeqAlignPtr PNTR salp);

extern SeqAlignPtr SortSeqAlignFromList (SeqAlignPtr salp, Int2Ptr sortlst);

/**SeqId**/
extern SeqAlignPtr SeqAlignIdReplace (SeqAlignPtr salp, Int2 index, 
	SeqIdPtr newsip);

extern SeqAlignPtr SeqAlignIDDelete (SeqAlignPtr salphead, SeqIdPtr sip);

extern Boolean SeqAlignIDCache (SeqAlignPtr salphead, SeqIdPtr sip);

extern SeqAlignPtr SeqAlignIDUncache (SeqAlignPtr salphead, SeqIdPtr sip);

extern SeqAlignPtr SeqAlignIDUncacheAll (SeqAlignPtr salphead);

/**Offset positions**/
extern void SeqAlignStartUpdate (SeqAlignPtr salp, SeqIdPtr target_sip, 
	Int4 offset, Uint1 strand);

/**Extend**/
extern SeqAlignPtr SeqAlignEndExtend (SeqAlignPtr sap, Int4 start1, Int4 start2, 
	Int4 stop1, Int4 stop2, Int4 x1, Int4 y1, Int4 x2, Int4 y2, 
	Uint1 strand1, Uint1 strand2);


/**Delete, Truncate**/
extern SeqAlignPtr SeqAlignDeleteByLoc (SeqLocPtr slp, SeqAlignPtr salp);

extern SeqAlignPtr SeqAlignTrunc (SeqAlignPtr salp, Int4 from, Int4 to);
extern SeqAlignPtr SeqAlignTrunc2 (SeqAlignPtr salp, Int4 from, Int4 to);

extern SeqAlignPtr SeqAlignMapOnFirstSeq (SeqAlignPtr salp);
	
/**Clean**/
extern SeqAlignPtr CleanStrandsSeqAlign (SeqAlignPtr salp);
	
/**
*Functions taking a SeqAlign and returning a new SeqAlign
**/
/**Duplicate**/
extern SeqAlignPtr SeqAlignDup (SeqAlignPtr salp);

extern DenseDiagPtr DenseDiagDup (DenseDiagPtr ddp);

/**Transfer of Format**/
extern SeqAlignPtr  DenseSegToDenseDiag (SeqAlignPtr salp);

/**Nuc/AA**/
extern SeqAlignPtr aaSeqAlign_to_dnaSeqAlign (SeqAlignPtr salp, ValNodePtr vnp, 
	ValNodePtr framep);
	
extern SeqAnnotPtr aaSeqAnnot_to_dnaSeqAnnot (SeqAnnotPtr sap, ValNodePtr vnp, 
	ValNodePtr framep);

/**misc**/
extern SeqAnnotPtr SeqAnnotForSeqAlign (SeqAlignPtr salp);

/**
*SeqAlign and SeqEntry
**/
extern void ReplaceSeqAlignInSeqEntry (Uint2 entityID, Uint2 itemID, 
	SeqAlignPtr salp);

/*******************************************************/

extern SeqAnnotPtr CompSeqAnnotFree (SeqAnnotPtr sap);
extern SeqAlignPtr CompSeqAlignFree (SeqAlignPtr salp);
extern SeqAnnotPtr SeqAlignBoolSegCpy (SeqAnnotPtr sap, Int4 from, Int4 to);
extern SeqAlignPtr SeqAlignDenseSegToBoolSeg (SeqAlignPtr salp);
extern SeqAnnotPtr SeqAnnotDenseSegToBoolSeg (SeqAnnotPtr sap);
extern SeqAlignPtr SeqAlignBoolSegToDenseSeg (SeqAlignPtr salp);
extern SeqAnnotPtr SeqAnnotBoolSegToDenseSeg (SeqAnnotPtr sap);
extern void CompSeqAlignPrint (SeqAlignPtr salp);

extern SeqAlignPtr SeqAlignDupRegion (SeqAlignPtr salp, Int2 to_numseg, Int4 subseg, Boolean first_part);
extern SeqAlignPtr SeqAlignDupAdd (SeqAlignPtr *salp_head, SeqAlignPtr salp, Int2 to_numseg, Int4 subseg, Boolean first_part);

extern SeqAlignPtr SeqAlignExtend (SeqAlignPtr salp1, SeqAlignPtr salp2);
extern SeqAlignPtr DeleteRegion (SeqIntPtr sip, SeqAlignPtr salp);

extern SeqAlignPtr DenseDiagToDenseSegFunc (SeqAlignPtr salp, Boolean add_ends);
extern SeqAlignPtr DenseDiagToDenseSeg (SeqAlignPtr salp, Boolean add_ends);

#endif
