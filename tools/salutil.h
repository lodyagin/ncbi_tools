/*   salutil.h
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
* File Name:  salutil.h
*
* Author:  Colombe Chappey, H. Sicotte
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.32 $
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
#ifndef _SALUTIL_
#define _SALUTIL_

#include <salsa.h>
#include <txalign.h>

/********************************************************
*** BandAlign Data Structure 
*********************************************************/
typedef struct mash {
   Uint1 band_method;
   Int2 reward;
   Int2 penalty;
   Int4 gap_open;
   Int4 gap_extend;
   Int4 wordsize;
   Int4 filter;
   Int4 gap_x_dropoff;
   Int4 gap_x_dropoff_final;
   Boolean is_prot;
   Boolean translate_prot;
   SeqAlignPtr transalp;
   Boolean use_gapped_blast;
   CharPtr matrixname;
   Int4 lg1_ext;
   Int4 rg1_ext;
   Int4 lg2_ext;
   Int4 rg2_ext;
   Int4 lg1_open; 
   Int4 lg2_open; 
   Int4 rg1_open; 
   Int4 rg2_open;
   Int4    blast_threshold;
   Uint1   choice_blastfilter;
   Boolean splicing;
   Boolean multidim;
   Boolean map_align;
} Mash, PNTR MashPtr;

/*********************************************************
***
***   Read SeqPort.bsp from SeqPort.start to stop
***   in : spp, from + to in seq coordinates 
***   out: length of buffer + buffer
***
**********************************************************/
extern Int4 ReadBufferFromSep (SeqPortPtr spp, CharPtr buffer, Int4 from, Int4 to, Int4 buffsegstart);
extern CharPtr ReadBufferFromSap (CharPtr str, CharPtr buffer, SeqAlignPtr salp, SeqIdPtr sip, Int4 from, Int4 to, Int4 *ncar);

/*******************************************
***
***   String Proc
***
********************************************/
extern Boolean CCStrToInt (CharPtr str, Int2Ptr intval);
extern Boolean CCStrToLong (CharPtr str, Int4Ptr longval);

extern CharPtr dashedstring (CharPtr buffer, Int4 maxbufflength);
extern CharPtr emptystring (CharPtr buffer, Int4 maxbufflength);
extern Boolean not_empty_string (CharPtr str, Int4 lens);
extern Int1    getgapsfromstring (CharPtr str,Int4 from,Int4 to, BoolPtr *gapline);
extern Boolean stringhasnochar (CharPtr str, Int4 from, Int4 to);
extern Boolean stringhasnocharplus (CharPtr str);
extern CharPtr purge_string (CharPtr *str);
extern CharPtr reverse_string (CharPtr str);
extern CharPtr to_lower (CharPtr str);
extern CharPtr complement_string (CharPtr str);
extern Int4    compare_string (CharPtr str1, CharPtr str2, Int4 *bestscorep);
extern CharPtr load_seq_data (SeqIdPtr sip, Int4 from, Int4 to, Boolean is_prot, Int4 *lenp);
extern Boolean IS_protSeqLoc (SeqLocPtr slp);
extern SeqEntryPtr StringToSeqEntry (CharPtr str, SeqIdPtr sip, Int4 length_align, Uint1 mol_type);

/*******************************************
***
***   ValNode Proc
***
********************************************/
extern ValNodePtr   ValNodeFind (ValNodePtr head, Int2 start, Int2 index);
extern ValNodePtr   ValNodeFreeType (ValNodePtr *head, Int2 seqtype);

extern SeqLocPtr    seqloc2fuzzloc(SeqLocPtr slp, Boolean is_from, Boolean is_to);
extern TextAlignBufPtr TextAlignBufFind (ValNodePtr anpvnp, Uint2 entityID, Uint2 itemID, Uint2 itemtype);
extern CharPtr PNTR buf2array (ValNodePtr list, Int2 seq1, Int2 seq2);
extern AlignNodePtr AlignNodeFind (ValNodePtr anpvnp, Uint2 entityID, Uint2 itemID, Uint2 itemtype);

extern Int2         AlignNodeIndex (ValNodePtr anpvnp, Uint2 entityID, Uint2 itemID, Uint2 itemtype);
extern void         OrderFeatProc (ValNodePtr vnpanp);
extern ValNodePtr   SeqfeatlistFree (ValNodePtr feathead);
extern ValNodePtr   SeqfeatlistFree_fromID (ValNodePtr feathead, Uint2 entityID);
extern SelEdStructPtr get_feat_fromid (ValNodePtr feat_vnp, Uint2 feattype, Uint2 ei, Uint2 ii, Int4 pos, SelEdStructPtr *prec);

extern SeqLocPtr    CollectSeqLocFromAlignNode (AlignNodePtr anp);
extern Int4         GetAlignLengthFromAlignNode (AlignNodePtr anp);

extern SeqIdPtr     SeqIdFromAlignNode (ValNodePtr anp_lst, Uint2 entityID, Uint2 itemID, Uint2 itemtype);
extern Uint1 StrandFromAlignNode (ValNodePtr anp_lst, Uint2 entityID, Uint2 itemID, Uint2 itemtype);

/*********************************************************
***
***  SeqIdPtr procedures
***    AddSeqId  : create a new seqid and add at the end
***                of the list starting with sip_head
***
***    SeqIdDupList : duplicate a list of SeqIdPtr
***
**********************************************************/
extern CharPtr      matching_seqid (SeqIdPtr sip1);
extern CharPtr      check_seqid (Uint2 choice, CharPtr ptr);
extern SeqIdPtr     AddSeqId (SeqIdPtr *sip_head, SeqIdPtr sip);
extern SeqIdPtr     SeqIdDupList (SeqIdPtr id_list);
extern SeqIdPtr     SeqIdDupBestList (SeqIdPtr id_list);
extern SeqIdPtr     SeqIdListfromSeqLoc (ValNodePtr vnpslp);
extern SeqIdPtr     ValNodeSeqIdListDup (ValNodePtr id_list);
extern CharPtr PNTR SeqIdListToCharArray (SeqIdPtr id_list, Int2 n);

extern SeqIdPtr     SeqIdReplaceID (SeqIdPtr head, SeqIdPtr pre, SeqIdPtr sip, SeqIdPtr next);
extern BioseqPtr    BioseqReplaceID (BioseqPtr bsp, SeqIdPtr newsip);
extern SeqEntryPtr  SeqEntryReplaceSeqID (SeqEntryPtr source_sep, SeqIdPtr sip);
/*********************************************************
***
***  ScorePtr procedures
***
**********************************************************/
extern ScorePtr ScoreDup (ScorePtr sp);
extern ScorePtr ScoreDupAdd (ScorePtr *sp_head, ScorePtr sp);
extern ScorePtr ScoreAdd (ScorePtr *sp_head, ScorePtr sp);

/*********************************************************
***
***  SeqLocPtr procedures
***
**********************************************************/
extern Int2      chkloc (SeqIdPtr sip, Int4 position, ValNodePtr vnp, Int4 *newpos);
extern SeqLocPtr expand_seq_loc(Int4 start, Int4 stop, Uint1 strand, SeqLocPtr loc);
extern Int4      MaxLengthSeqLoc (ValNodePtr sqloc_list);
extern Boolean   SeqLocListMatch (ValNodePtr vnp1, ValNodePtr vnp2, Boolean *Fp, Boolean *Tp);

/***********************************************************************
***    SeqEntryToSeqLoc
***      read SeqEntry (1->Bioseq or 2->BioseqSetPtr)
***      return list of ValNodePtr->SeqLocPtr
************************************************************************/
extern ValNodePtr SeqEntryToSeqLoc (SeqEntryPtr sep, Int2 *n, Uint1 bsp_mol);

/***********************************************************************
***     SeqLocToSeqId
***        read a list of SeqIdPtr, check if each SIP is NOT already open
***        built a list of SIP and open SeqPort on the SIP
***********************************************************************/
extern SeqIdPtr SeqLocToSeqId (ValNodePtr sqloc_list);

/*********************************************************
***
***   SelStruct Procedures
***
**********************************************************/
extern SelStructPtr SelStructNew (Uint2 ssp_ed, Uint2 ssp_id, Uint2 ssp_it, Int4 from, Int4 to, SeqIdPtr sip, Uint1 strand, Boolean is_fuzz);
extern SelStructPtr SelStructCpy (SelStructPtr ssp, SelStructPtr ssp2);
extern SelStructPtr SelStructDup (SelStructPtr ssp);
extern SelStructPtr SelStructAdd (SelStructPtr head, SelStructPtr ssp);
extern SelStructPtr SelStructDel (SelStructPtr ssp);
extern SelStructPtr SelStructDelList (SelStructPtr ssp);
extern void setposition_tossp (SelStructPtr ssp, Int4 from, Int4 to);

extern Boolean is_samessp (SelStructPtr ssp1, SelStructPtr ssp2);
extern Boolean is_sameId (Uint2 sei, Uint2 sii, Uint2 sit, Uint2 sist, Uint2 ei, Uint2 ii, Uint2 it, Uint2 ist);
extern Boolean is_samepos (SelStructPtr ssp1, SelStructPtr ssp2);
extern ValNodePtr del_ssp_fromid (ValNodePtr headp, Uint2 itemsubtype, SelEdStructPtr target);
extern Boolean include_ssp (SeqLocPtr slp1, SeqLocPtr slp2);
extern Int4    overlapp_startssp (SeqLocPtr slp1, SeqLocPtr slp2);
extern Boolean overlapp_ssp (SeqLocPtr slp1, SeqLocPtr slp2);
extern Boolean precede_ssp (SeqLocPtr slp1, SeqLocPtr slp2);
extern Boolean succeed_ssp (SeqLocPtr slp1, SeqLocPtr slp2);
extern SelStructPtr addssp (SelStructPtr *ssphead, Uint2 choice, Pointer pt, Uint2 iID);

/*********************************************************
***
***   SelEdStruct Procedures
***
**********************************************************/
extern SelEdStructPtr new_seledstruct (Uint2 ssp_ed, Uint2 ssp_id, Uint2 ssp_it, Uint2 bspiID, Int4 from, Int4 to, SeqIdPtr sip, Uint1 strand, Boolean is_fuzz, CharPtr label, Pointer data, Int4 offset, Uint1 codonstart);
extern SelEdStructPtr new_seledstruct_fromseqloc (Uint2 ssp_ed, Uint2 ssp_id, Uint2 ssp_it, Uint2 bspiID, SeqLocPtr slp, CharPtr label, Pointer data, Int4 offset, Uint1 codonstart);
extern SeqLocPtr      sesp_to_slp (SelEdStructPtr ses, SeqAlignPtr salp, ValNodePtr sqlocs, Boolean partial);

extern SelEdStructPtr SelEdStructDup (SelEdStructPtr ssp);
extern SelEdStructPtr SelEdStructAdd (SelEdStructPtr head, SelEdStructPtr ssp);
extern SelEdStructPtr SelEdStructDel (SelEdStructPtr ssp);
extern SelEdStructPtr SelEdStructListDel (SelEdStructPtr ssp);
extern void           setposition_toses (SelEdStructPtr ssp, Int4 from, Int4 to);
extern void           set_seqnot_visible (Uint2 eID, Uint2 iID, SelEdStructPtr sesp);
extern void           set_seqvisible (Uint2 eID, Uint2 iID, SelEdStructPtr sesp);
extern Boolean        is_seqvisible (Uint2 eID, Uint2 iID, SelEdStructPtr sesp);
extern SelEdStructPtr ss_to_ses (SelStructPtr ssp);
extern SelStructPtr   ses_to_ss (SelEdStructPtr sesp);
extern Boolean        is_samess_ses (SelStructPtr ssp1, SelEdStructPtr ssp2);
extern Boolean        is_sameses (SelEdStructPtr ssp1, SelEdStructPtr ssp2);

/*********************************************************
***
***   ObjMgr Procedures
***
**********************************************************/
extern Boolean        AfterAlsoSelect (void);
extern void           ObjMgrSelectPrint (void);
extern void           SelectType (EditAlignDataPtr adp, Uint2 feattype, Int4 slpto);
extern Int2           GetNumberObjMgrSelect (void);
extern Int2           checkOMss_for_itemtype (Uint2 itemtype);
extern SelStructPtr   getOMselect_for_itemtype (Uint2 itemtype);
extern SelStructPtr   is_selectedbyID (Uint2 entityID, Uint2 itemID, Uint2 itemtype);
extern SelEdStructPtr getCDSselect (ValNodePtr seqfeathead, ValNodePtr feathead);
extern Int2           checkCDSselect_forprotein (ValNodePtr seqfeathead, ValNodePtr feathead, Boolean with_prot);
extern Boolean        checkssp_for_editor (SelStructPtr ssp);
extern SeqLocPtr      checkseqlocfeature_for_editor (Uint2 entityID, Uint2 itemID, ValNodePtr headfeat);
extern void           checkselectsequinfeature_for_editor (ValNodePtr headfeat);
extern Int4           getminpos_fromOMselect (Uint2 itemsubtype);

/*********************************************************
*** * * * * * * * * * * * * * * * * * * * * * * * * * * * 
***
***  SeqCoordToAlignCoord procedures
***  AlignCoordToSeqCoord procedures
***
*** * * * * * * * * * * * * * * * * * * * * * * * * * * * 
**********************************************************/
extern Boolean  locate_in_seqalignds  (Int4 pos, Int2 dim, Int2 dspnumseg, Int4Ptr *dspstart, Int4Ptr *dsplens, Int2 *numseg_before, Int2 *subdsplens);
extern Boolean  locate_in_seqalign (Int4 pos, Int2 dim, Int2 dspnumseg, BoolPtr *dspstart, Int4Ptr *dsplens, Int2 *numseg_before, Int2 *subdsplens, Int4 *sumdsplens_before);
extern Int4 SeqCoordToAlignCoord (Int4 position, SeqIdPtr sip, SeqAlignPtr salp, Int2 intersalpwidth, Int2 is_end);
extern Int4 AlignCoordToSeqCoord (Int4 position, SeqIdPtr sip, SeqAlignPtr salp, ValNodePtr sqloc_list, Int2 intersalpwidth);
extern Int4 AlignCoordToSeqCoord2 (Int4 position, SeqIdPtr sip, SeqAlignPtr salp,ValNodePtr sqloc_list, Int2 intersalpwidth);
extern Boolean GetAlignCoordFromSeqLoc (SeqLocPtr slp, SeqAlignPtr salp, Int4 *start, Int4 *stop);
extern Boolean  SeqPosToLineColumn (Uint2 itemID, Uint2 entityID, Uint2 itemtype, Int4 pos, Int2 *line, Int2 *column, Int4 hoffset, EditAlignDataPtr adp);
extern Boolean  SeqPosInLineColumn (Int4 pos, Int2 alignline, Int2 *column, Int4 hoffset, EditAlignDataPtr adp);
extern Boolean  LocateInSeqAlign (Int4 pos, Int2 dim, Int2 dspnumseg, BoolPtr *dspstart, Int4Ptr *dsplens, Int2 *numseg_before, Int4 *subdsplens, Int4 *sumdsplens_before);
extern Boolean LocateInSeqAlignDenSeg (Int4 pos, Int2 dim, Int2 dspnumseg, Int4Ptr *dspstart, Int4Ptr *dsplens, Int2 *numseg_before, Int4 *subdsplens);

/*********************************************************
***
***   Find Matching Pattern Procedures
***
**********************************************************/
extern SelStructPtr SetupMatchPatList (ValNodePtr match_pat, ValNodePtr anp_list);
extern SelStructPtr ShowNextPattern (SelStructPtr match_pat, SelStructPtr cur_pat, Int4 start);
extern SelStructPtr ShowPrecPattern (SelStructPtr match_pat, SelStructPtr cur_pat, Int4 start);
/*
extern ValNodePtr EditFindPattern (CharPtr str, ValNodePtr sqloc_list, Uint1 mol_type);
*/

/*********************************************************
***
***   Editing Procedures
***
**********************************************************/
extern CharPtr char_to_insert (Char *ch, Int4 lens, Uint1 mol_type);
extern Boolean insertchar (CharPtr str, Int4 pos, SeqIdPtr target, Uint1 mol_type, Boolean spliteditmode);
extern Boolean insertchar_atcaret (CharPtr str, EditAlignDataPtr adp);

extern SeqEntryPtr getfirst_sep(SeqEntryPtr sep, Uint1 bsp_mol);

/*********************************************
*** SeqLocListToSeqAlign
***    aligns the sequences from the SeqLocs list (sqloc_list) 
***    returns a SeqAlign
***    Alignment methods:  
***      1: FASTA 
***      5: BandAlign (GlobalAlignTwoSeqLocs) 
**********************************************/
extern MashPtr MashNew (Boolean is_prot);
extern SeqAlignPtr SeqLocListToSeqAlign (ValNodePtr sqloc_list, Int2 choice, Pointer param);
extern SeqLocPtr AlignmRNA2genomic (BioseqPtr bsp1, BioseqPtr bsp2);

#endif

