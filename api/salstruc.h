/*   salstruc.h
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
* File Name:  salstruc.h
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.2 $
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

#ifndef _SALSTRUC_
#define _SALSTRUC_

#include <salsa.h>
#include <txalign.h>

/**********************************************************************/
extern SelStructPtr     BufferFree (SelStructPtr ssp);
extern EditAlignDataPtr SetupDataBuffer (EditAlignDataPtr adp);
extern EditAlignDataPtr SetupDataPanel (EditAlignDataPtr adp);

extern CharPtr      next_notemptyline (ValNodePtr anp_list, ValNodePtr linebuff, Int2 numberalignline, Int2 *index, Int4 start, Int4 *drw_width, TextAlignBufPtr *tdp, AlignNodePtr *anp);

extern SelEdStructPtr is_feature_to_buffer (ValNodePtr vnphead, Uint2 bspitemID, Uint2 entityID, Int4 from, Int4 drw_width, SeqAlignPtr salp, Uint2 seqedit, ValNodePtr sqloc_list);
extern ByteStorePtr cds_to_pept (SeqLocPtr slp, Uint1 frame, Int2 gencode, Boolean include_stop);
extern void         data_collect_arrange (EditAlignDataPtr adp, Boolean recollect);
 
extern CharPtr      get_master (ValNodePtr linebuff, Uint2 entityID, Uint2 itemID, Uint2 itemtype);

extern Boolean      read_buffer_fromalignnode (EditAlignDataPtr adp, ValNodePtr *linebuff, Int4 bufferstart, Int4 minbufferlength, Int2 *numberalignline);

/**********************************
***
*** BioseqTrimN 
***   truncates the nnn's at the extremities of a bioseq bsp
***   providing TopSeqEntry sep allows to modify the SeqAlign if any
***
***********************************/
extern Boolean BioseqTrimN (BioseqPtr bsp, SeqEntryPtr sep);

/*******************************************************************
*** AddFeatFunc
******************************************************************/
extern ValNodePtr AddFeatFunc (SelEdStructPtr feat, ValNodePtr *featlist, Uint2 itemsubtype);
/******************************************************************
*
***    CollectFeatureForEditor
***      slp: the target Seq-loc
***      anp: the AlignNode belong to the target Seq-loc
***      csop: the option for gathering the features
******************************************************************/
extern ValNodePtr CollectFeatureForEditor (SeqLocPtr slp, ValNodePtr seqfeat, Uint2 seq_entityID, Uint2 bsp_itemID, Uint1 *featOrder, Boolean all_feat);

/**********************************************************************
***  Display alignment in several formats
*** 
**********************************************************************/
extern void ShowAlignmentText (FILE *fout, EditAlignDataPtr adp, SelStructPtr ssp, Int2 leftmargin, Int4 printfrom, Int4 printto, Boolean html_option);

extern void showfastagap_fromalign (SeqAlignPtr salp, Int4 line, FILE *f);

extern SeqAnnotPtr multseqalign_from_pairseqalign (SeqAlignPtr salp);
extern SeqAlignPtr PairSeqAlign2MultiSeqAlign (SeqAlignPtr salp);
extern Int2 LIBCALLBACK MultSeqAlignFromPairSeqAlign (Pointer data);



#endif
