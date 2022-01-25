/* ===========================================================================
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
* ===========================================================================*/
/*****************************************************************************

File name: mblast.h

Author: Ilya Dondoshansky

Contents: prototypes for "public" Mega BLAST functions (ones that other utilitiles
	may safely call).

******************************************************************************/

/* $Revision: 6.33 $ 
* $Log: mblast.h,v $
* Revision 6.33  2000/12/21 22:29:29  dondosha
* Added prototype for MegaBlastGetHspPercentIdentity
*
* Revision 6.32  2000/10/31 15:06:16  dondosha
* Added Boolean parameter to PrintMaskedSequence function
*
* Revision 6.31  2000/10/24 18:56:45  dondosha
* Moved prototype of UniqueLocalId() to blast.h
*
* Revision 6.30  2000/10/19 16:02:29  dondosha
* Changed MegaBlastGetPercentIdentity to MegaBlastGetNumIdentical
*
* Revision 6.29  2000/10/05 19:47:51  dondosha
* Added prototype for MegaBlastSaveCurrentHitlist
*
* Revision 6.28  2000/08/30 22:08:45  dondosha
* Added function BioseqMegaBlastEngineByLoc
*
* Revision 6.27  2000/08/25 17:00:17  dondosha
* Added prototype for MegaBlastReevaluateWithAmbiguities
*
* Revision 6.26  2000/07/11 16:45:39  dondosha
* Removed mask SeqLocPtr array from all related prototypes
*
* Revision 6.25  2000/07/08 20:44:12  vakatov
* Get all "#include" out of the 'extern "C" { }' scope;  other cleanup...
*
* Revision 6.24  2000/06/30 17:52:45  madden
* Move AWAKE_THR_MIN_SIZE to blastdef.h
*
* Revision 6.23  2000/06/27 22:21:02  dondosha
* Added 3 prototypes for masked query output
*
* Revision 6.22  2000/06/27 14:48:04  dondosha
* Changed and added macros for testing HSP inclusion
*
* Revision 6.21  2000/05/26 19:21:39  dondosha
* Added two defines for neighboring
*
* Revision 6.20  2000/05/12 19:40:17  dondosha
* Added prototype for BinarySearchInt4; macro MB_HSP_CONTAINED
*
* Revision 6.19  2000/05/03 20:20:15  dondosha
* Added prototype for MegaBlastGappedAlign
*
* Revision 6.18  2000/04/12 21:10:05  dondosha
* Added MegaBlastGetPercentIdentity prototype
*
* Revision 6.17  2000/04/10 17:43:44  dondosha
* Changed prototype for BlastNeedHumanRepeatFiltering
*
* Revision 6.16  2000/04/07 16:52:44  dondosha
* Changed prototype for BioseqMegaBlastEngine, removed BioseqMegaBlastEngineByLoc, BlastSearchHandleResults
*
* Revision 6.15  2000/04/05 18:12:33  dondosha
* Moved SeqIdSetDup to objloc.h
*
* Revision 6.14  2000/04/04 20:51:57  dondosha
* No need for BlastAdjustHitOffsets prototype anymore
*
* Revision 6.13  2000/04/04 16:11:16  dondosha
* Added prototype for BlastNeedHumanRepeatFiltering
*
* Revision 6.12  2000/03/31 19:10:26  dondosha
* Changed some names related to MegaBlast
*
* Revision 6.11  2000/03/29 22:09:11  dondosha
* Added prototypes for gap info related functions
*
* Revision 6.10  2000/03/16 18:12:25  dondosha
* Added prototype for SeqIdSetDup
*
* Revision 6.9  2000/03/13 21:09:46  dondosha
* Added prototype for BlastSortUniqHspArray
*
* Revision 6.8  2000/03/08 20:50:17  madden
* Remove prototype for BlastGetAllowedGis
*
* Revision 6.7  2000/03/03 18:08:52  dondosha
* Added prototype for MegaBlastWordFinderDeallocate
*
* Revision 6.6  2000/03/02 17:22:26  dondosha
* Added array of lower case mask SeqLocPtrs as parameter to several routines
*
* Revision 6.5  2000/02/24 17:49:17  dondosha
* Added prototype for BlastAdjustHitOffsets
*
* Revision 6.4  2000/02/11 20:55:38  dondosha
* Added prototypes for new word finder and extension routines
*
* Revision 6.3  2000/02/03 22:12:26  dondosha
* Added header comments
*
* Revision 6.2  2000/02/02 15:04:14  dondosha
* Removed unused routine ReapHitlistByContext
*
* */

#ifndef __MBLAST__
#define __MBLAST__

#include <ncbi.h>
#include <blastpri.h>
#include <lookup.h>
#include <objcode.h>
#include <objseq.h>
#include <sequtil.h>
#include <tofasta.h>
#include <seqport.h>
#include <readdb.h>
#include <ncbithr.h>
#include <gapxdrop.h>
#include <dust.h>
#include <mbutils.h>
#include <mbalign.h>


#ifdef __cplusplus
extern "C" {
#endif

#define MB_DIAG_CLOSE 10
#define MB_DIAG_NEAR 30

#define MB_HSP_CLOSE(q1, q2, s1, s2, c) \
(ABS((q1-s1) - (q2-s2)) < c)

#define MB_HSP_CONTAINED(qo1,qo2,qe2,so1,so2,se2,c) \
(qo1>=qo2 && qo1<=qe2 && so1>=so2 && so1<=se2 && \
MB_HSP_CLOSE(qo1,qo2,so1,so2,c))

#define MIN_NEIGHBOR_PERC_IDENTITY 96
#define MIN_NEIGHBOR_HSP_LENGTH 100

Int2 MegaBlastSetUpSearchInternalByLoc PROTO((BlastSearchBlkPtr search,
					      SeqLocPtr query_slp, BioseqPtr
					      query_bsp, CharPtr prog_name, Int4
					      qlen, BLAST_OptionsBlkPtr options,
					      int (LIBCALLBACK
						   *callback)PROTO((Int4 done,
								    Int4
								    positives))));

SeqAlignPtr PNTR 
BioseqMegaBlastEngine PROTO((BioseqPtr PNTR bspp, CharPtr progname, CharPtr
			     database, BLAST_OptionsBlkPtr options, ValNodePtr
			     *other_returns, ValNodePtr *error_returns, int
			     (LIBCALLBACK *callback)PROTO((Int4 done, Int4
							   positives)), 
			     SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, 
			     Int4 gi_list_total, 
			     int (LIBCALLBACK *results_callback)PROTO((VoidPtr Ptr))));
SeqAlignPtr PNTR
BioseqMegaBlastEngineByLoc PROTO((SeqLocPtr slp, CharPtr progname, CharPtr database,
                            BLAST_OptionsBlkPtr options, ValNodePtr *other_returns,
                            ValNodePtr *error_returns, 
                            int (LIBCALLBACK *callback)(Int4 done, Int4 positives),
                            SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, 
                            Int4 gi_list_total, 
                            int (LIBCALLBACK *results_callback)PROTO((VoidPtr Ptr))));

SeqAlignPtr PNTR
BioseqMegaBlastEngineCore PROTO((BlastSearchBlkPtr search, BLAST_OptionsBlkPtr options, Int4Ptr *pos_matrix));

BlastSearchBlkPtr BlastFillQueryOffsets PROTO((BlastSearchBlkPtr search, SeqLocPtr
					query_slp));

Boolean MegaBlastGetFirstAndLastContext PROTO((CharPtr prog_name, SeqLocPtr query_slp,
					Int2Ptr first_context, Int2Ptr
					last_context, Uint1 strand_options));

Int4 SeqLocTotalLen PROTO((CharPtr prog_name, SeqLocPtr slp));

BlastSearchBlkPtr MegaBlastSetUpSearchWithReadDbInternal PROTO((SeqLocPtr
								query_slp,
								BioseqPtr
								query_bsp,
								CharPtr
								prog_name, Int4
								qlen, CharPtr
								dbname,
								BLAST_OptionsBlkPtr options, int (LIBCALLBACK *callback)PROTO((Int4 done, Int4 positives)), SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, Int4 gi_list_total, ReadDBFILEPtr rdfp));

SeqAlignPtr PNTR 
MegaBlastPackAlignmentsByQuery PROTO((BlastSearchBlkPtr search,
						SeqAlignPtr seqalign));
Int2 LIBCALL
MegaBlastSequenceAddSequence PROTO((BlastSequenceBlkPtr sequence_blk, Uint1Ptr sequence, Uint1Ptr sequence_start, Int4 length, Int4 original_length, Int4 effective_length));

BlastDoubleInt4Ptr GetGisFromFile PROTO((CharPtr file_name, Int4Ptr gi_list_size));
int LIBCALLBACK compare PROTO((VoidPtr v1, VoidPtr v2));
Uint1Ptr GetPrivatTranslationTable PROTO((CharPtr genetic_code, Boolean
					  reverse_complement));

BioseqPtr BlastMakeTempProteinBioseq PROTO((Uint1Ptr sequence, Int4 length, Uint1
				      alphabet));

VoidPtr index_proc PROTO((VoidPtr dummy));

CharPtr BlastConstructFilterString PROTO((Int4 filter_value));

int LIBCALLBACK evalue_compare_hits PROTO((VoidPtr v1, VoidPtr v2));

Int2 blast_set_parameters PROTO((BlastSearchBlkPtr search, Int4
				 dropoff_number_of_bits_1st_pass, Int4
				 dropoff_number_of_bits_2nd_pass, Nlm_FloatHi
				 avglen, Nlm_FloatHi searchsp, Int4 window));
void HackSeqLocId PROTO((SeqLocPtr slp, SeqIdPtr id));

Uint1 FrameToDefine PROTO((Int2 frame));

Boolean 
BlastGetFirstAndLastContext PROTO((CharPtr prog_name, SeqLocPtr query_slp, Int2Ptr first_context, Int2Ptr last_context, Uint1 strand_options));

Int4
MegaBlastExtendHit PROTO((BlastSearchBlkPtr search, LookupTablePtr lookup, 
		   Int4 s_off, Int4 q_off));
Int2
MegaBlastNtWordExtend PROTO((BlastSearchBlkPtr search, Uint1Ptr subject0, 
			     Int4 q_off, Int4 s_off));

Int4 
MegaBlastWordFinder PROTO((BlastSearchBlkPtr search, LookupTablePtr lookup));


   /*void BlastAdjustHitOffsets PROTO((BlastSearchBlkPtr search));*/

void MegaBlastMaskTheResidues PROTO((Uint1Ptr buffer, Int4 max_length, Uint1
			      mask_residue, SeqLocPtr mask_slp, Boolean reverse,
			      Int4 offset, Boolean lowercase_mask));
BLAST_WordFinderPtr
MegaBlastWordFinderDeallocate PROTO((BLAST_WordFinderPtr wfp));

void
BlastSortUniqHspArray PROTO((BLAST_HitListPtr hitlist));

void
MegaBlastFillHspGapInfo PROTO((BLAST_HSPPtr hsp, edit_script_t PNTR ed_script));

GapXEditScriptPtr
MBToGapXEditScript PROTO((edit_script_t PNTR ed_script));

SeqAlignPtr
MegaBlastSeqAlignFromResultHitlist PROTO((BlastSearchBlkPtr search,
				   BLASTResultHitlistPtr result_hitlist,
				   SeqIdPtr subject_id));
Boolean 
BlastNeedHumanRepeatFiltering PROTO((BlastDoubleInt4Ptr gi_list, 
				     Int4 gi_list_size, Uint4 query_gi));

Int4
MegaBlastGetNumIdentical PROTO((Uint1Ptr query, Uint1Ptr subject, 
                                Int4 q_start, Int4 s_start, Int4 length, 
                                Boolean reverse));

Int2
MegaBlastReevaluateWithAmbiguities PROTO((BlastSearchBlkPtr search, Int4 sequence_number));

Int4 MegaBlastGappedAlign PROTO((BlastSearchBlkPtr search));

Int4 BinarySearchInt4 PROTO((Int4 n, Int4Ptr A, Int4 size));

void BlastLCaseMaskTheResidues PROTO((Uint1Ptr buffer, Int4 max_length,
				      SeqLocPtr mask_slp, Boolean reverse, 
				      Int4 offset));

SeqLocPtr MaskSeqLocFromSeqAlign PROTO((SeqAlignPtr seqalign));
   
void PrintMaskedSequence PROTO((BioseqPtr query_bsp, SeqLocPtr mask_slp, 
				CharPtr file_name, Boolean first));
Int4 LIBCALL 
MegaBlastSaveCurrentHitlist PROTO((BlastSearchBlkPtr search));

FloatLo 
MegaBlastGetHspPercentIdentity PROTO((BlastSearchBlkPtr search, BLAST_HSPPtr hsp));

#ifdef __cplusplus
}
#endif
#endif /* !__MBLAST__ */
