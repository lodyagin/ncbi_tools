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

File name: blastpri.h

Author: Tom Madden

Contents: prototypes for "private" BLAST functions, these should not be called
	by outside utilities. 

******************************************************************************/

/* $Revision: 6.40 $ */
/* $Log: blastpri.h,v $
/* Revision 6.40  1999/01/19 13:26:47  madden
/* Change to HspArrayPurge
/*
 * Revision 6.39  1999/01/08 22:08:43  madden
 * BlastScaleMatrix returns factor as FloatHi
 *
 * Revision 6.38  1998/12/18 16:19:58  madden
 * Make BLASTSetUpSearchWithReadDbInternal public, add BlastSearchBlkNewExtra
 *
 * Revision 6.37  1998/09/17 19:53:03  madden
 * Added fillCandLambda
 *
 * Revision 6.36  1998/09/16 18:59:55  madden
 * Added subset Boolean
 *
 * Revision 6.35  1998/09/14 15:11:16  egorov
 * Add support for Int8 length databases; remove unused variables
 *
 * Revision 6.34  1998/09/10 22:36:10  madden
 * Added convertSeqAlignListToValNodeList and convertValNodeListToSeqAlignList
 *
 * Revision 6.33  1998/09/09 21:18:10  madden
 * Added PrintKAParametersExtra
 *
 * Revision 6.32  1998/09/04 14:45:43  madden
 * Moved code from blast.c blastool.c
 *
 * Revision 6.31  1998/07/21 20:58:07  madden
 * Changes to allow masking at hash only
 *
 * Revision 6.30  1998/06/12 16:08:50  madden
 * BlastHitRange stuff
 *
 * Revision 6.29  1998/06/04 16:23:04  madden
 * BioseqSeg to MyBioseqSeg
 *
 * Revision 6.28  1998/05/17 16:28:44  madden
 * Allow changes to filter options and cc filtering.
 *
 * Revision 6.27  1998/05/05 14:05:36  madden
 * Added functions BlastStartAwakeThread and BlastStopAwakeThread
 *
 * Revision 6.26  1998/04/24 19:28:48  madden
 * Added BlastScaleMatrix (and other rescaling code moved from posit.c)
 *
 * Revision 6.25  1998/04/15 20:23:49  madden
 * offset arg removed from BlastMaskTheResidues
 *
 * Revision 6.24  1998/03/26 14:21:35  madden
 * Added GetScoreSetFromBlastResultHsp prototype
 *
 * Revision 6.23  1998/03/25 22:27:27  madden
 * Remove GetScoreSetFromBlastResultHsp prototype
 *
 * Revision 6.22  1998/03/24 15:38:24  madden
 * Use BlastDoubleInt4Ptr to keep track of gis and ordinal_ids
 *
 * Revision 6.21  1998/03/18 14:14:16  madden
 * Support random access by gi list
 *
 * Revision 6.20  1998/03/16 14:02:13  madden
 * Changed call to BlastSeqIdListDestruct
 *
 * Revision 6.19  1998/02/27 16:52:07  madden
 * Added BlastGetSequenceFromBioseq
 *
 * Revision 6.18  1998/02/27 14:34:51  madden
 * Added error message prototypes
 *
 * Revision 6.17  1998/02/26 22:34:25  madden
 * Changes for 16 bit windows
 *
 * Revision 6.16  1998/02/26 19:11:38  madden
 * Added prototypes for BlastNtFindWords BlastPopulateAllWordArrays
 *
 * Revision 6.15  1998/02/19 17:17:12  madden
 * Use of Int4 rather than Int2 when pruning SeqAlign
 *
 * Revision 6.14  1998/02/11 17:18:17  madden
 * Made BlastGetGappedAlignmentTraceback functions to BlastGetGapAlgnTbck (shorter than 32 chars)
 *
 * Revision 6.13  1998/01/05 16:46:54  madden
 * One or both strands can be searched, as opposed to only both, changes to number of contexts
 *
 * Revision 6.12  1997/12/31 17:53:11  madden
 * Removed BLAST_WordFinderNew and BLAST_WordFinderDestruct prototypes
 *
 * Revision 6.11  1997/12/23 18:12:38  madden
 * Changes for range-dependent blast
 *
 * Revision 6.10  1997/12/12 20:38:36  madden
 * ContextToFrame lost last parameter
 *
 * Revision 6.9  1997/11/28 18:19:37  madden
 * Changes to TxDfDbInfoNew
 *
 * Revision 6.8  1997/11/07 00:48:29  madden
 * Added TXMATRIX defintion
 *
 * Revision 6.7  1997/10/24 20:46:55  madden
 * Removed BLASTResultsStructNew prototype
 *
 * Revision 6.6  1997/10/24 19:09:19  madden
 * Removed BlastSetReadDB and BlastGetReadDB_ID, changed to ReadDBGetDb and ReadDBGetDbId
 *
 * Revision 6.5  1997/10/03 21:27:34  madden
 * Added BlastGetTypes
 *
 * Revision 6.4  1997/10/02 17:29:27  madden
 * Added PrintDbInformationBasic
 *
 * Revision 6.3  1997/09/18 22:22:09  madden
 * Added prune functions
 *
 * Revision 6.2  1997/09/16 16:31:33  madden
 * More changes for multiple db runs
 *
 * Revision 6.1  1997/09/11 18:49:28  madden
 * Changes to enable searches against multiple databases.
 *
 * Revision 6.0  1997/08/25 18:52:44  madden
 * Revision changed to 6.0
 *
 * Revision 1.48  1997/08/22 18:37:49  madden
 * Added function BlastOtherReturnsPrepare
 *
 * Revision 1.47  1997/07/24 20:34:50  madden
 * define change for masking
 *
 * Revision 1.46  1997/07/18 14:26:42  madden
 * call to AcknowledgeBlastQuery changed
 *
 * Revision 1.45  1997/07/17 20:27:11  madden
 * Changed defines to indicate frame
 *
 * Revision 1.44  1997/07/16 20:34:48  madden
 * Added function BlastConvertProteinSeqLoc
 *
 * Revision 1.43  1997/07/15 20:36:11  madden
 * Added BioseqSeg and SeqLocSeg
 *
 * Revision 1.42  1997/07/14 15:32:25  madden
 * prototype for BlastConstructErrorMessage
 *
 * Revision 1.41  1997/06/27 14:30:42  madden
 * prototypes for BlastAddSeqIdToList and BlastSeqIdListDestruct
 *
 * Revision 1.40  1997/06/06 21:29:36  madden
 * Added Boolean html to AcknowledgeBlastQuery and PrintDbInformation
 *
 * Revision 1.39  1997/06/06 19:50:58  madden
 * Added BlastMakeFakeBioseq and BlastDeleteFakeBioseq
 *
 * Revision 1.38  1997/05/27 20:20:08  madden
 * Added function BlastMaskTheResidues
 *
 * Revision 1.37  1997/04/23 21:56:07  madden
 * Changes in BlastGetGappedAlignmentTraceback for in-frame gapping tblastn.
 *
 * Revision 1.36  1997/04/11  21:18:45  madden
 * Added GetSequenceWithDenseSeg.
 *
 * Revision 1.35  1997/04/07  18:17:09  madden
 * Added prototype for BioseqBlastEngineCore
 *
 * Revision 1.34  1997/03/06  21:47:10  madden
 * Added FormatBlastParameters.
 *
 * Revision 1.33  1997/03/05  14:29:46  madden
 * Added prototype for BlastSaveCurrentHsp.
 *
 * Revision 1.32  1997/03/04  21:34:59  madden
 * Added in HspArrayPurge.
 *
 * Revision 1.31  1997/03/04  20:36:51  madden
 * *** empty log message ***
 *
 * Revision 1.30  1997/03/03  22:39:45  madden
 * Moved code from blast.c to blastutl.c.
 *
 * Revision 1.29  1997/03/03  21:48:52  madden
 * *** empty log message ***
 *
 * Revision 1.28  1997/03/01  18:25:33  madden
 * reverse flag added to BlastGetGappedAlignmentTraceback functions.
 *
 * Revision 1.27  1997/02/26  23:39:54  madden
 * Added Txdfline stuff.
 *
 * Revision 1.26  1997/02/12  22:19:08  madden
 * Added prototype for BlastNewFindWords.
 *
 * Revision 1.25  1997/02/11  19:30:54  madden
 * Added prototypes for BlastGetGappedScoreWithReaddb and BlastGetGapped
 *
 * Revision 1.24  1997/02/10  20:03:58  madden
 * Added specific to BlastAllWordNew.
 *
 * Revision 1.23  1997/02/07  22:43:03  madden
 * Moved BLAST_WordFinderNew and Destruct from blast.c to blastutl.c, made
 * non-static.
 *
 * Revision 1.22  1997/02/07  22:32:40  madden
 * Changed prototypes for BlastGetSubjectId and GetSeqAlignForResultHitList.
 *
 * Revision 1.21  1997/01/30  21:41:17  madden
 * Prototype for FormatBlastParameters added.
 *
 * Revision 1.20  1997/01/11  18:58:29  madden
 * Removed defunct PerformBlastSearch... functions.
 *
 * Revision 1.19  1997/01/07  20:40:29  madden
 * Added reverse Boolean to GetSeqAlignForResultHitList.
 *
 * Revision 1.18  1997/01/06  22:41:46  madden
 * Added prototype BlastGetSubjectId.
 *
 * Revision 1.17  1996/12/23  22:02:05  madden
 * Changes to allow two sequences to be compared.
 *
 * Revision 1.16  1996/12/20  15:31:39  madden
 * Removed prototype for Perform2PassBlastSearchWithReadDb.
 *
 * Revision 1.15  1996/12/20  14:22:48  madden
 * Added discontinuous Boolean to GetSeqAlignForResultHitList.
 *
 * Revision 1.14  1996/12/12  16:46:25  madden
 * Changed CONTAINED_IN_HSP.
 *
 * Revision 1.13  1996/12/08  15:19:59  madden
 * Added defines and prototypes for gapped alignments.
 *
 * Revision 1.12  1996/11/14  16:21:55  madden
 * changed CharPtr to Uint1Ptr in GetTranslation.
 *
 * Revision 1.11  1996/11/13  22:35:18  madden
 * Added prototype for GetTranslation.
 *
 * Revision 1.10  1996/11/05  23:19:08  madden
 * Changed BlastTranslateUnambiguousSequence prototype.
 *
 * Revision 1.9  1996/09/26  20:18:43  madden
 * Changed prototype for GetSeqAlignForResultHitList
 *
 * Revision 1.8  1996/09/12  21:12:23  madden
 * Removed prototypes for BLAST_WordFinderNew and BLAST_WordFinderDestruct.
 *
 * Revision 1.7  1996/09/11  22:22:12  madden
 * Added prototpe for BLASTPerformSearchWithReadDb.
 *
 * Revision 1.6  1996/08/26  17:24:19  shavirin
 * Added definition of function Win32TimeFill()
 *
 * Revision 1.5  1996/08/23  16:30:11  shavirin
 * Fixed NT compiler warnings type mismatch
 *
 * Revision 1.4  1996/08/15  18:58:49  madden
 * Changed context from Int2 to Int1
 *
 * Revision 1.3  1996/08/14  15:20:37  madden
 * Added prototype for BlastTranslateUnambiguousSequence.
 *
 * Revision 1.2  1996/08/07  14:24:15  madden
 * Removed functions that depend on BLAST0 structures.
 *
 * Revision 1.1  1996/08/05  19:46:53  madden
 * Initial revision
 *
 * Revision 1.32  1996/07/31  13:10:53  madden
 * Added BlastSearchBlkDuplicate prototype.
 *
 * Revision 1.31  1996/07/25  20:47:49  madden
 * Change to arguments of Perform2PassBlastSearchWithReadDb.
 *
 * Revision 1.30  1996/07/18  22:01:35  madden
 * Changed call to BlastFindWords
 *
 * Revision 1.29  1996/06/21  15:15:21  madden
 * Removed unused prototype
 *
 * Revision 1.28  1996/06/20  16:15:57  madden
 * Replaced int's with Int4's.
 *
 * Revision 1.27  1996/06/17  21:07:32  madden
 * Added prototype for GetSeqAlignForResultHitList and some LIBCALL's.
 *
 * Revision 1.26  1996/06/17  19:02:35  madden
 * Removed unused prototypes.
 *
 * Revision 1.25  1996/06/13  21:16:33  madden
 * removed BLAST_ExtendWordNew prototype.
 *
 * Revision 1.24  1996/06/11  17:58:31  madden
 * removed prototype for BLAST_ExtendWordDiagResize.
 *
 * Revision 1.23  1996/06/06  14:09:22  madden
 * Removed blast_set_parameters prototype (it became static).
 *
 * Revision 1.22  1996/06/06  13:54:51  madden
 * Removed defunct function BLAST_ParameterBlkFill
 *
 * Revision 1.21  1996/06/04  15:33:55  madden
 * Changed BlastHitList function prototypes.
 *
 * Revision 1.20  1996/05/29  12:44:40  madden
 * Added prototype for BlastTimeFillStructure.
 *
 * Revision 1.19  1996/05/16  19:51:09  madden
 * Added documentation block.
 *
 * Revision 1.18  1996/05/14  16:15:59  madden
 * Added protoytpe for BLASTResultsStruc functions.
 *
 * Revision 1.17  1996/05/01  14:59:41  madden
 * *** empty log message ***
 *
 * Revision 1.16  1996/04/03  19:14:08  madden
 * added functions PerformBlastSearchWithReadDb and Perform2PassBlastSearchWithReadDb.
 *
 * Revision 1.15  1996/03/29  21:28:20  madden
 * *** empty log message ***
 *
 * Revision 1.14  1996/03/29  14:09:37  madden
 * prototype for GetSeqAlignForSparseHitList added.
 *
 * Revision 1.13  1996/02/28  21:38:36  madden
 * Changed prototypes for discontiguous words.
 *
 * Revision 1.12  1996/02/05  18:46:57  madden
 * *** empty log message ***
 *
 * Revision 1.11  1996/02/02  19:25:32  madden
 * Changed BlastFindWords prototype.
 *
 * Revision 1.10  1996/01/31  22:28:46  madden
 * Added prototype for BlastReapHitlistByEvalue.
 *
 * Revision 1.9  1996/01/17  13:47:13  madden
 * *** empty log message ***
 *
 * Revision 1.8  1996/01/11  15:17:58  madden
 * Added prototype for do_MPblast_search.
 *
 * Revision 1.7  1996/01/10  17:51:09  madden
 * Added SortHitListByPvalue.
 *
 * Revision 1.6  1996/01/06  18:58:45  madden
 * Added prototype for BlastLinkHsps.
 *
 * Revision 1.5  1995/12/30  18:39:27  madden
 * Added prototype for GetBLAST0KABlk.
 *
 * Revision 1.4  1995/12/28  21:26:30  madden
 * Added in prototype for do_blast_search.
 *
 * Revision 1.3  1995/12/26  23:05:29  madden
 * Added prototype for blast_set_parameters.
 *
 * Revision 1.2  1995/12/21  23:11:11  madden
 *  BLAST_Score prototypes moved to blastkar.h.
 *
 * Revision 1.1  1995/12/19  22:31:17  madden
 * Initial revision
 *
 * Revision 1.1  1995/12/08  15:48:23  madden
 * Initial revision
 *
 * */
#ifndef __BLASTPRI__
#define __BLASTPRI__
#ifdef __cplusplus
extern "C" {
#endif

#include <blast.h>
#include <blastkar.h>
#include <posit.h>
#include <seed.h>
#include <ffprint.h>


typedef struct _txdbinfo {
   struct _txdbinfo PNTR next;
   Boolean   is_protein;
   CharPtr   name;
   CharPtr   definition;
   CharPtr   date;
   Int8   total_length;
   Int4   number_seqs;
   Boolean subset;	/* Print the subset message. */
} TxDfDbInfo, *TxDfDbInfoPtr;
	

/*
	Defines for the return values in "other_returns".
*/

#define SEQLOC_MASKING_NOTSET 0
#define SEQLOC_MASKING_PLUS1 1
#define SEQLOC_MASKING_PLUS2 2
#define SEQLOC_MASKING_PLUS3 3
#define SEQLOC_MASKING_MINUS1 4
#define SEQLOC_MASKING_MINUS2 5
#define SEQLOC_MASKING_MINUS3 6
#define TXDBINFO 10
#define TXKABLK_NOGAP 12
#define TXKABLK_GAP 13
#define TXPARAMETERS 14
#define TXMATRIX 15

	
/*
	Allocates memory for TxDfDbInfoPtr.
	Link up new (returned) value to 'old', if non-NULL.
*/
TxDfDbInfoPtr LIBCALL TxDfDbInfoNew PROTO((TxDfDbInfoPtr old));

/*
	Deallocates memory (including strings for name, definition, and date).
*/
TxDfDbInfoPtr LIBCALL TxDfDbInfoDestruct PROTO((TxDfDbInfoPtr dbinfo));

/*
Print a summary of the query.
*/
Boolean LIBCALL AcknowledgeBlastQuery PROTO((BioseqPtr bsp, Int4 line_length, FILE *outfp, Boolean believe_query, Boolean html));

/*
	Print a report of the database used.
*/
Boolean LIBCALL PrintDbReport PROTO((TxDfDbInfoPtr dbinfo, Int4 line_length, FILE *outfp));

/*
	print out some of the Karlin-Altschul parameters.
*/
Boolean LIBCALL PrintKAParameters PROTO((Nlm_FloatHi Lambda, Nlm_FloatHi K, Nlm_FloatHi H, Int4 line_length, FILE *outfp, Boolean gapped));
Boolean LIBCALL PrintKAParametersExtra PROTO((Nlm_FloatHi Lambda, Nlm_FloatHi K, Nlm_FloatHi H, Nlm_FloatHi C, Int4 line_length, FILE *outfp, Boolean gapped));

/*
	Print a CharPtr (VisibleString), printing a new line every time
	a tilde is encountered.
*/
Boolean LIBCALL PrintTildeSepLines PROTO((CharPtr buffer, Int4 line_length, FILE *outfp));

/* How many interations should be done in the bisection. */
#define BLAST_SAVE_ITER_MAX 20


/*
	TRUE if c is between a and b; f between d and f.  Determines if the
	coordinates are already in an HSP that has been evaluated. 
*/
#define CONTAINED_IN_HSP(a,b,c,d,e,f) (((a <= c && b >= c) && (d <= f && e >= f)) ? TRUE : FALSE)

Int2 LIBCALL BlastFindWords PROTO((BlastSearchBlkPtr search, Int4 start, Int4 len, BLAST_Score threshold, Int1 context_index));

/*AAS*/
Int2 LIBCALL BlastNewFindWords PROTO((BlastSearchBlkPtr search, Int4 start, Int4 len, BLAST_Score threshold, Int1 context_index));

Int2 LIBCALL BlastLinkHsps PROTO ((BlastSearchBlkPtr search));

Int2 LIBCALL BlastReapHitlistByEvalue PROTO ((BlastSearchBlkPtr search));

Int2 LIBCALL BlastSaveCurrentHitlist PROTO((BlastSearchBlkPtr search));

BlastSearchBlkPtr LIBCALL BLASTPerformSearchWithReadDb PROTO((BlastSearchBlkPtr search, Int4 sequence_number));

BlastSearchBlkPtr LIBCALL BLASTPerformSearch PROTO((BlastSearchBlkPtr search, Int4 subject_length, Uint1Ptr subject_seq));

BLASTResultsStructPtr LIBCALL BLASTResultsStructDelete PROTO((BLASTResultsStructPtr result_struct));

SeqAlignPtr LIBCALL GetSeqAlignForResultHitList PROTO((BlastSearchBlkPtr search, Boolean getdensediag, Boolean ordinal_number, Boolean discontinuous, Boolean reverse, Boolean get_redundant_seq));

Int2 LIBCALL BlastTimeFillStructure PROTO((BlastTimeKeeperPtr btkp));

BlastSearchBlkPtr LIBCALL BlastSearchBlkDuplicate PROTO((BlastSearchBlkPtr search));

BlastSearchBlkPtr LIBCALL BlastSearchBlkNew PROTO((Int2 wordsize, Int4 qlen, CharPtr dbname, Boolean multiple_hits, BLAST_Score threshold_first, BLAST_Score threshold_second, Int4 result_size, CharPtr prog_name, BlastAllWordPtr all_words, Int2 first_context, Int2 last_context));

/* Allocates a search Block, except it only attaches to the rdfp, does not allocate it. */
BlastSearchBlkPtr LIBCALL BlastSearchBlkNewExtra PROTO((Int2 wordsize, Int4 qlen, CharPtr dbname, Boolean multiple_hits, BLAST_Score threshold_first, BLAST_Score threshold_second, Int4 result_size, CharPtr prog_name, BlastAllWordPtr all_words, Int2 first_context, Int2 last_context, ReadDBFILEPtr rdfp));

BlastSearchBlkPtr LIBCALL BlastSearchBlkDestruct PROTO((BlastSearchBlkPtr search));

BlastSearchBlkPtr BLASTSetUpSearchWithReadDbInternal PROTO((SeqLocPtr query_slp, BioseqPtr query_bsp, CharPtr prog_name, Int4 qlen, CharPtr dbname, BLAST_OptionsBlkPtr options, int (LIBCALLBACK *callback)PROTO((Int4 done, Int4 positives)), SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, Int4 gi_list_total, ReadDBFILEPtr rdfp));



Int4 LIBCALL BlastTranslateUnambiguousSequence PROTO((BlastSearchBlkPtr search, Int4 length, Uint1Ptr prot_seq, Uint1Ptr nt_seq, Int2 frame));

Uint1Ptr LIBCALL GetTranslation PROTO((Uint1Ptr query_seq, Int4 nt_length, Int2 frame, Int4Ptr length, CharPtr genetic_code));

SeqAlignPtr LIBCALL BlastGetGapAlgnTbck PROTO((BlastSearchBlkPtr search, Int4 hit_number, Boolean reverse, Boolean ordinal_number, Uint1Ptr subject, Int4 subject_length, Uint1Ptr rev_subject, Int4 rev_subject_length));

SeqAlignPtr LIBCALL BlastGetGapAlgnTbckWithReaddb PROTO((BlastSearchBlkPtr search, Int4 hit_number, Boolean ordinal_number));

Int2 LIBCALL BlastGetGappedScoreWithReaddb PROTO((BlastSearchBlkPtr search, Int4 sequence_number));

Int2 LIBCALL BlastGetGappedScore PROTO((BlastSearchBlkPtr search, Int4 subject_length, Uint1Ptr subject, Int2 frame));


SeqIdPtr LIBCALL BlastGetSubjectId PROTO((BlastSearchBlkPtr search, Int4 hit_number, Boolean ordinal_number, ValNodePtr *vnpp));



BlastAllWordPtr LIBCALL BlastAllWordNew PROTO((Int4 num_of_cols, Int4 wordsize, Boolean rows_allocated, Boolean specific));

BlastAllWordPtr LIBCALL BlastAllWordDestruct PROTO((BlastAllWordPtr all_words));


BLAST_HitListPtr LIBCALL BlastHitListDestruct PROTO((BLAST_HitListPtr hitlist));

BLAST_HitListPtr LIBCALL BlastHitListNew PROTO((BlastSearchBlkPtr search));



void LIBCALL BlastExtendWordExit PROTO((BlastSearchBlkPtr search));



Boolean LIBCALL FilterDNA PROTO((BioseqPtr bsp, Int4 filter));

Boolean LIBCALL FilterWithSeg PROTO((Uint1Ptr sequence, Int4 length, Uint1 alphabet));

BLASTResultHitlistPtr LIBCALL BLASTResultHitlistFree PROTO((BLASTResultHitlistPtr result));

BLASTResultHitlistPtr LIBCALL BLASTResultHitlistNew PROTO((Int4 hspcnt));

Nlm_FloatHi LIBCALL GetDbSubjRatio PROTO((BlastSearchBlkPtr search, Int4 subject_length));

Int2 LIBCALL BlastPreliminaryGappedScore PROTO((BlastSearchBlkPtr search, Uint1Ptr subject, Int4 subject_length, Int2 frame));

Int2 LIBCALL BlastHitListPurge PROTO((BLAST_HitListPtr hitlist));

Int4 LIBCALL HspArrayPurge PROTO((BLAST_HSPPtr PNTR hsp_array, Int4 hspcnt));


void BlastSaveCurrentHsp PROTO((BlastSearchBlkPtr search, BLAST_Score score, Int4 q_offset, Int4 s_offset, Int4 length, Int2 context));

CharPtr FormatBlastParameters PROTO((BlastSearchBlkPtr search));


SeqAlignPtr LIBCALL BioseqBlastEngineCore PROTO((BlastSearchBlkPtr search, BLAST_OptionsBlkPtr options, Int4Ptr *pos_matrix));

Uint1Ptr GetSequenceWithDenseSeg PROTO((DenseSegPtr dsp, Boolean query, Int4Ptr start, Int4Ptr length));

void BlastMaskTheResidues PROTO((Uint1Ptr buffer, Int4 max_length, Uint1 mask_residue, SeqLocPtr mask_slp, Boolean reverse));

BioseqPtr LIBCALL BlastMakeFakeBioseq PROTO((BioseqPtr bsp, CharPtr name));

BioseqPtr LIBCALL BlastDeleteFakeBioseq PROTO((BioseqPtr fake_bsp));

Boolean BlastAddSeqIdToList PROTO((BlastSearchBlkPtr search, Int4 ordinal_id, SeqIdPtr sip));

ValNodePtr BlastConstructErrorMessage PROTO((CharPtr function, CharPtr message, Uint1 level, ValNodePtr PNTR vnpp));

BlastErrorMsgPtr BlastDestroyErrorMessage PROTO((BlastErrorMsgPtr error_msg));

ValNodePtr BlastErrorChainDestroy PROTO((ValNodePtr vnp));

ValNodePtr LIBCALL BlastOtherReturnsPrepare PROTO((BlastSearchBlkPtr search));

SeqLocPtr BlastBioseqFilter PROTO((BioseqPtr bsp, CharPtr instructions));

SeqLocPtr BlastSeqLocFilter PROTO((SeqLocPtr slp, CharPtr instructions));

SeqLocPtr BlastBioseqFilterEx PROTO((BioseqPtr bsp, CharPtr instructions, BoolPtr mask_at_hash));

SeqLocPtr BlastSeqLocFilterEx PROTO((SeqLocPtr slp, CharPtr instructions, BoolPtr mask_at_hash));

SeqLocPtr MyBioseqSeg PROTO((BioseqPtr bsp_unfilter));

SeqLocPtr SeqLocSeg PROTO((SeqLocPtr slp));

Boolean BlastConvertProteinSeqLoc PROTO((SeqLocPtr slp, Int2 frame, Int4 full_length));

BlastPruneSapStructPtr BlastPruneSapStructDestruct PROTO((BlastPruneSapStructPtr prune));

BlastPruneSapStructPtr BlastPruneHitsFromSeqAlign PROTO((SeqAlignPtr sap, Int4 number, BlastPruneSapStructPtr prune));

Uint1 LIBCALL BlastGetTypes PROTO((CharPtr blast_program, Boolean PNTR query_is_na, Boolean PNTR db_is_na));

BLASTResultsStructPtr BLASTResultsStructNew PROTO((Int4 results_size, Int4 max_pieces, Int4 range_max));

Int2 BlastNtFindWords PROTO((BlastSearchBlkPtr search, Int4 start, Int4 len, Int1 context_index));

BlastAllWordPtr BlastPopulateAllWordArrays PROTO((Int4 wordsize, Int4 alphabet_size));

Uint1Ptr BlastGetSequenceFromBioseq PROTO((BioseqPtr bsp, Int4Ptr length));

BlastSeqIdListPtr BlastSeqIdListNew PROTO((void));
BlastSeqIdListPtr BlastSeqIdListDestruct PROTO((BlastSeqIdListPtr seqid_list));

Boolean BlastAdjustDbNumbers PROTO((ReadDBFILEPtr rdfp_list, Int8Ptr db_length, Int4Ptr db_number, SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, BlastDoubleInt4Ptr PNTR gi_list_pointers, Int4 gi_list_total));

BlastGiListPtr BlastGiListDestruct PROTO((BlastGiListPtr blast_gi_list));
BlastGiListPtr BlastGiListNew PROTO((BlastDoubleInt4Ptr gi_list, BlastDoubleInt4Ptr PNTR gi_list_pointers, Int4 total));


ScorePtr LIBCALL GetScoreSetFromBlastResultHsp PROTO((BLASTResultHspPtr hsp, SeqIdPtr gi_list));

Nlm_FloatHi BlastScaleMatrix PROTO((BlastMatrixRescalePtr matrix_rescale, Boolean position_dependent));

BlastMatrixRescalePtr BlastMatrixRescaleNew PROTO((Int4 alphabet_size, Int4 query_length, Uint1Ptr query,  Nlm_FloatHiPtr standardProb, Int4Ptr *matrix, Int4Ptr *private_matrix, BLAST_KarlinBlkPtr *kbp_std, BLAST_KarlinBlkPtr *kbp_psi, BLAST_KarlinBlkPtr *kbp_gap_std, BLAST_KarlinBlkPtr *kbp_gap_psi, Nlm_FloatHi lambda_ideal,  Nlm_FloatHi K_ideal));

BlastMatrixRescalePtr BlastMatrixRescaleDestruct PROTO((BlastMatrixRescalePtr matrix_rescale));

/*
        starts the awake thread using static variables in this file.
*/

void BlastStartAwakeThread PROTO((BlastSearchBlkPtr search));

/* Change the awake flag.  This thread will die in one second. */
void BlastStopAwakeThread PROTO((void));

SeqLocPtr HitRangeToSeqLoc PROTO((BlastHitRangePtr bhrp, Int4 link_value));


ValNodePtr BlastSeqLocFillDoubleInt PROTO((SeqLocPtr mask_slp, Int4 max_length, Boolean reverse));


Int2 BlastInsertList2Heap PROTO((BlastSearchBlkPtr search, BLASTResultHitlistPtr result_hitlist));

void BlastFreeHeap PROTO((BlastSearchBlkPtr search, BLASTResultHitlistPtr result_hitlist));

ValNodePtr convertSeqAlignListToValNodeList(SeqAlignPtr seqAlignList, SeqAlignPtr * lastSeqAligns, Int4 numLastSeqAligns);

SeqAlignPtr convertValNodeListToSeqAlignList(ValNodePtr seqAlignDoubleList, SeqAlignPtr ** lastSeqAligns, Int4 * numLastSeqAligns);

void LIBCALL fillCandLambda(seedSearchItems * seedSearch, Char *matrixName, BLAST_OptionsBlkPtr options);

Int2 RealBlastGetGappedAlignmentTraceback(BlastSearchBlkPtr search, Uint1Ptr subject, Int4 subject_length, Uint1Ptr rev_subject, Int4 rev_subject_length, SeqIdPtr subject_id, BLAST_HSPPtr *hsp_array, Int4 hspcnt, SeqAlignPtr *head, BlastHitRangePtr bhrp, Int4 min_score_to_keep, Boolean reverse, Int4 ordinal_id, Boolean do_traceback);

#ifdef __cplusplus
}
#endif
#endif /* !__BLASTPRI__ */