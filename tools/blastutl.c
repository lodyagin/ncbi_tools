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

File name: blastutl.c

Author: Tom Madden

Contents: Utilities for BLAST

******************************************************************************/
/* $Revision: 6.89 $ */
/* $Log: blastutl.c,v $
 * Revision 6.89  1998/09/24 15:26:38  egorov
 * Fix lint complaints
 *
 * Revision 6.88  1998/09/16 19:00:16  madden
 * Added subset Boolean
 *
 * Revision 6.87  1998/09/15 13:12:29  madden
 * Fixed memory leak
 *
 * Revision 6.86  1998/09/14 15:11:18  egorov
 * Add support for Int8 length databases; remove unused variables
 *
 * Revision 6.85  1998/09/04 20:48:48  madden
 * typo fix (= instead of ==)
 *
 * Revision 6.84  1998/09/03 20:23:42  madden
 * Copied seq_ext and seq_ext_type in MakeFakeBioseq
 *
 * Revision 6.83  1998/09/03 19:41:09  madden
 * do not switch sequences for Blast2Sequences if filtering is performed
 *
 * Revision 6.82  1998/08/24 14:59:59  madden
 * readdb_get_sequence_ex function
 *
 * Revision 6.81  1998/07/30 19:00:56  madden
 * Fix memory leak
 *
 * Revision 6.80  1998/07/29 21:29:45  madden
 * Fixed UMR with longest_db_seq that showed up in Blast 2 sequences
 *
 * Revision 6.79  1998/07/28 21:18:35  madden
 * Change to BLAST_ExtendWordParamsNew saves memory
 *
 * Revision 6.78  1998/07/24 14:58:53  madden
 * Jinqhuis call to SeqLocRevCmp put back
 *
 * Revision 6.77  1998/07/22 20:31:51  madden
 * Replaced cutvalue of 1000000 with INT4_MAX
 *
 * Revision 6.76  1998/07/22 12:17:03  madden
 * Added BioseqHitRange call for repeat filtering
 *
 * Revision 6.75  1998/07/21 20:58:10  madden
 * Changes to allow masking at hash only
 *
 * Revision 6.74  1998/07/20 15:51:28  zjing
 * add a check for plus-minus before SeqLocRevCmp
 *
 * Revision 6.73  1998/07/17 15:39:59  madden
 * Changes for Effective search space.
 *
 * Revision 6.72  1998/07/14 21:31:43  madden
 * Fix for incorrectly sorted HSP bug and speed-up of CheckHspOverlap
 *
 * Revision 6.71  1998/07/06 13:39:04  madden
 * Fixed improper use of Int4 in parse_seg_options
 *
 * Revision 6.70  1998/07/02 21:00:39  egorov
 * Remove memory leak in threaded version
 *
 * Revision 6.69  1998/06/12 22:09:14  madden
 * Added call to SegParamsFree
 *
 * Revision 6.68  1998/06/12 16:08:51  madden
 * BlastHitRange stuff
 *
 * Revision 6.67  1998/06/08 15:07:32  madden
 * Fixed bug in BlastConvertProteinSeqLoc
 *
 * Revision 6.66  1998/06/04 16:23:17  madden
 * Use new seg
 *
 * Revision 6.65  1998/05/28 19:59:58  madden
 * Zhengs new culling code
 *
 * Revision 6.64  1998/05/22 20:20:38  madden
 * Added BlastTwoSequencesByLocEx and BlastTwoSequencesEx
 *
 * Revision 6.63  1998/05/18 17:58:31  madden
 * fixed parsing of coil-coil options, added parsing of dust options
 *
 * Revision 6.62  1998/05/17 16:28:41  madden
 * Allow changes to filter options and cc filtering.
 *
 * Revision 6.61  1998/05/05 14:05:35  madden
 * Added functions BlastStartAwakeThread and BlastStopAwakeThread
 *
 * Revision 6.60  1998/04/28 21:04:19  madden
 * Reset number of HSPs to zero if relinking
 *
 * Revision 6.59  1998/04/24 21:52:09  madden
 * Protection against NULL pointers
 *
 * Revision 6.58  1998/04/24 19:10:59  egorov
 * Fix bug when if wordsize == 2 blastall produces extra alignments
 *
 * Revision 6.57  1998/04/23 21:15:09  egorov
 * Show exact matching even if score is below threshold (case of two sequences)
 *
 * Revision 6.56  1998/04/15 20:24:54  madden
 * BlastMaskTheResidues optimized
 *
 * Revision 6.55  1998/04/10 17:46:58  madden
 * Changed FALSE to NULL in BioseqSeg
 *
 * Revision 6.54  1998/04/02 21:12:55  madden
 * Properly set value for linking HSPs in blastx and tblastn
 *
 * Revision 6.53  1998/04/01 22:47:35  madden
 * Check for query_invalid flag
 *
 * Revision 6.52  1998/03/26 14:20:20  madden
 * Changed GetScoreSetFromBlastResultHsp1 from static to LIBCALL
 *
 * Revision 6.51  1998/03/25 22:28:16  madden
 * Changes to allow random access BLAST by gi
 *
 * Revision 6.50  1998/03/24 15:38:25  madden
 * Use BlastDoubleInt4Ptr to keep track of gis and ordinal_ids
 *
 * Revision 6.49  1998/03/19 22:16:24  madden
 * Changes to allow blasting by gi list
 *
 * Revision 6.48  1998/03/18 14:14:11  madden
 * Support random access by gi list
 *
 * Revision 6.47  1998/03/16 17:41:59  madden
 * Fixed leaks
 *
 * Revision 6.46  1998/03/14 18:28:10  madden
 * Added BioseqBlastEngineEx
 *
 * Revision 6.45  1998/03/09 16:35:10  madden
 * Fixed bug with tblastn and blastx gapped searches
 *
 * Revision 6.44  1998/02/27 14:32:33  madden
 * Functions moved to blastool.c
 *
 * Revision 6.43  1998/02/26 22:34:27  madden
 * Changes for 16 bit windows
 *
 * Revision 6.42  1998/02/26 19:12:39  madden
 *  Removed AdjustOffSetsInSeqAlign, added BlastNtFindWords BlastPopulateAllWordArrays BlastFindWords and BlastNewFindWords
 *
 * Revision 6.41  1998/02/24 22:47:06  madden
 * Fixed problem with Option validation
 *
 * Revision 6.40  1998/02/23 16:09:57  madden
 * Corrected from offset for subject in tblastx search
 *
 * Revision 6.39  1998/02/19 17:17:05  madden
 * Use of Int4 rather than Int2 when pruning SeqAlign
 *
 * Revision 6.38  1998/02/12 21:50:39  madden
 * protection against NULL hitlist in blastx and tblastn
 *
 * Revision 6.37  1998/02/11 17:18:19  madden
 * Made BlastGetGappedAlignmentTraceback functions to BlastGetGapAlgnTbck (shorter than 32 chars)
 *
 * Revision 6.36  1998/01/31 21:34:09  madden
 * Fix to SeqAlign pruning
 *
 * Revision 6.35  1998/01/06 18:26:22  madden
 * Use SeqLocLen rather than bsp->length, wordsize done properly for nucl
 *
 * Revision 6.34  1998/01/05 22:41:40  madden
 * Added seqalign_reverse_strand
 *
 * Revision 6.33  1998/01/05 20:53:16  madden
 * Added ability to align minus-minus or plus-minus in BlastTwoSeqsByLoc
 *
 * Revision 6.32  1998/01/05 16:46:55  madden
 * One or both strands can be searched, as opposed to only both, changes to number of contexts
 *
 * Revision 6.31  1997/12/31 17:52:09  madden
 * Change to BLAST_WordFinderNew
 *
 * Revision 6.30  1997/12/23 19:16:52  madden
 * Minor efficiency in ExtendWordExit
 *
 * Revision 6.29  1997/12/23 18:12:34  madden
 * Changes for range-dependent blast
 *
 * Revision 6.28  1997/12/12 20:38:55  madden
 * ContextToFrame lost last parameter, fix to sprintf
 *
 * Revision 6.27  1997/12/11 22:22:24  madden
 * Proper casting of variables
 *
 * Revision 6.26  1997/12/10 22:43:09  madden
 * proper casting
 *
 * Revision 6.25  1997/12/01 22:07:10  madden
 * Changed call to BLASTOptionValidateEx
 *
 * Revision 6.24  1997/11/28 18:19:33  madden
 * Changes to TxDfDbInfoNew
 *
 * Revision 6.23  1997/11/18 22:23:20  madden
 * Added BLASTOptionSetGapParams
 *
 * Revision 6.22  1997/11/14 17:15:29  madden
 * Realign matches when they contain ambiguities in blastx/tblastn
 *
 * Revision 6.21  1997/11/07 00:49:02  madden
 * Added call to BLAST_MatrixFill
 *
 * Revision 6.20  1997/10/29 22:11:13  madden
 * ABS value of frames
 *
 * Revision 6.19  1997/10/24 20:44:52  madden
 * Removed BlastSetReadDB and BlastGetReadDB_ID
 *
 * Revision 6.18  1997/10/22 21:46:34  madden
 * Changed default values
 *
 * Revision 6.17  1997/10/21 20:39:18  madden
 * Fix for more alignments than descriptions.
 *
 * Revision 6.16  1997/10/21 19:50:00  madden
 * Fix for no valid query sequence and hitlist_max of 1
 *
 * Revision 6.15  1997/10/03 21:27:28  madden
 * Added BlastGetTypes
 *
 * Revision 6.14  1997/10/02 17:29:29  madden
 * Added PrintDbInformationBasic
 *
 * Revision 6.13  1997/10/01 13:35:31  madden
 * Changed BLAST_VERSION to BLAST_ENGINE_VERSION
 *
 * Revision 6.12  1997/09/30 20:03:07  madden
 * Saved db filename in dbinfo
 *
 * Revision 6.11  1997/09/24 22:36:35  madden
 * Fixes for MT multidb searches
 *
 * Revision 6.10  1997/09/23 16:43:41  madden
 * removed unneeded DenseSegPtr
 *
 * Revision 6.9  1997/09/22 18:18:35  madden
 * Added umlaut to Schaffer in reference
 *
 * Revision 6.8  1997/09/18 22:22:03  madden
 * Added prune functions
 *
 * Revision 6.7  1997/09/16 16:54:09  kans
 * return FASLE instead of NULL for Boolean value
 *
 * Revision 6.6  1997/09/16 16:31:28  madden
 * More changes for multiple db runs
 *
 * Revision 6.5  1997/09/11 18:49:31  madden
 * Changes to enable searches against multiple databases.
 *
 * Revision 6.4  1997/09/10 21:28:00  madden
 * Changes to set CPU limits
 *
 * Revision 6.3  1997/09/08 16:25:32  madden
 * Fixed bug that did not mask low-complexity regions at the end of a query
 *
 * Revision 6.2  1997/08/27 14:46:51  madden
 * Changes to enable multiple DB searches
 *
 * Revision 6.1  1997/08/26 15:05:26  madden
 * Fix for negative effective search space
 *
 * Revision 6.0  1997/08/25 18:52:49  madden
 * Revision changed to 6.0
 *
 * Revision 1.105  1997/08/22 18:37:43  madden
 * Added function BlastOtherReturnsPrepare
 *
 * Revision 1.104  1997/08/20 21:43:34  madden
 * Added page numbers
 *
 * Revision 1.103  1997/08/14 21:07:08  madden
 * ignored gapped for tblastx
 *
 * Revision 1.102  1997/08/14 14:30:35  madden
 * BlastNewFindWords called with range set for ranged blast
 *
 * Revision 1.101  1997/07/31 21:18:11  madden
 * Removed left-over file from seg
 *
 * Revision 1.100  1997/07/30 16:39:30  madden
 * Print gap existence and extension parameters for blastn
 *
 * Revision 1.99  1997/07/30 16:31:37  madden
 * tblastx prepares StdSeg
 *
 * Revision 1.98  1997/07/29 17:07:27  madden
 * better tblastx error messages.
 *
 * Revision 1.97  1997/07/25 15:39:49  madden
 * Corrected citation
 *
 * Revision 1.96  1997/07/25 13:47:46  madden
 * Made buffer longer to avoid ABR
 *
 * Revision 1.95  1997/07/23 20:59:02  madden
 * Changed blastn defaults for gap opening and extension
 *
 * Revision 1.94  1997/07/22 17:22:41  madden
 * Added NULL arg (for index callback) to BLASTSetUpSearch funcs
 *
 * Revision 1.93  1997/07/21 17:36:42  madden
 * Added BlastGetReleaseDate
 *
 * Revision 1.92  1997/07/18 20:57:02  madden
 * Added functions BlastGetVersionNumber and BlastGetReference
 *
 * Revision 1.91  1997/07/18 14:26:20  madden
 * call to AcknowledgeBlastQuery changed, SeqId no longer deleted there.
 *
 * Revision 1.90  1997/07/16 20:34:35  madden
 * Added function BlastConvertProteinSeqLoc
 *
 * Revision 1.89  1997/07/15 20:36:14  madden
 * Added BioseqSeg and SeqLocSeg
 *
 * Revision 1.88  1997/07/14 20:11:10  madden
 * Removed unused variables
 *
 * Revision 1.87  1997/07/14 16:15:41  madden
 * call to BLASTOptionValidateEx in BlastBioseqEngine
 *
 * Revision 1.86  1997/07/14 15:31:49  madden
 * Added BlastErrorMessage functions
 *
 * Revision 1.85  1997/07/11 19:29:37  madden
 * Added function BioseqBlastEngineByLoc
 *
 * Revision 1.84  1997/07/10 20:35:43  madden
 * Changed parameter output
 *
 * Revision 1.83  1997/07/02 20:18:39  madden
 * Made continuous SeqAlign the default
 *
 * Revision 1.82  1997/07/02 18:31:39  madden
 * changed defaults
 *
 * Revision 1.81  1997/07/01 19:15:44  madden
 * More changes to FormatBlastParameters
 *
 * Revision 1.80  1997/07/01 17:51:36  madden
 * changed gap_decay rate, gap_prob
 *
 * Revision 1.79  1997/07/01 15:44:44  madden
 * Changes to FormatBlastParameters per S. Altschul
 *
 * Revision 1.78  1997/06/30 15:50:06  madden
 * Changes to FormatBlastParameters
 *
 * Revision 1.77  1997/06/27 22:18:51  madden
 * Updated default parameters
 *
 * Revision 1.76  1997/06/27 14:31:08  madden
 * Added functions BlastAddSeqIdToList and BlastSeqIdListDestruct
 *
 * Revision 1.75  1997/06/24 13:51:27  madden
 * Fixed SeqLoc leak
 *
 * Revision 1.74  1997/06/23 20:49:31  madden
 * BLASTOptionValidate checks for proper gapping parameters
 *
 * Revision 1.73  1997/06/20 13:11:33  madden
 * Made AdjustOffSetsInSeqAlign non-static, Fixed purify error
 *
 * Revision 1.72  1997/06/06 21:29:48  madden
 * Added Boolean html to AcknowledgeBlastQuery and PrintDbInformation
 *
 * Revision 1.71  1997/06/06 19:49:46  madden
 * Added BlastMakeFakeBioseq and BlastDeleteFakeBioseq
 *
 * Revision 1.70  1997/05/30 21:05:59  madden
 * corrected call to readdb_new
 *
 * Revision 1.69  1997/05/27 20:20:02  madden
 * Added function BlastMaskTheResidues
 *
 * Revision 1.68  1997/05/22 21:24:55  madden
 * Added support for final gapX dropoff value
 *
 * Revision 1.67  1997/05/20 17:52:58  madden
 * Added functions BlastTwoSequencesByLoc and BlastSequencesOnTheFlyByLoc
 *
 * Revision 1.66  1997/05/12 21:34:16  madden
 * readdb_new allows indeterminate database type
 *
 * Revision 1.65  1997/05/06 22:17:59  madden
 * Duplicate dblen_eff, dbseq_num, and length_adjustment
 *
 * Revision 1.64  1997/05/01  15:53:19  madden
 * Addition of extra KarlinBlk's for psi-blast
 *
 * Revision 1.63  1997/04/29  14:07:45  madden
 * Fixed problem with hits failing PreliminaryGapping; fixed UMR.
 *
 * Revision 1.62  1997/04/25  20:23:06  madden
 * Freed SeqPort to clear mem leak.
 *
 * Revision 1.61  1997/04/24  14:43:07  madden
 * Fix for minus strand (ungapped) tblastn runs.
 *
 * Revision 1.60  1997/04/23  21:56:07  madden
 * Changes in BlastGetGappedAlignmentTraceback for in-frame gapping tblastn.
 *
 * Revision 1.59  1997/04/22  14:00:14  madden
 * Removed unused variables.
 *
 * Revision 1.58  1997/04/22  13:04:19  madden
 * Changes for in-frame blastx gapping.
 *
 * Revision 1.57  1997/04/21  15:35:26  madden
 * Fixes for 'gapped' StdSegs.
 *
 * Revision 1.56  1997/04/18  17:08:35  madden
 * Corrected printing of threshold values.
 *
 * Revision 1.55  1997/04/17  22:12:43  madden
 * Fix for offset in GetStartForGappedAlignment.
 *
 * Revision 1.54  1997/04/17  22:07:48  madden
 * Changes to allow in-frame gapped tblastn.
 *
 * Revision 1.53  1997/04/15  22:02:59  madden
 * Set original_length1 for translating searches.
 *
 * Revision 1.52  1997/04/14  21:31:58  madden
 * Checking for NULL pointer.
 *
 * Revision 1.51  1997/04/14  15:59:47  madden
 * Changes for ungapped psi-blast.
 *
 * Revision 1.50  1997/04/11  21:18:45  madden
 * Added GetSequenceWithDenseSeg.
 *
 * Revision 1.49  1997/04/11  19:02:49  madden
 * Changes for in-frame blastx, tblastn gapping.
 *
 * Revision 1.48  1997/04/09  20:01:53  madden
 * Copied seqid_list from search structure to duplicate, for use on threads.
 *
 * Revision 1.47  1997/04/08  16:27:28  madden
 * Fixed leaks; fix for blastn formatting of parameters.
 *
 * Revision 1.46  1997/04/07  21:42:56  madden
 * Freed SeqLocPtr used for dust.
 *
 * Revision 1.45  1997/04/07  18:17:09  madden
 * Formatted parameters for Stephen.
 *
 * Revision 1.44  1997/04/04  20:44:09  madden
 * Added check for NULL return.
 *
 * Revision 1.43  1997/04/04  20:42:35  madden
 * Added function BioseqBlastEngineCore.
 *
 * Revision 1.42  1997/04/03  19:50:56  madden
 * Changes to use effective database length instead of the length of each
 * sequence in statistical calculations.
 *
 * Revision 1.41  1997/03/27  22:30:51  madden
 * Correctly checked for overlapping HSP's.
 *
 * Revision 1.40  1997/03/20  22:56:24  madden
 * Added gap_info to hsp.
 *
 * Revision 1.39  1997/03/20  21:52:10  madden
 * Fix for segmented query BioseqPtr when gapped alignment is performed.
 *
 * Revision 1.39  1997/03/20  21:52:10  madden
 * Fix for segmented query BioseqPtr when gapped alignment is performed.
 *
 * Revision 1.38  1997/03/14  22:06:11  madden
 * fixed MT bug in BlastReevaluateWithAmbiguities.
 *
 * Revision 1.37  1997/03/14  15:57:23  madden
 * Removed superfluous call to SeqAlignNew
 *
 * Revision 1.36  1997/03/14  15:22:11  madden
 * Fixed UMR of seqalign in BlastTwoSequencesCore.
 *
 * Revision 1.35  1997/03/11  14:38:40  madden
 * Added BlastSequencesOnTheFly and BlastTwoSequencesCore.
 *
 * Revision 1.34  1997/03/07  22:35:54  madden
 * Fix for BLASTOptionNew.
 *
 * Revision 1.33  1997/03/07  21:58:36  madden
 * Added Boolean gapped argument to BLASTOptionNew.
 *
 * Revision 1.32  1997/03/07  21:11:22  madden
 * Added in check for blastn on gapped calculations.
 *
 * Revision 1.31  1997/03/06  21:47:27  madden
 * Made FormatBlastParameters non-static.
 *
 * Revision 1.30  1997/03/05  18:16:16  madden
 * SeqIdFree replaced by SeqIdSetFree, fixed memory leak.
 *
 * Revision 1.29  1997/03/05  14:29:46  madden
 * Moved BlastSaveCurrentHsp from blast.c; Added function CheckHspOverlap.
 *
 * Revision 1.28  1997/03/04  21:34:59  madden
 * Added in HspArrayPurge.
 *
 * Revision 1.27  1997/03/04  20:08:19  madden
 * Moved gapped alignment code from blast.c to blastutl.c
 *
 * Revision 1.26  1997/03/03  22:39:45  madden
 * Moved code from blast.c to blastutl.c.
 *
 * Revision 1.25  1997/03/03  21:47:22  madden
 * Moved functions from blast.c to blastutl.c for 16-bit windows.
 *
 * Revision 1.24  1997/03/03  20:58:09  madden
 * Fixed offsets for minus strands.
 *
 * Revision 1.23  1997/03/03  17:30:21  madden
 * Set SeqAlignPtr to NULL in BlastTwoSequences and BlastBioseqEngine, possible UMR.
 *
 * Revision 1.22  1997/03/01  18:25:33  madden
 * reverse flag added to BlastGetGappedAlignmentTraceback functions.
 *
 * Revision 1.21  1997/02/27  22:47:07  madden
 * Replaced tblastx with tblastn in BioseqBlastEngine.
 *
 * Revision 1.20  1997/02/26  23:39:54  madden
 * Added Txdfline stuff.
 *
 * Revision 1.19  1997/02/26  20:37:31  madden
 * Added *error_returns to BioseqBlastEngine.
 *
 * Revision 1.18  1997/02/25  19:17:05  madden
 * Changes to BioseqBlastEngine.
 *
 * Revision 1.17  1997/02/20  23:00:34  madden
 * Checked for NULL return in BlastTwoSequences.
 *
 * Revision 1.16  1997/02/20  18:38:34  madden
 * Set Default db_length to zero in Options.
 *
 * Revision 1.15  1997/02/19  16:25:22  madden
 * Reset gapped_calculation for blastn; returned proper SeqAlign for blastx, tblastn
 * in BioseqBlastEngine.
 *
 * Revision 1.14  1997/02/19  13:45:13  madden
 * replaced zero in call to BlastGetGappedAlignmentTraceback with FALSE.
 *
 * Revision 1.13  1997/02/18  22:09:02  madden
 * Removed unused variable.
 *
 * Revision 1.12  1997/02/18  21:03:00  madden
 * Changes to BioseqBlastEngine for gapped calculations.
 *
 * Revision 1.11  1997/02/18  18:31:34  madden
 * Used SeqIdFindBest in BlastTwoSequences.
 *
 * Revision 1.10  1997/02/18  17:58:52  madden
 * Added BioseqBlastEngine.
 *
 * Revision 1.9  1997/02/14  17:17:59  madden
 * Changes to default options and BlastTwoSequences for nucl.
 * sequences with ambiguites.
 *
 * Revision 1.8  1997/02/13  18:23:56  madden
 * Fixed ID type from BlastTwoSequences.
 *
 * Revision 1.7  1997/02/11  19:30:54  madden
 * Changes to BlastTwoSequences for gapped alignments.
 *
 * Revision 1.6  1997/02/10  20:03:58  madden
 * BlastTwoSequences indexes only the subject.
 *
 * Revision 1.5  1997/02/10  15:24:26  madden
 * Removed unused variable.
 *
 * Revision 1.4  1997/02/07  22:43:03  madden
 * Moved BLAST_WordFinderNew and Destruct from blast.c to blastutl.c, made
 * non-static.
 *
 * Revision 1.3  1997/02/07  22:32:40  madden
 * Changed prototypes for BlastGetSubjectId and GetSeqAlignForResultHitList.
 *
 * Revision 1.2  1997/02/05  13:36:48  madden
 * Removed Unused variable.
 *
 * Revision 1.1  1997/02/04  18:23:58  madden
 * Initial revision
 *
*/

#include <ncbi.h>
#include <blastpri.h>
#include <objcode.h>
#include <objseq.h>
#include <sequtil.h>
#include <tofasta.h>
#include <seqport.h>
#include <readdb.h>
#include <ncbithr.h>
#include <dust.h>
#include <urkpcc.h>
#include <txalign.h>
#include <seg.h>

/* Window size used to scan HSP for highest score region, where gapped
extension starts. */
#define HSP_MAX_WINDOW 11

#define BLASTFILTER_DIR "/usr/ncbi/blast/filter"

static Int2 ContextToFrame PROTO((BlastSearchBlkPtr search, Int2 context_number));
static Boolean SumBlastGetGappedAlignmentEx PROTO((BlastSearchBlkPtr search, Int4 hit_number, Boolean reverse, Boolean ordinal_number, Uint1Ptr subject, Int4 subject_length, Boolean do_traceback, SeqAlignPtr PNTR seqalignP, BlastHitRangePtr bhrp));

/*
	Goes through the list of gi's/ordinal id's looking for matches
	to the ordinal ID.  Returns those acceptable gi's as SeqIdPtr's.
*/
static SeqIdPtr
BlastGetAllowedGis (BlastSearchBlkPtr search, Int4 ordinal_id)

{
	BlastGiListPtr blast_gi_list;
	Boolean found=FALSE;
	Int4 index, total;
	ValNodePtr gi_list=NULL;

	gi_list = NULL;
	if (search->blast_gi_list)
	{
		blast_gi_list = search->blast_gi_list;
		total = blast_gi_list->total;
		found = FALSE;
		for (index=0; index<total; index++)
		{
			if (ordinal_id == blast_gi_list->gi_list_pointer[index]->ordinal_id)
			{
				found = TRUE;
				break;
			}
		}

		if (found)
		{
			while (index < total)
			{
				if (ordinal_id == blast_gi_list->gi_list_pointer[index]->ordinal_id)
				{
					ValNodeAddInt(&gi_list, SEQID_GI, blast_gi_list->gi_list_pointer[index]->gi);
				}
				else
				{
					break;
				}
				index++;
			}
		}
	}
	return (SeqIdPtr) gi_list;
}

/* 
	SOME FUNCTIONS TO PRODUCE A SeqAlign from the BLAST results.
*/

/*****************************************************************************

	Finds the best SeqId for the SeqAlign.  Looks for the GI, then takes
	anything if that's not found and makes up a local ID if no ID is
	found at all.
*****************************************************************************/

static SeqIdPtr
GetTheSeqAlignID(SeqIdPtr seq_id)
{
	SeqIdPtr new_id, ret_id;
	ObjectIdPtr obidp;
	
	ret_id = NULL;
	if (seq_id)
	{
		/* Get the gi from the chain, if it's there. */
		new_id = SeqIdFindBest(seq_id, SEQID_GI);
		if (new_id)
		{
			ret_id = SeqIdDup(new_id);
		}
		else
		{	/* No Gi was found, use any ID. */
			ret_id = SeqIdDup(seq_id);
		}
	}

	if (ret_id == NULL)
	{	/* make up an ID. */
		obidp = ObjectIdNew();
		obidp->str = StringSave("lcl|unknown");
		ValNodeAddPointer(&ret_id, SEQID_LOCAL, obidp);
	}

	return ret_id;
}
static SeqAlignPtr 
FillInSegsInfo(SeqAlignPtr sap_head, StdSegPtr ssp_head, DenseDiagPtr ddp_head)

{
	SeqAlignPtr sap;

	if (ddp_head || ssp_head)
	{
		if (sap_head)
		{
			sap = sap_head;
			while (sap->next)
				sap = sap->next;
			sap->next = SeqAlignNew();
			sap = sap->next;
		}
		else
		{
			sap_head = sap = SeqAlignNew();
		}

		if (ddp_head)
		{
			sap->type = 2;
			sap->segs = ddp_head;
			sap->segtype = 1;
		}
		else if (ssp_head)
		{
			sap->type = 2;
			sap->segs = ssp_head;
			sap->segtype = 3;
		}
	}
	return sap_head;
}


/*************************************************************************
*
*	This function fills in the DenseDiag Information from the variable
*	hsp.  On the first call to this function *old should be
*	NULL, after that pass in the head of the DenseDiagPtr chain.
*	The newest DenseDiagPtr is returned.
*
************************************************************************/

static DenseDiagPtr
FillInDenseDiagInfo(DenseDiagPtr PNTR old, BLASTResultHspPtr hsp, Boolean reverse, Int4 query_length, Int4 subject_length, SeqIdPtr gi_list)

{
	DenseDiagPtr		ddp, new;

	new = DenseDiagNew();
	
	new->dim = 2;	/* Only 2 is supported in spec. */
	new->len = hsp->query_length;
	new->starts = (Int4Ptr) MemNew(2 * sizeof(Int4));
	new->strands = (Uint1Ptr) MemNew(2 * sizeof(Uint1));
	if (reverse)
	{
		if (hsp->subject_frame >= 0)
		{
			new->strands[0] = Seq_strand_plus;
			new->starts[0] = hsp->subject_offset;
		}
		else
		{
			new->strands[0] = Seq_strand_minus;
			new->starts[0] = subject_length - hsp->subject_offset - hsp->subject_length;
		}
		if (hsp->query_frame >= 0)
		{
			new->strands[1] = Seq_strand_plus;
			new->starts[1] = hsp->query_offset;
		}
		else
		{
			new->strands[1] = Seq_strand_minus;
			new->starts[1] = query_length - hsp->query_offset - hsp->query_length;
		}
	}
	else
	{
		if (hsp->query_frame >= 0)
		{
			new->strands[0] = Seq_strand_plus;
			new->starts[0] = hsp->query_offset;
		}
		else
		{
			new->strands[0] = Seq_strand_minus;
			new->starts[0] = query_length - hsp->query_offset - hsp->query_length;
		}
		if (hsp->subject_frame >= 0)
		{
			new->strands[1] = Seq_strand_plus;
			new->starts[1] = hsp->subject_offset;
		}
		else
		{
			new->strands[1] = Seq_strand_minus;
			new->starts[1] = subject_length - hsp->subject_offset - hsp->subject_length;
		}
	}
	new->scores = GetScoreSetFromBlastResultHsp(hsp, gi_list);

/* Go to the end of the chain, and then attach "new" */
	if (*old)
	{
		ddp = *old;
		while (ddp->next)
			ddp = ddp->next;
		ddp->next = new;
	}
	else
	{
		*old = new;
	}

	new->next = NULL;

	return new;
}

/*************************************************************************
*
*	This function fills in the StdSeg Information from the variable
*	hsp.  On the first call to this function *old should be
*	NULL, after that pass in the head of the DenseDiagPtr chain.
*	The newest StdSegPtr is returned.
*
************************************************************************/

static StdSegPtr
FillInStdSegInfo(BlastSearchBlkPtr search, BLASTResultHitlistPtr hitlist, StdSegPtr PNTR old, BLASTResultHspPtr hsp, SeqIdPtr sip, Boolean reverse, SeqIdPtr gi_list)

{
	Int4			subject_length;
	StdSegPtr		ssp, new;
	SeqIdPtr		query_sip, subject_sip;
	SeqIntPtr		seq_int1, seq_int2;
	SeqLocPtr		slp=NULL;

	new = StdSegNew();
/* Duplicate the id and split it up into query and subject parts */
	query_sip = SeqIdDup(sip);
	subject_sip = SeqIdDup(sip->next);
	
	new->dim = 2;	/* Only 2 is supported in spec. */
	seq_int1 = SeqIntNew();
	if (hsp->query_frame == 0)
	{
		seq_int1->from = hsp->query_offset;
		seq_int1->to = hsp->query_offset + hsp->query_length - 1;
		seq_int1->strand = Seq_strand_unknown;
	}
	else if (hsp->query_frame < 0)
	{
		seq_int1->to = search->context[hsp->context].query->original_length - CODON_LENGTH*hsp->query_offset + hsp->query_frame;
		seq_int1->from = search->context[hsp->context].query->original_length - CODON_LENGTH*(hsp->query_offset+hsp->query_length) + hsp->query_frame + 1;
		seq_int1->strand = Seq_strand_minus;
	}
	else if (hsp->query_frame > 0)
	{
		seq_int1->from = CODON_LENGTH*(hsp->query_offset) + hsp->query_frame - 1;
		seq_int1->to = CODON_LENGTH*(hsp->query_offset+hsp->query_length) + hsp->query_frame - 2;
		seq_int1->strand = Seq_strand_plus;
	}
	seq_int1->id = query_sip;
	seq_int2 = SeqIntNew();
	if (hsp->subject_frame == 0)
	{
		seq_int2->from = hsp->subject_offset;
		seq_int2->to = hsp->subject_offset + hsp->subject_length - 1;
		seq_int2->strand = Seq_strand_unknown;
	} 
	else if (hsp->subject_frame < 0)
	{
	    	if (search->rdfp)
			subject_length = readdb_get_sequence_length(search->rdfp, hitlist->subject_id);
	    else if (hitlist->subject_info)
			subject_length = hitlist->subject_info->length;
	    else
			subject_length = 0;
		seq_int2->from = subject_length - CODON_LENGTH*(hsp->subject_offset + hsp->subject_length) + hsp->subject_frame + 1;
		seq_int2->to = subject_length - CODON_LENGTH*(hsp->subject_offset) + hsp->subject_frame;
		seq_int2->strand = Seq_strand_minus;
	}
	else if (hsp->subject_frame > 0)
	{
		seq_int2->from = CODON_LENGTH*(hsp->subject_offset) + hsp->subject_frame - 1;
		seq_int2->to = CODON_LENGTH*(hsp->subject_offset + hsp->subject_length) + hsp->subject_frame - 2;
		seq_int2->strand = Seq_strand_plus;
	}
	seq_int2->id = subject_sip;

	if (reverse)
	{
		ValNodeAddPointer(&slp, SEQLOC_INT, seq_int2); 
		ValNodeAddPointer(&slp, SEQLOC_INT, seq_int1); 
	}
	else
	{
		ValNodeAddPointer(&slp, SEQLOC_INT, seq_int1); 
		ValNodeAddPointer(&slp, SEQLOC_INT, seq_int2); 
	}
	new->loc = slp;
	new->scores = GetScoreSetFromBlastResultHsp(hsp, gi_list);

/* Go to the end of the chain, and then attach "new" */
	if (*old)
	{
		ssp = *old;
		while (ssp->next)
			ssp = ssp->next;
		ssp->next = new;
	}
	else
	{
		*old = new;
	}

	new->next = NULL;

	return new;
}


/************************************************************************
*
*	This function assembles all the components of the Seq-align from
*	a "sparse" BLAST HitList.  "sparse" means that the hitlist 
*	may contain no sequence and not even a descriptor.  It is only 
*	required to contain the sequence_number that readdb refers to
*	and scoring/alignment information.
*
*	If dbname is non-NULL, then only a general ("gnl") ID is 
*	issued, with the ordinal number of the subject sequence in
*	the ObjectIdPtr.
*
*	Boolean reverse: reverse the query and db order in SeqAlign.
*
************************************************************************/
SeqAlignPtr LIBCALL
GetSeqAlignForResultHitList(BlastSearchBlkPtr search, Boolean getdensediag, Boolean ordinal_number, Boolean discontinuous, Boolean reverse, Boolean get_redundant_seqs)

{
	BLASTResultHspPtr	hsp;
	BLASTResultHitlistPtr	results;
	BLASTResultsStructPtr	result_struct;
	DenseDiagPtr		ddp_head=NULL, ddp;
	SeqIdPtr		gi_list, sip, sip_subject, sip_subject_start;
	StdSegPtr		ssp_head=NULL, ssp;
	SeqAlignPtr		last, seqalign_head, seqalign, sap_head;
	Int4 			hsp_cnt, index, index2, hspset_cnt_old;
	Int4			hitlist_count;
	Int4			subject_length;
	ValNodePtr		vnp, vnp_start;

	ddp_head = NULL;
	ssp_head = NULL;
	sap_head = NULL;
	seqalign_head = NULL;

	result_struct = search->result_struct;
	hitlist_count = result_struct->hitlist_count;

	last = NULL;
	sip = NULL;
	sip_subject_start = NULL;
	for (index=0; index<hitlist_count; index++)
	{
	    results = result_struct->results[index];
	    sip_subject_start = NULL;
	    if (get_redundant_seqs)
	    {
		vnp = NULL;
	    	sip = BlastGetSubjectId(search, index, ordinal_number, &vnp);
		vnp_start = vnp;
		while (vnp)
		{
			sip = GetTheSeqAlignID(vnp->data.ptrvalue);
			SeqIdFree(vnp->data.ptrvalue);
			if (sip_subject_start == NULL)
			{
				sip_subject_start = sip;
			}
			else
			{
				sip_subject = sip_subject_start;
				while (sip_subject->next)
					sip_subject = sip_subject->next;
				sip_subject->next = sip;
			}
			vnp = vnp->next;
		}
		vnp_start = vnp = ValNodeFree(vnp_start);
	    }
	    else
	    {
	    	sip = BlastGetSubjectId(search, index, ordinal_number, NULL);
	    	sip_subject_start = sip_subject = GetTheSeqAlignID(sip);
	    	sip = SeqIdSetFree(sip);
	    }

	    results = result_struct->results[index];
	    if (search->rdfp)
		subject_length = readdb_get_sequence_length(search->rdfp, results->subject_id);
	    else if (results->subject_info)
			subject_length = results->subject_info->length;
	    else
			subject_length = 0;

	gi_list = BlastGetAllowedGis(search, results->subject_id);
	/* right now sip_subject should only contain one ID.  At some
	point it will contain multiple ID's for identical sequences. */
	    sip_subject = sip_subject_start;
	    while (sip_subject)
	    {
	    	seqalign = SeqAlignNew();
	    	seqalign->type = 2;		/* alignment is diags */
	    	if (last == NULL)	/* First sequence. */
			seqalign_head = seqalign;
	    	else
			last->next = seqalign;

	    	last = seqalign;
		
		hspset_cnt_old = -1;
		hsp_cnt = results->hspcnt;
		for (index2=0; index2<hsp_cnt; index2++)
		{
			hsp = &(results->hsp_array[index2]);
			if (discontinuous && hspset_cnt_old != hsp->hspset_cnt)
			{
			    hspset_cnt_old = hsp->hspset_cnt;
			    if (index2 != 0)
			    { /* nothing to save on first pass. */
				if (getdensediag)
				{
					sap_head = FillInSegsInfo(sap_head, NULL, ddp_head);
					ddp_head = NULL;
				}
				else
				{
					sap_head = FillInSegsInfo(sap_head, ssp_head, NULL);
					ssp_head = NULL;
				}
			    }
			}

			if (reverse)
			{
				sip = SeqIdDup(sip_subject);
	    			sip->next = GetTheSeqAlignID(search->query_id);
			}
			else
			{
	    			sip = GetTheSeqAlignID(search->query_id);
				sip->next = SeqIdDup(sip_subject);
			}

			if (getdensediag)
			{
		    		ddp = FillInDenseDiagInfo(&ddp_head, hsp, reverse, search->context[hsp->context].query->length, subject_length, gi_list);
		    		ddp->id = sip;
			}
			else
			{
		    		ssp = FillInStdSegInfo(search, results, &ssp_head, hsp, sip, reverse, gi_list);
		    		ssp->ids = sip;
			}
			sip = NULL; /* This SeqIdPtr is now on the SeqAlign. */
		}

		if (discontinuous)
		{
			if (getdensediag)
			{
				sap_head = FillInSegsInfo(sap_head, NULL, ddp_head);
				ddp_head = NULL;
			}
			else
			{
				sap_head = FillInSegsInfo(sap_head, ssp_head, NULL);
				ssp_head = NULL;
			}
	        	seqalign->segs = sap_head;
	        	seqalign->segtype = 5;	/* Discontinuous */
		}
		else
		{
			if (getdensediag)
			{
				seqalign->segs = ddp_head;
                                seqalign->segtype = 1;  /* DenseDiag */
				ddp_head = NULL;
			}
			else
			{
                                seqalign->segs = ssp_head;
                                seqalign->segtype = 3;  /* StdSeg */
				ssp_head = NULL;
			}
		}

		sap_head = NULL;

		sip_subject = sip_subject->next;
	     }
	     if (sip_subject_start)
			sip_subject_start = SeqIdFree(sip_subject_start);
	}

	return seqalign_head;
}

/* Reverses the strand of a SeqAlignPtr (currently only of type DenseSegPtr). */
SeqAlignPtr LIBCALL
seqalign_reverse_strand (SeqAlignPtr salp)
{
  SeqAlignPtr salptmp;
  DenseSegPtr dsp;
  Int4Ptr     lenp;
  Int4Ptr     startp;
  Uint1Ptr    strandp;
  Int4        numseg;
  Int2        dim;
  Int4        j, k, n, tmp;

  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next) 
  {
     if (salptmp->segtype == 2) 
     {  
        dsp = (DenseSegPtr) salptmp->segs;
        strandp = dsp->strands;
        numseg = dsp->numseg;
        dim = dsp->dim;
        lenp = dsp->lens;
        startp = dsp->starts;
     }
     for (j=0; j < numseg*dim && strandp!=NULL; j++, strandp++) 
     {
           if (*strandp == Seq_strand_minus) 
              *strandp = Seq_strand_plus;
           else if (*strandp == Seq_strand_plus)
              *strandp = Seq_strand_minus;
     }
     for (j=0, k=numseg-1; j<numseg/2; j++, k--) {
           tmp=lenp[j];
           lenp[j]=lenp[k];
           lenp[k]=tmp;
     }
     for (j=0, k=(dim*numseg-dim); j<(dim*numseg-1)/2; j+=dim, k-=dim) {
           for (n=0; n<dim; n++) {
              tmp=startp[j+n];
              startp[j+n]=startp[k+n];
              startp[k+n]=tmp;
           }
     }
  }
  return salp;
}
/*
	"Core" function to compare two sequences, for use by 
	BlastTwoSequences and BlastSequencesOnTheFly.

	The subject_bsp is redundant with the subject_seq_start and
	subject_length (or visa-versa), but the subject must be
	extracted from the subject_bsp for BlastTwoSequences anyway, while
	the title and ID are needed from subject_bsp.
*/

static SeqAlignPtr 
BlastTwoSequencesCore (BlastSearchBlkPtr search, SeqLocPtr slp, Uint1Ptr subject_seq, Int4 subject_length, Boolean reverse)

{
	BLASTResultsStructPtr result_struct;
	BioseqPtr subject_bsp;
	Int2 status;
	Int4 index, hitlist_count;
	SeqAlignPtr seqalign=NULL;
	SeqIdPtr sip;
	SeqPortPtr spp;
	Uint1 residue;
	Uint1Ptr sequence;

	sip = SeqLocId(slp);
	subject_bsp = BioseqLockById(sip);

	search->subject_info = BLASTSubjectInfoDestruct(search->subject_info);
	search->subject_info = BLASTSubjectInfoNew(SeqIdDup(SeqIdFindBest(subject_bsp->id, SEQID_GI)), StringSave(BioseqGetTitle(subject_bsp)), subject_length);

	search = BLASTPerformSearch(search, subject_length, subject_seq);
	if (StringCmp(search->prog_name, "blastn") == 0 || search->pbp->gapped_calculation == FALSE)
	{
		status = BlastLinkHsps(search);
	}
	status = BlastReapHitlistByEvalue(search);
	status = BlastSaveCurrentHitlist(search);

	if (StringCmp(search->prog_name, "blastn") == 0 &&
                search->pbp->gapped_calculation == TRUE)
	{
		/* kbp_gap used in Traceback function. */
		search->sbp->kbp_gap[search->first_context] = search->sbp->kbp[search->first_context];
                result_struct = search->result_struct;
                hitlist_count = result_struct->hitlist_count;
		if (hitlist_count > 0)
		{
    			spp = SeqPortNewByLoc(slp, Seq_code_ncbi4na);
			/* make one longer to "protect" ALIGN. */
			sequence = MemNew((1+SeqLocLen(slp))*sizeof(Uint1));
			index=0;
			while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
			{
				if (residue == SEQPORT_VIRT)
					continue;
				sequence[index] = ncbi4na_to_blastna[residue];
				index++;
			}
			/* Gap character in last space. */
			sequence[index] = ncbi4na_to_blastna[0];
			seqalign = SumBlastGetGappedAlignmentTraceback(search, 0, reverse, FALSE, sequence, SeqLocLen(slp));
			sequence = MemFree(sequence);
    			spp = SeqPortFree(spp);
		}
		search->sbp->kbp_gap[search->first_context] = NULL;
	}
	else if (StringCmp(search->prog_name, "blastp") == 0 && 
		search->pbp->gapped_calculation == TRUE)
	{
                result_struct = search->result_struct;
                hitlist_count = result_struct->hitlist_count;
		if (hitlist_count > 0)
		{
			seqalign = BlastGetGapAlgnTbck(search, 0, reverse, FALSE, subject_seq, subject_length, NULL, 0);
		}
	}
	else
	{
		seqalign = GetSeqAlignForResultHitList(search, TRUE, FALSE, search->pbp->discontinuous, reverse, FALSE);
	}

	BioseqUnlock(subject_bsp);

	return seqalign;
}

/*
	Runs blast between two sequence, only works for blastp and blastn right
	now.
*/


SeqAlignPtr LIBCALL
BlastTwoSequencesByLocEx(SeqLocPtr slp1, SeqLocPtr slp2, CharPtr progname, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns)
{
	BlastAllWordPtr all_words;
	BlastSearchBlkPtr search;
	Boolean complement=FALSE, reverse, reverse_forbidden, options_alloc;
	Int4 index, subject_length, num_of_cols;
	SeqAlignPtr seqalign=NULL;
	SeqLocPtr query_slp, subject_slp;
	SeqPortPtr spp;
	SPCompressPtr spc=NULL;
	Uint1 residue;
	Uint1Ptr subject_seq, subject_seq_start;
	Uint1Ptr *array;

	if (slp1 == NULL || slp2 == NULL)
		return NULL;

	if (error_returns)
	{
		*error_returns = NULL;
	}

	if (other_returns)
	{
		*other_returns = NULL;
	}

	/* If filtering is performed, do not reverse the sequence.  In this case
		the wrong sequence would be filtered. */
	reverse_forbidden = FALSE;
	if (options && options->filter_string && StringCmp(options->filter_string, "F"))
	{
		reverse_forbidden = TRUE;
	}

	/* Select the shortest sequence as the query. */
	if (reverse_forbidden || SeqLocLen(slp1) < SeqLocLen(slp2))
	{
		query_slp = slp1;
		subject_slp = slp2;
		reverse = FALSE;
	}
	else
	{
		query_slp = slp2;
		subject_slp = slp1;
		reverse = TRUE;
	}

	if (progname == NULL && options == NULL)
		return NULL;

	if (progname == NULL)
	{
		progname = options->program_name;
	}

	/* If the subject strand is minus, turn it into plus for blastn. */
	/* Complement the other strand to keep things straight. */
	if (StringCmp(progname, "blastn") == 0 && SeqLocStrand(subject_slp) == Seq_strand_minus)
	{
		complement = TRUE;
                if(SeqLocStrand(query_slp) == Seq_strand_plus ||
			SeqLocStrand(query_slp) == Seq_strand_minus)
				SeqLocRevCmp(query_slp);
		SeqLocRevCmp(subject_slp);
	}

	subject_seq_start = subject_seq = NULL;

        /* Allocate default options if none are allocated yet. */
        options_alloc = FALSE;
        if (options == NULL)
        {
                options = BLASTOptionNew(progname, FALSE);
                options_alloc = TRUE;
        }

	all_words = NULL;

	subject_length = SeqLocLen(subject_slp);

	if (StringCmp(progname, "blastp") == 0)
	{
		subject_seq_start = (Uint1Ptr) MemNew(((subject_length)+2)*sizeof(Uint1));
		/* The first residue is the sentinel. */
		subject_seq_start[0] = NULLB;
		subject_seq = subject_seq_start+1;
		index = 0;
		spp = SeqPortNewByLoc(subject_slp, Seq_code_ncbistdaa);
		while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
		{
			subject_seq[index] = residue;
			index++;
		}
		subject_seq[index] = NULLB;

		num_of_cols = subject_length+1-options->wordsize;
		all_words = BlastAllWordNew(num_of_cols, options->wordsize, FALSE, TRUE);
		array = (Uint1Ptr *) MemNew(num_of_cols*sizeof(Uint1Ptr));
		for (index=0; index<num_of_cols; index++)
		{
			array[index] = subject_seq+index;
		}
		all_words->array = array;
		spp = SeqPortFree(spp);
		if (options->gapped_calculation == TRUE)
		{ 
			options->two_pass_method  = FALSE;
			options->multiple_hits_only  = TRUE;
		}
	}
	else if (StringCmp(progname, "blastn") == 0)
	{
		spp = SeqPortNewByLoc(subject_slp, Seq_code_ncbi4na);
		spc = SPCompressDNA(spp);
		if (spc == NULL)
			return NULL;
		subject_seq_start = subject_seq = spc->buffer;
		spp = SeqPortFree(spp);
	}
	else
	{
		return NULL;
	}

	search = BLASTSetUpSearchByLoc(query_slp, progname, SeqLocLen(query_slp), subject_length, all_words, options, NULL);

	if (search == NULL)
		return NULL;

	seqalign = BlastTwoSequencesCore(search, subject_slp, subject_seq, subject_length, reverse);

	if (complement)
	{
		seqalign = seqalign_reverse_strand(seqalign);
		SeqLocRevCmp(query_slp);
		SeqLocRevCmp(subject_slp);
	}

	if (spc)
	{
		SPCompressFree(spc);
		spc = NULL;
	}
	else
	{
		subject_seq_start = MemFree(subject_seq_start);
	}
	
	if (search->error_return)
	{
		ValNodeLink(error_returns, search->error_return);
		search->error_return = NULL;
	}

	if (other_returns)
	{ /* format dbinfo etc.  */
		*other_returns = BlastOtherReturnsPrepare(search);
	}


        if (options_alloc)
                options = BLASTOptionDelete(options);

	search = BlastSearchBlkDestruct(search);

	AdjustOffSetsInSeqAlign(seqalign, slp1, slp2);

	return seqalign;
}

SeqAlignPtr LIBCALL
BlastTwoSequencesByLoc(SeqLocPtr slp1, SeqLocPtr slp2, CharPtr progname, BLAST_OptionsBlkPtr options)
{
	return BlastTwoSequencesByLocEx(slp1, slp2, progname, options, NULL, NULL);
}

SeqAlignPtr LIBCALL
BlastTwoSequencesEx(BioseqPtr bsp1, BioseqPtr bsp2, CharPtr progname, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns)
{
	SeqAlignPtr seqalign;
	SeqLocPtr slp1=NULL, slp2=NULL;

	if (bsp1 == NULL || bsp2 == NULL)
		return NULL;

	slp1 = NULL;
	ValNodeAddPointer(&slp1, SEQLOC_WHOLE, SeqIdDup(SeqIdFindBest(bsp1->id, SEQID_GI)));

	slp2 = NULL;
	ValNodeAddPointer(&slp2, SEQLOC_WHOLE, SeqIdDup(SeqIdFindBest(bsp2->id, SEQID_GI)));

	seqalign = BlastTwoSequencesByLocEx(slp1, slp2, progname, options, other_returns, error_returns);

	slp1 = SeqLocFree(slp1);
	slp2 = SeqLocFree(slp2);

	return seqalign;
}

SeqAlignPtr LIBCALL
BlastTwoSequences(BioseqPtr bsp1, BioseqPtr bsp2, CharPtr progname, BLAST_OptionsBlkPtr options)
{
	return BlastTwoSequencesEx(bsp1, bsp2, progname, options, NULL, NULL);
}

/*
	Runs blast on the fly between the query BioseqPtr (specified with a
	call to BLASTSetUpSearch) and the subject BioseqPtr.
*/


SeqAlignPtr LIBCALL
BlastSequencesOnTheFlyByLoc(BlastSearchBlkPtr search, SeqLocPtr subject_slp)
{
	Int4 index, subject_length;
	SeqAlignPtr seqalign=NULL;
	SeqPortPtr spp;
	SPCompressPtr spc=NULL;
	Uint1Ptr subject_seq, subject_seq_start;
	Uint1 residue;

	if (subject_slp == NULL)
		return NULL;

	if (search == NULL || search->query_invalid)
		return NULL;

	if (search->result_struct)
		search->result_struct = BLASTResultsStructDelete(search->result_struct);
	search->result_struct = 
		BLASTResultsStructNew(search->result_size, search->pbp->max_pieces, search->pbp->hsp_range_max);

	BlastHitListPurge(search->current_hitlist);

	subject_seq_start = subject_seq = NULL;

	subject_length = SeqLocLen(subject_slp);

	if (StringCmp(search->prog_name, "blastp") == 0)
	{
		subject_seq_start = (Uint1Ptr) MemNew(((subject_length)+2)*sizeof(Uint1));
		/* The first residue is the sentinel. */
		subject_seq_start[0] = NULLB;
		subject_seq = subject_seq_start+1;
		index = 0;
		spp = SeqPortNewByLoc(subject_slp, Seq_code_ncbistdaa);
		while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
		{
			subject_seq[index] = residue;
			index++;
		}
		subject_seq[index] = NULLB;
		spp = SeqPortFree(spp);
	}
	else if (StringCmp(search->prog_name, "blastn") == 0)
	{
		spp = SeqPortNewByLoc(subject_slp, Seq_code_ncbi4na);
		spc = SPCompressDNA(spp);
		subject_seq = spc->buffer;
		spp = SeqPortFree(spp);
	}
	else
	{
		return NULL;
	}

	seqalign = BlastTwoSequencesCore(search, subject_slp, subject_seq, subject_length, FALSE);

	if (spc)
	{
		SPCompressFree(spc);
		spc = NULL;
	}
	else
	{
		subject_seq_start = MemFree(subject_seq_start);
	}
	
	AdjustOffSetsInSeqAlign(seqalign, search->query_slp, subject_slp);

	return seqalign;
}

SeqAlignPtr LIBCALL
BlastSequencesOnTheFly(BlastSearchBlkPtr search, BioseqPtr subject_bsp)
{
	SeqLocPtr slp;

	slp = NULL;
	ValNodeAddPointer(&slp, SEQLOC_WHOLE, SeqIdDup(SeqIdFindBest(subject_bsp->id, SEQID_GI)));
	return BlastSequencesOnTheFlyByLoc(search, slp);
}
/*
	Translate a nucleotide sequence without ambiguity codes.
	This is used for the first-pass translation of the database.
	
	BlastSearchBlkPtr search: overall BLAST structure.
	Int4 length: length of the nucl. sequence
	Uint1Ptr prot_seq: the (translated) protein sequence, with NULLB
		sentinels on either end.  This array should be allocated
		with sufficient memory before the function is called.
	Uint1Ptr nt_seq: the original nucl. sequence.
	
	The genetic code to be used is determined by the translation_table
	on the BlastSearchBlkPtr.

	This function translates a packed (ncbi2na) nucl. alphabet.  It
	views a basepair as being in one of four sets of 2-bits:

	|0|1|2|3||0|1|2|3||0|1|2|3||...

	1st byte | 2 byte | 3rd byte...

	A codon that starts at the beginning of the above sequence starts in
	state "0" and includes basepairs 0, 1, and 2.  The next codon, in the
	same frame, after that starts in state "3" and includes 3, 0, and 1.
*/

Int4 LIBCALL
BlastTranslateUnambiguousSequence(BlastSearchBlkPtr search, Int4 length, Uint1Ptr prot_seq, Uint1Ptr nt_seq, Int2 frame)

{
	register int state;
	Int2 total_remainder;
	Int4 prot_length;
	register Uint1 byte_value, codon;
	Uint1 last_remainder, last_byte, remainder;
	register Uint1Ptr translation, nt_seq_end, nt_seq_start;
	Uint1Ptr prot_seq_start;
  
	prot_length=0;
	if (nt_seq == NULL || prot_seq == NULL || (length-ABS(frame)+1) < CODON_LENGTH)
	return prot_length;

	*prot_seq = NULLB;
	prot_seq++;

/* record to determine protein length. */
	prot_seq_start = prot_seq;
  
	if (frame > 0)
		translation = search->translation_table;
	else
		translation = search->translation_table_rc;

	remainder = length%4;

	if (frame > 0)
	{
		nt_seq_end = nt_seq + (length)/4 - 1;
		last_remainder = (4*(length/4) - frame + 1)%CODON_LENGTH;
		total_remainder = last_remainder+remainder;
			
		state = frame-1;
		byte_value = *nt_seq;
		while (nt_seq < nt_seq_end)
		{
			switch (state)
			{
				case 0:
					codon = (byte_value >> 2);
					*prot_seq = translation[codon];
					prot_seq++;
				/* do state = 3 now, break is NOT missing. */
				case 3:
					codon = ((byte_value & 3) << 4);
					nt_seq++;
					byte_value = *nt_seq;	
					codon += (byte_value >> 4);
					*prot_seq = translation[codon];
					prot_seq++;
					if (nt_seq >= nt_seq_end)
					{
						state = 2;
						break;
					}
				/* Go on to state = 2 if not at end. */
				case 2:
					codon = ((byte_value & 15) << 2);
					nt_seq++;
					byte_value = *nt_seq;	
					codon += (byte_value >> 6);
					*prot_seq = translation[codon];
					prot_seq++;
					if (nt_seq >= nt_seq_end)
					{
						state = 1;
						break;
					}
				/* Go on to state = 1 if not at end. */
				case 1:
					codon = byte_value & 63;
					*prot_seq = translation[codon];
					prot_seq++;
					nt_seq++;
					byte_value = *nt_seq;	
					state = 0;
					break;
			}
		}

		if (state == 1)
		{ 
		/* This doesn't get done above, DON't do the state = 0
		   below if this is done. */
			byte_value = *nt_seq;
			codon = byte_value & 63;
			state = 0;
			*prot_seq = translation[codon];
			prot_seq++;
		}
		else if (state == 0)
		{ /* This one doesn't get done above. */
			byte_value = *nt_seq;
			codon = ((byte_value) >> 2);
			state = 3;
			*prot_seq = translation[codon];
			prot_seq++;
		}

		if (total_remainder >= CODON_LENGTH)
		{
			byte_value = *(nt_seq_end);
			last_byte = *(nt_seq_end+1);
			if (state == 0)
			{
				codon = (last_byte >> 2);
			}
			else if (state == 2)
			{
				codon = ((byte_value & 15) << 2);
				codon += (last_byte >> 6);
			}
			else if (state == 3)
			{
				codon = ((byte_value & 3) << 4);
				codon += (last_byte >> 4);
			}
			*prot_seq = translation[codon];
			prot_seq++;
		}
		*prot_seq = NULLB;
	}
	else
	{
		nt_seq_start = nt_seq;
		nt_seq += length/4;
		state = remainder+frame;
	/* Do we start in the last byte?  This one has the lowest order
	bits set to represent the remainder, hence the odd coding here. */
		if (state >= 0)
		{
			last_byte = *nt_seq;
			nt_seq--;
			if (state == 0)
			{
				codon = (last_byte >> 6);
				byte_value = *nt_seq;
				codon += ((byte_value & 15) << 2);
				state = 1;
			}
			else if (state == 1)
			{
				codon = (last_byte >> 4);
				byte_value = *nt_seq;
				codon += ((byte_value & 3) << 4);
				state = 2;
			}
			else if (state == 2)
			{
				codon = (last_byte >> 2);
				state = 3;
			}
			*prot_seq = translation[codon];
			prot_seq++;

		}
		else
		{
			state = 3 + (remainder + frame + 1);
			nt_seq--;
		}

		byte_value = *nt_seq;	
		while (nt_seq > nt_seq_start)
		{
			switch (state)
			{
				case 3:
					codon = (byte_value & 63);
					*prot_seq = translation[codon];
					prot_seq++;
				/* do state = 0 now, break is NOT missing. */
				case 0:
					codon = (byte_value >> 6);
					nt_seq--;
					byte_value = *nt_seq;	
					codon += ((byte_value & 15) << 2);
					*prot_seq = translation[codon];
					prot_seq++;
					if (nt_seq <= nt_seq_start)
					{
						state = 1;
						break;
					}
				/* Go on to state = 2 if not at end. */
				case 1:
					codon = (byte_value >> 4);
					nt_seq--;
					byte_value = *nt_seq;
					codon += ((byte_value & 3) << 4);
					*prot_seq = translation[codon];
					prot_seq++;
					if (nt_seq <= nt_seq_start)
					{
						state = 2;
						break;
					}
				/* Go on to state = 2 if not at end. */
				case 2:
					codon = (byte_value >> 2);
					*prot_seq = translation[codon];
					prot_seq++;
					nt_seq--;
					byte_value = *nt_seq;	
					state = 3;
					break;
			}
		}
		
		byte_value = *nt_seq;
		if (state == 3)
		{
			codon = (byte_value & 63);
			*prot_seq = translation[codon];
			prot_seq++;
		}
		else if (state == 2)
		{
			codon = (byte_value >> 2);
			*prot_seq = translation[codon];
			prot_seq++;
		}
	}

	*prot_seq = NULLB;

	return (prot_seq - prot_seq_start);
}	/* BlastTranslateUnambiguousSequence */



/*
	Gets an appropriate ID for the database (subject) sequence.
	Int4 hit_number is the index into the BLASTResultHitlistPtr,
	Boolean ordinal_number specifies whether an ordinal number (the
	db sequence number) or a real ID should be used.
*/
SeqIdPtr LIBCALL
BlastGetSubjectId(BlastSearchBlkPtr search, Int4 hit_number, Boolean ordinal_number, ValNodePtr *vnpp)

{
	BLASTResultHitlistPtr   results;
	DbtagPtr dbtagptr;
	ObjectIdPtr obidp;
	SeqIdPtr subject_id=NULL, sip;
	Uint4	header;

        results = search->result_struct->results[hit_number];
	if (ordinal_number)
	{

		obidp = ObjectIdNew();
                obidp->str = NULL;
                obidp->id = results->subject_id;
		dbtagptr = DbtagNew();
		if (search->rdfp)
		{
			dbtagptr->db = StringSave(search->rdfp->filename);
		}
		dbtagptr->tag = obidp;
                ValNodeAddPointer(&subject_id, SEQID_GENERAL, dbtagptr);
	} 
	else if (search->rdfp)
	{
		if (vnpp == NULL)
		{
           		readdb_get_descriptor(search->rdfp, results->subject_id, &subject_id, NULL);
		}
		else
		{
			header = 0;
			sip = NULL;
			while (readdb_get_header(search->rdfp, results->subject_id, &header, &sip, NULL) == TRUE)
				ValNodeAddPointer(vnpp, 0, sip);
			/* Get the last one. */
			ValNodeAddPointer(vnpp, 0, sip);
		}
	}
	else
	{
		if (results->subject_info)
			subject_id = SeqIdDup(results->subject_info->sip);
	}

	return subject_id;
}


/*
	Use by HeapSort (in BioseqBlastEngine) to rank Hitlist's.
*/

static int LIBCALLBACK
evalue_compare_hits(VoidPtr v1, VoidPtr v2)

{
        BLASTResultHitlistPtr h1, h2;
        BLASTResultHitlistPtr *hp1, *hp2;

        hp1 = (BLASTResultHitlistPtr *) v1;
        hp2 = (BLASTResultHitlistPtr *) v2;
        h1 = *hp1;
        h2 = *hp2;

	/* Sort first by evalue, then by score in case all evalues are zero. */
        if (h1->best_evalue < h2->best_evalue)
                return -1;
        if (h1->best_evalue > h2->best_evalue)
                return 1;
        if (h1->high_score > h2->high_score)
                return -1;
        if (h1->high_score < h2->high_score)
                return 1;
        return 0;
}

SeqAlignPtr LIBCALL
BioseqBlastEngineCore(BlastSearchBlkPtr search, BLAST_OptionsBlkPtr options, Int4Ptr *pos_matrix)
{
	BLASTResultHspPtr hsp;
	BLASTResultsStructPtr result_struct;
	BLASTResultHitlistPtr   result_hitlist;
	GapXEditBlockPtr edit_block;
	Int4 hitlist_count, hitlist_max, hspcnt, index, index1, length, sequence_length;
	SeqAlignPtr head, seqalign, seqalign_var;
	SeqIdPtr gi_list, subject_id;
	Uint1Ptr sequence;

	head = seqalign = NULL;

	if (search == NULL || search->query_invalid)
		return NULL;


	/* If pos_matrix is not NULL, then psi-blast iterations are being 
	performed.  The first psi-blast iteration should be with normal
	blast. */
	if (pos_matrix)
	{
		search->sbp->posMatrix = pos_matrix;
		search->positionBased = TRUE;
                search->sbp->kbp = search->sbp->kbp_psi;
                search->sbp->kbp_gap = search->sbp->kbp_gap_psi;
		hitlist_max = search->result_struct->hitlist_max;
                search->result_struct = BLASTResultsStructDelete(search->result_struct);
		search->result_struct = BLASTResultsStructNew(hitlist_max, search->pbp->max_pieces, search->pbp->hsp_range_max);
                if (search->allocated & BLAST_SEARCH_ALLOC_WFP_FIRST)
		{
                       search->wfp_first = BLAST_WordFinderDestruct(search->wfp_first);
		       search->wfp_first = BLAST_WordFinderNew(search->sbp->alphabet_size,options->wordsize,1, FALSE);
		}

		if (search->allocated & BLAST_SEARCH_ALLOC_WFP_SECOND)
		{
		       search->wfp_second = BLAST_WordFinderDestruct(search->wfp_second);
		       search->wfp_second = BLAST_WordFinderNew(search->sbp->alphabet_size,options->wordsize,1, FALSE);
		}

		/* Only find words once if thresholds are the same. */
                 search->wfp = search->wfp_first;
		 if (search->whole_query == TRUE)
                 	BlastNewFindWords(search, 0, search->context[search->first_context].query->length, search->pbp->threshold_first, (Uint1) 0);
		 else
                 	BlastNewFindWords(search, search->required_start, search->required_end, search->pbp->threshold_first, (Uint1) 0);
		 if (search->pbp->threshold_first != search->pbp->threshold_second)
		 {
                 	search->wfp = search->wfp_second;
			if (search->whole_query == TRUE)
                    		BlastNewFindWords(search, 0, search->context[search->first_context].query->length, search->pbp->threshold_second, (Uint1) 0);
			else
                    		BlastNewFindWords(search, search->required_start, search->required_end, search->pbp->threshold_second, (Uint1) 0);
		 }
		 else
		 {
			search->wfp_second = search->wfp_first;
		 }
	}

	/* Starting awake thread if multithreaded. */
	BlastStartAwakeThread(search);

	do_the_blast_run(search);
	
	head = NULL;
	if (StringCmp(search->prog_name, "blastn") == 0 && 
		search->pbp->gapped_calculation)
        {
		search->sbp->kbp_gap[search->first_context] = search->sbp->kbp[search->first_context];

		search->pbp->gap_open = options->gap_open;
		search->pbp->gap_extend = options->gap_extend;
		search->pbp->gap_x_dropoff = (BLAST_Score) (options->gap_x_dropoff*NCBIMATH_LN2 / search->sbp->kbp_gap[search->first_context]->Lambda);
		search->pbp->gap_x_dropoff_final = (BLAST_Score) (options->gap_x_dropoff_final*NCBIMATH_LN2 / search->sbp->kbp_gap[search->first_context]->Lambda);


		result_struct = search->result_struct;
       		hitlist_count = result_struct->hitlist_count;

		sequence=NULL;
		sequence_length=0;
		for (index=0; index<hitlist_count; index++)
		{
			length = readdb_get_sequence_ex(search->rdfp, result_struct->results[index]->subject_id, &sequence, &sequence_length);
			for (index1=0; index1<length; index1++)
			{
				sequence[index1] = ncbi4na_to_blastna[sequence[index1]];
			}
			/* Gap character in last space. */
			sequence[length] = ncbi4na_to_blastna[0];
			seqalign = SumBlastGetGappedAlignmentTraceback(search, index, FALSE, FALSE, sequence, length);
			result_struct->results[index]->seqalign = seqalign;
		}
		sequence = MemFree(sequence);
		search->sbp->kbp_gap[search->first_context] = NULL;

		HeapSort(result_struct->results, hitlist_count, sizeof(BLASTResultHitlistPtr), evalue_compare_hits);

		/* 
		The next loop organizes the SeqAligns (and the alignments in the
		BLAST report) in the same order as the deflines.
		*/
		head = NULL;
		for (index=0; index<hitlist_count; index++)
		{
		    seqalign = result_struct->results[index]->seqalign;
		    if (seqalign)
		    {
			if (head == NULL)
			{
				head = seqalign;
			}
			else
			{
				for (seqalign_var=head; seqalign_var->next;)
               		                seqalign_var = seqalign_var->next;
               		        seqalign_var->next = seqalign;
               	        }
		    }
		}
	}
	else if (search->pbp->gapped_calculation)
	{
		result_struct = search->result_struct;
                hitlist_count = result_struct->hitlist_count;
                for (index=0; index<hitlist_count; index++)
                {
                     seqalign = BlastGetGapAlgnTbckWithReaddb(search, index, FALSE);
		     result_struct->results[index]->seqalign = seqalign;
                }

		HeapSort(result_struct->results, hitlist_count, sizeof(BLASTResultHitlistPtr), evalue_compare_hits);

		/* 
		The next loop organizes the SeqAligns (and the alignments in the
		BLAST report) in the same order as the deflines.
		*/
		head = NULL;
		for (index=0; index<hitlist_count; index++)
		{
		    seqalign = result_struct->results[index]->seqalign;
		    if (seqalign)
		    {
			if (head == NULL)
			{
				head = seqalign;
			}
			else
			{
				for (seqalign_var=head; seqalign_var->next;)
               		                seqalign_var = seqalign_var->next;
               		        seqalign_var->next = seqalign;
               	        }
		    }
		}
	}
	else
	{
		/* Ungapped psi-blast. */
		if (StringCmp("blastp", search->prog_name) == 0 && search->positionBased == TRUE)
		{
		      result_struct = search->result_struct;
                      hitlist_count = result_struct->hitlist_count;
		      for (index=0; index<hitlist_count; index++)
		      {
			  result_hitlist = result_struct->results[index];
		      	  gi_list = BlastGetAllowedGis(search, result_hitlist->subject_id);
			  hspcnt = result_hitlist->hspcnt;
			  subject_id = BlastGetSubjectId(search, index, TRUE, NULL);
			  for (index1=0; index1<hspcnt; index1++)
			    {
			      hsp = &(result_hitlist->hsp_array[index1]);
			      edit_block = SimpleIntervalToGapXEditBlock(hsp->query_offset, hsp->subject_offset, hsp->query_length);
			      seqalign = GapXEditBlockToSeqAlign(edit_block, SeqIdDup(subject_id), SeqIdDup(search->query_id));
                              seqalign->score = GetScoreSetFromBlastResultHsp(hsp, gi_list);
			      
			      if (head == NULL)
				{
				  head = seqalign;
				}
			      else
				{
				  for (seqalign_var=head; seqalign_var->next;)
				    seqalign_var = seqalign_var->next;
				  seqalign_var->next = seqalign;
				}
			    }
			}
		}
		else
		{
		    if (StringCmp("blastn", search->prog_name) == 0 || StringCmp("blastp", search->prog_name) == 0)
			head = GetSeqAlignForResultHitList(search, TRUE, FALSE, options->discontinuous, FALSE, FALSE);
		    else
			head = GetSeqAlignForResultHitList(search, FALSE, FALSE, options->discontinuous, FALSE, FALSE);
		}
	}

	/* Stop the awake thread. */
	BlastStopAwakeThread();

	return head;
}

/*
	Deallocates all memory involved with the BlastHitRangePtr.
*/

BlastHitRangePtr LIBCALL
BlastHitRangeDestruct(BlastHitRangePtr old)

{
	if (old == NULL)
		return NULL;

	MemFree(old->range_list);
	MemFree(old->range_list_pointer);

	return MemFree(old);
}

/*
	Allocates a a BlastHitRangePtr, with two 'total' 
	BlastDoubleInt4Ptr's.
*/

BlastHitRangePtr LIBCALL
BlastHitRangeNew(Int4 total)

{
	BlastHitRangePtr bhrp;
	Int4 index;

	bhrp = MemNew(sizeof(BlastHitRange));

	bhrp->range_list = (BlastDoubleInt4Ptr) MemNew(total*sizeof(BlastDoubleInt4));
	bhrp->range_list_pointer = (BlastDoubleInt4Ptr PNTR) MemNew(total*sizeof(BlastDoubleInt4Ptr));
	for (index=0; index<total; index++)
	{
		bhrp->range_list_pointer[index] = &(bhrp->range_list[index]);
	}

	bhrp->current = 0;
	bhrp->total = total;

	return bhrp;
}

static int LIBCALLBACK
bhrp_compare(VoidPtr v1, VoidPtr v2)

{
	BlastDoubleInt4Ptr h1, h2;
	BlastDoubleInt4Ptr *hp1, *hp2;

	hp1 = (BlastDoubleInt4Ptr PNTR) v1;
	hp2 = (BlastDoubleInt4Ptr PNTR) v2;
	h1 = *hp1;
	h2 = *hp2;

	if (h1->gi < h2->gi)
		return -1;
	if (h1->gi > h2->gi)
		return 1;

	return 0;
}

BlastHitRangePtr LIBCALL
BioseqHitRangeEngineCore(BlastSearchBlkPtr search, BLAST_OptionsBlkPtr options)

{
	BioseqPtr subject_bsp;
	BlastHitRangePtr bhrp=NULL;
	BLASTResultsStructPtr result_struct;
	Int4 hitlist_count, index, index1, total_hsps;
	SeqPortPtr spp;
	Uint1Ptr sequence;

	if (search == NULL || search->query_invalid)
		return NULL;

	/* Starting awake thread if multithreaded. */
	BlastStartAwakeThread(search);

	do_the_blast_run(search);
	
	if (StringCmp(search->prog_name, "blastn") == 0 && 
		search->pbp->gapped_calculation)
        {
		search->sbp->kbp_gap[search->first_context] = search->sbp->kbp[search->first_context];

		search->pbp->gap_open = options->gap_open;
		search->pbp->gap_extend = options->gap_extend;
		search->pbp->gap_x_dropoff = (BLAST_Score) (options->gap_x_dropoff*NCBIMATH_LN2 / search->sbp->kbp_gap[search->first_context]->Lambda);
		search->pbp->gap_x_dropoff_final = (BLAST_Score) (options->gap_x_dropoff_final*NCBIMATH_LN2 / search->sbp->kbp_gap[search->first_context]->Lambda);


		result_struct = search->result_struct;
                hitlist_count = result_struct->hitlist_count;
		total_hsps = 0;
		for (index=0; index<hitlist_count; index++)
		{
			total_hsps += result_struct->results[index]->hspcnt;
		}
		bhrp = BlastHitRangeNew(total_hsps);
		
		result_struct = search->result_struct;
       		hitlist_count = result_struct->hitlist_count;

		for (index=0; index<hitlist_count; index++)
		{
			subject_bsp = readdb_get_bioseq(search->rdfp, result_struct->results[index]->subject_id);
    			spp = SeqPortNew(subject_bsp, 0, -1, Seq_strand_plus, Seq_code_ncbi4na);
			/* make one longer to "protect" ALIGN. */
			sequence = MemNew((1+subject_bsp->length)*sizeof(Uint1));
			for (index1=0; index1<subject_bsp->length; index1++)
			{
				sequence[index1] = ncbi4na_to_blastna[SeqPortGetResidue(spp)];
			}
			/* Gap character in last space. */
			sequence[subject_bsp->length] = ncbi4na_to_blastna[0];
			SumBlastGetGappedAlignmentEx(search, index, FALSE, FALSE, sequence, subject_bsp->length, FALSE, NULL, bhrp);
			sequence = MemFree(sequence);
			subject_bsp = BioseqFree(subject_bsp);
    			spp = SeqPortFree(spp);
		}
		search->sbp->kbp_gap[search->first_context] = NULL;
	}
	else if (search->pbp->gapped_calculation)
	{
		result_struct = search->result_struct;
                hitlist_count = result_struct->hitlist_count;
		/* FIll in. */
	}
	else
	{
/*	Change all this stuff. 
	   if (StringCmp("blastn", search->prog_name) == 0 || StringCmp("blastp", search->prog_name) == 0)
		head = GetSeqAlignForResultHitList(search, TRUE, FALSE, options->discontinuous, FALSE, FALSE);
	    else
		head = GetSeqAlignForResultHitList(search, FALSE, FALSE, options->discontinuous, FALSE, FALSE);
*/
	}

	HeapSort(bhrp->range_list_pointer, bhrp->current, sizeof(BlastHitRangePtr PNTR), bhrp_compare);

	/* Stop the awake thread. */
	BlastStopAwakeThread();

	return bhrp;
}

SeqAlignPtr LIBCALL
BioseqBlastEngineEx(BioseqPtr bsp, CharPtr progname, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, int (LIBCALLBACK *callback)PROTO((Int4 done, Int4 positives)), SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, Int4 gi_list_total)

{
	SeqLocPtr slp;

	slp = NULL;
	ValNodeAddPointer(&slp, SEQLOC_WHOLE, SeqIdDup(SeqIdFindBest(bsp->id, SEQID_GI)));
	return BioseqBlastEngineByLocEx(slp, progname, database, options, other_returns, error_returns, callback, seqid_list, gi_list, gi_list_total);
}

SeqAlignPtr LIBCALL
BioseqBlastEngine(BioseqPtr bsp, CharPtr progname, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, int (LIBCALLBACK *callback)PROTO((Int4 done, Int4 positives)))

{
	SeqLocPtr slp;

	slp = NULL;
	ValNodeAddPointer(&slp, SEQLOC_WHOLE, SeqIdDup(SeqIdFindBest(bsp->id, SEQID_GI)));
	return BioseqBlastEngineByLoc(slp, progname, database, options, other_returns, error_returns, callback);
}



SeqAlignPtr LIBCALL
BioseqBlastEngineByLoc(SeqLocPtr slp, CharPtr progname, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, int (LIBCALLBACK *callback)PROTO((Int4 done, Int4 positives)))

{
	return BioseqBlastEngineByLocEx(slp, progname, database, options, other_returns, error_returns, callback, NULL, NULL, 0);

}

SeqAlignPtr LIBCALL
BioseqBlastEngineByLocEx(SeqLocPtr slp, CharPtr progname, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, int (LIBCALLBACK *callback)PROTO((Int4 done, Int4 positives)), SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, Int4 gi_list_total)

{
	Boolean options_allocated=FALSE;
	BlastSearchBlkPtr search;
	Int2 status;
	SeqAlignPtr head, seqalign;

	head = seqalign = NULL;

	if (error_returns)
	{
		*error_returns = NULL;
	}

	if (other_returns)
	{
		*other_returns = NULL;
	}

	if (progname == NULL)
		return NULL;

	/* If no options, use default. */
	if (options == NULL)
	{
		options = BLASTOptionNew(progname, FALSE);
		options_allocated = TRUE;
	}

	status = BLASTOptionValidateEx(options, progname, error_returns);
	if (status != 0)
	{	/* error messages in other_returns? */
		return NULL;
	}

	if (slp == NULL || database == NULL)
		return NULL;

	search = BLASTSetUpSearchByLocWithReadDbEx(slp, progname, SeqLocLen(slp), database, options, NULL, seqid_list, gi_list, gi_list_total);

	if (search == NULL)
	{
		return NULL;
	}

	search->tick_callback = callback;
	search->star_callback = callback;

	head = BioseqBlastEngineCore(search, options, NULL);
	
	if (search->error_return)
	{
		ValNodeLink(error_returns, search->error_return);
		search->error_return = NULL;
	}

	if (other_returns)
	{ /* format dbinfo etc.  */
		*other_returns = BlastOtherReturnsPrepare(search);
	}

	if (options_allocated)
	{
		options = BLASTOptionDelete(options);
	}
	search = BlastSearchBlkDestruct(search);

	return head;
}

SeqLocPtr LIBCALL
BioseqHitRangeEngine(BioseqPtr bsp, CharPtr progname, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, int (LIBCALLBACK *callback)PROTO((Int4 done, Int4 positives)), SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, Int4 gi_list_total)

{
	SeqLocPtr slp;

	slp = NULL;
	ValNodeAddPointer(&slp, SEQLOC_WHOLE, SeqIdDup(SeqIdFindBest(bsp->id, SEQID_GI)));
	return BioseqHitRangeEngineByLoc(slp, progname, database, options, other_returns, error_returns, callback, seqid_list, gi_list, gi_list_total);
}

SeqLocPtr 
HitRangeToSeqLoc(BlastHitRangePtr bhrp, Int4 link_value)

{
	Boolean make_seqloc, start=TRUE;
	Int4 index, total, start_pos, stop_pos, largest_stop_pos;
	SeqIntPtr sint;
	SeqLocPtr retval=NULL;

	if (bhrp == NULL)
		return NULL;

	total = bhrp->current;
	index=0;
	while (index < total)
	{
		if (start == TRUE)
		{
			start_pos = bhrp->range_list_pointer[index]->gi;
			start = FALSE;
			largest_stop_pos = 0;
		}
		else
		{
			/* Keep track of largest stop position. */
			largest_stop_pos = MAX(largest_stop_pos, bhrp->range_list_pointer[index]->ordinal_id);
			make_seqloc = FALSE;
			if (index == total-1)	/* Last one. */
			{
				stop_pos = bhrp->range_list_pointer[index]->ordinal_id;
				start = TRUE;
				make_seqloc = TRUE;
			}
			else if (largest_stop_pos+link_value < bhrp->range_list_pointer[index+1]->gi)
			{ /* Check overlap with next one. */
				stop_pos = bhrp->range_list_pointer[index]->ordinal_id;
				start = TRUE;
				make_seqloc = TRUE;
			}
			
			if (make_seqloc)
			{
				sint = SeqIntNew();
				sint->from = start_pos;
				sint->to = MAX(largest_stop_pos, stop_pos);
				sint->strand = Seq_strand_plus;
				ValNodeAddPointer(&retval, SEQLOC_INT, sint);
			}
			index++;
		}
	}

	return retval;
}

#define HITRANGE_LINKVALUE 5

SeqLocPtr LIBCALL
BioseqHitRangeEngineByLoc(SeqLocPtr slp, CharPtr progname, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, int (LIBCALLBACK *callback)PROTO((Int4 done, Int4 positives)), SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, Int4 gi_list_total)

{
	Boolean options_allocated=FALSE;
	BlastHitRangePtr bhrp;
	BlastSearchBlkPtr search;
	Int2 status;
	SeqLocPtr seqloc;

	if (error_returns)
	{
		*error_returns = NULL;
	}

	if (other_returns)
	{
		*other_returns = NULL;
	}

	if (progname == NULL)
		return NULL;

	/* If no options, use default. */
	if (options == NULL)
	{
		options = BLASTOptionNew(progname, FALSE);
		options_allocated = TRUE;
	}

	status = BLASTOptionValidateEx(options, progname, error_returns);
	if (status != 0)
	{	/* error messages in other_returns? */
		return NULL;
	}

	if (slp == NULL || database == NULL)
		return NULL;

	search = BLASTSetUpSearchByLocWithReadDbEx(slp, progname, SeqLocLen(slp), database, options, NULL, seqid_list, gi_list, gi_list_total);

	if (search == NULL)
	{
		return NULL;
	}

	search->tick_callback = callback;
	search->star_callback = callback;

	bhrp = BioseqHitRangeEngineCore(search, options);
	
	seqloc = HitRangeToSeqLoc(bhrp, HITRANGE_LINKVALUE);
	bhrp = BlastHitRangeDestruct(bhrp);
	if (search->error_return)
	{
		ValNodeLink(error_returns, search->error_return);
		search->error_return = NULL;
	}

	if (other_returns)
	{ /* format dbinfo etc.  */
		*other_returns = BlastOtherReturnsPrepare(search);
	}

	if (options_allocated)
	{
		options = BLASTOptionDelete(options);
	}
	search = BlastSearchBlkDestruct(search);

	return seqloc;
}


ValNodePtr LIBCALL
BlastOtherReturnsPrepare(BlastSearchBlkPtr search)

{
	BLAST_KarlinBlkPtr ka_params;
	BLAST_MatrixPtr blast_matrix;
	CharPtr parameters;
	ReadDBFILEPtr rdfp_var;
	TxDfDbInfoPtr dbinfo, head, dbinfo_var;
	ValNodePtr other_returns=NULL;
	
	head = NULL;
	if (search->blast_gi_list)
	{
		dbinfo = MemNew(sizeof(TxDfDbInfo));
		dbinfo->total_length = search->dblen;
		dbinfo->number_seqs = search->dbseq_num;
		dbinfo->subset = TRUE;
		head = dbinfo;
		dbinfo_var = dbinfo;
	}

	rdfp_var = search->rdfp;
	while (rdfp_var)
	{
		dbinfo = MemNew(sizeof(TxDfDbInfo));
		dbinfo->name = StringSave(readdb_get_filename(rdfp_var));	
		dbinfo->definition = StringSave(readdb_get_title(rdfp_var));	
		dbinfo->date = StringSave(readdb_get_date(rdfp_var));	
		dbinfo->total_length = readdb_get_dblen(rdfp_var);
		dbinfo->number_seqs = readdb_get_num_entries(rdfp_var);
		if (head == NULL)
		{
			head = dbinfo;
			dbinfo_var = dbinfo;
		}
		else
		{
			dbinfo_var->next = dbinfo;
			dbinfo_var = dbinfo_var->next;
		}
		rdfp_var = rdfp_var->next;
	}
	ValNodeAddPointer (&other_returns, TXDBINFO, head);

	if (search->sbp->kbp)
	{
		ka_params = BlastKarlinBlkCreate();
		ka_params->Lambda = search->sbp->kbp[search->first_context]->Lambda;
		ka_params->K = search->sbp->kbp[search->first_context]->K;
		ka_params->H = search->sbp->kbp[search->first_context]->H;
		ValNodeAddPointer (&other_returns, TXKABLK_NOGAP, ka_params);
	}

	if (search->pbp->gapped_calculation == TRUE)
	{
		if (StringCmp(search->prog_name, "blastn") == 0)
		{
			if (search->sbp->kbp)
			{
				ka_params = BlastKarlinBlkCreate();
				ka_params->Lambda = search->sbp->kbp[search->first_context]->Lambda;
				ka_params->K = search->sbp->kbp[search->first_context]->K;
				ka_params->H = search->sbp->kbp[search->first_context]->H;
				ValNodeAddPointer (&other_returns, TXKABLK_GAP, ka_params);
			}
		}
		else
		{
			if (search->sbp->kbp_gap)
			{
				ka_params = BlastKarlinBlkCreate();
				ka_params->Lambda = search->sbp->kbp_gap[search->first_context]->Lambda;
				ka_params->K = search->sbp->kbp_gap[search->first_context]->K;
				ka_params->H = search->sbp->kbp_gap[search->first_context]->H;
				ValNodeAddPointer (&other_returns, TXKABLK_GAP, ka_params);
			}
		}
	}

	parameters = FormatBlastParameters(search);
	ValNodeAddPointer (&other_returns, TXPARAMETERS, parameters);

	blast_matrix = BLAST_MatrixFill(search->sbp, search->positionBased);
	ValNodeAddPointer (&other_returns, TXMATRIX, blast_matrix);

	if (search->mask)
		ValNodeLink(&other_returns, search->mask);

	return other_returns;
}


/*
	Deallocates memory for BLAST_ExtendWordParamsPtr
	
*/

static BLAST_ExtendWordParamsPtr
BLAST_ExtendWordParamsDestruct (BLAST_ExtendWordParamsPtr ewp_params)

{
	ewp_params = MemFree(ewp_params);

	return ewp_params;
}


/*
	Allocates memory for the BLAST_ExtendWordParamsPtr.  

	This function also sets many of the parametes such as min_diag_length etc.

	Int4 qlen: length of the query.
	Boolean multiple_hits: specifies whether multiple hits method is used.
*/

static BLAST_ExtendWordParamsPtr
BLAST_ExtendWordParamsNew (Int4 qlen, Boolean multiple_hits)

{
	BLAST_ExtendWordParamsPtr ewp_params;
	Int4 min_diag_length, bits_to_shift;

	ewp_params= MemNew(sizeof(BLAST_ExtendWordParams));

	if (ewp_params)
	{
		min_diag_length = 1;
		bits_to_shift = 0;
		/* What power of 2 is just longer than the query? */
		while (min_diag_length < qlen)
		{
			min_diag_length = min_diag_length << 1;
			bits_to_shift++;
		}
		/* These are used in the word finders to shift and mask 
		rather than dividing and taking the remainder. */
		ewp_params->bits_to_shift = bits_to_shift;
		ewp_params->min_diag_length = min_diag_length;
		ewp_params->min_diag_mask = min_diag_length-1;
		ewp_params->multiple_hits = multiple_hits;
	}
	return ewp_params;
}

/*
	Deallocates memory for the BLAST_ExtendWordPtr.

*/
static BLAST_ExtendWordPtr LIBCALL 
BLAST_ExtendWordDestruct (BLAST_ExtendWordPtr ewp)

{
	if (ewp)
	{
		if (ewp->_buffer)
			ewp->_buffer = MemFree(ewp->_buffer);

		ewp = MemFree(ewp);
	}

	return ewp;

}

/*
	Allocates memory for the BLAST_ExtendWordPtr.  

	All of the memory for the arrays is allocated in one chunk
	called "_buffer".  If multiple_hits is specified them room
	for "diag_level, "last_hit", and "version" is allocated and
	pointers into the array for these are set.  If multiple_hits
	is not set, then only room for diag_level and version is allocated;
	last_hit is not needed.

	Int4 qlen, dblen: length of the query and the LONGEST subject sequence.
	Boolean multiple_hits: specifies whether multiple hits method is used.

*/
static BLAST_ExtendWordPtr
BLAST_ExtendWordNew (BLAST_ExtendWordParamsPtr ewp_params)

{
	BLAST_ExtendWordPtr ewp;

	ewp = MemNew(sizeof(BLAST_ExtendWord));

	if (ewp)
	{
		/* Allocate the buffer to be used for all 3 arrays. */
		if (ewp_params->multiple_hits)
		{
                    ewp->_buffer = (Int4Ptr) MemNew(3*ewp_params->min_diag_length*sizeof(Int4));
		}
		else
		{
                    ewp->_buffer = (Int4Ptr) MemNew(2*ewp_params->min_diag_length*sizeof(Int4));
		}
		if (ewp->_buffer == NULL)
		{
			ewp = BLAST_ExtendWordDestruct(ewp);
			return NULL;
		}

		/* Set up the pointers to the arrays. */
		ewp->diag_level = ewp->_buffer;
		if (ewp_params->multiple_hits)
		{
			ewp->last_hit = ewp->_buffer + ewp_params->min_diag_length;
			ewp->version = ewp->_buffer + 2*ewp_params->min_diag_length;
		}
		else
		{
			ewp->version = ewp->_buffer + ewp_params->min_diag_length;
		}
	}

	return ewp;
}

/*****************************************************************************
*
*	Zeroe's out the memory in the array _buffer, if offset is greater than
*	INT4_MAX/2.  The first "min_diag_length" spaces in the array are used 
*	by the array "diag_level", the second "min_diag_length" spaces are used 
*	by "last_hit".  All of these are zeroed out.  The last "min_diag_length" 
*	spaces are used by "version; these are not zeroed out.
*
*	If offset is not greater than INT4_MAX/2, then the memory is not
*	zeroed out.  Rather "offset" is used as a "zero-point" that is
*	always greater than the next possible value when the word finder
*	starts working on a new subject sequence.
*
******************************************************************************/
void LIBCALL
BlastExtendWordExit(BlastSearchBlkPtr search)

{
	BLAST_ExtendWordPtr ewp;
	BLAST_ExtendWordParamsPtr ewp_params;
	Int2 index;

	ewp_params = search->ewp_params;

	for (index=search->first_context; index<=search->last_context; index++)
	{

		if (ewp_params->offset >= INT4_MAX/2)
		{
			ewp = search->context[index].ewp;
			if (ewp_params->multiple_hits == TRUE)
			{
				Nlm_MemSet(ewp->_buffer, 0, 2*ewp_params->min_diag_length*sizeof(ewp->_buffer[0]));
			}
			else
			{
				Nlm_MemSet(ewp->_buffer, 0, ewp_params->min_diag_length*sizeof(ewp->_buffer[0]));
			}
		}
	}

	if (ewp_params->offset < INT4_MAX/2)
	{
		ewp_params->offset += search->subject->length + ewp_params->window;
	}
	else
	{
		ewp_params->offset = 0;
	}
}


static BlastSequenceBlkPtr
BlastSequenceBlkDestruct(BlastSequenceBlkPtr seq_blk)

{

	if (seq_blk == NULL)
		return NULL;

	/* Free from the start of sequence if it's filled in. */
	if (seq_blk->sequence_start != NULL)
	{
		seq_blk->sequence_start = MemFree(seq_blk->sequence_start);
	}
	else
	{	
		seq_blk->sequence = MemFree(seq_blk->sequence);
	}

	seq_blk = MemFree(seq_blk);

	return seq_blk;
}



static BLASTContextStructPtr 
BLASTContextFree(BLASTContextStructPtr context, Int2 number)

{
	Int2 index;
	
	for (index=0; index<number; index++)
	{
		context[index].ewp = BLAST_ExtendWordDestruct(context[index].ewp);
		if (context[index].query_allocated == TRUE)
		{
			context[index].query = BlastSequenceBlkDestruct(context[index].query);
		}
	}
	context = MemFree(context);

	return context;
}

/* 
	Allocates space for a copy of the BlastSearchBlk for use in
	multi-processing BLAST.
*/

BlastSearchBlkPtr LIBCALL
BlastSearchBlkDuplicate (BlastSearchBlkPtr search)

{

	BlastSearchBlkPtr new_search;
	Int2 index;

	if (search == NULL)
		return NULL;

	new_search = (BlastSearchBlkPtr) MemNew(sizeof(BlastSearchBlk));
	if (new_search == NULL)
		return NULL;

	/* What's allocated here? */
	new_search->allocated = 0;	
	new_search->allocated += BLAST_SEARCH_ALLOC_SUBJECT;
	new_search->allocated += BLAST_SEARCH_ALLOC_PBP;
	new_search->allocated += BLAST_SEARCH_ALLOC_CONTEXT;
	new_search->allocated += BLAST_SEARCH_ALLOC_READDB;
	new_search->allocated += BLAST_SEARCH_ALLOC_EWPPARAMS;
		
	/* Duplicate the rfdp struct, but not the contents. */
	new_search->rdfp = readdb_attach(search->rdfp);
	if (new_search->rdfp == NULL)
	{
		new_search = BlastSearchBlkDestruct(new_search);
		return NULL;
	}

	new_search->positionBased = search->positionBased;

	if (search->blast_gi_list)
		new_search->blast_gi_list = MemDup(search->blast_gi_list, sizeof(BlastGiList));

	/* Changes, need to allocate. */
	new_search->pbp = MemDup(search->pbp, sizeof(BLAST_ParameterBlk));
	new_search->sbp = search->sbp;
	new_search->wfp_first = search->wfp_first;
	new_search->wfp_second = search->wfp_second;
	new_search->prog_name = StringSave(search->prog_name);
	new_search->prog_number = search->prog_number;
	new_search->first_context = search->first_context;
	new_search->last_context = search->last_context;
	new_search->ewp_params = MemDup(search->ewp_params, sizeof(BLAST_ExtendWordParams));
	new_search->dblen = search->dblen;
	new_search->dblen_eff = search->dblen_eff;
	new_search->dblen_eff_real = search->dblen_eff_real;
	new_search->dbseq_num = search->dbseq_num;
	new_search->length_adjustment = search->length_adjustment;
	new_search->searchsp_eff = search->searchsp_eff;

	/* Allocate last_context+1 elements, even if there are only last_context-first_context
	being used. */
	new_search->context = (BLASTContextStructPtr) MemNew((search->last_context+1)*sizeof(BLASTContextStruct));
	for (index=new_search->first_context; index<=new_search->last_context; index++)
	{
		new_search->context[index].ewp = BLAST_ExtendWordNew(new_search->ewp_params);
		new_search->context[index].query = search->context[index].query;
		new_search->context[index].query->frame = ContextToFrame(new_search, index);
		new_search->context[index].query_allocated = FALSE;
	}

	new_search->context_factor = search->context_factor;

	new_search->subject = (BlastSequenceBlkPtr) MemNew(sizeof(BlastSequenceBlk));
	/* 100 is the size limit in the present BLAST for hsp's. */
	new_search->hsp_array_size = search->hsp_array_size;
	/* The results are held here. */
	new_search->result_struct = search->result_struct;

	new_search->worst_evalue = DBL_MAX;

	new_search->translation_table = search->translation_table;
	new_search->translation_table_rc = search->translation_table_rc;
	new_search->genetic_code = search->genetic_code;
	new_search->db_genetic_code = search->db_genetic_code;

	if (search->translation_buffer_size > 0)
	{	/* two extra for the NULLB's on end. */
		new_search->translation_buffer = MemNew((2+search->translation_buffer_size)*sizeof(Uint1));
		new_search->translation_buffer_size = search->translation_buffer_size;
	}

	new_search->gap_align = NULL;	/* Allocated automatically. */

	new_search->whole_query = search->whole_query;
	new_search->required_start = search->required_start;
	new_search->required_end = search->required_end;

	new_search->blast_seqid_list = search->blast_seqid_list;

#ifdef BLAST_COLLECT_STATS
	new_search->first_pass_hits = 0;
	new_search->second_pass_hits = 0;
	new_search->second_pass_trys = 0;
	new_search->first_pass_extends = 0;
	new_search->second_pass_extends = 0;
	new_search->first_pass_good_extends = 0;
	new_search->second_pass_good_extends = 0;
	new_search->number_of_seqs_better_E = 0;
	new_search->prelim_gap_no_contest = 0;
	new_search->prelim_gap_passed = 0;
	new_search->prelim_gap_attempts = 0;
	new_search->real_gap_number_of_hsps = 0;
#endif

	return new_search;
}

/* 
	Allocates space for the new BlastSearchBlk and some sturctures
	attached to it.
*/

BlastSearchBlkPtr LIBCALL
BlastSearchBlkNew (Int2 wordsize, Int4 qlen, CharPtr dbname, Boolean multiple_hits, BLAST_Score threshold_first, BLAST_Score threshold_second, Int4 result_size, CharPtr prog_name, BlastAllWordPtr all_words, Int2 first_context, Int2 last_context)

{

	BlastSearchBlkPtr search;
	BLASTContextStructPtr context;
	Uint1 is_prot;
	Int2 index;
	Uint1 alphabet;
	Int4 longest_db_seq=INT4_MAX;

	search = (BlastSearchBlkPtr) MemNew(sizeof(BlastSearchBlk));

	if (search != NULL)
	{
		search->allocated = 0;	/* everything's allocated here. */
		search->allocated += BLAST_SEARCH_ALLOC_QUERY;
		search->allocated += BLAST_SEARCH_ALLOC_SUBJECT;
		search->allocated += BLAST_SEARCH_ALLOC_PBP;
		search->allocated += BLAST_SEARCH_ALLOC_SBP;
		search->allocated += BLAST_SEARCH_ALLOC_EWPPARAMS;
		search->allocated += BLAST_SEARCH_ALLOC_CONTEXT;
		search->allocated += BLAST_SEARCH_ALLOC_RESULTS;
		search->allocated += BLAST_SEARCH_ALLOC_READDB;
		search->allocated += BLAST_SEARCH_ALLOC_ALL_WORDS;
		
		search->positionBased = FALSE;

		if (StringCmp(prog_name, "blastn") == 0)
		{
			alphabet = BLASTNA_SEQ_CODE;
		}
		else
		{
			alphabet = Seq_code_ncbistdaa;
		}

		if (dbname != NULL)
		{
			
			if (StringCmp(prog_name, "blastp") == 0 || StringCmp(prog_name, "blastx") == 0)
			{ /* Protein DB for blastp and blastx. */
				is_prot = READDB_DB_IS_PROT;
			}
			else
			{
				is_prot = READDB_DB_IS_NUC;
			}
			
			if ((search->rdfp=readdb_new(dbname, is_prot)) == NULL)
			{
				return NULL;
			}

			longest_db_seq = readdb_get_maxlen(search->rdfp);
		}

		search->first_context = first_context;
		search->last_context = last_context;

		search->pbp = 
		   (BLAST_ParameterBlkPtr) MemNew(sizeof(BLAST_ParameterBlk));

		search->sbp = BLAST_ScoreBlkNew(alphabet, last_context+1);

		/* Only allocate these if thresholds are above zero, i.e. they will be used. */
		if (StringCmp(prog_name, "blastn") != 0)
		{
			if (threshold_first > 0)
			{
				search->wfp_first = BLAST_WordFinderNew(search->sbp->alphabet_size, wordsize, 1, FALSE);
				search->allocated += BLAST_SEARCH_ALLOC_WFP_FIRST;
			}
		
			if (threshold_second > 0)
			{
		/* Only allocate a new WFP if 2nd th differs from 1st. */
				if (threshold_second != threshold_first)
				{
					search->wfp_second = BLAST_WordFinderNew(search->sbp->alphabet_size, wordsize, 1, FALSE);
					search->allocated += BLAST_SEARCH_ALLOC_WFP_SECOND;
				}
				else
				{
					search->wfp_second = search->wfp_first;
				}
			}
		}
		else
		{
			if (multiple_hits)
				search->wfp_second = BLAST_WordFinderNew(256, wordsize, READDB_COMPRESSION_RATIO, FALSE);
			else
				search->wfp_second = BLAST_WordFinderNew(256, wordsize, READDB_COMPRESSION_RATIO, TRUE);
			search->allocated += BLAST_SEARCH_ALLOC_WFP_SECOND;
		}

		search->prog_name = StringSave(prog_name);
		search->prog_number = BlastGetProgramNumber(prog_name);

		search->ewp_params = BLAST_ExtendWordParamsNew(MIN(qlen, longest_db_seq), multiple_hits);

		context = search->context = (BLASTContextStructPtr) MemNew((1+search->last_context)*sizeof(BLASTContextStruct));
		for (index=0; index<=last_context; index++)
		{
			context[index].ewp = BLAST_ExtendWordNew(search->ewp_params);
			context[index].query = (BlastSequenceBlkPtr) MemNew(sizeof(BlastSequenceBlk));
			context[index].query->frame = ContextToFrame(search, index);
			context[index].query_allocated = TRUE;
		}

		search->subject = (BlastSequenceBlkPtr) MemNew(sizeof(BlastSequenceBlk));
		/* 100 is the size limit in the present BLAST for hsp's. */
		search->hsp_array_size = 100;
		/* The results are held here. */
		search->result_size = result_size;
/*
		search->result_struct = BLASTResultsStructNew(result_size, search->pbp->max_pieces, search->pbp->hsp_range_max);
*/

		search->worst_evalue = DBL_MAX;

		search->whole_query = TRUE;
		search->required_start = 0;
		search->required_end = -1;

		search->all_words = all_words;

#ifdef BLAST_COLLECT_STATS
		search->first_pass_hits = 0;
		search->second_pass_hits = 0;
		search->second_pass_trys = 0;
		search->first_pass_extends = 0;
		search->second_pass_extends = 0;
		search->first_pass_good_extends = 0;
		search->second_pass_good_extends = 0;
		search->number_of_seqs_better_E = 0;
		search->prelim_gap_no_contest = 0;
		search->prelim_gap_passed = 0;
		search->prelim_gap_attempts = 0;
		search->real_gap_number_of_hsps = 0;
#endif
	}

	return search;
}

/*
	Deallocates memory associated with the BlastSearchBlkPtr.
*/

BlastSearchBlkPtr LIBCALL 
BlastSearchBlkDestruct (BlastSearchBlkPtr search)

{

	if (search != NULL)
	{
		if (search->allocated & BLAST_SEARCH_ALLOC_QUERY)
			search->original_seq = MemFree(search->original_seq);

		if (search->allocated & BLAST_SEARCH_ALLOC_SUBJECT)
			search->subject = BlastSequenceBlkDestruct(search->subject);

		if (search->allocated & BLAST_SEARCH_ALLOC_PBP)
			search->pbp = MemFree(search->pbp);

		if (search->allocated & BLAST_SEARCH_ALLOC_SBP)
			search->sbp = BLAST_ScoreBlkDestruct(search->sbp);

		if (search->allocated & BLAST_SEARCH_ALLOC_WFP_FIRST)
			search->wfp_first = BLAST_WordFinderDestruct(search->wfp_first);

		if (search->allocated & BLAST_SEARCH_ALLOC_WFP_SECOND)
		{
			search->wfp_second = BLAST_WordFinderDestruct(search->wfp_second);
		}

		if (search->allocated & BLAST_SEARCH_ALLOC_EWPPARAMS)
		{
			search->ewp_params = BLAST_ExtendWordParamsDestruct(search->ewp_params);
		}

		if (search->allocated & BLAST_SEARCH_ALLOC_CONTEXT)
		{
			search->context = BLASTContextFree(search->context, 1+search->last_context);
		}

		if (search->allocated & BLAST_SEARCH_ALLOC_RESULTS)
		{
			search->result_struct = BLASTResultsStructDelete(search->result_struct);
		}

		if (search->allocated & BLAST_SEARCH_ALLOC_READDB)
		{
			search->rdfp = readdb_destruct(search->rdfp);
		}

		if (search->current_hitlist)
		{
			search->current_hitlist = BlastHitListDestruct(search->current_hitlist);
		}
		search->subject_info = BLASTSubjectInfoDestruct(search->subject_info);


		if (search->prog_name)
		{
			search->prog_name = MemFree(search->prog_name);
		}

		if (search->query_id)
		{
			search->query_id = SeqIdFree(search->query_id);
		}

		if (search->translation_buffer_size > 0)
		{
			search->translation_buffer = MemFree(search->translation_buffer);
		}

		if (search->allocated & BLAST_SEARCH_ALLOC_TRANS_INFO)
		{

			if (search->translation_table)
			{
				search->translation_table = MemFree(search->translation_table);
			}

			if (search->translation_table_rc)
			{
				search->translation_table_rc = MemFree(search->translation_table_rc);
			}
		}

		if (search->allocated & BLAST_SEARCH_ALLOC_ALL_WORDS)
		{
			search->all_words = BlastAllWordDestruct(search->all_words);
		}

		if (search->blast_gi_list)
			MemFree(search->blast_gi_list);

		BlastSeqIdListDestruct(search->blast_seqid_list);

		search->gap_align = GapAlignBlkDelete(search->gap_align);

		if (search->allocated & BLAST_SEARCH_ALLOC_QUERY_SLP)
		{
			if (search->query_slp)
				search->query_slp = SeqLocFree(search->query_slp);
		}

		search = MemFree(search);
	}

	return search;
}


/* 
	Deallocates all the memory associated with the BlastAllWordPtr.
*/

BlastAllWordPtr LIBCALL
BlastAllWordDestruct(BlastAllWordPtr all_words)

{
	Int4 index1, num_of_cols;
	Uint1Ptr *array;

	if (all_words == NULL)
		return NULL;

	if (all_words->array)
	{
		array = all_words->array;
		num_of_cols = all_words->num_of_cols;

		if (all_words->rows_allocated)
		{
			for (index1=0; index1<num_of_cols; index1++)
			{
				array[index1] = MemFree(array[index1]);
			}
		}
		array = MemFree(array);
	}

	MemFree(all_words);

	return NULL;
}

/*
	Allocates the BlastAllWordPtr and sets some flags.
*/
BlastAllWordPtr LIBCALL
BlastAllWordNew(Int4 num_of_cols, Int4 wordsize, Boolean rows_allocated, Boolean specific)

{
	BlastAllWordPtr all_words;

	all_words = MemNew(sizeof(BlastAllWord));
	if (all_words)
	{
		all_words->rows_allocated = rows_allocated;
		all_words->specific = specific;
		all_words->num_of_cols = num_of_cols;
		all_words->wordsize = wordsize;
	}

	return all_words;
}

BLAST_HitListPtr LIBCALL
BlastHitListDestruct(BLAST_HitListPtr hitlist)
{
        BLAST_HSPPtr PNTR hsp_array;
        Int4 hspcnt_max, index;

        if (hitlist == NULL)
                return NULL;

        hspcnt_max = hitlist->hspcnt_max;
        hsp_array = hitlist->hsp_array;

        for (index=0; index<hspcnt_max; index++)
        {
                hsp_array[index] = MemFree(hsp_array[index]);
        }

        hitlist->hsp_array = MemFree(hsp_array);
        hitlist = MemFree(hitlist);

        return hitlist;
}

/****************************************************************

        Functions to allocate and destroy the BLAST_HitList.

***************************************************************/
BLAST_HitListPtr LIBCALL
BlastHitListNew(BlastSearchBlkPtr search)
{
        BLAST_HitListPtr hitlist;

        hitlist = (BLAST_HitListPtr) MemNew(sizeof(BLAST_HitList));

        if (hitlist == NULL)
                return hitlist;

        hitlist->hspmax = search->hsp_array_size;
        hitlist->hsp_array = (BLAST_HSPPtr PNTR) MemNew(hitlist->hspmax*sizeof
(BLAST_HSPPtr));

        if (hitlist->hsp_array == NULL)
        {
                hitlist = BlastHitListDestruct(hitlist);
                return NULL;
        }

        return hitlist;
}


/*
	This function translates the context number of a context into
	the frame of the sequence.

	Arguments:
	
	BlastSearchBlkPtr search: search structure,
	Int2 context_number: context number used by BLASTContextStruct array
	Boolean is_query: if TRUE, refers to query, otherwise the subject.
*/

static Int2
ContextToFrame(BlastSearchBlkPtr search, Int2 context_number)

{
	Int2 frame = 0;

	if (search->prog_number == blast_type_blastp || search->prog_number == blast_type_tblastn)
	{	/* Query and subject are protein, no frame. */
		frame = 0;
	}
	if (search->prog_number == blast_type_blastx || search->prog_number == blast_type_tblastx)
	{
		frame = context_number < 3 ? context_number-3 : context_number-2;
	}
	else if (search->prog_number == blast_type_blastn)
	{
		if (context_number == 0)
			frame = 1;
		else
			frame = -1;
	}

	return frame;
}

/*
	Allocates and fills in the BLASTSubjectInfo structure.
*/

BLASTSubjectInfoPtr LIBCALL
BLASTSubjectInfoNew(SeqIdPtr sip, CharPtr defline, Int4 length)

{
	BLASTSubjectInfoPtr subject_info;

	subject_info = (BLASTSubjectInfoPtr) MemNew(sizeof(BLASTSubjectInfo));	

	if (subject_info == NULL)
		return NULL;

	subject_info->sip = sip;
	subject_info->defline = defline;
	subject_info->length = length;

	return subject_info;
}

/*
	Deallocates the BLASTSubjectInfo structure and the
	SeqIdPtr, as well as the defline.
*/

BLASTSubjectInfoPtr LIBCALL
BLASTSubjectInfoDestruct(BLASTSubjectInfoPtr subject_info)

{

	if (subject_info == NULL)
		return NULL;

	SeqIdFree(subject_info->sip);
	MemFree(subject_info->defline);
	subject_info = MemFree(subject_info);

	return subject_info;
}



/*
	Destroys BLASTResultsStructure and associated memory.
*/

BLASTResultsStructPtr LIBCALL
BLASTResultsStructDelete(BLASTResultsStructPtr result_struct)

{
	Int4 index;
	BLASTResultHitlistPtr PNTR results;
	BLASTHeapPtr hp, hpt;
	
	results = result_struct->results;
	for (index=0; index<result_struct->hitlist_max; index++)
	{
		if (results[index])
		{
			results[index] = BLASTResultHitlistFree(results[index]);
		}
	}


	for (hp = result_struct->heap_ptr; hp; ) 
	{
	  hpt = hp->next;
	  hp->heap = MemFree(hp->heap);
	  hp = MemFree(hp);
	  hp = hpt;
	}
	result_struct->results = MemFree(result_struct->results);
	result_struct = MemFree(result_struct);

	return result_struct;
}

/*
	returns BLASTResultsStruct.
*/

BLASTResultsStructPtr
BLASTResultsStructNew(Int4 results_size, Int4 max_pieces, Int4 range_max)

{
	BLASTResultsStructPtr new;
	Int4 index;

	new = MemNew(sizeof(BLASTResultsStruct));
	new->results = (BLASTResultHitlistPtr PNTR) MemNew(results_size*sizeof(BLASTResultHitlistPtr));

	for (index=0; index<results_size; index++)
		new->results[index] = NULL;

	new->hitlist_max = results_size;
	new->hitlist_count = 0;
	new->max_pieces = max_pieces;
	new->heap_ptr = (BLASTHeapPtr) MemNew(sizeof(BLASTHeapStruct));
	new->heap_ptr->cutvalue = INT4_MAX;
	new->heap_ptr->num_in_heap = new->heap_ptr->num_of_ref = 0;
	new->heap_ptr->prev = new->heap_ptr->next = NULL;
	new->heap_ptr->heap = (BLASTResultHspPtr PNTR) MemNew(sizeof(BLASTResultHspPtr)*range_max);

	return new;
}


Uint1 AAForCodon (Uint1Ptr codon, CharPtr codes);

/*
	GetTranslation to get the translation of the nucl. sequence in the
	appropriate frame and with the appropriate GeneticCode.

	The function return an allocated CharPtr, the caller must delete this.
	The first and last spaces of this CharPtr contain NULLB's.
*/

Uint1Ptr LIBCALL
GetTranslation(Uint1Ptr query_seq, Int4 nt_length, Int2 frame, Int4Ptr length, CharPtr genetic_code)
{
	Uint1 codon[CODON_LENGTH];
	Int4 index, index_prot;
	SeqMapTablePtr smtp;
	Uint1 residue;
	Uint1Ptr prot_seq;

	smtp = SeqMapTableFind(Seq_code_ncbistdaa, Seq_code_ncbieaa);

	/* Allocate two extra spaces for NULLB's at beginning and end of seq. */
	prot_seq = (Uint1Ptr) MemNew((2+(nt_length+2)/CODON_LENGTH)*sizeof(Uint1));

	/* The first character in the protein is the NULLB sentinel. */
	prot_seq[0] = NULLB;
	index_prot = 1;
	for (index=ABS(frame)-1; index<nt_length-2; index += CODON_LENGTH)
	{
		codon[0] = query_seq[index];
		codon[1] = query_seq[index+1];
		codon[2] = query_seq[index+2];
		residue = AAForCodon(codon, genetic_code);
		prot_seq[index_prot] = SeqMapTableConvert(smtp, residue);
		index_prot++;
	}
	prot_seq[index_prot] = NULLB;
	*length = index_prot-1;
	
	return prot_seq;
}




/*************************************************************************
*
*	MaskTheResidues masks up to max_length residues in buffer.
*	The residue to be used for masking (generally 'N' for nucleotides
*	and 'X' for proteins) is mask_residue.  offset tells how far
*	along the sequence the first residue in buffer is.  mask_slp
*	specifies which parts of the sequence to mask.  'max_length is
*	the total length of the sequence.
*
*************************************************************************/

void
BlastMaskTheResidues(Uint1Ptr buffer, Int4 max_length, Uint1 mask_residue, SeqLocPtr mask_slp, Boolean reverse)

{
	SeqLocPtr slp=NULL;
        Int4 index, start, stop;
       
	while (mask_slp)
	{
		slp=NULL;
        	while((slp = SeqLocFindNext(mask_slp, slp))!=NULL)
        	{
			if (reverse)
			{
				start = max_length - 1 - SeqLocStop(slp);
				stop = max_length - 1 - SeqLocStart(slp);
			}
			else
			{
              			start = SeqLocStart(slp);
              			stop = SeqLocStop(slp);
			}

			for (index=start; index<=stop; index++)
			{
				buffer[index] = mask_residue;
			}
        	}
		mask_slp = mask_slp->next;
	}

}

Boolean LIBCALL FilterDNA(BioseqPtr bsp, Int4 filter)
{
  SeqLocPtr slp = NULL, dust_slp=NULL, dust_slp_start=NULL;
  Uint1 mask_residue = 'N';
  Int4 index, start, stop;
  Uint1 oldcode;
  CharPtr buffer;

  switch(filter) {
    
  case FILTER_DUST:
    oldcode = bsp->seq_data_type;
    dust_slp = BioseqDust(bsp, 0, -1, -1, -1, -1, -1); 
    dust_slp_start = dust_slp;
    BioseqConvert(bsp, Seq_code_iupacna);
    buffer = MemNew(bsp->length);
    BSSeek(bsp->seq_data, 0, SEEK_SET);    
    BSRead(bsp->seq_data, buffer, bsp->length);
    while (dust_slp) {
      slp=NULL;
      while((slp = SeqLocFindNext(dust_slp, slp))!=NULL) {
        start = SeqLocStart(slp);
        stop = SeqLocStop(slp);
        for (index=0; index < bsp->length; index++) {
          if (index >= start) {
            if (index <= stop)
              buffer[index] = mask_residue;
            else if (index > stop)
              break;
          }
        }
      }
      dust_slp = dust_slp->next;
    }
    dust_slp_start = SeqLocSetFree(dust_slp_start);
    BSSeek(bsp->seq_data, 0, SEEK_SET);
    BSWrite(bsp->seq_data, buffer, bsp->length);
    BioseqConvert(bsp, oldcode);
    MemFree(buffer);
    break;

  default:
    return FALSE;  /* wrong type of filter used */
  }
  return TRUE;
}

/*
	COnverts a protein (translated) SeqLocPtr from the protein
	coordinates to the nucl. coordinates.

	Only works on a SeqLocPtr of type SeqIntPtr right now.
*/

Boolean
BlastConvertProteinSeqLoc(SeqLocPtr slp, Int2 frame, Int4 full_length)

{
	SeqIntPtr seq_int;
	Int4 from, to;

	if (slp == NULL)
		return TRUE;

	if (slp->choice == SEQLOC_PACKED_INT)
		slp = slp->data.ptrvalue;

	while (slp)
	{
		if (slp->choice != SEQLOC_INT)
			return FALSE;

		seq_int = slp->data.ptrvalue;
		from = seq_int->from;
		to = seq_int->to;

		if (frame < 0)
		{
			seq_int->to = full_length - CODON_LENGTH*from + frame;
			seq_int->from = full_length - CODON_LENGTH*to + frame + 1;
			seq_int->strand = Seq_strand_minus;
		}
		else
		{
			seq_int->from = CODON_LENGTH*from + frame - 1;
			seq_int->to = CODON_LENGTH*to + frame - 1;
			seq_int->strand = Seq_strand_plus;
		}
		slp = slp->next;
	}
	
	return TRUE;
}

SeqLocPtr
BioseqSegEx(BioseqPtr bsp_unfilter, CharPtr options)

{
	BioseqPtr bsp_filter;
	Boolean mask_state;
	Char cmd_buf[2*PATH_MAX], temp_file[PATH_MAX];
	CharPtr filter_dir;
	Int4 index, mask_begin;
	SeqEntryPtr sep;
	SeqLocPtr slp_mask;
	SeqPortPtr spp_filter, spp_unfilter;
	Uint1 res_filter, res_unfilter;
	FILE *fp;


	if (bsp_unfilter == NULL)
		return NULL;

#ifdef OS_UNIX

	TmpNam(temp_file);
	fp = FileOpen(temp_file, "w");
	if (BioseqToFasta(bsp_unfilter, fp, FALSE) == FALSE)
	{
		BioseqUnlock(bsp_unfilter);
		FileClose(fp);
		return NULL;
	}
	FileClose(fp);

	filter_dir = getenv("BLASTFILTER");
	if (filter_dir == NULL)
		filter_dir = BLASTFILTER_DIR;

	if (options != NULL)
		sprintf(cmd_buf, "%s%s%s%s %s%s", filter_dir, DIRDELIMSTR, "seg ", temp_file, options, " -x");
	else
		sprintf(cmd_buf, "%s%s%s%s%s", filter_dir, DIRDELIMSTR, "seg ", temp_file, " -x");

	fp = popen(cmd_buf, "r");
	if (fp == NULL)
	{
		ErrPostEx(SEV_WARNING, 0, 0, "Call to seg failed.");
		return NULL;
	}
	
	sep = FastaToSeqEntry(fp, FALSE);
	FileClose(fp);
	if (sep == NULL)
	{
		ErrPostEx(SEV_WARNING, 0, 0, "Call to seg failed.");
		return NULL;
	}
	bsp_filter = sep->data.ptrvalue;

	spp_filter = SeqPortNew(bsp_filter, 0, -1, Seq_strand_plus, Seq_code_ncbistdaa);
	spp_unfilter = SeqPortNew(bsp_unfilter, 0, -1, Seq_strand_plus, Seq_code_ncbistdaa);

	mask_state = FALSE;
	index = 0;
	slp_mask = NULL;
	while ((res_filter=SeqPortGetResidue(spp_filter)) != SEQPORT_EOF)
	{
		res_unfilter=SeqPortGetResidue(spp_unfilter);
		if (res_filter != res_unfilter)
		{
			if (mask_state == FALSE)
			{
				mask_begin = index;
				mask_state = TRUE;
			}
		}
		else if (mask_state == TRUE)
		{
			ValNodeLink(&slp_mask, SeqLocIntNew(mask_begin, index-1, Seq_strand_plus, bsp_filter->id));
			mask_state = FALSE;
		}
		index++;
	}

	/* If the last portion of the sequence was masked. */
	if (mask_state == TRUE)
	{
		ValNodeLink(&slp_mask, SeqLocIntNew(mask_begin, index-1, Seq_strand_plus, bsp_filter->id));
	}

	sep = SeqEntryFree(sep);
	SeqPortFree(spp_filter);
	SeqPortFree(spp_unfilter);

	pclose(fp);
	FileRemove(temp_file);

	return slp_mask;
#else
	return NULL;
#endif
}

/*
	Runs seg and obtains a SeqLocPtr from it.
*/
static SeqLocPtr 
SeqLocSegEx(SeqLocPtr slp, CharPtr instructions)

{
	BioseqPtr bsp_unfilter;
	SeqLocPtr slp_mask;
	SeqIdPtr sip;


	if (slp == NULL)
		return NULL;

	sip = SeqIdFindBest(SeqLocId(slp), SEQID_GI);
	bsp_unfilter = BioseqLockById(sip);
	slp_mask = BioseqSegEx(bsp_unfilter, instructions);

	BioseqUnlock(bsp_unfilter);

	return slp_mask;
}

SeqLocPtr 
SeqLocSeg(SeqLocPtr slp)

{
	return SeqLocSegEx(slp, NULL);
}

SeqLocPtr
MyBioseqSeg(BioseqPtr bsp_unfilter)

{
	return BioseqSegEx(bsp_unfilter, NULL);
}

#define BLASTSEQLOC_BUFFER_SIZE 32

static Boolean
parse_dust_options(CharPtr ptr, Int4Ptr level, Int4Ptr window, Int4Ptr cutoff, Int4Ptr linker)

{
	Char buffer[BLASTSEQLOC_BUFFER_SIZE];
	Int4 arg, index, index1, window_pri, linker_pri, level_pri, cutoff_pri;

	arg = 0;
	index1 = 0;
	for (index=0; index<BLASTSEQLOC_BUFFER_SIZE; index++)
	{
		if (*ptr == ' ' || *ptr == NULLB)
		{
			buffer[index1] = NULLB;
			index1 = 0;
			switch(arg) {
				case 0:
					sscanf(buffer, "%ld", &level_pri);
					break;
				case 1:
					sscanf(buffer, "%ld", &window_pri);
					break;
				case 2:
					sscanf(buffer, "%ld", &cutoff_pri);
					break;
				case 3:
					sscanf(buffer, "%ld", &linker_pri);
					break;
				default:
					break;
			}

			arg++;
			while (*ptr == ' ')
				ptr++;

			/* end of the buffer. */
			if (*ptr == NULLB)
				break;
		}
		else
		{
			buffer[index1] = *ptr; ptr++;
			index1++;
		}
	}

	*level = level_pri; 
	*window = window_pri; 
	*cutoff = cutoff_pri; 
	*linker = linker_pri; 

	return TRUE;
}


static Boolean
parse_seg_options(CharPtr ptr, Int4Ptr window, FloatHiPtr locut, FloatHiPtr hicut)

{
	Char buffer[BLASTSEQLOC_BUFFER_SIZE];
	Int4 arg, index, index1, window_pri; 
	FloatHi locut_pri, hicut_pri;

	arg = 0;
	index1 = 0;
	for (index=0; index<BLASTSEQLOC_BUFFER_SIZE; index++)
	{
		if (*ptr == ' ' || *ptr == NULLB)
		{
			buffer[index1] = NULLB;
			index1 = 0;
			switch(arg) {
				case 0:
					sscanf(buffer, "%ld", &window_pri);
					break;
				case 1:
					sscanf(buffer, "%le", &locut_pri);
					break;
				case 2:
					sscanf(buffer, "%le", &hicut_pri);
					break;
				default:
					break;
			}

			arg++;
			while (*ptr == ' ')
				ptr++;

			/* end of the buffer. */
			if (*ptr == NULLB)
				break;
		}
		else
		{
			buffer[index1] = *ptr; ptr++;
			index1++;
		}
	}

	*window = window_pri; 
	*locut = locut_pri; 
	*hicut = hicut_pri; 

	return TRUE;
}

static Boolean
parse_cc_options(CharPtr ptr, Int4Ptr window, FloatHiPtr cutoff, Int4Ptr linker)

{
	Char buffer[BLASTSEQLOC_BUFFER_SIZE];
	Int4 arg, index, index1, window_pri, linker_pri;
	Nlm_FloatHi cutoff_pri;

	arg = 0;
	index1 = 0;
	for (index=0; index<BLASTSEQLOC_BUFFER_SIZE; index++)
	{
		if (*ptr == ' ' || *ptr == NULLB)
		{
			buffer[index1] = NULLB;
			index1 = 0;
			switch(arg) {
				case 0:
					sscanf(buffer, "%ld", &window_pri);
					break;
				case 1:
					sscanf(buffer, "%le", &cutoff_pri);
					break;
				case 2:
					sscanf(buffer, "%ld", &linker_pri);
					break;
				default:
					break;
			}

			arg++;
			while (*ptr == ' ')
				ptr++;

			/* end of the buffer. */
			if (*ptr == NULLB)
				break;
		}
		else
		{
			buffer[index1] = *ptr; ptr++;
			index1++;
		}
	}

	*window = window_pri; 
	*cutoff = cutoff_pri; 
	*linker = linker_pri; 

	return TRUE;
}

static CharPtr
load_options_to_buffer(CharPtr instructions, CharPtr buffer)

{
	Boolean not_started=TRUE;
	CharPtr buffer_ptr, ptr;
	Int4 index;

	ptr = instructions;
	buffer_ptr = buffer;
	for (index=0; index<BLASTSEQLOC_BUFFER_SIZE && *ptr != NULLB; index++)
	{
		if (*ptr == ';')
		{
			ptr++;
			break;
		}
		/* Remove blanks at the beginning. */
		if (not_started && *ptr == ' ')
		{
			ptr++;
		}
		else
		{
			not_started = FALSE;
			*buffer_ptr = *ptr;
			buffer_ptr++; ptr++;
		}
	}

	*buffer_ptr = NULLB;

	if (not_started == FALSE)
	{	/* Remove trailing blanks. */
		buffer_ptr--;
		while (*buffer_ptr == ' ' && buffer_ptr > buffer)
		{
			*buffer_ptr = NULLB;
			buffer_ptr--;
		}
	}

	return ptr;
}

#define CC_WINDOW 22
#define CC_CUTOFF 40.0
#define CC_LINKER 32

/*
	This function parses the 'instructions' string and then calls the appopriate
	filtering functions.
*/
SeqLocPtr
BlastBioseqFilter(BioseqPtr bsp, CharPtr instructions)

{
	return BlastBioseqFilterEx(bsp, instructions, NULL);
}

SeqLocPtr
BlastBioseqFilterEx(BioseqPtr bsp, CharPtr instructions, BoolPtr mask_at_hash)

{
	BLAST_OptionsBlkPtr options;
	Boolean do_all=FALSE, do_seg=FALSE, do_coil_coil=FALSE, do_dust=FALSE, do_repeats=FALSE;
	Char buffer[BLASTSEQLOC_BUFFER_SIZE];
	CharPtr ptr;
	Int2 seqloc_num;
	Int4 window_cc, linker_cc, window_dust, level_dust, minwin_dust, linker_dust;
        SeqLocPtr cc_slp=NULL, dust_slp=NULL, seg_slp=NULL, seqloc_head=NULL, repeat_slp=NULL;
        PccDatPtr pccp;
        Nlm_FloatHiPtr scores;
	Nlm_FloatHi cutoff_cc;
	SegParamsPtr sparamsp=NULL;

	cutoff_cc = CC_CUTOFF;

	if (instructions == NULL || StringICmp(instructions, "F") == 0)
		return NULL;

	/* FALSE is the default right now. */
	if (mask_at_hash)
		*mask_at_hash = FALSE;

	if (StringICmp(instructions, "T") == 0)
	{
		do_all = TRUE;
	}
	else
	{
		ptr = instructions;
		if (*ptr == 'm')
		{
			if (mask_at_hash)
				*mask_at_hash = TRUE;
			ptr += 2;
		}
		while (*ptr != NULLB)
		{
			if (*ptr == 'S')
			{
				sparamsp = SegParamsNewAa();
				ptr = load_options_to_buffer(ptr+1, buffer);
				if (buffer[0] != NULLB)
				{
					parse_seg_options(buffer, &sparamsp->window, &sparamsp->locut, &sparamsp->hicut);
				}
				do_seg = TRUE;
			}
			else if (*ptr == 'C')
			{
				ptr = load_options_to_buffer(ptr+1, buffer);
				window_cc = CC_WINDOW;
				cutoff_cc = CC_CUTOFF;
				linker_cc = CC_LINKER;
				if (buffer[0] != NULLB)
					parse_cc_options(buffer, &window_cc, &cutoff_cc, &linker_cc);
				do_coil_coil = TRUE;
			}
			else if (*ptr == 'D')
			{
				ptr = load_options_to_buffer(ptr+1, buffer);
				/* -1 indicates defaults. */
				level_dust = -1;
				window_dust = -1;
				minwin_dust = -1;
				linker_dust = -1;
				if (buffer[0] != NULLB)
					parse_dust_options(buffer, &level_dust, &window_dust, &minwin_dust, &linker_dust);
				do_dust = TRUE;
			}
			else if (*ptr == 'R')
			{
				ptr = load_options_to_buffer(ptr+1, buffer);
				do_repeats = TRUE;
			}
			else
			{	/* Nothing applied. */
				ptr++;
			}
		}
	}

	seqloc_num = 0;
	seqloc_head = NULL;
	if (ISA_aa(bsp->mol))
	{
		if (do_all || do_seg)
		{
			seg_slp = BioseqSeg(bsp, sparamsp);
			SegParamsFree(sparamsp);
			sparamsp=NULL;
			seqloc_num++;
		}
		if (do_coil_coil)
		{
        		pccp = PccDatNew ();
     		   	pccp->window = window_cc;
        		ReadPccData (pccp);
        		scores = PredictCCBioseq(bsp, 0, bsp->length-1, pccp);
        		cc_slp = FilterCC(scores, cutoff_cc, bsp->length, linker_cc, SeqIdDup(bsp->id));
        		MemFree(scores);
        		PccDatFree (pccp);
			seqloc_num++;
		}
	}
	else
	{
		if (do_dust)
		{
			dust_slp = BioseqDust(bsp, 0, -1, level_dust, window_dust, minwin_dust, linker_dust);
			seqloc_num++;
		}
		if (do_repeats)
		{
			options = BLASTOptionNew("blastn", TRUE);
			options->expect_value = 0.01;
			repeat_slp = BioseqHitRangeEngine(bsp, "blastn", "humlines.lib  humsines.lib  retrovir.lib", options, NULL, NULL, NULL, NULL, NULL, 0);
			options = BLASTOptionDelete(options);
			seqloc_num++;
		}
	}

	if (seqloc_num == 0)
	{ /* nothing. */
		;
	} 
	else if (seqloc_num == 1)
	{
		if (seg_slp)
			seqloc_head = seg_slp;
		if (cc_slp)
			seqloc_head = cc_slp;
		if (dust_slp)
			seqloc_head = dust_slp;
		if (repeat_slp)
                        seqloc_head = repeat_slp;
	}
	else
	{
		if (seg_slp)
			ValNodeAddPointer(&seqloc_head, SEQLOC_MIX, seg_slp);
		if (cc_slp)
			ValNodeAddPointer(&seqloc_head, SEQLOC_MIX, cc_slp);
		if (dust_slp)
			ValNodeAddPointer(&seqloc_head, SEQLOC_MIX, dust_slp);
		if (repeat_slp)
                        ValNodeAddPointer(&seqloc_head, SEQLOC_MIX, repeat_slp);
	}

	return seqloc_head;
}

SeqLocPtr
BlastSeqLocFilter(SeqLocPtr slp, CharPtr instructions)

{

	return BlastSeqLocFilterEx(slp, instructions, NULL);

}

SeqLocPtr
BlastSeqLocFilterEx(SeqLocPtr slp, CharPtr instructions, BoolPtr mask_at_hash)

{
	BioseqPtr bsp;
	SeqIdPtr sip;
	SeqLocPtr slp_mask;

	sip = SeqIdFindBest(SeqLocId(slp), SEQID_GI);
	bsp = BioseqLockById(sip);
	slp_mask = BlastBioseqFilterEx(bsp, instructions, mask_at_hash);
	BioseqUnlock(bsp);

	return slp_mask;
}

/*
	Program to run seg on a sequence.  Note that this program only
	really works in UNIX systems.
*/
Boolean LIBCALL
FilterWithSeg (Uint1Ptr sequence, Int4 length, Uint1 alphabet)

{

#ifdef OS_UNIX

	BioseqPtr bsp;
	Char cmd_buf[2*PATH_MAX], temp_file[PATH_MAX];
	CharPtr filter_dir;
	FILE PNTR fp;
	Int4 byte_store_length;
	Nlm_ByteStorePtr byte_store;
	SeqEntryPtr sep;

	if (sequence == NULL || length == 0)
		return FALSE;

	byte_store = Nlm_BSNew(length);

	byte_store_length = Nlm_BSWrite(byte_store, (VoidPtr) sequence, length);
	if (length != byte_store_length)
	{
		Nlm_BSDelete(byte_store, length);
		return FALSE;
	}

	bsp = BioseqNew();
	bsp->seq_data = byte_store;
	bsp->length = length;
	bsp->seq_data_type = alphabet;
	bsp->mol = Seq_mol_aa;
	bsp->repr = Seq_repr_raw;

	TmpNam(temp_file);
	fp = FileOpen(temp_file, "w");
	if (BioseqToFasta(bsp, fp, FALSE) == FALSE)
	{
		bsp = BioseqFree(bsp);
		return FALSE;
	}
	FileClose(fp);

	bsp = BioseqFree(bsp);

	filter_dir = getenv("BLASTFILTER");
	if (filter_dir != NULL)
		sprintf(cmd_buf, "%s%s%s%s%s", filter_dir, DIRDELIMSTR, "seg ", temp_file, " -x");
	else
		sprintf(cmd_buf, "%s%s%s%s%s", BLASTFILTER_DIR, DIRDELIMSTR, "seg ", temp_file, " -x");

	fp = popen(cmd_buf, "r");
	if (fp == NULL)
	{
		ErrPostEx(SEV_WARNING, 0, 0, "Call to seg failed.");
		return FALSE;
	}
	
	sep = FastaToSeqEntry(fp, FALSE);
	if (sep == NULL)
	{
		ErrPostEx(SEV_WARNING, 0, 0, "Call to seg failed.");
		return FALSE;
	}

	pclose(fp);

	bsp = sep->data.ptrvalue;
	BioseqRawConvert(bsp, Seq_code_ncbistdaa);

	BSSeek(bsp->seq_data, 0, SEEK_SET);
	Nlm_BSRead(bsp->seq_data, (VoidPtr) sequence, length);

	SeqEntryFree(sep);

	FileRemove(temp_file);

	return TRUE;
#else
	return FALSE;
#endif
}



BLASTResultHitlistPtr LIBCALL
BLASTResultHitlistFree(BLASTResultHitlistPtr result)

{
  BLASTResultHspPtr hsp;
  register Int2 i;
	
	if (result == NULL)
		return NULL;

        for(i=0; i < result->hspcnt; i++) {
          hsp = &result->hsp_array[i];
	  if (hsp)
                GapXEditBlockDelete(hsp->gap_info);
        }

	if (result->hspcnt != 0)
		result->hsp_array = MemFree(result->hsp_array);

	BLASTSubjectInfoDestruct(result->subject_info);	

	result = MemFree(result);

	return result;
}

/*
	Creates a new BLASTResultHitlist, with the an hsp-array of length hspcnt.  If the
	allocation fails, then NULL is returned.
*/

BLASTResultHitlistPtr LIBCALL
BLASTResultHitlistNew(Int4 hspcnt)

{

	BLASTResultHitlistPtr new;

	new = (BLASTResultHitlistPtr) MemNew(sizeof(BLASTResultHitlist));
	if (new == NULL)
		return NULL;

	new->hsp_array = (BLASTResultHspPtr) MemNew(hspcnt*sizeof(BLASTResultHsp));
	if (new->hsp_array == NULL)
	{
		new = BLASTResultHitlistFree(new);
		return NULL;
	}
	new->hspcnt = hspcnt;

	return new; 
}


static Boolean 
CopyHSPToResultHsp(BLAST_KarlinBlkPtr kbp, BLAST_HSPPtr hsp, BLASTResultHspPtr result_hsp)
{
	if (result_hsp == NULL || hsp == NULL)
		return FALSE;

	result_hsp->ordering_method = hsp->ordering_method;
	result_hsp->number = hsp->num;
	result_hsp->score = hsp->score;
	result_hsp->bit_score = ((hsp->score*kbp->Lambda) - kbp->logK)/NCBIMATH_LN2;
	result_hsp->e_value = hsp->evalue;
	result_hsp->p_value = hsp->pvalue;
	result_hsp->query_offset = hsp->query.offset;
	result_hsp->query_length = hsp->query.length;
	result_hsp->query_frame = hsp->query.frame;
	result_hsp->query_gapped_start = hsp->query.gapped_start;
	result_hsp->subject_offset = hsp->subject.offset;
	result_hsp->subject_length = hsp->subject.length;
	result_hsp->subject_frame = hsp->subject.frame;
	result_hsp->subject_gapped_start = hsp->subject.gapped_start;
	result_hsp->context = hsp->context;
	result_hsp->gap_info = hsp->gap_info;

	/* Not set in the other type of HSP? */
	result_hsp->hspset_cnt = 0;

	return TRUE;
}

static Boolean 
CopyResultHspToHSP(BLASTResultHspPtr result_hsp, BLAST_HSPPtr hsp)
{
	if (result_hsp == NULL || hsp == NULL)
		return FALSE;

	hsp->ordering_method = result_hsp->ordering_method;
	hsp->num = result_hsp->number;
	hsp->score = result_hsp->score;
	hsp->evalue = result_hsp->e_value;
	hsp->pvalue = result_hsp->p_value;
	hsp->query.offset = result_hsp->query_offset;
	hsp->query.length = result_hsp->query_length;
	hsp->query.end = result_hsp->query_offset + result_hsp->query_length;
	hsp->query.frame = result_hsp->query_frame;
	hsp->query.gapped_start = result_hsp->query_gapped_start;
	hsp->subject.offset = result_hsp->subject_offset;
	hsp->subject.length = result_hsp->subject_length;
	hsp->subject.end = result_hsp->subject_offset + result_hsp->subject_length;
	hsp->subject.frame = result_hsp->subject_frame;
	hsp->subject.gapped_start = result_hsp->subject_gapped_start;
	hsp->context = result_hsp->context;

	return TRUE;
}


/*
	Sort the HSP's by score.
*/

static int LIBCALLBACK
score_compare_hsps(VoidPtr v1, VoidPtr v2)

{
	BLAST_HSPPtr h1, h2;
	BLAST_HSPPtr PNTR hp1, PNTR hp2;

	hp1 = (BLAST_HSPPtr PNTR) v1;
	hp2 = (BLAST_HSPPtr PNTR) v2;
	h1 = *hp1;
	h2 = *hp2;

	if (h1->score < h2->score) 
		return 1;
	if (h1->score > h2->score)
		return -1;

	return 0;
}

/*
	Function to look for the highest scoring window (of size HSP_MAX_WINDOW)
	in an HSP and return the middle of this.  Used by the gapped-alignment
	functions to start the gapped alignments.
*/

static Int4
GetStartForGappedAlignment (BlastSearchBlkPtr search, BLAST_HSPPtr hsp, Uint1Ptr query, Uint1Ptr subject, Int4Ptr PNTR matrix)
{
	Int4 index1, max_offset, score, max_score, hsp_end;
	Uint1Ptr query_var, subject_var;

	if (hsp->query.length <= HSP_MAX_WINDOW)
	{
		max_offset = hsp->query.offset + hsp->query.length/2;
	}
	else
	{
		hsp_end = hsp->query.offset + HSP_MAX_WINDOW;
		query_var = query + hsp->query.offset;
		subject_var = subject + hsp->subject.offset;
		score=0;
		for (index1=hsp->query.offset; index1<hsp_end; index1++)
		{
		  if (!(search->positionBased))
		    score += matrix[*query_var][*subject_var];
		  else
		    score += search->sbp->posMatrix[index1][*subject_var];
		  query_var++; subject_var++;
		}
		max_score = score;
		max_offset = hsp_end - 1;
		hsp_end = hsp->query.end;
		for (index1=hsp->query.offset + HSP_MAX_WINDOW; index1<hsp_end; index1++)
		{
		  if (!(search->positionBased)) {
		    score -= matrix[*(query_var-HSP_MAX_WINDOW)][*(subject_var-HSP_MAX_WINDOW)];
		    score += matrix[*query_var][*subject_var];
		  }
		  else {
		    score -= search->sbp->posMatrix[index1-HSP_MAX_WINDOW][*(subject_var-HSP_MAX_WINDOW)];
		    score += search->sbp->posMatrix[index1][*subject_var];
                  }
		if (score > max_score)
		  {
		    max_score = score;
		    max_offset = index1;
		  }
		query_var++; subject_var++;
	      }
	max_offset -= HSP_MAX_WINDOW/2;
	}
	return max_offset;
}

/*
	Gets the ratio used to change an evalue calculated with the subject
	sequence length to one with a db length.
*/

Nlm_FloatHi LIBCALL
GetDbSubjRatio(BlastSearchBlkPtr search, Int4 subject_length)
{
	Nlm_FloatHi db_subj_ratio;

	db_subj_ratio = (search->context_factor)*(search->dblen)/(subject_length);
       	if (StringCmp(search->prog_name, "tblastn") == 0 || StringCmp(search->prog_name, "tblastx") == 0)
	{
		db_subj_ratio *= 3;
	}
	
	return db_subj_ratio;
}

SeqAlignPtr LIBCALL 
BlastGetGapAlgnTbckWithReaddb (BlastSearchBlkPtr search, Int4 hit_number, Boolean ordinal_number)

{
	BLASTResultHitlistPtr   result_hitlist;
	BioseqPtr subject_bsp;
	Boolean subject_allocated = FALSE;
	Int4 index1, subject_length, rev_subject_length;
	SeqAlignPtr seqalign;
	SeqPortPtr spp;
	Uint1Ptr subject, rev_subject;

        result_hitlist = search->result_struct->results[hit_number];

	if (StringCmp(search->prog_name, "tblastn") == 0)
	{
		subject_bsp = readdb_get_bioseq(search->rdfp, result_hitlist->subject_id);
    		spp = SeqPortNew(subject_bsp, 0, -1, Seq_strand_plus, Seq_code_ncbi4na);
		/* make one longer to "protect" ALIGN. */
		subject = MemNew((1+subject_bsp->length)*sizeof(Uint1));
		for (index1=0; index1<subject_bsp->length; index1++)
		{
			subject[index1] = SeqPortGetResidue(spp);
		}
		/* Gap character in last space. */
		subject[subject_bsp->length] = NULLB;
		subject_length = subject_bsp->length;
		spp = SeqPortFree(spp);

    		spp = SeqPortNew(subject_bsp, 0, -1, Seq_strand_minus, Seq_code_ncbi4na);
		/* make one longer to "protect" ALIGN. */
		rev_subject = MemNew((1+subject_bsp->length)*sizeof(Uint1));
		for (index1=0; index1<subject_bsp->length; index1++)
		{
			rev_subject[index1] = SeqPortGetResidue(spp);
		}
		/* Gap character in last space. */
		rev_subject[subject_bsp->length] = NULLB;
		rev_subject_length = subject_bsp->length;
		spp = SeqPortFree(spp);
		subject_bsp = BioseqFree(subject_bsp);
		subject_allocated = TRUE;
	}
	else
	{
		subject_length = readdb_get_sequence(search->rdfp, result_hitlist->subject_id, (Uint1Ptr PNTR) &subject);
		rev_subject = NULL;
		rev_subject_length = 0;
	}

	seqalign = BlastGetGapAlgnTbck (search, hit_number,  FALSE, ordinal_number, subject, subject_length, rev_subject, rev_subject_length);
		
	if (subject_allocated)
	{
		subject = MemFree(subject);
		rev_subject = MemFree(rev_subject);
	}

	return seqalign;
}

/*
	Check the gapped alignments for an overlap of two different alignments.
	A sufficient overlap is when two alignments have the same start values
	of have the same final values.

	The number of valid alignments remaining is returned.
*/

static Int2
CheckGappedAlignmentsForOverlap(BlastSearchBlkPtr search, BLAST_HSPPtr *hsp_array, Int4 hsp_count, Int2 frame)

{
	BLAST_HSPPtr hsp, hsp1;
	Int4 index1, index;

	if (search == NULL || *hsp_array == NULL || hsp_count == 0)
		return 0;

	for (index=0; index<hsp_count; index++)
	{
		hsp = hsp_array[index];
		if (hsp == NULL)
			continue;
		if (frame != 0 && hsp->subject.frame != frame)
			continue;

		for (index1=0; index1<index; index1++)
		{
			hsp1 = hsp_array[index1];
			if (hsp1 == NULL)
				continue;
			if (frame != 0 && hsp1->subject.frame != frame)
				continue;

	/* Check if both HSP's start on or end on the same digonal. */
			if ((hsp1->query.offset == hsp->query.offset &&
			  hsp1->subject.offset == hsp->subject.offset) ||
			    (hsp1->query.end == hsp->query.end &&
			      hsp1->subject.end == hsp->subject.end))
			{
				if (hsp1->score > hsp->score)
				{
					hsp_array[index] = MemFree(hsp_array[index]);
					break; /* break out of inner loop. */
				}
				else
				{
					hsp_array[index1] = MemFree(hsp_array[index1]);
					/* look at rest of inner loop. */
				}
			}
		}
	}

	index1 = 0;
	for (index=0; index<hsp_count; index++)
	{
		if (hsp_array[index] != NULL)
			index1++;
	}

	return index1;

}

/*
	Engine to get the gapped scores from an array of HSP's.
*/
static BLAST_HSPPtr PNTR
BlastGappedScoreInternal(BlastSearchBlkPtr search, Uint1Ptr subject, Int4 subject_length, GapAlignBlkPtr gap_align, BLAST_HSPPtr *hsp_array, Int2Ptr hspcnt, Int2Ptr hspcnt_max, Int2 hspmax, Int2 frame)

{
	BLAST_HSPPtr hsp, hsp1;
	BLAST_HSPPtr PNTR hsp_array_new;
	Boolean hsp_start_is_contained, hsp_end_is_contained;
	Int2 hsp_cnt=0, index, index1;
	Int4 max_offset;

	for (index=0; index<(*hspcnt); index++)
	{
		hsp_start_is_contained = FALSE;
		hsp_end_is_contained = FALSE;
		hsp = hsp_array[index];
		if (frame != 0 && hsp->subject.frame != frame)
			continue;

		for (index1=0; index1<index; index1++)
		{
			hsp_start_is_contained = FALSE;
			hsp_end_is_contained = FALSE;

			hsp1 = hsp_array[index1];
			if (hsp1 == NULL)
				continue;
			if (frame != 0 && hsp1->subject.frame != frame)
				continue;

			if (CONTAINED_IN_HSP(hsp1->query.offset, hsp1->query.end, hsp->query.offset, hsp1->subject.offset, hsp1->subject.end, hsp->subject.offset) == TRUE)
			{
				hsp_start_is_contained = TRUE;
			}
			if (CONTAINED_IN_HSP(hsp1->query.offset, hsp1->query.end, hsp->query.end, hsp1->subject.offset, hsp1->subject.end, hsp->subject.end) == TRUE)
			{
				hsp_end_is_contained = TRUE;
			}
			if (hsp_start_is_contained && hsp_end_is_contained && hsp->score <= hsp1->score)
			{
				break;
			}
		}

		if (hsp_start_is_contained == FALSE ||
			 hsp_end_is_contained == FALSE || 
				hsp->score > hsp1->score)
		{
			gap_align->include_query = 0;
			max_offset = GetStartForGappedAlignment(search, hsp, search->context[hsp->context].query->sequence, subject, search->sbp->matrix);
#ifdef BLAST_COLLECT_STATS
			search->real_gap_number_of_hsps++;
#endif
			Nlm_MemSet((VoidPtr) &(hsp_array[index]->hsp_link), 0, sizeof(BLAST_HSP_LINK));
			hsp_array[index]->linked_set = FALSE;
			hsp_array[index]->start_of_chain = FALSE;
			hsp_array[index]->num = 0;
			hsp_array[index]->sumscore = 0;
			gap_align->query = search->context[hsp->context].query->sequence;
			gap_align->subject = subject;
			gap_align->query_length = search->context[hsp->context].query->length;
			gap_align->subject_length = subject_length;
			gap_align->q_start = max_offset;
			gap_align->s_start = (hsp->subject.offset - hsp->query.offset) + max_offset;
			hsp->query.gapped_start = gap_align->q_start;
			hsp->subject.gapped_start = gap_align->s_start;
			PerformGappedAlignment(gap_align);
			hsp->query.offset = gap_align->query_start;
			hsp->subject.offset = gap_align->subject_start;
	/* The end is one further for BLAST than for the gapped align. */
			hsp->query.end = gap_align->query_stop + 1;
			hsp->subject.end = gap_align->subject_stop + 1;
			hsp->query.length = hsp->query.end - hsp->query.offset;
			hsp->subject.length = hsp->subject.end - hsp->subject.offset;
			hsp->score = gap_align->score;
/* TLM */
			hsp->evalue = BlastKarlinStoE_simple(hsp->score, search->sbp->kbp_gap[search->first_context], search->searchsp_eff);
			hsp_cnt++;
		}
		else
		{ /* Contained within another HSP, delete. */
			hsp_array[index] = MemFree(hsp_array[index]);
		}
	}

	hsp_cnt = CheckGappedAlignmentsForOverlap(search, hsp_array, *hspcnt, frame);

	if (hsp_cnt < (*hspcnt))
	{
/* Save HSP's again, discarding those that have been NULLed out. */
		hsp_array_new = MemNew(hspmax*sizeof(BLAST_HSPPtr));
		index1 = 0;
		for (index=0; index<(*hspcnt_max); index++)
		{
			if (hsp_array[index] != NULL)
			{
				hsp_array_new[index1] = hsp_array[index];
				index1++;
			}
		}

		hsp_array = MemFree(hsp_array);

		*hspcnt = index1;	
		*hspcnt_max = index1;	
		hsp_array = hsp_array_new;
	}
	*hspcnt = hsp_cnt;

	return hsp_array;
}

/*
	Loads the HSP's into the BlastHitRangePtr.
*/
static Boolean
BlastHitRangeLoad (BlastSearchBlkPtr search, BLAST_HSPPtr *hsp_array, Int4 hspcnt, BlastHitRangePtr bhrp)

{
	BlastDoubleInt4Ptr tmp;
	BLAST_HSPPtr hsp;
	Int4 index, query_length;

	if (bhrp->current+hspcnt > bhrp->total)
		return FALSE;

	if (hspcnt <= 0)
		return TRUE;

	tmp = bhrp->range_list;

	tmp += bhrp->current;
	
	for (index=0; index<hspcnt; index++)
	{
		hsp = hsp_array[index];
		query_length = search->context[hsp->context].query->length;
		if (hsp->query.frame >= 0)
		{
			tmp->gi = hsp->query.offset;
			tmp->ordinal_id = hsp->query.end;
		}
		else
		{
			tmp->gi = query_length - hsp->query.offset - hsp->query.length;
			tmp->ordinal_id = query_length - hsp->query.offset;
		}
		tmp++;
	}

	bhrp->current += hspcnt;

	return TRUE;
}

/*
	Take a BLAST_HSPPtr (array of HSP's) and get a traceback for them.
*/

static Int2
RealBlastGetGappedAlignmentTraceback(BlastSearchBlkPtr search, Uint1Ptr subject, Int4 subject_length, Uint1Ptr rev_subject, Int4 rev_subject_length, SeqIdPtr subject_id, BLAST_HSPPtr *hsp_array, Int2 hspcnt, SeqAlignPtr *head, BlastHitRangePtr bhrp, Int4 min_score_to_keep, Boolean reverse, Int4 ordinal_id, Boolean do_traceback)

{
	BLAST_HSPPtr hsp, hsp1, hsp2;
	BLAST_HSPPtr *old_hsp_array;
	BLAST_ParameterBlkPtr pbp;
        BLASTResultHsp       	result_hsp;
	Boolean hsp_start_is_contained, hsp_end_is_contained, keep;
	Boolean current_hitlist_created, do_not_do;
	GapAlignBlkPtr gap_align;
	Int2 new_hspcnt=0;
	Int4 index, index1, index2, query_length, max_offset;
	Int4 old_hspmax, old_hspcnt, total;
	Int4Ptr translated_subject_length=NULL; 
	Int4Ptr translated_subject_length_orig=NULL;
	SeqAlignPtr seqalign, seqalign_var, *seqalign_array;
	StdSegPtr ssp;
	Uint1Ptr query, PNTR translated_subject=NULL, PNTR translated_subject_orig=NULL;
	ValNodePtr gi_list;

	pbp = search->pbp;
	MemSet(&result_hsp, 0, sizeof(BLASTResultHsp));

	seqalign=NULL;
	if (do_traceback)
		seqalign_array = MemNew(hspcnt*sizeof(SeqAlignPtr));

	if (search->gap_align == NULL)
	{
		search->gap_align = GapAlignBlkNew(1, 1);
	}
	gap_align = search->gap_align;

	gi_list = BlastGetAllowedGis(search, ordinal_id);

	gap_align->positionBased = search->positionBased;
	gap_align->gap_open = pbp->gap_open;
	gap_align->gap_extend = pbp->gap_extend;
	gap_align->x_parameter = pbp->gap_x_dropoff_final;
	gap_align->matrix = search->sbp->matrix;
	gap_align->posMatrix = search->sbp->posMatrix;
	for (index=0; index<hspcnt; index++)
	{
		hsp_start_is_contained = FALSE;
		hsp_end_is_contained = FALSE;
		hsp = hsp_array[index];

		for (index1=0; index1<index; index1++)
		{
			hsp_start_is_contained = FALSE;
			hsp_end_is_contained = FALSE;

			hsp1 = hsp_array[index1];
			if (hsp1 == NULL)
				continue;

			if (ContextToFrame(search, hsp->context) != ContextToFrame(search, hsp1->context))
				continue;

			if (CONTAINED_IN_HSP(hsp1->query.offset, hsp1->query.end, hsp->query.offset, hsp1->subject.offset, hsp1->subject.end, hsp->subject.offset) == TRUE)
			{
				hsp_start_is_contained = TRUE;
			}
			if (CONTAINED_IN_HSP(hsp1->query.offset, hsp1->query.end, hsp->query.end, hsp1->subject.offset, hsp1->subject.end, hsp->subject.end) == TRUE)
			{
				hsp_end_is_contained = TRUE;
			}
			if (hsp_start_is_contained && hsp_end_is_contained && hsp->score <= hsp1->score)
			{
				break;
			}
		}

		do_not_do = FALSE;
		/* Check whether this part of query has already been covered. */
		if (bhrp)
		{
			total = bhrp->current;
			for (index1=0; index1<total; index1++)
			{
				if (hsp->query.offset >= bhrp->range_list_pointer[index1]->gi &&
					hsp->query.end <= bhrp->range_list_pointer[index1]->ordinal_id)
				{
					do_not_do = TRUE;
					break;
				}
			}
		}

		if (do_not_do == FALSE && (hsp_start_is_contained == FALSE || hsp_end_is_contained == FALSE || 
					hsp->score > hsp1->score))
		{
			query = (Uint1Ptr) search->context[hsp->context].query->sequence;
			query_length = search->context[hsp->context].query->length;

			gap_align->include_query = 0;
			/* these should never both be zero. */
			if (hsp->query.gapped_start == 0 && hsp->subject.gapped_start == 0)
			{
				max_offset = GetStartForGappedAlignment(search, hsp, query, subject, search->sbp->matrix);
				gap_align->q_start = max_offset;
				gap_align->s_start = (hsp->subject.offset - hsp->query.offset) + max_offset;
				hsp->query.gapped_start = gap_align->q_start;
				hsp->subject.gapped_start = gap_align->s_start;
			}
			else
			{
				gap_align->q_start = hsp->query.gapped_start;
				gap_align->s_start = hsp->subject.gapped_start;
			}

/*
			gap_align->q_start_limit = hsp->query.offset;
			gap_align->s_start_limit = hsp->subject.offset;
			gap_align->q_end_limit = hsp->query.end;
			gap_align->s_end_limit = hsp->subject.end;
*/

			gap_align->query_frame = ContextToFrame(search, hsp->context);
			gap_align->query = query;

			gap_align->subject_frame = hsp->subject.frame;
			gap_align->subject = subject;

			gap_align->query_length = query_length;
			gap_align->subject_length = subject_length;

			gap_align->translate1 = FALSE;
			gap_align->translate2 = FALSE;
			if (StringCmp(search->prog_name, "blastx") == 0)
			{
				gap_align->translate1 = TRUE;
				gap_align->translate2 = FALSE;
			}
			if (StringCmp(search->prog_name, "tblastn") == 0)
			{
				gap_align->translate1 = FALSE;
				gap_align->translate2 = TRUE;
				if (translated_subject == NULL)
				{
					translated_subject_orig = MemNew(8*sizeof(Uint1Ptr));
					translated_subject = translated_subject_orig + 3;
					translated_subject_length_orig = MemNew(8*sizeof(Int4));
					translated_subject_length = translated_subject_length_orig + 3;
				}
				if (translated_subject[hsp->subject.frame] == NULL)
				{
				   if (hsp->subject.frame > 0)
					translated_subject[hsp->subject.frame] =
						GetTranslation(subject, subject_length, hsp->subject.frame, &translated_subject_length[hsp->subject.frame], search->db_genetic_code);
				   else
					translated_subject[hsp->subject.frame] =
						GetTranslation(rev_subject, rev_subject_length, hsp->subject.frame, &translated_subject_length[hsp->subject.frame], search->db_genetic_code);
				}
				gap_align->subject = translated_subject[hsp->subject.frame] + 1;
				gap_align->subject_length = translated_subject_length[hsp->subject.frame];
			}

			if (do_traceback)
				PerformGappedAlignmentWithTraceback(gap_align);
			else
				PerformGappedAlignment(gap_align);

			if (gap_align->score >= min_score_to_keep)
			{
				hsp->query.offset = gap_align->query_start;
				hsp->subject.offset = gap_align->subject_start;
			/* The end is one further for BLAST than for the gapped align. */
				hsp->query.end = gap_align->query_stop + 1;
				hsp->subject.end = gap_align->subject_stop + 1;
				hsp->query.length = hsp->query.end - hsp->query.offset;
				hsp->subject.length = hsp->subject.end - hsp->subject.offset;
				hsp->score = gap_align->score;
				if (do_traceback)
					hsp->gap_info = gap_align->edit_block;
				if (StringCmp(search->prog_name, "blastp") == 0 ||
					StringCmp(search->prog_name, "blastn") == 0)
				{
					hsp->evalue = BlastKarlinStoE_simple(hsp->score, search->sbp->kbp_gap[search->first_context], search->searchsp_eff);
					hsp->pvalue = BlastKarlinEtoP(hsp->evalue);
				}
				/* only one alignment considered for blast[np]. */
				/* This may be changed by LinkHsps for blastx or tblastn. */
				hsp->num = 1;
		
				keep = TRUE;
				for (index2=0; index2<index; index2++)
				{
					hsp2 = hsp_array[index2];
					if (hsp2 == NULL)
						continue;
 
        	/* Check if both HSP's start on or end on the same digonal. */
               			        if ((hsp->query.offset == hsp2->query.offset &&
               			          hsp->subject.offset == hsp2->subject.offset) ||
               			            (hsp->query.end == hsp2->query.end &&
               			              hsp->subject.end == hsp2->subject.end))
					{
						if (hsp2->score > hsp->score)
						{
							keep = FALSE;
							break;
						}
						else
						{
						    new_hspcnt--;
						    if (do_traceback)
						    {
							seqalign_array[index2] = SeqAlignFree(seqalign_array[index2]);
							hsp_array[index2]->gap_info = GapXEditBlockDelete(hsp_array[index2]->gap_info);
						    }
						    hsp_array[index2] = MemFree(hsp_array[index2]);
						}
					}
				}
				
				if (keep)
				{
					new_hspcnt++;
				}
				else
				{
					hsp_array[index] = MemFree(hsp_array[index]);
					if (do_traceback)
						gap_align->edit_block = GapXEditBlockDelete(gap_align->edit_block);
				}
			}
			else
			{	/* Should be kept? */
				hsp_array[index] = MemFree(hsp_array[index]);
				if (do_traceback)
					gap_align->edit_block = GapXEditBlockDelete(gap_align->edit_block);
			}
		}
		else
		{ /* Contained within another HSP, delete. */
			hsp_array[index] = MemFree(hsp_array[index]);
		}
	}

	/* Make up fake hitlist, relink and rereap. */
	if (StringCmp(search->prog_name, "blastx") == 0 || 
		StringCmp(search->prog_name, "tblastn") == 0)
	{
		new_hspcnt = HspArrayPurge(hsp_array, hspcnt);
		current_hitlist_created = FALSE;
		if (search->current_hitlist == NULL)
		{
			current_hitlist_created = TRUE;
			search->current_hitlist = BlastHitListNew(search);
		}
		old_hsp_array = search->current_hitlist->hsp_array;
		old_hspmax = search->current_hitlist->hspmax;
		old_hspcnt = search->current_hitlist->hspcnt;
		search->current_hitlist->hsp_array = hsp_array;
		search->current_hitlist->hspcnt = new_hspcnt;
		search->current_hitlist->hspmax = hspcnt;
		/* Use translated lenght for tblastn, real for others. */
		if (StringCmp(search->prog_name, "tblastn") == 0)
		{
			search->subject->length = translated_subject_length[1];
		}
		else
		{
			search->subject->length = subject_length;
		}
		BlastLinkHsps(search);
		BlastReapHitlistByEvalue(search);
		new_hspcnt = search->current_hitlist->hspcnt;
		search->current_hitlist->hsp_array = old_hsp_array;
		search->current_hitlist->hspcnt = old_hspcnt;
		search->current_hitlist->hspmax = old_hspmax;
		if (current_hitlist_created)
		{
			search->current_hitlist = BlastHitListDestruct(search->current_hitlist);
		}
	}

	new_hspcnt = HspArrayPurge(hsp_array, hspcnt);

	HeapSort(hsp_array,new_hspcnt,sizeof(BLAST_HSPPtr), score_compare_hsps);

	if (do_traceback)
	{
	    for (index=0; index<new_hspcnt; index++)
	    {
		hsp = hsp_array[index];
		hsp->gap_info->reverse = reverse;
		hsp->gap_info->original_length1 = search->context[search->first_context].query->original_length;
		hsp->gap_info->original_length2 = subject_length;
		if (do_traceback)
			seqalign = GapXEditBlockToSeqAlign(hsp->gap_info, subject_id, search->query_id); 
		CopyHSPToResultHsp(search->sbp->kbp_gap[search->first_context], hsp, &result_hsp);
		seqalign->score = GetScoreSetFromBlastResultHsp(&result_hsp, gi_list);
		if (seqalign->segtype == 3)
		{
			ssp = seqalign->segs;
			while (ssp)
			{
				ssp->scores = GetScoreSetFromBlastResultHsp(&result_hsp, gi_list);
				ssp = ssp->next;
			}
		}
		seqalign_array[index] = seqalign;
	   }
					
	   *head = NULL;
	   for (index=0; index<new_hspcnt; index++)
	   {
		if (seqalign_array[index] != NULL)
		{
			if (*head == NULL)
			{
				*head = seqalign_array[index];
			}
			else
			{
				for (seqalign_var=*head; seqalign_var->next != NULL;)
				{
					seqalign_var = seqalign_var->next;
				}
				seqalign_var->next = seqalign_array[index];
			}
		}
	   }

	   seqalign_array = MemFree(seqalign_array);
	}
	else
	{
		BlastHitRangeLoad(search, hsp_array, new_hspcnt, bhrp);
	}

	if (StringCmp(search->prog_name, "tblastn") == 0 && 
		translated_subject_orig)
	{
		for (index=0; index<8; index++)
		{
			MemFree(translated_subject_orig[index]);
		}
		MemFree(translated_subject_orig);
		MemFree(translated_subject_length_orig);
	}

	return new_hspcnt;
}


/*
	find the traceback for a gapped alignment.  Do this by 
	organizing the list of HSP's by sum group, then order
	these groups by score.  Then attempt to perform the alignment
	by using the highest scoring HSP of every sum group, then the
	2nd highest scoring HSP, etc. until all the HSP's of a sum
	group have been examined.  Then move onto the next sum group.
*/
SeqAlignPtr LIBCALL
SumBlastGetGappedAlignmentTraceback (BlastSearchBlkPtr search, Int4 hit_number, Boolean reverse, Boolean ordinal_number, Uint1Ptr subject, Int4 subject_length)

{
	SeqAlignPtr seqalign;

	SumBlastGetGappedAlignmentEx(search, hit_number, reverse, ordinal_number, subject, subject_length, TRUE, &seqalign, NULL);

	return seqalign;
}

static Boolean
SumBlastGetGappedAlignmentEx (BlastSearchBlkPtr search, Int4 hit_number, Boolean reverse, Boolean ordinal_number, Uint1Ptr subject, Int4 subject_length, Boolean do_traceback, SeqAlignPtr PNTR seqalignP, BlastHitRangePtr bhrp)

{
	BLAST_HSPPtr PNTR hsp_array;
	BLASTResultHitlistPtr   result_hitlist;
        BLASTResultHspPtr       result_hsp_array, hsp;
	Boolean not_done;
	Int2 hspcnt=0, new_hspcnt=0, hspset_cnt_old;
	Int4 index, index1, high_score=0, ordinal_id, next_start, start, stop;
	SeqAlignPtr seqalign=NULL;
	SeqIdPtr subject_id=NULL;
	Nlm_FloatHi current_evalue=DBL_MAX;

	if (search == NULL)
		return FALSE;


        result_hitlist = search->result_struct->results[hit_number];
        hspcnt = result_hitlist->hspcnt;

	subject_id = BlastGetSubjectId(search, hit_number, ordinal_number, NULL);
	ordinal_id = result_hitlist->subject_id;

	hsp_array = MemNew(hspcnt*sizeof(BLAST_HSPPtr));
	not_done = TRUE;
	start=0;
	next_start=0;
	while (not_done)
	{
		hsp = &(result_hitlist->hsp_array[start]);
		hspset_cnt_old = hsp->hspset_cnt;
		for (index=start; index<hspcnt; index++)
		{
			hsp = &(result_hitlist->hsp_array[index]);
			if(hspset_cnt_old != hsp->hspset_cnt)
			{
				hspset_cnt_old = hsp->hspset_cnt;
				stop = index;
				next_start = stop;
				break;
			}
		}

		if (index == hspcnt)
		{
			stop = hspcnt;
			not_done = FALSE;
		}

		index1=0;
		for (index=start; index<stop; index++)
		{
			hsp_array[index] = MemNew(sizeof(BLAST_HSP));
			CopyResultHspToHSP(&(result_hitlist->hsp_array[index]), hsp_array[index]);
			index1++;
		}

		/* heap sort the last sum group */
		HeapSort(hsp_array+start,(stop-start),sizeof(BLAST_HSPPtr), score_compare_hsps);
		start = next_start;
	}

	new_hspcnt = RealBlastGetGappedAlignmentTraceback(search, subject, subject_length, NULL, 0, subject_id, hsp_array, hspcnt, &seqalign, bhrp, search->pbp->cutoff_s, reverse, ordinal_id, do_traceback);

/* Save HSP's again, discarding those that have been NULLed out. */
/* If no HSP's were valid, best_evalue is set to DBL_MAX. */
	index1 = 0;
	if (new_hspcnt > 0)
	{
		result_hsp_array = MemNew((new_hspcnt)*sizeof(BLASTResultHsp));
		index1 = 0;
		for (index=0; index<hspcnt; index++)
		{
			if (hsp_array[index] != NULL)
			{
				if (current_evalue > hsp_array[index]->evalue)
					current_evalue = hsp_array[index]->evalue;
				if (high_score < hsp_array[index]->score)
       	                        	high_score = hsp_array[index]->score;
				CopyHSPToResultHsp(search->sbp->kbp_gap[search->first_context], hsp_array[index], &(result_hsp_array[index1]));
				index1++;
				hsp_array[index] = MemFree(hsp_array[index]);
			}
		}
	}
	hsp_array = MemFree(hsp_array);

	result_hitlist->hspcnt = index1;	
	if (result_hitlist->hsp_array)
		MemFree(result_hitlist->hsp_array);
	result_hitlist->hsp_array = result_hsp_array;
	result_hitlist->best_evalue = current_evalue;
	result_hitlist->high_score = high_score;

	subject_id = SeqIdSetFree(subject_id);

	if (seqalignP)
	*	seqalignP = seqalign;

	return TRUE;
}

/*
	Performs a gapped alignment on the HSP's in a hitlist.
	Discards those that do not meet the standard.
*/

SeqAlignPtr LIBCALL 
BlastGetGapAlgnTbck (BlastSearchBlkPtr search, Int4 hit_number, Boolean reverse, Boolean ordinal_number, Uint1Ptr subject, Int4 subject_length, Uint1Ptr rev_subject, Int4 rev_subject_length)

{
	BLAST_HSPPtr PNTR hsp_array;
	BLASTResultHitlistPtr   result_hitlist;
        BLASTResultHspPtr       result_hsp_array;
	Int2 hspcnt=0, new_hspcnt=0;
	Int4 index, index1, high_score=0, ordinal_id;
	SeqAlignPtr seqalign, head, seqalign_var;
	SeqIdPtr subject_id=NULL;
	Nlm_FloatHi current_evalue=DBL_MAX;

	if (search == NULL)
		return NULL;


        result_hitlist = search->result_struct->results[hit_number];
        hspcnt = result_hitlist->hspcnt;
	ordinal_id = result_hitlist->subject_id;

	subject_id = BlastGetSubjectId(search, hit_number, ordinal_number, NULL);

	head = NULL;

	hsp_array = MemNew(hspcnt*sizeof(BLAST_HSPPtr));
	for (index=0; index<hspcnt; index++)
	{
		hsp_array[index] = MemNew(sizeof(BLAST_HSP));
		CopyResultHspToHSP(&(result_hitlist->hsp_array[index]), hsp_array[index]);
	}
	HeapSort(hsp_array,hspcnt,sizeof(BLAST_HSPPtr), score_compare_hsps);

	new_hspcnt = RealBlastGetGappedAlignmentTraceback(search, subject, subject_length, rev_subject, rev_subject_length, subject_id, hsp_array, hspcnt, &seqalign, NULL, 0, reverse, ordinal_id, TRUE);
	if (seqalign != NULL)
	{
		if (head == NULL)
		{
			head = seqalign;
		}
		else
		{
			for (seqalign_var=head; seqalign_var->next != NULL;)
			{
				seqalign_var = seqalign_var->next;
			}
			seqalign_var->next = seqalign;
		}
	}

/* Save HSP's again, discarding those that have been NULLed out. */
	result_hsp_array = MemNew((new_hspcnt)*sizeof(BLASTResultHsp));
	index1 = 0;
	for (index=0; index<hspcnt; index++)
	{
		if (hsp_array[index] != NULL)
		{
			if (current_evalue > hsp_array[index]->evalue)
				current_evalue = hsp_array[index]->evalue;
			if (high_score < hsp_array[index]->score)
                               	high_score = hsp_array[index]->score;
			CopyHSPToResultHsp(search->sbp->kbp_gap[search->first_context], hsp_array[index], &(result_hsp_array[index1]));
			index1++;
			hsp_array[index] = MemFree(hsp_array[index]);
		}
	}
	hsp_array = MemFree(hsp_array);

	result_hitlist->hspcnt = index1;	
	if (result_hitlist->hsp_array)
		MemFree(result_hitlist->hsp_array);
	result_hitlist->hsp_array = result_hsp_array;
	result_hitlist->best_evalue = current_evalue;
	result_hitlist->high_score = high_score;

	subject_id = SeqIdSetFree(subject_id);

	return head;
}

/*
	Performs a gapped alignment on the HSP's in a hitlist.
	Discards those that do not meet the standard.
*/

Int2 LIBCALL
BlastPreliminaryGappedScore (BlastSearchBlkPtr search, Uint1Ptr subject, Int4 subject_length, Int2 frame)

{
	BLAST_HitListPtr hitlist;
	BLAST_HSPPtr hsp;
	BLAST_HSPPtr PNTR hsp_array;
	GapAlignBlkPtr gap_align;
	Int2 hspcnt_max, status;
	Int4 index, max_offset, query_length;
	Nlm_FloatHi e_value;
	BLAST_ParameterBlkPtr pbp;

	if (search == NULL)
		return 1;

	pbp = search->pbp;

	if (search->gap_align == NULL)
	{
		search->gap_align = GapAlignBlkNew(1, 1);
	}
	gap_align = search->gap_align;

	status = 0;
	hitlist = search->current_hitlist;
	if (hitlist && hitlist->hspcnt > 0)
	{
		hspcnt_max = hitlist->hspcnt;
		query_length = search->context[search->first_context].query->length;

		hitlist->hspcnt_max = hitlist->hspcnt;
		hsp_array = hitlist->hsp_array;
		if (frame != 0)
		{
			for (index=0; index<hitlist->hspcnt; index++)
			{
				hsp = hsp_array[index];
				if (frame == hsp->subject.frame)
					break;
			}
			if (frame != hsp->subject.frame)
				return 0;
		}
		else
		{ /* The first HSP has the highest score. */
			hsp = hsp_array[0];
		}

		/* The first HSP has the highest score. */
		e_value = BlastKarlinStoE_simple(hsp->score, search->sbp->kbp_gap[search->first_context], search->searchsp_eff);
		if (e_value <= pbp->cutoff_e)
		{
#ifdef BLAST_COLLECT_STATS
			search->prelim_gap_no_contest++;
#endif
			hitlist->further_process = TRUE;
			return 1;
		}

		gap_align->positionBased = search->positionBased;
		gap_align->include_query = 0;
		gap_align->gap_open = pbp->gap_open;
		gap_align->gap_extend = pbp->gap_extend;
		gap_align->x_parameter = pbp->gap_x_dropoff;
		gap_align->matrix = search->sbp->matrix;
		gap_align->posMatrix = search->sbp->posMatrix;
		for (index=0; index<hitlist->hspcnt; index++)
		{
			hsp = hsp_array[index];
			if (frame != 0)
			{
				if (frame != hsp->subject.frame)
					continue;
			}

			if (hsp->score < search->pbp->gap_trigger)
			{	/* Stop looking, we're below the cutoff. */
				status = 0;
				break;
			}

#ifdef BLAST_COLLECT_STATS
				search->prelim_gap_attempts++;
#endif
			gap_align->score = 0;
			gap_align->query = search->context[hsp->context].query->sequence;
			gap_align->subject = subject;
			gap_align->query_length = search->context[hsp->context].query->length;
			gap_align->subject_length = subject_length;
			gap_align->include_query = 0;
			max_offset = GetStartForGappedAlignment(search, hsp, search->context[hsp->context].query->sequence, subject, search->sbp->matrix);
			gap_align->q_start = max_offset;
			gap_align->s_start = (hsp->subject.offset - hsp->query.offset) + max_offset;
	/* Perform only if the query's required start corresponds to a point after the start of the subject. */
			if (gap_align->s_start >= 0)
				PerformGappedAlignment(gap_align);
			e_value = BlastKarlinStoE_simple(gap_align->score, search->sbp->kbp_gap[search->first_context], search->searchsp_eff);
			if (e_value <= pbp->cutoff_e)
			{	/* Found one, stop looking. */
				hitlist->further_process = TRUE;
				status = 1;
#ifdef BLAST_COLLECT_STATS
				search->prelim_gap_passed++;
#endif
				break;
			}
		}
	}

	return status;
}

/*
	Performs a gapped alignment on the HSP's in a hitlist.
	Discards those that do not meet the standard.
	Do this by obtaining the sequence from readdb and calling
	BlastGetGappedScore.
*/

Int2 LIBCALL
BlastGetGappedScoreWithReaddb (BlastSearchBlkPtr search, Int4 sequence_number)

{
	BLAST_HitListPtr hitlist;
	Int2 retval;
	Int4 subject_length;
	Uint1Ptr subject;

	if (search == NULL)
		return 1;

	retval=0;
	hitlist = search->current_hitlist;
	if (hitlist && hitlist->hspcnt > 0)
	{
		if (hitlist->further_process == FALSE)
		{
			BlastHitListPurge(hitlist);
			return 0;
		}
		subject_length = readdb_get_sequence(search->rdfp, sequence_number, &subject);
		retval = BlastGetGappedScore(search, subject_length, subject, 0);
	}

	return retval;
}


/*
	Performs a gapped alignment on the HSP's in a hitlist.
	Discards those that do not meet the standard.
*/

Int2 LIBCALL
BlastGetGappedScore (BlastSearchBlkPtr search, Int4 subject_length, Uint1Ptr subject, Int2 frame)

{
	BLAST_HitListPtr hitlist;
	BLAST_HSPPtr PNTR hsp_array, PNTR hsp_array_new;
	BLAST_ParameterBlkPtr pbp;
	GapAlignBlkPtr gap_align;
	Int2 hsp_cnt=0, hspcnt_max, status=0;
	Int4 index, index1;

	if (search == NULL)
		return 1;

	pbp = search->pbp;


	if (search->gap_align == NULL)
	{
		search->gap_align = GapAlignBlkNew(1, 1);
	}
	gap_align = search->gap_align;

	hitlist = search->current_hitlist;
	if (hitlist && hitlist->hspcnt > 0)
	{
		if (hitlist->further_process == FALSE)
		{
			BlastHitListPurge(hitlist);
			return 0;
		}

		
		hsp_array = hitlist->hsp_array;
		if (hitlist->hspcnt != hitlist->hspcnt_max)
		{
/* Save HSP's again, discarding those that have been NULLed out. */
			hsp_array_new = MemNew((hitlist->hspmax)*sizeof(BLAST_HSPPtr));
			index1 = 0;
			for (index=0; index<hitlist->hspcnt_max; index++)
			{
				if (hsp_array[index] != NULL)
				{
					hsp_array_new[index1] = hsp_array[index];
					index1++;
				}
			}
			hsp_array = MemFree(hsp_array);
			hsp_array = hsp_array_new;
			hitlist->hsp_array = hsp_array_new;
			hitlist->hspcnt = index1;
			hitlist->hspcnt_max = index1;
		}

		gap_align->positionBased = search->positionBased;
		gap_align->include_query = 0;
		gap_align->gap_open = pbp->gap_open;
		gap_align->gap_extend = pbp->gap_extend;
		gap_align->x_parameter = pbp->gap_x_dropoff;
		gap_align->matrix = search->sbp->matrix;
		gap_align->posMatrix = search->sbp->posMatrix;

		HeapSort(hsp_array,hitlist->hspcnt,sizeof(BLAST_HSPPtr), score_compare_hsps);
		hitlist->hspcnt_max = hitlist->hspcnt;
		hsp_array = hitlist->hsp_array;
		hspcnt_max = hitlist->hspcnt;
		hsp_cnt = hitlist->hspcnt;

		hsp_array = BlastGappedScoreInternal(search, subject, subject_length, gap_align, hsp_array, &hsp_cnt, &hspcnt_max, hitlist->hspmax, frame);
		hitlist->hspcnt = hsp_cnt;	
		hitlist->hspcnt_max = hspcnt_max;	
		hitlist->hsp_array = hsp_array;
	}

	return status;
}


/******************************************************************

	Purges (i.e., cleans) the HitList for reuse.

*******************************************************************/

Int2 LIBCALL
BlastHitListPurge(BLAST_HitListPtr hitlist)

{
	BLAST_HSPPtr PNTR hsp_array;
	Int4 hspcnt_max, index;

	if (hitlist == NULL)
		return 1;

	hsp_array = hitlist->hsp_array;

	if (hitlist->hspcnt > hitlist->hspcnt_max)
		hspcnt_max = hitlist->hspcnt;
	else
		hspcnt_max = hitlist->hspcnt_max;

	for (index=0; index<hspcnt_max; index++)
	{
		hsp_array[index] = MemFree(hsp_array[index]);
	}

	hitlist->hspcnt = 0;
	hitlist->hspcnt_max = 0;
	hitlist->further_process = FALSE;

	return 0;
}

/*
	Cleans out the NULLed out HSP's from the HSP array,
	moving the BLAST_HSPPtr's up to fill in the gaps.

	returns the number of valid HSP's.
*/

Int2 LIBCALL
HspArrayPurge (BLAST_HSPPtr PNTR hsp_array, Int2 hspcnt)

{
	Int2 index, index1;

	if (hspcnt == 0 || hsp_array == NULL)
		return 0;

	index1 = 0;
	for (index=0; index<hspcnt; index++)
	{
		if (hsp_array[index] != NULL)
		{
			hsp_array[index1] = hsp_array[index];
			hsp_array[index1]->num = 0;
			index1++;
		}
	}

	for (index=index1; index<hspcnt; index++)
	{
		hsp_array[index] = NULL;
	}

	hspcnt = index1;

	return index1;
}


/*
	Compares two HSP's looking for overlap in both the
	query and the subject sequence.  

	if hsp1 should be deleted, then -1 is returned,
	if hsp2 should be deleted, then +1 is returned,
	if neither should be deleted, then 0 is returned.
*/

static Int2 
CheckHspOverlap (BLAST_HSPPtr PNTR hsp_array, BLAST_HSPPtr hsp2, Int4 hspcnt, Boolean PNTR hsp_deleted)

{
	BLAST_HSPPtr hsp1;
	BLAST_SegPtr query, subject;
	Int4 index;

	if (hsp_array == NULL || hsp2 == NULL)
		return 0;


	*hsp_deleted = FALSE;
	query = &(hsp2->query);
	subject = &(hsp2->subject);

	for (index=0; index<hspcnt; index++)
	{
		hsp1 = hsp_array[index];
		if (hsp1 == NULL)
			continue;

		if (SIGN(hsp1->query.frame) != SIGN(query->frame))
			continue;

		if (SIGN(hsp1->subject.frame) != SIGN(subject->frame))
			continue;

		if (hsp1->query.offset > query->offset && 
			hsp1->query.end > query->end) 
			continue;

		if (hsp1->query.offset < query->offset && 
			hsp1->query.end < query->end) 
			continue;

		if (hsp1->subject.offset > subject->offset && 
			hsp1->subject.end > subject->end) 
			continue;

		if (hsp1->subject.offset < subject->offset && 
			hsp1->subject.end < subject->end) 
			continue;

		if (hsp1->score > hsp2->score)
		{
			if (*hsp_deleted == FALSE)
			{
				return 1;
			}
		}
		else
		{
			hsp_array[index] = MemFree(hsp_array[index]);
			*hsp_deleted = TRUE;
		}
	}

	return 0;
}
	
/**************************************************************************
*
*	Save the current HSP in the appropriate ranking.
*
**************************************************************************/

void 
BlastSaveCurrentHsp(BlastSearchBlkPtr search, BLAST_Score score, Int4 q_offset, Int4 s_offset, Int4 length, Int2 context)

{
	BLAST_HitListPtr current_hitlist;
	BLAST_HSPPtr PNTR hsp_array, new_hsp;
	BLAST_Score highscore, lowscore;
	Boolean hsp_deleted;
	Int2 status;
	Int4 hspcnt, hspmax, index, new_index, high_index, old_index, low_index;

	current_hitlist = search->current_hitlist;
	hsp_array = current_hitlist->hsp_array;
	hspcnt = current_hitlist->hspcnt;
	hspmax = current_hitlist->hspmax;


	if (hspcnt != 0)
	{
		highscore = hsp_array[0]->score;
		lowscore = hsp_array[hspcnt-1]->score;
	}
	else
	{
		highscore = 0;
		lowscore = 0;
	}

	/* Check if list is already full and this is a lower score. */
	if (score <= lowscore && hspcnt >= hspmax-1)
		return;

	new_hsp = (BLAST_HSPPtr) MemNew(sizeof(BLAST_HSP));
	new_hsp->score = score;
	new_hsp->query.offset = q_offset;
	new_hsp->subject.offset = s_offset;
	new_hsp->query.length = length;
	new_hsp->subject.length = length;
	new_hsp->query.end = q_offset + length;
	new_hsp->subject.end = s_offset + length;
	new_hsp->context = context;
	new_hsp->query.frame = ContextToFrame(search, context);
	new_hsp->subject.frame = search->subject->frame;

	hsp_deleted = FALSE;
	status = CheckHspOverlap(hsp_array, new_hsp, hspcnt, &hsp_deleted);
	if (status == 1 && hsp_deleted == FALSE)
	{ /* keep if this HSP covers another one. */
		new_hsp = MemFree(new_hsp);
		return;
	}

	if (hsp_deleted)
	{
		hspcnt = HspArrayPurge(hsp_array, hspcnt);	
		current_hitlist->hspcnt = hspcnt;
	}

	if (score >= highscore)
	{
		new_index = 0;
	}
	else if (score <= lowscore)
	{
		new_index = hspcnt;
	}
	else
	{
		low_index = 0;
		high_index = hspcnt-1;
		new_index = (low_index+high_index)/2;
		old_index = new_index;

		for (index=0; index<BLAST_SAVE_ITER_MAX; index++)
		{
			if (score > hsp_array[new_index]->score)
			{
				high_index = new_index;
			}
			else
			{
				low_index = new_index;
			}
			new_index = (low_index+high_index)/2;
                        if (new_index == old_index)
                        { /* Perform this check as new_index get rounded DOWN a
bove.*/
                                if (score < hsp_array[new_index]->score)
                                {
                                        new_index++;
                                }
                                break;
                        }
                        old_index = new_index;
		}
	}

	if (hspcnt >= hspmax-1)
	{
		if (new_index >= hspcnt)
		{ /* this HSP is less significant than others on a full list.*/
			return;
		}
		else
		{ /* Delete the last HPS on the list. */
			hspcnt = current_hitlist->hspcnt--;
			hsp_array[hspcnt-1] = MemFree(hsp_array[hspcnt-1]);
		}
	}
	current_hitlist->hspcnt++;
	Nlm_MemMove((hsp_array+new_index+1), (hsp_array+new_index), (hspcnt-new_index)*sizeof(hsp_array[0]));
	hsp_array[new_index] = new_hsp;

	return;
}

Uint1Ptr
GetSequenceWithDenseSeg(DenseSegPtr dsp, Boolean query, Int4Ptr start, Int4Ptr length)

{
	BioseqPtr bsp;
	Int4 index, offset;
	SeqIdPtr id;
	SeqPortPtr spp;
	Uint1Ptr buffer;

	if (dsp == NULL)
		return NULL;

	if (query == TRUE)
	{
		offset = 0;
		id = dsp->ids;
	}
	else
	{
		offset = 1;
		id = dsp->ids->next;
	}

	*start = dsp->starts[offset];
	*length = 0;
	for (index=0; index<dsp->numseg; index++)
	{
		if (dsp->starts[offset+2*index] != -1)
			*length += dsp->lens[index];
	}

	bsp = BioseqLockById(id);

	spp = SeqPortNew(bsp, *start, (*start)+(*length)-1, Seq_strand_unknown, Seq_code_ncbistdaa);

	buffer = MemNew((*length)*sizeof(Uint1));

	for (index=0; index<*length; index++)
		buffer[index] = SeqPortGetResidue(spp);

	spp = SeqPortFree(spp);
	BioseqUnlock(bsp);

	return buffer;
}

/* 
Produces a 'fake' BioseqPtr, for use with BLAST when the
ID of the original BioseqPtr cannot be trusted.  Note that
the ID of the original BioseqPtr is removed. 
*/

BioseqPtr LIBCALL
BlastMakeFakeBioseq(BioseqPtr bsp, CharPtr name)

{
	BioseqPtr fake_bsp;
	ObjectIdPtr obidp;

	if (bsp == NULL)
		return NULL;

        fake_bsp = BioseqNew();
        fake_bsp->descr = bsp->descr;
        fake_bsp->repr = bsp->repr;
        fake_bsp->mol = bsp->mol;
        fake_bsp->length = bsp->length;
        fake_bsp->strand = bsp->strand;
        fake_bsp->seq_data_type = bsp->seq_data_type;
        fake_bsp->seq_ext_type = bsp->seq_ext_type;
        fake_bsp->seq_data = bsp->seq_data;
        fake_bsp->seq_ext = bsp->seq_ext;

        obidp = ObjectIdNew();
	if (name)
   	     obidp->str = StringSave(name);
	else
   	     obidp->str = StringSave("QUERY");
        ValNodeAddPointer(&(fake_bsp->id), SEQID_LOCAL, obidp);

	return fake_bsp;
}

BioseqPtr LIBCALL
BlastDeleteFakeBioseq(BioseqPtr fake_bsp)

{
	if (fake_bsp == NULL)
		return NULL;

         fake_bsp->descr = NULL;
         fake_bsp->length = 0;
         fake_bsp->seq_data = NULL;
         fake_bsp->seq_ext = NULL;

	return BioseqFree(fake_bsp);
}

/*
	Makes a new BlastSeqIdListPtr, with the mode set ot BLAST_NOT_OWN.
*/
BlastSeqIdListPtr BlastSeqIdListNew (void)

{
	BlastSeqIdListPtr blast_seqid_list;

	blast_seqid_list = (BlastSeqIdListPtr) MemNew(sizeof(BlastSeqIdList));
	blast_seqid_list->mode = BLAST_NOT_OWN;

	return blast_seqid_list;
}

/*
	Deletes the BlastSeqIdListPtr, deletes the SeqIdPtr's only
	if they are non-NULL.
*/
BlastSeqIdListPtr BlastSeqIdListDestruct (BlastSeqIdListPtr blast_seqid_list)

{
	if (blast_seqid_list == NULL)
		return NULL;

	if (blast_seqid_list->mode == BLAST_OWN)
	{
		blast_seqid_list->seqid_list = SeqIdSetFree(blast_seqid_list->seqid_list);
	}

	return MemFree(blast_seqid_list);
}

/*
	This function sets the 'seqid_list' member of the BlastSearchBlkPtr,
	which limits the BLAST search to the specified members of the
	database.

	If the ordinal_id should be used, then it should be greater 
	than or equal to zero; the SeqIdPtr should be NULL.   If the 
	SeqIdPtr will be used, then it should non-NULL, ordinal_id
	should be '-1'.

	If a SeqIdPtr is used, then it will not be deleted by BLAST, it is up to 
	the user to delete it.  The SeqIdPtr created for the ordinal ID will
	be deleted.

	This function may be called multiple times.  It should only be called with
	either ordianl_id's or seqIdPtr's.
*/

Boolean BlastAddSeqIdToList(BlastSearchBlkPtr search, Int4 ordinal_id, SeqIdPtr new_sip)

{
	BlastSeqIdListPtr	blast_seqid_list;
        DbtagPtr        	dbtagptr;


	if (search == NULL)
		return FALSE;

	blast_seqid_list = search->blast_seqid_list;
	if (blast_seqid_list == NULL)
	{
		blast_seqid_list = search->blast_seqid_list = BlastSeqIdListNew();
	}

	if (ordinal_id >= 0)
	{
		if (blast_seqid_list->mode == BLAST_NOT_OWN)
			return FALSE;
		dbtagptr = DbtagNew();
		dbtagptr->db = StringSave("BL_ORD_ID");
		dbtagptr->tag = ObjectIdNew();
		dbtagptr->tag->id = ordinal_id;
		ValNodeAddPointer(&blast_seqid_list->seqid_list, SEQID_GENERAL, dbtagptr);
		blast_seqid_list->mode = BLAST_OWN;
		return TRUE;
	}
	else if (new_sip)
	{
		ValNodeLink(&(blast_seqid_list->seqid_list), new_sip);
		blast_seqid_list->mode = BLAST_NOT_OWN;
		return TRUE;
	}
	return FALSE;
}

/*
	Counts the number of SeqAligns present.
*/

static Int4 
GetSeqAlignCount(SeqAlignPtr sap)

{
	Int4 count = 1;
	SeqIdPtr last_id=NULL, id;

	while (sap)
	{
		id = TxGetSubjectIdFromSeqAlign(sap);
		if (last_id)
		{
			if(SeqIdComp(id, last_id) != SIC_YES)
				count++;
		}
		else
		{
			count = 1;
		}
		last_id = id;
		sap = sap->next;
	}

	return count;

}

/*
	Duplicates a SeqAlignPtr, up to the number of unique
	records specified.
*/

static SeqAlignPtr 
GetPrivateSeqAlign(SeqAlignPtr sap, Int4 number, Int4Ptr number_returned)

{
	Int2 count=0;
	SeqIdPtr last_id=NULL, id;
	SeqAlignPtr new_head=NULL, var;

	last_id = TxGetSubjectIdFromSeqAlign(sap);

	while (count<number && sap)
	{
		count++;
		while (sap)
		{
			id = TxGetSubjectIdFromSeqAlign(sap);
			if(SeqIdComp(id, last_id) != SIC_YES)
			{
				last_id = id;
				break;
			}
			if (new_head == NULL)
			{
				new_head = AsnIoMemCopy(sap, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);
				var = new_head;
			}
			else
			{
				var->next = AsnIoMemCopy(sap, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);
				var = var->next;
			}
			last_id = id;
			sap = sap->next;
		}
	}

	*number_returned = count;

	return new_head;
}


/*
	Duplicate a SeqAlignPtr, keeping on the number of unique db
	hits specified.
*/

BlastPruneSapStructPtr
BlastPruneHitsFromSeqAlign(SeqAlignPtr sap, Int4 number, BlastPruneSapStructPtr prune)

{
	if (prune == NULL)
	{
		prune = MemNew(sizeof(BlastPruneSapStruct));
	}
	else
	{
		if (prune->number == number)
			return prune;
		if (prune->allocated)
			prune->sap = SeqAlignSetFree(prune->sap);
		prune->sap = NULL;
		prune->allocated = FALSE;
		prune->original_number = 0;
		prune->number = 0;
	}

	prune->original_number = GetSeqAlignCount(sap);

	if (prune->original_number < number)
	{
		prune->number = prune->original_number;
		prune->sap = sap;
		prune->allocated = FALSE;
	}
	else
	{
		prune->sap = GetPrivateSeqAlign(sap, number, &(prune->number));
		prune->allocated = TRUE;
	}
	
	return prune;
}


BlastPruneSapStructPtr
BlastPruneSapStructDestruct(BlastPruneSapStructPtr prune)

{
	if (prune == NULL)
		return NULL;

	if (prune->allocated)
	{
		prune->sap = SeqAlignSetFree(prune->sap);
	}
	prune = MemFree(prune);

	return prune;
}

/*
	Returns the program number for a string containing the
	program name.
*/

Uint1 LIBCALL
BlastGetProgramNumber(CharPtr blast_program)

{
	if (blast_program == NULL)
		return blast_type_undefined;

	if (StringICmp("blastn", blast_program) == 0)
	{
		return blast_type_blastn;
	}
	else if (StringICmp("blastp", blast_program) == 0)
	{
		return blast_type_blastp;
	}
	else if (StringICmp("blastx", blast_program) == 0)
	{
		return blast_type_blastx;
	}
	else if (StringICmp("tblastn", blast_program) == 0)
	{
		return blast_type_tblastn;
	}
	else if (StringICmp("tblastx", blast_program) == 0)
	{
		return blast_type_tblastx;
	}

	return blast_type_undefined;
}

/*
	returns information aobut the db and query types (protein or dna)
	as well as the 'align_type' that should be attached to the SeqAnnot
	for formatting.

	If an invalid program is entered, then 0 is returned.
*/

Uint1 LIBCALL
BlastGetTypes(CharPtr blast_program, Boolean PNTR query_is_na, Boolean PNTR db_is_na)

{
	Uint1 align_type=0;

	align_type = BlastGetProgramNumber(blast_program);

	if (align_type == blast_type_blastn)
	{
		*query_is_na = TRUE;
		*db_is_na = TRUE;
	}
	else if (align_type == blast_type_blastp)
	{
		*query_is_na = FALSE;
		*db_is_na = FALSE;
	}
	else if (align_type == blast_type_blastx)
	{
		*query_is_na = TRUE;
		*db_is_na = FALSE;
	}
	else if (align_type == blast_type_tblastn)
	{
		*query_is_na = FALSE;
		*db_is_na = TRUE;
	}
	else if (align_type == blast_type_tblastx)
	{
		*query_is_na = TRUE;
		*db_is_na = TRUE;
	}

	return align_type;
}
 

/*
	Find the word hits for a nucl. query.  No neighbors are found here.
	If no indices are saved, then return 1, indicating that the
	search should not be performed.

*/

Int2
BlastNtFindWords(BlastSearchBlkPtr search, Int4 start, Int4 len, Int1 context_index) 
{
	register Int4 offset, initial_wordsize;
	Boolean found_ambig, saved_index=FALSE;
	BLAST_ScorePtr PNTR	matrix;
	BLAST_WordFinderPtr	wfp;
	Int4 end, index, index_addition, lookup_index, stop;
	LookupTablePtr		lookup;
	Uint1Ptr str;
	ValNodePtr              vnp, vnp_start=NULL;

	matrix = search->sbp->matrix;

	if (search == NULL)
	{
		return -1;
	}

	wfp = search->wfp;
	if (wfp == NULL)
	{
		return -2;
	}

	lookup = wfp->lookup;
	if (lookup == NULL)
	{
		return -3;
	}

	initial_wordsize = (lookup->wordsize)*(wfp->compression_ratio);

	vnp = search->context[context_index].location;
	if (vnp == NULL)
	{
		ValNodeAddInt(&vnp, 1, -1);
		vnp_start = vnp;
		ValNodeAddInt(&vnp, 0, len);
	}

        while (vnp)
	{
	    if (vnp->choice == 1)
	    {
	    	start = vnp->data.intvalue + 1;
	    	vnp = vnp->next;
	    	if (vnp == NULL)
	    		end = len;
	    }
	    if (vnp && vnp->choice == 0)
	    {
	    	end = vnp->data.intvalue - initial_wordsize;
	    	vnp = vnp->next;
	    }

	    end = MIN(end, len-initial_wordsize);

	    str = (Uint1Ptr) search->context[context_index].query->sequence + start;
		
	    for (offset=start; offset<end; offset++, str++)
	    {
		found_ambig= FALSE;
		lookup_index = 0;
		stop = initial_wordsize/wfp->compression_ratio;
		index_addition = 0;
		for (index=0; index<stop; index++)
		{
			if (*(str+index_addition) > 3 || *(str+index_addition+1) > 3 || *(str+index_addition+2) > 3 || *(str+index_addition+3) > 3)
			{
				found_ambig = TRUE;
				break;
			}

			lookup_index += (*(str+index_addition)   << 6);
			lookup_index += (*(str+1+index_addition) << 4);
			lookup_index += (*(str+2+index_addition) << 2);
			lookup_index += *(str+3+index_addition);

			if (index != stop-1)
			{	/* 8 bits/byte */
					lookup_index <<= 8;
			  		index_addition += 4;
			}
		}
			
		if (found_ambig == FALSE)
		{
			lookup_add_index(lookup, (Int4) lookup_index, offset+initial_wordsize, context_index);
			saved_index = TRUE;
		}
	    }
	}

	if (vnp_start)
	{
		vnp_start = ValNodeFree(vnp_start);
	}

	if (saved_index == FALSE)
		return 1;

	return 0;
}

/*
	This functions finds the words.
	return values:
		0: success, words saved
		1: no words saved, no error
		-1: error

*/
		    
	
Int2 LIBCALL
BlastFindWords(BlastSearchBlkPtr search, Int4 start, Int4 len, BLAST_Score threshold, Int1 context_index) 

{
	register Uint1		last_char, last_char2;
	Uint1Ptr		words, PNTR array;
	Uint1Ptr		s_string_start, s_string;
	register Uint1Ptr	str;
	BLAST_Score		best_total, delta_score, diff, diff2, first_score;
	BLAST_Score		second_score, start_score, start_score2, score;
	BLAST_ScoreBlkPtr 	sbp;
	register BLAST_ScorePtr PNTR	matrix;
	BLAST_WordFinderPtr	wfp;
	Boolean			exact_match, saved_index=FALSE;
	LookupTablePtr		lookup;
	register Int4		index1, index3, offset;
	register Int1		index2;
	Int4			num_of_cols, alphabet_size, wordsize;
	Int4 			loop_increment, loop_increment2, stop;
	SeqCodeTablePtr 	sctp;
	ValNodePtr		vnp, vnp_start;

	sbp = search->sbp;
	matrix = sbp->matrix;
	str = (Uint1Ptr) search->context[context_index].query->sequence + start;
	wfp = search->wfp;
	if (wfp == NULL)
		return -2;
	lookup = wfp->lookup;
	if (lookup == NULL)
		return -3;
	wordsize = wfp->wordsize;


	sctp = SeqCodeTableFindObj(sbp->alphabet_code);
	alphabet_size=sctp->num;
	if (search->all_words == NULL)
	{
		search->all_words = BlastPopulateAllWordArrays(wordsize, alphabet_size);
		if (search->all_words == NULL)
		{
			return -1;
		}
		num_of_cols = search->all_words->num_of_cols;
		array = search->all_words->array;
	}
	else
	{
		num_of_cols = search->all_words->num_of_cols;
		array = search->all_words->array;
	}

	/* Index a specific small set, such as one db sequence. */
	if (search->all_words->specific)
	{
		len -= (wordsize-1);
		for (offset=start; offset<len; offset++, str++)
		{
			for (index1=0; index1<num_of_cols; index1++)
			{
				Boolean		ExactMatch = TRUE;
				words = array[index1];
				score = 0;
				for (index2=0; index2<wordsize; index2++)
				{
					score += matrix[*(str+index2)][*(words+index2)];
					if (*(str+index2) != *(words+index2))
					    ExactMatch = FALSE;
				}	
		     	     	if (score >= threshold || ExactMatch)
			     	{
					lookup_add(lookup, (CharPtr) words, offset+wordsize-1, context_index);
					saved_index = TRUE;
				}
			}
		}
		if (saved_index)
		    return 0;
		else
		    return 1;
	}

	s_string_start = s_string = MemNew((wordsize+2)*sizeof(Uint1));

	if (s_string_start == NULL)
		return -1;

/* Amounts to advance loops if the same character is to be checked again. */
	loop_increment=(long) (Nlm_Powi((Nlm_FloatHi)alphabet_size,(wordsize-2)));
	loop_increment2=loop_increment/alphabet_size;
/* Shorten len so up to the last complete word is checked. */
	len -= (wordsize-1);

	vnp_start = NULL;
	vnp = search->context[context_index].location;
	if (vnp == NULL)
	{
		ValNodeAddInt(&vnp, 1, -1);
		ValNodeAddInt(&vnp, 0, len+wordsize);
		vnp_start = vnp;
	}

        while (vnp)
	{
		if (vnp->choice == 1)
		{
			start = vnp->data.intvalue + 1;
			vnp = vnp->next;
			if (vnp == NULL)
				stop = len;
		}
		if (vnp && vnp->choice == 0)
		{
			stop = vnp->data.intvalue - (wordsize-1);
			vnp = vnp->next;
		}

		stop = MIN(stop, len);

		str = (Uint1Ptr) search->context[context_index].query->sequence + start;
		
	for (offset=start; offset<stop; offset++, str++)
	{
/* Put query into the lookup table, after checking that word would give
a positive value. */
		best_total=0;
		for (index1=0; index1<wordsize; index1++)
		{
		    best_total += matrix[(Int4) *(str+index1)][(Int4) *(str+index1)];
		}
		if (best_total > 0)
		{
			lookup_add(lookup, (CharPtr) str, offset+wordsize-1, context_index);
			saved_index = TRUE;
		}

/* Check if a match with a non-identical word could give a score above T. */
		best_total=0;
		for (index1=0; index1<wordsize; index1++)
		{
			best_total += sbp->maxscore[str[index1]];
		}

		if (best_total < threshold)
		{	/* no chance of a match! */
			continue;
		}

		delta_score = best_total-threshold;

/* pick a last_char that is at end of the array, could this be improved? */
		last_char=array[num_of_cols-1][wordsize-2];
		last_char2=array[num_of_cols-1][wordsize-2];

		for (index1=0; index1<num_of_cols; index1++)
		{
			words = array[index1];

/* 
only do this check if the letter has changed from last time.  See if
the new letter, matched with the first letter of the word, changes the
total possible score to below threshold.  If so, move ahead to the next letter.
This is repeated with the second letter in the word.

The order of the letters in the first and second columns of array is
important here!								
*/
			if (last_char != *words)
			{
				last_char = *words;
				first_score = matrix[(Int4) *str][(Int4) *words];
				diff = delta_score + first_score - sbp->maxscore[*str];
				if (diff < 0)
				{ 
/* index1 should be advanced by loop_increment, decrement by one as the "for" 
loop above increments by one.	*/
					index1 += loop_increment;
					index1--;
					continue;
				}
				start_score = first_score;
			}

			if (wordsize > 2 && last_char2 != *(words+1) && wordsize != 1)
			{
				last_char2 = *(words+1);
				second_score = matrix[(Int4) *(str+1)][(Int4) *(words+1)];
				diff2 =  second_score - sbp->maxscore[*(str+1)];
				diff2 += diff;
				if (diff2 < 0)
				{
/* index1 should be advanced by loop_increment2, decrement by one as the "for" 
loop above increments by one.	*/
					index1 += loop_increment2;
					index1--;
					continue;
				}
				start_score = second_score+first_score;
			}

			start_score2 = start_score;

			for (index2=2; index2<wordsize-1; index2++)
			{
				start_score2 += matrix[(Int4) *(str+index2)][*(words+index2)];
			}

			for (index2=0; index2<alphabet_size; index2++)
			{
			     score = start_score2;
			     score += matrix[(Int4) *(str+wordsize-1)][index2];

		     	     if (score >= threshold)
			     {
				exact_match=TRUE;
				for (index3=0; index3<wordsize-1; index3++)
				{
					if (*(str+index3) != *(words+index3))
					{
						exact_match=FALSE;
						break;
					}
				}
				if (*(str+wordsize-1) != index2)
				{
					exact_match=FALSE;
				}
						
				if (exact_match == FALSE)
				{
					s_string = s_string_start;
					for (index3=0; index3<wordsize-1; index3++)
					{
						*s_string = *(words+index3);
						s_string++;
					}
					*s_string = index2;
					lookup_add(lookup, (CharPtr) s_string_start, offset+wordsize-1, context_index);
					saved_index = TRUE;
				}
			     }
			}
		}
	}
	}

	if (vnp_start)
	{
		vnp_start = ValNodeFree(vnp_start);
	}

	s_string_start = MemFree(s_string_start);

	if (saved_index == FALSE)
		return 1;

	return 0;
}

/*AAS position-based version of BlastFindWords*/
Int2 LIBCALL
BlastNewFindWords(BlastSearchBlkPtr search, Int4 start, Int4 len, BLAST_Score threshold, Int1 context_index) 

{
	register Uint1		last_char, last_char2;
	Uint1Ptr		words, PNTR array;
	Uint1Ptr		s_string_start, s_string;
	register Uint1Ptr	str;
	BLAST_Score		best_total, delta_score, diff, diff2, first_score;
	BLAST_Score		second_score, start_score, start_score2, score;
	BLAST_ScoreBlkPtr 	sbp;
	register BLAST_ScorePtr PNTR	matrix;
	BLAST_WordFinderPtr	wfp;
	Boolean			exact_match;
	LookupTablePtr		lookup;
	register Int4		index1, index3, offset;
	register Int1		index2;
	Int4			num_of_cols, alphabet_size, wordsize;
	Int4 			loop_increment, loop_increment2;
	SeqCodeTablePtr 	sctp;

	sbp = search->sbp;
	matrix = sbp->matrix;
	str = (Uint1Ptr) search->context[context_index].query->sequence + start;
	wfp = search->wfp;
	lookup = wfp->lookup;
	wordsize = wfp->wordsize;

	sctp = SeqCodeTableFindObj(sbp->alphabet_code);
	alphabet_size=sctp->num;
	if (search->all_words == NULL)
	{
		search->all_words = BlastPopulateAllWordArrays(wordsize, alphabet_size);
		if (search->all_words == NULL)
		{
			return -1;
		}
		num_of_cols = search->all_words->num_of_cols;
		array = search->all_words->array;
	}
	else
	{
		num_of_cols = search->all_words->num_of_cols;
		array = search->all_words->array;
	}

	/* Index a specific small set, such as one db sequence. */
	if (search->all_words->specific)
	{
		len -= (wordsize-1);
		for (offset=start; offset<len; offset++, str++)
		{
			for (index1=0; index1<num_of_cols; index1++)
			{
				words = array[index1];
				score = 0;
				for (index2=0; index2<wordsize; index2++)
				{
					score += MtrxScorePosSearch(search->sbp,
							offset + index2,
							*(words+index2));
				}	
		     	     	if (score >= threshold)
			     	{
					lookup_add(lookup, (CharPtr) words, offset+wordsize-1, context_index);
				}
			}
		}
		return 0;
	}

	s_string_start = s_string = MemNew((wordsize+2)*sizeof(Uint1));

	if (s_string_start == NULL)
		return 1;

/* Amounts to advance loops if the same character is to be checked again. */
	loop_increment=(long) (Nlm_Powi((Nlm_FloatHi)alphabet_size,(wordsize-2)));
	loop_increment2=loop_increment/alphabet_size;
/* Shorten len so up to the last complete word is checked. */
	len -= (wordsize-1);
	for (offset=start; offset<len; offset++, str++)
	{
/* Put query into the lookup table, after checking that word would give
a positive value. */
		best_total=0;
		for (index1=0; index1<wordsize; index1++)
		{
		    best_total += MtrxScorePosSearch(search->sbp,
					offset+index1,(Int4) *(str+index1));
		}
		if (best_total > 0)
			lookup_add(lookup, (CharPtr) str, offset+wordsize-1, context_index);

/* Check if a match with a non-identical word could give a score above T. */
		best_total=0;
		for (index1=0; index1<wordsize; index1++)
		{
			best_total += sbp->maxscore[str[index1]];
		}

		if (best_total < threshold)
		{	/* no chance of a match! */
			continue;
		}

		delta_score = best_total-threshold;

/* pick a last_char that is at end of the array, could this be improved? */
		last_char=array[num_of_cols-1][wordsize-2];
		last_char2=array[num_of_cols-1][wordsize-2];

		for (index1=0; index1<num_of_cols; index1++)
		{
			words = array[index1];

/* 
only do this check if the letter has changed from last time.  See if
the new letter, matched with the first letter of the word, changes the
total possible score to below threshold.  If so, move ahead to the next letter.
This is repeated with the second letter in the word.

The order of the letters in the first and second columns of array is
important here!								
*/
			if (last_char != *words)
			{
				last_char = *words;
				first_score = MtrxScorePosSearch(search->sbp,
					offset,(Int4) *words);
				diff = delta_score + first_score - sbp->maxscore[*str];
				if (diff < 0)
				{ 
/* index1 should be advanced by loop_increment, decrement by one as the "for" 
loop above increments by one.	*/
					index1 += loop_increment;
					index1--;
					continue;
				}
				start_score = first_score;
			}

			if (last_char2 != *(words+1) && wordsize != 1)
			{
				last_char2 = *(words+1);
				second_score = MtrxScorePosSearch(search->sbp,
						offset+1,(Int4) *(words+1));
				diff2 =  second_score - sbp->maxscore[*(str+1)];
				diff2 += diff;
				if (diff2 < 0)
				{
/* index1 should be advanced by loop_increment2, decrement by one as the "for" 
loop above increments by one.	*/
					index1 += loop_increment2;
					index1--;
					continue;
				}
				start_score = second_score+first_score;
			}

			start_score2 = start_score;

			for (index2=2; index2<wordsize-1; index2++)
			{
				start_score2 += MtrxScorePosSearch(search->sbp,
						offset+index2,*(words+index2));
			}

			for (index2=0; index2<alphabet_size; index2++)
			{
			     score = start_score2;
			     score += MtrxScorePosSearch(search->sbp,
					offset+wordsize-1,index2);

		     	     if (score >= threshold)
			     {
				exact_match=TRUE;
				for (index3=0; index3<wordsize-1; index3++)
				{
					if (*(str+index3) != *(words+index3))
					{
						exact_match=FALSE;
						break;
					}
				}
				if (*(str+wordsize-1) != index2)
				{
					exact_match=FALSE;
				}
						
				if (exact_match == FALSE)
				{
					s_string = s_string_start;
					for (index3=0; index3<wordsize-1; index3++)
					{
						*s_string = *(words+index3);
						s_string++;
					}
					*s_string = index2;
					lookup_add(lookup, (CharPtr) s_string_start, offset+wordsize-1, context_index);

				}
			     }
			}
		}
	}

	s_string_start = MemFree(s_string_start);

	return 0;
}

/*******************************************************************

This function allocates and populates an array containing every possible
letter combination for an alphabet of size alphabet_size for the first
wordsize-1 letters.   The last letter of the word is done on the fly.
The approach is best described with a table that demonstrates
the strategy with a two-letter alphabet and a wordsize of three:

	index   1	2
	col.	1	0

		A	A
		A	B
		B	A
		B	B
		A	A
		A	B
		B	A
		B	B

col. 0: basic pattern repeated N**(W-1) times, where N is the size of the
	alphabet and W is the wordsize.
col. 1: AABB repeated N**(W-2) times.

Each pattern is repeated N**(W-1-C) times, where "C" is the column number.
The number of repeats of a given letter is N**C.
The total number of rows in the array is then N**(W-1-C) * N**C * N = N**W,
as we expect.

NOTE:  The order of the columns is important, as it is used in
BlastWFContextNeighborhoodAdd above.  In particular it is useful to have
all the letters grouped together.
*********************************************************************/

BlastAllWordPtr
BlastPopulateAllWordArrays(Int4 wordsize, Int4 alphabet_size)

{
	BlastAllWordPtr all_words;
	Uint1Ptr *array_ptr, *array;
	register Int4 index, index1, index3, num_of_cols, times, repeat_num;
	register Int1 index2;
	num_of_cols = (Int4) Nlm_Powi((Nlm_FloatHi)alphabet_size, wordsize-1);
	array = (Uint1Ptr *) MemNew(num_of_cols*sizeof(Uint1Ptr));

	for (index=0; index<num_of_cols; index++)
	{
	    array[index] = (Uint1Ptr) MemNew((wordsize-1)*sizeof(Uint1));
	}

	for (index=0; index<wordsize-1; index++)
	{
	    array_ptr = array;
	    repeat_num= (Int4) Nlm_Powi((Nlm_FloatHi)alphabet_size,(wordsize-2-index));
	    times = (Int4) Nlm_Powi((Nlm_FloatHi)alphabet_size, index);
	    for (index1=0; index1<times; index1++)
	    {
	    	for (index2=0; index2<alphabet_size; index2++)
	    	{
	    	    for (index3=0; index3<repeat_num; index3++)
		    {
			(*array_ptr)[index] = index2;
			array_ptr++;
		    }
		 }
	     }
	}

	all_words = BlastAllWordNew(num_of_cols, wordsize, TRUE, FALSE);
	if (all_words)
	{
		all_words->array = array;
	}

	return all_words;
}

/**************************************************************************
*
*	Get the "ScoreSet" ScorePtr from the BLAST data, which is provided
*	by "hsp".  "score_set" should be NULL, a chain of scores is added.
**************************************************************************/

ScorePtr LIBCALL
GetScoreSetFromBlastResultHsp(BLASTResultHspPtr hsp, SeqIdPtr gi_list)

{
	ScorePtr	score_set=NULL;
	Nlm_FloatHi	prob;
	Int4		score;
	CharPtr		scoretype;
	
	score = hsp->score;
	if (score > 0)
		MakeBlastScore(&score_set, "score", 0.0, score);

	score = hsp->number;
	scoretype = "sum_n";

	if (score > 1)
		MakeBlastScore(&score_set, scoretype, 0.0, score);

	prob = hsp->p_value;
	if (hsp->number == 1)
	{
		scoretype = "p_value";
	}
	else 
	{
		scoretype = "sum_p";
	}

	if (prob >= 0.)
		MakeBlastScore(&score_set, scoretype, prob, 0);

	prob = hsp->e_value;
	if (hsp->number == 1)
	{
		scoretype = "e_value";
	}
	else 
	{
		scoretype = "sum_e";
	}
	if (prob >= 0.)
		MakeBlastScore(&score_set, scoretype, prob, 0);

	prob = hsp->bit_score;
	if (prob >= 0.)
		MakeBlastScore(&score_set, "bit_score", prob, 0);

	if (hsp->ordering_method == BLAST_SMALL_GAPS)
	{
		MakeBlastScore(&score_set, "small_gap", 0.0, 1);
	}

	while (gi_list)
	{
		MakeBlastScore(&score_set, "use_this_gi", 0.0, gi_list->data.intvalue);
		gi_list = gi_list->next;
	}

	return score_set;
}

