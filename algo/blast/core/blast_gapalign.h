/* $Id: blast_gapalign.h,v 1.28 2003/12/10 23:14:34 dondosha Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
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

File name: blast_gapalign.h

Author: Ilya Dondoshansky

Contents: Structures and functions prototypes used for BLAST gapped extension

******************************************************************************
 * $Revision: 1.28 $
 * */

#ifndef __BLAST_GAPALIGN__
#define __BLAST_GAPALIGN__

#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_extend.h>
#include <algo/blast/core/gapinfo.h>
#include <algo/blast/core/greedy_align.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Defines extension algorithm types */
typedef enum {
    EXTEND_DYN_PROG = 1,
    EXTEND_GREEDY,
    EXTEND_GREEDY_NO_TRACEBACK,
    EXTEND_ALGO_MAX
} ExtensionAlgorithmType;

#define MB_DIAG_NEAR 30
#define MB_DIAG_CLOSE 6
#define MIN_NEIGHBOR_HSP_LENGTH 100
#define MIN_NEIGHBOR_PERC_IDENTITY 96

#define MAX_DBSEQ_LEN 5000000

/* One sequence segment within an HSP */
typedef struct BlastSeg {
   Int2 frame;  /**< Translation frame */
   Int4 offset; /**< Start of hsp */
   Int4 length; /**< Length of hsp */
   Int4 end;    /**< End of HSP */
   Int4 offset_trim; /**< Start of trimmed hsp */
   Int4 end_trim;    /**< End of trimmed HSP */
   Int4 gapped_start;/**< Where the gapped extension started. */
} BlastSeg;

/** BLAST_NUMBER_OF_ORDERING_METHODS tells how many methods are used
 * to "order" the HSP's.
*/
#define BLAST_NUMBER_OF_ORDERING_METHODS 2

/** The following structure is used in "link_hsps" to decide between
 * two different "gapping" models.  Here link is used to hook up
 * a chain of HSP's (this is a void* as _blast_hsp is not yet
 * defined), num is the number of links, and sum is the sum score.
 * Once the best gapping model has been found, this information is
 * transferred up to the BlastHSP.  This structure should not be
 * used outside of the function link_hsps.
*/
typedef struct BlastHSPLink {
   void* link[BLAST_NUMBER_OF_ORDERING_METHODS]; /**< Used to order the HSPs
                                           (i.e., hook-up w/o overlapping). */
   Int2 num[BLAST_NUMBER_OF_ORDERING_METHODS]; /**< number of HSP in the
                                                  ordering. */
   Int4 sum[BLAST_NUMBER_OF_ORDERING_METHODS]; /**< Sum-Score of HSP. */
   double xsum[BLAST_NUMBER_OF_ORDERING_METHODS]; /**< Sum-Score of HSP,
                                     multiplied by the appropriate Lambda. */
   Int4 changed;
} BlastHSPLink;

/** Structure holding all information about an HSP */
typedef struct BlastHSP {
   struct BlastHSP* next; /**< The next HSP */
   struct BlastHSP* prev; /**< The previous HSP. */
   BlastHSPLink  hsp_link;
   Boolean linked_set;        /**< Is this HSp part of a linked set? */
   Int2 ordering_method;/**< Which method (max or no max for gaps) was used? */
   Int4 num;            /**< How many HSP's make up this (sum) segment */
   Int4 sumscore;/**< Sumscore of a set of "linked" HSP's. */
   Boolean start_of_chain; /**< If TRUE, this HSP starts a chain along the
                              "link" pointer. */
   Int4 score;         /**< This HSP's raw score */
   Int4 num_ident;         /**< Number of identical base pairs in this HSP */
   double evalue;        /**< This HSP's e-value */
   BlastSeg query;            /**< Query sequence info. */
   BlastSeg subject;          /**< Subject sequence info. */
   Int2     context;          /**< Context number of query */
   GapEditBlock* gap_info; /**< ALL gapped alignment is here */
   Int4 num_ref;              /**< Number of references in the linked set */
   Int4 linked_to;            /**< Where this HSP is linked to? */
} BlastHSP;

/** Auxiliary structure for dynamic programming gapped extension */
typedef struct BlastGapDP {
  Int4 CC, DD, FF;      /**< Values for gap opening and extensions. */
} BlastGapDP;

/** Structure supporting the gapped alignment */
typedef struct BlastGapAlignStruct {
   Boolean positionBased; /**< Is this PSI-BLAST? */
   GapStateArrayStruct* state_struct; /**< Structure to keep extension 
                                                state information */
   GapEditBlock* edit_block; /**< The traceback (gap) information */
   BlastGapDP* dyn_prog; /**< Preallocated memory for the dynamic 
                              programming extension */
   GreedyAlignMem* greedy_align_mem;/**< Preallocated memory for the greedy 
                                         gapped extension */
   BlastScoreBlk* sbp; /**< Pointer to the scoring information block */
   Int4 gap_x_dropoff; /**< X-dropoff parameter to use */
   Int4 query_start, query_stop;/**< Return values: query offsets */
   Int4 subject_start, subject_stop;/**< Return values: subject offsets */
   Int4 score;   /**< Return value: alignment score */
   double percent_identity;/**< Return value: percent identity - filled only 
                               by the greedy non-affine alignment algorithm */
} BlastGapAlignStruct;

/** Initializes the BlastGapAlignStruct structure 
 * @param score_options Options related to scoring alignments [in]
 * @param ext_params Options and parameters related to gapped extension [in]
 * @param max_subject_length Maximum length of any subject sequence (needed 
 *        for greedy extension allocation only) [in]
 * @param query_length The length of the query sequence [in]
 * @param sbp The scoring information block [in]
 * @param gap_align_ptr The BlastGapAlignStruct structure [out]
*/
Int2
BLAST_GapAlignStructNew(const BlastScoringOptions* score_options, 
   BlastExtensionParameters* ext_params, 
   Uint4 max_subject_length, Int4 query_length, 
   BlastScoreBlk* sbp, BlastGapAlignStruct** gap_align_ptr);

/** Deallocates memory in the BlastGapAlignStruct structure */
BlastGapAlignStruct* 
BLAST_GapAlignStructFree(BlastGapAlignStruct* gap_align);

/** The structure to hold all HSPs for a given sequence after the gapped 
 *  alignment.
 */
typedef struct BlastHSPList {
   Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
   BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
   Int4 hspcnt; /**< Number of HSPs saved */
   Int4 allocated; /**< The allocated size of the hsp_array */
   Int4 hsp_max; /**< The maximal number of HSPs allowed to be saved */
   Boolean do_not_reallocate; /**< Is reallocation of the hsp_array allowed? */
   Boolean traceback_done; /**< Has the traceback already been done on HSPs in
                              this list? */
} BlastHSPList;

/** Creates HSP list structure with a default size HSP array */
BlastHSPList* BlastHSPListNew(void);

/** Deallocate memory for the HSP list */
BlastHSPList* BlastHSPListDestruct(BlastHSPList* hsp_list);

/** Mega BLAST function performing gapped alignment: 
 *  Sorts initial HSPs by diagonal; 
 *  For each HSP:
 *    - Removes HSP if it is included in another already extended HSP;
 *    - If required, performs ungapped extension;
 *    - Performs greedy gapped extension;
 *    - Saves qualifying HSPs with gapped information into an HSP list 
 *      structure.
 * @param program_number Not needed: added for prototype consistency.
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param gap_align A placeholder for gapped alignment information and 
 *        score block. [in] [out]
 * @param score_options Options related to scoring alignments [in]
 * @param ext_params Options related to alignment extension [in]
 * @param hit_params Options related to saving HSPs [in]
 * @param init_hitlist Contains all the initial hits [in]
 * @param hsp_list_ptr List of HSPs with full extension information [out]
*/
Int2 BLAST_MbGetGappedScore(Uint1 program_number, 
        BLAST_SequenceBlk* query, 
			    BLAST_SequenceBlk* subject,
			    BlastGapAlignStruct* gap_align,
			    const BlastScoringOptions* score_options, 
			    BlastExtensionParameters* ext_params,
			    BlastHitSavingParameters* hit_params,
			    BlastInitHitList* init_hitlist,
			    BlastHSPList** hsp_list_ptr);



/** Performs gapped extension for all non-Mega BLAST programs, given
 * that ungapped extension has been done earlier.
 * Sorts initial HSPs by score (from ungapped extension);
 * Deletes HSPs that are included in already extended HSPs;
 * Performs gapped extension;
 * Saves HSPs into an HSP list.
 * @param program_number Type of BLAST program [in]
 * @param query The query sequence block [in]
 * @param subject The subject sequence block [in]
 * @param gap_align The auxiliary structure for gapped alignment [in]
 * @param score_options Options related to scoring [in]
 * @param ext_params Options and parameters related to extensions [in]
 * @param hit_params Options related to saving hits [in]
 * @param init_hitlist List of initial HSPs (offset pairs with additional 
 *        information from the ungapped alignment performed earlier) [in]
 * @param hsp_list_ptr Structure containing all saved HSPs [out]
 */
Int2 BLAST_GetGappedScore (Uint1 program_number, BLAST_SequenceBlk* query, 
		      BLAST_SequenceBlk* subject,
		      BlastGapAlignStruct* gap_align,
		      const BlastScoringOptions* score_options, 
		      BlastExtensionParameters* ext_params,
		      BlastHitSavingParameters* hit_params,
		      BlastInitHitList* init_hitlist,
		      BlastHSPList** hsp_list_ptr);

/** Perform a gapped alignment with traceback
 * @param program Type of BLAST program [in]
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param gap_align The gapped alignment structure [in] [out]
 * @param score_options Scoring parameters [in]
 * @param q_start Offset in query where to start alignment [in]
 * @param s_start Offset in subject where to start alignment [in]
 * @param subject_length Maximal allowed extension in subject [in]
 */
Int2 BLAST_GappedAlignmentWithTraceback(Uint1 program, 
        Uint1* query, Uint1* subject, 
        BlastGapAlignStruct* gap_align, 
        const BlastScoringOptions* score_options,
        Int4 q_start, Int4 s_start, Int4 query_length, Int4 subject_length);

/** Greedy gapped alignment, with or without traceback.
 * Given two sequences, relevant options and an offset pair, fills the
 * gap_align structure with alignment endpoints and, if traceback is 
 * performed, gap information.
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param query_length The query sequence length [in]
 * @param subject The subject sequence length [in]
 * @param gap_align The structure holding various information and memory 
 *        needed for gapped alignment [in] [out]
 * @param score_options Options related to scoring alignments [in]
 * @param q_off Starting offset in query [in]
 * @param s_off Starting offset in subject [in]
 * @param compressed_subject Is subject sequence compressed? [in]
 * @param do_traceback Should traceback be saved? [in]
 */
Int2 
BLAST_GreedyGappedAlignment(Uint1* query, Uint1* subject, 
   Int4 query_length, Int4 subject_length, BlastGapAlignStruct* gap_align,
   const BlastScoringOptions* score_options, 
   Int4 q_off, Int4 s_off, Boolean compressed_subject, Boolean do_traceback);

/** Perform a gapped alignment with traceback for PHI BLAST
 * @param program Type of BLAST program [in]
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param gap_align The gapped alignment structure [in] [out]
 * @param score_options Scoring parameters [in]
 * @param q_start Offset in query where to start alignment [in]
 * @param s_start Offset in subject where to start alignment [in]
 * @param subject_length Maximal allowed extension in subject [in]
 */
Int2 PHIGappedAlignmentWithTraceback(Uint1 program, 
        Uint1* query, Uint1* subject, 
        BlastGapAlignStruct* gap_align, 
        const BlastScoringOptions* score_options,
        Int4 q_start, Int4 s_start, Int4 query_length, Int4 subject_length);

/** Convert initial HSP list to an HSP list: to be used in ungapped search.
 * Ungapped data must be available in the initial HSP list for this function 
 * to work.
 * @param init_hitlist List of initial HSPs with ungapped extension 
 *                     information [in]
 * @param subject Subject sequence block containing frame information [in]
 * @param hit_options Hit saving options [in]
 * @param hsp_list_ptr HSPs in the final form [out]
 */
Int2 BLAST_GetUngappedHSPList(BlastInitHitList* init_hitlist, 
        BLAST_SequenceBlk* subject, 
        const BlastHitSavingOptions* hit_options, 
        BlastHSPList** hsp_list_ptr);

/** Preliminary gapped alignment for PHI BLAST.
 * @param program_number Type of BLAST program [in]
 * @param query The query sequence block [in]
 * @param subject The subject sequence block [in]
 * @param gap_align The auxiliary structure for gapped alignment [in]
 * @param score_options Options related to scoring [in]
 * @param ext_params Options and parameters related to extensions [in]
 * @param hit_params Options related to saving hits [in]
 * @param init_hitlist List of initial HSPs, including offset pairs and
 *                     pattern match lengths [in]
 * @param hsp_list_ptr Structure containing all saved HSPs [out]
 */
Int2 PHIGetGappedScore (Uint1 program_number, 
        BLAST_SequenceBlk* query, 
        BLAST_SequenceBlk* subject, 
        BlastGapAlignStruct* gap_align,
        const BlastScoringOptions* score_options,
        BlastExtensionParameters* ext_params,
        BlastHitSavingParameters* hit_params,
        BlastInitHitList* init_hitlist,
        BlastHSPList** hsp_list_ptr);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_GAPALIGN__ */