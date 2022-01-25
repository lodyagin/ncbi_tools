/* $Id: twoseq_api.h,v 1.3 2004/05/14 17:24:03 dondosha Exp $
***************************************************************************
*                                                                         *
*                             COPYRIGHT NOTICE                            *
*                                                                         *
* This software/database is categorized as "United States Government      *
* Work" under the terms of the United States Copyright Act.  It was       *
* produced as part of the author's official duties as a Government        *
* employee and thus can not be copyrighted.  This software/database is    *
* freely available to the public for use without a copyright notice.      *
* Restrictions can not be placed on its present or future use.            *
*                                                                         *
* Although all reasonable efforts have been taken to ensure the accuracy  *
* and reliability of the software and data, the National Library of       *
* Medicine (NLM) and the U.S. Government do not and can not warrant the   *
* performance or results that may be obtained by using this software,     *
* data, or derivative works thereof.  The NLM and the U.S. Government     *
* disclaim any and all warranties, expressed or implied, as to the        *
* performance, merchantability or fitness for any particular purpose or   *
* use.                                                                    *
*                                                                         *
* In any work or product derived from this material, proper attribution   *
* of the author(s) as the source of the software or data would be         *
* appreciated.                                                            *
***************************************************************************/
/**************************************************************************

  File name: twoseq_api.h

  Author: Jason Papadopoulos

  Contents: Functions for C toolkit applications to compare two sequences
                using the rewritten blast engine 

***************************************************************************/

#ifndef _TWOSEQ_API_H_
#define _TWOSEQ_API_H_

#include <ncbi.h>
#include <objseq.h>
#include <tofasta.h>
#include <sqnutils.h>

#define MEGABLAST_CUTOFF 10000

/**
 * The type of blast search to perform. For nucleotide searches,
 * the blastn algorithm is used unless at least one input sequence
 * exceeds MEGABLAST_CUTOFF bases in size. In that case, megablast
 * is used instead. If the blast_hint is eSensitive, discontiguous
 * megablast with word size 12 is used.
 */
enum blast_type {
    eChoose = 0,        /**< Choose type of search by sequences molecule type:
                           n-n=blastn, p-p=blastp, n-p=blastx, p-n=tblastn */
    eBlastn = 1,        /**< blastn or megablast (determined automatically) */
    eBlastp = 2,        /**< blastp search between protein sequences */
    eBlastx = 3,        /**< blastx for nucleotide vs protein sequences */
    eTblastn = 4,       /**< tblastn for protein vs nucleotide sequences */
    eTblastx = 5        /**< tblastx for translated nucleotide sequences */
};

/**
 * Provide a hint on how the search is to be set up. At
 * present this only applies to nucleotide searches
 */
enum blast_hint {
    eSensitive = 0,     /**< trade off speed for sensitivity */
    eFast = 1           /**< trade off sensitivity for speed */
};

typedef enum seed_type {
   eDefaultSeedType = 0, /**< BLAST will decide which method to use based on 
                            program and other information. */
   eOneHit = 1,          /**< Require only one initial hit for extension */
   eTwoHits = 2           /**< Require more than one hit within a window 
                            for extension */
} seed_type;

/**
  * The main user-visible setup structure for the API. This
  * only makes a (small) subset of the complete options available
  */
typedef struct {
    enum blast_hint hint;       /**< for nucleotide searches, how should
                                     the search be set up? 
                                     Default = eSensitive */
    enum blast_type program;    /**< the BLAST program to use.
                                     Default = eChoose */
    char strand;                /**< For nucleotide searches, the strand
                                     of the first sequence to check: 
                                     choices are Seq_strand_{plus|minus|both}
                                     Default is Seq_strand_both */
    double cutoff_evalue;       /**< Alignments whose E value is larger than
                                     this number are discarded. Default 10.0 */
    char* matrix;               /**< The scoring matrix to use (protein
                                     searches only). NULL means "BLOSUM62".
                                     Default is NULL */
    char* filter_string;        /**< Specifies filtering to apply to the
                                     first of the two input sequences. 
                                     NULL or "T" implies DUST/SEG, "F"
                                     turns off filtering. Default = NULL */
    Int4 word_size;             /**< The word size to use. 0 chooses the
                                     default for the specified program
                                     (i.e. 3 for blastp, 11 for blastn, 
                                     28 for blastn with large sequences).
                                     Default = 0 */
    Boolean gapped_calculation; /**< Perform gapped alignments. Default = TRUE*/
    Int4 nucleotide_match;      /**< For nucleotide searches, the reward
                                     for matching letters (default 1) */
    Int4 nucleotide_mismatch;   /**< For nucleotide searches, the penalty
                                     for mismatching letters (default -3) */
    Int4 gap_open;              /**< Cost of opening a gap. Default=0, invokes 
                                     default values: 5 for nucleotide; 
                                     depends on matrix for protein search.*/
    Int4 gap_extend;            /**< Cost of extending a gap. Default=0, 
                                     invokes default values: 2 for nucleotide; 
                                     depends on matrix for protein search.*/
    Int4 gap_x_dropoff;         /**< Dropoff value for the gapped extension.
                                     Default=0, invokes default values. */
    double db_length;           /**< Database length to use in statistical 
                                     calculations. 
                                     Default=0 means "database length" is set
                                     to the subject sequence length for each
                                     subject sequence. */
    Int4 word_threshold;        /**< Threshold for finding neighboring words
                                     in protein searches. */
    seed_type init_seed_method; /**< Single-hit or multiple-hit choice of 
                                     initial seeds for extension. */
} BLAST_SummaryOptions;

/** Small structure containing the just those Karlin-Altschul parameters needed
 * for the BLAST formatting */
typedef struct BLAST_KAParameters {
   double Lambda;
   double K;
   double H;
} BLAST_KAParameters;

/** Structure holding all calculated data returned from a BLAST search other
 * than the alignment.
 */
typedef struct BLAST_SummaryReturn {
   BLAST_KAParameters* ka_params; /**< Ungapped Karlin-Altschul parameters */
   BLAST_KAParameters* ka_params_gap;/**< Gapped Karlin-Altschul parameters */
   char* params_buffer; /**< Buffer holding the bottom of BLAST report. */
} BLAST_SummaryReturn;


/**
  * Allocate storage for an API setup structure and set the
  * default options for it.
  *
  * @param options pointer to be updated with newly allocated structure [out]
  * @return 0 for successful allocation, -1 otherwise
  */
Int2 BLAST_SummaryOptionsInit(BLAST_SummaryOptions **options);

/**
  * Free the storage previously allocated for an API setup structure
  *
  * @param options pointer tothe structure to be freed [in]
  * @return always NULL
  */
BLAST_SummaryOptions* BLAST_SummaryOptionsFree(BLAST_SummaryOptions *options);

/**
  * Perform a BLAST search on the two input sequences and return
  * the list of alignments the search generates
  * @param options structure describing how the search will be configured [in]
  * @param bsp1 the first sequence to be compared. Filtering and selection
  *             of nucleotide strand apply only to this sequence [in]
  * @param bsp2 the second sequence to be compared [in]
  * @param seqalign_out the list of alignments generated by the search.
  *             If search failed or no alignments were found, set to NULL [out]
  * @return 0 for a successful search, nonzero if search failed
  */
Int2 BLAST_TwoSequencesSearch(BLAST_SummaryOptions *options,
                              Bioseq *bsp1, 
                              Bioseq *bsp2,
                              SeqAlign **seqalign_out);

/**
  * Perform a BLAST search on the two input sequences and return
  * the list of alignments the search generates
  * @param options Structure describing how the search will be configured [in]
  * @param seqloc1 The first list of sequences (queries) to be compared. 
  *                Filtering is applied only to these sequences [in]
  * @param seqloc2 The second list of sequences (subjects) to be compared [in]
  * @param seqalign_out The list of alignments generated by the search.
  *                     Alignments are sorted by query; then by subject among 
  *                     same query alignments.
  *                     If search failed or no alignments were found, 
  *                     set to NULL [out]
  * @param filter_out Masking locations [out]
  * @param extra_returns Data needed to print the bottom of BLAST report [out]
  * @return 0 for a successful search, nonzero if search failed
  */
Int2 BLAST_TwoSeqLocSets(const BLAST_SummaryOptions *options,
                         SeqLoc* seqloc1, SeqLoc* seqloc2,
                         SeqAlign **seqalign_out, 
                         SeqLoc** filter_out,
                         BLAST_SummaryReturn* *extra_returns);

#endif  /* !_TWOSEQ_API_H_ */
