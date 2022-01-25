/*  $Id: blast_stat.h,v 1.47 2004/06/10 13:22:33 madden Exp $
 * ===========================================================================
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
 * Author:  Tom Madden
 *
 */

/** @file blast_stat.h
 * Definitions and prototypes used by blast_stat.c to calculate BLAST
 * statistics. @todo FIXME: needs doxygen comments
 */

#ifndef __BLAST_STAT__
#define __BLAST_STAT__

#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_message.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
	Defines for the matrix 'preferences' (as specified by S. Altschul).
*/
#define BLAST_MATRIX_NOMINAL 0
#define BLAST_MATRIX_PREFERRED 1
#define BLAST_MATRIX_BEST 2


/* 
	Where are the BLAST matrices located?
*/
#define BLASTMAT_DIR "/usr/ncbi/blast/matrix"

/**
  Structure to hold the Karlin-Altschul parameters.
*/

typedef struct Blast_KarlinBlk {
		double	Lambda; /**< Lambda value used in statistics */
		double	K, logK; /**< K value used in statistics */
		double	H; /**< H value used in statistics */
		double	paramC;	/**< for use in seed. */
	} Blast_KarlinBlk;




/********************************************************************

	Structures relating to scoring or the BlastScoreBlk

********************************************************************/

/*
SCORE_MIN is (-2**31 + 1)/2 because it has been observed more than once that
a compiler did not properly calculate the quantity (-2**31)/2.  The list
of compilers which failed this operation have at least at some time included:
NeXT and a version of AIX/370's MetaWare High C R2.1r.
For this reason, SCORE_MIN is not simply defined to be LONG_MIN/2.
*/
#define BLAST_SCORE_MIN	INT2_MIN
#define BLAST_SCORE_MAX	INT2_MAX


/** Holds score frequencies used in calculation
of Karlin-Altschul parameters for an ungapped search.
*/
typedef struct Blast_ScoreFreq {
    Int4	score_min, score_max; /**< lowest and highest allowed scores */
    Int4	obs_min, obs_max; /**< lowest and highest observed (actual) scores */
    double	score_avg; /**< average score, must be negative for local alignment. */
    double* sprob0,* sprob; /**< arrays for frequency of given score, sprob is shifted down by score_min. */
} Blast_ScoreFreq;

#define BLAST_MATRIX_SIZE 32

/* Remove me */
typedef struct SBLASTMatrixStructure {
    Int4 *matrix[BLAST_MATRIX_SIZE];
    Int4 long_matrix[BLAST_MATRIX_SIZE*BLAST_MATRIX_SIZE]; /* not used */
} SBLASTMatrixStructure;

/** Structure used for scoring calculations.
*/
typedef struct BlastScoreBlk {
	Boolean		protein_alphabet; /**< TRUE if alphabet_code is for a 
protein alphabet (e.g., ncbistdaa etc.), FALSE for nt. alphabets. */
	Uint1		alphabet_code;	/**< NCBI alphabet code. */
	Int2 		alphabet_size;  /**< size of alphabet. */
	Int2 		alphabet_start;  /**< numerical value of 1st letter. */
	SBLASTMatrixStructure* matrix_struct;	/**< Holds info about matrix. */
	Int4 **matrix;  /**< Substitution matrix */
	Int4 **posMatrix;  /**< Sub matrix for position depend BLAST. */
   double karlinK; /**< Karlin-Altschul parameter associated with posMatrix */
	Int2	mat_dim1, mat_dim2;	/**< dimensions of matrix. */
	Int4 *maxscore; /**< Max. score for each letter */
	Int4	loscore, hiscore; /**< Min. & max. substitution scores */
	Int4	penalty, reward; /**< penalty and reward for blastn. */
        double  scale_factor; /**< multiplier for all cutoff and dropoff scores */
	Boolean		read_in_matrix; /**< If TRUE, matrix is read in, otherwise
					produce one from penalty and reward above. */
	Blast_ScoreFreq** sfp;	/**< score frequencies. */
	double **posFreqs; /**<matrix of position specific frequencies*/
	/* kbp & kbp_gap are ptrs that should be set to kbp_std, kbp_psi, etc. */
	Blast_KarlinBlk** kbp; 	/**< Karlin-Altschul parameters. */
	Blast_KarlinBlk** kbp_gap; /**< K-A parameters for gapped alignments. */
	/* Below are the Karlin-Altschul parameters for non-position based ('std')
	and position based ('psi') searches. */
	Blast_KarlinBlk **kbp_std,	
                    **kbp_psi,	
                    **kbp_gap_std,
                    **kbp_gap_psi;
	Blast_KarlinBlk* 	kbp_ideal;	/**< Ideal values (for query with average database composition). */
	Int4 number_of_contexts;	/**< Used by sfp and kbp, how large are these*/
	char* 	name;		/**< name of matrix. */
	Uint1* 	ambiguous_res;	/**< Array of ambiguous res. (e.g, 'X', 'N')*/
	Int2		ambig_size,	/**< size of array above. FIXME: not needed here? */
			ambig_occupy;	/**< How many occupied? */
	ListNode*	comments;	/**< Comments about matrix. */
	Int4    	query_length;   /**< the length of the query. */
	Int8	effective_search_sp;	/**< product of above two */
} BlastScoreBlk;

/** 
Stores the letter frequency of a sequence or database.
*/
typedef struct Blast_ResFreq {
    Uint1   alphabet_code;    /**< indicates alphabet. */
    double* prob;	      /**< letter probs, (possible) non-zero offset. */
    double* prob0;            /**< probs, zero offset. */
} Blast_ResFreq;

/**
 *   Allocates and initializes BlastScoreBlk
 * @param alphabet either BLASTAA_SEQ_CODE or BLASTNA_SEQ_CODE [in]
 * @param number_of_contexts how many strands or sequences [in]
 * @return BlastScoreBlk*
*/
BlastScoreBlk* BlastScoreBlkNew (Uint1 alphabet, Int4 number_of_contexts);

/** Deallocates BlastScoreBlk as well as all associated structures.
 * @param sbp BlastScoreBlk to be deallocated [in]
 * @return NULL pointer.
 */
BlastScoreBlk* BlastScoreBlkFree (BlastScoreBlk* sbp);

/* FIXME make private? */
Int2 BlastScoreBlkMatrixLoad(BlastScoreBlk* sbp);

/** Set the ambiguous residue (e.g, 'N', 'X') in the BlastScoreBlk*.
 *  Convert from ncbieaa to sbp->alphabet_code (i.e., ncbistdaa) first.
 *
 * @param sbp the object to be modified [in|out]
 * @param ambiguous_res the residue to be set on the BlastScoreBlk
 * @return zero on success, others on error
 */
Int2 BLAST_ScoreSetAmbigRes (BlastScoreBlk* sbp, char ambiguous_res);


/** Calculate the Karlin parameters.  This function should be called once
 *  for each context, or frame translated.
 * @param sbp the object to be modified [in|out]
 * @param string the query sequence [in]
 * @param length length of above sequence [in]
 * @param context_number which element in various arrays [in]
 * @return zero on success.
 */
Int2 BLAST_ScoreBlkFill (BlastScoreBlk* sbp, char* string, Int4 length, Int4 context_number);
 
/** This function fills in the BlastScoreBlk structure.  
 * Tasks are:
 *	-read in the matrix
 *	-set maxscore
 * @param sbp Scoring block [in] [out]
 * @param matrix Full path to the matrix in the directory structure [in]
*/
Int2 BLAST_ScoreBlkMatFill (BlastScoreBlk* sbp, char* matrix);
 
/** Callocs a Blast_KarlinBlk
 * @return pointer to the Blast_KarlinBlk
*/
Blast_KarlinBlk* Blast_KarlinBlkCreate (void);

/** Deallocates the KarlinBlk
 * @param kbp KarlinBlk to be deallocated [in]
 * @return NULL
*/
Blast_KarlinBlk* Blast_KarlinBlkDestruct(Blast_KarlinBlk* kbp);


Int2 Blast_KarlinBlkGappedCalc (Blast_KarlinBlk* kbp, Int4 gap_open, 
        Int4 gap_extend, Int4 decline_align, char* matrix_name, 
        Blast_Message** error_return);


/** Calculates the standard Karlin parameters.  This is used
 *       if the query is translated and the calculated (real) Karlin
 *       parameters are bad, as they're calculated for non-coding regions.
 * @param sbp ScoreBlk used to calculate "ideal" values. [in]
*/
Blast_KarlinBlk* Blast_KarlinBlkIdealCalc(BlastScoreBlk* sbp);


/** Fills KarlinBlk pointers in BlastScoreBlk with "ideal" values if the 
 * ideal Lambda is less than the actual Lambda.  This happens if 
 * if the query is translated and the calculated (real) Karlin
 *
 * parameters are bad, as they're calculated for non-coding regions.
 * @param sbp the object to be modified [in|out]
 * @param context_start first context to start with [in]
 * @param context_end last context to work on [in]
 * @return zero on success
 */
Int2 Blast_KarlinBlkStandardCalc(BlastScoreBlk* sbp, Int4 context_start, 
                                 Int4 context_end);

/*
        Attempts to fill KarlinBlk for given gap opening, extensions etc.
        Will return non-zero status if that fails.

        return values:  -1 if matrix_name is NULL;
                        1 if matrix not found
                        2 if matrix found, but open, extend etc. values not supported.
*/
Int2 Blast_KarlinkGapBlkFill(Blast_KarlinBlk* kbp, Int4 gap_open, Int4 gap_extend, Int4 decline_align, char* matrix_name);

/** Prints a messages about the allowed matrices, BlastKarlinkGapBlkFill should return 1 before this is called. 
 * @param matrix the matrix to print a message about [in]
 * @return the message
 */
char* BLAST_PrintMatrixMessage(const char *matrix);

/** Prints a messages about the allowed open etc values for the given matrix, 
 * BlastKarlinkGapBlkFill should return 2 before this is called. 
 * @param matrix name of the matrix [in]
 * @param gap_open gap existence cost [in]
 * @param gap_extend cost to extend a gap by one [in]
 * @param decline_align cost of declining to align [in]
 * @return message
 */
char* BLAST_PrintAllowedValues(const char *matrix, Int4 gap_open, Int4 gap_extend, Int4 decline_align);

/** Calculates the parameter Lambda given an initial guess for its value */
double
Blast_KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess);

/** Calculates the Expect value based upon the search space and some Karlin-Altschul 
 * parameters.  It is "simple" as it does not use sum-statistics.
 * @param S the score of the alignment. [in]
 * @param kbp the Karlin-Altschul parameters. [in]
 * @param searchsp total search space to be used [in]
 * @return the expect value
 */

double BLAST_KarlinStoE_simple (Int4 S, Blast_KarlinBlk* kbp, Int8  searchsp);
double BLAST_GapDecayDivisor(double decayrate, unsigned nsegs );

/** Calculate the cutoff score from the expected number of HSPs or vice versa.
 * @param S The calculated score [in] [out]
 * @param E The calculated e-value [in] [out]
 * @param kbp The Karlin-Altschul statistical parameters [in]
 * @param searchsp The effective search space [in]
 * @param dodecay Use gap decay feature? [in]
 * @param gap_decay_rate Gap decay rate to use, if dodecay is set [in]
 */
Int2 BLAST_Cutoffs (Int4 *S, double* E, Blast_KarlinBlk* kbp, 
                    Int8 searchsp, Boolean dodecay, double gap_decay_rate);

/* Functions to calculate SumE (for large and small gaps). */
double BLAST_SmallGapSumE (Blast_KarlinBlk* kbp, Int4 gap, Int2 num,  double xsum, Int4 query_length, Int4 subject_length, double weight_divisor);

double BLAST_UnevenGapSumE (Blast_KarlinBlk* kbp, Int4 p_gap, Int4 n_gap, Int2 num,  double xsum, Int4 query_length, Int4 subject_length, double weight_divisor);

double BLAST_LargeGapSumE (Blast_KarlinBlk* kbp, Int2 num,  double xsum, Int4 query_length, Int4 subject_length, double weight_divisor );

/*
Obtains arrays of the allowed opening and extension penalties for gapped BLAST for
the given matrix.  Also obtains arrays of Lambda, K, and H.  The pref_flags field is
used for display purposes, with the defines: BLAST_MATRIX_NOMINAL, BLAST_MATRIX_PREFERRED, and
BLAST_MATRIX_BEST.

Any of these fields that
are not required should be set to NULL.  The Int2 return value is the length of the
arrays.
*/

void BLAST_GetAlphaBeta (char* matrixName, double *alpha,
double *beta, Boolean gapped, Int4 gap_open, Int4 gap_extend);

Int4 ** RPSCalculatePSSM(double scalingFactor, Int4 rps_query_length, 
                   Uint1 * rps_query_seq, Int4 db_seq_length, Int4 **posMatrix);

Int4
BLAST_ComputeLengthAdjustment(double K,
                              double logK,
                              double alpha_d_lambda,
                              double beta,
                              Int4 query_length,
                              Int8 db_length,
                              Int4 db_num_seqs,
                              Int4 * length_adjustment);


/** Allocates a new Blast_ResFreq structure and fills in the prob element
    based upon the contents of sbp.
 * @param sbp The BlastScoreBlk* used to init prob [in]
*/
Blast_ResFreq* Blast_ResFreqNew(const BlastScoreBlk* sbp);

/** Deallocates Blast_ResFreq and prob0 element.
 * @param rfp the Blast_ResFreq to be deallocated.
*/
Blast_ResFreq* Blast_ResFreqDestruct(Blast_ResFreq* rfp);


/** Calculates residues frequencies given a standard distribution.
 * @param sbp the BlastScoreBlk provides information on alphabet.
 * @param rfp the prob element on this Blast_ResFreq is used.
*/
Int2 Blast_ResFreqStdComp(const BlastScoreBlk* sbp, Blast_ResFreq* rfp);

/** Creates a new structure to keep track of score frequencies for a scoring
 * system.
 * @param score_min Minimum score [in]
 * @param score_max Maximum score [in]
 */
Blast_ScoreFreq*
Blast_ScoreFreqNew(Int4 score_min, Int4 score_max);

/** Deallocates the score frequencies structure 
 * @param sfp the structure to deallocate [in]
 * @return NULL
 */
Blast_ScoreFreq*
Blast_ScoreFreqDestruct(Blast_ScoreFreq* sfp);

/** Fills a buffer with the 'standard' alphabet 
 * (given by STD_AMINO_ACID_FREQS[index].ch).
 *
 * @return Number of residues in alphabet or negative returns upon error.
 */
Int2
Blast_GetStdAlphabet(Uint1 alphabet_code, Uint1* residues, 
                     Uint4 residues_size);

/* Please see comment on blast_stat.c  */
Int2
Blast_KarlinBlkCalc(Blast_KarlinBlk* kbp, Blast_ScoreFreq* sfp);

/**  Given a sequence of 'length' amino acid residues, compute the
 *   probability of each residue and put that in the array resProb
 *
 * @param sequence the sequence to be computed upon [in]
 * @param length the length of the sequence [in]
 * @param resProb the object to be filled in [in|out]
 */
void
Blast_FillResidueProbability(const Uint1* sequence, Int4 length, double * resProb);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_STAT__ */
