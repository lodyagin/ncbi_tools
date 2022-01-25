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

File name: blastkar.c

Author: Tom Madden

Contents: Functions to calculate BLAST probabilities etc.

Detailed Contents:

	- allocate and deallocate structures used by BLAST to calculate
	probabilities etc.

	- calculate residue frequencies for query and "average" database.

	- read in matrix.

        - calculate sum-p from a collection of HSP's, for both the case
	of a "small" gap and a "large" gap, when give a total score and the
	number of HSP's.

	- calculate expect values for p-values.

	- calculate pseuod-scores from p-values.

******************************************************************************/

/* $Revision: 6.32 $ */
/* $Log: blastkar.c,v $
/* Revision 6.32  1999/08/20 14:42:17  madden
/* Changed Robinson frequencies per Stephen A. request
/*
/* Revision 6.31  1999/08/05 13:13:34  madden
/* Improved help messages for permitted matrices and gap values
/*
/* Revision 6.30  1999/07/30 13:25:51  shavirin
/* Fixed bug in the function BLAST_MatrixFill which created wrong matrix size.
/*
 * Revision 6.29  1998/12/31 18:17:04  madden
 * Added strand option
 *
 * Revision 6.28  1998/09/28 12:28:50  madden
 * Protection for a DEC Alpha problem with HUGE_VAL
 *
 * Revision 6.27  1998/09/24 15:26:36  egorov
 * Fix lint complaints
 *
 * Revision 6.26  1998/08/11 12:57:23  madden
 * Changed importance of BlastKarlinBlkGappedCalc error message
 *
 * Revision 6.25  1998/07/17 15:39:57  madden
 * Changes for Effective search space.
 *
 * Revision 6.24  1998/06/12 15:52:50  madden
 * Fixed warnings
 *
 * Revision 6.23  1998/06/02 21:21:20  madden
 * Changes for DNA matrices
 *
 * Revision 6.22  1998/05/03 17:56:45  madden
 * Fixed memory leaks
 *
 * Revision 6.21  1998/04/24 21:51:40  madden
 * Check return value on BlastKarlinBlkCalc
 *
 * Revision 6.20  1998/04/24 19:27:54  madden
 * Added BlastKarlinBlkStandardCalcEx for ideal KarlinBlk
 *
 * Revision 6.19  1998/04/10 20:57:55  madden
 * Added data for BLOSUM80 and BLOSUM45
 *
 * Revision 6.18  1998/04/10 15:05:48  madden
 * Added pref_flags return by value to BlastKarlinGetMatrixValues
 *
 * Revision 6.17  1998/04/02 21:12:29  madden
 * Added infinite values to the arrays returned by BlastKarlinGetMatrixValues
 *
 * Revision 6.16  1998/03/16 17:41:56  madden
 * Fixed leaks
 *
 * Revision 6.15  1998/03/09 17:15:01  madden
 * Added BlastKarlinGetMatrixValues
 *
 * Revision 6.14  1998/02/26 22:34:32  madden
 * Changes for 16 bit windows
 *
 * Revision 6.13  1998/02/06 18:28:05  madden
 * Added functions to produce pseudo-scores from p and e values
 *
 * Revision 6.12  1997/12/15 21:55:44  madden
 * Added new parameters for BLOSUM62_20
 *
 * Revision 6.11  1997/12/11 22:21:02  madden
 * Removed unused variables
 *
 * Revision 6.10  1997/12/10 22:41:58  madden
 * proper casting done
 *
 * Revision 6.9  1997/12/01 22:07:55  madden
 * Added missing values, fixed UMR
 *
 * Revision 6.8  1997/11/14 21:29:58  madden
 * Stephens new parameters
 *
 * Revision 6.7  1997/11/14 17:14:52  madden
 * Added Karlin parameter to matrix structure
 *
 * Revision 6.6  1997/11/07 22:27:29  madden
 * Fix in BLAST_MatrixFill
 *
 * Revision 6.5  1997/11/07 00:48:07  madden
 * Added defintitions and functions for BLAST_Matrix
 *
 * Revision 6.4  1997/10/30 15:41:01  madden
 * Casts and fixes for DEC alpha
 *
 * Revision 6.3  1997/10/22 21:46:59  madden
 * Added BLOSUM90, changed to better BLOSUM62 values
 *
 * Revision 6.2  1997/10/09 19:47:35  madden
 * changes for 1/20th bit BLOSUM62 matrices
 *
 * Revision 6.1  1997/10/06 17:59:17  madden
 * Added checks to BLAST_ScoreBlkDestruct
 *
 * Revision 6.0  1997/08/25 18:52:37  madden
 * Revision changed to 6.0
 *
 * Revision 1.39  1997/07/18 20:56:13  madden
 * Added more BLOSUM50 values
 *
 * Revision 1.38  1997/07/14 20:11:17  madden
 * Removed unused variables
 *
 * Revision 1.37  1997/07/14 15:30:59  madden
 * Changed call to BlastKarlinBlkGappedCalc
 *
 * Revision 1.36  1997/06/23 20:48:37  madden
 * Added support for BLOSUM50 and PAM250 matrices
 *
 * Revision 1.35  1997/06/23 17:53:07  madden
 * BlastKarlinBlkGappedCalc checks for valid values if KarlinBlkPtr NULL
 *
 * Revision 1.34  1997/05/01 15:53:13  madden
 * Addition of extra KarlinBlk's for psi-blast
 *
 * Revision 1.33  1997/04/22  16:36:49  madden
 * Added Robinson amino acid frequencies.
 *
 * Revision 1.32  1997/04/03  19:48:13  madden
 * Changes to use effective database length instead of the length of each
 * sequence in statistical calculations.
 *
 * Revision 1.31  1997/03/07  22:24:43  madden
 * Closed File* for matrix after use.
 *
 * Revision 1.30  1997/02/20  21:33:51  madden
 * Added FindPath call to blastkar.c to look for matrix.
 *
 * Revision 1.29  1997/02/06  23:15:43  madden
 * Added additional gap extension and opening penalties.
 *
 * Revision 1.28  1997/02/05  19:54:59  madden
 * *** empty log message ***
 *
 * Revision 1.27  1997/01/16  22:58:08  madden
 * Secured BlastScoreFreqDestruct against NULL pointers.
 *
 * Revision 1.26  1997/01/13  17:15:29  madden
 * Save matrix name for blastn type matrices.
 *
 * Revision 1.25  1996/12/11  17:59:42  madden
 * Fixed purify leaks.
 *
 * Revision 1.24  1996/12/10  17:30:59  madden
 * Changed statistics for gapped blastp
 *
 * Revision 1.23  1996/12/03  19:13:47  madden
 * Made BlastRes functions non-static.
 *
 * Revision 1.22  1996/11/27  16:39:11  madden
 * Added NULLB to end of array, purify nit.
 *
 * Revision 1.21  1996/11/26  20:38:01  madden
 * Removed unused variable.
 *
 * Revision 1.20  1996/11/26  19:55:25  madden
 * Added check for matrices in standard places.
 *
 * Revision 1.19  1996/11/18  19:32:09  madden
 * Fixed uninitialized variable found by CodeWarrior.
 *
 * Revision 1.18  1996/11/18  14:49:26  madden
 * Rewrote BlastScoreSetAmbigRes to take multiple ambig. chars.
 *
 * Revision 1.17  1996/11/08  21:45:03  madden
 * Fix for blastn matrices.
 *
 * Revision 1.16  1996/11/07  22:31:15  madden
 * Corrected Nucl. matrix for ambiguity characters.
 *
 * Revision 1.15  1996/10/04  20:12:26  madden
 * Fixed memory leaks found by purify.
 *
 * Revision 1.14  1996/10/03  20:49:29  madden
 * Added function BlastKarlinBlkStandardCalc to calculate standard
 * Karlin parameters for blastx and tblastx.
 *
 * Revision 1.13  1996/10/03  18:04:48  madden
 * Each context (or frame) now uses "proper" Karlin parameters.
 *
 * Revision 1.12  1996/10/03  13:07:55  madden
 * Added return statement.
 *
 * Revision 1.11  1996/10/02  20:00:38  madden
 * Calc. different Karlin parameters for each frame.
 *
 * Revision 1.10  1996/10/01  21:24:02  madden
 * Corrected statistics for partial blastn matching.
 *
 * Revision 1.9  1996/09/30  21:56:12  madden
 * Replaced query alphabet of ncbi2na with blastna alphabet.
 *
 * Revision 1.8  1996/09/11  20:08:55  madden
 * *** empty log message ***
 *
 * Revision 1.6  1996/09/11  19:14:09  madden
 * Changes to blastn matrix.
 *
 * Revision 1.5  1996/09/10  19:40:35  madden
 * Added function BlastScoreBlkMatCreate for blastn matrices.
 *
 * Revision 1.4  1996/08/23  19:07:01  madden
 * Adjust changes for NT warning to give correct results.
 *
 * Revision 1.3  1996/08/23  15:32:30  shavirin
 * Fixed a lot of NT compiler warnings about type mismatch
 *
 * Revision 1.2  1996/08/22  20:38:11  madden
 * Changed all doubles and floats to Nlm_FloatHi.
 *
 * Revision 1.1  1996/08/05  19:47:42  madden
 * Initial revision
 *
 * Revision 1.26  1996/06/21  15:23:54  madden
 * Corelibed BlastSumP.
 *
 * Revision 1.25  1996/06/21  15:14:55  madden
 * Made some functions static, added LIBCALL to others.
 *
 * Revision 1.24  1996/06/20  16:52:23  madden
 * Changed "pow" to "Nlm_Powi".
 *
 * Revision 1.23  1996/06/11  15:42:46  madden
 * Replaced strncasempc with toolbox function.
 *
 * Revision 1.22  1996/05/29  12:44:04  madden
 * Added function BlastRepresentativeResidues.
 *
 * Revision 1.21  1996/05/22  20:38:05  madden
 * *** empty log message ***
 *
 * Revision 1.20  1996/05/22  20:21:23  madden
 * removed unused variables.
 *
 * Revision 1.19  1996/05/20  21:18:51  madden
 * Added BLASTMatrixStructure.
 *
 * Revision 1.18  1996/05/16  19:50:15  madden
 * Added documentation block.
 *
 * Revision 1.17  1996/05/14  21:30:37  madden
 * *** empty log message ***
 *
 * Revision 1.16  1996/04/16  15:33:41  madden
 * Changes for backward-compatability with old blast.
 *
 * Revision 1.15  1996/03/25  16:34:19  madden
 * Changes to mimic old statistics.
 *
 * Revision 1.14  1996/03/01  19:40:37  madden
 * Added factorial correction to SmallGapSum.
 *
 * Revision 1.13  1996/02/28  21:38:08  madden
 * *** empty log message ***
 *
 * Revision 1.12  1996/02/15  23:20:01  madden
 * Added query length to a number of calls.
 *
 * Revision 1.12  1996/02/15  23:20:01  madden
 * Added query length to a number of calls.
 *
 * Revision 1.11  1996/02/09  13:51:08  madden
 * Change to prevent log sign error.
 *
 * Revision 1.10  1996/01/22  22:33:06  madden
 * Changed effective length, fixed (improper) calculation of sum_e.
 *
 * Revision 1.9  1996/01/17  23:18:41  madden
 * Added BlastScoreSetAmbigRes.
 * ci -l blastkar.h
 *
 * Revision 1.8  1996/01/17  17:01:19  madden
 * Fixed BlastSmallGapSumE  and BlastLargeGapSumE.
 *
 * Revision 1.7  1996/01/17  13:46:36  madden
 * Added function BlastKarlinPtoE.
 *
 * Revision 1.6  1996/01/06  18:58:07  madden
 * Added BlastSumP.
 *
 * Revision 1.5  1995/12/28  21:26:16  madden
 * *** empty log message ***
 *
 * Revision 1.4  1995/12/26  23:04:46  madden
 * removed GetBLAST0Matrix, Added "cutoff" functions.
 *
 * Revision 1.3  1995/12/26  20:27:47  madden
 * Replaced BLAST_ScoreMat wiht BLAST_ScorePtr PNTR.
 *
 * Revision 1.2  1995/12/26  14:25:56  madden
 * *** empty log message ***
 *
 * Revision 1.1  1995/12/21  23:07:35  madden
 * Initial revision
 *
 * */
#include <ncbi.h>
#include <ncbimath.h>
#include <objcode.h>
#include <objseq.h>
#include <sequtil.h>
#include <blastkar.h>
#include <blastpri.h>

static Int2 BlastScoreFreqCalc PROTO((BLAST_ScoreBlkPtr sbp, BLAST_ScoreFreqPtr sfp, BLAST_ResFreqPtr rfp1, BLAST_ResFreqPtr rfp2));

/* OSF1 apparently doesn't like this. */
#if defined(HUGE_VAL) && !defined(OS_UNIX_OSF1)
#define BLASTKAR_HUGE_VAL HUGE_VAL
#else
#define BLASTKAR_HUGE_VAL 1.e30
#endif


/* Allocates and Deallocates the two-dimensional matrix. */
static BLASTMatrixStructurePtr BlastMatrixAllocate PROTO((Int2 alphabet_size));
static BLASTMatrixStructurePtr BlastMatrixDestruct PROTO((BLASTMatrixStructurePtr matrix_struct));

/* performs sump calculation, used by BlastSumPStd */
static Nlm_FloatHi BlastSumPCalc PROTO((int r, Nlm_FloatHi s));

#define COMMENT_CHR	'#'
#define TOKSTR	" \t\n\r"
#define BLAST_MAX_ALPHABET 40 /* ncbistdaa is only 26, this should be enough */

/*
	How many of the first bases are not ambiguous
	(it's four, of course).
*/
#define NUMBER_NON_AMBIG_BP 4
/*
	translates between ncbi4na and blastna.  the first four elements
	of this array match ncbi2na.
*/

Uint1 ncbi4na_to_blastna[] = {
15,/* Gap, 0 */
0, /* A,   1 */
1, /* C,   2 */
6, /* M,   3 */
2, /* G,   4 */
4, /* R,   5 */
9, /* S,   6 */
13, /* V,   7 */
3, /* T,   8 */
8, /* W,   9 */
 5, /* Y,  10 */
 12, /* H,  11 */
 7, /* K,  12 */
 11, /* D,  13 */
 10, /* B,  14 */
 14  /* N,  15 */
};

Uint1 blastna_to_ncbi4na[] = {	 1, /* A, 0 */
				 2, /* C, 1 */
				 4, /* G, 2 */
				 8, /* T, 3 */
				 5, /* R, 4 */
				10, /* Y, 5 */
				 3, /* M, 6 */
				12, /* K, 7 */
				 9, /* W, 8 */
				 6, /* S, 9 */
				14, /* B, 10 */
				13, /* D, 11 */
				11, /* H, 12 */
				 7, /* V, 13 */
				15, /* N, 14 */
				 0, /* Gap, 15 */
			};

/* Used in BlastKarlinBlkGappedCalc */
typedef FloatHi array_of_5[5];

/* Used to temporarily store matrix values for retrieval. */

typedef struct _matrix_info {
	CharPtr		name;			/* name of matrix (e.g., BLOSUM90). */
	array_of_5 	*values;		/* The values (below). */
	Int4		*prefs;			/* Preferences for display. */
	Int4		max_number_values;	/* number of values (e.g., BLOSUM90_VALUES_MAX). */
} MatrixInfo, PNTR MatrixInfoPtr;


/**************************************************************************************

		NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!

	The arrays below list all the matrices allowed for gapped BLAST.
	If new ones are added, remember to edit the function BlastLoadMatrixValues!!!

***************************************************************************************/
	
	
#define BLOSUM90_VALUES_MAX 7
static Nlm_FloatHi blosum90_values[BLOSUM90_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX,  (Nlm_FloatHi) INT2_MAX,	0.3346,     0.190,      0.75},     
	{  8.0,  2.0,	0.297,     0.088,      0.50},     
	{  7.0,  2.0,   0.285,     0.077,      0.42}, 
	{  6.0,  2.0,   0.269,     0.072,      0.31},
	{ 11.0,  1.0,   0.304,     0.098,      0.52},
	{ 10.0,  1.0,   0.289,     0.072,      0.42},
	{  9.0,  1.0,   0.263,     0.040,      0.31},
};

static Int4 blosum90_prefs[BLOSUM90_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL
};


#define BLOSUM80_VALUES_MAX 7
static Nlm_FloatHi blosum80_values[BLOSUM80_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX,  (Nlm_FloatHi) INT2_MAX,	0.3430,     0.177,      0.66},     
	{  8.0,  2.0,	0.308,     0.089,      0.46},     
	{  7.0,  2.0,   0.295,     0.077,      0.38},
	{  6.0,  2.0,   0.271,     0.051,      0.28},
	{ 11.0,  1.0,   0.314,     0.096,      0.48},
	{ 10.0,  1.0,   0.300,     0.072,      0.39},
	{  9.0,  1.0,   0.277,     0.046,      0.30}
};

static Int4 blosum80_prefs[BLOSUM80_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_BEST,
BLAST_MATRIX_PREFERRED
};


#define BLOSUM62_20_VALUES_MAX 8
static Nlm_FloatHi blosum62_20_values[BLOSUM62_20_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX, 0.03391, 0.125, 0.45},
	{100.0, 10.0, 0.0293, 0.054, 0.26},
	{97.0, 10.0, 0.0286, 0.046, 0.24},
	{95.0, 10.0, 0.0280, 0.041, 0.22},
	{93.0, 10.0, 0.0273, 0.035, 0.21},
	{90.0, 10.0, 0.0266, 0.031, 0.19},
	{88.0, 10.0, 0.0261, 0.029, 0.17},
	{87.0, 10.0, 0.0256, 0.026, 0.16}
}; 

static Int4 blosum62_20_prefs[BLOSUM62_20_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL
};

#define BLOSUM62_VALUES_MAX 7
static Nlm_FloatHi blosum62_values[BLOSUM62_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX,     0.3176,     0.134,    0.40},
	{  9.0,  2.0,     0.285,     0.075,    0.27},
	{  8.0,  2.0,     0.265,     0.046,    0.22},
	{  7.0,  2.0,     0.243,     0.032,    0.17},
	{ 12.0,  1.0,     0.281,     0.057,    0.27},
	{ 11.0,  1.0,     0.270,     0.047,    0.23},
	{ 10.0,  1.0,     0.250,     0.033,    0.17},
}; 
				
static Int4 blosum62_prefs[BLOSUM62_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_BEST,
BLAST_MATRIX_PREFERRED
};

#define BLOSUM50_VALUES_MAX 13
static Nlm_FloatHi blosum50_values[BLOSUM50_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX,     0.232,     0.11,      0.34},
	{12.0,  3.0,     0.206,     0.055,      0.23},
	{11.0,  3.0,     0.198,     0.046,      0.20},
	{10.0,  3.0,     0.189,     0.038,      0.17},
	{ 9.0,  3.0,     0.177,     0.030,      0.14},
	{15.0,  2.0,     0.211,     0.062,      0.25},
	{14.0,  2.0,     0.205,     0.053,      0.22},
	{13.0,  2.0,    0.197,     0.043,      0.19},
	{12.0,  2.0,     0.183,     0.028,      0.15},
	{18.0,  1.0,     0.208,     0.055,      0.23},
	{17.0,  1.0,     0.200,     0.042,      0.20},
	{16.0,  1.0,     0.189,     0.030,      0.17},
	{15.0,  1.0,     0.175,     0.020,      0.13},
};

static Int4 blosum50_prefs[BLOSUM50_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL
};

#define BLOSUM45_VALUES_MAX 13
static Nlm_FloatHi blosum45_values[BLOSUM45_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX,     0.2291,     0.092,      0.25},
	{13.0,  3.0,     0.209,     0.057,      0.19},
	{12.0,  3.0,     0.203,     0.049,      0.17},
	{11.0,  3.0,     0.193,     0.037,      0.15},
	{10.0,  3.0,     0.182,     0.029,      0.12},
	{15.0,  2.0,     0.206,     0.049,      0.18},
	{14.0,  2.0, 	 0.199,     0.040,      0.16},
	{13.0,  2.0,	 0.190,     0.032,      0.14},
	{12.0,  2.0,     0.177,     0.023,      0.11},
	{19.0,  1.0,     0.209,     0.049,      0.19},
	{18.0,  1.0,     0.202,     0.041,      0.17},
	{17.0,  1.0,     0.195,     0.034,      0.14},
	{16.0,  1.0,     0.183,     0.024,      0.12}
};

static Int4 blosum45_prefs[BLOSUM45_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_BEST,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED
};

#define PAM30_VALUES_MAX 10
static Nlm_FloatHi pam30_values[PAM30_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX,     0.340,     0.283,       1.75},
	{5.0,  3.0,     0.301,     0.14,       1.06},
	{4.0,  3.0,     0.286,     0.12,       0.85},
	{3.0,  3.0,     0.259,     0.081,      0.61},
	{7.0,  2.0,     0.306,     0.15,       1.07 },
	{6.0,  2.0,     0.292,     0.13,       0.86},
	{5.0,  2.0,     0.263,     0.077,      0.60},
	{10.0,  1.0,     0.309,     0.15,       1.07 },
	{9.0,  1.0,    0.295,     0.12,       0.85  },
	{8.0,  1.0,     0.270,     0.070,      0.59},
};


static Int4 pam30_prefs[PAM30_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_BEST,
BLAST_MATRIX_PREFERRED
};

#define PAM70_VALUES_MAX 10
static Nlm_FloatHi pam70_values[PAM70_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX,     0.3345,     0.229,      1.03},
	{6.0,  3.0,     0.297,     0.11,      0.67},
	{5.0,  3.0,     0.285,     0.10,       0.55},
	{4.0,  3.0,     0.258,     0.062,      0.40},
	{8.0,  2.0,     0.303,     0.13,       0.67}, 
	{7.0,  2.0,     0.287,     0.095,      0.56},
	{6.0,  2.0,     0.269,     0.079,      0.42},
	{11.0,  1.0,     0.307,     0.13,       0.70},
	{10.0,  1.0,    0.291,     0.089,      0.57},
	{9.0,  1.0,     0.269,     0.058,      0.42},
};

static Int4 pam70_prefs[PAM70_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_PREFERRED,
BLAST_MATRIX_BEST,
BLAST_MATRIX_PREFERRED
};


#define PAM250_VALUES_MAX 13
static Nlm_FloatHi pam250_values[PAM250_VALUES_MAX][5] = {
	{(Nlm_FloatHi) INT2_MAX, (Nlm_FloatHi) INT2_MAX,     0.229,     0.09,      0.23},
	{13.0,  3.0,     0.207,     0.051,      0.17},
	{12.0,  3.0,     0.200,     0.043,      0.15},
	{11.0,  3.0,     0.191,     0.034,      0.13},
	{10.0,  3.0,     0.181,     0.028,      0.11},
	{15.0,  2.0,     0.203,     0.043,      0.16},
	{14.0,  2.0,    0.196,     0.036,      0.14},
	{13.0,  2.0,     0.188,     0.030,      0.12},
	{12.0,  2.0,     0.175,     0.022,      0.10},
	{19.0,  1.0,     0.208,     0.049,      0.17 },
	{18.0,  1.0,     0.202,     0.040,      0.15},
	{17.0,  1.0,     0.194,     0.034,      0.13},
	{16.0,  1.0,     0.180,     0.021,      0.10}
};

static Int4 pam250_prefs[PAM250_VALUES_MAX] = {
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL,
BLAST_MATRIX_NOMINAL
};

/*
	Allocates memory for the BLAST_ScoreBlkPtr.
*/

BLAST_ScoreBlkPtr LIBCALL
BLAST_ScoreBlkNew(Uint1 alphabet, Int2 number_of_contexts)

{
	BLAST_ScoreBlkPtr sbp;
	SeqCodeTablePtr sctp;

	if (alphabet != BLASTNA_SEQ_CODE)
	{
		sctp = SeqCodeTableFindObj(alphabet);
		if (sctp == NULL)
			return NULL;
	}

	sbp = (BLAST_ScoreBlkPtr) MemNew(sizeof(BLAST_ScoreBlk));

	if (sbp != NULL)
	{
		sbp->alphabet_code = alphabet;
		if (alphabet != BLASTNA_SEQ_CODE)
			sbp->alphabet_size = sctp->num;
		else
			sbp->alphabet_size = BLASTNA_SIZE;

/* set the alphabet type (protein or not). */
		switch (alphabet) {
			case Seq_code_iupacaa:
			case Seq_code_ncbi8aa:
			case Seq_code_ncbieaa:
			case Seq_code_ncbipaa:
			case Seq_code_iupacaa3:
				return NULL;

			case Seq_code_ncbistdaa:
				sbp->protein_alphabet = TRUE;
				break;
			case BLASTNA_SEQ_CODE:
				sbp->protein_alphabet = FALSE;
				break;
			default:
				break;
		}

		sbp->matrix_struct = BlastMatrixAllocate(sbp->alphabet_size);
		if (sbp->matrix_struct == NULL)
		{
			sbp = BLAST_ScoreBlkDestruct(sbp);
			return sbp;
		}
		sbp->matrix = sbp->matrix_struct->matrix;
		sbp->maxscore = (BLAST_ScorePtr) MemNew(BLAST_MATRIX_SIZE*sizeof(BLAST_Score));
		sbp->number_of_contexts = number_of_contexts;
		sbp->sfp = MemNew(sbp->number_of_contexts*sizeof(BLAST_ScoreFreqPtr));
		sbp->kbp_std = MemNew(sbp->number_of_contexts*sizeof(BLAST_KarlinBlkPtr));
		sbp->kbp_gap_std = MemNew(sbp->number_of_contexts*sizeof(BLAST_KarlinBlkPtr));
		sbp->kbp_psi = MemNew(sbp->number_of_contexts*sizeof(BLAST_KarlinBlkPtr));
		sbp->kbp_gap_psi = MemNew(sbp->number_of_contexts*sizeof(BLAST_KarlinBlkPtr));
	}

	return sbp;
}

BLAST_ScoreBlkPtr LIBCALL
BLAST_ScoreBlkDestruct(BLAST_ScoreBlkPtr sbp)

{
	Int2 index;

	if (sbp == NULL)
		return NULL;

	for (index=0; index<sbp->number_of_contexts; index++)
	{
		if (sbp->sfp)
			sbp->sfp[index] = BlastScoreFreqDestruct(sbp->sfp[index]);
		if (sbp->kbp_std)
			sbp->kbp_std[index] = BlastKarlinBlkDestruct(sbp->kbp_std[index]);
		if (sbp->kbp_gap_std)
			sbp->kbp_gap_std[index] = BlastKarlinBlkDestruct(sbp->kbp_gap_std[index]);
		if (sbp->kbp_psi)
			sbp->kbp_psi[index] = BlastKarlinBlkDestruct(sbp->kbp_psi[index]);
		if (sbp->kbp_gap_psi)
			sbp->kbp_gap_psi[index] = BlastKarlinBlkDestruct(sbp->kbp_gap_psi[index]);
	}
	if (sbp->kbp_ideal)
		sbp->kbp_ideal = BlastKarlinBlkDestruct(sbp->kbp_ideal);
	sbp->sfp = MemFree(sbp->sfp);
	sbp->kbp_std = MemFree(sbp->kbp_std);
	sbp->kbp_psi = MemFree(sbp->kbp_psi);
	sbp->kbp_gap_std = MemFree(sbp->kbp_gap_std);
	sbp->kbp_gap_psi = MemFree(sbp->kbp_gap_psi);
	sbp->matrix_struct = BlastMatrixDestruct(sbp->matrix_struct);
	sbp->maxscore = MemFree(sbp->maxscore);
	sbp->comments = ValNodeFreeData(sbp->comments);
	sbp->name = MemFree(sbp->name);
	sbp->ambiguous_res = MemFree(sbp->ambiguous_res);

	sbp = MemFree(sbp);

	return sbp;
}

/* 
	Set the ambiguous residue (e.g, 'N', 'X') in the BLAST_ScoreBlkPtr.
	Convert from ncbieaa to sbp->alphabet_code (i.e., ncbistdaa) first.
*/
Int2 LIBCALL
BlastScoreSetAmbigRes(BLAST_ScoreBlkPtr sbp, Char ambiguous_res)

{
	Int2 index;
	SeqMapTablePtr smtp;
	Uint1Ptr ambig_buffer;

	if (sbp == NULL)
		return 1;

	if (sbp->ambig_occupy >= sbp->ambig_size)
	{
		sbp->ambig_size += 5;
		ambig_buffer = MemNew(sbp->ambig_size*sizeof(Uint1));
		for (index=0; index<sbp->ambig_occupy; index++)
		{
			ambig_buffer[index] = sbp->ambiguous_res[index];
		}
		sbp->ambiguous_res = MemFree(sbp->ambiguous_res);
		sbp->ambiguous_res = ambig_buffer;
	}

	if (sbp->alphabet_code == Seq_code_ncbistdaa)
	{
		smtp = SeqMapTableFind(sbp->alphabet_code, Seq_code_ncbieaa);
		sbp->ambiguous_res[sbp->ambig_occupy] = SeqMapTableConvert(smtp, ambiguous_res);
	}
	else if (sbp->alphabet_code == BLASTNA_SEQ_CODE)
	{
		smtp = SeqMapTableFind(Seq_code_ncbi4na, Seq_code_iupacna);
		sbp->ambiguous_res[sbp->ambig_occupy] = ncbi4na_to_blastna[SeqMapTableConvert(smtp, ambiguous_res)];
	}
	else if (sbp->alphabet_code == Seq_code_ncbi4na)
	{
		smtp = SeqMapTableFind(Seq_code_ncbi4na, Seq_code_iupacna);
		sbp->ambiguous_res[sbp->ambig_occupy] = SeqMapTableConvert(smtp, ambiguous_res);
	}

	(sbp->ambig_occupy)++;
	

	return 0;
}

/*
	Deallocates all data associated with the BLAST_MatrixPtr.
*/
BLAST_MatrixPtr LIBCALL
BLAST_MatrixDestruct(BLAST_MatrixPtr blast_matrix)

{
	Int4 index;
	if (blast_matrix == NULL)
		return NULL;

	blast_matrix->name = MemFree(blast_matrix->name);
	for (index=0; index<blast_matrix->rows; index++)
	{
		MemFree(blast_matrix->matrix[index]);
	}
	MemFree(blast_matrix->matrix);
	MemFree(blast_matrix);

	return NULL;
}


/*
	Allocates and fills the BLAST_MatrixPtr.  
	positionBased indicates that the posMatrix
	on BLAST_ScoreBlkPtr should be used rather
	than the 'normal' matrix.
*/
BLAST_MatrixPtr LIBCALL
BLAST_MatrixFill(BLAST_ScoreBlkPtr sbp, Boolean positionBased)

{
	BLAST_MatrixPtr blast_matrix;
	FloatHi karlinK = 0.0;
	Int4 index1, index2, dim1, dim2;
	Int4Ptr PNTR matrix, PNTR original_matrix;

	if (sbp == NULL)
		return NULL;

	blast_matrix = (BLAST_MatrixPtr) MemNew(sizeof(BLAST_Matrix));

	if (sbp->posMatrix)
	{
		dim1 = sbp->query_length + 1;
		dim2 = sbp->alphabet_size;
		original_matrix = sbp->posMatrix;
	}
	else
	{
		dim1 = sbp->alphabet_size;
		dim2 = sbp->alphabet_size;
		original_matrix = sbp->matrix;
	}
	if (sbp->kbp_gap_psi[0])
		karlinK = sbp->kbp_gap_psi[0]->K;

	matrix = MemNew(dim1*sizeof(Int4Ptr));
	for (index1=0; index1<dim1; index1++)
	{
		matrix[index1] = MemNew(dim2*sizeof(Int4));
		for (index2=0; index2<dim2; index2++)
		{
			matrix[index1][index2] = original_matrix[index1][index2];
		}
	}

	blast_matrix->is_prot = sbp->protein_alphabet;
	blast_matrix->name = StringSave(sbp->name);
	blast_matrix->rows = dim1;
	blast_matrix->columns = dim2;
	blast_matrix->matrix = matrix;
	blast_matrix->karlinK = karlinK;

	return blast_matrix;
}

/* 
	Fill in the matrix for blastn using the penaly and rewards on
	the BLAST_ScoreBlkPtr.

	The query sequence alphabet is blastna, the subject sequence
	is ncbi2na.  The alphabet blastna is defined in blastkar.h
	and the first four elements of blastna are identical to ncbi2na.

	The query is in the first index, the subject is the second.
*/

static Int2
BlastScoreBlkMatCreate(BLAST_ScoreBlkPtr sbp)
{

	BLAST_Score     penalty, reward;
	BLAST_ScorePtr PNTR	matrix;
	Char buffer[25];
	Int2	index1, index2, degen;
	Int2 degeneracy[BLASTNA_SIZE+1];

	matrix = sbp->matrix;
	
	for (index1 = 0; index1<BLASTNA_SIZE; index1++) /* blastna */
		for (index2 = 0; index2<BLASTNA_SIZE; index2++) /* blastna */
			matrix[index1][index2] = 0;

	/* In blastna the 1st four bases are A, C, G, and T, exactly as it is ncbi2na. */
	/* ncbi4na gives them the value 1, 2, 4, and 8.  */

	/* Set the first four bases to degen. one */
	for (index1=0; index1<NUMBER_NON_AMBIG_BP; index1++)
		degeneracy[index1] = 1;

	for (index1=NUMBER_NON_AMBIG_BP; index1<BLASTNA_SIZE; index1++) /* blastna */
	{
		degen=0;
		for (index2=0; index2<NUMBER_NON_AMBIG_BP; index2++) /* ncbi2na */
		{
			if (blastna_to_ncbi4na[index1] & blastna_to_ncbi4na[index2])
				degen++;
		}
		degeneracy[index1] = degen;
	}

	reward = sbp->reward;
	penalty = sbp->penalty;

	for (index1=0; index1<BLASTNA_SIZE; index1++) /* blastna */
	{
		for (index2=index1; index2<BLASTNA_SIZE; index2++) /* blastna */
		{
			if (blastna_to_ncbi4na[index1] & blastna_to_ncbi4na[index2])
			{ /* round up for positive scores, down for negatives. */
				matrix[index1][index2] = Nlm_Nint( (double) ((degeneracy[index2]-1)*penalty + reward))/degeneracy[index2];
				if (index1 != index2)
				{
				      matrix[index2][index1] = matrix[index1][index2];
				}
			}
			else
			{
				matrix[index1][index2] = penalty;
				matrix[index2][index1] = penalty;
			}
		}
	}

	sbp->mat_dim1 = BLASTNA_SIZE;
	sbp->mat_dim2 = BLASTNA_SIZE;

	sprintf(buffer, "blastn matrix:%ld %ld", (long) sbp->reward, (long) sbp->penalty);
	sbp->name = StringSave(buffer);

	return 0;
}



/*
	This function fills in the BLAST_ScoreBlk structure.  
	Tasks are:
		-read in the matrix
		-set maxscore
*/

Int2 LIBCALL
BlastScoreBlkMatFill(BLAST_ScoreBlkPtr sbp, CharPtr matrix)

{
	Char string[PATH_MAX], alphabet_type[3];
	CharPtr matrix_dir;
	Int2 status;
	FILE *fp;

	status=0;

	fp = NULL;
	if (sbp->read_in_matrix)
	{
		/* Convert matrix name to upper case. */
		matrix = Nlm_StrUpper(matrix);

		sbp->name = StringSave(matrix);	/* Save the name of the matrix. */

		fp = FileOpen(matrix, "r");

		if (fp == NULL)
		{
			if(FindPath("ncbi", "ncbi", "data", string, PATH_MAX))
			{
				StringCat(string, matrix);
				fp = FileOpen(string, "r");
			}
		}
			
#ifdef OS_UNIX
		/* Get the matrix locations from the environment for UNIX. */
		if (fp == NULL)
		{
			if (sbp->protein_alphabet)
				Nlm_StringNCpy(alphabet_type, "aa", 2);
			else
				Nlm_StringNCpy(alphabet_type, "nt", 2);
			alphabet_type[2] = NULLB;

			matrix_dir = getenv("BLASTMAT");
			if (matrix_dir != NULL)
			{
				sprintf(string, "%s%s%s%s%s", matrix_dir, DIRDELIMSTR,alphabet_type, DIRDELIMSTR, matrix); 
			}
			else
			{
				sprintf(string, "%s%s%s%s%s", BLASTMAT_DIR, DIRDELIMSTR, alphabet_type, DIRDELIMSTR, matrix); 
			}
			fp = FileOpen(string, "r");
			/* Try again without "aa" or "nt" */
			if (fp == NULL)
			{
				if (matrix_dir != NULL)
				{
					sprintf(string, "%s%s%s", matrix_dir, DIRDELIMSTR, matrix); 
				}
				else
				{
					sprintf(string, "%s%s%s", BLASTMAT_DIR, DIRDELIMSTR, matrix); 
				}
				fp = FileOpen(string, "r");
			}
		}
#endif
		if (fp == NULL)
		{
			ErrPostEx(SEV_WARNING, 0, 0, "Unable to open %s", matrix);
			return 4;
		}

		if((status=BlastScoreBlkMatRead(sbp, fp)) != 0)
		{
			FileClose(fp);
			return status;
		}
		FileClose(fp);
	}
	else
	{
		if((status=BlastScoreBlkMatCreate(sbp)) != 0)
			return status;
	}
	
	if((status=BlastScoreBlkMaxScoreSet(sbp)) != 0)
		return status;

	return status;
}

/*
	Calculate the Karlin parameters.  This function should be called once
	for each context, or frame translated.

	The rfp and stdrfp are calculated for each context, this should be
	fixed. 
*/

Int2 LIBCALL
BlastScoreBlkFill(BLAST_ScoreBlkPtr sbp, CharPtr query, Int4 query_length, Int2 context_number)
{
	BLAST_ResFreqPtr rfp, stdrfp;
	Int2 retval=0;

	rfp = BlastResFreqNew(sbp);
	stdrfp = BlastResFreqNew(sbp);
	BlastResFreqStdComp(sbp, stdrfp);
	BlastResFreqString(sbp, rfp, query, query_length);
	sbp->sfp[context_number] = BlastScoreFreqNew(sbp->loscore, sbp->hiscore);
	BlastScoreFreqCalc(sbp, sbp->sfp[context_number], rfp, stdrfp);
	sbp->kbp_std[context_number] = BlastKarlinBlkCreate();
	retval = BlastKarlinBlkCalc(sbp->kbp_std[context_number], sbp->sfp[context_number]);
	if (retval)
	{
		rfp = BlastResFreqDestruct(rfp);
		stdrfp = BlastResFreqDestruct(stdrfp);
		return retval;
	}
	sbp->kbp_psi[context_number] = BlastKarlinBlkCreate();
	retval = BlastKarlinBlkCalc(sbp->kbp_psi[context_number], sbp->sfp[context_number]);
	rfp = BlastResFreqDestruct(rfp);
	stdrfp = BlastResFreqDestruct(stdrfp);

	return retval;
}

/*
	Calculates the standard Karlin parameters.  This is used
	if the query is translated and the calculated (real) Karlin
	parameters are bad, as they're calculated for non-coding regions.
*/

BLAST_KarlinBlkPtr LIBCALL
BlastKarlinBlkStandardCalcEx(BLAST_ScoreBlkPtr sbp)

{
	BLAST_KarlinBlkPtr kbp_ideal;
	BLAST_ResFreqPtr stdrfp;
	BLAST_ScoreFreqPtr sfp;

	stdrfp = BlastResFreqNew(sbp);
	BlastResFreqStdComp(sbp, stdrfp);
	sfp = BlastScoreFreqNew(sbp->loscore, sbp->hiscore);
	BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
	kbp_ideal = BlastKarlinBlkCreate();
	BlastKarlinBlkCalc(kbp_ideal, sfp);
	stdrfp = BlastResFreqDestruct(stdrfp);

	sfp = BlastScoreFreqDestruct(sfp);

	return kbp_ideal;
}

Int2 LIBCALL
BlastKarlinBlkStandardCalc(BLAST_ScoreBlkPtr sbp, Int2 context_start, Int2 context_end)
{

	BLAST_KarlinBlkPtr kbp_ideal, kbp;
	Int2 index;

	kbp_ideal = BlastKarlinBlkStandardCalcEx(sbp);
/* Replace the calculated values with ideal ones for blastx, tblastx. */
	for (index=context_start; index<=context_end; index++)
	{
		kbp = sbp->kbp[index];	
		if (kbp->Lambda >= kbp_ideal->Lambda)
		{
			kbp->Lambda = kbp_ideal->Lambda;
			kbp->K = kbp_ideal->K;
			kbp->logK = kbp_ideal->logK;
			kbp->H = kbp_ideal->H;
		}
	}
	kbp_ideal = BlastKarlinBlkDestruct(kbp_ideal);

	return 0;
}

/*
	Creates the Karlin Block.
*/

BLAST_KarlinBlkPtr LIBCALL
BlastKarlinBlkCreate(void)

{
	BLAST_KarlinBlkPtr kbp;

	kbp = (BLAST_KarlinBlkPtr) MemNew(sizeof(BLAST_KarlinBlk));

	return kbp;
}

/*
	Deallocates the Karlin Block.
*/

BLAST_KarlinBlkPtr LIBCALL
BlastKarlinBlkDestruct(BLAST_KarlinBlkPtr kbp)

{
	kbp = MemFree(kbp);

	return kbp;
}


/* 
	Read in the matrix from the FILE *fp.

	This function ASSUMES that the matrices are in the ncbistdaa
	format.  BLAST should be able to use any alphabet, though it
	is expected that ncbistdaa will be used.
*/

static Char ASCII_TO_BLASTNA_CONVERT[128]={
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};

Int2 LIBCALL
BlastScoreBlkMatRead(BLAST_ScoreBlkPtr sbp, FILE *fp)
{
	Char	buf[512+3];
	Char	temp[512];
	CharPtr	cp, lp;
	Char		ch;
	BLAST_ScorePtr PNTR	matrix;
	BLAST_ScorePtr	m;
	BLAST_Score	score;
	Int2		a1cnt = 0, a2cnt = 0;
	Char    a1chars[BLAST_MAX_ALPHABET], a2chars[BLAST_MAX_ALPHABET];
	long	lineno = 0;
	Nlm_FloatHi	xscore;
	register int	index1, index2, total;
	SeqCodeTablePtr sctp;
	SeqMapTablePtr smtp=NULL;
	Int2 status;

	matrix = sbp->matrix;	

	if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
	  sctp = SeqCodeTableFindObj(sbp->alphabet_code);
	  if(sctp==NULL) {
	    return 1;
	  }
	  total = sctp->start_at + sctp->num;
	  for (index1 = sctp->start_at; index1 < total; index1++)
	        for (index2 = sctp->start_at; index2 < total; index2++)
			matrix[index1][index2] = BLAST_SCORE_MIN;

	  if (sbp->alphabet_code != Seq_code_ncbieaa)
	    {
	      smtp = SeqMapTableFind(sbp->alphabet_code, Seq_code_ncbieaa);
	      if (smtp == NULL)
		{
			return 1;
		}
	    }
	} else {
	  /* Fill-in all the defaults ambiguity and normal codes */
	  status=BlastScoreBlkMatCreate(sbp); 
	  if(status != 0);
	    return status;
	}

	/* Read the residue names for the second alphabet */
	while (Nlm_FileGets(buf, sizeof(buf), fp) != NULL) {
		++lineno;
		if (Nlm_StrChr(buf, '\n') == NULL)
			return 2;
		if (buf[0] == COMMENT_CHR) {
			/* save the comment line in a linked list */
			*Nlm_StrChr(buf, '\n') = NULLB;
			ValNodeCopyStr(&sbp->comments, 0, buf+1);
			continue;
		}
		if ((cp = Nlm_StrChr(buf, COMMENT_CHR)) != NULL)
			*cp = NULLB;
		lp = (CharPtr)Nlm_StrTok(buf, TOKSTR);
		if (lp == NULL) /* skip blank lines */
			continue;
		while (lp != NULL)
		{
		
			if (smtp)
				ch = SeqMapTableConvert(smtp,  *lp);
			else {
			  if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
				ch = *lp;
                          } else {
				ch = ASCII_TO_BLASTNA_CONVERT[toupper(*lp)];
			  }
			}
			a2chars[a2cnt++] = ch;
			lp = (CharPtr)Nlm_StrTok(NULL, TOKSTR);
		}

		break;	/* Exit loop after reading one line. */
	}

	if (a2cnt <= 1)
	{
		return 2;
	}
	if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
	  sbp->mat_dim2 = a2cnt;
	}
	while (Nlm_FileGets(buf, sizeof(buf), fp) != NULL) 
	{
		++lineno;
		if ((cp = Nlm_StrChr(buf, '\n')) == NULL)
			return 2;
		if ((cp = Nlm_StrChr(buf, COMMENT_CHR)) != NULL)
			*cp = NULLB;
		if ((lp = (CharPtr)Nlm_StrTok(buf, TOKSTR)) == NULL)
			continue;
		ch = *lp;
		cp = (CharPtr) lp;
		if ((cp = Nlm_StrTok(NULL, TOKSTR)) == NULL)
			return 2;
		if (a1cnt >= DIM(a1chars))
			return 2;

		if (smtp) {
			ch = SeqMapTableConvert(smtp,  ch);

		} else {
			  if (sbp->alphabet_code == BLASTNA_SEQ_CODE) {
				ch = ASCII_TO_BLASTNA_CONVERT[toupper(ch)];
			  }
		}
		a1chars[a1cnt++] = ch;
		m = &matrix[ch][0];
		index2 = 0;
		while (cp != NULL)
		{
			if (index2 >= a2cnt)
				return 2;
			Nlm_StrCpy(temp, cp);

			if (Nlm_StrICmp(temp, "na") == 0) 
			{
				score = BLAST_SCORE_1MIN;
			}
			else 
			{
				if (sscanf(temp, "%lg", &xscore) != 1)
					return 2;
				/*xscore = MAX(xscore, BLAST_SCORE_1MIN);*/
				if (xscore > BLAST_SCORE_1MAX || xscore < BLAST_SCORE_1MIN)
					return 2;
				xscore += (xscore >= 0. ? 0.5 : -0.5);
				score = (BLAST_Score)xscore;
			}

			m[a2chars[index2++]] = score;

			cp = Nlm_StrTok(NULL, TOKSTR);
		}
	}

	if (a1cnt <= 1)
		return 2;
	if (sbp->alphabet_code != BLASTNA_SEQ_CODE) {
	        sbp->mat_dim1 = a1cnt;
	}
	return 0;
}

Int2 LIBCALL
BlastScoreBlkMaxScoreSet(BLAST_ScoreBlkPtr sbp)
{
	BLAST_Score score, maxscore;
	BLAST_ScorePtr PNTR  matrix; 
	Int2 index1, index2;

	sbp->loscore = BLAST_SCORE_1MAX;
        sbp->hiscore = BLAST_SCORE_1MIN;
	matrix = sbp->matrix;
	for (index1=0; index1<sbp->mat_dim1; index1++)
	{
		maxscore=BLAST_SCORE_MIN;
		for (index2=0; index2<sbp->mat_dim2; index2++)
		{
			score = matrix[index1][index2];
			if (score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX)
				continue;
			if (score > maxscore)
			{
				maxscore = score;
			}
			if (sbp->loscore > score)
				sbp->loscore = score;
			if (sbp->hiscore < score)
				sbp->hiscore = score;
		}
		sbp->maxscore[index1] = maxscore;
	}
/* If the lo/hi-scores are BLAST_SCORE_MIN/BLAST_SCORE_MAX, (i.e., for
gaps), then use other scores. */

	if (sbp->loscore < BLAST_SCORE_1MIN)
		sbp->loscore = BLAST_SCORE_1MIN;
	if (sbp->hiscore > BLAST_SCORE_1MAX)
		sbp->hiscore = BLAST_SCORE_1MAX;

	return 0;
}

static BLASTMatrixStructurePtr
BlastMatrixAllocate(Int2 alphabet_size)

{
	BLASTMatrixStructurePtr matrix_struct;
	Int2 index;

	if (alphabet_size <= 0 || alphabet_size >= BLAST_MATRIX_SIZE)
		return NULL;

	matrix_struct =	(BLASTMatrixStructurePtr) MemNew(sizeof(BLASTMatrixStructure));

	if (matrix_struct == NULL)
		return NULL;

	for (index=0; index<BLAST_MATRIX_SIZE-1; index++)
	{
		matrix_struct->matrix[index] = matrix_struct->long_matrix + index*BLAST_MATRIX_SIZE;
	}

	return matrix_struct;
}

static BLASTMatrixStructurePtr
BlastMatrixDestruct(BLASTMatrixStructurePtr matrix_struct)

{

	if (matrix_struct == NULL)
		return NULL;

	matrix_struct = MemFree(matrix_struct);

	return matrix_struct;
}

/* 
	Allocated the BLAST_ResCompPtr for a given alphabet.  Only the
	alphabets ncbistdaa and ncbi4na should be used by BLAST.
*/

BLAST_ResCompPtr LIBCALL
BlastResCompNew(BLAST_ScoreBlkPtr sbp)
{
	BLAST_ResCompPtr	rcp;

	rcp = (BLAST_ResCompPtr) MemNew(sizeof(BLAST_ResComp));
	if (rcp == NULL)
		return NULL;

	rcp->alphabet_code = sbp->alphabet_code;

/* comp0 has zero offset, comp starts at sctp->start_at, only one 
array is allocated.  */
	rcp->comp0 = (Int4Ptr) MemNew(BLAST_MATRIX_SIZE*sizeof(Int4));
	if (rcp->comp0 == NULL) 
	{
		rcp = BlastResCompDestruct(rcp);
		return rcp;
	}

	rcp->comp = rcp->comp0 - sbp->alphabet_start;

	return rcp;
}


BLAST_ResCompPtr LIBCALL
BlastResCompDestruct(BLAST_ResCompPtr rcp)
{
	if (rcp == NULL)
		return NULL;

	if (rcp->comp0 != NULL)
		rcp->comp0 = MemFree(rcp->comp0);

	return MemFree(rcp);
}

/*
	Store the composition of a (query) string.  
*/
Int2 LIBCALL
BlastResCompStr(BLAST_ScoreBlkPtr sbp, BLAST_ResCompPtr rcp, CharPtr str, Int4 length)
{
	CharPtr	lp, lpmax;
	Int2 index;

	if (sbp == NULL || rcp == NULL || str == NULL)
		return 1;

	if (rcp->alphabet_code != sbp->alphabet_code)
		return 1;

/* comp0 starts at zero and extends for "num", comp is the same array, but 
"start_at" offset. */
	for (index=0; index<(sbp->alphabet_size); index++)
		rcp->comp0[index] = 0;

	for (lp = str, lpmax = lp+length; lp < lpmax; lp++)
	{
		++rcp->comp[*lp];
	}

	/* Don't count ambig. residues. */
	for (index=0; index<sbp->ambig_occupy; index++)
	{
		rcp->comp[sbp->ambiguous_res[index]] = 0;
	}

	return 0;
}

static Int2
BlastScoreChk(BLAST_Score lo, BLAST_Score hi)
{
	if (lo >= 0 || hi <= 0 ||
			lo < BLAST_SCORE_1MIN || hi > BLAST_SCORE_1MAX)
		return 1;

	if (hi - lo > BLAST_SCORE_RANGE_MAX)
		return 1;

	return 0;
}

BLAST_ScoreFreqPtr
BlastScoreFreqNew(BLAST_Score score_min, BLAST_Score score_max)
{
	BLAST_ScoreFreqPtr	sfp;
	BLAST_Score	range;

	if (BlastScoreChk(score_min, score_max) != 0)
		return NULL;

	sfp = (BLAST_ScoreFreqPtr) MemNew(sizeof(BLAST_ScoreFreq));
	if (sfp == NULL)
		return NULL;

	range = score_max - score_min + 1;
	sfp->sprob = (Nlm_FloatHi PNTR) MemNew(sizeof(Nlm_FloatHi) * range);
	if (sfp->sprob == NULL) 
	{
		BlastScoreFreqDestruct(sfp);
		return NULL;
	}

	sfp->sprob0 = sfp->sprob;
	sfp->sprob -= score_min;
	sfp->score_min = score_min;
	sfp->score_max = score_max;
	sfp->obs_min = sfp->obs_max = 0;
	sfp->score_avg = 0.0;
	return sfp;
}

BLAST_ScoreFreqPtr
BlastScoreFreqDestruct(BLAST_ScoreFreqPtr sfp)
{
	if (sfp == NULL)
		return NULL;

	if (sfp->sprob0 != NULL)
		sfp->sprob0 = MemFree(sfp->sprob0);
	sfp = MemFree(sfp);
	return sfp;
}

static Int2
BlastScoreFreqCalc(BLAST_ScoreBlkPtr sbp, BLAST_ScoreFreqPtr sfp, BLAST_ResFreqPtr rfp1, BLAST_ResFreqPtr rfp2)
{
	BLAST_ScorePtr PNTR	matrix;
	BLAST_Score	score, obs_min, obs_max;
	Nlm_FloatHi		score_sum, score_avg;
	Int2 		alphabet_start, alphabet_end, index1, index2;

	if (sbp == NULL || sfp == NULL)
		return 1;

	if (sbp->loscore < sfp->score_min || sbp->hiscore > sfp->score_max)
		return 1;

	for (score = sfp->score_min; score <= sfp->score_max; score++)
		sfp->sprob[score] = 0.0;

	matrix = sbp->matrix;

	alphabet_start = sbp->alphabet_start;
	alphabet_end = alphabet_start + sbp->alphabet_size;
	for (index1=alphabet_start; index1<alphabet_end; index1++)
	{
		for (index2=alphabet_start; index2<alphabet_end; index2++)
		{
			score = matrix[index1][index2];
			if (score >= sbp->loscore) 
			{
				sfp->sprob[score] += rfp1->prob[index1] * rfp2->prob[index2];
			}
		}
	}

	score_sum = 0.;
	obs_min = obs_max = BLAST_SCORE_MIN;
	for (score = sfp->score_min; score <= sfp->score_max; score++)
	{
		if (sfp->sprob[score] > 0.) 
		{
			score_sum += sfp->sprob[score];
			obs_max = score;
			if (obs_min == BLAST_SCORE_MIN)
				obs_min = score;
		}
	}
	sfp->obs_min = obs_min;
	sfp->obs_max = obs_max;

	score_avg = 0.0;
	if (score_sum > 0.0001 || score_sum < -0.0001)
	{
		for (score = obs_min; score <= obs_max; score++) 
		{
			sfp->sprob[score] /= score_sum;
			score_avg += score * sfp->sprob[score];
		}
	}
	sfp->score_avg = score_avg;

	return 0;
}

typedef struct {
		Char	ch;
		Nlm_FloatHi	p;
	} BLAST_LetterProb;

#define STD_AMINO_ACID_FREQS Robinson_prob

#if STD_AMINO_ACID_FREQS == Dayhoff_prob
/*  M. O. Dayhoff amino acid background frequencies   */
static BLAST_LetterProb	Dayhoff_prob[] = {
		{ 'A', 87.13 },
		{ 'C', 33.47 },
		{ 'D', 46.87 },
		{ 'E', 49.53 },
		{ 'F', 39.77 },
		{ 'G', 88.61 },
		{ 'H', 33.62 },
		{ 'I', 36.89 },
		{ 'K', 80.48 },
		{ 'L', 85.36 },
		{ 'M', 14.75 },
		{ 'N', 40.43 },
		{ 'P', 50.68 },
		{ 'Q', 38.26 },
		{ 'R', 40.90 },
		{ 'S', 69.58 },
		{ 'T', 58.54 },
		{ 'V', 64.72 },
		{ 'W', 10.49 },
		{ 'Y', 29.92 }
	};
#endif

#if STD_AMINO_ACID_FREQS == Altschul_prob
/* Stephen Altschul amino acid background frequencies */
static BLAST_LetterProb Altschul_prob[] = {
		{ 'A', 81.00 },
		{ 'C', 15.00 },
		{ 'D', 54.00 },
		{ 'E', 61.00 },
		{ 'F', 40.00 },
		{ 'G', 68.00 },
		{ 'H', 22.00 },
		{ 'I', 57.00 },
		{ 'K', 56.00 },
		{ 'L', 93.00 },
		{ 'M', 25.00 },
		{ 'N', 45.00 },
		{ 'P', 49.00 },
		{ 'Q', 39.00 },
		{ 'R', 57.00 },
		{ 'S', 68.00 },
		{ 'T', 58.00 },
		{ 'V', 67.00 },
		{ 'W', 13.00 },
		{ 'Y', 32.00 }
	};
#endif

#if STD_AMINO_ACID_FREQS == Robinson_prob
/* amino acid background frequencies from Robinson and Robinson */
static BLAST_LetterProb Robinson_prob[] = {
		{ 'A', 78.05 },
		{ 'C', 19.25 },
		{ 'D', 53.64 },
		{ 'E', 62.95 },
		{ 'F', 38.56 },
		{ 'G', 73.77 },
		{ 'H', 21.99 },
		{ 'I', 51.42 },
		{ 'K', 57.44 },
		{ 'L', 90.19 },
		{ 'M', 22.43 },
		{ 'N', 44.87 },
		{ 'P', 52.03 },
		{ 'Q', 42.64 },
		{ 'R', 51.29 },
		{ 'S', 71.20 },
		{ 'T', 58.41 },
		{ 'V', 64.41 },
		{ 'W', 13.30 },
		{ 'Y', 32.16 }
	};
#endif

static BLAST_LetterProb	nt_prob[] = {
		{ 'A', 25.00 },
		{ 'C', 25.00 },
		{ 'G', 25.00 },
		{ 'T', 25.00 }
	};

/*
	Allocates the BLAST_ResFreqPtr and then fills in the frequencies
	in the probabilities.
*/ 

BLAST_ResFreqPtr LIBCALL
BlastResFreqNew(BLAST_ScoreBlkPtr sbp)
{
	BLAST_ResFreqPtr	rfp;

	if (sbp == NULL)
	{
		return NULL;
	}

	rfp = (BLAST_ResFreqPtr) MemNew(sizeof(BLAST_ResFreq));
	if (rfp == NULL)
		return NULL;

	rfp->alphabet_code = sbp->alphabet_code;

	rfp->prob0 = (Nlm_FloatHi PNTR) MemNew(sizeof(Nlm_FloatHi) * sbp->alphabet_size);
	if (rfp->prob0 == NULL) 
	{
		rfp = BlastResFreqDestruct(rfp);
		return rfp;
	}
	rfp->prob = rfp->prob0 - sbp->alphabet_start;

	return rfp;
}

/*
	Normalize the frequencies to "norm".
*/
Int2 LIBCALL
BlastResFreqNormalize(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp, Nlm_FloatHi norm)
{
	Int2	alphabet_stop, index;
	Nlm_FloatHi	sum = 0., p;

	if (norm == 0.)
		return 1;

	alphabet_stop = sbp->alphabet_start + sbp->alphabet_size;
	for (index=sbp->alphabet_start; index<alphabet_stop; index++)
	{
		p = rfp->prob[index];
		if (p < 0.)
			return 1;
		sum += p;
	}
	if (sum <= 0.)
		return 0;

	for (index=sbp->alphabet_start; index<alphabet_stop; index++)
	{
		rfp->prob[index] /= sum;
		rfp->prob[index] *= norm;
	}
	return 0;
}

/*
	Fills a buffer with the 'standard' alphabet (given by 
	STD_AMINO_ACID_FREQS[index].ch).

	Return value is the number of residues in alphabet.  
	Negative returns upon error.
*/

Int2 LIBCALL
BlastGetStdAlphabet (Uint1 alphabet_code, Uint1Ptr residues, Int4 residues_size)

{
	Int2 index;
	SeqMapTablePtr smtp=NULL;

	if (residues_size < DIM(STD_AMINO_ACID_FREQS))
		return -2;

	if (alphabet_code != Seq_code_ncbieaa)
	{
		smtp = SeqMapTableFind(alphabet_code, Seq_code_ncbieaa);
		if (smtp == NULL)
		{
			return -1;
		}
	}

	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) 
	{
		if (smtp)
		{
	 		residues[index] = SeqMapTableConvert(smtp,  STD_AMINO_ACID_FREQS[index].ch);
		}
		else
		{
			residues[index] = STD_AMINO_ACID_FREQS[index].ch;
		}
	}

	return index;
}

Int2 LIBCALL
BlastResFreqStdComp(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp)
{
	Int2 index, retval;
	Uint1Ptr residues;

	if (sbp->protein_alphabet == TRUE)
	{
		residues = (Uint1Ptr) MemNew(DIM(STD_AMINO_ACID_FREQS)*sizeof(Uint1));
		retval = BlastGetStdAlphabet(sbp->alphabet_code, residues, DIM(STD_AMINO_ACID_FREQS));
		if (retval < 1)
			return retval;

		for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) 
		{
			rfp->prob[residues[index]] = STD_AMINO_ACID_FREQS[index].p;
		}
		residues = MemFree(residues);
	}
	else
	{	/* beginning of blastna and ncbi2na are the same. */
		/* Only blastna used  for nucleotides. */
		for (index=0; index<DIM(nt_prob); index++) 
		{
			rfp->prob[index] = nt_prob[index].p;
		}
	}

	BlastResFreqNormalize(sbp, rfp, 1.0);

	return 0;
}

CharPtr  LIBCALL
BlastRepresentativeResidues(Int2 length)

{
	CharPtr buffer, ptr;
	Int2 index, total;
	Int4 number;
	total=0;

	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) 
	{
		total += (Int2) STD_AMINO_ACID_FREQS[index].p;
	}

	buffer = (CharPtr) MemNew((length+1)*sizeof(Char));

	ptr = buffer;
	for (index=0; index<DIM(STD_AMINO_ACID_FREQS); index++) 
	{
		number = Nint((STD_AMINO_ACID_FREQS[index].p)*((Nlm_FloatHi) length)/((Nlm_FloatHi) total));
		while (number > 0)
		{
			*ptr = STD_AMINO_ACID_FREQS[index].ch;
			ptr++;
			number--;
		}
	}

	return buffer;
}

Int2 LIBCALL
BlastResFreqString(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp, CharPtr string, Int4 length)
{
	BLAST_ResCompPtr rcp;
	
	rcp = BlastResCompNew(sbp);

	BlastResCompStr(sbp, rcp, string, length);

	BlastResFreqResComp(sbp, rfp, rcp);

	rcp = BlastResCompDestruct(rcp);

	return 0;
}

BLAST_ResFreqPtr LIBCALL
BlastResFreqDestruct(BLAST_ResFreqPtr rfp)
{
	if (rfp == NULL)
		return NULL;

	if (rfp->prob0 != NULL)
		MemFree(rfp->prob0);

	rfp = MemFree(rfp);

	return rfp;
}

/*
	Calculate the residue frequencies associated with the provided ResComp
*/
Int2 LIBCALL
BlastResFreqResComp(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp, BLAST_ResCompPtr rcp)
{
	Int2	alphabet_max, index;
	Nlm_FloatHi	sum = 0.;

	if (rfp == NULL || rcp == NULL)
		return 1;

	if (rfp->alphabet_code != rcp->alphabet_code)
		return 1;

	alphabet_max = sbp->alphabet_start + sbp->alphabet_size;
	for (index=sbp->alphabet_start; index<alphabet_max; index++)
		sum += rcp->comp[index];

	if (sum == 0.) {
		BlastResFreqClr(sbp, rfp);
		return 0;
	}

	for (index=sbp->alphabet_start; index<alphabet_max; index++)
		rfp->prob[index] = rcp->comp[index] / sum;

	return 0;
}

Int2 LIBCALL
BlastResFreqClr(BLAST_ScoreBlkPtr sbp, BLAST_ResFreqPtr rfp)
{
	Int2	alphabet_max, index;
 
	alphabet_max = sbp->alphabet_start + sbp->alphabet_size;
	for (index=sbp->alphabet_start; index<alphabet_max; index++)
                rfp->prob[index] = 0.0;
 
        return 0;
}


/*
	Deallocates MatrixInfoPtr
*/

static MatrixInfoPtr
MatrixInfoDestruct(MatrixInfoPtr matrix_info)

{
	if (matrix_info == NULL)
		return NULL;

	MemFree(matrix_info->name);
	return MemFree(matrix_info);
}

/*
	Makes New MatrixInfoPtr
*/

static MatrixInfoPtr
MatrixInfoNew(CharPtr name, array_of_5 *values, Int4Ptr prefs, Int4 max_number)

{
	MatrixInfoPtr matrix_info;

	matrix_info = (MatrixInfoPtr) MemNew(sizeof(MatrixInfo));
	matrix_info->name = StringSave(name);
	matrix_info->values = values;
	matrix_info->prefs = prefs;
	matrix_info->max_number_values = max_number;

	return matrix_info;
}

static ValNodePtr
BlastMatrixValuesDestruct(ValNodePtr vnp)

{
	ValNodePtr head;

	head = vnp;
	while (vnp)
	{
		MatrixInfoDestruct((MatrixInfoPtr) vnp->data.ptrvalue);
		vnp = vnp->next;
	}

	head = ValNodeFree(head);

	return head;
}

/*
	ValNodePtr BlastLoadMatrixValues (void)
	
	Loads all the matrix values, returns a ValNodePtr chain that contains
	MatrixInfoPtr's.

*/
static ValNodePtr 
BlastLoadMatrixValues (void)

{
	MatrixInfoPtr matrix_info;
	ValNodePtr retval=NULL;

	matrix_info = MatrixInfoNew("BLOSUM80", blosum80_values, blosum80_prefs, BLOSUM80_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);
	matrix_info = MatrixInfoNew("BLOSUM62", blosum62_values, blosum62_prefs, BLOSUM62_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);
	matrix_info = MatrixInfoNew("BLOSUM50", blosum50_values, blosum50_prefs, BLOSUM50_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);
	matrix_info = MatrixInfoNew("BLOSUM45", blosum45_values, blosum45_prefs, BLOSUM45_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);
	matrix_info = MatrixInfoNew("PAM250", pam250_values, pam250_prefs, PAM250_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);
	matrix_info = MatrixInfoNew("BLOSUM62_20", blosum62_20_values, blosum62_20_prefs, BLOSUM62_20_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);
	matrix_info = MatrixInfoNew("BLOSUM90", blosum90_values, blosum90_prefs, BLOSUM90_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);
	matrix_info = MatrixInfoNew("PAM30", pam30_values, pam30_prefs, PAM30_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);
	matrix_info = MatrixInfoNew("PAM70", pam70_values, pam70_prefs, PAM70_VALUES_MAX);
	ValNodeAddPointer(&retval, 0, matrix_info);

	return retval;
}

/*
Int2 LIBCALL
BlastKarlinGetMatrixValues(CharPtr matrix, Int4Ptr open, Int4Ptr extension, FloatHiPtr lambda, FloatHiPtr K, FloatHiPtr H)
	
Obtains arrays of the allowed opening and extension penalties for gapped BLAST for
the given matrix.  Also obtains arrays of Lambda, K, and H.  Any of these fields that
are not required should be set to NULL.  The Int2 return value is the length of the
arrays.
*/

Int2 LIBCALL
BlastKarlinGetMatrixValues(CharPtr matrix, Int4Ptr PNTR open, Int4Ptr PNTR extension, FloatHiPtr PNTR lambda, FloatHiPtr PNTR K, FloatHiPtr PNTR H, Int4Ptr PNTR pref_flags)

{
	array_of_5 *values;
	Boolean found_matrix=FALSE;
	Int4 index, max_number_values=0;
	Int4Ptr open_array=NULL, extension_array=NULL, pref_flags_array=NULL, prefs;
	Nlm_FloatHiPtr lambda_array=NULL, K_array=NULL, H_array=NULL;
	MatrixInfoPtr matrix_info;
	ValNodePtr vnp, head;

	if (matrix == NULL)
		return 0;

	vnp = head = BlastLoadMatrixValues();

	while (vnp)
	{
		matrix_info = vnp->data.ptrvalue;
		if (StringICmp(matrix_info->name, matrix) == 0)
		{
			values = matrix_info->values;
			max_number_values = matrix_info->max_number_values;
			prefs = matrix_info->prefs;
			found_matrix = TRUE;
			break;
		}
		vnp = vnp->next;
	}

	if (found_matrix)
	{
		if (open)
			*open = open_array = MemNew(max_number_values*sizeof(Int4));
		if (extension)
			*extension = extension_array = MemNew(max_number_values*sizeof(Int4));
		if (lambda)
			*lambda = lambda_array = (FloatHiPtr) MemNew(max_number_values*sizeof(FloatHi));
		if (K)
			*K = K_array = (FloatHiPtr) MemNew(max_number_values*sizeof(FloatHi));
		if (H)
			*H = H_array = (FloatHiPtr) MemNew(max_number_values*sizeof(FloatHi));
		if (pref_flags)
			*pref_flags = pref_flags_array = MemNew(max_number_values*sizeof(Int4));

		for (index=0; index<max_number_values; index++)
		{
			if (open)
				open_array[index] = values[index][0];
			if (extension)
				extension_array[index] = values[index][1];
			if (lambda)
				lambda_array[index] = values[index][2];
			if (K)
				K_array[index] = values[index][3];
			if (H)
				H_array[index] = values[index][4];
			if (pref_flags)
				pref_flags_array[index] = prefs[index];
		}
	}

	BlastMatrixValuesDestruct(head);

	return max_number_values;
}
	
/*
	Supplies lambda, H, and K values, as calcualted by Stephen Altschul 
	in Methods in Enzy. (vol 266, page 474).

	if kbp is NULL, then a validation is perfomed.
*/

Int2 LIBCALL
BlastKarlinBlkGappedCalc(BLAST_KarlinBlkPtr kbp, Int4 gap_open, Int4 gap_extend, CharPtr matrix_name, ValNodePtr PNTR error_return)

{
	Boolean found_matrix, found_values;
	array_of_5 *values;
	Char buffer[256];
	Int4 index, max_number_values=0;
	MatrixInfoPtr matrix_info;
	ValNodePtr vnp, head;
	
	if (matrix_name == NULL)
		return 1;

	values = NULL;
	found_matrix = FALSE;
	found_values = FALSE;

	vnp = head = BlastLoadMatrixValues();
	while (vnp)
	{
		matrix_info = vnp->data.ptrvalue;
		if (StringICmp(matrix_info->name, matrix_name) == 0)
		{
			values = matrix_info->values;
			max_number_values = matrix_info->max_number_values;
			found_matrix = TRUE;
			break;
		}
		vnp = vnp->next;
	}


	if (found_matrix)
	{
		for (index=0; index<max_number_values; index++)
		{
			if (Nint(values[index][0]) == gap_open &&
				Nint(values[index][1]) == gap_extend)
			{
				if (kbp)
				{
					kbp->Lambda_real = kbp->Lambda = values[index][2];
					kbp->K_real = kbp->K = values[index][3];
					kbp->logK_real = kbp->logK = log(kbp->K);
					kbp->H_real = kbp->H = values[index][4];
				}
				found_values = TRUE;
				break;
			}
		}

		if (found_values == TRUE)
		{
			BlastMatrixValuesDestruct(head);
			return 0;
		}
	}
	else
	{
		sprintf(buffer, "%s is not a supported matrix", matrix_name);
		BlastConstructErrorMessage("BlastKarlinBlkGappedCalc", buffer, 2, error_return);
		vnp = head;
		while (vnp)
		{
			matrix_info = vnp->data.ptrvalue;
			sprintf(buffer, "%s is a supported matrix", matrix_info->name);
			BlastConstructErrorMessage("BlastKarlinBlkGappedCalc", buffer, 2, error_return);
			vnp = vnp->next;
		}
		BlastMatrixValuesDestruct(head);
		return 2;
	}

	sprintf(buffer, "Gap existence and extension values of %ld and %ld not supported for %s", (long) gap_open, (long) gap_extend, matrix_name);
	BlastConstructErrorMessage("BlastKarlinBlkGappedCalc", buffer, 2, error_return);
	for (index=0; index<max_number_values; index++)
	{
		sprintf(buffer, "Gap existence and extension values of %ld and %ld are supported", (long) Nint(values[index][0]), (long) Nint(values[index][1]));
		BlastConstructErrorMessage("BlastKarlinBlkGappedCalc", buffer, 2, error_return);
	}

	BlastMatrixValuesDestruct(head);

	return 1;
}

/* 
        Everything below here was (more or less) copied from the old 
        karlin.c and could work separately from the stuff above. 
*/ 

/**************** Statistical Significance Parameter Subroutine ****************

    Version 1.0     February 2, 1990
    Version 1.2     July 6,     1990

    Program by:     Stephen Altschul

    Address:        National Center for Biotechnology Information
                    National Library of Medicine
                    National Institutes of Health
                    Bethesda, MD  20894

    Internet:       altschul@ncbi.nlm.nih.gov

See:  Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
    Significance of Molecular Sequence Features by Using General Scoring
    Schemes,"  Proc. Natl. Acad. Sci. USA 87 (1990), 2264-2268.

    Computes the parameters lambda and K for use in calculating the
    statistical significance of high-scoring segments or subalignments.

    The scoring scheme must be integer valued.  A positive score must be
    possible, but the expected (mean) score must be negative.

    A program that calls this routine must provide the value of the lowest
    possible score, the value of the greatest possible score, and a pointer
    to an array of probabilities for the occurence of all scores between
    these two extreme scores.  For example, if score -2 occurs with
    probability 0.7, score 0 occurs with probability 0.1, and score 3
    occurs with probability 0.2, then the subroutine must be called with
    low = -2, high = 3, and pr pointing to the array of values
    { 0.7, 0.0, 0.1, 0.0, 0.0, 0.2 }.  The calling program must also provide
    pointers to lambda and K; the subroutine will then calculate the values
    of these two parameters.  In this example, lambda=0.330 and K=0.154.

    The parameters lambda and K can be used as follows.  Suppose we are
    given a length N random sequence of independent letters.  Associated
    with each letter is a score, and the probabilities of the letters
    determine the probability for each score.  Let S be the aggregate score
    of the highest scoring contiguous segment of this sequence.  Then if N
    is sufficiently large (greater than 100), the following bound on the
    probability that S is greater than or equal to x applies:

            P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ].

    In other words, the p-value for this segment can be written as
    1-exp[-KN*exp(-lambda*S)].

    This formula can be applied to pairwise sequence comparison by assigning
    scores to pairs of letters (e.g. amino acids), and by replacing N in the
    formula with N*M, where N and M are the lengths of the two sequences
    being compared.

    In addition, letting y = KN*exp(-lambda*S), the p-value for finding m
    distinct segments all with score >= S is given by:

                           2             m-1           -y
            1 - [ 1 + y + y /2! + ... + y   /(m-1)! ] e

    Notice that for m=1 this formula reduces to 1-exp(-y), which is the same
    as the previous formula.

*******************************************************************************/

Int2 LIBCALL
BlastKarlinBlkCalc(BLAST_KarlinBlkPtr kbp, BLAST_ScoreFreqPtr sfp)
{
	

	if (kbp == NULL || sfp == NULL)
		return 1;

	/* Calculate the parameter Lambda */

	kbp->Lambda_real = kbp->Lambda = BlastKarlinLambdaNR(sfp);
	if (kbp->Lambda < 0.)
		goto ErrExit;


	/* Calculate H */

	kbp->H_real = kbp->H = BlastKarlinLtoH(sfp, kbp->Lambda);
	if (kbp->H < 0.)
		goto ErrExit;


	/* Calculate K and log(K) */

	kbp->K_real = kbp->K = BlastKarlinLHtoK(sfp, kbp->Lambda, kbp->H);
	if (kbp->K < 0.)
		goto ErrExit;
	kbp->logK_real = kbp->logK = log(kbp->K);

	/* Normal return */
	return 0;

ErrExit:
	kbp->Lambda = kbp->H = kbp->K
		= kbp->Lambda_real = kbp->H_real = kbp->K_real = -1.;
#ifdef BLASTKAR_HUGE_VAL
	kbp->logK_real = kbp->logK = BLASTKAR_HUGE_VAL;
#else
	kbp->logK_real = kbp->logK = 1.e30;
#endif
	return 1;
}
#define DIMOFP0	(iter*range + 1)
#define DIMOFP0_MAX (BLAST_KARLIN_K_ITER_MAX*BLAST_SCORE_RANGE_MAX+1)

Nlm_FloatHi
BlastKarlinLHtoK(BLAST_ScoreFreqPtr sfp, Nlm_FloatHi	lambda, Nlm_FloatHi H)
{
#ifndef BLAST_KARLIN_STACKP
	Nlm_FloatHi	PNTR P0 = NULL;
#else
	Nlm_FloatHi	P0 [DIMOFP0_MAX];
#endif
	BLAST_Score	low;	/* Lowest score (must be negative) */
	BLAST_Score	high;	/* Highest score (must be positive) */
	Nlm_FloatHi	K;			/* local copy of K */
	Nlm_FloatHi	ratio;
	int		i, j;
	BLAST_Score	range, lo, hi, first, last;
	register Nlm_FloatHi	sum;
	Nlm_FloatHi	Sum, av, oldsum, oldsum2;
	int		iter;
	Nlm_FloatHi	sumlimit;
	Nlm_FloatHi	PNTR p, PNTR ptrP, PNTR ptr1, PNTR ptr2, PNTR ptr1e;
	Nlm_FloatHi	etolami, etolam;

	if (lambda <= 0. || H <= 0.) {
		return -1.;
	}

	if (sfp->score_avg >= 0.0) {
		return -1.;
	}

	low = sfp->obs_min;
	high = sfp->obs_max;
	range = high - low;

	av = H/lambda;
	etolam = exp((Nlm_FloatHi)lambda);
	if (low == -1 || high == 1) {
		if (high == 1)
			K = av;
		else
			K = (sfp->score_avg * sfp->score_avg) / av;
		return K * (1.0 - 1./etolam);
	}

	sumlimit = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;

	iter = BLAST_KARLIN_K_ITER_MAX;

	if (DIMOFP0 > DIMOFP0_MAX) {
		return -1.;
	}
#ifndef BLAST_KARLIN_STACKP
	P0 = (Nlm_FloatHi PNTR)MemNew(DIMOFP0 * sizeof(*P0));
	if (P0 == NULL)
		return -1.;
#else
	Nlm_MemSet((CharPtr)P0, 0, DIMOFP0*sizeof(P0[0]));
#endif

	Sum = 0.;
	lo = hi = 0;
	p = &sfp->sprob[low];
	P0[0] = sum = oldsum = oldsum2 = 1.;
    for (j = 0; j < iter && sum > sumlimit; Sum += sum /= ++j) {
        first = last = range;
		lo += low;
		hi += high;
        for (ptrP = P0+(hi-lo); ptrP >= P0; *ptrP-- =sum) {
            ptr1 = ptrP - first;
            ptr1e = ptrP - last;
            ptr2 = p + first;
            for (sum = 0.; ptr1 >= ptr1e; )
                sum += *ptr1--  *  *ptr2++;
            if (first)
                --first;
            if (ptrP - P0 <= range)
                --last;
        }
		etolami = Nlm_Powi((Nlm_FloatHi)etolam, lo - 1);
        for (sum = 0., i = lo; i != 0; ++i) {
			etolami *= etolam;
            sum += *++ptrP * etolami;
		}
        for (; i <= hi; ++i)
            sum += *++ptrP;
		oldsum2 = oldsum;
		oldsum = sum;
    }

	/* Terms of geometric progression added for correction */
	ratio = oldsum / oldsum2;
	if (ratio >= (1.0 - sumlimit*0.001)) {
		K = -1.;
		goto CleanUp;
	}
	sumlimit *= 0.01;
	while (sum > sumlimit) {
		oldsum *= ratio;
		Sum += sum = oldsum / ++j;
	}

/* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
Karlin&Altschul (1990) */
    	for (i = 1, j = -low; i <= range && j > 1; ++i)
        	if (p[i])
            		j = Nlm_Gcd(j, i);

	if (j*etolam > 0.05) 
	{
		etolami = Nlm_Powi((Nlm_FloatHi)etolam, -j);
    		K = j*exp((Nlm_FloatHi)-2.0*Sum) / (av*(1.0 - etolami));
	}
	else
	    K = -j*exp((Nlm_FloatHi)-2.0*Sum) / (av*Nlm_Expm1(-j*(Nlm_FloatHi)lambda));

CleanUp:
#ifndef BLAST_KARLIN_K_STACKP
	if (P0 != NULL)
		MemFree(P0);
#endif
	return K;
}

/*
	BlastKarlinLambdaBis

	Calculate Lambda using the bisection method (slow).
*/
Nlm_FloatHi
BlastKarlinLambdaBis(BLAST_ScoreFreqPtr sfp)
{
	register Nlm_FloatHi	PNTR sprob;
	Nlm_FloatHi	lambda, up, newval;
	BLAST_Score	i, low, high;
	int		j;
	register Nlm_FloatHi	sum, x0, x1;

	if (sfp->score_avg >= 0.) {
		return -1.;
	}
	low = sfp->obs_min;
	high = sfp->obs_max;
	if (BlastScoreChk(low, high) != 0)
		return -1.;

	sprob = sfp->sprob;
	up = BLAST_KARLIN_LAMBDA0_DEFAULT;
	for (lambda=0.; ; ) {
		up *= 2;
		x0 = exp((Nlm_FloatHi)up);
		x1 = Nlm_Powi((Nlm_FloatHi)x0, low - 1);
		if (x1 > 0.) {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * (x1 *= x0);
		}
		else {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * exp(up * i);
		}
		if (sum >= 1.0)
			break;
		lambda = up;
	}

	for (j=0; j<BLAST_KARLIN_LAMBDA_ITER_DEFAULT; ++j) {
		newval = (lambda + up) / 2.;
		x0 = exp((Nlm_FloatHi)newval);
		x1 = Nlm_Powi((Nlm_FloatHi)x0, low - 1);
		if (x1 > 0.) {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * (x1 *= x0);
		}
		else {
			for (sum=0., i=low; i<=high; ++i)
				sum += sprob[i] * exp(newval * i);
		}
		if (sum > 1.0)
			up = newval;
		else
			lambda = newval;
	}
	return (lambda + up) / 2.;
}

/******************* Fast Lambda Calculation Subroutine ************************
	Version 1.0	May 16, 1991
	Program by:	Stephen Altschul

	Uses Newton-Raphson method (fast) to solve for Lambda, given an initial
	guess (lambda0) obtained perhaps by the bisection method.
*******************************************************************************/

Nlm_FloatHi
BlastKarlinLambdaNR(BLAST_ScoreFreqPtr sfp)
{
	BLAST_Score	low;			/* Lowest score (must be negative)  */
	BLAST_Score	high;			/* Highest score (must be positive) */
	int		j;
	BLAST_Score	i;
	Nlm_FloatHi PNTR	sprob;
	Nlm_FloatHi	lambda0, sum, slope, temp, x0, x1, amt;

	low = sfp->obs_min;
	high = sfp->obs_max;
	if (sfp->score_avg >= 0.) {	/* Expected score must be negative */
		return -1.0;
	}
	if (BlastScoreChk(low, high) != 0)
		return -1.;

	lambda0 = BLAST_KARLIN_LAMBDA0_DEFAULT;

	/* Calculate lambda */

	sprob = sfp->sprob;
	for (j=0; j<20; ++j) { /* limit of 20 should never be close-approached */
		sum = -1.0;
		slope = 0.0;
		if (lambda0 < 0.01)
			break;
		x0 = exp((Nlm_FloatHi)lambda0);
		x1 = Nlm_Powi((Nlm_FloatHi)x0, low - 1);
		if (x1 == 0.)
			break;
		for (i=low; i<=high; ++i) {
			sum += (temp = sprob[i] * (x1 *= x0));
			slope += temp * i;
		}
		lambda0 -= (amt = sum/slope);
		if (ABS(amt/lambda0) < BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT) {
			/*
			Does it appear that we may be on the verge of converging
			to the ever-present, zero-valued solution?
			*/
			if (lambda0 > BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT)
				return lambda0;
			break;
		}
	}
	return BlastKarlinLambdaBis(sfp);
}

/*
	BlastKarlinLtoH

	Calculate H, the relative entropy of the p's and q's
*/
Nlm_FloatHi LIBCALL
BlastKarlinLtoH(BLAST_ScoreFreqPtr sfp, Nlm_FloatHi	lambda)
{
	BLAST_Score	score;
	Nlm_FloatHi	av, etolam, etolami;

	if (lambda < 0.) {
		return -1.;
	}
	if (BlastScoreChk(sfp->obs_min, sfp->obs_max) != 0)
		return -1.;

	etolam = exp((Nlm_FloatHi)lambda);
	etolami = Nlm_Powi((Nlm_FloatHi)etolam, sfp->obs_min - 1);
	if (etolami > 0.) 
	{
	    av = 0.0;
	    for (score=sfp->obs_min; score<=sfp->obs_max; score++)
   			av += sfp->sprob[score] * score * (etolami *= etolam);
	}
	else 
	{
	    av = 0.0;
	    for (score=sfp->obs_min; score<=sfp->obs_max; score++)
   			av += sfp->sprob[score] * score * exp(lambda * score);
	}

    	return lambda * av;
}


static Nlm_FloatHi
BlastGapDecayInverse(Nlm_FloatHi pvalue, unsigned nsegs, Nlm_FloatHi decayrate)
{
	if (decayrate <= 0. || decayrate >= 1. || nsegs == 0)
		return pvalue;

	return pvalue * (1. - decayrate) * Nlm_Powi(decayrate, nsegs - 1);
}

static Nlm_FloatHi
BlastGapDecay(Nlm_FloatHi pvalue, unsigned nsegs, Nlm_FloatHi decayrate)
{
	if (decayrate <= 0. || decayrate >= 1. || nsegs == 0)
		return pvalue;

	return pvalue / ((1. - decayrate) * Nlm_Powi(decayrate, nsegs - 1));
}

/*
	BlastCutoffs

	Calculate the cutoff score, S, and the highest expected score.

	WRG (later modified by TLM).
*/
Int2 LIBCALL
BlastCutoffs(BLAST_ScorePtr S, /* cutoff score */
	Nlm_FloatHi PNTR E, /* expected no. of HSPs scoring at or above S */
	BLAST_KarlinBlkPtr kbp,
	Nlm_FloatHi qlen, /* length of query sequence */
	Nlm_FloatHi dblen, /* length of database or database sequence */
	Nlm_Boolean dodecay) /* TRUE ==> use gapdecay feature */
{
	return BlastCutoffs_simple(S, E, kbp, qlen*dblen, dodecay);
}

Int2 LIBCALL
BlastCutoffs_simple(BLAST_ScorePtr S, /* cutoff score */
	Nlm_FloatHi PNTR E, /* expected no. of HSPs scoring at or above S */
	BLAST_KarlinBlkPtr kbp,
	Nlm_FloatHi searchsp, /* size of search space. */
	Nlm_Boolean dodecay) /* TRUE ==> use gapdecay feature */
{
	BLAST_Score	s = *S, es;
	Nlm_FloatHi	e = *E, esave;
	Boolean		s_changed = FALSE;

	if (kbp->Lambda == -1. || kbp->K == -1. || kbp->H == -1.)
		return 1;

	/*
	Calculate a cutoff score, S, from the Expected
	(or desired) number of reported HSPs, E.
	*/
	es = 1;
	esave = e;
	if (e > 0.) 
	{
		if (dodecay)
			e = BlastGapDecayInverse(e, 1, 0.5);
		es = BlastKarlinEtoS_simple(e, kbp, searchsp);
	}
	/*
	Pick the larger cutoff score between the user's choice
	and that calculated from the value of E.
	*/
	if (es > s) {
		s_changed = TRUE;
		*S = s = es;
	}

	/*
		Re-calculate E from the cutoff score, if E going in was too high
	*/
	if (esave <= 0. || !s_changed) 
	{
		e = BlastKarlinStoE_simple(s, kbp, searchsp);
		if (dodecay)
			e = BlastGapDecay(e, 1, 0.5);
		*E = e;
	}

	return 0;
}

/*
	BlastKarlinEtoS() -- given an Expect value, return the associated cutoff score

	Error return value is BLAST_SCORE_MIN
*/
BLAST_Score LIBCALL
BlastKarlinEtoS(Nlm_FloatHi	E,	/* Expect value */
	BLAST_KarlinBlkPtr	kbp,
	Nlm_FloatHi	qlen,	/* length of query sequence */
	Nlm_FloatHi	dblen)	/* length of database */
{
	return BlastKarlinEtoS_simple(E, kbp, qlen*dblen);
}


BLAST_Score LIBCALL
BlastKarlinEtoS_simple(Nlm_FloatHi	E,	/* Expect value */
	BLAST_KarlinBlkPtr	kbp,
	Nlm_FloatHi	searchsp)	/* size of search space */
{

	Nlm_FloatHi	Lambda, K, H; /* parameters for Karlin statistics */
	BLAST_Score	S;

	Lambda = kbp->Lambda;
	K = kbp->K;
	H = kbp->H;
	if (Lambda < 0. || K < 0. || H < 0.) 
	{
		return BLAST_SCORE_MIN;
	}

	S = (BLAST_Score) (ceil( log((Nlm_FloatHi)(K * searchsp / E)) / Lambda ));
	return S;
}

/*
	BlastKarlinStoP
	
	Calculate the probability (as opposed to expectation)
	of achieving a particular score.

	On error, return value is -1. (same as BlastKarlinStoE()).
*/
Nlm_FloatHi LIBCALL
BlastKarlinStoP(BLAST_Score S,
		BLAST_KarlinBlkPtr kbp,
		Nlm_FloatHi	qlen,	/* length of query sequence */
		Nlm_FloatHi	dblen)	/* length of database */
{
	return BlastKarlinStoP_simple(S, kbp, qlen*dblen);
}

Nlm_FloatHi LIBCALL
BlastKarlinStoP_simple(BLAST_Score S,
		BLAST_KarlinBlkPtr kbp,
		Nlm_FloatHi	searchsp)	/* size of search space. */
{
	Nlm_FloatHi	x, p;

	x = BlastKarlinStoE_simple(S, kbp, searchsp);
	if (x == -1.)
		return x;
	p = BlastKarlinEtoP(x);

	return p;
}

/*
	BlastKarlinStoE() -- given a score, return the associated Expect value

	Error return value is -1.
*/
Nlm_FloatHi LIBCALL
BlastKarlinStoE(BLAST_Score S,
		BLAST_KarlinBlkPtr kbp,
		Nlm_FloatHi	qlen,	/* length of query sequence */
		Nlm_FloatHi	dblen)	/* length of database */
{
	return BlastKarlinStoE_simple(S, kbp, qlen*dblen);
}

Nlm_FloatHi LIBCALL
BlastKarlinStoE_simple(BLAST_Score S,
		BLAST_KarlinBlkPtr kbp,
		Nlm_FloatHi	searchsp)	/* size of search space. */
{
	Nlm_FloatHi	Lambda, K, H; /* parameters for Karlin statistics */

	Lambda = kbp->Lambda;
	K = kbp->K;
	H = kbp->H;
	if (Lambda < 0. || K < 0. || H < 0.) {
		return -1.;
	}

	return searchsp * exp((Nlm_FloatHi)(-Lambda * S) + kbp->logK);
}

/*
BlastKarlinPtoE -- convert a P-value to an Expect value
*/
Nlm_FloatHi LIBCALL
BlastKarlinPtoE(Nlm_FloatHi p)
{
        if (p < 0. || p > 1.0) 
	{
                return INT4_MIN;
        }

	if (p == 1)
		return INT4_MAX;

        return -Nlm_Log1p(-p);
}

/*
BlastKarlinEtoP -- convert an Expect value to a P-value
*/
Nlm_FloatHi LIBCALL
BlastKarlinEtoP(Nlm_FloatHi x)
{
	return -Nlm_Expm1((Nlm_FloatHi)-x);
}

/*
	BlastKarlinStoLen()

	Given a score, return the length expected for an HSP of that score
*/
Nlm_FloatHi LIBCALL
BlastKarlinStoLen(BLAST_KarlinBlkPtr kbp, BLAST_Score S)
{
	return kbp->Lambda * S / kbp->H;
}

static Nlm_FloatHi	tab2[] = { /* table for r == 2 */
0.01669,  0.0249,   0.03683,  0.05390,  0.07794,  0.1111,   0.1559,   0.2146,   
0.2890,   0.3794,   0.4836,   0.5965,   0.7092,   0.8114,   0.8931,   0.9490,   
0.9806,   0.9944,   0.9989
		};

static Nlm_FloatHi	tab3[] = { /* table for r == 3 */
0.9806,   0.9944,   0.9989,   0.0001682,0.0002542,0.0003829,0.0005745,0.0008587,
0.001278, 0.001893, 0.002789, 0.004088, 0.005958, 0.008627, 0.01240,  0.01770,  
0.02505,  0.03514,  0.04880,  0.06704,  0.09103,  0.1220,   0.1612,   0.2097,   
0.2682,   0.3368,   0.4145,   0.4994,   0.5881,   0.6765,   0.7596,   0.8326,   
0.8922,   0.9367,   0.9667,   0.9846,   0.9939,   0.9980
		};

static Nlm_FloatHi	tab4[] = { /* table for r == 4 */
2.658e-07,4.064e-07,6.203e-07,9.450e-07,1.437e-06,2.181e-06,3.302e-06,4.990e-06,
7.524e-06,1.132e-05,1.698e-05,2.541e-05,3.791e-05,5.641e-05,8.368e-05,0.0001237,
0.0001823,0.0002677,0.0003915,0.0005704,0.0008275,0.001195, 0.001718, 0.002457,
0.003494, 0.004942, 0.006948, 0.009702, 0.01346,  0.01853,  0.02532,  0.03431,
0.04607,  0.06128,  0.08068,  0.1051,   0.1352,   0.1719,   0.2157,   0.2669,
0.3254,   0.3906,   0.4612,   0.5355,   0.6110,   0.6849,   0.7544,   0.8168,
0.8699,   0.9127,   0.9451,   0.9679,   0.9827,   0.9915,   0.9963
		};

static Nlm_FloatHi PNTR table[] = { tab2, tab3, tab4 };
static short tabsize[] = { DIM(tab2)-1, DIM(tab3)-1, DIM(tab4)-1 };

static Nlm_FloatHi LIBCALL f PROTO((Nlm_FloatHi,Nlm_VoidPtr));
static Nlm_FloatHi LIBCALL g PROTO((Nlm_FloatHi,Nlm_VoidPtr));


/*
    Estimate the Sum P-value by calculation or interpolation, as appropriate.
	Approx. 2-1/2 digits accuracy minimum throughout the range of r, s.
	r = number of segments
	s = total score (in nats), adjusted by -r*log(KN)
*/
Nlm_FloatHi LIBCALL
BlastSumP(Int4 r, Nlm_FloatHi s)
{
	Int4		i, r1, r2;
	Nlm_FloatHi	a;

	if (r == 1)
		return -Nlm_Expm1(-exp(-s));

	if (r <= 4) {
		if (r < 1)
			return 0.;
		r1 = r - 1;
		if (s >= r*r+r1) {
			a = Nlm_LnGammaInt(r+1);
			return r * exp(r1*log(s)-s-a-a);
		}
		if (s > -2*r) {
			/* interpolate */
			i = (Int4) (a = s+s+(4*r));
			a -= i;
			i = tabsize[r2 = r - 2] - i;
			return a*table[r2][i-1] + (1.-a)*table[r2][i];
		}
		return 1.;
	}

	return BlastSumPCalc(r, s);
}

/*
    BlastSumPCalc

    Evaluate the following Nlm_FloatHi integral, where r = number of segments
    and s = the adjusted score in nats:

                    (r-2)         oo           oo
     Prob(r,s) =   r              -            -   (r-2)
                 -------------   |   exp(-y)  |   x   exp(-exp(x - y/r)) dx dy
                 (r-1)! (r-2)!  U            U
                                s            0
*/
static Nlm_FloatHi
BlastSumPCalc(int r, Nlm_FloatHi s)
{
	int		r1, itmin;
	Nlm_FloatHi	t, d, epsilon;
	Nlm_FloatHi	est_mean, mean, stddev, stddev4;
	Nlm_FloatHi	xr, xr1, xr2, logr;
	Nlm_FloatHi	args[6];

	epsilon = BLAST_SUMP_EPSILON_DEFAULT; /* accuracy for SumP calcs. */

	if (r == 1) {
		if (s > 8.)
			return exp(-s);
		return -Nlm_Expm1(-exp(-s));
	}
	if (r < 1)
		return 0.;

	xr = r;

	if (r < 8) {
		if (s <= -2.3*xr)
			return 1.;
	}
	else if (r < 15) {
			if (s <= -2.5*xr)
				return 1.;
	}
	else if (r < 27) {
			if (s <= -3.0*xr)
				return 1.;
	}
	else if (r < 51) {
			if (s <= -3.4*xr)
				return 1.;
	}
	else if (r < 101) {
			if (s <= -4.0*xr)
				return 1.;
	}

	/* stddev in the limit of infinite r, but quite good for even small r */
	stddev = sqrt(xr);
	stddev4 = 4.*stddev;
	xr1 = r1 = r - 1;

	if (r > 100) {
		/* Calculate lower bound on the mean using inequality log(r) <= r */
		est_mean = -xr * xr1;
		if (s <= est_mean - stddev4)
			return 1.;
	}

	/* mean is rather close to the mode, and the mean is readily calculated */
	/* mean in the limit of infinite r, but quite good for even small r */
	logr = log(xr);
	mean = xr * (1. - logr) - 0.5;
	if (s <= mean - stddev4)
		return 1.;

	if (s >= mean) {
		t = s + 6.*stddev;
		itmin = 1;
	}
	else {
		t = mean + 6.*stddev;
		itmin = 2;
	}

#define ARG_R args[0]
#define ARG_R2 args[1]
#define ARG_ADJ1 args[2]
#define ARG_ADJ2 args[3]
#define ARG_SDIVR args[4]
#define ARG_EPS args[5]

	ARG_R = xr;
	ARG_R2 = xr2 = r - 2;
	ARG_ADJ1 = xr2*logr - Nlm_LnGammaInt(r1) - Nlm_LnGammaInt(r);
	ARG_EPS = epsilon;

	do {
		d = Nlm_RombergIntegrate(g, args, s, t, epsilon, 0, itmin);
#ifdef BLASTKAR_HUGE_VAL
		if (d == BLASTKAR_HUGE_VAL)
			return d;
#endif
	} while (s < mean && d < 0.4 && itmin++ < 4);

	return (d < 1. ? d : 1.);
}

static Nlm_FloatHi LIBCALL
g(Nlm_FloatHi	s, Nlm_VoidPtr	vp)
{
	register Nlm_FloatHi PNTR	args = vp;
	Nlm_FloatHi	mx;
	
	ARG_ADJ2 = ARG_ADJ1 - s;
	ARG_SDIVR = s / ARG_R;	/* = s / r */
	mx = (s > 0. ? ARG_SDIVR + 3. : 3.);
	return Nlm_RombergIntegrate(f, vp, 0., mx, ARG_EPS, 0, 1);
}

static Nlm_FloatHi LIBCALL
f(Nlm_FloatHi	x, Nlm_VoidPtr	vp)
{
	register Nlm_FloatHi PNTR	args = vp;
	register Nlm_FloatHi	y;

	y = exp(x - ARG_SDIVR);
#ifdef BLASTKAR_HUGE_VAL
	if (y == BLASTKAR_HUGE_VAL)
		return 0.;
#endif
	if (ARG_R2 == 0.)
		return exp(ARG_ADJ2 - y);
	if (x == 0.)
		return 0.;
	return exp(ARG_R2*log(x) + ARG_ADJ2 - y);
}

/*
	Calculates the p-value for alignments with "small" gaps (typically
	under fifty residues/basepairs) following ideas of Stephen Altschul's.
	"gap" gives the size of this gap, "gap_prob" is the probability
	of this model of gapping being correct (it's thought that gap_prob
	will generally be 0.5).  "num" is the number of HSP's involved, "sum" 
	is the "raw" sum-score of these HSP's. "subject_len" is the (effective)
	length of the database sequence, "query_length" is the (effective) 
	length of the query sequence.  min_length_one specifies whether one
	or 1/K will be used as the minimum expected length.

*/

Nlm_FloatHi LIBCALL
BlastSmallGapSumE(BLAST_KarlinBlkPtr kbp, Int4 gap, Nlm_FloatHi gap_prob, Nlm_FloatHi gap_decay_rate, Int2 num, Int4 sum, Nlm_FloatHi score_prime, Int4 query_length, Int4 subject_length, Boolean min_length_one)

{

	Nlm_FloatHi sum_p, sum_e;
		
	score_prime -= kbp->logK + log((Nlm_FloatHi)subject_length*(Nlm_FloatHi)query_length) + (num-1)*(kbp->logK + 2*log((Nlm_FloatHi)gap));
	score_prime -= log(Nlm_Factorial((int) num)); 

	sum_p = BlastSumP(num, score_prime);

	sum_e = BlastKarlinPtoE(sum_p);

	sum_e = sum_e/((1.0-gap_decay_rate)*Nlm_Powi(gap_decay_rate, (num-1)));

	if (num > 1)
	{
		if (gap_prob == 0.0)
			sum_e = INT4_MAX;
		else
			sum_e = sum_e/gap_prob;
	}

	return sum_e;
}

/*
	Calculates the p-values for alignments with "large" gaps (i.e., 
	infinite) followings an idea of Stephen Altschul's.
*/

Nlm_FloatHi LIBCALL
BlastLargeGapSumE(BLAST_KarlinBlkPtr kbp, Nlm_FloatHi gap_prob, Nlm_FloatHi gap_decay_rate, Int2 num, Int4 sum, Nlm_FloatHi score_prime, Int4 query_length, Int4 subject_length, Boolean old_stats)

{

	Nlm_FloatHi sum_p, sum_e;
/* The next two variables are for compatability with Warren's code. */
	Nlm_FloatHi lcl_subject_length, lcl_query_length;

        lcl_query_length = (Nlm_FloatHi) query_length;
        lcl_subject_length = (Nlm_FloatHi) subject_length;

	score_prime -= num*(kbp->logK + log(lcl_subject_length*lcl_query_length)) - log(Nlm_Factorial((int) num)); 

	sum_p = BlastSumP(num, score_prime);

	sum_e = BlastKarlinPtoE(sum_p);

	sum_e = sum_e/((1.0-gap_decay_rate)*Nlm_Powi(gap_decay_rate, (num-1)));

	if (num > 1)
	{
		if (gap_prob == 1.0)
			sum_e = INT4_MAX;
		else
			sum_e = sum_e/(1.0 - gap_prob);
	}

	return sum_e;
}

/********************************************************************
*
*	The following function, from Stephen Altschul, calculates a 
*	pseudoscore from a p-vlaue and n, the product of the database
*	and query lengths.
*	double	p;		 p-value	
*	double	n;		 search space 
*********************************************************************/
/* Move the following constant into blast.h??, or only the last one. */
#define		PSCALE 20.0
#define 	PSEUDO_SCORE_MAX 32767
#define 	SYBASE_MIN 1.0e-300

Int2 
ConvertPtoPseudoS(Nlm_FloatHi p, Nlm_FloatHi n)
{
	Int2	s;
	Nlm_FloatHi	E;

/* If p is 1.0, then E is very large and E/n is about one. */
	if (p > 0.99)
		return 0.5;
/* If p is very small, the highest score should be returned. */
	else if (p < SYBASE_MIN)
		return PSEUDO_SCORE_MAX;

/* E = -ln(1-p); the following calculates it. */
        E = -Nlm_Log1p(-p);

	s= ConvertEtoPseudoS (E, n);

	return s;
}

/*******************************************************************
*
*	This function calculates a pseudo-score from an E value.
*	The searchsp is the product of the query and database
*	lengths.  As the E value is directly related to the search
*	space, this effectively scales the database size out of
*	the calculation of the pseudo-score.
*******************************************************************/
Int2 
ConvertEtoPseudoS(Nlm_FloatHi E, Nlm_FloatHi searchsp)
{
	Int2	s;

/* If E is very small, the highest score should be returned. */
	if (E < SYBASE_MIN)
		return PSEUDO_SCORE_MAX;

	if (E > searchsp)
		return 0;

/* 0.5 is added to make sure this is rounded up, not down. */
	s= 0.5-PSCALE*log(E/searchsp);

	return s;
}


 
/*
Given a pseudoscore, a subroutine for calculating an E-value for a comparison
of size n (the product of the sequence length) is:
*/

Nlm_FloatHi 
ConvertPseudoStoE(Int2 s, Nlm_FloatHi n)
{
	return n*exp(-s/PSCALE);
}

/****************************************************************************
*
*	This function attaches a new ScorePtr to the "*old" one passed
*	in.  If *old is NULL, then *old is set equal to the new ptr.
*	The new pointer (NOT the head of the chain) is returned.
*
*	If the value is an int, then set prob to zero and pass the value in
*	as score; if the value is a Nlm_FloatHi, pass it in as prob.
*
*	The type of score is stored in an Object-id.str
****************************************************************************/

ScorePtr
MakeBlastScore (ScorePtr PNTR old, CharPtr scoretype, Nlm_FloatHi prob, Int4 score)

{
	ScorePtr	scrp, scrp1;

	scrp = ScoreNew();

	if (score) 
	{
		scrp->choice = 1;
		scrp->value.intvalue = score;
	}
	else 
	{
		scrp->choice = 2;
		scrp->value.realvalue = (Nlm_FloatHi) prob;
	}

	scrp->id = ObjectIdNew();
	scrp->id->str = StringSave(scoretype);


	if (*old == NULL)
		*old = scrp;
	else
	{
		for (scrp1=*old; scrp1->next; scrp1=scrp1->next)
			;
		scrp1->next = scrp;	
	}
		

	return scrp;
}

