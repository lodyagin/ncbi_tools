/*
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
*/

/*****************************************************************************

File name: posit2.c

Author: Alejandro Schaffer

Contents: utilities for makematrices.


*****************************************************************************/


#include<ncbi.h>
#include <blastpri.h>
#include<objcode.h>
#include<objseq.h>
#include<objsset.h>
#include<sequtil.h>
#include <posit.h>
#include <txalign.h>
#include <profiles.h>

/*Used to check that diagnostics printing routine will work*/
#define PROTEIN_ALPHABET 26
#define EFFECTIVE_ALPHABET 20

#define POSIT_PERCENT 0.05

#define POSIT_NUM_ITERATIONS 10

#define POSIT_SCALE_FACTOR 1000


static BLAST_ScoreFreqPtr fillSfp(BLAST_Score **matrix, Int4 matrixLength, Nlm_FloatHi *queryProbArray, Nlm_FloatHi *scoreArray,  BLAST_ScoreFreqPtr return_sfp)
{
  Int4 minScore, maxScore; /*observed minimum and maximum scores*/
  Int4 i,j; /* indices */
  Nlm_FloatHi onePosFrac; /*1/matrix length as a double*/

  minScore = maxScore = 0;

  for(i = 0; i < matrixLength; i++) {
    for(j = 0 ; j < PROTEIN_ALPHABET; j++) {
      if ((matrix[i][j] != BLAST_SCORE_MIN) && (matrix[i][j] < minScore))
	minScore = matrix[i][j];
      if (matrix[i][j] > maxScore)
        maxScore = matrix[i][j];
    }
  }
  return_sfp->obs_min = minScore;
  return_sfp->obs_max = maxScore;
  for (i = 0; i < scoreRange; i++)
    scoreArray[i] = 0.0;
  return_sfp->sprob = &(scoreArray[-minScore]); /*center around 0*/
  onePosFrac = 1.0/ ((Nlm_FloatHi) matrixLength);
  for(i = 0; i < matrixLength; i++) {
    for (j = 0; j < PROTEIN_ALPHABET; j++) {
      if(matrix[i][j] >= minScore) {
        return_sfp->sprob[matrix[i][j]] += (onePosFrac * queryProbArray[j]);
      }
    }
  }
  return_sfp->score_avg = 0;
  for(i = minScore; i <= maxScore; i++)
    return_sfp->score_avg += i * return_sfp->sprob[i];
  return(return_sfp);
}


static void 
impalaScaleMatrix(BlastMatrixRescalePtr matrix_rescale, Nlm_FloatHi scalingFactor)
{
   Int4 dim1, dim2; /*number of rows and number of columns*/
   Int4 a,c; /*loop indices*/
   Boolean too_high=TRUE, done, first_time; /*control variables for binary search*/
   Nlm_FloatHi factor, factor_low=1.0, factor_high=1.0; /*multiplicative factors in binary search*/
   Nlm_FloatHi lambda, new_lambda; /*Karlin-Altschul parameter*/
   Int4 index; /*loop index for binary search*/
   Int4Ptr *private_matrix; /*pointer to locally manipulated version of the
                              matrix*/
   Int4Ptr *matrix;
   Nlm_FloatHi scalefactor; /*local version of amount of scaling*/
   Nlm_FloatHi divFactor; /*1/scalefacor*/
   BLAST_ScoreFreqPtr this_sfp, return_sfp; /*score frequency pointers to
                                              compute lambda*/
   Nlm_FloatHi *scoreArray; /*array of score probabilities*/

   scoreArray = (Nlm_FloatHi *) MemNew(scoreRange * sizeof(Nlm_FloatHi));
   return_sfp = (BLAST_ScoreFreqPtr) MemNew(1 * sizeof(BLAST_ScoreFreq));


   private_matrix = matrix_rescale->private_matrix;
   matrix = matrix_rescale->matrix;

/* Bracket the values. */
   dim1 = matrix_rescale->query_length;
   dim2 = matrix_rescale->alphabet_size;

   lambda = matrix_rescale->lambda_ideal/scalingFactor;
   divFactor = ((Nlm_FloatHi) POSIT_SCALE_FACTOR)/scalingFactor;

   done = FALSE;
   first_time = TRUE;
   factor = 1.0;
   while (done != TRUE)
   {
   	for(c = 0; c < dim1; c++)
   	{
       	    for(a = 0; a < dim2; a++)
	    {
		if (private_matrix[c][a] == BLAST_SCORE_MIN)
		{
			matrix[c][a] = BLAST_SCORE_MIN;
		}
		else
		{
			matrix[c][a] = (factor*private_matrix[c][a]) /divFactor;
		}
	    }
        }

	this_sfp =  fillSfp(matrix, dim1, matrix_rescale->standardProb, scoreArray, return_sfp);
	new_lambda = impalaKarlinLambdaNR(this_sfp, matrix_rescale->kbp_psi[0]->Lambda/scalingFactor);

	if (new_lambda > lambda)
	{
		if (first_time)
		{
			factor_high = 1.0 + POSIT_PERCENT;
			factor = factor_high;
			factor_low = 1.0;
			too_high = TRUE;
			first_time = FALSE;
		}
		else
		{
			if (too_high == FALSE)
				break;
			factor_high += (factor_high-1.0);
			factor = factor_high;
		}
	}
	else 
	{
		if (first_time)
		{
			factor_high = 1.0;
			factor_low = 1.0 - POSIT_PERCENT;
			factor = factor_low;
			too_high = FALSE;
			first_time = FALSE;
		}
		else
		{
			if (too_high == TRUE)
				break;
			factor_low += (factor_low-1.0);
			factor = factor_low;
		}
	}
   }

/* binary search for ten times. */
   for (index=0; index<POSIT_NUM_ITERATIONS; index++)
   {
        factor = 0.5*(factor_high+factor_low);
   	for(c = 0; c < dim1; c++)
   	{
       	    for(a = 0; a < dim2; a++)
	    {
		if (private_matrix[c][a] == BLAST_SCORE_MIN)
		{
			matrix[c][a] = BLAST_SCORE_MIN;
		}
		else
		{
			matrix[c][a] = (factor*private_matrix[c][a])/divFactor;
		}
	    }
   	}

	this_sfp =  fillSfp(matrix, dim1, matrix_rescale->standardProb, scoreArray, return_sfp);
	new_lambda = impalaKarlinLambdaNR(this_sfp, matrix_rescale->kbp_psi[0]->Lambda/scalingFactor);

	if (new_lambda > lambda)
	{
		factor_low = factor;
	}
	else
	{
		factor_high = factor;
	}
    }

   for(c = 0; c < dim1; c++)
     for(a = 0; a < dim2; a++) {
       if (BLAST_SCORE_MIN != private_matrix[c][a]) {
         matrix[c][a] = Nlm_Nint((Nlm_FloatHi) private_matrix[c][a] *
             factor / POSIT_SCALE_FACTOR);
       }
     }

   updateLambdaK(matrix_rescale, TRUE);

   scalefactor = ((Nlm_FloatHi) scalingFactor) / POSIT_SCALE_FACTOR;
   for(c = 0; c < dim1; c++)
     for(a = 0; a < dim2; a++) {
       if (BLAST_SCORE_MIN != private_matrix[c][a]) {
         private_matrix[c][a] = Nlm_Nint((Nlm_FloatHi) private_matrix[c][a] *
             factor * scalefactor);
       }
     }

   for(a = 0; a < matrix_rescale->alphabet_size; a++) {
     matrix[dim1][a] = BLAST_SCORE_MIN;
     private_matrix[dim1][a] = BLAST_SCORE_MIN;
   }

   MemFree(scoreArray);
   MemFree(return_sfp);
}

void LIBCALL impalaScaling(posSearchItems *posSearch, compactSearchItems * compactSearch, Nlm_FloatHi scalingFactor)
{
	BlastMatrixRescalePtr matrix_rescale;

	matrix_rescale = BlastMatrixRescaleNew(compactSearch->alphabetSize, 
						compactSearch->qlength,
						compactSearch->query,
						compactSearch->standardProb,
						posSearch->posMatrix,
						posSearch->posPrivateMatrix,
						compactSearch->kbp_std,
						compactSearch->kbp_psi,
						compactSearch->kbp_gap_std,
						compactSearch->kbp_gap_psi,
						compactSearch->lambda_ideal,
						compactSearch->K_ideal);

	impalaScaleMatrix(matrix_rescale, scalingFactor);

	matrix_rescale = BlastMatrixRescaleDestruct(matrix_rescale);

	return;
}

/*Some of the following checkpointing code is taken and adapted from
code written by K. Shriram for FASTLINK.
Reference:
 A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. 
 Avoiding Recomputation in Linkage Analysis,
 Human Heredity 44(1994), pp. 225-237. */


/* Front-ends to retrieve numbers. */

#define  getCkptNlm_FloatHi(d, ckptFile)  (getCkptNumber(&(d),sizeof(Nlm_FloatHi),ckptFile))
#define  getCkptInt4(i, ckptFile)         (getCkptNumber(&(i),sizeof(Int4),ckptFile))
#define  getCkptChar(c, ckptFile)         (getCkptNumber(&(c),sizeof(Char),ckptFile))
 

/*Code to get a vector of type Nlm_FloatHi*/

static void    getCkptNlmFloat_HiVector (Nlm_FloatHi * theVector, Int4 length, FILE * ckptFile)
{
  int  vectorRef ;
 
  for(vectorRef = 0; vectorRef < length; vectorRef++)
    getCkptNlm_FloatHi(theVector[vectorRef],ckptFile) ;
}

static void    getFreqVector (Nlm_FloatHi * theVector, Int4 length, FILE * ckptFile)
{
   int  vectorRef ;

   Int4 charOrder[EFFECTIVE_ALPHABET]; /*standard order of letters according to S. Altschul*/


   charOrder[0] =  1;  /*A*/
   charOrder[1] =  16; /*R*/
   charOrder[2] =  13; /*N*/  
   charOrder[3] =  4;  /*D*/ 
   charOrder[4] =  3;  /*C*/
   charOrder[5] =  15; /*Q*/
   charOrder[6] =  5; /*E*/ 
   charOrder[7] =  7;  /*G*/
   charOrder[8] =  8;  /*H*/
   charOrder[9] =  9;  /*I*/
   charOrder[10] = 11; /*L*/
   charOrder[11] = 10; /*K*/
   charOrder[12] = 12; /*M*/  
   charOrder[13] =  6; /*F*/
   charOrder[14] = 14; /*P*/
   charOrder[15] = 17; /*S*/
   charOrder[16] = 18; /*T*/
   charOrder[17] = 20; /*W*/
   charOrder[18] = 22; /*Y*/
   charOrder[19] = 19; /*V*/
 
  for(vectorRef = 0; vectorRef < length; vectorRef++)
    theVector[vectorRef] = 0;
  for(vectorRef = 0; vectorRef < EFFECTIVE_ALPHABET; vectorRef++)
    getCkptNlm_FloatHi(theVector[charOrder[vectorRef]],ckptFile) ;
}

/*Read a checkpoint from the end of a previous PSI-BLAST round, get
query length, query, and position-specific target frequencies.
Returns TRUE if checkpoint was read sucessfully and FALSE otherwise.
scalingFactor records by how much to scale output scores */

Boolean LIBCALL  impalaReadCheckpoint(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, ValNodePtr * error_return,
Nlm_FloatHi scalingFactor)
{
  FILE * checkFile; /*file in which to take the checkpoint*/
  Int4 length1, length2, c; /*length of query sequence, and an index for it*/
  Char  nextRes; /*next residue in stored copy of the query sequence*/
  Uint1Ptr oldQuery; /*array to hold the query sequence*/
  Int4 i,j; /*indices for position and character in alphabet*/

  BlastConstructErrorMessage("posReadCheckpoint", "Attempting to recover data from previous checkpoint\n", 1, error_return);
  checkFile = FileOpen(fileName, "r");  
  if (NULL == checkFile) {
    BlastConstructErrorMessage("posReadCheckpoint", "Could not open checkpoint file\n", 1, error_return);
    return(FALSE);
  }
  length1 = compactSearch->qlength;
  getCkptInt4(length2,checkFile);
  if (length1 != length2) {
    ErrPostEx(SEV_WARNING, 0, 0, "Invalid usage of checkpoint recovery; old query has length %ld, new query has length %ld", (long) length2,  (long) length1);
    BlastConstructErrorMessage("posReadCheckpoint", "Failed to recover data\n", 1, error_return);
    FileClose(checkFile);
    return(FALSE);
  }
  oldQuery = (Uint1Ptr) MemNew(length1 * sizeof(Uint1));
  if (NULL == oldQuery) {
    BlastConstructErrorMessage("posReadCheckpoint", "Failed to reconstruct previous query\n", 1, error_return);
    BlastConstructErrorMessage("posReadCheckpoint", "Failed to recover data\n", 1, error_return);
    FileClose(checkFile);
    return(FALSE);
  }  
  for(c = 0; c < length1; c++) {
    getCkptChar(nextRes, checkFile);
    oldQuery[c] = ResToInt(nextRes);
    if ((oldQuery[c] != compactSearch->query[c]) && (oldQuery[c] != Xchar)) {
      ErrPostEx(SEV_WARNING, 0, 0, "Stored query has a %c at position %ld, while new query has a %c there",getRes(oldQuery[c]), (long) c, getRes(compactSearch->query[c]));
      BlastConstructErrorMessage("posReadCheckpoint", "Failed to recover data\n", 1, error_return);
      MemFree(oldQuery);
      FileClose(checkFile);
      return(FALSE);
    }
  }
  posSearch->posMatrix = (BLAST_Score **) MemNew((length1 + 1) * sizeof(BLAST_Score *));
  posSearch->posPrivateMatrix = (BLAST_Score **) MemNew((length1 + 1) * sizeof(BLAST_Score *));
  posSearch->posFreqs = (Nlm_FloatHi **) MemNew((length1 + 1) * sizeof(Nlm_FloatHi *));
  if ((NULL == posSearch->posMatrix) || (NULL == posSearch->posPrivateMatrix) || (NULL == posSearch->posFreqs)) {

    BlastConstructErrorMessage("posReadCheckpoint", "Failed to allocate position-specific score matrix", 1, error_return);
    BlastConstructErrorMessage("posReadCheckpoint", "Failed to recover data\n", 1, error_return);
    MemFree(oldQuery);
    FileClose(checkFile);
    return(FALSE);
  }
  for(i = 0; i <= length1; i++) {
    posSearch->posMatrix[i] = (BLAST_Score *) MemNew(compactSearch->alphabetSize * sizeof(BLAST_Score));
    posSearch->posPrivateMatrix[i] = (BLAST_Score *) MemNew(compactSearch->alphabetSize * sizeof(BLAST_Score));
    posSearch->posFreqs[i] = (Nlm_FloatHi *) MemNew(compactSearch->alphabetSize * sizeof(Nlm_FloatHi));

    if ((NULL == posSearch->posMatrix[i]) || (NULL == posSearch->posPrivateMatrix[i]) || (NULL == posSearch->posFreqs[i])) {
      BlastConstructErrorMessage("posReadCheckpoint", "Failed to allocate position-specific score matrix", 1, error_return);
      BlastConstructErrorMessage("posReadCheckpoint", "Failed to recover data\n", 1, error_return);
      MemFree(oldQuery);
      FileClose(checkFile);
      return(FALSE);
    }
    for(j = 0; j < compactSearch->alphabetSize; j++) {
      posSearch->posFreqs[i][j] = 0.0;
    }
  }
  getCkptFreqMatrix(posSearch->posFreqs,length1,compactSearch->alphabetSize,checkFile);
  posFreqsToMatrix(posSearch,compactSearch);
  impalaScaling(posSearch, compactSearch, scalingFactor);
  BlastConstructErrorMessage("posReadCheckpoint", "Data recovered successfully\n", 1, error_return);
  MemFree(oldQuery);
  FileClose(checkFile);
  return(TRUE);
}




