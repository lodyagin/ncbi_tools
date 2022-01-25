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

File name: posit.c

Author: Alejandro Schaffer

Contents: utilities for position-based BLAST.


*****************************************************************************/
/* $Revision: 6.22 $ */
/* $Log: posit.c,v $
/* Revision 6.22  1998/12/09 18:51:51  madden
/* fixed counting bug in posCancel
/*
 * Revision 6.21  1998/09/28 12:31:31  madden
 * Used BlastConstructErrorMessage
 *
 * Revision 6.20  1998/09/09 21:18:33  madden
 * AS fixed warnings
 *
 * Revision 6.19  1998/09/09 16:09:20  madden
 * Changes for PHI-BLAST
 *
 * Revision 6.18  1998/08/26 18:07:00  kans
 * fixed -v -fd warnings (AS)
 *
 * Revision 6.17  1998/06/18 18:20:22  madden
 * Fixed typo in posConvergenceTest
 *
 * Revision 6.16  1998/06/14 19:43:02  madden
 * Added function posFreqsToInformation
 *
 * Revision 6.15  1998/06/12 20:38:48  madden
 * Fix for no hits to build model situation
 *
 * Revision 6.14  1998/06/09 19:38:16  madden
 * Changes rounddown to posit_rounddown to avoid conflict
 *
 * Revision 6.13  1998/04/24 19:29:30  madden
 * Moved rescaling code to blastool.c
 *
 * Revision 6.12  1998/03/25 22:36:17  egorov
 * Change type of posRepeatSequences
 *
 * Revision 6.11  1998/03/23 18:32:30  madden
 * Fix for zero/zero problem
 *
 * Revision 6.10  1998/02/06 18:34:17  madden
 * Added check that residue was not masked in posReadCheckpoint
 *
 * Revision 6.9  1998/02/03 15:57:28  madden
 * Cpos arg in posComputePseudoFreqs set to FALSE for WposComputation call
 *
 * Revision 6.8  1998/01/02 22:19:46  madden
 * Replaced printf by ErrPostEx of SEV_WARNING
 *
 * Revision 6.7  1997/12/23 21:07:06  madden
 * Changes for checkpointing
 *
 * Revision 6.6  1997/12/12 22:14:35  kans
 * changed round to rounddown to avoid CodeWarrior 68K collision
 *
 * Revision 6.5  1997/11/19 15:29:31  madden
 * Changed OS_UNIX ifdef to POSIT_DEBUG
 *
 * Revision 6.4  1997/09/08 13:33:29  madden
 * Added posEpsilon2 to check for small numbers
 *
 * Revision 6.3  1997/09/05 22:29:13  madden
 * Check for zero denominator and replace log(2) by NCBIMATH_LN2
 *
 * Revision 6.2  1997/09/02 22:23:01  madden
 * Removed redundant calls to updateLambdaK
 *
 * Revision 6.1  1997/08/27 21:18:18  madden
 * Fixed problem with deleted matrix
 *
 * Revision 6.0  1997/08/25 18:53:48  madden
 * Revision changed to 6.0
 *
 * Revision 1.22  1997/08/20 21:35:01  madden
 * ALL_ROUNDS replaced by Boolean
 *
 * Revision 1.21  1997/08/11 15:45:24  madden
 * eliminated obsolete fields
 *
 * Revision 1.20  1997/07/28 18:35:06  madden
 * Removed and ifdefed printf
 *
 * Revision 1.19  1997/06/27 19:14:01  madden
 * Fixed two bugs in posComputeSequenceWeights for the special case where all participating sequences are identical in a block
 *
 * Revision 1.17  1997/06/23 14:42:46  madden
 *  Made posComputeSequenceWeights faster by catching the special case where the set of participating sequences does not change from one column to the next.
 *
 * Revision 1.16  1997/05/29 20:35:47  madden
 * Eliminate duplicate sequences and alignments that are 98 perc. identical and ignore columns with all identical sequence weights.
 *
 * Revision 1.15  1997/05/27 20:26:09  madden
 * Fixed problem with matrix
 *
 * Revision 1.14  1997/05/23 20:52:50  madden
 * Fixed bug in setting of matrix for psi-blast
 *
 * Revision 1.13  1997/05/22 21:25:28  madden
 * fixed memory leaks
 *
 * Revision 1.12  1997/05/16 20:56:35  madden
 * replace hard coded numbers by defines
 *
 * Revision 1.11  1997/05/16 20:09:42  madden
 * Fixes for statistical problems
 *
 * Revision 1.10  1997/05/07 21:00:03  madden
 * Call to SeqId2OrdinalId replaces call to readdb_gi2seq
 *
 * Revision 1.9  1997/05/01 15:53:25  madden
 * Addition of extra KarlinBlk's for psi-blast
 *
 * Revision 1.8  1997/04/23  13:31:20  madden
 * Changed diagnostic output.
 *
 * Revision 1.7  1997/04/22  16:36:49  madden
 * Changes for use of psi-blast with www.
 *
 * Revision 1.6  1997/04/10  19:25:53  madden
 * Added casts, COMMAND_LINE replaced by ALL_ROUNDS.
 *
 * Revision 1.5  1997/04/09  20:01:53  madden
 * Functions CposComputation and WposComputation replace posComputations.
 *
 * Revision 1.4  1997/04/04  20:44:55  madden
 * Changed posComputation to return Int4Ptr *.
 *
 * Revision 1.3  1997/03/27  22:30:51  madden
 * Fix for Array Bounds Read.
 *
 * Revision 1.2  1997/03/11  14:38:40  madden
 * Fixes for GI's instead of ordinal numbers.
 *
 * Revision 1.1  1997/02/13  15:22:13  madden
 * Initial revision
 *
*/


#include<ncbi.h>
#include <blastpri.h>
#include<objcode.h>
#include<objseq.h>
#include<objsset.h>
#include<sequtil.h>
#include <posit.h>
#include <txalign.h>

extern BLAST_ScoreFreqPtr BlastScoreFreqNew (BLAST_Score score_min, BLAST_Score score_max);

extern BLAST_ScoreFreqPtr BlastScoreFreqDestruct (BLAST_ScoreFreqPtr sfp);

/*small constants to test against 0*/
#define posEpsilon 0.0001
#define posEpsilon2 0.0000001
/*Representation of a gap in a motif*/
#define GAP_CHAR 0
/*Used inside a seqAlign to reprsent the presence of a gap*/
#define GAP_HERE (-1)
/*Used to check that diagnostics printing routine will work*/
#define PROTEIN_ALPHABET 26
#define EFFECTIVE_ALPHABET 20

#define POSIT_PERCENT 0.05
#define POSIT_NUM_ITERATIONS 10

#define POSIT_SCALE_FACTOR 1000

#define POS_RESTING 0
#define POS_COUNTING 1

#define IDENTITY_RATIO 0.98

/*Allocate memory for  data structures inside posSearch used in
* position-specific caculations
* posSearch -- to be filled in 
* alphabetSize -- number of distinct characters used in the sequences
* querySize -- number of characters in the query sequence
* numSequences -- number of matching sequences potentially in the model */
static void posAllocateMemory(posSearchItems * posSearch, 
		       Int4 alphabetSize, Int4 querySize, Int4 numSequences)
{
  Int4 i, j;  /*loop indices*/

  posSearch->posCount = (Int4 *) MemNew(querySize * sizeof(Int4));
  if (NULL == posSearch->posCount)
    exit(EXIT_FAILURE);
  for(i = 0; i < querySize; i++)
    posSearch->posCount[i] = 0;

  posSearch->posC = (Int4 **) MemNew((querySize + 1) * sizeof(Int4 *));
  if (NULL == posSearch->posC)
    exit(EXIT_FAILURE);
  for(i = 0; i <= querySize; i++) {
    posSearch->posC[i] = (Int4 *) MemNew(alphabetSize * sizeof(Int4));
    if (NULL == posSearch->posC[i])
      exit(EXIT_FAILURE);   
    for(j = 0; j < alphabetSize; j++)
      posSearch->posC[i][j] = 0;
 
  }


  posSearch->posMatchWeights = (Nlm_FloatHi **) MemNew((querySize+1) * sizeof(Nlm_FloatHi *));
  if (NULL == posSearch->posMatchWeights)
    exit(EXIT_FAILURE);
  for (i = 0; i <= querySize ; i++) {
    posSearch->posMatchWeights[i] = (Nlm_FloatHi *) MemNew(alphabetSize * sizeof(Nlm_FloatHi));
    if (NULL == posSearch->posMatchWeights[i])
      exit(EXIT_FAILURE);
    for(j = 0; j < alphabetSize; j++) 
      posSearch->posMatchWeights[i][j] = 0.0;
  }  

  posSearch->posMatrix = (BLAST_Score **) MemNew((querySize + 1) * sizeof(BLAST_Score *));
  posSearch->posPrivateMatrix = (BLAST_Score **) MemNew((querySize + 1) * sizeof(BLAST_Score *));
  if (NULL == posSearch->posMatrix)
    exit(EXIT_FAILURE);
  for(i = 0; i <= querySize; i++) {
    posSearch->posMatrix[i] = (BLAST_Score *) MemNew(alphabetSize * sizeof(BLAST_Score));
    posSearch->posPrivateMatrix[i] = (BLAST_Score *) MemNew(alphabetSize * sizeof(BLAST_Score));
    if (NULL == posSearch->posMatrix[i])
      exit(EXIT_FAILURE);   
    for(j = 0; j < alphabetSize; j++)
      posSearch->posMatrix[i][j] = 0;
 
  }

  posSearch->posSigma = (Nlm_FloatHi *) MemNew((querySize) * sizeof(Nlm_FloatHi));
  if (NULL == posSearch->posSigma)
    exit(EXIT_FAILURE);
  for(i = 0; i < querySize; i++) {
    posSearch->posSigma[i] = 0.0;
  }

  posSearch->posIntervalSizes = (Int4 *) MemNew((querySize) * sizeof(Int4));
  if (NULL == posSearch->posIntervalSizes)
    exit(EXIT_FAILURE);
  for(i=0; i < querySize; i++)
    posSearch->posIntervalSizes[i] = 0;

  posSearch->posDescMatrixLength = numSequences;
  posSearch->posDescMatrix = (posDesc **) MemNew((numSequences + 1) * sizeof(posDesc *));
  if (NULL == posSearch->posDescMatrix)
    exit(EXIT_FAILURE);
  for (i = 0; i <= numSequences; i++) {
    posSearch->posDescMatrix[i] = (posDesc *) MemNew(querySize * sizeof(posDesc));
    if (NULL == posSearch->posDescMatrix[i])
      exit(EXIT_FAILURE);
    for(j = 0; j < querySize; j++) {
      posSearch->posDescMatrix[i][j].letter = UNUSED;
      posSearch->posDescMatrix[i][j].used = FALSE;
      posSearch->posDescMatrix[i][j].e_value = 1.0;
      posSearch->posDescMatrix[i][j].leftExtent = -1;
      posSearch->posDescMatrix[i][j].rightExtent = querySize;
    }
  }  
  posSearch->posExtents = (posDesc *) MemNew(querySize * sizeof(posDesc));
  if (NULL == posSearch->posExtents)
    exit(EXIT_FAILURE);
  for(j = 0; j < querySize; j++) {
    posSearch->posExtents[j].used = FALSE;
    posSearch->posExtents[j].leftExtent = -1;
    posSearch->posExtents[j].rightExtent = querySize;
  }
   posSearch->posA = (Nlm_FloatHi *) MemNew((numSequences+ 1) * sizeof(Nlm_FloatHi));
   if (NULL == posSearch->posA)
     exit(EXIT_FAILURE);
   posSearch->posRowSigma = (Nlm_FloatHi *) MemNew((numSequences + 1) * sizeof(Nlm_FloatHi));
   if (NULL == posSearch->posRowSigma)
     exit(EXIT_FAILURE);
}

static void freePosFreqs(Nlm_FloatHi ** posFreqs, Int4 length)
{
  Int4 i;

  for (i = 0; i <= length; i++)
    MemFree(posFreqs[i]);
  MemFree(posFreqs); 
}

/*Deallocate memory allocated in posReadCheckpoint
* posSearch -- pointer to record used in building the position-specific model
* querySize -- number of characters in the query sequence
*/
void LIBCALL posCheckpointFreeMemory(posSearchItems *posSearch, Int4 querySize)
{
  Int4 i;  /*loop index*/

  freePosFreqs(posSearch->posFreqs, querySize);
  for(i = 0; i <= querySize; i++){
    MemFree(posSearch->posMatrix[i]);
    MemFree(posSearch->posPrivateMatrix[i]);
  }
  MemFree(posSearch->posMatrix);
  MemFree(posSearch->posPrivateMatrix);
}

/*Deallocate memory allocated in posAllocateMemory
* posSearch -- pointer to record used in building the position-specific model
* querySize -- number of characters in the query sequence
*/
static void posFreeMemory(posSearchItems *posSearch, Int4 querySize)
{
  Int4 i;  /*loop index*/

  MemFree(posSearch->posCount);
  MemFree(posSearch->posExtents);
  MemFree(posSearch->posSigma);
  MemFree(posSearch->threshSequences);
  for(i = 0; i <= querySize; i++){
    MemFree(posSearch->posC[i]);
    MemFree(posSearch->posMatrix[i]);
    MemFree(posSearch->posPrivateMatrix[i]);
    MemFree(posSearch->posMatchWeights[i]);
  }
  MemFree(posSearch->posC);
  for(i = 0; i <= posSearch->posDescMatrixLength; i++)
    MemFree(posSearch->posDescMatrix[i]);
  MemFree(posSearch->posMatrix);
  MemFree(posSearch->posPrivateMatrix);
  MemFree(posSearch->posDescMatrix);
  MemFree(posSearch->posMatchWeights);
  MemFree(posSearch->posA);
  MemFree(posSearch->posRowSigma);
  MemFree(posSearch->posIntervalSizes);
  MemFree(posSearch->posUseSequences);
  freePosFreqs(posSearch->posFreqs,querySize);
}

/*Cleanup position-specific  data structures after one pass*/
void LIBCALL posCleanup(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
  posFreeMemory(posSearch, compactSearch->qlength);
}

/*Find which sequences that match in the i-th pass have an e-value below
  the specified threshold. These sequences will be used to make the
  score matrix for the next pass*/
static void findThreshSequences(posSearchItems *posSearch, BlastSearchBlkPtr search)
{

   Int4 numseq, seqIndex; /*number of sequences, index for them*/
   Int4 hspIndex, hspBnd; /*Index for hsps, number of hsps for one sequence*/


   numseq = search->result_struct->hitlist_count;
     /*Allocate boolean array to store values*/
   posSearch->threshSequences = (Int4 *) MemNew(numseq * sizeof(Int4)); 
   posSearch->posResultSequences = (Int4 *) MemNew(numseq * sizeof(Int4));
   posSearch->posResultsCounter = 0;
   posSearch->posMaxThresh = 0;
   if (NULL == posSearch->threshSequences)
     exit(EXIT_FAILURE);
   for(seqIndex = 0; seqIndex < numseq; seqIndex++)
    posSearch->threshSequences[seqIndex] = 0;
   for(seqIndex = 0; seqIndex < numseq; seqIndex++) {
     hspBnd = search->result_struct->results[seqIndex]->hspcnt;
      /*e-value for a sequence is the smallest e-value among the
        hsps matching a region of the sequence to the query*/
     for(hspIndex = 0; hspIndex < hspBnd; hspIndex++) {
       if(search->result_struct->results[seqIndex]->hsp_array[hspIndex].e_value < (search->pbp->ethresh)) {
         posSearch->threshSequences[seqIndex] = 1;
         posSearch->posMaxThresh = seqIndex + 1;
         posSearch->posResultSequences[posSearch->posResultsCounter] =
           search->result_struct->results[seqIndex]->subject_id;
         posSearch->posResultsCounter++;
         break;
       }
     }
   }
}

/*extract the e-value that applies to an entire dense
  diagonal alignment from its ScorePtr, based on similar
  code from Tom Madden*/

static Nlm_FloatHi getEvalueFromSeqAlign(SeqAlignPtr thisSeqAlign)
{
  ScorePtr thisScorePtr;

  thisScorePtr = thisSeqAlign->score;
  while ((thisScorePtr != NULL) &&
         (StringICmp(thisScorePtr->id->str, "e_value") != 0) &&
         (StringICmp(thisScorePtr->id->str, "sum_e") != 0))
    thisScorePtr = thisScorePtr->next;
  if(NULL == thisScorePtr)
    return(10.0);
  else
    return((Nlm_FloatHi) (thisScorePtr->value.realvalue));
}

/*Find the lowest e-value among all seqAligns for the sequence represented by
curSeqAlign*/
static Nlm_FloatHi minEvalueForSequence(SeqAlignPtr curSeqAlign, SeqAlignPtr listOfSeqAligns) 
{
   SeqAlignPtr testSeqAlign; /*Index into listOfSeqALigns*/
   DenseSegPtr curSegs, testSegs; /*Used to extract ids from curSeqAlign, testSeqAlign*/
   SeqIdPtr curId, testId; /*Ids of target sequences in testSeqAlign*/
   Nlm_FloatHi  returnValue; /*stores current best e-value*/
   Nlm_FloatHi  testEvalue; /*temporary e-value for one seqAlign*/
   Boolean seen;   /*have we seen a seqAlign matching the sequence yet*/

   returnValue = getEvalueFromSeqAlign(curSeqAlign);
   curSegs = (DenseSegPtr) curSeqAlign->segs;
   curId = curSegs->ids->next; 
   seen = FALSE;

   testSeqAlign = listOfSeqAligns;
   while (NULL != testSeqAlign) {
     testSegs = (DenseSegPtr) testSeqAlign->segs;
     testId = testSegs->ids->next; 
     if (SeqIdMatch(curId, testId)) {
         seen = TRUE;
	 if ((testEvalue = getEvalueFromSeqAlign(testSeqAlign)) < returnValue)
	   returnValue = testEvalue;
       }
     else
      /*if we have already seen a match and this one doesn't match,
        then stop looking*/
       if (seen)  
	 break;
     testSeqAlign = testSeqAlign->next;
   }
   return(returnValue);
}


/*Count the number of seqAligns in a list (returned) and count the number of
distinct target sequences represented (passed back in numSequences);
if useThreshold is TRUE, only those sequences with e-values below the threshold are counted.
Important assumption: All SeqAligns with  the same target sequence
are consecutive in the list*/
static Int4 countSeqAligns(SeqAlignPtr listOfSeqAligns, Int4 * numSequences, Boolean useThreshold, Nlm_FloatHi threshold)
{
   SeqAlignPtr curSeqAlign, prevSeqAlign;
   Int4 seqAlignCounter;
   DenseSegPtr curSegs;
   SeqIdPtr curId, prevId; /*Ids of target sequences in current and previous
                             SeqAlign*/

   seqAlignCounter = 0;
   *numSequences = 0;
   curSeqAlign = listOfSeqAligns;
   prevSeqAlign = NULL;
   while (NULL != curSeqAlign) {
     seqAlignCounter++;
     curSegs = (DenseSegPtr) curSeqAlign->segs;
     curId = curSegs->ids->next; 
     if ((NULL == prevSeqAlign) ||  (!(SeqIdMatch(curId, prevId))))
       if (!useThreshold || (threshold > minEvalueForSequence(curSeqAlign, listOfSeqAligns)))
	 (*numSequences)++;
     prevSeqAlign = curSeqAlign;
     prevId = curId;
     curSeqAlign = curSeqAlign->next;
   }
   return(seqAlignCounter);
}



/* Determines if  the search has converged from round to round.
*  Checks whether every new sequence found is in posSearch->posResultSequences.
*  Also sets up posSearch->posRepeatSequences, a boolean array that
*  indicates whether the sequence represented by the i-th new seqAlign is a repeat.
*  This is used in printing the sequences where they are
*  subdivided into two categories: sequences that were found previously
*  and new sequences.
*  posSearch is the data structure representing the parameters of the position-specific part
*  search represents the overall  BLAST search
*  listOfSeqAligns is one representation of the results of the current round.
*  If thissPassNum is 1, then it checks only to see that some sequence
* distinct from the query was found */
void LIBCALL posConvergenceTest(posSearchItems *posSearch, BlastSearchBlkPtr search, SeqAlignPtr listOfSeqAligns, Int4 thisPassNum)
{
  Int4 numseq;   /*Number of sequences found*/
  Int4 numalign; /* Number of items in listOfSeqAligns*/
  Int4 oldSeqIndex; /*Ordinal number of a sequence in old results (previous round)*/
  Int4 alignIndex; /*index into the list of seqAligns*/
  Boolean found;  /*Have we found the new sequence on the old list?*/
  SeqAlignPtr curSeqAlign, prevSeqAlign, startAlign; /* pointers into list of seqAligns*/
  DenseSegPtr curSegs;  /*Item in list of seqAligns*/
  SeqIdPtr thisId, prevId; /*Ids of target sequences in current and previous
			   SeqAlign*/
  Int4 ordinalNumber; /*Ordinal number of a sequence in the database*/
  Nlm_FloatHi thisEvalue; /*lowest evalue from all seqAligns for a sequence*/
  Int4 queryOffset, subjectOffset, retrievalOffset; /*offsets needed
                                                    to make a match align*/
  Int4 qplace, splace; /*index into query string and matching string*/
  Uint1Ptr q,s; /*Pointers into query and matching string*/
  Int4 queryLength; /*length of query*/
  Int4 matchLength; /* length of match*/
  Int4 subjectLength; /* length of a matching string*/
  Int4 c;  /*index into a string*/
  Int4 numsegs; /*number of segments in an alignment*/
  Int4 startQ,startS; /*Indices into array of starting positions*/


  numalign = countSeqAligns(listOfSeqAligns, &numseq, FALSE, 0.0);
  search->posConverged = TRUE;
  curSeqAlign = listOfSeqAligns;
  if (thisPassNum > 1) {
    posSearch->posRepeatSequences = (Int2Ptr) MemNew(numalign * sizeof(Int2));
    prevSeqAlign = NULL;
    for(alignIndex = 0; alignIndex < numalign; alignIndex++) {
      posSearch->posRepeatSequences[alignIndex] = (Int2) 0;
      curSegs = (DenseSegPtr) curSeqAlign->segs;
      thisId = curSegs->ids->next; 
      if ((NULL == prevSeqAlign) ||  (!(SeqIdMatch(thisId, prevId)))) {
	startAlign = curSeqAlign;
	thisEvalue = minEvalueForSequence(curSeqAlign, startAlign);
	if (thisEvalue < search->pbp->ethresh) {    
	  /*Extract the ordinal number from the SeqAlign*/
	  curSegs = (DenseSegPtr) curSeqAlign->segs;
	  thisId = curSegs->ids->next;  /*id of target sequence is second*/
	  /*Get ordinal ids of sequences in result*/
	  ordinalNumber = SeqId2OrdinalId(search->rdfp, thisId);
	  found = FALSE;
	  for(oldSeqIndex = 0; oldSeqIndex < posSearch->posResultsCounter; oldSeqIndex++)
	    if(ordinalNumber ==  posSearch->posResultSequences[oldSeqIndex]) {
	      posSearch->posRepeatSequences[alignIndex] = SEQ_ALIGN_MARK_REPEAT;
	      found = TRUE;
	      break;
	    }      
	  if (!found) 
	    search->posConverged = FALSE;
	}
      }
      else  /*both alignments come from the same sequence*/
	posSearch->posRepeatSequences[alignIndex] = posSearch->posRepeatSequences[alignIndex - 1];
      prevSeqAlign = curSeqAlign;
      prevId = thisId;
      curSeqAlign = curSeqAlign->next;
    }
    MemFree(posSearch->posResultSequences);
  }
  else {
    q = search->context[0].query->sequence;
    queryLength = search->context[0].query->length;
    prevSeqAlign = NULL;
    while (curSeqAlign != NULL) {
      curSegs = (DenseSegPtr) curSeqAlign->segs;
      s = GetSequenceWithDenseSeg(curSegs, FALSE, &retrievalOffset, &subjectLength);
      numsegs = curSegs->numseg;
      thisId = curSegs->ids->next; 
      if ((NULL == prevSeqAlign) ||  (!(SeqIdMatch(thisId, prevId)))) {
	startAlign = curSeqAlign;
	thisEvalue = minEvalueForSequence(curSeqAlign, startAlign);
	if (thisEvalue < search->pbp->ethresh) {    
	  if (numsegs > 1) {
	    search->posConverged = FALSE;
	    return;
	  }
	  startQ = 0;
	  startS = 1;
	  queryOffset = curSegs->starts[startQ];
	  if (curSegs->starts[startS] != GAP_HERE)
	    subjectOffset = curSegs->starts[startS] - retrievalOffset;
	  else
	    subjectOffset = GAP_HERE;
	  matchLength = curSegs->lens[0];
	  if ((queryOffset != 0) || (subjectOffset != 0) ||
	      (matchLength != queryLength) || (matchLength != subjectLength)) {
	    search->posConverged = FALSE;
	    return;
	  }
	  for (c = 0, qplace = queryOffset, splace = subjectOffset;
	       c < matchLength; c++, qplace++, splace++)
	    if (s[splace] != q[qplace]) {
	      search->posConverged = FALSE;
	      return;
	    }        
	}
      }
      prevSeqAlign = curSeqAlign;
      prevId = thisId;
      curSeqAlign = curSeqAlign->next;
    }
  }
}


/*Eliminate the matches from sequence second starting at position
matchStart and extending for intervalLength characters */
static void posCancel(posSearchItems *posSearch, compactSearchItems * compactSearch, Int4 first, Int4 second, Int4 matchStart, Int4 intervalLength)
{
  Int4 c, i;
  Boolean stillNeeded;

  for(c = matchStart, i = 0; i < intervalLength; i++, c++) {
    posSearch->posDescMatrix[second][c].used = FALSE;
    posSearch->posDescMatrix[second][c].letter = 0;
  }
  stillNeeded = FALSE;
  for (c = 0; c < compactSearch->qlength; c++)
    if (posSearch->posDescMatrix[second][c].used) {
      stillNeeded = TRUE;
      break;
    }
   if (!stillNeeded)
      posSearch->posUseSequences[second] = FALSE;
}

/*Eliminate sequences that are identical to the query and partial alignments
  that are identical in two matching sequences*/
static void  posPurgeMatches(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
  Int4 i, j; /*index over sequences*/
  Boolean matchesQuery; /*Is a matching sequence identical to the query?*/
  Int4 c; /*index over demographics of matching sequence*/
  Int4 state; /*state of checking for a match*/
  Int4 intervalLength, matchStart; /*Length and start of a matching region*/
  Int4 matchNumber; /*number of characters matching*/

  posSearch->posUseSequences =  (Boolean *) MemNew((posSearch->posNumSequences + 1) * sizeof(Boolean));
   if (NULL == posSearch->posUseSequences)
     exit(EXIT_FAILURE);
  for(i = 0; i <= posSearch->posNumSequences; i++)
    posSearch->posUseSequences[i] = TRUE;
  for(i = 1; i <= posSearch->posNumSequences; i++) {
    matchesQuery = TRUE;
    for (c = 0; c < compactSearch->qlength; c++) {
      if ((!posSearch->posDescMatrix[i][c].used) ||
          (posSearch->posDescMatrix[i][c].letter !=
           posSearch->posDescMatrix[0][c].letter)) {
        matchesQuery = FALSE;
        break;
      }
    }
    if (matchesQuery)
      posSearch->posUseSequences[i] = FALSE;
  }
  for(j = 1; j <= posSearch->posNumSequences; j++) {
    if (!posSearch->posUseSequences[j])
      continue;
    state = POS_COUNTING;
    c = 0;
    matchStart = 0;
    intervalLength = 0;
    matchNumber = 0;
    while (c < compactSearch->qlength) {
      if (posSearch->posDescMatrix[j][c].used) {
	if ((posSearch->posDescMatrix[0][c].letter != Xchar) &&
	    (posSearch->posDescMatrix[j][c].letter != Xchar)) { 
	  if (state == POS_RESTING) {
	    matchStart = c;
	    intervalLength = 1;
	    state = POS_COUNTING;
	    matchNumber = 0;
	  }
	  else 
	    intervalLength++;
	  if (posSearch->posDescMatrix[j][c].used &&
	      (posSearch->posDescMatrix[0][c].letter == posSearch->posDescMatrix[j][c].letter))
	    matchNumber++;
	  }
      }
      else {
	if (state == POS_COUNTING) {
	  if ((intervalLength > 0) && (matchNumber == intervalLength))
	    posCancel(posSearch,compactSearch,0,j,matchStart,intervalLength);
	  state = POS_RESTING;
	}
      }
      c++;
    }
    if (state == POS_COUNTING) /*at end of sequence i*/
      if ((intervalLength > 0) && (matchNumber == intervalLength))
	posCancel(posSearch,compactSearch,0,j,matchStart,intervalLength);
  }
  for (i = 1; i < posSearch->posNumSequences; i++) {
    if (!posSearch->posUseSequences[i])
      continue;
    for(j = i+1; j <= posSearch->posNumSequences; j++) {
      if (!posSearch->posUseSequences[j])
	continue;
      state = POS_COUNTING;
      c = 0;
      matchStart = 0;
      intervalLength = 0;
      matchNumber = 0;
      while (c < compactSearch->qlength) {
	if (posSearch->posDescMatrix[i][c].used ||
	    posSearch->posDescMatrix[j][c].used) {
	  if ((posSearch->posDescMatrix[i][c].letter != Xchar) &&
	      (posSearch->posDescMatrix[j][c].letter != Xchar)) { 
	    if (state == POS_RESTING) {
	      matchStart = c;
	      intervalLength = 1;
	      state = POS_COUNTING;
	      matchNumber = 0;
	    }
	    else 
	      intervalLength++;
	    if (posSearch->posDescMatrix[i][c].used &&
		posSearch->posDescMatrix[j][c].used &&
		(posSearch->posDescMatrix[i][c].letter == posSearch->posDescMatrix[j][c].letter))
	      matchNumber++;
	    }
	}
	else {
	  if (state == POS_COUNTING) {
	    if ((intervalLength > 0) && ((((Nlm_FloatHi) matchNumber)/intervalLength) >= IDENTITY_RATIO))
	      posCancel(posSearch,compactSearch,i,j,matchStart,intervalLength);
	    state = POS_RESTING;
	  }
	}
	c++;
      }
      if (state == POS_COUNTING) /*at end of sequence i*/
	if ((intervalLength > 0) && ((((Nlm_FloatHi) matchNumber)/intervalLength) >= IDENTITY_RATIO))
	  posCancel(posSearch,compactSearch,i,j,matchStart,intervalLength);
    }
  }
}

/*Compute general information about the sequences that matched on the
  i-th pass such as how many matched at each query position and what letter
  matched*/
static void posDemographics(posSearchItems *posSearch, compactSearchItems * compactSearch, SeqAlignPtr listOfSeqAligns)
{
   Uint1Ptr q; /*pointers into query */
   Uint1Ptr s; /*pointer into a matching string */
   Int4 length, subjectLength;  /*length of query and subject*/
   Int4 c; /*index into a string*/
   Int4 numseq, numSeqAligns;  /*number of matching sequences and SeqAligns*/
   Int4 seqIndex;  /*index for the array of matching sequences*/
   Int4 matchLength; /*length of a match*/
   Int4  queryOffset, subjectOffset, retrievalOffset;  /*offsets needed to make a match align*/
   Int4 qplace, splace; /*index into query string and matching string*/
   SeqAlignPtr curSeqAlign, prevSeqAlign; /*pointers into listOfSeqAligns*/
   DenseSegPtr curSegs, prevSegs;  /*used to extract alignments from curSeqAlign*/
   SeqIdPtr curId, prevId;  /*Used to compare sequences that come from different SeqAligns*/
   Int4 startQ, startS; /*Indices into array of starting positions*/
   Int4 numsegs; /*Number of pieces in the gapped alignment*/
   Int4 segIndex; /*Index for which piece we are at*/
   Nlm_FloatHi thisEvalue;  /*evalue of current partial alignment*/

   q = compactSearch->query;
   length = compactSearch->qlength;
   for(c = 0; c < length; c++) {
     posSearch->posDescMatrix[0][c].letter = (Int1) q[c];
     posSearch->posDescMatrix[0][c].used = TRUE;
     posSearch->posDescMatrix[0][c].leftExtent = 0;
     posSearch->posDescMatrix[0][c].rightExtent = length;
     posSearch->posDescMatrix[0][c].e_value = compactSearch->ethresh/2;
     posSearch->posC[c][q[c]]++;
     posSearch->posCount[c]++;
   }
   numSeqAligns = countSeqAligns(listOfSeqAligns, &numseq, TRUE, compactSearch->ethresh);
   posSearch->posNumSequences = numseq;
   /*use only those sequences below e-value threshold*/
   seqIndex = 0;
   curSeqAlign = listOfSeqAligns;
   prevSeqAlign = NULL;
   while (curSeqAlign != NULL) {
     if ((thisEvalue = getEvalueFromSeqAlign(curSeqAlign)) < compactSearch->ethresh) {
       curSegs = (DenseSegPtr) curSeqAlign->segs;
       if (NULL != prevSeqAlign) {
	 prevSegs = (DenseSegPtr) prevSeqAlign->segs;
         curId = curSegs->ids->next; 
         prevId = prevSegs->ids->next;
	 if (!(SeqIdMatch(curId, prevId)))
	   seqIndex++;
       }
       s = GetSequenceWithDenseSeg(curSegs, FALSE, &retrievalOffset, &subjectLength);
       startQ = 0;
       startS = 1;
       numsegs = curSegs->numseg;
       for(segIndex = 0; segIndex < numsegs; segIndex++) {
	 queryOffset = curSegs->starts[startQ];
         if (curSegs->starts[startS] != GAP_HERE)
	   subjectOffset = curSegs->starts[startS] - retrievalOffset;
	 else
           subjectOffset = GAP_HERE;
	 matchLength = curSegs->lens[segIndex];
	 if ((GAP_HERE ) == queryOffset) {
	   ; /*do nothing, gap in query*/
	 }
	 else
	   if ((GAP_HERE) == subjectOffset) {
	     for(c = 0, qplace = queryOffset;
		 c < matchLength; c++, qplace++) {
	       posSearch->posDescMatrix[seqIndex + 1][qplace].used = TRUE;
	       posSearch->posDescMatrix[seqIndex + 1][qplace].letter = GAP_CHAR;
	       posSearch->posDescMatrix[seqIndex + 1][qplace].e_value = 1.0;
	     }
	   }
	   else {  /*no gap*/
	     for(c = 0, qplace = queryOffset, splace = subjectOffset;
		 c < matchLength; c++, qplace++, splace++) {
	       if ((!posSearch->posDescMatrix[seqIndex+1][qplace].used) ||
		   (thisEvalue 
		    < posSearch->posDescMatrix[seqIndex+1][qplace].e_value)) 
		 if (!posSearch->posDescMatrix[seqIndex+1][qplace].used)
		   {
		     posSearch->posDescMatrix[seqIndex+1][qplace].letter = (Int1) s[splace]; 
		     posSearch->posDescMatrix[seqIndex+1][qplace].used = TRUE; 
		     posSearch->posDescMatrix[seqIndex+1][qplace].e_value = 
		       thisEvalue;
		   }
	     }
	   }
	 startQ += 2;
	 startS += 2;
       }
       prevSeqAlign = curSeqAlign;
       s = MemFree(s);
     }
     curSeqAlign = curSeqAlign->next;
   } /*closes the while loop over seqAligns*/
 }

static void posComputeExtents(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
   Int4 seqIndex; /*index of sequence*/
   Int4 length; /*length of query*/
   Int4 qplace, qplace2; /*place in query*/
   Int4 numseq; /*number of sequences including query*/
   Uint1Ptr q; /*pointers into query */

   length = compactSearch->qlength;
   numseq = posSearch->posNumSequences;
   q = compactSearch->query;
   for(seqIndex = 0; seqIndex < numseq; seqIndex++) {
     if (!posSearch->posUseSequences[seqIndex+1])
       continue;
     if ((posSearch->posDescMatrix[seqIndex+1][0].used)
	 && (posSearch->posDescMatrix[seqIndex+1][length-1].letter != GAP_CHAR))
       posSearch->posDescMatrix[seqIndex+1][0].leftExtent = 0;
     for(qplace = 1; qplace < length; qplace++)
       if(posSearch->posDescMatrix[seqIndex+1][qplace].used)
	 if(posSearch->posDescMatrix[seqIndex+1][qplace-1].used)
	   posSearch->posDescMatrix[seqIndex+1][qplace].leftExtent =
	     posSearch->posDescMatrix[seqIndex+1][qplace -1].leftExtent;
	 else
	   posSearch->posDescMatrix[seqIndex+1][qplace].leftExtent = qplace;
     if ((posSearch->posDescMatrix[seqIndex+1][length-1].used)
	 && (posSearch->posDescMatrix[seqIndex+1][length-1].letter != GAP_CHAR))
       posSearch->posDescMatrix[seqIndex+1][length-1].rightExtent = length -1;
     for(qplace = length -2; qplace >= 0; qplace--)
       if(posSearch->posDescMatrix[seqIndex+1][qplace].used)
	 if(posSearch->posDescMatrix[seqIndex+1][qplace+1].used)
	   posSearch->posDescMatrix[seqIndex+1][qplace].rightExtent =
	     posSearch->posDescMatrix[seqIndex+1][qplace + 1].rightExtent;
	 else
	   posSearch->posDescMatrix[seqIndex+1][qplace].rightExtent = qplace;
     for(qplace = 0; qplace < length; qplace++) 
       if (posSearch->posDescMatrix[seqIndex+1][qplace].used) {
	 posSearch->posExtents[qplace].leftExtent = MAX(posSearch->posExtents[qplace].leftExtent,
							posSearch->posDescMatrix[seqIndex+1][qplace].leftExtent);
	 posSearch->posExtents[qplace].rightExtent = MIN(posSearch->posExtents[qplace].rightExtent,
							 posSearch->posDescMatrix[seqIndex+1][qplace].rightExtent);
	 
       }

     for(qplace = 0; qplace < length; qplace++) 
     /*used to check qplace for GAP_CHAR here*/
       if (posSearch->posDescMatrix[seqIndex+1][qplace].used) {
	 posSearch->posC[qplace][posSearch->posDescMatrix[seqIndex+1][qplace].letter]++;
	 posSearch->posCount[qplace]++; /*Add to number of matches in this query position*/
       }
   }
   for(qplace = 0; qplace < length; qplace++)
     posSearch->posIntervalSizes[qplace] = posSearch->posExtents[qplace].rightExtent - 
       posSearch->posExtents[qplace].leftExtent + 1;
   for(qplace =0; qplace < length; qplace++) {
     if(Xchar == q[qplace]) {
       posSearch->posIntervalSizes[qplace] = 0;
       for(qplace2 = 0; qplace2 <qplace; qplace2++) {
	 if(posSearch->posExtents[qplace2].rightExtent >= qplace)
	   posSearch->posIntervalSizes[qplace2]--;
       }
     }
   }
 }
 
/*Compute weight of each sequence and letter in each position*/
static void posComputeSequenceWeights(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
   Int4 length; /*length of query*/
   Int4 numseq, seqIndex; /*number of matches, index for them*/
   Int4  i; /*index over a multi-alignment block*/
   Int4 qplace; /*index into query*/
   Nlm_FloatHi Sigma; /*Number of different characters occurring in matches within
                   a multi-alignment block, excluding identical columns*/
   Nlm_FloatHi intervalSigma; /*Same as Sigma but includes identical columns*/
   Int4 alphabetSize; /*number of characters in alphabet*/
   Int4 *participatingSequences; /*array of participating sequences at a position*/
   Int4 *oldParticipatingSequences; /*array of participating sequences at a position*/
   Int4 posLocalVariety;  /*number of different characters at a position*/
   Int4 *posLocalC; /*counts of how many of each letter in this column*/
   Int4 c;
   Int4 thisSeq;
   Int4 numParticipating; /*number of sequences in this alignment block*/
   Int4 oldNumParticipating; /*number of sequences in this alignment block*/
   Boolean newSequenceSet; 
   Int4 p; /*index on sequences*/

   alphabetSize = compactSearch->alphabetSize;
   length = compactSearch->qlength;
   numseq = posSearch->posNumSequences;
   participatingSequences = (Int4 *) MemNew((numseq+1) * sizeof(Int4));
   if (NULL == participatingSequences)
     exit(EXIT_FAILURE);
   oldParticipatingSequences = (Int4 *) MemNew((numseq+1) * sizeof(Int4));
   if (NULL == oldParticipatingSequences)
     exit(EXIT_FAILURE);
   posLocalC = (Int4 *) MemNew(alphabetSize * sizeof(Int4));
   if (NULL == posLocalC)
     exit(EXIT_FAILURE);
   for (qplace = 0; qplace < length; qplace++) {
     posSearch->posSigma[qplace] = 0.0;
   }
   numParticipating = 0;
   for(qplace = 0; qplace < length; qplace++) {
     if ((posSearch->posCount[qplace] > 1) && (posSearch->posIntervalSizes[qplace] > 0)) {
       oldNumParticipating = numParticipating;
       for(p =0; p < numParticipating; p++)
         oldParticipatingSequences[p] = participatingSequences[p];
       numParticipating = 0;
       for (seqIndex = 0; seqIndex <= numseq; seqIndex++) {
         if (!posSearch->posUseSequences[seqIndex])
           continue;
	 if ((posSearch->posDescMatrix[seqIndex][qplace].used) &&
	     (posSearch->posDescMatrix[seqIndex][qplace].letter != GAP_CHAR)) {
	   participatingSequences[numParticipating] = seqIndex;
	   numParticipating++;
	 }
       }
       newSequenceSet = TRUE;
       if (numParticipating == oldNumParticipating) {
         for(p = 0; p < numParticipating; p++)
           if (oldParticipatingSequences[p] != participatingSequences[p])
             break;
         if (p == numParticipating)
           newSequenceSet = FALSE;
       }
         
       if (newSequenceSet) {
	 Sigma = 0;
	 intervalSigma = 0;
	 for (seqIndex = 0; seqIndex <= numseq; seqIndex++) {
	   if (!posSearch->posUseSequences[seqIndex])
	     continue;
	   posSearch->posRowSigma[seqIndex] = 0.0;
	   posSearch->posA[seqIndex] = 0.0;
	 }
	 for (i = posSearch->posExtents[qplace].leftExtent;
	      i <= posSearch->posExtents[qplace].rightExtent; i++) {
	   posLocalVariety = 0;
	   for(c = 0; c < alphabetSize; c++)
	     posLocalC[c] = 0;
	   for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	     thisSeq = participatingSequences[seqIndex];
	     /*used to check for GAP here*/
	     if (0 == posLocalC[posSearch->posDescMatrix[thisSeq][i].letter])
	       /*letter (not a gap) not seen before in this query pos.*/
	       posLocalVariety++;  
	     posLocalC[posSearch->posDescMatrix[thisSeq][i].letter]++;
	   }
	   intervalSigma += posLocalVariety;
	   if (posLocalVariety > 1) {
	     Sigma += posLocalVariety;
	     for(seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	       thisSeq = participatingSequences[seqIndex];
	       /*used to check for gap here*/
	       posSearch->posRowSigma[thisSeq] += 
		 ( 1.0 / (Nlm_FloatHi) posLocalC[posSearch->posDescMatrix[thisSeq][i].letter]);
	     }
	   }
	 }
       }
       if (Sigma > 0) {
	 for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = posSearch->posRowSigma[thisSeq]/Sigma;
	 }
       }
       else {
         for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	   thisSeq = participatingSequences[seqIndex];
	   posSearch->posA[thisSeq] = ((Nlm_FloatHi) 1 / (Nlm_FloatHi) numParticipating);
         }
       }
       posSearch->posSigma[qplace] = intervalSigma;
       for (seqIndex = 0; seqIndex < numParticipating; seqIndex++) {
	 thisSeq = participatingSequences[seqIndex];
	 posSearch->posMatchWeights[qplace][posSearch->posDescMatrix[thisSeq][qplace].letter] += posSearch->posA[thisSeq];
       }
     }
   }
   MemFree(participatingSequences);
   MemFree(oldParticipatingSequences);
   MemFree(posLocalC);
}

static Nlm_FloatHi countsFunction(Nlm_FloatHi Sigma, Int4 intervalLength)
{
  return(Sigma / intervalLength - 1);
}

static Nlm_FloatHi posit_rounddown(Nlm_FloatHi value)
{
  return (Nlm_FloatHi) Nlm_Nint(value);
}

/*check that weights add to 1 in each column */
static void posCheckWeights(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
   Uint1Ptr q;  /*pointer to query*/
   Int4 length, alphabetSize; /*length of query and number of characters in alphabet*/
   Int4  a, c; /*loop indices*/
   Nlm_FloatHi runningSum; /*partial total for a column*/

   length = compactSearch->qlength;
   alphabetSize = compactSearch->alphabetSize;

   q = compactSearch->query;
   for(c = 0; c < length; c++) {
     if ((posSearch->posCount[c] > 1) && (q[c] != Xchar)) {
       runningSum = 0;
       for(a = 0; a < alphabetSize; a++) 
           runningSum += posSearch->posMatchWeights[c][a];
       if((runningSum < 0.99) || (runningSum > 1.01))
         ErrPostEx(SEV_ERROR, 0, 0, "\nERROR IN WEIGHTS, column %d, value %lf\n",c, runningSum);
     }
   }
}

/*Fill in information content per position from pseudo-count frequencies*/
static void  posFreqsToInformation(posSearchItems * posSearch, compactSearchItems * compactSearch)
{
   Int4 length;  /*length of the query*/
   Int4 c; /*loop index*/
   Int4 a, alphabetSize; /*loop index and size of alphabet*/
   Nlm_FloatHi  qOverPEstimate; /*intermediate term*/
   Nlm_FloatHi  infoSum; /*information content sum for this position*/
  
   length = compactSearch->qlength;
   alphabetSize = compactSearch->alphabetSize;
   for (c = 0; c < length; c++) {
     infoSum = 0;
     for(a = 0; a < alphabetSize; a++) {
       if (compactSearch->standardProb[a] > posEpsilon) {
         qOverPEstimate = posSearch->posFreqs[c][a] / compactSearch->standardProb[a];
         if (qOverPEstimate > posEpsilon)
	   infoSum += posSearch->posFreqs[c][a] * log(qOverPEstimate)/
                    NCBIMATH_LN2;
       }
     }
     posSearch->posInformation[c] = infoSum;
   }
}




/*Convert pseudo-count frequencies to a score matrix */
static void posFreqsToMatrix(posSearchItems *posSearch, compactSearchItems * compactSearch)
{
   Uint1Ptr q;  /*pointer to the query*/
   Int4 length;  /*length of the query*/
   Int4 c; /*loop index*/
   Int4 a, alphabetSize; /*loop index and size of alphabet*/
   Nlm_FloatHi lambda; /*Karlin-Altschul parameter*/
   Nlm_FloatHi  qOverPEstimate, value; /*intermediate terms*/
   Boolean allZeros; /*are all frequencies in a column 0?*/

   q = compactSearch->query;
   length = compactSearch->qlength;

   alphabetSize = compactSearch->alphabetSize;
   lambda = compactSearch->lambda_ideal;

   for(c = 0; c < length; c++) {
     allZeros = TRUE;
     for(a = 0; a < alphabetSize; a++) {
       /*Division compensates for multiplication in posComputePsedoFreqs*/
       if (compactSearch->standardProb[a] > posEpsilon)
	 qOverPEstimate = posSearch->posFreqs[c][a]/compactSearch->standardProb[a];
       else
        qOverPEstimate = 0.0;
       if (qOverPEstimate != 0.0)
         allZeros = FALSE;
       if (0.0 == qOverPEstimate || (compactSearch->standardProb[a] < posEpsilon))
	 posSearch->posPrivateMatrix[c][a] = BLAST_SCORE_MIN;
       else {
	 value = log(qOverPEstimate)/lambda;
	 posSearch->posPrivateMatrix[c][a] = (BLAST_Score) posit_rounddown(POSIT_SCALE_FACTOR*value);
       }
     }    
     if (allZeros) {
       for(a = 0; a < alphabetSize; a++) {
         posSearch->posMatrix[c][a] = compactSearch->matrix[q[c]][a];
	 if (compactSearch->matrix[q[c]][a] == BLAST_SCORE_MIN)
		posSearch->posPrivateMatrix[c][a] = BLAST_SCORE_MIN;
	 else
         	posSearch->posPrivateMatrix[c][a] = POSIT_SCALE_FACTOR*compactSearch->matrix[q[c]][a];
       }
     }
   }
   for(a = 0; a < alphabetSize; a++) {
     posSearch->posPrivateMatrix[length][a] = posSearch->posPrivateMatrix[length][a] = BLAST_SCORE_MIN;
   }
}

static Nlm_FloatHi ** allocatePosFreqs(Int4 length, Int4 alphabetSize)
{
  Int4 c, i; /*loop indices*/
  Nlm_FloatHi ** returnArray;

  returnArray = (Nlm_FloatHi **) MemNew((length + 1) * sizeof(Nlm_FloatHi *));
  if (NULL == returnArray)
    exit(EXIT_FAILURE);
  for(i = 0; i <= length; i++) {
    returnArray[i] = (Nlm_FloatHi *) MemNew(alphabetSize * sizeof(Nlm_FloatHi));
    if (NULL == returnArray[i])
      exit(EXIT_FAILURE);   
    for(c = 0; c < alphabetSize; c++)
      returnArray[i][c] = 0.0;
  }
  return(returnArray); 
}

static Nlm_FloatHi ** posComputePseudoFreqs(posSearchItems *posSearch, compactSearchItems * compactSearch, Boolean Cpos)
{
   Uint1Ptr q;  /*pointer to the query*/
   Int4 length;  /*length of the query*/
   Int4 c; /*loop index*/
   Int4 a, aSub, alphabetSize; /*loop indices and size of alphabet*/
   Nlm_FloatHi lambda; /*Karlin-Altschul parameter*/
   Nlm_FloatHi Sigma;  /*number of characters in an interval*/ 
   Int4 intervalLength;  /*length of a block*/
   Nlm_FloatHi pseudo, numerator, denominator, qOverPEstimate; /*intermediate terms*/
   Nlm_FloatHi infoSum; /*sum used for information content*/
   Nlm_FloatHi **posFreqs; /*store frequencies*/

   q = compactSearch->query;
   length = compactSearch->qlength;

   alphabetSize = compactSearch->alphabetSize;
   lambda = compactSearch->lambda_ideal;
   posFreqs = allocatePosFreqs(length, alphabetSize);

   for(c = 0; c < length; c++) {
     if ((posSearch->posCount[c] > 1) && (Xchar != q[c])) {
       infoSum = 0;
       for(a = 0; a < alphabetSize; a++) {
         if (compactSearch->standardProb[a] > posEpsilon) {
	   pseudo = 0;
	   for (aSub = 0; aSub < alphabetSize; aSub++)
	     if(compactSearch->matrix[a][aSub] != BLAST_SCORE_MIN) 
	       pseudo += (posSearch->posMatchWeights[c][aSub] *
			exp(lambda * compactSearch->matrix[a][aSub]));
	   pseudo *= (compactSearch->pseudoCountConst);
           Sigma = posSearch->posSigma[c];
           intervalLength = posSearch->posIntervalSizes[c];
	   numerator = pseudo + 
             (countsFunction(Sigma, intervalLength) * posSearch->posMatchWeights[c][a]/
                compactSearch->standardProb[a]);
	   denominator = countsFunction(Sigma, intervalLength) + (compactSearch->pseudoCountConst); 
	   qOverPEstimate = numerator / denominator;
	   /*Note artificial multiplication by standard probability to
             normalize*/
           posFreqs[c][a] = qOverPEstimate * compactSearch->standardProb[a];
	 if (0.0 != qOverPEstimate && (compactSearch->standardProb[a] > posEpsilon))
	   infoSum += qOverPEstimate * compactSearch->standardProb[a] * log(qOverPEstimate)/ NCBIMATH_LN2;
	 }
        else
          posFreqs[c][a] = 0.0;
       }
       if (Cpos)
	 posSearch->posInformation[c] = infoSum;
     }
     else
       for(a = 0; a < alphabetSize; a++) {
         posFreqs[c][a] = 0;
       }
   }
  return(posFreqs);
}

static void posScaling(posSearchItems *posSearch, compactSearchItems * compactSearch)
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

	BlastScaleMatrix(matrix_rescale, TRUE);

	matrix_rescale = BlastMatrixRescaleDestruct(matrix_rescale);

	return;
}


Int4Ptr * LIBCALL CposComputation(posSearchItems *posSearch, BlastSearchBlkPtr search, compactSearchItems * compactSearch, SeqAlignPtr listOfSeqAligns, Char *ckptFileName, Boolean patternSearchStart, ValNodePtr * error_return)
{
  search->posConverged = FALSE;
  if (patternSearchStart)
    posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, posSearch->posNumSequences);
  else
    posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, search->result_struct->hitlist_count);
  if (!patternSearchStart)
    findThreshSequences(posSearch, search);
  posDemographics(posSearch, compactSearch, listOfSeqAligns);
  posPurgeMatches(posSearch, compactSearch);
  posComputeExtents(posSearch, compactSearch);
  posComputeSequenceWeights(posSearch, compactSearch);
  posCheckWeights(posSearch, compactSearch);
  posSearch->posFreqs = posComputePseudoFreqs(posSearch, compactSearch, TRUE);
  if (NULL != ckptFileName)
    posTakeCheckpoint(posSearch, compactSearch, ckptFileName, error_return);
  posFreqsToMatrix(posSearch,compactSearch);
  posScaling(posSearch, compactSearch);
  return posSearch->posMatrix;
}

/* Top-level routine to compute position-specific matrix, when used through
the Web, one round at a time*/
Int4Ptr * LIBCALL WposComputation(compactSearchItems * compactSearch, SeqAlignPtr listOfSeqAligns)
{
  posSearchItems *posSearch;
  Int4 numSeqAligns, numseq;
  Int4Ptr *posMatrix;

  posSearch = (posSearchItems *) MemNew(1 * sizeof(posSearchItems));
  numSeqAligns = countSeqAligns(listOfSeqAligns, &numseq, FALSE, 0.0);
  posAllocateMemory(posSearch, compactSearch->alphabetSize, compactSearch->qlength, numseq);
  posDemographics(posSearch, compactSearch, listOfSeqAligns);
  posPurgeMatches(posSearch, compactSearch);
  posComputeExtents(posSearch, compactSearch);
  posComputeSequenceWeights(posSearch, compactSearch);
  posCheckWeights(posSearch, compactSearch);
  posSearch->posFreqs = posComputePseudoFreqs(posSearch, compactSearch, FALSE);
  posFreqsToMatrix(posSearch,compactSearch);
  posScaling(posSearch, compactSearch);
  posMatrix = posSearch->posMatrix;
  MemFree(posSearch);
  return posMatrix;
}


static Char getRes(Char input)
{
    switch(input) 
      {
      case 0: 
	return('-');
      case 1: 
	return('A');
      case 2: 
	return('B');
      case 3: 
	return('C');
      case 4: 
	return('D');
      case 5: 
	return('E');
      case 6: 
	return('F');
      case 7: 
	return('G');
      case 8: 
	return('H');
      case 9: 
	return('I');
      case 10: 
	return('K');
      case 11: 
	return('L');
      case 12: 
	return('M');
      case 13: 
	return('N');
      case 14: 
	return('P');
      case 15: 
	return('Q');
      case 16: 
	return('R');
      case 17: 
	return('S');
      case 18: 
	return('T');
      case 19: 
	return('V');
      case 20: 
	return('W');
      case 21: 
	return('X');
      case 22: 
	return('Y');
      case 23: 
	return('Z');
      case 24: 
	return('U');
      case 25: 
	return('*');
      default:
        return('?');
    }
}
static Uint1 ResToInt(Char input)
{
    switch(input) 
      {
      case '-': 
        return(0);
      case 'A': 
        return(1);
      case 'B': 
        return(2);
      case 'C': 
        return(3);
      case 'D': 
        return(4);
      case 'E': 
        return(5);
      case 'F': 
        return(6);
      case 'G': 
        return(7);
      case 'H': 
        return(8); 
      case 'I': 
        return(9);
      case 'K': 
        return(10);
      case 'L': 
        return(11);
      case 'M': 
        return(12);
      case 'N': 
        return(13);
      case 'P': 
        return(14);
      case 'Q': 
        return(15);
      case 'R': 
        return(16);
      case 'S': 
        return(17);
      case 'T': 
        return(18);
      case 'V': 
        return(19);
      case 'W': 
        return(20);
      case 'X': 
        return(21);
      case 'Y': 
        return(22);
      case 'Z': 
        return(23);
      case 'U': 
        return(24);
      case '*': 
        return(25);
      default:
        return(-1);
    }
}


/*Print out the position-specific matrix*/
void LIBCALL outputPosMatrix(posSearchItems *posSearch, compactSearchItems *compactSearch)
{
   Uint1Ptr q; /*query sequence*/
   Int4 i, index; /*loop indices*/
   Int4 c; /*index over alphabet*/
   Int4 length; /*length of query*/
   Int4 charOrder[EFFECTIVE_ALPHABET]; /*standard order of letters according to S. Altschul*/

   if (compactSearch->alphabetSize != PROTEIN_ALPHABET){
     ErrPostEx(SEV_ERROR, 0, 0, "\nCannot print diagnostic information because alphabet size is not %ld", (long) compactSearch->alphabetSize);
     return;
   }
   
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

   q = compactSearch->query;
   length = compactSearch->qlength;
   index = 0;   
   
/* Used ifdef until final decision is made on output. */
#ifdef POSIT_DEBUG

   printf("\nCharacter Frequencies by positon\n");
   printf("         ");
   for (c = 0; c< EFFECTIVE_ALPHABET; c++)
      printf("  %c",getRes((Char) charOrder[c]));
   for(i=0; i < length; i++) {
     printf("\n%5d %c   ", i + 1, getRes(q[i]));
     for (c = 0; c < EFFECTIVE_ALPHABET; c++) 
       printf("%2d ", posSearch->posC[i][charOrder[c]]);
   }
   printf("\n\n");
   printf("\nposition counts used. multiplied by 10 and rounded and");
   printf("\nposition character weights used, multiplied by 10 and rounded\n");
   printf("        Counts");
   for (c = 0; c< EFFECTIVE_ALPHABET; c++)
      printf("  %c",getRes((Char) charOrder[c]));
   printf(" Extent ");
   for(i=0; i < length; i++) {
     printf("\n%5d %c   ", i + 1, getRes(q[i]));
     if ((posSearch->posCount[index+i] > 1) && (Xchar != q[index+i]))
       printf("%4d ", (Int4) posit_rounddown(10 * countsFunction
				     (posSearch->posSigma[i],posSearch->posIntervalSizes[i])));
     else
       printf("     ");
     for (c = 0; c< EFFECTIVE_ALPHABET; c++)
       if((posSearch->posMatrix[i][charOrder[c]] == BLAST_SCORE_MIN) ||
             (0.0 == posSearch->posMatchWeights[i][charOrder[c]]))
           printf(" - ");
         else
	   printf("%2d ", (Int4) posit_rounddown(10 * posSearch->posMatchWeights[i][charOrder[c]]));
     printf(" %4d",posSearch->posExtents[i].rightExtent - posSearch->posExtents[i].leftExtent +1);
   }
   printf("\n\n");
   printf("\nPosition-based scoring matrix used\n");
   printf("         ");
   for (c = 0; c< EFFECTIVE_ALPHABET; c++)
      printf("  %c",getRes((Char) charOrder[c]));
   for(i=0; i < length; i++) {
    /* printf("\n%5d %c   ", i + 1, getRes(q[i])); */
     printf("\n          ");
     for (c = 0; c < EFFECTIVE_ALPHABET; c++) 
       if(posSearch->posMatrix[i][charOrder[c]] == BLAST_SCORE_MIN)
	 printf("-I ");
       else
	 printf("%2d ", (Int4) posSearch->posMatrix[i][charOrder[c]]);
   }
   printf("\n\n");
   printf("                      K         Lambda\n");
   printf("Standard Ungapped    %6.4f     %6.4f\n",compactSearch->kbp_std[0]->K,compactSearch->kbp_std[0]->Lambda);
   printf("Standard Gapped      %6.4f     %6.4f\n",compactSearch->kbp_gap_std[0]->K,compactSearch->kbp_gap_std[0]->Lambda);
   printf("PSI Ungapped         %6.4f     %6.4f\n",compactSearch->kbp_psi[0]->K,compactSearch->kbp_psi[0]->Lambda);
   printf("PSI Gapped           %6.4f     %6.4f\n",compactSearch->kbp_gap_psi[0]->K,compactSearch->kbp_gap_psi[0]->Lambda);
   printf("\n\n");
#endif
}


void LIBCALL posPrintInformation(posSearchItems *posSearch, BlastSearchBlkPtr search, Int4 passNum)
{
  Int4 querySize;
  Int4 c;

  querySize = search->context[0].query->length;

/* Used ifdef until final decision is made on output. */
#ifdef POSIT_DEBUG
  printf("\nInformation content by position for pass %d\n", passNum);
  for(c = 0; c < querySize; c++)
    printf(" %5d", c); 
  printf("\n");
  for(c = 0; c < querySize; c++)
    printf(" %5.2lf", posSearch->posInformation[c]); 
  printf("\n");
#endif
}   
 
void LIBCALL posInitializeInformation(posSearchItems *posSearch, BlastSearchBlkPtr search)
{
  Uint1Ptr query;
  Int4 querySize;
  Int4 c, a, alphabetSize;
  BLAST_ScoreBlkPtr sbp;
  BLAST_ResFreqPtr stdrfp; /*standard frequencies*/
  Nlm_FloatHi lambda;
  Nlm_FloatHi term1, term2, term3, term4;
  Nlm_FloatHi infoSum;
 
  querySize = search->context[0].query->length;
  query = search->context[0].query->sequence;
  posSearch->posInformation = (Nlm_FloatHi *) MemNew(querySize * sizeof(Nlm_FloatHi));
  if (NULL == posSearch->posInformation)
    exit(EXIT_FAILURE);
  for(c = 0; c < querySize; c++)
    posSearch->posInformation[c] = 0.0;
  alphabetSize = search->sbp->alphabet_size;
  /*Compute standard frequencies as in BlastScoreBlkFill in blastkar.c*/
  sbp = search->sbp;
  stdrfp = BlastResFreqNew(sbp);
  BlastResFreqStdComp(sbp,stdrfp); 
  lambda = search->sbp->kbp[0]->Lambda;
  for(c = 0; c < querySize; c++) {
    infoSum = 0;
    for(a = 0; a < alphabetSize; a++)
      if (stdrfp->prob[a] > posEpsilon) {
        term1 = search->sbp->matrix[query[c]][a];
	term2 = term1 * lambda;
	term3 = exp(term2);
	term4 = stdrfp->prob[a] * term3;
	infoSum += term4 * log(term4/stdrfp->prob[a])/NCBIMATH_LN2;
      }
    posSearch->posInformation[c] = infoSum;
  }
}

/*
	Is this function used?
*/

void LIBCALL posFreeInformation(posSearchItems *posSearch)
{
  MemFree(posSearch->posInformation);
}

/*Copy a few fields from the lasrge record search into the small record
  compactSearch, so that a small amount of information
  is passed into posit.c*/
void LIBCALL copySearchItems(compactSearchItems * compactSearch, BlastSearchBlkPtr search)
{
   BLAST_ResFreqPtr stdrfp; /* gets standard frequencies in prob field */
   Int4 a; /*index over characters*/

   compactSearch->query = search->context[0].query->sequence;
   compactSearch->qlength = search->context[0].query->length;
   compactSearch->gapped_calculation = search->pbp->gapped_calculation;
   compactSearch->alphabetSize = search->sbp->alphabet_size;
   compactSearch->pseudoCountConst = search->pbp->pseudoCountConst;
   compactSearch->ethresh = search->pbp->ethresh;
   compactSearch->lambda =  search->sbp->kbp[0]->Lambda;
   compactSearch->matrix = search->sbp->matrix;
   compactSearch->kbp_psi = search->sbp->kbp_psi;
   compactSearch->kbp_gap_psi = search->sbp->kbp_gap_psi;
   compactSearch->kbp_std = search->sbp->kbp_std;
   compactSearch->kbp_gap_std = search->sbp->kbp_gap_std;
   compactSearch->lambda_ideal = search->sbp->kbp_ideal->Lambda;
   compactSearch->K_ideal = search->sbp->kbp_ideal->K;

   stdrfp = BlastResFreqNew(search->sbp);
   BlastResFreqStdComp(search->sbp,stdrfp); 
   compactSearch->standardProb = MemNew(compactSearch->alphabetSize * sizeof(Nlm_FloatHi));
   if (NULL == compactSearch->standardProb)
     exit(EXIT_FAILURE);
   for(a = 0; a < compactSearch->alphabetSize; a++)
     compactSearch->standardProb[a] = stdrfp->prob[a];
   stdrfp = BlastResFreqDestruct(stdrfp);
}

/*allocate memory for a record of type compactSearchItems*/
compactSearchItems * LIBCALL  compactSearchNew(compactSearchItems * compactSearch)
{
   compactSearch = MemNew(1 * sizeof(compactSearchItems));
   if (NULL == compactSearch)
     exit(EXIT_FAILURE);
   return(compactSearch);
}

/*De-allocate memory for a record of type compactSearchItems*/
void LIBCALL compactSearchDestruct(compactSearchItems * compactSearch)
{

   MemFree(compactSearch->standardProb);
   MemFree(compactSearch);
}

/*Some of the following checkpointing code is taken and adapted from
code written by K. Shriram for FASTLINK.
Reference:
 A. A. Schaffer, S. K. Gupta, K. Shriram, and R. W. Cottingham, Jr. 
 Avoiding Recomputation in Linkage Analysis,
 Human Heredity 44(1994), pp. 225-237. */


#define  putCkptNlm_FloatHi(d, ckptFile)  (putCkptNumber(&(d),sizeof(Nlm_FloatHi),ckptFile))
#define  putCkptInt4(i, ckptFile)         (putCkptNumber(&(i),sizeof(Int4),ckptFile))
#define  putCkptChar(c, ckptFile)         (putCkptNumber(&(c),sizeof(Char),ckptFile))
 
/* General routine for putting the internal representation of a number. */
 
static void  putCkptNumber(void * numberPtr, Int4 numberSize, FILE * ckptFile )
{
  FileWrite(numberPtr,numberSize,1,ckptFile) ;
}

static void    putCkptNlmFloat_HiVector (Nlm_FloatHi * theVector, Int4 length, FILE * ckptFile)
{
  int  vectorRef;
 
  for(vectorRef = 0; vectorRef < length; vectorRef++)
    putCkptNlm_FloatHi(theVector[vectorRef],ckptFile) ;
}


/*Code to put a vector of frequencies; put only the interesting
  entries*/
static void  putFreqVector(Nlm_FloatHi * theVector, Int4 length, FILE * ckptFile)
{
   int  vectorRef;
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

 
   for(vectorRef = 0; vectorRef < EFFECTIVE_ALPHABET; vectorRef++)
     putCkptNlm_FloatHi(theVector[charOrder[vectorRef]],ckptFile);
}

 
/* Code to put a matrix, vector-by-vector. */
static void    putCkptFreqMatrix (Nlm_FloatHi **theMatrix, Int4 length, Int4 width, FILE * ckptFile)
{
  int  matrixRef;  /*loop index*/
 
  for (matrixRef = 0; matrixRef < length ; matrixRef++ )
    putFreqVector(theMatrix[matrixRef], width, ckptFile);
}

 
/* Front-ends to retrieve numbers. */

#define  getCkptNlm_FloatHi(d, ckptFile)  (getCkptNumber(&(d),sizeof(Nlm_FloatHi),ckptFile))
#define  getCkptInt4(i, ckptFile)         (getCkptNumber(&(i),sizeof(Int4),ckptFile))
#define  getCkptChar(c, ckptFile)         (getCkptNumber(&(c),sizeof(Char),ckptFile))
 
/* General routine for getting the internal representation of a number. */
 
static void  getCkptNumber(void * numberPtr, Int4 numberSize, FILE * ckptFile )
{
  FileRead(numberPtr,numberSize,1,ckptFile) ;
}

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

/* Code to frequency matrix, vector-by-vector. */
 
static void    getCkptFreqMatrix (Nlm_FloatHi ** theMatrix, Int4 length, Int4 width, FILE * ckptFile)
{
  Int4  matrixRef;  /*loop index*/
 
  for (matrixRef = 0; matrixRef < length ; matrixRef++ )
    getFreqVector(theMatrix[matrixRef], width, ckptFile);
}

/*Take a checkpoint at the end of the current PSI-BLAST round, stores
query length, query, and position-specific target frequencies.
Returns TRUE if checkpoint was sucessful and FALSE otherwise. */
Boolean LIBCALL posTakeCheckpoint(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, ValNodePtr *error_return)
{
  FILE * checkFile; /*file in which to take the checkpoint*/
  Int4 length; /*length of query sequence, and an index for it*/
  Int4 i; /*indices to position and alphabet */
  Char localChar; /*temporary character*/

  checkFile = FileOpen(fileName, "w");
  if (NULL == checkFile) {
    BlastConstructErrorMessage("posTakeCheckpoint", "Could not open checkpoint file", 1, error_return);
    return(FALSE);
  }
  length = compactSearch->qlength;
  putCkptInt4(length, checkFile);
  for(i = 0; i < length; i++) {
    localChar = getRes(compactSearch->query[i]);
    putCkptChar(localChar, checkFile);
  }  
  putCkptFreqMatrix(posSearch->posFreqs,length,compactSearch->alphabetSize, checkFile);
  FileClose(checkFile);
  return(TRUE);
}

/*Read a checkpoint from the end of a previous PSI-BLAST round, get
query length, query, and position-specific target frequencies.
Returns TRUE if checkpoint was read sucessfully and FALSE otherwise. */
Boolean LIBCALL  posReadCheckpoint(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, ValNodePtr * error_return)
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
  posFreqsToInformation(posSearch,compactSearch);
  posFreqsToMatrix(posSearch,compactSearch);
  posScaling(posSearch, compactSearch);
  BlastConstructErrorMessage("posReadCheckpoint", "Data recovered successfully\n", 1, error_return);
  MemFree(oldQuery);
  FileClose(checkFile);
  return(TRUE);
}

