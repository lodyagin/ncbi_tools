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

File name: posit.h

Author: Alejandro Schaffer

Contents: header file for position-based BLAST.


*****************************************************************************/
/* $Revision: 6.5 $ */
/* $Log: posit.h,v $
/* Revision 6.5  1998/09/28 12:31:32  madden
/* Used BlastConstructErrorMessage
/*
 * Revision 6.3  1998/04/24 19:29:50  madden
 * Added ideal values to compactSearch
 *
 * Revision 6.2  1998/03/25 22:36:19  egorov
 * Change type of posRepeatSequences
 *
 * Revision 6.1  1997/12/23 21:07:08  madden
 * Changes for checkpointing
 *
 * Revision 6.0  1997/08/25 18:53:52  madden
 * Revision changed to 6.0
 *
 * Revision 1.11  1997/08/11 15:45:28  madden
 * eliminated obsolete fields
 *
 * Revision 1.10  1997/06/25 14:04:51  madden
 * prototype change
 *
 * Revision 1.9  1997/05/29 20:36:23  madden
 * Add Boolean *posUseSequences
 *
 * Revision 1.8  1997/05/22 21:25:30  madden
 * fixed memory leaks
 *
 * Revision 1.7  1997/05/16 20:10:10  madden
 * Added BLAST_Score **posPrivateMatrix
 *
 * Revision 1.6  1997/05/01 15:53:27  madden
 * Addition of extra KarlinBlk's for psi-blast
 *
 * Revision 1.5  1997/04/22  16:36:49  madden
 * Changes for use of psi-blast with www.
 *
 * Revision 1.4  1997/04/10  19:25:53  madden
 * COMMAND_LINE replaced by ALL_ROUNDS, Char to Int1.
 *
 * Revision 1.3  1997/04/09  20:01:53  madden
 * Functions CposComputation and WposComputation replace posComputations.
 *
 * Revision 1.2  1997/04/04  20:44:55  madden
 * Changed posComputation to return Int4Ptr *.
 *
 * Revision 1.1  1997/02/13  15:22:13  madden
 * Initial revision
 *
*/

#ifndef __POSIT__
#define __POSIT__

#ifdef __cplusplus
extern "C" {
#endif


#include <ncbi.h>
#include <math.h>
#include <blast.h>
#include <blastdef.h>

#define charsPerLine 20 /*Number of characters of a sequence to print
                          per line for score matrix*/

#define UNUSED (-1)

#define Xchar   21    /*character for low-complexity columns*/

#define ALL_ROUNDS 1 /*do all rounds without interruption*/

typedef struct posDesc {
  Int1 letter;  /*what is the preferred letter here*/
  Boolean used;  /*is there any letter here */
  Nlm_FloatHi e_value; /*score of highest hsp including this position */
  Int4 leftExtent; /*How far left do same sequences match?*/
  Int4 rightExtent; /*How far right do same sequences match?*/
} posDesc;

typedef struct posSearchItems {
  Int4 *posCount; /*count of how many sequences match at
                  each query position, default value is 1 to
                  include query*/
  Int4 **posC; /*position-sepcific occurrence counts*/
  Nlm_FloatHi **posMatchWeights;
  BLAST_Score **posMatrix;
  BLAST_Score **posPrivateMatrix;
  Nlm_FloatHi **posFreqs;
  Int4 *threshSequences;   /*Which sequences are below p-value threshold*/
  Int4 posMaxThresh;  /*Highest index of a sequence below p-value threshold*/
  Int4 posNumSequences;
  Int4 posResultsCounter;
  Int4 *posResultSequences;
  Nlm_FloatHi *posA;
  Nlm_FloatHi *posRowSigma;
  Int4 posDescMatrixLength;	/* Length of posDescMatrix, for deallocation. */
  posDesc **posDescMatrix;
  posDesc *posExtents;
  Nlm_FloatHi *posSigma;
  Int4 *posIntervalSizes;  /*interval size used for this column*/
  Int2Ptr posRepeatSequences;
  Boolean *posUseSequences;
  Nlm_FloatHi *posInformation;
} posSearchItems;

typedef struct compactSearchItems {
  Uint1Ptr  query;
  Int4 qlength;
  Boolean gapped_calculation;
  Int4 alphabetSize;
  Int4 pseudoCountConst;
  Nlm_FloatHi ethresh;
  Nlm_FloatHi lambda;
  Nlm_FloatHi *standardProb;
  Int4Ptr  *matrix;
  BLAST_KarlinBlkPtr *kbp_std, *kbp_psi, *kbp_gap_std, *kbp_gap_psi;
  Nlm_FloatHi	lambda_ideal,
		K_ideal;
} compactSearchItems;
  

void LIBCALL outputPosMatrix PROTO((posSearchItems *posSearch, compactSearchItems * compactSearch));

Int4Ptr * LIBCALL CposComputation PROTO((posSearchItems *posSearch, BlastSearchBlkPtr search, compactSearchItems * compactSearch, SeqAlignPtr listOfSeqAligns, Char *ckptFileName, Boolean patternSearchStart, ValNodePtr * error_return));

Int4Ptr * LIBCALL WposComputation PROTO((compactSearchItems *compactSearch, SeqAlignPtr listOfSeqAligns));

void LIBCALL posPrintInformation PROTO((posSearchItems *posSearch, BlastSearchBlkPtr search, Int4 passNum));

void LIBCALL posInitializeInformation PROTO((posSearchItems *posSearch, BlastSearchBlkPtr search));

void LIBCALL posFreeInformation PROTO((posSearchItems *posSearch));

void LIBCALL posConvergenceTest PROTO((posSearchItems *posSearch, BlastSearchBlkPtr search, SeqAlignPtr listOfSeqAligns, Int4 thisPassNum));

/*Cleanup position-specific  data structures after one pass*/
void LIBCALL posCleanup PROTO((posSearchItems *posSearch, compactSearchItems * compactSearch));

void LIBCALL copySearchItems(compactSearchItems * compactSearch, BlastSearchBlkPtr search);

compactSearchItems * LIBCALL compactSearchNew(compactSearchItems * compactSearch);

void LIBCALL compactSearchDestruct(compactSearchItems * compactSearch);

Boolean LIBCALL posTakeCheckpoint(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, ValNodePtr * error_return);

Boolean LIBCALL posReadCheckpoint(posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, ValNodePtr * error_return);

void LIBCALL posCheckpointFreeMemory(posSearchItems *posSearch, Int4 querySize);

#ifdef __cplusplus

}
#endif

#endif /* __POSIT__ */