/**************************************************************************
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
*                                                                         *
**************************************************************************/

/**************************************************************************
File name: profiles.c

Author: Alejandro Schaffer

Contents: utilities for IMPALA

**************************************************************************/



#define MAXLINELEN 2000
#define MAX_NAME_LENGTH 100

#define   ARG_DB 0
#define   ARG_QUERY_FILE 1
#define   ARG_E_VALUE 2
#define   ARG_ALIGN_VIEW 3
#define   ARG_OUTPUT_FILE 4
#define   ARG_FILTER 5
#define   ARG_GAP_OPEN 6
#define   ARG_GAP_EXTEND 7
#define   ARG_X_DROP  8
#define   ARG_PROCESSORS 9
#define   ARG_SHOW_GI 10
#define   ARG_E_MULTIPASS 11
#define   ARG_PSEUDO 12
#define   ARG_NUM_PASS 13
#define   ARG_BELIEVE_DEF 14
#define   ARG_SEQALIGN_FILE 15
#define   ARG_MATRIX 16
#define   ARG_MATRICES_DB 17
#define   ARG_DB_LENGTH  18


#define PRO_VERSION "0.2"
#define PRO_DATE "11-Jan-1998"
#define PRO_MAX_HIT_LIST 250
#define PRO_NUM_TICKS 50

#define PRO_ALPHABET_SIZE  26

#define PRO_DEFAULT_SCALING_UP  100
#define PRO_DEFAULT_SCALING_DOWN 0.01

#define scoreRange 10000

/*factor used to multiply the gapped K parameter to make it
  more accurate in most cases*/
#define PRO_K_MULTIPLIER 1.2

typedef BLAST_Score ScoreRow[PRO_ALPHABET_SIZE];

typedef struct SWpairs {
  Int4 noGap;
  Int4 gapExists;
} SWpairs;

typedef struct SWResults {
  Uint1Ptr seq;
  Int4 seqStart;
  Int4 seqEnd;
  Int4 queryStart;
  Int4 queryEnd;
  Int4 *reverseAlignScript;
  BLAST_Score score;
  Nlm_FloatHi eValue;
  Nlm_FloatHi eValueThisAlign;
  Nlm_FloatHi Lambda;
  Nlm_FloatHi logK;
  SeqIdPtr subject_id;  /*used to display the sequence in alignment*/
  struct SWResults *next;
} SWResults;

typedef struct proDemographicsItems {
  Uint1 matrixName[MAX_NAME_LENGTH];
  Int4  numSequencesTested;
  Int4  numBetterThanEthresh;
  Int4  numCallsALIGN;
  Int4  queryLength;
  Int4  dbLength;
  Nlm_FloatHi effDbLength;
  Nlm_FloatHi effSearchSpace;
  Nlm_FloatHi X;
  Nlm_FloatHi XinBits;
} proDemographicsItems;

Char * LIBCALL addSuffixToName PROTO((Char *prefix, Char *suffix));


Char LIBCALL getRes PROTO((Char input));

Boolean LIBCALL IMPALAPrintReference PROTO((Boolean html, Int4 line_length, FILE *outfp));

Nlm_FloatHi LIBCALL
impalaKarlinLambdaNR PROTO((BLAST_ScoreFreqPtr sfp, Nlm_FloatHi initialLambda));

void LIBCALL impalaScaling PROTO((posSearchItems *posSearch, compactSearchItems * compactSearch, Nlm_FloatHi scalingFactor));

Boolean LIBCALL  impalaReadCheckpoint PROTO((posSearchItems * posSearch, compactSearchItems * compactSearch, CharPtr fileName, ValNodePtr * error_return,
Nlm_FloatHi scalingFactor));
