
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
/*****************************************************************************
 
File name: pseed3.c
 
Original Author: Zheng Zhang
Maintainer: Alejandro Schaffer
 
Contents: high-level routines for PHI-BLAST and pseed3
 
 
*****************************************************************************/



#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <txalign.h>
#include <simutil.h>
#include <gapxdrop.h>
#include <posit.h>
#include <readdb.h>
#include <ncbithr.h>
#include <seed.h> 

#if defined(OS_UNIX) && !defined(OS_UNIX_SUN) && !defined(OS_UNIX_LINUX)
#include <sys/resource.h>
#include <signal.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif


#define SEED_NTICKS 50



static Int4 get_pat PROTO((FILE *fp, Char *stringForPattern, Char *pname));
static Uint1Ptr reverseSequence PROTO((Uint1Ptr seqFromDb, Int4 lenSeqFromDb));

static void initThreadInfo PROTO((threadInfoItems *threadInfo, BlastSearchBlkPtr search));


/*creates duplicate data structure for each thread of parallel process*/
static seedParallelItems * seedParallelFill(ReadDBFILEPtr rdpt, qseq_ptr query_seq, Int4 lenPatMatch, Boolean is_dna, GapAlignBlkPtr gap_align,
patternSearchItems * patternSearch, seedSearchItems * seedSearch,
threadInfoItems * threadParallelCreate);

static Boolean
reentrant_get_db_chunk PROTO((ReadDBFILEPtr rdpt, Int4Ptr start, Int4Ptr stop, Int4Ptr id_list, Int4Ptr id_list_number, threadInfoItems *threadInfo));


static void do_the_seed_search PROTO((BlastSearchBlkPtr search,
  ReadDBFILEPtr rdfp, Int4 num_seq, 
  qseq_ptr query_seq, Int4 lenPatMatch, Boolean is_dna, 
  GapAlignBlkPtr gap_align, patternSearchItems * patternSearch, 
  seedSearchItems * seedSearch, Int4 * matchIndex, Int4 * totalOccurrences,
  seedResultItems * seedResults));

static void reentrant_tick_proc PROTO((threadInfoItems * threadInfo,
                 Int4 sequence_number));

static Boolean BlastSetLimits PROTO((Int4 cpu_limit, Int2 num_cpu));



/*seedReturn is a list of lists of SeqAligns in which each list is shortened
    to at most number_of_descriptions elements, the rest are deallocated
   and the list of shortened lists is returned*/
ValNodePtr  LIBCALL SeedPruneHitsFromSeedReturn(ValNodePtr seedReturn, Int4 number_of_descriptions)
{
  ValNodePtr returnList; /*list of SeqAlignPtrs to return*/
  ValNodePtr thisValNode; /*scan down input list of ValNodes*/
  SeqAlignPtr thisList; /*one list of seqAlignPtrs*/
  Int4 counter; /*counter for SeqAligns*/
  SeqAlignPtr prevSeqAlign = NULL, currentSeqAlign = NULL , nextSeqAlign = NULL; 
  /*used to walk down list*/
  DenseSegPtr curSegs, testSegs; /*Used to extract ids from curSeqAlign, testSeqAlign*/
  SeqIdPtr curId, testId; /*Ids of target sequences in testSeqAlign*/

  curId = NULL;
  returnList = seedReturn;   
  if (number_of_descriptions > 1) {
    thisValNode = returnList;
    while (NULL != thisValNode) {
      thisList = thisValNode->data.ptrvalue;
      counter = 0;
      while ((thisList != NULL) && (counter < (number_of_descriptions + 1))) {
	if (NULL != currentSeqAlign)
	  prevSeqAlign = currentSeqAlign;
	currentSeqAlign = thisList;
	testSegs = (DenseSegPtr) currentSeqAlign->segs;
	testId = testSegs->ids->next; 
	if ((NULL == curId) || (!(SeqIdMatch(curId, testId)))) {
	  counter++;
	  curSegs = testSegs;
	  curId = testId;
	}
	thisList = thisList->next;
      }
      if (counter == (number_of_descriptions + 1)) {
	prevSeqAlign->next = NULL;
	while(NULL != currentSeqAlign) {
	  nextSeqAlign = currentSeqAlign->next;
	  MemFree(currentSeqAlign);
	  currentSeqAlign = nextSeqAlign;
	}
      }
      thisValNode = thisValNode->next;
    }
  }
  return(returnList);
}

ValNodePtr LIBCALL seedEngineCore(BlastSearchBlkPtr search, BLAST_OptionsBlkPtr options, 
 Uint1Ptr query, Uint1Ptr unfilter_query,
 CharPtr database, CharPtr patfile, Int4 program_flag, 
 FILE * patfp, FILE *outfp, 
 Boolean is_dna, Boolean reverseDb, seedSearchItems *seedSearch, 
 Nlm_FloatHi posEThresh, Nlm_FloatLo searchSpEff,
 posSearchItems *posSearch, SeqLocPtr *seed_seq_loc, Boolean showDiagnostics)
{
        Nlm_FloatHi  dbLength, adjustdbLength;  /*total number of characters in database*/
        qseq_ptr query_seq; /*query sequence in record format*/
        Char  *pattern; /*string description of a pettern*/
        Char *pname; /*name of pattern*/
        Int4 seed; /*position in sequence where pattern match starts*/
        Int4 lenPatMatch; /*number of positions taken by pattern match*/
        Int4  num_seq; /*number of sequences in database*/
	Int4 list[MAX_HIT]; /*list of lengths and start positions where
                              pattern matches sequence*/
        Int4  *occurArray; /*places in query where pattern occurs*/
        Int4  *hitArray; /* beginning and end of pattern in query. */
        Int4  numPatOccur; /*number of pattern occurrences in query string*/
        Int4  effectiveOccurrences; /*number of occurrences not overlapping
                                      in more than half the pattern*/
        Int4  occurIndex;  /*index over pattern ocuurences*/
        Int4  twiceNumMatches; /*stores return value from find_hits*/
        Int4  matchIndex;  /*index for matches to a single sequence*/
        Int4 totalOccurrences = 0; /*total occurrences of pattern in database*/
        Int4 totalBelowEThresh = 0;
        seedResultItems *seedResults = NULL; /*holds list of matching sequences*/
        patternSearchItems *patternSearch = NULL; /*holds parameters
                                                     related to pattern1.c*/
        SeqAlign *thisSeqAlign = NULL; /*return value from output_hits for
					 list of matches to one pattern occurrence*/
        ValNodePtr seqAlignList = NULL;  /*list of SeqAlign lists to pass back*/
        ValNodePtr nextValNodePtr = NULL, lastValNodePtr = NULL;  
            /*pointers into list of SeqAlign lists to pass back*/

	SeqIntPtr seq_int;

	ReadDBFILEPtr rdpt=NULL;  /*holds result of attempt to read database*/
        GapAlignBlkPtr gap_align; /*local holder for gap_align*/

        gap_align = search->gap_align;
        seedResults = (seedResultItems *) ckalloc(sizeof(seedResultItems));
        patternSearch = (patternSearchItems *) ckalloc(sizeof(patternSearchItems));
	rdpt = readdb_new(database, !is_dna);
	if (program_flag == PATTERN_FLAG) {
	    search_pat(rdpt, patfile, is_dna, seedSearch, patternSearch, &(search->error_return), outfp);
	    rdpt = readdb_destruct(rdpt);
	    exit(1);
	} 
	occurArray = (Int4 *) ckalloc(sizeof(Int4)*search->sbp->query_length*2);
	hitArray = (Int4 *) MemNew(sizeof(Int4)*search->sbp->query_length*2);

	dbLength = 0;
	num_seq = readdb_get_num_entries(rdpt);
	dbLength = (Nlm_FloatHi) readdb_get_dblen(rdpt);              
	/*correct the effective size of the database to take out
	  the number of positions at the end of each sequence where
	  pattern match cannot start*/

		 
	while (pattern = get_a_pat(patfp, &pname, occurArray, hitArray, &numPatOccur, 
                                   &effectiveOccurrences,
				   program_flag, unfilter_query, query, search->sbp->query_length, is_dna,
                                   patternSearch, seedSearch, outfp, showDiagnostics, &(search->error_return))) {
          if (patternSearch->patternProbability > PAT_PROB_THRESH &&
	      (patternSearch->patternProbability * dbLength > EXPECT_MATCH_THRESH)) {
             fprintf(outfp, "Pattern %s is too likely to occur in the database to be informative\n",pname);
          }
	  else {
            if (patternSearch->wildcardProduct > WILDCARD_THRESH) {
              fprintf(outfp, "Due to variable wildcards pattern %s is likely to occur too many times in a single sequence\n",pname);
            }
	    else {
	      *seed_seq_loc = NULL;
	      adjustdbLength = dbLength - (num_seq * patternSearch->minPatternMatchLength);
	      if (0.0 < searchSpEff)
		adjustdbLength = searchSpEff;

	      for (occurIndex = 0; occurIndex < numPatOccur; occurIndex++) {
		seed = occurArray[occurIndex];
		if (showDiagnostics)
		  fprintf(outfp, "%s pattern %s at position %d of query sequence\n",
			  pname, pattern, seed);
		if ((twiceNumMatches=find_hits(list, &query[seed-1], search->sbp->query_length-seed+1, FALSE, patternSearch)) < 2 || 
		    list[1] != 0) {
		  fprintf(outfp, "twiceNumMatches=%ld list[1]=%ld\n", (long) twiceNumMatches, (long) list[1]);
		  BlastConstructErrorMessage("seedEngineCore", "pattern does not match the query at the place", 1, &(search->error_return));
		  exit(1);
		}
		seq_int = SeqIntNew();
		seq_int->from = occurArray[occurIndex] - 1;
		seq_int->to = list[2*occurIndex] - list[2*occurIndex+1]
                             + seq_int->from;
		seq_int->id = SeqIdDup(search->query_id);
		ValNodeAddPointer(seed_seq_loc, SEQLOC_INT, seq_int);
		if (program_flag != PAT_MATCH_FLAG) {
		  lenPatMatch = list[0]+1;
		  matchIndex = 0;
		  query_seq = split_target_seq(query, seed, lenPatMatch, search->sbp->query_length);
		  if (showDiagnostics)
		    fprintf(outfp, "effective database length=%.1e\n pattern probability=%.1e\nlengthXprobability=%.1e\n", adjustdbLength, 
			    patternSearch->patternProbability, 
			    patternSearch->patternProbability * adjustdbLength);
		  if (!is_dna) {
		    /*extra caution about what values are tolerated*/
		    seedSearch->cutoffScore = eValueFit(
							MIN(MAX_EVALUE, 10 * options->expect_value), adjustdbLength, 
							seedSearch, effectiveOccurrences, patternSearch->patternProbability);
		  }
		}
		totalOccurrences = 0;

		if (program_flag != PAT_MATCH_FLAG) {

		  do_the_seed_search(search, rdpt, num_seq,
				   query_seq, lenPatMatch, is_dna, 
				   gap_align, patternSearch, seedSearch,
                                   &matchIndex, &totalOccurrences,
                                   seedResults);

		  if (matchIndex > 0) 
		    quicksort_hits(matchIndex, seedResults);
		  if (showDiagnostics)
		    fprintf(outfp,"\nNumber of occurrences of pattern in the database is %d\n", totalOccurrences);
		  search->second_pass_hits += totalOccurrences;
		  search->second_pass_extends += totalOccurrences;
		  thisSeqAlign = output_hits(rdpt, FALSE, query, query_seq, 
					     lenPatMatch, adjustdbLength, gap_align, is_dna, 
					     effectiveOccurrences, seedSearch, seedResults,
					     patternSearch, reverseDb, totalOccurrences, 
					     options->expect_value,  
					     search->query_id, posEThresh, posSearch, 
					     matchIndex, &totalBelowEThresh, showDiagnostics,
					     outfp);

		  nextValNodePtr = (ValNodePtr) MemNew(sizeof(ValNode));
		  nextValNodePtr->data.ptrvalue = thisSeqAlign;
		  if (NULL == seqAlignList)
		    seqAlignList = nextValNodePtr;
		  else
		    lastValNodePtr->next = nextValNodePtr;
		  lastValNodePtr = nextValNodePtr;
		  
		  search->number_of_seqs_better_E += totalBelowEThresh;
		  search->second_pass_good_extends += totalBelowEThresh;
		  seed_free_all(seedResults);
		} /*if (program_flag...) */
	      } /*for (occurIndex ...)*/
	    }  /*else*/
	  } /*else*/
	}
        MemFree(occurArray);
	rdpt = readdb_destruct(rdpt);
        if (seedResults)
	  MemFree(seedResults);
        MemFree(patternSearch);
	return(seqAlignList);
}


/*creates duplicate data structure for each thread of parallel process*/
static seedParallelItems * seedParallelFill(
     ReadDBFILEPtr rdpt, qseq_ptr query_seq,
     Int4 lenPatMatch, Boolean is_dna, GapAlignBlkPtr gap_align,
     patternSearchItems * patternSearch, seedSearchItems * seedSearch,
     threadInfoItems *threadInfo)
{
    seedParallelItems *returnData;

    returnData = (seedParallelItems *) ckalloc(sizeof(seedParallelItems));
    returnData->rdpt = readdb_attach(rdpt);
    if (returnData->rdpt == NULL) {
      MemFree(returnData);
      return NULL;
    }
    returnData->query_seq = query_seq;
    returnData->lenPatMatch = lenPatMatch;
    returnData->gap_align = GapAlignBlkNew(1,1); 
    returnData->gap_align->gap_open = gap_align->gap_open;
    returnData->gap_align->gap_extend = gap_align->gap_extend;
    returnData->gap_align->decline_align = gap_align->decline_align;
    returnData->gap_align->x_parameter = gap_align->x_parameter;
    returnData->gap_align->matrix = gap_align->matrix;

    returnData->is_dna = is_dna;
    returnData->patternSearch = patternSearch;
    returnData->seedSearch = seedSearch;
    returnData->seedResults = 
      (seedResultItems *) ckalloc(sizeof(seedResultItems));
    returnData->totalOccurrences = 0;
    returnData->matchIndex = 0;
    returnData->threadInfo = threadInfo;
    return(returnData);
    /*NEED to free seedResults when concatenating*/
}



/*Handles the time intensive part of the seed search in parallel*/
static VoidPtr
parallel_seed_search(VoidPtr seedParallel)
{
   Int4 seqIndex;  /*index over sequences*/
   Int4 seqIndex1; /*index over array of sequence indices*/
   Int4 lenSeqFromDb;  /* length of seqFromDb */
   Uint1Ptr seqFromDb; /*newly read sequence from database*/
   Int4  newOccurrences = 0; /*number of new occurrences of pattern*/
   hit_ptr hit_list = NULL; /*list of matches to one database sequence*/
   seedParallelItems *localSeedParallel; /*local copy of seedParallel
                                   coerced to be of the desired type*/
   Int4  start=0, stop=0; /*starting and stopping indices for
                            database search*/
   Int4 id_list_length; /*number of ids of strings to test against*/
   Int4Ptr id_list=NULL; /*list of ids of sequences to test against*/

   localSeedParallel = (seedParallelItems * ) seedParallel;
   id_list = NULL;
   if (localSeedParallel->threadInfo->global_gi_being_used)
     id_list = MemNew(localSeedParallel->threadInfo->db_chunk_size*sizeof(Int4));

   while (reentrant_get_db_chunk(localSeedParallel->rdpt, &start, &stop, id_list, &id_list_length,
                   localSeedParallel->threadInfo) 
	  != TRUE) {
     if (id_list) {
       for (seqIndex1 = 0; seqIndex1 < id_list_length; seqIndex1++) {	    	        seqIndex = id_list[seqIndex1];
	 lenSeqFromDb = readdb_get_sequence(localSeedParallel->rdpt, 
					seqIndex, &seqFromDb);
	 hit_list = get_hits(localSeedParallel->query_seq, 
                         localSeedParallel->lenPatMatch, seqFromDb, 
 			 lenSeqFromDb, 
			 localSeedParallel->gap_align, 
			 localSeedParallel->is_dna, 
			 localSeedParallel->patternSearch, 
			 localSeedParallel->seedSearch, &newOccurrences);
	 if (newOccurrences > 0)
	   localSeedParallel->totalOccurrences += newOccurrences;
	 if (hit_list) {
	   storeOneMatch(hit_list, seqIndex, seqFromDb, 
			 localSeedParallel->seedResults); 
	   localSeedParallel->matchIndex++;
	 }
	 /* Emit a tick if needed and we're not MT. */
	 if (localSeedParallel->threadInfo->db_mutex == NULL)
	   reentrant_tick_proc(localSeedParallel->threadInfo,seqIndex);
       }
     }
     else {
       for (seqIndex = start; seqIndex < stop; seqIndex++) {	    	
	 lenSeqFromDb = readdb_get_sequence(localSeedParallel->rdpt, 
					seqIndex, &seqFromDb);
	 hit_list = get_hits(localSeedParallel->query_seq, 
                         localSeedParallel->lenPatMatch, seqFromDb, 
 			 lenSeqFromDb, 
			 localSeedParallel->gap_align, 
			 localSeedParallel->is_dna, 
			 localSeedParallel->patternSearch, 
			 localSeedParallel->seedSearch, &newOccurrences);
	 if (newOccurrences > 0)
	   localSeedParallel->totalOccurrences += newOccurrences;
	 if (hit_list) {
	   storeOneMatch(hit_list, seqIndex, seqFromDb, 
			 localSeedParallel->seedResults); 
	   localSeedParallel->matchIndex++;
	 }
	 /* Emit a tick if needed and we're not MT. */
	 if (localSeedParallel->threadInfo->db_mutex == NULL)
	   reentrant_tick_proc(localSeedParallel->threadInfo,seqIndex);
       }
     }
   }
   return ((VoidPtr) localSeedParallel);
}


/*Converts the string in program to a number; sets is_dna to TRUE,
  if program name indicates DNA sequences will be used*/
Int4 LIBCALL convertProgramToFlag(Char * program, Boolean * is_dna)
{
  Int4 program_flag;  /*flag for program name to return*/

  if (strsame(program, "seed")) { 
    program_flag = SEED_FLAG; 
    *is_dna = TRUE;
  } 
  else 
    if (strsame(program, "pattern")) {
      program_flag = PATTERN_FLAG; 
      *is_dna = TRUE;
    } 
    else 
      if (strsame(program, "patseed")) {
	program_flag = PAT_SEED_FLAG; 
	*is_dna = TRUE;
      } 
      else 
	if (strsame(program, "patmatch")) {
	  program_flag = PAT_MATCH_FLAG; 
	  *is_dna = TRUE;
	} 
	else 
	  if (strsame(program, "seedp")) 
	    program_flag = SEED_FLAG; 
	  else 
	    if (strsame(program, "patternp")) 
	      program_flag = PATTERN_FLAG;
	    else 
	      if (strsame(program, "patseedp")) 
		program_flag = PAT_SEED_FLAG;
	      else 
		if (strsame(program, "patmatchp")) 
		  program_flag = PAT_MATCH_FLAG;
		else {
		  ErrPostEx(SEV_FATAL, 0, 0, "name of program not recognized %s \n", program);
		  return(1);
		}
  return(program_flag);
}


/*count the number of occurrences of pattern in sequences that
  do not overlap by more than half the pattern match length*/
static Int4 countEffectiveOccurrences(Int4 numOccur, Int4 *hitArray, patternSearchItems *patternSearch)
{
    Int4 j; /*loop index*/
    Int4 lastEffectiveOccurrence; /*last nonoverlapping occurrence*/
    Int4 count; /*what to return*/

    count = 0;
    if (numOccur > 0) {
      count = 1;
      lastEffectiveOccurrence = hitArray[0];
      for(j = 1; j < numOccur; j++) {
        if (((hitArray[j] - lastEffectiveOccurrence) * 2) >
	    patternSearch->minPatternMatchLength) {
          lastEffectiveOccurrence = hitArray[j];
          count++;
	}
      }
    }
    return(count);
}

/*Retrieve a pattern and return its string description and pass back the
  hits to query
  fp is the file descriptor to read the pattern from
  name is used to pass back the name of pattern, if any
  hitArray is the start position for matches
  numPatOccur is to pass back the number of pattern occurrences
  program_flag tells in which mode the program is being used
  seq is the query sequence
  unfilter_seq is the unfiltered version of seq, if filtering is on
  len is the length of seq
  is_dna says whether seq is DNA or proteins*/
Char * LIBCALL get_a_pat(FILE *fp, Char **name, Int4Ptr hitArray, Int4Ptr fullHitArray, 
   Int4 * numPatOccur, Int4 *numEffectiveOccurrences, Int4 program_flag, 
   Uint1Ptr unfilter_seq, Uint1Ptr seq, Int4 len, Boolean is_dna,
   patternSearchItems *patternSearch, seedSearchItems * seedSearch, FILE * outfp, Boolean showDiagnostics, ValNodePtr * error_return)
{

    Char line[BUF_SIZE]; /*line of pattern description read in*/
    Char  c; /*single character in the line*/
    Int4 i; /*number of pattern occurrences when occurrences 
             read in from file*/
    Int4  hitIndex; /*index over hits between pattern and seq*/
    Int4 twiceNumHits; /*twice the number of matches*/
    Char tmp[BUF_SIZE]; /*buffer for copying*/
    Int4 *unfilterHitArray = NULL ; /*hitArray for unfiltered hits*/
    Int4 twiceUnfilteredHits; /*twice the number of hits to unfiltered
                                sequence*/
    Int4 linePlace; /*index for characters on a line*/


    if ((program_flag == PAT_SEED_FLAG) || (program_flag == PAT_MATCH_FLAG)) {
      while (get_pat(fp, seedSearch->pat_space, seedSearch->name_space)) {
	if (init_pattern((Uint1 *) seedSearch->pat_space, is_dna, patternSearch, seedSearch, error_return)>=0) {
          if (NULL == unfilter_seq)
	    twiceNumHits = find_hits(hitArray, seq, len, FALSE, patternSearch);
          else {
            unfilterHitArray = (Int4 *) ckalloc(len * sizeof(Int4));
	    twiceNumHits = find_hits(hitArray, seq, len, FALSE, patternSearch);
            twiceUnfilteredHits = find_hits(unfilterHitArray, unfilter_seq, len, FALSE, patternSearch);
            if ((twiceUnfilteredHits > twiceNumHits) && showDiagnostics )
              fprintf(outfp,"\nWARNING: SEG filtering has wiped out %ld occurrence(s) of this pattern\n", (long) ((twiceUnfilteredHits - twiceNumHits)/2));

          }
          if (program_flag == PAT_SEED_FLAG)
            if (showDiagnostics)
	      fprintf(outfp,"\n%ld occurrence(s) of pattern in query\n", (long) twiceNumHits/2);
	  if (twiceNumHits >0) {
	    /* copy start and stop positions. */
	    for (hitIndex=0; hitIndex<twiceNumHits; hitIndex++)
		fullHitArray[hitIndex] = hitArray[hitIndex];
		
	    /*Keep starting positions and add 1 to each starting position*/
	    for (hitIndex = 0; hitIndex < twiceNumHits; hitIndex+=2) {
	      hitArray[hitIndex/2] = hitArray[hitIndex+1]+1;
	    }
            *numEffectiveOccurrences = countEffectiveOccurrences(twiceNumHits/2,hitArray,patternSearch);
	    *numPatOccur = twiceNumHits/2;
	    *name = seedSearch->name_space;
            if (unfilterHitArray)
              MemFree(unfilterHitArray);
	    return seedSearch->pat_space;
	  }
	}
      }
      FileClose(fp);
      if (unfilterHitArray)
        MemFree(unfilterHitArray);
      return NULL;
    }
    i = 0;
    seedSearch->pat_space[0] = '\0';
    while (FileGets(line, BUF_SIZE, fp)) {
      if ((c=line[0]) == '\n') 
	continue;
      if (c == 'I') {
	if (i > 0) {
	  strcpy(tmp, seedSearch->name_space);
	  *name = tmp;
	  strcpy(seedSearch->name_space, &line[4]);
	  *numPatOccur = i;
	  init_pattern((Uint1 *) seedSearch->pat_space, is_dna, patternSearch, seedSearch, error_return);

          *numEffectiveOccurrences = countEffectiveOccurrences(i,hitArray,patternSearch);
	  return seedSearch->pat_space;
	} 
	else {	      
	  strcpy(seedSearch->name_space, &line[4]);
	}
      }
      else 
	if (c != 'H') {
	  if (c == 'P') strcat(seedSearch->pat_space, &line[2]);
	  else 
	    strcat(seedSearch->pat_space, line);
	} 
	else {
          linePlace = 2;
          while ((line[linePlace] < '0') || (line[linePlace] > '9'))
            linePlace++;
   	  hitArray[i++] = atoi(&line[linePlace]);
	}
    }
    if (i > 0) {
      *name = seedSearch->name_space;
      *numPatOccur = i;
      init_pattern((Uint1 *) seedSearch->pat_space, is_dna, patternSearch, seedSearch, error_return);
      *numEffectiveOccurrences = countEffectiveOccurrences(i,hitArray,patternSearch);
      return seedSearch->pat_space;
    }
    return NULL;
}

/*Free all the memory on a list of hits associated with one matching sequence*/
static void free_list(hit_ptr hp)
{
    hit_ptr thisHit, nextHit;  /*placeholders for two items on the list*/
    for (thisHit = hp; thisHit; ) {
	nextHit = thisHit->next;
	MemFree(thisHit);
	thisHit = nextHit;
    }
}

/*Frees all the memory for items on the list of results*/
void LIBCALL seed_free_all(seedResultItems *seedResults)
{
    store_ptr sp, ap; /*this item on list, and next item on lits*/
    for (sp = seedResults->listOfMatchingSequences; sp; sp = ap) {
      free_list(sp->hit_list);
      ap = sp->next;
      MemFree(sp);
    }
    seedResults->listOfMatchingSequences = NULL;
}

/*Set up a structure in which the query sequence is split into three
  pieces
    left piece: from position 0 until just before the seed position
        reversed
    middle piece: from start of seed through end of pattern
    right_piece: from 1 past end of pattern through end of sequence
  put sequences and their lengths into a structure that is returned
  seq is the entire sequence
  seed is the position where the pattern match starts
  len_pat is the number of positions taken by the pattern match
  len_query is the length of seq*/
qseq_ptr LIBCALL split_target_seq(Uint1 *seq, Int4 seed, Int4 len_pat, Int4 len_query)
{
    Int4 i; /*index over seq*/
    Uint1 *str; /*local string used to construct left piece*/
    qseq_ptr qp; /*structure to return*/

    qp = (qseq_ptr) ckalloc(sizeof(query_seq));

    qp->rseq = &seq[seed+len_pat - 1];
    qp->rlen = len_query-seed-len_pat+1;

    qp->sseq = &seq[seed - 1]; 
    qp->slen = len_pat;
    str = (Uint1Ptr) ckalloc(seed);
    /*copy over left piece reversed*/
    for (i = 0; i < seed-1; i++)
      str[i] = seq[seed-2-i];
    str[i] = '\0';
    qp->lseq = str; 
    qp->llen = seed-1;
    return qp;
}

/*insert a hit into hitlist, score, l_score, startPos,
  endPos, bi,bj,ei,ej, mul, describe the hit*/
void insert_hit(Int4 score, Int4 l_score, Int4 startPos, Int4 endPos, 
    Int4 bi, Int4 bj, Int4 ei, Int4 ej, hit_ptr *hit_list, Nlm_FloatHi mul)
{
    hit_ptr newHit, hitListItem; /*new item to initialize and add to list,
                                    and pointer into hit_list */
    
    /*set up all the fields for the new hit*/
    newHit = (hit_ptr) ckalloc(sizeof(hit_node));
    newHit->score = score;
    newHit->mul = mul;
    newHit->l_score = l_score;
    newHit->hit_pos = startPos;
    newHit->hit_end = endPos;
    newHit->bi = bi;
    newHit->bj = bj;
    newHit->ei = ei;
    newHit->ej = ej;
    /*start a new list*/
    if (*hit_list == NULL) {
      *hit_list = newHit; 
      newHit->next = NULL; 
      return;
    }
    /*insert at begining of list*/
    if ((*hit_list)->l_score < l_score) {
      newHit->next = *hit_list; 
      *hit_list = newHit;
    } 
    /*insert in middle or at end of list*/
    else {
      for (hitListItem = *hit_list; hitListItem->next; 
	   hitListItem = hitListItem->next) 
	if (hitListItem->next->l_score <= l_score) 
	  break;
      newHit->next = hitListItem->next;
      hitListItem->next = newHit;
    }
}

/* pos is assumed to be the address of a chracter in the array seq
   if so, copy from pos backwards to the start of seq
   into target.
   return the number of characters copied*/
Int4 reverse_seq(Uint1 *seq, Uint1 *pos, Uint1 *target) 
{
    Uint1 *c; /*index over addresses of characters in seq*/
    Int4 numCopied; /*number of characters copied*/

    for (c = pos, numCopied = 0; c >= seq; c--, numCopied++) 
	*target++ = *c;
    *target = '\0';
    return numCopied;
}

/*seq is a packed DNA sequence with 4 DNA letters packed into one
  item of type Uint1, len is the number of DNA characters, dseq
  is the unpacked sequence*/
void unpack_dna_seq(Uint1Ptr seq, Int4 len, Uint1Ptr dseq)
{
  Int4 remainder; /*number of DNA characters packed into last entry
                    of seq*/
  Int4 i, j; /*indices into seq and dseq respectively*/
  Uint1 ch; /*placeholder for one item in seq*/

  remainder = len % 4;

  for (i = 0, j = 0; i < len/4; i++) {
    dseq[j++] = READDB_UNPACK_BASE_1(seq[i]);
    dseq[j++] = READDB_UNPACK_BASE_2(seq[i]);
    dseq[j++] = READDB_UNPACK_BASE_3(seq[i]);
    dseq[j++] = READDB_UNPACK_BASE_4(seq[i]);
  }
  ch = seq[i];
  while (remainder > 0) { 
    dseq[j++] = READDB_UNPACK_BASE_1(ch);
    ch <<= 2; 
    remainder--;
  } 
}

/* Find matches of input sequence qp against database
   sequence seq_db at those places where qp matches the pattern
   qp is a structured representation of the query sequence
   len-of_pat is the length of the pattern
   seq_db is a sequence from the database
   len_seq_db is the length of seq_db
   gap_align is a structure to keep track of a gapped
   alignment between the query sequence and the database sequence
   is_dna indicates whether the sequences are DNA or protein
   return the list of matches in the form of a hit_ptr list*/
hit_ptr LIBCALL get_hits(qseq_ptr qp, Int4 len_of_pat, 
		 Uint1Ptr seq_db, Int4 len_seq_db, GapAlignBlkPtr gap_align, 
		 Boolean is_dna, patternSearchItems * patternSearch,
                 seedSearchItems * seedSearch, Int4 * newOccurrences) 
{
    Int4 matchIndex; /*Index over matches between seq_db and pattern*/
    Int4  twiceNumMatches;  /*twice the number of matches between the
                              databse sequence and the pattern*/
    Uint1  *matchStart, *matchEnd; /*for matches between database sequence
                                     and pattern*/
    Int4 endPoslseq;  /*ending position in qp->lseq for optimal alignment*/
    Int4 endPoslseq_db;  /*ending position of left part of seq_db reversed in
                 optimal alignment with qp->lseq*/
    Int4  endPosrseq; /*ending position in qp->rseq for optimal alignment*/
    Int4 endPosrseq_db;  /*ending position of right part of seq_db  in
                 optimal alignment with qp->lseq*/
    Nlm_FloatHi mul; /*pattern probability multiplier to get 
				  passed back from align_of_pattern*/
    Int4 patWildcardScore; /*score from wildcards */
    Int4 score; /*score for matching pattern against piec of sequence in qp*/
    Int4  scoreLeft, scoreRight, scoreOutside; /*scores for gapped alignments
                                   of left and right parts*/
    Uint1 *leftPartseq_dbReversed;  /*Holds the left part of seq_db
                                     before the match, in reverse order*/
    Int4 lenLeft;  /*length of leftPartseq_dbReversed*/
    Int4 lenRight; /*length of the part of the database sequence
                     to the right of where the match occurs*/
    Uint1 * buffer = NULL; /*used to unpack the database sequence if
                           is a DNA sequence*/
    Uint1  *seq_db_local; /*local copy of seq_db, unpacked if DNA*/
    Uint1Ptr lseq, rseq, sseq; /*three piece of the original input sequence
                                lseq is reversed, sseq is where the
                                pattern matches*/
    Int4 llen, rlen; /*lengths of lseq, rseq*/
    hit_ptr hit_list = NULL;
    Int4 hitL[MAX_HIT]; /*array of hit pairs between seq_db and pattern
                         two entries are used per hit*/
    
    (*newOccurrences) = 0;
    lseq = qp->lseq; 
    llen = qp->llen;
    rseq = qp->rseq; 
    rlen = qp->rlen;
    sseq = qp->sseq;
    leftPartseq_dbReversed = (Uint1Ptr) ckalloc(len_seq_db+1);
    twiceNumMatches = find_hits(hitL, seq_db, len_seq_db, is_dna, patternSearch);
    if (twiceNumMatches > 0 && is_dna) {
      if (len_seq_db > MAXDNA) {
	ErrPostEx(SEV_FATAL, 0, 0, "MAX DNA Sequence length exceeded %d \n",
MAXDNA);
        exit(1);
      }
      buffer = (Uint1 *) MemNew(MAXDNA * sizeof(Uint1));
      unpack_dna_seq(seq_db, len_seq_db, buffer);
      seq_db_local = buffer;
    } 
    else 
      seq_db_local =  seq_db;
    (*newOccurrences) += (twiceNumMatches/2);
    for (matchIndex = 0; matchIndex < twiceNumMatches; matchIndex+= 2) {
      matchStart = &seq_db_local[hitL[matchIndex+1]];
      matchEnd = &seq_db_local[hitL[matchIndex]];
      score = align_of_pattern(sseq, matchStart, len_of_pat,
			   matchEnd-matchStart+1, NULL, NULL, gap_align, 
			   &patWildcardScore, &mul, patternSearch, seedSearch);
      lenLeft = reverse_seq(seq_db_local, matchStart-1, leftPartseq_dbReversed);
      scoreLeft = SEMI_G_ALIGN((Uchar *) (lseq-1), (Uchar *) (leftPartseq_dbReversed-1), llen, lenLeft, NULL,
                       &endPoslseq, &endPoslseq_db, TRUE, NULL, 	gap_align,0, TRUE);
      lenRight = len_seq_db - hitL[matchIndex] /*strlen(matchEnd)+1*/;
      scoreRight = SEMI_G_ALIGN((Uchar *) (rseq-1), (Uchar *) matchEnd, rlen, lenRight, 
			       NULL, &endPosrseq, &endPosrseq_db, TRUE,NULL, gap_align,0, TRUE);
      scoreOutside = scoreRight + scoreLeft;
      if (scoreOutside > seedSearch->cutoffScore) { 
        /*No longer include patWildCardScore in 2nd argument*/
	insert_hit(scoreLeft+scoreRight+score, scoreLeft+scoreRight,
		   lenLeft+1,  hitL[matchIndex], llen-endPoslseq,lenLeft-endPoslseq_db,
		   llen+len_of_pat+endPosrseq, lenLeft+(matchEnd-matchStart)+1+endPosrseq_db, &hit_list, mul);
      }
    }
    free(leftPartseq_dbReversed);
    if (buffer)
      MemFree(buffer);
    return hit_list;
}

/*Return an e_value for getting raw score rawScore in a database of size 
  databaseSize, paramP is pattern probability*/
static Nlm_FloatHi e_value(Int4 rawScore, Nlm_FloatHi databaseSize, 
    Nlm_FloatHi paramC, Nlm_FloatHi paramLambda, Nlm_FloatHi paramP,
    Int4 effectiveOccurrences )
{
    Nlm_FloatHi E; /*e value to return*/
    Nlm_FloatHi localparamC; /*local version of the "constant" C*/

    if (rawScore < 0) 
      return 0.0;
    if (rawScore == 0)
      localparamC = 1.0;
    else
      localparamC = paramC;
    E = effectiveOccurrences*paramP*databaseSize
         *localparamC*(1.0+paramLambda*rawScore)*exp(-paramLambda*rawScore);
    return E;
}

/*return the integer score that just exceeds the e value threshold eThresh*/
Int4 LIBCALL eValueFit(Nlm_FloatHi eThresh, Nlm_FloatHi dbLength, seedSearchItems *seedSearch, 
          Int4 numOccurrences, Nlm_FloatHi patternProbability)
{
  Int4 targetScore;
  Int4 lowScore, highScore; /*for binary search*/

  lowScore = 0;
  highScore = 20;
  while (e_value(highScore, dbLength, seedSearch->paramC, 
            seedSearch->paramLambda, patternProbability, 
            numOccurrences) > eThresh)
      highScore *= 2;
  targetScore = highScore / 2;
  while (lowScore < targetScore) {
    if (e_value(targetScore, dbLength, seedSearch->paramC, 
		seedSearch->paramLambda, patternProbability, 
		numOccurrences) >  eThresh) {
      lowScore = targetScore;
      targetScore = (lowScore + highScore) / 2;
    }
    else {
      highScore = targetScore;
      targetScore = (lowScore + highScore) / 2;
    }
  }
  return(lowScore);
}

/*Gets the 3 scores of an alignment together into a ScorePtr*/
ScorePtr putScoresInSeqAlign(Int4 rawScore, Nlm_FloatHi eValue, Nlm_FloatHi lambda, Nlm_FloatHi C)
{
  Nlm_FloatHi bitScoreUnrounded; /*conversion for raw score to bit score*/
  ScorePtr returnScore = NULL;

  MakeBlastScore(&returnScore,"score",0.0, rawScore);
  MakeBlastScore(&returnScore,"e_value",eValue,0);
  bitScoreUnrounded = (rawScore*lambda - log(C) -log(lambda*rawScore+1.0))/NCBIMATH_LN2;
  MakeBlastScore(&returnScore,"bit_score",bitScoreUnrounded, 0);
  
  return(returnScore);
}

/*output the matches
    rdpt is a pointer to the sequence database
    score_only determines how much output to give per match
    seq1 is the query sequence
    qp is the split reprsentation of the query sequence
    len is the length of the query sequence
    dbLength is the total length of the database
    gap_align is structure to keep track of a gapped alignment
    is_dna tells whether the sequences are DNA or proteins 
    effectiveOccurrences counts number of not significantly overlapping
    occurrences of pattern in seq1*/
    
SeqAlignPtr LIBCALL output_hits(ReadDBFILEPtr rdpt,
	    Boolean score_only, Uint1 *seq1, qseq_ptr qp, 
	    Int4 len, Nlm_FloatHi dbLength, GapAlignBlkPtr gap_align, Boolean is_dna,
            Int4 effectiveOccurrences,
            seedSearchItems *seedSearch, seedResultItems *seedResults, 
            patternSearchItems * patternSearch, Boolean reverse, 
            Int4 numOccurences, Nlm_FloatHi eThresh,
            SeqIdPtr query_id, Nlm_FloatHi posEthresh, 
            posSearchItems *posSearch, Int4 numMatches,
            Int4 * totalBelowEThresh, Boolean showDiagnostics, FILE * outfp)
{
    store_ptr oneMatch;  /*one matching sequence, used as loop index*/
    hit_ptr oneHit; /*one hit reprsenting one place query matches 
                      database sequence; the list of matching
                      places is in oneMatch->hit_list*/
    Uint1 *dbSeqPrefixReversed; /*prefix of database sequence before the
                                 position where the hit occurs, used
                                 in reverse form for alignment purposes*/
    Int4 lenPrefix; /*length of dbSeqPrefixReversed*/
    Uint1  *buffer; /*buffer to hold an unpacked matching
                            DNA sequence*/
    Uint1  *reverseSeqFromDb;
    Uint1  *lseq, *sseq,  *rseq; /*three pieces of seq1, left of pattern match
                                  at pattern match, right of pattern match*/
    Int4  llen, rlen; /*lengths of lseq, rseq*/
    Char buff[SeqIdBufferSize]; /*buffer to hold descriptor of
                                  database sequence that matches*/
    Nlm_FloatHi mul; /*pattern probability multiplier to get 
				  passed back from align_of_pattern*/
    Int4 patWildcardScore; /*score on wildcards for Altschul stats*/
    Int4  endPosQuery; /*ending position in part of query for optimal alignment*/
    Int4  endPosDbSeq; /*ending position in part of database sequence 
                         for optimal alignment*/
    Int4 *alignScript; /*edit script that describes pairwise alignment*/
    Int4 *reverseAlignScript; /*used to revrse align script for left parts
                              of sequence that are fed in reversed*/
    Int4 *ASptr, *revASptr, temp; /*pointers/indices to the 2 alignment scripts*/
    Int4  lenDbSequence; /*length of sequence from database that matches seq1*/
    SeqIdPtr sip; /*used to extract sequence from database*/
    Nlm_FloatHi eValueForMatch; /*e-value for sequence match*/
    Nlm_FloatHi probability; /*probability of hitting pattern*/
    Boolean  eThreshWarning = TRUE; /*warn about possible missing
                                      sequences*/
    Int4 oldNumber; /*sequence number of previous match*/
    SeqAlignPtr seqAlignList =NULL; /*list of SeqAligns to return*/
    SeqAlignPtr nextSeqAlign; /*new one to add*/
    SeqAlignPtr lastSeqAlign = NULL; /*last one on list*/
    GapXEditBlockPtr nextEditBlock;  /*intermediate structure towards seqlign*/
    Int4 p; /*loop index*/
    Int4 *tempPosThreshSequences, *tempPosResultSequences; /*temporary
                       arrays for copying over old posSearch structure*/
    Int4 totalNumMatches; /*number of <threshold matches from previous calls
                            + numMatches*/
    Int4 posBaseIndex; /*loop index for copying pos results arrays*/

    probability = ((Nlm_FloatHi) numOccurences /  dbLength);
    oldNumber = -1;
    if (score_only) {
      for (oneMatch = seedResults->listOfMatchingSequences; 
	   oneMatch; oneMatch = oneMatch->next) {
        /*retrieve description of sequence*/
	sip = NULL;
	readdb_get_descriptor(rdpt, oneMatch->seqno, &sip, NULL);
	if (sip)
	  SeqIdWrite(sip, buff, PRINTID_FASTA_LONG, sizeof(buff));
	sip = SeqIdSetFree(sip);
        /*print information about the matching sequence*/
	for (oneHit = oneMatch->hit_list; oneHit && eThreshWarning; oneHit = oneHit->next) {
	  eValueForMatch = e_value(oneHit->l_score, dbLength, 
			 seedSearch->paramC, seedSearch->paramLambda,
                         probability, effectiveOccurrences)
	                 *oneHit->mul;

	  if (eValueForMatch <= eThresh) {
	    if (oneMatch->seqno != oldNumber) {
	      fprintf(outfp,"%d\t%s\n", oneMatch->seqno, buff);
	    }
	    fprintf(outfp,"%.3g Total Score %ld Outside Pattern Score %ld Match start in db seq %ld\n       Extent in query seq %ld %ld Extent in db seq %ld %ld\n", eValueForMatch, 
		   oneHit->score, (long) oneHit->l_score, (long) oneHit->hit_pos,
		   (long) oneHit->bi+1, (long) oneHit->ei, (long) oneHit->bj+1, (long) oneHit->ej);
	  }
          else
	    if (oneMatch->seqno != oldNumber)
	      eThreshWarning = FALSE;
	  oldNumber = oneMatch->seqno;
	}
      }
    } 
    else {
      lseq = qp->lseq; 
      rseq = qp->rseq; 
      sseq = qp->sseq;
      llen = qp->llen; 
      rlen = qp->rlen;
      if (posEthresh > 0.0) {
        if (posSearch->threshSequences) {
	  totalNumMatches = posSearch->posNumSequences + numMatches;
          posBaseIndex = posSearch->posNumSequences;
          tempPosThreshSequences = (Int4 *) 
	    MemNew(totalNumMatches * sizeof(Int4));
          for (p = 0; p < posSearch->posNumSequences; p++)
	    tempPosThreshSequences[p] = posSearch->threshSequences[p];
          MemFree(posSearch->threshSequences);
          posSearch->threshSequences = tempPosThreshSequences;
          totalNumMatches = posSearch->posResultsCounter + numMatches;
          tempPosResultSequences = (Int4 *) 
	    MemNew(totalNumMatches * sizeof(Int4));
          for (p = 0; p < posSearch->posResultsCounter; p++)
	    tempPosResultSequences[p] = posSearch->posResultSequences[p];
          MemFree(posSearch->posResultSequences);
          posSearch->posResultSequences = tempPosResultSequences;

        }
        else {
          posSearch->posNumSequences = 0;
          posSearch->threshSequences = (Int4 *) MemNew(numMatches * 
							  sizeof(Int4));
          posSearch->posResultSequences = (Int4 *) MemNew(numMatches * 
							  sizeof(Int4));
          posSearch->posResultsCounter = 0;
          posBaseIndex = 0;
	}
	for (p = posBaseIndex; p < posBaseIndex + numMatches; p++) 
	  posSearch->threshSequences[p] = 0;
      }
      for (oneMatch = seedResults->listOfMatchingSequences; 
	   oneMatch; oneMatch = oneMatch->next) {
        /*retrieve description of sequence*/
	sip = NULL;
	readdb_get_descriptor(rdpt, oneMatch->seqno, &sip, NULL);
	/* if (sip)
	  SeqIdWrite(sip, buff, PRINTID_FASTA_LONG, sizeof(buff)); */
	lenDbSequence = readdb_get_sequence(rdpt, oneMatch->seqno, 
					    (Uint1Ptr *) &(oneMatch->seq));
        if (reverse) {
	  reverseSeqFromDb = reverseSequence(oneMatch->seq,lenDbSequence);
          /*MemFree(oneMatch->seq);*/
          oneMatch->seq = reverseSeqFromDb;
	}
	reverseAlignScript =  (Int4 *) ckalloc(sizeof(Int4)*(lenDbSequence+llen+rlen+3));
	dbSeqPrefixReversed = (Uint1Ptr) ckalloc(lenDbSequence + 1);
        buffer = NULL;
	if (is_dna) {
	  buffer = (Uint1 *) MemNew(MAXDNA * sizeof(Uint1));
	  unpack_dna_seq((Uint1Ptr) oneMatch->seq, lenDbSequence, buffer);
	  oneMatch->seq = buffer;
	}
	for (oneHit = oneMatch->hit_list; oneHit && eThreshWarning; oneHit = oneHit->next) {
          eValueForMatch = e_value(oneHit->l_score, 
           dbLength, seedSearch->paramC, seedSearch->paramLambda, probability, effectiveOccurrences);
          eValueForMatch*=oneHit->mul;
          if (eValueForMatch <= eThresh) {
	    if ((posEthresh > 0.0) && (oneMatch->seqno != oldNumber)) { 
	      if (eValueForMatch < posEthresh) {
                posSearch->threshSequences[posSearch->posNumSequences] = 1;
                posSearch->posResultSequences[posSearch->posResultsCounter] =
                  oneMatch->seqno;
                posSearch->posResultsCounter++;
              }
              posSearch->posNumSequences++;
	    }
	    /* if (oneMatch->seqno != oldNumber)
	       printf("\n%d\t%s\n", oneMatch->seqno, buff);
	       printf("\n E-value = %.3g: \n", eValueForMatch);
	       printf("\n Total Score = %d, Outside Pattern Score = %d: \n", oneHit->score, oneHit->l_score); */
	    lenPrefix = reverse_seq(oneMatch->seq, &oneMatch->seq[oneHit->hit_pos-2], dbSeqPrefixReversed);
	    SEMI_G_ALIGN((Uchar *) (lseq-1), (Uchar *) (dbSeqPrefixReversed-1), llen, lenPrefix,  
			 reverseAlignScript,  &endPosQuery, &endPosDbSeq, FALSE,
			 &alignScript, gap_align,0, TRUE);
	    for (revASptr = reverseAlignScript, ASptr = alignScript-1; 
		 revASptr < ASptr; revASptr++, ASptr--) {
	      temp = *revASptr; 
	      *revASptr = *ASptr; 
	      *ASptr = temp;
	    }
	    align_of_pattern(sseq, &oneMatch->seq[oneHit->hit_pos-1], len,
			     oneHit->hit_end-oneHit->hit_pos+2, alignScript, &alignScript, 
			     gap_align, &patWildcardScore, &mul, patternSearch,
			     seedSearch);
	    SEMI_G_ALIGN((Uchar *) (rseq-1), (Uchar *) &oneMatch->seq[oneHit->hit_end], rlen, 
			 lenDbSequence-oneHit->hit_end-1, 
			 alignScript,  &endPosQuery, &endPosDbSeq, FALSE, &alignScript,gap_align,0, FALSE);
	    /* DISPLAY(seq1-1+oneHit->bi, oneMatch->seq-1+oneHit->bj, 
	       oneHit->ei-oneHit->bi,
	       oneHit->ej-oneHit->bj, reverseAlignScript, oneHit->bi+1, 
	       oneHit->bj+1, is_dna); */
	    
	    nextEditBlock = TracebackToGapXEditBlock(seq1 - 1 + oneHit->bi,
						     oneMatch->seq-1 + oneHit->bj, oneHit->ei - oneHit->bi,
						     oneHit->ej - oneHit->bj, reverseAlignScript, oneHit->bi,
						     oneHit->bj);  
	    nextSeqAlign = GapXEditBlockToSeqAlign(nextEditBlock, sip, 
						   query_id);
	    nextSeqAlign->score = putScoresInSeqAlign(oneHit->l_score,
						      eValueForMatch,seedSearch->paramLambda, seedSearch->paramC);
	    if (NULL == seqAlignList)
	      seqAlignList = nextSeqAlign;
	    else
	      lastSeqAlign->next = nextSeqAlign;
	    lastSeqAlign = nextSeqAlign;
	    if (oneMatch->seqno != oldNumber)
              (*totalBelowEThresh)++;
	  }
	  else {
	    if (oneMatch->seqno != oldNumber)
	      eThreshWarning = FALSE;
	  }
	  oldNumber = oneMatch->seqno;
	}
	sip = SeqIdSetFree(sip);
	MemFree(dbSeqPrefixReversed); 
	MemFree(reverseAlignScript); 
      }
    }
    if (eThreshWarning && showDiagnostics)
      fprintf(outfp, "WARNING: There may be more matching sequences with e-values below the threshold of %lf\n", eThresh);
    if (buffer)
      MemFree(buffer);
    return(seqAlignList);
}

/* Put a new match to seq described by hit_list onto the
   list of matching sequences in seedResults
   seqno is the number of seq in the database indexing*/ 
void LIBCALL storeOneMatch(hit_ptr hit_list, Int4 seqno, Uint1Ptr seq, 
	      seedResultItems *seedResults)
{
    store_ptr sp; /*new node to store this match*/

    sp = (store_ptr) ckalloc(sizeof(store_node));
    sp->seq = seq;
    sp->hit_list = hit_list;
    sp->l_score = hit_list->l_score;
    sp->seqno = seqno;
    sp->next = seedResults->listOfMatchingSequences; 
    seedResults->listOfMatchingSequences = sp;
}



/*Bubble sort the entries in qs from index i through j*/
void bbsort(store_ptr *qs, Int4 i, Int4 j)
{
    Int4 x, y; /*loop bounds for the two ends of the array to be sorted*/
    store_ptr sp; /*temporary pointer for swapping*/
    for (x = j; x > i; x--) {
      for (y = i; y < x; y++) {
	if (qs[y]->l_score > qs[y+1]->l_score) {
	  /*swap pointers for inverted adjacent elements*/
	  sp = qs[y];
	  qs[y] = qs[y+1]; 
	  qs[y+1] = sp;
	}
      }
    }
}

/*quicksort the entries in qs from qs[i] through qs[j] */
void quicksort(store_ptr *qs, Int4 i, Int4 j)
{
    Int4 lf, rt;  /*left and right fingers into the array*/
    Int4 partitionScore; /*score to partition around*/
    store_ptr tp; /*temporary pointer for swapping*/
    if (j-i <= SORT_THRESHOLD) {
      bbsort(qs, i,j);
      return;
    }
    /*implicitly choose qs[i] as the partition element*/
    lf = i+1; 
    rt = j; 
    partitionScore = qs[i]->l_score;
    /*partititon around partitionScore = qs[i]*/
    while (lf <= rt) {
      while (qs[lf]->l_score <= partitionScore) 
	lf++;
      while (qs[rt]->l_score > partitionScore) 
	rt--;
      if (lf < rt) {
	/*swap elements on wrong side of partition*/
	tp = qs[lf];
	qs[lf] = qs[rt];
	qs[rt] = tp;
	rt--;
	lf++;
      } 
      else 
	break;
    }
    /*swap partition element into middle position*/
    tp = qs[i];
    qs[i] = qs[rt]; 
    qs[rt] = tp;
    /*call recursively on two parts*/
    quicksort(qs, i,rt-1); 
    quicksort(qs, rt+1, j);
}

/*Sort the sequences that hit the query by increasing score;
  This routine converts the list starting at 
  seedResults->listOfMatchingSequences to
  an array for sorting and then converts back to a singly-linked list*/
void LIBCALL quicksort_hits(Int4 no_of_seq, seedResultItems *seedResults)
{
    Int4 i;
    store_ptr sp;
    store_node sentinel;
    store_ptr *qs; /*local array for sorting*/

    /*Copy the list starting from seedResults->listOfMatchingSequences 
      into the array qs*/
    qs = (store_ptr *) ckalloc(sizeof(store_ptr)*(no_of_seq+1));
    for (i = 0, sp = seedResults->listOfMatchingSequences; 
	 i < no_of_seq; i++, sp = sp->next) 
      qs[i] = sp;
    /*Put sentinel at the end of the array*/
    qs[i] = &sentinel; 
    sentinel.l_score = SEED_INFINITY;
    quicksort(qs, 0, no_of_seq-1);
    /*Copy back to the list starting at seedResults->listOfMatchingSequences*/
    for (i = no_of_seq-1; i > 0; i--)
      qs[i]->next = qs[i-1];
    qs[0]->next = NULL;
    seedResults->listOfMatchingSequences = qs[no_of_seq-1];
    free(qs);
}

/* ckopen -------------------------------------- open file; check for success */
/* FILE *ckopen(Char *name, Char *mode)
{
  FILE *fp;

  if ((fp = FileOpen(name, mode)) == NULL)
    fatalf("Cannot open %s.", name);
  return fp;
} */

/* ckalloc ----------------------------- allocate space; check for success 
  amount is the number of bytes to allocate*/
/* Char *ckalloc(Int4 amount)
{
  char *p; 

  if ((p = MemNew( (unsigned) amount)) == NULL)
    fatal("Ran out of memory.");
  return p;
} */

/* strsame --------------------------- tell whether two strings are identical */
/* Int4 strsame(char *s, char *t)
{
  return (strcmp(s, t) == 0);
} */

/* strsave ----------------------- save string s somewhere; return address */
/*Char *strsave(Char *s)
{
  Char *copy;
  
  copy = (Char *) ckalloc(strlen(s)+1);	comment: +1 to hold '\0' 
  return strcpy(copy, s);
} */

/*Initialize the order of letters in the alphabet, the score matrix,
and the row sums of the score matrix. matrixToFill is the
score matrix, program_flag says which variant of the program is
used; is_dna tells whether the strings are DNA or protein*/
void LIBCALL init_order(Int4 **matrix, Int4 program_flag, Boolean is_dna, seedSearchItems *seedSearch)
{
    Char i, j; /*loop indices for matrix*/ 
    Int4 *matrixRow; /*row of matrixToFill*/ 
    Nlm_FloatHi rowSum; /*sum of scaled substitution probabilities on matrixRow*/
    
    if (is_dna) {
      seedSearch->order['A'] = 0; 
      seedSearch->order['C'] = 1;
      seedSearch->order['G'] = 2; 
      seedSearch->order['T'] = 3;
    } 
    else {
      for (i = 0; i < ALPHABET_SIZE; i++) 
	seedSearch->order[seedSearch->pchars[i]] = i;
    }
    if (program_flag == SEED_FLAG) {
      for (i = 0; i < ALPHABET_SIZE; i++) 
        seedSearch->charMultiple[i] = 0;
      for (i = 0; i < ALPHABET_SIZE; i++) {
	matrixRow = matrix[i];
        rowSum= 0;
	for (j = 0; j < ALPHABET_SIZE; j++) {
	  rowSum += seedSearch->standardProb[j]*exp(-(seedSearch->paramLambda)*matrixRow[j]);
	}
	seedSearch->charMultiple[i] = rowSum;
      }
    }
}

/*Read in a pattern from a file; fp is the input file descriptor
  stringForPattern stores a string representation of the pattern
  pname is the name of the pattern if there is a name, NULL otherwise
  return value is 1 if pname is not NULL, 0 if pname is NULL*/
static Int4 get_pat(FILE *fp, Char *stringForPattern, Char *pname)
{
    Char line[BUF_SIZE];  /*Line read in from the file*/
    Char *name;  /*name of pattern, if there is a name*/
    Char  *rp;   /*pattern string read in from file*/

    name = NULL; 
    stringForPattern[0] = '\0';
    while (FileGets(line, BUF_SIZE, fp)) {
      /*recognize a pattern name by the string 'ID */
      if (line[0] == 'I' && line[1] == 'D') {
	/*rest of the line is the pattern name*/
	strcpy(pname, &line[3]);
	name = pname;
      } 
      else {
        if (NULL == name) {
          strcpy(pname,"No pattern name given\n");
          name = pname;
        }
	/*recognize pattern content by line starting 'PA'*/
	if (line[0] == 'P' && line[1] == 'A') {
	  /*skip over spaces*/
	  for (rp = &line[2]; *rp == ' '; rp++);
	  /*copy rest of pattern to return parameter stringForPattern*/
	  strcat(stringForPattern, rp);
	  while (FileGets(line, BUF_SIZE, fp)) {
	    if (line[0] == 'P' && line[1] == 'A') {
	      for (rp = &line[2]; *rp == ' '; rp++);
	      strcat(stringForPattern, rp); 
	    } 
	    else {
	      if (line[0] == '/') 
		break;
	      /* PROSITE definition of uniformative*/ 
	      /* if (strncmp(line, "CC   /SKIP-FLAG=TRUE", 20)==0) {
		name = NULL; 
		stringForPattern[0] = '\0'; 
		break;
	      } */
	    }
	  }
	  if (name) {
	    return 1;
	  }
	}
      }
    }
    return 0;
}

/*write out a line with pattern ID pname and pattern content pat*/
static void pat_write_head(char *pat, char *pname, FILE *outfp)
{
    fprintf(outfp,"\nID  %sPA  %s", pname, pat);
}


/*Search for occurrences of all patterns in the file
  patternFileName occurring the database pointed to by rdpt
  is_dna is true if and only if the sequences are DNA sequences*/
void LIBCALL search_pat(ReadDBFILEPtr rdpt, Char *patternFileName, Boolean is_dna, seedSearchItems *seedSearch, patternSearchItems *patternSearch,
ValNodePtr * error_return, FILE *outfp)
{
    FILE *fp;  /*file descriptor for pattern file*/
    Char *pat; /*description of current pattern*/
    Char *pname;  /*name of current paatern*/
    Char buff[SeqIdBufferSize];  /*buffer for description of database sequence*/
    Uint1Ptr seq; /*sequence retrieved from the database*/
    Int4  len; /*length of sequence retrieved from the database*/
    Int4 *hitArray; /*array of hit pairs describing positions of matches
                  between the pattern and a database sequence*/
    Int4  i; /*loop index over matches*/
    Int4  numMatches; /*number of matches to current sequence*/
    Int4  seqno; /*index over all sequences*/
    SeqIdPtr sip; /*Description of the sequence from the database*/

    /*FIXED bug here on amount allocated*/
    hitArray = (Int4 *) ckalloc(sizeof(Int4)*MAX_HIT*2);
    fp = FileOpen(patternFileName, "r");
    pat = ckalloc(PATTERN_BUF_SIZE);
    pname = ckalloc(PATTERN_NAME_SIZE);
    while (get_pat(fp, pat, pname)) {
      if (init_pattern((Uint1 *) pat, is_dna, patternSearch, seedSearch, error_return)>=0) {  
	for (seqno = 0; seqno < rdpt->num_seqs; seqno++) {
          /*if (30 == seqno)
            printf("\n stopping at 30\n"); */
	  len = readdb_get_sequence(rdpt, seqno, &seq);
	  numMatches = find_hits(hitArray, seq, len, is_dna, patternSearch);
	  if (numMatches >0) {
	    sip = NULL;
	    readdb_get_descriptor(rdpt, seqno, &sip, NULL);
	    if (sip)
	      SeqIdWrite(sip, buff, PRINTID_FASTA_LONG, sizeof(buff));
	    sip = SeqIdSetFree(sip);
	    fprintf(outfp,"seqno=%d\t%s\n", seqno, buff);
	    pat_write_head(pat, pname,outfp);
	  }
	  for (i = 0; i < numMatches; i+=2) {
            if (is_dna)
              fprintf(outfp,"HI (%ld %ld)\n", (long) hitArray[i+1], (long) hitArray[i]);
	    else
	      pat_output(seq, hitArray[i+1], hitArray[i], patternSearch, outfp);
	  }
	}
      }
    }
    MemFree(hitArray);
    MemFree(pat); 
    MemFree(pname);
    FileClose(fp);
}

/*reverse the sequence seqFromDb of length lenSeqFromDb*/
static Uint1Ptr reverseSequence(Uint1Ptr seqFromDb, Int4 lenSeqFromDb)
{
  Int4 indexBack, indexFront;
  Uint1Ptr returnSeq;

  returnSeq = (Uint1Ptr) ckalloc(lenSeqFromDb * sizeof(Uint1));
  indexFront = 0;
  indexBack = lenSeqFromDb - 1;
  while (indexFront < lenSeqFromDb) {
    returnSeq[indexFront] = seqFromDb[indexBack];
    indexFront++;
    indexBack--;
  }
  return(returnSeq);
}

static void initThreadInfo(threadInfoItems *threadInfo, BlastSearchBlkPtr search)
{
   threadInfo->global_gi_current = 0;
   threadInfo->global_gi_being_used = FALSE;
   threadInfo->db_mutex = NULL;
   threadInfo->results_mutex = NULL;
   threadInfo->callback_mutex = NULL;
   threadInfo->global_seqid_list = NULL;
   threadInfo->global_seqid_ptr = NULL;
   threadInfo->global_gi_current = 0;
   threadInfo->db_chunk_last = 0;
   threadInfo->db_chunk_size = BLAST_DB_CHUNK_SIZE;  
   threadInfo->final_db_seq = 0;  
   threadInfo->db_incr = 0;
   threadInfo->number_seqs_done = 0;
   threadInfo->blast_gi_list = search->blast_gi_list;
   threadInfo->tick_callback = search->tick_callback;
}


/*
	Function to assign chunks of the database to a thread.  
	The "start" and "stop" points are returned by the arguments.
	Note that this is a half-closed interval (stop is not searched).

	The Int4 "db_chunk_last" (a global variable) keeps track of the last 
	database number assigned and is only changed if the db_mutex has been acquired.

	The Boolean done specifies that the search has already been
	completed.
*/

static Boolean
reentrant_get_db_chunk(ReadDBFILEPtr rdpt, Int4Ptr start, Int4Ptr stop, Int4Ptr id_list, Int4Ptr id_list_number, threadInfoItems *threadInfo)

{
	BlastGiListPtr blast_gi_list;
	Boolean done=FALSE;
        Int4 index, index1, number, gi_list_start, gi_list_end, ordinal_id, last_ordinal_id;
	SeqIdPtr sip;


	if (threadInfo->global_gi_being_used)
	{
		blast_gi_list = threadInfo->blast_gi_list;
		if (threadInfo->global_gi_current < blast_gi_list->total)
		{
			*id_list_number = 0;
			number = 0;
			while (*id_list_number == 0)
			{
				NlmMutexLock(threadInfo->db_mutex);
				number = MIN(threadInfo->db_chunk_size, ((blast_gi_list->total) - threadInfo->global_gi_current));
				if (number <= 0)
				{
					NlmMutexUnlock(threadInfo->db_mutex);
					break;
				}

				gi_list_start = threadInfo->global_gi_current;
				threadInfo->global_gi_current += number;
				gi_list_end = threadInfo->global_gi_current;
				
				NlmMutexUnlock(threadInfo->db_mutex);
				index1 = 0;
				last_ordinal_id = -1;
				for (index=gi_list_start; index<gi_list_end; index++)
				{
					ordinal_id = blast_gi_list->gi_list_pointer[index]->ordinal_id;
					if (ordinal_id > 0 && last_ordinal_id != ordinal_id)
					{
						id_list[index1] = ordinal_id;
						last_ordinal_id = ordinal_id;
						index1++;
					}
				}
				*id_list_number = index1;
			}
		}
		else
		{
			done = TRUE;
		}
	}
	else if (threadInfo->global_seqid_list)
	{    /* global_seqid_list indicates there is a list of selected ID's, 
	     global_seqid_ptr is the position in the list. */
		NlmMutexLock(threadInfo->db_mutex);
		sip = threadInfo->global_seqid_ptr;
		if (sip == NULL)
		{
			NlmMutexUnlock(threadInfo->db_mutex);
			return TRUE;
		}

		ordinal_id = SeqId2OrdinalId(rdpt, sip);
		while (ordinal_id < 0 && sip)
		{
			threadInfo->global_seqid_ptr = threadInfo->global_seqid_ptr->next;
			sip = threadInfo->global_seqid_ptr;
			ordinal_id = SeqId2OrdinalId(rdpt, sip);
		}

		if (threadInfo->global_seqid_ptr)
		{
			*start = ordinal_id;
			*stop = ordinal_id+1;
			done = FALSE;
			threadInfo->global_seqid_ptr = threadInfo->global_seqid_ptr->next;
		}
		else
		{
			done = TRUE;
		}

		NlmMutexUnlock(threadInfo->db_mutex);
	}
	else
	{
		if (threadInfo->db_mutex)
		{
			NlmMutexLock(threadInfo->db_mutex);

			/* Emit a tick if needed. */
			reentrant_tick_proc(threadInfo, threadInfo->db_chunk_last);
			*start = threadInfo->db_chunk_last;
			if (threadInfo->db_chunk_last < threadInfo->final_db_seq)
			{
				*stop = MIN((threadInfo->db_chunk_last+ threadInfo->db_chunk_size), threadInfo->final_db_seq);
			}
			else 
			{/* Already finished. */
				*stop = threadInfo->db_chunk_last;
				done = TRUE;
			}
			threadInfo->db_chunk_last = *stop;
			NlmMutexUnlock(threadInfo->db_mutex);
		}
		else
		{
			if (*stop != threadInfo->final_db_seq)
			{
				done = FALSE;
				*start = 0;
				*stop = threadInfo->final_db_seq;
			}
			else
			{
				done = TRUE;
			}
		}
	}

	return done;
}

static void do_the_seed_search(BlastSearchBlkPtr search,
  ReadDBFILEPtr rdfp, Int4 num_seq, 
  qseq_ptr query_seq, Int4 lenPatMatch, Boolean is_dna, 
  GapAlignBlkPtr gap_align, patternSearchItems * patternSearch, 
  seedSearchItems * seedSearch, Int4 * matchIndex, Int4 * totalOccurrences,
  seedResultItems * seedResults)
{
  Int2 threadIndex, secondIndex; /*indices over threads*/
  Int4 number_of_entries;  /*number of entries in the database*/
  threadInfoItems threadInfo;  /*keeps some information that is manipulated
                                       by the different threads*/

  seedParallelItems **seedParallelArray; /*keeps search information for each thread*/
  TNlmThread PNTR thread_array;
  VoidPtr status=NULL;
  store_ptr endOfList; /*end of list of matching seuqnecs*/
  ReadDBFILEPtr local_rdfp; /*local copy for indexing down a list*/
  if (search == NULL)
    return;

  initThreadInfo(&threadInfo, search);
  threadInfo.number_seqs_done = search->pbp->first_db_seq;  /* The 1st sequence to compare against. */
  threadInfo.db_incr = num_seq / BLAST_NTICKS;

  /*readjustment of final_db_seq needs to be after initThreadInfo*/
  if (search->pbp->final_db_seq > 0)
    threadInfo.final_db_seq = search->pbp->final_db_seq;
  else
    threadInfo.final_db_seq = readdb_get_num_entries_total(search->rdfp);

  if (search->blast_gi_list) {
    threadInfo.global_gi_current = 0;
    threadInfo.global_gi_being_used = TRUE;
  }
  else 
    if (search->blast_seqid_list)	{
      threadInfo.global_seqid_list = search->blast_seqid_list->seqid_list;
      threadInfo.global_seqid_ptr = search->blast_seqid_list->seqid_list;
    }
    else {
      threadInfo.global_seqid_list = NULL;
      threadInfo.global_seqid_ptr = NULL;
    }

  BlastSetLimits(search->pbp->cpu_limit, search->pbp->process_num);	

  if (NlmThreadsAvailable() && search->pbp->process_num > 1) {
    number_of_entries = INT4_MAX;
    /* Look for smallest database. */
    local_rdfp = rdfp;
    while (local_rdfp) {
      number_of_entries = MIN(number_of_entries, readdb_get_num_entries(local_rdfp));
      local_rdfp = local_rdfp->next;
    }
    /* Divide up the chunks differently if db is small. */
    if (threadInfo.db_chunk_size*(search->pbp->process_num) > number_of_entries) {
      /* check that it is at least one. */
      threadInfo.db_chunk_size = MAX(number_of_entries/(search->pbp->process_num), 1);
    }
    NlmMutexInit(&(threadInfo.db_mutex));
    NlmMutexInit(&(threadInfo.results_mutex));

    seedParallelArray = (seedParallelItems **) 
      MemNew((search->pbp->process_num)*sizeof(seedParallelItems *));
    for (threadIndex = 0; threadIndex<search->pbp->process_num; threadIndex++) 
      seedParallelArray[threadIndex] =  seedParallelFill(
	  rdfp, query_seq, lenPatMatch, is_dna, gap_align,
	   patternSearch,  seedSearch, &threadInfo);

    thread_array = (TNlmThread PNTR) MemNew((search->pbp->process_num)*sizeof(TNlmThread));
    for (threadIndex=0; threadIndex<search->pbp->process_num; threadIndex++) {
      thread_array[threadIndex] = NlmThreadCreateEx(parallel_seed_search, (VoidPtr) seedParallelArray[threadIndex], THREAD_RUN|THREAD_BOUND, eTP_Default, NULL, NULL);

      if (NlmThreadCompare(thread_array[threadIndex], NULL_thread)) {
	ErrPostEx(SEV_ERROR, 0, 0, "Unable to open thread.");
      }
    }
    for (threadIndex=0; threadIndex<search->pbp->process_num; threadIndex++) 
      {
	NlmThreadJoin(thread_array[threadIndex], &status);
      }
    
    seedResults->listOfMatchingSequences = NULL;
    for (threadIndex=0; threadIndex< search->pbp->process_num; threadIndex++) {
      if (NULL != seedParallelArray[threadIndex]->seedResults->listOfMatchingSequences) {
        seedResults->listOfMatchingSequences = seedParallelArray[threadIndex]->seedResults->listOfMatchingSequences;
        break;
      }
    }

    endOfList = seedResults->listOfMatchingSequences;
    for (secondIndex= threadIndex; secondIndex< search->pbp->process_num; secondIndex++) {
      while(NULL != endOfList->next)
        endOfList = endOfList->next;
      if (secondIndex < (search->pbp->process_num - 1))
	endOfList->next = seedParallelArray[secondIndex+1]->seedResults->listOfMatchingSequences;
    }
   
    thread_array = MemFree(thread_array);

    NlmMutexDestroy(threadInfo.db_mutex);
    threadInfo.db_mutex = NULL;
    NlmMutexDestroy(threadInfo.results_mutex);
    threadInfo.results_mutex = NULL;
    *matchIndex = 0;
    *totalOccurrences = 0;
    for (threadIndex=0; threadIndex<search->pbp->process_num; threadIndex++) {
      (*matchIndex) += seedParallelArray[threadIndex]->matchIndex;
      (*totalOccurrences) += seedParallelArray[threadIndex]->totalOccurrences;
    }
  }
  else {
    seedParallelArray = (seedParallelItems **) 
      MemNew(1*sizeof(seedParallelItems *));
    seedParallelArray[0] =  seedParallelFill( rdfp, query_seq, 
						lenPatMatch, is_dna, gap_align,
						patternSearch, seedSearch,
						&threadInfo);
    parallel_seed_search((VoidPtr) seedParallelArray[0]);
    (*matchIndex) += seedParallelArray[0]->matchIndex;
    (*totalOccurrences) += seedParallelArray[0]->totalOccurrences;
    seedResults->listOfMatchingSequences = 
      seedParallelArray[0]->seedResults->listOfMatchingSequences;
  }
    for (threadIndex=0; threadIndex < MAX(1,threadIndex<search->pbp->process_num); threadIndex++)
    MemFree(seedParallelArray[threadIndex]);
  
  MemFree(seedParallelArray); 

#ifdef RLIMIT_CPU
  signal(SIGXCPU, SIG_IGN);
#endif
  return;
}


/*
	The function that decides whether or not a tick should be
	emitted.  This is performed through the callback function
	("tick_callback") that is set in "do_the_blast_run".  This
	function is called from "parallel_seed_search" for single processing
	machines and "reentrant_get_db_chunk" 
        for MT machines, after the db_mutex
	has been obtained in "get_db_chunk".
*/

static void reentrant_tick_proc(threadInfoItems * threadInfo,
                 Int4 sequence_number)

{
	if(sequence_number > (threadInfo->number_seqs_done+ threadInfo->db_incr))
	{
		NlmMutexLockEx(&(threadInfo->callback_mutex));
		threadInfo->number_seqs_done += threadInfo->db_incr;
		if (threadInfo->tick_callback)
		{
			threadInfo->tick_callback(sequence_number, 0);
		}
		NlmMutexUnlock(threadInfo->callback_mutex);
	}
	return;
}

Boolean  my_time_out_boolean;

/* Called by UNIX signal for timeouts. */
#ifdef RLIMIT_CPU
static void 
timeout_shutdown(int flag)

{
	my_time_out_boolean = TRUE;
	signal(SIGXCPU, SIG_IGN);
}
#endif

static Boolean BlastSetLimits(Int4 cpu_limit, Int2 num_cpu)

{

#ifdef RLIMIT_CPU
	struct rlimit   rl;

	if (cpu_limit <= 0)
		return TRUE;

	if (getrlimit(RLIMIT_CPU, &rl) != -1 )
	{
		if (rl.rlim_cur == RLIM_INFINITY)
			rl.rlim_cur = 0;
		rl.rlim_cur += cpu_limit/num_cpu;
		setrlimit(RLIMIT_CPU, &rl);
#ifdef SIGXCPU
		sigset(SIGXCPU, timeout_shutdown);
#endif
	}
#endif
	my_time_out_boolean = FALSE;

	return TRUE;

}

#ifdef __cplusplus
}
#endif


