/* $Id: blast.c,v 1.9 2003/10/23 20:15:37 dondosha Exp $
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

File name: blast_driver.c

Author: Ilya Dondoshansky

Contents: Main function for running BLAST

******************************************************************************
 * $Revision: 1.9 $
 * */

static char const rcsid[] = "$Id: blast.c,v 1.9 2003/10/23 20:15:37 dondosha Exp $";

#include <ncbi.h>
#include <sqnutils.h>
#include <readdb.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_message.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_engine.h>
#include "blast_seqalign.h"
#include "blast_format.h"
#include "seqsrc_readdb.h"

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

typedef enum {
   ARG_PROGRAM = 0,
   ARG_DB,
   ARG_QUERY,
   ARG_SUBJECT,
   ARG_STRAND,
   ARG_GENCODE,
   ARG_DBGENCODE,
   ARG_FILTER,
   ARG_LCASE,
   ARG_LOOKUP,
   ARG_MATRIX,
   ARG_MISMATCH,
   ARG_MATCH,
   ARG_WORDSIZE,
   ARG_TEMPL_LEN,
   ARG_TEMPL_TYPE,
   ARG_PHI,
   ARG_THRESHOLD,
   ARG_WINDOW,
   ARG_AG,
   ARG_VARIABLE_WORD,
   ARG_STRIDE,
   ARG_XDROP_UNGAPPED,
   ARG_UNGAPPED,
   ARG_GREEDY,
   ARG_GAPOPEN,
   ARG_GAPEXT,
   ARG_FRAMESHIFT,
   ARG_XDROP,
   ARG_XDROP_FINAL,
   ARG_EVALUE,
   ARG_SEARCHSP,
   ARG_PERC_IDENT,
   ARG_DESCRIPTIONS,
   ARG_ALIGNMENTS,
   ARG_OUT,
   ARG_FORMAT,
   ARG_HTML,
   ARG_ASNOUT
} BlastArguments;

static Args myargs[] = {
   { "Program Name",           /* ARG_PROGRAM */
      NULL, NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
   { "Database name (if not set, second sequence FASTA must be provided)",
     NULL, NULL, NULL, TRUE, 'd', ARG_STRING, 0.0, 0, NULL}, /* ARG_DB */
   { "Query File",               /* ARG_QUERY */
     "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
   { "Subject File (for two sequences comparison", /* ARG_SUBJECT */
     "stdin", NULL, NULL, FALSE, 'j', ARG_FILE_IN, 0.0, 0, NULL},
   { "Query strands to search against database: 0 or 3 is both, 1 is top, "
     "2 is bottom", /* ARG_STRAND */
     "0", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},
   { "Genetic code for translation of the query sequence", /* ARG_GENCODE */
     "0", NULL, NULL, FALSE, 'Q', ARG_INT, 0.0, 0, NULL},
   { "Genetic code for translation of the database", /* ARG_DBGENCODE */
     "0", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},
   { "Filter query sequence (DUST with blastn, SEG with others)", /* ARG_FILTER */
     "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
   { "Mask lower case",
     "F", NULL, NULL, FALSE, 'c', ARG_BOOLEAN, 0.0, 0, NULL}, /* ARG_LCASE */
   { "Use (classical Mega BLAST) lookup table with width 12", 
     "F", NULL, NULL, FALSE, 'L', ARG_BOOLEAN, 0.0, 0, NULL},/* ARG_LOOKUP */
   { "Matrix",                 /* ARG_MATRIX */
     "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
   { "Penalty for a nucleotide mismatch (blastn only)", /* ARG_MISMATCH */
     "-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},
   { "Reward for a nucleotide match (blastn only)", /* ARG_MATCH */
     "1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},
   { "Word size, default if 0 (blastn 11, others 3) ", /* ARG_WORDSIZE */  
     "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL}, 
   { "Length of a discontiguous word template (contiguous word if 0)",
     "0", NULL, NULL, FALSE, 't', ARG_INT, 0.0, 0, NULL}, /* ARG_TEMPL_LEN */
   { "Type of a discontiguous word template (0 - coding, 1 - optimal, "
     "2 - two simultaneous", /* ARG_TEMPL_TYPE */
     "0", NULL, NULL, FALSE, 'T', ARG_INT, 0.0, 0, NULL},
   { "Pattern for PHI BLAST",
     NULL, NULL, NULL, TRUE, 'k', ARG_STRING, 0.0, 0, NULL}, /* ARG_PHI */
   { "Threshold for extending hits, default if zero\n" /* ARG_THRESHOLD */
     "      blastp 11, blastn 0, blastx 12, tblastn 13\n"
     "      tblastx 13, megablast 0",
     "0", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},
   { "Window size (max. allowed distance between a pair of initial hits)", 
     "0", NULL, NULL, FALSE, 'w', ARG_INT, 0.0, 0, NULL}, /* ARG_WINDOW */
   { "Use AG BLAST approach to database scanning", /* ARG_AG */
     "T", NULL, NULL, FALSE, 'A', ARG_BOOLEAN, 0.0, 0, NULL},
   { "Use variable word size approach to database scanning",/* ARG_VARIABLE_WORD */ 
     "F", NULL, NULL, FALSE, 'V', ARG_BOOLEAN, 0.0, 0, NULL},
   { "Database scanning stride", 
     "0", NULL, NULL, FALSE, 's', ARG_INT, 0.0, 0, NULL}, /* ARG_STRIDE */
   { "X dropoff value for ungapped extensions in bits (0 invokes default "
     "behavior)\n      blastn 20, others 7",/*ARG_XDROP_UNGAPPED*/
      "0", NULL, NULL, FALSE, 'y', ARG_INT, 0.0, 0, NULL},
   { "Do only ungapped alignment (always TRUE for tblastx)",/*ARG_UNGAPPED*/
     "F", NULL, NULL, FALSE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
   { "Use greedy algorithm for gapped extensions", /* ARG_GREEDY */
     "F", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
   { "Gap open penalty (default: non-affine if greedy; 5 if dyn. prog.)", 
     "0", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL}, /* ARG_GAPOPEN */
   { "Gap extension penalty (default: non-affine if greedy; 2 otherwise)",
     "0", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL}, /* ARG_GAPEXT */
   { "Frame shift penalty for out-of-frame gapping (blastx, tblastn only)",
     "0", NULL, NULL, FALSE, 'h', ARG_INT, 0.0, 0, NULL}, /* ARG_FRAMESHIFT */
   { "X dropoff value for gapped alignment (in bits) (zero invokes default "
     "behavior)\n      blastn 30, tblastx 0, others 15", /* ARG_XDROP */
     "0", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
   { "X dropoff value for final gapped alignment in bits "
     "(0 invokes default behavior)\n"
     "      blastn 50, tblastx 0, others 25",  /* ARG_XDROP_FINAL */
     "0", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},
   { "Expected value",                         /* ARG_EVALUE */
     "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
   { "Effective length of the search space (use zero for the real size)", 
     "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL}, /* ARG_SEARCHSP */
   {"Identity percentage cut-off",  /* ARG_PERC_IDENT */
    "0", NULL, NULL, FALSE, 'P', ARG_FLOAT, 0.0, 0, NULL},
   { "Number of database sequences to show one-line descriptions for (V)",
     "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL}, /* ARG_DESCRIPTIONS */
   { "Number of database sequence to show alignments for (B)", /* ARG_ALIGNMENTS */
     "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
   { "Final output file name",             /* ARG_OUT */
     "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL}, 
   { "alignment view options:\n0 = pairwise,\n1 = query-anchored showing "
     "identities,\n2 = query-anchored no identities,\n3 = flat "
     "query-anchored, show identities,\n4 = flat query-anchored, no "
     "identities,\n5 = query-anchored no identities and blunt ends,\n6 = "
     "flat query-anchored, no identities and blunt ends,\n7 = XML Blast "
     "output,\n8 = tabular, \n9 tabular with comment lines\n10 ASN, text\n"
     "11 ASN, binary",                         /* ARG_FORMAT */
     "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},
   { "Produce HTML output",                    /* ARG_HTML */
     "F", NULL, NULL, FALSE, 'H', ARG_BOOLEAN, 0.0, 0, NULL},
   { "File name for output in ASN.1 format",   /* ARG_ASNOUT */
     NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL} 
};

/** Fills all the options structures with user defined values. Uses the 
 * myargs global structure obtained from GetArgs.
 * @param lookup_options Lookup table options [in]
 * @param query_setup_options Query options [in]
 * @param word_options Initial word processing options [in]
 * @param ext_options Extension options [in]
 * @param hit_options Hit saving options [out]
 * @param score_options Scoring options [out]
 * @param eff_len_options Effective length options [out]
 * @param psi_options Protein BLAST options [out]
 * @param db_options BLAST database options [out]
 * @param subject_length Length of the subject sequence in 2 sequences 
 *                       case [in]
 */
static Int2 
BLAST_FillOptions(LookupTableOptions* lookup_options,
   QuerySetUpOptions* query_setup_options, 
   BlastInitialWordOptions* word_options,
   BlastExtensionOptions* ext_options,
   BlastHitSavingOptions* hit_options,
   BlastScoringOptions* score_options,
   BlastEffectiveLengthsOptions* eff_len_options,
   PSIBlastOptions* psi_options,
   BlastDatabaseOptions* db_options, 
   BlastSeqSrc* bssp, Int4 subject_length)
{
   char* blast_program;
   Boolean ag_blast = TRUE, variable_wordsize = FALSE, mb_lookup = FALSE;
   Boolean greedy_extension = FALSE;
   Boolean is_gapped = FALSE;
   Int8 totlen = 0;
   Int4 numseqs = 0;
   Uint1 program_number;
   Int2 status;
   Boolean use_pssm = FALSE;

   blast_program = myargs[ARG_PROGRAM].strvalue;
   BlastProgram2Number(blast_program, &program_number);

   /* The following options are for blastn only */
   if (program_number == blast_type_blastn) {
      if (myargs[ARG_TEMPL_LEN].intvalue == 0) {
         ag_blast = (Boolean) myargs[ARG_AG].intvalue;
         mb_lookup = (Boolean) myargs[ARG_LOOKUP].intvalue;
         /* Variable word size can only be used for word sizes divisible 
            by 4 */
         if (myargs[ARG_WORDSIZE].intvalue % COMPRESSION_RATIO == 0)
            variable_wordsize = (Boolean) myargs[ARG_VARIABLE_WORD].intvalue;
      } else {
         /* Discontiguous words */
         ag_blast = FALSE;
         mb_lookup = TRUE;
         variable_wordsize = FALSE;
      }
      greedy_extension = (Boolean) myargs[ARG_GREEDY].intvalue;
   }

   BLAST_FillLookupTableOptions(lookup_options, program_number, mb_lookup,
      myargs[ARG_THRESHOLD].intvalue, myargs[ARG_WORDSIZE].intvalue, 
      ag_blast, variable_wordsize, use_pssm);
   /* Fill the rest of the lookup table options */
   lookup_options->mb_template_length = 
      (Uint1) myargs[ARG_TEMPL_LEN].intvalue;
   lookup_options->mb_template_type = 
      (Uint1) myargs[ARG_TEMPL_TYPE].intvalue;

   if (myargs[ARG_STRIDE].intvalue)
      lookup_options->scan_step = myargs[ARG_STRIDE].intvalue;

   BLAST_FillQuerySetUpOptions(query_setup_options, program_number, 
      myargs[ARG_FILTER].strvalue, myargs[ARG_STRAND].intvalue);

   if (myargs[ARG_GENCODE].intvalue &&
       (program_number == blast_type_blastx || 
        program_number == blast_type_tblastx))
      query_setup_options->genetic_code = myargs[ARG_GENCODE].intvalue;

   BLAST_FillInitialWordOptions(word_options, program_number, 
      greedy_extension, myargs[ARG_WINDOW].intvalue, variable_wordsize, 
      ag_blast, mb_lookup, myargs[ARG_XDROP_UNGAPPED].intvalue);

   BLAST_FillExtensionOptions(ext_options, program_number, greedy_extension, 
      myargs[ARG_XDROP].intvalue, myargs[ARG_XDROP_FINAL].intvalue);

   BLAST_FillScoringOptions(score_options, program_number, greedy_extension, 
      myargs[ARG_MISMATCH].intvalue, myargs[ARG_MATCH].intvalue,
      myargs[ARG_MATRIX].strvalue, myargs[ARG_GAPOPEN].intvalue,
      myargs[ARG_GAPEXT].intvalue);

   if (program_number != blast_type_tblastx)
      is_gapped = !myargs[ARG_UNGAPPED].intvalue;

   score_options->gapped_calculation = is_gapped;
   if (myargs[ARG_FRAMESHIFT].intvalue) {
      score_options->shift_pen = myargs[ARG_FRAMESHIFT].intvalue;
      score_options->is_ooframe = TRUE;
   }

   BLAST_FillHitSavingOptions(hit_options, is_gapped, 
      myargs[ARG_EVALUE].floatvalue, 
      MAX(myargs[ARG_DESCRIPTIONS].intvalue, 
          myargs[ARG_ALIGNMENTS].intvalue));
 
   hit_options->percent_identity = myargs[ARG_PERC_IDENT].floatvalue;
   
   if (bssp) {
      totlen =  BLASTSeqSrcGetTotLen(bssp);
      numseqs = BLASTSeqSrcGetNumSeqs(bssp);
   } else {
      totlen = (Int8) subject_length;
      numseqs = 1;
   }

   BLAST_FillEffectiveLengthsOptions(eff_len_options, 
      numseqs, totlen, (Int8) myargs[ARG_SEARCHSP].floatvalue);

   if (myargs[ARG_DBGENCODE].intvalue && db_options &&
       (program_number == blast_type_tblastn ||
        program_number == blast_type_tblastx)) {
      db_options->genetic_code = myargs[ARG_DBGENCODE].intvalue;
      if ((status = BLAST_GeneticCodeFind(db_options->genetic_code, 
                       &db_options->gen_code_string)))
         return status;
   }

   return 0;
}

Int2 Nlm_Main(void)
{
   BLAST_SequenceBlk* subject = NULL, *query = NULL;
   SeqLoc* subject_slp = NULL; /* SeqLoc for the subject sequence in two
                                    sequences case */
   Boolean query_is_na, db_is_na;
   LookupTableOptions* lookup_options;
   char buf[256] = { '\0' };
   char* blast_program;
   Uint1 program_number;
   BlastInitialWordOptions* word_options;
   BlastScoringOptions* score_options;
   BlastExtensionOptions* ext_options;
   BlastHitSavingOptions* hit_options;
   BlastSeqSrcNewInfo bssn_info;
   ReaddbNewArgs* readdb_args = NULL;
   char* dbname = NULL;
   LookupTableWrap* lookup_wrap;
   Int4 subject_length = 0;
   Int2 status;
   QuerySetUpOptions* query_options=NULL;	
   BlastEffectiveLengthsOptions* eff_len_options=NULL;
   BlastMask* lcase_mask = NULL;
   BlastMask* filter_loc=NULL;	/* All masking locations */
   SeqLoc* query_slp = NULL;
   BlastScoreBlk* sbp = NULL;
   FILE *infp, *outfp;
   BlastQueryInfo* query_info;
   BlastResults* results = NULL;
   Blast_Message* blast_message = NULL;
   SeqAlign* seqalign;
   BlastFormattingOptions* format_options;
   Boolean done;
   BlastReturnStat* return_stats;
   Int4 ctr = 0;
   PSIBlastOptions* psi_options = NULL;
   BlastDatabaseOptions* db_options = NULL;
   ListNode* lookup_segments = NULL;
   Boolean translated_query;
   Int4 num_queries;
   BlastSeqSrc* bssp = NULL;

   if (! GetArgs (buf, NUMARG, myargs))
      return (1);
   
   UseLocalAsnloadDataAndErrMsg ();
   
   if (! SeqEntryLoad())
      return 1;
   
   ErrSetMessageLevel(SEV_WARNING);
   
   if ((outfp = fopen(myargs[ARG_OUT].strvalue, "w")) == NULL) {
      ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open output file %s\n", 
                myargs[ARG_OUT].strvalue);
     return (1);
   }
   
   blast_program = strdup(myargs[ARG_PROGRAM].strvalue);
   BlastProgram2Number(myargs[ARG_PROGRAM].strvalue, &program_number);

   db_is_na = (program_number == blast_type_blastn || 
               program_number == blast_type_tblastn || 
               program_number == blast_type_tblastx);
   query_is_na = (program_number == blast_type_blastn || 
                  program_number == blast_type_blastx || 
                  program_number == blast_type_tblastx);

   BlastGetTypes(myargs[ARG_PROGRAM].strvalue, &query_is_na, &db_is_na);

   BLAST_InitDefaultOptions(program_number, &lookup_options,
      &query_options, &word_options, &ext_options, &hit_options,
      &score_options, &eff_len_options, &psi_options,
      &db_options);

   if (!myargs[ARG_DB].strvalue) {
      FILE *infp2;
      char *subject_file = strdup(myargs[ARG_SUBJECT].strvalue);
      if ((infp2 = fopen(subject_file, "r")) == NULL) {
         ErrPostEx(SEV_FATAL, 1, 0, 
                   "blast: Unable to open second input file %s\n", 
                   subject_file);
         return (1);
      }
      sfree(subject_file);

      BLAST_GetQuerySeqLoc(infp2, db_is_na, 0, 0, 0, NULL, &subject_slp, 
                           0, NULL);
      fclose(infp2);
      
      BLAST_SetUpSubject(program_number, subject_slp, &subject);
      subject_length = subject->length;
      ctr = 1;
   } else {
      readdb_args = (ReaddbNewArgs*) malloc(sizeof(ReaddbNewArgs));
      readdb_args->dbname = dbname = myargs[ARG_DB].strvalue;
      readdb_args->is_protein = !db_is_na;
      
      bssn_info.constructor = &ReaddbSeqSrcNew;
      bssn_info.ctor_argument = (void*) readdb_args;
      bssp = BlastSeqSrcNew(&bssn_info);
   }

   BLAST_FillOptions(lookup_options, query_options, word_options, 
      ext_options, hit_options, score_options, eff_len_options, 
      psi_options, db_options, bssp, subject_length);

   if ((status = BlastFormattingOptionsNew(program_number, 
                    myargs[ARG_OUT].strvalue, 
                    myargs[ARG_DESCRIPTIONS].intvalue, 
                    myargs[ARG_ALIGNMENTS].intvalue, 
                    myargs[ARG_FORMAT].intvalue, &format_options)) != 0)
      return status;
   format_options->html = (Boolean) myargs[ARG_HTML].intvalue;

   BLAST_PrintOutputHeader(format_options, 
      myargs[ARG_GREEDY].intvalue, readdb_args);

   if ((infp = fopen(myargs[ARG_QUERY].strvalue, "r")) == NULL) {
      ErrPostEx(SEV_FATAL, 1, 0, "blast: Unable to open input file %s\n", 
                myargs[ARG_QUERY].strvalue);
      return (1);
   }
   
   return_stats = (BlastReturnStat*) calloc(1, sizeof(BlastReturnStat));

   translated_query = (program_number == blast_type_blastx || 
                       program_number == blast_type_tblastx);

   done = FALSE;
   /* Get the query (queries), loop if necessary. */
   while (!done) {
      if ((Boolean)myargs[ARG_LCASE].intvalue) {
         done = BLAST_GetQuerySeqLoc(infp, query_is_na, 
                   myargs[ARG_STRAND].intvalue, 0, 0,
                   &lcase_mask, &query_slp, ctr, 
                   &num_queries);
      } else {
         done = BLAST_GetQuerySeqLoc(infp, query_is_na,
                   myargs[ARG_STRAND].intvalue, 0, 0, NULL, &query_slp,
                   ctr, &num_queries);
      }

      if (translated_query) {
         BlastMaskDNAToProtein(&lcase_mask, query_slp);
      }

      status = BLAST_SetUpQuery(program_number, query_slp, 
                  query_options, &query_info, &query);

      query->lcase_mask = lcase_mask;

      status = 
         BLAST_MainSetUp(program_number, query_options, score_options, 
            lookup_options, hit_options, query, query_info, &lookup_segments, 
            &filter_loc, &sbp, &blast_message);

      if (translated_query) {
         /* Filter locations were returned in protein coordinates; convert them
            back to nucleotide here */
         BlastMaskProteinToDNA(&filter_loc, query_slp);
      }

      if (status) {
         fprintf(stderr, "BlastSetUp_Standard returned non-zero status: %d\n", 
                 status);
         Blast_MessagePost(blast_message);
         return status;
      }

      BLAST_ResultsInit(query_info->num_queries, &results);
      LookupTableWrapInit(query, lookup_options, 
                          lookup_segments, sbp, &lookup_wrap);
    
      return_stats = (BlastReturnStat*) calloc(1, sizeof(BlastReturnStat));

      if (bssp) {
         BLAST_DatabaseSearchEngine(program_number, query, query_info, 
            bssp, sbp, score_options, lookup_wrap, 
            word_options, ext_options, hit_options, eff_len_options, 
            psi_options, db_options, results, return_stats);
      } else {
         BLAST_TwoSequencesEngine(program_number, query, query_info, 
            subject, sbp, score_options, lookup_wrap, 
            word_options, ext_options, hit_options, eff_len_options, 
            psi_options, db_options, results, return_stats);
      }

      lookup_wrap = LookupTableWrapFree(lookup_wrap);

      /* The following works because the ListNodes' data point to simple
         double-integer structures */
      lookup_segments = ListNodeFreeData(lookup_segments);

      /* Convert results to the SeqAlign form */
      BLAST_ResultsToSeqAlign(program_number, results, query_slp, bssp, 
         subject_slp, score_options, sbp, hit_options->gapped_calculation,
         &seqalign);

      results = BLAST_ResultsFree(results);
      bssp = BlastSeqSrcFree(bssp);
      
      if (myargs[ARG_ASNOUT].strvalue) {
         AsnIoPtr asnout = AsnIoOpen(myargs[ARG_ASNOUT].strvalue, (char*)"w");
         GenericSeqAlignSetAsnWrite(seqalign, asnout);
         asnout = AsnIoClose(asnout);
      }
 
      /* Format the results; note that seqalign and filter locations 
         are freed inside. */
      status = BLAST_FormatResults(seqalign, dbname, 
                  blast_program, query_info->num_queries, query_slp,
                  filter_loc, format_options, score_options->is_ooframe);

      PrintOutputFooter(program_number, format_options, score_options, sbp, 
         lookup_options, word_options, ext_options, hit_options, query_info, 
         readdb_args, return_stats);

      query = BlastSequenceBlkFree(query);
      query_info = BlastQueryInfoFree(query_info);
      BlastScoreBlkFree(sbp);
      query_slp = SeqLocSetFree(query_slp);
   } /* End loop on sets of queries */
   
   subject = BlastSequenceBlkFree(subject);
   sfree(return_stats);
   LookupTableOptionsFree(lookup_options);
   BlastQuerySetUpOptionsFree(query_options);
   BlastExtensionOptionsFree(ext_options);
   BlastHitSavingOptionsFree(hit_options);
   BlastInitialWordOptionsFree(word_options);
   BlastScoringOptionsFree(score_options);
   BlastEffectiveLengthsOptionsFree(eff_len_options);
   PSIBlastOptionsFree(psi_options);
   BlastDatabaseOptionsFree(db_options);
   
   BlastFormattingOptionsFree(format_options);

   if (infp)
      fclose(infp);
   
   sfree(blast_program);

   return status;
}
