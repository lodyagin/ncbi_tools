/* $Id: blast_returns.c,v 1.1 2004/05/14 17:19:03 dondosha Exp $
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

File name: blast_returns.c

Author: Ilya Dondoshansky

Contents: Manipulating data returned from BLAST other than Seq-aligns

Detailed Contents: 

******************************************************************************
 * $Revision: 1.1 $
 * */

static char const rcsid[] = "$Id: blast_returns.c,v 1.1 2004/05/14 17:19:03 dondosha Exp $";

#include <algo/blast/api/blast_returns.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_seqsrc.h>

BLAST_DbSummary* LIBCALL
Blast_DbSummaryFree (BLAST_DbSummary* dbinfo)
{
        BLAST_DbSummary* next;

        if (dbinfo == NULL)
                return NULL;

        while (dbinfo)
        {
                sfree(dbinfo->name);
                sfree(dbinfo->definition);
                sfree(dbinfo->date);
                next = dbinfo->next;
                sfree(dbinfo);
                dbinfo = next;
        }

        return dbinfo;
}

BLAST_DbSummary* Blast_GetDbSummary(const BlastSeqSrc* seq_src)
{
   BLAST_DbSummary* dbinfo = NULL;
   char* chptr = NULL;

   dbinfo = calloc(1, sizeof(BLAST_DbSummary));
   dbinfo->name = BLASTSeqSrcGetName(seq_src);	
      
   if((chptr = BLASTSeqSrcGetDefinition(seq_src)) == NULL)
      chptr = dbinfo->name;

   if (chptr)
      dbinfo->definition = strdup(chptr);	
      
   dbinfo->date = BLASTSeqSrcGetDate(seq_src);

   dbinfo->is_protein = BLASTSeqSrcGetIsProt(seq_src);
     
   if ((dbinfo->total_length = BLASTSeqSrcGetTotLen(seq_src)) == 0)
      dbinfo->total_length = BLASTSeqSrcGetMaxSeqLen(seq_src); 
   dbinfo->number_seqs = BLASTSeqSrcGetNumSeqs(seq_src);

   return dbinfo;
}

/*
        adds the new string to the buffer, separating by a tilde.
        Checks the size of the buffer for Blast_GetParametersBuffer and
        allocates longer replacement if needed.
*/

static Boolean
add_string_to_buffer(char* buffer, char* *old, Int2* old_length)

{
        char* new,* ptr;
        Int2 length = 0, new_length;

        if (!buffer) 
           return FALSE;

        if (*old)
           length = (strlen(*old));

        if((Int2)(strlen(buffer)+length+3) > *old_length)
        {
                new_length = *old_length + 255;
                new = calloc(new_length, sizeof(char));
                if (*old_length > 0 && *old != NULL)
                {
                        memcpy(new, *old, *old_length);
                        sfree(*old);
                }
                *old = new;
                *old_length = new_length;
        }

        ptr = *old;
        ptr += length;

        /* Add a tilde */
        *ptr = '~';
        ptr++;

        while (*buffer != NULLB)
        {
                *ptr = *buffer;
                buffer++; ptr++;
        }

        return TRUE;
}

char*
Blast_GetParametersBuffer(Uint1 program_number, 
   const BlastScoringOptions* score_options,
   const BlastScoreBlk* sbp, const LookupTableOptions* lookup_options,
   const BlastInitialWordOptions* word_options,
   const BlastExtensionOptions* ext_options,
   const BlastHitSavingOptions* hit_options,
   const BlastEffectiveLengthsOptions* eff_len_options,
   const BlastQueryInfo* query_info, const BlastSeqSrc* seq_src, 
   const BlastDiagnostics* diagnostics)
{
   Int4 cutoff = 0;
   char buffer[128];
   char* ret_buffer;
   Int2 ret_buffer_length;
   Int4 num_entries = 0;
   Int8 total_length = 0;
   Int4 qlen = 0;
   double evalue;
   Int2 num_frames;
   Boolean single_query;
   Blast_KarlinBlk* kbp;
   BlastUngappedStats* ungapped_stats = NULL;
   BlastGappedStats* gapped_stats = NULL;
   BlastRawCutoffs* raw_cutoffs = NULL;
   
   ret_buffer = NULL;
   ret_buffer_length = 0;
   
   if (program_number == blast_type_blastx ||
       program_number == blast_type_tblastx)
      num_frames = NUM_FRAMES;
   else if (program_number == blast_type_blastn)
      num_frames = 2;
   else
      num_frames = 1;

   single_query = (query_info->last_context < num_frames);

   if (diagnostics) {
      ungapped_stats = diagnostics->ungapped_stat;
      gapped_stats = diagnostics->gapped_stat;
      raw_cutoffs = diagnostics->cutoffs;
   }

   if (score_options->matrix) {
      sprintf(buffer, "Matrix: %s", score_options->matrix);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }

   if (score_options->gapped_calculation) {
      sprintf(buffer, "Gap Penalties: Existence: %ld, Extension: %ld",
              (long) score_options->gap_open, 
              (long) score_options->gap_extend);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }
   
   if (eff_len_options->db_length)
      total_length = eff_len_options->db_length;
   else if (seq_src) {
      if ((total_length = BLASTSeqSrcGetTotLen(seq_src)) == 0)
         total_length = BLASTSeqSrcGetMaxSeqLen(seq_src);
   }
      
   if (program_number == blast_type_tblastn || 
       program_number == blast_type_rpstblastn ||
       program_number == blast_type_tblastx)
      total_length /= 3;

   if (eff_len_options->dbseq_num)
      num_entries = eff_len_options->dbseq_num;
   else if (seq_src)
      num_entries = BLASTSeqSrcGetNumSeqs(seq_src);

   sprintf(buffer, "Number of Sequences: %ld", (long) num_entries);
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   if (ungapped_stats) {
      sprintf(buffer, "Number of Hits to DB: %s", 
              Nlm_Int8tostr(ungapped_stats->lookup_hits, 1));
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   
      sprintf(buffer, "Number of extensions: %ld", 
              (long) ungapped_stats->init_extends);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      sprintf(buffer, "Number of successful extensions: %ld", 
              (long) ungapped_stats->good_init_extends);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }

   if (gapped_stats) {
      if (hit_options->expect_value > 0.1) {
         sprintf(buffer, "Number of sequences better than %4.1f: %ld", 
                 hit_options->expect_value, 
                 (long) gapped_stats->num_seqs_passed);
      } else {
         sprintf(buffer, "Number of sequences better than %3.1e: %ld", 
                 hit_options->expect_value, 
                 (long) gapped_stats->num_seqs_passed);
      }
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   
      if (score_options->gapped_calculation) {
         sprintf(buffer, 
                 "Number of HSP's better than %4.1f without gapping: %ld", 
                 hit_options->expect_value, 
                 (long) gapped_stats->seqs_ungapped_passed);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
         sprintf(buffer, 
                 "Number of HSP's gapped: %ld", 
                 (long) gapped_stats->extensions);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
         sprintf(buffer, 
                 "Number of HSP's successfully gapped: %ld", 
                 (long) gapped_stats->good_extensions);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
         sprintf(buffer, 
               "Number of extra gapped extensions for HSPs above %4.1f: %ld", 
                 hit_options->expect_value,
                 (long) gapped_stats->extra_extensions);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      }
   }

   /* Query length makes sense only for single query sequence. */
   if (single_query) {
      qlen = BLAST_GetQueryLength(query_info, query_info->first_context);
      sprintf(buffer, "Length of query: %ld", (long)qlen);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }

   sprintf(buffer, "Length of database: %s", Nlm_Int8tostr (total_length, 1));
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   
   if (single_query) {
      Int4 length_adjustment = 
         query_info->length_adjustments[query_info->first_context];
      Int4 eff_qlen;
      Int8 eff_dblen;
      sprintf(buffer, "Length adjustment: %ld", (long) length_adjustment);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);

      /** FIXME: Should this be different for RPS BLAST? */
      eff_qlen = qlen - length_adjustment;
      sprintf(buffer, "Effective length of query: %ld", (long)eff_qlen);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      
      eff_dblen = total_length - num_entries*length_adjustment;
      sprintf(buffer, "Effective length of database: %s", 
              Nlm_Int8tostr (eff_dblen , 1));
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      sprintf(buffer, "Effective search space: %8.0f", 
              ((double) eff_dblen)*((double) eff_qlen));
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      sprintf(buffer, "Effective search space used: %8.0f",
         (double)query_info->eff_searchsp_array[query_info->first_context]);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }
   sprintf(buffer, "Neighboring words threshold: %ld", 
           (long) lookup_options->threshold);
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   sprintf(buffer, "Window for multiple hits: %ld", 
           (long) word_options->window_size);
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   
   if (raw_cutoffs) {
      kbp = sbp->kbp[query_info->first_context];
      sprintf(buffer, "X1: %ld (%4.1f bits)", 
              (long)raw_cutoffs->x_drop_ungapped, 
              raw_cutoffs->x_drop_ungapped*kbp->Lambda/NCBIMATH_LN2);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      if (score_options->gapped_calculation) {
         sprintf(buffer, "X2: %ld (%4.1f bits)", 
                 (long)raw_cutoffs->x_drop_gap, ext_options->gap_x_dropoff);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
         sprintf(buffer, "X3: %ld (%4.1f bits)", 
                 (long)raw_cutoffs->x_drop_gap_final,
                 ext_options->gap_x_dropoff_final);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      
         sprintf(buffer, "S1: %ld (%4.1f bits)", 
                 (long)raw_cutoffs->gap_trigger, ext_options->gap_trigger);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      }
   }
   cutoff = 0;
   if (single_query) {
      Int4 context = query_info->first_context;
      double searchsp = (double) query_info->eff_searchsp_array[context];

      /* For translated RPS blast the search space must be scaled down */
      if (program_number == blast_type_rpstblastn)
         searchsp = searchsp / NUM_FRAMES;

      evalue = hit_options->expect_value;
      if (!score_options->gapped_calculation)
         kbp = sbp->kbp[query_info->first_context];
      else
         kbp = sbp->kbp_gap[query_info->first_context];
   
      BLAST_Cutoffs(&cutoff, &evalue, kbp, searchsp, FALSE, 0);
      sprintf(buffer, "S2: %ld (%4.1f bits)", (long) cutoff, 
              (((cutoff)*(kbp->Lambda))-(kbp->logK))/NCBIMATH_LN2);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }
   return ret_buffer;
}

/** Save the Karlin-Altschul parameters calculated in the BLAST search.
 * @param sbp Internal scoring block structure [in]
 * @param context Index into the array of structures containing
 *                Karlin-Altschul parameters [in]
 * @param sum_returns Returns summary structure [out]
*/
static void 
Blast_SummaryFillKAParameters(const BlastScoreBlk* sbp, Int4 context, 
                              BLAST_SummaryReturn* sum_returns)
{
   Blast_KarlinBlk* kbp;

   if (!sbp)
      return;

   if (sbp->kbp) {
      kbp = sbp->kbp[context];
      sum_returns->ka_params = 
         (BLAST_KAParameters*) malloc(sizeof(BLAST_KAParameters));
      sum_returns->ka_params->Lambda = kbp->Lambda;
      sum_returns->ka_params->K = kbp->K;
      sum_returns->ka_params->H = kbp->H;
   }

   if (sbp->kbp_gap) {
      kbp = sbp->kbp_gap[context];
      sum_returns->ka_params_gap = 
         (BLAST_KAParameters*) malloc(sizeof(BLAST_KAParameters));
      sum_returns->ka_params_gap->Lambda = kbp->Lambda;
      sum_returns->ka_params_gap->K = kbp->K;
      sum_returns->ka_params_gap->H = kbp->H;
   }
}

void Blast_SummaryReturnFill(Uint1 program_number, 
        const BlastScoringOptions* score_options, 
        const BlastScoreBlk* sbp,
        const LookupTableOptions* lookup_options,
        const BlastInitialWordOptions* word_options,
        const BlastExtensionOptions* ext_options,
        const BlastHitSavingOptions* hit_options,
        const BlastEffectiveLengthsOptions* eff_len_options,
        const BlastQueryInfo* query_info,
        const BlastSeqSrc* seq_src,
        const BlastDiagnostics* diagnostics, 
        BLAST_SummaryReturn** sum_returns_out)
{
   BLAST_SummaryReturn* sum_returns = 
      (BLAST_SummaryReturn*) calloc(1, sizeof(BLAST_SummaryReturn));
   Blast_SummaryFillKAParameters(sbp, query_info->first_context, sum_returns);
   sum_returns->params_buffer = 
      Blast_GetParametersBuffer(program_number, score_options, sbp, 
                                lookup_options, word_options, ext_options,
                                hit_options, eff_len_options, query_info, 
                                seq_src, diagnostics);
   *sum_returns_out = sum_returns;
}
