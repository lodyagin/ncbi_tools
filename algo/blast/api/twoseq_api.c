/* $Id: twoseq_api.c,v 1.13 2004/06/08 17:47:24 dondosha Exp $
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

  File name: twoseq_api.c

  Author: Jason Papadopoulos

  Contents: Functions for C toolkit applications to compare two sequences
                using the rewritten blast engine 

***************************************************************************/

#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_message.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_engine.h>
#include <algo/blast/core/mb_lookup.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/hspstream_collector.h>
#include <algo/blast/api/multiseq_src.h>
#include <algo/blast/api/blast_seqalign.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/api/twoseq_api.h>
#include <algo/blast/api/blast_returns.h>
/* For the AdjustOffSetsInSeqalign function */
#include <sequtil.h>

Int2 BLAST_SummaryOptionsInit(BLAST_SummaryOptions **options)
{
    BLAST_SummaryOptions *new_options = (BLAST_SummaryOptions *)malloc(
                                                 sizeof(BLAST_SummaryOptions));
    if (new_options == NULL) {
        *options = NULL;
        return -1;
    }
    
    new_options->hint = eSensitive;
    new_options->program = eChoose;
    new_options->strand = Seq_strand_both;
    new_options->cutoff_evalue = 10.0;
    new_options->matrix = NULL;
    new_options->filter_string = NULL;
    new_options->word_size = 0;
    new_options->gapped_calculation = TRUE;
    new_options->nucleotide_match = 1;
    new_options->nucleotide_mismatch = -3;
    new_options->word_threshold = 0;
    new_options->init_seed_method = eDefaultSeedType;

    *options = new_options;
    return 0;
} 

BLAST_SummaryOptions*
BLAST_SummaryOptionsFree(BLAST_SummaryOptions *options)
{
    free(options);
    return NULL;
}

static Int2 
BLAST_FillOptions(Uint1 program_number,
    const BLAST_SummaryOptions* basic_options,
    LookupTableOptions* lookup_options,
    QuerySetUpOptions* query_setup_options, 
    BlastInitialWordOptions* word_options,
    BlastExtensionOptions* ext_options,
    BlastHitSavingOptions* hit_options,
    BlastScoringOptions* score_options,
    BlastEffectiveLengthsOptions* eff_len_options,
    PSIBlastOptions* psi_options,
    BlastDatabaseOptions* db_options, 
    Int4 query_length,
    Int4 subject_length, RPSInfo *rps_info)
{
    Boolean do_megablast = FALSE;
    Boolean do_ag_blast = FALSE;
    Boolean do_discontig = FALSE;
    Int4 greedy_align = 0;
    Int4 word_size = basic_options->word_size;
    char *matrix;

    if (program_number == blast_type_blastn) {
        if (basic_options->strand != Seq_strand_plus &&
            basic_options->strand != Seq_strand_minus &&
            basic_options->strand != Seq_strand_both) {
            return -2;
        }

        /* If one of the sequences is large enough,
           set up a megablast search */

        if (basic_options->hint == eFast ||
            query_length > MEGABLAST_CUTOFF || 
            subject_length > MEGABLAST_CUTOFF) {
            do_megablast = TRUE;
            if (basic_options->gapped_calculation)
                greedy_align = 1;       /* one-pass, no ungapped */
        }

        /* For a megablast search or a blastn search with
           a larger-than-default word size, turn on striding */

        if (word_size > 13 || do_megablast)
            do_ag_blast = TRUE;

        /* If megablast was turned on but the input indicates
           a sensitive search is desired, switch to discontiguous
           megablast (hardwired to the 11-of-21 optimal template) */

        if (do_megablast && basic_options->hint == eSensitive) {
           if (word_size == 0)
              word_size = 11;
            do_discontig = TRUE;
            do_ag_blast = FALSE;
        }
    }
    

    BLAST_FillLookupTableOptions(lookup_options, 
                                 program_number, 
                                 do_megablast,
                                 basic_options->word_threshold,
                                 word_size,
                                 do_ag_blast,
                                 0,              /* no variable wordsize */ 
                                 FALSE);         /* no PSSM */ 
 
    /* Fill the rest of the lookup table options */
    if (do_discontig) {
        lookup_options->mb_template_length = 21; 
        lookup_options->mb_template_type = MB_WORD_OPTIMAL;
    }
    else {
        lookup_options->mb_template_length = 0; 
        lookup_options->mb_template_type = 0;
    }
    
    BLAST_FillQuerySetUpOptions(query_setup_options, 
                                program_number, 
                                basic_options->filter_string,
                                basic_options->strand);
 
    BLAST_FillInitialWordOptions(word_options, 
                                 program_number, 
                                 (Boolean)greedy_align, 
                                 0,      /* default for ungapped extensions */
                                 0,              /* no variable wordsize */ 
                                 do_ag_blast,
                                 do_megablast,
                                 0);     /* default ungapped X dropoff */
 
    /* If we need to enforce a single-hit method, reset window size to 0. 
       To enforce two-hit method, set window size to a default non-zero 
       value */
    if (basic_options->init_seed_method == eOneHit)
       word_options->window_size = 0;
    else if (basic_options->init_seed_method == eTwoHits)
       word_options->window_size = BLAST_WINDOW_SIZE_PROT;

    BLAST_FillExtensionOptions(ext_options, 
                               program_number, 
                               greedy_align, 
                               basic_options->gap_x_dropoff,
                               0);       /* default final X dropoff */
 
    if (basic_options->matrix == NULL)
        matrix = "BLOSUM62";
    else
        matrix = basic_options->matrix;

    BLAST_FillScoringOptions(score_options, 
                             program_number, 
                             (Boolean)greedy_align, 
                             basic_options->nucleotide_mismatch,
                             basic_options->nucleotide_match,
                             matrix,
                             basic_options->gap_open,
                             basic_options->gap_extend);
 
    score_options->gapped_calculation = basic_options->gapped_calculation;
 
    BLAST_FillHitSavingOptions(hit_options, 
                               basic_options->cutoff_evalue,
                               0);    /* default number of alignments saved */
  
    hit_options->percent_identity = 0;   /* no percent identity cutoff */
  
    eff_len_options->db_length = basic_options->db_length;

    return 0;
}

Int2 
BLAST_TwoSequencesSearch(BLAST_SummaryOptions *basic_options,
                         BioseqPtr bsp1, BioseqPtr bsp2, 
                         SeqAlign **seqalign_out)
{
    Uint1 program_type = blast_type_undefined;
    SeqLocPtr query_slp = NULL;      /* sequence variables */
    SeqLocPtr subject_slp = NULL;
    Boolean seq1_is_aa, seq2_is_aa;
    Int2 status = 0;

    /* sanity checks */

    *seqalign_out = NULL;
    if (bsp1 == NULL || bsp2 == NULL)
        return 0;

    seq1_is_aa = ISA_aa(bsp1->mol);
    seq2_is_aa = ISA_aa(bsp2->mol);

    /* Find program type consistent with the sequences. */
    if (!seq1_is_aa && !seq2_is_aa) {
       if (basic_options->program == eTblastx)
          program_type = blast_type_tblastx;
       else
          program_type = blast_type_blastn;
    } else if (seq1_is_aa && seq2_is_aa) {
       program_type = blast_type_blastp;
    } else if (!seq1_is_aa && seq2_is_aa) {
       program_type = blast_type_blastx;
    } else if (seq1_is_aa && !seq2_is_aa) {
       program_type = blast_type_tblastn;
    }

    /* Check if program type in options is consistent with the one determined
       from sequences. */
    if (basic_options->program == eChoose)
       basic_options->program = program_type;
    else if (basic_options->program != program_type)
       return -1;

    /* Convert the bioseqs into seqlocs. */

    ValNodeAddPointer(&query_slp, SEQLOC_WHOLE,
                      SeqIdDup(SeqIdFindBest(bsp1->id, SEQID_GI)));
    if (!query_slp)
       return -1;
    ValNodeAddPointer(&subject_slp, SEQLOC_WHOLE,
                      SeqIdDup(SeqIdFindBest(bsp2->id, SEQID_GI)));
    if (!subject_slp)
       return -1;

    status = BLAST_TwoSeqLocSets(basic_options, query_slp, subject_slp, 
                                 seqalign_out, NULL, NULL);
    SeqLocFree(query_slp);
    SeqLocFree(subject_slp);

    return status;
}

static Int4 SeqLocListLen(SeqLoc* seqloc)
{
   Int4 length = 0;

   for ( ; seqloc; seqloc = seqloc->next)
      length += SeqLocLen(seqloc);

   return length;
}

/** Compares one list of SeqLoc's against another list of SeqLoc's using the
 * BLAST algorithm.
 */
Int2 
BLAST_TwoSeqLocSets(const BLAST_SummaryOptions *basic_options,
                    SeqLoc* query_seqloc, SeqLoc* subject_seqloc, 
                    SeqAlign **seqalign_out,
                    SeqLoc** filter_out,
                    BLAST_SummaryReturn* *extra_returns)
{
    Uint1 program_type;
    BlastSeqSrc *seq_src = NULL;
    BLAST_SequenceBlk *query = NULL;
    BlastQueryInfo* query_info = NULL;

    ListNode* lookup_segments = NULL;      /* query filtering structures */
    BlastMaskLoc* filter_loc = NULL;

    LookupTableOptions* lookup_options = NULL;  /* new engine options */
    BlastInitialWordOptions* word_options = NULL;
    BlastScoringOptions* score_options = NULL;
    BlastExtensionOptions* ext_options = NULL;
    BlastHitSavingOptions* hit_options = NULL;
    LookupTableWrap* lookup_wrap = NULL;
    QuerySetUpOptions* query_options = NULL;	
    BlastEffectiveLengthsOptions* eff_len_options = NULL;
    BlastScoreBlk* sbp = NULL;
    BlastHSPResults* results = NULL;
    BlastDiagnostics* diagnostics = NULL;
    PSIBlastOptions* psi_options = NULL;
    BlastDatabaseOptions* db_options = NULL;
    RPSInfo *rps_info = NULL;
    Int2 status = 0;
    BlastHSPStream* hsp_stream;

    switch(basic_options->program) {
    case eBlastn:
       program_type = blast_type_blastn;
       break;
    case eBlastp:
       program_type = blast_type_blastp;
       break;
    case eBlastx:
       program_type = blast_type_blastx;
       break;
    case eTblastn:
       program_type = blast_type_tblastn;
       break;
    case eTblastx:
       program_type = blast_type_tblastx;
       break;
    default:
       return -1;
    }

    /* fill in the engine-specific options */

    status = BLAST_InitDefaultOptions(program_type, &lookup_options,
                &query_options, &word_options, &ext_options, 
                &hit_options, &score_options, &eff_len_options, 
                NULL, &db_options);
    if (status != 0)
        goto bail_out;

    if (program_type == blast_type_tblastn || 
        program_type == blast_type_tblastx) {
       if ((status = BLAST_GeneticCodeFind(db_options->genetic_code,
                                           &db_options->gen_code_string)))
          return status;
    }

    seq_src = MultiSeqSrcInit(subject_seqloc, program_type);
    if (seq_src == NULL)
        goto bail_out;

    status = BLAST_FillOptions(program_type, basic_options, 
                 lookup_options, query_options, word_options, 
                 ext_options, hit_options, score_options, 
                 eff_len_options, psi_options, db_options, 
                 SeqLocListLen(query_seqloc), 
                 BLASTSeqSrcGetMaxSeqLen(seq_src), 
                 rps_info);
    if (status != 0)
        goto bail_out;


    /* convert the seqlocs into SequenceBlks, and fill in query_info */

    status = BLAST_SetUpQuery(program_type, query_seqloc, query_options,
                     &query_info, &query);
    if (status != 0)
        goto bail_out;

    /* perform final setup */

    status = BLAST_ValidateOptions(program_type, ext_options,
                score_options, lookup_options, hit_options, NULL);
    if (status != 0)
       goto bail_out;

    status = BLAST_MainSetUp(program_type, query_options, score_options,
                    hit_options, query, query_info, 1.0, &lookup_segments,
                    &filter_loc, &sbp, NULL);
    if (status != 0)
       goto bail_out;

    if (extra_returns) {
       if ((diagnostics = Blast_DiagnosticsInit()) == NULL)
          goto bail_out;
    }

    LookupTableWrapInit(query, lookup_options,
                        lookup_segments, sbp, &lookup_wrap, NULL);

    /* Initialize the HSPList collector stream. Results should not be sorted
       before reading from it. */
    hsp_stream = 
       Blast_HSPListCollectorInit(program_type, hit_options, 
                                  query_info->num_queries, FALSE);
    /* finally, do the search */

    status = BLAST_SearchEngine(program_type, query, query_info,
                seq_src, sbp, score_options, lookup_wrap, word_options, 
                ext_options, hit_options, eff_len_options,
                psi_options, db_options, hsp_stream, diagnostics, &results);

    hsp_stream = BlastHSPStreamFree(hsp_stream);
    if (status != 0)
       goto bail_out;

    /* Convert results to the SeqAlign form */

    status = BLAST_ResultsToSeqAlign(program_type, results, query_seqloc, 
                seq_src, score_options->gapped_calculation,
                score_options->is_ooframe, seqalign_out);

bail_out:
    AdjustOffSetsInSeqAlign(*seqalign_out, query_seqloc, subject_seqloc);

    if (!status && extra_returns) {
       Blast_SummaryReturnFill(program_type, score_options, sbp, 
                               lookup_options, word_options, ext_options, 
                               hit_options, eff_len_options, query_info, 
                               seq_src, diagnostics, extra_returns);
    }

    if (filter_out) {
       *filter_out = 
          BlastMaskLocToSeqLoc(program_type, filter_loc, query_seqloc);
    }

    Blast_DiagnosticsFree(diagnostics);
    BlastMaskLocFree(filter_loc);
    BlastSeqSrcFree(seq_src);
    LookupTableWrapFree(lookup_wrap);
    ListNodeFreeData(lookup_segments);
    Blast_HSPResultsFree(results);
    BlastSequenceBlkFree(query);
    BlastQueryInfoFree(query_info);
    BlastScoreBlkFree(sbp);
    LookupTableOptionsFree(lookup_options);
    BlastQuerySetUpOptionsFree(query_options);
    BlastExtensionOptionsFree(ext_options);
    BlastHitSavingOptionsFree(hit_options);
    BlastInitialWordOptionsFree(word_options);
    BlastScoringOptionsFree(score_options);
    BlastEffectiveLengthsOptionsFree(eff_len_options);
    PSIBlastOptionsFree(psi_options);
    BlastDatabaseOptionsFree(db_options);

    return status;
}


