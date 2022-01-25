/* $Id: blast_tabular.c,v 1.3 2004/06/14 20:43:30 dondosha Exp $
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

File name: blast_tabular.c

Author: Ilya Dondoshansky

Contents: On-the-fly tabular formatting of BLAST results

Detailed Contents: 

******************************************************************************
 * $Revision: 1.3 $
 * */

static char const rcsid[] = "$Id: blast_tabular.c,v 1.3 2004/06/14 20:43:30 dondosha Exp $";

#include <algo/blast/api/blast_tabular.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_engine.h>
#include <algo/blast/core/blast_traceback.h>
#include <algo/blast/api/blast_format.h>
#include <txalign.h>

BlastTabularFormatData* 
Blast_TabularFormatDataInit(Uint1 program, BlastHSPStream* hsp_stream,
   BlastSeqSrc* seq_src, BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
   const BlastScoringOptions* score_options, BlastScoreBlk* sbp,
   const BlastEffectiveLengthsOptions* eff_len_options,
   const BlastExtensionOptions* ext_options,
   const BlastHitSavingOptions* hit_options,
   const BlastDatabaseOptions* db_options, 
   SeqLoc* query_slp, FILE* outfp)
{
   BlastTabularFormatData* tf_data = 
      (BlastTabularFormatData*) calloc(1, sizeof(BlastTabularFormatData));
   tf_data->perform_traceback =
      (score_options->gapped_calculation && 
       ext_options->eTbackExt != eSkipTbck);

   tf_data->program = program;
   tf_data->hsp_stream = hsp_stream;
   tf_data->query = query;
   tf_data->gen_code_string = db_options->gen_code_string;
   tf_data->query_slp = query_slp;
   tf_data->outfp = outfp;
   /* Sequence source must be copied, to guarantee multi-thread safety. */
   tf_data->seq_src = BlastSeqSrcCopy(seq_src);
   /* Effective lengths must be duplicated in query info structure, because
      they might be changing in the preliminary search. */
   tf_data->query_info = BlastQueryInfoDup(query_info);

   /* If traceback will have to be performed before tabular output, 
      do the preparation for it here. */
   if (tf_data->perform_traceback) {
      BLAST_GapAlignSetUp(program, seq_src, 
         score_options, eff_len_options, ext_options, hit_options,
         query_info, sbp, &tf_data->score_params, &tf_data->ext_params, 
         &tf_data->hit_params, &tf_data->eff_len_params, &tf_data->gap_align);
   }

   return tf_data;
}

void BlastTabularFormatDataFree(BlastTabularFormatData* tf_data)
{
   /* Free only the structures that have been initialized internally */
   tf_data->query_info = BlastQueryInfoFree(tf_data->query_info);
   tf_data->score_params = BlastScoringParametersFree(tf_data->score_params);
   tf_data->ext_params = BlastExtensionParametersFree(tf_data->ext_params);
   tf_data->hit_params = BlastHitSavingParametersFree(tf_data->hit_params);
   tf_data->eff_len_params = 
      BlastEffectiveLengthsParametersFree(tf_data->eff_len_params);
   tf_data->gap_align = BLAST_GapAlignStructFree(tf_data->gap_align);
   tf_data->seq_src = BlastSeqSrcFree(tf_data->seq_src);
   sfree(tf_data);
}

/* This function might serve as a starting point for a callback function 
 * that prints results before the traceback stage, e.g. the on-the-fly 
 * tabular output, a la the -D3 option of the old megablast.
 */
void* Blast_TabularFormatThread(void* data) 
{
   BlastTabularFormatData* tf_data;
   Uint1 program;
   BlastHSPList* hsp_list = NULL;
   BlastSeqSrc* seq_src;
   BLAST_SequenceBlk* query = NULL; 
   BlastQueryInfo* query_info = NULL;
   BlastScoringParameters* score_params = NULL;
   BlastExtensionParameters* ext_params = NULL;
   BlastHitSavingParameters* hit_params = NULL;
   BlastEffectiveLengthsParameters* eff_len_params = NULL;
   Uint1* gen_code_string = NULL;
   BlastGapAlignStruct* gap_align = NULL;
   Int4 query_index, index;
   char* query_buffer = NULL;
   char* subject_buffer = NULL;
   Int4 q_start=0, q_end=0, s_start=0, s_end=0;
   SeqLocPtr slp;
   char bit_score_buff[10], eval_buff[10];
   char* eval_buff_ptr = NULL;
   BlastHSP* hsp;
   SeqId** query_id_array = NULL;
   Int4 align_length = 0;
   Int4 num_gaps = 0, num_gap_opens = 0, num_mismatches = 0;
   double perc_ident = 0;
   GetSeqArg seq_arg;
   Boolean one_seq_update_params;
 
   tf_data = (BlastTabularFormatData*) data;
   if (!tf_data || !tf_data->query_slp || !tf_data->hsp_stream ||
       !tf_data->seq_src || !tf_data->outfp) 
      return NULL;

   program = tf_data->program;
   seq_src = tf_data->seq_src;

   if (tf_data->perform_traceback) {
      query = tf_data->query;
      query_info = tf_data->query_info;
      score_params = tf_data->score_params;
      ext_params = tf_data->ext_params;
      hit_params = tf_data->hit_params;
      eff_len_params = tf_data->eff_len_params;
      gap_align = tf_data->gap_align;
      gen_code_string = tf_data->gen_code_string;

      memset((void*) &seq_arg, 0, sizeof(seq_arg));
      seq_arg.encoding = Blast_TracebackGetEncoding(program);
   }

   query_id_array = 
      (SeqId**) malloc(ValNodeLen(tf_data->query_slp)*sizeof(SeqId*));

   for (index = 0, slp = tf_data->query_slp; slp; ++index, slp = slp->next) {
      query_id_array[index] = SeqLocId(slp);
   }

   one_seq_update_params = (BLASTSeqSrcGetTotLen(seq_src) == 0);

   while (BlastHSPStreamRead(tf_data->hsp_stream, &hsp_list) 
          != kBlastHSPStream_Eof) {
      if (!hsp_list) {
         /* This should not happen, but just in case */
         continue;
      }

      /* Perform traceback if necessary */
      if (tf_data->perform_traceback) {
         BlastSequenceBlkClean(seq_arg.seq);
         seq_arg.oid = hsp_list->oid;
         if (BLASTSeqSrcGetSequence(seq_src, (void*) &seq_arg) < 0)
             continue;
         
         if (one_seq_update_params) {
            Int2 status;
            /* This is not a database search, so effective search spaces
               need to be recalculated based on this subject sequence length */
            if ((status = BLAST_OneSubjectUpdateParameters(program, 
                             seq_arg.seq->length, 
                             score_params->options, 
                             query_info, gap_align->sbp, 
                             ext_params, hit_params, NULL, 
                             eff_len_params)) != 0) {
               continue;
            }
         }

         Blast_TracebackFromHSPList(program, hsp_list, query,
            seq_arg.seq, query_info, gap_align, gap_align->sbp, score_params,
            ext_params->options, hit_params, gen_code_string);
         BLASTSeqSrcRetSequence(seq_src, (void*)&seq_arg);
         /* Recalculate the bit scores, since they might have changed. */
         Blast_HSPListGetBitScores(hsp_list, 
            score_params->options->gapped_calculation, gap_align->sbp);
      }
      subject_buffer = 
         BLASTSeqSrcGetSeqIdStr(seq_src, (void*) &hsp_list->oid); 

      for (index = 0; index < hsp_list->hspcnt; ++index) {
         hsp = hsp_list->hsp_array[index];
         query_index = 
            Blast_GetQueryIndexFromContext(hsp->context, program);
         Blast_SeqIdGetDefLine(query_id_array[query_index], NULL, 
                               &query_buffer, TRUE, FALSE, TRUE, FALSE);
         
         eval_buff_ptr = eval_buff;
         ScoreAndEvalueToBuffers(hsp->bit_score, hsp->evalue, 
                                 bit_score_buff, &eval_buff_ptr, FALSE);
         
         /* Calculate percentage of identities */
         Blast_HSPCalcLengthAndGaps(hsp, &align_length, &num_gaps, 
                                    &num_gap_opens);
         perc_ident = ((double)hsp->num_ident)/align_length * 100;
         num_mismatches = align_length - hsp->num_ident - num_gaps;
         
         Blast_HSPGetAdjustedOffsets(hsp, &q_start, &q_end, &s_start, &s_end);
         
         fprintf(tf_data->outfp, 
                 "%s\t%s\t%.2f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%s\t%s\n",
                 query_buffer, subject_buffer, perc_ident, 
                 (long) align_length, (long) num_mismatches, 
                 (long) num_gap_opens, (long) q_start, (long) q_end, 
                 (long) s_start, (long) s_end, eval_buff, bit_score_buff);
      }
      fflush(tf_data->outfp);
   }
   /* Deallocate the formatting thread data structure */
   BlastTabularFormatDataFree(tf_data);
   return NULL;
}
