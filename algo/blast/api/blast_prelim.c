/* $Id: blast_prelim.c,v 1.2 2004/10/06 15:00:44 dondosha Exp $
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

File name: blast_prelim.c

Author: Ilya Dondoshansky

Contents: Preliminary stage of a BLAST search performed by one of the threads
          in a multi-threaded search.

Detailed Contents: 

******************************************************************************
 * $Revision: 1.2 $
 * */

static char const rcsid[] = "$Id: blast_prelim.c,v 1.2 2004/10/06 15:00:44 dondosha Exp $";

#include <algo/blast/api/blast_prelim.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_gapalign.h>
#include <algo/blast/core/blast_setup.h>
#include <algo/blast/core/blast_engine.h>

BlastPrelimSearchThreadData* 
BlastPrelimSearchThreadDataInit(EBlastProgramType program,
   BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
   BlastSeqSrc* seq_src, LookupTableWrap* lut, 
   BlastScoringOptions* score_options, BlastInitialWordOptions* word_options,
   BlastExtensionOptions* ext_options, BlastHitSavingOptions* hit_options,
   BlastEffectiveLengthsOptions* eff_len_options,
   PSIBlastOptions* psi_options, BlastDatabaseOptions* db_options,
   BlastScoreBlk* sbp, BlastDiagnostics* diagnostics,
   BlastHSPStream* hsp_stream)
{
   BlastPrelimSearchThreadData* data = (BlastPrelimSearchThreadData*)
      calloc(1, sizeof(BlastPrelimSearchThreadData));
   
   data->program = program;
   data->query = query;
   data->query_info = BlastQueryInfoDup(query_info);
   data->seq_src = BlastSeqSrcCopy(seq_src);
   data->lut = lut;
   data->score_options = score_options;
   data->word_options = word_options;
   data->ext_options = ext_options;
   data->hit_options = hit_options;
   data->eff_len_options = eff_len_options;
   data->psi_options = psi_options;
   data->db_options = db_options;
   data->sbp = sbp;
   data->diagnostics = diagnostics;
   data->hsp_stream = hsp_stream;
   
   return data;
}


BlastPrelimSearchThreadData* 
BlastPrelimSearchThreadDataFree(BlastPrelimSearchThreadData* data)
{
   if (!data)
      return NULL;

   BlastSeqSrcFree(data->seq_src);
   BlastQueryInfoFree(data->query_info);
   sfree(data);
   return NULL;
}

void* Blast_PrelimSearchThreadRun(void* data)
{
   void* ret_status = NULL;
   Int2 status;
   BlastPrelimSearchThreadData* search_data = 
      (BlastPrelimSearchThreadData*) data;
   
   BlastScoringParameters* score_params = NULL; /**< Scoring parameters */
   BlastExtensionParameters* ext_params = NULL; /**< Gapped extension 
                                                   parameters */
   BlastHitSavingParameters* hit_params = NULL; /**< Hit saving parameters*/
   BlastEffectiveLengthsParameters* eff_len_params = NULL; /**< Parameters
                                        for effective lengths calculations */
   BlastGapAlignStruct* gap_align = NULL; /**< Gapped alignment structure */
   BlastDiagnostics* local_diagnostics = NULL;

   /* No need for a mutex in the local diagnostics structure. */
   local_diagnostics = Blast_DiagnosticsInit();

   status = 
      BLAST_GapAlignSetUp(search_data->program, search_data->seq_src, 
         search_data->score_options, search_data->eff_len_options, 
         search_data->ext_options, search_data->hit_options,
         search_data->query_info, search_data->sbp, &score_params, 
         &ext_params, &hit_params, &eff_len_params, &gap_align);
   if (status)
      return ret_status;

   BLAST_PreliminarySearchEngine(search_data->program,
      search_data->query, search_data->query_info, search_data->seq_src, 
      gap_align, score_params, search_data->lut, search_data->word_options, 
      ext_params, hit_params, eff_len_params, search_data->psi_options, 
      search_data->db_options, search_data->hsp_stream, 
      local_diagnostics);

   /* Do not destruct score block here */
   gap_align->sbp = NULL;
   gap_align = BLAST_GapAlignStructFree(gap_align);
   
   score_params = BlastScoringParametersFree(score_params);
   hit_params = BlastHitSavingParametersFree(hit_params);
   ext_params = BlastExtensionParametersFree(ext_params);
   eff_len_params = BlastEffectiveLengthsParametersFree(eff_len_params);
   
   Blast_DiagnosticsUpdate(search_data->diagnostics, local_diagnostics);
   
   Blast_DiagnosticsFree(local_diagnostics);

   BlastPrelimSearchThreadDataFree(search_data);

   return ret_status;
}
