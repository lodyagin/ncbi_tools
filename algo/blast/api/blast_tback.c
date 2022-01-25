/* $Id: blast_tback.c,v 1.1 2004/10/06 19:03:15 dondosha Exp $
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
 * ===========================================================================
 *
 * Author: Ilya Dondoshansky
 *
 */

/** @file blast_tback.c
 * API level functions to perform the traceback stage of the BLAST algorithm
 */

static char const rcsid[] = 
    "$Id: blast_tback.c,v 1.1 2004/10/06 19:03:15 dondosha Exp $";

#include <algo/blast/api/blast_tback.h>
#include <algo/blast/core/blast_traceback.h>
#include <algo/blast/core/blast_setup.h>

Int2 
Blast_RunTracebackSearch(EBlastProgramType program, 
   BLAST_SequenceBlk* query, BlastQueryInfo* query_info, 
   const BlastSeqSrc* seq_src, const BlastScoringOptions* score_options,
   const BlastExtensionOptions* ext_options,
   const BlastHitSavingOptions* hit_options,
   const BlastEffectiveLengthsOptions* eff_len_options,
   const BlastDatabaseOptions* db_options, 
   const PSIBlastOptions* psi_options, BlastScoreBlk* sbp,
   BlastHSPStream* hsp_stream, BlastHSPResults** results)
{
   Int2 status = 0;
   BlastScoringParameters* score_params = NULL; /**< Scoring parameters */
   BlastExtensionParameters* ext_params = NULL; /**< Gapped extension 
                                                   parameters */
   BlastHitSavingParameters* hit_params = NULL; /**< Hit saving parameters*/
   BlastEffectiveLengthsParameters* eff_len_params = NULL; /**< Parameters
                                        for effective lengths calculations */
   BlastGapAlignStruct* gap_align = NULL; /**< Gapped alignment structure */

   status = 
      BLAST_GapAlignSetUp(program, seq_src, score_options, eff_len_options, 
         ext_options, hit_options, query_info, sbp, &score_params, 
         &ext_params, &hit_params, &eff_len_params, &gap_align);
   if (status)
      return status;

   /* Prohibit any subsequent writing to the HSP stream. */
   BlastHSPStreamClose(hsp_stream);

   status = 
      BLAST_ComputeTraceback(program, hsp_stream, query, query_info,
         seq_src, gap_align, score_params, ext_params, hit_params,
         eff_len_params, db_options, psi_options, results);

   /* Do not destruct score block here */
   gap_align->sbp = NULL;
   BLAST_GapAlignStructFree(gap_align);

   score_params = BlastScoringParametersFree(score_params);
   hit_params = BlastHitSavingParametersFree(hit_params);
   ext_params = BlastExtensionParametersFree(ext_params);
   eff_len_params = BlastEffectiveLengthsParametersFree(eff_len_params);
   return status;
}
