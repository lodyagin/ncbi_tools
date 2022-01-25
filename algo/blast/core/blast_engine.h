/* $Id: blast_engine.h,v 1.41 2004/06/16 14:53:03 dondosha Exp $
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
 *
 * Author:  Ilya Dondoshansky
 *
 */

/** @file blast_engine.h
 * High level BLAST functions
 */

#ifndef __BLAST_ENGINE__
#define __BLAST_ENGINE__

#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_extend.h>
#include <algo/blast/core/blast_gapalign.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_diagnostics.h>   
#include <algo/blast/core/blast_hspstream.h>

#ifdef __cplusplus
extern "C" {
#endif

/** How many subject sequences to process in one database chunk. */
#define BLAST_DB_CHUNK_SIZE 1024

/** The high level function performing the BLAST search against a BLAST 
 * database after all the setup has been done.
 * @param program_number Type of BLAST program [in]
 * @param query The query sequence [in]
 * @param query_info Additional query information [in]
 * @param bssp Structure containing BLAST database [in]
 * @param sbp Scoring and statistical parameters [in]
 * @param score_options Hit scoring options [in]
 * @param lookup_wrap The lookup table, constructed earlier [in] 
 * @param word_options Options for processing initial word hits [in]
 * @param ext_options Options and parameters for the gapped extension [in]
 * @param hit_options Options for saving the HSPs [in]
 * @param eff_len_options Options for setting effective lengths [in]
 * @param psi_options Options specific to PSI-BLAST [in]
 * @param db_options Options for handling BLAST database [in]
 * @param results Structure holding all saved results [in] [out]
 * @param diagnostics Return statistics containing numbers of hits on 
 *                    different stages of the search [out]
 * @param results Results of the BLAST search [out]
 */
Int4 
BLAST_SearchEngine(Uint1 program_number, 
   BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
   const BlastSeqSrc* bssp, BlastScoreBlk* sbp, 
   const BlastScoringOptions* score_options, 
   LookupTableWrap* lookup_wrap, 
   const BlastInitialWordOptions* word_options, 
   const BlastExtensionOptions* ext_options, 
   const BlastHitSavingOptions* hit_options,
   const BlastEffectiveLengthsOptions* eff_len_options,
   const PSIBlastOptions* psi_options, 
   const BlastDatabaseOptions* db_options,
   BlastHSPStream* hsp_stream, BlastDiagnostics* diagnostics,
   BlastHSPResults** results);

/** The high level function performing an RPS BLAST search 
 * @param program_number Type of BLAST program [in]
 * @param query The query sequence [in]
 * @param query_info Additional query information [in]
 * @param bssp Structure containing BLAST database [in]
 * @param sbp Scoring and statistical parameters [in]
 * @param score_options Hit scoring options [in]
 * @param lookup_wrap The lookup table, constructed earlier [in] 
 * @param word_options Options for processing initial word hits [in]
 * @param ext_options Options and parameters for the gapped extension [in]
 * @param hit_options Options for saving the HSPs [in]
 * @param eff_len_options Options for setting effective lengths [in]
 * @param psi_options Options specific to PSI-BLAST [in]
 * @param db_options Options for handling BLAST database [in]
 * @param hsp_stream Placeholder for saving results [in]
 * @param diagnostics Return statistics containing numbers of hits on 
 *                    different stages of the search [out]
 * @param results Structure holding all saved results [in] [out]
 */
Int4 
BLAST_RPSSearchEngine(Uint1 program_number, 
   BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
   const BlastSeqSrc* bssp, BlastScoreBlk* sbp, 
   const BlastScoringOptions* score_options, 
   LookupTableWrap* lookup_wrap, 
   const BlastInitialWordOptions* word_options, 
   const BlastExtensionOptions* ext_options, 
   const BlastHitSavingOptions* hit_options,
   const BlastEffectiveLengthsOptions* eff_len_options,
   const PSIBlastOptions* psi_options, 
   const BlastDatabaseOptions* db_options,
   BlastHSPStream* hsp_stream, BlastDiagnostics* diagnostics, 
   BlastHSPResults** results);

/** Gapped extension function pointer type */
typedef Int2 (*BlastGetGappedScoreType) 
     (Uint1, BLAST_SequenceBlk*, BlastQueryInfo* query_info,
      BLAST_SequenceBlk*, BlastGapAlignStruct*, const BlastScoringParameters*,
      const BlastExtensionParameters*, const BlastHitSavingParameters*,
      BlastInitHitList*, BlastHSPList**, BlastGappedStats*);
     
/** Word finder function pointer type */
typedef Int2 (*BlastWordFinderType) 
     (BLAST_SequenceBlk*, BLAST_SequenceBlk*,
      LookupTableWrap*, Int4**, const BlastInitialWordParameters*,
      Blast_ExtendWord*, Uint4*, Uint4*, Int4, BlastInitHitList*,
      BlastUngappedStats*);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_ENGINE__ */
