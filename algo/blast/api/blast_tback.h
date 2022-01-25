/* $Id: blast_tback.h,v 1.1 2004/10/06 19:03:15 dondosha Exp $
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

/** @file blast_tback.h
 * API level functions to do perform traceback stage of the BLAST algorithm
 */

#ifndef __BLAST_TBACK__
#define __BLAST_TBACK__

#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_hspstream.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Perform traceback stage of the BLAST search, given the source of HSP lists,
 * obtained from the preliminary stage. The parameters internal to the engine
 * are calculated here independently of the similar calculation in the
 * preliminary stage, effectively making the two stages independent of each 
 * other.
 * @param program BLAST program type [in]
 * @param query Query sequence(s) structure [in]
 * @param query_info Additional query information [in]
 * @param seq_src Source of subject sequences [in]
 * @param score_options Scoring options [in]
 * @param ext_options Word extension options, needed for cutoff scores 
 *                    calculation only [in]
 * @param hit_options Hit saving options [in]
 * @param eff_len_options Options for calculating effective lengths [in]
 * @param db_options Database options (database genetic code) [in]
 * @param psi_options PSI BLAST options [in]
 * @param sbp Scoring block with statistical parameters and matrix [in]
 * @param hsp_stream Source of HSP lists. [in]
 * @param results Where to save the results after traceback. [out]
 */
Int2 
Blast_RunTracebackSearch(EBlastProgramType program, 
   BLAST_SequenceBlk* query, BlastQueryInfo* query_info, 
   const BlastSeqSrc* seq_src, const BlastScoringOptions* score_options,
   const BlastExtensionOptions* ext_options,
   const BlastHitSavingOptions* hit_options,
   const BlastEffectiveLengthsOptions* eff_len_options,
   const BlastDatabaseOptions* db_options, 
   const PSIBlastOptions* psi_options, BlastScoreBlk* sbp,
   BlastHSPStream* hsp_stream, BlastHSPResults** results);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_TBACK__ */
