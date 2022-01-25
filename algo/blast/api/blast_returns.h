/* $Id: blast_returns.h,v 1.1 2004/05/14 17:19:03 dondosha Exp $
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

File name: blast_returns.h

Author: Ilya Dondoshansky

Contents: Manipulation of data returned from BLAST other than Seq-aligns

******************************************************************************
 * $Revision: 1.1 $
 * */
#ifndef __BLAST_RETURNS__
#define __BLAST_RETURNS__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_C_TOOLKIT
#define NCBI_C_TOOLKIT
#endif

#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_diagnostics.h>   
#include <algo/blast/api/twoseq_api.h>

typedef struct BLAST_DbSummary {
   struct BLAST_DbSummary* next;
   Boolean   is_protein;
   char*   name;
   char*   definition;
   char*   date;
   Int8   total_length;
   Int4   number_seqs;
   Boolean subset;      /* Print the subset message. */
} BLAST_DbSummary;

BLAST_DbSummary* LIBCALL
Blast_DbSummaryFree (BLAST_DbSummary* dbinfo);

/** Retrieves necessary information from a sequence source and fills the 
 * BLAST_DbSummary structure.
 */
BLAST_DbSummary* Blast_GetDbSummary(const BlastSeqSrc* seq_src);

/** Formats the BLAST parameters for the BLAST report.
 *	One char* is returned, newlines are indicated by tildes ('~').
 */	
char*
Blast_GetParametersBuffer(Uint1 program_number, 
   const BlastScoringOptions* score_options,
   const BlastScoreBlk* sbp, const LookupTableOptions* lookup_options,
   const BlastInitialWordOptions* word_options,
   const BlastExtensionOptions* ext_options,
   const BlastHitSavingOptions* hit_options,
   const BlastEffectiveLengthsOptions* eff_len_options,
   const BlastQueryInfo* query_info, const BlastSeqSrc* seq_src, 
   const BlastDiagnostics* diagnostics);

/** Fills the summary returns */
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
        BLAST_SummaryReturn** sum_returns_out);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_FORMAT__ */

