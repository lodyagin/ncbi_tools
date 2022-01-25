/* $Id: blast_tabular.h,v 1.5 2004/08/31 16:59:15 dondosha Exp $
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

File name: blast_tabular.h

Author: Ilya Dondoshansky

Contents: Functions needed for formatting of BLAST results

******************************************************************************
 * $Revision: 1.5 $
 * */
#ifndef __BLAST_TABULAR__
#define __BLAST_TABULAR__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_C_TOOLKIT
#define NCBI_C_TOOLKIT
#endif

#include <ncbi.h>
#include <asn.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/lookup_wrap.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_hspstream.h>
#include <algo/blast/core/blast_gapalign.h>
#include <objloc.h>

/** Tabular formatting options. */
typedef enum {
   eBlastTabularDefault=1,
   eBlastTabularAddSequences
} EBlastTabularFormatOptions;

/** Data structure containing all information necessary for production of the
 * tabular output.
 */
typedef struct BlastTabularFormatData {
   EBlastProgramType program; /**< Type of BLAST program */
   BlastHSPStream* hsp_stream; /**< Source of the BLAST results */
   BlastSeqSrc* seq_src; /**< Source of the subject sequences */
   BLAST_SequenceBlk* query; /**< Query sequence */
   BlastQueryInfo* query_info; /**< Query information, including context
                                  offsets and effective lengths. */
   BlastScoringParameters* score_params;
   BlastExtensionParameters* ext_params;
   BlastHitSavingParameters* hit_params;
   BlastEffectiveLengthsParameters* eff_len_params;
   Uint1* gen_code_string;
   BlastGapAlignStruct* gap_align;
   SeqLoc* query_slp; /**< Source of query sequences identifiers */
   FILE* outfp; /**< Output stream */
   Boolean perform_traceback; /**< Must gapped extension with traceback be
                                 performed before formatting? */
   Boolean show_gi; /**< Show gi's instead of full ids in output, if 
                       possible. */
   Boolean show_accession; /**< Show accessions instead of full ids in output,
                              if possible. This option has lower priority than
                              show_gi. */
   EBlastTabularFormatOptions format_options; /**< Tabular formatting options. */
} BlastTabularFormatData;

/** Function initializing the BlastTabularFormatData data structure fields. */
BlastTabularFormatData* 
Blast_TabularFormatDataInit(EBlastProgramType program, 
   BlastHSPStream* hsp_stream, BlastSeqSrc* seq_src, 
   BLAST_SequenceBlk* query, BlastQueryInfo* query_info,
   const BlastScoringOptions* scoring_options, BlastScoreBlk* sbp,
   const BlastEffectiveLengthsOptions* eff_len_options,
   const BlastExtensionOptions* ext_options,
   const BlastHitSavingOptions* hit_options,
   const BlastDatabaseOptions* db_options, SeqLoc* query_slp, FILE* outfp);

/** Free the tabular formatting data structure and all its internally 
 * allocated substructures. 
 */
void BlastTabularFormatDataFree(BlastTabularFormatData* tf_data);

/** Driver for the thread producing tabular output. */
void* Blast_TabularFormatThread(void* data);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_TABULAR__ */

