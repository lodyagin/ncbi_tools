/* $Id: blast_seqalign.h,v 1.19 2004/09/08 16:07:07 dondosha Exp $
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

File name: blast_seqalign.h

Author: Ilya Dondoshansky

Contents: Functions to convert BLAST results to the SeqAlign form

******************************************************************************
 * $Revision: 1.19 $
 * */
#ifndef __BLAST_SEQALIGN__
#define __BLAST_SEQALIGN__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_C_TOOLKIT
#define NCBI_C_TOOLKIT
#endif

#include <readdb.h>
#include <algo/blast/core/blast_hits.h>

/** This should be defined somewhere else!!!!!!!!!!!!! */
#define BLAST_SMALL_GAPS 0

/** Convert BLAST results structure to a list of SeqAlign's.
 * @param program_number Type of BLAST program [in]
 * @param results The BLAST results [in]
 * @param query_slp List of query SeqLoc's [in]
 * @param rdfp Pointer to a BLAST database structure [in]
 * @param subject_slp List of subject sequences locations [in]
 * @param is_gapped Is this a gapped alignment search? [in]
 * @param is_ooframe Is this a search with out-of-frame gapping? [in]
 * @param head_seqalign List of SeqAlign's [out]
 */
Int2 BLAST_ResultsToSeqAlign(EBlastProgramType program_number, 
        BlastHSPResults* results, SeqLocPtr query_slp, 
        ReadDBFILE* rdfp, SeqLoc* subject_slp, 
        Boolean is_gapped, Boolean is_ooframe, SeqAlignPtr* head_seqalign);

Boolean GapCollectDataForSeqalign(GapEditBlock* edit_block,
                                  GapEditScript* curr_in, Int4 numseg,
                                  Int4** start_out,
                                  Int4** length_out,
                                  Uint1** strands_out,
                                  Int4* start1, Int4* start2);

SeqAlignPtr LIBCALL GapEditBlockToSeqAlign(GapEditBlock* edit_block,
                                           SeqIdPtr subject_id,
                                           SeqIdPtr query_id);

/** Adjusts the offsets in a Seq-align list produced by a BLAST search
 * of multiple queries against multiple subjects. Seq-aligns in a list are
 * assumed to be sorted by query, and then by subject, i.e. all Seq-aligns for
 * the first query are at the beginning of the list, then all for the second
 * query, etc. Within a list for any single query, all Seq-aligns for one 
 * subject are grouped together. The order of queries in the Seq-align list 
 * must be the same as in the query Seq-loc list; the order of subjects
 * within a Seq-align list for a given query is the same as in the subject
 * Seq-loc list. Some or all queries or subjects might be skipped in the 
 * Seq-align list.
 */
void Blast_AdjustOffsetsInSeqAlign(SeqAlign* head, SeqLoc* query_slp,
				   SeqLoc* subject_slp);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_SEQALIGN__ */

