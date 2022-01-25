/* $Id: blast_seq.h,v 1.17 2004/10/06 15:00:23 dondosha Exp $
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

File name: blast_seq.h

Author: Ilya Dondoshansky

Contents: Functions converting from SeqLocs to structures used in BLAST and 
          back.

Detailed Contents: 

******************************************************************************
 * $Revision: 1.17 $
 * */
#ifndef __BLAST_SEQ__
#define __BLAST_SEQ__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_C_TOOLKIT
#define NCBI_C_TOOLKIT
#endif

#include <objseq.h>
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_options.h>

#define NUM_FRAMES 6

/** Convert a BlastMaskLoc list to a list of SeqLocs, used for formatting 
 * BLAST results.
 * @param program_number identifies blastn, blastp, etc. [in]
 * @param mask_loc internal mask structure [in]
 * @param slp SeqLoc of query [in]
 * @return Pointer to SeqLoc
 */
SeqLocPtr 
BlastMaskLocToSeqLoc(EBlastProgramType program_number, 
                     const BlastMaskLoc* mask_loc, 
                     const SeqLoc* slp);

/* Converts a SeqLocPtr to a BlastSeqLoc, used for formatting.
 * @param mask_slp SeqLocPtr to be converted [in]
 * @return pointer to BlastSeqLoc
 */
BlastSeqLoc* BlastSeqLocFromSeqLoc(SeqLocPtr mask_slp);

/** Duplicate masks in 6 frames for each nucleotide sequence, converting
 * all masks' coordinates from DNA to protein.
 * @param mask_loc_ptr Masks list to be modified [in] [out]
 * @param slp List of nucleotide query SeqLoc's [in]
 */
Int2 BlastMaskLocDNAToProtein(BlastSeqLoc* dna_seqloc, BlastMaskLoc* prot_maskloc, Int4 start, SeqLocPtr slp);
/*
Int2 BlastMaskLocDNAToProtein(BlastMaskLoc** mask_loc_ptr, SeqLocPtr slp);
*/

/** Convert all masks' protein coordinates to nucleotide.
 * @param mask_loc_ptr Masks list to be modified [in] [out]
 * @param slp List of nucleotide query SeqLoc's [in]
 */
Int2 BlastMaskLocProteinToDNA(BlastMaskLoc** mask_loc_ptr, SeqLocPtr slp);

/** Given a list of query SeqLoc's, create the sequence block and the query
 * info structure. This is the last time SeqLoc is needed before formatting.
 * @param query_slp List of query SeqLoc's [in]
 * @param query_options Query setup options, containing genetic code for
 *                      translation [in]
 * @param program_number Type of BLAST program [in]
 * @param query_info Query information structure, containing offsets into 
 *                   the concatenated sequence [out]
 * @param query_blk Query block, containing (concatenated) sequence [out]
 */
Int2 BLAST_SetUpQuery(EBlastProgramType program_number, SeqLocPtr query_slp, 
        const QuerySetUpOptions* query_options, 
        BlastQueryInfo** query_info, BLAST_SequenceBlk* *query_blk);

/** Set up the subject sequence block in case of two sequences BLAST.
 * @param program_number Type of BLAST program [in]
 * @param subject_slp SeqLoc for the subject sequence [in]
 * @param subject Subject sequence block [out]
 */
Int2 BLAST_SetUpSubject(EBlastProgramType program_number, 
        SeqLocPtr subject_slp, BLAST_SequenceBlk** subject);

/** Find a genetic code string in ncbistdaa encoding, given an integer
 *  genetic code value.
 * @param gc genetic code value [in]
 * @param genetic_code genetic code string [out]
 */

Int2 BLAST_GeneticCodeFind(Int4 gc, Uint1** genetic_code);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_SEQ__ */
