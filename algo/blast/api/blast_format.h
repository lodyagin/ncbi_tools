/* $Id: blast_format.h,v 1.31 2004/10/04 14:00:18 madden Exp $
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

File name: blast_format.h

Author: Ilya Dondoshansky

Contents: Functions needed for formatting of BLAST results

******************************************************************************
 * $Revision: 1.31 $
 * */
#ifndef __BLAST_FORMAT__
#define __BLAST_FORMAT__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NCBI_C_TOOLKIT
#define NCBI_C_TOOLKIT
#endif

#include <ncbi.h>
#include <asn.h>
#include <bxmlobj.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_diagnostics.h>   
#include <algo/blast/api/twoseq_api.h>

/** Options for formatting BLAST results 
 */
typedef struct BlastFormattingOptions {
   Uint1 align_view;
   Uint4 align_options;
   Uint4 print_options;
   Boolean believe_query;
   Boolean html;
   void* outfp; /**< Output file: either FILE* or AsnIoPtr */
   Int4 number_of_descriptions;
   Int4 number_of_alignments;
   Boolean ungapped; /**< Should this be here????? */
} BlastFormattingOptions;

/** Initialize the formatting options structure.
 * @param program Number of the BLAST program [in]
 * @param file_out Name of the file where output is printed [in]
 * @param num_descriptions Number of definition lines to report per query [in]
 * @param num_alignments Number of alignments to show per query [in]
 * @param align_view What kind of formatted output to show? [in]
 * @param format_options_ptr The initialized structure [out]
*/
Int2 BlastFormattingOptionsNew(EBlastProgramType program, char* file_out, 
        Int4 num_descriptions, Int4 num_alignments, Int4 align_view,
        BlastFormattingOptions** format_options_ptr);

/** Deallocate memory for formatting options. In particular,
 * close the output file.
 * @param format_options Structure to free [in]
 */
BlastFormattingOptions* 
BlastFormattingOptionsFree(BlastFormattingOptions* format_options);

/** Print formatted output.
 * @param head Results in the SeqAlign form [in]
 * @param blast_database BLAST database name [in]
 * @param blast_program BLAST program name [in]
 * @param num_queries Number of query sequences [in]
 * @param query_slp Linked list of query SeqLocs [in]
 * @param mask_loc Masking locations for all queries [in]
 * @param format_options Formatting options [in]
 * @param is_ooframe Are frame shifts allowed in alignments? [in]
 */
Int2 BLAST_FormatResults(SeqAlignPtr head, char* blast_database,
        char* blast_program, Int4 num_queries, 
        SeqLocPtr query_slp, BlastMaskLoc* mask_loc, 
        const BlastFormattingOptions* format_options, Boolean is_ooframe, Int4** matrix,
        Blast_SummaryReturn* sum_returns);

/** Print the summary at the end of the BLAST report.
 * @param program_number Type of BLAST program [in]
 * @param format_options Formatting options [in]
 * @param rdfp Pointer to the BLAST database [in]
 * @param sum_returns infor from inside blast engine [in]
 */
Int2 Blast_PrintOutputFooter(EBlastProgramType program_number,
        const BlastFormattingOptions* format_options,
        ReadDBFILE* rdfp, const Blast_SummaryReturn* sum_returns);

/** Prints the top part of the traditional BLAST output, including version, 
 * reference(s) and database information.
 * @param format_options Formatting options [in]
 * @param is_megablast Is this a megablast search (i.e. greedy gapped
 *                     extension used)? [in]
 * @param program_name blastn, blastp, blastx, etc. [in]
 * @param dbname BLAST database name [in]
 * @param is_protein Is the database protein or nucleotide? [in]
 */
Int2 BLAST_PrintOutputHeader(const BlastFormattingOptions* format_options,
        Boolean is_megablast, const char* program_name, char* dbname, Boolean is_protein);

/** Prints the "log info" at the bottom of the megablast output, looks something like:
 * "Mega BLAST run finished, processed 2 queries"
 *
 * @param format_options Formatting options [in]
 * @param count number of searches completed [in]
 * @return 0 on success 
 */
Int2 BlastPrintLogReport(const BlastFormattingOptions* format_options, Int4 count);

void 
Blast_SeqIdGetDefLine(SeqIdPtr sip, char** buffer_ptr, Boolean ncbi_gi, Boolean accession_only);


#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_FORMAT__ */

