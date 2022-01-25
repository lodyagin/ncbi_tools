/* $Id: blast_format.c,v 1.29 2003/10/23 20:15:37 dondosha Exp $
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

File name: blast_format.c

Author: Ilya Dondoshansky

Contents: Formatting of BLAST results (SeqAlign)

Detailed Contents: 

******************************************************************************
 * $Revision: 1.29 $
 * */

static char const rcsid[] = "$Id: blast_format.c,v 1.29 2003/10/23 20:15:37 dondosha Exp $";

#include "blast_format.h"
#include "blast_seq.h"
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <sequtil.h>
#include <txalign.h>
#include <readdb.h>

extern Uint1 LIBCALL 
BlastGetTypes PROTO((char* blast_program, Boolean* query_is_na, 
                     Boolean* db_is_na));
extern Boolean LIBCALL 
AcknowledgeBlastQuery PROTO((BioseqPtr bsp, Int4 line_length, FILE *outfp, 
                             Boolean believe_query, Boolean html));
extern void PrintTabularOutputHeader PROTO((char* blast_database, 
               BioseqPtr query_bsp, SeqLocPtr query_slp, char* blast_program, 
               Int4 iteration, Boolean believe_query, FILE *outfp));
extern void BlastPrintTabulatedResults PROTO((SeqAlignPtr seqalign, 
               BioseqPtr query_bsp, SeqLocPtr query_slp, Int4 num_alignments, 
               char* blast_program, Boolean is_ungapped, 
               Boolean believe_query, Int4 q_shift, Int4 s_shift, FILE *fp, 
               Boolean print_query_info));
extern IterationPtr 
BXMLBuildOneQueryIteration (SeqAlignPtr seqalign, ValNodePtr other_returns, 
   Boolean is_ooframe, Boolean ungapped, Int4 iter_num, char* message,
   BioseqPtr query);
extern BlastPruneSapStruct* 
BlastPruneHitsFromSeqAlign PROTO((SeqAlignPtr sap, Int4 number, 
   BlastPruneSapStruct* prune));
extern BlastPruneSapStruct* 
BlastPruneSapStructDestruct PROTO((BlastPruneSapStruct* prune));
extern Int4 SeqLocTotalLen (char* prog_name, SeqLocPtr slp);
extern Boolean LIBCALL 
PrintKAParameters PROTO((double Lambda, double K, double H, 
   Int4 line_length, FILE *outfp, Boolean gapped));
extern Boolean LIBCALL 
PrintTildeSepLines PROTO((char* buffer, Int4 line_length, FILE *outfp));
extern Boolean 
BlastPrintVersionInfo PROTO((char* program, Boolean html, FILE *outfp));
extern Boolean LIBCALL 
BlastPrintReference PROTO((Boolean html, Int4 line_length, FILE *outfp));
extern Boolean LIBCALL 
MegaBlastPrintReference PROTO((Boolean html, Int4 line_length, FILE *outfp));

typedef struct TxDfDbInfo {
   struct TxDfDbInfo* next;
   Boolean   is_protein;
   char*   name;
   char*   definition;
   char*   date;
   Int8   total_length;
   Int4   number_seqs;
   Boolean subset;      /* Print the subset message. */
} TxDfDbInfo, *TxDfDbInfoPtr;

extern Boolean LIBCALL 
PrintDbReport PROTO((TxDfDbInfo* dbinfo, Int4 line_length, FILE *outfp));

Int2 BlastFormattingOptionsNew(Uint1 program_number, char* file_out, 
        Int4 num_descriptions, Int4 num_alignments, Int4 align_view, 
        BlastFormattingOptions** format_options_ptr)
{
   BlastFormattingOptions* format_options; 
   if ((format_options = (BlastFormattingOptions*)
        calloc(1, sizeof(BlastFormattingOptions))) == NULL)
      return -1;
   format_options->believe_query = FALSE;
   if ((format_options->outfp = FileOpen(file_out, "w")) == NULL)
      return -1;
   format_options->number_of_descriptions = num_descriptions;
   format_options->number_of_alignments = num_alignments;
   
   format_options->print_options = 0;
   format_options->align_options = 0;
   format_options->align_options += TXALIGN_COMPRESS;
   format_options->align_options += TXALIGN_END_NUM;
   format_options->align_options += TXALIGN_SHOW_GI;
   format_options->print_options += TXALIGN_SHOW_GI;
   
   format_options->align_options += TXALIGN_MATRIX_VAL;
   format_options->align_options += TXALIGN_SHOW_QS;
   if (program_number == blast_type_blastx)
      format_options->align_options += TXALIGN_BLASTX_SPECIAL;
   
   format_options->align_view = align_view;

   *format_options_ptr = format_options;
   return 0;
}

BlastFormattingOptions* 
BlastFormattingOptionsFree(BlastFormattingOptions* format_options)
{
   FileClose(format_options->outfp);
   sfree(format_options);
   return NULL;
}

Int2 BLAST_FormatResults(SeqAlignPtr head, char* blast_database,
        char* blast_program, Int4 num_queries, 
        SeqLocPtr query_slp, BlastMask* blast_mask, 
        BlastFormattingOptions* format_options, Boolean is_ooframe)
{  
   SeqAlignPtr seqalign = head, sap, next_seqalign;
   SeqLocPtr mask_loc, next_mask_loc = NULL, tmp_loc = NULL, mask_loc_head;
   Uint1 align_type, program_number;
   Boolean db_is_na, query_is_na;
   Int4 query_index;
   SeqLocPtr slp, mask_slp;
   SeqIdPtr query_id;
   BioseqPtr bsp;
   AsnIoPtr aip = NULL;
   BlastPruneSapStruct* prune;
   SeqAnnotPtr seqannot;
   MBXml* mbxp = NULL;
   FILE *outfp;

   if (!format_options || !format_options->outfp || !query_slp)
      return -1;

   outfp = format_options->outfp;
   align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);
   BlastProgram2Number(blast_program, &program_number);
   
   if (blast_database)
      ReadDBBioseqFetchEnable ("blast", blast_database, db_is_na, TRUE);

   query_index = 0;
   slp = query_slp;
   mask_loc_head = mask_loc = 
      BlastMaskToSeqLoc(program_number, blast_mask, slp);

   if (!seqalign && format_options->align_view < 7) {
      fprintf(outfp, "\n\n ***** No hits found ******\n\n");
   }
      
   while (seqalign) {
      /* Find which query the current SeqAlign is for */
      query_id = TxGetQueryIdFromSeqAlign(seqalign);
      for ( ; slp && query_index < num_queries; 
            ++query_index, slp = slp->next) {
         if (SeqIdComp(query_id, SeqLocId(slp)) == SIC_YES)
            break;
      }
      /* The following should never happen! We mustn't have a SeqAlign left
         without a corresponding query */
      if (!slp || query_index >= num_queries)
         break;

      slp = slp->next;

      /* Find the masking location for this query */
      for ( ; mask_loc; mask_loc = mask_loc->next) {
         mask_slp = (SeqLocPtr) mask_loc->data.ptrvalue;
         if (SeqIdComp(query_id, SeqLocId(mask_slp)) == SIC_YES) {
            break;
         }
      }
      /* Unlink the masking location for this query and save the next one */
      if (mask_loc) {
         for (next_mask_loc = mask_loc; next_mask_loc->next; 
              next_mask_loc = next_mask_loc->next) {
            mask_slp = (SeqLocPtr) next_mask_loc->next->data.ptrvalue;
            if (SeqIdComp(query_id, SeqLocId(mask_slp))
                != SIC_YES) {
               break;
            }
         }
         tmp_loc = next_mask_loc;
         next_mask_loc = next_mask_loc->next;
         tmp_loc->next = NULL;
      }

      /* On the next iteration we can start from the next query */
      ++query_index;

      /* Find where this query's SeqAlign chain ends */
      for (sap = seqalign; sap->next; sap = sap->next) {
         if (SeqIdComp(query_id, TxGetQueryIdFromSeqAlign(sap->next)) != 
             SIC_YES)
            break;
      }
      /* Unlink this query's SeqAlign chain and save the start of the 
         next chain*/
      next_seqalign = sap->next;
      sap->next = NULL;

      /* Retrieve this query's Bioseq */
      bsp = BioseqLockById(query_id);

      if (format_options->align_view < 7) {
         init_buff_ex(70);
         AcknowledgeBlastQuery(bsp, 70, outfp, 
            format_options->believe_query, format_options->html);
         free_buff();
      }
      if (format_options->align_view == 8 || format_options->align_view == 9) {
         if (format_options->align_view == 9)
            PrintTabularOutputHeader(blast_database, bsp, NULL, blast_program,
               0, format_options->believe_query, outfp);
         
         BlastPrintTabulatedResults(seqalign, bsp, NULL, 
            format_options->number_of_alignments, blast_program, 
            format_options->ungapped, format_options->believe_query, 0, 0, 
            outfp, (format_options->align_view == 9));
         
         SeqAlignSetFree(seqalign);
      } else if(format_options->align_view == 7 || 
                format_options->align_view == 10 || 
                format_options->align_view == 11) {
         IterationPtr iterp;
         
         iterp = BXMLBuildOneQueryIteration(seqalign, NULL, FALSE, 
                    format_options->ungapped, query_index, NULL, bsp);
         IterationAsnWrite(iterp, mbxp->aip, mbxp->atp);
         AsnIoFlush(mbxp->aip);
         IterationFree(iterp);
         SeqAlignSetFree(seqalign);
      } else {
         seqannot = SeqAnnotNew();
         seqannot->type = 2;
         AddAlignInfoToSeqAnnot(seqannot, align_type);
         seqannot->data = seqalign;
         if (aip) {
            SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, aip, NULL);
            AsnIoReset(aip);
         }
         /* Uncacheing causes problems with ordinal nos. vs. gi's. */
         prune = BlastPruneHitsFromSeqAlign(seqalign, 
                    format_options->number_of_descriptions, NULL);
         ObjMgrSetHold();
         init_buff_ex(85);
         PrintDefLinesFromSeqAlign(prune->sap, 80, outfp, 
            format_options->print_options, FIRST_PASS, NULL);
         free_buff();
         
         prune = BlastPruneHitsFromSeqAlign(seqalign, 
                    format_options->number_of_alignments, prune);
         seqannot->data = prune->sap;
         if(is_ooframe) {
            OOFShowBlastAlignment(prune->sap, mask_loc, outfp, 
                                  format_options->align_options, NULL);
         } else if (format_options->align_view != 0) {
            ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, 
               format_options->align_options, NULL, mask_loc, NULL);
         } else {
            ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, 
               format_options->align_options, NULL, mask_loc, 
               FormatScoreFunc);
         }
         seqannot->data = seqalign;
         prune = BlastPruneSapStructDestruct(prune);
         ObjMgrClearHold();
         ObjMgrFreeCache(0);
         seqannot = SeqAnnotFree(seqannot);
      }

      seqalign = next_seqalign;
      /* Relink the mask locations so chain can be freed in the end.
       The 'tmp_loc' variable points to the location that was unlinked. */
      if (tmp_loc)
          tmp_loc->next = next_mask_loc;
      
      mask_loc = next_mask_loc;

   } /* End loop on seqaligns for different queries */

   /* Free the mask locations. Note that mask_loc is not a real SeqLoc, 
      but a ValNode wrapper around the actual masking locations. */
   for (mask_loc = mask_loc_head; mask_loc; ) {
      SeqLocSetFree((SeqLocPtr)mask_loc->data.ptrvalue);
      next_mask_loc = mask_loc->next;
      sfree(mask_loc);
      mask_loc = next_mask_loc;
   }
   
   if (blast_database)
      ReadDBBioseqFetchDisable();

   return 0;
}

static TxDfDbInfo* LIBCALL
TxDfDbInfoDestruct (TxDfDbInfo* dbinfo)
{
        TxDfDbInfo* next;

        if (dbinfo == NULL)
                return NULL;

        while (dbinfo)
        {
                sfree(dbinfo->name);
                sfree(dbinfo->definition);
                sfree(dbinfo->date);
                next = dbinfo->next;
                sfree(dbinfo);
                dbinfo = next;
        }

        return dbinfo;
}

static TxDfDbInfo* BLAST_GetDbInfo(ReadDBFILEPtr rdfp)
{
   TxDfDbInfo* dbinfo,* head = NULL,* dbinfo_var = NULL;
   char* chptr;

   while (rdfp) {
      dbinfo = calloc(1, sizeof(TxDfDbInfo));
      dbinfo->name = strdup(readdb_get_filename(rdfp));	
      
      if((chptr = readdb_get_title(rdfp)) == NULL)
         chptr = readdb_get_filename(rdfp);
      dbinfo->definition = strdup(chptr);	
      
      dbinfo->date = strdup(readdb_get_date(rdfp));	

      dbinfo->is_protein = readdb_is_prot(rdfp);
      
      if (rdfp->aliaslen)
         dbinfo->total_length = rdfp->aliaslen;
      else
         dbinfo->total_length = readdb_get_dblen(rdfp);
      if (rdfp->aliasnseq)
         dbinfo->number_seqs = rdfp->aliasnseq;
      else
         dbinfo->number_seqs = readdb_get_num_entries(rdfp);
      if (head == NULL) {
         head = dbinfo;
         dbinfo_var = dbinfo;
      } else {
         dbinfo_var->next = dbinfo;
         dbinfo_var = dbinfo_var->next;
      }
      rdfp = rdfp->next;
   }
   return head;
}
/*
        adds the new string to the buffer, separating by a tilde.
        Checks the size of the buffer for FormatBlastParameters and
        allocates longer replacement if needed.
*/

static Boolean
add_string_to_bufferEx(char* buffer, char* *old, Int2* old_length, Boolean add_tilde)

{
        char* new,* ptr;
        Int2 length, new_length;

        length = (StringLen(*old));

        if((Int2)(StringLen(buffer)+length+3) > *old_length)
        {
                new_length = *old_length + 255;
                new = calloc(new_length, sizeof(char));
                if (*old_length > 0 && *old != NULL)
                {
                        memcpy(new, *old, *old_length);
                        sfree(*old);
                }
                *old = new;
                *old_length = new_length;
        }

        ptr = *old;
        ptr += length;
        if (add_tilde)
        {
                *ptr = '~';
                ptr++;
        }

        while (*buffer != NULLB)
        {
                *ptr = *buffer;
                buffer++; ptr++;
        }

        return TRUE;
}

static Boolean
add_string_to_buffer(char* buffer, char* *old, Int2* old_length)

{
        return add_string_to_bufferEx(buffer, old, old_length, TRUE);
}



/*
	Formats the BLAST parameters for the BLAST report.
	One char* is returned, newlines are indicated by tildes ('~').
*/	


static char*
FormatBlastParameters(Uint1 program_number, 
   BlastScoringOptions* score_options,
   BlastScoreBlk* sbp, LookupTableOptions* lookup_options,
   BlastInitialWordOptions* word_options,
   BlastExtensionOptions* ext_options,
   BlastHitSavingOptions* hit_options,
   BlastQueryInfo* query_info, ReadDBFILEPtr rdfp, 
   BlastReturnStat* return_stats)
{
   Int4 cutoff = 0;
   char buffer[128];
   char* ret_buffer;
   Int2 ret_buffer_length;
   Int4 num_entries;
   Int8 total_length;
   Int4 qlen;
   double evalue;
   Boolean single_query = (query_info->last_context <= 1);
   BLAST_KarlinBlk* kbp;
   char* prog_name = NULL;
   
   ret_buffer = NULL;
   ret_buffer_length = 0;
   
   
   if (score_options->matrix) {
      sprintf(buffer, "Matrix: %s", score_options->matrix);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }

   if (hit_options->gapped_calculation) {
      sprintf(buffer, "Gap Penalties: Existence: %ld, Extension: %ld",
              (long) score_options->gap_open, 
              (long) score_options->gap_extend);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }
   
   if (rdfp) {
      readdb_get_totals_ex(rdfp, &total_length, &num_entries, TRUE);
   } else {
      num_entries = 1;
      total_length = 1000;
   }

   sprintf(buffer, "Number of Sequences: %ld", (long) num_entries);
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   if (return_stats) {
      sprintf(buffer, "Number of Hits to DB: %s", 
              Nlm_Int8tostr((Int8) return_stats->db_hits, 1));
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   
      sprintf(buffer, "Number of extensions: %ld", 
              (long) return_stats->init_extends);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      sprintf(buffer, "Number of successful extensions: %ld", 
              (long) return_stats->good_init_extends);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);

      if (hit_options->expect_value > 0.1) {
         sprintf(buffer, "Number of sequences better than %4.1f: %ld", 
                 hit_options->expect_value, 
                 (long) return_stats->number_of_seqs_better_E);
      } else {
         sprintf(buffer, "Number of sequences better than %3.1e: %ld", 
                 hit_options->expect_value, 
                 (long) return_stats->number_of_seqs_better_E);
      }
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   
      if (hit_options->gapped_calculation) {
         sprintf(buffer, 
                 "Number of HSP's better than %4.1f without gapping: %ld", 
                 hit_options->expect_value, 
                 (long) return_stats->prelim_gap_no_contest);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
         sprintf(buffer, 
                 "Number of HSP's successfully gapped in prelim test: %ld", 
                 (long) return_stats->prelim_gap_passed);
         add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      }
   }
   BlastNumber2Program(program_number, &prog_name);
   qlen = query_info->total_length;
   sfree(prog_name);
   sprintf(buffer, "length of query: %ld", (long)qlen);
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   sprintf(buffer, "length of database: %s", Nlm_Int8tostr (total_length, 1));
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   
   if (single_query) {
      sprintf(buffer, "effective search space used: %s",
         Nlm_Int8tostr(
            query_info->eff_searchsp_array[query_info->first_context], 1));
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }
   sprintf(buffer, "T: %ld", (long) lookup_options->threshold);
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   sprintf(buffer, "A: %ld", (long) word_options->window_size);
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   
   if (program_number == blast_type_blastn || 
       !hit_options->gapped_calculation)
      kbp = sbp->kbp[query_info->first_context];
   else
      kbp = sbp->kbp_gap[query_info->first_context];
   
   sprintf(buffer, "X1: %ld (%4.1f bits)", 
           (long)return_stats->x_drop_ungapped, word_options->x_dropoff);
   add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   if (hit_options->gapped_calculation) {
      sprintf(buffer, "X2: %ld (%4.1f bits)", 
              (long)return_stats->x_drop_gap, ext_options->gap_x_dropoff);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      sprintf(buffer, "X3: %ld (%4.1f bits)", 
              (long)return_stats->x_drop_gap_final,
              ext_options->gap_x_dropoff_final);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
      
      sprintf(buffer, "S1: %4.1f (%4.1f bits)", 
              return_stats->gap_trigger, ext_options->gap_trigger);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }
   
   cutoff = 0;
   if (single_query) {
      evalue = hit_options->expect_value;
      BLAST_Cutoffs_simple(&cutoff, &evalue, kbp,
         (double) query_info->eff_searchsp_array[query_info->first_context], 
         FALSE);
      sprintf(buffer, "S2: %ld (%4.1f bits)", (long) cutoff, 
              (((cutoff)*(kbp->Lambda))-(kbp->logK))/NCBIMATH_LN2);
      add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
   }
   return ret_buffer;
}

Int2 PrintOutputFooter(Uint1 program_number, 
        BlastFormattingOptions* format_options, 
        BlastScoringOptions* score_options, BlastScoreBlk* sbp,
        LookupTableOptions* lookup_options,
        BlastInitialWordOptions* word_options,
        BlastExtensionOptions* ext_options,
        BlastHitSavingOptions* hit_options,
        BlastQueryInfo* query_info, ReaddbNewArgs* readdb_args, 
        BlastReturnStat* return_stats) 
{
   FILE *outfp;
   TxDfDbInfo* dbinfo_head,* dbinfo;
   char* params_buffer;
   ReadDBFILEPtr rdfp = NULL;

   if (!format_options || !format_options->outfp)
      return -1;

   if (format_options->align_view >= 7)
      return 0;

   outfp = format_options->outfp;

   if(format_options->html) 
      fprintf(outfp, "<PRE>\n");
   init_buff_ex(85);

   if (readdb_args && 
       (rdfp = readdb_new(readdb_args->dbname, readdb_args->is_protein))) {
      dbinfo_head = BLAST_GetDbInfo(rdfp);

      for (dbinfo = dbinfo_head; dbinfo; dbinfo = dbinfo->next) {
         PrintDbReport(dbinfo, 70, outfp);
      }
      dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
   }
   
   if (sbp && sbp->kbp) {
      BLAST_KarlinBlk* ka_params;
      
      if (sbp->kbp[query_info->first_context]) {
         ka_params = sbp->kbp[query_info->first_context];
         PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, 
                           outfp, FALSE);
      }

      if (sbp->kbp_gap[query_info->first_context]) {
         ka_params = sbp->kbp_gap[query_info->first_context];
         PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, 
                           outfp, TRUE);
      }
   }
   
   params_buffer = 
      FormatBlastParameters(program_number, score_options, sbp, 
         lookup_options, word_options, ext_options, hit_options, query_info, 
         rdfp, return_stats);

   PrintTildeSepLines(params_buffer, 70, outfp);
   sfree(params_buffer);
   free_buff();
   
   readdb_destruct(rdfp);

   return 0;
}

/** Prints the top part of the traditional BLAST output, including version, 
 * reference(s) and database information.
 */
Int2 BLAST_PrintOutputHeader(BlastFormattingOptions* format_options,
        Boolean greedy_extension, ReaddbNewArgs* readdb_args)
{
   if (format_options->align_view < 7) {
      if (format_options->html) {
         fprintf(format_options->outfp, 
                 "<HTML>\n<TITLE>MEGABLAST Search Results</TITLE>\n");
         fprintf(format_options->outfp, 
                 "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
                 "VLINK=\"#660099\" ALINK=\"#660099\">\n");
         fprintf(format_options->outfp, "<PRE>\n");
      }
      init_buff_ex(90);
      BlastPrintVersionInfo("BLAST", format_options->html, 
                            format_options->outfp);
      fprintf(format_options->outfp, "\n");
      BlastPrintReference(format_options->html, 90, format_options->outfp);
      fprintf(format_options->outfp, "\n");
      if (greedy_extension) {
         MegaBlastPrintReference(format_options->html, 90, 
                                 format_options->outfp);
         fprintf(format_options->outfp, "\n");
      }

      if(readdb_args && 
         !PrintDbInformation(readdb_args->dbname, readdb_args->is_protein, 70, 
             format_options->outfp, format_options->html))
         return 1;
              
      free_buff();
	
#ifdef OS_UNIX
      fprintf(format_options->outfp, "%s", "Searching\n\n");
#endif
   } else if (readdb_args && format_options->align_view == 9) {
      PrintTabularOutputHeader(readdb_args->dbname, NULL,
         NULL, "BLAST", 0, format_options->believe_query,
         format_options->outfp);
   }
   return 0;
}

#define BUFFER_LENGTH 256

static void 
PrintSeqDefline(SeqIdPtr sip, char* descr, char** buffer_ptr, 
   Boolean ncbi_gi, Boolean accession_only, Boolean seqid_only,
   Boolean believe_local_id)
{
   char* seqid_buffer = NULL, *title = NULL, *defline_buffer;
   Int4 gi;
   Boolean numeric_sip_type = FALSE;

   if (!descr) {
      /* Get the title */
      BioseqPtr bsp = BioseqLockById(sip);
      title = strdup(BioseqGetTitle(bsp));
      BioseqUnlock(bsp);
   } else {
      title = descr;
   }

   if ((believe_local_id || sip->choice != SEQID_LOCAL) && 
       (sip->choice != SEQID_GENERAL ||
       StringCmp(((DbtagPtr)sip->data.ptrvalue)->db, "BL_ORD_ID"))) {
      if (!accession_only && sip->choice != SEQID_LOCAL) {
         seqid_buffer = (char*) malloc(BUFFER_LENGTH + 1);
         SeqIdWrite(sip, seqid_buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
      } else if (ncbi_gi) {
         numeric_sip_type = 
            GetAccessionFromSeqId(SeqIdFindBest(sip, SEQID_GI), 
                                  &gi, &seqid_buffer);
      } else {
         numeric_sip_type = 
            GetAccessionFromSeqId(SeqIdFindBestAccession(sip), 
                                  &gi, &seqid_buffer);
      }
   } else if (seqid_only && title) {
      seqid_buffer = StringTokMT(title, " \t\n\r", &title); 
   }
   
   if (numeric_sip_type && !seqid_buffer) {
      seqid_buffer = (char*) malloc(16);
      sprintf(seqid_buffer, "%ld", (long)gi);
   }

   if (!title || seqid_only) {
      defline_buffer = seqid_buffer;
   } else if (!seqid_buffer) {
      defline_buffer = title;
   } else {
      defline_buffer = (char *) malloc((StrLen(seqid_buffer) + StrLen(title) + 2));
      sprintf(defline_buffer, "%s %s", seqid_buffer, title);
      sfree(seqid_buffer);
      sfree(title);
   }
      
   *buffer_ptr = defline_buffer;
}

static void ScoreAndEvalueToBuffers(double bit_score, double evalue, 
               char* *bit_score_buf, char* *evalue_buf)
{
   if (evalue < 1.0e-180)
      sprintf(*evalue_buf, "0.0");
   else if (evalue < 1.0e-99)
      sprintf(*evalue_buf, "%2.0e", evalue);
   else if (evalue < 0.0009) 
         sprintf(*evalue_buf, "%3.1e", evalue);
   else if (evalue < 1.0) 
      sprintf(*evalue_buf, "%4.3f", evalue);
   else 
      sprintf(*evalue_buf, "%5.1f", evalue);
   
   if (bit_score > 9999)
      sprintf(*bit_score_buf, "%4.3e", bit_score);
   else if (bit_score > 99.9)
      sprintf(*bit_score_buf, "%4.1f", bit_score);
   else
      sprintf(*bit_score_buf, "%4.2f", bit_score);
}

/* This function might serve as a starting point for a callback function 
 * that prints results before the traceback stage, e.g. the on-the-fly 
 * tabular output, a la the -D3 option of the old megablast.
 */
void BLAST_PrintIntermediateResults(BlastResults* results, 
        BlastQueryInfo* query_info, SeqLocPtr query_slp, 
        ReadDBFILEPtr rdfp, SeqIdPtr seqid, BlastScoreBlk* sbp, 
        char* filename)
{
   Int4 query_index, hit_index, hsp_index;
   SeqIdPtr query_id, subject_id;
   char* query_buffer,* subject_buffer;
   BlastHitList* hit_list;
   BlastHSPList* hsp_list;
   Int4 query_length;
   Int4 q_start, q_end, s_start, s_end;
   FILE *outfp;
   SeqLocPtr slp;
   char* subject_descr;
   BLAST_KarlinBlk* kbp;
   char* bit_score_buff,* eval_buff;
   double bit_score;
   BlastHSP* hsp;

   if (!results || !query_info || !query_slp || (!rdfp && !seqid) || 
       !sbp || !filename)
      return;

   outfp = FileOpen(filename, "w");

   eval_buff = (char *) malloc(10);
   bit_score_buff = (char *) malloc(10);
   
   if (!rdfp) {
      /* Two sequences case */
      subject_id = seqid;
      subject_descr = NULL;
   }

   for (query_index = 0, slp = query_slp; 
        query_index < results->num_queries && slp; 
        ++query_index, slp = slp->next) {
      hit_list = results->hitlist_array[query_index];
      query_id = SeqLocId(slp);
      PrintSeqDefline(query_id, NULL, &query_buffer, TRUE, FALSE, FALSE,
                      FALSE);
      fprintf(outfp, "#Query = %s\n\n", query_buffer);
      sfree(query_buffer);

      if (!hit_list || hit_list->hsplist_count == 0) {
         fprintf(outfp, "# No hits found.\n\n");
      } else {
         query_length = SeqLocLen(query_slp);
         
         for (hit_index = 0; hit_index < hit_list->hsplist_count; 
              ++hit_index) {
            hsp_list = hit_list->hsplist_array[hit_index];
            if (rdfp) {
               readdb_get_descriptor(rdfp, hsp_list->oid, &subject_id,
                                     &subject_descr);
            }
            PrintSeqDefline(subject_id, subject_descr, &subject_buffer, 
                            TRUE, FALSE, FALSE, TRUE);
            fprintf(outfp, ">%s\n\n", subject_buffer);
            sfree(subject_buffer);

            for (hsp_index = 0; hsp_index < hsp_list->hspcnt; ++hsp_index) {
               hsp = hsp_list->hsp_array[hsp_index];
               if (hsp->query.offset > query_length) {
                  q_end = 2*query_length - hsp->query.offset + 1;
                  q_start = 2*query_length - hsp->query.end + 2;
                  s_start = hsp->subject.end;
                  s_end = hsp->subject.offset + 1;
               } else {
                  q_start = hsp->query.offset + 1;
                  q_end = hsp->query.end;
                  s_start = hsp->subject.offset + 1;
                  s_end = hsp->subject.end;
               }

               kbp = sbp->kbp[hsp->context];
               bit_score = (hsp->score*kbp->Lambda - kbp->logK) / 
                  NCBIMATH_LN2;

               ScoreAndEvalueToBuffers(bit_score, hsp->evalue, 
                  &bit_score_buff, &eval_buff);

               fprintf(outfp, "[%ld %ld] [%ld %ld] %s %s\n", 
                       (long)q_start, (long)q_end, (long)s_start, (long)s_end,
                       bit_score_buff, eval_buff);
            }
         }
      }
   }
   FileClose(outfp);
}

