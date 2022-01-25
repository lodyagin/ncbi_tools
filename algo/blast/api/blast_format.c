/* $Id: blast_format.c,v 1.51 2004/06/07 18:40:34 dondosha Exp $
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
 * $Revision: 1.51 $
 * */

static char const rcsid[] = "$Id: blast_format.c,v 1.51 2004/06/07 18:40:34 dondosha Exp $";

#include <algo/blast/api/blast_format.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/api/blast_returns.h>
#include <sequtil.h>
#include <txalign.h>

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
   BioseqPtr query, ValNodePtr mask_loc);
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

extern void MBXmlClose(MBXml* mbxp, ValNodePtr other_returns, 
                       Boolean ungapped);

extern CharPtr LIBCALL BlastGetReference (Boolean html);

extern CharPtr LIBCALL BlastGetVersionNumber PROTO((void));
extern CharPtr LIBCALL BlastGetReleaseDate PROTO((void));

/** This function is defined in distrib/tools/blastool.c; 1st argument is
 * a pointer to TxDfDbInfo, defined in blastpri.h, which is identical to 
 * BLAST_DbSummary, defined in algo/blast/api/blast_returns.h
 */
extern Boolean LIBCALL 
PrintDbReport PROTO((BLAST_DbSummary* dbinfo, Int4 line_length, FILE *outfp));

/** The following 3 functions are from distrib/tools/readdb.h */
extern Boolean LIBCALL 
ReadDBBioseqFetchEnable PROTO((CharPtr program, CharPtr dbname, Boolean is_na,
                               Boolean now));

extern void LIBCALL ReadDBBioseqFetchDisable PROTO((void));

extern Boolean LIBCALL 
PrintDbInformation PROTO((CharPtr database, Boolean is_aa, Int4 line_length,
                          FILE *outfp, Boolean html));

Int2 BlastFormattingOptionsNew(Uint1 program_number, char* file_out, 
        Int4 num_descriptions, Int4 num_alignments, Int4 align_view, 
        BlastFormattingOptions** format_options_ptr)
{
   BlastFormattingOptions* format_options; 

   /* For translated RPS blast, the query will be a
      nucleotide sequence and subject sequences will
      be protein; thus the formatting must behave as
      if a blastx search was performed */

   if (program_number == blast_type_rpstblastn)
      program_number = blast_type_blastx;

   if ((format_options = (BlastFormattingOptions*)
        calloc(1, sizeof(BlastFormattingOptions))) == NULL)
      return -1;
   format_options->believe_query = FALSE;
   if (align_view != 7 && align_view < 10) {
      FILE* outfp;
      if ((outfp = FileOpen(file_out, "w")) == NULL)
         return -1;
      format_options->outfp = (void *) outfp;
   } else {
      AsnIoPtr aip;
      char write_mode[3];

      switch (align_view) {
      case 7:
         strcpy(write_mode, "wx"); break;
      case 10:
         strcpy(write_mode, "w"); break;
      case 11:
         strcpy(write_mode, "wb"); break;
      default:
         return -1;
      }
      
      if ((aip = AsnIoOpen(file_out, write_mode)) == NULL)
         return -1;
      format_options->outfp = (void *) aip;
   }

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
   if (format_options->align_view != 7 && format_options->align_view < 10)
      FileClose( (FILE *)format_options->outfp);
   else
      AsnIoClose( (AsnIoPtr)format_options->outfp);

   sfree(format_options);
   return NULL;
}

#define BXML_INCLUDE_QUERY 0x1

static BlastOutputPtr 
CreateBlastOutputHead(CharPtr program, CharPtr database, 
                      SeqLocPtr query_loc, Int4 flags)
{
    BlastOutputPtr boutp;
    Char buffer[1024];
    SeqPortPtr spp;
    Boolean is_aa = FALSE;
    Int4 i;
    
    if((boutp = BlastOutputNew()) == NULL)
        return FALSE;
    
    if (strcmp(program, "rpsblast") == 0)
       program = "blastp";
    else if (strcmp(program, "rpstblastn") == 0)
       program = "blastx";

    /* For optimization BLOSUM62 may be loaded ones */
    if (query_loc) {
       SeqIdPtr sip = SeqLocId(query_loc);
       BioseqPtr bsp;
       SeqIdWrite(sip, buffer, PRINTID_FASTA_LONG, sizeof(buffer));
       boutp->query_ID = StringSave(buffer);

       bsp = BioseqLockById(sip);
       if(!bsp ||
          (boutp->query_def = StringSave(BioseqGetTitle(bsp))) == NULL) {
          boutp->query_def = StringSave("No definition line found");
       }
       BioseqUnlock(bsp);

       boutp->query_len = SeqLocLen(query_loc);

       if(flags & BXML_INCLUDE_QUERY) {
          boutp->query_seq = MemNew(boutp->query_len+1);
          is_aa = (strcmp(program, "blastp") || strcmp(program, "tblastn"));
          spp = SeqPortNewByLoc(query_loc, 
                   (is_aa) ? Seq_code_ncbieaa : Seq_code_iupacna);
          
          for (i = 0; i < boutp->query_len; i++) {
             boutp->query_seq[i] = SeqPortGetResidue(spp);
          }
          spp = SeqPortFree(spp);
       } else {
          boutp->query_seq = NULL;    /* Do we need sequence here??? */
       }
    }
    /* Program name */
    boutp->program = StringSave(program);

    /* Database name */
    boutp->db = StringSave(database);

    /* Version text */
    sprintf(buffer, "%s %s [%s]", program, BlastGetVersionNumber(), 
            BlastGetReleaseDate());
    boutp->version = StringSave(buffer);

    /* Reference */
    boutp->reference = BlastGetReference(FALSE);

    /* Filling parameters */
    
    boutp->param = ParametersNew();    
    /* THE FOLLOWING IS NOT IMPLEMENTED. JUST FILL DEFAULT NUMBERS */
    boutp->param->expect = 10;
    boutp->param->gap_open = 0;
    boutp->param->gap_extend = 0;
    
    return boutp;
}

#define MACRO_atp_find(atp,name)\
        if((atp = AsnTypeFind(amp, #name))==NULL){\
                ErrPostEx(SEV_ERROR,0,0,\
                        "Could not find type <%s>", #name);\
                return NULL; \
        }

static MBXml* MBXmlInit(AsnIoPtr aip, CharPtr program, CharPtr database, 
                          SeqLocPtr query_loc, Int4 flags)
{
    MBXml* mbxp;
    AsnModulePtr amp;
    DataVal av;
    AsnTypePtr atp;
    Boolean retval = FALSE;

    AsnTypePtr       BLASTOUTPUT;
    AsnTypePtr       BLASTOUTPUT_program;
    AsnTypePtr       BLASTOUTPUT_version;
    AsnTypePtr       BLASTOUTPUT_reference;
    AsnTypePtr       BLASTOUTPUT_db;
    AsnTypePtr       BLASTOUTPUT_query_ID;
    AsnTypePtr       BLASTOUTPUT_query_def;
    AsnTypePtr       BLASTOUTPUT_query_len;
    AsnTypePtr       BLASTOUTPUT_query_seq;
    AsnTypePtr       BLASTOUTPUT_param;
    AsnTypePtr       BLASTOUTPUT_iterations;
    AsnTypePtr       BLASTOUTPUT_iterations_E;
    AsnTypePtr       BLASTOUTPUT_mbstat;

    if (strcmp(program, "rpsblast") == 0)
       program = "blastp";
    else if (strcmp(program, "rpstblastn") == 0)
       program = "blastx";

    mbxp = (MBXml*) MemNew(sizeof(MBXml));
    
    mbxp->aip = aip;
    
    if (! bxmlobjAsnLoad()) {
        return NULL;
    }
    
    amp = AsnAllModPtr();

    MACRO_atp_find(BLASTOUTPUT,BlastOutput);
    MACRO_atp_find(BLASTOUTPUT_program,BlastOutput.program);
    MACRO_atp_find(BLASTOUTPUT_version,BlastOutput.version);
    MACRO_atp_find(BLASTOUTPUT_reference,BlastOutput.reference);
    MACRO_atp_find(BLASTOUTPUT_db,BlastOutput.db);
    MACRO_atp_find(BLASTOUTPUT_query_ID,BlastOutput.query-ID);
    MACRO_atp_find(BLASTOUTPUT_query_def,BlastOutput.query-def);
    MACRO_atp_find(BLASTOUTPUT_query_len,BlastOutput.query-len);
    MACRO_atp_find(BLASTOUTPUT_query_seq,BlastOutput.query-seq);
    MACRO_atp_find(BLASTOUTPUT_param,BlastOutput.param);
    MACRO_atp_find(BLASTOUTPUT_iterations,BlastOutput.iterations);
    MACRO_atp_find(BLASTOUTPUT_iterations_E,BlastOutput.iterations.E);
    MACRO_atp_find(BLASTOUTPUT_mbstat,BlastOutput.mbstat);

    /* Start of iterations structure */
    mbxp->atp = BLASTOUTPUT_iterations_E;

    /* Head of all BlastOutput structure */
    mbxp->BlastOutput = BLASTOUTPUT;
    
    /* Head of iterations strucure */
    mbxp->BlastOutput_iterations = BLASTOUTPUT_iterations;

    /* Head of the final statistics for Mega BLAST */
    mbxp->BlastOutput_mbstat = BLASTOUTPUT_mbstat;
    
    mbxp->boutp = CreateBlastOutputHead(program, database, query_loc, flags);
    
    atp = AsnLinkType(NULL, BLASTOUTPUT);   /* link local tree */
    
    if (atp == NULL) {
        return NULL;
    }
    
    if (! AsnOpenStruct(mbxp->aip, atp, (Pointer) mbxp->boutp)) {
        return NULL;
    }

    if (mbxp->boutp->program != NULL) {
        av.ptrvalue = mbxp->boutp -> program;
        retval = AsnWrite(mbxp->aip, BLASTOUTPUT_program,  &av);
    }
    
    if (mbxp->boutp->version != NULL) {
        av.ptrvalue = mbxp->boutp->version;
        retval = AsnWrite(mbxp->aip, BLASTOUTPUT_version,  &av);
    }
    
    if (mbxp->boutp->reference != NULL) {
        av.ptrvalue = mbxp->boutp->reference;
        retval = AsnWrite(mbxp->aip, BLASTOUTPUT_reference,  &av);
    }

    if (mbxp->boutp -> db != NULL) {
        av.ptrvalue = mbxp->boutp->db;
        retval = AsnWrite(mbxp->aip, BLASTOUTPUT_db,  &av);
    }

    if (mbxp->boutp -> query_ID != NULL) {
        av.ptrvalue = mbxp->boutp->query_ID;
        retval = AsnWrite(mbxp->aip, BLASTOUTPUT_query_ID,  &av);
    }

    if (mbxp->boutp->query_def != NULL) {
        av.ptrvalue = mbxp->boutp->query_def;
        retval = AsnWrite(mbxp->aip, BLASTOUTPUT_query_def,  &av);
    }

    av.intvalue = mbxp->boutp->query_len;
    retval = AsnWrite(mbxp->aip, BLASTOUTPUT_query_len,  &av);

    if (mbxp->boutp->query_seq != NULL) {
        av.ptrvalue = mbxp->boutp->query_seq;
        retval = AsnWrite(mbxp->aip, BLASTOUTPUT_query_seq,  &av);
    }

    if (mbxp->boutp->param != NULL) {
        if (!ParametersAsnWrite(mbxp->boutp->param, 
                                mbxp->aip, BLASTOUTPUT_param)) {
            return NULL;
        }
    }

    if(!AsnOpenStruct(mbxp->aip, BLASTOUTPUT_iterations, NULL))
        return NULL;
    
    AsnIoFlush(mbxp->aip);
    
    return mbxp;
}

Int2 BLAST_FormatResults(SeqAlignPtr head, char* blast_database,
        char* blast_program, Int4 num_queries, 
        SeqLocPtr query_slp, BlastMaskLoc* blast_mask, 
        const BlastFormattingOptions* format_options, Boolean is_ooframe)
{  
   SeqAlignPtr seqalign = head, sap, next_seqalign;
   SeqLocPtr mask_loc, next_mask_loc = NULL, tmp_loc = NULL, mask_loc_head;
   Uint1 program_number;
   Uint1 align_type;
   Boolean query_is_na, db_is_na;
   Int4 query_index;
   SeqLocPtr slp, mask_slp;
   SeqIdPtr query_id;
   BioseqPtr bsp;
   AsnIoPtr aip = NULL;
   BlastPruneSapStruct* prune;
   SeqAnnotPtr seqannot;
   MBXml* mbxp = NULL;
   FILE *outfp = NULL;

   if (!format_options || !format_options->outfp || !query_slp)
      return -1;

   if (format_options->align_view != 7 && format_options->align_view < 10)
      outfp = (FILE *) format_options->outfp;
   else 
      aip = (AsnIoPtr) format_options->outfp; 

   if (strcmp(blast_program, "rpsblast") == 0)
      blast_program = "blastp";
   else if (strcmp(blast_program, "rpstblastn") == 0)
      blast_program = "blastx";

   BlastProgram2Number(blast_program, &program_number);

   /* Old toolkit might have different values for program numbers, so 
      use old toolkit function to determine alignment type. */

   align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);

   if (blast_database)
      ReadDBBioseqFetchEnable ("blast", blast_database, db_is_na, TRUE);

   query_index = 0;
   slp = query_slp;
   mask_loc_head = mask_loc = 
      BlastMaskLocToSeqLoc(program_number, blast_mask, slp);

   if (!seqalign && format_options->align_view < 7) {
      fprintf(outfp, "\n\n ***** No hits found ******\n\n");
   }
      
   if (format_options->align_view == 7) {
      mbxp = MBXmlInit(aip, "megablast", blast_database, slp, 0);
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
      } else if(format_options->align_view == 7) {
         /* XML printout - NOT IMPLEMENTED YET */
         IterationPtr iterp;
         
         iterp = BXMLBuildOneQueryIteration(seqalign, NULL, FALSE, 
                    format_options->ungapped, query_index, NULL, bsp, 
                    mask_loc);
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
         } else if (outfp) {
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
         }
         seqannot = SeqAnnotFree(seqannot);
      }
      seqalign = next_seqalign;
      /* Relink the mask locations so chain can be freed in the end.
       The 'tmp_loc' variable points to the location that was unlinked. */
      if (tmp_loc)
          tmp_loc->next = next_mask_loc;
      
      mask_loc = next_mask_loc;

   } /* End loop on seqaligns for different queries */

   if (format_options->align_view == 7 && mbxp) {
      MBXmlClose(mbxp, NULL, format_options->ungapped);
   }

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

Int2 PrintOutputFooter(Uint1 program_number, 
        const BlastFormattingOptions* format_options, 
        const BlastScoringOptions* score_options, const BlastScoreBlk* sbp,
        const LookupTableOptions* lookup_options,
        const BlastInitialWordOptions* word_options,
        const BlastExtensionOptions* ext_options,
        const BlastHitSavingOptions* hit_options,
        const BlastEffectiveLengthsOptions* eff_len_options,
        const BlastQueryInfo* query_info, const BlastSeqSrc* seq_src,
        const BlastDiagnostics* diagnostics) 
{
   FILE *outfp;
   BLAST_DbSummary* dbinfo_head,* dbinfo;
   char* params_buffer;

   if (!format_options || !format_options->outfp)
      return -1;

   if (format_options->align_view >= 7)
      return 0;

   outfp = format_options->outfp;

   if(format_options->html) 
      fprintf(outfp, "<PRE>\n");
   init_buff_ex(85);

   dbinfo_head = Blast_GetDbSummary(seq_src);

   for (dbinfo = dbinfo_head; dbinfo; dbinfo = dbinfo->next) {
      PrintDbReport(dbinfo, 70, outfp);
   }
   dbinfo_head = Blast_DbSummaryFree(dbinfo_head);
   
   if (sbp && sbp->kbp) {
      Blast_KarlinBlk* ka_params;
      
      if (sbp->kbp[query_info->first_context]) {
         ka_params = sbp->kbp[query_info->first_context];
         PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, 
                           outfp, FALSE);
      }

      if (score_options->gapped_calculation) {
         ka_params = sbp->kbp_gap[query_info->first_context];
         PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, 
                           outfp, TRUE);
      }
   }
   
   params_buffer = 
      Blast_GetParametersBuffer(program_number, score_options, sbp, 
                                lookup_options, word_options, ext_options, 
                                hit_options, eff_len_options, query_info, 
                                seq_src, diagnostics);

   PrintTildeSepLines(params_buffer, 70, outfp);
   sfree(params_buffer);
   free_buff();
   
   return 0;
}

Int2 BLAST_PrintOutputHeader(const BlastFormattingOptions* format_options,
        Boolean is_megablast, char* dbname, Boolean is_protein)
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
      if (is_megablast) {
         MegaBlastPrintReference(format_options->html, 90, 
                                 format_options->outfp);
         fprintf(format_options->outfp, "\n");
      }

      if(dbname && 
         !PrintDbInformation(dbname, is_protein, 70, 
             format_options->outfp, format_options->html))
         return 1;
              
      free_buff();
	
#ifdef OS_UNIX
      fprintf(format_options->outfp, "%s", "Searching\n\n");
#endif
   } else if (dbname && format_options->align_view == 9) {
      PrintTabularOutputHeader(dbname, NULL, NULL, "BLAST", 0, 
                               format_options->believe_query,
                               format_options->outfp);
   }
   return 0;
}

#define BUFFER_LENGTH 256

void 
Blast_SeqIdGetDefLine(SeqIdPtr sip, char* descr, char** buffer_ptr, 
   Boolean ncbi_gi, Boolean accession_only, Boolean seqid_only,
   Boolean believe_local_id)
{
   char* seqid_buffer = NULL, *title = NULL, *defline_buffer = NULL;
   Int4 gi;
   Boolean numeric_sip_type = FALSE;

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
   } else {
      if (!descr) {
         /* Get the title */
         BioseqPtr bsp = BioseqLockById(sip);
         if (bsp) {
            title = strdup(BioseqGetTitle(bsp));
            BioseqUnlock(bsp);
         }
      } else {
         title = descr;
      }

      if (seqid_only) 
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

/* This function might serve as a starting point for a callback function 
 * that prints results before the traceback stage, e.g. the on-the-fly 
 * tabular output, a la the -D3 option of the old megablast.
 */
void BLAST_PrintIntermediateResults(BlastHSPResults* results, 
        BlastQueryInfo* query_info, SeqLocPtr query_slp, 
        BlastSeqSrc* seq_src, BlastScoreBlk* sbp, 
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
   Blast_KarlinBlk* kbp;
   char bit_score_buff[10], eval_buff[10];
   char* eval_buff_ptr = NULL;
   BlastHSP* hsp;
   ListNode* seqid_wrap = NULL;

   if (!results || !query_info || !query_slp || !seq_src || 
       !sbp || !filename)
      return;

   outfp = FileOpen(filename, "w");

   for (query_index = 0, slp = query_slp; 
        query_index < results->num_queries && slp; 
        ++query_index, slp = slp->next) {
      hit_list = results->hitlist_array[query_index];
      query_id = SeqLocId(slp);
      Blast_SeqIdGetDefLine(query_id, NULL, &query_buffer, TRUE, FALSE, FALSE,
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
            BLASTSeqSrcGetSeqId(seq_src, (void*) &hsp_list->oid);
            if (seqid_wrap->choice == BLAST_SEQSRC_C_SEQID) {
               subject_id = (SeqId*) seqid_wrap->ptr;
               ListNodeFree(seqid_wrap);
            } else {
               /* Could not retrieve id for this subject sequence. This should
                  never happen, but if it does, skip this HSP. */
               continue;
            }
            Blast_SeqIdGetDefLine(subject_id, NULL, &subject_buffer, 
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

               eval_buff_ptr = eval_buff;

               ScoreAndEvalueToBuffers(hsp->bit_score, hsp->evalue, 
                  bit_score_buff, &eval_buff_ptr, TRUE);

               fprintf(outfp, "[%ld %ld] [%ld %ld] %s %s\n", 
                       (long)q_start, (long)q_end, (long)s_start, (long)s_end,
                       bit_score_buff, eval_buff_ptr);
            }
         }
      }
   }
   FileClose(outfp);
}
