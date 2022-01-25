static char const rcsid[] = "$Id: blast_seq.c,v 1.52 2004/10/06 18:16:17 dondosha Exp $";
/*
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

File name: blast_seq.c

Author: Ilya Dondoshansky

Contents: Functions converting between SeqLocs and structures used in BLAST.

******************************************************************************
 * $Revision: 1.52 $
 * */

#include <seqport.h>
#include <sequtil.h>
#include <objloc.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/core/blast_filter.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_encoding.h>

BlastSeqLoc* BlastSeqLocFromSeqLoc(SeqLocPtr mask_slp)
{
   BlastSeqLoc* last_loc = NULL,* head_loc = NULL;

   if (mask_slp == NULL)
      return NULL;

   if (mask_slp->choice == SEQLOC_PACKED_INT)
      mask_slp = (SeqLocPtr) mask_slp->data.ptrvalue;

   for ( ; mask_slp; mask_slp = mask_slp->next) {
      SeqIntPtr si = (SeqIntPtr) mask_slp->data.ptrvalue;
      if (!head_loc) {
         last_loc = head_loc = BlastSeqLocNew(&last_loc, si->from, si->to);
      } else {
         last_loc = BlastSeqLocNew(&last_loc, si->from, si->to);
      }
   }
   return head_loc;
}


SeqLocPtr BlastMaskLocToSeqLoc(EBlastProgramType program_number, 
                               const BlastMaskLoc* mask_loc, 
                               const SeqLoc* seq_loc)
{
   SeqLocPtr mask_head = NULL, last_mask = NULL;
   Int4 index;
   const Boolean k_translate = (program_number == eBlastTypeBlastx || 
                                program_number == eBlastTypeTblastx ||
                                program_number == eBlastTypeRpsTblastn);
   const Uint1 k_num_frames = (k_translate ? NUM_FRAMES : 1);
   SeqLoc* slp;

   if (mask_loc == NULL || mask_loc->seqloc_array == NULL)
      return NULL;

   for (index=0, slp = (SeqLoc*)seq_loc; slp; ++index, slp = slp->next)
   {
      Int4 frame_index = index*k_num_frames;
      Int4 tmp_index;
      for (tmp_index=frame_index; tmp_index<(frame_index+k_num_frames); tmp_index++)
      {
         BlastSeqLoc* loc = NULL;
         SeqIdPtr seqid = SeqLocId(slp);
         SeqLocPtr mask_slp_head = NULL, mask_slp_last = NULL;
         for (loc = mask_loc->seqloc_array[tmp_index]; loc; loc = loc->next)
         {
            SSeqRange* di = loc->ssr;
            SeqIntPtr si = SeqIntNew();
            si->from = di->left;
            si->to = di->right;
            si->id = SeqIdDup(seqid);
            if (!mask_slp_last)
               mask_slp_last = 
                  ValNodeAddPointer(&mask_slp_head, SEQLOC_INT, si);
            else 
               mask_slp_last = 
                  ValNodeAddPointer(&mask_slp_last, SEQLOC_INT, si);
         }

         if (mask_slp_head) {
            SeqLocPtr new_mask_slp = ValNodeAddPointer(NULL, SEQLOC_PACKED_INT, 
                                             mask_slp_head);
            /* Uint1 tmp_choice = (k_translate ? (frame_index+1) : 0); */
            Uint1 tmp_choice = (k_translate ? (tmp_index+1) : 0);
            /* The 'choice' of the SeqLoc in masks should show the frame,
               with values 1..6 when queries are translated; otherwise
               it does not matter. */
            if (!last_mask) {
               last_mask = ValNodeAddPointer(&mask_head, tmp_choice, new_mask_slp);
            } else {
               last_mask = ValNodeAddPointer(&last_mask, tmp_choice, new_mask_slp);
            }
        }
      }
   }
   return mask_head;
}

Int2 BlastMaskLocDNAToProtein(BlastSeqLoc* dna_seqloc, BlastMaskLoc* prot_maskloc, Int4 start, SeqLocPtr slp)
{

   Int4 dna_length;
   BlastSeqLoc** prot_seqloc;
   Int4 context;

   if (!dna_seqloc)
      return 1;

   prot_seqloc = &(prot_maskloc->seqloc_array[NUM_FRAMES*(start)]);

   dna_length = SeqLocLen(slp);
   /* Reproduce this mask for all 6 frames, with translated 
      coordinates */
   for (context = 0; context < NUM_FRAMES; ++context) {
       BlastSeqLoc* prot_head=NULL;
       BlastSeqLoc* seqloc_var;
       Int2 frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);

       prot_head = NULL;
       for (seqloc_var = dna_seqloc; seqloc_var; seqloc_var = seqloc_var->next) {
           Int4 from, to;
           SSeqRange* dip = seqloc_var->ssr;
           if (frame < 0) {
               from = (dna_length + frame - dip->right)/CODON_LENGTH;
               to = (dna_length + frame - dip->left)/CODON_LENGTH;
           } else {
               from = (dip->left - frame + 1)/CODON_LENGTH;
               to = (dip->right - frame + 1)/CODON_LENGTH;
           }
           BlastSeqLocNew(&prot_head, from, to);
       }
       prot_seqloc[context] = prot_head;
   }
   return 0;
}


Int2 BlastMaskLocProteinToDNA(BlastMaskLoc** mask_loc_ptr, SeqLocPtr slp)
{
   Int2 status = 0;
   Int4 index;
   BlastMaskLoc* mask_loc;
   BlastSeqLoc* loc;
   SSeqRange* dip;
   Int4 dna_length;
   Int2 frame;
   Int4 from, to;
   Int4 total;  /* total number of BlastSeqLoc's in array. */

   if (!mask_loc_ptr) 
      return -1;

   mask_loc = *mask_loc_ptr;
   total = mask_loc->total_size;

   index = NUM_FRAMES;

   for (index=0; index<total; index++)
   {
         dna_length = SeqLocLen(slp);
         frame = BLAST_ContextToFrame(eBlastTypeBlastx, index % NUM_FRAMES);

         for (loc = mask_loc->seqloc_array[index]; loc; loc = loc->next) {
            dip = loc->ssr;
            if (frame < 0) {
               to = dna_length - CODON_LENGTH*dip->left + frame;
               from = dna_length - CODON_LENGTH*dip->right + frame + 1;
            } else {
               from = CODON_LENGTH*dip->left + frame - 1;
               to = CODON_LENGTH*dip->right + frame - 1;
            }
            dip->left = from;
            dip->right = to;
         }
   }
   return status;
}

static Int4 BLAST_SetUpQueryInfo(SeqLocPtr slp, EBlastProgramType program, 
               BlastQueryInfo** query_info_ptr)
{
   Uint4 length, protein_length;
   Boolean translate = 
      (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
       program == eBlastTypeRpsTblastn);
   Boolean is_na = (program == eBlastTypeBlastn);
   Int2 num_frames, frame;
   Uint1 strand;
   BlastQueryInfo* query_info;
   Int4* context_offsets;
   Int4 index;
   Int4 total_contexts;
   Uint4 max_length = 0;

   if (translate)
      num_frames = NUM_FRAMES;
   else if (is_na)
      num_frames = 2;
   else
      num_frames = 1;

   if ((query_info = (BlastQueryInfo*) malloc(sizeof(BlastQueryInfo)))
       == NULL)
      return -1;

   query_info->first_context = 0;
   query_info->num_queries = ValNodeLen(slp);
   query_info->last_context = query_info->num_queries*num_frames - 1;
   total_contexts = query_info->last_context + 1;

   if ((strand = SeqLocStrand(slp)) == Seq_strand_minus) {
      if (translate)
         query_info->first_context = 3;
      else
         query_info->first_context = 1;
   }

   if ((context_offsets = (Int4*) 
      calloc((total_contexts+1), sizeof(Int4))) == NULL)
      return -1;

   if ((query_info->eff_searchsp_array = 
      (Int8*) calloc(total_contexts, sizeof(Int8))) == NULL)
      return -1;
   if ((query_info->length_adjustments =
        (Int4*) calloc(total_contexts, sizeof(Int4))) == NULL)
       return -1;

   context_offsets[0] = 0;

   query_info->context_offsets = context_offsets;
   
   /* Fill the context offsets */
   for (index = 0; slp; slp = slp->next, index += num_frames) {
      length = SeqLocLen(slp);  /* FIXME: could return -1 */
      strand = SeqLocStrand(slp);
      if (translate) {
         Int2 first_frame, last_frame;
         if (strand == Seq_strand_plus) {
            first_frame = 0;
            last_frame = 2;
         } else if (strand == Seq_strand_minus) {
            first_frame = 3;
            last_frame = 5;
         } else {
            first_frame = 0;
            last_frame = 5;
         }
         for (frame = 0; frame < first_frame; ++frame)
            context_offsets[index+frame+1] = context_offsets[index+frame];
         for (frame = first_frame; frame <= last_frame; ++frame) {
            protein_length = (length - frame%CODON_LENGTH)/CODON_LENGTH;
            max_length = MAX(max_length, protein_length);
            context_offsets[index+frame+1] = 
               context_offsets[index+frame] + protein_length + 1;
         }
         for ( ; frame < num_frames; ++frame)
            context_offsets[index+frame+1] = context_offsets[index+frame];
      } else {
         max_length = MAX(max_length, length);
         
         if (is_na) {
            if (strand == Seq_strand_plus) {
               context_offsets[index+1] = context_offsets[index] + length + 1;
               context_offsets[index+2] = context_offsets[index+1];
            } else if (strand == Seq_strand_minus) {
               context_offsets[index+1] = context_offsets[index];
               context_offsets[index+2] = 
                  context_offsets[index+1] + length + 1;
            } else {
               context_offsets[index+1] = context_offsets[index] + length + 1;
               context_offsets[index+2] = 
                  context_offsets[index+1] + length + 1;
            }
         } else {
            context_offsets[index+1] = context_offsets[index] + length + 1;
         }
      }
   }
   query_info->max_length = max_length;

   *query_info_ptr = query_info;
   return 0;
}

/** Given a SeqPort, fills a preallocated sequence buffer 
 * in the correct encoding.
 */
static Int2 SeqPortToSequenceBuffer(SeqPortPtr spp, Uint1 encoding, 
                                    Uint1** buffer)
{
   Uint1* buffer_var = *buffer;
   Uint1 residue;

   if (!buffer_var || !spp)
      return -1;

   SeqPortSet_do_virtual(spp, TRUE);

   switch (encoding) {
   case BLASTP_ENCODING: case NCBI4NA_ENCODING:
      while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
         if (IS_residue(residue)) {
           *(buffer_var++) = residue;
         }
      }
      break;
   case BLASTNA_ENCODING:
      while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
         if (IS_residue(residue)) {
            *(buffer_var++) = NCBI4NA_TO_BLASTNA[residue];
         }
      }
      break;
   default:
      break;
   }

   *buffer = buffer_var;
   return 0;
}

/** Fills sequence buffer for a single SeqLoc; 
 * fills both strands if necessary.
 */
static Int2 SeqLocFillSequenceBuffer(SeqLocPtr slp, Uint1 encoding, 
        Boolean add_sentinel_bytes, Boolean both_strands, Uint1* buffer)
{
   Uint1* buffer_var;
   SeqPortPtr spp = NULL;
   Uint1 sentinel = 
      (encoding == BLASTNA_ENCODING ? NCBI4NA_TO_BLASTNA[NULLB] : NULLB);
   Uint1 seq_code, strand;

   buffer_var = buffer;

   if (add_sentinel_bytes) {
      *buffer_var = sentinel;
      ++buffer_var;
   }

   if (encoding == BLASTP_ENCODING) {
      seq_code = Seq_code_ncbistdaa;
      strand = Seq_strand_unknown;
   } else {
      seq_code = Seq_code_ncbi4na;
      strand = SeqLocStrand(slp);
   }

   spp = SeqPortNewByLoc(slp, seq_code);
   SeqPortSet_do_virtual(spp, TRUE);

   SeqPortToSequenceBuffer(spp, encoding, &buffer_var);
   spp = SeqPortFree(spp);

   if (add_sentinel_bytes)
      *buffer_var = sentinel;

   if (both_strands && strand == Seq_strand_both) {
      SeqLocPtr tmp_slp=NULL;

      ++buffer_var;

      tmp_slp = SeqLocIntNew(SeqLocStart(slp), SeqLocStop(slp),
                             Seq_strand_minus, SeqLocId(slp));
            
      spp = SeqPortNewByLoc(tmp_slp, Seq_code_ncbi4na);
      SeqPortToSequenceBuffer(spp, encoding, &buffer_var);
      if (add_sentinel_bytes)
         *buffer_var = sentinel;

      spp = SeqPortFree(spp);
      SeqLocFree(tmp_slp);
   }

   return 0;
}

/** Find a genetic code string in ncbistdaa encoding, given an integer 
 * genetic code value.
 */
Int2 BLAST_GeneticCodeFind(Int4 gc, Uint1** genetic_code)
{
   ValNodePtr vnp;
   GeneticCodePtr gcp;
   char* gen_code_eaa = NULL;
   Uint1* gen_code_stdaa = NULL;
   Int4 gen_code_length = 0, index;
   SeqMapTablePtr smtp;

   gcp = GeneticCodeFind(gc, NULL);
   for (vnp = (ValNodePtr)gcp->data.ptrvalue; vnp != NULL; 
        vnp = vnp->next) {
      if (vnp->choice == 3) {  /* ncbieaa */
         gen_code_eaa = (char*)vnp->data.ptrvalue;
         break;
      }
   }

   if (!gen_code_eaa)
      return -1;
   smtp = SeqMapTableFind(Seq_code_ncbistdaa, Seq_code_ncbieaa);
   gen_code_length = StrLen(gen_code_eaa);
   *genetic_code = gen_code_stdaa = (Uint1*) calloc(gen_code_length+1, 1);

   if (!gen_code_stdaa)
      return -2;

   for (index = 0; index < gen_code_length; ++index) {
      gen_code_stdaa[index] = 
         SeqMapTableConvert(smtp, gen_code_eaa[index]);
   }
   
   return 0;
}

/** BLAST_GetSequence
 * Purpose:     Get the sequence for the BLAST engine, put in a Uint1 buffer
 * @param slp SeqLoc to extract sequence for [in]
 * @param query_info The query information structure, pre-initialized,
 *                   but filled here [in]
 * @param query_options Query setup options, containing the genetic code for
 *                      translation [in]
 * @param num_frames How many frames to get for this sequence? [in]
 * @param encoding In what encoding to retrieve the sequence? [in]
 * @param buffer_out Buffer to hold plus strand or protein [out]
 * @param buffer_length Length of buffer allocated [out]
 */
static Int2 
BLAST_GetSequence(SeqLocPtr slp, BlastQueryInfo* query_info, 
   const QuerySetUpOptions* query_options, Uint1 num_frames, Uint1 encoding, 
   Uint1* *buffer_out, Int4 *buffer_length)
{
   Int2		status=0; /* return value. */
   Int4 total_length; /* Total length of all queries/frames/strands */
   Int4		index; /* Loop counter */
   SeqLocPtr	slp_var; /* loop variable */
   Uint1*	buffer; /* buffer to fill. */
   Boolean add_sentinel_bytes = TRUE;
   Uint1* genetic_code=NULL;
   Boolean translate = FALSE;
   Int4 offset = 0;

   if (query_info) {
      *buffer_length = total_length = 
         query_info->context_offsets[query_info->last_context+1] + 1;
   } else {
      /* Subject sequence in 2 sequences comparison */
      *buffer_length = SeqLocLen(slp);
      if (encoding == NCBI4NA_ENCODING) {
         /* Searches with translated subjects (tblastn, tblastx) */
         add_sentinel_bytes = FALSE;
         total_length = *buffer_length;
      } else {
         total_length = (*buffer_length) + 2;
      }
   }

   if (num_frames == NUM_FRAMES) {
      /* Sequence must be translated in 6 frames. This can only happen
         for query - subject sequences are translated later. */
      Int4 gc;
      
      translate = TRUE;
      gc = (query_options ? query_options->genetic_code : 1);

      if ((status = BLAST_GeneticCodeFind(gc, &genetic_code)) != 0)
         return status;
   }

   *buffer_out = buffer = (Uint1 *) malloc((total_length)*sizeof(Uint1));
   
   for (index = 0, slp_var = slp; slp_var; 
        slp_var = slp_var->next, index += num_frames)
   {
      if (translate) {
         Uint1* na_buffer, *buffer_rev = NULL;
         Int4 context, context_start, context_end;
         Int4 na_length;
         Uint1 strand;
         

         na_length = SeqLocLen(slp_var);
         strand = SeqLocStrand(slp_var);
         /* Retrieve nucleotide sequence in an auxiliary buffer; 
            then translate into the appropriate place in the 
            preallocated buffer */
         if (strand == Seq_strand_plus) {
            na_buffer = (Uint1 *) malloc(na_length + 2);
            context_start = 0;
            context_end = 2;
         } else if (strand == Seq_strand_minus) {
            na_buffer = (Uint1 *) malloc(na_length + 2);
            context_start = 3;
            context_end = 5;
         } else {
            na_buffer = (Uint1*) malloc(2*na_length + 3);
            context_start = 0;
            context_end = 5;
         }
         SeqLocFillSequenceBuffer(slp_var, encoding, TRUE, TRUE, na_buffer);
         if (strand == Seq_strand_both)
            buffer_rev = na_buffer + na_length + 1;
	 else if (strand == Seq_strand_minus)
	    buffer_rev = na_buffer;

         for (context = context_start; context <= context_end; context++) {
            offset = query_info->context_offsets[index+context];
           
            BLAST_GetTranslation(na_buffer+1, buffer_rev, na_length, 
                BLAST_ContextToFrame(eBlastTypeBlastx, context), 
                &buffer[offset], genetic_code);
         }
         sfree(na_buffer);
      } else {
         /* This can happen both for query and subject, so query_info 
            might not be initialized here. */
         if (query_info)
            offset = query_info->context_offsets[index];
         SeqLocFillSequenceBuffer(slp_var, encoding, add_sentinel_bytes, 
                              (Boolean)(num_frames == 2), &buffer[offset]);
      }
      /* For subjects, do only one SeqLoc at a time */
      if (!query_info)
         break;
   }

   sfree(genetic_code);

   return status;
}

Int2 BLAST_SetUpQuery(EBlastProgramType program_number, 
        SeqLocPtr query_slp, const QuerySetUpOptions* query_options, 
        BlastQueryInfo** query_info, BLAST_SequenceBlk* *query_blk)
{
   Uint1* buffer;	/* holds sequence for plus strand or protein. */
   Int4 buffer_length;
   Int2 status;
   Uint1 num_frames;
   Uint1 encoding;

   if (query_slp == NULL || query_options == NULL ||
       query_info == NULL || query_blk == NULL)
      return -1;

   if ((status = BLAST_SetUpQueryInfo(query_slp, program_number, query_info)))
      return status;

   if (program_number == eBlastTypeBlastn) {
      encoding = BLASTNA_ENCODING;
      num_frames = 2;
   } else if (program_number == eBlastTypeBlastp ||
              program_number == eBlastTypeRpsBlast ||
              program_number == eBlastTypeTblastn) {
      encoding = BLASTP_ENCODING;
      num_frames = 1;
   } else { /* blastx or rpstblastn, which is also essentially blastx */
      encoding = NCBI4NA_ENCODING;
      num_frames = NUM_FRAMES;
   }

   if ((status=BLAST_GetSequence(query_slp, *query_info, query_options,
                  num_frames, encoding, &buffer, &buffer_length)))
      return status; 
        
   /* Do not count the first and last sentinel bytes in the 
      query length */
   if ((status=BlastSetUp_SeqBlkNew(buffer, buffer_length-2, 
                                    0, query_blk, TRUE)))
      return status;

   return 0;
}

Int2 BLAST_SetUpSubject(EBlastProgramType program_number, 
        SeqLocPtr subject_slp, BLAST_SequenceBlk** subject)
{
   Int2 status = 0;
   Uint1* subject_buffer = NULL; /* Buffer for the compressed subject 
                                      sequence in two sequences case */
   Int4 buffer_length=0; /* Length of subject sequence for two sequences 
                            case */
   Uint1 encoding;

   if (program_number == eBlastTypeBlastn)
      encoding = BLASTNA_ENCODING;
   else if (program_number == eBlastTypeTblastn ||
            program_number == eBlastTypeTblastx) {
      encoding = NCBI4NA_ENCODING;
   } else {
      encoding = BLASTP_ENCODING;
   }

   if ((status = BLAST_GetSequence(subject_slp, NULL, NULL, 1, encoding,
                                   &subject_buffer, &buffer_length)))
      return status;
   
   /* Initialize the sequence block, saving the sequence buffer in 
      'sequence_start'. */
   if ((status=BlastSetUp_SeqBlkNew(subject_buffer, buffer_length,
                                    0, subject, TRUE)))
      return status;

   /* If subject sequence is nucleotide, create compressed sequence buffer
      and save it in 'sequence'. For blastn, the sentinel bytes should not 
      be included in the packed sequence. */
   if (program_number == eBlastTypeBlastn)
      ++subject_buffer;

   if (program_number == eBlastTypeBlastn || 
       program_number == eBlastTypeTblastn || 
       program_number == eBlastTypeTblastx) {
      BLAST_PackDNA(subject_buffer, buffer_length, encoding, 
                    &((*subject)->sequence));
      (*subject)->sequence_allocated = TRUE;
   }

   return 0;
}
