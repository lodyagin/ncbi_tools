static char const rcsid[] = "$Id: blast_input.c,v 1.17 2004/09/13 16:36:09 madden Exp $";
/* ===========================================================================
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
* ===========================================================================*/

/*****************************************************************************

File name: blast_input.c

Author: Ilya Dondoshansky

Contents: Reading FASTA sequences for BLAST

$Revision: 1.17 $

******************************************************************************/
#include <objloc.h>
#include <tofasta.h>
#include <algo/blast/api/blast_input.h>
#include <algo/blast/api/blast_seq.h>
#include <algo/blast/core/blast_filter.h>

#define MAX_NUM_QUERIES 16383 /* == 1/2 INT2_MAX */
#define MAX_TOTAL_LENGTH 2000000

Int4
BLAST_GetQuerySeqLoc(FILE *infp, Boolean query_is_na, Uint1 strand, Int4 max_total_length,
   Int4 start, Int4 end, BlastMaskLoc** lcase_mask, SeqLocPtr* query_slp, 
   Int2Ptr ctr, Int4* num_queries, Boolean believe_query)
{
   Int4 total_length=0; /* total number of letters read this call, also final return value. */
   SeqEntryPtr sep;
   SeqLocPtr mask_slp, last_slp;
   char prefix[2];     /* for FastaToSeqEntryForDb */
   Int4 query_index = 0;  /* total number of sequence read. */
   ValNodePtr vnp=NULL; /* used to keep lower-case masking SeqLoc's */

   if (!query_slp)
   {
      ErrPostEx(SEV_FATAL, 0, 0, "NULL query_slp obtained in BLAST_GetQuerySeqLoc");
      return -1;
   }

   if (!ctr)  /* Not providing this can cause problems if multiple calls to this function are made. */
   {
      ErrPostEx(SEV_FATAL, 0, 0, "ctr must be non-NULL in BLAST_GetQuerySeqLoc");
      return -1;
   }

   if (max_total_length <= 0)
     max_total_length = MAX_TOTAL_LENGTH;

   *query_slp = NULL;
   last_slp = NULL;

   if (query_is_na && strand == Seq_strand_unknown)
      strand = Seq_strand_both;

   SeqMgrHoldIndexing(TRUE);
   mask_slp = NULL;
   if (lcase_mask) /* Make sure we don't get old (possibly freed) locations. */
     *lcase_mask = NULL;
   
   StrCpy(prefix, "");
   
   while ((sep=FastaToSeqEntryForDb(infp, query_is_na, NULL, believe_query, prefix, 
                                    ctr, (lcase_mask ? &mask_slp : NULL))) != NULL)
   {
      BioseqPtr query_bsp;
      Int4 from, to;

      if (lcase_mask)  /* Only keep if lcase masking is being read in. */
         ValNodeAddPointer(&vnp, 0, mask_slp);
      ++query_index;
      mask_slp = NULL;
      
      query_bsp = NULL;
      if (query_is_na) {
         SeqEntryExplore(sep, &query_bsp, FindNuc);
      } else {
         SeqEntryExplore(sep, &query_bsp, FindProt);
      }
      
      if (query_bsp == NULL) {
         ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
         return -1;
      }

      
      /* Original from and to are 1-offsets, except when they are 0's,
         in which case they are start and end of sequence respectively */
      from = ((start > 0) ? start - 1 : 0);
      to = ((end > 0) ? end - 1 : query_bsp->length - 1);

      from = MIN(from, query_bsp->length - 1);
      to = MIN(to, query_bsp->length - 1);

      if ((strand == Seq_strand_plus) || (strand == Seq_strand_minus) ||
          (from > 0) || (to < query_bsp->length - 1))
      {
         SeqLocPtr new_slp = SeqLocIntNew(from, to, strand, 
                                SeqIdFindBest(query_bsp->id, SEQID_GI)); 
         if (last_slp) {
            last_slp->next = new_slp;
            last_slp = last_slp->next;
         } else {
            *query_slp = last_slp = new_slp;
         }
      } else {
         last_slp = ValNodeAddPointer(&last_slp, SEQLOC_WHOLE, 
                       SeqIdSetDup(query_bsp->id));
         if (*query_slp == NULL)
            *query_slp = last_slp;
      }

      total_length += query_bsp->length;
      if (total_length > max_total_length || query_index >= MAX_NUM_QUERIES) {
         break;  /* Read maximum allowed amount of data. */
      }
   }

   if (lcase_mask)
   {
       BlastMaskLoc* mask_loc = BlastMaskLocNew(query_index);
       Int4 tmp_index=0;
       SeqLocPtr tmp_slp=NULL;
       ValNodePtr vnp_var = vnp;
       while (vnp_var)
       {
           tmp_slp = (SeqLocPtr) vnp_var->data.ptrvalue;
           if (tmp_slp)
               mask_loc->seqloc_array[tmp_index] = BlastSeqLocFromSeqLoc(tmp_slp);
           tmp_slp = SeqLocSetFree(tmp_slp);
           tmp_index++;
           vnp_var = vnp_var->next;
       }
       vnp = ValNodeFree(vnp);
       *lcase_mask = mask_loc;
   }
   
   SeqMgrHoldIndexing(FALSE);

   if (num_queries)
      *num_queries = query_index;

   return total_length;
}

/** Probably can be removed */
Int4 MakeBlastSequenceBlkFromFasta(FILE *fasta_fp, BLAST_SequenceBlk* seq)
{
   BioseqPtr query_bsp;
   SeqEntryPtr query_sep;
   Boolean is_na = FALSE;
   Boolean believe_defline = TRUE;
   SeqLocPtr mask_slp = NULL;
   Int2 ctr=1;
   
   Uint1* sequence = NULL;
   
   query_sep=FastaToSeqEntryForDb(fasta_fp,
                               is_na, /* query is nucleotide? */
                               NULL, /* error message */
                               believe_defline, /* believe query defline? */
                               "", /* prefix for localid if not parseable */
                               &ctr, /* starting point for constructing a unique id */
                               &mask_slp);

   if (query_sep == NULL)
      return 1;
   
   SeqEntryExplore(query_sep, &query_bsp, FindProt);
   
   if (query_bsp == NULL)
      return 1;
   
   /* allocate contiguous space for the sequence */
   sequence = (Uint1 *) malloc(query_bsp->length);
   
   if (sequence == NULL)
      return 1;
   
   /* convert to ncbistdaa encoding */
   BioseqRawConvert(query_bsp, Seq_code_ncbistdaa);
   
   /* read the sequence */
   BSSeek(query_bsp->seq_data, 0, SEEK_SET);
   BSRead(query_bsp->seq_data, sequence, query_bsp->length);
   
   seq->length = query_bsp->length;
   seq->sequence = sequence;
   return 0;
}

