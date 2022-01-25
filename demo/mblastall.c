/**************************************************************************
*                                                                         *
*                             COPYRIGHT NOTICE                            *
*                                                                         *
* This software/database is categorized as "United States Government      *
* Work" under the terms of the United States Copyright Act.  It was       *
* produced as part of the author's official duties as a Government        *
* employee and thus can not be copyrighted.  This software/database is    *
* freely available to the public for use without a copyright notice.      *
* Restrictions can not be placed on its present or future use.            *
*                                                                         *
* Although all reasonable efforts have been taken to ensure the accuracy  *
* and reliability of the software and data, the National Library of       *
* Medicine (NLM) and the U.S. Government do not and can not warrant the   *
* performance or results that may be obtained by using this software,     *
* data, or derivative works thereof.  The NLM and the U.S. Government     *
* disclaim any and all warranties, expressed or implied, as to the        *
* performance, merchantability or fitness for any particular purpose or   *
* use.                                                                    *
*                                                                         *
* In any work or product derived from this material, proper attribution   *
* of the author(s) as the source of the software or data would be         *
* appreciated.                                                            *
*                                                                         *
************************************************************************** 
 * $Revision 6.9$ *  
 * $Log: mblastall.c,v $
 * Revision 6.11  2000/04/25 22:04:51  dondosha
 * Report number of differences instead of score in -D [01] mode
 *
 * Revision 6.10  2000/04/24 16:45:58  dondosha
 * Check evalues of hsps in the callbacks; default expect value 0 means it is ignored
 *
 * Revision 6.9  2000/04/21 19:43:29  dondosha
 * Added command line option: total length of queries for a single search
 *
 * Revision 6.8  2000/04/12 22:48:47  dondosha
 * Fixed synchronization problem in writing segments information
 *
 * Revision 6.7  2000/04/12 21:19:20  dondosha
 * Cleaned argument list
 *
 * Revision 6.6  2000/04/12 18:29:52  dondosha
 * Added callback MegaBlastPrintSegments, made BlastSearchHandleResults static
 *
 * Revision 6.5  2000/04/07 16:54:58  dondosha
 * Added option for segments only output; moved result handling callbacks from 
 * mblast.c
 *
 * Revision 6.4  2000/03/31 19:13:41  dondosha
 * Changed some names related to MegaBlast
 *
 * Revision 6.3  2000/02/17 17:42:58  dondosha
 * Added extra parameter to FastaToSeqEntryForDb call due to prototype change;
 * do not fill mask_loc to avoid sequence substitution when filtering
 *
 * Revision 6.2  2000/02/03 22:02:53  dondosha
 * Fixed some memory leaks, added option for one-line results printout
 *
 * Revision 6.1  2000/02/01 21:41:14  dondosha
 * Initial revision
 *
 *
*/

#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <simutil.h>
#include <txalign.h>
#include <gapxdrop.h>
#include <sqnutils.h>
#include <mbalign.h>
#include <mblast.h>

#define DEFLINE_BUF 255


/* Used by the callback function. */
FILE *global_fp=NULL;
/*
	Callback to print out ticks, in UNIX only due to file systems
	portability issues.
*/
static int LIBCALLBACK
dummy_callback(Int4 sequence_number, Int4 number_of_positive_hits)
{
   return 0;
}

static int LIBCALLBACK
tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{

#ifdef OS_UNIX
   
   fprintf(global_fp, "%s", ".");
   fflush(global_fp);
#endif
	return 0;
}

enum {
   MBLAST_ENDPOINTS = 0,
   MBLAST_SEGMENTS,
   MBLAST_ALIGNMENTS
};

#define BUFFER_LENGTH 255

static int LIBCALLBACK
BlastSearchHandleResults(VoidPtr ptr)
{
   BlastSearchBlkPtr search = (BlastSearchBlkPtr) ptr;
   CharPtr subject_descr;
   SeqIdPtr sip, query_id, best_sip;
   Char subject_buffer[BUFFER_LENGTH], query_buffer[BUFFER_LENGTH];
   Int4 query_offset, query_end, query_length, s_start, s_end;
   Int4 index, hsp_index;
   BLAST_HSPPtr hsp; 
   Int2 context, descr_len;
   Char context_sign;
   BioseqPtr bsp;

   Uint4 header_index = 0;
   Int4 subject_gi, n_diff;
   
   if (search->current_hitlist == NULL || search->current_hitlist->hspcnt <= 0) {
      search->subject_info = BLASTSubjectInfoDestruct(search->subject_info);
      return 0;
   }
   readdb_get_descriptor(search->rdfp, search->subject_id, &sip,
			 &subject_descr);
   for (descr_len=0; subject_descr[descr_len]!=' ' && 
	   descr_len<StrLen(subject_descr); descr_len++);
   /* Find the correct query by comparing the hit start with query offsets */
   if (search->current_hitlist->hspcnt > 1)
      BlastSortUniqHspArray(search->current_hitlist);
   for (hsp_index=0; hsp_index<search->current_hitlist->hspcnt; hsp_index++) {
      hsp = search->current_hitlist->hsp_array[hsp_index];
      if (search->pbp->cutoff_e > 0 && 
	  hsp->evalue > search->pbp->cutoff_e)
	 continue;
      query_id = search->query_id;
      context = 1;
      while (search->query_context_offsets[context]<=
	     hsp->query.offset && query_id != NULL) {
	 context++;
	 if (context & 1)
	    query_id = query_id->next;
      }      
      if (query_id == NULL) /* Bad hsp, something wrong */
	 continue; 
      hsp->query.offset -= search->query_context_offsets[context-1];
      hsp->context = (context - 1) & 1;
      if (hsp->context) {
	 query_length = search->query_context_offsets[context] -
	    search->query_context_offsets[context-1] - 1;
	 hsp->query.end = query_length - hsp->query.offset;
	 hsp->query.offset = 
	    hsp->query.end - hsp->query.length + 1;
	 context_sign = '-'; 
      } else {
	 hsp->query.end = (++hsp->query.offset) + hsp->query.length - 1;
	 context_sign = '+';  
      }
      hsp->subject.offset++;
      SeqIdWrite(query_id, query_buffer, PRINTID_TEXTID_ACCESSION, BUFFER_LENGTH);
      if (context_sign == '+') {
	 s_start = hsp->subject.offset;
	 s_end = hsp->subject.offset+hsp->subject.length-1;
      } else {
	 s_start = hsp->subject.offset+hsp->subject.length-1;
	 s_end = hsp->subject.offset;
      }
      
      n_diff = ((hsp->subject.length + hsp->query.length)*
		search->sbp->reward / 2	- hsp->score) / 
	 (search->sbp->reward - search->sbp->penalty);

      fprintf(global_fp, "'%.*s'=='%c%s' (%d %d %d %d) %d\n", descr_len, 
	      subject_descr, context_sign, query_buffer, s_start, 
	      hsp->query.offset, s_end, hsp->query.end, n_diff);
   }
   sip = SeqIdSetFree(sip);
   fflush(global_fp);
   return 0;
}

#define LARGE_BUFFER_LENGTH 1024
static int LIBCALLBACK
MegaBlastPrintSegments(VoidPtr ptr)
{
   BlastSearchBlkPtr search = (BlastSearchBlkPtr) ptr;
   ReadDBFILEPtr rdfp = search->rdfp;
   BioseqPtr bsp;
   BLAST_HSPPtr hsp; 
   Int4 i, j, subject_gi;
   Int2 context;
   Char query_buffer[BUFFER_LENGTH];
   SeqIdPtr sip, query_id; 
   Int4 oid, descr_len;
   SeqPortPtr spp;
   Int4 index, hsp_index, n_diff;
   Uint1 residue;
   Uint1Ptr query_seq, subject_seq;
   FloatHi perc_ident;
   Boolean reverse;
   Char strand;
   GapXEditScriptPtr esp;
   Int4 q_start, q_end, s_start, s_end, query_length, subj_length, numseg;
   Int4Ptr length, start;
   Uint1Ptr strands;
   CharPtr subject_descr;
   CharPtr buffer;
   Char tmp_buffer[BUFFER_LENGTH];
   Int4 buffer_size, max_buffer_size = LARGE_BUFFER_LENGTH;

   buffer = (CharPtr) Malloc(LARGE_BUFFER_LENGTH);
   if (search->current_hitlist == NULL || search->current_hitlist->hspcnt <= 0) {
      search->subject_info = BLASTSubjectInfoDestruct(search->subject_info);
      return 0;
   }

   subj_length = readdb_get_sequence(rdfp, search->subject_id, &subject_seq);
   readdb_get_descriptor(search->rdfp, search->subject_id, &sip, &subject_descr);
   for (descr_len=0; subject_descr[descr_len]!=' ' && 
	   descr_len<StrLen(subject_descr); descr_len++);
   /* Find the correct query by comparing the hit start with query offsets */
   if (search->current_hitlist->hspcnt > 1)
      BlastSortUniqHspArray(search->current_hitlist);
   for (hsp_index=0; hsp_index<search->current_hitlist->hspcnt; hsp_index++) {
      hsp = search->current_hitlist->hsp_array[hsp_index];
      if (search->pbp->cutoff_e > 0 && 
	  hsp->evalue > search->pbp->cutoff_e)
	 continue;
      query_id = search->query_id;
      context = 1;
      while (search->query_context_offsets[context]<=
	     hsp->query.offset && query_id != NULL) {
	 context++;
	 if (context & 1)
	    query_id = query_id->next;
      }      
      if (query_id == NULL) /* Bad hsp, something wrong */
	 continue; 
      hsp->query.offset -= search->query_context_offsets[context-1];
      hsp->context = (context - 1) & 1;

      if (hsp->context) {
	 query_length = search->query_context_offsets[context] -
	    search->query_context_offsets[context-1] - 1;
	 s_start = hsp->subject.offset+hsp->subject.length;
	 s_end = hsp->subject.offset + 1;
	 q_end = query_length - hsp->query.offset;
	 q_start = q_end - hsp->query.length + 1;
	 strand = '-'; 
      } else {
	 strand = '+';  
	 s_start = hsp->subject.offset + 1;
	 s_end = hsp->subject.offset+hsp->subject.length;
	 q_start = hsp->query.offset + 1;
	 q_end = hsp->query.offset + hsp->query.length;
      }
      SeqIdWrite(query_id, query_buffer, PRINTID_TEXTID_ACCESSION,
		 BUFFER_LENGTH);
      n_diff = ((hsp->subject.length + hsp->query.length) *
		search->sbp->reward / 2	- hsp->score) / 
	 (search->sbp->reward - search->sbp->penalty);

      sprintf(buffer, "\n#'>%.*s'=='%c%s' (%d %d %d %d) %d\na {\n  s %d\n  b %d %d\n  e %d %d\n", 
	      descr_len, subject_descr, strand, query_buffer, 
	      s_start, q_start, s_end, q_end, n_diff, n_diff, 
	      s_start, q_start, s_end, q_end);
      buffer_size = StringLen(buffer);

      query_seq = search->context[context-1].query->sequence;

      esp = hsp->gap_info->esp;
        
      for (numseg=0; esp; esp = esp->next, numseg++);

      GXECollectDataForSeqalign(hsp->gap_info, hsp->gap_info->esp, numseg,
				&start, &length, &strands, 
				&hsp->query.offset, &hsp->subject.offset);

      for (index=0; index<numseg; index++) {
	 if (strand == '+') {
	    i = index;
	    q_start = start[2*i] + 1;
	    q_end = q_start + length[i] - 1;
	 } else {
	    i = numseg - 1 - index;
	    q_start = query_length - start[2*i];
	    q_end = q_start - length[i] + 1;
	 }
	 if (start[2*i] != -1 && start[2*i+1] != -1) {
	    perc_ident = MegaBlastGetPercentIdentity(query_seq, subject_seq, 
						     start[2*i], 
						     start[2*i+1], length[i],
						     FALSE);
	    sprintf(tmp_buffer, "  l %d %d %d %d (%.0f)\n", start[2*i+1]+1, 
		    q_start, start[2*i+1]+length[i],
		    q_end, perc_ident);	 
	    if ((buffer_size += StringLen(tmp_buffer)) > max_buffer_size - 1) {
	       max_buffer_size *= 2;
	       buffer = (CharPtr) Realloc(buffer, max_buffer_size);
	    }
	    StringCat(buffer, tmp_buffer);
	 }
      }
      StringCat(buffer, "}\n");
      fprintf(global_fp, buffer);
   } /* End loop on hsp's */
   MemFree(buffer);
   fflush(global_fp);
}

#define NUMARG 26

static Args myargs [NUMARG] = {
 { "Program Name",
   "blastn", NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL}, /* 0 */
  { "Database", 
    "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL}, /* 1 */
  { "Query File", 
	"stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},/* 2 */
  { "Expectation value", 
	"1000000.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},/* 3 */
  { "alignment view options:\n0 = pairwise,\n1 = master-slave showing identities,\n2 = master-slave no identities,\n3 = flat master-slave, show identities,\n4 = flat master-slave, no identities,\n5 = master-slave no identities and blunt ends,\n6 = flat master-slave, no identities and blunt ends", 
        "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},/* 4 */
  { "BLAST report Output File", 
	"stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},/* 5 */
  { "Filter query sequence (DUST with blastn, SEG with others)",
        "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},/* 6 */
  { "X dropoff value for gapped alignment (in bits)",
	"20", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},/* 7 */
  { "Show GI's in deflines",
        "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},/* 8 */
  { "Penalty for a nucleotide mismatch",
	"-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},/* 9 */
  { "Reward for a nucleotide match",
	"1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},/* 10 */
  { "Number of one-line descriptions (V)",
        "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},/* 11 */
  { "Number of alignments to show (B)",
        "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},/* 12 */
  { "Level of details in alignment output (0 - alignment endpoints, 1 - all ungapped segment endpoints, 2 - full alignments)",
        "0", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},/* 13 */
  { "Number of processors to use",
        "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},/* 14 */
  { "SeqAlign file", 
	NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},/* 15 */
  { "Believe the query defline",
        "T", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},/* 16 */
  { "Maximal total length of queries for a single search iteration",
        "2000000", NULL, NULL, FALSE, 'M', ARG_INT, 0.0, 0, NULL},/* 17 */
  { "Word size (length of best perfect match)", 
        "32", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},/* 18 */
  { "Effective length of the database (use zero for the real size)", 
        "0", NULL, NULL, FALSE, 'z', ARG_FLOAT, 0.0, 0, NULL},/* 19 */
  { "Number of best hits from a region to keep",
        "0", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},/* 20 */
  { "Length of region used to judge hits",
        "20", NULL, NULL, FALSE, 'L', ARG_INT, 0.0, 0, NULL},/* 21 */
  { "Effective length of the search space (use zero for the real size)",
        "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},/* 22 */
  { "Query strands to search against database: 3 is both, 1 is top, 2 is bottom",
        "3", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},/* 23 */
  { "Produce HTML output",
        "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},/* 24 */
  { "Restrict search of database to list of GI's",
	NULL, NULL, NULL, TRUE, 'l', ARG_STRING, 0.0, 0, NULL}/* 25 */
};

#define MAX_NUM_QUERIES 20000

Int2 Main (void)
 
{
	AsnIoPtr aip;
	BioseqPtr query_bsp, PNTR query_bsp_array;
	BioSourcePtr source;
	BLAST_MatrixPtr matrix;
	BLAST_OptionsBlkPtr options;
	BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
	BlastPruneSapStructPtr prune;
	Boolean db_is_na, query_is_na, show_gi, believe_query=FALSE;
	Boolean html=FALSE;
	CharPtr params_buffer=NULL;
	Int4 number_of_descriptions, number_of_alignments;
	SeqAlignPtr  seqalign, PNTR seqalign_array;
        SeqAnnotPtr seqannot;
	SeqEntryPtr PNTR sepp;
	TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
	Uint1 align_type, align_view;
	Uint4 align_options, print_options;
	ValNodePtr  mask_loc, mask_loc_start, vnp, other_returns, error_returns;

	CharPtr blast_program, blast_database, blast_inputfile, blast_outputfile;
	FILE *infp, *outfp;
	Int4 index, num_bsps, total_length;
	Int2 ctr = 1;
	Char prefix[2];
	SeqLocPtr PNTR mask_slpp;
	Boolean done;

        if (! GetArgs ("mblastall", NUMARG, myargs))
	   return (1);

	UseLocalAsnloadDataAndErrMsg ();

	if (! SeqEntryLoad())
		return 1;

	ErrSetMessageLevel(SEV_WARNING);

	blast_program = myargs [0].strvalue;
        blast_database = myargs [1].strvalue;
        blast_inputfile = myargs [2].strvalue;
        blast_outputfile = myargs [5].strvalue;
	if (myargs[24].intvalue)
		html = TRUE;

	if ((infp = FileOpen(blast_inputfile, "r")) == NULL) {
	   ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile);
	   return (1);
	}

	outfp = NULL;
	if (blast_outputfile != NULL) {
	   if ((outfp = FileOpen(blast_outputfile, "w")) == NULL) {
	      ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", blast_outputfile);
	      return (1);
	   }
	}

	align_view = (Int1) myargs[4].intvalue;

	align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);
 
        believe_query = FALSE;
        if (myargs[16].intvalue != 0)
                believe_query = TRUE;

	if (believe_query == FALSE && myargs[15].strvalue) 
	   ErrPostEx(SEV_FATAL, 0, 0, "-J option must be TRUE to produce a SeqAlign file");

	options = BLASTOptionNew(blast_program, TRUE);
	if (options == NULL)
		return 3;

	options->do_sum_stats = FALSE;
	options->is_neighboring = FALSE;
        options->expect_value  = (Nlm_FloatHi) myargs [3].floatvalue;
	number_of_descriptions = myargs[11].intvalue;	
	number_of_alignments = myargs[12].intvalue;	
	options->hitlist_size = MAX(number_of_descriptions, number_of_alignments);

	if (myargs[7].intvalue != 0)
		options->gap_x_dropoff = myargs[7].intvalue;
	if (StringICmp(myargs[6].strvalue, "T") == 0)
	   options->filter_string = StringSave("D");
	else
	   options->filter_string = StringSave(myargs[6].strvalue);
	
	show_gi = (Boolean) myargs[8].intvalue;
	options->penalty = myargs[9].intvalue;
	options->reward = myargs[10].intvalue;

	options->genetic_code = 1;
	options->db_genetic_code = 1; /* Default; it's not needed here anyway */
	options->number_of_cpus = myargs[14].intvalue;
	if (myargs[18].intvalue != 0)
		options->wordsize = myargs[18].intvalue;

	if (myargs[19].floatvalue != 0)
		options->db_length = (Int8) myargs[19].floatvalue;

        options->hsp_range_max  = myargs[20].intvalue;
        if (options->hsp_range_max != 0)
                options->perform_culling = TRUE;
        options->block_width  = myargs[21].intvalue;
        if (myargs[22].floatvalue)
                 options->searchsp_eff = (Nlm_FloatHi) myargs[22].floatvalue;

	options->strand_option = myargs[23].intvalue;


        print_options = 0;
        align_options = 0;
        align_options += TXALIGN_COMPRESS;
        align_options += TXALIGN_END_NUM;
        if (show_gi) {
	   align_options += TXALIGN_SHOW_GI;
	   print_options += TXALIGN_SHOW_GI;
        }
			
        if (align_view) {
	   align_options += TXALIGN_MASTER;
	   if (align_view == 1 || align_view == 3)
	      align_options += TXALIGN_MISMATCH;
	   if (align_view == 3 || align_view == 4 || align_view == 6)
	      align_options += TXALIGN_FLAT_INS;
	   if (align_view == 5 || align_view == 6)
	      align_options += TXALIGN_BLUNT_END;
        } else {
	   align_options += TXALIGN_MATRIX_VAL;
	   align_options += TXALIGN_SHOW_QS;
	}

	if (html) {
	   align_options += TXALIGN_HTML;
	   print_options += TXALIGN_HTML;
	}

	if (myargs[25].strvalue)
	   options->gifile = StringSave(myargs[25].strvalue);
   
	options->is_megablast_search = TRUE;
	if (myargs[13].intvalue == MBLAST_ENDPOINTS)
	   options->one_line_results = TRUE;
	else
	   options->one_line_results = FALSE;

	query_bsp_array = (BioseqPtr PNTR) MemNew(MAX_NUM_QUERIES*sizeof(BioseqPtr));
	sepp = (SeqEntryPtr PNTR) MemNew(MAX_NUM_QUERIES*sizeof(SeqEntryPtr));
	mask_slpp = (SeqLocPtr PNTR) MemNew(MAX_NUM_QUERIES*sizeof(SeqLocPtr));

	StrCpy(prefix, "Q");

	global_fp = outfp;
	if (myargs[13].intvalue==MBLAST_ALIGNMENTS) {
	   if (html)
	      fprintf(outfp, "<PRE>\n");
	   init_buff_ex(90);
	   BlastPrintVersionInfo(blast_program, html, outfp);
	   fprintf(outfp, "\n");
	   BlastPrintReference(html, 90, outfp);
	   fprintf(outfp, "\n");
	   PrintDbInformation(blast_database, !db_is_na, 70, outfp, html);
	   free_buff();
	
#ifdef OS_UNIX
	fprintf(global_fp, "%s", "Searching");
#endif
	}
	
	done = FALSE;
	while (!done) {
	   num_bsps = 0;
	   total_length = 0;
	   done = TRUE;
	   while ((sepp[num_bsps]=FastaToSeqEntryForDb(infp, query_is_na, NULL,
						       believe_query, prefix, &ctr, 
						       &mask_slpp[num_bsps])) != NULL) {
	      query_bsp = NULL;
	      if (query_is_na) 
		 SeqEntryExplore(sepp[num_bsps], &query_bsp, FindNuc);
	      else
		 SeqEntryExplore(sepp[num_bsps], &query_bsp, FindProt);
	      
	      if (query_bsp == NULL) {
		 ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
		 return 2;
	      }
	      
	      source = BioSourceNew();
	      source->org = OrgRefNew();
	      source->org->orgname = OrgNameNew();
	      source->org->orgname->gcode = options->genetic_code;
	      ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);
	      
	      query_bsp_array[num_bsps++] = query_bsp;
	      
	      total_length += query_bsp->length;
	      if (total_length > myargs[17].intvalue) {
		 done = FALSE;
		 break;
	      }
	   }
	   
	   other_returns = NULL;
	   error_returns = NULL;
	   
	   if (myargs[13].intvalue==MBLAST_ENDPOINTS) 
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
						     blast_database, options,
						     &other_returns, &error_returns,
						     dummy_callback, NULL, NULL, 0, 
						     mask_slpp, BlastSearchHandleResults);
	   else if (myargs[13].intvalue==MBLAST_SEGMENTS) 
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
						     blast_database, options,
						     &other_returns, &error_returns,
						     dummy_callback, NULL, NULL, 0,
						     mask_slpp, MegaBlastPrintSegments);
	   else /* if (myargs[13].intvalue==MBLAST_ALIGNMENTS) */
	      seqalign_array = BioseqMegaBlastEngine(query_bsp_array, blast_program,
						     blast_database, options,
						     &other_returns, &error_returns,
						     tick_callback, NULL, NULL, 0,
						     mask_slpp, NULL);
	   
	   if (myargs[13].intvalue==MBLAST_ALIGNMENTS) {
	      BlastErrorPrint(error_returns);
	      dbinfo = NULL;
	      ka_params = NULL;
	      ka_params_gap = NULL;
	      params_buffer = NULL;
	      mask_loc = NULL;
	      matrix = NULL;
	      for (vnp=other_returns; vnp; vnp = vnp->next) {
		 switch (vnp->choice) {
		 case TXDBINFO:
		    dbinfo = vnp->data.ptrvalue;
		    break;
		 case TXKABLK_NOGAP:
		    ka_params = vnp->data.ptrvalue;
		    break;
		 case TXKABLK_GAP:
		    ka_params_gap = vnp->data.ptrvalue;
		    break;
		 case TXPARAMETERS:
		    params_buffer = vnp->data.ptrvalue;
		    break;
		 case TXMATRIX:
		    matrix = vnp->data.ptrvalue;
		    break;
		 case SEQLOC_MASKING_NOTSET:
		 case SEQLOC_MASKING_PLUS1:
		 case SEQLOC_MASKING_PLUS2:
		 case SEQLOC_MASKING_PLUS3:
		 case SEQLOC_MASKING_MINUS1:
		 case SEQLOC_MASKING_MINUS2:
		 case SEQLOC_MASKING_MINUS3:
		    /*ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);*/
		    break;
		 default:
		    break;
		 }
	      }	
	      
	      
#ifdef OS_UNIX
	      fflush(global_fp);
#endif
	      
#ifdef OS_UNIX
	      fprintf(global_fp, "%s\n", "done");
#endif
	      
	      aip = NULL;
	      if (myargs[15].strvalue != NULL) {
		 if ((aip = AsnIoOpen (myargs[15].strvalue,"w")) == NULL) {
		    ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", "blastngp.sat");
		    return 1;
		 }
	      }
	      
	      ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);
	      
	      for (index=0; index<num_bsps; index++) {
		 if (seqalign_array[index]==NULL) continue;
		 init_buff_ex(70);
		 AcknowledgeBlastQuery(query_bsp_array[index], 70, outfp, believe_query, html);
		 free_buff();
		 
		 seqalign = seqalign_array[index];
		 if (seqalign) {
		    seqannot = SeqAnnotNew();
		    seqannot->type = 2;
		    AddAlignInfoToSeqAnnot(seqannot, align_type);
		    seqannot->data = seqalign;
		    if (aip) {
		       SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, aip, NULL);
		       AsnIoReset(aip);
		       aip = AsnIoClose(aip);
		    }
		    if (outfp) { /* Uncacheing causes problems with ordinal nos. vs. gi's. */
		       prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_descriptions, NULL);
		       ObjMgrSetHold();
		       init_buff_ex(85);
		       PrintDefLinesFromSeqAlign(prune->sap, 80,
						 outfp, print_options, FIRST_PASS, NULL);
		       free_buff();
		       
		       prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_alignments, prune);
		       seqannot->data = prune->sap;
		       if (align_view != 0)
			  ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL,
						 NULL, align_options, NULL, 
						 mask_loc, NULL);
		       else
			  ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, align_options, NULL, mask_loc, FormatScoreFunc);
		       seqannot->data = seqalign;
		       prune = BlastPruneSapStructDestruct(prune);
		       ObjMgrClearHold();
		       ObjMgrFreeCache(0);
		    }
		    seqannot = SeqAnnotFree(seqannot);
		 } else 
		    fprintf(outfp, "\n\n ***** No hits found ******\n\n");
	      } /* End loop on seqaligns for different queries */
	      
	      matrix = BLAST_MatrixDestruct(matrix);
	      
	      init_buff_ex(85);
	      dbinfo_head = dbinfo;
	      while (dbinfo) {
		 PrintDbReport(dbinfo, 70, outfp);
		 dbinfo = dbinfo->next;
	      }
	      dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
	      
	      if (ka_params) {
		 PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
		 MemFree(ka_params);
	      }
	      
	      if (ka_params_gap) {
		 PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
		 MemFree(ka_params_gap);
	      }
	      
	      PrintTildeSepLines(params_buffer, 70, outfp);
	      MemFree(params_buffer);
	      free_buff();
	      mask_loc_start = mask_loc;
	      while (mask_loc) {
		 SeqLocSetFree(mask_loc->data.ptrvalue);
		 mask_loc = mask_loc->next;
	      }
	      ValNodeFree(mask_loc_start);
	      
	      ReadDBBioseqFetchDisable();
	      other_returns = ValNodeFree(other_returns);
	      MemFree(seqalign_array);
	   } /* End producing alignment output */
	   
	   
	   for (index=0; index<num_bsps; index++) {
	      sepp[index] = SeqEntryFree(sepp[index]);
	      query_bsp_array[index] = NULL;
	      mask_slpp[index] = NULL;
	   }	   
	} /* End of loop on complete searches */
	MemFree(query_bsp_array);
	MemFree(sepp);
	MemFree(mask_slpp);
	options = BLASTOptionDelete(options);
	FileClose(infp);
	
	return 0;
}
