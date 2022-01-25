/* $Id: blastall.c,v 6.45 2000/05/26 19:28:44 shavirin Exp $
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
 * 
 * $Log: blastall.c,v $
 * Revision 6.45  2000/05/26 19:28:44  shavirin
 * Added adjustment of dropoff_1st_pass if dropoff_1st_pass > dropoff_2nd_pass
 *
 * Revision 6.44  2000/05/26 18:48:23  shavirin
 * Added two new parameters; '-y' and '-Z'
 *
 * Revision 6.43  2000/05/09 15:57:26  shavirin
 * Added call to the function ReadDBBioseqSetDbGeneticCode().
 *
 * Revision 6.42  2000/04/25 20:50:45  dondosha
 * Removed unavailable option to use greedy algorithm
 *
 * Revision 6.41  2000/04/13 13:34:19  shavirin
 * Added call to ObjMgrFreeCache() back after fixes in API.
 *
 * Revision 6.40  2000/04/04 18:29:13  shavirin
 * Added some missing HTML tags.
 *
 * Revision 6.39  2000/03/31 19:13:33  dondosha
 * Changed some names related to MegaBlast
 *
 * Revision 6.38  2000/03/24 21:49:30  madden
 * Comment out ObjMgrFreeCache
 *
 * Revision 6.37  2000/03/02 21:06:09  shavirin
 * Added -U option, that allows to consider low characters in FASTA files
 * as filtered regions (for blastn, blastp and tblastn).
 *
 * Revision 6.36  2000/02/01 20:05:31  dondosha
 * Added option -B: use greedy basic alignment search if set to T
 *
 * Revision 6.35  2000/01/28 16:46:54  madden
 * Added function BlastGetMaskingLoc
 *
 * Revision 6.34  1999/12/17 20:48:53  egorov
 * Fix 'gcc -Wall' warnings and remove old stuff.
 *
 * Revision 6.33  1999/10/12 19:35:26  madden
 * Deallocate Mask information
 *
 * Revision 6.32  1999/08/26 14:58:06  madden
 * Use float for db length
 *
 * Revision 6.31  1999/05/26 13:12:56  madden
 * Initialized matrix to NULL
 *
 * Revision 6.30  1999/03/31 16:58:04  madden
 * Removed static FindProt and FindNuc
 *
 * Revision 6.29  1999/02/10 21:12:26  madden
 * Added HTML and GI list option, fixed filtering
 *
 * Revision 6.28  1999/01/22 17:24:51  madden
 * added line breaks for alignment views
 *
 * Revision 6.27  1998/12/31 18:18:27  madden
 * Added strand option
 *
 * Revision 6.26  1998/12/29 20:03:14  kans
 * calls UseLocalAsnloadDataAndErrMsg at startup
 *
 * Revision 6.25  1998/11/19 14:04:34  madden
 * Changed message level to SEV_WARNING
 *
 * Revision 6.24  1998/11/16 16:29:19  madden
 * Added ErrSetMessageLevel(SEV_INFO)
 *
 * Revision 6.23  1998/07/17 15:41:36  madden
 * Added effective search space flag
 *
 * Revision 6.22  1998/06/29 13:02:01  madden
 * Deallocate matrix
 *
 * Revision 6.21  1998/06/10 13:33:14  madden
 * Change -K from zero to 100
 *
 * Revision 6.20  1998/06/05 21:48:42  madden
 * Added -K and -L options
 *
 * Revision 6.19  1998/05/18 18:01:04  madden
 * Changed args to allow filter options to be changed
 *
 * Revision 6.18  1998/05/01 18:31:02  egorov
 * Add new parametes to BLASTOptionSetGapParam()
 *
 * Revision 6.17  1998/04/30 14:32:32  madden
 * init_buff_ex arg changed to 90 for reference
 *
 * Revision 6.16  1998/04/29 14:29:30  madden
 * Made reference line longer
 *
 * Revision 6.15  1998/04/01 22:49:12  madden
 * Print No hits found message
 *
 * Revision 6.14  1998/02/25 20:50:48  madden
 * Added arg for db length
 *
 * Revision 6.13  1998/02/24 22:48:34  madden
 * Removed options for culling
 *
 * Revision 6.12  1998/01/31 21:35:17  madden
 * zeroed out values between searches
 *
 * Revision 6.11  1997/12/31 17:48:52  madden
 * Added wordsize option
 *
 * Revision 6.10  1997/12/23 21:09:47  madden
 * Added -K and -L for range-dependent blast
 *
 * Revision 6.9  1997/11/19 14:26:43  madden
 * Removed extra break statement
 *
 * Revision 6.8  1997/11/18 22:24:22  madden
 * Added call to BLASTOptionSetGapParams
 *
 * Revision 6.7  1997/10/27 22:26:52  madden
 * Added call to ObjMgrFreeCache(0)
 *
 * Revision 6.6  1997/10/23 20:26:12  madden
 * Use of init_buff_ex rather than init_buff
 *
 * Revision 6.5  1997/10/22 21:56:04  madden
 * Added matrix option
 *
 * Revision 6.3  1997/10/07 21:33:38  madden
 * Added BLUNT option
 *
 * Revision 6.2  1997/09/23 22:13:19  madden
 * enabled descriptions and alignment options
 *
 * Revision 6.1  1997/09/16 16:34:32  madden
 * Dbinfo printing changed for multiple db searches
 *
 * Revision 6.0  1997/08/25 18:19:14  madden
 * Revision changed to 6.0
 *
 * Revision 1.16  1997/07/29 19:33:02  madden
 * Added TXALIGN_SHOW_QS flag
 *
 * Revision 1.15  1997/07/28 17:01:23  madden
 * Added include for simutil.h
 *
 * Revision 1.14  1997/07/28 14:31:09  madden
 * Changes for masking alignments.
 *
 * Revision 1.13  1997/07/22 19:06:35  madden
 * Option changes, Printing of verison info
 *
 * Revision 1.12  1997/07/18 20:09:22  madden
 * Conversion from blast2 output to new output
 *
 * Revision 1.3  1997/02/24  22:08:38  madden
 * Added reward and penalty for match and mismatch.
 *
 * Revision 1.2  1997/02/23  16:48:52  madden
 * Call to AcknowledgeBlastQuery added.
 *
 * Revision 1.1  1997/02/19  21:44:28  madden
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


#define DEFLINE_BUF 255


/* Used by the callback function. */
FILE *global_fp=NULL;
/*
	Callback to print out ticks, in UNIX only due to file systems
	portability issues.
*/

static int LIBCALLBACK
tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{

#ifdef OS_UNIX

	fprintf(global_fp, "%s", ".");
	fflush(global_fp);
#endif
	return 0;
}

static Int2
BlastGetMaskingLoc(FILE *infp, FILE *outfp, CharPtr instructions)
{
	BioseqPtr bsp;
	Char buffer[50];
	SeqEntryPtr sep;
	SeqLocPtr slp, slp_start, tmp_slp;

	if (infp == NULL || outfp == NULL || instructions == NULL)
		return 1;

	while ((sep=FastaToSeqEntryEx(infp, TRUE, NULL, TRUE)) != NULL) 
	{
		bsp = NULL;
		SeqEntryExplore(sep, &bsp, FindNuc);

		if (bsp == NULL)
		{
	  	 	ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
	   		return 2;
		}
		SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, 50);
		fprintf(outfp, ">%s\n", buffer);
		slp_start = slp = BlastBioseqFilter(bsp, instructions);
        	while (slp)
        	{
               		tmp_slp=NULL;
               		while((tmp_slp = SeqLocFindNext(slp, tmp_slp))!=NULL)
               	 	{
				fprintf(outfp, "%ld %ld\n", (long) (1+SeqLocStart(tmp_slp)), (long) (1+SeqLocStop(tmp_slp)));
                 	}
                	slp = slp->next;
        	}

/* used for debugging. */
#if 0
{{
	BioseqPtr bsp_tmp;
	ByteStorePtr byte_sp;
	Int4 index;
	SeqLocPtr tmp_slp_1, tmp_filter_slp;
	SeqPortPtr spp;
	Uint1Ptr tmp_query_seq, tmp_query_seq_start;
	Uint1 residue;
	FILE *tmp_fp;

		spp = SeqPortNew(bsp, 0, -1, 0, Seq_code_iupacna);
                SeqPortSet_do_virtual(spp, TRUE);
		tmp_query_seq_start = (Uint1Ptr) MemNew(((BioseqGetLen(bsp))+2)*sizeof(Uint1));
		tmp_query_seq_start[0] = NULLB;
		tmp_query_seq = tmp_query_seq_start+1;
		index=0;
		while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
		{

			if (IS_residue(residue))
			{
				tmp_query_seq[index] = residue;
				index++;
			}
		}
		BlastMaskTheResidues(tmp_query_seq, BioseqGetLen(bsp), 78, slp_start, FALSE, 0);
		bsp_tmp = BioseqNew();
		bsp_tmp->length = BioseqGetLen(bsp);
		byte_sp = BSNew(1);
		BSWrite(byte_sp, tmp_query_seq, bsp->length);
		bsp_tmp->seq_data = byte_sp;
		bsp_tmp->repr = Seq_repr_raw;
		bsp_tmp->seq_data_type = Seq_code_iupacna;
		bsp_tmp->mol = 1;

		bsp_tmp->id = bsp->id;
		bsp_tmp->descr = bsp->descr;

		tmp_fp = FileOpen("masked.fsa", "w");
		BioseqRawToFastaExtra(bsp_tmp, tmp_fp, 50);

		bsp_tmp->id = NULL;
		bsp_tmp->descr = NULL;

		spp = SeqPortFree(spp);
		bsp_tmp = BioseqFree(bsp_tmp);
		tmp_query_seq_start = MemFree(tmp_query_seq_start);
		FileClose(tmp_fp);

		tmp_filter_slp = slp_start;
		tmp_fp = FileOpen("locations.msk", "w");
        	while (tmp_filter_slp)
        	{
               	 tmp_slp_1=NULL;
               	 while((tmp_slp_1 = SeqLocFindNext(tmp_filter_slp, tmp_slp_1))!=NULL)
               	 {
			fprintf(tmp_fp, "%ld %ld\n", (long) (1+SeqLocStart(tmp_slp_1)), (long) (1+SeqLocStop(tmp_slp_1)));

                 }
                	tmp_filter_slp = tmp_filter_slp->next;
        	}


		FileClose(tmp_fp);
}}
#endif
		slp_start = SeqLocSetFree(slp_start);
		sep = SeqEntryFree(sep);
	}

	return 0;
}

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs [] = {
    { "Program Name",           /* 0 */
      NULL, NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
    { "Database",               /* 1 */
      "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
    { "Query File",             /* 2 */
      "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
    { "Expectation value (E)",  /* 3 */
      "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
    { "alignment view options:\n0 = pairwise,\n1 = master-slave showing identities,\n2 = master-slave no identities,\n3 = flat master-slave, show identities,\n4 = flat master-slave, no identities,\n5 = master-slave no identities and blunt ends,\n6 = flat master-slave, no identities and blunt ends", /* 4 */
      "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},
    { "BLAST report Output File", /* 5 */
      "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Filter query sequence (DUST with blastn, SEG with others)", /* 6 */
      "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
    { "Cost to open a gap (zero invokes default behavior)", /* 7 */
      "0", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},
    { "Cost to extend a gap (zero invokes default behavior)", /* 8 */
      "0", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},
    { "X dropoff value for gapped alignment (in bits) (zero invokes default behavior)", /* 9 */
      "0", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
    { "Show GI's in deflines",  /* 10 */
      "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Penalty for a nucleotide mismatch (blastn only)", /* 11 */
      "-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},
    { "Reward for a nucleotide match (blastn only)", /* 12 */
      "1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},
    { "Number of one-line descriptions (V)", /* 13 */
      "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},
    { "Number of alignments to show (B)", /* 14 */
      "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
    { "Threshold for extending hits, default if zero", /* 15 */
      "0", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},
    { "Perfom gapped alignment (not available with tblastx)", /* 16 */
        "T", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Query Genetic code to use", /* 17 */
      "1", NULL, NULL, FALSE, 'Q', ARG_INT, 0.0, 0, NULL},
    { "DB Genetic code (for tblast[nx] only)", /* 18 */
      "1", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},
    { "Number of processors to use", /* 19 */
      "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},
    { "SeqAlign file",          /* 20 */
      NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Believe the query defline", /* 21 */
      "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Matrix",                 /* 22 */
      "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
    { "Word size, default if zero", /* 23 */
      "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},
    { "Effective length of the database (use zero for the real size)", /* 24 */
      "0", NULL, NULL, FALSE, 'z', ARG_FLOAT, 0.0, 0, NULL},
    { "Number of best hits from a region to keep", /* 25 */
      "100", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},
    { "Length of region used to judge hits", /* 26 */
      "20", NULL, NULL, FALSE, 'L', ARG_INT, 0.0, 0, NULL},
    { "Effective length of the search space (use zero for the real size)", /* 27 */
      "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},
    { "Query strands to search against database (for blast[nx], and tblastx).  3 is both, 1 is top, 2 is bottom", /* 28 */
      "3", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},
    { "Produce HTML output",    /* 29 */
      "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Restrict search of database to list of GI's", /* 30 */
      NULL, NULL, NULL, TRUE, 'l', ARG_STRING, 0.0, 0, NULL},
    {"Use lower case filtering of FASTA sequence", /* 31 */
     "F", NULL,NULL,TRUE,'U',ARG_BOOLEAN, 0.0,0,NULL},
    { "Dropoff (X) for blast extensions in bits (default if zero)", /* 32 */
      "0.0", NULL, NULL, FALSE, 'y', ARG_FLOAT, 0.0, 0, NULL},
    { "X dropoff value for final gapped alignment (in bits)", /* 33 */
      "0", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},
};

Int2 Main (void)
 
{
    AsnIoPtr aip;
    BioseqPtr fake_bsp, query_bsp;
    BioSourcePtr source;
    BLAST_MatrixPtr matrix;
    BLAST_OptionsBlkPtr options;
    BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
    BlastPruneSapStructPtr prune;
    Boolean db_is_na, query_is_na, show_gi, believe_query=FALSE;
    Boolean html = FALSE;
    CharPtr params_buffer=NULL;
    Int4 number_of_descriptions, number_of_alignments;
    SeqAlignPtr  seqalign;
    SeqAnnotPtr seqannot;
    SeqEntryPtr sep;
    TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
    Uint1 align_type, align_view;
    Uint4 align_options, print_options;
    ValNodePtr  mask_loc, mask_loc_start, vnp, other_returns, error_returns;
    
    CharPtr blast_program, blast_database, blast_inputfile, blast_outputfile;
    FILE *infp, *outfp;
    
    if (! GetArgs ("blastall", NUMARG, myargs)) {
        return (1);
    }
    
    UseLocalAsnloadDataAndErrMsg ();
    
    if (! SeqEntryLoad())
        return 1;
    
    ErrSetMessageLevel(SEV_WARNING);
    
    blast_program = myargs [0].strvalue;
    blast_database = myargs [1].strvalue;
    blast_inputfile = myargs [2].strvalue;
    blast_outputfile = myargs [5].strvalue;
    if (myargs[29].intvalue)
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
    
    if (StringCmp("filter", blast_program) == 0) {
        BlastGetMaskingLoc(infp, outfp, myargs[6].strvalue);
        FileClose(outfp);
        FileClose(infp);	
        return 0;
    }
    align_view = (Int1) myargs[4].intvalue;
    
    align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);
    if (StringICmp("blastx", blast_program) == 0) {
        if (align_view != 0) {
            ErrPostEx(SEV_FATAL, 0, 0, "This option is not available with blastx");
            return 1;
        }
    } else if (StringICmp("tblastx", blast_program) == 0) {
        if (align_view != 0) {
            ErrPostEx(SEV_FATAL, 0, 0, "This option is not available with tblastx");
            return 1;
        }
    }
    
    believe_query = FALSE;
    if (myargs[21].intvalue != 0)
        believe_query = TRUE;
    
    if (believe_query == FALSE && myargs[20].strvalue) {
        ErrPostEx(SEV_FATAL, 0, 0, "-J option must be TRUE to produce a SeqAlign file");
    }
    
    options = BLASTOptionNew(blast_program, (Boolean) myargs [16].intvalue);
    if (options == NULL)
        return 3;
    
    BLASTOptionSetGapParams(options, myargs[22].strvalue, 0, 0); 
    options->expect_value  = (Nlm_FloatHi) myargs [3].floatvalue;
    number_of_descriptions = myargs[13].intvalue;	
    number_of_alignments = myargs[14].intvalue;	
    options->hitlist_size = MAX(number_of_descriptions, number_of_alignments);
    if (myargs[7].intvalue != 0)
        options->gap_open = myargs[7].intvalue;
    if (myargs[8].intvalue != 0)
        options->gap_extend = myargs[8].intvalue;
    if (myargs[9].intvalue != 0)
        options->gap_x_dropoff = myargs[9].intvalue;

    if(myargs[33].intvalue != 0) 
        options->gap_x_dropoff_final = myargs[33].intvalue;

    if (StringICmp(myargs[6].strvalue, "T") == 0) {
        if (StringICmp("blastn", blast_program) == 0)
            options->filter_string = StringSave("D");
        else
            options->filter_string = StringSave("S");
    } else {
        options->filter_string = StringSave(myargs[6].strvalue);
    }
    
    show_gi = (Boolean) myargs[10].intvalue;
    if (StringICmp("blastn", blast_program) == 0) {
        options->penalty = myargs[11].intvalue;
        options->reward = myargs[12].intvalue;
    } else {
        if (myargs[15].intvalue != 0) {
            options->threshold_first = myargs[15].intvalue;
            options->threshold_second = myargs[15].intvalue;
        }
    }
    
    options->genetic_code = myargs[17].intvalue;
    options->db_genetic_code = myargs[18].intvalue;
    options->number_of_cpus = myargs[19].intvalue;
    if (myargs[23].intvalue != 0)
        options->wordsize = myargs[23].intvalue;
    if (myargs[24].floatvalue != 0)
        options->db_length = (Int8) myargs[24].floatvalue;
    
    options->hsp_range_max  = myargs[25].intvalue;
    if (options->hsp_range_max != 0)
        options->perform_culling = TRUE;
    options->block_width  = myargs[26].intvalue;
    if (myargs[27].floatvalue)
        options->searchsp_eff = (Nlm_FloatHi) myargs[27].floatvalue;
    
    options->strand_option = myargs[28].intvalue;

    if(myargs [32].floatvalue != 0.0) {
        options->dropoff_2nd_pass  = myargs [32].floatvalue;
        if(options->dropoff_1st_pass > options->dropoff_2nd_pass)
            options->dropoff_1st_pass = options->dropoff_2nd_pass;
    }
    
    print_options = 0;
    align_options = 0;
    align_options += TXALIGN_COMPRESS;
    align_options += TXALIGN_END_NUM;
    if (StringICmp("blastx", blast_program) == 0) {
        align_options += TXALIGN_BLASTX_SPECIAL;
    }
    if (show_gi) {
        align_options += TXALIGN_SHOW_GI;
        print_options += TXALIGN_SHOW_GI;
    }
    if (myargs[16].intvalue == 0)
        print_options += TXALIGN_SHOW_NO_OF_SEGS;
    
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
    
    if (myargs[30].strvalue) {
        options->gifile = StringSave(myargs[30].strvalue);
    }
    
    options->is_megablast_search = FALSE;
    
    while (TRUE) {
        if(myargs[31].intvalue) {
            sep = FastaToSeqEntryForDb (infp, query_is_na, NULL, believe_query, NULL, NULL, &options->query_lcase_mask);
       
        } else {
            sep = FastaToSeqEntryEx(infp, query_is_na, NULL, believe_query);
        }
        
        if(sep == NULL)
            break;
        
        query_bsp = NULL;
        if (query_is_na) {
            SeqEntryExplore(sep, &query_bsp, FindNuc);
        } else {
            SeqEntryExplore(sep, &query_bsp, FindProt);
        }
        
        if (query_bsp == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
            return 2;
        }
        
        fake_bsp = BlastMakeFakeBioseq(query_bsp, NULL);
        
        /* If fake_bsp created mask should be updated to use it's id */
        BLASTUpdateSeqIdInSeqInt(options->query_lcase_mask, fake_bsp->id);
        
        source = BioSourceNew();
        source->org = OrgRefNew();
        source->org->orgname = OrgNameNew();
        source->org->orgname->gcode = options->genetic_code;
        ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);
        
        global_fp = outfp;
        
        if (html) {
            fprintf(outfp, "<HTML>\n<TITLE>BLAST Search Results</TITLE>\n");
            fprintf(outfp, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
                    "VLINK=\"#660099\" ALINK=\"#660099\">\n");
            fprintf(outfp, "<PRE>\n");
        }

        init_buff_ex(90);
        BlastPrintVersionInfo(blast_program, html, outfp);
        fprintf(outfp, "\n");
        BlastPrintReference(html, 90, outfp);
        fprintf(outfp, "\n");
        AcknowledgeBlastQuery(query_bsp, 70, outfp, believe_query, html);
        PrintDbInformation(blast_database, !db_is_na, 70, outfp, html);
        free_buff();
        
#ifdef OS_UNIX
        fprintf(global_fp, "%s", "Searching");
#endif
        other_returns = NULL;
        error_returns = NULL;
        
        seqalign = BioseqBlastEngine(believe_query ? query_bsp : fake_bsp, blast_program, blast_database, options, &other_returns, &error_returns, tick_callback);
        
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
                ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
                break;
            default:
                break;
            }
        }	
        
        
#ifdef OS_UNIX
        fflush(global_fp);
#endif
        
#ifdef OS_UNIX
        fprintf(global_fp, "%s", "done");
#endif
        
        aip = NULL;
        if (myargs[20].strvalue != NULL) {
            if ((aip = AsnIoOpen (myargs[20].strvalue,"w")) == NULL) {
                ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", "blastngp.sat");
                return 1;
            }
        }
        
        ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);
        ReadDBBioseqSetDbGeneticCode(options->db_genetic_code);

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
                PrintDefLinesFromSeqAlign(prune->sap, 80, outfp, print_options, FIRST_PASS, NULL);
                free_buff();
                
                prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_alignments, prune);
                seqannot->data = prune->sap;
                if (align_view != 0)
                    ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, align_options, NULL, mask_loc, NULL);
                else
                    ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, align_options, NULL, mask_loc, FormatScoreFunc);
                seqannot->data = seqalign;
                prune = BlastPruneSapStructDestruct(prune);
                ObjMgrClearHold();

                ObjMgrFreeCache(0);

            }
            seqannot = SeqAnnotFree(seqannot);
        } else {
            fprintf(outfp, "\n\n ***** No hits found ******\n\n");
        }
        
        matrix = BLAST_MatrixDestruct(matrix);

        if(html) {
            fprintf(outfp, "<PRE>\n");
        }
        
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
        
        fake_bsp = BlastDeleteFakeBioseq(fake_bsp);
        
        ReadDBBioseqFetchDisable();
        other_returns = ValNodeFree(other_returns);
        sep = SeqEntryFree(sep);
    }

    if (html) {
        fprintf(outfp, "</PRE>\n</BODY>\n</HTML>\n");
    }

    options = BLASTOptionDelete(options);
    FileClose(infp);
    
    return 0;
}
	

