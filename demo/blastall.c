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
**************************************************************************/
/* $Revision 1.0$ */ 
/* $Log: blastall.c,v $
/* Revision 6.28  1999/01/22 17:24:51  madden
/* added line breaks for alignment views
/*
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

static int LIBCALLBACK
star_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{

#ifdef OS_UNIX

	fprintf(global_fp, "%s", "*");
	fflush(global_fp);
#endif
	return 0;
}

/* find the last nucleotide bioseq in the bioseqset */
static void FindNuc(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
    BioseqPtr PNTR bp;
    BioseqPtr local_bsp;

    bp = (BioseqPtr PNTR) data;
    if (IS_Bioseq(sep))
    {
        local_bsp = (BioseqPtr) sep->data.ptrvalue;
        if (ISA_na(local_bsp->mol))
          *bp = local_bsp;
    }
}

/* find the last protein bioseq in the bioseqset */
static void FindProt(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
    BioseqPtr PNTR bp;
    BioseqPtr local_bsp;

    bp = (BioseqPtr PNTR) data;
    if (IS_Bioseq(sep))
    {
        local_bsp = (BioseqPtr) sep->data.ptrvalue;
        if (ISA_aa(local_bsp->mol))
          *bp = local_bsp;
    }
}

#define NUMARG 29

static Args myargs [NUMARG] = {
 { "Program Name",
        NULL, NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
  { "Database", 
	"nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
  { "Query File", 
	"stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  { "Expectation value (E)", 
	"10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
  { "alignment view options:\n0 = pairwise,\n1 = master-slave showing identities,\n2 = master-slave no identities,\n3 = flat master-slave, show identities,\n4 = flat master-slave, no identities,\n5 = master-slave no identities and blunt ends,\n6 = flat master-slave, no identities and blunt ends", 
        "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},
  { "BLAST report Output File", 
	"stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  { "Filter query sequence (DUST with blastn, SEG with others)",
        "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
  { "Cost to open a gap (zero invokes default behavior)",
	"0", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},
  { "Cost to extend a gap (zero invokes default behavior)",
	"0", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},
  { "X dropoff value for gapped alignment (in bits) (zero invokes default behavior)",
	"0", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
  { "Show GI's in deflines",
        "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},
  { "Penalty for a nucleotide mismatch (blastn only)",
	"-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},
  { "Reward for a nucleotide match (blastn only)",
	"1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},
  { "Number of one-line descriptions (V)",
        "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},
  { "Number of alignments to show (B)",
        "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
  { "Threshold for extending hits, default if zero",
        "0", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},
  { "Perfom gapped alignment (not available with tblastx)",
        "T", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
  { "Query Genetic code to use",
        "1", NULL, NULL, FALSE, 'Q', ARG_INT, 0.0, 0, NULL},
  { "DB Genetic code (for tblast[nx] only)",
        "1", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},
  { "Number of processors to use",
        "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},
  { "SeqAlign file", 
	NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},
  { "Believe the query defline",
        "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
  { "Matrix",
        "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
  { "Word size, default if zero", 
        "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},
  { "Effective length of the database (use zero for the real size)", 
        "0", NULL, NULL, FALSE, 'z', ARG_INT, 0.0, 0, NULL},
  { "Number of best hits from a region to keep",
        "100", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},
  { "Length of region used to judge hits",
        "20", NULL, NULL, FALSE, 'L', ARG_INT, 0.0, 0, NULL},
  { "Effective length of the search space (use zero for the real size)",
        "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},
  { "Query strands to search against database (for blast[nx], and tblastx).  3 is both, 1 is top, 2 is bottom",
        "3", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},
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
	Char buffer[256];
	CharPtr ret_buffer=NULL, params_buffer=NULL;
	Int4 number_of_descriptions, number_of_alignments;
	ObjectIdPtr obidp;
	SeqAlignPtr  seqalign;
        SeqAnnotPtr seqannot;
	SeqEntryPtr sep;
	SeqIdPtr seqid_list=NULL;
	TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
	Uint1 align_type, align_view;
	Uint4 align_options, print_options;
	ValNodePtr  mask_loc, vnp, other_returns, error_returns;

	CharPtr blast_program, blast_database, blast_inputfile, blast_outputfile;
	FILE *infp, *outfp;

        if (! GetArgs ("blastall", NUMARG, myargs))
        {
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

	if ((infp = FileOpen(blast_inputfile, "r")) == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile);
		return (1);
	}

	outfp = NULL;
	if (blast_outputfile != NULL)
	{
		if ((outfp = FileOpen(blast_outputfile, "w")) == NULL)
		{
			ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", blast_outputfile);
			return (1);
		}
	}

	align_view = (Int1) myargs[4].intvalue;

	align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);
	if (StringICmp("blastx", blast_program) == 0)
	{
		if (align_view != 0)
		{
	  	 	ErrPostEx(SEV_FATAL, 0, 0, "This option is not available with blastx");
			return 1;
		}
	}
	else if (StringICmp("tblastx", blast_program) == 0)
	{
		if (align_view != 0)
		{
	  	 	ErrPostEx(SEV_FATAL, 0, 0, "This option is not available with tblastx");
			return 1;
		}
	}
 
        believe_query = FALSE;
        if (myargs[21].intvalue != 0)
                believe_query = TRUE;

	if (believe_query == FALSE && myargs[20].strvalue)
	{
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
	if (StringICmp(myargs[6].strvalue, "T") == 0)
	{
		if (StringICmp("blastn", blast_program) == 0)
			options->filter_string = StringSave("D");
		else
			options->filter_string = StringSave("S");
	}
	else
	{
		options->filter_string = StringSave(myargs[6].strvalue);
	}

	show_gi = (Boolean) myargs[10].intvalue;
	if (StringICmp("blastn", blast_program) == 0)
	{
		options->penalty = myargs[11].intvalue;
		options->reward = myargs[12].intvalue;
	}
	else
	{
		if (myargs[15].intvalue != 0)
		{
			options->threshold_first = myargs[15].intvalue;
			options->threshold_second = myargs[15].intvalue;
		}
	}

	options->genetic_code = myargs[17].intvalue;
	options->db_genetic_code = myargs[18].intvalue;
	options->number_of_cpus = myargs[19].intvalue;
	if (myargs[23].intvalue != 0)
		options->wordsize = myargs[23].intvalue;
	if (myargs[24].intvalue != 0)
		options->db_length = myargs[24].intvalue;

        options->hsp_range_max  = myargs[25].intvalue;
        if (options->hsp_range_max != 0)
                options->perform_culling = TRUE;
        options->block_width  = myargs[26].intvalue;
        if (myargs[27].floatvalue)
                 options->searchsp_eff = (Nlm_FloatHi) myargs[27].floatvalue;

	options->strand_option = myargs[28].intvalue;


        print_options = 0;
        align_options = 0;
        align_options += TXALIGN_COMPRESS;
        align_options += TXALIGN_END_NUM;
	if (StringICmp("blastx", blast_program) == 0)
	{
		align_options += TXALIGN_BLASTX_SPECIAL;
	}
        if (show_gi) {
                align_options += TXALIGN_SHOW_GI;
                print_options += TXALIGN_SHOW_GI;
        }
	if (myargs[16].intvalue == 0)
                 print_options += TXALIGN_SHOW_NO_OF_SEGS;
			
        if (align_view)
        {
                        align_options += TXALIGN_MASTER;
                        if (align_view == 1 || align_view == 3)
                                align_options += TXALIGN_MISMATCH;
                        if (align_view == 3 || align_view == 4 || align_view == 6)
                                align_options += TXALIGN_FLAT_INS;
                        if (align_view == 5 || align_view == 6)
                                align_options += TXALIGN_BLUNT_END;
        }
        else
        {
                        align_options += TXALIGN_MATRIX_VAL;
                        align_options += TXALIGN_SHOW_QS;
	}

	while ((sep=FastaToSeqEntryEx(infp, query_is_na, NULL, believe_query)) != NULL) 
	{
		query_bsp = NULL;
		if (query_is_na)
		{
			SeqEntryExplore(sep, &query_bsp, FindNuc);
		}
		else
		{
			SeqEntryExplore(sep, &query_bsp, FindProt);
		}

		if (query_bsp == NULL)
		{
	  	 	ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
	   		return 2;
		}

		fake_bsp = BlastMakeFakeBioseq(query_bsp, NULL);
		
		source = BioSourceNew();
		source->org = OrgRefNew();
		source->org->orgname = OrgNameNew();
		source->org->orgname->gcode = options->genetic_code;
		ValNodeAddPointer(&(query_bsp->descr), Seq_descr_source, source);

		global_fp = outfp;

		init_buff_ex(90);
		BlastPrintVersionInfo(blast_program, FALSE, outfp);
		fprintf(outfp, "\n");
		BlastPrintReference(FALSE, 90, outfp);
		fprintf(outfp, "\n");
		AcknowledgeBlastQuery(query_bsp, 70, outfp, believe_query, FALSE);
                PrintDbInformation(blast_database, !db_is_na, 70, outfp, FALSE);
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
		for (vnp=other_returns; vnp; vnp = vnp->next)
		{
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
		if (myargs[20].strvalue != NULL)
		{
       			if ((aip = AsnIoOpen (myargs[20].strvalue,"w")) == NULL)
        		{
               		 ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", "blastngp.sat");
               		 return 1;
        		}
		}

       		ReadDBBioseqFetchEnable ("blastall", blast_database, db_is_na, TRUE);

        	if (seqalign)
        	{
               	 	seqannot = SeqAnnotNew();
               	 	seqannot->type = 2;
			AddAlignInfoToSeqAnnot(seqannot, align_type);
                	seqannot->data = seqalign;
			if (aip)
			{
		 		SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, aip, NULL);
		 		AsnIoReset(aip);
		 		aip = AsnIoClose(aip);
			}
			if (outfp)
			{	/* Uncacheing causes problems with ordinal nos. vs. gi's. */
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
		}
		else
		{
			fprintf(outfp, "\n\n ***** No hits found ******\n\n");
		}

		matrix = BLAST_MatrixDestruct(matrix);

		init_buff_ex(85);
		dbinfo_head = dbinfo;
		while (dbinfo)
		{
                	PrintDbReport(dbinfo, 70, outfp);
			dbinfo = dbinfo->next;
		}
		dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);

		if (ka_params)
		{
                	PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
			MemFree(ka_params);
		}

		if (ka_params_gap)
		{
                	PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, TRUE);
			MemFree(ka_params_gap);
		}

                PrintTildeSepLines(params_buffer, 70, outfp);
                MemFree(params_buffer);
                free_buff();

		fake_bsp = BlastDeleteFakeBioseq(fake_bsp);

		ReadDBBioseqFetchDisable();
		other_returns = ValNodeFree(other_returns);
		sep = SeqEntryFree(sep);
	}
	options = BLASTOptionDelete(options);
	FileClose(infp);

	return 0;
}
	

