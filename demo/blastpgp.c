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
/* $Log: blastpgp.c,v $
/* Revision 6.48  1999/08/27 18:17:42  shavirin
/* Added new parameter to command line - decline_align.
/*
/* Revision 6.47  1999/08/26 14:58:06  madden
/* Use float for db length
/*
/* Revision 6.46  1999/08/04 13:26:49  madden
/* Added -B option
/*
/* Revision 6.45  1999/03/31 16:58:04  madden
/* Removed static FindProt and FindNuc
/*
/* Revision 6.44  1999/03/24 18:16:33  madden
/* zero-out groupOrder and featureOrder
/*
/* Revision 6.43  1999/03/21 19:43:23  madden
/* Added -Q option to store ASCII version of last position-specific matrix in a file
/*
/* Revision 6.42  1999/03/04 14:17:20  egorov
/* Added new parameter to BlastMaskTheResidues() function for correct masking
/* when query is seqloc.  The paramter is not used in this file and is 0.
/*
/* Revision 6.41  1999/02/18 21:13:11  madden
/* Check for non-pattern search before call to BlastPruneSapStructDestruct
/*
/* Revision 6.40  1999/02/10 21:12:27  madden
/* Added HTML and GI list option, fixed filtering
/*
/* Revision 6.39  1999/01/22 17:24:51  madden
/* added line breaks for alignment views
/*
/* Revision 6.38  1998/12/29 20:03:14  kans
/* calls UseLocalAsnloadDataAndErrMsg at startup
/*
 * Revision 6.37  1998/12/17 21:54:39  madden
 * Added call to BlastPruneHitsFromSeqAlign for other than first round
 *
 * Revision 6.36  1998/12/16 14:13:36  madden
 * Fix for more than one pattern in query
 *
 * Revision 6.35  1998/11/19 14:04:34  madden
 * Changed message level to SEV_WARNING
 *
 * Revision 6.34  1998/11/16 16:29:41  madden
 * Added ErrSetMessageLevel(SEV_INFO) and fixed call to PrintKAParametersExtra
 *
 * Revision 6.33  1998/09/28 12:33:02  madden
 * Used BlastErrorPrint
 *
 * Revision 6.31  1998/09/17 19:54:31  madden
 * Removed fillCandLambda
 *
 * Revision 6.29  1998/09/10 22:38:24  madden
 * Moved convertSeqAlignListToValNodeList and convertValNodeListToSeqAlignList
 *
 * Revision 6.28  1998/09/09 21:20:02  madden
 * AS fixed warnings, removed stderr fprintf, added call to PrintKAParametersExtra
 *
 * Revision 6.27  1998/09/09 16:10:49  madden
 * PHI-BLAST changes
 *
 * Revision 6.26  1998/07/17 15:41:36  madden
 * Added effective search space flag
 *
 * Revision 6.24  1998/06/14 19:44:46  madden
 * Checkpointing fix
 *
 * Revision 6.23  1998/06/12 21:32:27  madden
 * Removed extra FnPtr cast
 *
 * Revision 6.22  1998/06/10 13:33:39  madden
 * change -h from 0.01 to 0.001
 *
 * Revision 6.21  1998/06/05 21:48:43  madden
 * Added -K and -L options
 *
 * Revision 6.20  1998/05/18 18:01:05  madden
 * Changed args to allow filter options to be changed
 *
 * Revision 6.19  1998/05/01 18:31:03  egorov
 * Add new parametes to BLASTOptionSetGapParam()
 *
 * Revision 6.18  1998/04/30 14:32:33  madden
 * init_buff_ex arg changed to 90 for reference
 *
 * Revision 6.17  1998/04/29 14:29:31  madden
 * Made reference line longer
 *
 * Revision 6.16  1998/04/01 22:49:13  madden
 * Print No hits found message
 *
 * Revision 6.15  1998/02/25 20:50:50  madden
 * Added arg for db length
 *
 * Revision 6.14  1998/02/24 22:48:36  madden
 * Removed options for culling
 *
 * Revision 6.13  1998/01/02 14:33:37  madden
 * Added comma
 *
 * Revision 6.12  1997/12/31 17:48:56  madden
 * Added wordsize option
 *
 * Revision 6.11  1997/12/23 21:09:44  madden
 * Added -K and -L for range-dependent blast
 *
 * Revision 6.10  1997/12/23 20:57:44  madden
 * Changes for checkpointing
 *
 * Revision 6.9  1997/11/19 14:26:46  madden
 * Removed extra break statement
 *
 * Revision 6.8  1997/11/18 22:24:24  madden
 * Added call to BLASTOptionSetGapParams
 *
 * Revision 6.7  1997/11/08 22:04:43  madden
 * Called BlastOtherReturnsPrepare earlier to before posMatrix is deleted
 *
 * Revision 6.6  1997/10/27 22:26:49  madden
 * Added call to ObjMgrFreeCache(0)
 *
 * Revision 6.5  1997/10/23 20:26:15  madden
 * Use of init_buff_ex rather than init_buff
 *
 * Revision 6.4  1997/10/22 21:56:49  madden
 * Changed default values
 *
 * Revision 6.3  1997/10/07 21:33:34  madden
 * Added BLUNT option
 *
 * Revision 6.2  1997/09/18 22:25:02  madden
 * b and v options now work
 *
 * Revision 6.1  1997/09/16 16:40:50  madden
 * Dbinfo printing changed for multiple db searches
 *
 * Revision 6.0  1997/08/25 18:19:19  madden
 * Revision changed to 6.0
 *
 * Revision 1.24  1997/08/24 19:38:23  madden
 * Used function BlastOtherReturnsPrepare
 *
 * Revision 1.23  1997/08/14 21:48:57  madden
 * Added descriptions and alignments options
 *
 * Revision 1.22  1997/07/29 19:33:05  madden
 * Added TXALIGN_SHOW_QS flag
 *
 * Revision 1.21  1997/07/28 18:36:45  madden
 * Replaced printf with ErrPostEx and fprintf
 *
 * Revision 1.20  1997/07/28 16:59:21  madden
 * Added include for simutil.h
 *
 * Revision 1.19  1997/07/28 16:50:51  madden
 * Changed call to ShowTextAlignFromAnnot
 *
 * Revision 1.18  1997/07/22 19:18:40  madden
 * printing version, etc.
 *
 * Revision 1.17  1997/06/25 14:06:21  madden
 * minor changes to check convergence
 *
 * Revision 1.16  1997/06/23 20:51:29  madden
 * Added matrix option
 *
 * Revision 1.15  1997/06/20 19:30:00  madden
 * Added align_type info, support for SeqAligns
 *
 * Revision 1.14  1997/05/23 20:54:48  madden
 * Added -Z option for final gapped alignment
 *
 * Revision 1.13  1997/05/07 15:06:01  madden
 * replacement of search by compactSearch
 *
 * Revision 1.12  1997/04/21  14:09:27  madden
 * Added four new master-slave alignment types.
 *
 * Revision 1.11  1997/04/11  19:03:47  madden
 * Changes to ignore query ID and show master-slave alignments.
 *
 * Revision 1.10  1997/04/10  19:28:28  madden
 * COMMAND_LINE replaced by ALL_ROUNDS.
 *
 * Revision 1.9  1997/04/10  13:27:32  madden
 * Added COMMAND_LINE define, option for multiple alignments.
 *
 * Revision 1.8  1997/04/07  21:44:50  madden
 * Removed unused variable stats.
 *
 * Revision 1.7  1997/04/04  21:13:33  madden
 * Used function BioseqBlastEngineCore, remove PrivateScoreFunc.
 *
 * Revision 1.6  1997/04/04  16:38:04  madden
 * Changed filtering option, ObjMgrHold.
 *
 * Revision 1.5  1997/03/21  20:30:02  madden
 * Expect value arg made a float.
 *
 * Revision 1.4  1997/03/13  21:18:51  madden
 * Changed default extension value to 1 from 2.
 *
 * Revision 1.3  1997/02/19  21:43:04  madden
 * Extensive changes for psi-blast.  blastp runs now occur multiple times.
 *
 * Revision 1.2  1997/02/12  15:16:58  madden
 * Changed from blast2 to new formatting.
 *
 * Revision 1.1  1997/01/16  21:35:42  madden
 * Initial revision
 *
 * Revision 1.1  1997/01/16  21:34:23  madden
 * Initial revision
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
#include <txalign.h>
#include <simutil.h>
#include <gapxdrop.h>
#include <posit.h>
#include <seed.h>
#include <sqnutils.h>

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



/* #define NUMARG 41 */
#define NUMARG 40

static Args myargs [NUMARG] = {
  { "Database", 
	"nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
  { "Query File", 
	"stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  { "Multiple Hits window size (zero for single hit algorithm)", 
	"40", NULL, NULL, FALSE, 'A', ARG_INT, 0.0, 0, NULL},
  { "Threshold for extending hits", 
	"11", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},
  { "Expectation value (E)", 
	"10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
  { "alignment view options:\n0 = pairwise,\n1 = master-slave showing identities,\n2 = master-slave no identities,\n3 = flat master-slave, show identities,\n4 = flat master-slave, no identities,\n5 = master-slave no identities and blunt ends,\n6 = flat master-slave, no identities and blunt ends", 
        "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},
  { "Output File for Alignment", 
	"stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  { "Dropoff (X) for blast extensions in bits (default if zero)",
        "7.0", NULL, NULL, FALSE, 'y', ARG_FLOAT, 0.0, 0, NULL},
  { "0 for multiple hits 1-pass, 1 for single hit 1-pass, 2 for 2-pass",
	"0", NULL, NULL, FALSE, 'P', ARG_INT, 0.0, 0, NULL},
  { "Filter query sequence with SEG",
    "F", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
  { "Cost to open a gap",
	"11", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},
  { "Cost to extend a gap",
	"1", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},
  { "X dropoff value for gapped alignment (in bits)",
	"15", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
  { "Number of bits to trigger gapping",
        "22.0", NULL, NULL, FALSE, 'N', ARG_FLOAT, 0.0, 0, NULL},
  { "Gapped",
        "T", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
  { "Start of required region in query",
        "1", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},
  { "End of required region in query (-1 indicates end of query)",
        "-1", NULL, NULL, FALSE, 'H', ARG_INT, 0.0, 0, NULL},
  { "Number of processors to use",
	  "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},
  { "Show GI's in deflines", 
        "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},
  { "e-value threshold for inclusion in multipass model",
        "0.001", NULL, NULL, FALSE, 'h', ARG_FLOAT, 0.0, 0, NULL},
  { "Constant in pseudocounts for multipass version",
        "10", NULL, NULL, FALSE, 'c', ARG_INT, 0.0, 0, NULL},
  { "Maximum number of passes to use in  multipass version",
        "1", NULL, NULL, FALSE, 'j', ARG_INT, 0.0, 0, NULL},
  { "Believe the query defline", 
        "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
  { "X dropoff value for final gapped alignment (in bits)",
	"25", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},
  { "SeqAlign file ('Believe the query defline' must be TRUE)",
	NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},
  { "Matrix", 
	"BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
  { "Number of one-line descriptions (V)",
        "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},
  { "Number of alignments to show (B)",
        "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
 { "Output File for PSI-BLAST Checkpointing", 
	NULL, NULL, NULL, TRUE, 'C', ARG_FILE_OUT, 0.0, 0, NULL},
  { "Input File for PSI-BLAST Restart", 
	NULL, NULL, NULL, TRUE, 'R', ARG_FILE_IN, 0.0, 0, NULL},
  { "Word size, default if zero",
        "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},
  { "Effective length of the database (use zero for the real size)",
        "0", NULL, NULL, FALSE, 'z', ARG_INT, 0.0, 0, NULL},
  { "Number of best hits from a region to keep",
        "0", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},
  { "Length of region used to judge hits",
        "20", NULL, NULL, FALSE, 'L', ARG_INT, 0.0, 0, NULL},
  { "Effective length of the search space (use zero for the real size)",
        "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},
  { "program option for PHI-BLAST",
        "blastpgp", NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
  { "Hit File for PHI-BLAST", 
	"hit_file", NULL, NULL, FALSE, 'k', ARG_FILE_IN, 0.0, 0, NULL},
  { "Produce HTML output",
        "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  { "Output File for PSI-BLAST Matrix in ASCII", 
	NULL, NULL, NULL, TRUE, 'Q', ARG_FILE_OUT, 0.0, 0, NULL},
  { "Input Alignment File for PSI-BLAST Restart", 
	NULL, NULL, NULL, TRUE, 'B', ARG_FILE_IN, 0.0, 0, NULL}

  /*  { "Cost to decline alignment",
        "2", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL} */
}; 

Int2 Main (void)
 
{
	AsnIoPtr aip=NULL;
	BioseqPtr fake_bsp, query_bsp;
	BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
	BlastSearchBlkPtr search;
	BLAST_OptionsBlkPtr options;
	BlastPruneSapStructPtr prune;
	Boolean query_is_na, show_gi, believe_query=FALSE;
	Boolean html=FALSE;
	CharPtr params_buffer=NULL;
	Int4 number_of_descriptions, number_of_alignments;
	ObjectIdPtr obidp;
	SeqAlignPtr  head;
        SeqAnnotPtr annot, seqannot = NULL;
	SeqEntryPtr sep;
	SeqFeatPtr sfp;
	SeqLocPtr next, seqloc, seqloc_duplicate;
	TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
	Uint1 featureOrder[FEATDEF_ANY];
	Uint1 groupOrder[FEATDEF_ANY];
	Uint4 align_options, print_options;
	ValNodePtr  mask_loc, vnp, other_returns;
	/* used for psi-blast */
        Int4 thisPassNum;       
        posSearchItems *posSearch;
	compactSearchItems *compactSearch;
	Boolean  recoverCheckpoint = FALSE;
        Boolean  alreadyRecovered = FALSE;
	Boolean  freqCheckpoint = FALSE;
        Boolean  alignCheckpoint = FALSE;
        Boolean  checkReturn = FALSE;


	CharPtr blast_database, blast_inputfile, blast_outputfile;
	FILE *infp, *outfp;
        FILE *matrixfp = NULL; /*output file for ASCII version of PSI-BLAST matrix*/
 
        /*following declarations are for PHI-BLAST*/
        CharPtr patfile;
        FILE *patfp; 
        seedSearchItems *seedSearch;
        Int4 program_flag;
        Int4 queryLength; /*length of query sequence*/
        Boolean is_dna = FALSE; /*cannot use DNA queries in blastpgp*/
        Int4 i; /*index over characters*/
        Uint1Ptr query = NULL; /*query sequence read in*/
        Uint1Ptr unfilter_query = NULL; /*needed if seg will filter query*/
        SeqLocPtr seg_slp; /*pointer to structure for seg filtering*/
        ValNodePtr seedReturn; /*return value from seedEngineCore, which
                                 is a list of lists of SeqAligns, one
                                 list of SeqAligns per pattern occurrence*/
        ValNodePtr pruneSeed;  /*possibly reduced version of seedReturn*/
        SeqAlignPtr * lastSeqAligns; /*keeps track of the last SeqAlign in
                                     each list of seedReturn so that the
                                     2-level list can be converted to a 1-level
                                     list and then back to 2-level*/
        Int4 numLastSeqAligns;
        Int4 alignCount;
        


        if (! GetArgs ("blastpgp", NUMARG, myargs))
        {
                return (1);
        }


	UseLocalAsnloadDataAndErrMsg ();

	if (! SeqEntryLoad())
		return 1;

	ErrSetMessageLevel(SEV_WARNING);

        blast_database = myargs [0].strvalue;
        blast_inputfile = myargs [1].strvalue;
        blast_outputfile = myargs [6].strvalue;
	if (myargs[37].intvalue)
		html = TRUE;

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


	believe_query = FALSE;
	if (myargs[22].intvalue != 0)
		believe_query = TRUE;
	query_is_na = FALSE;

	if (myargs[24].strvalue != NULL)
	{
		if (believe_query == FALSE)
		{
			ErrPostEx(SEV_FATAL, 0, 0, "-J option must be TRUE to use this option");
			return (1);
		}
       		if ((aip = AsnIoOpen (myargs[24].strvalue,"w")) == NULL)
        	{
               	 ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", myargs[24].strvalue);
               	 return 1;
        	}
	}

	sep = FastaToSeqEntryEx(infp, query_is_na, NULL, believe_query);
	if (sep != NULL)
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
		
		options = BLASTOptionNew("blastp", (Boolean)myargs[14].intvalue);

		options->isPatternSearch = FALSE;
		if (0 != (StringCmp("blastpgp",myargs[35].strvalue))) {
		  options->isPatternSearch = TRUE;
		  program_flag = convertProgramToFlag(myargs[35].strvalue, &is_dna);
		}

		if (options->isPatternSearch) {
		  patfile = myargs[36].strvalue;
		  if ((patfp = FileOpen(patfile, "r")) == NULL)
		    {
		      ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open pattern file %s\n", patfile);
		      return (1);
		    }
		  seedSearch = (seedSearchItems *) ckalloc(sizeof(seedSearchItems));
		}


		/* Set default gap params for matrix. */
		BLASTOptionSetGapParams(options, myargs[25].strvalue, 0, 0);

		/* decrement by one to agree with program values. */
                options->required_start = myargs[15].intvalue - 1;
                options->required_end = myargs[16].intvalue;
                if (options->required_end != -1)
		{
			options->required_end--;
		}

                options->window_size = myargs [2].intvalue;

		options->threshold_second = (Int4) myargs [3].intvalue;


                options->dropoff_2nd_pass  = myargs [7].floatvalue;
                options->expect_value  = (Nlm_FloatHi) myargs [4].floatvalue;
		number_of_descriptions = myargs[26].intvalue;
		number_of_alignments = myargs[27].intvalue;
		options->hitlist_size = MAX(number_of_descriptions, number_of_alignments);
		if (myargs[14].intvalue != 0)
		{
			if (myargs[8].intvalue == 0)
			{
                		options->two_pass_method  = FALSE;
               			options->multiple_hits_only  = TRUE;
			}
			else if (myargs[8].intvalue == 1)
			{
               	 		options->two_pass_method  = FALSE;
               			options->multiple_hits_only  = FALSE;
			}
			else
			{
                		options->two_pass_method  = TRUE;
               			options->multiple_hits_only  = FALSE;
			}
			options->gap_open = myargs[10].intvalue;
			options->gap_extend = myargs[11].intvalue;
                        options->decline_align = INT2_MAX;
                        /* options->decline_align = myargs[40].intvalue; */
			options->gap_x_dropoff = myargs[12].intvalue;
			options->gap_x_dropoff_final = myargs[23].intvalue;
			options->gap_trigger = myargs[13].floatvalue;
		}
		if(options->isPatternSearch)
                  fillCandLambda(seedSearch, myargs[25].strvalue, options);

		if (StringICmp(myargs[9].strvalue, "T") == 0)
		{
			options->filter_string = StringSave("S");
		}
		else
		{
			options->filter_string = StringSave(myargs[9].strvalue);
		}
		options->number_of_cpus = (Int2) myargs[17].intvalue;
		show_gi = (Boolean) myargs[18].intvalue;
                options->ethresh = (Nlm_FloatHi) myargs[19].floatvalue;
                options->pseudoCountConst = myargs[20].intvalue;
                options->maxNumPasses = myargs[21].intvalue;
                /*zero out e-value threshold if it will not be used*/
                if (options->maxNumPasses == 1)
                  options->ethresh = 0.0;
		if (myargs[30].intvalue)
			options->wordsize = myargs[30].intvalue;
		if (myargs[31].floatvalue)
			options->db_length = (Int8) myargs[31].floatvalue;

		options->hsp_range_max  = myargs[32].intvalue;
                if (options->hsp_range_max != 0)
                        options->perform_culling = TRUE;
                options->block_width  = myargs[33].intvalue;
		if (myargs[34].floatvalue)
			options->searchsp_eff = (Nlm_FloatHi) myargs[34].floatvalue;

		options = BLASTOptionValidate(options, "blastp");

		if (options == NULL)
		{
			return 1;
		}

		if (believe_query == TRUE)
		{
			fake_bsp = query_bsp;
		}
		else
		{
			fake_bsp = BioseqNew();
			fake_bsp->descr = query_bsp->descr;
			fake_bsp->repr = query_bsp->repr;
			fake_bsp->mol = query_bsp->mol;
			fake_bsp->length = query_bsp->length;
			fake_bsp->seq_data_type = query_bsp->seq_data_type;
			fake_bsp->seq_data = query_bsp->seq_data;
	
			obidp = ObjectIdNew();
               		obidp->str = StringSave("QUERY");
                	ValNodeAddPointer(&(fake_bsp->id), SEQID_LOCAL, obidp);
			
			/* FASTA defline not parsed, ignore the "lcl|tempseq" ID. */
			query_bsp->id = SeqIdSetFree(query_bsp->id);
		}
 

		if (options->isPatternSearch) {
		  query = BlastGetSequenceFromBioseq(fake_bsp,&queryLength);
		  seg_slp = BlastBioseqFilter(fake_bsp, options->filter_string);
		  if (seg_slp)
		    {
		      unfilter_query = MemNew((queryLength + 1) * sizeof(Uint1));
		      for (i = 0; i < queryLength; i++)
			unfilter_query[i] = query[i];
		      BlastMaskTheResidues(query,queryLength,21,seg_slp,FALSE, 0);
		    }
		}



		search = BLASTSetUpSearchWithReadDb(fake_bsp, "blastp", query_bsp->length, blast_database, options, NULL);

		if (search == NULL)
			return 1;

                /*AAS*/
 		if ((NULL != myargs[29].strvalue) || (NULL != myargs[39].strvalue)) {
		  recoverCheckpoint = TRUE;
                  if (NULL != myargs[29].strvalue) {
		    freqCheckpoint = TRUE;
		    alignCheckpoint = FALSE;
		  }
		  else {
		    freqCheckpoint = FALSE;
		    alignCheckpoint = TRUE;
		  }
		}
                if (recoverCheckpoint)
		  search->positionBased = TRUE;
		else
		  search->positionBased = FALSE;

		global_fp = outfp;

		if (html)
			fprintf(outfp, "<PRE>\n");
                init_buff_ex(90);
		BlastPrintVersionInfo("blastp", html, outfp);
		fprintf(outfp, "\n");
		BlastPrintReference(html, 90, outfp);
		fprintf(outfp, "\n");
		AcknowledgeBlastQuery(query_bsp, 70, outfp, believe_query, html);
		PrintDbInformation(blast_database, TRUE, 70, outfp, html);
		free_buff();

		MemSet((Pointer)(groupOrder), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
		MemSet((Pointer)(featureOrder), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));

		posSearch = NULL;
                thisPassNum = 0;
		compactSearch = NULL;
                search->posConverged = FALSE;
		global_fp = outfp;
                search->error_return = NULL;
                /*AAS*/
		if (recoverCheckpoint) {
		  posSearch = (posSearchItems *) MemNew(1 * sizeof(posSearchItems));
                  compactSearch = compactSearchNew(compactSearch);
                  copySearchItems(compactSearch, search);
                  posInitializeInformation(posSearch,search);
		  /*AAS*/
		  if (freqCheckpoint) {
		    checkReturn = posReadCheckpoint(posSearch, compactSearch, myargs[29].strvalue, &(search->error_return));
		    search->sbp->posMatrix = posSearch->posMatrix;
		  }
		  else {
		    search->sbp->posMatrix = BposComputation(posSearch, search, compactSearch, myargs[39].strvalue, myargs[28].strvalue, &(search->error_return)); 
                    if (NULL == search->sbp->posMatrix)
		      checkReturn = FALSE;
		    else
                      checkReturn = TRUE;
		  }
		  BlastErrorPrint(search->error_return);
                  if (!checkReturn) {
		    ErrPostEx(SEV_FATAL, 0, 0, "blast: Error recovering from checkpoint");
		    return 1;
		  }
		  if (myargs[38].strvalue != NULL) {
		    if ((matrixfp = FileOpen(myargs[38].strvalue, "w")) == NULL)
		      {
			ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open matrix output file %s\n", myargs[38].strvalue);
			return (1);
		      }
		    outputPosMatrix(posSearch, compactSearch, matrixfp); /* a diagnostic, partly an option with -Q. */
                    if (NULL != matrixfp)
		      FileClose(matrixfp);
		  }
		}

                do {  /*AAS*/

                  thisPassNum++;
                  if (thisPassNum > 1)
                    options->isPatternSearch = FALSE;

#ifdef OS_UNIX
		  search->tick_callback =  tick_callback;
		  fprintf(global_fp, "%s", "Searching");
		  fflush(global_fp);
#endif
                  if (ALL_ROUNDS && 1 == thisPassNum && (!recoverCheckpoint)) 
                    posSearch = (posSearchItems *) MemNew(1 * sizeof(posSearchItems));


                  posSearch->threshSequences = NULL;
		  if (options->isPatternSearch) {
                    search->gap_align = GapAlignBlkNew(1,1);
                    search->gap_align->gap_open = options->gap_open;
                    search->gap_align->gap_extend = options->gap_extend;
                    search->gap_align->decline_align = (-(BLAST_SCORE_MIN));
                    search->gap_align->x_parameter = options->gap_x_dropoff
                          *NCBIMATH_LN2/seedSearch->paramLambda;
                    search->gap_align->matrix = search->sbp->matrix;
                    initProbs(seedSearch);
                    init_order(search->gap_align->matrix,program_flag,is_dna,seedSearch);
                    for(i = 0; i < queryLength; i++)
		      query[i] = seedSearch->order[query[i]];
                    if (unfilter_query)
		      for(i = 0; i < queryLength; i++)
			unfilter_query[i] = seedSearch->order[unfilter_query[i]];
		    if ((1 == thisPassNum) && (!recoverCheckpoint))
		    {
		        seqloc = NULL;
		        seedReturn = seedEngineCore(search, options, query, 
                          unfilter_query, blast_database, patfile, program_flag, 
                          patfp,  outfp, is_dna, FALSE, seedSearch, 
                          options->ethresh, myargs[34].floatvalue, posSearch, &seqloc, TRUE);
			BlastErrorPrint(search->error_return);
                        seqloc_duplicate = seqloc;
			head = convertValNodeListToSeqAlignList(seedReturn, &lastSeqAligns, &numLastSeqAligns);
			MemSet((Pointer)(groupOrder), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
			MemSet((Pointer)(featureOrder), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
			featureOrder[FEATDEF_REGION] = 1;
			groupOrder[FEATDEF_REGION] = 1;
			annot = fake_bsp->annot = SeqAnnotNew();
			fake_bsp->annot->type = 1;	/* ftable. */
			while (seqloc)
			{
				next = seqloc->next;
				sfp = SeqFeatNew();
				sfp->location = seqloc;
				sfp->data.choice = SEQFEAT_REGION;
				sfp->data.value.ptrvalue = StringSave("pattern");
				annot->data = sfp;
				if (next)
				{
					annot->next = SeqAnnotNew();
					annot = annot->next;
					annot->type = 1;
				}
				seqloc = next;
			}
		    }
		    else
		    {
		      ErrPostEx(SEV_FATAL, 0, 0, "blast: should not be using PHI-BLAST in this way");
		    }
		  }
		  else {
		    if ((1 == thisPassNum) && (!recoverCheckpoint))
		      head = BioseqBlastEngineCore(search, options, NULL);
		    else
		      head = BioseqBlastEngineCore(search, options, search->sbp->posMatrix);  
		  }
		  if (recoverCheckpoint) {
                    compactSearchDestruct(compactSearch);
                    recoverCheckpoint = FALSE;
                    alreadyRecovered = TRUE;
                  }


                  compactSearch = compactSearchNew(compactSearch);
                  copySearchItems(compactSearch, search);

/* The next two calls (after the "if") are diagnostics for Stephen. */
/* Don't perform this if only one pass will be made (i.e., standard BLAST) */
		  if (ALL_ROUNDS && 1 != search->pbp->maxNumPasses)
		  {
                  	if ((1 == thisPassNum)  && (!alreadyRecovered))
		    		posInitializeInformation(posSearch, search);
		        posPrintInformation(posSearch, search, thisPassNum);
		  }



#ifdef OS_UNIX
                fprintf(global_fp, "%s", "done");
#endif

/* AAS */
                if (1 == thisPassNum) {
		  ReadDBBioseqFetchEnable ("blastpgp", blast_database, FALSE, TRUE);
		}
		/* Have things converged? */
		if (ALL_ROUNDS && 1 != search->pbp->maxNumPasses)
		{
		  posConvergenceTest(posSearch, search, head, thisPassNum);
                }


                
                print_options = 0;
		if (search->pbp->gapped_calculation == FALSE)
		  print_options += TXALIGN_SHOW_NO_OF_SEGS;
                align_options = 0;
                align_options += TXALIGN_COMPRESS;
                align_options += TXALIGN_END_NUM;
                if (show_gi) {
		   align_options += TXALIGN_SHOW_GI;
		   print_options += TXALIGN_SHOW_GI;
		 } 

		if (html)
		{
			align_options += TXALIGN_HTML;
			print_options += TXALIGN_HTML;
		}

		if (myargs[5].intvalue != 0)
		{
		   	align_options += TXALIGN_MASTER;
			if (myargs[5].intvalue == 1 || myargs[5].intvalue == 3)
		   		align_options += TXALIGN_MISMATCH;
			if (myargs[5].intvalue == 3 || myargs[5].intvalue == 4 || myargs[5].intvalue == 6)
		   		align_options += TXALIGN_FLAT_INS;
			if (myargs[5].intvalue == 5 || myargs[5].intvalue == 6)
                                align_options += TXALIGN_BLUNT_END;
		}
		else
		{
			align_options += TXALIGN_MATRIX_VAL;
			align_options += TXALIGN_SHOW_QS;
		}

/*AAS*/
                search->positionBased = TRUE;
		other_returns = BlastOtherReturnsPrepare(search);
                if (alreadyRecovered) {
		  posCheckpointFreeMemory(posSearch, compactSearch->qlength);
                  alreadyRecovered = FALSE;
                }
                if (ALL_ROUNDS && thisPassNum > 1) {
		  posCleanup(posSearch, compactSearch);
                }
                if (!search->posConverged && (search->pbp->maxNumPasses == 0 || (thisPassNum < search->pbp->maxNumPasses)))
                  if (ALL_ROUNDS) {
		    search->sbp->posMatrix = CposComputation(posSearch, search, compactSearch, head, myargs[28].strvalue, (options->isPatternSearch && (1== thisPassNum)), &(search->error_return)); /*AAS*/
		    BlastErrorPrint(search->error_return);
		  }
                  else
		    search->sbp->posMatrix = WposComputation(compactSearch, head); /*AAS*/
        	if (head)
        	{
			if (seqannot)
				seqannot = SeqAnnotFree(seqannot);
			seqannot = SeqAnnotNew();
               	 	seqannot->type = 2;
			AddAlignInfoToSeqAnnot(seqannot, 2);
                	seqannot->data = head;
			if (outfp)
			{
			        if (1 != search->pbp->maxNumPasses)
				  fprintf(outfp, "\nResults from round %d\n", thisPassNum);


				ObjMgrSetHold();
                                init_buff_ex(85);
				if (1 == thisPassNum)
				  search->positionBased = FALSE;
				if (1 == thisPassNum) {
                                  if (!options->isPatternSearch) {
				    prune = BlastPruneHitsFromSeqAlign(head, number_of_descriptions, NULL);
				    PrintDefLinesFromSeqAlign(prune->sap, 80, outfp, print_options, FIRST_PASS, NULL);
				  }
                                  else {
				    seedReturn = convertSeqAlignListToValNodeList(head,lastSeqAligns, numLastSeqAligns);

                                    pruneSeed = SeedPruneHitsFromSeedReturn(seedReturn, number_of_descriptions);
				    PrintDefLinesExtra(pruneSeed, 80, outfp, print_options, FIRST_PASS, NULL, seqloc_duplicate);
				  }
				}
                                else {
				  prune = BlastPruneHitsFromSeqAlign(head, number_of_descriptions, NULL);
                                  if (ALL_ROUNDS) {
				    PrintDefLinesFromSeqAlign(prune->sap, 80, outfp, print_options, NOT_FIRST_PASS_REPEATS, posSearch->posRepeatSequences);
				    PrintDefLinesFromSeqAlign(prune->sap, 80, outfp, print_options, NOT_FIRST_PASS_NEW, posSearch->posRepeatSequences);
				  }
				  else
				  PrintDefLinesFromSeqAlign(prune->sap, 80, outfp, print_options, FIRST_PASS, NULL);
				}
				if (ALL_ROUNDS && search->posConverged)
				  fprintf(outfp, "\nCONVERGED!\n");
				free_buff();

                                if (!(options->isPatternSearch)) {
				  prune = BlastPruneHitsFromSeqAlign(head, number_of_alignments, prune);

				  seqannot->data = prune->sap;
				  if (myargs[5].intvalue != 0)
				    ShowTextAlignFromAnnot(seqannot, 60, outfp, featureOrder, groupOrder, align_options, NULL, search->mask, NULL);
				  else
				    ShowTextAlignFromAnnot(seqannot, 60, outfp, featureOrder, groupOrder, align_options, NULL, search->mask, FormatScoreFunc);
				}
                                else {
                                  if (number_of_alignments < number_of_descriptions)
                                    pruneSeed = SeedPruneHitsFromSeedReturn(pruneSeed, number_of_alignments);
				  if (myargs[5].intvalue != 0) 
				    ShowTextAlignFromAnnotExtra(query_bsp, pruneSeed, seqloc_duplicate, 60, outfp, featureOrder, groupOrder, align_options, NULL, search->mask, NULL);
                                  else
				    ShowTextAlignFromAnnotExtra(query_bsp, pruneSeed, seqloc_duplicate, 60, outfp, featureOrder, groupOrder, align_options, NULL, search->mask, FormatScoreFunc);
				    }
				seqannot->data = head;
                                if (!(options->isPatternSearch)) 
					prune = BlastPruneSapStructDestruct(prune);
				search->positionBased = TRUE;
				ObjMgrClearHold();
				ObjMgrFreeCache(0);
			}
		}
                else
                {
                        fprintf(outfp, "\n\n ***** No hits found ******\n\n");
                }


                if (ALL_ROUNDS && thisPassNum > 1) {
		  MemFree(posSearch->posRepeatSequences);
                }
		if (!search->posConverged && (0 == search->pbp->maxNumPasses || thisPassNum < search->pbp->maxNumPasses)) {
		  if (ALL_ROUNDS && (myargs[38].strvalue != NULL)) {
		    if ((matrixfp = FileOpen(myargs[38].strvalue, "w")) == NULL)
		      {
			ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open matrix output file %s\n", myargs[38].strvalue);
			return (1);
		      }
		  }
		  outputPosMatrix(posSearch, compactSearch, matrixfp); /* a diagnostic, partly an option with -Q. */
                    if (NULL != matrixfp)
		      FileClose(matrixfp);
		}
		  
		} while (( 0 == search->pbp->maxNumPasses || thisPassNum < (search->pbp->maxNumPasses)) && (! (search->posConverged)));

		if (aip && seqannot)
		{
			 SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, aip, NULL);
			 AsnIoReset(aip);
                         aip = AsnIoClose(aip);
		}


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

                init_buff_ex(85);
                dbinfo_head = dbinfo;
                while (dbinfo)
                {
                        PrintDbReport(dbinfo, 70, outfp);
                        dbinfo = dbinfo->next;
                }
                dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);

		if (ka_params && !options->isPatternSearch)
		{
                	PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, outfp, FALSE);
			MemFree(ka_params);
		}

		if (ka_params_gap)
		{
			if (options->isPatternSearch)
                		PrintKAParametersExtra(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, seedSearch->paramC, 70, outfp, FALSE);
			else
                		PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, outfp, FALSE);
			MemFree(ka_params_gap);
		}

                PrintTildeSepLines(params_buffer, 70, outfp);
                MemFree(params_buffer);
                free_buff();

		ReadDBBioseqFetchDisable();
		if (options->isPatternSearch) {
		  if (NULL != query)
		    MemFree(query);
		  if (NULL != unfilter_query)
		    unfilter_query = MemFree(unfilter_query);
		  seedSearch = MemFree(seedSearch);
		  FileClose(patfp);
		}

        }
	if (ALL_ROUNDS)
	  posSearch = MemFree(posSearch);
        compactSearchDestruct(compactSearch);
	options = BLASTOptionDelete(options);
	sep = SeqEntryFree(sep);
	search = BlastSearchBlkDestruct(search);
	if (believe_query == FALSE)
	{
		fake_bsp->descr = NULL;
		fake_bsp->length = 0;
		fake_bsp->seq_data = NULL;
	}
	FileClose(infp);

	return 0;
}
	


