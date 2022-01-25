/* $Id: blastpgp.c,v 6.58 2000/01/07 22:01:04 shavirin Exp $ */
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
/* $Revision: 6.58 $ */ 
/* $Log: blastpgp.c,v $
/* Revision 6.58  2000/01/07 22:01:04  shavirin
/* Fixed problem with printing alignment.
/*
/* Revision 6.57  1999/12/21 21:34:05  shavirin
/* Fixed some memory leaks.
/*
/* Revision 6.56  1999/11/08 19:12:41  shavirin
/* Fixed minor SGI warning.
/*
/* Revision 6.55  1999/10/21 20:27:52  shavirin
/* Fixed bug resulted in failue to print out seqannot. (-O)
/*
/* Revision 6.53  1999/10/14 15:52:51  shavirin
/* Added possibility to make search by list of gis. Fixed some bugs.
/*
/* Revision 6.52  1999/10/05 17:41:08  shavirin
/* Corrected in accordance to blast.c changes.
/*
/* Revision 6.51  1999/09/24 16:06:15  shavirin
/* Matrix is set to NULL if no matrix calculation produced.
/*
/* Revision 6.50  1999/09/22 17:53:25  shavirin
/* Now functions will collect messages in ValNodePtr before printing out.
/*
/* Revision 6.49  1999/09/21 16:01:46  shavirin
/* Rearanged all file. Main function was disintegrated into few small
/* functions.
/*
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

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs [] = {
    { "Database",               /* 0 */
      "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
    { "Query File",             /* 1 */
      "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
    { "Multiple Hits window size (zero for single hit algorithm)", /* 2 */
      "40", NULL, NULL, FALSE, 'A', ARG_INT, 0.0, 0, NULL},
    { "Threshold for extending hits", /* 3 */
      "11", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},
    { "Expectation value (E)",  /* 4 */
      "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
    { "alignment view options:\n0 = pairwise,\n1 = master-slave showing identities,\n2 = master-slave no identities,\n3 = flat master-slave, show identities,\n4 = flat master-slave, no identities,\n5 = master-slave no identities and blunt ends,\n6 = flat master-slave, no identities and blunt ends", /* 5 */
      "0", NULL, NULL, FALSE, 'm', ARG_INT, 0.0, 0, NULL},
    { "Output File for Alignment", /* 6 */
      "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Dropoff (X) for blast extensions in bits (default if zero)", /* 7 */
      "7.0", NULL, NULL, FALSE, 'y', ARG_FLOAT, 0.0, 0, NULL},
    { "0 for multiple hits 1-pass, 1 for single hit 1-pass, 2 for 2-pass", /* 8 */
      "0", NULL, NULL, FALSE, 'P', ARG_INT, 0.0, 0, NULL},
    { "Filter query sequence with SEG", /* 9 */
      "F", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
    { "Cost to open a gap",     /* 10 */
      "11", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},
    { "Cost to extend a gap",   /* 11 */
      "1", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},
    { "X dropoff value for gapped alignment (in bits)", /* 12 */
      "15", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
    { "Number of bits to trigger gapping", /* 13 */
      "22.0", NULL, NULL, FALSE, 'N', ARG_FLOAT, 0.0, 0, NULL},
    { "Gapped",                 /* 14 */
      "T", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Start of required region in query", /* 15 */
      "1", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},
    { "End of required region in query (-1 indicates end of query)", /* 16 */
      "-1", NULL, NULL, FALSE, 'H', ARG_INT, 0.0, 0, NULL},
    { "Number of processors to use", /* 17 */
      "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},
    { "Show GI's in deflines",  /* 18 */
      "F", NULL, NULL, FALSE, 'I', ARG_BOOLEAN, 0.0, 0, NULL},
    { "e-value threshold for inclusion in multipass model", /* 19 */
      "0.001", NULL, NULL, FALSE, 'h', ARG_FLOAT, 0.0, 0, NULL},
    { "Constant in pseudocounts for multipass version", /* 20 */
      "10", NULL, NULL, FALSE, 'c', ARG_INT, 0.0, 0, NULL},
    { "Maximum number of passes to use in  multipass version", /* 21 */
      "1", NULL, NULL, FALSE, 'j', ARG_INT, 0.0, 0, NULL},
    { "Believe the query defline", /* 22 */
      "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
    { "X dropoff value for final gapped alignment (in bits)", /* 23 */
      "25", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},
    { "SeqAlign file ('Believe the query defline' must be TRUE)", /* 24 */
      NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Matrix",                 /* 25 */
      "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
    { "Number of one-line descriptions (V)", /* 26 */
      "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},
    { "Number of alignments to show (B)", /* 27 */
      "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
    { "Output File for PSI-BLAST Checkpointing", /* 28 */
      NULL, NULL, NULL, TRUE, 'C', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Input File for PSI-BLAST Restart", /* 29 */
      NULL, NULL, NULL, TRUE, 'R', ARG_FILE_IN, 0.0, 0, NULL},
    { "Word size, default if zero", /* 30 */
      "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},
    { "Effective length of the database (use zero for the real size)", /* 31 */
      "0", NULL, NULL, FALSE, 'z', ARG_INT, 0.0, 0, NULL},
    { "Number of best hits from a region to keep", /* 32 */
      "0", NULL, NULL, FALSE, 'K', ARG_INT, 0.0, 0, NULL},
    { "Length of region used to judge hits", /* 33 */
      "20", NULL, NULL, FALSE, 'L', ARG_INT, 0.0, 0, NULL},
    { "Effective length of the search space (use zero for the real size)", /* 34 */
      "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},
    { "program option for PHI-BLAST", /* 35 */
      "blastpgp", NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
    { "Hit File for PHI-BLAST", /* 36 */
      "hit_file", NULL, NULL, FALSE, 'k', ARG_FILE_IN, 0.0, 0, NULL},
    { "Produce HTML output",    /* 37 */
      "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Output File for PSI-BLAST Matrix in ASCII", /* 38 */
      NULL, NULL, NULL, TRUE, 'Q', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Input Alignment File for PSI-BLAST Restart", /* 39 */
      NULL, NULL, NULL, TRUE, 'B', ARG_FILE_IN, 0.0, 0, NULL},
    { "Restrict search of database to list of GI's", /* 40 */
	NULL, NULL, NULL, TRUE, 'l', ARG_STRING, 0.0, 0, NULL},
    /*    { "Cost to decline alignment",  41 
          "10000", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL}  */
};

typedef struct _pgp_blast_options {
    BLAST_OptionsBlkPtr options;
    CharPtr blast_database;
    BioseqPtr query_bsp, fake_bsp;
    Int4 number_of_descriptions, number_of_alignments;
    FILE *infp, *outfp;
    AsnIoPtr aip_out;
    Boolean html;
    Boolean believe_query;
    Uint4 align_options, print_options;
    /* PHI-PSI Blast variables */
    Uint1 featureOrder[FEATDEF_ANY];
    Uint1 groupOrder[FEATDEF_ANY];
    Int4 program_flag;
    CharPtr patfile;
    FILE *patfp; 
    seedSearchItems *seedSearch;
} PGPBlastOptions, PNTR PGPBlastOptionsPtr;

void PGPGetPrintOptions(Boolean gapped, Uint4Ptr align_options_out, 
                        Uint4Ptr print_options_out)
{
    Uint4 print_options, align_options;

    print_options = 0;
    if (gapped == FALSE)
        print_options += TXALIGN_SHOW_NO_OF_SEGS;
    
    align_options = 0;
    align_options += TXALIGN_COMPRESS;
    align_options += TXALIGN_END_NUM;

    if (myargs[18].intvalue) {
        align_options += TXALIGN_SHOW_GI;
        print_options += TXALIGN_SHOW_GI;
    } 
    
    if (myargs[37].intvalue) {
        align_options += TXALIGN_HTML;
        print_options += TXALIGN_HTML;
    }
    
    if (myargs[5].intvalue != 0) {
        align_options += TXALIGN_MASTER;
        if (myargs[5].intvalue == 1 || myargs[5].intvalue == 3)
            align_options += TXALIGN_MISMATCH;
        if (myargs[5].intvalue == 3 || myargs[5].intvalue == 4 || myargs[5].intvalue == 6)
            align_options += TXALIGN_FLAT_INS;
        if (myargs[5].intvalue == 5 || myargs[5].intvalue == 6)
            align_options += TXALIGN_BLUNT_END;
    } else {
        align_options += TXALIGN_MATRIX_VAL;
        align_options += TXALIGN_SHOW_QS;
    }

    *align_options_out = align_options;
    *print_options_out = print_options;

    return;
}
void PGPFreeBlastOptions(PGPBlastOptionsPtr bop)
{
    bop->options = BLASTOptionDelete(bop->options);
    MemFree(bop->blast_database);
    MemFree(bop->patfile);
    MemFree(bop->seedSearch);
    MemFree(bop);

    return;
} 
    
PGPBlastOptionsPtr PGPReadBlastOptions(void)
{
    PGPBlastOptionsPtr bop;
    BLAST_OptionsBlkPtr options;
    SeqEntryPtr sep;
    Boolean is_dna;
    ObjectIdPtr obidp;

    bop = MemNew(sizeof(PGPBlastOptions));
    
    bop->blast_database   = StringSave(myargs [0].strvalue);

    if ((bop->infp = FileOpen(myargs [1].strvalue, "r")) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", 
                  myargs [1].strvalue);
        return NULL;
    }
    
    if (myargs [6].strvalue != NULL) {
        if ((bop->outfp = FileOpen(myargs [6].strvalue, "w")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output "
                      "file %s\n", myargs [6].strvalue);
            return NULL;
        }
    }
    
    bop->number_of_descriptions = myargs[26].intvalue;
    bop->number_of_alignments = myargs[27].intvalue;
    
    if (myargs[22].intvalue != 0)
        bop->believe_query = TRUE;
    
    if (myargs[24].strvalue != NULL) {
        
        if (bop->believe_query == FALSE) {
            ErrPostEx(SEV_FATAL, 0, 0, 
                      "-J option must be TRUE to use this option");
            return NULL;
        }
        
        if ((bop->aip_out = AsnIoOpen (myargs[24].strvalue,"w")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output "
                      "file %s\n", myargs[24].strvalue);
            return NULL;
        }
    }
    
    if((sep = FastaToSeqEntryEx(bop->infp, FALSE, NULL, 
                                bop->believe_query)) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "Unable to read input FASTA file\n");
        return NULL;
    }
    
    SeqEntryExplore(sep, &bop->query_bsp, FindProt);    
    sep->data.ptrvalue = NULL;
    SeqEntryFree(sep);
    
    if (bop->query_bsp == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
        return NULL;
    }    
    
    options = BLASTOptionNew("blastp", (Boolean)myargs[14].intvalue);
    bop->options = options;

    /* Set default gap params for matrix. */
    BLASTOptionSetGapParams(options, myargs[25].strvalue, 0, 0);

    PGPGetPrintOptions(options->gapped_calculation, &bop->align_options, 
                       &bop->print_options);

    /* decrement by one to agree with program values. */
    options->required_start = myargs[15].intvalue - 1;
    options->required_end = myargs[16].intvalue;
    if (options->required_end != -1) {
        options->required_end--;
    }
    
    options->window_size = myargs [2].intvalue;
    
    options->threshold_second = (Int4) myargs [3].intvalue;
    
    options->dropoff_2nd_pass  = myargs [7].floatvalue;
    options->expect_value  = (Nlm_FloatHi) myargs [4].floatvalue;
    options->hitlist_size = MAX(bop->number_of_descriptions, 
                                bop->number_of_alignments);
    
    if (myargs[14].intvalue != 0) {
        if (myargs[8].intvalue == 0) {
            options->two_pass_method  = FALSE;
            options->multiple_hits_only  = TRUE;
        } else if (myargs[8].intvalue == 1) {
            options->two_pass_method  = FALSE;
            options->multiple_hits_only  = FALSE;
        } else {
            options->two_pass_method  = TRUE;
            options->multiple_hits_only  = FALSE;
        }
        options->gap_open = myargs[10].intvalue;
        options->gap_extend = myargs[11].intvalue;

        options->decline_align = INT2_MAX;
        /* options->decline_align = myargs[41].intvalue; */

        options->gap_x_dropoff = myargs[12].intvalue;
        options->gap_x_dropoff_final = myargs[23].intvalue;
        options->gap_trigger = myargs[13].floatvalue;
    }
    
    if (StringICmp(myargs[9].strvalue, "T") == 0) {
        options->filter_string = StringSave("S");
    } else {
        options->filter_string = StringSave(myargs[9].strvalue);
    }
    
    options->number_of_cpus = (Int2) myargs[17].intvalue;
    
    
    options->isPatternSearch = FALSE;
    
    if (0 != (StringCmp("blastpgp",myargs[35].strvalue))) {
        options->isPatternSearch = TRUE;
        bop->program_flag = convertProgramToFlag(myargs[35].strvalue, 
                                                 &is_dna);
    }
    
    if (options->isPatternSearch) {
        bop->patfile = StringSave(myargs[36].strvalue);
        if ((bop->patfp = FileOpen(bop->patfile, "r")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open pattern "
                      "file %s\n", bop->patfile);
            return NULL;
        }
        
        bop->seedSearch = (seedSearchItems *) 
            ckalloc(sizeof(seedSearchItems));
    }
    
    if(options->isPatternSearch)
        fillCandLambda(bop->seedSearch, myargs[25].strvalue, options);
    
    options->ethresh = (Nlm_FloatHi) myargs[19].floatvalue;
    options->pseudoCountConst = myargs[20].intvalue;
    options->maxNumPasses = myargs[21].intvalue;
    /*zero out e-value threshold if it will not be used*/
    if (options->maxNumPasses == 1)
        options->ethresh = 0.0;
    if (myargs[30].intvalue)
        options->wordsize = myargs[30].intvalue;
    if (myargs[31].intvalue)
        options->db_length = (Int8) myargs[31].intvalue;
    
    options->hsp_range_max  = myargs[32].intvalue;
    if (options->hsp_range_max != 0)
        options->perform_culling = TRUE;
    options->block_width  = myargs[33].intvalue;
    
    if (myargs[34].floatvalue)
        options->searchsp_eff = (Nlm_FloatHi) myargs[34].floatvalue;

    /* Seting list of gis to restrict search */
    
    if (myargs[40].strvalue) {
        options->gifile = StringSave(myargs[40].strvalue);
    }
    
    options = BLASTOptionValidate(options, "blastp");
    
    if (options == NULL)
        return NULL;

    if (bop->believe_query == TRUE) {
        bop->fake_bsp = bop->query_bsp;
    } else {
        bop->fake_bsp = BioseqNew();
        bop->fake_bsp->descr = bop->query_bsp->descr;
        bop->fake_bsp->repr = bop->query_bsp->repr;
        bop->fake_bsp->mol = bop->query_bsp->mol;
        bop->fake_bsp->length = bop->query_bsp->length;
        bop->fake_bsp->seq_data_type = bop->query_bsp->seq_data_type;
        bop->fake_bsp->seq_data = bop->query_bsp->seq_data;
        
        obidp = ObjectIdNew();
        obidp->str = StringSave("QUERY");
        ValNodeAddPointer(&(bop->fake_bsp->id), SEQID_LOCAL, obidp);
        
        /* FASTA defline not parsed, ignore the "lcl|tempseq" ID. */
        bop->query_bsp->id = SeqIdSetFree(bop->query_bsp->id);
    }
    
    return bop;
}
Boolean PGPFormatHeader(PGPBlastOptionsPtr bop)
{
    Boolean html = myargs[37].intvalue;

    if (html)
        fprintf(bop->outfp, "<PRE>\n");
    
    init_buff_ex(90);
    BlastPrintVersionInfo("blastp", html, bop->outfp);
    fprintf(bop->outfp, "\n");
    BlastPrintReference(html, 90, bop->outfp);
    fprintf(bop->outfp, "\n");
    AcknowledgeBlastQuery(bop->query_bsp, 70, 
                          bop->outfp, bop->believe_query, html);
    PrintDbInformation(bop->blast_database, TRUE, 70, bop->outfp, html);
    free_buff();

    return TRUE;
}
Boolean  PGPFormatFooter(PGPBlastOptionsPtr bop, BlastSearchBlkPtr search)
{
    
    ValNodePtr  mask_loc, vnp;
    BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
    TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
    CharPtr params_buffer=NULL;
    ValNodePtr other_returns;
    BLAST_MatrixPtr blast_matrix;

    other_returns = BlastOtherReturnsPrepare(search);

    mask_loc = NULL;
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
            blast_matrix = vnp->data.ptrvalue;
            BLAST_MatrixDestruct(blast_matrix);
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
    while (dbinfo) {
        PrintDbReport(dbinfo, 70, bop->outfp);
        dbinfo = dbinfo->next;
    }
    dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
    
    if (ka_params && !bop->options->isPatternSearch) {
        PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, bop->outfp, FALSE);
    }
    
    if (ka_params_gap) {
        if (bop->options->isPatternSearch)
            PrintKAParametersExtra(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, bop->seedSearch->paramC, 70, bop->outfp, FALSE);
        else
            PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, bop->outfp, FALSE);
    }
    
    MemFree(ka_params);
    MemFree(ka_params_gap);

    PrintTildeSepLines(params_buffer, 70, bop->outfp);
    MemFree(params_buffer);
    free_buff();

    return TRUE;
}

Boolean PGPrintPosMatrix(CharPtr filename, posSearchItems *posSearch, 
                         compactSearchItems *compactSearch)
{
    FILE *fp;
    
    if ((fp = FileOpen(filename, "w")) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "Unable to open matrix output file %s\n", 
                  filename);
        return FALSE;
    }

    /* a diagnostic, partly an option with -Q. */
    outputPosMatrix(posSearch, compactSearch, fp); 
    FileClose(fp);

    return TRUE;
}

SeqAlignPtr PGPSeedSearch(PGPBlastOptionsPtr bop, BlastSearchBlkPtr search,
                          posSearchItems *posSearch, 
                          SeqLocPtr PNTR seqloc_duplicate,
                          SeqAlignPtr PNTR PNTR lastSeqAligns,
                          Int4Ptr numLastSeqAligns)
{
    Uint1Ptr query = NULL; /*query sequence read in*/
    Uint1Ptr unfilter_query = NULL; /*needed if seg will filter query*/
    SeqLocPtr seg_slp;  /*pointer to structure for seg filtering*/
    Int4 i, queryLength;   /*length of query sequence*/
    SeqAlignPtr head;
    SeqAnnotPtr annot;
    SeqFeatPtr sfp;
    SeqLocPtr seqloc, next;
    ValNodePtr seedReturn; /*return value from seedEngineCore, which
                             is a list of lists of SeqAligns, one
                             list of SeqAligns per pattern occurrence*/
    
    ValNodePtr info_vnp;  /* Output text messages from seedEngineCore() */
    
    query = BlastGetSequenceFromBioseq(bop->fake_bsp, &queryLength);
    seg_slp = BlastBioseqFilter(bop->fake_bsp, 
                                bop->options->filter_string);
    if (seg_slp) {
        unfilter_query = MemNew((queryLength + 1) * sizeof(Uint1));
        for (i = 0; i < queryLength; i++)
            unfilter_query[i] = query[i];
        BlastMaskTheResidues(query,queryLength,21,seg_slp,FALSE, 0);
    }
    
    search->gap_align = GapAlignBlkNew(1,1);
    search->gap_align->gap_open = bop->options->gap_open;
    search->gap_align->gap_extend = bop->options->gap_extend;

    search->gap_align->decline_align = (-(BLAST_SCORE_MIN));
    /* search->gap_align->decline_align = myargs[41].intvalue; */

    search->gap_align->x_parameter = bop->options->gap_x_dropoff
        *NCBIMATH_LN2/bop->seedSearch->paramLambda;
    search->gap_align->matrix = search->sbp->matrix;
    initProbs(bop->seedSearch);
    init_order(search->gap_align->matrix, bop->program_flag, 
               FALSE, bop->seedSearch);
    
    for(i = 0; i < queryLength; i++)
        query[i] = bop->seedSearch->order[query[i]];
    
    if (unfilter_query) {
        for(i = 0; i < queryLength; i++)
            unfilter_query[i] = bop->seedSearch->order[unfilter_query[i]];
    }
            
    seqloc = NULL;
    seedReturn = seedEngineCore(search, bop->options, query, 
                                unfilter_query,
                                bop->blast_database, bop->patfile, 
                                bop->program_flag, 
                                bop->patfp, FALSE, FALSE, 
                                bop->seedSearch, bop->options->ethresh, 
                                myargs[34].floatvalue, posSearch, 
                                &seqloc, TRUE, &info_vnp);
    
    PGPOutTextMessages(info_vnp, bop->outfp);
    ValNodeFreeData(info_vnp);
    
    BlastErrorPrint(search->error_return);
    *seqloc_duplicate = seqloc;
    head = convertValNodeListToSeqAlignList(seedReturn, lastSeqAligns, 
                                            numLastSeqAligns);

    bop->featureOrder[FEATDEF_REGION] = 1;
    bop->groupOrder[FEATDEF_REGION] = 1;
    annot = bop->fake_bsp->annot = SeqAnnotNew();
    bop->fake_bsp->annot->type = 1;	/* ftable. */
    while (seqloc) {
        next = seqloc->next;
        sfp = SeqFeatNew();
        sfp->location = seqloc;
        sfp->data.choice = SEQFEAT_REGION;
        sfp->data.value.ptrvalue = StringSave("pattern");
        annot->data = sfp;
        if (next) {
            annot->next = SeqAnnotNew();
            annot = annot->next;
            annot->type = 1;
        }
        seqloc = next;
    }

    if (query != NULL)
        MemFree(query);
    if (unfilter_query != NULL)
        unfilter_query = MemFree(unfilter_query);
    
    return head;
}

void PGPFormatMainOutput(SeqAlignPtr head, PGPBlastOptionsPtr bop,
                         BlastSearchBlkPtr search, Int4 thisPassNum,
                         SeqAlignPtr PNTR lastSeqAligns, 
                         Int4 numLastSeqAligns, SeqLocPtr seed_seqloc,
                         Int2Ptr posRepeatSequences)
{
    SeqAnnotPtr seqannot;
    ValNodePtr pruneSeed, seedReturn;  
    BlastPruneSapStructPtr prune;

    if(head == NULL) {
        fprintf(bop->outfp, "\n\n ***** No hits found ******\n\n");
        return;
    }
    
    seqannot = SeqAnnotNew();
    seqannot->type = 2;
    AddAlignInfoToSeqAnnot(seqannot, 2);
    seqannot->data = head;

    if (search->pbp->maxNumPasses != 1) {
        fprintf(bop->outfp, "\nResults from round %d\n", 
                thisPassNum);
    }
    
    ObjMgrSetHold();
    init_buff_ex(85);

    /* ------- Printing deflines for the BLAST output ------- */

    if (thisPassNum == 1) {
        search->positionBased = FALSE;
        if (!bop->options->isPatternSearch) {
            prune = BlastPruneHitsFromSeqAlign(head, 
                             bop->number_of_descriptions, NULL);
            PrintDefLinesFromSeqAlign(prune->sap, 80, bop->outfp, 
                                      bop->print_options, FIRST_PASS, NULL);
        } else {
            seedReturn = convertSeqAlignListToValNodeList(head,lastSeqAligns, 
                                                          numLastSeqAligns);
                        
            pruneSeed = SeedPruneHitsFromSeedReturn(seedReturn, 
                                      bop->number_of_descriptions);
            PrintDefLinesExtra(pruneSeed, 80, bop->outfp, bop->print_options, 
                               FIRST_PASS, NULL, seed_seqloc);
        }
    } else {
        prune = BlastPruneHitsFromSeqAlign(head, 
                              bop->number_of_descriptions, NULL);
        if (ALL_ROUNDS) {
            PrintDefLinesFromSeqAlign(prune->sap, 80, bop->outfp, 
                                      bop->print_options, 
                                      NOT_FIRST_PASS_REPEATS, 
                                      posRepeatSequences);
            PrintDefLinesFromSeqAlign(prune->sap, 80, bop->outfp, 
                                      bop->print_options, NOT_FIRST_PASS_NEW, 
                                      posRepeatSequences);
        } else
            PrintDefLinesFromSeqAlign(prune->sap, 80, bop->outfp, 
                                      bop->print_options, FIRST_PASS, NULL);
    } /*thisPassNum == 1 */

    /* ------- --------------------------------------- ------- */

    if (ALL_ROUNDS && search->posConverged) {
        fprintf(bop->outfp, "\nCONVERGED!\n");
    }

    free_buff();


    
#ifdef OLD_ALIGNMENT
        prune = BlastPruneHitsFromSeqAlign(head, bop->number_of_alignments, 
                                           prune);

        if(!DDV_DisplayBlastSAP(prune->sap, bop->outfp, 
                                FALSE, bop->align_options)) {
            fprintf(bop->outfp, 
                    "\n\n!!!\n   "
                    "    --------  Failure to print alignment...  --------"
                    "\n!!!\n\n");
            fflush(bop->outfp);
        }
#else
    if (!(bop->options->isPatternSearch)) {
        prune = BlastPruneHitsFromSeqAlign(head, bop->number_of_alignments, 
                                           prune);
        seqannot->data = prune->sap;

        if (myargs[5].intvalue != 0) {
            ShowTextAlignFromAnnot(seqannot, 60, bop->outfp, 
                                   bop->featureOrder, bop->groupOrder, 
                                   bop->align_options, NULL, 
                                   search->mask, NULL);
        } else {
            ShowTextAlignFromAnnot(seqannot, 60, bop->outfp, 
                                   bop->featureOrder, bop->groupOrder, 
                                   bop->align_options, NULL, 
                                   search->mask, FormatScoreFunc);
        }

        /* seqannot->data = head; */

    } else {
        if (bop->number_of_alignments < bop->number_of_descriptions) {
            pruneSeed = SeedPruneHitsFromSeedReturn(pruneSeed, 
                                                    bop->number_of_alignments);
        }

        if (myargs[5].intvalue != 0) {
            ShowTextAlignFromAnnotExtra(bop->query_bsp, pruneSeed, 
                                        seed_seqloc, 60, bop->outfp, 
                                        bop->featureOrder, bop->groupOrder, 
                                        bop->align_options, NULL, 
                                        search->mask, NULL);
        } else {
            ShowTextAlignFromAnnotExtra(bop->query_bsp, pruneSeed, 
                                        seed_seqloc, 60, bop->outfp, 
                                        bop->featureOrder, bop->groupOrder, 
                                        bop->align_options, NULL, 
                                        search->mask, FormatScoreFunc);
        }
    }

#endif

    if (!(bop->options->isPatternSearch)) {
        prune = BlastPruneSapStructDestruct(prune);
    }

    search->positionBased = TRUE;
    ObjMgrClearHold();
    ObjMgrFreeCache(0);

    seqannot->data = NULL;
    seqannot = SeqAnnotFree(seqannot);

    return;
}

void PGPSeqAlignOut(PGPBlastOptionsPtr bop, SeqAlignPtr head)
{
    SeqAnnotPtr seqannot;

    if (!bop->aip_out || !head)
        return;
    
    seqannot = SeqAnnotNew();
    seqannot->type = 2;
    AddAlignInfoToSeqAnnot(seqannot, 2);
    seqannot->data = head;
    SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, bop->aip_out, NULL);
    AsnIoReset(bop->aip_out);
    bop->aip_out = AsnIoClose(bop->aip_out);
    seqannot->data = NULL;
    seqannot = SeqAnnotFree(seqannot);
}

Int2 Main (void)
     
{
    PGPBlastOptionsPtr bop;
    BlastSearchBlkPtr search;
    SeqAlignPtr  head = NULL;
    SeqLocPtr seed_seqloc = NULL;
    
    /* used for psi-blast */
    Int4 thisPassNum;       
    posSearchItems *posSearch;
    compactSearchItems *compactSearch;
    Boolean  recoverCheckpoint = FALSE;
    Boolean  alreadyRecovered = FALSE;
    Boolean  freqCheckpoint = FALSE;
    Boolean  alignCheckpoint = FALSE;
    Boolean  checkReturn = FALSE;
    
    SeqAlignPtr PNTR lastSeqAligns = NULL; 
                                /*keeps track of the last SeqAlign in
                                  each list of seedReturn so that the
                                  2-level list can be converted to a 1-level
                                  list and then back to 2-level*/
    Int4 numLastSeqAligns = 0;

    /* ----- End of definitions ----- */
    
    if (! GetArgs ("blastpgp", NUMARG, myargs))
        return (1);
    ErrSetMessageLevel(SEV_WARNING);
    
    UseLocalAsnloadDataAndErrMsg ();

    if (! SeqEntryLoad())
        return 1;    
    
    bop = PGPReadBlastOptions();
    
    search = BLASTSetUpSearchWithReadDb(bop->fake_bsp, "blastp", 
                                        bop->query_bsp->length, 
                                        bop->blast_database,
                                        bop->options, NULL);
    
    if (search == NULL)
        return 1;
    
    /*AAS*/
    if ((NULL != myargs[29].strvalue) || (NULL != myargs[39].strvalue)) {
        recoverCheckpoint = TRUE;
        if (NULL != myargs[29].strvalue) {
            freqCheckpoint = TRUE;
            alignCheckpoint = FALSE;
        } else {
            freqCheckpoint = FALSE;
            alignCheckpoint = TRUE;
        }
    }
    
    if (recoverCheckpoint)
        search->positionBased = TRUE;
    else
        search->positionBased = FALSE;
    
    global_fp = bop->outfp;
    
    PGPFormatHeader(bop);

    posSearch = NULL;
    thisPassNum = 0;
    compactSearch = NULL;
    search->posConverged = FALSE;
    global_fp = bop->outfp;
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
        } else {
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
        
        /* Print out Pos matrix if necessary */
        if (myargs[38].strvalue != NULL)
            PGPrintPosMatrix(myargs[38].strvalue, posSearch, compactSearch);
    }
    
    do {  /*AAS*/
        thisPassNum++;
        if (thisPassNum > 1)
            bop->options->isPatternSearch = FALSE;

        if(head != NULL)
            SeqAlignSetFree(head);
        
#ifdef OS_UNIX
        search->thr_info->tick_callback =  tick_callback;
        fprintf(global_fp, "%s", "Searching");
        fflush(global_fp);
#endif
        if (1 == thisPassNum && (!recoverCheckpoint)) {
            
            posSearch = (posSearchItems *) 
                MemNew(1 * sizeof(posSearchItems));
        }

        /* ----- Here is real BLAST search done ------- */

        if (bop->options->isPatternSearch && 
            (1 == thisPassNum && (!recoverCheckpoint))) {
            head = PGPSeedSearch(bop, search, posSearch, 
                                 &seed_seqloc,
                                 &lastSeqAligns, &numLastSeqAligns);
        } else {
            if ((1 == thisPassNum) && (!recoverCheckpoint))
                head = BioseqBlastEngineCore(search, bop->options, NULL);
            else
                head = BioseqBlastEngineCore(search, bop->options, search->sbp->posMatrix);  
        }
        /* ---------------------------------------------- */

        if (recoverCheckpoint) {
            compactSearchDestruct(compactSearch);
            recoverCheckpoint = FALSE;
            alreadyRecovered = TRUE;
        }
        
        
        compactSearch = compactSearchNew(compactSearch);
        copySearchItems(compactSearch, search);
        
        /* The next two calls (after the "if") are diagnostics 
           for Stephen. Don't perform this if only one pass will 
           be made (i.e., standard BLAST) */
        
        if (ALL_ROUNDS && 1 != search->pbp->maxNumPasses) {
            if ((1 == thisPassNum)  && (!alreadyRecovered))
                posInitializeInformation(posSearch, search);
            posPrintInformation(posSearch, search, thisPassNum);
        }
        
#ifdef OS_UNIX
        fprintf(global_fp, "%s", "done\n\n");
#endif
        
        /* AAS */
        if (thisPassNum == 1) {
            ReadDBBioseqFetchEnable ("blastpgp", bop->blast_database, 
                                     FALSE, TRUE);
        } else {
        
            /* Have things converged? */
            if (ALL_ROUNDS && search->pbp->maxNumPasses != 1) {
                posConvergenceTest(posSearch, search, head, thisPassNum);
            }
        }

        /*AAS*/
        search->positionBased = TRUE;
        if (alreadyRecovered) {
            posCheckpointFreeMemory(posSearch, compactSearch->qlength);
            alreadyRecovered = FALSE;
        }
        
        if (ALL_ROUNDS && thisPassNum > 1) {
            posCleanup(posSearch, compactSearch);
        }
        
        if (!search->posConverged && (search->pbp->maxNumPasses == 0 || 
                            (thisPassNum < search->pbp->maxNumPasses))) {
            if (ALL_ROUNDS) {
                search->sbp->posMatrix = 
                    CposComputation(posSearch, search, compactSearch, 
                                    head, myargs[28].strvalue, 
                                    (bop->options->isPatternSearch && 
                                     (1== thisPassNum)), 
                                    &(search->error_return)); /*AAS*/
                BlastErrorPrint(search->error_return);
            } else {
                search->sbp->posMatrix = 
                    WposComputation(compactSearch, head); 
            }
#if 0
            /* DEBUG Printing of the matrix */
            {{
                Int4 i, j;
                FILE *fd;
                
                fd = FileOpen("omatrix.out", "w");
                for(i = 0; i < bop->query_bsp->length; i++) {
                    for(j = 0; j < 26; j++) {
                        fprintf(fd, "%d ", search->sbp->posMatrix[i][j]);
                    } 
                    fprintf(fd, "\n");
                } 
                FileClose(fd);
            }}
#endif
        } else {
            search->sbp->posMatrix = NULL;
        }

        /* Here is all BLAST formating of the main output done */
        PGPFormatMainOutput(head, bop, search, thisPassNum,
                            lastSeqAligns, numLastSeqAligns, 
                            seed_seqloc, posSearch->posRepeatSequences);
        
        if (ALL_ROUNDS && thisPassNum > 1) {
            MemFree(posSearch->posRepeatSequences);
        }

        if (!search->posConverged && (0 == search->pbp->maxNumPasses || 
                              thisPassNum < search->pbp->maxNumPasses)) {
            
            /* Print out pos matrix if necessary */
            if (ALL_ROUNDS && (myargs[38].strvalue != NULL))
                PGPrintPosMatrix(myargs[38].strvalue, posSearch, 
                                 compactSearch);
        }
        
    } while (( 0 == search->pbp->maxNumPasses || thisPassNum < (search->pbp->maxNumPasses)) && (! (search->posConverged)));

    if (bop->aip_out != NULL && head != NULL)
        PGPSeqAlignOut(bop, head);

    SeqAlignSetFree(head);
    
    /* Here we will print out footer of BLAST output */
    PGPFormatFooter(bop, search);
    
    ReadDBBioseqFetchDisable();
    if (bop->options->isPatternSearch) {
        bop->seedSearch = MemFree(bop->seedSearch);
        FileClose(bop->patfp);
    }
    
    if (ALL_ROUNDS)
        posSearch = MemFree(posSearch);
    
    compactSearchDestruct(compactSearch);
    bop->options = BLASTOptionDelete(bop->options);
    search = BlastSearchBlkDestruct(search);
    
    ObjMgrFreeCache(0);
    FileClose(bop->infp);
    PGPFreeBlastOptions(bop);
    
    return 0;
}
	

/* Nothing below this line is executable code */

#ifdef PRINT_ONLY_ALIGNMENT
        {{
            AsnIoPtr aip;
            
            if (seqannot)
                seqannot = SeqAnnotFree(seqannot);
            
            seqannot = SeqAnnotNew();
            seqannot->type = 2;
            AddAlignInfoToSeqAnnot(seqannot, 2);
            seqannot->data = head;
            aip = AsnIoOpen("stdout", "w");
            SeqAnnotAsnWrite(seqannot, aip, NULL);
            AsnIoReset(aip);
            AsnIoClose(aip);
            
            seqannot->data = NULL;
            seqannot = SeqAnnotFree(seqannot);
            
            return 0;
        }}
#endif
