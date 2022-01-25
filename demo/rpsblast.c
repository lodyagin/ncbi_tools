/* $Id: rpsblast.c,v 6.3 2000/01/07 22:34:05 shavirin Exp $
* ===========================================================================
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
* ===========================================================================
*
* File Name:  $RCSfile: rpsblast.c,v $
*
* Author:  Sergei Shavirin
*
* Initial Version Creation Date: 12/14/1999
*
* $Revision: 6.3 $
*
* File Description:
*         Main file for RPS BLAST program
*
* $Log: rpsblast.c,v $
* Revision 6.3  2000/01/07 22:34:05  shavirin
* Added printing of SeqAlignment if necessary.
*
* Revision 6.2  1999/12/30 18:37:53  shavirin
* Added NCBI header and Log information.
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <lookup.h>
#include <tofasta.h>
#include <txalign.h>

#include <rpsutil.h>

typedef struct _rps_blast_options {
    BLAST_OptionsBlkPtr options;
    CharPtr blast_database;
    BioseqPtr query_bsp, fake_bsp;
    Int4 number_of_descriptions, number_of_alignments;
    FILE *outfp;
    AsnIoPtr aip_out;
    Boolean html;
    Boolean believe_query;
    Uint4 align_options, print_options;
    /* RPS Blast variables */
    CharPtr rps_database;
    CharPtr rps_matrix;
    CharPtr rps_lookup;
    RPSInfoPtr rpsinfo;
} RPSBlastOptions, PNTR RPSBlastOptionsPtr;

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs [] = {
    {"Input query sequence (this parameter must be set)",  /* 0 */
     "stdin", NULL,NULL,FALSE,'i',ARG_FILE_IN, 0.0,0,NULL},
    {"RPS BLAST Database",            /* 1 */
     NULL, NULL,NULL,FALSE,'d',ARG_FILE_IN, 0.0,0,NULL},
    {"Logfile name ",     /* 2 */
     "rpsblast.log", NULL,NULL,TRUE,'l',ARG_FILE_OUT, 0.0,0,NULL},
    { "Threshold for extending hits", /* 3 */
      "11", NULL, NULL, FALSE, 'f', ARG_INT, 0.0, 0, NULL},
    { "Expectation value (E)",        /* 4 */
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
    { "Believe the query defline", /* 19 */
      "F", NULL, NULL, FALSE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
    { "X dropoff value for final gapped alignment (in bits)", /* 20 */
      "25", NULL, NULL, FALSE, 'Z', ARG_INT, 0.0, 0, NULL},
    { "SeqAlign file ('Believe the query defline' must be TRUE)", /*21*/
      NULL, NULL, NULL, TRUE, 'O', ARG_FILE_OUT, 0.0, 0, NULL},
    { "Matrix",                 /* 22 */
      "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
    { "Number of one-line descriptions (V)", /* 23 */
      "500", NULL, NULL, FALSE, 'v', ARG_INT, 0.0, 0, NULL},
    { "Number of alignments to show (B)", /* 24 */
      "250", NULL, NULL, FALSE, 'b', ARG_INT, 0.0, 0, NULL},
    { "Word size, default if zero", /* 25 */
      "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},
    { "Effective length of the database (use zero for the real size)", /* 26 */
      "0", NULL, NULL, FALSE, 'z', ARG_INT, 0.0, 0, NULL},
    { "Effective length of the search space (use zero for the real size)", /* 27 */
      "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL},
    { "Produce HTML output",  /* 28 */
      "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Cost to decline alignment", /* 29 */
      "10000", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL}  
};

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
    
    if (myargs[28].intvalue) {
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

void RPSBlastOptionsFree(RPSBlastOptionsPtr rpsbop)
{

    BLASTOptionDelete(rpsbop->options);
    RPSClose(rpsbop->rpsinfo);
    
    /*    BioseqFree(rpsbop->query_bsp);
          BioseqFree(rpsbop->fake_bsp); */
    
    MemFree(rpsbop->rps_database);
    MemFree(rpsbop->rps_matrix);
    
    MemFree(rpsbop);
    
    return;
}

RPSBlastOptionsPtr RPSReadBlastOptions(void)
{
    RPSBlastOptionsPtr rpsbop;
    FILE *infp;
    BLAST_OptionsBlkPtr options;
    SeqEntryPtr sep;
    Char buffer[512];
    rpsbop = MemNew(sizeof(RPSBlastOptions));

    rpsbop->rps_database = StringSave(myargs [1].strvalue);

    sprintf(buffer, "%s.rps", rpsbop->rps_database);
    rpsbop->rps_matrix = StringSave(buffer);
    
    sprintf(buffer, "%s.loo", rpsbop->rps_database);
    rpsbop->rps_lookup = StringSave(buffer);
    
    /* Initializing RPS Blast database */

    /* sprintf(buffer, "%s.mat", myargs[1].strvalue); */
    if((rpsbop->rpsinfo = RPSInit(rpsbop->rps_database, 
                                  rpsbop->rps_matrix,
                                  rpsbop->rps_lookup)) == NULL) {
        ErrPostEx(SEV_ERROR, 0,0, "Failure to initialize RPS Blast database");
        return NULL;
    }
    
    /* Reading query sequence */
    if ((infp = FileOpen(myargs [0].strvalue, "r")) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "rpsblast: Unable to open input file %s\n", 
                  myargs [0].strvalue);
        return NULL;
    }

    if (myargs [6].strvalue != NULL) {
        if ((rpsbop->outfp = FileOpen(myargs [6].strvalue, "w")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "rpsblast: Unable to open output "
                      "file %s\n", myargs [6].strvalue);
            return NULL;
        }
    }

    if (myargs[19].intvalue != 0)
        rpsbop->believe_query = TRUE;
    
    if (myargs[21].strvalue != NULL) {
#if 0        
        if (rpsbop->believe_query == FALSE) {
            ErrPostEx(SEV_FATAL, 0, 0, 
                      "-J option must be TRUE to use this option");
            return NULL;
        }
#endif        
        if ((rpsbop->aip_out = AsnIoOpen (myargs[21].strvalue,"w")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output "
                      "file %s\n", myargs[21].strvalue);
            return NULL;
        }
    }

    if((sep = FastaToSeqEntryEx(infp, FALSE, NULL, 
                                rpsbop->believe_query)) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "Unable to read input FASTA file\n");
        return NULL;
    }

    FileClose(infp);
    SeqEntryExplore(sep, &rpsbop->query_bsp, FindProt);    
    sep->data.ptrvalue = NULL;
    SeqEntryFree(sep);
    
    if (rpsbop->query_bsp == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
        return NULL;
    }

    rpsbop->number_of_descriptions = myargs[23].intvalue;
    rpsbop->number_of_alignments = myargs[24].intvalue;
    
    options = BLASTOptionNew("blastp", (Boolean)myargs[14].intvalue);
    rpsbop->options = options;

    /* Update size of the database in accordance with RPS Database size */
    RPSUpdateDbSize(rpsbop->options, rpsbop->rpsinfo);

    /* Set default gap params for matrix. */
    BLASTOptionSetGapParams(options, myargs[22].strvalue, 0, 0);
    
    PGPGetPrintOptions(options->gapped_calculation, &rpsbop->align_options, 
                       &rpsbop->print_options);
    
    /* decrement by one to agree with program values. */
    options->required_start = myargs[15].intvalue - 1;
    options->required_end = myargs[16].intvalue;
    if (options->required_end != -1) {
        options->required_end--;
    }

    /* This is not used */
    options->threshold_second = (Int4) myargs [3].intvalue;
    
    options->dropoff_2nd_pass  = myargs [7].floatvalue;
    options->expect_value  = (Nlm_FloatHi) myargs [4].floatvalue;
    options->hitlist_size = MAX(rpsbop->number_of_descriptions, 
                                rpsbop->number_of_alignments);
    
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

        /* options->decline_align = INT2_MAX; */
        options->decline_align = myargs[29].intvalue;

        options->gap_x_dropoff = myargs[12].intvalue;
        options->gap_x_dropoff_final = myargs[20].intvalue;
        options->gap_trigger = myargs[13].floatvalue;
    }
    
    if (StringICmp(myargs[9].strvalue, "T") == 0) {
        options->filter_string = StringSave("S");
    } else {
        options->filter_string = StringSave(myargs[9].strvalue);
    }

    /* Only one CPU may be used at this time */    
    options->number_of_cpus = (Int2) myargs[17].intvalue;
    
    options->isPatternSearch = FALSE;

    if (myargs[25].intvalue)
        options->wordsize = myargs[25].intvalue;
    if (myargs[26].intvalue)
        options->db_length = (Int8) myargs[26].intvalue;
    
    if (myargs[27].floatvalue)
        options->searchsp_eff = (Nlm_FloatHi) myargs[27].floatvalue;
    
    options = BLASTOptionValidate(options, "blastp");
    
    if (options == NULL)
        return NULL;

    if (options == NULL)
        return NULL;

    if (rpsbop->believe_query == TRUE) {
        rpsbop->fake_bsp = rpsbop->query_bsp;
    } else {
        ObjectIdPtr obidp;
        rpsbop->fake_bsp = BioseqNew();
        rpsbop->fake_bsp->descr = rpsbop->query_bsp->descr;
        rpsbop->fake_bsp->repr = rpsbop->query_bsp->repr;
        rpsbop->fake_bsp->mol = rpsbop->query_bsp->mol;
        rpsbop->fake_bsp->length = rpsbop->query_bsp->length;
        rpsbop->fake_bsp->seq_data_type = rpsbop->query_bsp->seq_data_type;
        rpsbop->fake_bsp->seq_data = rpsbop->query_bsp->seq_data;
        
        obidp = ObjectIdNew();
        obidp->str = StringSave("QUERY");
        ValNodeAddPointer(&(rpsbop->fake_bsp->id), SEQID_LOCAL, obidp);
        
        /* FASTA defline not parsed, ignore the "lcl|tempseq" ID. */
        rpsbop->query_bsp->id = SeqIdSetFree(rpsbop->query_bsp->id);
    }
    
    return rpsbop;
}
void RPSViewSeqAlign(SeqAlignPtr seqalign, RPSBlastOptionsPtr rpsbop, 
                     ValNodePtr mask)
{
    SeqAnnotPtr seqannot;
    AsnIoPtr aip;
    BlastPruneSapStructPtr prune;
    
    ObjMgrSetHold();
    init_buff_ex(128);

    BlastPrintReference(FALSE, 90, stdout);
    fprintf(stdout, "\n");
    AcknowledgeBlastQuery(rpsbop->fake_bsp, 70, stdout, FALSE, FALSE);
    
    prune = BlastPruneHitsFromSeqAlign(seqalign, 
                                       rpsbop->number_of_descriptions, NULL);
    PrintDefLinesFromSeqAlign(prune->sap, 80, rpsbop->outfp, 
                              rpsbop->print_options, FIRST_PASS, NULL);
    
    seqannot = SeqAnnotNew();
    seqannot->type = 2;
    AddAlignInfoToSeqAnnot(seqannot, 2); /* blastp */
    
    prune = BlastPruneHitsFromSeqAlign(seqalign, rpsbop->number_of_alignments, 
                                       prune);
    seqannot->data = prune->sap;

    if(rpsbop->aip_out != NULL) {     
        SeqAnnotAsnWrite(seqannot, rpsbop->aip_out, NULL);
        AsnIoClose(rpsbop->aip_out); 
    }

    if (myargs[5].intvalue != 0) {
        ShowTextAlignFromAnnot(seqannot, 60, rpsbop->outfp, 
                               NULL, NULL, rpsbop->align_options, NULL, 
                               mask, NULL);
    } else {
        ShowTextAlignFromAnnot(seqannot, 60, rpsbop->outfp, 
                               NULL, NULL, rpsbop->align_options, NULL, 
                               mask, FormatScoreFunc);
    }
    
    prune = BlastPruneSapStructDestruct(prune);
    ObjMgrClearHold();
    ObjMgrFreeCache(0);

    seqannot->data = NULL;
    seqannot = SeqAnnotFree(seqannot);

    free_buff();
    return;
}
Int2 Main(void)
{
    SeqAlignPtr seqalign;
    ValNodePtr other_returns, error_returns;
    BlastSearchBlkPtr search;
    RPSBlastOptionsPtr rpsbop;
    
    if (!GetArgs("rpsblast", NUMARG, myargs))
	return 1;
    
    if ( !ErrSetLog (myargs[2].strvalue) ) { /* Logfile */
        ErrShow();
        return 1;
    } else {
        ErrSetOpts (ERR_CONTINUE, ERR_LOG_ON);
    }

    UseLocalAsnloadDataAndErrMsg ();
    
    if (!SeqEntryLoad())
        return 1;
    
    if((rpsbop = RPSReadBlastOptions()) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "Unable to initialize RPS Blast.");
        return 1;
    }
    
    search = BLASTSetUpSearch (rpsbop->fake_bsp, 
                               "blastp", rpsbop->fake_bsp->length, 0, 
                               NULL/*all_words*/, rpsbop->options, NULL);
    
    seqalign = RPSBlastSearch(search, rpsbop->fake_bsp, rpsbop->rpsinfo);
    
    RPSViewSeqAlign(seqalign, rpsbop, search->mask);
    
    /* Final cleanup */
    
    SeqAlignSetFree(seqalign);
    search = BlastSearchBlkDestruct(search);
    
    RPSBlastOptionsFree(rpsbop);
    
    return 0;
}
