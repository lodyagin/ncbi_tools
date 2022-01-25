/* $Id: rpsblast.c,v 6.15 2000/04/13 18:50:55 shavirin Exp $
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
* $Revision: 6.15 $
*
* File Description:
*         Main file for RPS BLAST program
*
* $Log: rpsblast.c,v $
* Revision 6.15  2000/04/13 18:50:55  shavirin
* Fixed serious memory leaks.
*
* Revision 6.14  2000/04/12 14:15:41  shavirin
* Added back ObjMgrFreeCache together wirg seqmgr changes, those removed
* deadlock.
*
* Revision 6.13  2000/03/28 20:33:44  shavirin
* Changed logic of processing MT - multiple FASTA files.
*
* Revision 6.12  2000/03/10 20:00:15  shavirin
* Added multi-thread support for multi-FASTA files.X
*
* Revision 6.11  2000/03/02 21:06:09  shavirin
* Added -U option, that allows to consider low characters in FASTA files
* as filtered regions (for blastn, blastp and tblastn).
*
* Revision 6.10  2000/02/23 21:03:36  shavirin
* Fixed -z and -Y options in rpsblast.
*
* Revision 6.9  2000/02/17 21:28:55  shavirin
* Added option is_rps_blast = TRUE.
*
* Revision 6.8  2000/02/15 16:14:04  shavirin
* Minor changes.
*
* Revision 6.7  2000/02/11 22:05:01  shavirin
* Oprion do_sum_stats set to FALSE.
*
* Revision 6.6  2000/02/11 20:51:01  shavirin
* Added possibility to search PSSM database against DNA sequences.
*
* Revision 6.5  2000/02/08 17:39:08  shavirin
* Empty log message.
*
* Revision 6.4  2000/02/01 17:22:48  shavirin
* Updated function RPSViewSeqAlign().
*
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
#include <ncbithr.h>
#include <rpsutil.h>

#if PURIFY
#include "/am/purew/solaris2/new/../purify/purify-4.5-solaris2/purify.h"
#endif

typedef struct _rps_blast_options {
    BLAST_OptionsBlkPtr options;
    CharPtr blast_database;
    SeqEntryPtr sep;
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

typedef struct _rps_thr_data {
    RPSBlastOptionsPtr rpsbop;
    RPSInfoPtr rpsinfo;
} RPSThrData, PNTR RPSThrDataPtr;

#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs [] = {
    {"Input query sequence (this parameter must be set)",  /* 0 */
     "stdin", NULL,NULL,FALSE,'i',ARG_FILE_IN, 0.0,0,NULL},
    {"RPS BLAST Database",            /* 1 */
     NULL, NULL,NULL,FALSE,'d',ARG_FILE_IN, 0.0,0,NULL},
    {"Query sequence is protein ",     /* 2 */
     "T", NULL,NULL,TRUE, 'p', ARG_BOOLEAN, 0.0,0,NULL},
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
      "10000", NULL, NULL, FALSE, 'D', ARG_INT, 0.0, 0, NULL},
    {"Logfile name ",     /* 30 */
     "rpsblast.log", NULL,NULL,TRUE,'l',ARG_FILE_OUT, 0.0,0,NULL},
    {"Use lower case filtering of FASTA sequence",    /* 31 */
     "F", NULL,NULL,TRUE,'U',ARG_BOOLEAN, 0.0,0,NULL},
};

static TNlmSemaphore MaxThreadsSem;

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

    FileClose(rpsbop->outfp);
    BLASTOptionDelete(rpsbop->options);
    RPSInfoDetach(rpsbop->rpsinfo);
    
    SeqEntryFree(rpsbop->sep);
    
    if(!rpsbop->believe_query) {
        rpsbop->fake_bsp->descr = NULL;
        rpsbop->fake_bsp->seq_data = NULL;
        BioseqFree(rpsbop->fake_bsp);
    }

    MemFree(rpsbop->rps_lookup);
    MemFree(rpsbop->rps_database);
    MemFree(rpsbop->rps_matrix);
    
    MemFree(rpsbop);
    
    return;
}

RPSBlastOptionsPtr RPSReadBlastOptions(RPSInfoPtr rpsinfo_main, 
                                       SeqEntryPtr sep, SeqLocPtr slp)
{
    RPSBlastOptionsPtr rpsbop;
    BLAST_OptionsBlkPtr options;
    Char buffer[512];
    BioseqPtr bsp;
    static Int4 count;
    
    rpsbop = MemNew(sizeof(RPSBlastOptions));
    
    rpsbop->rps_database = StringSave(myargs [1].strvalue);
    
    sprintf(buffer, "%s.rps", rpsbop->rps_database);
    rpsbop->rps_matrix = StringSave(buffer);
    
    sprintf(buffer, "%s.loo", rpsbop->rps_database);
    rpsbop->rps_lookup = StringSave(buffer);
    
    rpsbop->rpsinfo = RPSInfoAttach(rpsinfo_main);
    
    if (myargs [6].strvalue != NULL) {
        if ((rpsbop->outfp = FileOpen(myargs [6].strvalue, "a")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "rpsblast: Unable to open output "
                      "file %s\n", myargs [6].strvalue);
            return NULL;
        }
    }
    
    rpsbop->sep = sep;
    
    SeqEntryExplore(sep, &rpsbop->query_bsp, 
                    myargs[2].intvalue ? FindProt : FindNuc); 
    
    if (rpsbop->query_bsp == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
        return NULL;
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
    
    options = BLASTOptionNew(rpsbop->rpsinfo->query_is_prot ?
                             "blastp" : "tblastn", 
                             (Boolean)myargs[14].intvalue);
    rpsbop->options = options;

    rpsbop->options->query_lcase_mask = slp; /* External filtering */
    
   if (myargs[26].intvalue)
        options->db_length = (Int8) myargs[26].intvalue;
    
    if (myargs[27].floatvalue)
        options->searchsp_eff = (Nlm_FloatHi) myargs[27].floatvalue;
    
    /* Necessary options for RPS Blast */
    options->do_sum_stats = FALSE;
    options->is_rps_blast = TRUE; 
    
    rpsbop->number_of_descriptions = myargs[23].intvalue;
    rpsbop->number_of_alignments = myargs[24].intvalue;
    
    /* Update size of the database in accordance with RPS Database size */
    RPSUpdateDbSize(rpsbop->options, rpsbop->rpsinfo, rpsbop->query_bsp->length);
    
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
    options->threshold_first = (Int4) myargs [3].intvalue;
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
    
    options = BLASTOptionValidate(options, rpsbop->rpsinfo->query_is_prot ? 
                                  "blastp" : "tblastn");
    
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
        sprintf(buffer, "QUERY_%d", count++);
        obidp->str = StringSave(buffer);
        ValNodeAddPointer(&(rpsbop->fake_bsp->id), SEQID_LOCAL, obidp);
        
        /* FASTA defline not parsed, ignore the "lcl|tempseq" ID. */
        rpsbop->query_bsp->id = SeqIdSetFree(rpsbop->query_bsp->id);

        BLASTUpdateSeqIdInSeqInt(options->query_lcase_mask, 
                                 rpsbop->fake_bsp->id);
        
    }
    
    return rpsbop;
}
void RPSViewSeqAlign(SeqAlignPtr seqalign, RPSBlastOptionsPtr rpsbop, 
                     ValNodePtr mask)
{
    SeqAnnotPtr seqannot;
    AsnIoPtr aip;
    BlastPruneSapStructPtr prune;
    Uint1 align_type;

    free_buff();    
    init_buff_ex(128);

    BlastPrintReference(FALSE, 90, rpsbop->outfp);
    fprintf(rpsbop->outfp, "\n");
    AcknowledgeBlastQuery(rpsbop->fake_bsp, 70, rpsbop->outfp, FALSE, FALSE);

    if(seqalign == NULL) {
        fprintf(rpsbop->outfp, "\nNo hits found for the sequence...\n\n");
        return;
    }
    
    prune = BlastPruneHitsFromSeqAlign(seqalign, 
                                       rpsbop->number_of_descriptions, NULL);
    PrintDefLinesFromSeqAlign(prune->sap, 80, rpsbop->outfp, 
                              rpsbop->print_options, FIRST_PASS, NULL);

    seqannot = SeqAnnotNew();
    seqannot->type = 2;
    align_type = BlastGetProgramNumber(rpsbop->rpsinfo->query_is_prot ?
                                       "blastp" : "blastx");
    
    AddAlignInfoToSeqAnnot(seqannot, align_type); /* blastp or tblastn */
    
    if(!rpsbop->rpsinfo->query_is_prot)
        rpsbop->align_options += TXALIGN_BLASTX_SPECIAL;
    
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

    seqannot->data = NULL;
    seqannot = SeqAnnotFree(seqannot);

    free_buff();
    return;
}
static BioseqPtr createFakeProtein(void)
{
    BioseqPtr bsp;
    CharPtr sequence = "THEFAKEPROTEIN"; 
    
    bsp = BioseqNew();
    
    bsp->mol = Seq_mol_aa;
    bsp->seq_data_type = Seq_code_iupacaa;
    bsp->repr = Seq_repr_raw;
    bsp->length = StringLen(sequence);
    
    bsp->seq_data = BSNew(64);
    BSWrite(bsp->seq_data, sequence, StringLen(sequence));
    
    bsp->id = MakeNewProteinSeqIdEx(NULL, NULL, "ssh_seq", NULL);

    return bsp;
}

Boolean RPSFormatFooter(RPSBlastOptionsPtr rpsbop, BlastSearchBlkPtr search)
{
    ValNodePtr  mask_loc, mask_loc_start, vnp;
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
    
    free_buff();    
    init_buff_ex(85);
    dbinfo_head = dbinfo;
    while (dbinfo) {
        PrintDbReport(dbinfo, 70, rpsbop->outfp);
        dbinfo = dbinfo->next;
    }
    dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
    
    if (ka_params) {
        PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 
                          70, rpsbop->outfp, FALSE);
    }
    
    if (ka_params_gap) {
        PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K,
                          ka_params_gap->H, 70, rpsbop->outfp, TRUE);
    }
    
    MemFree(ka_params);
    MemFree(ka_params_gap);
    
    PrintTildeSepLines(params_buffer, 70, rpsbop->outfp);
    MemFree(params_buffer);
    free_buff();


    mask_loc_start = mask_loc;
    while (mask_loc) {
        SeqLocSetFree(mask_loc->data.ptrvalue);
        mask_loc = mask_loc->next;
    }
    ValNodeFree(mask_loc_start);
    search->mask = NULL;
    
    other_returns = ValNodeFree(other_returns);

    fflush(rpsbop->outfp);

    /* ----- */    
    return TRUE;
}

SeqEntryPtr RPSGetNextSeqEntry(SeqLocPtr PNTR slp)
{
    SeqEntryPtr sep;
    static TNlmMutex read_mutex;
    static FILE *infp;
    static Boolean end_of_data = FALSE;

    NlmMutexInit(&read_mutex);
    NlmMutexLock(read_mutex);
    
    if(end_of_data) {
        NlmMutexUnlock(read_mutex);
        return NULL;
    }

    /* Opening file with input sequences */
    if(infp == NULL) {
        if ((infp = FileOpen(myargs [0].strvalue, "r")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, 
                      "rpsblast: Unable to open input file %s\n", 
                      myargs [0].strvalue);
            end_of_data = TRUE;
            NlmMutexUnlock(read_mutex);
        }
    }
    
    if(myargs[31].intvalue) {
        sep = FastaToSeqEntryForDb (infp, !myargs[2].intvalue, NULL, myargs[19].intvalue, NULL, NULL, slp);
    } else {
        sep = FastaToSeqEntryEx (infp, !myargs[2].intvalue, NULL, myargs[19].intvalue);
    }
    
    if(sep == NULL) {            /* Probably last FASTA entry */
        end_of_data = TRUE;
        FileClose(infp);
        NlmMutexUnlock(read_mutex);
        return NULL;
    }
    
    NlmMutexUnlock(read_mutex);
    
    return sep;
}

void RPSThreadOnExit(VoidPtr data)
{
    /* NlmSemaPost(MaxThreadsSem); */
    return;
}

VoidPtr RPSEngineThread(VoidPtr data)
{
    SeqAlignPtr seqalign;
    ValNodePtr other_returns, error_returns;
    BlastSearchBlkPtr search;
    RPSBlastOptionsPtr rpsbop;
    BioseqPtr bsp;
    SeqEntryPtr sep;
    Char buffer[64];
    static TNlmMutex print_mutex;
    SeqLocPtr slp = NULL;
    RPSInfoPtr rpsinfo_main;

    if((rpsinfo_main = (RPSInfoPtr) data) == NULL) {
        return NULL;
    }
    
    NlmThreadAddOnExit(RPSThreadOnExit, NULL);
    
    ReadDBBioseqFetchEnable("rpsblast",  myargs [1].strvalue, FALSE, TRUE);
    
    NlmMutexInit(&print_mutex);

    while(TRUE) {               /* Main loop */

        if((sep = RPSGetNextSeqEntry(&slp)) == NULL)
            break;
        
        if((rpsbop = RPSReadBlastOptions(rpsinfo_main, 
                                         sep, slp)) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "Unable to initialize RPS Blast.");
            break;
        }
        
        if(rpsbop->rpsinfo->query_is_prot)
            bsp = rpsbop->fake_bsp;
        else
            bsp = createFakeProtein();
        
        search = BLASTSetUpSearch (bsp, rpsbop->rpsinfo->query_is_prot ? 
                                   "blastp" : "tblastn", 
                                   bsp->length, 0, 
                                   NULL, rpsbop->options, NULL);

        seqalign = RPSBlastSearch(search, rpsbop->fake_bsp, rpsbop->rpsinfo);

        NlmMutexLock(print_mutex);

        RPSViewSeqAlign(seqalign, rpsbop, search->mask);    

        RPSFormatFooter(rpsbop, search);

        ObjMgrFreeCache(0);

        /* Final cleanup */
        
        if(!rpsbop->rpsinfo->query_is_prot) 
            BioseqFree(bsp);
        
        SeqAlignSetFree(seqalign);
        search = BlastSearchBlkDestruct(search);

        RPSBlastOptionsFree(rpsbop);

        NlmMutexUnlock(print_mutex); 
    }

    return NULL;
}

Int2 Main(void)
{
    FILE *fd;
    Char rps_matrix[128], rps_lookup[128];
    RPSInfoPtr rpsinfo_main;
    Int4 i;

    if (!GetArgs("rpsblast", NUMARG, myargs))
	return 1;
    
    if ( !ErrSetLog (myargs[30].strvalue) ) { /* Logfile */
        ErrShow();
        return 1;
    } else {
        ErrSetOpts (ERR_CONTINUE, ERR_LOG_ON);
    }
    
    UseLocalAsnloadDataAndErrMsg ();
    
    if (!SeqEntryLoad())
        return 1;
        
    if((MaxThreadsSem = NlmSemaInit(myargs[17].intvalue)) == NULL) {
        ErrPostEx(SEV_ERROR, 0,0, "NlmSemaInit failed with errno %d", errno);
        return 1;
    }

    /* Truncate output file */
    if((fd = FileOpen(myargs [6].strvalue, "w")) != NULL)
        FileClose(fd);
    
    /* Initializing RPS Database */

    sprintf(rps_matrix, "%s.rps", myargs [1].strvalue);
    sprintf(rps_lookup, "%s.loo", myargs [1].strvalue);
    
    /* Initializing RPS Blast database */
    
    /* sprintf(buffer, "%s.mat", myargs[1].strvalue); */
    /*  myargs[2].intvalue == query is protein, default = 1 */
    if((rpsinfo_main = RPSInit(myargs [1].strvalue, 
                               rps_matrix, rps_lookup,
                               myargs[2].intvalue)) == NULL) {
        ErrPostEx(SEV_ERROR, 0,0, 
                  "Failure to initialize RPS Blast database");
        return NULL;
    }

    if(myargs[17].intvalue > 1) {
        for(i = 0; i < myargs[17].intvalue; i++) {
            NlmThreadCreate(RPSEngineThread, (VoidPtr) rpsinfo_main);
        }
    } else {
        RPSEngineThread((VoidPtr) rpsinfo_main);
    }
    
    NlmThreadJoinAll();
    RPSClose(rpsinfo_main);
    
    return 0;
}
