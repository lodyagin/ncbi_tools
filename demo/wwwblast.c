/* $Id: wwwblast.c,v 6.5 2000/05/17 15:50:51 shavirin Exp $
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
* File Name:  $RCSfile: wwwblast.c,v $
*
* Author:  Sergei Shavirin
*
* Initial Creation Date: 03/15/2000
*
* $Revision: 6.5 $
*
* File Description:
*        Standalone WWW Blast CGI program.
*
* $Log: wwwblast.c,v $
* Revision 6.5  2000/05/17 15:50:51  shavirin
* Moved many functions to the wwwbutl.c file.
*
* Revision 6.4  2000/04/21 18:10:59  shavirin
* Added possibility to print Patrick's alignment.
*
* Revision 6.3  2000/03/28 14:44:20  shavirin
* Changed function ctime_r to ctime() for compatibility.
*
* Revision 6.2  2000/03/24 16:05:37  shavirin
* Added option to be used as NCBI client/server.
*
* Revision 6.1  2000/03/20 19:01:00  shavirin
* Initial revision.
*
*
* ==========================================================================
*/

#include <wwwblast.h>

static int LIBCALLBACK
tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{

    /* do nothing */
    return 0;
}

static Int4 get_number_alignment(SeqAlignPtr align)
{
    Int4 num = 0;

    while(align)
    {
	++num;
	align = align->next;
    }

    return num;
}
static void
PrintMotd(CharPtr string, FILE *fp, Boolean html_format)

{
    Char buffer[100];
    CharPtr ptr;
    
    if (string == NULL)
        return;
    
    buffer[0] = NULLB;
    ptr = buffer;
    
    if (html_format) {
        fprintf(fp, "<PRE>\n");
    }
    
    while (*string != NULLB) {
        if (*string == '~') {
            *ptr = NULLB;
            fprintf(fp, "%s\n", buffer);
            buffer[0] = NULLB;
            ptr = buffer;
            string++;
            if (*string == NULLB)
                break;
        } else {
            *ptr=*string;
            ptr++;  string++;
        }
    }
    *ptr = NULLB;
    fprintf(fp, "%s\n", buffer);
    
    if (html_format) {
        fprintf(fp, "</PRE>\n");
    }
    
    fflush(fp);

    return;
}
#ifdef NCBI_CLIENT_SERVER
static Boolean
TraditionalBlastReportWithImage(BioseqPtr bsp, BLAST_OptionsBlkPtr options, BlastNet3Hptr bl3hp, CharPtr program, CharPtr database, Uint4 print_options, Uint4 align_options, Int4 number_of_descriptions, Int4 number_of_alignments, Boolean show_overview, CharPtr blast_type, Int4 color_schema)

{
    BlastDbinfoPtr dbinfo;
    BlastKABlkPtr ka_params=NULL, ka_params_gap=NULL;
    BlastPruneSapStructPtr prune;
    BLAST_MatrixPtr matrix;
    Boolean query_is_na, db_is_na;
    Boolean status;
    CharPtr params_buffer=NULL;
    Int4 number_of_hits_private=0, length;
    SeqAlignPtr seqalign;
    SeqAnnotPtr seqannot=NULL;
    TxDfDbInfoPtr tx_dbinfo=NULL, tx_dbinfo_head;
    ValNodePtr mask_loc, mask_loc_start, other_returns, error_returns, vnp, vnp1=NULL;
    Uint1 align_type;
    Uint1 f_order[FEATDEF_ANY], g_order[FEATDEF_ANY];
    
    MemSet((Pointer)(g_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    MemSet((Pointer)(f_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    
    if (bsp == NULL)
        return FALSE;
    
    if (bl3hp == NULL || program == NULL || database == NULL)
        return FALSE;
    
    align_type = BlastGetTypes(program, &query_is_na, &db_is_na);

    init_buff_ex(85);
    dbinfo = BlastRequestDbInfo(bl3hp, database, !db_is_na);

    if (dbinfo)
        PrintDbInformationBasic(database, !db_is_na, 70, dbinfo->definition, dbinfo->number_seqs, dbinfo->total_length, stdout, TRUE);
    dbinfo = BlastDbinfoFree(dbinfo);
    free_buff();

    if (bsp) {
        seqalign = BlastBioseqNetCore(bl3hp, bsp, program, database, options, &other_returns, &error_returns, NULL, NULL, &status);
    }
    
    BlastErrorPrintExtra(error_returns, TRUE, stdout);
    
    mask_loc = NULL;
    for (vnp=other_returns; vnp; vnp = vnp->next) {
        switch (vnp->choice) {
        case TXDBINFO:
            tx_dbinfo = (TxDfDbInfoPtr) vnp->data.ptrvalue;
            break;
        case TXKABLK_NOGAP:
            ka_params = (BlastKABlkPtr) vnp->data.ptrvalue;
            break;
        case TXKABLK_GAP:
            ka_params_gap = (BlastKABlkPtr) vnp->data.ptrvalue;
            break;
        case TXPARAMETERS:
            params_buffer = (CharPtr) vnp->data.ptrvalue;
            break;
        case TXMATRIX:
            /* This function should not use matrix */
            matrix = (BLAST_MatrixPtr) vnp->data.ptrvalue;
            BLAST_MatrixDestruct(matrix);
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
    
    if (seqalign) {
        seqannot = SeqAnnotNew();
        seqannot->type = 2;
        AddAlignInfoToSeqAnnot(seqannot, align_type);
        seqannot->data = seqalign;
        
        if(show_overview) {
            Char f_name[64], title[1024], href[64];
            Int4 align_num;       

            sprintf(f_name, "%ld%ld.gif", (long)random(), (long)getpid());
            sprintf(href, "nph-viewgif.cgi?");
        
            align_num = get_number_alignment(seqalign); 
            sprintf(title, 
                    "<H3><a href=\"docs/newoptions.html#graphical-overview\"> "
                    "Distribution of %ld Blast Hits on the Query Sequence</a> "
                    "</H3>\n", (long)align_num);  
            
            
            /* Open HTML form */
            fprintf(stdout, "<FORM NAME=\"BLASTFORM\">\n");
            fflush(stdout);
            
            PrintAlignmentOverview(seqannot, stdout, 
                                   "BLASTFORM", href, f_name, title); 
        }

        prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_descriptions, NULL);
        ObjMgrSetHold();
        init_buff_ex(85);

        PrintDefLinesFromSeqAlignEx2(prune->sap, 80, stdout, print_options, FIRST_PASS, NULL, number_of_descriptions, database, blast_type);
        free_buff();
        
        prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_alignments, prune);
        seqannot->data = prune->sap;

        if(theInfo->color_schema != 0 && 
           (!StringICmp(program, "blastn") || 
            !StringICmp(program, "blastp"))) {
            
            if(!DDV_DisplayBlastPairList(prune->sap, mask_loc, stdout, 
                                         query_is_na, align_options, 
                                         color_schema)) { 
                fprintf(stdout, 
                        "\n\n!!!\n   "
                        "    --------  Failure to print alignment...  --------"
                        "\n!!!\n\n");
                fflush(stdout);
            }
        } else {
            
            if (align_options & TXALIGN_MASTER)
                ShowTextAlignFromAnnot(seqannot, 60, stdout, f_order, g_order, align_options, NULL, mask_loc, NULL);
            else
                ShowTextAlignFromAnnot(seqannot, 60, stdout, f_order, g_order, align_options, NULL, mask_loc, FormatScoreFunc);
        }
        
        printf("<PRE>\n");
        
        seqannot->data = seqalign;
        number_of_hits_private = prune->original_number; 
        prune = BlastPruneSapStructDestruct(prune);
        ObjMgrClearHold();
    }
    
    mask_loc_start = mask_loc;
    while (mask_loc) {
        SeqLocSetFree(mask_loc->data.ptrvalue);
        mask_loc = mask_loc->next;
    }
    ValNodeFree(mask_loc_start);
    
    init_buff_ex(85);
    tx_dbinfo_head = tx_dbinfo;
    while (tx_dbinfo) {
        PrintDbReport(tx_dbinfo, 70, stdout);
        tx_dbinfo = tx_dbinfo->next;
    }
    tx_dbinfo_head = TxDfDbInfoDestruct(tx_dbinfo_head);
    
    if (ka_params) {
        PrintKAParameters(ka_params->lambda, ka_params->k, ka_params->h, 
                          70, stdout, FALSE);
        MemFree(ka_params);
    }
    
    if (ka_params_gap) {
        PrintKAParameters(ka_params_gap->lambda, ka_params_gap->k, ka_params_gap->h, 70, stdout, TRUE);
        MemFree(ka_params_gap);
    }
    
    PrintTildeSepLines(params_buffer, 70, stdout);
    MemFree(params_buffer);
    free_buff();
    
    other_returns = ValNodeFree(other_returns);
    
    if (seqannot)
        seqannot = SeqAnnotFree(seqannot);
    
    return status;
}

Boolean WWWBlastDoClientSearch(WWWBlastInfoPtr theInfo)
{
    BlastNet3Hptr	bl3hp;
    BlastResponsePtr	response;
    BlastVersionPtr	blast_version;
    CharPtr		date, motd, version;
    Boolean status;
    
    if(theInfo == NULL)
	return FALSE;
        
    if (!BlastInit("blastcl3", &bl3hp, &response)) {
        WWWBlastErrMessage(BLASTErrClient, NULL);                
        return FALSE;
    }
    
    if (response && response->choice == BlastResponse_init) {
        blast_version = response->data.ptrvalue;
        version = blast_version->version;
        date = blast_version->date;
    } else {
        WWWBlastErrMessage(BLASTErrClient, NULL);                
        return FALSE;
    }
    
    BlastNetBioseqFetchEnable(bl3hp, theInfo->database, 
                              theInfo->db_is_na, TRUE);
    
#ifdef BLAST_PRINT_MOTD
    motd = Blast3GetMotd(bl3hp);
    PrintMotd(motd, stdout, TRUE);
    motd = MemFree(motd);
#endif
    
    status = TraditionalBlastReportWithImage(theInfo->fake_bsp, 
                                             theInfo->options, 
                                             bl3hp, theInfo->program, 
                                             theInfo->database,
                                             theInfo->print_options, 
                                             theInfo->align_options, 
                                             theInfo->number_of_descriptions, 
                                             theInfo->number_of_alignments, 
                                             theInfo->show_overview,
                                             theInfo->blast_type,
                                             theInfo->color_schema);
    if (status == FALSE) {
        WWWBlastErrMessage(BLASTErrServer, NULL);
        return FALSE;
    }
    
    return TRUE;
}
#endif

Boolean WWWBlastDoSearch(WWWBlastInfoPtr theInfo)
{
    SeqAlignPtr  seqalign;
    ValNodePtr  mask_loc, mask_loc_start, vnp, other_returns, error_returns;
    TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
    BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
    BLAST_MatrixPtr matrix;
    CharPtr params_buffer=NULL;
    SeqAnnotPtr seqannot;
    BlastPruneSapStructPtr prune;
    
    if(theInfo == NULL)
	return FALSE;
    
    PrintDbInformation(theInfo->database, !theInfo->db_is_na, 
                       70, stdout, TRUE);

    other_returns = NULL;
    error_returns = NULL;
    
    seqalign = BioseqBlastEngine(theInfo->fake_bsp, theInfo->program, 
                                 theInfo->database, theInfo->options, 
                                 &other_returns, &error_returns, 
                                 tick_callback);
    
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
    
    
    fflush(stdout);
    
    ReadDBBioseqFetchEnable ("blastall", theInfo->database, 
                             theInfo->db_is_na, TRUE);
    
    if (seqalign) {
        seqannot = SeqAnnotNew();
        seqannot->type = 2;
        AddAlignInfoToSeqAnnot(seqannot, theInfo->align_type);
        seqannot->data = seqalign;

#if !defined(NCBI_CLIENT_SERVER) && defined (NCBI_ENTREZ_CLIENT)
        if(theInfo->show_tax_blast) {
            fprintf(stdout, "<a href=#taxblast>Taxonomy reports</a><BR>");
        }
#endif
        
        /* Now printing nice gif with alignment overview */

        if(theInfo->show_overview) {
            Char f_name[64], title[1024], href[64];
            Int4 align_num;
            
            sprintf(f_name, "%ld%ld.gif", (long)random(), (long)getpid());
            sprintf(href, "nph-viewgif.cgi?");
            
            align_num = get_number_alignment(seqalign); 
            sprintf(title, 
                    "<H3><a href=\"docs/newoptions.html#graphical-overview\"> "
                    "Distribution of %ld Blast Hits on the Query Sequence</a> "
                    "</H3>\n", (long)align_num);  
	

            /* Open HTML form */
            fprintf(stdout, "<FORM NAME=\"BLASTFORM\">\n");
            fflush(stdout);
            
            PrintAlignmentOverview(seqannot, stdout, 
                                   "BLASTFORM", href, f_name, title); 
        }

        prune = BlastPruneHitsFromSeqAlign(seqalign, theInfo->number_of_descriptions, NULL);
        ObjMgrSetHold();
        init_buff_ex(85);
        PrintDefLinesFromSeqAlign(prune->sap, 80, stdout, 
                                  theInfo->print_options, FIRST_PASS, NULL);
        free_buff();
        
        prune = BlastPruneHitsFromSeqAlign(seqalign, theInfo->number_of_alignments, prune);
        seqannot->data = prune->sap;
        
        if(theInfo->color_schema != 0 && 
           (!StringICmp(theInfo->program, "blastn") || 
            !StringICmp(theInfo->program, "blastp"))) {
            
            if(!DDV_DisplayBlastPairList(prune->sap, mask_loc, stdout, 
                                         theInfo->query_is_na, 
                                         theInfo->align_options, 
                                         theInfo->color_schema)) { 
                fprintf(stdout, 
                        "\n\n!!!\n   "
                        "    --------  Failure to print alignment...  --------"
                        "\n!!!\n\n");
                fflush(stdout);
            }
            printf("<PRE>\n");
        } else {
            
            if (theInfo->align_view != 0)
                ShowTextAlignFromAnnot(seqannot, 60, stdout, NULL, NULL, 
                                       theInfo->align_options, NULL, 
                                       mask_loc, NULL);
            else
                ShowTextAlignFromAnnot(seqannot, 60, stdout, NULL, NULL, 
                                       theInfo->align_options, NULL, mask_loc, 
                                       FormatScoreFunc);
        }
        
        seqannot->data = seqalign;
        prune = BlastPruneSapStructDestruct(prune);
        ObjMgrClearHold();
        ObjMgrFreeCache(0);
        
    } else {
        fprintf(stdout, "\n\n ***** No hits found ******\n\n");
    }
    
    matrix = BLAST_MatrixDestruct(matrix);
    
    fprintf(stdout, "<PRE>");
    init_buff_ex(85);
    dbinfo_head = dbinfo;
    while (dbinfo) {
        PrintDbReport(dbinfo, 70, stdout);
        dbinfo = dbinfo->next;
    }
    dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
    
    if (ka_params) {
        PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, stdout, FALSE);
        MemFree(ka_params);
    }
    
    if (ka_params_gap) {
        PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, stdout, TRUE);
        MemFree(ka_params_gap);
    }
    
    PrintTildeSepLines(params_buffer, 70, stdout);

    fprintf(stdout, "</PRE>\n");

#if !defined(NCBI_CLIENT_SERVER) && defined (NCBI_ENTREZ_CLIENT)

    fprintf(stdout, "<HR>\n");

    if(theInfo->show_tax_blast) {
        TXBHtmlReport(seqalign, stdout, theInfo->query_is_na,
                      theInfo->db_is_na, theInfo->database, 
                      NULL, NULL, theInfo->show_gi);
    }
    
#endif
    fprintf(stdout, "</BODY>\n</HTML>\n");

    seqannot = SeqAnnotFree(seqannot);

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

    return TRUE;
}

void WWWBlastPrintHeader(WWWBlastInfoPtr theInfo)
{

    fprintf(stdout, "<HTML>\n");
    fprintf(stdout, "<TITLE>BLAST Search Results</TITLE>\n"); 
    fflush(stdout);
    
    fprintf(stdout, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
            "VLINK=\"#660099\" ALINK=\"#660099\">\n");
    fprintf(stdout, "<A HREF=\"blast_form.map\"> \r"
            "<IMG SRC=\"images/blast_results.gif\" "
            "BORDER=0 ISMAP>\r</A><P>\n");    
    fprintf(stdout, "<PRE>\n");
    init_buff_ex(90);
    BlastPrintVersionInfo(theInfo->program, TRUE, stdout);
    fprintf(stdout, "\n");
    BlastPrintReference(TRUE, 90, stdout);
    fprintf(stdout, "\n");
    AcknowledgeBlastQuery(theInfo->query_bsp, 70, stdout, 
                          theInfo->believe_query, TRUE);
    free_buff();

    return;
}

Int2 Main(void)
{
    WWWBlastInfoPtr theInfo;
    
    UseLocalAsnloadDataAndErrMsg ();

    /* This function will read posting data, set-up config file and
       write small message into logfile (if it exists) */
    
    if((theInfo = WWWBlastReadArgs(NULL)) == NULL)
        return 1;
    
    /* Read options into structure */
    if(!WWWCreateSearchOptions(theInfo)) {
        return 1;
    }
    
    /* validate them */
    if(!WWWValidateOptions(theInfo)) {
        return 1;
    }

    /* Print BLAST Header */

    WWWBlastPrintHeader(theInfo);

    /* Do the search and Format output */

#ifdef NCBI_CLIENT_SERVER
    WWWBlastDoClientSearch(theInfo);
#else
    WWWBlastDoSearch(theInfo);
#endif

    WWWBlastInfoFree(theInfo);

    return 0;
}
