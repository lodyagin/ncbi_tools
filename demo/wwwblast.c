/* $Id: wwwblast.c,v 6.25 2000/10/23 20:19:57 dondosha Exp $
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
* $Revision: 6.25 $
*
* File Description:
*        Standalone WWW Blast CGI program.
*
* $Log: wwwblast.c,v $
* Revision 6.25  2000/10/23 20:19:57  dondosha
* Open and close AsnIo outside calls to BXMLPrintOutput function
*
* Revision 6.24  2000/10/18 20:19:31  shavirin
* Added title for OOF Blastx.
*
* Revision 6.23  2000/10/16 22:18:17  shavirin
* Added possibility to perform OOF blastx
*
* Revision 6.22  2000/10/16 20:27:03  shavirin
* Added possibility to run RPS Blast.
*
* Revision 6.21  2000/09/28 16:48:20  dondosha
* Changed MegaBlast related code to get a single SeqAlignPtr from server
*
* Revision 6.20  2000/09/28 15:16:55  shavirin
* Added message if request was limited to results of Entrez query.
*
* Revision 6.19  2000/09/27 22:17:04  shavirin
* Added possibility to limit search to results of entrez query.
*
* Revision 6.18  2000/09/13 22:28:10  dondosha
* Removed extra </PRE> that is now printed in PrintDefLinesFromSeqAlign
*
* Revision 6.17  2000/09/13 20:47:51  dondosha
* Small cleanup with closures of html blocks
*
* Revision 6.16  2000/09/12 21:57:28  dondosha
* Pass the correct scoring matrix to ShowTextAlignFromAnnot
*
* Revision 6.15  2000/09/11 17:51:07  shavirin
* Removed redundant <PRE> tag.
*
* Revision 6.14  2000/09/08 20:18:12  dondosha
* Print the title, background and GIF image before creating options
*
* Revision 6.13  2000/09/08 14:49:28  dondosha
* Allow graphical overview with multiple queries
*
* Revision 6.12  2000/09/07 18:01:38  dondosha
* Pass a callback to the server from TraditionalBlastReportEngineWithImage; allow multiple queries for all BLAST searches
*
* Revision 6.11  2000/09/05 18:00:25  dondosha
* Added query acknowledgement for each query in WWWBlastDoSearch for megablast
*
* Revision 6.10  2000/09/01 17:55:14  dondosha
* Call SeqEntryLoad at the beginning; corrections for megablast page
*
* Revision 6.9  2000/08/30 22:22:35  dondosha
* Enhance function WWWBlastDoSearch to handle megablast search
*
* Revision 6.8  2000/08/28 20:20:59  dondosha
* Added functionality for megablast web page
*
* Revision 6.7  2000/08/09 20:30:30  shavirin
* Added possibility to print XML output.
*
* Revision 6.6  2000/07/31 20:44:12  shavirin
* Minor change (initialized variable) from Haruna Cofer (haruna@detroit.sgi.com)
*
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

static  Boolean LIBCALLBACK
callback (BlastResponsePtr brp, Boolean PNTR cancel)

{
        fprintf(stdout, "</PRE>\n<PRE>");
	return TRUE;
}

static Boolean
TraditionalBlastReportEngineWithImage(SeqLocPtr slp, BioseqPtr bsp, BlastNet3Hptr bl3hp, WWWBlastInfoPtr theInfo)

{
    BlastDbinfoPtr dbinfo;
    BlastKABlkPtr ka_params=NULL, ka_params_gap=NULL;
    BlastPruneSapStructPtr prune;
    BLAST_MatrixPtr matrix;
    Int4Ptr PNTR txmatrix;
    Boolean query_is_na, db_is_na;
    Boolean status;
    CharPtr params_buffer=NULL;
    Int4 number_of_hits_private=0, length;
    SeqAlignPtr seqalign = NULL, sap, next_seqalign;
    SeqAnnotPtr seqannot=NULL;
    TxDfDbInfoPtr tx_dbinfo=NULL, tx_dbinfo_head;
    ValNodePtr mask_loc, mask_loc_start, other_returns, error_returns, vnp, vnp1=NULL;
    Uint1 align_type;
    Uint1 f_order[FEATDEF_ANY], g_order[FEATDEF_ANY];
    /* Variables for multiple query output */
    SeqLocPtr tmp_slp;
    Boolean done = FALSE;
    BLAST_OptionsBlkPtr options = theInfo->options;
    CharPtr program = theInfo->program, database = theInfo->database;
    Uint4 align_options = theInfo->align_options;
    DenseSegPtr dsp, next_dsp;
    AsnIoPtr xml_aip;

    MemSet((Pointer)(g_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    MemSet((Pointer)(f_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    
    if (bsp == NULL && slp == NULL)
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

    if (bsp) 
        seqalign = BlastBioseqNetCore(bl3hp, bsp, program, database, options, &other_returns, &error_returns, callback, NULL, &status);
    else if (options->is_megablast_search) {
       seqalign = MegaBlastSeqLocNetCore(bl3hp, slp, program, database, options, &other_returns, &error_returns, callback, &status);
    } else
        seqalign = BlastSeqLocNetCore(bl3hp, slp, program, database, options, &other_returns, &error_returns, callback, NULL, &status);

    
    BlastErrorPrintExtra(error_returns, TRUE, stdout);
    
    mask_loc = NULL;
    matrix = NULL;
    txmatrix = NULL;
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
            matrix = (BLAST_MatrixPtr) vnp->data.ptrvalue;
            if (matrix)
               txmatrix = BlastMatrixToTxMatrix(matrix);
            /*BLAST_MatrixDestruct(matrix);*/
            break;
        case SEQLOC_MASKING_NOTSET:
        case SEQLOC_MASKING_PLUS1:
        case SEQLOC_MASKING_PLUS2:
        case SEQLOC_MASKING_PLUS3:
        case SEQLOC_MASKING_MINUS1:
        case SEQLOC_MASKING_MINUS2:
        case SEQLOC_MASKING_MINUS3:
	    if (!options->is_megablast_search)
	       ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
	    break;
        default:
            break;
        }
    }	
    
    
    ReadDBBioseqFetchEnable ("blastall", database, db_is_na, TRUE);

    tmp_slp = slp;

    if(theInfo->xml_output)
       xml_aip = AsnIoOpen("stdout", "wx");

    while (seqalign) {
       if (!options->is_megablast_search)
          next_seqalign = NULL;
       else {
          sap = seqalign;
          while (sap != NULL) { 
             if (sap->next != NULL) {
                dsp = (DenseSegPtr) (sap->segs);
                next_dsp = (DenseSegPtr) (sap->next->segs);
                
                if (SeqIdComp(dsp->ids, next_dsp->ids) != SIC_YES) {
                   next_seqalign = sap->next;
                   sap->next = NULL;
                }
             } else
                next_seqalign = NULL;
             sap = sap->next;
          }
      
          dsp = (DenseSegPtr) (seqalign->segs);
          while (tmp_slp && SeqIdComp(dsp->ids, SeqLocId(tmp_slp)) != SIC_YES)
             tmp_slp = tmp_slp->next;
          if (tmp_slp == NULL) /* Should never happen */
             break;
          bsp = BioseqLockById(SeqLocId(tmp_slp));
          init_buff_ex(85);
          fprintf(stdout, "<HR><BR>");
          AcknowledgeBlastQuery(bsp, 70, stdout, FALSE, TRUE);
          free_buff();
          BioseqUnlock(bsp);
       }

       if(theInfo->xml_output) {
          printf("<PRE>");
          BXMLPrintOutput(xml_aip, seqalign, options, 
                          program, database, 
                          bsp, other_returns, 0, 0, NULL);
          AsnIoReset(xml_aip);
          printf("</PRE>");
          
       } else {
          
          seqannot = SeqAnnotNew();
          seqannot->type = 2;
          AddAlignInfoToSeqAnnot(seqannot, align_type);
          seqannot->data = seqalign;
          
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
          
          PrintDefLinesFromSeqAlignEx2(prune->sap, 80, stdout, theInfo->print_options, FIRST_PASS, NULL, theInfo->number_of_descriptions, database, theInfo->blast_type);
          free_buff();
          
          prune = BlastPruneHitsFromSeqAlign(seqalign, theInfo->number_of_alignments, prune);
          seqannot->data = prune->sap;
          
          if(theInfo->color_schema != 0 && 
             (!StringICmp(program, "blastn") || 
              !StringICmp(program, "blastp"))) {
             
             if(!DDV_DisplayBlastPairList(prune->sap, mask_loc, stdout, 
                                          query_is_na, align_options, 
                                          theInfo->color_schema)) { 
                fprintf(stdout, 
                        "\n\n!!!\n   "
                        "    --------  Failure to print alignment...  --------"
                        "\n!!!\n\n");
                fflush(stdout);
             }
          } else {
             
              if(options->is_ooframe) {
                  printf("<PRE>");
                  OOFShowBlastAlignment(seqalign, /*mask*/ NULL,
                                        stdout, align_options, txmatrix);
              } else {
                  if (align_options & TXALIGN_MASTER)
                      ShowTextAlignFromAnnot(seqannot, 60, stdout, f_order,
                                             g_order, align_options, txmatrix, 
                                             mask_loc, NULL);
                  else
                      ShowTextAlignFromAnnot(seqannot, 60, stdout, f_order, 
                                             g_order, align_options, txmatrix, 
                                             mask_loc, FormatScoreFunc);
              }
          }
          
          seqannot->data = seqalign;
          number_of_hits_private = prune->original_number; 
          prune = BlastPruneSapStructDestruct(prune);
          ObjMgrClearHold();
          ObjMgrFreeCache(0);
       }

       if (options->is_megablast_search)
          tmp_slp = tmp_slp->next;
       if (seqannot)
           seqannot = SeqAnnotFree(seqannot);
       seqalign = next_seqalign;
       fprintf(stdout, "<PRE>\n");
    }
    
    if(theInfo->xml_output)
       xml_aip = AsnIoClose(xml_aip);

    BLAST_MatrixDestruct(matrix);
    if (txmatrix)
       txmatrix = TxMatrixDestruct(txmatrix);

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

    printf("</PRE></BODY></HTML>\n");
    MemFree(params_buffer);
    free_buff();
    
    other_returns = ValNodeFree(other_returns);
    
    
    return status;
}

static Boolean
TraditionalBlastReportWithImage(BioseqPtr bsp, BlastNet3Hptr bl3hp, 
                                WWWBlastInfoPtr theInfo)
{
   return TraditionalBlastReportEngineWithImage(NULL, bsp, bl3hp, theInfo);
}
static Boolean
TraditionalBlastReportLocWithImage(SeqLocPtr slp, BlastNet3Hptr bl3hp,
                                   WWWBlastInfoPtr theInfo)
{
   return TraditionalBlastReportEngineWithImage(slp, NULL, bl3hp, theInfo);
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
    
    if (!theInfo->options->is_megablast_search)
       status = TraditionalBlastReportWithImage(theInfo->fake_bsp, bl3hp, theInfo);
    else 
       status = TraditionalBlastReportLocWithImage(theInfo->query_slp, bl3hp, theInfo);
    if (status == FALSE) {
        WWWBlastErrMessage(BLASTErrServer, NULL);
        return FALSE;
    }
    
    return TRUE;
}
#endif

Boolean WWWBlastDoSearch(WWWBlastInfoPtr theInfo)
{
    SeqAlignPtr  seqalign = NULL;
    ValNodePtr  mask_loc, mask_loc_start, vnp, other_returns, error_returns;
    TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
    BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
    BLAST_MatrixPtr matrix;
    Int4Ptr PNTR txmatrix;
    CharPtr params_buffer=NULL;
    SeqAnnotPtr seqannot = NULL;
    BlastPruneSapStructPtr prune;
    SeqAlignPtr PNTR seqalignp = NULL;
    Boolean is_megablast, done = FALSE, all_done = FALSE;
    SeqLocPtr query_slp;
    Int4 index;
    AsnIoPtr xml_aip;

    if(theInfo == NULL)
	return FALSE;
    
    is_megablast = theInfo->options->is_megablast_search;

    query_slp = theInfo->query_slp;

    if(!theInfo->xml_output) {
        PrintDbInformation(theInfo->database, !theInfo->db_is_na, 
                           70, stdout, TRUE);
    }

    if(theInfo->options->entrez_query != NULL && 
       theInfo->gi_list_total > 0) {

        printf("Your search was limited by an Entrez query: '%s'\n"
               "<!-- %d sequences-->\n<P>\n",
               theInfo->options->entrez_query, 
               theInfo->gi_list_total);
    }
    
    ReadDBBioseqFetchEnable ("blastall", theInfo->database, 
                             theInfo->db_is_na, TRUE);
    
    if(theInfo->xml_output)
       xml_aip = AsnIoOpen("stdout", "wx");

    while (!all_done) { /* Loop on complete BLAST searches */
        if (!query_slp || is_megablast)
            all_done = TRUE;
        
        other_returns = NULL;
        error_returns = NULL;
        
        if (!is_megablast) {
            if (query_slp) {
                seqalign = BioseqBlastEngineByLocEx(query_slp, theInfo->program, theInfo->database, theInfo->options, &other_returns, &error_returns, tick_callback, NULL, theInfo->gi_list, theInfo->gi_list_total);
            } else {
                seqalign = BioseqBlastEngineEx(theInfo->fake_bsp, theInfo->program, theInfo->database, theInfo->options, &other_returns, &error_returns, tick_callback, NULL, theInfo->gi_list, theInfo->gi_list_total);
            }
        } else {
            seqalignp = BioseqMegaBlastEngineByLoc(query_slp, theInfo->program, theInfo->database, theInfo->options, &other_returns, &error_returns, tick_callback, NULL, theInfo->gi_list, theInfo->gi_list_total, NULL);
        }
        
        BlastErrorPrint(error_returns);
        
        dbinfo = NULL;
        ka_params = NULL;
        ka_params_gap = NULL;
        params_buffer = NULL;
        mask_loc = NULL;
        matrix = NULL;
        txmatrix = NULL;
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
                if (matrix)
                    txmatrix = (Int4Ptr PNTR) BlastMatrixToTxMatrix(matrix);
                break;
            case SEQLOC_MASKING_NOTSET:
            case SEQLOC_MASKING_PLUS1:
            case SEQLOC_MASKING_PLUS2:
            case SEQLOC_MASKING_PLUS3:
            case SEQLOC_MASKING_MINUS1:
            case SEQLOC_MASKING_MINUS2:
            case SEQLOC_MASKING_MINUS3:
                if (!is_megablast)
                    ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
                break;
            default:
                break;
            }
        }	
        
        
        fflush(stdout);
        
        done = FALSE;
        
        for (index=0; !done; index++) {
            if (is_megablast)
                seqalign = seqalignp[index];
            else 
                done = TRUE;
            if (seqalign) {
                if(theInfo->xml_output) {
                    /* printf("<XM>"); */
                   BXMLPrintOutput(xml_aip, seqalign, theInfo->options, 
                                   theInfo->program, theInfo->database, 
                                   theInfo->fake_bsp, other_returns, 0, 0, NULL);
                   AsnIoReset(xml_aip);
                    /* printf("</XM>"); */
                } else {
                    seqannot = SeqAnnotNew();
                    seqannot->type = 2;
                    AddAlignInfoToSeqAnnot(seqannot, theInfo->align_type);
                    seqannot->data = seqalign;
                    
#if !defined(NCBI_CLIENT_SERVER) && defined (NCBI_ENTREZ_CLIENT)
                    if(theInfo->show_tax_blast) {
                        fprintf(stdout, "<a href=#taxblast>Taxonomy reports</a><BR>");
                    }
#endif
                    
                    if (is_megablast || query_slp) {
                        BioseqPtr bsp = BioseqLockById(SeqLocId(query_slp));
                        init_buff_ex(85);
                        fprintf(stdout, "<HR><BR>");
                        AcknowledgeBlastQuery(bsp, 70, stdout, FALSE, TRUE);
                        free_buff();
                        BioseqUnlock(bsp);
                    }
                    
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
                    } else {

                        if(theInfo->options->is_ooframe) {
                            printf("<PRE>");
                            OOFShowBlastAlignment(seqalign, /*mask*/ NULL,
                                                  stdout, theInfo->align_options, txmatrix);
                        } else {                        
                            if (theInfo->align_view != 0)
                                ShowTextAlignFromAnnot(seqannot, 60, stdout, NULL, NULL, 
                                                       theInfo->align_options, txmatrix, 
                                                       mask_loc, NULL);
                            else
                                ShowTextAlignFromAnnot(seqannot, 60, stdout, NULL, NULL, 
                                                       theInfo->align_options, txmatrix, mask_loc, 
                                                       FormatScoreFunc);
                        }
                    }
                    seqannot->data = NULL; /* Don't want to delete seqalign yet */
                    prune = BlastPruneSapStructDestruct(prune);
                    ObjMgrClearHold();
                    ObjMgrFreeCache(0);
                }
            } else if (!is_megablast) {
                fprintf(stdout, "\n\n ***** No hits found ******\n\n");
            }
            
            if (is_megablast) { 
                query_slp = query_slp->next;
                if (!query_slp)
                    done = TRUE;
            }
            if(!theInfo->xml_output && seqannot != NULL)
                seqannot = SeqAnnotFree(seqannot);
        }
        
        matrix = BLAST_MatrixDestruct(matrix);
        if (txmatrix)
            txmatrix = (Int4Ptr PNTR) TxMatrixDestruct(txmatrix);
        
        dbinfo_head = dbinfo;
        if(!theInfo->xml_output) {        
            fprintf(stdout, "<PRE>");
            init_buff_ex(85);
            while (dbinfo) {
                PrintDbReport(dbinfo, 70, stdout);
                dbinfo = dbinfo->next;
            }
        }
        
        dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
        
        if (ka_params) {
            if(!theInfo->xml_output) {        
                PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70,
                                  stdout, FALSE);
            }
            MemFree(ka_params);
        }
        
        if (ka_params_gap) {
            if(!theInfo->xml_output) {        
                PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K,
                                  ka_params_gap->H, 70, stdout, TRUE);
            }
            MemFree(ka_params_gap);
        }
        
        if(!theInfo->xml_output) {        
            
            PrintTildeSepLines(params_buffer, 70, stdout);
            
#if !defined(NCBI_CLIENT_SERVER) && defined (NCBI_ENTREZ_CLIENT)
            
            fprintf(stdout, "<HR>\n");
            
            if(theInfo->show_tax_blast) {
                TXBHtmlReport(seqalign, stdout, theInfo->query_is_na,
                              theInfo->db_is_na, theInfo->database, 
                              NULL, NULL, theInfo->show_gi);
            }
            
#endif
        }
        
        if (is_megablast) {
            Int4 i;
            for (i=0; i<index; i++)
                SeqAlignSetFree(seqalignp[i]);
            MemFree(seqalignp);
        } else if(seqalign != NULL)
            seqalign = SeqAlignSetFree(seqalign);
        
        MemFree(params_buffer);
        free_buff();
        
        mask_loc_start = mask_loc;
        while (mask_loc) {
            SeqLocSetFree(mask_loc->data.ptrvalue);
            mask_loc = mask_loc->next;
        }
        ValNodeFree(mask_loc_start);
        
        other_returns = ValNodeFree(other_returns);
        
        if (!is_megablast && query_slp) {
            query_slp = query_slp->next;
            if (!query_slp)
                all_done = TRUE;
        }
    } /* End of loop over all searches */

    if (theInfo->xml_output)
       AsnIoClose(xml_aip);
        
    ReadDBBioseqFetchDisable();

    if (!theInfo->xml_output)
       fprintf(stdout, "</PRE>\n</BODY>\n</HTML>\n");

    return TRUE;
}

void WWWBlastPrintTopHeader()
{
    fprintf(stdout, "<HTML>\n");
    fprintf(stdout, "<TITLE>BLAST Search Results</TITLE>\n"); 
    fflush(stdout);
    
    fprintf(stdout, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
            "VLINK=\"#660099\" ALINK=\"#660099\">\n");
    fprintf(stdout, "<A HREF=\"blast_form.map\"> \r"
            "<IMG SRC=\"images/blast_results.gif\" "
            "BORDER=0 ISMAP>\r</A><P>\n");    
}

void WWWBlastPrintHeader(WWWBlastInfoPtr theInfo)
{

    fprintf(stdout, "<PRE>\n");
    init_buff_ex(90);
    if (theInfo->options->is_megablast_search) {
        BlastPrintVersionInfo("megablast", TRUE, stdout);
    } else if (theInfo->options->is_rps_blast){
        BlastPrintVersionInfo("rps-blast", TRUE, stdout);
    } else if (theInfo->options->is_ooframe){
        BlastPrintVersionInfo("OOF BLASTX", TRUE, stdout);
    } else {
        BlastPrintVersionInfo(theInfo->program, TRUE, stdout);
        fprintf(stdout, "\n");
        BlastPrintReference(TRUE, 90, stdout);
    }
    fprintf(stdout, "\n");
    if (!theInfo->options->is_megablast_search && !theInfo->query_slp)
       AcknowledgeBlastQuery(theInfo->query_bsp, 70, stdout, 
                             theInfo->believe_query, TRUE);
    free_buff();

    return;
}

Int2 Main(void)
{
    WWWBlastInfoPtr theInfo;
    FILE *fp;

    UseLocalAsnloadDataAndErrMsg ();
            
    if (! SeqEntryLoad())
        return 1;
    
    ErrSetMessageLevel(SEV_WARNING);

    /* This function will read posting data, set-up config file and
       write small message into logfile (if it exists) */
    
    if((theInfo = WWWBlastReadArgs(NULL)) == NULL)
        return 1;
    
    if (!theInfo->xml_output) 
       WWWBlastPrintTopHeader();

    /* Read options into structure */
    if(!WWWCreateSearchOptions(theInfo)) {
        return 1;
    }
    
    /* validate them */
    if(!WWWValidateOptions(theInfo)) {
        return 1;
    }

    /* Print BLAST Header */

    if(!theInfo->xml_output) {        
        WWWBlastPrintHeader(theInfo);
    }

    /* Do the search and Format output */

#ifdef NCBI_CLIENT_SERVER
    WWWBlastDoClientSearch(theInfo);
#else
    WWWBlastDoSearch(theInfo);
#endif

    WWWBlastInfoFree(theInfo);

    return 0;
}
