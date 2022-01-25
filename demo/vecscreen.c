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
/* $Revision 1.0 $  
* $Log: vecscreen.c,v $
* Revision 6.2  2000/01/24 19:09:07  vakatov
* + #include <vecscrn.h>
*
* Revision 6.1  2000/01/20 18:58:30  madden
* Main file for vector screening
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
#include <salsap.h>
#include <vecscrn.h>


#define NUMARG 3

static Args myargs [NUMARG] = {
  { "Query File", 
	"stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  { "BLAST report Output File", 
	"stdout", NULL, NULL, FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  { "Database", 
	"UniVec", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL}
};

Int2 Main (void)
 
{
	AsnIoPtr aip;
	BioseqPtr query_bsp;
	BLAST_MatrixPtr matrix=NULL;
	BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
	Boolean db_is_na=TRUE, query_is_na=TRUE, show_gi, believe_query=FALSE;
	CharPtr ret_buffer=NULL, params_buffer=NULL;
	CharPtr database=NULL;
	Int4 number_of_descriptions, number_of_alignments;
	SeqAlignPtr  seqalign;
        SeqAnnotPtr seqannot;
	SeqEntryPtr sep;
	TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
	Uint1 align_type, align_view;
	Uint4 align_options, print_options;
	ValNodePtr  mask_loc, vnp, vnp1, other_returns, error_returns;

	CharPtr blast_inputfile, blast_outputfile;
	FILE *infp, *outfp;

        if (! GetArgs ("vecscreen", NUMARG, myargs))
        {
                return (1);
        }

	UseLocalAsnloadDataAndErrMsg ();

	if (! SeqEntryLoad())
		return 1;

	ErrSetMessageLevel(SEV_WARNING);

        blast_inputfile = myargs [0].strvalue;
        blast_outputfile = myargs [1].strvalue;
	database = myargs[2].strvalue;

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

	align_type = BlastGetTypes("blastn", &query_is_na, &db_is_na);

        believe_query = FALSE;

        print_options = 0;
        align_options = 0;
	print_options += TXALIGN_HTML;
	align_options += TXALIGN_HTML;
        align_options += TXALIGN_COMPRESS;
        align_options += TXALIGN_END_NUM;
        align_options += TXALIGN_MATRIX_VAL;
        align_options += TXALIGN_SHOW_QS;

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


		fprintf(outfp, "<PRE>");
		init_buff_ex(90);
		BlastPrintVersionInfo("blastn", TRUE, outfp);
		fprintf(outfp, "\n");
		BlastPrintReference(FALSE, 90, outfp);
		fprintf(outfp, "\n");
		AcknowledgeBlastQuery(query_bsp, 70, outfp, believe_query, TRUE);
                PrintDbInformation(database, !db_is_na, 70, outfp, TRUE);
                free_buff();

		error_returns = NULL;
		VSScreenSequence(query_bsp, NULL, database, &seqalign, &vnp1, &other_returns, &error_returns);
		VSPrintOverviewFromSeqLocs(vnp1, query_bsp->length, outfp);

		BlastErrorPrint(error_returns);

		for (vnp=vnp1; vnp; vnp = vnp->next)
		{
			SeqLocFree(vnp->data.ptrvalue);
		}
		vnp1 = ValNodeFree(vnp1);

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



       		ReadDBBioseqFetchEnable ("vecscreen", database, db_is_na, TRUE);

        	if (seqalign)
        	{
               	 	seqannot = SeqAnnotNew();
               	 	seqannot->type = 2;
			AddAlignInfoToSeqAnnot(seqannot, align_type);
                	seqannot->data = seqalign;
			if (outfp)
			{	/* Uncacheing causes problems with ordinal nos. vs. gi's. */
				ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, align_options, NULL, mask_loc, FormatScoreFunc);
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

		ReadDBBioseqFetchDisable();
		other_returns = ValNodeFree(other_returns);
		sep = SeqEntryFree(sep);
	}
	FileClose(infp);

	return 0;
}
	

