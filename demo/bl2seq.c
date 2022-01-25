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
***************************************************************************
*
* $Log: bl2seq.c,v $
* Revision 6.15  2000/10/30 19:05:44  madden
* Added -Y search space option
*
* Revision 6.14  2000/10/20 20:10:01  dondosha
* Revert previous change - would fail on some platforms
*
* Revision 6.13  2000/10/20 19:45:16  dondosha
* Changed -d option type to string to allow 8-byte integer db sizes
*
* Revision 6.12  2000/10/13 20:32:32  madden
* Added call to AcknowledgeBlastQuery
*
* Revision 6.11  2000/10/06 21:37:08  dondosha
* Set Mega BlAST default parameters correctly before calling the engine
*
* Revision 6.10  2000/10/05 22:44:28  dondosha
* Set block_width to 0 for Mega BLAST - it has different meaning
*
* Revision 6.9  2000/10/03 15:16:12  madden
* Set program name in call to BlastTwoSequencesEx
*
* Revision 6.8  2000/09/12 21:48:41  dondosha
* Pass the correct scoring matrix to ShowTextAlignFromAnnot
*
* Revision 6.7  2000/06/15 15:29:57  dondosha
* Fixed several memory leaks; tblastn and tblastx enabled
*
* Revision 6.6  2000/05/24 20:36:08  dondosha
* If megablast is used, set cutoff_s and cutoff_s2 parameters
*
* Revision 6.5  2000/04/10 15:23:33  dondosha
* Added option to use MegaBlast for search
*
* Revision 6.2  1999/11/26 20:16:11  vakatov
* Added <sqnutils.h> to pick up proto of 'UseLocalAsnloadDataAndErrMsg()'
*
* Revision 6.1  1999/07/06 18:48:20  madden
* Compares two sequences
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
#include <sqnutils.h>

		
#define NUMARG 20

static Args myargs [NUMARG] = {
  { "First sequence",
	NULL, NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  { "Second sequence",
	NULL, NULL, NULL, FALSE, 'j', ARG_FILE_IN, 0.0, 0, NULL},
  { "Program name: blastp, blastn, blastx, tblastn, tblastx. For blastx 1st sequence should be nucleotide, tblastn 2nd sequence nucleotide",
        "blastp", NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
  { "Gapped",
        "T", NULL, NULL, FALSE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
  { "alignment output file",
	"stdout", NULL, NULL, FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  { "theor. db size (zero is real size)", 
	"0", NULL, NULL, FALSE, 'd', ARG_INT, 0.0, 0, NULL},
  { "SeqAnnot output file",
	NULL, NULL, NULL, TRUE, 'a', ARG_FILE_OUT, 0.0, 0, NULL},
  { "Cost to open a gap (zero invokes default behavior)",
        "0", NULL, NULL, FALSE, 'G', ARG_INT, 0.0, 0, NULL},
  { "Cost to extend a gap (zero invokes default behavior)",
        "0", NULL, NULL, FALSE, 'E', ARG_INT, 0.0, 0, NULL},
  { "X dropoff value for gapped alignment (in bits) (zero invokes default behavior)",
        "0", NULL, NULL, FALSE, 'X', ARG_INT, 0.0, 0, NULL},
  { "Wordsize (zero invokes default behavior)",
        "0", NULL, NULL, FALSE, 'W', ARG_INT, 0.0, 0, NULL},
  { "Matrix",
        "BLOSUM62", NULL, NULL, FALSE, 'M', ARG_STRING, 0.0, 0, NULL},
  { "Penalty for a nucleotide mismatch (blastn only)",
        "-3", NULL, NULL, FALSE, 'q', ARG_INT, 0.0, 0, NULL},
  { "Reward for a nucleotide match (blastn only)",
        "1", NULL, NULL, FALSE, 'r', ARG_INT, 0.0, 0, NULL},
  { "Filter query sequence (DUST with blastn, SEG with others)",
        "T", NULL, NULL, FALSE, 'F', ARG_STRING, 0.0, 0, NULL},
  { "Expectation value (E)",
        "10.0", NULL, NULL, FALSE, 'e', ARG_FLOAT, 0.0, 0, NULL},
  { "Query strands to search against database (blastn only).  3 is both, 1 is top, 2 is bottom",
        "3", NULL, NULL, FALSE, 'S', ARG_INT, 0.0, 0, NULL},
  { "Produce HTML output",
        "F", NULL, NULL, FALSE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  { "Use Mega Blast for search",
        "F", NULL, NULL, FALSE, 'm', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Effective length of the search space (use zero for the real size)",
      "0", NULL, NULL, FALSE, 'Y', ARG_FLOAT, 0.0, 0, NULL}
};

Int2 Main (void)
 
{
	
	AsnIoPtr aip;
	BioseqPtr fake_bsp, fake_subject_bsp, query_bsp, subject_bsp;
	BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
	BLAST_OptionsBlkPtr options;
	Boolean seq1_is_na, seq2_is_na;
	CharPtr ret_buffer=NULL, params_buffer=NULL;
        DbtagPtr        dbtagptr;
	Uint1 align_type;
	Uint4 align_options;
	SeqAlignPtr  seqalign;
        SeqAnnotPtr seqannot;
	SeqEntryPtr sep, sep1;
	TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
	CharPtr blast_inputfile, blast_inputfile1, program_name, blast_outputfile;
	FILE *infp, *infp1, *outfp;
	ValNodePtr  mask_loc, mask_loc_start, vnp, other_returns, error_returns;
	BLAST_MatrixPtr matrix;
        Int4Ptr PNTR txmatrix;

        if (! GetArgs ("bl2seq", NUMARG, myargs))
        {
                return (1);
        }

	UseLocalAsnloadDataAndErrMsg ();

        if (! SeqEntryLoad())
                return 1;

	ErrSetMessageLevel(SEV_WARNING);

        blast_inputfile = myargs [0].strvalue;
        blast_inputfile1 = myargs [1].strvalue;
        blast_outputfile = myargs [4].strvalue;

	program_name = StringSave(myargs[2].strvalue);
	if (StringCmp(program_name, "blastn") && 
	    StringCmp(program_name, "blastp") && 
	    StringCmp(program_name, "blastx") && 
	    StringCmp(program_name, "tblastn") && 
	    StringCmp(program_name, "tblastx")) {
		ErrPostEx(SEV_FATAL, 0, 0, "Program name must be blastn, blastp, blastx, tblastn or tblastx\n");
		return (1);
	}
	   

	align_type = BlastGetTypes(program_name, &seq1_is_na, &seq2_is_na);

	if ((infp = FileOpen(blast_inputfile, "r")) == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile);
		return (1);
	}

	if ((infp1 = FileOpen(blast_inputfile1, "r")) == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile1);
		return (1);
	}

	if ((outfp = FileOpen(blast_outputfile, "w")) == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", blast_outputfile);
		return (1);
	}


	sep = FastaToSeqEntry(infp, seq1_is_na);
	if (sep != NULL)
	{
		query_bsp = NULL;
		if (seq1_is_na)
		{
			SeqEntryExplore(sep, &query_bsp, FindNuc);
		}
		else
		{
			SeqEntryExplore(sep, &query_bsp, FindProt);
		}

		fake_bsp = BlastMakeFakeBioseq(query_bsp, NULL);

		if (query_bsp == NULL)
		{
	  	 	ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
	   		return 2;
		}
	}

	if (myargs[6].strvalue != NULL || myargs[17].intvalue != 0)
		sep1 = FastaToSeqEntry(infp1, seq2_is_na);
	else
		sep1 = FastaToSeqEntryEx(infp1, seq2_is_na, NULL, FALSE);
	
	if (sep1 != NULL)
	{
		subject_bsp = NULL;
		if (seq2_is_na)
		{
			SeqEntryExplore(sep1, &subject_bsp, FindNuc);
		}
		else
		{
			SeqEntryExplore(sep1, &subject_bsp, FindProt);
		}

		if (subject_bsp == NULL)
		{
	  	 	ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
	   		return 2;
		}

		if (myargs[6].strvalue == NULL && myargs[17].intvalue == 0)
		{
			fake_subject_bsp = BioseqNew();
			fake_subject_bsp->descr = subject_bsp->descr;
			fake_subject_bsp->repr = subject_bsp->repr;
			fake_subject_bsp->mol = subject_bsp->mol;
			fake_subject_bsp->length = subject_bsp->length;
			fake_subject_bsp->seq_data = subject_bsp->seq_data;
			fake_subject_bsp->seq_data_type = subject_bsp->seq_data_type;
                	dbtagptr = DbtagNew();
                	dbtagptr->db = StringSave("BL_ORD_ID");
                	dbtagptr->tag = ObjectIdNew();
                	dbtagptr->tag->id = 0;
                	ValNodeAddPointer(&fake_subject_bsp->id, SEQID_GENERAL, dbtagptr);
		}
	}
		
	options = BLASTOptionNew(program_name, (Boolean) myargs[3].intvalue);
	options->is_megablast_search = (Boolean) myargs[18].intvalue;
	if (options->is_megablast_search) {
           options->gap_open = options->gap_extend = 0;
           options->wordsize = 32;
           options->block_width = 0;
	}

    	if (myargs[19].floatvalue)
        	options->searchsp_eff = (Nlm_FloatHi) myargs[19].floatvalue;


	options->filter_string = StringSave(myargs[14].strvalue);
	options->expect_value  = (Nlm_FloatHi) myargs [15].floatvalue;

        if (StringICmp("blastn", program_name) == 0)
        {
                options->penalty = myargs[12].intvalue;
                options->reward = myargs[13].intvalue;
        }

	options->db_length = myargs[5].intvalue;

	options->discontinuous = FALSE;

        if (myargs[7].intvalue != 0)
              options->gap_open = myargs[7].intvalue;
        if (myargs[8].intvalue != 0)
               options->gap_extend = myargs[8].intvalue;
        if (myargs[9].intvalue != 0)
	{
               options->gap_x_dropoff = myargs[9].intvalue;
	}
        if (myargs[10].intvalue != 0)
               options->wordsize = (Int2) myargs[10].intvalue;

	if (options->is_megablast_search) {
	   options->cutoff_s2 = options->wordsize - 4;
	   options->cutoff_s = options->wordsize;
        }
	MemFree(options->matrix);
	options->matrix = myargs[11].strvalue;

	options->strand_option = myargs[16].intvalue;
	if (myargs[6].strvalue || myargs[17].intvalue)
		seqalign = BlastTwoSequencesEx(query_bsp, subject_bsp, program_name, options, &other_returns, &error_returns);
	else
		seqalign = BlastTwoSequencesEx(fake_bsp, fake_subject_bsp, program_name, options, &other_returns, &error_returns);

		BlastErrorPrint(error_returns);

		dbinfo = NULL;
		ka_params = NULL;
		ka_params_gap = NULL;
		params_buffer = NULL;
		mask_loc = NULL;
		matrix = NULL;
                txmatrix = NULL;
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
                                        if (matrix)
                                           txmatrix = BlastMatrixToTxMatrix(matrix);
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


    align_options = 0;
    align_options += TXALIGN_MATRIX_VAL;
    align_options += TXALIGN_SHOW_QS;
    align_options += TXALIGN_COMPRESS;
    align_options += TXALIGN_END_NUM;
    if (StringICmp("blastx", program_name) == 0) {
        align_options += TXALIGN_BLASTX_SPECIAL;
    }

    if (myargs[17].intvalue)
       align_options += TXALIGN_HTML;
	
	AcknowledgeBlastQuery(query_bsp, 70, outfp, FALSE, myargs[17].intvalue);

      	seqannot = SeqAnnotNew();
        seqannot->type = 2;
	AddAlignInfoToSeqAnnot(seqannot, align_type);
        seqannot->data = seqalign;
	ShowTextAlignFromAnnot(seqannot, 60, outfp, NULL, NULL, align_options, txmatrix, mask_loc, FormatScoreFunc);

	aip = NULL;
	if (myargs[6].strvalue)
	{
		aip = AsnIoOpen (myargs[6].strvalue,"w");
	}

        if (aip && seqannot)
        {
                  SeqAnnotAsnWrite((SeqAnnotPtr) seqannot, aip, NULL);
                  AsnIoReset(aip);
                  aip = AsnIoClose(aip);
        }
        seqannot = SeqAnnotFree(seqannot);
		
        matrix = BLAST_MatrixDestruct(matrix);
        if (txmatrix)
           txmatrix = TxMatrixDestruct(txmatrix);
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
        mask_loc_start = mask_loc;
        while (mask_loc) {
            SeqLocSetFree(mask_loc->data.ptrvalue);
            mask_loc = mask_loc->next;
        }
        ValNodeFree(mask_loc_start);

	fake_bsp = BlastDeleteFakeBioseq(fake_bsp);
	other_returns = ValNodeFree(other_returns);
	options = BLASTOptionDelete(options);
	MemFree(program_name);
	FileClose(outfp);

	SeqEntryFree(sep);
	SeqEntryFree(sep1);

	return 0;
}
	

