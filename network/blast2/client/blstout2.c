/*
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
* File Name: blstout2.c
*
* Author:  Roman L. Tatusov, Warren Gish, Jonathan Epstein, Tom Madden, Yuri Sadykov
*
* Version Creation Date:   06/16/95
*
* $Revision: 6.1 $
*
* File Description: 
*       Creating BLAST report
*       Warren Gish's utilities needed for producing BLAST output
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*
*
* RCS Modification History:
* $Log: blstout2.c,v $
* Revision 6.1  1997/11/28 18:21:59  madden
* fprintf fixes
*
* Revision 6.0  1997/08/25 18:34:08  madden
* Revision changed to 6.0
*
* Revision 5.4  1997/07/02 17:58:50  madden
* Corrected pointer dereferenced to find gi
*
* Revision 5.3  1997/05/13 01:10:09  shavirin
* Added ability to print absolute and relative URLs on BLAST output
*
 * Revision 5.2  1996/11/27  14:01:16  madden
 * Replaced strlen by Nlm_StringLen.
 *
 * Revision 5.1  1996/06/06  21:08:57  madden
 * Removed unused variable.
 *
 * Revision 5.0  1996/05/28  14:09:11  ostell
 * Set to revision 5.0
 *
 * Revision 4.43  1996/05/28  13:02:21  madden
 * Changed max ID length (BLASTCLI_ID_MAX) from 30 to 35.
 *
 * Revision 4.42  1996/05/28  12:35:02  madden
 * Changed maximum ID length (BLASTCLI_ID_MAX) from 25 to 30.
 *
 * Revision 4.41  1996/05/27  14:10:41  madden
 * Added function ResidueSymbolTranslate that protects against
 * invalid residues.
 *
 * Revision 4.40  1996/05/22  15:29:43  vakatov
 * Modified "sort_letters()" function prototype to suit the MSVC++ compiler
 *
 * Revision 4.39  1996/05/17  22:02:09  madden
 * Optimized lookup of Residue for a given Symbol.
 *
 * Revision 4.38  1996/04/16  14:31:32  shavirin
 * Fixed reference problem in Mosaic for HTML output
 *
 * Revision 4.36  1996/04/12  15:40:23  madden
 * Changed "%d" to "%ld" in two places.
 *
 * Revision 4.35  1996/03/31  03:15:59  shavirin
 * Added HTML output in function TraditionalBlastOutputHTML()
 *
 * Revision 4.32  1996/03/22  17:29:52  madden
 * Removed include of netblap2.h
 *
 * Revision 4.31  1996/03/11  22:02:09  madden
 * portnew now calls BLAST0SeqData2Bioseq
 *
 * Revision 4.30  1995/12/04  20:22:34  madden
 * fixed bug when "-gi" option is used.
 *
 * Revision 4.29  1995/11/28  17:11:10  madden
 * Checked for defline before dereferencing, replaced strlen with StringLen.
 *
 * Revision 4.28  1995/11/07  17:56:27  madden
 * removed PrintTraditionalBlastPreface.
 *
 * Revision 4.27  1995/11/06  15:38:02  madden
 * Fixed "strand" for blastn.
 *
 * Revision 4.26  1995/11/03  18:36:56  madden
 * Fixed formatting of deflines.
 *
 * Revision 4.25  1995/11/03  16:03:40  madden
 * Replaced some fprintf calls in TradtionalTopBlas.. with BlastPrint...
 *
 * Revision 4.24  1995/11/02  23:15:18  madden
 * Used BlastPrintTabToColumn.
 *
 * Revision 4.23  1995/11/02  21:48:45  madden
 * Used BlastPrintIntegerAdd and BlastPrintDoubleAdd.
 *
 * Revision 4.22  1995/11/02  21:13:37  madden
 * Converted to using BlastPrint functions in some parts.
 *
 * Revision 4.21  1995/11/02  14:24:19  madden
 * Changed format_an_id to calculate real number of digits for a gi ID
 * rather than just using nine.
 *
 * Revision 4.20  1995/11/01  15:40:56  madden
 * removed TraditionalBlastReportSetUp and TraditionalBlastReportCleanUp.
 *
 * Revision 4.19  1995/10/27  21:16:31  madden
 * Removed TraditionalHistBlastOutput (moved to blast2.c).
 *
 * Revision 4.18  1995/10/26  17:10:01  madden
 * moved wrap, PrintBlastPreface, acknowledge_blast_request, and
 * BlastPrintValNodeStack to blast2.c.
 *
 * Revision 4.17  1995/09/05  20:50:01  madden
 * Check that definition is not just a NULLB.
 *
 * Revision 4.16  1995/09/05  20:03:48  madden
 * Newline added if only id was present.
 *
 * Revision 4.15  1995/08/29  21:51:04  madden
 * removed unused variables that were lint complaints.
 *
 * Revision 4.14  1995/08/29  21:28:12  madden
 * Histogram corrected to agree with stand-alone histogram.
 *
 * Revision 4.13  1995/08/29  18:26:36  madden
 * Sequence ID now printed out even if FASTA definition not found.
 *
 * Revision 4.12  1995/08/28  17:33:49  madden
 * replaced write to stderr with call to ErrPostEx.
 *
 * Revision 4.11  1995/08/24  17:23:01  madden
 * Added support for printing gi's in the FASTA id line.
 *
 * Revision 4.10  1995/08/23  22:14:33  madden
 * removed "qoffset" in the alignment functions.
 *
 * Revision 4.9  1995/08/23  13:32:55  madden
 * Added support for BLAST "V" and "B" options.
 *
 * Revision 4.8  1995/08/04  15:06:41  madden
 * Added function BlastPrintValNodeStack.
 *
 * Revision 4.7  1995/08/04  12:44:04  kans
 * changed return NULL to return FALSE
 *
 * Revision 4.6  1995/08/03  13:17:39  madden
 * Added check to PrintBlastPreface that the BLAST0PrefacePtr is not NULL.
 *
 * Revision 4.5  1995/08/02  16:26:24  madden
 * fixed lint nits (unused variable, functions that could be static).
 *
 * Revision 4.4  1995/08/02  16:18:41  madden
 * Added new preface and acknowledgment functions, enabled the "-qoffset"
 * option.
 *
 * Revision 4.3  1995/08/01  20:44:33  madden
 * traditional functions have BlastReportStructPtr as argument,
 * fixed purify errors.
 *
 * Revision 4.2  1995/08/01  14:51:32  madden
 * replaced the array ralpha with a call to GetResidueForSymbol.
 *
 * Revision 4.1  1995/07/28  16:16:37  madden
 * Made functions that format individual parts of the traditional
 * BLAST report non-static.
 *
 * Revision 4.0  1995/07/26  13:55:34  ostell
 * force revision to 4.0
 *
 * Revision 1.4  1995/07/19  15:00:44  madden
 * Changed some spacings in TraditionalHistBlastOutput.
 *
 * Revision 1.3  1995/07/18  21:06:51  madden
 * changed p-value etc. to p_value etc.
 *
 * Revision 1.2  1995/06/22  17:03:44  madden
 * HitDataPtr replaced by BLAST0ResultPtr.
 *
 * Revision 1.1  1995/06/16  11:26:47  epstein
 * Initial revision
 *
*/

#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#define NLM_GENERATED_CODE_PROTO
#include <objblst2.h>
#include <blast2.h>

#define BLASTCLI_MATRIX_SIZE 26
#define QUERY_HREF "http://www.ncbi.nlm.nih.gov"

static Pointer find PROTO((BLAST0ResponsePtr p, Uint1 v));
static void	doindent PROTO((FILE *fp, int ncols));
static CharPtr	gish_str_strand PROTO ((Uint2 strand));
static Int2 format_an_id PROTO ((CharPtr id, ValNodePtr vnp, Int2 max_id_length, Boolean get_gi));
static Int2 get_max_ID_length PROTO ((BLAST0HitListPtr hlp, Boolean get_gi));
static void get_frame PROTO ((BLAST0SeqIntervalPtr loc, CharPtr array, Int4 total_length));
static Boolean LIBCALL TraditionalBlastOutputInternal(BLAST0ResultPtr hdp, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp, Int4 type);

/**** MISC  ****/

static CharPtr
gish_str_strand(Uint2 strand)
{
    switch (strand)
    {
    case BLAST0_Seq_interval_strand_plus:
    case BLAST0_Seq_interval_strand_plus_rf:
        return "Plus"; /* plus */
    case BLAST0_Seq_interval_strand_minus:
    case BLAST0_Seq_interval_strand_minus_rf:
        return "Minus"; /* minus */
    case BLAST0_Seq_interval_strand_both:
    default: /* unknown, or unable to convert */
        return "Undefined";
    }
}

static BLAST0ResponsePtr _find(BLAST0ResponsePtr p, Uint1 v)
{
	BLAST0ResponsePtr b;

	for (b = p; b != NULL && b->choice != v; b = b->next) 
		;
	return b;
}

static Pointer find(BLAST0ResponsePtr p, Uint1 v)
{
	BLAST0ResponsePtr b;

	for (b = p; b != NULL && b->choice != v; b = b->next) 
		;

	return (b == NULL) ? NULL : b->data.ptrvalue;
}

/*
	static Int2 PrivateGetNumberOfDigits(Int4 number)

	This function calculates the lenght of "number" when printed.
*/

static Int2 
PrivateGetNumberOfDigits(Int4 number)
{

	Char temp[20];
	int length;

	length = sprintf(temp, "%ld", (long) number);
	return (Int2) length;
}

/***************************************************************************
*
*	Format the query and database acknowledgement of BLAST
*
*Query=  195747|Mitochondrion Rupicapra rupi
*        (646 letters)
*
*Database:  Non-redundant PDB+GBupdate+GenBank+EMBLupdate+EMBL
*           248,447 sequences; 236,652,708 total letters.
*Searching..................................................done
*
**************************************************************************/
Boolean LIBCALL
TraditionalHeadBlastOutput(BlastReportStructPtr blrp)
{
	BLAST0ResponsePtr a, blresp;
	BLAST0StatusPtr status;
	BLAST0SequencePtr query;
	BLAST0DbDescPtr db;
	BLAST0JobDescPtr jobd;
	Boolean found_top, found_bottom;
	Int2 strand_value;
	Int4 progress;
	CharPtr program, string;
	ValNodePtr parms;
	FILE *fp;

	if (blrp == NULL)
		return FALSE;

	blresp = blrp->blresp;	
	fp = blrp->bpsp->fp;
	program = blrp->program;
	
/* The next section looks through the "parms" (i.e., BLAST command-line options
to set values needed by the function acknowledge_blast_request */
	found_top = FALSE;
	found_bottom = FALSE;
	if ((parms = find(blresp, BLAST0Response_parms)) != NULL)
	{
		while (parms != NULL)
		{
			string = parms->data.ptrvalue;
			if (StringCmp(string, "-top") == 0)
				found_top = TRUE;
			else if (StringCmp(string, "-bottom") == 0)
				found_bottom = TRUE;
			parms = parms->next;
		}
	}

	strand_value=0;
	if (StringCmp(program, "blastx") == 0 || 
		StringCmp(program, "tblastx") == 0)
	{
		strand_value = 3; /* if neither top or bottom are found
				then both are TRUE! */
		if (found_top == TRUE && found_bottom == TRUE)
			strand_value = 3;
		else if (found_top == TRUE && found_bottom == FALSE)
			strand_value = 1;
		else if (found_top == FALSE && found_bottom == TRUE)
			strand_value = 2;
	}

	if ((query = find(blresp, BLAST0Response_query)) != NULL) 
	{
		acknowledge_blast_request(program, query, strand_value, blrp->qoffset, fp, blrp->type);
	}

	if ((db = find(blresp, BLAST0Response_dbdesc)) != NULL) {
		string = db->def;
                if(blrp->type & HTML_OUTPUT) 
                  blast2_wrap(fp, "\n<b>Database:</b>  ", string, Nlm_StringLen(string), 79, 11);
                else  /* blrp->type == TEXT_OUTPUT */
                  blast2_wrap(fp, "\nDatabase:  ", string, Nlm_StringLen(string), 79, 11);

		fprintf(fp, "           %s sequences; %s total letters.\n",
			Nlm_Ultostr((long)(db->count), 1), Nlm_Ultostr(db->totlen, 1));
	}

	a = blresp;
	while ((a = _find(a, BLAST0Response_job_start)) != NULL) {
		jobd = a->data.ptrvalue;
		if (jobd->jid != BLAST0_Job_desc_jid_search) {
			a = a->next;
			continue;
		}
		fprintf(fp, "%s", jobd->desc);
		fflush(fp);
		progress = BLAST0Response_job_progress;
		for (a = a->next; a != NULL && a->choice == progress; a = a->next) {
			fprintf(fp, ".");
			fflush(fp);
		}
		if (a != NULL && a->choice == BLAST0Response_job_done) {
			fprintf(fp, "done\n");
		} else {
			fprintf(fp, "Interrupted!!!\n");
		}
		fflush(fp);
	}
	
	if ((status = find(blresp, BLAST0Response_status)) != NULL) {
		if (status->code != 0) {
			fprintf(fp, "ERROR: %ld\n", (long) status->code);
			fprintf(fp, "reason: %s\n", status->reason);
		}
	}
	fflush(fp);

	return TRUE;
}

#define DEF_AND_ID_CUTOFF 59
#define BLASTCLI_ID_MAX 35

/**************************************************************************
*	
*	This function formats the One-line summaries of the traditional
*	BLAST output.  This includes the High-Scoring Segment Pair (HSP) 
*	id, the defline of the HSP, the high score of that HSP, the P 
*	value, and the number of HSP's found.  If a translation 
*	was performed on the database, the frame is also printed.  An
*	example is:
*
*                                                                     Smallest
*                                                                       Sum
*                                                              High  Probability
*Sequences producing High-scoring Segment Pairs:              Score  P(N)      N
*
*dbj|D32197|RURMTCB27 Mitochondrion Rupicapra rupicapra ge...  3230  7.0e-264  1
*dbj|D32191|CCRMTCB21 Mitochondrion Capricornis crispus ge...  2609  1.8e-212  1
*gb|U17862|OMU17862   Ovibos moschatus moschatus cytochrom...  2600  1.0e-211  1
*dbj|D32195|CCRMTCB25 Mitochondrion Capricornis sumatrensi...  2582  5.6e-210  1
*gb|U17861|NCU17861   Nemorhaedus caudatus cytochrome b ge...  2573  1.8e-209  1
*
*	and so on.
******************************************************************************/

Boolean LIBCALL
TraditionalTopBlastOutput(BlastReportStructPtr blrp)
{
    BlastPrintStructurePtr bpsp;
    BLAST0HitListPtr hlp;
    BLAST0ResultPtr hdp;
    BLAST0SegmentPtr segp;
    /* Int4 dim; */
    BLAST0SeqDescPtr sdp, sdp1;
    CharPtr defline=NULL, ptr;
    Char frame_array[3], buffer[DEF_AND_ID_CUTOFF+5];
    ScorePtr score; 
    CharPtr scoreDescr;
    BLAST0HSPPtr hsp;
    Int4 highScore, string_length;
    FloatHi poissonScore;
    Int2 max_id_length, max_def_length, length, length_left, old_length, total_length;
    Int4 nScore, query_length;
    Int4 TemphighScore;
    FloatHi TemppoissonScore;
    Int4 TempnScore, max_defline_number=blrp->num_of_defline, defline_num;
    CharPtr str1, program;
    Char str2[DEF_AND_ID_CUTOFF+4];
    Char id[BLASTCLI_ID_MAX+1];

    /*   HTML stuff  */

    ValNodePtr vnp1;
    Int4 HTML_gi;
    Char HTML_tmp[512];
    Char HTML_tmp1[512];
    CharPtr HTML_database;
    register Int4 i, n;
    Boolean use_text_id = FALSE;
    CharPtr HTML_text_id;
    Int4 HTML_query_len;
    Int4 num_of_alignments = 0, max_alignment_number=blrp->num_of_align;
    /* ------------- */

    if (blrp == NULL)
	return FALSE;

    bpsp = blrp->bpsp;
 
    hdp = blrp->result;
    query_length = blrp->query_length;
    program = blrp->program;

    
    BlastPrintStart(bpsp, 0, 0);
    BlastPrintNewLine(bpsp);
    BlastPrintTabToColumn(bpsp, 70);
    BlastPrintStringAdd(bpsp, "Smallest");
    BlastPrintNewLine(bpsp);
    BlastPrintTabToColumn(bpsp, 72);
    BlastPrintStringAdd(bpsp, "Sum");
    if (StringCmp("blastn", program) == 0 || 
	StringCmp("blastp", program) == 0)
    {
    	BlastPrintNewLine(bpsp);
    	BlastPrintTabToColumn(bpsp, 63);
    	BlastPrintStringAdd(bpsp, "High  Probability");
    	BlastPrintNewLine(bpsp);
    	BlastPrintStringAdd(bpsp, "Sequences producing High-scoring Segment Pairs:");
    	BlastPrintTabToColumn(bpsp, 62);
    	BlastPrintStringAdd(bpsp, "Score  P(N)      N");
    	BlastPrintNewLine(bpsp);
    }
    else
    {
    	BlastPrintNewLine(bpsp);
    	BlastPrintTabToColumn(bpsp, 54);
    	BlastPrintStringAdd(bpsp, "Reading  High  Probability");
    	BlastPrintNewLine(bpsp);
    	BlastPrintStringAdd(bpsp, "Sequences producing High-scoring Segment Pairs:");
    	BlastPrintTabToColumn(bpsp, 56);
    	BlastPrintStringAdd(bpsp, "Frame Score  P(N)      N");
    	BlastPrintNewLine(bpsp);
    }

    if (hdp == NULL || hdp->hitlists == NULL) {
	    if (hdp->hitlists == NULL) {
    		BlastPrintStart(bpsp, 0, 0);
    		BlastPrintNewLine(bpsp);
    		BlastPrintStringAdd(bpsp, "      *** NONE ***");
	    }
        return FALSE;
    }
    BlastPrintEnd(bpsp);
 
/*    dim = hdp->dim;*/ /* should be 2, or 3 for BLAST3 */

    max_id_length = get_max_ID_length(hdp->hitlists, blrp->get_gi);
    if (StringCmp("blastn", program) == 0 || 
	StringCmp("blastp", program) == 0)
    	max_def_length = DEF_AND_ID_CUTOFF - max_id_length;
    else
    	max_def_length = DEF_AND_ID_CUTOFF - max_id_length - 3;
 
    defline_num=0;
    for (hlp = hdp->hitlists; hlp != NULL; hlp = hlp->next)
    {
    	defline_num++;
	if (defline_num > max_defline_number)
		break;
	
        if (hlp->hsps != NULL)
        {
            hsp = hlp->hsps; /* first HSP */
            segp = hsp->segs;
            if (segp != NULL && segp->next != NULL &&
                segp->next->loc != NULL)
            {
                if (hlp->seqs != NULL)
                {
                    defline = NULL;
                    sdp = hlp->seqs->desc;
                    if (sdp->defline != NULL) {
                    	defline = sdp->defline;
                    }
                }

                highScore = -1;
                TemphighScore = highScore;
                poissonScore = 99999;
                TemppoissonScore = poissonScore;
                nScore = 0;
                for (; hsp != NULL; hsp = hsp->next)
                {   /* If value associated with BLAST0_Score_sid_sum_n is 1,
			then this is not set in the ASN.1 */
		    TempnScore = 1; 
            	    segp = hsp->segs;
                    for (score = (ScorePtr) hsp->scores; score != NULL; score = score->next)
                    {
			if (score->id->str == NULL)
			    scoreDescr = "";
			else
			    scoreDescr = score->id->str;
			if (StrCmp(scoreDescr, "score") == 0)
			{
                            TemphighScore = score->value.intvalue;
			} else if (StrCmp(scoreDescr, "p_value") == 0 ||
			           StrCmp(scoreDescr, "poisson_p") == 0 ||
			           StrCmp(scoreDescr, "sum_p") == 0)
			{
                            TemppoissonScore = score->value.realvalue;
			} else if (StrCmp(scoreDescr, "poisson_n") == 0 ||
			           StrCmp(scoreDescr, "sum_n") == 0)
			{
                            TempnScore = score->value.intvalue;
			} else { /* default */
			}
                    }

                    if (TemphighScore > highScore)
                    {
                        highScore = TemphighScore;
			if (segp != NULL)
			{
			   if (StringCmp(program, "blastx") == 0)
				get_frame(segp->loc, frame_array, query_length);
			   else if (StringCmp(program, "tblastn") == 0)
				get_frame(segp->next->loc, frame_array, hlp->seqs->length);
			   else if (StringCmp(program, "tblastx") == 0)
				get_frame(segp->next->loc, frame_array, hlp->seqs->length);
			}
		    }
		    if (TemppoissonScore < poissonScore)
		    {
                        poissonScore = TemppoissonScore;
                        nScore = TempnScore;
                    }
                }

/* Add secondary id's (and their deflines) onto the first defline if it isn't 
too long already.  */
		str1 = defline;
		string_length = 0;
                if (str1 != NULL)
                {
		   	string_length = StringLen(str1);
/* secondary ID's get added on if the first defline doesn't exceed the max.*/

                   	if (string_length <= max_def_length)
			{
				if (sdp && sdp->next)
				{	
				    ptr = &buffer[0];
				    StrNCpy(ptr, str1, string_length);
				    ptr += string_length;
				    *ptr = NULLB;
				    total_length = string_length;
				    for (sdp1=sdp->next; sdp1; sdp1=sdp1->next)
				    {
					length = format_an_id(&id[0], sdp1->id, 0, blrp->get_gi);
					old_length = total_length;
					total_length += length;
					total_length += 2;
					StringCpy(ptr, " >");
					ptr += 2;
					if (total_length > max_def_length)
					{
/* Add on part of last (too long) id to agree with Traditional Output.
There must be some slack (extra room) in buffer for this. */
					    length_left = max_def_length-old_length;
					    if (length_left > 0)
					    {
					    	StrNCpy(ptr, id, length_left);
						ptr += length_left;
						*ptr = '\0';
					    }
					    else
					    {
					    	StrNCpy(ptr, "...", 3);
						ptr += 3;
						*ptr = '\0';
					    }
					    break;
					}
					StrNCpy(ptr, id, length);
					ptr += length;
					*ptr = ' '; ptr++;
					*ptr = NULLB;
					total_length++;
					if (sdp1->defline)
					{
						length = StringLen(sdp1->defline);
					}
					else
					{
						length = 0;
					}
					if ((total_length+length) > max_def_length)
					{
/* Add on part of last (too long) id to agree with Traditional Output.
There must be some slack (extra room) in buffer for this. */
					    length = max_def_length-total_length;
					    length += 3;
					    if (length > 0)
					    {
					    	StrNCpy(ptr, sdp1->defline, length);
						ptr += length;
					    }
					    else
					    {
					    	StrNCpy(ptr, sdp1->defline, 1);
						ptr++;
					    }
					    *ptr = NULLB;
					    break;
					}
					else
						total_length += length; 
					StrNCpy(ptr, sdp1->defline, length);
					ptr += length;
				    	*ptr = '\0';
				    }
				    str1 = &buffer[0];
				}
			}
		}

                if (str1 != NULL)
		   	string_length = StringLen(str1);
		else
		   	string_length = 0;
                if (string_length > max_def_length)
                {
                      	StrNCpy (str2, str1, max_def_length);
                        str2[max_def_length-3] = '.';
                        str2[max_def_length-2] = '.';
                        str2[max_def_length-1] = '.';
                        str2[max_def_length] = '\0';
		}
		else
		{
                	if (str1 != NULL)
                      		StrNCpy (str2, str1, string_length);
			str2[string_length] = NULLB;
		}
		str1 = str2;

		

		format_an_id(&id[0], sdp->id, max_id_length, blrp->get_gi);

/* Cutoff the id's after max_id_length characters to agree with Traditional 
Output */

		id[max_id_length] = '\0';

		BlastPrintStart(bpsp, 0, 0);


                /* HTML stuff */ 

                if(blrp->type & HTML_OUTPUT) {

                  /* setting database to use in Entrez query */
                  HTML_database = malloc(3);
                  if(!StringICmp(program, "blastx") || 
                     !StringICmp(program, "blastp"))
                    HTML_database = "p";
                  else
                    HTML_database = "n";  

                  /* determine HTML_gi variable - main reference point */
                  HTML_gi = 0;
                  use_text_id = FALSE;
                  for (vnp1=sdp->id; vnp1 != NULL; vnp1=vnp1->next) {
                    if (vnp1->choice == BLAST0SeqId_giid) {
                      HTML_gi=vnp1->data.intvalue;
                      break;
                    }
                  }
                  if (HTML_gi == 0) { /* use text id as main reference point */
                    use_text_id = TRUE; 
                    /* Look for text id's */
                    for (vnp1=sdp->id; vnp1 != NULL; vnp1=vnp1->next) {
                      if (vnp1->choice == BLAST0SeqId_textid) {
                        HTML_text_id = StringSave(vnp1->data.ptrvalue);
                        break;
                      }
                    }
                  }

                  /* filtering ">" and "<" characters to "/" */

                  for(i = 0; i < StringLen(str1); i++)
                    if ((str1[i] == '<') || (str1[i] == '>'))
                      str1[i] = '/';

                  /* preparing Entrez reference */

                  if (use_text_id) {
                    sprintf(HTML_tmp, "<a href=\""
                            "%s/htbin-post/Entrez/query?form=6&dopt=g&db=%s&"
                            "uid=%s\">", 
                            blrp->type & ABSOLUTE_LINKS ? QUERY_HREF : "",
                            HTML_database,  HTML_text_id);

                    HTML_query_len = StringLen(HTML_tmp) + 4;
                    sprintf(HTML_tmp1, "%s</a>", id);
                    strcat(HTML_tmp, HTML_tmp1);
                  } else {
                      sprintf(HTML_tmp, "<a href=\""
                            "%s/htbin-post/Entrez/query?form=6&dopt=g&db=%s&"
                            "uid=%08ld\">", 
                              blrp->type & ABSOLUTE_LINKS ? QUERY_HREF : "",
                              HTML_database, (long) HTML_gi);
                    HTML_query_len = StringLen(HTML_tmp) + 4;
                    sprintf(HTML_tmp1, "%s</a>", id);
                    strcat(HTML_tmp, HTML_tmp1);
                  }
                  /* some formatting to follow Traditional output */
                  for ( i =StringLen(HTML_tmp); 
                        i < (HTML_query_len + max_id_length +1) ; i++)
                    strcat(HTML_tmp, " ");

                  strcat(HTML_tmp, str1);


                  if (StringCmp("blastx", program) == 0 || 
                      StringCmp("tblastn", program) == 0 || 
                      StringCmp("tblastx", program) == 0)  {
                    
                    for ( i =StringLen(HTML_tmp); i < (HTML_query_len + 58) ; i++)
                      strcat(HTML_tmp, " ");
                    strcat(HTML_tmp, frame_array);
                  }
                  
                  fprintf(bpsp->fp, "%s",  HTML_tmp); /* actual printing */

                  for ( i =StringLen(HTML_tmp); i < (HTML_query_len +61) ; i++)
                    BlastPrintCharAdd(bpsp, ' ');
                  
                } else {   /* if(blrp->type == HTML_OUTPUT) */

                  BlastPrintStringAdd(bpsp, id);
                  BlastPrintTabToColumn(bpsp, max_id_length+1);
                  BlastPrintCharAdd(bpsp, ' ');
                  BlastPrintStringAdd(bpsp, str1);
                  BlastPrintCharAdd(bpsp, ' ');

                  if (StringCmp("blastx", program) == 0 || 
                      StringCmp("tblastn", program) == 0 || 
                      StringCmp("tblastx", program) == 0)
                    {
                      BlastPrintTabToColumn(bpsp, DEF_AND_ID_CUTOFF);
                      BlastPrintStringAdd(bpsp, frame_array);
                      BlastPrintCharAdd(bpsp, ' ');
                    } else {
                      BlastPrintTabToColumn(bpsp, DEF_AND_ID_CUTOFF+3);
                    }
		} /* blrp->type == TEXT_OUTPUT */

                if(blrp->type & HTML_OUTPUT &&
                   num_of_alignments < max_alignment_number) {

		  sprintf(HTML_tmp,"%ld", (long) highScore);
		  n = StringLen(HTML_tmp);
		  for(i = 0; i < (5 - n); i++)
		    BlastPrintCharAdd(bpsp, ' ');

                  if(!use_text_id)
                    sprintf(HTML_tmp, 
                            "<a href = #%08ld>%ld</a>", 
                            HTML_gi, (long) highScore);
                  else
                    sprintf(HTML_tmp, 
                            "<a href = #%s>%ld</a>", 
                            HTML_text_id, (long) highScore);

                  BlastPrintStringAdd(bpsp, HTML_tmp);

                } else /* blrp->type == TEXT_OUTPUT */

                  BlastPrintIntegerAdd(bpsp, "%5ld", highScore);

                num_of_alignments++;

/* The following code was lifted from "print_header" in blastapp/lib/prt_hdr.c to match the traditional output. */
        	if (poissonScore <= 0.99) 
		{
               	    if (poissonScore != 0.)
                        BlastPrintDoubleAdd(bpsp, "  %#-8.2lg", (double) poissonScore);
               	    else
                        BlastPrintDoubleAdd(bpsp, "  %#-8.2lg", 0.);
        	}
        	else if (poissonScore <= 0.999)
               		BlastPrintDoubleAdd(bpsp, "  %#-8.3lg", (double) poissonScore);
        	else if (poissonScore <= 0.9999)
                	BlastPrintDoubleAdd(bpsp, "  %#-8.4lg", (double) poissonScore);
        	else if (poissonScore <= 0.99999)
                	BlastPrintDoubleAdd(bpsp, "  %#-8.5lg", (double) poissonScore);
        	else if (poissonScore <  0.9999995)
                	BlastPrintDoubleAdd(bpsp, "  %#-8.6lg", (double) poissonScore);
        	else
                	BlastPrintDoubleAdd(bpsp, "  %#-8.7lg", (double) 1.0);
		BlastPrintIntegerAdd(bpsp, "%3ld", (int) nScore);
		BlastPrintEnd(bpsp);
            }
        }
    }
    return TRUE;
}

/* The traditional BLAST output uses fasta style id's and not gi's.  Look 
for a textid (i.e., fasta) and use this, otherwise use the first id. 
This function operates in two modes: if id is NULL, then the length of 
the id is returned; if id is not NULL, then the formatted id is returned in 
"id", which should already contain space for this operation.  The ValNodePtr 
vnp is actually BLAST0SeqDesc.id 
The parameter max_id_length should be greater than zero if id is not NULL
and one wishes to pad the end of the id with white spaces.
If the gi should be prepended to the FASTA output then set get_gi should be
TRUE.
*/

static Int2
format_an_id (CharPtr id, ValNodePtr vnp, Int2 max_id_length, Boolean get_gi)

{
	Boolean found_id;
	Char temp[100];
	Int2 index, length;
	Int4 gi=0;
	ValNodePtr vnp1;

	if (id != NULL)
		id[0] = '\0';

	found_id = FALSE;
	length=0;

/* IF gi was requested, then look for that first. */
	if (get_gi)
	{
		for (vnp1=vnp; vnp1 != NULL; vnp1=vnp1->next)
		{
			if (vnp1->choice == BLAST0SeqId_giid)
			{
				gi=vnp->data.intvalue;
				length = PrivateGetNumberOfDigits(gi);
				length += 3;	/* 3 extra chars for "gi|" */
			}
		}
	}

/* Look for text id's */
	for (vnp1=vnp; vnp1 != NULL; vnp1=vnp1->next)
	{
		if (vnp1->choice == BLAST0SeqId_textid)
		{
		    if (length > 0)
			length++; 	/* add in room for dividing "|" */ 
		    length += StringLen(vnp1->data.ptrvalue);
		    if (id != NULL)
		    {
			if (get_gi)
	  	    		sprintf(temp, "gi|%ld|%s", (long) gi, vnp1->data.ptrvalue);
			else
	  	    		sprintf(temp, "%s", vnp1->data.ptrvalue);
		    }
		    found_id = TRUE;
		    break;
		}
	}

/* If no textid, but a gi, then use that. */
	if (found_id == FALSE && gi != 0)
	{
		sprintf(temp, "gi|%ld", (long) gi);
		length = PrivateGetNumberOfDigits(gi);
		length += 3;	/* 3 extra chars for "gi|" */
		found_id=TRUE;
	}


/* Id's of last resort. */
	if (found_id == FALSE)
	    switch (vnp->choice) 
	    {
		default:
		    	if (id != NULL)
				sprintf(temp, "Unknown");
			length = 7;
			break;
		case BLAST0SeqId_giid:
		    	length = PrivateGetNumberOfDigits(vnp->data.intvalue);
			length += 3;
		    	if (id != NULL)
			{
				sprintf(temp, "gi|%ld", (long) vnp->data.intvalue);
				temp[length] = NULLB;
			}
			break;
		case BLAST0SeqId_textid:
		    	length = StringLen(vnp->data.ptrvalue);
		    	if (id != NULL)
			{
	  	    		sprintf(temp, "%s", vnp->data.ptrvalue);
				temp[length] = NULLB;
			}
		    	break;
	    }

	if (max_id_length > 0)
	{
		if (max_id_length < length)
			length = max_id_length;
	}
	else
	{
		if (BLASTCLI_ID_MAX < length)
			length = BLASTCLI_ID_MAX;
	}

	if (id)
	{
		for (index=0; index<length; index++)
			id[index] = temp[index];
		id[index] = '\0';	
	}

	return length;
}

/***************************************************************************
*	Get the length of an id, truncate if necessary.
***************************************************************************/
static Int2 
get_max_ID_length(BLAST0HitListPtr hlp, Boolean get_gi)

{
	BLAST0SeqDescPtr sdp;
	Int2 length, id_length=0;

	for (; hlp != NULL; hlp=hlp->next)
	{
		if (hlp->seqs != NULL)
		{
			sdp = hlp->seqs->desc;
			if (sdp)
			{
				length = format_an_id(NULL, sdp->id, 0, get_gi);
				if (length > BLASTCLI_ID_MAX)
					id_length = BLASTCLI_ID_MAX;
				if (length > id_length)
					id_length = length;
			}
		}
	}
	return id_length;
}

static Int4 get_pos(BLAST0SeqIntervalPtr loc, Int4 interval)
{

    Int4 position=0;

    position = loc->from + interval + 1; /* plus, default */

    switch (loc->strand)
    {
    	case BLAST0_Seq_interval_strand_minus:
    	case BLAST0_Seq_interval_strand_minus_rf:
        	position = loc->to - interval + 1; /* minus */
		break;
	default:
		break;
    }

    return position;
}

/* 
"get_frame" gets the frame and puts it into an array, which the user
should provide.  The array should have room for three characters.
*/

static void get_frame(BLAST0SeqIntervalPtr loc, CharPtr array, Int4 total_length)
{

    Int2 frame;
    Int4 from;

    switch (loc->strand)
    {
    	case BLAST0_Seq_interval_strand_minus:
    	case BLAST0_Seq_interval_strand_minus_rf:
		from = total_length-(loc->to + 1);
    		frame = (Int2) (from%3) + 1;
		array[0] = '-';
		if (frame == 1)
			array[1] = '1';
		else if (frame == 2)
			array[1] = '2';
		else if (frame == 3)
			array[1] = '3';
		break;
	default:
    		from = loc->from;
    		frame = (Int2) (from%3) + 1;
		array[0] = '+';
		if (frame == 1)
			array[1] = '1';
		else if (frame == 2)
			array[1] = '2';
		else if (frame == 3)
			array[1] = '3';
		break;
    }
    array[2] = '\0';
}

/************************************************************************
*	Set up a SeqPort to allow the user to retrive residues/basepairs
*	one at a time.  Note that the BioseqPtr is really just there
*	to allow the opening of the SeqPortPtr, and is just a skeletal
*	Bioseq.
*********************************************************************/
static SeqPortPtr portnew(BioseqPtr bsp, ValNodePtr str, Int4 len)
{
	Uint1 code;

	BLAST0SeqData2Bioseq(bsp, str, len);
	switch (bsp->mol)
	{
		case Seq_mol_aa:
			code = Seq_code_ncbieaa;
			break;
		case Seq_mol_na: 
			code = Seq_code_iupacna;
			break;
		default:
			return NULL;
	}
	return SeqPortNew(bsp, 0, -1, 0, code);
}

void LIBCALL
TraditionalBlastWarning(BlastReportStructPtr blrp)

{
	BLAST0ResponsePtr blresp, var;
	BLAST0WarningPtr bwp;
	CharPtr warning;

	if (blrp == NULL)
		return;

	blresp = blrp->blresp;

	for (var=blresp; var; var=var->next)
	{
		if (var->choice == BLAST0Response_warning)
		{
			bwp = var->data.ptrvalue;
			warning = bwp->reason;
			ErrPostEx(SEV_WARNING, 0, 0, "%s", warning);
		}
	}
}

typedef struct _residue_symbol_array {
		Uint1Ptr _letters,  /* The actual array that was allocated
					and should be deallocated. */
			letters;    /* Offset from _letter by start of alphabet.*/
		Uint1	first_letter,	/* 1st valid letter in alphabet. */
			last_letter;	/* last valid letter in alphabet. */
        } ResidueSymbolArray, PNTR ResidueSymbolArrayPtr;

static int LIBCALLBACK
sort_letters(VoidPtr vp1, VoidPtr vp2)

{
	Uint1Ptr letter1=vp1, letter2=vp2;

	if (*letter1 < *letter2)
		return -1;
	if (*letter1 > *letter2)
		return 1;
	return 0;
}

/******************************************************************************
*
*	Sorts all the residues for an alphabet into an array so that they
*	can be quickly retrieved.
*
*******************************************************************************/

static ResidueSymbolArrayPtr
ResidueSymbolArrayNew(SeqCodeTablePtr sctp)

{
	ResidueSymbolArrayPtr rsa;
	Int2 index;
	Uint1 first_letter, last_letter;
	Uint1Ptr sctp_letters, letters, _letters;

	if (sctp == NULL)
		return NULL;

	if (! sctp->one_letter) 
		return NULL;

	rsa = (ResidueSymbolArrayPtr) MemNew(sizeof(ResidueSymbolArray));
	if (rsa == NULL)
		return NULL;

	sctp_letters = MemNew((sctp->num + 1)*sizeof(Uint1));

        for (index=0; index < (Int2)sctp->num; index++)
        {
		sctp_letters[index] = sctp->letters[index];
        }

	HeapSort((Nlm_VoidPtr)sctp_letters, (size_t)sctp->num, sizeof(Uint1), sort_letters);
	
	first_letter = (Uint1) sctp_letters[0];
	last_letter = (Uint1) sctp_letters[sctp->num - 1];

	sctp_letters = MemFree(sctp_letters);

	_letters = MemNew((last_letter-first_letter+1)*sizeof(Uint1));

	MemSet(_letters, INVALID_RESIDUE, (last_letter-first_letter+1));

	letters = _letters - first_letter;

        for (index=0; index < (Int2)sctp->num; index++)
	{
		letters[sctp->letters[index]] = (Uint1) index + sctp->start_at;
	}

	rsa->letters = letters;
	rsa->_letters = _letters;
	rsa->first_letter = first_letter;
	rsa->last_letter = last_letter;
	
	return rsa;
}

static ResidueSymbolArrayPtr
ResidueSymbolArrayDestroy(ResidueSymbolArrayPtr rsa)

{
	if (rsa == NULL)
		return NULL;

	rsa->_letters = MemFree(rsa->_letters);
	rsa->letters = NULL;
	rsa = MemFree(rsa);

	return rsa;
}

/***********************************************************************
*
*	Translates letter chh, after checking that a valid translation 
*	exists.  If no valid translation exists, return INVALID_RESIDUE.
*
*	A ResidueSymbolArrayPtr must be initialized first. 
*
*************************************************************************/

static Uint1 ResidueSymbolTranslate(ResidueSymbolArrayPtr rsa, Uint1 chh)

{
	if (rsa == NULL)
		return INVALID_RESIDUE;

	if (chh < rsa->first_letter || chh > rsa->last_letter)
		return INVALID_RESIDUE;

	return rsa->letters[chh];
}

#define BLSTOUT_PRINT_LEN 60

/************************************************************************
*	
*	Print out the alignment portion of the traditional BLAST output,
*	that is:
*
*>dbj|D32197|RURMTCB27 Mitochondrion Rupicapra rupicapra gene for cytochrome b.
*            >emb|D32197|MIRRCB27 Mitochondrion Rupicapra rupicapra gene for
*            cytochrome b.
*            Length = 646
*
*  Plus Strand HSPs:
*
* Score = 3230 (892.5 bits), Expect = 7.0e-264, P = 7.0e-264
* Identities = 646/646 (100%), Positives = 646/646 (100%), Strand = Plus / Plus
*
*Query:     1 AATACACTACACATCCGATACAGCAACAGCATTCTCCTCTGTAACCCACATTTGCCGAGA 60
*             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*Sbjct:     1 AATACACTACACATCCGATACAGCAACAGCATTCTCCTCTGTAACCCACATTTGCCGAGA 60
*
*Query:    61 TGTAAACTACGGCTGAATCATCCGATACATACATGCAAATGGAGCATCAATATTTTTCAT 120
*             ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*Sbjct:    61 TGTAAACTACGGCTGAATCATCCGATACATACATGCAAATGGAGCATCAATATTTTTCAT 120
*
* and so on.
***************************************************************************/
Boolean LIBCALL
TraditionalBottomBlastOutput(BlastReportStructPtr blrp)
{
	SeqPortPtr spps, sppq, spp;
	SeqCodeTablePtr sctp;
	Int2 num_open, num, sum_n;
	static Bioseq bios, bioq;
	BLAST0ResultPtr hdp;
	BLAST0ResponsePtr blresp;
	BLAST0SegmentPtr seg;
	BLAST0SeqIntervalPtr loc, loc1, loc2, loc_next, loc_s;
	BLAST0HitListPtr hlp;
	Int4 nbr_ident, nbr_pos, pos_temp;
	Int2 index1, index2;
	Scores_scaled_intsPtr sco;
   	ValNodePtr   ints, mat_sc;
	BLAST0KABlkPtr ka;
	BLAST0MatrixPtr ma;
    BLAST0SeqDescPtr sdp;
    Char ch, chh, qchh, chh_std, qchh_std;
    Char frame_array[3];
    CharPtr defline, cp, head, s;
    CharPtr program, strand1, strand2;
    ResidueSymbolArrayPtr rsa=NULL;
    ScorePtr score;
    CharPtr scoreDescr;
    BLAST0HSPPtr hsp;
    int i, sc=0, ldef;	/* Left as "int" for use in Warren's blast2_wrap */
    FloatHi pv=0, ev=0;
    Int4 j, pos, range, spos, query_length; 
    Int4 first, width;
    Int4 hsp_length, qposition, sposition;
    Int4 matr[BLASTCLI_MATRIX_SIZE+1][BLASTCLI_MATRIX_SIZE+1];
    Int4 num_of_alignments, max_alignment_number=blrp->num_of_align;
    Uint2 prev_strand;
    double la2;
    FILE *fp;

    /* HTML stuff */

    ValNodePtr vnp1;
    Int4 HTML_gi;
    CharPtr HTML_database;
    Int4 count;
    Boolean use_text_id = FALSE;
    CharPtr HTML_text_id;
    /* ---------- */

    if (blrp == NULL)
	return FALSE;

	hdp = blrp->result;
	query_length = blrp->query_length;
	program = blrp->program;
	fp = blrp->bpsp->fp;
	blresp = blrp->blresp;

/* Extra work is required for alignment of (translated) proteins */
	if (blrp->is_prot)
	{
/* Find the SeqCodeTablePtr to convert for the matrices */
		sctp = SeqCodeTableFind(Seq_code_ncbistdaa);
/* Read in the matrix */
		if ((ma = find(blresp, BLAST0Response_matrix)) != NULL) 
		{
			mat_sc = ma->Scores_scores;
			sco = mat_sc->data.ptrvalue;
			ints = sco->ints;
			if (ints)
			{
			   for (index1=0; index1<BLASTCLI_MATRIX_SIZE; index1++)
			       for (index2=0; index2<BLASTCLI_MATRIX_SIZE; index2++)
			       {
			   	   matr[index1][index2] = ints->data.intvalue;
				   if (ints->next == NULL)
				   {	/* Get out of both loops! */
					index1 = BLASTCLI_MATRIX_SIZE;
					index2 = BLASTCLI_MATRIX_SIZE;
					break;
				   }
				   ints = ints->next;
			       }
			}
		}
		rsa = ResidueSymbolArrayNew(sctp);
	}
	else
	{
		sctp = NULL;
		rsa = NULL;
	}

/* Get lambda over ln(2) */
	if ((ka = find(blresp, BLAST0Response_kablk)) != NULL)
		la2 = ka->lambda / log(2.);  /* lambda / ln2  */
	
	num_of_alignments=0;
/* Start of main loop that runs throughout the function. */
	for (hlp = hdp->hitlists; hlp != NULL; hlp = hlp->next) 	{
		num_of_alignments++;
		if (num_of_alignments > max_alignment_number)
			break;
		if ((hsp = hlp->hsps) == NULL) {
			continue;
		}
		if (hlp->seqs == NULL) {
			continue;
		}
		if ((sdp = hlp->seqs->desc) == NULL) {
			continue;
		}
		ldef = 0;
		for (sdp = hlp->seqs->desc; sdp != NULL; sdp = sdp->next) {
			if (sdp->defline)
			{
				ldef += StringLen(sdp->defline) + 1;
			}
			else
			{
				ldef++;
			}
			ldef += format_an_id(NULL, sdp->id, 0, blrp->get_gi);
			ldef += 2;
		}
		s = defline = MemNew(ldef + 201);
		count = 0;
		for (sdp = hlp->seqs->desc; sdp != NULL; sdp = sdp->next) 
		{
			if (s != defline) {
				sprintf(s++, " ");
			}

                        if(blrp->type == TEXT_OUTPUT)
                          sprintf(s++, "%c", '>');

                        /* HTML stuff */

                        if(blrp->type & HTML_OUTPUT && (count < 1)) {

                          /* HTML database for link to Entrez */
                          
                          if(!StringICmp(program, "blastx") || 
                             !StringICmp(program, "blastp"))
                            HTML_database = "p";
                          else
                            HTML_database = "n";  
                          
                          /* Calculating HTML_gi */
                          
                          use_text_id = FALSE;
                          
                          for (vnp1=sdp->id; vnp1 != NULL; vnp1=vnp1->next) {
                            if (vnp1->choice == BLAST0SeqId_giid) {
                              HTML_gi=sdp->id->data.intvalue;
                            }
                          }
                          if(HTML_gi == 0)
                            use_text_id = TRUE;
                          
                 if(!use_text_id) {
                     fprintf(fp, "\n<a name =%08ld> "
                             "</a><a href=\""
                             "%s/htbin-post/Entrez/query?form=6&dopt=g&db=%s&"
                             "uid=%08ld\"><b>",(long) HTML_gi, 
                             blrp->type & ABSOLUTE_LINKS ? QUERY_HREF : "",
                             
                             HTML_database,  (long) HTML_gi);
                 } else { /* use text id as main reference point */
                     /* Look for text id's */
                     for (vnp1=sdp->id; vnp1 != NULL; vnp1=vnp1->next) {
                         if (vnp1->choice == BLAST0SeqId_textid) {
                             HTML_text_id = StringSave(vnp1->data.ptrvalue);
                             break;
                         }
                     }
                     fprintf(fp, "\n<a name =%s></a> <a href=\""
                             "%s/htbin-post/Entrez/query?form=6&dopt=g&db=%s&"
                             "uid=%s\"><b>", 
                             HTML_text_id, 
                             blrp->type & ABSOLUTE_LINKS ? QUERY_HREF : "",
                             HTML_database,
                             HTML_text_id);
                 }
                        } /* if(blrp->type == HTML_OUTPUT) */

/* Don't pad the id here, as it's secondary and goes into the def line */
			format_an_id(s, sdp->id, 0, blrp->get_gi);
			s += StringLen(s);
                        if (blrp->type & HTML_OUTPUT && (count < 1)) {
                          strcat(s, "</b></a>");
                          s += StringLen(s);
                        }
			if (sdp->defline)
				sprintf(s, " %s", sdp->defline);
			s += StringLen(s);
			count++;
		}
		
		if (defline != NULL) {
			for (cp = defline; ((ch = *cp) != NULLB && !IS_WHITESP(ch)); cp++);
			for (; IS_WHITESP(*cp); cp++) ;
			i = cp - defline;
			i = MIN(i, 12);
			fprintf(fp, "\n\n");
			blast2_wrap(fp, "", defline, ldef, 79, i);
			defline = MemFree(defline);
		}
		fprintf(fp, "%*sLength = %s\n", i, "", 
			Ltostr(hlp->seqs->length, 1));
		prev_strand = UINT2_MAX;	/* Value that will not be used*/
		for (; hsp != NULL; hsp = hsp->next) {
			sum_n = 1;	/* Default (one) not in structure */
			for (score = (ScorePtr) hsp->scores; score != NULL; score = score->next) {
			    if (score->id->str == NULL)
			        scoreDescr = "";
			    else
			        scoreDescr = score->id->str;
			    if (StrCmp(scoreDescr, "score") == 0)
			    {
                                sc = score->value.intvalue;
			    } else if (StrCmp(scoreDescr, "p_value") == 0 ||
			               StrCmp(scoreDescr, "poisson_p") == 0 ||
			               StrCmp(scoreDescr, "sum_p") == 0)
			    {
                                pv = score->value.realvalue;
			    } else if (StrCmp(scoreDescr, "e_value") == 0 ||
			               StrCmp(scoreDescr, "poisson_e") == 0 ||
			               StrCmp(scoreDescr, "sum_e") == 0)
			    {
                                ev = score->value.realvalue;
			    } else if (StrCmp(scoreDescr, "poisson_n") == 0 ||
			               StrCmp(scoreDescr, "sum_n") == 0)
			    {
                                sum_n = score->value.intvalue;
			    } else { /* default */
			    }
			}
			hsp_length = hsp->len;
/* loc1 is the query; loc2 is the subject. */
			loc1 = hsp->segs->loc;
			loc2 = hsp->segs->next->loc;
/* loc and loc_s are set to loc1 or loc2, depending on which "Strand" should
be printed out below. */
			if (StringCmp(program, "blastn") == 0 ||
			    StringCmp(program, "blastp") == 0)
			{
				loc = loc1;
				loc_s = loc2;
			}
			else if (StringCmp(program, "blastx") == 0) 
			{
				loc = loc1;
				loc_s = loc2;
			}
			else if (StringCmp(program, "tblastn") == 0) 
			{
				loc_s = loc1;
				loc = loc2;
			}
			else if (StringCmp(program, "tblastx") == 0) 
			{
				loc_s = loc1;
				loc = loc2;
			}

			if (loc->strand != 0 && loc->strand != prev_strand) 
			{
				prev_strand = loc->strand;
				fprintf(fp, "\n  %s Strand HSPs:\n", 
					gish_str_strand(loc->strand));
			}
			fprintf(fp, "\n Score = %ld (%.1lf bits), Expect = %#0.2lg, ",
				(long)sc, (double) sc * la2, (double) ev);

			if (sum_n < 2) {
				fprintf(fp, "P = %#0.2lg\n", (double) pv);
			} else {
				fprintf(fp, "Sum P(%u) = %#0.2lg\n", sum_n, (double) pv);
			}

			nbr_ident = nbr_pos = 0;
			if ((sppq = portnew(&bioq, hsp->segs->str, hsp_length)) == NULL) {
				return FALSE;
			}
			if ((spps = portnew(&bios, hsp->segs->next->str, hsp_length)) == NULL) {
				return FALSE;
			}
			
			range = loc_s->to - loc_s->from;
			if (StringCmp(program, "tblastx") == 0) 
				range /= 3;
			for (pos=0; pos <= range; pos++) 
			{
				chh = SeqPortGetResidue(spps);
				qchh = SeqPortGetResidue(sppq);
				if (chh == qchh) 
				{
					nbr_ident++;
					if (StringCmp(program, "blastn") == 0) 
						nbr_pos++;
				}
				if (StringCmp(program, "blastn") != 0) 
				{
				     chh_std = ResidueSymbolTranslate(rsa, chh);
				     if (chh_std == INVALID_RESIDUE)
					chh_std = 0;
				     qchh_std = ResidueSymbolTranslate(rsa, qchh);
				     if (qchh_std == INVALID_RESIDUE)
					qchh_std = 0;
				    if (matr[chh_std][qchh_std] > 0) 
					    nbr_pos++;
				}
			}

			fprintf(fp,
			" Identities = %ld/%ld (%ld%%), Positives = %ld/%ld (%ld%%)",
			(long) nbr_ident, (long) pos, (long)((100*nbr_ident)/pos),
			(long) nbr_pos, (long) pos, (long)((100*nbr_pos)/pos) );
			if (StringCmp( program, "blastn") == 0)
			{
			    strand1 = gish_str_strand(loc1->strand);
			    strand2 = gish_str_strand(loc2->strand);
			    if (StringCmp(strand1, "Undefined") != 0 &&
			        StringCmp(strand2, "Undefined") != 0)
				    fprintf(fp, ", Strand = %s / %s", strand1, strand2);
			}
			else if (StringCmp( program, "blastx") == 0)
			{
				loc = hsp->segs->loc;
				get_frame(loc, frame_array, query_length);
				fprintf(fp, ", Frame = %s", frame_array);
			}
			else if (StringCmp( program, "tblastn") == 0)
			{
				loc = hsp->segs->next->loc;
				get_frame(loc, frame_array, hlp->seqs->length);
				fprintf(fp, ", Frame = %s", frame_array);
			}
			else if (StringCmp( program, "tblastx") == 0)
			{
				loc = hsp->segs->loc;
				get_frame(loc, frame_array, query_length);
				fprintf(fp, ", Frame = %s", frame_array);
				loc = hsp->segs->next->loc;
				get_frame(loc, frame_array, hlp->seqs->length);
				fprintf(fp, " / %s", frame_array);
			}
			putc('\n', fp);
			num_open = 1;
/* Print out the actual alignment. */ 
			for (pos=0; pos < hsp_length; pos += BLSTOUT_PRINT_LEN) {
				head = "Query";
				first = TRUE;
				num = 0;
				for (seg = hsp->segs; seg != NULL; seg = seg->next, num++) {
					loc = seg->loc;
					if (first)
					{
					     qposition = get_pos(loc, pos);
					     if (StringCmp(head, "Query") == 0)
					     width = MAX(5, Lwidth(qposition, 1));
					     loc_next = seg->next->loc;
					     sposition = get_pos(loc_next, pos);
					     width = MAX(width, Lwidth(sposition, 1));
					}
					if (!first) {
						SeqPortSeek(sppq, pos, 0);
						if (num_open != num) {
							spps = portnew(&bios, seg->str, hsp_length);
							if (spps == NULL) {
								return FALSE;
							}
							num_open = num;
						}
						SeqPortSeek(spps, pos, 0);
						for (j=0; j<(width+8); j++)
							putc(' ', fp);
/*
						fprintf(fp, "%s", " ");
*/
						for (j=0; j+pos < hsp_length && j < BLSTOUT_PRINT_LEN; j++) {
							 chh = SeqPortGetResidue(spps);
							qchh = SeqPortGetResidue(sppq);
							if (blrp->is_prot) {
								if (chh == qchh) {
									putc(chh, fp);
								} else {
				     					chh_std = ResidueSymbolTranslate(rsa, chh);
				     					if (chh_std == INVALID_RESIDUE)
										chh_std = 0;
				     					qchh_std = ResidueSymbolTranslate(rsa, qchh);
				     					if (qchh_std == INVALID_RESIDUE)
										qchh_std = 0;
									putc((matr[chh_std][qchh_std]>0) ? '+' : ' ', fp);
								}
							} else {
								putc((chh == qchh) ? '|' : ' ', fp);
							}
						}
						spp = spps;
					} else {
						spp = sppq;
					}
					fprintf(fp, "\n");
					SeqPortSeek(spp, pos, 0);
					if (pos+BLSTOUT_PRINT_LEN < hsp_length)
						spos = BLSTOUT_PRINT_LEN - 1;
					else
						spos = hsp_length - pos - 1;
					pos_temp = pos;
				        if (StringCmp( program, "tblastx") == 0)
					{
						pos_temp *= 3;
					}
					else if (StringCmp(head, "Query") == 0)
					{
					     if (StringCmp(program, "blastx") == 0)
					     {
							pos_temp *= 3;
					     }
					}
					else if (StringCmp(head, "Sbjct") == 0)
					{
					     if (StringCmp(program, "tblastn") == 0)
					     {
							pos_temp *= 3;
					     }
					}
					qposition = get_pos(loc, pos_temp);
					fprintf(fp, "%s: %*ld ", head, (int) width, (long) qposition);
					for (j=0; j <= spos; j++) {
						chh = SeqPortGetResidue(spp);
						fputc(chh ,fp);
					}
					pos_temp = j-1+pos;
				        if (StringCmp( program, "tblastx") == 0)
					{
						pos_temp *= 3;
						pos_temp += 2;
					}
					else if (StringCmp(head, "Query") == 0)
					{
					     if (StringCmp(program, "blastx") == 0 ||
				                  StringCmp( program, "tblastx") == 0)
						{
							pos_temp *= 3;
							pos_temp += 2;
						}
					}
					else if (StringCmp(head, "Sbjct") == 0)
					{
					     if (StringCmp(program, "tblastn") == 0)
					     {
							pos_temp *= 3;
							pos_temp += 2;
					     }
					}
					qposition = get_pos(loc, pos_temp);
					fprintf(fp, " %ld\n", (long) qposition);
					head = "Sbjct";
					first = FALSE;
				}
			}
			SeqPortFree(sppq);
			SeqPortFree(spps);
		}
	}
	fprintf(fp, "\n");
	rsa = ResidueSymbolArrayDestroy(rsa);
    return TRUE;
}

/********************************************************************
*	
*	Formats the parameters (e.g., command line options) sent 
*	back by the server and the statistics found at the bottom
*	of the output.
*********************************************************************/
Boolean LIBCALL
TraditionalTailBlastOutput(BlastReportStructPtr blrp)
{
	BLAST0ResponsePtr blresp;
	FILE *fp;
	ValNodePtr stack;

	blresp = blrp->blresp;
	fp = blrp->bpsp->fp;

	if ((stack = find(blresp, BLAST0Response_parms)) != NULL)
			BlastPrintValNodeStack(stack, "Parameters", fp);

	if ((stack = find(blresp, BLAST0Response_stats)) != NULL)
			BlastPrintValNodeStack(stack, "Statistics", fp);

    	return TRUE;
}

Boolean LIBCALL TraditionalBlastOutput(BLAST0ResultPtr hdp, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp) {
  return  TraditionalBlastOutputInternal(hdp, blresp, 
                                         program, fp, 
                                         TEXT_OUTPUT);
}
Boolean LIBCALL TraditionalBlastOutputHTML(BLAST0ResultPtr hdp, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp) {
  return  TraditionalBlastOutputInternal(hdp, blresp, 
                                         program, fp, 
                                         HTML_OUTPUT|ABSOLUTE_LINKS);
}
Boolean LIBCALL TraditionalBlastOutputHTML2(BLAST0ResultPtr hdp, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp) {
  return  TraditionalBlastOutputInternal(hdp, blresp, 
                                         program, fp, 
                                         HTML_OUTPUT);
}

static Boolean LIBCALL TraditionalBlastOutputInternal(BLAST0ResultPtr hdp, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp, Int4 type)
{

	BlastReportStructPtr blrp;
	
	blrp = TraditionalBlastReportSetUp(hdp, blresp, program, fp, type); 

	if (blrp == NULL)
		return FALSE;

	if (PrintTraditionalBlastPreface(blrp) == FALSE)
		return FALSE;

	if (!TraditionalHeadBlastOutput(blrp)) {
		return FALSE;
	}
   
   if (!TraditionalHistBlastOutput(blrp)) {
      return FALSE;
   }
   
  
	if (TraditionalTopBlastOutput(blrp))
	{
		TraditionalBlastWarning(blrp);
		if (!TraditionalBottomBlastOutput(blrp)) {
			return FALSE;
		}
    }
    if (!TraditionalTailBlastOutput(blrp)) {
    	return FALSE;
    }

    blrp = TraditionalBlastReportCleanUp(blrp);

    return TRUE;
}
