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
* File Name: blast2.c
*
* Author: Tom Madden
*
* Version Creation Date:   10/26/95
*
* $Revision: 6.1 $
*
* File Description: 
*       Functions to format parts of the traditional BLAST output. 
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
* $Log: blast2.c,v $
* Revision 6.1  1997/11/28 18:21:57  madden
* fprintf fixes
*
* Revision 6.0  1997/08/25 18:33:53  madden
* Revision changed to 6.0
*
* Revision 5.2  1996/12/13 18:10:58  madden
* Added TraditionalBlastReportAddResult.
*
 * Revision 5.1  1996/11/27  14:01:16  madden
 * Replaced strlen by Nlm_StringLen.
 *
 * Revision 5.0  1996/05/28  14:09:11  ostell
 * Set to revision 5.0
 *
 * Revision 1.12  1996/05/22  18:04:17  madden
 * Fixed formatting of Query ID if only gi is present.
 *
 * Revision 1.11  1996/03/25  17:30:33  shavirin
 * Changes due to new function TraditionalBlastOutputHTML(),
 * that produce HTML output from the BLAST server
 *
 * Revision 1.10  1996/03/11  22:02:09  madden
 * Added function BLAST0SeqData2Bioseq.
 *
 * Revision 1.9  1995/11/22  22:43:08  madden
 * Added return value for PrintTraditionalBlastPreface.
 *
 * Revision 1.8  1995/11/07  17:56:51  madden
 * Added PrintTraditionalBlastPreface and BlastPrintCheckBufferStatus.
 *
 * Revision 1.7  1995/11/03  16:02:43  madden
 * Added BlastPrintNewLine, renames BlastPrintInit to BlastPrintStructInit.
 *
 * Revision 1.6  1995/11/02  23:14:49  madden
 * Fixed bugs with BlastPrintTabToColumn.
 *
 * Revision 1.5  1995/11/02  21:46:44  madden
 * Added functions BlastPrintDoubleAdd and BlastPrintIntegerAdd.
 *
 * Revision 1.4  1995/11/02  21:12:25  madden
 * Added "BlastPrint" functions to print out lines.
 *
 * Revision 1.3  1995/11/01  15:40:18  madden
 * Addition of TraditionalBlastReportSetUp and TraditionalBlastReportCleanUp.
 *
 * Revision 1.2  1995/10/27  21:17:10  madden
 * Added TraditionalHistBlastOutput.
 *
 * Revision 1.1  1995/10/26  17:06:44  madden
 * Initial revision
 *
*
*/

#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#define NLM_GENERATED_CODE_PROTO
#include <blast2.h>
#include <objblst2.h>


/*
	This function is used, by TraditionalBlastReportSetUp, to find
	a node specified by the Uint1 v.
*/

static Pointer find(BLAST0ResponsePtr p, Uint1 v)
{
	BLAST0ResponsePtr b;

	for (b = p; b != NULL && b->choice != v; b = b->next) 
		;

	return (b == NULL) ? NULL : b->data.ptrvalue;
}

/*
	Adds the BLAST0ResultPtr result, BLAST0ResponsePtr blresp
	to a BlastReportStructPtr that has already been allocated
	and (possibly) used to print out the BLAST header before
	the full results were available.

	The BlastReportStructPtr should be allocated with the
	function TraditionalBlastReportSetUp first.
*/

Boolean LIBCALL
TraditionalBlastReportAddResult(BlastReportStructPtr blrp, BLAST0ResultPtr result, BLAST0ResponsePtr blresp)

{
	if (blrp == NULL)
		return FALSE;

	blrp->result = result;
	blrp->blresp = blresp;

	return TRUE;
}

/*****************************************************************************
*
*	Allocates and fills in the BlastReportStruct.  This function should
*	be called before any of the other "traditional" formatter functions
*	are called.
****************************************************************************/ 
BlastReportStructPtr
TraditionalBlastReportSetUp(BLAST0ResultPtr result, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp, Int4 type)

{
	Boolean get_gi=FALSE;
	BlastReportStructPtr blrp;
	BLAST0SequencePtr query;
	CharPtr string;
	Int4 offset, num_deflines, num_of_aligns;
	ValNodePtr parms;

	if ((blrp=MemNew(sizeof(BlastReportStruct))) == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, 
		  "Couldn't allocate memory in SetUpForTraditionalBlastReport");
		return NULL;
	}
	
	blrp->result = result;
	blrp->blresp = blresp;
	blrp->program = StringSave(program);
	blrp->bpsp = BlastPrintStructInit(BLAST2_LINE_LENGTH, fp);
	
/* Set is_prot if alignments done with proteins (regardless of query type) */
	if (StringCmp(program, "blastn") != 0)
		blrp->is_prot = TRUE;
	else
		blrp->is_prot = FALSE;

	if ((query = find(blresp, BLAST0Response_query)) != NULL)
		blrp->query_length = query->length;
	else
		blrp->query_length = 0;

	offset = 0;
	num_deflines = 500;
	num_of_aligns = 250;
	if ((parms = find(blresp, BLAST0Response_parms)) != NULL)
	{
		while (parms != NULL)
		{
			string = parms->data.ptrvalue;
			if (StringNCmp(string, "-qoffset", 8) == 0)
				sscanf(string+9, "%ld", &offset);
			else if (StringNCmp(string, "V=", 2) == 0)
				sscanf(string+2, "%ld", &num_deflines);
			else if (StringNCmp(string, "B=", 2) == 0)
				sscanf(string+2, "%ld", &num_of_aligns);
			else if (StringNCmp(string, "-gi", 3) == 0)
				get_gi=TRUE;
			parms = parms->next;
		}
	}

	blrp->qoffset = offset;
	blrp->num_of_defline = num_deflines;
	blrp->num_of_align = num_of_aligns;
	blrp->get_gi = get_gi;
	blrp->type = type; /*   TEXT_OUTPUT or HTML_OUTPUT */

	return blrp;
}

/**************************************************************************
*
*	Deallocate the BlastReportStructPtr and all allocated memory
*	within it.
************************************************************************/

BlastReportStructPtr
TraditionalBlastReportCleanUp(BlastReportStructPtr blrp)

{
	if (blrp == NULL)
		return NULL;

	blrp->program = MemFree(blrp->program);
	blrp->bpsp = BlastPrintStructClose(blrp->bpsp);

	blrp = MemFree(blrp);

	return blrp;
}

/***************************************************************************
*
*
*	Format the "preface" of the BLAST output.  This includes the
*	program, version, compile date: 
*
*BLASTN 1.4.7MP [16-Oct-94] [Build 16:01:48 Dec 16 1994]
*
*Reference:  Altschul, Stephen F., Warren Gish, Webb Miller, Eugene W. Myers,
*and David J. Lipman (1990).  Basic local alignment search tool.  J. Mol. Biol.
*215:403-10.
*
*Notice:  this program and its default parameter settings are optimized to find
*nearly identical sequences rapidly.  To identify weak similarities encoded in
*nucleic acid, use BLASTX, TBLASTN or TBLASTX.
*
**************************************************************************/

Boolean LIBCALL
PrintTraditionalBlastPreface(BlastReportStructPtr blrp)
{
	BlastPrintStructurePtr bpsp;
	BLAST0ResponsePtr blresp; 
	BLAST0PrefacePtr preface;
	Char temp[BLAST2_LINE_LENGTH];
	CharPtr string;
	ValNodePtr node;

	if (blrp == NULL)
		return FALSE;

	blresp = blrp->blresp; 
	bpsp = blrp->bpsp;

	preface = find(blresp, BLAST0Response_preface);

	if (preface == NULL)
		return FALSE;

	BlastPrintStart(bpsp, 0, 0);
        if(blrp->type == HTML_OUTPUT)
          sprintf(temp, "<PRE><b>%s %s [%s]", preface->program, preface->version, preface->dev_date);
        else  /* blrp->type == TEXT_OUTPUT */
          sprintf(temp, "%s %s [%s]", preface->program, preface->version, preface->dev_date);

	BlastPrintStringAdd(bpsp, temp);

	if (preface->bld_date != NULL) 
	{
          if(blrp->type == HTML_OUTPUT) 
            sprintf(temp, " [Build %s]</b>", preface->bld_date);
          else  /* blrp->type == TEXT_OUTPUT */
            sprintf(temp, " [Build %s]", preface->bld_date);
          
          BlastPrintStringAdd(bpsp, temp);
	}
	BlastPrintNewLine(bpsp);
	BlastPrintNewLine(bpsp);

	for (node = preface->cit; node != NULL; node = node->next) 
	{
		string = (CharPtr) node->data.ptrvalue;
		if (node == preface->cit)
		{
                  if(blrp->type == HTML_OUTPUT) 
			BlastPrintStringAdd(bpsp, "<b>Reference:</b>  ");
                  else  /* blrp->type == TEXT_OUTPUT */
			BlastPrintStringAdd(bpsp, "Reference:  ");
		}
		BlastPrintStringAdd(bpsp, string);
		BlastPrintNewLine(bpsp);
	}

	for (node = preface->notice; node != NULL; node = node->next) 
	{
		BlastPrintNewLine(bpsp);
		string = (CharPtr) node->data.ptrvalue;
                if(blrp->type == HTML_OUTPUT) 
                  BlastPrintStringAdd(bpsp, "<b>Notice:</b>  ");
                else  /* blrp->type == TEXT_OUTPUT */
                  BlastPrintStringAdd(bpsp, "Notice:  ");

		BlastPrintStringAdd(bpsp, string);
	}
	BlastPrintEnd(bpsp);

	return TRUE;
}

/*************************************************************************
*
*	This function acknowledges the BLAST request.  The information
*	comes in the BLAST0SequencePtr, which is then formatted as something
*	of the form:
*
*	Query=  195747|Mitochondrion Rupicapra rupi
*        	(646 letters)
*
*	The query length is returned.
*	
*	This code was orignally part of TraditionalHeadBlastOutput in
*	blstout1.c and is still called in TraditionalHeadBlastOutput.
*	
*	"strands" tells how which strands are being translated (set by 
*	"-top" and "-bottom").
*	"offset" reports the offset of the search (set by "-qoffset" option).
*	
*
*************************************************************************/

#define BLASTCOM_BUFFER_SIZE 50

Int4 LIBCALL
acknowledge_blast_request(CharPtr program, BLAST0SequencePtr query, Int2 strands, Int4 offset, FILE *fp, Int4 type)

{

	Char buf[BLASTCOM_BUFFER_SIZE];
	CharPtr string;
	Int4 query_length=0;
	ValNodePtr node;

        if(type == HTML_OUTPUT) 
          fprintf(fp, "<b>\nQuery=</b>  ");
        else  /* type == TEXT_OUTPUT */
          fprintf(fp, "\nQuery=  ");
	if (query != NULL && query->desc != NULL)
	{
	    node = query->desc->id;
	    while (node)
	    {
		if (node->choice == BLAST0SeqId_giid)
		{
			if (node->next)
				sprintf(buf, "gi|%ld|", (long) node->data.intvalue);
			else
				sprintf(buf, "gi|%ld", (long) node->data.intvalue);
			string = &buf[0];
		}
		else if (node->choice == BLAST0SeqId_textid)
		{
	    		string = node->data.ptrvalue;
	    		if (strncmp(string, "lcl|", 4) == 0)
			{
				StringNCpy(buf, string+4, BLASTCOM_BUFFER_SIZE);
				string = &buf[0];
			}
		}
	    	fprintf(fp, "%s", string);
		node = node->next;
	    }
	    if((string=query->desc->defline) != NULL && *string != NULLB)
	    	blast2_wrap(fp, " ", string, Nlm_StringLen(string), 79, 8);
	    else
		fprintf(fp, "\n");
	} 
	else 
	{
		fprintf(fp, "Unknown\n");
	}

	if (query)
	{
		query_length = query->length;
		fprintf(fp, "        (%s letters", 
			Nlm_Ultostr(query->length, 1));
	}

	if (offset > 0)
		fprintf(fp, ", %s offset", Ltostr(offset,1));

	fprintf(fp, ")\n");

	if (StringCmp(program, "blastx") == 0 ||
		StringCmp(program, "tblastx") == 0)
	{
		fprintf(fp, "\n  Translating ");
				switch (strands) {
			case 0:
				break;
			case 1:
				fprintf(fp, "top strand of query sequence in 3 reading frames\n");
				break;
			case 2:
				fprintf(fp, "bottom strand of query sequence in 3 reading frames\n");
				break;
			case 3:
				fprintf(fp, "both strands of query sequence in all 6 reading frames\n");
				break;
			default:
				break;
			}
	}

	return query_length;
}

/*************************************************************************
*
*	Produces an NCBI SeqId from a BLAST0 SeqId.
*	
**************************************************************************/

static SeqIdPtr 
BLAST0_ID2NCBI_ID(ValNodePtr vnp)

{
	ValNodePtr seqid=NULL, text_seqid;

	while (vnp)
	{
		if (vnp->choice == BLAST0SeqId_giid)
		{
			ValNodeAddInt(&seqid, SEQID_GI, vnp->data.intvalue);
		}
		else if (vnp->choice == BLAST0SeqId_textid)
		{
			text_seqid = SeqIdParse(vnp->data.ptrvalue);
			ValNodeLink(&seqid, text_seqid);
		}
	}
	return seqid;
}

/***************************************************************************
*
*	Produces a BioseqPtr from the data in BLAST0SequencePtr.
*
****************************************************************************/

BioseqPtr LIBCALL
BLAST0Sequence2Bioseq(BLAST0SequencePtr sequence)

{

	BioseqPtr bsp;

	bsp = BioseqNew();
	if (bsp == NULL)
		return NULL;

	BLAST0SeqData2Bioseq(bsp, sequence->seq, sequence->length);
	bsp->id = BLAST0_ID2NCBI_ID(sequence->desc->id);
	if (sequence->desc && sequence->desc->defline)
		ValNodeCopyStr(&(bsp->descr), Seq_descr_title, sequence->desc->defline);

	return bsp;
}


/************************************************************************
*
*	Transfers the sequence data from the BLAST0SequencePtr sequence
*	(in the ValNodePtr str) to the BioseqPtr bsp.  Note that the 
*	BioseqPtr is not allocated here.
*
************************************************************************/
	
Int2 LIBCALL
BLAST0SeqData2Bioseq(BioseqPtr bsp, ValNodePtr sequence, Int4 length)

{
	if (bsp == NULL || sequence == NULL)
		return 1;

	bsp->length = length;
	switch (sequence->choice) 
	{
		case BLAST0SeqData_ncbistdaa: 
			bsp->seq_data_type = Seq_code_ncbistdaa;
			bsp->mol = Seq_mol_aa;
			break;
		case BLAST0SeqData_ncbi4na: 
			bsp->seq_data_type = Seq_code_ncbi4na;
			bsp->mol = Seq_mol_na;
			break;
		default:
			return 1;
	}
	bsp->seq_data = sequence->data.ptrvalue;
	bsp->repr = Seq_repr_raw;

	return 0;
}


/*************************************************************************
*
*	Print out the strings stored in the ValNodePtr stack.  "title"
*	is printed first.
*************************************************************************/

Int2 LIBCALL
BlastPrintValNodeStack(ValNodePtr stack, CharPtr title, FILE *fp)

{
	Boolean new;
	CharPtr	string;
	ValNodePtr vnp;

	if (fp == NULL || stack == NULL)
		return 0;

	if (stack != NULL) 
	{
	     if (title != NULL)
                     fprintf(fp, "\n%s:\n", title);
             new = TRUE;
             for (vnp=stack; vnp != NULL; vnp = vnp->next) 
	     {
                     string = vnp->data.ptrvalue;
                     if (*string == NULLB) 
		     {
                          fprintf(fp, "\n");
                           new = TRUE;
                     } 
		     else 
		     {
                          if (new)
                                  fprintf(fp, "  ");
                           fprintf(fp, "%s", string);
                           new = FALSE;
                     }
               }
        }
	return 0;
}


static void
doindent(FILE *fp, int ncols)
{
      while (ncols-- > 0)
              fputc(' ', fp);
}

/*************************************************************
*                                                             *
*  Warren Gish's utilities needed for producing BLAST output  *
*                                                             *
 *************************************************************/

/* blast2_wrap -- wordwrap lines of output */

void
blast2_wrap(FILE *fp, /*  the output stream */
	CharPtr	title,	/* title string */
	CharPtr	s,	/* pointer to null-terminated string for output */
	int	slen,	/* strlen(s), or -1 */
	int	linelen, /* max. length of output line before each '\n' */
	int	indent	/* no. of columns to indent any continuation lines */
	)
{
	register char	*savep, *savep2;
	register char	*cp, ch;
	int		outlen, len, olinelen;
	int		titlelen;
	Boolean	once = TRUE;

	if (title != NULL) {
		cp = title;
		len = 0;
		while ((ch = *cp++) != NULLB) {
			fputc(ch, fp);
			if (ch == '\n')
				len = cp - title;
		}
		titlelen = (cp - title) - len - 1;
	}
	else
		titlelen = 0;

	if (slen < 0)
		slen = Nlm_StringLen(s);

	if (indent >= linelen) {
		indent = MIN(1 + titlelen, linelen-5);
		indent = MAX(indent, 0);
	}
	olinelen = linelen;
	linelen -= titlelen;
	for (cp = s; cp < s + slen;) {
		/* Skip leading white space */
		for (;; ++cp) {
			if ((ch = *cp) == NULLB)
				return;
			if (!IS_WHITESP(ch))
				break;
		}

		outlen = cp - s;
		if (slen - outlen <= linelen) {
			/* Remainder is short enough to fit on one line */
			if (!once)
				doindent(fp, indent);
			fprintf(fp, "%s\n", cp);
			break;
		}
		else {
			savep2 = cp + linelen;
			ch = *savep2;
			if (IS_WHITESP(ch)) {
Phase2:
				for (savep = savep2; savep >= cp; --savep)
					if (!IS_WHITESP(*savep)) {
						++savep;
						break;
					}
			}
			else {
				for (savep = savep2; savep >= cp; --savep)
					if (IS_WHITESP(*savep)) {
						savep2 = savep;
						goto Phase2;
					}
				/* a _very_ long word here! */
				savep = savep2;
			}
		}
		if (!once)
			doindent(fp, indent);
		once = FALSE;
		while (cp < savep)
			fputc(*cp++, fp);
		putc('\n', fp);
		cp = savep2;
		linelen = olinelen - indent;
	}

}

/* 
	This function is used by TraditionalHistBlastOutput
	(below).
*/
static CharPtr
print_double(CharPtr pCh, double dX, int iWidth, int iPrecision)
{
   double dY;
   int    i;
   
   dY = log(dX * 1.000000001) / NCBIMATH_LN10;
   i = dY;
   dY = Nlm_Powi(10., i - iPrecision + 1);
   dX = Nlm_Nint(dX / dY) * dY;
   dY = Nlm_Powi(10., iPrecision - 1) - 1.e-7;
   if (dX >= dY)
      sprintf(pCh, "%*.0lf", iWidth, dX);
   else
      sprintf(pCh, "%*.*lf", iWidth, iPrecision - i - 1, dX);
   
   return pCh;
}

#define PAGE_W 80
#define EXPECT_W 6
#define EXPECT_PRCSN 3
#define SEPARATOR_SYMBOL '|'
#define DRAWING_SYMBOL '='
#define SMALLDRAWING_SYMBOL ':'
/**************************************************************************
 *    The function formats histgram data and is part of the "traditional"
 * BLAST output.
 * 
 * 
 **************************************************************************/
Boolean LIBCALL
TraditionalHistBlastOutput(BlastReportStructPtr blrp)
{
   BLAST0HistogramPtr pHist;
   BLAST0ResultPtr result;
   Int4 nBars;     /* Number of histogram bars */
   Int4 i;
   Int4 nMaxObserved;
   Int4 nMaxHist;
   Int4 nHistValue;
   Int2 nObsWidth; /* Width (in symbols) of field wich represents observed value */
   Int2 nHistWidth;/* Width (in symbols) of field wich represents histogram value */
   Int2 nCols;     /* Width of space (in symbols) where histogram is painted */
   Int2 nPerCol;   /* Number of sequences per unit */
   Int2 nUnits;    /* Number units needed to represent histogram bar */
   Int2 j;
   FloatLo fDelta;
   Char aCh[20];
   Boolean bNeedCheck = TRUE;
   BLAST0HistogramBarPtr pBar;
   FILE *fp;
      
   
	if (blrp == NULL)
		return FALSE;

	fp = blrp->bpsp->fp;
	result = blrp->result;
	pHist = result->hist;

   if (pHist == NULL)
      return TRUE;
   
   
   putc('\n', fp);
   putc('\n', fp);
   fprintf(fp, "     Observed Numbers of Database Sequences Satisfying\n");
   fprintf(fp, "    Various EXPECTation Thresholds (E parameter values)\n");
   putc('\n', fp);

   nBars  = pHist->nbars;
      
   /* Find out max observed and histogram data values
    */
   for (i = 0, nMaxObserved = 0, nMaxHist = 0, pBar = pHist->bar;
	i < nBars; i++) {
      nMaxObserved = MAX(nMaxObserved, pBar->n);
       
      if (i == nBars - 1) /* the last one */
	 nMaxHist = MAX(nMaxHist, pBar->n - pHist->base);
      else
	 nMaxHist = MAX(nMaxHist, pBar->n - pBar->next->n);
       
      pBar = pBar->next;
   }
   
   nObsWidth = Nlm_Ulwidth(nMaxObserved, 0);
   nObsWidth = MAX(nObsWidth, 2);
   nHistWidth = Nlm_Ulwidth(nMaxHist, 0);
   nHistWidth = MAX(nHistWidth, 2);
   
/* Gets the number of hits per "=" */
   nCols = PAGE_W - (EXPECT_W+1 + nObsWidth+1 + nHistWidth+1 + 2);
   fDelta = (float)nMaxHist / (float)nCols;
   nPerCol = ceil(fDelta);
   nPerCol = MAX(nPerCol, 1);
   
   fprintf(fp, "        Histogram units:      %c %d Sequence%s",
	   DRAWING_SYMBOL, nPerCol, (nPerCol == 1 ? "" : "s"));

    if (fDelta > 1.0)
      fprintf(fp, "     %c less than %d sequences\n", SMALLDRAWING_SYMBOL, nPerCol);
    else
      putc('\n', fp);
    
   putc('\n', fp);
   fprintf(fp, " EXPECTation Threshold\n");
   fprintf(fp, " (E parameter)\n");
   fprintf(fp, "    |\n");
   fprintf(fp, "    V   Observed Counts-->\n");
   
/* If there were no hits, exit. */
   if (nMaxHist == 0) 
   {
      fprintf(fp, "\n\n***  histogram slots are all empty ***\n\n");
      return TRUE;
   }

   /* Draw the histogram
    */
   for (i = 0, pBar = pHist->bar; i < nBars; i++) 
   {
      if ( bNeedCheck && (pBar->x <= pHist->expect + 1e-6)) {
	 /* This line must appear only once
	  */
	 fprintf(fp,
		 " >>>>>>>>>>>>>>>>>>>>>  Expect = %#0.3lg, Observed = %lu  <<<<<<<<<<<<<<<<<\n",
		 pHist->expect, pHist->observed);
	 bNeedCheck = FALSE;
      }

      nHistValue = (i < nBars - 1) ? (pBar->n - pBar->next->n) : (pBar->n - pHist->base);
      fprintf(fp, " %s %*lu %*lu %c", print_double(aCh, pBar->x, 6, 3),
	      nObsWidth, pBar->n, nHistWidth, nHistValue, SEPARATOR_SYMBOL);

      /* Draw long line ("=")
	nHistValue is the number of hits observed, nPerCol is the number of
	hits that each symbol represents. 
       */
	nUnits = (Int2) (nHistValue/nPerCol);
	for (j = 0; j < nUnits; j++)
	    putc(DRAWING_SYMBOL, fp);
          
      /* Draw tiny line if nothing was drawn above (":").
       */
      if (j == 0 && nHistValue > 0 && nPerCol > 1)
	putc(SMALLDRAWING_SYMBOL, fp);
      
      putc('\n', fp);
      
      pBar = pBar->next;
   }
   
   return TRUE;
}

/************************************************************************

	The following functions are used to print the output to a file:

	BlastPrintStructInit: call once at beginning of program, performs
		initialization.
	BlastPrintStart: call when a newline should be started or the
		initial or continuation line indentation changes.
	BlastPrintAddString: Adds String to the print "buffer".
	BlastPrintAddChar: Adds a Char to the print "buffer".
	BlastPrintEnd: finishes section, newline inserted.
	BlastPrintStructClose: call once at end of program, does deallocation.

*************************************************************************/

/*
	Call this function ONCE. If failure, NULL is returned,
	otherwise the BlastPrintStructurePtr is returned. 
*/

BlastPrintStructurePtr
BlastPrintStructInit(Int2 line_length, FILE *fp)

{
	BlastPrintStructurePtr bpsp;
	CharPtr buffer;

	bpsp = (BlastPrintStructurePtr) MemNew(sizeof(BlastPrintStructure));

	if (bpsp == NULL)
	{
		return NULL;
	}
/* Make the buffer 2 longer for safety */
	buffer = (CharPtr) MemNew((line_length+2)*sizeof(Char));
	if (buffer == NULL)
	{
		bpsp = MemFree(bpsp);
		return NULL;
	}

	bpsp->fp = fp;
	bpsp->buffer = buffer;
	bpsp->line_length = line_length;
	bpsp->position = 0;

	return bpsp;
}


/*
	Call at the start of every "paragraph".  TRUE is returned on
	success, otherwise FALSE.
*/
Boolean
BlastPrintStart(BlastPrintStructurePtr bpsp, Int2 init_indent, Int2 cont_indent)

{
	if (bpsp == NULL)
		return FALSE;

	bpsp->init_indent = init_indent;
	bpsp->cont_indent = cont_indent;
	bpsp->position = 0;	/* set position of buffer to zero. */
	bpsp->first_line = TRUE;
	bpsp->buffer[0] = NULLB;

	return TRUE;
}

/*
	Adds a string to the buffer; if buffer is full, print buffer,
	add line return, and reset position to zero.

	The total number of characters printed is returned.
*/

Int2
BlastPrintStringAdd(BlastPrintStructurePtr bpsp, CharPtr string)

{
	CharPtr buffer;
	Int2 indent, line_length, num_of_chars=0, position;

	if (bpsp == NULL)
		return 0;

	buffer = bpsp->buffer;
	position = bpsp->position;
	buffer += position;
	line_length = bpsp->line_length;

/* if position is zero we're starting a new line. */
	if (position == 0)
	{
		if (bpsp->first_line == TRUE)
		{
			indent = bpsp->init_indent;
			bpsp->first_line = FALSE;
		}
		else
		{
			indent = bpsp->cont_indent;
		}	

		while (indent > 0)
		{
			*buffer = ' ';
			buffer++;
			position++;
		}
	}
	
	while (*string != NULLB)
	{
		*buffer = *string;
		string++;
		buffer++;
		position++;
		if (position >= line_length && *string != NULLB)
		{
			bpsp->position = position;
			num_of_chars += BlastPrintCheckBufferStatus(bpsp, *(string));
			if (*string == ' ')
				string++;
			num_of_chars += BlastPrintStringAdd(bpsp, string);
			break;
		}
	}


/* If num_of_chars is non-zero, then position was set. */
	if (num_of_chars == 0)
	{
		bpsp->position = position;
		return position;
	}
	else
	{
		return num_of_chars;
	}
}

#define NUM_TO_BACKTRACK 10	/* How far to look for whitespace. */
Int2
BlastPrintCheckBufferStatus (BlastPrintStructurePtr bpsp, Char next_char)

{
	Boolean found_whitespace=FALSE;
	CharPtr buffer, ptr;
	Char temp[NUM_TO_BACKTRACK];
	Int2 index, num_of_chars=0;

	if (next_char != NULLB && next_char != ' ')
	{
		buffer = bpsp->buffer;

		buffer += bpsp->position;

		/* Look for a whitespace to break on. */
		for (index=0; index<NUM_TO_BACKTRACK; index++)
		{
			if (*buffer == ' ')
			{
				found_whitespace=TRUE;
				break;
			}
			buffer--;
		}

		if (found_whitespace == TRUE)
		{
			*buffer = NULLB;
			buffer++;
			bpsp->position -= index;

			ptr = &temp[0];
			while (index > 0)
			{
				*ptr = *buffer;
				ptr++; buffer++;
				index--;
			}
			*ptr = NULLB;
		}
		num_of_chars = bpsp->position;
		BlastPrintFlush(bpsp);
		buffer = bpsp->buffer;
		ptr = temp;
		index=0;
		while (*ptr != NULLB)
		{
			*buffer = *ptr;
			buffer++; ptr++;
			index++;
		}
		bpsp->position = index;
	}
	else
	{
		num_of_chars = bpsp->position;
		BlastPrintFlush(bpsp);
	}

	return num_of_chars;

}

Int2
BlastPrintDoubleAdd(BlastPrintStructurePtr bpsp, CharPtr format, double number)

{
	Char temp[BLAST2_LINE_LENGTH];

	sprintf(temp, format, (double) number);
	number = BlastPrintStringAdd(bpsp, temp);
	return number;
}

Int2
BlastPrintIntegerAdd(BlastPrintStructurePtr bpsp, CharPtr format, Int4 number)

{
	Char temp[BLAST2_LINE_LENGTH];

	sprintf(temp, format, (long) number);
	number = BlastPrintStringAdd(bpsp, temp);
	return number;
}


Int2
BlastPrintCharAdd(BlastPrintStructurePtr bpsp, Char character)

{
	CharPtr buffer;
	Int2 indent, line_length, num_of_chars=0, position;

	if (bpsp == NULL)
		return 0;

	buffer = bpsp->buffer;
	position = bpsp->position;
	buffer += position;
	line_length = bpsp->line_length;

/* if position is zero we're starting a new line. */
	if (position == 0)
	{
		if (bpsp->first_line == TRUE)
		{
			indent = bpsp->init_indent;
			bpsp->first_line = FALSE;
		}
		else
		{
			indent = bpsp->cont_indent;
		}	

		while (indent > 0)
		{
			*buffer = ' ';
			buffer++;
			position++;
		}
	}

	if (position > line_length)
	{
		bpsp->position = position;
		BlastPrintFlush(bpsp);
		position = bpsp->position;
		buffer = bpsp->buffer;
	}

	position++;
	bpsp->position = position;
	*buffer = character;

	return 1;
}

/*
	Tab to column given by column.  If column is less than position
	a negative number is returned, otherwise the number of blank
	spaces added is returned.

*/
Int2 
BlastPrintTabToColumn(BlastPrintStructurePtr bpsp, Int2 column)

{
	CharPtr buffer;
	Int2 diff, position;

	if (bpsp == NULL)
		return -1;

	column -= 1;	/* Change column so it's zero offset, like position. */

	if (column > bpsp->line_length)
		return -1;

	position = bpsp->position;
	diff = column - position; 
	if (diff <= 0)
		return diff;

	buffer = bpsp->buffer;
	buffer += position;

	while (diff > 0)
	{
		*buffer = ' ';
		buffer++;
		diff--;
	}

	diff = column - position; 
	bpsp->position = column;

	return diff;
}

/* 
	Starts a new-line.
*/
void
BlastPrintNewLine(BlastPrintStructurePtr bpsp)

{
	BlastPrintFlush(bpsp);
	bpsp->position = 0;
	bpsp->buffer[0] = NULLB;
}


/*
	Ends printing.
*/
void
BlastPrintEnd(BlastPrintStructurePtr bpsp)

{
	BlastPrintFlush(bpsp);
	bpsp->position = 0;
	bpsp->buffer[0] = NULLB;
}

/*
	Prints the buffer to the File*.
	This funciton could also store the buffer in a
	ByteStorePtr for later use.
*/
void
BlastPrintFlush(BlastPrintStructurePtr bpsp)

{
	bpsp->buffer[bpsp->position] = NULLB;
	fprintf(bpsp->fp, "%s\n", bpsp->buffer);
	bpsp->position = 0;
}

/*
	Deallocates the CharPtr buffer and the BlastPrintStructurePtr.
	This function does NOT close the file for FILE* fp!
*/
BlastPrintStructurePtr
BlastPrintStructClose(BlastPrintStructurePtr bpsp)

{
	bpsp->buffer = MemFree(bpsp->buffer);
	bpsp = MemFree(bpsp);
	return bpsp;
}
