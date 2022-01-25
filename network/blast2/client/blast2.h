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
* File Name: blast2.h
*
* Author: Tom Madden
*
* Version Creation Date:   10/26/95
*
* $Revision: 6.0 $
*
* File Description: 
*       Functions that format traditional BLAST output.
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
* $Log: blast2.h,v $
* Revision 6.0  1997/08/25 18:33:55  madden
* Revision changed to 6.0
*
* Revision 5.2  1997/05/14 14:19:15  shavirin
* Added define for absolute links
*
 * Revision 5.1  1996/12/13  18:10:58  madden
 * Added TraditionalBlastReportAddResult.
 *
 * Revision 5.0  1996/05/28  14:09:11  ostell
 * Set to revision 5.0
 *
 * Revision 1.11  1996/03/25  17:28:26  shavirin
 * Changes due to new function TraditionalBlastOutputHTML(),
 * that produce HTML output from the BLAST server
 *
 * Revision 1.10  1996/03/22  17:28:49  madden
 * Moved formatting prototypes from netblap2.h
 *
 * Revision 1.9  1996/03/12  19:10:32  madden
 * Added prototype fro BLAST0Sequence2Bioseq.
 *
 * Revision 1.8  1996/03/11  22:02:09  madden
 * Added function BLAST0SeqData2Bioseq.
 *
 * Revision 1.7  1995/11/07  17:55:47  madden
 * Added prototypes for PrintTraditionalBlastPreface and BlastPrintCheckBufferStatus.
 *
 * Revision 1.6  1995/11/03  16:02:43  madden
 * Added BlastPrintNewLine, renames BlastPrintInit to BlastPrintStructInit.
 *
 * Revision 1.5  1995/11/02  21:46:44  madden
 * Added functions BlastPrintDoubleAdd and BlastPrintIntegerAdd.
 *
 * Revision 1.4  1995/11/02  21:12:25  madden
 * Added prototypes for "BlastPrint" functions and typedef for blast_print_struct.
 *
 * Revision 1.3  1995/11/01  15:40:18  madden
 * Addition of TraditionalBlastReportSetUp and TraditionalBlastReportCleanUp.
 *
 * Revision 1.2  1995/10/27  21:17:10  madden
 * Added prototype for TraditionalHistBlastOutput and struct blast_report_struct.
 *
 * Revision 1.1  1995/10/26  17:06:44  madden
 * Initial revision
 *
 * 
*/
#ifndef _BLAST2_
#define _BLAST2_

#include <ncbi.h>
#include <objblst2.h>
#include <objseq.h>
#define BLAST2_LINE_LENGTH 79

#define TEXT_OUTPUT    0
#define HTML_OUTPUT    1
#define ABSOLUTE_LINKS 2

/*
The following structure contains all the info for the "BlastPrint"
functions (e.g., BlastPrintStructInit, BlastPrintStart etc.)	
*/

typedef struct blast_print_struct {
	CharPtr buffer;		/* contains line to be printed. */
	FILE *fp;		/* File to print to */
	Int2 init_indent,	/* Indentation of first line */
	     cont_indent,	/* indentation of continuation lines. */
	     line_length,	/* size of allocated buffer (above) */
	     position;		/* position in buffer. */
	Boolean first_line;	/* TRUE if first call after BlastStartPrint. */
} BlastPrintStructure, PNTR BlastPrintStructurePtr;

BlastPrintStructurePtr BlastPrintStructInit PROTO((Int2 line_length, FILE *fp));

Boolean BlastPrintStart PROTO((BlastPrintStructurePtr bpsp, Int2 init_indent, Int2 cont_indent));

Int2 BlastPrintStringAdd PROTO((BlastPrintStructurePtr bpsp, CharPtr string));

Int2 BlastPrintCharAdd PROTO((BlastPrintStructurePtr bpsp, Char character));

Int2 BlastPrintDoubleAdd PROTO((BlastPrintStructurePtr bpsp, CharPtr format, double number));

Int2 BlastPrintIntegerAdd PROTO((BlastPrintStructurePtr bpsp, CharPtr format, Int4 number));

Int2 BlastPrintCheckBufferStatus PROTO((BlastPrintStructurePtr bpsp, Char next_char));
void BlastPrintNewLine PROTO((BlastPrintStructurePtr bpsp));

Int2 BlastPrintTabToColumn PROTO((BlastPrintStructurePtr bpsp, Int2 column));

void BlastPrintEnd PROTO((BlastPrintStructurePtr bpsp));

void BlastPrintFlush PROTO((BlastPrintStructurePtr bpsp));

BlastPrintStructurePtr BlastPrintStructClose PROTO((BlastPrintStructurePtr bpsp));

typedef struct blast_report_struct {
	BLAST0ResultPtr	result;	/* contains histogram & hitlist */
	BLAST0ResponsePtr blresp;/* contains preface, matrix, etc. */
	Int4 query_length,	/* number of residues/basepairs in query. */ 
		qoffset,	/* offset for command-line option "-qoffset" */
		num_of_defline,	/* number of (one-line) descriptions to show*/
		num_of_align;	/* number of alignments to show*/
	CharPtr program;	/* name of the program (e.g., blastn) */
	Boolean is_prot,	/* TRUE if the protein alignments are shown */ 
		get_gi;		/* TRUE if gi's are shown in FASTA id */
	Int4    type;           /* TEXT_OUTPUT or HTML_OUTPUT         */
	BlastPrintStructurePtr bpsp;	/* Printing info in this structure. */
	Int2Ptr	window;		/* Array to save actual window sizes for
				Stephen A's development of BLAST. */
	Int2Ptr	right_dropoff_sizes;
	Int2Ptr	left_dropoff_sizes;
} BlastReportStruct, PNTR BlastReportStructPtr;

Boolean LIBCALL TraditionalBlastReportAddResult PROTO((BlastReportStructPtr blrp, BLAST0ResultPtr result, BLAST0ResponsePtr blresp));

BlastReportStructPtr TraditionalBlastReportSetUp PROTO((BLAST0ResultPtr result, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp, Int4 type));

BlastReportStructPtr TraditionalBlastReportCleanUp PROTO((BlastReportStructPtr));

Boolean LIBCALL TraditionalHistBlastOutput PROTO((BlastReportStructPtr blrp));

Boolean LIBCALL PrintTraditionalBlastPreface PROTO((BlastReportStructPtr blrp));

Int4 LIBCALL acknowledge_blast_request PROTO((CharPtr program, BLAST0SequencePtr query, Int2 strands, Int4 offset, FILE *fp, Int4 type));

Int2 LIBCALL BlastPrintValNodeStack PROTO((ValNodePtr stack, CharPtr title, FILE *fp));

void blast2_wrap PROTO((FILE *fp, CharPtr title, CharPtr s, int slen, int linelen, int indent));

BioseqPtr LIBCALL BLAST0Sequence2Bioseq PROTO((BLAST0SequencePtr sequence));

Int2 LIBCALL BLAST0SeqData2Bioseq PROTO((BioseqPtr bsp, ValNodePtr sequence, Int4 length));

/* This function produces the traditional BLAST output. */
Boolean LIBCALL TraditionalBlastOutput PROTO((BLAST0ResultPtr hdp, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp));

/* This function produces the traditional BLAST output in HTML format */
Boolean LIBCALL TraditionalBlastOutputHTML PROTO((BLAST0ResultPtr hdp, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp));

/* This function produces the BLAST output in short HTML format */
Boolean LIBCALL TraditionalBlastOutputHTML2 PROTO((BLAST0ResultPtr hdp, BLAST0ResponsePtr blresp, CharPtr program, FILE *fp));

/* These functions produce parts of the traditional output and are used
by TraditionalBlastOutput */

Boolean LIBCALL TraditionalHeadBlastOutput PROTO((BlastReportStructPtr blrp));

Boolean LIBCALL TraditionalHistBlastOutput PROTO((BlastReportStructPtr blrp));

void LIBCALL TraditionalBlastWarning PROTO((BlastReportStructPtr blrp));


Boolean LIBCALL TraditionalTopBlastOutput PROTO((BlastReportStructPtr blrp));

Boolean LIBCALL TraditionalBottomBlastOutput PROTO((BlastReportStructPtr blrp));

Boolean LIBCALL TraditionalTailBlastOutput PROTO((BlastReportStructPtr blrp));


#endif /* _BLAST2_ */
