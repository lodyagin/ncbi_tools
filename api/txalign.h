/* $Id: txalign.h,v 6.25 2000/06/09 19:00:06 shavirin Exp $
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
* File Name:  $RCSfile: txalign.h,v $
*
* Author:  Jinghui Zhang
*
* Initial Version Creation Date: 03/13/94
*
* $Revision: 6.25 $
*
* File Description:
*         External include file for various alignments
* Revision 5.13  1997/06/05 20:55:34  madden
 * Added PrintDefLinesFromSeqAlign prototype
*
*
* $Log: txalign.h,v $
* Revision 6.25  2000/06/09 19:00:06  shavirin
* Function GetGeneticCodeFromSeqId() made external and added to header file.
*
* Revision 6.24  2000/06/08 20:44:50  shavirin
* Added calculation of start/stop values in the function find_score_in_align().
*
* Revision 6.23  2000/03/07 21:58:41  shavirin
* Now will use PSSM Matrix to show positives in PSI Blast
*
* Revision 6.22  1999/11/24 21:24:33  vakatov
* Fixed for the C++ and/or MSVC DLL compilation
*
* Revision 6.21  1999/11/09 22:15:08  shavirin
* Added parameter follower to the Blast score printing function
*
* Revision 6.20  1999/10/07 16:08:05  shavirin
* Passed matrix to the function FormatScoreFromSeqAlign().
*
* Revision 6.19  1999/09/29 17:15:38  shavirin
* Added new funtion FormatScoreFromSeqAlign()
*
* Revision 6.18  1999/06/07 18:43:17  madden
* added TXALIGN_NO_DUMPGNL if dumpgnl is not desired
*
* Revision 6.17  1999/04/15 20:57:23  madden
* overview printing for vector stuff
*
* Revision 6.16  1999/04/06 15:13:25  madden
* Add support for non-gnl queries with dumpgnl syntax
*
* Revision 6.15  1999/02/26 21:28:06  victorov
* taking different sections of config file depending on WWW_BLAST_TYPE
*
* Revision 6.14  1999/02/19 20:51:07  victorov
* changed URL to the tool reporting incomplete
* sequences. URL now includes starts/stops for all hits
*
* Revision 6.13  1999/01/13 21:52:43  victorov
* added links to incomplete genomes in hit details
*
* Revision 6.12  1999/01/06 22:51:00  victorov
* added hyperlinks for incomplete sequences
*
* Revision 6.11  1999/01/05 14:52:04  madden
* Add frame and strand information
*
* Revision 6.10  1998/11/09 19:06:47  vakatov
* Added "NLM_EXTERN" to the ShowTextAlignFromAnnotExtra() prototype
*
* Revision 6.9  1998/09/01 13:27:02  madden
* PrintDefLinesExtra function
*
* Revision 6.8  1998/08/26 21:32:23  madden
* Added ShowTextAlignFromAnnotExtra for PHI-BLAST
*
* Revision 6.7  1998/07/23 13:35:29  egorov
* Allow print specified number of descriptions in PrintDefLinesFromSeqAlign()
*
* Revision 6.6  1998/07/02 21:22:56  madden
* Changes for random-access BLAST
*
* Revision 6.5  1998/03/25 22:38:50  egorov
* Change prototypes for PrintDefLinesFromAnnot and PrintDefLinesFromSeqAlign
*
* Revision 6.4  1997/10/06 14:01:11  zjing
* move TxGetSubjectId, GetScoreAndEvalue to sequtil.ch
*
* Revision 6.3  1997/09/25 17:17:37  zjing
* Add the option for showing blunt-end alignment
*
* Revision 6.2  1997/09/25 02:00:27  vakatov
* Added NLM_EXTERN specifier to some functions(necessary for MS-Win DLLs)
*
* Revision 6.1  1997/09/18 22:24:23  madden
* Made TxGetSubjectIdFromSeqAlign public
*
* Revision 6.0  1997/08/25 18:08:14  madden
* Revision changed to 6.0
*
* Revision 5.20  1997/08/14 17:55:49  zjing
* minor changes
*
* Revision 5.18  1997/07/28 13:55:46  madden
* Added mask_loc to prototypes.
*
* Revision 5.17  1997/07/11 15:28:13  madden
* Added TXALIGN_HTML_RELATIVE define
*
* Revision 5.16  1997/07/07 20:22:26  madden
* changes to show the results as query-subect
*
* Revision 5.15  1997/06/19 18:39:42  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.14  1997/06/09 21:47:25  madden
* Added Boolean follower to AlignStatOptions
*
 * Revision 5.12  1997/03/13  21:43:13  shavirin
 * Added protection for C++ compiler
 *
*
* ==========================================================================
*/
#ifndef _TXALIGN_
#define _TXALIGN_

/****************************************************************************/
/* INCLUDES */
/****************************************************************************/

#include <jzcoll.h>
#include <ffprint.h>

/****************************************************************************/
/* DEFINES */
/****************************************************************************/

#define WEBB_asize 23		/* webb's matrix */
#define TX_MATRIX_SIZE 128	/*size of the matrix for showing the 
                                  text alignment*/

#define TXALIGN_LOCUS_NAME	((Uint4)256)	/*display the locus name*/
#define TXALIGN_MASTER		((Uint4)2)	/*display the alignment as multiple pairwise alignment*/
#define TXALIGN_MISMATCH	((Uint4)4)	/*display the mismatched residue of the sequence */
#define TXALIGN_MATRIX_VAL	((Uint4)8)	/*display the matrix of the alignment */
#define TXALIGN_HTML		((Uint4)16)	/*display the format in a HTML page*/
#define TXALIGN_HTML_RELATIVE	((Uint4)8192)	/*the HTML (if enabled by TXALIGN_HTML) should be relative*/
#define TXALIGN_SHOW_RULER	((Uint4)32)	/*display the ruler for the text alignment*/
#define TXALIGN_COMPRESS	((Uint4)64)	/*make the space for label smaller*/
#define TXALIGN_END_NUM		((Uint4)128)	/*show the number at the end */
#define TXALIGN_FLAT_INS	((Uint4)1)	/*flat the insertions in multiple pairwise alignment */
#define TXALIGN_SHOW_GI         ((Uint4)512)    /*show the gi in the defline. */
#define TXALIGN_SHOW_NO_OF_SEGS ((Uint4)1024)    /*show the number of (sum statistics) segments in the one-line descriptions? */

#define TXALIGN_BLASTX_SPECIAL  ((Uint4)2048)	/*display the BLASTX results 
						  as protein alignment */
#define	TXALIGN_SHOW_QS		((Uint4)4096)	/*show the results as query-subect*/
#define TXALIGN_SPLIT_ANNOT	((Uint4)16384)	/*for Seq-annot from the same alignment, split the 
												the display into individual panel*/
#define TXALIGN_SHOW_STRAND	((Uint4)32768)	/*for displaying the stradn even in the compact form*/
#define TXALIGN_BLUNT_END	((Uint4)65536)	/*showing the blunt-end for the end gaps*/
#define TXALIGN_DO_NOT_PRINT_TITLE	((Uint4)131072)	/* do not print title before list of deflines */
#define TXALIGN_CHECK_BOX	((Uint4)262144)	/* place checkbox before the line (HTML only) */
#define TXALIGN_CHECK_BOX_CHECKED	((Uint4)524288)	/* make default value for checkboxes ON (HTML only) */
#define TXALIGN_NEW_GIF		((Uint4)1048576)	/* print new.gif near new alignments (HTML only) */
#define TXALIGN_NO_ENTREZ	((Uint4)2097152)	/* Use dumpgnl syntax instead of ENTREZ. */
#define TXALIGN_NO_DUMPGNL	((Uint4)4194304)	/* No dumpgnl output, even if GNL. */
/*
	Used by psi-blast to distinguish first from subsequent passes.
*/

#define FIRST_PASS 1
#define NOT_FIRST_PASS_REPEATS 2
#define NOT_FIRST_PASS_NEW  3

/****************************************************************************/
/* TYPEDEFS */
/****************************************************************************/

typedef struct text_buf{	/*for a generic feature comment*/
	Int4 pos;	        /*position for label*/
	Uint1 strand;	        /*the orientation*/
	CharPtr label;	        /*label for the feature*/
	CharPtr buf;	   /*the buffer for features other than cds for aa*/
	Int2Ptr matrix_val;	/*the value of each residue from the matrix */
	CharPtr codon[3];	/*for features such as cds for aa*/
	Int2 frame;	        /*for cds for feature*/
	Int4 f_pos;	        /*position of the current buf*/
	Uint2 exonCount;	/*count the number of exons, useded in 
                                  cds for aa*/
        Uint2 itemID;	/*feature's itemID. It is used to check identity*/
	Uint2 feattype;
	Uint2 subtype;
	Uint2 entityID;
	Uint2 seqEntityID;	/*the entityID for the sequence*/
	Uint2 bsp_itemID;	/*itemID for the Bioseqs*/
	Boolean extra_space;
}TextAlignBuf, PNTR TextAlignBufPtr;

typedef struct align_summary {
    Int4 positive;	        /*number of positive residues*/
    Int4 identical;	        /*number of identical residues*/
    Int4 gaps;		        /*number of the gaps*/
    Int4 totlen;	        /*total length of the alignemtns*/
    Int4Ptr PNTR matrix;	/*matrix for protein alignments*/
    Int4Ptr PNTR posMatrix;	/*matrix for PSSM protein alignments*/
    SeqIdPtr master_sip;	/*the Seq-id of the master sequence*/
    SeqIdPtr target_sip;	/*the Seq-id for the target sequence*/
    Boolean is_aa;              /*are the sequences nucleotide or protein?*/
    Uint1 m_strand,	        /* strand of the query. */
          t_strand;	        /* strand of the database sequence. */
    Int4  m_frame,	        /* Frame of the query. */
        t_frame;	        /* Frame of the database sequence. */
    Boolean	m_frame_set,	/* query frame was set. */
        t_frame_set;	        /* database sequence frame was set. */
    Int4 master_from;           /* from for master sequence */
    Int4 master_to;             /* to for master sequence */
    Int4 target_from;           /* from for target sequence */
    Int4 target_to;             /* to region for master sequence */
}AlignSum, PNTR AlignSumPtr;

typedef struct align_stat_option { /*options for printing the statistics*/
	Int2 line_len;
	Int2 indent_len;
	Boolean html_hot_link;			/* Prepare HTML output. */
	Boolean html_hot_link_relative;		/* Make the HTML link relative. */
	Boolean show_gi;
	Boolean no_entrez;			/* Do not use Entrez format for HTML links. */
	Boolean no_dumpgnl;			/* Do not use dumpgnl format even if GNL. */
	FILE *fp;
	CharPtr buf;
	BioseqPtr bsp;
	ScorePtr sp;
	Int4 identical;         /*number of identical residues*/
	Int4 gaps;		/*number of the gaps*/
	Int4 positive;	        /*number of the positive residues*/
	Int4 align_len;	  /*the length of the alignment. EXCLUDE the GAPS*/
	Boolean follower; /* If TRUE, this is NOT the first alignment for this sequences. */
	Uint1 	m_strand,	/* strand of the query. */
		t_strand;	/* strand of the database sequence. */
	Int2	m_frame,	/* Frame of the query. */
		t_frame;	/* Frame of the database sequence. */
	CharPtr segs; /* <start> "-" <stop> ("," <start> "-" <stop>)* */
	CharPtr db_name; /* searched databases list */
	CharPtr blast_type; /* string used to choose proper config parms */
}AlignStatOption, PNTR AlignStatOptionPtr;

/****************************************************************************/
/* FINCTION DEFINITIONS */
/****************************************************************************/

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************
*
*	find_score_in_align(align, chain, asp)
*	align: the Seq-align point
*	chain: for multiple segment Seq-aligns, such as DenseDiag and 
*	StdSeg, the order within the Seq-align
*	asp:   the structure that records and stores the positive, 
*			identical residues
*	the function only works for DenseDiag and Stdseg for now
*
*****************************************************************/
NLM_EXTERN ScorePtr find_score_in_align PROTO((SeqAlignPtr align, 
                                               Uint2 chain, AlignSumPtr asp));

/*the default formatting function for printing the scores*/

NLM_EXTERN int LIBCALLBACK FormatScoreFunc PROTO((AlignStatOptionPtr asop));

/**********************************************************************************
*
*       Given a chain of annots (ValNodePtrs) they are all printed out, one pattern
*       at a time.
*
*************************************************************************************/

NLM_EXTERN Boolean LIBCALL ShowTextAlignFromAnnotExtra PROTO((BioseqPtr bsp, ValNodePtr vnp, SeqLocPtr seqloc,
        Int4 line_len, FILE *fp,
        Uint1Ptr featureOrder, Uint1Ptr groupOrder, Uint4 option, Int4Ptr PNTR matrix,
        ValNodePtr mask_loc, int (LIBCALLBACK *fmt_score_func)PROTO((AlignStatOptionPtr))));

/*****************************************************************************
*
*	ShowTextAlignFromAnnot(annot, locus, line_len, fp, master, f_order)
*	display the alignment stored in a Seq-annot in a text file
*	annot: the Seq-annot pointer
*	locus: if TRUE, show the locus name as the sequence label, otherwise, 
*		use the accession
*	line_len: the number of sequence char per line
*	fp: The file pointer to store the text output
*	master: if TRUE, show the result as a master-slave type multiple pair
*	wise alignment. if FALSE, display one alignment after the other
*	f_order: the user selected feature type and order to be shown together
*	with the alignment
*	return TRUE for success, FALSE for fail
*
*****************************************************************************/
NLM_EXTERN Boolean ShowTextAlignFromAnnot PROTO((
                    SeqAnnotPtr hannot, Int4 line_len, 
                    FILE *fp, Uint1Ptr featureOrder, 
                    Uint1Ptr groupOrder, Uint4 option, 
                    Int4Ptr PNTR matrix, ValNodePtr mask_loc,
                    int (LIBCALLBACK *fmt_score_func)
                    PROTO((AlignStatOptionPtr))
                    ));
/**
 * same as ShowTextAlignFromAnnot
 * the db_name argument is used to make links to
 * incomplete genomes
 */
NLM_EXTERN Boolean ShowTextAlignFromAnnot2 PROTO((
                    SeqAnnotPtr hannot, Int4 line_len, 
                    FILE *fp, Uint1Ptr featureOrder, 
                    Uint1Ptr groupOrder, Uint4 option, 
                    Int4Ptr PNTR matrix, ValNodePtr mask_loc,
                    int (LIBCALLBACK *fmt_score_func)
                    PROTO((AlignStatOptionPtr)),
                    CharPtr db_name,
                    CharPtr blast_type
                    ));
/**
 * same as ShowTextAlignFromAnnot
 * the posMatrix used to show alignments using PSSM
 */
NLM_EXTERN Boolean ShowTextAlignFromAnnot3 PROTO((
                    SeqAnnotPtr hannot, Int4 line_len, 
                    FILE *fp, Uint1Ptr featureOrder, 
                    Uint1Ptr groupOrder, Uint4 option, 
                    Int4Ptr PNTR matrix, ValNodePtr mask_loc,
                    int (LIBCALLBACK *fmt_score_func)
                    PROTO((AlignStatOptionPtr)),
                    CharPtr db_name,
                    CharPtr blast_type,
                    Int4Ptr PNTR posMatrix
                    ));


/***********************************************************************
*
*	ShowAlignNodeText(anp_list, num_node, line_len, locus,
*	fp)
*	convert the alignment data in the list of AlignNode into text written
*	to a file
*	anp_list: a list (ValNodePtr) of AlignNode processed from Seq-aligns
*	num_node: the number of AlignNode to be processed currently. It can
*	be used in the cases where only the top num_node in the anp_list is 
*	going to be processed. This can be useful to make vertically cashed
*	buffer
*	line_len: the length of sequence char per line
*	locus: if TRUE, show the locus name
*	fp: the file Pointer
*	left: the leftmost position for display
*	right: the rightmost position for display
*	align_type:	the type of alignment. DNA-protein alignment?
*
*	return TRUE for success, FALSE for fail
*
************************************************************************/

NLM_EXTERN Boolean ShowAlignNodeText PROTO((
                     ValNodePtr anp_list, Int2 num_node, 
                     Int4 line_len, FILE *fp, Int4 left, 
                     Int4 right, Uint4 option, 
                     Int4Ptr PNTR u_matrix, 
                     int (LIBCALLBACK *fmt_score_func)
                     PROTO((AlignStatOptionPtr))
                     ));

NLM_EXTERN Boolean ShowAlignNodeText2 PROTO((
                     ValNodePtr anp_list, Int2 num_node, 
                     Int4 line_len, FILE *fp, Int4 left, 
                     Int4 right, Uint4 option, 
                     Int4Ptr PNTR u_matrix, 
                     int (LIBCALLBACK *fmt_score_func)
                     PROTO((AlignStatOptionPtr)),
                     CharPtr db_name,
                     CharPtr blast_type,
                     Int4Ptr PNTR posMatrix
                     ));

/***********************************************************************
*
*	ProcessTextAlignNode(anp, left, right, p_stop, m_buf, locus)
*	process an AlignNode to generate a list of text buffer
*
*	anp: the AlignNode
*	left, right: the range of alignment in process. mapped to 
*	anp->extremes.left, and anp->extremes.right
*	p_stop: the previous stop position in the sequence. It is used 
*	to label the position of line which is a gap
*	m_buf: the buffer of the master sequence. Can be used to compare
*	mismatches
*	locus: if TRUE, use the locus name for sequence
*
*
*
************************************************************************/
NLM_EXTERN ValNodePtr ProcessTextAlignNode PROTO((
                    AlignNodePtr anp, Int4 m_left, 
                    Int4 m_right, Int4Ptr p_stop, 
                    CharPtr m_buf, Int4 line_len, 
                    Int1 m_frame, 
                    Uint4 option, Int4Ptr PNTR matrix
                    ));
NLM_EXTERN ValNodePtr ProcessTextAlignNode2 PROTO((
                    AlignNodePtr anp, Int4 m_left, 
                    Int4 m_right, Int4Ptr p_stop, 
                    CharPtr m_buf, Int4 line_len, 
                    Int1 m_frame, 
                    Uint4 option, Int4Ptr PNTR matrix,
                    Int4Ptr PNTR posMatrix, Int4 q_start
                    ));


NLM_EXTERN ValNodePtr FreeTextAlignList PROTO((ValNodePtr tdp_list));

/*
  Print a summary of the Sequences producing significant alignments.
*/

NLM_EXTERN Boolean LIBCALL PrintDefLinesExtra PROTO((
		ValNodePtr vnp, Int4 line_length, FILE *outfp, Uint4 options, 
		Int4 mode, Int2Ptr marks, SeqLocPtr seqloc));


NLM_EXTERN Boolean LIBCALL PrintDefLinesFromAnnot PROTO((
                    SeqAnnotPtr seqannot, 
                    Int4 line_length, FILE *fp, 
                    Uint4 options, Int4 mode, 
                    Int2Ptr marks
                    ));

NLM_EXTERN Boolean LIBCALL PrintDefLinesFromSeqAlign PROTO((
                    SeqAlignPtr seqalign, 
                    Int4 line_length, FILE *fp, 
                    Uint4 options, Int4 mode, 
                    Int2Ptr marks
                    ));

NLM_EXTERN Boolean LIBCALL PrintDefLinesFromSeqAlignEx PROTO((
		    SeqAlignPtr seqalign, 
		    Int4 line_length, 
		    FILE *outfp, 
		    Uint4 options, 
		    Int4 mode, 
		    Int2Ptr marks, 
		    Int4 number_of_descriptions
		    ));

NLM_EXTERN Boolean LIBCALL PrintDefLinesFromSeqAlignEx2 PROTO((
		    SeqAlignPtr seqalign, 
		    Int4 line_length, 
		    FILE *outfp, 
		    Uint4 options, 
		    Int4 mode, 
		    Int2Ptr marks, 
		    Int4 number_of_descriptions,
		    CharPtr db_name,
		    CharPtr blast_type
		    ));

/*
	Fills in the slots with score, bit_score, etc. from the SeqAlign.
*/


/* setting up the matrix for the positive residue of the alignment */

NLM_EXTERN Int4Ptr PNTR load_default_matrix PROTO((void));
NLM_EXTERN void free_default_matrix PROTO((Int4Ptr PNTR matrix));


/*options for display of the text alignment*/
#define TEXT_MP_MISMATCH	1	/*multiple pairwise alignment with mismatch*/
#define TEXT_MP			2	/*multiple pairwise without mismatch*/
#define TEXT_MPFLAT_MISMATCH	3	/*flat multile with mismatch*/
#define TEXT_MPFLAT		4	/*flat multiple without mismatch*/
#define TEXT_BLAST		5	/*traditional blast output*/


/*can the current alignnode be printed for text view*/
NLM_EXTERN Boolean PrintAlignForText PROTO((AnnotInfoPtr info, AlignNodePtr anp));

/*
*
*	determine the option for alignment based on the named tx_option
*
*/
NLM_EXTERN Uint4 GetTxAlignOptionValue PROTO((Uint1 tx_option, BoolPtr hide_feature, 
	BoolPtr print_score, BoolPtr split_display));

/*
        Gets the SeqIdPtr for the subject sequence from the first SeqAlign.
        The SeqIdPtr is not saved  and should not be deleted.
*/

/* Marks structure is used for PSI Blast to print .gif marsk 
   near alignments and to check for convergence */

#define	SEQ_ALIGN_MARK_PREVGOOD		1
#define	SEQ_ALIGN_MARK_PREVCHECKED	2
/* the following serves only for old stuff which uses posRepeat... */
#define	SEQ_ALIGN_MARK_REPEAT		4


typedef struct MarkSeqAlign {
    Int4		kind;	/* bitmask for the mark */
    struct MarkSeqAlign	*next;
} MarkSeqAlign, PNTR MarkSeqAlignPtr;


NLM_EXTERN SeqIdPtr LIBCALL GetUseThisGi PROTO((SeqAlignPtr seqalign));
NLM_EXTERN Boolean LIBCALL FilterTheDefline PROTO((BioseqPtr bsp, SeqIdPtr gi_list_head, CharPtr buffer_id, Int4 buffer_id_length, CharPtr PNTR titlepp));


/* Printoverview stuff. */
NLM_EXTERN Boolean LIBCALL MakeDisplaySeqLoc PROTO((SeqAlignPtr PNTR seqalign_ptr, ValNodePtr PNTR vnp, Int4 length));
NLM_EXTERN Boolean LIBCALL PrintOverviewFromSeqLocs PROTO((ValNodePtr vnp, Int4 query_length, FILE *outfp));

NLM_EXTERN Boolean FormatScoreFromSeqAlign
(SeqAlignPtr sap, Uint4 option, FILE *fp,
Int4Ptr PNTR matrix, Boolean follower);    

/*
  Obtains the genetic code from a BioseqPtr, assuming that a fetch function
  has been enabled.
*/
NLM_EXTERN CharPtr GetGeneticCodeFromSeqId (SeqIdPtr sip);

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
