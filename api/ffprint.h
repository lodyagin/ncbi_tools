/*   ffprint.h
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  ffprint.h
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov
*
* Version Creation Date:   7/15/95
*
* $Revision: 6.14 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#ifndef _FFPRINT_
#define _FFPRINT_

#include <asn2ffg.h>

/*--------------the structure for the buffered printing-----*/

#define LINKS 20

#define LOCUS_line 1
#define DEF_line 2
#define ACC_line 3
#define NID_line 4
#define KW_line 5
#define SOURCE_line 6
#define BASECOUNT_line 7
#define ORIGIN_line 8
#define SEQ_line 9

typedef void (*HeadTailProc) PROTO((Pointer, FILE*));

typedef struct buffstruct {

/* The next eight variables are used by the "printing" utilities of asn2ff6.c
(StartPrint, AddChar, CheckBufferState, NewContLine) to perform the "buffered"
printing */ 
	CharPtr buffer;		/* buffer to hold line */
	Int2 init_indent;  /*indentation of the first line, set by StartPrint */
	Int2 cont_indent;  /*indentation of continuation lines */
	Int2 line_max;	/* maximum allowable length of line, set in StartPrint*/
	CharPtr line_prefix; /* prefix, such as "ID" on EMBL id lines */
	Char newline;		/* newline character */
	FILE *fp;		/* file to print to. */
	ByteStorePtr byte_sp;	/* Used to save paragraph (i.e., several lines)
				until printing. */
	CharPtr line_return;
/*  next three variables are used for HTML URLs   */
	Int4  PNTR	pos_links;
	CharPtr	PNTR links;
	Int2	n_links;
	Int2	buf_n_links;
} BuffStruct, PNTR BuffStructPtr;


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

NLM_EXTERN void LIBCALL asn2ff_set_output PROTO((FILE *fp, CharPtr line_return));
NLM_EXTERN CharPtr LIBCALL ff_MergeString PROTO((void));
NLM_EXTERN CharPtr LIBCALL FFPrint PROTO((FFPrintArrayPtr pap, Int4 index, Int4 pap_size));
NLM_EXTERN void LIBCALL FFBSPrint PROTO((FFPrintArrayPtr pap, Int4 index, Int4 pap_size));
NLM_EXTERN void LIBCALL ff_print_string PROTO((FILE *fp, CharPtr string, CharPtr line_return));
NLM_EXTERN CharPtr LIBCALL ff_print_string_mem PROTO((CharPtr string));
NLM_EXTERN Int2 LIBCALL ff_StartPrint PROTO((Int2 init_indent, Int2 cont_indent, Int2 line_max, CharPtr line_prefix));
NLM_EXTERN void LIBCALL ff_AddString PROTO((CharPtr string));
NLM_EXTERN void LIBCALL ff_AddInteger PROTO((CharPtr fmt, long integer));
NLM_EXTERN void LIBCALL AddLink PROTO((CharPtr str));
NLM_EXTERN void LIBCALL AddLinkLater PROTO((CharPtr str, Int2 prevlen));
NLM_EXTERN void LIBCALL ff_AddChar PROTO((Char character));
NLM_EXTERN void LIBCALL PrintXX PROTO((void));
NLM_EXTERN void LIBCALL ff_AddStringWithTildes PROTO ((CharPtr string));
NLM_EXTERN void LIBCALL ChangeStringWithTildes PROTO ((CharPtr string));
NLM_EXTERN CharPtr LIBCALL CheckBufferState  PROTO((Int2Ptr increment_string, Char next_char));
NLM_EXTERN Int2 LIBCALL NewContLine PROTO((void));
NLM_EXTERN Int2 LIBCALL TabToColumn PROTO((Int2 column));
NLM_EXTERN void LIBCALL ff_EndPrint PROTO((void));
NLM_EXTERN void LIBCALL FlushBuffer PROTO((void));
NLM_EXTERN CharPtr LIBCALL CheckEndPunctuation PROTO((CharPtr string, Char end));
NLM_EXTERN CharPtr PrintDate PROTO((NCBI_DatePtr date));
NLM_EXTERN void LIBCALL BuffFree PROTO((void));
NLM_EXTERN void LIBCALL init_buff PROTO((void));
NLM_EXTERN void LIBCALL init_buff_ex PROTO((Int2 init_size));
NLM_EXTERN void LIBCALL free_buff PROTO((void));
NLM_EXTERN void LIBCALL init_www PROTO((void));
NLM_EXTERN void LIBCALL fini_www PROTO((void));
NLM_EXTERN void LIBCALL head_tail_ff PROTO((Pointer mydata, HeadTailProc headfun, HeadTailProc tailfun));
NLM_EXTERN void LIBCALL head_www PROTO((FILE *fp, SeqEntryPtr sep));
NLM_EXTERN void LIBCALL tail_www PROTO((FILE *fp));
NLM_EXTERN Boolean LIBCALL get_www PROTO((void));
NLM_EXTERN Boolean LIBCALL www_muid PROTO((Int4 muid));
NLM_EXTERN Boolean LIBCALL www_gcode PROTO((CharPtr gcode));
NLM_EXTERN Boolean LIBCALL www_source PROTO((CharPtr orgname, OrgRefPtr orp));
NLM_EXTERN Boolean LIBCALL www_organism PROTO((CharPtr orgname, Int4 id));
NLM_EXTERN Boolean LIBCALL www_taxid PROTO((CharPtr orgname, Int4 id));
NLM_EXTERN Boolean LIBCALL www_extra_acc PROTO((CharPtr acc, Boolean ncbi));
NLM_EXTERN Boolean LIBCALL www_note_gi PROTO((CharPtr str));
NLM_EXTERN Boolean LIBCALL www_db_xref PROTO((CharPtr str));
NLM_EXTERN Boolean LIBCALL www_protein_id PROTO((CharPtr str));
NLM_EXTERN Boolean LIBCALL www_map PROTO((CharPtr str));
NLM_EXTERN Boolean LIBCALL www_genpept_gi PROTO((CharPtr str));
NLM_EXTERN Boolean LIBCALL www_dbsource PROTO((CharPtr str, Boolean first, Uint1 choice));
NLM_EXTERN Boolean LIBCALL www_xref PROTO((CharPtr str, Uint1 link));
NLM_EXTERN Boolean LIBCALL www_xref_button PROTO((FILE *fp, CharPtr str, Uint1 link, Uint1 db));
NLM_EXTERN CharPtr LIBCALL ReportPrint PROTO((FFPrintArrayPtr pap, Int4 index, Int4 pap_size));
NLM_EXTERN Boolean LIBCALL PrintSPBlock PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN Boolean LIBCALL ff_PrintLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp, Int2 type));
NLM_EXTERN CharPtr LIBCALL www_featloc PROTO((CharPtr loc));
NLM_EXTERN void LIBCALL GetHelpMsg PROTO((SeqEntryPtr sep));
NLM_EXTERN void LIBCALL www_PrintComment  PROTO((CharPtr string, Boolean identifier, Uint1 format));
NLM_EXTERN Boolean LIBCALL www_featkey PROTO((CharPtr key, BIG_ID gi, Int2 entityID, Uint4 itemID));
NLM_EXTERN void LIBCALL www_accession PROTO((CharPtr string));
NLM_EXTERN void LIBCALL ff_RecalculateLinks(Int4 indent);

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* _FFPRINT_ */
