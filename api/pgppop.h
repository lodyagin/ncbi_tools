/*  $Id: pgppop.h,v 6.15 1999/09/16 18:52:27 durand Exp $
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
* File Name:  pgppop.h
*
* Author:  Patrick Durand
*
* Version Creation Date:   05/03/99
*
* $Revision: 6.15 $
*
* File Description: 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: pgppop.h,v $
* Revision 6.15  1999/09/16 18:52:27  durand
* redesign the PopSet viewer toolbar
*
* Revision 6.14  1999/09/10 14:27:10  durand
* add RulerDescr in MsaParaGPopList struct to be used for complex ruler types
*
* Revision 6.13  1999/09/02 17:32:23  durand
* add some defines for new display types in DDV
*
* Revision 6.12  1999/08/31 16:45:39  durand
* add sap in MsaParaGPopList structure
*
* Revision 6.11  1999/08/31 13:51:31  durand
* add PubMed link for PopSet Viewer
*
* Revision 6.10  1999/08/30 20:29:56  durand
* make PrintSeqAlignCallback estern function
*
* Revision 6.9  1999/08/30 14:18:09  durand
* use ByteStore to format the SeqAlign output
*
* Revision 6.8  1999/08/16 18:57:00  lewisg
* made DDV_GetSeqAlign extern and added prototype to header
*
* Revision 6.7  1999/08/06 21:43:16  chappey
* SeqAlignToBS new function to save in ByteStore structure the text output of the SeqAlign(s) packaged in a SeqEntry
*
* Revision 6.6  1999/07/22 13:23:14  durand
* made DDV_SearchAli external function
*
* Revision 6.5  1999/07/21 21:52:13  durand
* add some functions to display a summary for a PopSet entry
*
* Revision 6.4  1999/07/20 17:01:48  durand
* add eID field in URLData structure for PopSet Viewer
*
* Revision 6.3  1999/07/19 21:16:02  durand
* add DDV_ResetParaGSeqAlignCoord to reset the seqalign coord in the display data structures of DDV
*
* Revision 6.2  1999/07/15 18:20:51  durand
* add display options to support BLAST outputs
*
* Revision 6.1  1999/07/09 13:59:58  durand
* move pgppop from desktop to api
*
* Revision 6.2  1999/07/06 22:31:24  kans
* removed whitespace before #define directives
*
* Revision 6.1  1999/07/06 20:18:07  kans
* initial public checkin
*
* Revision 1.30  1999/07/06 18:54:27  durand
* add new features for the display of PopSet viewer
*
* Revision 1.29  1999/07/02 13:22:18  durand
* fix bugs for the display of minus strand sequences
*
* Revision 1.28  1999/06/29 16:48:11  shavirin
* Changed definition of function DDV_ShowSeqAlign()
*
* Revision 1.27  1999/06/28 22:07:22  durand
* add loader functions and clean the code with Lint and Purify
*
* Revision 1.26  1999/06/25 14:21:07  durand
* add a new command-line to wwwddv.cgi
*
* Revision 1.25  1999/06/24 20:49:42  shavirin
* Added new function DDV_ShowSeqAlign().
*
* Revision 1.24  1999/06/23 17:24:24  durand
* use a binary encoding to manage the display styles
*
* Revision 1.23  1999/06/21 18:37:57  durand
* update DDV_DisplayDefaultAlign to produce full text output
*
* Revision 1.22  1999/06/19 18:36:13  durand
* new display procedure
*
* Revision 1.21  1999/06/14 23:49:44  durand
* add function for Vibrant DDV
*
* Revision 1.20  1999/06/11 22:33:02  durand
* add new functions for Vibrant DDV
*
* Revision 1.19  1999/06/11 18:03:14  durand
* add declarations
*
* Revision 1.17  1999/06/11 14:01:07  durand
* add a define _PGPPOP_ to avoid compilation error when mutilple include
*
* Revision 1.16  1999/06/09 21:35:30  durand
* add constructors/destructors for BspInfo struct as well as read seq function
*
*
*
*
* ==========================================================================
*/
#ifndef _PGPPOP_
#define _PGPPOP_

#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>
#include <objseq.h>

#include <blocks.h>
#include <udvseq.h>

/*#define _min_(a,b)        ((a)>=(b)?(b):(a))
#define _max_(a,b)        ((a)>=(b)?(a):(b))
*/
#define WWW_SCRIPT_SIZE 1000
#define LETTER_BLOCK_WIDTH 10
#define LETTER_HALF_BLOCK_WIDTH 5
#define SCALE_TICK 1
#define SCALE_NUM 2
#define TABLE_LEFT_MARGIN 16

#define INFO_SIZE 255
#define ParaG_Size 70

	
/*****************************************************************************
	Populated lists of the elements ready for the display
	
    Note : use the following data structure to display in the 
           Vibrant viewer/editor
	       if (DisplayType==DDV_DISP_HORZ) use 'TableHead '
		   else use 'DisplayVert'
	
*****************************************************************************/
	typedef struct msaparagpoplist {
		Int4 LengthAli;				/*Alignment size after population from index*/
		Int4 nBsp;					/*number of Bioseq after...*/
		SABlockPtr sabp;/*indexed seqalign*/
		SeqAlignPtr sap;
		ValNodePtr PNTR TableHead;	/*Table list of the heads of the ValNode
									ParaG list for each bioseq*/
		ValNodePtr DisplayVert;		/*List of ParaG Ready for display*/
		BspInfoPtr bspp; /*info for each bsp of the seqalign*/
		Uint1 DisplayType;/*see defines above*/
		ValNodePtr RulerDescr; 
	} MsaParaGPopList, PNTR MsaParaGPopListPtr;

/*****************************************************************************
	Command-Line or URL content
*****************************************************************************/
typedef struct urldata{	/*structure containing the decoded URL*/
	Int4		gi;		/*gi to retrieve in ID1*/
	Int4		from;	/*want to display from...*/
	Int4		to;		/*... to (SeqAlign coord)*/
	Uint4		disp_format;/*display format (seebelow)*/
	CharPtr		szFileNameIn;/*used only for Unix command-line*/
	CharPtr		szFileNameOut;/*used only for Unix command-line*/
	SeqEntryPtr sep; /*used by DDV_ToolsScript*/
	CharPtr     szGIs;/*used of the realign stuff*/
	Uint2        eID;/*entityID; used by PopSet viewer to display onSeqAlign
			among others*/
	} URLData, PNTR URLDataPtr;

/*****************************************************************************
	data block to get BioSource(s) 
*****************************************************************************/
typedef struct popsource{
	CharPtr szTaxName;
	CharPtr szCommonName;
} PopSource, PNTR PopSourcePtr;
	
/*****************************************************************************
	data block for enhanced display options
*****************************************************************************/
typedef struct ddvoptionsblock{	
	Int2 LineSize;/*size of a sequence line*/
	ByteStorePtr PNTR byteSpp;/*use this to store the formatted align as a
	    ByteStore*/
	} DDVOptionsBlock, PNTR DDVOptionsBlockPtr;
	
/*************************************************************************
	DISPLAY FORMATS. the following values can be putted in a single Uint4
		using the binary operator '|='
*************************************************************************/
		/*ruler styles : left, right, top (all), ticks*/
#define RULER_TOP ((Uint4)1)
#define RULER_LEFT ((Uint4)2)
#define RULER_RIGHT ((Uint4)4)
#define RULER_TICK ((Uint4)8)
		/*display enhancement : color, show by block of 10 letters*/
#define DISPE_COLOR ((Uint4)16)
#define DISPE_SHOWBLOCK ((Uint4)32)
		/*seq. letters to show : all, variations*/
#define VIEW_FULLSEQ ((Uint4)64)
#define VIEW_VARIA ((Uint4)128)
		/*output format */
#define DISP_FULL_HTML ((Uint4)256)
#define DISP_PHYLIP_TXT ((Uint4)512)
#define DISP_FASTA_NOGAP ((Uint4)1024)
#define DISP_FASTA_GAP ((Uint4)2048)
#define DISP_FULL_TXT ((Uint4)4096)
		/*content of a line */
		/* a line is like that : [...] Name < 4125 AGCT AGCT ... 4254*/
		/*Name is always displayed*/
#define DISP_ORDER_NUM ((Uint4)8192) /*[...]*/
#define DISP_STRAND ((Uint4)16384)   /* < or > */
#define DISP_BSP_COORD ((Uint4)32768)
		/*use this flag to force a rebuild a new align*/
#define ALIGN_REBUILD ((Uint4)65536)
		/*use this flag to display lower case letters; otherwise upper case is used*/
#define TEXT_LOWERCASE ((Uint4)131072)
		/*use this flag to display 'QUERY:' and 'SBJCT:' as seq. name
		should be used for BLAST output only*/
#define SEQNAME_BLAST_STD ((Uint4)262144)

extern void DDV_DeleteTxtList(ValNodePtr PNTR vnp);
extern void DDV_DeleteParaGList(ValNodePtr PNTR vnp);
extern void DDV_DeleteDisplayList(MsaParaGPopListPtr mpplp);
extern void DDV_ResetParaGSeqAlignCoord(MsaParaGPopListPtr mpplp,Int2 LineSize);
extern void DDV_AffichageParaG(MsaParaGPopListPtr mpplp,Int4 gi,Int4 from,Int4 to,
		Int4 TotalAliLength,Int4 numBlockAffich,Uint4 disp_format,Int2 LineSize,
		FILE *fp,ByteStorePtr PNTR bspp);
extern Boolean DDV_PopDisplay(SABlockPtr sabp,MsaParaGPopListPtr mpplp,
			Int4 nBSP,Uint1 MasterScaleStyle,Boolean ShowTick,
			Int2 cbl,Int4 nBlockToPop);
extern Boolean DDV_PopDisplay2(SABlockPtr sabp,MsaParaGPopListPtr mpplp,
			Int4 nBSP,Uint1 MasterScaleStyle,Boolean ShowTick,
			Int2 cbl,Int4 nBlockToPop,Int4 ali_from,Int4 ali_to);
extern Boolean DDV_GetSequenceFromParaG(ParaGPtr pgp,
		CharPtr PNTR szSequence,Int4 bspLength,Boolean IsAA,Uint1 PNTR bsp_strand,
		Int4 PNTR bsp_start,Int4 PNTR bsp_stop);
extern Boolean DDV_GetFullFASTAforIdxAli(SABlockPtr sabp,Uint4 disp_format,
	Int4 SAsize,FILE *fp);
extern void DDV_PrintStudyName(CharPtr szPopSetName,CharPtr szPopSetAuth,
		Int4 pmid,FILE *fp);
extern void PrintHTML_Header(Uint4 type,FILE *fp);
extern void PrintHTML_Tail(Uint4 type,FILE *fp);
extern void PrintPopContHeader(URLDataPtr udp,FILE *fp);
extern void PrintPopSizeAliTable(Int4 sapLength,Int4 sapNumBioseqs,FILE *fp);
extern void PrintPopRightContent(URLDataPtr udp,Int4 sapLength,Int4 sapNumBioseqs,FILE *fp,
		Int4 numBlockAffich,Int2 LineSize);
extern void PrintPopulationTail(Uint4 type,FILE *fp);
extern void PrintPopulationHeadG(Uint4 type,FILE *fp);
extern void PrintPopulationBar(CharPtr db,FILE *fp);
extern void	DDV_GetArticleInfo(SeqEntryPtr sep,CharPtr szPopSetName,
	CharPtr szPopSetAuth,Int4Ptr pmid);
extern void DDV_GetBspList(SeqEntryPtr sep,FILE *fp);
extern Boolean DDV_DisplayDefaultAlign(SeqAlignPtr sap,Int4 gi,Int4 from,Int4 to,
		Uint4 disp_format,Pointer disp_data,FILE *fp);
extern CharPtr DDV_ReadSeqGivenSpp(SeqPortPtr spp,Int4 from,Int4 to,Uint1 strand,
	Boolean IsAA,BoolPtr bError);
extern SeqPortPtr DDV_OpenBspFullSeqPort(BioseqPtr bsp,Uint1 strand);
extern BspInfoPtr DDV_BspInfoNew(Boolean InitInfo,Boolean InitSeqPort,SeqIdPtr sip,
		Uint1 strand);
extern BspInfoPtr DDV_BspInfoDelete(BspInfoPtr bip);
extern BspInfoPtr DDV_BspInfoDeleteList(BspInfoPtr bip);
extern void DDV_LocateParaG(ValNodePtr PNTR TableHead,Int4 nBsp);
extern Boolean DDV_CreateDisplay(SeqAlignPtr sap,MsaParaGPopListPtr mpplp);
Boolean DDV_ShowSeqAlign(SeqAlignPtr seqalign, Int4 gi, Int4 from, Int4 to,
                         Uint4 disp_format);
extern void	DDV_GetEntryBioSource(SeqEntryPtr sep,ValNodePtr PNTR vnpp);
extern void DDV_DisplayBlastSAP(SeqAlignPtr sap);
extern void DDV_PrintPopSetSummary(SeqEntryPtr sep, Int4 gi, FILE *FileOut);
extern void DDV_SearchAli(SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);
extern void DDV_GetSeqAlign(Pointer dataptr, Uint2 datatype,
		ValNodePtr PNTR vnp);


extern void PrintSeqAlignCallback (SeqEntryPtr sep, Pointer mydata,
                                   Int4 index, Int2 indent);
extern ByteStorePtr SeqAlignToBS (Uint2 entity);

#ifdef __cplusplus
}
#endif

#endif /* ndef _PGPPOP_ */

