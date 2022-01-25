/*  $Id: ddvcreate.h,v 1.15 1999/12/20 20:20:41 lewisg Exp $
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
* File Name:  ddvcreate.h
*
* Author:  Sarah Wheelan
*
* Version Creation Date:   08/99
*
* $Revision: 1.15 $
*
* File Description: 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ddvcreate.h,v $
* Revision 1.15  1999/12/20 20:20:41  lewisg
* allow cn3d to do color and ddv to do case when both are running
*
* Revision 1.14  1999/12/20 14:45:37  durand
* add new functions for the new BLAST outputs
*
* Revision 1.13  1999/12/08 22:42:17  durand
* add the code to produce colored BLAST outputs
*
* Revision 1.12  1999/12/03 23:17:22  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 1.11  1999/12/03 15:11:49  durand
* add DDV_AllGapInLowerCase declaration
*
* Revision 1.10  1999/11/24 21:26:30  vakatov
* Fixed for the C++ and/or MSVC DLL compilation
*
* Revision 1.9  1999/11/17 22:43:57  durand
* speed up the selection manager for large SeqAlign
*
* Revision 1.8  1999/11/09 17:08:59  durand
* transfer some functions from ddvgraph to ddvcreate, so that ddvcreate remains Vibrant free and can be compiled with BLAST
*
* Revision 1.7  1999/11/03 21:29:48  durand
* add CTRL and SHFT keys for mouse selection. redesign the loader functions of DDV to properly register the message callbacks
*
* Revision 1.6  1999/10/29 14:15:40  durand
* add simple mouse selection functions
*
* Revision 1.5  1999/10/22 20:12:47  durand
* add Export command (text, HTML and Phylip formats)
*
* Revision 1.4  1999/10/20 13:17:19  durand
* add display for disc. SeqAlign tails
*
* Revision 1.3  1999/10/15 21:57:36  durand
* add a UI for display options
*
* Revision 1.2  1999/10/12 15:01:29  lewisg
* resolve confict with internal/ddv
*
* Revision 1.1  1999/09/30 14:10:26  durand
* add ddv to toolkit
*
* Revision 1.5  1999/09/30 13:38:09  durand
* DDV_CreateDisplayFromIndex takes ParaG_Size as an argument
*
*
*
* ==========================================================================
*/
#ifndef _DDVCREATE_
#define _DDVCREATE_

#include <ncbi.h>
#include <alignmgr.h>
#include <udvseq.h>
#include <pgppop.h>
#include <udvdef.h>
#include <ddvcreate.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif


/*******************************************************************************

  structure

*******************************************************************************/
typedef struct uamsg{ /*UnAligned Messaging structure*/
	Int4    from_ua;/*UA coord system (always 0->UA max lenght*/
	Int4    to_ua;
	Int4    from_r; /*BSP coord.*/
	Int4    to_r;
	Int4    CurrentPos; /*for internal use only*/
	Uint1   DispType;
	Uint1   strand;
	Boolean IsGap;
} UAMsg, PNTR UAMsgPtr;

typedef struct descridisp {
        Int4            from;
        Int4            to;
        Int4            UAnum;
        Int4            UAMaxlength;
        Boolean         IsGap;
        Uint1           TextStyle;
        Uint1           strand;
} DescriDisp, PNTR DescriDispPtr;

typedef struct ddv_disp_opt{
	/*color display ?*/
	Boolean bUseColors;
	/*styles for a disc. seqalign*/
	Boolean ShowLeftTail;
	Boolean ShowRightTail;
	Uint1   DispDiscStyle;
	Uint1   SpacerSize;
    Uint1   DiscJustification;
	/*this field is here because it's closely related to the SeqAlign*/
	Uint1   RulerStyle; /*for bioseq coord system; used only when rebuild the display*/
	}DDV_Disp_Opt, PNTR DDV_Disp_OptPtr;

/*******************************************************************************

  defines

*******************************************************************************/
#define RulerStyle_One_Start      ((Uint1)1)
#define RulerStyle_Continue_Start ((Uint1)2)


/*******************************************************************************

  functions

*******************************************************************************/
NLM_EXTERN Boolean DDV_CreateDisplayFromIndex(SeqAlignPtr sap, MsaParaGPopListPtr mpplp, 
		Int2 LineSize, DDV_Disp_OptPtr ddop);
NLM_EXTERN Boolean DDV_CreateDisplay_DiscAlign(SeqAlignPtr sap, 
		MsaParaGPopListPtr mpplp, Int2 LineSize,DDV_Disp_OptPtr ddop);
NLM_EXTERN UAMsgPtr UAMgrIntUAMsg(Int4 from, Int4 to, Uint1 DispType, Uint1 strand);
NLM_EXTERN UAMsgPtr UAMgrDelUAMsg(UAMsgPtr PNTR uamp);
NLM_EXTERN void UAMgrUAMsgGetInfo(UAMsgPtr uamp,Int4Ptr from_bsp, Int4Ptr to_bsp, 
		BoolPtr IsGap);
NLM_EXTERN Boolean UAMgrGetNextUAbit(UAMsgPtr uamp, Int4 MaxLength, Int4 BspLength,
		Int4 BspStart,Int4 BspStop);
NLM_EXTERN ValNodePtr UABuildDescriptor(SeqAlignPtr sap, Int4 nBsp, Int2 LineSize,
	DDV_Disp_OptPtr ddop, Int4Ptr TotLength,Boolean AddLeftUAPart,
	Boolean AddRightUAPart);
NLM_EXTERN void DDV_InitDefSAPdispStyles(DDV_Disp_OptPtr ddop);
NLM_EXTERN void DDV_InitCn3DSAPdispStyles(DDV_Disp_OptPtr ddop);
NLM_EXTERN Uint4Ptr DDV_BuildBspEntitiesTbl(ValNodePtr PNTR TableHead,Int4 nBsp);

/*export functions*/
NLM_EXTERN void DDV_DumpSAPInAFile(MsaParaGPopListPtr mpplp,DDVOptionsBlockPtr dobp, 
		FILE *hFile,Uint4 option, DDV_ColorGlobal *gclr);
NLM_EXTERN void DDV_DumpSAPInFastaFile(MsaParaGPopListPtr mpplp,DDVOptionsBlockPtr dobp, 
		FILE *hFile,Boolean bPrintGap);

/*BLAST stuffs; function to build new outputs*/
NLM_EXTERN Boolean DDV_DisplayBlastSAP(SeqAlignPtr sap, ValNodePtr mask,
	FILE *fp, Boolean is_na, Uint4 tx_option);

/*layout manager stuff*/
NLM_EXTERN Boolean DDV_InitColour_When_Start(SeqAlignPtr sap,MsaParaGPopListPtr mpplp,
	DDV_ColorGlobal * * gclr, Boolean NoColor);
NLM_EXTERN Boolean DDV_LayoutUAregion(SeqAlignPtr sap,MsaParaGPopListPtr mpplp,
	DDV_ColorGlobal * * gclr,Uint1 * pRGB,Uint1 style);
NLM_EXTERN void DDV_LayoutGap(DDV_ColorGlobal *pcg,
	ValNodePtr vnp_seq1,ValNodePtr vnp_seq2,Uint1 style,Uint1 *rgb);
NLM_EXTERN void DDV_LayoutMaskRegions(DDV_ColorGlobal *pcg,
	ValNodePtr vnp_query,ValNodePtr mask,Uint1 style,Uint1 *rgb);
NLM_EXTERN void DDV_LayoutISOColors(DDV_ColorGlobal *pcg,ValNodePtr * row_list,Int4 nRow,
	Int4 Master,Boolean bSetMaster,Int4Ptr * matrix,Uint1 IdentClr,Uint1 SimilClr,
	Uint1 OtherClr);

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* ndef _DDVCREATE_ */

