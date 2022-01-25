/*   undviewer.h
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
* File Name:  undviewer.h
*
* Author:  Patrick Durand
*
* Version Creation Date:   5/3/99
*
* $Revision: 6.3 $
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

#ifndef _UNDVIEWER_
#define _UNDVIEWER_

#ifdef __cplusplus
extern "C" {
#endif


/*******************************************************************************

  INCLUDE SECTION 

*******************************************************************************/
#include <explore.h>
#include <gather.h>
#include <ncbi.h>
#include <objfdef.h>
#include <objseq.h>
#include <objsub.h>
#include <seqmgr.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <vibrant.h>

#include <odlbox.h>
#include <udvseq.h>

/*******************************************************************************

  The major elements of the single sequence viewer : window and paragraph 

*******************************************************************************/

	/*****************************************************************************

	Structure of the UnDViewer area

	-Panel--------------------------------------------------    -
	|         |             |                              |    |
	| Name    | Scale area  |                              |    |
	| area    | (if left    |    Sequence area             |    |
	|         |  scale)     |    (paragraph)               |    cyClient
	|         |             |                              |    |
	|         |             |                              |    |
	--------------------------------------------------------    -
	<-cxName->
        	  <-cxLeftScale->
	<-------------------------cxClient--------------------->

	*****************************************************************************/

	/*****************************************************************************

	Structure of a ParaG : see module udvseq.h

	*****************************************************************************/


/*******************************************************************************

  DEFINE SECTION 

*******************************************************************************/

	/*number of letter by block: default is ten-letter blocks*/
#define LETTER_BLOCK_WIDTH 10

	/*margins around the UnDViewer*/
#define VIEWER_HORZ_MARGIN 10
#define VIEWER_VERT_MARGIN 10

	/*Numerical scale */
#define SCALE_POS_LEFT 1	/*scale position*/
#define SCALE_POS_TOP  2
#define SCALE_POS_BOTH  3
#define SCALE_MAJOR_TICK_EVERY 10	/*position of the ticks*/
#define SCALE_MINOR_TICK_EVERY 5
#define SCALE_LETTER_COLOR (GetColorRGB(0,0,255))	/*scale colours*/
#define SCALE_MAJTICK_COLOR (GetColorRGB(0,0,0))
#define SCALE_MINTICK_COLOR (GetColorRGB(128,128,128))
#define SCALE_WIDTH_IF_LEFT 9 /*nb. of cxChar (see Font struct); used to
                                   the scale on the left*/

	/*Font*/	
#define FONT_DEF_COLOR (GetColorRGB(0,0,0))/*default font colour*/

	/*Panel*/	
#define PANEL_NAME_WIDTH 30 /*nb. of cxChar (see Font struct) MAX: 50 !!;
	                            used to put names pn the left on the UDV window*/

	/*Info Panel*/	
#define BUFFER_SIZE  255 /*nb. of cxChar for Info Buffer*/

	/*NA and AA color layout*/	
#define NA_A_LAYOUT	0
#define NA_G_LAYOUT	1
#define NA_C_LAYOUT	2
#define NA_T_LAYOUT	3
#define NA_U_LAYOUT	4
#define LAYOUT_UPPER_CASE	1
#define LAYOUT_LOWER_CASE	2
	
	/*Mouse actions*/
#define MS_ACTION_FEAT_NOTHING	1 /*no action*/
#define MS_ACTION_FEAT_CURSOR	2 /* double cursor for features*/
#define MS_ACTION_RESIZE_WIN	3 /*resize cxName region*/	

	/*************************************************************************

	  I defined _min_ & _max_ because I needed the test '>=' in both 
	  _max_ and _min_

	*************************************************************************/
#define _min_(a,b)        ((a)>=(b)?(b):(a))
#define _max_(a,b)        ((a)>=(b)?(a):(b))

	/*used to draw features*/
		/*     |===>......===>.....===>|   : a big arrow feature*/
#define FEATURE_START_BOX		1	/*used to draw big arrows*/
#define FEATURE_START_ARROW		2
#define FEATURE_START_ARROW_END	3
#define FEATURE_START_NOTHING	4
		/* /\/\/\/\   : a helix feature (2D structures)*/
#define DRAW_HELIX_DOWN		1  /*used to draw helix faeture*/
#define DRAW_HELIX_MIDDLE	2
#define DRAW_HELIX_UP		3
		/*  >------<  : a bond feature*/
#define BOND_RIGHT 1	/*used to draw a bond feature*/
#define BOND_LEFT  2
#define SZBUF_SIZE 250
	
	/*sequence buffer management ; depending on VScroll up/down*/
#define BUFFER_REPOP_VCRL_LUP	1   /*line up*/
#define BUFFER_REPOP_VCRL_LDN	2   /*line down*/
#define BUFFER_REPOP_VCRL_PUP	3   /*page up*/
#define BUFFER_REPOP_VCRL_PDN	4   /*page down*/
	
/*******************************************************************************

  DATA STRUCTURE SECTION 

*******************************************************************************/
	typedef struct udvscaledata {
		Int2 ScalePosition;	/*position of the scale, see #define above*/
		Int2 cxLeftScale;	/*width of the 'SCALE_POS_LEFT' area, pixel unit
							NULL if ScalePosition== _TOP or if ShowScale==FALSE*/
		Boolean ShowMajorTick;	/*display tick marks if TRUE*/
		Boolean ShowMMinorTick;	/*display tick marks if TRUE*/
		Int2 MajorTickEvery;	/*put a major tick every*/
		Int2 MinorTickEvery;	/*put a minor tick every*/
		Uint4 ScaleColor;	/*letter color*/
		Uint4 TickMajColor;	/*tick color*/
		Uint4 TickMinColor;	/*tick color*/
		} UDVScaleData, PNTR UDVScaleDataPtr;

	typedef struct udvletterlayout {
		Uint4 LetClr;		/*color of a letter*/
		Uint4 bkClr;		/*bk color of a letter*/
		Uint1 AspLet;		/*aspect (uppercase,...)*/
		} UDVLetterLayout, PNTR UDVLetterLayoutPtr;

	typedef struct udvfontdata {
		FonT hFnt;			/*Current Font*/
		Int2 cxChar;		/*width of a char, pixel unit*/
		Int2 cyChar;		/*height of a char, pixel unit*/
		Int2 LineHeight;	/*height of a single data line*/
		Boolean UseDefaultColorLetter; /*if true, use the following color for 
							letter (NA or AA); draft mode*/
		Uint4 LetterColor;
		} UDVFontData, PNTR UDVFontDataPtr;

	typedef struct udvpaneldata {
		/*panel size*/
		Int2 cxClient;		/*width of the client panel area, pixel unit*/
		Int2 cyClient;		/*height of the client panel area, pixel unit*/
		Int2 cxName;		/*width of the 'Name' area, pixel unit*/

		/*Letter lines*/
		Int2 nCharByLine;	/*number of letter displayed on one line of the viewer*/
		Int2 nBlockByLine;	/*number Ten-letter blocks displayed on one line of 
							the viewer*/
		Int4 nTotLines;		/*total number of line to display*/

		/*Display options*/
		Boolean ShowFeatures; /*TRUE: show features*/
		Boolean ShowScale;	/*TRUE: show num. scale*/
		} UDVPanelData, PNTR UDVPanelDataPtr;

	typedef struct udvvscrolldata {
		/*vertical Scroll Bar*/
		Int4 ScrollMax;		/*Size of the vertical Scroll Bar*/
							/*also equal number of data lines displayed within
							UnDViewer*/
		Int4 ScrollPos;		/*Current position within the vertical Scroll Bar*/
		Int4 ScrollPage;	/*Page decal of the vertical Scroll Bar*/
		} UDVVScrollData, PNTR UDVVScrollDataPtr;

	typedef struct undviewergraphdata {
		UDVFontData 	udv_font;	/*Font*/
		UDVPanelData 	udv_panel;	/*Panel*/		
		UDVScaleData	udv_scale;	/*Numerical scale*/
		UDVVScrollData	udv_vscrl;	/*vertical Scroll Bar*/
		Uint4Ptr 		pClr;		/*Color table for Features*/
		UDVLetterLayout NA_LayoutPal[5];
		UDVLetterLayout AA_LayoutPal[26];
		} UnDViewerGraphData, PNTR UnDViewerGraphDataPtr;

	typedef struct udv_item_select {
		Uint2		eIDsel;			/*-|               */
		Uint2		iIDsel;			/* |--selected item*/
		Uint2		iTypeSel;		/*-|               */
		} UDV_Item_Select, PNTR UDV_Item_SelectPtr;
		
	typedef struct udv_mouse_select {
		Int1		Action_type;	/*action with the mouse*/
		RecT		rcClip;			/*limit mvt of the mouse to this rc*/
		PoinT		oldPos;			/*old pos of the mouse*/
		PoinT		newPos;			/*new pos of the mouse*/
		ParaGPtr	pgp;			/*current ParaG data*/
		Char        szPos[20];		/*current position within the bsp*/
		} UDV_mouse_select,PNTR UDV_mouse_selectPtr;

	typedef  struct  scanfeat {/*use to scan feature when the user clicks
						on a feature*/
		Uint2 eID;		/*identity of the feature*/
		Uint2 iID;		
		Int4 index;		/*position within the Feature List Box*/
		} ScanFeat, PNTR ScanFeatPtr;

	typedef struct viewerdialogdata {
		Uint2 				procid;		/*identification of ObjMgr User Data*/
		Uint2 				proctype;	/* idem */
		Uint2 				userkey;	/* idem */
		BspInfo				bsp_i;		/*bsp information*/
		UnDViewerGraphData 	udv_graph;	/*UnDviewer graphical panel info*/
		ValNodePtr 			ParaG;		/*ParaG data of the entire display*/
		PaneL				UnDViewer;	/*handle of the seq viewer panel*/
		PaneL				InfoPanel;	/*handle of the info panel*/
		UDV_mouse_select	UDV_ms;		/*used for mouse manipulation*/
		UDV_Item_Select		Item_select;/*selected Item*/
		UDV_Item_Select		Old_Item_select;/*old selected Item*/
		Boolean				ClickFeatFromDlg;/*click from Feat List Dlg==TRUE*/
		Boolean				AlreadyInit;	/*to avoid multiple init*/
		} ViewerDialogData, PNTR ViewerDialogDataPtr;

	typedef struct udvmainmenu {/*handles of the main menu commands*/
		MenU File;
		MenU Options;
		MenU Help;
		IteM FileOpen;		/*open file command*/
		IteM EntrezOpen;	/*open from Entrez command*/
		IteM FileClose;		/*close file command*/
		IteM ListSequence;	/*list of available sequence*/
		IteM QuitViewer;	/*close the viewer*/
		IteM ShowFeature;	/*show feature on/off command*/
		IteM ShowFeatureList;	/*show feature on/off command*/
		IteM SearchForFeature;/*search features by keyword*/
		MenU ScalePos;		/*scale position command*/
		IteM RefreshScreen;	/*redraw the entire screen*/
		IteM HelpAbout;	/*what's that ?*/
	} UDVMainMenu, PNTR UDVMainMenuPtr;
		
	typedef struct viewermain {
			/*main win menu*/
		UDVMainMenu			MainMenu;		/*menu command list*/
			/*Features List DlgBox*/
		WindoW				hFeatDlg;
			/*Logo Panel*/
		PaneL				Logo_Panel;
		Boolean				Show_logo;
		FonT				f1;
		FonT				f2;
		FonT				f3;
			/*viewer data used only when AutonomeViewer is TRUE */
		Boolean				AutonomeViewer;
		Pointer				dataptr;		/*identification of the object*/
		Uint2				datatype;		/*idem */
		Uint2				entityID;		/*idem */
		ValNodePtr			BspTable;		/*list of bsp*/
		Uint2				BspChoice;		/*selection in menu (ONE-based)*/
		WindoW				hWndMain;
			/*BSP displayed in the viewer*/
		ViewerDialogDataPtr vdp;		
		} ViewerMain, PNTR ViewerMainPtr;

#define REGISTER_UDV_AUTONOMOUS ObjMgrProcLoad(OMPROC_VIEW, \
		"UnD-Viewer", "SingleSeqViewer", OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0, \
		NULL, UDV_ObjRegAutonomous, 0)	
	
	/*struture passed to the FEAT List dialog box*/
	typedef  struct  flmdata {
		WindoW 		hWndMain;	/*UnDviewer*/
		WindoW 		lbox;		/*list of features (OwnerDraw Lbox)*/
		GrouP		gFeatClass;	/*Hidden group of the Feat Class popup*/
		PopuP		pop;	/*available class list*/
		Int1Ptr		SeqFeatClass;/*classes reported for a sequence*/
		} FLMData, PNTR FLMDataPtr;

	/*use to initialize the Feature List Dialog box*/
	typedef  struct  scanfeatforselect {
		Uint2 	eID;	/*if !=-1, UnDviewer slected feature*/		
		Uint2 	iID;		
		Int4 	index;	/*return the index of the UnDviewer slected feature*/
		Int4 	compteur;/*count the number of entries in Feature List DlgBox*/
		WindoW 	lbox;	/*the list itself*/
		Boolean SeqFeatListAvail[SEQFEAT_MAX+1];/*use to retrieve the feature
				classes*/
		Boolean SearchFeat;
		Char    szKeyword[50];
		} ScanFeatForSelect, PNTR ScanFeatForSelectPtr;

	/*
	*  The UdvGlobalsPtr may be registered with a call to SetAppProperty
	*  e.g., SetAppProperty ("UdvGlobals", &udvglobals), where udvglobals
	*  is a persistent structure filled with Vibrant objects or callback
	*  function pointers specific for a given application.
	*/

	typedef SeqEntryPtr (*UdvFetchSeqEntryProc) (Int4 uid, Int2 retcode);

	typedef  struct  udvglobals {
		ViewerMainPtr         vmp;
		UdvFetchSeqEntryProc  fetchSepProc;
		} UdvGlobals, PNTR UdvGlobalsPtr;

/*******************************************************************************

  EXTERNAL FUNCTIONS SECTION (see the .c files of UnDViewer for a complete
	  description of each function)

*******************************************************************************/
	/*drawing and UDV graphical management*/
	NLM_EXTERN void LogoFontCreate(FonT PNTR f1,FonT PNTR f2,FonT PNTR f3);
	NLM_EXTERN FonT UDV_InitFont(UDVFontDataPtr udvfp);
	NLM_EXTERN Int2 UDV_ComputeLineHeight(Int2 cyChar);
	NLM_EXTERN void UDV_FontDim(Int2Ptr cxChar,Int2Ptr cyChar);
	NLM_EXTERN void UDV_ComputePanelSize(RecT rc,Int2Ptr cxClient,
				Int2Ptr cyClient);
	NLM_EXTERN void UDV_ComputeBlockByLine(Int2 cxClient,Int2 cxName,
				Int2 cxLeftScale,Int2 cxChar,Int2Ptr nCharByLine,
				Int2Ptr nBlockByLine);
	NLM_EXTERN void UDV_Init_GraphData(PaneL p, 
											UnDViewerGraphDataPtr GrDataPtr);
	NLM_EXTERN Uint4Ptr UDV_BuildFeatColorTable(void);
	NLM_EXTERN void UDV_Build_NA_LayoutPalette(UnDViewerGraphDataPtr GrData);
	NLM_EXTERN void UDV_Build_AA_LayoutPalette(UnDViewerGraphDataPtr GrData);
	NLM_EXTERN void UDV_select_feature(PaneL p,ViewerDialogDataPtr vdp,
			Uint2 entityID,Uint2 itemID,Boolean bRepos);
	NLM_EXTERN void UDV_ClickMouse(PaneL p,ViewerDialogDataPtr vdp, PoinT pt);
	NLM_EXTERN void UDV_ClickProc(PaneL p, PoinT pt);
	NLM_EXTERN void UDV_DragMouse(PaneL p,ViewerDialogDataPtr vdp, PoinT pt);
	NLM_EXTERN void UDV_DragProc(PaneL p, PoinT pt);
	NLM_EXTERN void UDV_ReleaseMouse(PaneL p,ViewerDialogDataPtr vdp, PoinT pt);
	NLM_EXTERN void UDV_ReleaseProc(PaneL p, PoinT pt);
	NLM_EXTERN void UDV_draw_viewer (PaneL p);
	NLM_EXTERN void UDV_InfoPanelDrawProc(PaneL p);
	NLM_EXTERN void UDV_Logo_onDraw (PaneL p);

	/*Features List Dialog Box*/
	NLM_EXTERN Boolean LIBCALLBACK UDV_FeaturesListBoxFind (SeqFeatPtr sfp, 
			SeqMgrFeatContextPtr context);

	/*UDV main window management*/
	NLM_EXTERN Int2 LIBCALLBACK UDV_ObjRegAutonomous (Pointer data);
	NLM_EXTERN void UDV_WinMainCleanupExtraProc (GraphiC g, VoidPtr data);
	NLM_EXTERN void UDV_Resize_Logo_Panel (WindoW Parent,RecT PNTR rcL);
	NLM_EXTERN void * UDV_FreeVDPstruct(ViewerDialogDataPtr PNTR vdp);
	NLM_EXTERN void UDV_set_PullMenus(UDVMainMenuPtr mmp,Boolean enable);
	NLM_EXTERN void UDV_set_MainMenus(UDVMainMenuPtr mmp,Boolean enable);
	NLM_EXTERN void UDV_SetupMenus(WindoW w,Boolean isEntrezOk);
	NLM_EXTERN void	UDV_resize_viewer(PaneL p,ViewerDialogDataPtr vdp);
	NLM_EXTERN void UnDViewerVScrlUpdate(PaneL p,Boolean bInit,Int4 CurPos);
	NLM_EXTERN void UDV_WinMainResize(WindoW w);
	NLM_EXTERN Boolean  UDV_init_bsp_forViewer(PaneL p,BioseqPtr bsp,Uint2 eID,
			Uint2 iID,Uint2 itype,ViewerDialogDataPtr vdp);
	NLM_EXTERN Boolean UDV_Init_NonAutonomous(PaneL p,
				ViewerDialogDataPtr PNTR vdp,FonT f);

	/*sequence buffer management*/
	NLM_EXTERN CharPtr UDV_Read_Sequence (SeqIdPtr sip, Int4 from, Int4 to, 
		Boolean IsProt,Int2 len);
	NLM_EXTERN void UDV_create_buffer(UnDViewerGraphDataPtr GrData,
		ValNodePtr ParaG_list,BspInfoPtr bsp_i,ValNodePtr StartScan);
		
	/*open a SeqEntry management*/
	NLM_EXTERN Boolean  UDV_analyze_SEP_for_open(FILE *fp,SeqEntryPtr the_set,
		ViewerMainPtr vmp,WindoW w);
	NLM_EXTERN void UDV_Init_vdp_struct(PaneL p,ViewerDialogDataPtr vdp, 
		Boolean EraseParaG,Boolean EraseMainTitle,Boolean EraseInfoPanel);
	NLM_EXTERN void UDV_NetOpen (IteM i);
	NLM_EXTERN void UDV_FileOpen(IteM i);
	NLM_EXTERN void UDV_FileClose(IteM i);
	NLM_EXTERN void UDV_CreateListBioseqDlg(IteM i);

	/*clean quit*/
	NLM_EXTERN void UDV_WinMainProgQuit(WindoW w);




#ifdef __cplusplus
}
#endif

#endif /* ndef _UNDVIEWER_ */

