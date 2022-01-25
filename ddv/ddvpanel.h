/*  $Id: ddvpanel.h,v 1.10 2000/01/10 15:09:45 durand Exp $
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
* File Name:  ddvpanel.h
*
* Author:  Patrick Durand
*
* Version Creation Date:   06/19/99
*
* $Revision: 1.10 $
*
* File Description: 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ddvpanel.h,v $
* Revision 1.10  2000/01/10 15:09:45  durand
* Use Entrez instead of ID1
*
* Revision 1.9  1999/12/07 21:40:14  durand
* add mouse modes menu and caret facility for the editor
*
* Revision 1.8  1999/12/06 16:19:20  durand
* add GoTo facility to DDV
*
* Revision 1.7  1999/12/03 23:17:23  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 1.6  1999/11/30 17:29:00  durand
* fix a problem of redeclaration of the function DDV_CloseData
*
* Revision 1.5  1999/11/29 15:26:26  durand
* designed a new GUI to fix problems under MacOS, Linux and SGI
*
* Revision 1.4  1999/10/22 20:12:47  durand
* add Export command (text, HTML and Phylip formats)
*
* Revision 1.3  1999/10/20 13:17:19  durand
* add display for disc. SeqAlign tails
*
* Revision 1.2  1999/10/15 21:57:37  durand
* add a UI for display options
*
* Revision 1.1  1999/09/30 14:10:29  durand
* add ddv to toolkit
*
* Revision 1.5  1999/09/09 21:55:07  durand
* instantiate the Fle|Close command of DDV
*
* Revision 1.4  1999/07/01 14:08:08  durand
* the loader functions of DDV or now in ddvopen.c
*
* Revision 1.2  1999/06/28 22:07:22  durand
* add loader functions and clean the code with Lint and Purify
*
* Revision 1.1  1999/06/19 17:21:08  durand
* add Vibrant DDV code
*
*
*
* ==========================================================================
*/

#ifndef _DDVPANEL_
#define _DDVPANEL_

#ifdef __cplusplus
extern "C" {
#endif

#include <udviewer.h>
#include <ddvmain.h>

/******************************************************************************

	Data structures

******************************************************************************/
typedef struct ddvdispstylesmsg {/*used by the Display Styles Dlg box*/
	WindoW hWinMain; /*DDV main window handle*/
	ButtoN  chk2;/*use color display chk box*/
	ButtoN  chk3;/*show left tail chk box*/
	ButtoN  chk4;/*show right tail chk box*/
	GrouP  g1;   /*main group for Highlight unaligned regions"*/
	GrouP  g2;   /*"Spacer" group for Highlight unaligned regions"*/
	GrouP  g5;   /*"seq. justification" group for Highlight unaligned regions"*/
	GrouP  g12;  /*BSP ruler style group*/
	TexT   edit1; /*size TexT */
	LisT   BspNames;/*list of BSP names*/
	} DdvDispStylesMSG, PNTR DdvDispStylesMSGPtr;

typedef struct ddvexporttextemsg {/*used by the Display Styles Dlg box*/
	WindoW hWinMain; /*DDV main window handle*/
	ButtoN  chk2;/*show number*/
	ButtoN  chk3;/*show ticks*/
	ButtoN  chk4;/*Use block of 10 letters*/
	ButtoN  chk5;/*Display strand orientation*/
	ButtoN  chk6;/*Display BioSeq coordinates*/
	ButtoN  ok;
	PopuP   pop; /*list of available formats*/
	TexT    edit1; /*filename */
	} DdvExportTexteMSG, PNTR DdvExportTexteMSGPtr;


/******************************************************************************

	Defines

******************************************************************************/

#define DDV_DEFAULT_PARAG_SIZE 70

/******************************************************************************

	Global varaibles

******************************************************************************/
extern Char szAppName[];


/******************************************************************************

	Exported functions

******************************************************************************/
extern void DDV_EnableGotoTBItems(WindoW hParent,Boolean bEnable);
extern void DDV_WhatSize(DdvMainPtr dmp);
extern void DDV_InitPanelData(UDVPanelDataPtr pdp);
extern void DDV_SetupWin (PaneL p,Boolean bInit,RecT PNTR rcPp);
extern void DDV_Resize_DDV (PaneL Viewer,RecT PNTR rcP);
extern void DDV_WinMainResize (WindoW w);
extern void DDV_VHScrl(PaneL p,UnDViewerGraphDataPtr gdp, Int4 newval, Int4 oldval,
	Boolean IsVscroll);
extern void DDV_VScrlProc (BaR sb, SlatE s, Int4 newval, Int4 oldval);
extern void DDV_HScrlProc (BaR sb, SlatE s, Int4 newval, Int4 oldval);
extern Boolean DDV_CreateViewerPanel(WindoW w,DdvMainWinPtr dmwp,
                                     SAM_ViewGlobal *vgp);
extern void DDV_SetupMenus(WindoW w,Boolean isID1Ok);
extern void DDV_WinMainResize (WindoW w);
extern void DDV_WinMainCleanup (GraphiC g, VoidPtr data);
extern void DDV_WinMainProgQuit(WindoW w);
extern void DDV_InitGraphGlobal(DdvMainPtr dmp);
extern void DDV_CloseData(DdvMainWinPtr mWin_d,Boolean bFinalExit);
extern void DDV_SetRulerAttribInPGP(ValNodePtr ParaG_Head, Uint1 RulerStyle);
extern void DDV_SortPGPLineNum(ValNodePtr PNTR Head, Int4 nBsp);
extern void DDV_CleanupDDVPdata_g (DdvMainPtr dmp);
extern void DDV_TimerProc (WindoW w);

#ifdef __cplusplus
}
#endif

#endif /* ndef _DDVPANEL_ */

