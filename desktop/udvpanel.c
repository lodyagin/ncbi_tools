/*   udvpanel.c
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
* File Name:  udvpanel.c
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

#include <udviewer.h>

/*local text*/
static Char szAppName[]="UnD-Viewer";

/*local text*/
static Char szAbout[]="UnD-Viewer : A sequence viewer for GenBank\n\
Version beta\n\nInformation Engineering Branch\n\
NCBI - NLM - NIH, Bldg 38A\n\
8600 Rockville Pike\n\
Bethesda, MD 20894 - USA\n\n\
info@ncbi.nlm.nih.gov";

/*struture passed to the FEAT List dialog box; one struct by entry*/
#define COMMENT_SIZE 100
typedef  struct  fldata {
	Int4 num;	/*order num*/		
	Int4 from;	/*position of the feature on the sequence*/
	Int4 to;
	Int2 segments;	/*feature constituted by several seg. (ex: RNA)*/
	Char szType[COMMENT_SIZE];	/*names of the feature*/
	Char szComment[COMMENT_SIZE];
	Uint2 eID;	/*entity ID of this feature*/
	Uint2 iID;	/*item ID of this feature*/
	Uint1 strand;	/*minus, plus, ...*/
	} FLData, PNTR FLDataPtr;

typedef  struct  fksdata {
	WindoW 	hWinMain;
	TexT	keyField;
	ButtoN  bOk;
	}FKSData, PNTR FKSDataPtr;

/*use to initialize the Feature List Dialog box*/
static char *szSeqFeatClassName[]={
				"All",
				"Gene",
				"Organism",
				"CD region",
				"Protein",
				"RNA",
				"Publication",
				"SEQ",
				"IMP",
				"Region",
				"Comment",
				"Bond",
				"Site",
				"RSITE",
				"USER",
				"TXINIT",
				"Numbering",
				"Secondary structure",
				"Non standard residue",
				"Heteroatom",
				"Biosource",
				};

/*local functions*/
static void  UnDViewerVScrlUpdatePage(Int4Ptr PageUpDn,Int2 cyClient,
				Int2 LineHeight);
static void  PanelOrgChange(PaneL p,ViewerDialogDataPtr vdp);
static Boolean  CreateUDVpanel(WindoW w,ViewerMainPtr PNTR vmp);
static void UDV_ReLocateRcParaGList(RecT rcP,Boolean ShowTop,Boolean ShowTick,
			Boolean ShowSequence, Boolean ShowFeatures,
			Boolean ShowBlank,Int4Ptr nTotL,ValNodePtr ParaG_head);
static void ShowFeaturesListDlg(IteM i);
static void UDV_SearchFeatForKey(IteM i);
static void  UDV_analyse_buffer(UnDViewerGraphDataPtr GrData,
		ValNodePtr ParaG_list,BspInfoPtr bsp_i,Uint2 ActionType);

/*******************************************************************************

  Function : UDV_OM_MsgFunc()
  
  Purpose : ObjMgr message loop of the viewer
  
  Parameters : see Toolkit
  
  Return value : see Toolkit 

*******************************************************************************/
static Int2 LIBCALLBACK UDV_OM_MsgFunc (OMMsgStructPtr ommsp)
{
OMUserDataPtr omudp;
/*BioseqPtr bsp;*/
Char buf[41];
ViewerDialogDataPtr vdp;
   
	omudp = (OMUserDataPtr)(ommsp->omuserdata);
	vdp = (ViewerDialogDataPtr)(omudp->userdata.ptrvalue);

	if (!vdp) return(OM_MSG_RET_OK);
	
	BioseqLabel(vdp->bsp_i.bsp, buf, 40, OM_LABEL_BOTH); 
	
	switch (ommsp->message)
	{
		case OM_MSG_DEL:
			break;
		case OM_MSG_CREATE:
			break;
		case OM_MSG_UPDATE:
			/*just for test, when I first implemented ObjMgr*/
			/*Message(MSG_OK, "Got an update message on [%s]", buf);*/
			break;
		case OM_MSG_SELECT:          /* add highlight code */
			/*Click Feat ?*/
			if (ommsp->itemtype==OBJ_SEQFEAT) {
				/*just for test, when I first implemented ObjMgr*/
				/*Message(MSG_OK, "Got a select message on [%s]", buf);*/
				UDV_select_feature(vdp->UnDViewer,vdp,ommsp->entityID,
					ommsp->itemID,
					(Boolean)(vdp->ClickFeatFromDlg ? TRUE : FALSE));
			}
			break;
		case OM_MSG_DESELECT:        /* add deselect code */
			break;
		case OM_MSG_CACHED:
			break;
		case OM_MSG_UNCACHED:
			break;
		case OM_MSG_TO_CLIPBOARD:  /* this is just because no clipboard now */
			break;
		case OM_MSG_SETCOLOR:
			if (vdp->UnDViewer){
				RecT rc;
				ObjectRect(vdp->UnDViewer,&rc);
				InvalRect(&rc);
				Update();
			}
			break;
		default:
			break;
	}

	return OM_MSG_RET_OK;
}

/*******************************************************************************

  Function : UDV_ObjRegAutonomous()
  
  Purpose : call by ObjMgr to start the viewer; Autonomous version
  
  Parameters : see Toolkit
  
  Note: only for the Autonomous viewer purpose
  
  Return value : see Toolkit 

*******************************************************************************/
NLM_EXTERN Int2 LIBCALLBACK UDV_ObjRegAutonomous (Pointer data)
{
WindoW              w=NULL;
OMProcControlPtr    ompcp;
OMUserDataPtr       omudp;
BioseqPtr			bsp=NULL;
ViewerMainPtr 		vmp=NULL;
ViewerDialogDataPtr vdp=NULL;
Char 				szBuf[255]={""};
UdvGlobalsPtr       ugp=NULL;

	/*retrieve data*/
	ompcp = (OMProcControlPtr) data;
	
	if (ompcp == NULL || ompcp->proc == NULL) {
		ErrPostEx (SEV_ERROR, 0, 0, "Data NULL [1]");
		return OM_MSG_RET_ERROR;
	}

	switch (ompcp->input_itemtype) {
		case OBJ_BIOSEQ :
			bsp = (BioseqPtr) ompcp->input_data;
			break;
		case 0 :
			return OM_MSG_RET_ERROR;
		default :
			return OM_MSG_RET_ERROR;
	}
	
	if (bsp == NULL) {
		ErrPostEx (SEV_ERROR, 0, 0, "Data NULL [2]");
		return OM_MSG_RET_ERROR;
	}

	/*initialize the viewer windows (viewer and info panels)*/
	/*autonomous viewer ?*/
	ugp=(UdvGlobalsPtr)GetAppProperty("UdvGlobals");
	if (ugp) vmp=ugp->vmp;

	/*if w and/or vmp are NULL, they'll be initialized in this function*/
	if (vmp) {/*reuse the viewer if already there*/
		w=vmp->hWndMain;
		if (!vmp->vdp || !vmp->vdp->AlreadyInit){
			if (!CreateUDVpanel(w,&vmp)) return OM_MSG_RET_ERROR;
		}
	}
	else{/*for each call, create a new viewer window*/
		/*non-autonomous viewer*/
		if (!CreateUDVpanel(w,&vmp)) return OM_MSG_RET_ERROR;
	}
	
	/*init specific data for the ObjMgr msg loop*/	
	vdp=vmp->vdp;
	if (!vdp) return OM_MSG_RET_ERROR;

	vdp->procid=ompcp->proc->procid;
	vdp->proctype=ompcp->proc->proctype;
	
	if (UDV_init_bsp_forViewer(vmp->vdp->UnDViewer,bsp,ompcp->input_entityID,
		ompcp->input_itemID,ompcp->input_itemtype,vdp)){
		UDV_set_MainMenus(&vmp->MainMenu,TRUE);
	}
	else{
		UDV_set_MainMenus(&vmp->MainMenu,FALSE);
		return OM_MSG_RET_ERROR;
	}

	/*Update dialog box title*/
	sprintf(szBuf,"%s - [%s - %d letters]",szAppName,
				vdp->bsp_i.bspAccNum,
				vdp->bsp_i.bspLength);
	SetTitle(w,szBuf);

	vdp->userkey=OMGetNextUserKey ();
	omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid, 
			OMPROC_VIEW, vdp->userkey);
	if (omudp != NULL) {
		omudp->userdata.ptrvalue = (Pointer) vdp;
		omudp->messagefunc = UDV_OM_MsgFunc;
	}

/*	ObjMgrSendMsg(OM_MSG_UPDATE,ompcp->input_entityID,ompcp->input_itemID,
			ompcp->input_itemtype);*/

	return(OM_MSG_RET_DONE);
}


/*******************************************************************************

  Function : UDV_WinMainProgQuit()
  
  Purpose : end of prog 
  
  Parameters : w; main window handle
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void UDV_WinMainProgQuit(WindoW w)
{
	QuitProgram();
}

/*******************************************************************************

  Function : QuitProc()
  
  Purpose : bye-bye procedure ! 
  
  Parameters : i; menu item which has called this function
  
  Return value : none 

*******************************************************************************/
static void QuitProc(IteM i)
{
ViewerMainPtr 		vmp;
WindoW				hWinMain;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;
	
	Remove(hWinMain);
	if (vmp->AutonomeViewer) QuitProgram();
	
}


/*******************************************************************************

  Function : UDV_FreeVDPstruct()
  
  Purpose : free VDP struct before leaving viewer
  
  Parameters : vdp; data to delete
  
  Return value : NULL

*******************************************************************************/
NLM_EXTERN void * UDV_FreeVDPstruct(ViewerDialogDataPtr PNTR vdp)
{

	if ((*vdp)){
		/*free some user data of Obj Mgr*/
		if ((*vdp)->userkey>0) {
			if ((*vdp)->bsp_i.bsp_entityID>0){
				ObjMgrFreeUserData((*vdp)->bsp_i.bsp_entityID,
						(*vdp)->procid,(*vdp)->proctype,(*vdp)->userkey);
			}
		}

		if ((*vdp)->udv_graph.pClr) MemFree((*vdp)->udv_graph.pClr);
		UDV_FreeListParaG(&(*vdp)->ParaG);
		MemFree ((*vdp));
	}

	return(NULL);
}

/*******************************************************************************

  Function : UDV_WinMainCleanupExtraProc()
  
  Purpose : free memory before leaving viewer
  
  Parameters : vmp; main data
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void UDV_WinMainCleanupExtraProc (GraphiC g, VoidPtr data)
{
ViewerMainPtr vmp=(ViewerMainPtr)data;

	if (!vmp) return;
	/*if (vmp->the_set) SeqEntryFree(vmp->the_set);	*/
	if (vmp->BspTable) ValNodeFree(vmp->BspTable);
 	
	/*done only for Autonomous viewer*/
	if (vmp->dataptr) ObjMgrFree(vmp->datatype,vmp->dataptr);

	UDV_FreeVDPstruct(&vmp->vdp);
	MemFree(vmp);
}


/*******************************************************************************

  Function : UDV_AboutProc()
  
  Purpose : about dialog box
  
  Parameters : -
  
  Return value : none 

*******************************************************************************/
static void  UDV_AboutProc(IteM i)
{
  MsgAlert(KEY_OK, SEV_INFO, "About UnD-Viewer",szAbout);

}

/*******************************************************************************

  Function : UDV_set_PullMenus()
  
  Purpose : manage main window pulldown menus 
  
  Parameters : 	mmp; menu data
  				enable; if TRUE, enable the command menus listed in mmp
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void  UDV_set_PullMenus(UDVMainMenuPtr mmp,Boolean enable)
{
	if (enable){
		Enable(mmp->FileOpen);
		Enable(mmp->EntrezOpen);
		Enable(mmp->FileClose);
		Enable(mmp->ListSequence);
		Enable(mmp->QuitViewer);
		Enable(mmp->ShowFeature);
		Enable(mmp->ShowFeatureList);
		Enable(mmp->SearchForFeature);
		Enable(mmp->RefreshScreen);
		Enable(mmp->HelpAbout);
	}
	else{
		Disable(mmp->FileOpen);
		Disable(mmp->EntrezOpen);
		Disable(mmp->FileClose);
		Disable(mmp->ListSequence);
		Disable(mmp->QuitViewer);
		Disable(mmp->ShowFeature);
		Disable(mmp->ShowFeatureList);
		Disable(mmp->SearchForFeature);
		Disable(mmp->RefreshScreen);
		Disable(mmp->HelpAbout);
	}
}


/*******************************************************************************

  Function : UDV_set_MainMenus()
  
  Purpose : manage main window command menus 
  
  Parameters : 	mmp; menu data
  				enable; if TRUE, enable the command menus listed in mmp
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void  UDV_set_MainMenus(UDVMainMenuPtr mmp,Boolean enable)
{
	if (enable){
		if (mmp->FileClose) Enable(mmp->FileClose);
		Enable(mmp->ShowFeature);
		Enable(mmp->ShowFeatureList);
		Enable(mmp->SearchForFeature);
		Enable(mmp->ScalePos);
		Enable(mmp->RefreshScreen);
		if (mmp->ListSequence) Enable(mmp->ListSequence);
	}
	else{
		if (mmp->FileClose) Disable(mmp->FileClose);
		Disable(mmp->ShowFeature);
		Disable(mmp->ShowFeatureList);
		Disable(mmp->SearchForFeature);
		Disable(mmp->ScalePos);
		Disable(mmp->RefreshScreen);
		if (mmp->ListSequence) Disable(mmp->ListSequence);
	}
}


/*******************************************************************************

  Function : RefreshScreenProc()
  
  Purpose : redraw the screen 
  
  Parameters :  i; menu item
  
  Return value : none 

*******************************************************************************/
static void RefreshScreenProc(IteM i)
{
ViewerMainPtr 		vmp;
WindoW				hWinMain,temport;
ViewerDialogDataPtr vdp;
RecT				rcP;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;

	vdp = vmp->vdp;
	if (vdp == NULL) return;

	if (vdp->UnDViewer==NULL) return;	

	temport=SavePort((WindoW)vdp->UnDViewer);
	Select (vdp->UnDViewer);
	ObjectRect(vdp->UnDViewer,&rcP);
	InvalRect(&rcP);
	RestorePort(temport);
}

/*******************************************************************************

  Function : ScalePosChoiceProc()
  
  Purpose : modification of the scale position 
  
  Parameters :  i; menu item
  
  Return value : none 

*******************************************************************************/
static void ScalePosChoiceProc(ChoicE i)
{
Int2 				value;
ViewerMainPtr 		vmp;
WindoW				hWinMain;

ViewerDialogDataPtr vdp;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;

	vdp = vmp->vdp;
	if (vdp == NULL) return;
	
	value=GetValue(i);
	switch(value){
		case 1:/*Left only*/
			vdp->udv_graph.udv_scale.ScalePosition=SCALE_POS_LEFT;
			vdp->udv_graph.udv_scale.ShowMajorTick=FALSE;
					break;
		case 2:/*Top only*/
			vdp->udv_graph.udv_scale.ScalePosition=SCALE_POS_TOP;
			vdp->udv_graph.udv_scale.ShowMajorTick=TRUE;
			break;
		case 3:/*Top and left*/
			vdp->udv_graph.udv_scale.ScalePosition=SCALE_POS_BOTH;
			vdp->udv_graph.udv_scale.ShowMajorTick=TRUE;
			break;
	}

	PanelOrgChange(vdp->UnDViewer,vdp);
}

/*******************************************************************************

  Function : ShowFeatProc()
  
  Purpose : show/hide features 
  
  Parameters : i; menu item
  
  Return value : none 

*******************************************************************************/
static void ShowFeatProc(IteM i)
{
ViewerMainPtr 		vmp;
WindoW				hWinMain;
Boolean 			ShowFeatures;
ViewerDialogDataPtr vdp;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;

	vdp = vmp->vdp;
	if (vdp == NULL) return;

	ShowFeatures=GetStatus(i);	
	vdp->udv_graph.udv_panel.ShowFeatures=ShowFeatures;

	if (ShowFeatures==FALSE && vmp->hFeatDlg){
		Remove(vmp->hFeatDlg);
		vmp->hFeatDlg=NULL;
		SetStatus(vmp->MainMenu.ShowFeatureList,FALSE);
	}
	if (ShowFeatures==FALSE){
		Disable(vmp->MainMenu.ShowFeatureList);
		Disable(vmp->MainMenu.SearchForFeature);
	}
	else{
		Enable(vmp->MainMenu.ShowFeatureList);	
		Enable(vmp->MainMenu.SearchForFeature);	
	}
	
	PanelOrgChange(vdp->UnDViewer,vdp);
}

/*******************************************************************************

  Function : PanelOrgChange()
  
  Purpose : redraw window after some graph changes (show features, scale,...) 
  
  Parameters : w; handle of the viewer window
  
  Return value : none 

*******************************************************************************/
static void  PanelOrgChange(PaneL p,ViewerDialogDataPtr vdp)
{
/*ViewerDialogDataPtr vdp;*/
RecT 				rcP;
BaR 				vsb;
ValNodePtr 			vnp;
ParaGPtr 			pgp;
Int4 				nLettre=-1;
Int4 				CurPos=0;
Boolean				ShowTop=FALSE,
					ShowTick=FALSE;
WindoW				temport;
/*ViewerMainPtr		vmp;*/

/*	vmp = (ViewerMainPtr) GetObjectExtra (w);
	if (vmp==NULL) return;
	
	vdp=vmp->vdp;

	if (vdp == NULL) return;
	if (vdp->UnDViewer == NULL) return;*/
	if (vdp->ParaG == NULL) return;
	
	temport=SavePort((WindoW)p);

	Select (p);
	ObjectRect(p,&rcP);

	/*current scroll status*/
	vsb = GetSlateVScrollBar ((SlatE) p);
	CurPos=GetBarValue(vsb);
	if (vdp->ParaG){
		for(vnp=vdp->ParaG ; vnp!=NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				if (CurPos>=pgp->StartLine && CurPos<=(pgp->StartLine+
					pgp->nLines)){
					nLettre=pgp->StartLetter;
					break;
				}
			}
		}		
	}

	/*repos. ParaG - compute TotalLines*/
	if (vdp->udv_graph.udv_panel.ShowScale){
		if (vdp->udv_graph.udv_scale.ScalePosition==SCALE_POS_LEFT)
			ShowTop=FALSE;
		else ShowTop=TRUE;

		if (vdp->udv_graph.udv_scale.ShowMajorTick)
			ShowTick=TRUE;
		else ShowTick=FALSE;
	}
	else{
		ShowTop=FALSE;
		ShowTick=FALSE;
	}

		
	UDV_ReLocateRcParaGList(/*vdp->udv_graph,*/rcP,
		ShowTop,ShowTick,TRUE,vdp->udv_graph.udv_panel.ShowFeatures,FALSE,
		&vdp->udv_graph.udv_panel.nTotLines,vdp->ParaG);

	/*compute new Vscroll pos */
	if (vdp->ParaG){
		for(vnp=vdp->ParaG ; vnp!=NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				if (nLettre>=pgp->StartLetter && nLettre<=pgp->StopLetter){
					CurPos=pgp->StartLine;
					break;
				}
			}
		}		
	}

	/*update the Vscroll Bar*/
	UnDViewerVScrlUpdatePage(&(vdp->udv_graph.udv_vscrl.ScrollPage),
		vdp->udv_graph.udv_panel.cyClient,
		vdp->udv_graph.udv_font.LineHeight);
	UnDViewerVScrlUpdate(p,FALSE,CurPos);
	
	/*create new buffer*/
	UDV_create_buffer(&vdp->udv_graph,vdp->ParaG,&vdp->bsp_i,NULL);

	/*Redraw*/
	InvalRect(&rcP);
	RestorePort(temport);
}

/*******************************************************************************

  Function : UDV_SetupMenus()
  
  Purpose : create the menu of the main dialog box 
  
  Parameters : w; handle of the dialog box
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void  UDV_SetupMenus(WindoW w,Boolean isEntrezOk)
{
MenU			m,s;	/*temp variable*/
ChoicE			c;	/* " */
ViewerMainPtr 	vmp;/*program data*/

	vmp = (ViewerMainPtr) GetObjectExtra (w);

	/*File menu*/
	m=PulldownMenu(w,"File");
	vmp->MainMenu.File=m;

	if (vmp->AutonomeViewer){/*available only for the Auntonomous viewer*/
		s=SubMenu(m,"Open from ");
		vmp->MainMenu.FileOpen=CommandItem(s,"a local file...",
				UDV_FileOpen);
		vmp->MainMenu.EntrezOpen=CommandItem(s,"the network...",
				UDV_NetOpen);
		if (isEntrezOk==FALSE) Disable(vmp->MainMenu.EntrezOpen);

		vmp->MainMenu.FileClose=CommandItem(m,"Close",UDV_FileClose);
		SeparatorItem(m);
		vmp->MainMenu.ListSequence=
			CommandItem(m,"Sequence list...",UDV_CreateListBioseqDlg);
		SeparatorItem(m);
	}
	vmp->MainMenu.QuitViewer=CommandItem(m,"Quit/Q",QuitProc);

	/*Options menu*/
	m=PulldownMenu(w,"Options");
	vmp->MainMenu.Options=m;
	vmp->MainMenu.ShowFeature=StatusItem(m,"Show features",
			ShowFeatProc);
	SetStatus(vmp->MainMenu.ShowFeature,TRUE);
	vmp->MainMenu.ShowFeatureList=StatusItem(m,"Features List",
			ShowFeaturesListDlg);
	vmp->MainMenu.SearchForFeature=CommandItem(m,"Search for features",
			UDV_SearchFeatForKey);
	SetStatus(vmp->MainMenu.ShowFeatureList,FALSE);
	SeparatorItem(m);
	vmp->MainMenu.ScalePos=SubMenu(m,"Scale position");

	c=ChoiceGroup(vmp->MainMenu.ScalePos,ScalePosChoiceProc);
	ChoiceItem(c,"Left only");
	ChoiceItem(c,"Top only");
	ChoiceItem(c,"Top and left");
	SetValue(c,3);

	SeparatorItem(m);
	vmp->MainMenu.RefreshScreen=CommandItem(m,"Refresh screen",
			RefreshScreenProc);
	
	/*Help menu*/
	m=PulldownMenu(w,"Help");
	vmp->MainMenu.Help=m;
	vmp->MainMenu.HelpAbout=CommandItem(m,"About.../Q",UDV_AboutProc);
}

/*******************************************************************************

  Function : WinMainMsgProc()
  
  Purpose : Main dialog box answering loop 
  
  Parameters : d; handle of the dialog box
  				mssg; message
  
  Return value : none 

*******************************************************************************/
/*static void WinMainMsgProc (DialoG d, Int2 mssg)

{
ViewerDialogDataPtr vdp;
ViewerMainPtr 		vmp;

	vmp = (ViewerMainPtr) GetObjectExtra (d);
	if (vmp==NULL) return;
	
	vdp=vmp->vdp;

	if (vdp != NULL) {
		switch (mssg) {
			case VIB_MSG_INIT :
				break;
			case VIB_MSG_CLOSE :
				Remove (d);
				break;
			default :
				break;
		}
	}
}
*/
/*******************************************************************************

  Function : UDV_Resize_Logo_Panel()
  
  Purpose : resize logo panel 
  
  Parameters : p; handle of the main window
  				rcP; new rc of the logo panel
  
  Return value : none (see rcP)

*******************************************************************************/
NLM_EXTERN void  UDV_Resize_Logo_Panel (WindoW Parent,RecT PNTR rcL)
{
RecT rcDlg;

	/*size of the parent Window of the Panel*/
	ObjectRect(Parent,&rcDlg);

	/*new size of the Logo Panel*/
	LoadRect (rcL, (Int2)4*VIEWER_HORZ_MARGIN, (Int2)4*VIEWER_VERT_MARGIN, 
					(Int2)(rcDlg.right - rcDlg.left - 4*VIEWER_HORZ_MARGIN), 
					(Int2)(rcDlg.bottom - rcDlg.top - 4*VIEWER_VERT_MARGIN));
}

/*******************************************************************************

  Function : Resize_Panels()
  
  Purpose : resize UnD viewer 
  
  Parameters : p; handle of the UnD viewer window
  				rcP; new rc of the viewer
  
  Return value : none (see rcP)

*******************************************************************************/
static void  Resize_Panels (PaneL Viewer,PaneL Info,RecT PNTR rcP,
		RecT PNTR rcI,Int2 LineHeight)
{
WindoW w;
RecT rcDlg;

	/*size of the parent Window of the Viewer*/
	w=ParentWindow(Viewer);
	ObjectRect(w,&rcDlg);

#ifdef WIN_MAC
	/*New size of the Info panel*/
	if (Info){
		GetPosition(Info,rcI);	
	
		rcI->right=rcDlg.right-rcDlg.left-VIEWER_HORZ_MARGIN;
		rcI->bottom=rcI->top+3*LineHeight/2;
	}
	else rcI->bottom=rcDlg.top+VIEWER_VERT_MARGIN;

	/*new size of the Viewer Panel*/
	GetPosition(Viewer,rcP);	
	rcP->top=rcI->bottom+4;
	rcP->right=rcDlg.right-rcDlg.left-VIEWER_HORZ_MARGIN;
	rcP->bottom=rcDlg.bottom-rcP->top-2*VIEWER_VERT_MARGIN;
#endif

#ifndef WIN_MAC	
	/*New size of the Info panel*/
	if (Info){
		LoadRect(rcI,(Int2)VIEWER_HORZ_MARGIN, 4,
			(Int2)(rcDlg.right - rcDlg.left - VIEWER_HORZ_MARGIN),
			(Int2)(3*LineHeight/2));
	}
	else rcI->bottom=rcDlg.top+VIEWER_VERT_MARGIN;

	/*new size of the Viewer Panel*/
	LoadRect (rcP, (Int2)VIEWER_HORZ_MARGIN, (Int2)(rcI->bottom+3), 
					(Int2)(rcDlg.right - rcDlg.left - VIEWER_HORZ_MARGIN), 
					(Int2)(rcDlg.bottom - rcDlg.top - 2*VIEWER_VERT_MARGIN));
#endif
}

/*******************************************************************************

  Function : UDV_VSscroll()
  
  Purpose : UnD Viewer Vertical Scroll Bar Callback 
  
  Parameters : p; viewer
  				vdp; viewer data structure
				newval; new value of the thumb
				oldval; old value of the thumb
  
  Note : this function MUST to be used by external software using the viewer.
  
  Return value : none 

*******************************************************************************/
static void UDV_VScroll(PaneL p,ViewerDialogDataPtr vdp, Int4 newval, 
			Int4 oldval)
{
RecT 				rcP;
WindoW 				temport;
Int4 				n;
Boolean 			IsGoUp=FALSE;

	if (!vdp) return;
	
	vdp->udv_graph.udv_vscrl.ScrollPos=newval;
	
	temport=SavePort(ParentWindow((WindoW)p));
	Select(p);
	ObjectRect(p,&rcP);
	InsetRect(&rcP,4,4);
	ClipRect(&rcP);

	if (oldval>newval){
		n=oldval-newval;
		IsGoUp=TRUE;
	}
	else {
		n=newval-oldval;
		IsGoUp=FALSE;
	}
	
	if (n<vdp->udv_graph.udv_vscrl.ScrollPage){/*Line UP/Down*/
		/*Create a new Buffer ?*/
		UDV_analyse_buffer(&vdp->udv_graph,vdp->ParaG,&vdp->bsp_i,
			(Uint2)(IsGoUp==TRUE ? BUFFER_REPOP_VCRL_LUP : 
				BUFFER_REPOP_VCRL_LDN));	
		
		/*move and redraw*/
		ScrollRect (&rcP, 0, (Int2)((oldval-newval)*
				vdp->udv_graph.udv_font.LineHeight));
	}
	else{/*Page Up/Down*/
		/*Create a new Buffer ?*/
		UDV_analyse_buffer(&vdp->udv_graph,vdp->ParaG,&vdp->bsp_i,
			(Uint2)(IsGoUp==TRUE ? BUFFER_REPOP_VCRL_PUP : 
				BUFFER_REPOP_VCRL_PDN));	

		/*redraw*/
		/*rcP.left=0;
		rcP.top=0;*/
		InvalRect(&rcP);
	}
	ResetClip();
	RestorePort(temport);
	Update ();
}

/*******************************************************************************

  Function : UnDViewerVScrlProc()
  
  Purpose : UnD Viewer Vertical Scroll Bar Callback 
  
  Parameters : sb; handle of the scroll
  				s; pane;
				newval; new value of the thumb
				oldval; old value of the thumb
  
  Return value : none 

*******************************************************************************/
static void UnDViewerVScrlProc (BaR sb, SlatE s, Int4 newval, 
			Int4 oldval)
{
WindoW 				w;
ViewerDialogDataPtr vdp;
ViewerMainPtr 		vmp;	

	/*get the parent Window of the Panel*/
	w=ParentWindow(s);

	if (w==NULL) return;	
	
	/*get some usefull data...*/
	vmp = (ViewerMainPtr) GetObjectExtra (w);
	if (vmp==NULL) return;
	
	vdp=vmp->vdp;

	UDV_VScroll((PaneL)s,vdp,newval,oldval);
}

/*******************************************************************************

  Function : UnDViewerVScrlUpdatePage()
  
  Purpose : compute PageUp & Down (VScroll) 
  
  Parameters : PageUpDn; initialize here
  				cyClient; height of UnDviewer panel
				LineHeight; height of one data line
  
  Return value : none 

*******************************************************************************/
static void  UnDViewerVScrlUpdatePage(Int4Ptr PageUpDn,Int2 cyClient,
				Int2 LineHeight)
{
	*PageUpDn=(Int4)(cyClient/LineHeight);
}

/*******************************************************************************

  Function : UnDViewerVScrlUpdate()
  
  Purpose : update the UnDviewer VScroll bar 
  
  Parameters : p; handel of the UnDviewer panel
  				bInit; TRUE only when the UnDviewer is created (soft. start up)
  				CurPos ; new position
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void  UnDViewerVScrlUpdate(PaneL p,Boolean bInit,Int4 CurPos)
{
WindoW 				hWnd;
BaR 				vsb;
ViewerDialogDataPtr vdp;
ViewerMainPtr 		vmp;

	/*get the parent Window of the Panel*/
	hWnd=ParentWindow(p);

	if (hWnd==NULL) return;

	/*get some usefull data...*/
	vmp = (ViewerMainPtr) GetObjectExtra (hWnd);
	if (vmp==NULL) return;
	
	vdp=vmp->vdp;

	if (vdp == NULL) return;

	/*current scroll status*/
	vsb = GetSlateVScrollBar ((SlatE) p);
	
	if (bInit) CurPos=0;
	
	/*given curpos (lineNumber) retrieve letter number*/
	
	/*compute new values*/	
	UnDViewerVScrlUpdatePage(&(vdp->udv_graph.udv_vscrl.ScrollPage),
				vdp->udv_graph.udv_panel.cyClient,
				vdp->udv_graph.udv_font.LineHeight);
	if (vdp->udv_graph.udv_panel.nTotLines>vdp->udv_graph.udv_vscrl.ScrollPage){
		vdp->udv_graph.udv_vscrl.ScrollMax=vdp->udv_graph.udv_panel.nTotLines-
							vdp->udv_graph.udv_vscrl.ScrollPage;
	}
	else{
		vdp->udv_graph.udv_vscrl.ScrollMax=0;
		vdp->udv_graph.udv_vscrl.ScrollPos=0;
		vdp->udv_graph.udv_vscrl.ScrollPage=0;
	}
	
	if (CurPos<0) CurPos=0;

	if (CurPos>=vdp->udv_graph.udv_vscrl.ScrollMax)
		vdp->udv_graph.udv_vscrl.ScrollPos=vdp->udv_graph.udv_vscrl.ScrollMax;
	else vdp->udv_graph.udv_vscrl.ScrollPos=CurPos;
				
	/*update scroll*/
	CorrectBarMax(vsb,vdp->udv_graph.udv_vscrl.ScrollMax);
	CorrectBarValue(vsb,vdp->udv_graph.udv_vscrl.ScrollPos);
	CorrectBarPage(vsb,vdp->udv_graph.udv_vscrl.ScrollPage,
			vdp->udv_graph.udv_vscrl.ScrollPage);
	
}

/*******************************************************************************

  Function : UDV_resize_viewer()
  
  Purpose : viewer Resize function 
  
  Parameters : p; handle of the viewer panel
  				vdp; viewer data structure
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void	UDV_resize_viewer(PaneL p,ViewerDialogDataPtr vdp)
{
BaR					vsb;
ValNodePtr			vnp;
ParaGPtr			pgp;
Int4				nLettre=-1;
Int4				CurPos=0;
Boolean				ShowTop=TRUE;
Boolean				ShowTick=TRUE;
RecT 				rcP;
WindoW				temport;

	if (vdp->ParaG == NULL) return;
	
	ObjectRect(p,&rcP);
	temport=SavePort(ParentWindow(p));
	Select(p);
	
	/*Update UnDviewer's panel size values*/
	UDV_ComputePanelSize(rcP,&(vdp->udv_graph.udv_panel.cxClient),
			&(vdp->udv_graph.udv_panel.cyClient));
	
	/*block lines*/
	UDV_ComputeBlockByLine(vdp->udv_graph.udv_panel.cxClient,
			vdp->udv_graph.udv_panel.cxName,
			vdp->udv_graph.udv_scale.cxLeftScale,
			vdp->udv_graph.udv_font.cxChar,
			&vdp->udv_graph.udv_panel.nCharByLine,
			&vdp->udv_graph.udv_panel.nBlockByLine);	

	/*current scroll status*/
	vsb = GetSlateVScrollBar ((SlatE) p);
	CurPos=GetBarValue(vsb);
	if (vdp->ParaG){
		for(vnp=vdp->ParaG ; vnp!=NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				if (CurPos>=pgp->StartLine && CurPos<=(pgp->StartLine+
					pgp->nLines)){
					nLettre=pgp->StartLetter;
					break;
				}
			}
		}		
	}	
	
	/*securite*/
	if (vdp->udv_graph.udv_panel.nBlockByLine<1){
		/*WindoW temport = SavePort(ParentWindow(vdp->UnDViewer));
		Select (vdp->UnDViewer);*/
		InvalRect(&rcP);
		RestorePort (temport);
		Update();
		return;
	}

	/*repos. ParaG - compute TotalLines*/
	if (vdp->udv_graph.udv_panel.ShowScale){
		if (vdp->udv_graph.udv_scale.ScalePosition==SCALE_POS_LEFT)
			ShowTop=FALSE;
		else ShowTop=TRUE;

		if (vdp->udv_graph.udv_scale.ShowMajorTick)
			ShowTick=TRUE;
		else ShowTick=FALSE;
	}
	else{
		ShowTop=FALSE;
		ShowTick=FALSE;
	}

	/*total lines in the viewer*/
	vdp->ParaG=UDV_CreateParaGList(vdp->udv_graph.udv_panel.nCharByLine,
				vdp->bsp_i.bspLength,0,vdp->bsp_i.bspLength-1,
				ShowTop,ShowTick,TRUE,FALSE,&vdp->udv_graph.udv_panel.nTotLines,
				vdp->ParaG);

	if (vdp->ParaG) UDV_PopulateParaGFeatures(/*vdp->udv_graph,*/vdp->bsp_i.bsp,
			vdp->ParaG,vdp->udv_graph.udv_panel.ShowFeatures,
			/*rcP,*/&vdp->udv_graph.udv_panel.nTotLines);

	/*compute new Vscroll pos */
	if (vdp->ParaG){
		for(vnp=vdp->ParaG ; vnp!=NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				if (nLettre>=pgp->StartLetter && nLettre<=pgp->StopLetter){
					CurPos=pgp->StartLine;
					break;
				}
			}
		}		
	}

	/*update the Vscroll Bar*/
	UnDViewerVScrlUpdate(p,FALSE,CurPos);
	
	/*create new buffer*/
	UDV_create_buffer(&vdp->udv_graph,vdp->ParaG,&vdp->bsp_i,NULL);
	InvalRect(&rcP);
	RestorePort(temport);
	Update();
}


/*******************************************************************************

  Function : UDV_WinMainResize()
  
  Purpose : Main Dialog Box Resize function 
  
  Parameters : w; handle of the dialog box
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void  UDV_WinMainResize (WindoW w)
{
ViewerDialogDataPtr vdp;
ViewerMainPtr 		vmp;
RecT				rcP,rcI,rcL;
	
	vmp = (ViewerMainPtr) GetObjectExtra (w);
	if (vmp==NULL) return;

	if (vmp->Show_logo){
		UDV_Resize_Logo_Panel (w,&rcL);
		SetPosition (vmp->Logo_Panel, &rcL );
		AdjustPrnt (vmp->Logo_Panel, &rcL, FALSE);
		return;
	}

	vdp=vmp->vdp;

	if (vdp == NULL) return;

	/*viewer not initialized; typically at the startup of the soft
	 just draw the logo*/
	if (vdp->UnDViewer == NULL) return;

	/*change UnDviewer's panel size*/
	Resize_Panels (vdp->UnDViewer,vdp->InfoPanel,&rcP,&rcI,
		vdp->udv_graph.udv_font.LineHeight);

	if (vdp->ParaG == NULL) {
		SetPosition (vdp->UnDViewer, &rcP );
		AdjustPrnt (vdp->UnDViewer, &rcP, FALSE);
		if (vdp->InfoPanel){
			SetPosition (vdp->InfoPanel, &rcI );
			AdjustPrnt (vdp->InfoPanel, &rcI, FALSE);
		}
		return;
	}

	SetPosition (vdp->UnDViewer, &rcP );
	AdjustPrnt (vdp->UnDViewer, &rcP, FALSE);
	if (vdp->InfoPanel){
		SetPosition (vdp->InfoPanel, &rcI );
		AdjustPrnt (vdp->InfoPanel, &rcI, FALSE);
	}

	UDV_resize_viewer(vdp->UnDViewer,vdp);
}

/*******************************************************************************

  Function : UDV_init_bsp_forViewer()
  
  Purpose : initialize a bsp to display it in the viewer 
  
  Parameters : 	w ; dialog box (parent window of the viewer)
  				bsp; the bsp
				eID, iID, itype; identification of the BSP
				vdp; general data to store ParaG
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN Boolean  UDV_init_bsp_forViewer(PaneL p,BioseqPtr bsp,Uint2 eID,
	Uint2 iID,Uint2 itype,ViewerDialogDataPtr vdp)
{
Boolean		ShowTop=TRUE;
Boolean		ShowTick=TRUE;
Boolean		nRet=FALSE;
RecT 		rcP;
MonitorPtr  mon;
WindoW		temport;

	/*size of UDV_VIEWER*/
	temport=SavePort(ParentWindow(p));
	Select(p);
	ObjectRect(p,&rcP);

	/*store BSP info*/
	vdp->bsp_i.bsp=bsp;
	UDV_ReadBspDataForViewer(&vdp->bsp_i);
	vdp->bsp_i.bsp_entityID=eID;
	vdp->bsp_i.bsp_itemID=iID;
	vdp->bsp_i.bsp_itemType=itype;

	WatchCursor();
	mon = MonitorStrNewEx (szAppName, 30, FALSE);
	MonitorStrValue (mon, "Preparing display...");
	Update ();

	/*Populate ParaG - compute TotalLines*/
	if (vdp->udv_graph.udv_panel.ShowScale){
		if (vdp->udv_graph.udv_scale.ScalePosition==SCALE_POS_LEFT)
			ShowTop=FALSE;
		else ShowTop=TRUE;

		if (vdp->udv_graph.udv_scale.ShowMajorTick)
			ShowTick=TRUE;
		else ShowTick=FALSE;
	}
	else{
		ShowTop=FALSE;
		ShowTick=FALSE;
	}

	vdp->ParaG=UDV_CreateParaGList(vdp->udv_graph.udv_panel.nCharByLine,
		vdp->bsp_i.bspLength,0,vdp->bsp_i.bspLength-1,
		ShowTop,ShowTick,TRUE,FALSE,&vdp->udv_graph.udv_panel.nTotLines,
		vdp->ParaG);

	if (vdp->ParaG) {
		UDV_PopulateParaGFeatures(vdp->bsp_i.bsp,vdp->ParaG,
			vdp->udv_graph.udv_panel.ShowFeatures,
			&vdp->udv_graph.udv_panel.nTotLines);
	}
	else {
		UDV_Init_vdp_struct(p,vdp,FALSE,TRUE,TRUE);
		nRet=FALSE;
		goto fin;
	}

	/*update scroll*/
	UnDViewerVScrlUpdatePage(&(vdp->udv_graph.udv_vscrl.ScrollPage),
		vdp->udv_graph.udv_panel.cyClient,
		vdp->udv_graph.udv_font.LineHeight);
	UnDViewerVScrlUpdate(p,TRUE,0);

	/*create sequence buffer*/
	UDV_create_buffer(&vdp->udv_graph,vdp->ParaG,&vdp->bsp_i,NULL);
	
	/*redraw the viewer*/
	rcP.top=0;rcP.left=0;
	InvalRect(&rcP);
	Update();
	nRet=TRUE;
	
fin:
	RestorePort(temport);
	ArrowCursor();
	MonitorFree (mon);
	return(nRet);
}

/*******************************************************************************

  Function : UDV_init_NonAutonomous()
  
  Purpose : viewer panel creation for external software using the viewer
  
  Parameters : 	vdp; viewer data block - allocated & initialized here
				p; viewer panel - MUST be not NULL
				f; font - MUST be not NULL
				
  Note : 1. the external software MUST initialize p and f before entering this
	  function
		 2. the PaneL p MUST have a correct size : UDV_Init_GraphData() computes
	  cxClient and cyClient used by the viewer to populate ParaG data blocks.
		 		  
  Return value : FALSE if failure 

*******************************************************************************/
NLM_EXTERN Boolean UDV_Init_NonAutonomous(PaneL p,ViewerDialogDataPtr PNTR vdp,
					FonT f)
{
	/*init data structure*/
	(*vdp)=(ViewerDialogDataPtr)MemNew(sizeof(ViewerDialogData));
	if ((*vdp)){
		MemSet((*vdp),0,sizeof(ViewerDialogData));
	}
	else{
		Message (MSG_ERROR, "Viewer creation failed.");
		return(FALSE);
	}

	/*store important data for the viewer*/
	/*(*vdp)->UnDViewer=p;*/
	(*vdp)->udv_graph.udv_font.hFnt=f;

	/*init font size*/
	Select(p);
	SelectFont((*vdp)->udv_graph.udv_font.hFnt);
	UDV_FontDim(&((*vdp)->udv_graph.udv_font.cxChar),
			&((*vdp)->udv_graph.udv_font.cyChar));

	/*compute Feature Color Palette*/
	(*vdp)->udv_graph.pClr=UDV_BuildFeatColorTable();

	/*compute layout letters for NA and AA*/
	UDV_Build_NA_LayoutPalette(&(*vdp)->udv_graph);
	UDV_Build_AA_LayoutPalette(&(*vdp)->udv_graph);			

	(*vdp)->udv_graph.udv_font.LineHeight=
				UDV_ComputeLineHeight((*vdp)->udv_graph.udv_font.cyChar);
	
	/*general graph data init; scale pos, etc.*/
	UDV_Init_GraphData((*vdp)->UnDViewer,&(*vdp)->udv_graph);
	
	(*vdp)->AlreadyInit=TRUE;
	
	return(TRUE);
}

/*******************************************************************************

  Function : CreateUDVpanel()
  
  Purpose : viewer panel creation ; only graphic strctures and panel itself
  
  Parameters : w; parent main window
				vmp; main data structure
				  
  Return value : FALSE if failure 

*******************************************************************************/
static Boolean  CreateUDVpanel(WindoW w,ViewerMainPtr PNTR vmp)
{
ViewerDialogDataPtr vdp;
RecT				rcP,rcI;
CharPtr				szTexte;
WindoW				hParent;



	/*w==NULL; we have to build a child viewer*/
	/*w!=NULL; we deal with an autonomous viewer*/
	if (w==NULL){
		Int2 Margins;
		
		*vmp=(ViewerMainPtr)MemNew(sizeof(ViewerMain));
		if (*vmp){
			MemSet(*vmp,0,sizeof(ViewerMain));
		}
		else{
			Message (MSG_ERROR, "Viewer creation failed.");
			return(FALSE);
		}

		Margins=4*stdCharWidth;
		hParent=DocumentWindow(Margins,Margins ,
				(Int2)((screenRect.right-screenRect.left)-2*Margins), 
				(Int2)((screenRect.bottom-screenRect.top)-2*Margins), 
				szAppName, 
				NULL,
				UDV_WinMainResize);

		if (hParent==NULL){
			Message (MSG_ERROR, "Viewer creation failed.");
			MemFree((*vmp));
			return(FALSE);
		}

		(*vmp)->hWndMain=hParent;
		SetObjectExtra (hParent, (Pointer) *vmp, (FreeProc)UDV_WinMainCleanupExtraProc);

		(*vmp)->AutonomeViewer=FALSE;

		UDV_SetupMenus(hParent,FALSE);

		UDV_set_MainMenus(&((*vmp)->MainMenu),FALSE);
	}
	else hParent=w;


/*	if (w) hParent=w;
	else return(FALSE);
*/
	if (!(*vmp)->vdp){
		(*vmp)->vdp=(ViewerDialogDataPtr)MemNew(sizeof(ViewerDialogData));
		if ((*vmp)->vdp){
			MemSet((*vmp)->vdp,0,sizeof(ViewerDialogData));
		}
		else{
			Message (MSG_ERROR, "Viewer creation failed.");
			Remove(hParent);			
			return(FALSE);
		}
	}
	
	vdp=(*vmp)->vdp;

	/*create info panel*/
	vdp->InfoPanel=AutonomousPanel(hParent,10,10,UDV_InfoPanelDrawProc,
			NULL,NULL,0,NULL,NULL);

	if (vdp->InfoPanel==NULL) return(FALSE);

	szTexte=(CharPtr)MemNew((BUFFER_SIZE+1)*sizeof(char));

	if (szTexte) MemSet(szTexte,0,BUFFER_SIZE+1);
	else return(FALSE);

	SetObjectExtra (vdp->InfoPanel, (Pointer) szTexte, StdCleanupExtraProc);

	/*create viewer panel*/
	vdp->UnDViewer=AutonomousPanel4(hParent,10,10,UDV_draw_viewer,
			UnDViewerVScrlProc,
			NULL,0,NULL,NULL);
	if (vdp->UnDViewer==NULL) return(FALSE);
	SetPanelClick(vdp->UnDViewer,UDV_ClickProc,
		UDV_DragProc,NULL,UDV_ReleaseProc);
	/*Dialog Box CallBack*/
	/*vdp->dialogmessage=(DialogMessageFunc)WinMainMsgProc;*/

	/*FonT & FonT size*/
	if (!UDV_InitFont(&(vdp->udv_graph.udv_font))) return(FALSE);
	Select(hParent);
	SelectFont(vdp->udv_graph.udv_font.hFnt);
	UDV_FontDim(&(vdp->udv_graph.udv_font.cxChar),
			&(vdp->udv_graph.udv_font.cyChar));

	/*compute Feature Color Palette*/
	vdp->udv_graph.pClr=UDV_BuildFeatColorTable();

	/*compute layout letters for NA and AA*/
	UDV_Build_NA_LayoutPalette(&vdp->udv_graph);
	UDV_Build_AA_LayoutPalette(&vdp->udv_graph);			

	vdp->udv_graph.udv_font.LineHeight=
				UDV_ComputeLineHeight(vdp->udv_graph.udv_font.cyChar);
	
	/*adjust panel size*/	
	Resize_Panels (vdp->UnDViewer,vdp->InfoPanel,&rcP,&rcI,
		vdp->udv_graph.udv_font.LineHeight);
	if (vdp->InfoPanel){
		SetPosition (vdp->InfoPanel, &rcI );
		AdjustPrnt (vdp->InfoPanel, &rcI, FALSE);
	}
	SetPosition (vdp->UnDViewer, &rcP );
	AdjustPrnt (vdp->UnDViewer, &rcP, FALSE);

	/*hide the logo panel; only fot the autonomous viewer variant*/
	if ((*vmp)->AutonomeViewer) Hide((*vmp)->Logo_Panel);

	/*general graph data init; scale pos, etc.*/
	UDV_Init_GraphData(vdp->UnDViewer,&vdp->udv_graph);
	
	vdp->AlreadyInit=TRUE;

	RealizeWindow(hParent);
	Show(hParent);
	Update();
	
	return(TRUE);
}



/*******************************************************************************

  Function : UDV_create_buffer()
  
  Purpose : create a sequence buffer for a specific bioseq 
  
  Parameters : vdp; general data
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void  UDV_create_buffer(UnDViewerGraphDataPtr GrData,
		ValNodePtr ParaG_list,BspInfoPtr bsp_i,ValNodePtr StartScan)
{
ParaGPtr 	pgp=NULL;
ValNodePtr 	vnp=NULL,vnp2=NULL;
Int4 		start_seq=0,stop_seq=0,stop=0,nLength=0,StartBuf=0,StopBuf=0;
SeqIdPtr 	sip;
	
	if (StartScan==NULL) vnp=ParaG_list;
	else vnp=StartScan;
	
	/*find the ParaG where ScrollPos is Located*/
	/*Find the end of the display; last PGP currently displayed*/
	if (GrData->udv_vscrl.ScrollMax){/*Vscroll ?*/
		for( ; vnp != NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				if ((GrData->udv_vscrl.ScrollPos>=pgp->StartLine) &&
					(GrData->udv_vscrl.ScrollPos<=(pgp->StartLine+pgp->nLines))){
						start_seq=pgp->StartLetter;
						break;
				}
			}		
		}
	}	
	else{/*no Vscroll: buffer starts at the beginning of the Bioseq*/
		start_seq=0;
	}

	/*Find the end of the display; last PGP currently displayed*/
	if (GrData->udv_vscrl.ScrollMax){/*Vscroll ?*/
		stop=MIN(GrData->udv_panel.nTotLines,
			GrData->udv_vscrl.ScrollPos+GrData->udv_vscrl.ScrollPage);
		for(vnp2=vnp ; vnp2 != NULL ; vnp2=vnp2->next){
			pgp=(ParaGPtr)vnp2->data.ptrvalue;
			if (pgp->StartLine+pgp->nLines>=stop){
				stop_seq=pgp->StopLetter;
				break;
			}
		}
	}
	else{/*no Vscroll: buffer ends at the end of the Bioseq*/
		stop_seq=bsp_i->bspLength;
	}

	/*nb of char within a panel*/
	nLength=stop_seq-start_seq+1;

	/*read sequence; buffer size: three panels*/
	StartBuf=MAX(0,start_seq-nLength);
	StopBuf=MIN(bsp_i->bspLength-1,StartBuf+3*nLength);
	
	/*if same buffer, return*/
	if (bsp_i->StartBuf==StartBuf && bsp_i->StopBuf==StopBuf) return;

	bsp_i->StartBuf=StartBuf;
	bsp_i->StopBuf=StopBuf;
	bsp_i->LengthBuf=bsp_i->StopBuf-bsp_i->StartBuf+1;

	if (bsp_i->SeqBuf) MemFree(bsp_i->SeqBuf);

	sip=SeqIdFindBest(bsp_i->bsp->id,0);
	if (sip==NULL) sip=bsp_i->bsp->id;

	bsp_i->SeqBuf=UDV_Read_Sequence (sip,bsp_i->StartBuf,bsp_i->StopBuf, 
		(Boolean)(!bsp_i->bspMolNuc),bsp_i->LengthBuf);

/*printf("%d; %d ; %d ; %d ; %s\n\n",bsp_i->bspMolNuc,bsp_i->StartBuf,
bsp_i->StopBuf,bsp_i->LengthBuf,bsp_i->SeqBuf);*/

	/*find the ParaG where bsp_i->StartBuf is Located*/
	if (StartScan==NULL) vnp=ParaG_list;
	else vnp=StartScan;
	for( ; vnp != NULL ; vnp=vnp->next){
		if (vnp->data.ptrvalue){
			pgp=(ParaGPtr)vnp->data.ptrvalue;
			if (pgp->StartLetter==bsp_i->StartBuf){
					bsp_i->PgpStartBuf=vnp;
					break;
			}
		}		
	}
}


/*******************************************************************************

  Function : UDV_analyse_buffer()
  
  Purpose : adjust a buffer in response to VScroll action 
  
  Parameters : vdp; general data
  
  Return value : none 

*******************************************************************************/
static void  UDV_analyse_buffer(UnDViewerGraphDataPtr GrData,
		ValNodePtr ParaG_list,BspInfoPtr bsp_i,Uint2 ActionType)
{
ParaGPtr pgp=NULL;
ValNodePtr vnp=NULL,vnp2=NULL;
Int4 start_seq=0,stop_seq=0,stop=0;

	if (ActionType==BUFFER_REPOP_VCRL_LUP ||
		ActionType==BUFFER_REPOP_VCRL_PUP){/*Vscroll  Up*/
	
		if (bsp_i->PgpStartBuf==NULL) vnp=ParaG_list;
		else vnp=bsp_i->PgpStartBuf;

		/*find the ParaG where ScrollPos is Located*/
		/*Find the end of the display; last PGP currently displayed*/
		for( ; vnp != NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				if ((GrData->udv_vscrl.ScrollPos>=pgp->StartLine) &&
					(GrData->udv_vscrl.ScrollPos<=(pgp->StartLine+pgp->nLines))){
						start_seq=pgp->StartLetter;
						break;
				}
			}		
		}
		
		if (start_seq<=bsp_i->StartBuf){
			UDV_create_buffer(GrData,ParaG_list,bsp_i,NULL);
			return;
		}
	}

	if (ActionType==BUFFER_REPOP_VCRL_LDN ||
		ActionType==BUFFER_REPOP_VCRL_PDN){/*Vscroll  Down*/
	
		if (bsp_i->PgpStartBuf==NULL) vnp=ParaG_list;
		else vnp=bsp_i->PgpStartBuf;

		/*find the last ParaG located in the panel*/
		
		stop=MIN(GrData->udv_panel.nTotLines,
			GrData->udv_vscrl.ScrollPos+GrData->udv_vscrl.ScrollPage);
		for(vnp2=vnp; vnp2 != NULL ; vnp2=vnp2->next){
			pgp=(ParaGPtr)vnp2->data.ptrvalue;
			if (pgp->StartLine+pgp->nLines>=stop){
				stop_seq=pgp->StopLetter;
				break;
			}
		}
		if (stop_seq>=bsp_i->StopBuf-GrData->udv_panel.nCharByLine){
			UDV_create_buffer(GrData,ParaG_list,bsp_i,vnp);
		}
	}
}


/*******************************************************************************

  Function : UDV_ReLocateRcParaGList()
  
  Purpose : relocate RecT (position) of each ParaG in the list. This function
  			is generally called in response to a Scale change (top, left, etc)
			or to Show Features change (on/off)
  
  Parameters : 	GrData; graphical data
				rcP; size and position of the UnD viewer's panel
				ShowTop; show scale on top if TRUE
				ShowTick; show scale's ticks if TRUE
				ShowSequence; show sequence if TRUE
				ShowFeature; show features if TRUE
				ShowBlank; add blank line at the bottom of the ParaG
				nTotL; total number of single lines for the whole ParaG list
				ParaG_head; head of the ParaG list

  Note : This function
  		is generally called in response to a Scale change (top, left, etc)
		or to Show Features change (on/off).		

  Return value : none (see nTotL)

*******************************************************************************/
static void UDV_ReLocateRcParaGList(/*UnDViewerGraphData GrData,*/
			RecT rcP,
			Boolean ShowTop,Boolean ShowTick,Boolean ShowSequence, 
			Boolean ShowFeatures,
			Boolean ShowBlank,Int4Ptr nTotL,ValNodePtr ParaG_head)
{
ValNodePtr 	vnp;					/*just for ParaG list scanning*/
Int1		minLineByParaG=0;		/*min height of a paraG, nb. of lines*/			
Int4		nTotLines=0;			/*Total lines in viewer*/
ParaGPtr 	pgp=NULL;				/*ParaG data*/

	if (ParaG_head == NULL) {*nTotL=0;return;}

	/*Minimum ParaG height*/
	if (ShowTop) minLineByParaG++;
	if (ShowTick) minLineByParaG++;
	if (ShowSequence) minLineByParaG++;
	if (ShowBlank) minLineByParaG++;

	/*modify the RecT of each ParaG*/
	for(vnp=ParaG_head ; vnp!=NULL ; vnp=vnp->next){
		pgp=(ParaGPtr)vnp->data.ptrvalue;
		/*modify in each ParaG just the needed*/
		pgp->StartLine=nTotLines;
		pgp->nLines=minLineByParaG;

		if (ShowFeatures) {
			/*nLines=pgp->nFeat+pgp->nTrans;*/
			pgp->nLines+=pgp->nFeatLines;/*nLines;*/
		}

		/*modify values*/
		if (ShowFeatures) nTotLines+=pgp->nLines;
		else nTotLines+=minLineByParaG;
	}

	*nTotL=nTotLines;
}

/*******************************************************************************

Function: FeaturesListBoxScan()

Purpose: scan features to retrieve the index of a seleted feature. This callback
		is used when the user clicks a feature on the UnDviewer panel while the
		Feature List Dialog Box is opened.

Parameters: see explore.h

*******************************************************************************/
static Boolean LIBCALLBACK FeaturesListBoxScan (SeqFeatPtr sfp, 
			SeqMgrFeatContextPtr context)
{
ScanFeatPtr sf=(ScanFeatPtr)context->userdata;

	if(sf){
		sf->index++;
		if (sf->eID==context->entityID && sf->iID==context->itemID){
			return(FALSE);
		}
	}
	return(TRUE);
}



/*******************************************************************************

Function: FLB_OnDraw()

Purpose: Owner Draw Features List box drawing procedure.

Parameter: mid; entry data to draw.

Return value: -

*******************************************************************************/
static void FLB_OnDraw(MIDataPtr mid)
{
RecT 		rc;
Char 		szTmp[255]={""};
FLDataPtr 	fldp=(FLDataPtr)mid->Data;

	if (fldp==NULL) return;
	
	rc=mid->rcItem;
	InsetRect(&rc,3,1);
	rc.left+=mid->cxChar/2;
	/*mid->rcItem.right=rc.right;*/
	/*selection ?*/	
	if (mid->Selected){
		SetColor(GetColorRGB(0,0,128));
		PaintRect(&rc);
		White();
	}
	else{
		White();
		PaintRect(&rc);
		Black();
	}
	mid->rcItem.bottom-=(mid->cyChar/3);
	/*draw features elements*/
	SelectFont(mid->f);
	MoveTo(mid->rcItem.left,(Int2)(mid->rcItem.bottom-mid->cyChar));
	sprintf(szTmp,"%5d.",fldp->num);
	PaintString(szTmp);
	if (mid->Selected) White();
	else Blue();
	MoveTo((Int2)(mid->rcItem.left+7*mid->cxChar),
			(Int2)(mid->rcItem.bottom-mid->cyChar));
	PaintString("Name :");
	if (!mid->Selected) Black();
	sprintf(szTmp,"%s",fldp->szComment);
	MoveTo((Int2)(mid->rcItem.left+14*mid->cxChar),
			(Int2)(mid->rcItem.bottom-mid->cyChar));
	PaintString(szTmp);

	if (fldp->segments>1){
		sprintf(szTmp,"%s [from %d to %d in %d segments]. Length: %d.",
			fldp->szType,fldp->from+1,fldp->to+1,
			fldp->segments,fldp->to-fldp->from+1);
	}
	else{
		sprintf(szTmp,"%s [from %d to %d]. Length: %d.",
			fldp->szType,fldp->from+1,fldp->to+1,fldp->to-fldp->from+1);
	}
	if (!mid->Selected) SetColor(GetColorRGB(0,128,0));
	MoveTo((Int2)(mid->rcItem.left+8*mid->cxChar),mid->rcItem.bottom);
	PaintString("Info :");
	if (!mid->Selected) Black();
	
	MoveTo((Int2)(mid->rcItem.left+15*mid->cxChar),mid->rcItem.bottom);
	PaintString(szTmp);
	LtGray();
	mid->rcItem.bottom+=(mid->cyChar/3);
	MoveTo(rc.left,mid->rcItem.bottom);	
	LineTo(rc.right,mid->rcItem.bottom);
	Black();
}

/*******************************************************************************

Function: FeatListlBoxProc()

Purpose: Owner Draw Features List box callback function.

Parameter: mid; entry data to draw.

Return value: -

*******************************************************************************/
static Boolean FeatListlBoxProc(WindoW lbox,Int2 msg,Pointer data)
{
MIDataPtr mid=(MIDataPtr)data;

	if (mid==NULL) return(TRUE);
	
	switch(msg){
		case ODLB_DrawItem:{
			/*WindoW temport;
			temport=SavePort(ParentWindow(lbox));
			Select(lbox);
			ClipRect(&mid->rcP);*/
			FLB_OnDraw(mid);
			/*ResetClip();
			RestorePort(temport);*/
			break;
		}
		case ODLB_MeasureItem:{
			
			mid->rcItem.left=mid->rcP.left;
			mid->rcItem.right=mid->rcP.right;
			mid->rcItem.top=0;
			mid->rcItem.bottom=5*mid->cyChar/2;
			
			break;
		}
		case ODLB_SelectItem:{
			ViewerMainPtr 	vmp;
			FLMDataPtr		pflm;
			FLDataPtr 		fldp=(FLDataPtr)mid->Data;

			pflm=(FLMDataPtr)GetObjectExtra(mid->Parent);
			if (pflm==NULL) break;
			vmp=(ViewerMainPtr)GetObjectExtra(pflm->hWndMain);
			if (vmp==NULL) break;
			if (fldp==NULL) break;
			vmp->vdp->ClickFeatFromDlg=TRUE;
			ObjMgrSendMsg(OM_MSG_SELECT,fldp->eID,fldp->iID,OBJ_SEQFEAT);
			/*UDV_select_feature(vmp->vdp,fldp->eID,fldp->iID,TRUE);*/
			break;
		}
	}
	return(TRUE);
}

/*******************************************************************************

  Function : FeatListDlgQuit()
  
  Purpose : end of prog Features List Dialog Box
  
  Parameters : w; Dialog Box window handle
  
  Return value : none 

*******************************************************************************/

static void FeatListDlgQuit(WindoW w)
{
FLMDataPtr		pflm;
ViewerMainPtr 	vmp;
	
	pflm = (FLMDataPtr) GetObjectExtra (w);
	if (!pflm) return;

	vmp = (ViewerMainPtr) GetObjectExtra (pflm->hWndMain);

	if (vmp==NULL) return;
	
	SetStatus(vmp->MainMenu.ShowFeatureList,FALSE);
	Remove(w);
	vmp->hFeatDlg=NULL;
}


/*******************************************************************************

  Function : FeatListDlgProc()
  
  Purpose : Features List Dialog Box Resize function 
  
  Parameters : w; handle of the dialog box
  
  Return value : none 

*******************************************************************************/
static void FeatListDlgProc (WindoW w)
{
RecT			rcDlg,rc;
FLMDataPtr		pflm;

	Select(w);
	ObjectRect(w,&rcDlg);
	
	pflm = (FLMDataPtr) GetObjectExtra (w);
	if (!pflm) return;
	if (pflm->gFeatClass){
		ObjectRect(pflm->gFeatClass,&rc);
		rc.bottom+=3;
	}
	else rc.bottom=10;
	
	OwnerDrawLbox_MoveWindow(pflm->lbox,10,rc.bottom,
		(Int2)(rcDlg.right-rcDlg.left-20),
		(Int2)(rcDlg.bottom-rcDlg.top-10-rc.bottom));
}

/*******************************************************************************

Function: UDV_FeaturesListBoxFind()

Purpose: feature selection

Parameters: see explore.h

*******************************************************************************/
NLM_EXTERN Boolean LIBCALLBACK UDV_FeaturesListBoxFind (SeqFeatPtr sfp, 
			SeqMgrFeatContextPtr context)
{
ScanFeatForSelectPtr sffs;


	sffs=(ScanFeatForSelectPtr)context->userdata;
	
	sffs->compteur++;

	if (sffs->eID==context->entityID && sffs->iID==context->itemID){
		sffs->index=sffs->compteur;
		return(FALSE);
	}
		
	return(TRUE);
}

/*******************************************************************************

Function: FeaturesListBoxInit()

Purpose: retrieve features for Features List Box

Parameters: see explore.h

*******************************************************************************/
static Boolean LIBCALLBACK FeaturesListBoxInit (SeqFeatPtr sfp, 
			SeqMgrFeatContextPtr context)
{
FLDataPtr fldp;
ScanFeatForSelectPtr sffs;

/* for test only
if (context->seqfeattype==SEQFEAT_HET){
 Char  str [256];
 int i;

  if (FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
    printf ( "   Feature item %d index %d (%d) (%d - %d) (%d) %s\n",
             (int) context->itemID, (int) context->index, 
                         (int) context->numivals,(int) context->left,
                         (int) context->right,(int) context->strand,str);
        if (context->numivals>1){
                for (i=0;i<context->numivals*2;i+=2){
                        printf ( "   Ivals: (%d - %d)\n", context->ivals[i],
                                                context->ivals[i+1]);
                }
        }
  }
	
}
*/
	sffs=(ScanFeatForSelectPtr)context->userdata;
	
	/*keyword search only*/
	if (sffs->SearchFeat==TRUE){
		Char szType[COMMENT_SIZE];
		
		FeatDefLabel(sfp, szType, sizeof (szType) - 1, OM_LABEL_BOTH);
		StringUpper(szType);
		if (!StringStr(szType,sffs->szKeyword)){
			if (context->sfp->comment){
				MemCopy(szType,context->sfp->comment,
					MIN(COMMENT_SIZE-1,StringLen(context->sfp->comment)));
				StringUpper(szType);
				if (!StringStr(szType,sffs->szKeyword)) return(TRUE);
			}
			else return(TRUE);
		}
	}

	sffs->compteur++;
	fldp=(FLDataPtr)MemNew(sizeof(FLData));
	if (fldp==NULL) return(TRUE);
	fldp->num=sffs->compteur;/*context->index;*/
	fldp->from=context->left;
	fldp->to=context->right;
	fldp->segments=context->numivals;
	fldp->eID=context->entityID;
	fldp->iID=context->itemID;
	fldp->strand=context->strand;
	FeatDefLabel(sfp, fldp->szType, sizeof (fldp->szType) - 1, OM_LABEL_BOTH);
	if (context->sfp->comment){
		MemCopy(fldp->szComment,context->sfp->comment,
		MIN(COMMENT_SIZE-1,StringLen(context->sfp->comment)));
	}
	OwnerDrawLbox_AddItem(sffs->lbox,(Pointer) fldp);

	if (sffs->eID==context->entityID && sffs->iID==context->itemID){
		sffs->index=sffs->compteur;
	}
	
	sffs->SeqFeatListAvail[context->seqfeattype]=TRUE;
	
	return(TRUE);
}

/*******************************************************************************

  Function : FeatClassChooseProc()
  
  Purpose : repopulate the Feature list box when the nice user choose various
  			feature classes
  
  Parameters : pop; popup list of feat classes
  
  Return value : none 

*******************************************************************************/
static void FeatClassChooseProc (PopuP pop)
{
WindoW 				w,temport;
Int1 				Choice;
ViewerMainPtr 		vmp;
ScanFeatForSelect 	sffs;
Boolean 			SeqFeatListAvail[SEQFEAT_MAX];
FLMDataPtr			pflm;
RecT				rcP;

	w=(WindoW)ParentWindow(pop);

	if (w==NULL) return;
	pflm = (FLMDataPtr) GetObjectExtra (w);
	if (!pflm) return;
	vmp = (ViewerMainPtr) GetObjectExtra (pflm->hWndMain);
	if (!vmp) return;
	if (!pflm->SeqFeatClass) return;
	Choice=pflm->SeqFeatClass[GetValue(pflm->pop)-1];
	/*reset lbox Content*/
	OwnerDrawLbox_ResetLBContent(pflm->lbox);
	
	/*repopulate lbox with the user Choice*/
	sffs.eID=vmp->vdp->Item_select.eIDsel;		
	sffs.iID=vmp->vdp->Item_select.iIDsel;		
	sffs.index=0;
	sffs.compteur=0;
	sffs.lbox=pflm->lbox;
	sffs.SearchFeat=FALSE;
	sffs.szKeyword[0]='\0';
	MemSet(sffs.SeqFeatListAvail,0,sizeof(sffs.SeqFeatListAvail));
	if (Choice==0){
		MemSet(SeqFeatListAvail,1,sizeof(SeqFeatListAvail));
	}
	else{
		MemSet(SeqFeatListAvail,0,sizeof(SeqFeatListAvail));
		SeqFeatListAvail[Choice]=TRUE;
	}	
	SeqMgrExploreFeatures (vmp->vdp->bsp_i.bsp, (Pointer) &sffs, 
				FeaturesListBoxInit, NULL, SeqFeatListAvail, NULL);
	if (sffs.index>0)
		OwnerDrawLbox_SelectItem(pflm->lbox,sffs.index);
	temport=SavePort(w);
	Select(pflm->lbox);
	ObjectRect(pflm->lbox,&rcP);
	InvalRect(&rcP);
	RestorePort(temport);
}

/*******************************************************************************

  Function : FeatListDlgCleanupProc()
  
  Purpose : clean data structure of Features DlgBox
  
  Parameters : data; something to delete
  
  Return value : none 

*******************************************************************************/
static void FeatListDlgCleanupProc (GraphiC g, VoidPtr data)
{
FLMDataPtr		pflm=(FLMDataPtr)data;

	if (pflm){
		if (pflm->SeqFeatClass) Free(pflm->SeqFeatClass);
		Free(pflm);
	}
}


/*******************************************************************************

  Function : InitFeaturesListDlg()
  
  Purpose : initialize Features DlgBox
  
  Parameters : vmp; main data 
  				hWndMain; handle of the viewer main window
				Search; TRUE is keyword search
				szKeyword; used only for keyword search
  
  Return value : none 

*******************************************************************************/
static void	InitFeaturesListDlg(ViewerMainPtr vmp,WindoW hWinMain,
	Boolean Search,CharPtr szKeyword)
{	
WindoW			w,lbox; 
Int2			Margins,nFeatClass,nClass=SEQFEAT_MAX+1;
Int1			k;
RecT			rcDlg,rc;
FLMDataPtr		pflm;
ScanFeatForSelect sffs;
GrouP g2;
PopuP pop;
PrompT txt;

	Margins=4*stdCharWidth;
	w=DocumentWindow(Margins,Margins ,
			(Int2)((screenRect.right-screenRect.left)/2-2*Margins), 
			(Int2)((screenRect.bottom-screenRect.top)/3-2*Margins), 
			"Features List", 
			FeatListDlgQuit,
			FeatListDlgProc);

	if (w==NULL){
		Message (MSG_ERROR, "Dialog creation failed.");
		return;
	}
	
	Select(w);
	ObjectRect(w,&rcDlg);
	
	if (!Search){
		g2=HiddenGroup(w,2,0,NULL);
		txt=StaticPrompt(g2,"Feature classes :",0,0,NULL,'l');

#ifdef WIN_MAC
		pop=PopupList(g2,TRUE,FeatClassChooseProc);
#endif

#ifndef WIN_MAC
		pop=PopupList(g2,FALSE,FeatClassChooseProc);
#endif
		AlignObjects(ALIGN_MIDDLE,(HANDLE) txt, (HANDLE) pop, NULL);

		ObjectRect(g2,&rc);
		rc.bottom+=3;
	}
	else rc.bottom=10;

	lbox=OwnerDrawLbox_Create(w,10,rc.bottom,
		(Int2)(rcDlg.right-rcDlg.left-20),
		(Int2)(rcDlg.bottom-rcDlg.top-10-rc.bottom), 
		FeatListlBoxProc,
		vmp->vdp->udv_graph.udv_font.hFnt,
		vmp->vdp->udv_graph.udv_font.cxChar,
		vmp->vdp->udv_graph.udv_font.cyChar);
		
	if (lbox==NULL) {
		Remove(w);
		return;
	}

	pflm=(FLMDataPtr)MemNew(sizeof(FLMData));
	if (pflm==NULL){
		Remove(w);
		return;
	}
	MemSet(pflm,0,sizeof(FLMData));
	pflm->hWndMain=hWinMain;
	pflm->lbox=lbox;
	if (!Search){
		pflm->gFeatClass=g2;
		pflm->pop=pop;
	}
	else{
		pflm->gFeatClass=NULL;
		pflm->pop=NULL;
	}
	/*add features to the listbox*/
	WatchCursor();
	sffs.eID=vmp->vdp->Item_select.eIDsel;		
	sffs.iID=vmp->vdp->Item_select.iIDsel;		
	sffs.index=0;
	sffs.compteur=0;
	sffs.lbox=lbox;
	sffs.SearchFeat=Search;
	StringCpy(sffs.szKeyword,szKeyword);
	MemSet(sffs.SeqFeatListAvail,0,sizeof(sffs.SeqFeatListAvail));

	SeqMgrExploreFeatures (vmp->vdp->bsp_i.bsp, (Pointer) &sffs, 
				FeaturesListBoxInit, NULL, NULL, NULL);

	
	/*no features found*/
	if (sffs.compteur==0){
		RealizeWindow(w);
		goto error;
	}
	/*Select the currently selected feature of UNDViewer*/
	if (sffs.index>0)
		OwnerDrawLbox_SelectItem(pflm->lbox,sffs.index);

	/*Fill the popup list with Feature classes*/
	if (!Search){
		nFeatClass=0;
		for(k=0;k<nClass;k++){
			if (sffs.SeqFeatListAvail[k]) nFeatClass++;
		}
		PopupItem(pop,szSeqFeatClassName[0]);/*ALL*/
		if (nFeatClass>0){
			nFeatClass++;
			pflm->SeqFeatClass=(Int1Ptr)MemNew(nFeatClass*sizeof(Int1));
			if (pflm->SeqFeatClass){
				Int1 j=1;
				pflm->SeqFeatClass[0]=0;
				for(k=1;k<nClass;k++){
					if (sffs.SeqFeatListAvail[k]) {
						pflm->SeqFeatClass[j++]=k;
						PopupItem(pop,szSeqFeatClassName[k]);
					}
				}
			}
		}
		SetValue(pop,1);
	}
	/*show the dialog box*/
	ArrowCursor();
	SetObjectExtra (w, (Pointer) pflm, FeatListDlgCleanupProc);
	RealizeWindow(w);
	Show(w);
	Update();
	vmp->hFeatDlg=w;
	/*if true means this function has been called for a string search*/
	if (Search) SetStatus(vmp->MainMenu.ShowFeatureList,TRUE);
	return;

error:
	ArrowCursor();
	/*delete unused data*/
	if (pflm) Free(pflm);
	/*delete the dlg box*/
	Remove(w);
	/*update menu*/
	SetStatus(vmp->MainMenu.ShowFeatureList,FALSE);
	MsgAlert(KEY_OK, SEV_INFO, "UnD-Viewer","Nothing found.");
}

/*******************************************************************************

  Function : ShowFeaturesListDlg()
  
  Purpose : show/hide Features DlgBox
  
  Parameters : i; menu item
  
  Return value : none 

*******************************************************************************/
static void ShowFeaturesListDlg(IteM i)
{
ViewerMainPtr 	vmp;
WindoW hWinMain;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;

	if (GetStatus(i)==FALSE){
		Remove(vmp->hFeatDlg);
		vmp->hFeatDlg=NULL;
		return;
	}
	InitFeaturesListDlg(vmp, hWinMain, FALSE, NULL);
}

/*******************************************************************************

  Function : UDV_FeatKeySearchCancelProc()
  
  Purpose : manage cancel button of the Feature keyword search dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_FeatKeySearchCancelProc(ButtoN g)
{
WindoW				w;
FKSDataPtr			fksd;
ViewerMainPtr 		vmp;

	w=(WindoW)ParentWindow(g);

	if (!w) return;
	fksd = (FKSDataPtr) GetObjectExtra (w);
	
	if (fksd){
		vmp = (ViewerMainPtr) GetObjectExtra (fksd->hWinMain);
	    UDV_set_PullMenus(&vmp->MainMenu,TRUE);
	}
	
	Remove(w);
}

/*******************************************************************************

  Function : UDV_FeatKeySearchAcceptProc()
  
  Purpose : manage search button of the Feature keyword search dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_FeatKeySearchAcceptProc(ButtoN g)
{
Char 				szFName[PATH_MAX]={""};
WindoW				w;
FKSDataPtr			fksd;
ViewerMainPtr 		vmp;

	w=(WindoW)ParentWindow(g);

	if (!w) return;
	fksd = (FKSDataPtr) GetObjectExtra (w);
	
	if (fksd){
		vmp = (ViewerMainPtr) GetObjectExtra (fksd->hWinMain);
		/*retrieve the keyword*/
		GetTitle(fksd->keyField, szFName, sizeof(szFName)-1);
		/*find it; if yes: Features List will be created*/
		StringUpper(szFName);
		InitFeaturesListDlg(vmp, fksd->hWinMain, TRUE, szFName);
	    UDV_set_PullMenus(&vmp->MainMenu,TRUE);
	}
	
	Remove(w);
}

/*******************************************************************************

  Function : UDV_FeatKeySearchTextProc()
  
  Purpose : manage keyword edit control of the  Feature keyword search  
  			dialog box 
  
  Parameters : t; edit control
  
  Return value : none 

*******************************************************************************/
static void UDV_FeatKeySearchTextProc(TexT t)
{
Char 				szFName[PATH_MAX]={""};
WindoW				w;
FKSDataPtr			fksd;

	w=(WindoW)ParentWindow(t);

	if (!w) return;
	fksd = (FKSDataPtr) GetObjectExtra (w);

	if (fksd==NULL) return;
	
	GetTitle(t, szFName, sizeof(szFName)-1);

	if (StringLen(szFName) == 0)
		Disable(fksd->bOk);
	else Enable(fksd->bOk);

	return;
}

/*******************************************************************************

  Function : UDV_SearchFeatForKey()
  
  Purpose : callback of the Options|Search for features command 
  
  Parameters : i; command item
  
  Return value : none 

*******************************************************************************/
static void UDV_SearchFeatForKey(IteM i)
{
ViewerMainPtr 		vmp;
WindoW				hWinMain,w;
GrouP				g1,g2,g3;
ButtoN				b1,b2;
PrompT				t1;
TexT				t2;    
FKSDataPtr			fksd;
MsgAnswer			nRet;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;
	
	/*if hFeatDlg opens, ask to close it*/
	if (vmp->hFeatDlg){
		nRet=MsgAlert(KEY_OKC, SEV_INFO, "UnD-Viewer",
			"Features List is already opened.\n Would you like to close it?");
		if (nRet==ANS_OK){
			Remove(vmp->hFeatDlg);
			vmp->hFeatDlg=NULL;
			SetStatus(vmp->MainMenu.ShowFeatureList,FALSE);
		}
		else return;
	}

	fksd=(FKSDataPtr)MemNew(sizeof(FKSData));
	if (!fksd) return;
	MemSet(fksd,0,sizeof(FKSData));		

    w = FixedWindow(-30, -20,  -10,  -10, 
				"UnD-Viewer - Find Features",  NULL);
	g1=HiddenGroup(w,0,2,NULL);

	SetGroupSpacing(g1,10,10);
	g2=HiddenGroup(g1,2,0,NULL);
	t1=StaticPrompt(g2,"Keyword :",0,0,NULL,'l');
    t2=DialogText(g2,"",20, UDV_FeatKeySearchTextProc);

	AlignObjects(ALIGN_MIDDLE,(HANDLE) t1, (HANDLE) t2, NULL);

	g3=HiddenGroup(g1,2,0,NULL);

    b1 = DefaultButton(g3, "Search",  UDV_FeatKeySearchAcceptProc);
    b2 = PushButton(g3, "Cancel",  UDV_FeatKeySearchCancelProc);
   
	SetObjectExtra (w, (Pointer) fksd, StdCleanupExtraProc);
	fksd->hWinMain=hWinMain;
	fksd->keyField=t2;
	fksd->bOk=b1;
	RealizeWindow(w);
    Show(w);
    Select(w);
	Select(t2);

	/*disable all menus*/
    UDV_set_PullMenus(&vmp->MainMenu,FALSE);
}

