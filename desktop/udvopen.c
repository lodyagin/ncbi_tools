/*   udvopen.c
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
* File Name:  udvopen.c
*
* Author:  Patrick Durand
*
* Version Creation Date:   5/3/99
*
* $Revision: 6.2 $
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

	typedef struct dlgfileopendata {/*use to manage the FileOpen dialog box*/
		WindoW 		parent;			/*main window of the application*/
		TexT		FNameEditCtrl;	/*handle of the file text control*/
		ButtoN		Ok;				/*handle of the Ok button*/
		GrouP		ReadMode;		/*handle of the file type control*/
	} DlgFileOpenData, PNTR DlgFileOpenDataPtr;


/*local struct used only by the Download sequence dialog box*/
	typedef struct udvnetopen {
		WindoW 			hWndMain;	/*main window*/
		TexT			Entry;		/*Entrez entry*/
		PopuP			AccessType;	/*database type*/
		ButtoN  		ok;			/*ok button*/
		UDVMainMenuPtr 	mmp;		/*main menu*/
	} UDVNetOpen, PNTR UDVNetOpenPtr;

/*local struct used only by the Choose sequence dialog box*/
	typedef struct udvchooseseq {
		WindoW 	hWndMain;		/*main window*/
		LisT	bsp_list;		/*listbox of bsps*/
	} UDVChooseSeq, PNTR UDVChooseSeqPtr;

	static Char szAppName[]="UnD-Viewer";

	static Uint1 DataBaseID[]={
				SEQID_GENBANK,
				SEQID_EMBL,
				SEQID_PIR,
				SEQID_SWISSPROT,
				SEQID_GI,
				SEQID_DDBJ,
				SEQID_PRF,
				SEQID_PDB};


/*******************************************************************************

  Function : UDV_Init_vdp_struct()
  
  Purpose : free memory of a bsp
  
  Parameters : vdp; contains data to delete
  				EraseParaG; TRUE means delete vdp->ParaG
				EraseMainTitle; TRUE means restore default viewer window name
				EraseInfoPanel; TRUE means delete content of InfoPanel
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void  UDV_Init_vdp_struct(PaneL p,ViewerDialogDataPtr vdp, 
	Boolean EraseParaG,Boolean EraseMainTitle,Boolean EraseInfoPanel)
{
BaR 	vsb;
WindoW 	w;

	if (!vdp) return;
	
	if (EraseParaG) UDV_FreeListParaG(&vdp->ParaG);

	MemSet(&vdp->bsp_i,0,sizeof(BspInfo));

	/*invalid Item selection identifier*/
	vdp->Item_select.eIDsel=(Uint2)-1;
	vdp->Item_select.iIDsel=(Uint2)-1;
	vdp->Item_select.iTypeSel=(Uint2)-1;
	vdp->Old_Item_select=vdp->Item_select;

	/*reinit scrollbar*/
	vsb = GetSlateVScrollBar ((SlatE)p);
	CorrectBarMax(vsb,0);
	CorrectBarValue(vsb,0);
	CorrectBarPage(vsb,0,0);
	
	/*reset winmain title*/
	if (EraseMainTitle && p){
		w=ParentWindow(p);
		if (w) SetTitle(w,szAppName);
	}
	
	/*reset infopanel title*/
	if(EraseInfoPanel && vdp->InfoPanel){
		CharPtr szFeatName=(CharPtr) GetObjectExtra (vdp->InfoPanel);
		RecT rcI;
		WindoW temport;
		
		if (szFeatName) {
			MemSet(szFeatName,0,sizeof(szFeatName));
			temport = SavePort(vdp->InfoPanel);
			Select(vdp->InfoPanel);
			ObjectRect(vdp->InfoPanel,&rcI);
			rcI.left=0;rcI.top=0;
			InvalRect(&rcI);
			Update();
			RestorePort(temport);
		}
	}
}

/*******************************************************************************

  Function : SearchBioseq()
  
  Purpose : store bsps of a SeqEntry 
  
  Parameters : see Toolbox doc !
  
  Return value : none 

*******************************************************************************/
static Boolean LIBCALLBACK SearchBioseq (BioseqPtr bsp, 
			SeqMgrBioseqContextPtr context)
{
ValNodePtr PNTR	BspTable;
	
	BspTable=(ValNodePtr PNTR)context->userdata;

	if (bsp) {
		ValNodeAddInt(BspTable,1,
				UDV_EncodeIdxFeat(context->entityID,context->itemID));
	}

	return(TRUE);
}


/*******************************************************************************

  Function : UDV_analyze_SEP_for_open()
  
  Purpose : analyse a SEP coming form a local SeqEntry file or from Entrez,
  			and find the BSPs. 
  
  Parameters : the_set; SeqEntry to analyze
  				vmp; general data software
				w; window which has called this function
				mon; monitor
  
  Return value : If error, return FALSE.  

*******************************************************************************/
NLM_EXTERN Boolean  UDV_analyze_SEP_for_open(FILE *fp,SeqEntryPtr the_set,
	ViewerMainPtr vmp,WindoW w)
{	
BioseqPtr 	bsp;
Boolean 	bRet=TRUE,bReadOk=FALSE;
Pointer     dataptr;
Uint2       datatype,entityID,eID,iID;
MonitorPtr  mon;
PaneL		p1,p2;
	
	if (vmp->hFeatDlg){/*Features List Dlg Box*/
		SetStatus(vmp->MainMenu.ShowFeatureList,FALSE);
		Remove(vmp->hFeatDlg);
		vmp->hFeatDlg=NULL;
	}
	/*free data, remove viewer*/
	if (vmp->vdp){
		p1=vmp->vdp->UnDViewer;
		p2=vmp->vdp->InfoPanel;
		vmp->vdp=UDV_FreeVDPstruct(&vmp->vdp);

		UDV_set_MainMenus(&vmp->MainMenu,FALSE);
		
		Remove(p1);
		if (p2) Remove(p2);
	}
	if (vmp->BspTable){
		ValNodeFree(vmp->BspTable);
		vmp->BspTable=NULL;
		vmp->BspChoice=0;
	}
	if (vmp->dataptr && vmp->datatype){
		ObjMgrFree(vmp->datatype,vmp->dataptr);
		vmp->datatype=0;
		vmp->dataptr=NULL;
	}

	mon = MonitorStrNewEx (szAppName, 30, FALSE);
	WatchCursor();
	
	if (fp!=NULL){/*open a file*/
		/*Read the file*/
		MonitorStrValue (mon, "Reading ASN.1 file");
		Update ();
		dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, 
			FALSE, FALSE, FALSE, FALSE);
		/*register for ObjMgr*/
		entityID = ObjMgrRegister (datatype, dataptr);
		if(datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        	datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET){
			bReadOk=TRUE;
		}
	}
	else if (the_set!=NULL){/*retrieve a Network Entry*/
		/*a complete verification of sep is done elsewhere; here
		sep can be only Bioseq or BioseqSet*/
		if (IS_Bioseq (the_set)) datatype = OBJ_BIOSEQ;
		else datatype = OBJ_BIOSEQSET;

		dataptr=(Pointer)the_set->data.ptrvalue;
		/*register for ObjMgr*/
		entityID = ObjMgrRegister (datatype, dataptr);
		bReadOk=TRUE;
		
	}
	
	if (bReadOk){
		/*compute the global Feature Index*/
		MonitorStrValue (mon, "Building Feature index...");
		Update ();
		UDV_CreateOneFeatureIndex(entityID,NULL);						
		/*scan the BSPs*/
		SeqMgrExploreBioseqs (entityID, NULL, (Pointer) &vmp->BspTable, 
			SearchBioseq, TRUE, TRUE, TRUE);
	}

	if (vmp->BspTable==NULL){
		MonitorFree (mon);
		ArrowCursor();
		Message (MSG_ERROR, "No sequence found in your entry.");
		bRet=FALSE;
		goto error;
	}

	vmp->entityID=entityID;
	vmp->datatype=datatype;
	vmp->dataptr=dataptr;
	
	/*choose the first sequence by default*/
	UDV_DecodeIdxFeat ((Uint4)vmp->BspTable->data.intvalue,&eID,&iID);
	bsp=GetBioseqGivenIDs (eID, iID, OBJ_BIOSEQ);
	if (bsp){
		vmp->BspChoice=1;
		/*temport=SavePort(w);
		Select(vmp->vdp->UnDViewer);
		ObjectRect(vmp->vdp->UnDViewer,&rcP);
		UDV_init_bsp_forViewer(w,bsp,eID,iID,OBJ_BIOSEQ,vmp->vdp);

		if (vmp->vdp->ParaG){
			UDV_set_MainMenus(&vmp->MainMenu,TRUE);
			rcP.top=0;rcP.left=0;
			InvalRect(&rcP);
			Update();
		}
		else{
			if (vmp->BspTable) ValNodeFree(vmp->BspTable);
			vmp->BspTable=NULL;
			UDV_set_MainMenus(&vmp->MainMenu,FALSE);
			UDV_Init_vdp_struct(vmp->vdp,FALSE,TRUE,TRUE);
		}
		RestorePort(temport);*/
		/*Ask ObjMgr for the viewer*/
		ArrowCursor();
		MonitorFree (mon);
		vmp->Show_logo=FALSE;	
		Hide(vmp->Logo_Panel);
		/*ask ObjMgr to load the viewer panel*/
		GatherProcLaunch(OMPROC_VIEW, FALSE, eID, iID, OBJ_BIOSEQ, 
			OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0);
	}
	else {
		MonitorFree (mon);
		ArrowCursor();
		Message (MSG_ERROR, "No sequence found in your file.");
		bRet=FALSE;
	}

error:
	if (!bRet){/*show the logo in case of error*/
		RecT rcL;
		UDV_Resize_Logo_Panel (vmp->hWndMain,&rcL);
		SetPosition (vmp->Logo_Panel, &rcL );
		AdjustPrnt (vmp->Logo_Panel, &rcL, FALSE);
		vmp->Show_logo=TRUE;	
		Show(vmp->Logo_Panel);
		UDV_set_MainMenus(&vmp->MainMenu,FALSE);
	}
	return (bRet);
}

/*******************************************************************************

  Function : AccessionToGi_ID1()
  
  Purpose : given s string, retrieve the UID 
  
  Parameters : string, either a name or an accession number	
  				type; entry choosen by the nice user in the dataBase type
					popup list
  
  Return value : none 

*******************************************************************************/
static Int4 AccessionToGi_ID1 (CharPtr string,Int2 type)
{
SeqIdPtr 	sip;
Int4 		uid=0;
Uint1 		Choice;

	sip=(SeqIdPtr)ValNodeNew(NULL);
	if (!sip) return(0);

	Choice=DataBaseID[type-1];

	switch(Choice){
		case SEQID_GENBANK:
		case SEQID_EMBL:
		case SEQID_PIR:
		case SEQID_DDBJ:
		case SEQID_PRF:
		case SEQID_SWISSPROT:{
			TextSeqIdPtr tsip;
			
			tsip=TextSeqIdNew();
			if (!tsip) break;
			sip->data.ptrvalue=(VoidPtr)tsip;
			sip->choice = Choice;
			/*try to retrieve an uid; to avoid the user to know whether string
			is an accession or a name, this function tests both*/
			tsip->name=StringSave(string);
			tsip->accession=StringSave(string);
			uid=GetGIForSeqId(sip);

			/*MemFree(tsip->name);
			tsip->name=NULL;
			if (uid==0){
				tsip->accession=StringSave(string);
				uid=GetGIForSeqId(sip);
			}*/
			
			TextSeqIdFree(tsip);
			break;
		}
		case SEQID_GI:
			if (!StrToLong(string,&uid)) uid=0;
			break;
		case SEQID_PDB:{
			PDBSeqIdPtr psip;
			
			psip=PDBSeqIdNew();
			if (!psip) break;
			sip->data.ptrvalue=(VoidPtr)psip;
			/*type is one of the SEQID_... reported in objloc.h*/
			sip->choice = Choice;
			/*try to retrieve an uid; to avoid the user to know whether string
			is an accession or a name, this function tests both*/
			psip->mol=StringSave(string);
			psip->chain=0;        /* 0 = no chain set.  default = 32 */
        	psip->rel=NULL;
			uid=GetGIForSeqId(sip);
			PDBSeqIdFree(psip);
			break;
		}
	}

	ValNodeFree(sip);
	return(uid);
}


/*******************************************************************************

  Function : UDV_NetOpen_okProc()
  
  Purpose : manage ok button of the Network open dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_NetOpen_okProc(ButtoN g)
{
WindoW			hOpenDlg;
UDVNetOpenPtr 	unop;
Char 			szAccess[50]={""};
Int2			AccessType;
Int4			uid=0;
ViewerMainPtr 	vmp;
MonitorPtr  	mon;
SeqEntryPtr		sep=NULL;
Boolean			bReadok=FALSE;
UdvGlobalsPtr   ugp=NULL;

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;
	
	unop = (UDVNetOpenPtr) GetObjectExtra (hOpenDlg);	

	if (unop==NULL) return;

	vmp = (ViewerMainPtr) GetObjectExtra (unop->hWndMain);
	if (vmp==NULL) return;

	Hide(hOpenDlg);
	
	/*retrieve and analyze Accession or Gi*/
	GetTitle(unop->Entry, szAccess, sizeof(szAccess));

	if (StringHasNoText (szAccess)) {
		Message (MSG_OK, "Please enter an accession number or an entry name");
		Show (hOpenDlg);
		Select (hOpenDlg);
		Select (unop->Entry);
		return;
	}

	AccessType=GetValue(unop->AccessType);
	WatchCursor();
	uid=AccessionToGi_ID1(szAccess,AccessType);	
	/*Connect ID1 with a valid UID*/
	if (uid==0){
		ArrowCursor();
		Message (MSG_OK, "Unable to find your record in the database.");
		Show (hOpenDlg);
		Select (hOpenDlg);
		Select (unop->Entry);
		return;
	}
	else{
		mon = MonitorStrNewEx (szAppName, 30, FALSE);
		MonitorStrValue (mon, "Retrieving your record...");
		Update ();
	
		ugp=(UdvGlobalsPtr)GetAppProperty("UdvGlobals");
		if (ugp && ugp->fetchSepProc) {
				sep=ugp->fetchSepProc(uid, 0);
		}
		
		/*sep=ID1SeqEntryGet(uid, 0);*/
		
		ArrowCursor();
		MonitorFree (mon);
		if (sep){
			if (IS_Bioseq (sep) || IS_Bioseq_set (sep))	bReadok=TRUE;
			else bReadok=FALSE;
		}
		else bReadok=FALSE;
		
		if (bReadok){
			UDV_analyze_SEP_for_open(NULL,sep,vmp,unop->hWndMain);
		}
		else{
			Message (MSG_OK, "Unable to retrieve your record in the database.");
		}
	}

	/*enable main menus and delete Download dlg box*/
	UDV_set_PullMenus(unop->mmp,TRUE);
	Remove(hOpenDlg);
}

/*******************************************************************************

  Function : UDV_NetOpen_cancelProc()
  
  Purpose : manage ok button of the Network open dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_NetOpen_cancelProc(ButtoN g)
{
WindoW			hOpenDlg;
UDVNetOpenPtr 	unop;

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;
	
	unop = (UDVNetOpenPtr) GetObjectExtra (hOpenDlg);

	if (unop==NULL) return;

	/*enable main menus and delete Download dlg box*/
	UDV_set_PullMenus(unop->mmp,TRUE);
	Remove(hOpenDlg);
}

/*******************************************************************************

  Function : UDV_NetOpen_access()
  
  Purpose : manage Accession edit control of the Download dialog box 
  
  Parameters : t; edit control
  
  Return value : none 

*******************************************************************************/
static void UDV_NetOpen_access(TexT t)
{
Char 			szAccess[50]={""};
WindoW			hOpenDlg;
UDVNetOpenPtr 	unop;

	hOpenDlg=(WindoW)ParentWindow(t);

	if (!hOpenDlg) return;
	
	unop = (UDVNetOpenPtr) GetObjectExtra (hOpenDlg);

	if (unop==NULL) return;
	
	GetTitle(t, szAccess, sizeof(szAccess));

	if (StringLen(szAccess) == 0)
		Disable(unop->ok);
	else Enable(unop->ok);

	return;
}

/*******************************************************************************

  Function : UDV_NetOpen()
  
  Purpose : retrieve a file from Download dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void UDV_NetOpen (IteM i)
{
GrouP			c,g;
WindoW			w,hWinMain;
UDVNetOpenPtr 	unop;
ViewerMainPtr 	vmp;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);
	if (vmp==NULL) return;

	w = FixedWindow (-50, -33, -10, -10, "Download From NCBI", NULL);

	if (w==NULL) return;

	unop=(UDVNetOpenPtr)MemNew(sizeof(UDVNetOpen));
	if (unop==NULL) {
		Remove(w);
		return;
	}
	unop->hWndMain=hWinMain;
	
	SetGroupSpacing (w, 10, 10);

	/*accesion*/
	g = NormalGroup (w, 3, 0, "Entry", systemFont,NULL);

#ifdef WIN_MAC
	unop->AccessType = PopupList (g, TRUE, NULL);
#endif

#ifndef WIN_MAC
	unop->AccessType = PopupList (g, FALSE, NULL);
#endif

	PopupItem(unop->AccessType,"GENBANK");
	PopupItem(unop->AccessType,"EMBL");
	PopupItem(unop->AccessType,"PIR");
	PopupItem(unop->AccessType,"SWISSPROT");
	PopupItem(unop->AccessType,"GI number");
	PopupItem(unop->AccessType,"DDBJ");
	PopupItem(unop->AccessType,"PRF");
	PopupItem(unop->AccessType,"PDB");
	
	SetValue (unop->AccessType, 6);

	unop->Entry = DialogText (g, "", 10, UDV_NetOpen_access);

	/*retrieve - cancel*/
	c = HiddenGroup (w, 4, 0, NULL);
	SetGroupSpacing (c, 10, 2);
	unop->ok=DefaultButton (c, "Retrieve", UDV_NetOpen_okProc);
	Disable (unop->ok);
	PushButton (c, "Cancel", UDV_NetOpen_cancelProc);
	
	SetObjectExtra (w, (Pointer) unop, StdCleanupExtraProc);

	/*display*/
	AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
	RealizeWindow (w);
	Show (w);
	Select (w);
	Update ();
	
	/*disable all menus*/
	unop->mmp=&vmp->MainMenu;
	UDV_set_PullMenus(unop->mmp,FALSE);
}

/*******************************************************************************

  Function : UDV_FOpenTextProc()
  
  Purpose : manage File name edit control of the FileOpen dialog box 
  
  Parameters : t; edit control
  
  Return value : none 

*******************************************************************************/
static void UDV_FOpenTextProc(TexT t)
{
Char 				szFName[PATH_MAX]={""};
WindoW				hOpenDlg;
DlgFileOpenDataPtr  dfodp;

	hOpenDlg=(WindoW)ParentWindow(t);

	if (!hOpenDlg) return;
	
	dfodp = (DlgFileOpenDataPtr) GetObjectExtra (hOpenDlg);

	if (dfodp==NULL) return;
	
	GetTitle(t, szFName, sizeof(szFName)-1);

	if (StringLen(szFName) == 0)
		Disable(dfodp->Ok);
	else Enable(dfodp->Ok);

	return;
}

/*******************************************************************************

  Function : UDV_FOpenBrowseProc()
  
  Purpose : manage browse button of the FileOpen dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_FOpenBrowseProc(ButtoN g)
{
DlgFileOpenDataPtr  dfodp;
WindoW				hOpenDlg;
Char 				path[PATH_MAX]={""};

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;

	dfodp = (DlgFileOpenDataPtr) GetObjectExtra (hOpenDlg);

	if (dfodp==NULL) return;
	if (!dfodp->FNameEditCtrl) return;

	if (GetInputFileName (path, sizeof(path), NULL, NULL)){ 
		SetTitle(dfodp->FNameEditCtrl, path);
		UDV_FOpenTextProc(dfodp->FNameEditCtrl);
	}

	return;   
}

/*******************************************************************************

  Function : UDV_FOpenCancelProc()
  
  Purpose : manage cancel button of the FileOpen dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_FOpenCancelProc(ButtoN g)
{
WindoW				hOpenDlg;
DlgFileOpenDataPtr  dfodp;
ViewerMainPtr 		vmp;

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;
	dfodp = (DlgFileOpenDataPtr) GetObjectExtra (hOpenDlg);
	
	if (dfodp){
		vmp = (ViewerMainPtr) GetObjectExtra (dfodp->parent);
		Enable(vmp->MainMenu.FileOpen);
	}
	
	Remove(hOpenDlg);
    UDV_set_PullMenus(&vmp->MainMenu,TRUE);
}


/*******************************************************************************

  Function : UDV_FOpenAcceptProc()
  
  Purpose : manage ok button of the FileOpen dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_FOpenAcceptProc(ButtoN g)
{
WindoW				hOpenDlg,w;
DlgFileOpenDataPtr  dfodp;
ViewerMainPtr 		vmp;
Char 				szFName[PATH_MAX]={""};
Boolean				isBinary;
FILE				*fp;

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;
	dfodp = (DlgFileOpenDataPtr) GetObjectExtra (hOpenDlg);
	
	if (dfodp){
		w=dfodp->parent;
		vmp = (ViewerMainPtr) GetObjectExtra (w);
		if (vmp!=NULL){
			Enable(vmp->MainMenu.FileOpen);
			/*file name*/
			GetTitle(dfodp->FNameEditCtrl, szFName, sizeof(szFName));

			/*file reading mode*/
			if (GetValue(dfodp->ReadMode)==2) isBinary=TRUE;
			else isBinary=FALSE;
	
			Remove(hOpenDlg);

			/* open the i/o files in the correct mode */
			if ((fp = FileOpen (szFName, isBinary?"rb":"r")) == NULL){
				Message (MSG_ERROR, "Unable to open your file.");
				goto error;
			}
			
			UDV_analyze_SEP_for_open(fp,NULL,vmp,w);
			FileClose(fp);
		}
	}
	
error:	
    UDV_set_PullMenus(&vmp->MainMenu,TRUE);
}

/*******************************************************************************

  Function : UDV_FileOpen()
  
  Purpose : callback of the File|Open command 
  
  Parameters : i; command item
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void UDV_FileOpen(IteM i)
{
DlgFileOpenDataPtr  dfodp;
ViewerMainPtr 		vmp;
WindoW				hWinMain,hOpenDlg;
GrouP				g,g4;
ButtoN				b,b1,b2;
TexT				t1;    

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;
		
	dfodp=(DlgFileOpenDataPtr)MemNew(sizeof(DlgFileOpenData));
	if (!dfodp) return;
	MemSet(dfodp,0,sizeof(DlgFileOpenData));

    hOpenDlg = FixedWindow(-30, -20,  -10,  -10, 
				"UnD-Viewer - Open a local file",  NULL);
    g = NormalGroup(hOpenDlg, 2, 1, "File name:",  systemFont, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 10, 20);  
    t1 = DialogText(g,"",25, UDV_FOpenTextProc);
    b = PushButton(g, " browse...", UDV_FOpenBrowseProc);
   
    g = HiddenGroup(hOpenDlg, 3, 1, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 30, 30);
 
    g4 = NormalGroup(g, 2, 1, "File type", systemFont,  NULL);
    SetGroupMargins(g4, 10, 10);
    RadioButton(g4, "Ascii");
    RadioButton(g4, "Binary");
    SetValue(g4, 1);

    b1 = DefaultButton(g, "OK",  UDV_FOpenAcceptProc);
    b2 = PushButton(g, "Cancel",  UDV_FOpenCancelProc);
   
    Disable(b1);

	dfodp->parent=hWinMain;	
	dfodp->FNameEditCtrl=t1;
	dfodp->Ok=b1;
	dfodp->ReadMode=g4;

	SetObjectExtra (hOpenDlg, (Pointer) dfodp, StdCleanupExtraProc);
	
    Select(hOpenDlg);
    Show(hOpenDlg);

	/*allow only one FOpen dlg at a time*/
    UDV_set_PullMenus(&vmp->MainMenu,FALSE);
}

/*******************************************************************************

  Function : UDV_FileClose()
  
  Purpose : callback of the File|Close command 
  
  Parameters : i; command item
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void UDV_FileClose(IteM i)
{
ViewerMainPtr 		vmp;
WindoW				hWinMain;
RecT				rcL;
PaneL				p1,p2;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;
	
	if (vmp->vdp==NULL) return;
		
	if (vmp->hFeatDlg){/*Features List Dlg Box*/
		SetStatus(vmp->MainMenu.ShowFeatureList,FALSE);
		Remove(vmp->hFeatDlg);
		vmp->hFeatDlg=NULL;
	}

	/*free data, remove viewer*/
	p1=vmp->vdp->UnDViewer;
	p2=vmp->vdp->InfoPanel;

	vmp->vdp=UDV_FreeVDPstruct(&vmp->vdp);

	if (vmp->BspTable){
		ValNodeFree(vmp->BspTable);
		vmp->BspTable=NULL;
		vmp->BspChoice=0;
	}
	if (vmp->dataptr && vmp->datatype){
		ObjMgrFree(vmp->datatype,vmp->dataptr);
		vmp->datatype=0;
		vmp->dataptr=NULL;
		vmp->entityID=0;
	}

	UDV_set_MainMenus(&vmp->MainMenu,FALSE);
	
	Remove(p1);
	if (p2) Remove(p2);
	
	UDV_Resize_Logo_Panel (hWinMain,&rcL);
	SetPosition (vmp->Logo_Panel, &rcL );
	AdjustPrnt (vmp->Logo_Panel, &rcL, FALSE);
	vmp->Show_logo=TRUE;	
	Show(vmp->Logo_Panel);

}

/*******************************************************************************

  Function : UDV_CSeq_okProc()
  
  Purpose : manage ok button of the Choose sequence dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_CSeq_okProc(ButtoN g)
{
ViewerDialogDataPtr vdp;
UDVChooseSeqPtr 	ucsp;
ViewerMainPtr 		vmp;
Int2 				value;
WindoW				hOpenDlg;
ValNodePtr			vnp;
Uint2				nCompt=1,eID,iID;
BioseqPtr			bsp=NULL;
Char 				szBuf[255]={""};
PaneL				p1,p2;

	hOpenDlg=(WindoW)ParentWindow(g);
	if (!hOpenDlg) return;
	
	ucsp=(UDVChooseSeqPtr)GetObjectExtra (hOpenDlg);
	if (ucsp==NULL) goto error;
	
	vmp = (ViewerMainPtr) GetObjectExtra (ucsp->hWndMain);
	if (vmp==NULL) goto error;
	
	UDV_set_PullMenus(&vmp->MainMenu,TRUE);
	
	vdp = vmp->vdp;
	if (vdp == NULL) return;
	
	value=GetValue(ucsp->bsp_list);/*one-base value*/
	if (vmp->BspChoice==value) goto error;

	/*look for a new bsp*/
	for (vnp=vmp->BspTable ; vnp!=NULL ;vnp=vnp->next){
		if (nCompt==value){
			UDV_DecodeIdxFeat ((Uint4)vnp->data.ptrvalue,&eID,&iID);
			bsp=GetBioseqGivenIDs (eID, iID, OBJ_BIOSEQ);
			vmp->BspChoice=value;
			break;
		}
		nCompt++;
	}
	
	if (bsp){
		if (vmp->hFeatDlg){/*Features List Dlg Box*/
			SetStatus(vmp->MainMenu.ShowFeatureList,FALSE);
			Remove(vmp->hFeatDlg);
			vmp->hFeatDlg=NULL;
		}
		/*free data, remove viewer*/
		p1=vmp->vdp->UnDViewer;
		p2=vmp->vdp->InfoPanel;
		vmp->vdp=UDV_FreeVDPstruct(&vmp->vdp);

		UDV_set_MainMenus(&vmp->MainMenu,FALSE);
		
		Remove(p1);
		if (p2) Remove(p2);

		GatherProcLaunch(OMPROC_VIEW, FALSE, eID, iID, OBJ_BIOSEQ, 
			OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0);
	
/*		UDV_Init_vdp_struct(vdp, TRUE,TRUE,TRUE);
		temport=SavePort(ucsp->hWndMain);
		Select(vdp->UnDViewer);
		ObjectRect(vdp->UnDViewer,&rcP);
		UDV_init_bsp_forViewer(ucsp->hWndMain,bsp,eID,iID,OBJ_BIOSEQ,vdp);	
		sprintf(szBuf,"%s - [%s - %d letters]",szAppName,
					vdp->bsp_i.bspAccNum,
					vdp->bsp_i.bspLength);
		SetTitle(ucsp->hWndMain,szBuf);
		rcP.left=0;rcP.top=0;
		InvalRect(&rcP);
		Update();
		RestorePort(temport);*/
	}
error:

	Remove(hOpenDlg);
}


/*******************************************************************************

  Function : UDV_CSeq_cancelProc()
  
  Purpose : manage cancel button of the Choose sequence dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void UDV_CSeq_cancelProc(ButtoN g)
{	
WindoW				hOpenDlg;
UDVChooseSeqPtr 	ucsp;
ViewerMainPtr 		vmp;

	hOpenDlg=(WindoW)ParentWindow(g);
	if (!hOpenDlg) return;
	
	ucsp=(UDVChooseSeqPtr)GetObjectExtra (hOpenDlg);
	if (ucsp==NULL) goto error;
	
	vmp = (ViewerMainPtr) GetObjectExtra (ucsp->hWndMain);
	if (vmp==NULL) goto error;
	
	UDV_set_PullMenus(&vmp->MainMenu,TRUE);

error:

	Remove(hOpenDlg);
}

/*******************************************************************************

  Function : UDV_CreateListBioseqDlg()
  
  Purpose : let the user choose between various bsp 
  
  Parameters : i; menu
  
  Return value : none 

*******************************************************************************/
NLM_EXTERN void UDV_CreateListBioseqDlg(IteM i)
{
GrouP 			g1,h,h1;
LisT 			lBox;
WindoW 			d,hWinMain;
ValNodePtr 		vnp;
BioseqPtr 		bsp;
Char 			szName[21]={""};
Char 			szBuf[50]={""};
ViewerMainPtr 	vmp;
UDVChooseSeqPtr ucsp;
Uint2			eID,iID;

	ucsp=(UDVChooseSeqPtr)MemNew(sizeof(UDVChooseSeq));
	if (ucsp==NULL) return;
	
	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	vmp = (ViewerMainPtr) GetObjectExtra (hWinMain);

	if (vmp==NULL) return;

	if (vmp->BspTable==NULL) return;
	
	d=DocumentWindow(-50, -33 ,-10, -10, szAppName, NULL,NULL);
	if (d!=NULL){
		UDV_set_PullMenus(&vmp->MainMenu,FALSE);
		/*create some controls*/
		h=HiddenGroup(d, 3, 0,  NULL);
		h1=HiddenGroup(h, 0, 2,  NULL);
		StaticPrompt(h1,"Choose a sequence :",0,0,systemFont,'l');

		lBox=SingleList(h1,20,6,NULL/*(LstActnProc)UDV_CSeq_LB_Proc*/);

		g1=HiddenGroup(h, 0, 2, NULL);
		PushButton(g1, "Ok",UDV_CSeq_okProc);
		PushButton(g1, "Cancel",UDV_CSeq_cancelProc);
		/*fill in the list box*/
		for (vnp=vmp->BspTable ; vnp!=NULL ;vnp=vnp->next){
			UDV_DecodeIdxFeat ((Uint4)vnp->data.ptrvalue,&eID,&iID);
			bsp=GetBioseqGivenIDs (eID, iID, OBJ_BIOSEQ);
			if (bsp){
				SeqIdWrite(bsp->id,szName,
						PRINTID_TEXTID_ACCESSION,20);
				if (szName){
					sprintf(szBuf,"%s (%s)",
						szName,
						(ISA_na(bsp->mol) ? "nucleic seq.":"protein"));
				}else{
					sprintf(szBuf,"Unknown name (%s)",
						(ISA_na(bsp->mol) ? "nucleic seq.":"protein"));
				}
				ListItem(lBox,szBuf);
			}
		}
		SetValue(lBox,vmp->BspChoice);

		ucsp->bsp_list=lBox;
		ucsp->hWndMain=hWinMain;
		SetObjectExtra (d, (Pointer) ucsp, StdCleanupExtraProc);
		RealizeWindow(d);
		Show(d);
	}
}
