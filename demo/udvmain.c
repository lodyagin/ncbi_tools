/*   udvmain.c
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
* File Name:  udvmain.c
*
* Author:  Patrick Durand
*
* Version Creation Date:   5/3/99
*
* $Revision: 6.5 $
*
* File Description: 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: udvmain.c,v $
* Revision 6.5  1999/07/30 20:10:15  durand
* updates for the new Entrez graphical viewer
*
* Revision 6.4  1999/06/07 15:40:30  durand
* add LOG line to keep track of the history
*
*
*
* ==========================================================================
*/

#include <accentr.h>
#include <accid1.h>
#include <udviewer.h>

/*local text*/
static Char szAppName[]="UnD-Viewer";

/*******************************************************************************

  Function : UDV_ID1Init()
  
  Purpose : Init the connexion to ID1 
  
  Parameters : -
  
  Return value : TRUE if success 

*******************************************************************************/
static Boolean  UDV_ID1Init(void)
{
MonitorPtr  mon;
Boolean     rsult;

	mon = MonitorStrNewEx (szAppName, 30, FALSE);
	MonitorStrValue (mon, "Connecting to ID1 service");
	Update ();
	rsult = ID1BioseqFetchEnable(szAppName,TRUE);
	MonitorFree (mon);
	Update ();
	return rsult;
}

/*******************************************************************************

  Function : Main()
  
  Purpose : Entry point of the software 
  
  Parameters : 
  
  Return value : 

*******************************************************************************/
Int2 Main(void)
{
ViewerMainPtr 	vmp;
WindoW			w; 
Int2			Margins;
Boolean 		isID1Ok;
RecT			rcL;
UdvGlobals      ug;
UDVLogoData		ldp;
	
	ErrSetMessageLevel(SEV_NONE);
	ErrSetOptFlags(EO_SHOW_CODES);
	ErrSetOptFlags(EO_XLATE_CODES);
	MemSet(&ug,0,sizeof(UdvGlobals));
		
	/*init some important stuffs*/
	UseLocalAsnloadDataAndErrMsg();

	isID1Ok=UDV_ID1Init();
	if (!isID1Ok){
		Message (MSG_ERROR, "Unable to connect Network.");
	}

	if (! AllObjLoad()){
		Message (MSG_ERROR, "AsnObjLoad() failed.");
		return(1);
	}

	if (! SubmitAsnLoad()){
		Message (MSG_ERROR, "SeqSubmitLoad() failed.");
		return(1);
	}

	if (!SeqCodeSetLoad ()){
		Message (MSG_ERROR, "SeqCodeSetLoad () failed.");
		return(1);
	}

	if (!GeneticCodeTableLoad()){
		Message (MSG_ERROR, "GeneticCodeTableLoad() failed.");
		return(1);
	}

	if (!FeatDefSetLoad()){
		Message (MSG_ERROR, "FeatDefSeqLoad() failed.");
		return(1);
	}
	
	/*init data blocks*/
	vmp=(ViewerMainPtr)MemNew(sizeof(ViewerMain));
	if (vmp){
		MemSet(vmp,0,sizeof(ViewerMain));
	}
	else{
		Message (MSG_ERROR, "Viewer creation failed.");
		return(1);
	}

	/*OBjMgr Callback Function declaration*/
	REGISTER_UDV_AUTONOMOUS;
	
	/*main window*/
	Margins=4*stdCharWidth;
	w=DocumentWindow(Margins,Margins ,
			(Int2)((screenRect.right-screenRect.left)-2*Margins), 
			(Int2)((screenRect.bottom-screenRect.top)-2*Margins), 
			szAppName, 
			UDV_WinMainProgQuit,
			UDV_WinMainResize);

	if (w==NULL){
		Message (MSG_ERROR, "Viewer creation failed.");
		return(1);
	}
	vmp->hWndMain=w;
	SetObjectExtra (w, (Pointer) vmp, (FreeProc)UDV_WinMainCleanupExtraProc);

	/*this is an autonomous viewer*/
	vmp->AutonomeViewer=TRUE;

	/*init menu*/
	UDV_SetupMenus(w,isID1Ok);

	UDV_set_MainMenus(&vmp->MainMenu,FALSE);
	/*init logo_panel*/
	LogoFontCreate(&ldp.f1,&ldp.f2,&ldp.f3);
	StringCpy(ldp.szTitle,"UnD-Viewer");
	StringCpy(ldp.szDesc,", a sequence viewer for GenBank");
	SetAppProperty("UDVLogoData",(Pointer)&ldp);	
	vmp->Logo_Panel=AutonomousPanel4(w,10,10,UDV_Logo_onDraw,
			NULL,NULL,0,NULL,NULL);
	UDV_Resize_Logo_Panel (w,&rcL);
	SetPosition (vmp->Logo_Panel, &rcL );
	AdjustPrnt (vmp->Logo_Panel, &rcL, FALSE);
	vmp->Show_logo=TRUE;	
	SetAppProperty("AutonomousUDVViewer",(Pointer)vmp);	
	
	/*ProcessUpdatesFirst(FALSE);*/
	
	RealizeWindow(w);
	Show(w);

	/*is there a file to open on the command line ?*/
	if (GetArgc()>1){
		/*is GetArgv()[1] a file name ?*/
		FILE *f;
		
		f=FileOpen(GetArgv()[1],"r");
		if (f){
			UDV_analyze_SEP_for_open(f,NULL,vmp,w);
			FileClose(f);
		}
	}

	ug.fetchSepProc=ID1SeqEntryGet;
	ug.vmp=vmp;
	SetAppProperty("UdvGlobals",(Pointer)&ug);
	
	ProcessEvents();

	ID1BioseqFetchDisable();

	RemoveAppProperty("AutonomousUDVViewer");	
	RemoveAppProperty("UDVLogoData");	
	RemoveAppProperty("UdvGlobals");

	return(0);
}

