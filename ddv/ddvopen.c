/*  $Id: ddvopen.c,v 1.31 2000/01/18 22:49:16 lewisg Exp $
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
* File Name:  ddvopen.c
*
* Author:  Patrick Durand
*
* Version Creation Date:   06/19/99
*
* $Revision: 1.31 $
*
* File Description: code to open a SeqAlign (file & Net) and code of the
* message callback for DeuxD-Viewer (DDV).
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ddvopen.c,v $
* Revision 1.31  2000/01/18 22:49:16  lewisg
* send OM_MSG_FLUSH to ddv/udv, tweak CPK coloration, misc bugs
*
* Revision 1.30  2000/01/12 21:52:16  durand
* add import function; update menus when DDV is loaded from Cn3D
*
* Revision 1.29  2000/01/11 15:05:23  durand
* remove network stuff
*
* Revision 1.28  2000/01/10 15:09:46  durand
* Use Entrez instead of ID1
*
* Revision 1.27  2000/01/05 21:11:14  durand
* update mouse click actions and DrawSequence function for a better use from ddv and cn3d
*
* Revision 1.26  1999/12/30 21:08:45  lewisg
* bioseq import dialog
*
* Revision 1.25  1999/12/30 13:45:50  beloslyu
* fix comment
*
* Revision 1.24  1999/12/29 22:55:03  lewisg
* get rid of seqalign id
*
* Revision 1.23  1999/12/23 19:22:06  durand
* modify default options for DDV when loaded from Cn3D
*
* Revision 1.22  1999/12/21 15:27:24  durand
* avoid to quit Cn3D when closing DDV
*
* Revision 1.21  1999/12/20 20:20:41  lewisg
* allow cn3d to do color and ddv to do case when both are running
*
* Revision 1.20  1999/12/07 21:40:13  durand
* add mouse modes menu and caret facility for the editor
*
* Revision 1.19  1999/12/06 16:19:19  durand
* add GoTo facility to DDV
*
* Revision 1.18  1999/12/03 23:17:23  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 1.17  1999/11/30 18:19:47  durand
* fix a problem with function declaration DDV_CloseData
*
* Revision 1.16  1999/11/29 15:26:26  durand
* designed a new GUI to fix problems under MacOS, Linux and SGI
*
* Revision 1.15  1999/11/17 22:43:58  durand
* speed up the selection manager for large SeqAlign
*
* Revision 1.14  1999/11/09 17:09:00  durand
* transfer some functions from ddvgraph to ddvcreate, so that ddvcreate remains Vibrant free and can be compiled with BLAST
*
* Revision 1.13  1999/11/04 22:11:38  durand
* add the Desktop to DDV. Add a better set of cleanup functions when closing DDV. Before creating color tables, try to get them from the SeqAlign
*
* Revision 1.12  1999/11/03 21:29:47  durand
* add CTRL and SHFT keys for mouse selection. redesign the loader functions of DDV to properly register the message callbacks
*
* Revision 1.11  1999/10/29 19:04:21  durand
* move DDVUpdateMSG in objmgr.h
*
* Revision 1.10  1999/10/29 14:15:39  durand
* add simple mouse selection functions
*
* Revision 1.9  1999/10/23 21:20:45  lewisg
* move g_hParent to ddvopen.c
*
* Revision 1.8  1999/10/23 14:54:34  durand
* resolve external symbol g_hParent
*
* Revision 1.7  1999/10/22 14:19:43  durand
* update the code for the startup functions of DDV drawing panel
*
* Revision 1.6  1999/10/20 18:37:31  lewisg
* add messagefunc for slave mode
*
* Revision 1.5  1999/10/20 13:17:18  durand
* add display for disc. SeqAlign tails
*
* Revision 1.4  1999/10/16 15:02:25  durand
* fixes due to toolkit build failed
*
* Revision 1.3  1999/10/15 21:57:36  durand
* add a UI for display options
*
* Revision 1.2  1999/10/12 15:01:29  lewisg
* resolve confict with internal/ddv
*
* Revision 1.1  1999/09/30 14:10:28  durand
* add ddv to toolkit
*
* Revision 1.14  1999/09/30 13:38:10  durand
* DDV_CreateDisplayFromIndex takes ParaG_Size as an argument
*
* Revision 1.13  1999/09/16 13:07:52  durand
* add File|Close and File|Open|Network commands
*
* Revision 1.12  1999/09/09 21:55:06  durand
* instantiate the Fle|Close command of DDV
*
* Revision 1.11  1999/09/02 17:36:07  durand
* reconcile ddvopen.c
*
* Revision 1.10  1999/08/19 17:15:27  wheelan
* added messages to indicate timing of index building vs color tables
*
* Revision 1.9  1999/08/04 18:02:12  wheelan
* changes to support new seqalign indexing
*
* Revision 1.8  1999/07/29 12:43:07  durand
* update DDV_GetAndCheckSeqAlign
*
* Revision 1.7  1999/07/20 17:18:23  durand
* update DDV_GetAndCheckSeqAlign for PopSet Viewer
*
* Revision 1.6  1999/07/20 14:58:01  durand
* use the Color Manager to display colored MSA
*
* Revision 1.5  1999/07/01 15:28:29  durand
* validate function loaders of DDV
*
* Revision 1.2  1999/06/28 22:07:19  durand
* add loader functions and clean the code with Lint and Purify
*
* Revision 1.1  1999/06/19 17:21:05  durand
* add Vibrant DDV code
*
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <vibrant.h>
#include <salfiles.h>
#include <pgppop.h>
#include <udviewer.h>
#include <ddvopen.h>
#include <ddvpanel.h>
#include <ddvmain.h>
#include <ddvgraph.h>
#include <ddvcreate.h>
#include <objmgr.h>
#include <gather.h>
#include <ddvcolor.h>
#ifdef _DDVEDIT
#include <samedit.h>
#endif /* _DDVEDIT */
#include <samutil.h>
#include <accutils.h>
#include <blast.h>


extern WindoW g_hParent = NULL;

/*local struct used only by the Download sequence dialog box*/
	typedef struct ddvnetopen {
		WindoW 			hWndMain;	/*main window*/
		TexT			Entry;		/*Entrez entry*/
		PopuP			AccessType;	/*database type*/
		ButtoN  		ok;			/*ok button*/
		UDVMainMenuPtr 	mmp;		/*main menu*/
	} DDVNetOpen, PNTR DDVNetOpenPtr;


/*******************************************************************************

  Function : DDV_GetSelectedRegions()
  
  Purpose :  get the selected region(s) of one bioseq.
  
  Parameters :	om_ssp;list of selected regions (usually this field points to
                  ObjMgr data)
				bsp_eID,bsp_iID; bioseq identifiers.
				    
  Return value : list of selected regions on the bioseq bsp_eID,bsp_iID

*******************************************************************************/
extern ValNodePtr DDV_GetSelectedRegions(SelStructPtr om_ssp, Uint2 bsp_eID,
	Uint2 bsp_iID)
{
SelStructPtr ssp;
SeqLocPtr    slp,slp2;
ValNodePtr   bsp_vnp=NULL,vnp;

	if (om_ssp==NULL || bsp_eID==0 || bsp_iID==0) return(NULL);

	ssp=om_ssp;
	
	while(ssp){
		if (ssp->entityID==bsp_eID && ssp->itemID==bsp_iID && 
			ssp->itemtype==OBJ_BIOSEQ && ssp->regiontype==OM_REGION_SEQLOC){
			slp=(SeqLocPtr)ssp->region;
			while(slp){
				slp2=slp->next;
				slp->next=NULL;
				if (!bsp_vnp){
					vnp=ValNodeAddPointer(NULL,0,(Pointer)slp);
					if (!vnp) return(NULL);
					bsp_vnp=vnp;
				}
				else{
					vnp=ValNodeAddPointer(&vnp,0,(Pointer)slp);
					if (!vnp){
						if (bsp_vnp) ValNodeFree(bsp_vnp);
						return(NULL);
					}
				}
				slp->next=slp2;
				slp=slp->next;
			}
		}
		ssp=ssp->next;
	}
	return(bsp_vnp);
}

/*******************************************************************************

  Function : DDV_IsLetterSelected()
  
  Purpose : check if a bsp_pos is selected (vnp_bsp is usually built with 
      DDV_GetSelectedRegions() function)
  
  Return value : TRUE if bsp_pos is selected

*******************************************************************************/
extern Boolean DDV_IsLetterSelected(ValNodePtr vnp_bsp, Int4 bsp_pos)
{
Boolean    bSelected=FALSE;
ValNodePtr vnp;
SeqLocPtr  slp;
Int4 bsp_start,bsp_stop;
	if (vnp_bsp==NULL || bsp_pos==(Int4)-1) return(FALSE);

	vnp=vnp_bsp;

	while(vnp){
		slp=(SeqLocPtr)vnp->data.ptrvalue;
		bsp_start=SeqLocStart(slp);
		bsp_stop=SeqLocStop(slp);
		if (bsp_pos>=bsp_start && bsp_pos<=bsp_stop){
			bSelected=TRUE;
			break;
		}
		vnp=vnp->next;
	}
	
	return(bSelected);
}

/*******************************************************************************

  Function : DDV_CheckBSPEntitiesInPgp()
  
  Purpose : check if entityID,itemID Object is in the current paragraph
  
  Return value : TRUE if success

*******************************************************************************/
/*static Boolean DDV_CheckBSPEntities(Uint4Ptr entitiesTbl,Int4 nRowBsp,
	Int4 nBsp,Uint2 eID,Uint2 iID)
{
Uint2     bsp_eID,bsp_iID;
Int4      i;
Boolean   bRet=FALSE;

	for (i=0;i<nBsp;i++){
		UDV_DecodeIdxFeat (entitiesTbl[i], &bsp_eID,&bsp_iID);
	}
	bsp=BioseqLockById(pgp->sip);
	if (bsp){
		bsp_eID=ObjMgrGetEntityIDForPointer((Pointer)bsp);
		bsp_iID = GetItemIDGivenPointer (bsp_eID, 
				OBJ_BIOSEQ, (Pointer) bsp);	
	
		if (bsp_eID==eID && bsp_iID==iID) bRet=TRUE;
		BioseqUnlockById(pgp->sip);
	}
	return(bRet);
}
*/
/*******************************************************************************

  Function : DDV_CheckBSPEntitiesInDisplay()
  
  Purpose : check if entityID,itemID Object is in the current DDV display
  
  Return value : TRUE if object found

*******************************************************************************/
static Boolean DDV_CheckBSPEntitiesInDisplay(MsaParaGPopListPtr mpplp,
		Uint2 entityID,Uint2 itemID)
{
Uint2     bsp_eID,bsp_iID;
Int4      i;
Boolean   bRet=FALSE;

	/*scan the ParaG list to find the Bioseq*/
	for (i=0;i<mpplp->nBsp;i++){
		UDV_DecodeIdxFeat (mpplp->entitiesTbl[i], &bsp_eID,&bsp_iID);
		if (bsp_eID==entityID && bsp_iID==itemID) {
			bRet=TRUE;
			break;
		}
	}
	return(bRet);
}

/*******************************************************************************

  Function : DDV_GetBspListGivenIDs()
  
  Purpose : analyse the ParaG struct and get a list of row(s) containing the
         BSP identified by eID and iID
  
  Return value : list of row number (one-based values)

*******************************************************************************/
static Int4Ptr DDV_GetBspListGivenIDs(MsaParaGPopListPtr mpplp,
		Uint2 entityID,Uint2 itemID,Int4Ptr nFound)
{
Int4       i;
Int4Ptr    pList=NULL;
Uint2     bsp_eID,bsp_iID;

	*nFound=0;
	
	/*scan the ParaG list to find all the ref. of the Bioseq
	(entityID,itemID) in the ParaG list*/
	for (i=0;i<mpplp->nBsp;i++){
		UDV_DecodeIdxFeat (mpplp->entitiesTbl[i], &bsp_eID,&bsp_iID);
		if (bsp_eID==entityID && bsp_iID==itemID) {
			if (!pList){
				pList=(Int4Ptr)MemNew(sizeof(Int4));
			}
			else{
				pList=(Int4Ptr)MemExtend((void *)pList, 
					((*nFound)+1)*sizeof(Int4), (*nFound)*sizeof(Int4));
			}
			if (!pList) {*nFound=0;return(NULL);}
			pList[*nFound]=i+1;
			(*nFound)++;
		}
	}
	return(pList);
}

/*******************************************************************************

  Function : DDV_MSG_SELECT()
  
  Purpose : manage OM_MSG_SELECT & OM_MSG_DESELECT message
  
  Return value : OM_MSG_RET_OK if success

*******************************************************************************/
static Int2 DDV_MSG_SELECT(OMMsgStructPtr ommsp,Boolean IsSelectMsg)
{
DdvMainPtr    dmp;  /*DDV panel data*/
OMUserDataPtr omudp;/*user data set when registering DDV panel*/
Int4Ptr       bsp_List=NULL;

Boolean       bRet;
SeqLocPtr     slp;
Int4          from_col,to_col,from_row,to_row,bsp_start,bsp_stop,nBspInstances=0,i,
              row_num;

	omudp = (OMUserDataPtr)(ommsp->omuserdata);
	if (omudp == NULL) return(OM_MSG_RET_ERROR);

	dmp=(DdvMainPtr)omudp->userdata.ptrvalue;
	if (dmp == NULL) return(OM_MSG_RET_ERROR);

	if (ommsp->itemtype != OBJ_BIOSEQ) return(OM_MSG_RET_OK);
	if (ommsp->regiontype!=OM_REGION_SEQLOC) return(OM_MSG_RET_OK);

	/*Am I concern by that Bioseq ?*/
	bRet=DDV_CheckBSPEntitiesInDisplay(&(dmp->MSA_d.pgp_l),ommsp->entityID,
		ommsp->itemID);
	
	if (bRet==FALSE) return(OM_MSG_RET_OK);  /* some other viewer */
	
	/*get the disp coord range*/
	DDV_GetCurrentDispRange(dmp->hWndDDV,&(dmp->GrData),dmp->MSA_d.pgp_l.LengthAli,
			&from_col,&to_col,&from_row,&to_row);

	/*get the occurence(s) (row(s) #) for that bioseq*/
	bsp_List=DDV_GetBspListGivenIDs(&(dmp->MSA_d.pgp_l),ommsp->entityID,
		ommsp->itemID,&nBspInstances);

	if (bsp_List==NULL || nBspInstances==0) return(OM_MSG_RET_ERROR);

	/*scan the SeqLoc list to figure out if it's currently on the screen*/
	slp=(SeqLocPtr)ommsp->region;
	while(slp){/*for each SeqLoc, I try to see what region(s) has (have) to be
		uptdated*/
		for(i=0;i<nBspInstances;i++){
			row_num=bsp_List[i]-1;
			bsp_start=DDV_GetDispCoordGivenBspCoord(dmp->MSA_d.pgp_l.TableHead[row_num],
				SeqLocStart(slp));
			bsp_stop=DDV_GetDispCoordGivenBspCoord(dmp->MSA_d.pgp_l.TableHead[row_num],
				SeqLocStop(slp));
			if (bsp_start!=(Int4)-1 && bsp_stop!=(Int4)-1)
				DDV_InvalRegion(dmp->hWndDDV,&(dmp->GrData),
					_max_(bsp_start,from_col),_min_(bsp_stop,to_col),
					bsp_List[i],IsSelectMsg);
		}
		slp=slp->next;
	}
	return(OM_MSG_RET_OK);
}

/*******************************************************************************

  Function : DDV_MSG_UPDATE_CaretPos()
  
  Purpose : update the position of the caret (edit mode only)
  
  Return value : OM_MSG_RET_OK if success

*******************************************************************************/
static Int2 DDV_MSG_UPDATE_CaretPos(OMMsgStructPtr ommsp)
{
DdvMainPtr    dmp;  /*DDV panel data*/
OMUserDataPtr omudp;/*user data set when registering DDV panel*/
Boolean       bRet;
Int4          from_col,to_col,from_row,to_row;

	omudp = (OMUserDataPtr)(ommsp->omuserdata);
	if (omudp == NULL) return(OM_MSG_RET_ERROR);

	dmp=(DdvMainPtr)omudp->userdata.ptrvalue;
	if (dmp == NULL) return(OM_MSG_RET_ERROR);
	
	/*Am I concern by that Bioseq ?*/
	bRet=DDV_CheckBSPEntitiesInDisplay(&(dmp->MSA_d.pgp_l),ommsp->entityID,
		ommsp->itemID);
	
	if (bRet==FALSE) return(OM_MSG_RET_ERROR);
	
	/*get the disp coord range*/
	DDV_GetCurrentDispRange(dmp->hWndDDV,&(dmp->GrData),dmp->MSA_d.pgp_l.LengthAli,
			&from_col,&to_col,&from_row,&to_row);
	
	/*switch to one-bsed value because dmp->dci values are one-based values*/
	from_row++;
	to_row++;
	/*hide the caret from old coordinates*/
	if (dmp->dci.old_col>=from_col && dmp->dci.old_col<=to_col &&
		dmp->dci.old_row>=from_row && dmp->dci.old_row<=to_row){
		DDV_InvalRegion(dmp->hWndDDV,&(dmp->GrData),
			dmp->dci.old_col,dmp->dci.old_col+1,
			dmp->dci.old_row,FALSE);
	}
	/*show caret on new coordinates*/
	if (dmp->dci.new_col>=from_col && dmp->dci.new_col<=to_col &&
		dmp->dci.new_row>=from_row && dmp->dci.new_row<=to_row){
		DDV_InvalRegion(dmp->hWndDDV,&(dmp->GrData),
			dmp->dci.new_col,dmp->dci.new_col+1,
			dmp->dci.new_row,FALSE);
	}
    return(OM_MSG_RET_OK);
}

/*******************************************************************************

  Function : DDV_MSG_UPDATE_DelBSP()
  
  Purpose : rebuilt a display after deletion of one BSP
  
  Return value : OM_MSG_RET_OK if success

*******************************************************************************/
static Int2 DDV_MSG_UPDATE_DelBSP(OMMsgStructPtr ommsp)
{
DdvMainPtr    dmp;  /*DDV panel data*/
OMUserDataPtr omudp;/*user data set when registering DDV panel*/

	omudp = (OMUserDataPtr)(ommsp->omuserdata);
	if (omudp == NULL) return(OM_MSG_RET_ERROR);

	dmp=(DdvMainPtr)omudp->userdata.ptrvalue;
	if (dmp == NULL) return(OM_MSG_RET_ERROR);
	
    /*delete the current display*/
	if (dmp->MSA_d.pgp_l.TableHead) 
		DDV_DeleteDisplayList(&dmp->MSA_d.pgp_l);
	if (dmp->MSA_d.pgp_l.RulerDescr) 
		ValNodeFreeData(dmp->MSA_d.pgp_l.RulerDescr);
	dmp->MSA_d.pgp_l.TableHead=NULL;
	dmp->MSA_d.pgp_l.RulerDescr=NULL;
	/*rebuild a new one*/
	if (!DDV_CreateDisplayFromIndex(dmp->MSA_d.pgp_l.sap,
		&(dmp->MSA_d.pgp_l),ParaG_Size,&(dmp->ddo))) 
		return(OM_MSG_RET_ERROR);

	/*build the Master Ruler descriptor*/
	dmp->MSA_d.pgp_l.RulerDescr=DDV_ComputeRuler(dmp->MSA_d.pgp_l.sap,&(dmp->ddo));

	/*delete old tables*/
    if(dmp->MasterViewer == SAMVIEWCN3D) {
        if (!DDV_InitColour_When_Start(dmp->MSA_d.pgp_l.sap,
            &(dmp->MSA_d.pgp_l),&(dmp->Globals.colorp), TRUE)){
            dmp->ddo.bUseColors=FALSE;
        }
    }
    else {
        DDV_ClearColor(dmp->Globals.colorp);
        /*build new tables*/
        if (!DDV_InitColour_When_Start(dmp->MSA_d.pgp_l.sap,
            &(dmp->MSA_d.pgp_l),&(dmp->Globals.colorp), FALSE)){
            dmp->ddo.bUseColors=FALSE;
        }
    }
    
    DDV_SortPGPLineNum(dmp->MSA_d.pgp_l.TableHead,dmp->MSA_d.pgp_l.nBsp);
	DDV_WhatSize(dmp);
	DDV_SetupWin(dmp->hWndDDV,TRUE,NULL);
	Update();

    return(OM_MSG_RET_OK);
}

/*******************************************************************************

  Function : DDV_MSG_UPDATE_Layout()
  
  Purpose : manage the OM_MSG_UPDATE message
  
  Return value : OM_MSG_RET_OK if success

*******************************************************************************/
static Int2 DDV_MSG_UPDATE_Layout(PaneL hWndDDV, DDVUpdateLayoutDataPtr dumdp)
{
DdvMainPtr  dmp;
Int4        j;
Boolean     bResetScrolls=FALSE,
			bRebuildDisplay=FALSE,
			bSwitchColors=FALSE,
			bUpdateColors=FALSE;

	if (!hWndDDV || !dumdp)
		return(OM_MSG_RET_OK);
	
	dmp = (DdvMainPtr) GetObjectExtra(hWndDDV);
	
	if (!dmp)
		return(OM_MSG_RET_OK);

	if (!Visible(hWndDDV))
		return(OM_MSG_RET_OK);
	
	
	dmp->ddo.SpacerSize=dumdp->SpacerSize;

	/*rebuild the display, if needed*/	
	if ((dmp->ddo.DispDiscStyle!=dumdp->DispDiscStyle) ||
		(dmp->ddo.SpacerSize!=dumdp->SpacerSize) ||
		(dmp->ddo.DiscJustification!=dumdp->DiscJustification) ||
		(dmp->ddo.ShowLeftTail!=dumdp->ShowLeftTail)||
		(dmp->ddo.ShowRightTail!=dumdp->ShowRightTail)){
		
		if((dmp->ddo.ShowLeftTail!=dumdp->ShowLeftTail)||
			(dmp->ddo.ShowRightTail!=dumdp->ShowRightTail)||
			(dmp->ddo.DispDiscStyle!=dumdp->DispDiscStyle)){
			bUpdateColors=TRUE;
		}
		/*update the styles*/
		dmp->ddo.DispDiscStyle=dumdp->DispDiscStyle;
		dmp->ddo.DiscJustification=dumdp->DiscJustification;
		dmp->ddo.ShowLeftTail=dumdp->ShowLeftTail;
		dmp->ddo.ShowRightTail=dumdp->ShowRightTail;
		
		/*delete the current display*/
		if (dmp->MSA_d.pgp_l.TableHead) 
			DDV_DeleteDisplayList(&dmp->MSA_d.pgp_l);
		if (dmp->MSA_d.pgp_l.RulerDescr) 
			ValNodeFreeData(dmp->MSA_d.pgp_l.RulerDescr);
		dmp->MSA_d.pgp_l.TableHead=NULL;
		dmp->MSA_d.pgp_l.RulerDescr=NULL;
		/*rebuild a new one*/
		if (!DDV_CreateDisplayFromIndex(dmp->MSA_d.pgp_l.sap,
			&(dmp->MSA_d.pgp_l),ParaG_Size,&(dmp->ddo))) 
			return(OM_MSG_RET_ERROR);

		/*build the Master Ruler descriptor*/
		dmp->MSA_d.pgp_l.RulerDescr=DDV_ComputeRuler(dmp->MSA_d.pgp_l.sap,&(dmp->ddo));

		/*init the colours*/
		if(bUpdateColors){
            if( dmp->MasterViewer == SAMVIEWCN3D ) {
                if (!DDV_InitColour_When_Start(dmp->MSA_d.pgp_l.sap,
                    &(dmp->MSA_d.pgp_l),&(dmp->Globals.colorp),TRUE)){
                    dmp->ddo.bUseColors=FALSE;
                }
            }
            else {
                /*delete old tables*/
                DDV_ClearColor(dmp->Globals.colorp);
                /*build new tables*/
                if (!DDV_InitColour_When_Start(dmp->MSA_d.pgp_l.sap,
                    &(dmp->MSA_d.pgp_l),&(dmp->Globals.colorp),FALSE)){
                    dmp->ddo.bUseColors=FALSE;
                }
            }
        }
        bResetScrolls=TRUE;
		bRebuildDisplay=TRUE;
	}

	/*set the ruler styles*/
	if(dumdp->nSeq && dumdp->SeqList){
		for(j=0;j<dumdp->nSeq;j++){
			switch(dumdp->RulerStyle){
				case SCALE_POS_TOP:
					DDV_SetRulerAttribInPGP(dmp->MSA_d.pgp_l.TableHead[dumdp->SeqList[j]-1], 
							SCALE_POS_TOP);
					break;
				case SCALE_POS_NONE:
					DDV_SetRulerAttribInPGP(dmp->MSA_d.pgp_l.TableHead[dumdp->SeqList[j]-1], 
							SCALE_POS_NONE);
					break;
				default:
					DDV_SetRulerAttribInPGP(dmp->MSA_d.pgp_l.TableHead[dumdp->SeqList[j]-1], 
							SCALE_POS_NONE);
					break;
			}			
		}
		bRebuildDisplay=TRUE;
	}

	/*colors*/
	if(dmp->ddo.bUseColors!=dumdp->bUseColors){
		dmp->ddo.bUseColors=dumdp->bUseColors;
		bSwitchColors=TRUE;
	}

	/*redraw if needed*/
	if (bRebuildDisplay){
	    DDV_SortPGPLineNum(dmp->MSA_d.pgp_l.TableHead,dmp->MSA_d.pgp_l.nBsp);
	    DDV_WhatSize(dmp);
	    DDV_SetupWin(hWndDDV,bResetScrolls,NULL);
	    Update();
	}
	else if (bSwitchColors){
		RecT   rcP;
		WindoW temport;
		
		temport=SavePort(ParentWindow(hWndDDV));
		Select(hWndDDV);
		
		ObjectRect(hWndDDV,&rcP);
		InvalRect(&rcP);
		RestorePort(temport);
		/*Update();*/
	}


	return(OM_MSG_RET_OK);
}

/*******************************************************************************

  Function : DDV_MSG_FLUSH()
  
  Purpose : kill the viewer in response to a OM_MSG_FLUSH message
  
  Return value : OM_MSG_RET_OK if success

*******************************************************************************/
static Int2 DDV_MSG_FLUSH(OMMsgStructPtr ommsp)
{
DdvMainPtr    dmp;  /*DDV panel data*/
OMUserDataPtr omudp;/*user data set when registering DDV panel*/

	omudp = (OMUserDataPtr)(ommsp->omuserdata);
	if (omudp == NULL) return(OM_MSG_RET_ERROR);

	dmp=(DdvMainPtr)omudp->userdata.ptrvalue;
	if (dmp == NULL) return(OM_MSG_RET_ERROR);
	
	if (dmp->MSA_d.entityID==ommsp->entityID &&
		((dmp->MSA_d.itemID==ommsp->itemID && ommsp->itemtype==OBJ_SEQALIGN)
         ||(ommsp->itemID == 0 && ommsp->itemtype==0))){
		Remove(dmp->hParent);
	}
	return(OM_MSG_RET_OK);
}


/*******************************************************************************

  Function : DDV_OM_MsgFunc()
  
  Purpose : ObjMgr message loop of the DDV viewer/editor
  
  Parameters : see Toolkit
  
  Return value : see Toolkit 

*******************************************************************************/
static Int2 LIBCALLBACK DDV_OM_MsgFunc (OMMsgStructPtr ommsp)
{
Int2 nRet=OM_MSG_RET_OK;
  
	switch (ommsp->message)
	{
		case OM_MSG_DEL:
			break;
		case OM_MSG_CREATE:
			break;
		case OM_MSG_UPDATE:{
			DDVUpdateMSGPtr        dump;
			
			dump = (DDVUpdateMSGPtr)(ommsp->procmsgdata);
			
			if (!dump)
				break;
            switch(dump->type){
                case UPDATE_TYPE_LAYOUT:{
                    DDVUpdateLayoutDataPtr dumdp;
                
                    dumdp=(DDVUpdateLayoutDataPtr)dump->data;
                
                    if (!dumdp)
                        break;
                
                    nRet=DDV_MSG_UPDATE_Layout(dumdp->ddv_panel,dumdp);
                    break;
                }
                case UPDATE_TYPE_EDIT_DELBSP:
                    nRet=DDV_MSG_UPDATE_DelBSP(ommsp);
                    break;
				case UPDATE_TYPE_CARETPOS:
					nRet=DDV_MSG_UPDATE_CaretPos(ommsp);
					break;
                default :{
		            RecT   rcP;
		            WindoW temport;
                    DdvMainPtr    dmp;  /*DDV panel data*/
                    OMUserDataPtr omudp;/*user data set when registering DDV panel*/

	                omudp = (OMUserDataPtr)(ommsp->omuserdata);
	                if (omudp == NULL) return(OM_MSG_RET_ERROR);

	                dmp=(DdvMainPtr)omudp->userdata.ptrvalue;
	                if (dmp == NULL) return(OM_MSG_RET_ERROR);
		            
		            temport=SavePort(ParentWindow(dmp->hWndDDV));
		            Select(dmp->hWndDDV);
		            
		            ObjectRect(dmp->hWndDDV,&rcP);
		            InvalRect(&rcP);
		            RestorePort(temport);
                }
            }
			
			break;
		}
		case OM_MSG_SELECT:
			nRet=DDV_MSG_SELECT(ommsp,TRUE);
			break;
		case OM_MSG_DESELECT:
			nRet=DDV_MSG_SELECT(ommsp,FALSE);
			break;
		case OM_MSG_CACHED:
			break;
		case OM_MSG_UNCACHED:
			break;
		case OM_MSG_TO_CLIPBOARD:
			break;
		case OM_MSG_SETCOLOR:
			break;
		case OM_MSG_FLUSH:
			DDV_MSG_FLUSH(ommsp);
			break;
		default:
			break;
	}

	return (nRet);
}

/*******************************************************************************

  Function : DDV_RegMsgFuncForBsp()
  
  Purpose : given the display data structure (nBsp & TableHead), register a
            MsgFunc callback for each bioseq of a SeqAlign.
  
  Parameters : TableHead; entry point of the display data structure (ParaG list)
			   nBsp; number of sequences in the SeqAlign
			   master_eID; entityID of the object containing the SeqAlign
			   procID; identifier of the viewer (usually DDV)
			   userkey; userkey of the viewer (usually DDV)
			   
  Return value : -

*******************************************************************************/
static void DDV_RegMsgFuncForBsp(ValNodePtr PNTR TableHead,DdvMainPtr dmp,Int4 nBsp, 
	Uint2 master_eID, Uint2 procID, Uint2 procType, Uint2 userkey)
{
OMUserDataPtr omudp,omudp_tmp;
ParaGPtr      pgp;
BioseqPtr     bsp;
Int4          i;
Boolean       bFound;
Uint2         bsp_eID;

	if (!TableHead || nBsp==0 || master_eID==0 || procID==0 
		|| procType==0 ||userkey==0) return;
		
	for(i=0;i<nBsp;i++){
		pgp=(ParaGPtr)(dmp->MSA_d.pgp_l.TableHead[i]->data.ptrvalue);
		bsp=BioseqLockById(pgp->sip);
		if (bsp){
			bsp_eID=ObjMgrGetEntityIDForPointer((Pointer)bsp);
			if (bsp_eID>0 && bsp_eID!=master_eID){
				omudp=ObjMgrGetUserData(bsp_eID,procID,procType,userkey);
				if (omudp){
					/*scan the user data associated with that bioseq and try to
					find if a Msg Func is already attached for DDV*/
					omudp_tmp=omudp;
					bFound=FALSE;
					while(omudp_tmp){
						if (omudp_tmp->procid==procID && 
							omudp_tmp->proctype==procType && 
							omudp_tmp->userkey==userkey && 
							omudp_tmp->messagefunc==DDV_OM_MsgFunc){
							
							bFound=TRUE;
							break;
						}
						omudp_tmp=omudp_tmp->next;
					}
				}
				else{/*add the Msg Callaback and DDV main data block to the bioseq*/
					bFound=FALSE;
				}
				if (bFound==FALSE){
    				omudp = ObjMgrAddUserData (bsp_eID,procID,procType,userkey);
    				if (omudp != NULL) {
    					omudp->messagefunc = DDV_OM_MsgFunc;
						omudp->userdata.ptrvalue = (Pointer)dmp;
					}
				}
			}
			BioseqUnlockById(pgp->sip);
		}
	}
}

/*******************************************************************************

  Function : DDV_StartPanel_Slave()
  
  Purpose : start DDV main window as a slave module: init 'dmp' data structure
  			(if needed), create the DDV panel, prepare a display from 'sap'.
  
  Parameters : 	sap; pointer to a SeqAlign
  				input_entityID,input_itemID; identity of the SeqAlign
  				dmp; main data of DDV. May be NULL
				hParent; parent window for DDV. Must not be NULL
  				bSetAppProp; if TRUE, use SetAppProperty().
				bInitGraph; if TRUE, initialize the graphical data for DDV
				bEditor; if TRUE, start the editor mode. Otherwise, a simple
						viewer will ba started
                rcPp; the RecT that defines the size of the panel.  can be NULL.
						
  Note :  This function assumes that DDV will be placed inside an existing 
     parent window provided by the host application.
	 
  Return value : handle to the DDV panel

*******************************************************************************/
static Boolean DDV_StartPanel_Slave(DdvMainPtr dmp,SeqAlignPtr sap,
		Uint2 input_entityID,Uint2 input_itemID,Boolean bEditor)
{
DDV_ColorGlobal * dcgp;
	
	WatchCursor();

	/*get the SeqAlign and build a default display structure*/
	if (!AlnMgrIndexSeqAlign(sap))
		return (FALSE);

#ifdef _DDVEDIT
    dmp->editGlobal = SAM_NewEditGlobal(TRUE, TRUE, TRUE, sap);
#endif /* _DDVEDIT */

	if (!DDV_CreateDisplayFromIndex(sap,&(dmp->MSA_d.pgp_l),ParaG_Size,
		&(dmp->ddo))) return(FALSE);
	dmp->MSA_d.pgp_l.entitiesTbl=
		DDV_BuildBspEntitiesTbl(dmp->MSA_d.pgp_l.TableHead,
		dmp->MSA_d.pgp_l.nBsp);
	if (!dmp->MSA_d.pgp_l.entitiesTbl) return(FALSE);
	
	/*build the Master Ruler descriptor*/
	dmp->MSA_d.pgp_l.RulerDescr=DDV_ComputeRuler(sap,&(dmp->ddo));
	
	dmp->MSA_d.entityID=input_entityID;
	dmp->MSA_d.itemID=input_itemID;
		
	/*get the color stuff from the sap object. Usually I can do that
	when DDV is loaded from Cn3D*/
	dcgp=DDV_GetColorGlobalEx((Pointer) sap);

	/*init the colours*/
	if (dcgp==NULL){
		if (!DDV_InitColour_When_Start(dmp->MSA_d.pgp_l.sap,&(dmp->MSA_d.pgp_l),
			&(dmp->Globals.colorp),FALSE)){
			dmp->ddo.bUseColors=FALSE;
		}
	}
    else dmp->Globals.colorp = dcgp;

    if(dmp->MasterViewer == SAMVIEWCN3D) {
        if (!DDV_InitColour_When_Start(dmp->MSA_d.pgp_l.sap,&(dmp->MSA_d.pgp_l),
            &(dmp->Globals.colorp),TRUE)) dmp->ddo.bUseColors=FALSE;
    }
        

	/*set the initial position of the caret (edit mode)*/
	dmp->dci.old_row=1;
	dmp->dci.new_row=1;
	dmp->dci.old_col=0;
	dmp->dci.new_col=0;
	
	/*relation between align size and display type; compute nTotLines for example*/
	DDV_WhatSize(dmp);

	ArrowCursor();
	DDV_SetupWin (dmp->hWndDDV,TRUE,NULL);
	/*Update();*/
	return(TRUE);
}


/*******************************************************************************

  Function : DDV_ObjRegMasterDDV()
  
  Purpose : call by ObjMgr to start DDV panel
  
  Parameters : see Toolkit
  
  Note : this callback is only used by the standalone DDV. Never use this
         callback for other purpose.
  
  Return value : see Toolkit 

*******************************************************************************/
extern Int2 LIBCALLBACK DDV_ObjRegMasterDDV (Pointer data)
{
OMProcControlPtr    ompcp;
SeqAlignPtr			sap=NULL;
OMUserDataPtr       omudp;
DdvMainWinPtr 		mWin_d;
DdvMainPtr          dmp;

	/*retrieve data*/
	ompcp = (OMProcControlPtr) data;

	if (ompcp == NULL || ompcp->proc == NULL) {
		ErrPostEx (SEV_ERROR, 0, 0, "Data NULL [1]");
		return OM_MSG_RET_ERROR;
	}
	/*only accept a SeqALign*/
	switch (ompcp->input_itemtype) {
		case OBJ_SEQALIGN :
			sap = (SeqAlignPtr) ompcp->input_data;
			break;
		default :
			return OM_MSG_RET_ERROR;
	}
	/*oups... nothing to deal with*/
	if (sap == NULL) {
		ErrPostEx (SEV_ERROR, 0, 0, "Data NULL [2]");
		return OM_MSG_RET_ERROR;
	}

	/*get the DDV main window data*/
	if (!g_hParent) goto error;
	mWin_d = (DdvMainWinPtr) GetObjectExtra (g_hParent);
	
	/*build the DDV panel (its parent is g_hParent)*/
	dmp = (DdvMainPtr) GetObjectExtra(mWin_d->hWndDDV);
	if (!DDV_StartPanel_Slave(dmp,sap,ompcp->input_entityID,
		ompcp->input_itemID,TRUE)) goto error;
	
	/*attach a Msg Func on the current seqalign*/
	dmp->userkey=OMGetNextUserKey();
	dmp->procid=ompcp->proc->procid;
	dmp->proctype=ompcp->proc->proctype;
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid, 
		OMPROC_EDIT,dmp->userkey);
    if (omudp == NULL) 
		goto error;
    omudp->messagefunc =  DDV_OM_MsgFunc;
	omudp->userdata.ptrvalue = (Pointer)dmp;

	/*attach a Msg Func on each Bioseq of the current seqalign*/
	DDV_RegMsgFuncForBsp(dmp->MSA_d.pgp_l.TableHead,dmp,dmp->MSA_d.pgp_l.nBsp, 
		ompcp->input_entityID, dmp->procid,dmp->proctype,dmp->userkey);

	DDV_EnableGotoTBItems(dmp->hParent,TRUE);
	return(OM_MSG_RET_DONE);

error:
	/*show the logo window*/
	mWin_d->Show_logo=TRUE;
	Update ();
	ArrowCursor();
	return(OM_MSG_RET_ERROR);
}


/*******************************************************************************

  Function : DDV_ObjRegSlaveDDV()
  
  Purpose : call by DDV_ObjRegSlaveEditDDV() & DDV_ObjRegSlaveViewDDV(), below
  
*******************************************************************************/
static Int2 DDV_ObjRegSlaveDDV(Pointer data,Boolean bEditor)
{
OMProcControlPtr    ompcp;
SeqAlignPtr			sap=NULL;
WindoW              hParent=NULL;
DdvMainWinPtr 		mWin_d;
OMUserDataPtr       omudp = NULL;
DdvMainPtr          dmp;
SAM_ViewGlobal      *vgp;


	/*retrieve data*/
	ompcp = (OMProcControlPtr) data;
    if (ompcp == NULL || ompcp->proc == NULL) {
		ErrPostEx (SEV_ERROR, 0, 0, "Data NULL [1]");
		return OM_MSG_RET_ERROR;
	}

	/*only accept a SeqALign*/
	switch (ompcp->input_itemtype) {
		case OBJ_SEQALIGN :
			sap = (SeqAlignPtr) ompcp->input_data;
			break;
		default :
			return OM_MSG_RET_ERROR;
	}

	/*oups... nothing to deal with*/
	if (sap == NULL) {
		ErrPostEx (SEV_ERROR, 0, 0, "Data NULL [2]");
		return OM_MSG_RET_ERROR;
	}
        
	/*start the slave module*/
    /* get parameters from launching program */
    vgp = GetAppProperty(SAM_ViewString);
    if(vgp == NULL) hParent=DDV_StartMainWin_Slave(NULL);
    else hParent=DDV_StartMainWin_Slave(vgp);

	if (hParent){
		RealizeWindow(hParent);
		Show(hParent);
		mWin_d=(DdvMainWinPtr)GetObjectExtra(hParent);
		dmp = (DdvMainPtr) GetObjectExtra(mWin_d->hWndDDV);
        if(vgp != NULL) dmp->MasterViewer = vgp->MasterViewer;
        if (!DDV_StartPanel_Slave(dmp,sap,ompcp->input_entityID,
            ompcp->input_itemID,bEditor)) 
            goto error;
        
		/*when the user will close DDV, DDV won't delete any data excepts
		its internal data structure*/
		mWin_d->dod.choice=DDV_OPENTYPE_NOTRESP;
		
		/*attach a Msg Func on the current seqalign*/
		dmp->userkey=OMGetNextUserKey();
		dmp->procid=ompcp->proc->procid;
		dmp->proctype=ompcp->proc->proctype;
    	omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid, 
			(bEditor==TRUE ? OMPROC_EDIT : OMPROC_VIEW),dmp->userkey);
    	if (omudp == NULL) 
			goto error;
    	omudp->messagefunc =  DDV_OM_MsgFunc;
		omudp->userdata.ptrvalue = (Pointer)dmp;
		/*attach a Msg Func on each Bioseq of the current seqalign*/
		DDV_RegMsgFuncForBsp(dmp->MSA_d.pgp_l.TableHead,dmp,dmp->MSA_d.pgp_l.nBsp, 
			ompcp->input_entityID, dmp->procid,dmp->proctype,dmp->userkey);
		DDV_EnableGotoTBItems(hParent,TRUE);
	} 
	else {
		goto error;
	}
	
	return(OM_MSG_RET_DONE);

error:
	if (hParent) Remove(hParent);
	return(OM_MSG_RET_ERROR);
}

/*******************************************************************************

  Function : DDV_ObjRegSlaveEditDDV()
  
  Purpose : call by ObjMgr to start DDV as a slave, editor mode on a SeqAlign
  
  Parameters : see Toolkit
  
  Return value : see Toolkit 

*******************************************************************************/
extern Int2 LIBCALLBACK DDV_ObjRegSlaveEditDDV (Pointer data)
{
	return(DDV_ObjRegSlaveDDV(data,TRUE));
}

/*******************************************************************************

  Function : DDV_ObjRegSlaveViewDDV()
  
  Purpose : call by ObjMgr to start DDV as a slave, editor mode on a SeqAlign
  
  Parameters : see Toolkit
  
  Return value : see Toolkit 

*******************************************************************************/
extern Int2 LIBCALLBACK DDV_ObjRegSlaveViewDDV (Pointer data)
{
	return(DDV_ObjRegSlaveDDV(data,FALSE));
}

/*******************************************************************************

  Function : DDV_GetAndCheckSeqAlign()
  
  Purpose : given a file OR a GI OR a SeqEntry , retrieve SeqAlign(s) 
  
  Parameters :  fp; a handle of an open file
  				gi; a GI numger
				sep2; a SeqEntry
				fetchSepProc; fetch procedure of the user (only if GI!=0)
				
  Return value : a list of nodes. For each node, data.intvalue is an entityID
           of a SeqAlign 

*******************************************************************************/
extern ValNodePtr DDV_GetAndCheckSeqAlign(FILE *fp,Int4 gi,SeqEntryPtr sep2,
	UdvFetchSeqEntryProc fetchSepProc,DdvOpenDataPtr dodp,Uint2Ptr entityID)
{
SeqEntryPtr sep=NULL;
ValNodePtr  head = NULL,vnp_ali=NULL;
Pointer     dataptr=NULL;
Uint2       datatype=0;
Uint2		nRet=DVV_MSG_M_OK;

	*entityID=0;
	if (fp!=NULL){/*open a file*/
		/*Read the file (it can contain several objects)*/
		while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, 
			FALSE, TRUE, FALSE)) != NULL) {
				ValNodeAddPointer (&head, datatype, dataptr);
		}
		/*scan the first node of the head list*/
		datatype = head->choice;
		dataptr = head->data.ptrvalue;
        *entityID=ObjMgrRegister(datatype,dataptr);
		AssignIDsInEntity (*entityID,0,NULL);
		switch(datatype){
			case OBJ_SEQENTRY:
			case OBJ_SEQALIGN:
			case OBJ_SEQANNOT:
			case OBJ_BIOSEQ:
			case OBJ_BIOSEQSET:
				DDV_GetSeqAlign(dataptr,datatype,&vnp_ali);
				break;
			default:
				break;
		}
		if (dodp){
			dodp->choice=DDV_OPENTYPE_FILE;
			dodp->vnp=head;
		}
	} else if (gi!=0){/*retrieve a Network Entry given a GI*/
		if (fetchSepProc) {/*fetch the SeqEntry*/
			sep=fetchSepProc(gi, 0);
			if (!sep) {
				nRet=DVV_MSG_O_E_READGI;
				goto error;
			}
			*entityID=ObjMgrGetEntityIDForPointer((Pointer)sep);
			AssignIDsInEntity (*entityID,0,NULL);
			SeqEntryExplore(sep,(Pointer)&vnp_ali,DDV_SearchAli);
			if (dodp){
				dodp->choice=DDV_OPENTYPE_GI;
				dodp->sep=sep;
			}
		} else {
			nRet=DVV_MSG_O_E_NOFETCHFUNC;
			goto error;
		}
	} else if (sep2!=NULL){/*analyse a SeqEntry*/
		*entityID=ObjMgrGetEntityIDForPointer((Pointer)sep2);
		AssignIDsInEntity (*entityID,0,NULL);
		SeqEntryExplore(sep2,(Pointer)&vnp_ali,DDV_SearchAli);
		if (dodp){
			dodp->choice=DDV_OPENTYPE_SEP;
			dodp->sep=sep2;
		}
	} else {
		nRet=DVV_MSG_O_E_NOTHINGTODO;
		goto error;
	}
	
	/*register for ObjMgr*/
	if (!vnp_ali) {
		nRet=DVV_MSG_O_E_BADTYPE;
		goto error;
	}

	
	return(vnp_ali);

error:
	/*to add : Msg Handler with nRet*/
	return(NULL);
}

/*******************************************************************************

  Function : DDV_FOpenTextProc()
  
  Purpose : manage File name edit control of the FileOpen dialog box 
  
  Parameters : t; edit control
  
  Return value : none 

*******************************************************************************/
static void DDV_FOpenTextProc(TexT t)
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

  Function : DDV_FOpenBrowseProc()
  
  Purpose : manage browse button of the FileOpen dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void DDV_FOpenBrowseProc(ButtoN g)
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
		DDV_FOpenTextProc(dfodp->FNameEditCtrl);
	}

	return;   
}

/*******************************************************************************

  Function : DDV_FOpenCancelProc()
  
  Purpose : manage cancel button of the FileOpen dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void DDV_FOpenCancelProc(ButtoN g)
{
WindoW				hOpenDlg;
DlgFileOpenDataPtr  dfodp;
DdvMainWinPtr 		mWin_d;

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;
	dfodp = (DlgFileOpenDataPtr) GetObjectExtra (hOpenDlg);
	
	if (dfodp){
		mWin_d = (DdvMainWinPtr) GetObjectExtra (dfodp->parent);
		if (mWin_d) Enable(mWin_d->MainMenu.FileOpen);
	}
	
	Remove(hOpenDlg);
}

/*******************************************************************************

  Function : DDV_FOpenAcceptProc()
  
  Purpose : manage ok button of the FileOpen dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void DDV_FOpenAcceptProc(ButtoN g)
{
WindoW				hOpenDlg;
DlgFileOpenDataPtr  dfodp;
DdvMainPtr          dmp;
DdvMainWinPtr 		mWin_d;
Char 				szFName[PATH_MAX]={""};
Boolean				isBinary;
FILE				*fp;
ValNodePtr 			vnp_ali=NULL;
Uint2 				nRet=DVV_MSG_M_OK,The_entityID;

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;
	dfodp = (DlgFileOpenDataPtr) GetObjectExtra (hOpenDlg);
	
	if (dfodp){

		mWin_d = (DdvMainWinPtr) GetObjectExtra (dfodp->parent);
		if (mWin_d!=NULL){
			/*try to connect Entrez; if failure, still ok : the user can
			use DDV with a local file. Otherwise he/she will be in trouble ;-)*/
			if (mWin_d->AutonomeViewer && mWin_d->NetStartProc){
				mWin_d->NetStartProc(mWin_d->UseNetwork);
			}

			Enable(mWin_d->MainMenu.FileOpen);
			/*file name*/
			GetTitle(dfodp->FNameEditCtrl, szFName, sizeof(szFName));

			/*file reading mode*/
			if (GetValue(dfodp->ReadMode)==2) isBinary=TRUE;
			else isBinary=FALSE;
	
			Remove(hOpenDlg);

			/* open the i/o files in the correct mode */
			if ((fp = FileOpen (szFName, isBinary?"rb":"r")) == NULL){
				nRet=DVV_MSG_O_E_OPENFILEFAIL;
				goto error;
			}
			/*restrieve a list of registered seqalign(s)*/
			WatchCursor();
			/*delete old data, if needed*/
			dmp=(DdvMainPtr)GetObjectExtra(mWin_d->hWndDDV);
			DDV_CleanupDDVPdata_g(dmp);
			DDV_CloseData(mWin_d,FALSE);
			
			vnp_ali=DDV_GetAndCheckSeqAlign(fp,0,NULL,mWin_d->fetchSepProc,&(mWin_d->dod),
				&The_entityID);
			ArrowCursor();
			if (vnp_ali && The_entityID!=0){
				Uint2 entityID,itemID;
				SeqAlignPtr sap;
								
				/*get the first SeqAlign in the list and load the viewer panel*/
				mWin_d->vnp_ali=vnp_ali;
				
				if(mWin_d->Show_logo){
					mWin_d->Show_logo=FALSE;
				}
				sap=(SeqAlignPtr)vnp_ali->data.ptrvalue;
				entityID=sap->idx.entityID;
				itemID=sap->idx.itemID;
				GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, itemID, 
					OBJ_SEQALIGN, 0, 0, OBJ_SEQALIGN, 0);
			}
			FileClose(fp);
		}
	}
	
error:	
	return;
}

/*******************************************************************************

  Function : DDV_OpenFile()
  
  Purpose : callback of the File|Open command 
  
  Parameters : i; command item
  
  Return value : none 

*******************************************************************************/
extern void DDV_OpenFile(IteM i)
{
DlgFileOpenDataPtr  dfodp;
DdvMainWinPtr 		mWin_d;
WindoW				hWinMain,hOpenDlg;
GrouP				g,g4;
ButtoN				b1;
TexT				t1;    

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	mWin_d = (DdvMainWinPtr) GetObjectExtra (hWinMain);

	if (mWin_d==NULL) return;
		
	dfodp=(DlgFileOpenDataPtr)MemNew(sizeof(DlgFileOpenData));
	if (!dfodp) return;
	MemSet(dfodp,0,sizeof(DlgFileOpenData));

    hOpenDlg = FixedWindow(-30, -20,  -10,  -10, 
				"DDV - Open a local file",  NULL);
    g = NormalGroup(hOpenDlg, 2, 1, "File name:",  systemFont, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 10, 20);  
    t1 = DialogText(g,"",25, DDV_FOpenTextProc);
    PushButton(g, " browse...", DDV_FOpenBrowseProc);
   
    g = HiddenGroup(hOpenDlg, 3, 1, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 30, 30);
 
    g4 = NormalGroup(g, 2, 1, "File type", systemFont,  NULL);
    SetGroupMargins(g4, 10, 10);
    RadioButton(g4, "Ascii");
    RadioButton(g4, "Binary");
    SetValue(g4, 1);

    b1 = DefaultButton(g, "OK",  DDV_FOpenAcceptProc);
    PushButton(g, "Cancel",  DDV_FOpenCancelProc);
   
    Disable(b1);

	dfodp->parent=hWinMain;	
	dfodp->FNameEditCtrl=t1;
	dfodp->Ok=b1;
	dfodp->ReadMode=g4;

	SetObjectExtra (hOpenDlg, (Pointer) dfodp, StdCleanupExtraProc);
	
    Select(hOpenDlg);
    Show(hOpenDlg);

	/*allow only one FOpen dlg at a time*/
    /*UDV_set_PullMenus(&vmp->MainMenu,FALSE);*/
}

/*******************************************************************************

  Function : DDV_NetOpen_cancelProc()
  
  Purpose : manage ok button of the Network open dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void DDV_NetOpen_cancelProc(ButtoN g)
{
WindoW			hOpenDlg;
DDVNetOpenPtr 	dnop;

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;
	
	dnop = (DDVNetOpenPtr) GetObjectExtra (hOpenDlg);

	if (dnop==NULL) return;

	Remove(hOpenDlg);
}

/*******************************************************************************

  Function : DDV_IsInt()
  
  Purpose : test if an entry access code is a GI number, i.e. an Int
  
  Parameters : t; edit control
  
  Return value : none 

*******************************************************************************/
static Boolean DDV_IsInt(CharPtr szAccess)
{
Int2 i=0;

	while(szAccess[i]){
		if (!IS_DIGIT(szAccess[i])){
			Message (MSG_OK, "Not a GI number.");
			return(FALSE);
		}
		i++;
	}
	return(TRUE);
}

/*******************************************************************************

  Function : DDV_NetOpen_okProc()
  
  Purpose : manage ok button of the Network open dialog box 
  
  Parameters : g; button
  
  Return value : none 

*******************************************************************************/
static void DDV_NetOpen_okProc(ButtoN g)
{
WindoW			hOpenDlg;
DDVNetOpenPtr 	dnop;
Char 			szAccess[50]={""};
Int2			AccessType;
Int4			uid=0;
DdvMainWinPtr 	mWin_d;
ValNodePtr      vnp_ali;
MonitorPtr      mon;
DdvMainPtr      dmp;
Uint2           The_entityID;

	hOpenDlg=(WindoW)ParentWindow(g);

	if (!hOpenDlg) return;
	
	dnop = (DDVNetOpenPtr) GetObjectExtra (hOpenDlg);	

	if (dnop==NULL) return;

	mWin_d = (DdvMainWinPtr) GetObjectExtra (dnop->hWndMain);
	if (mWin_d==NULL) return;

	dmp=(DdvMainPtr)GetObjectExtra(mWin_d->hWndDDV);	
	
	Hide(hOpenDlg);
	
	/*retrieve a Gi*/
	GetTitle(dnop->Entry, szAccess, sizeof(szAccess));

	if (StringHasNoText (szAccess)) {
		Message (MSG_OK, "Please enter a GI number");
		Show (hOpenDlg);
		Select (hOpenDlg);
		Select (dnop->Entry);
		return;
	}

	AccessType=GetValue(dnop->AccessType);
	if (!StrToLong(szAccess,&uid)) 
		uid=0;
	/*Connect ID1 with a valid UID*/
	if (uid==0){
		ArrowCursor();
		Message (MSG_OK, "Unable to find your record in the database.");
		Show (hOpenDlg);
		Select (hOpenDlg);
		Select (dnop->Entry);
		return;
	}
	else{
		/*delete old data, if needed*/
		DDV_CleanupDDVPdata_g(dmp);
		DDV_CloseData(mWin_d,FALSE);
		
		WatchCursor();

		mon = MonitorStrNewEx (szAppName, 30, FALSE);
		MonitorStrValue (mon, "Retrieving your record...");
		Update ();

		vnp_ali=DDV_GetAndCheckSeqAlign(NULL,uid,NULL,mWin_d->fetchSepProc,
			&(mWin_d->dod),&The_entityID);

		ArrowCursor();

		MonitorFree (mon);

		if (vnp_ali && The_entityID!=0){
			Uint2 entityID,itemID;
			SeqAlignPtr sap;
			
			/*get the first SeqAlign in the list and load the viewer panel*/
			mWin_d->vnp_ali=vnp_ali;

			if(mWin_d->Show_logo){
				/*Hide(mWin_d->Logo_Panel);*/
				mWin_d->Show_logo=FALSE;
			}
			sap=(SeqAlignPtr)vnp_ali->data.ptrvalue;
			entityID=sap->idx.entityID;
			itemID=sap->idx.itemID;
			GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, itemID, 
				OBJ_SEQALIGN, 0, 0, OBJ_SEQALIGN, 0);
		}
	}

	Remove(hOpenDlg);
}

/*******************************************************************************

  Function : DDV_NetOpen_access()
  
  Purpose : manage Accession edit control of the Download dialog box 
  
  Parameters : t; edit control
  
  Return value : none 

*******************************************************************************/
static void DDV_NetOpen_access(TexT t)
{
Char 			szAccess[50]={""};
WindoW			hOpenDlg;
DDVNetOpenPtr 	dnop;

	hOpenDlg=(WindoW)ParentWindow(t);

	if (!hOpenDlg) return;
	
	dnop = (DDVNetOpenPtr) GetObjectExtra (hOpenDlg);

	if (dnop==NULL) return;
	
	GetTitle(t, szAccess, sizeof(szAccess));

	if (StringLen(szAccess) == 0)
		Disable(dnop->ok);
	else Enable(dnop->ok);

	return;
}

/*******************************************************************************

  Function : DDV_OpenNetwork()
  
  Purpose : callback of the File|Network command 
  
  Parameters : i; command item
  
  Return value : none 

*******************************************************************************/
extern void DDV_OpenNetwork(IteM i)
{
DdvMainWinPtr 		mWin_d;
WindoW				hWinMain,hOpenDlg;
GrouP				g,c;
DDVNetOpenPtr       dnopp;

	hWinMain=(WindoW)ParentWindow(i);

	if (!hWinMain) return;
	
	mWin_d = (DdvMainWinPtr) GetObjectExtra (hWinMain);

	if (mWin_d==NULL) return;
		
	dnopp=(DDVNetOpenPtr)MemNew(sizeof(DDVNetOpen));
	if (!dnopp) return;

	/*init the Network, if needed*/
	if (mWin_d->AutonomeViewer && mWin_d->NetStartProc){
		if (!mWin_d->NetStartProc(mWin_d->UseNetwork)) return;
	}

    hOpenDlg = FixedWindow(-30, -20,  -10,  -10, 
				"DDV - Fetch a network entry",  NULL);
	SetGroupSpacing (hOpenDlg, 10, 10);
	
	dnopp->hWndMain=hWinMain;
	
	/*accesion*/
	g = NormalGroup (hOpenDlg, 3, 0, "Entry", systemFont,NULL);

#ifdef WIN_MAC
	dnopp->AccessType = PopupList (g, TRUE, NULL);
#endif

#ifndef WIN_MAC
	dnopp->AccessType = PopupList (g, FALSE, NULL);
#endif

	PopupItem(dnopp->AccessType,"PopSet");
	
	SetValue (dnopp->AccessType, 1);

	dnopp->Entry = DialogText (g, "", 10, DDV_NetOpen_access);

	c = HiddenGroup (hOpenDlg, 4, 0, NULL);
	SetGroupSpacing (c, 10, 2);
	dnopp->ok=DefaultButton (c, "Retrieve", DDV_NetOpen_okProc);
	Disable (dnopp->ok);
	PushButton (c, "Cancel", DDV_NetOpen_cancelProc);
	
	SetObjectExtra (hOpenDlg, (Pointer) dnopp, StdCleanupExtraProc);

	/*display*/
	AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
	RealizeWindow (hOpenDlg);
	Show (hOpenDlg);
	Select (hOpenDlg);
	Update ();
}

/*******************************************************************************

  Function : DDV_LaunchAlignEditor()
  
  Purpose : load DDV as an editor
  
  Parameters : entityID  OR  sap; identifier of a SeqAlign
  
  Note:  DO NOT use entityID and sap at a same time.
  
  Return value : none 

*******************************************************************************/
extern void DDV_LaunchAlignEditor (Uint2 entityID,SeqAlignPtr sap)
{
Uint2 itemID,options;
ObjMgrDataPtr   omdp;

	if (sap==NULL && entityID==0) return;
	
	if (sap) {/*given sap, get eID & iID*/
		entityID = ObjMgrRegister (OBJ_SEQALIGN, (Pointer) sap);
		if (entityID==0) return;
		itemID = GetItemIDGivenPointer (entityID, OBJ_SEQALIGN, (Pointer) sap);
	}
	else{/*given eID, get iID*/
		omdp = ObjMgrGetData (entityID);
		if (!omdp) return;
		itemID = GetItemIDGivenPointer (entityID, 
				OBJ_SEQALIGN, (Pointer) omdp->dataptr);

	}
	options = ObjMgrGetOptions(entityID);
	options |= OM_OPT_FREE_IF_NO_VIEW;
	ObjMgrSetOptions(options, entityID);
	GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, itemID, 
		OBJ_SEQALIGN, 0, 0, OBJ_SEQALIGN, 0);
}

/*******************************************************************************

  Function : DDV_LaunchAlignViewer()
  
  Purpose : load DDV as a simple viewer 
  
  Parameters : entityID  OR  sap; identifier of a SeqAlign
  
  Note:  DO NOT use entityID and sap at a same time.
  
  Return value : none 

*******************************************************************************/
extern void DDV_LaunchAlignViewer (Uint2 entityID,SeqAlignPtr sap)
{
Uint2 itemID,options;
ObjMgrDataPtr   omdp;

	if (sap==NULL && entityID==0) return;

	if (sap) {/*given sap, get eID & iID*/
		entityID = ObjMgrRegister (OBJ_SEQALIGN, (Pointer) sap);
		if (entityID==0) return;
		itemID = GetItemIDGivenPointer (entityID, OBJ_SEQALIGN, (Pointer) sap);
	}
	else{/*given eID, get iID*/
		omdp = ObjMgrGetData (entityID);
		if (!omdp) return;
		itemID = GetItemIDGivenPointer (entityID, 
				OBJ_SEQALIGN, (Pointer) omdp->dataptr);

	}
	options = ObjMgrGetOptions(entityID);
	options |= OM_OPT_FREE_IF_NO_VIEW;
	ObjMgrSetOptions(options, entityID);
	GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, itemID, 
		OBJ_SEQALIGN, 0, 0, OBJ_SEQALIGN, 0);
}

/*******************************************************************************

  Function : DDV_StartMainWin_Slave()
  
  Purpose : start DDV main window as a slave module: init 'dmp' data structure,
  			set 'DDV_MAIN_DATA' property, start a dlg (the parent window of DDV
			panel).
  
  Parameters : vgp; sets defaults for non-autonomous viewers
							 
  Return value : handle to the DDV main window

*******************************************************************************/
extern WindoW DDV_StartMainWin_Slave(SAM_ViewGlobal *vgp)
{
Int2	Margins;
DdvMainWinPtr mWin_d;
WindoW	w; 
Boolean bRet;
	
	/*init data blocks*/
	mWin_d=(DdvMainWinPtr)MemNew(sizeof(DdvMainWin));
	if (!mWin_d) return(NULL);
	
	/*main window*/
	Margins=10*stdCharWidth;
    if(vgp == NULL)
        w=DocumentWindow(Margins,Margins ,-10, -10,
           szAppName, 
           NULL,
           DDV_WinMainResize);
    else
        w=DocumentWindow(vgp->Rect.left,vgp->Rect.top ,-10, -10,
            szAppName, 
            NULL,
            DDV_WinMainResize);

	if (!w) return(NULL);

	SetObjectExtra (w, (Pointer) mWin_d, (FreeProc)DDV_WinMainCleanup);
	mWin_d->hWndMain=w;
	/*this is a slave  editor/viewer*/
	mWin_d->AutonomeViewer=FALSE;
	/*mWin_d->bEditor=bEditor;*//*editor or viewer ?*//*todo*/

	/*init menu*/
	DDV_SetupMenus(w,FALSE);
	/*build GUI*/
	bRet=DDV_CreateViewerPanel(w,mWin_d, vgp);

	if (!bRet) return(NULL);

	mWin_d->Show_logo=FALSE;

	return(w);
}


NLM_EXTERN Int4 DDV_Accession2Gi (CharPtr string, DocType type)
{
   CharPtr str;
   LinkSetPtr lsp;
   Int4 gi;

   if (!EntrezIsInited ()) {
       Message (MSG_OK, "Network connection to Entrez unavailable");
       return 0;
   }

   str = MemNew (StringLen (string) + 10);
   sprintf (str, "\"%s\" [ACCN]", string);
   lsp = EntrezTLEvalString (str, type, -1, NULL, NULL);
   MemFree (str);
   if (lsp == NULL) return 0;
   if (lsp->num <= 0) {
       LinkSetFree (lsp);
       return 0;
   }
   gi = lsp->uids [0];
   LinkSetFree (lsp);
   return gi;
}

static Boolean DDV_DigitsOnly(Char* str)
{
    while(*str) {
        if (!isdigit(*str)) 
            return FALSE;
        str++;
    }
    return TRUE;
}

static void DDV_ImportEnableProc(TexT t)
{
    Char str[32];
    WindoW hImportDlg;
    DDV_ImportDialog *idp;

	hImportDlg= (WindoW)ParentWindow(t);
	if(hImportDlg == NULL) return;	
	idp = (DDV_ImportDialog *) GetObjectExtra(hImportDlg);
	if(idp == NULL) return;

    GetTitle(idp->DDV_tAccession, str, sizeof(str));
    if (StringLen(str) == 0) {
        Disable(idp->DDV_bImportAccept);
    } else {
        Enable(idp->DDV_bImportAccept);
        if(DDV_DigitsOnly(str)) SetValue(idp->DDV_gAccType, 2);
        else SetValue(idp->DDV_gAccType, 1);
    }
    return;
}

static void DDV_ImportAcceptProc(ButtoN b)
{
    Char str[32];
    Char out[32];
    WindoW hImportDlg;
    DDV_ImportDialog *idp;
    Bioseq *bsp1, *bsp2;
    SeqId si;
    Uint1 align_type;
    SeqAlign *sap, *sap_tmp;
    Uint2 entityID, itemID;

	hImportDlg= (WindoW)ParentWindow(b);
	if(hImportDlg == NULL) goto out;	
	idp = (DDV_ImportDialog *) GetObjectExtra(hImportDlg);
	if(idp == NULL) goto out;

    WatchCursor();
    GetTitle(idp->DDV_tAccession, str, sizeof(str));

    si.choice = SEQID_GI;
    if(GetValue(idp->DDV_gAccType) == 1) si.data.intvalue = 
        DDV_Accession2Gi(str, idp->AAorNN);
    else si.data.intvalue = atoi(str);

    if (si.data.intvalue <= 0) goto out;

    bsp1 = BioseqLockById(idp->sip);
    bsp2 = BioseqLockById(&si);

    if(idp->AAorNN == TYP_AA) sap = DDV_Blast2Seqs(bsp1, bsp2, TRUE,
                         "blastp");
    else sap = DDV_Blast2Seqs(bsp1, bsp2, TRUE, "blastn");

    if(sap == NULL) goto out;

    for (sap_tmp = (SeqAlignPtr)idp->sap->segs; sap_tmp->next != NULL;
    sap_tmp = sap_tmp->next);
    
    sap_tmp->next = sap;
    sap->next = NULL;

    AlnMgrReIndexSeqAlign(idp->sap);
    entityID=ObjMgrGetEntityIDForPointer((Pointer)idp->sap);
    itemID=GetItemIDGivenPointer (entityID, 
        OBJ_SEQALIGN, (Pointer) idp->sap);	
    
    ObjMgrSendProcMsg(OM_MSG_UPDATE, entityID, itemID,
        OBJ_SEQALIGN,0,0,idp->userdata);
    MemFree(idp->userdata);
    
out:
    ArrowCursor();

    Remove(idp->DDV_wImport);

    return;
}

static void DDV_ImportCancelProc(ButtoN b)
{
    WindoW hImportDlg;
    DDV_ImportDialog *idp;

	hImportDlg= (WindoW)ParentWindow(b);
	if(hImportDlg == NULL) return;	
	idp = (DDV_ImportDialog *) GetObjectExtra(hImportDlg);
	if(idp == NULL) return;
    MemFree(idp->userdata);

    Remove(idp->DDV_wImport);
    return;
}


NLM_EXTERN void DDV_ImportBioseqDlg(DDV_ImportDialog *idp)
{
    GrouP g, hg;
    ButtoN b;

    if(idp == NULL) return;
    if (!EntrezIsInited ()) {
        Message (MSG_OK, "Network connection to Entrez unavailable");
        return;
    }

    if(idp->sap->saip == NULL) {
        Message (MSG_OK, "The seqalign must be indexed");
        return;
    }
    
    idp->DDV_wImport =
        FixedWindow(-30, -20, -10, -10, "Import sequence using BLAST",
                    NULL);
    Nlm_SetObjectExtra (idp->DDV_wImport, (Nlm_VoidPtr)idp, StdCleanupExtraProc);
    hg = HiddenGroup(idp->DDV_wImport, 3, 0, NULL);
    SetGroupSpacing(hg, 30, 30);
    g = NormalGroup(hg, 1, 0, " Enter sequence identifier:", systemFont, NULL);
    SetGroupMargins(g, 10, 15);
    idp->DDV_tAccession = DialogText(g, "", 10, (TxtActnProc) DDV_ImportEnableProc);
    idp->DDV_gAccType =
        NormalGroup(hg, 1, 2, " accession type", systemFont, NULL);
    SetGroupMargins(idp->DDV_gAccType, 10, 10);
    RadioButton(idp->DDV_gAccType, "Accession");
    RadioButton(idp->DDV_gAccType, "gi");

    g = HiddenGroup(hg, 1, 2, NULL);
    SetGroupSpacing(g, 15, 15);
    idp->DDV_bImportAccept =
        DefaultButton(g, "OK", (BtnActnProc) DDV_ImportAcceptProc);
    b = PushButton(g, "Cancel", (BtnActnProc) DDV_ImportCancelProc);


    SetValue(idp->DDV_gAccType, 1);
    Disable(idp->DDV_bImportAccept);
    Select(idp->DDV_tAccession);
    Show(idp->DDV_wImport);
    return;
}

NLM_EXTERN void DDV_ImportBioseq(IteM i)
{
    DDV_ImportDialog *idp;
    DdvMainPtr    dmp;  /*DDV panel data*/
    DdvMainWinPtr mWin_d;
    WindoW		  hWinMain;
    DDVUpdateMSG  *dump;

    idp = MemNew(sizeof(DDV_ImportDialog));
    if(idp == NULL) return;

	hWinMain=(WindoW)ParentWindow(i);
	if (hWinMain==NULL) return;
	mWin_d = (DdvMainWinPtr) GetObjectExtra (hWinMain);
	if (mWin_d==NULL) return;
	if (mWin_d->hWndDDV==NULL) return;
	dmp = (DdvMainPtr) GetObjectExtra(mWin_d->hWndDDV);
	if (dmp==NULL) return;

    if(dmp->MSA_d.pgp_l.sap == NULL) return;

    idp->sap = dmp->MSA_d.pgp_l.sap;

    idp->sip = AlnMgrFindMaster(idp->sap);
    if(idp->sip == NULL) return;

    idp->AAorNN = TYP_AA;
    idp->DDV_wImport = NULL;
    idp->DDV_gAccType = NULL;
    idp->DDV_bImportAccept = NULL;
    idp->DDV_tAccession = NULL;
    dump = MemNew(sizeof(DDVUpdateMSG));
    if(dump == NULL) {
        MemFree(idp);
        return;
    }
    idp->userdata = dump;
    dump->data = NULL;
	dump->type=UPDATE_TYPE_EDIT_DELBSP;
    DDV_ImportBioseqDlg(idp);
}
/*
blastp align_type = 2
blastn align_type = 1
*/

NLM_EXTERN SeqAlign *DDV_Blast2Seqs(Bioseq *bsp1, Bioseq *bsp2, Boolean gapped,
                         Char *progname)
{
    BLAST_OptionsBlkPtr options;
    SeqAlign  *salp;
    
    if(bsp1 == NULL || bsp2 == NULL || progname == NULL) return NULL;
    options = BLASTOptionNew(progname, gapped);
    if(options == NULL) return NULL;
    options->discontinuous = FALSE;
    salp = BlastTwoSequences(bsp1, bsp2, progname, options);
    BLASTOptionDelete(options);
    return salp;
}

/*******************************************************************************

  Function : DDV_ImportSeqAlign(), DDV_ImportNucSeqAlign() and DDV_ImportProtSeqAlign()
  
  Purpose : import a sequence alignment of text type
  
  Return value : none 

*******************************************************************************/
static void DDV_ImportSeqAlign(Boolean is_prot,IteM i)
{
WindoW        hWinMain;
DdvMainPtr    dmp;
DdvMainWinPtr mWin_d;
SeqEntryPtr   sep;
ValNodePtr    vnp_ali;
SeqAlignPtr   sap;
Uint2 		  The_entityID,entityID,itemID;

	hWinMain=(WindoW)ParentWindow(i);
	
	mWin_d = (DdvMainWinPtr) GetObjectExtra (hWinMain);

	sep=ReadAnyAlignment(is_prot,NULL);

	if (sep!=NULL && mWin_d!=NULL){
		/*try to connect Entrez; if failure, still ok : the user can
		use DDV with a local file. Otherwise he/she will be in trouble ;-)*/
		if (mWin_d->AutonomeViewer && mWin_d->NetStartProc){
			mWin_d->NetStartProc(mWin_d->UseNetwork);
		}
		/*delete old data, if needed*/
		dmp=(DdvMainPtr)GetObjectExtra(mWin_d->hWndDDV);
		DDV_CleanupDDVPdata_g(dmp);
		DDV_CloseData(mWin_d,FALSE);

		vnp_ali=DDV_GetAndCheckSeqAlign(NULL,0,sep,mWin_d->fetchSepProc,&(mWin_d->dod),
			&The_entityID);
		if (vnp_ali && The_entityID!=0){

			/*get the first SeqAlign in the list and load the viewer panel*/
			mWin_d->vnp_ali=vnp_ali;

			if(mWin_d->Show_logo){
				mWin_d->Show_logo=FALSE;
			}
			sap=(SeqAlignPtr)vnp_ali->data.ptrvalue;
			entityID=sap->idx.entityID;
			itemID=sap->idx.itemID;
			GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, itemID, 
				OBJ_SEQALIGN, 0, 0, OBJ_SEQALIGN, 0);
		}
	}
	else{
		ErrPostEx (SEV_ERROR, 0, 0, "We do not support this format yet");
	}
}

extern void DDV_ImportNucSeqAlign(IteM i)
{
	DDV_ImportSeqAlign(FALSE,i);
}

extern void DDV_ImportProtSeqAlign(IteM i)
{
	DDV_ImportSeqAlign(TRUE,i);
}
