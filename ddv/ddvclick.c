/*  $Id: ddvclick.c,v 1.15 2000/01/12 15:49:42 durand Exp $
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
* File Name:  ddvclick.c
*
* Author:  Patrick Durand
*
* Version Creation Date:   09/20/99
*
* $Revision: 1.15 $
*
* File Description: mouse management code for DeuxD-Viewer (DDV)
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ddvclick.c,v $
* Revision 1.15  2000/01/12 15:49:42  durand
* add About dlg box and fix a bug in selection
*
* Revision 1.14  2000/01/10 15:09:45  durand
* Use Entrez instead of ID1
*
* Revision 1.13  2000/01/05 21:11:13  durand
* update mouse click actions and DrawSequence function for a better use from ddv and cn3d
*
* Revision 1.12  1999/12/07 21:40:13  durand
* add mouse modes menu and caret facility for the editor
*
* Revision 1.11  1999/12/06 16:19:19  durand
* add GoTo facility to DDV
*
* Revision 1.10  1999/12/03 23:17:22  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 1.9  1999/11/29 15:26:25  durand
* designed a new GUI to fix problems under MacOS, Linux and SGI
*
* Revision 1.8  1999/11/18 00:21:42  lewisg
* draw speedups and selection on mouseup
*
* Revision 1.7  1999/11/17 22:44:01  durand
* speed up the selection manager for large SeqAlign
*
* Revision 1.6  1999/11/03 21:29:48  durand
* add CTRL and SHFT keys for mouse selection. redesign the loader functions of DDV to properly register the message callbacks
*
* Revision 1.5  1999/10/29 14:15:39  durand
* add simple mouse selection functions
*
* Revision 1.4  1999/10/22 14:19:05  durand
* add new code for mouse selection
*
* Revision 1.3  1999/10/21 13:12:43  durand
* update selection system
*
* Revision 1.2  1999/10/12 15:01:28  lewisg
* resolve confict with internal/ddv
*
* Revision 1.1  1999/09/30 14:10:24  durand
* add ddv to toolkit
*
* Revision 1.2  1999/09/22 20:41:17  durand
* add mouse-click and drag functions
*
* Revision 1.1  1999/09/21 14:19:59  durand
* add basic mouse management functions
*
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <vibrant.h>

#include <udviewer.h>
#include <ddvpanel.h>
#include <ddvclick.h>
#include <ddvgraph.h>
#include <alignmgr.h>
#include <sqnutils.h>
#include <seqmgr.h>
#include <pgppop.h>

/*static MonitorPtr  mon=NULL;*/

/*****************************************************************************

Function: DDV_GetRowGivenMousePos()

Purpose: given some graphical info, this function returns line_num, the row
		number (i.e. sequence).
					
Return value: line_num (-1 if not found ; otherwise a one-based value)

*****************************************************************************/
static Int4 DDV_GetRowGivenMousePos(DdvMainPtr dmp, RecT * rcP, PoinT * pt)
{
Int4 Line_num=(Int4)-1,mouse_row,i,n;
ParaGPtr pgp;

	/*size of the Panel; see also DDV_DrawPanelContent_H()*/
	InsetRect(rcP,4,4);
	DDV_AdjustDrawingRect(rcP,&(dmp->GrData.udv_font));
	rcP->left+=dmp->GrData.udv_panel.cxName+dmp->GrData.udv_scale.cxLeftScale;
	rcP->top+=3*dmp->GrData.udv_panel.cyScale/2;

	/*mouse position*/
	mouse_row=(pt->y-rcP->top)/dmp->GrData.udv_font.LineHeight+
		dmp->GrData.udv_vscrl.ScrollPos;

	/*use the first ParaG of each sequence and mouse_row to find 
	row (==sequence) number*/
	for(i=0;i<dmp->MSA_d.pgp_l.nBsp;i++){
		pgp=(ParaGPtr)(dmp->MSA_d.pgp_l.TableHead[i]->data.ptrvalue);
		if (pgp->ScaleStyle==SCALE_POS_NONE) n=0;
		else n=1;
		/*within each ParaG, the seq. in on line 0 or 1 (if ruler)*/
		if (mouse_row>=pgp->StartLine && mouse_row<=(pgp->StartLine+n)){
			Line_num=mouse_row+1;
			break;
		}
	}

	return(Line_num);
}

/*****************************************************************************

Function: DDV_GetClosetSeqLocGivenBspPos()

Purpose: given a position in the Bioseq, this function return the closest
         selected region on that bioseq.
		
Parameters:	sip, eID, iID; bsioseq identifiers
            bsp_pos; current position on that bioseq
			old_pos; returns old_pos of a selection (used only if bModify is TRUE)
			bModify; if TRUE, create a modified slp

Note : sip can be null ONLY if bModify if FALSE.

Return value: see Purpose.

*****************************************************************************/
static SeqLocPtr DDV_GetClosetSeqLocGivenBspPos(SeqIdPtr sip, Uint2 eID, 
		Uint2 iID, Int4 bsp_pos, Int4Ptr old_pos, Boolean bModify)
{
SelStructPtr ssp;
SeqLocPtr    slp=NULL;
ValNodePtr   vnp_bsp,vnp;
Int4         diff,old_diff,diff_l,diff_r,bsp_start,bsp_stop;
Uint1        strand;

	ssp=ObjMgrGetSelected();
	if (ssp==NULL) return(NULL);
	
	vnp_bsp=DDV_GetSelectedRegions(ssp,eID,iID);
	if (vnp_bsp==NULL) return(NULL);
	vnp=vnp_bsp;
	*old_pos=(Int4)-1;
	old_diff=(Int4)INT4_MAX;

	while(vnp){
		bsp_start=SeqLocStart((SeqLocPtr)vnp->data.ptrvalue);
		bsp_stop=SeqLocStop((SeqLocPtr)vnp->data.ptrvalue);
		strand=SeqLocStrand((SeqLocPtr)vnp->data.ptrvalue);
		if (bsp_pos<bsp_start){/*left of a selected region ?*/
			diff=bsp_start-bsp_pos;
			if (diff<old_diff){	
				if (slp) slp=SeqLocFree(slp);
				if (bModify){
					slp = SeqLocIntNew (bsp_pos, bsp_stop, strand/*Seq_strand_minus*/, sip);
					*old_pos=bsp_stop;
				}
				else{
					slp = SeqLocIntNew (bsp_start, bsp_stop, strand, sip);
					*old_pos=(Int4)-1;
				}
				old_diff=diff;
			}
		}
		else if (bsp_pos>bsp_stop){/*right of a selected region ?*/
			diff=bsp_pos-bsp_stop;
			if (diff<old_diff){	
				if (slp) slp=SeqLocFree(slp);
				if (bModify){
					slp = SeqLocIntNew (bsp_start, bsp_pos, strand/*Seq_strand_plus*/, sip);
					*old_pos=bsp_start;
				}
				else{
					slp = SeqLocIntNew (bsp_start, bsp_stop, strand, sip);
					*old_pos=(Int4)-1;
				}
				old_diff=diff;
			}
		}
		else{/*inside a selected region ?*/
			if(bModify){
				diff_l=bsp_pos-bsp_start;
				diff_r=bsp_stop-bsp_pos;
				if (diff_l<diff_r){
					bsp_start=bsp_pos;
					*old_pos=bsp_stop;
					/*strand=Seq_strand_plus;*/
				}
				else{
					bsp_stop=bsp_pos;
					*old_pos=bsp_start;
					/*strand=Seq_strand_minus;*/
				}
				slp = SeqLocIntNew (bsp_start, bsp_stop, strand, sip);
				break;
			}
		}
		vnp=vnp->next;
	}
	ValNodeFree(vnp_bsp);
	return(slp);
}

/*****************************************************************************

Function: DDV_IsFullBspAlreadySel()

Purpose: check if a row is already selected.
		
Parameters:	bsp_eID,bsp_iID; bsp identifiers
			bsp_start,bsp_stop; full sequence 
			
Return value: TRUE if the sequence is already fully selected

*****************************************************************************/
static Boolean DDV_IsFullBspAlreadySel(Uint2 bsp_eID,Uint2 bsp_iID,Int4 bsp_start,
		Int4 bsp_stop)
{
SelStructPtr ssp;
ValNodePtr   vnp_bsp,vnp;
Int4         sel_bsp_start,sel_bsp_stop;

	ssp=ObjMgrGetSelected();
	if (ssp==NULL) return(FALSE);
	
	vnp_bsp=DDV_GetSelectedRegions(ssp,bsp_eID,bsp_iID);
	if (vnp_bsp==NULL) return(FALSE);
	vnp=vnp_bsp;

	while(vnp){
		sel_bsp_start=SeqLocStart((SeqLocPtr)vnp->data.ptrvalue);
		sel_bsp_stop=SeqLocStop((SeqLocPtr)vnp->data.ptrvalue);
		if (bsp_start==sel_bsp_start && bsp_stop==sel_bsp_stop)
			return(TRUE);
		vnp=vnp->next;
	}
	return(FALSE);
}

/*****************************************************************************

Function: DDV_GetCoordsGivenAClick()

Purpose: given some graphical info, this function returns bsp_coord, SeqAlign_coord,
         Disp_coord. Future implementation : will also return eID, iID, idx of
		 a feature (if the nice user clicks on a feature).
		
Parameters:	dmp; panel main data block
            rcP; panel size
			pt; mouse position
			bsp_coord, SeqAlign_coord, Disp_coord,Line_num; return values
            uWhere; where the mouse is in a ParaG
			
Note : all but Line_num are zero-based values. Line_num is a one-based value.
			
Return value: TRUE if the user clicks within a ParaG (see also 
		bsp_coord, SeqAlign_coord, Disp_coord, uWhere). NOTE : if the function 
		returns FALSE then bsp_coord, SeqAlign_coord, Line_num, Disp_coord and 
		uWhere are undefined !

*****************************************************************************/
static Boolean DDV_GetCoordsGivenAClick(DdvMainPtr dmp, RecT * rcP, PoinT * pt,
		Int4Ptr bsp_coord, Int4Ptr SeqAlign_coord, Int4Ptr Disp_coord,
		Int4Ptr Line_num, Uint1Ptr uWhere, ParaGPtr PNTR cur_pgp)
{
ValNodePtr     vnp;
ParaGPtr       pgp,TheParaG;
/*MsaTxtDispPtr  mtdp;*/
DDVRulerDescrPtr drdp;
Int4           mouse_row,mouse_col,i,diff;
Boolean        bFound;

	/*size of the Panel; see also DDV_DrawPanelContent_H()*/
	InsetRect(rcP,4,4);
	DDV_AdjustDrawingRect(rcP,&(dmp->GrData.udv_font));
	rcP->left+=dmp->GrData.udv_panel.cxName+dmp->GrData.udv_scale.cxLeftScale;
	rcP->top+=3*dmp->GrData.udv_panel.cyScale/2;

	/*mouse position*/
	mouse_row=(pt->y-rcP->top)/dmp->GrData.udv_font.LineHeight+
		dmp->GrData.udv_vscrl.ScrollPos;
	mouse_col=(pt->x-rcP->left)/dmp->GrData.udv_font.ColWidth+
		dmp->GrData.udv_hscrl.ScrollPos;

	bFound=FALSE;

	for(i=0;i<dmp->MSA_d.pgp_l.nBsp;i++){
		/*use the first ParaG of each sequence and mouse_row to find 
		row (==sequence) number; then go though the ParaG list of that
		sequence to find The ParaG, using mouse_col*/
		pgp=(ParaGPtr)(dmp->MSA_d.pgp_l.TableHead[i]->data.ptrvalue);
		if (mouse_row>=pgp->StartLine && mouse_row<=(pgp->StartLine+pgp->nLines-1)){
			vnp=dmp->MSA_d.pgp_l.TableHead[i];
			while(vnp){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				if (mouse_col>=pgp->StartLetter && mouse_col<=pgp->StopLetter){
					bFound=TRUE;
					TheParaG=pgp;
					break;
				}
				vnp=vnp->next;
			}
		}
	}
	/*ok : the mouse is inside a ParaG*/
	if (bFound){
		*cur_pgp=TheParaG;
		/*where is the mouse in the ParaG ?*/
		*uWhere=PGP_REGION_NOWHERE;
		if (TheParaG->ScaleStyle==SCALE_POS_NONE)
			diff=0;
		else
			diff=1;
		if (mouse_row==(TheParaG->StartLine+diff)){
			*uWhere=PGP_REGION_SEQ;
		}
		else if (mouse_row>(TheParaG->StartLine+diff) && 
				mouse_row<=(pgp->StartLine+pgp->nLines-1)){
			*uWhere=PGP_REGION_FEAT;
		}

		/*a. get Disp Coord*/
		*Line_num=mouse_row+1;/*base 1 value*/
		*Disp_coord=mouse_col;

		/*on the Sequence ? Then, get the coords. (BSP, DISP, SEQALIGN)*/
		/*PGP_REGION_FEAT is not yet implemented...  ;-) */
		if (*uWhere==PGP_REGION_SEQ){
			/*b. get BSP coord*/
			*bsp_coord=DDV_GetBspCoordGivenDispCoord(TheParaG,mouse_col);
			
			/*c. get SeqAlign Coord*/
			vnp=dmp->MSA_d.pgp_l.RulerDescr;
			while(vnp){
				drdp=(DDVRulerDescrPtr)vnp->data.ptrvalue;
				if (!drdp->bUnAligned){
					diff=mouse_col-drdp->disp_start;
					if (mouse_col>=drdp->disp_start && mouse_col<=drdp->disp_stop){
						*SeqAlign_coord=drdp->align_start+diff;
						break;
					}
				}
				else{
					*SeqAlign_coord=(Int4)-1;
				}
				vnp=vnp->next;
			}
		}
	}

	return(bFound);
}

/*****************************************************************************

Function: DDV_AutoHScroll()

Purpose: used from drag/hold mouse to auto-scroll the panel horizontally when
         the mouse reaches the left/right borders of the panel


Return value: none

*****************************************************************************/
static void DDV_AutoHScroll(DdvMainPtr dmp,Int4 Disp_coord)
{
Int4 from_col,to_col,from_row,to_row,old_pos;
BaR hsb;

	DDV_GetCurrentDispRange(dmp->hWndDDV,&(dmp->GrData),dmp->MSA_d.pgp_l.LengthAli,
			&from_col,&to_col,&from_row,&to_row);
	
	hsb = GetSlateHScrollBar ((SlatE) dmp->hWndDDV);
	old_pos=GetBarValue(hsb);
	if (Disp_coord>=to_col){/*scroll to the right*/
		SetValue(hsb,old_pos+2);
	}
	else if (Disp_coord<=from_col){
		SetValue(hsb,old_pos-2);
	}

}

/*****************************************************************************

Function: DDV_DispSelRangeInStatus()

Purpose: display position while selection is in progress. (mouse mode : selection)

Return value: none

*****************************************************************************/
static void DDV_DispSelRangeInStatus(PrompT InfoPanel,Int4 bsp_start, Int4 bsp_stop,
		Int4 line_num,CharPtr szSeqName)
{
Char szText[250];

	sprintf(szText,"Selection : %s, [%i]:[%i..%i]",szSeqName,line_num,
		bsp_start,bsp_stop);
	SetTitle(InfoPanel,szText);
}

/*****************************************************************************

Function: DDV_DispPositionInStatus()

Purpose: display position (mouse mode : query).

Return value: none

*****************************************************************************/
NLM_EXTERN void DDV_DispPositionInStatus(PrompT InfoPanel,Int4 bsp_pos,
		Int4 line_num,CharPtr szSeqName)
{
Char szText[250];

	sprintf(szText,"Position : %s, [%i]:[%i]",szSeqName,line_num,bsp_pos);
	SetTitle(InfoPanel,szText);
}

/*****************************************************************************

Function: DDV_GetSeqNameGivenRow()

Purpose: get a sequence number, given a list of ParaG and a row number

Note: row_num is a one-based value

Return value: -

*****************************************************************************/
NLM_EXTERN void DDV_GetSeqNameGivenRow(ValNodePtr * ParaG_List, Int4 row_num,
	CharPtr szSeqName)
{
SeqIdPtr sip;
ParaGPtr pgp;

	pgp=(ParaGPtr)(ParaG_List[row_num-1]->data.ptrvalue);
	sip = SeqIdFindBest(pgp->sip, SEQID_GENBANK);
	if (!sip)
		sip = SeqIdFindBest(pgp->sip, 0);
	SeqIdWrite(sip, szSeqName,PRINTID_TEXTID_ACCESSION, 20);   
}

/*****************************************************************************

Function: DDV_SendBSPSelectMsg()

Purpose: send a message to select a range of letters on a bioseq.

Return value: none

*****************************************************************************/
static void DDV_SendBSPSelectMsg(DdvMainPtr dmp,Int4 bsp_coord,Int4 SeqAlign_coord,
		Int4 Disp_coord,Int4 Line_num,ParaGPtr cur_pgp,PoinT * pt)
{
Int4        first_bsp_coord, bsp_coord_old,bsp_pos;
Uint1       direction;/*use to tell ObjMgr the direction of the mouse (left,right)*/
SeqLocPtr   slp;/*to send an AlsoSelectMsg*/
Uint2       bsp_eID,bsp_iID;
Boolean     bDeselectAll=TRUE;
WindoW      hParent;
DdvMainWinPtr dmwp;

	/*analyse the bsp coords*/
	if (bsp_coord==(Int4)-1){/*this is a gap*/
		dmp->ms.oldPos=dmp->ms.newPos=*pt;
		dmp->ms.old_row=Line_num;
		dmp->ms.old_col=Disp_coord;
		dmp->ms.old_pgp=cur_pgp;
		return;
	}

	/*get the very first click (from DDV_ClickProc)*/
	first_bsp_coord=
		DDV_GetBspCoordGivenDispCoord(dmp->ms.first_pgp,
			dmp->ms.first_col);
	if (first_bsp_coord==(Int4)-1){/*the first time the user 
		clicked, it was on a gap*/
		/*if we are in this scope, that means we have a bsp coord*/
		first_bsp_coord=bsp_coord;
		dmp->ms.first_col=Disp_coord;
		dmp->ms.first_pgp=cur_pgp;
	}

	bsp_pos=bsp_coord;
	
	if (ctrlKey || shftKey)
		bDeselectAll=FALSE;

	/*the mouse is moving to the left or to the right ?*/
	/*don't forget to check the strand...  ;-) */
	bsp_coord_old=DDV_GetBspCoordGivenDispCoord(dmp->ms.old_pgp,
					dmp->ms.old_col);
	if (AlnMgrGetNthStrand(dmp->MSA_d.pgp_l.sap, Line_num)==
			Seq_strand_plus){
		if (dmp->ms.oldPos.x<pt->x)
			direction=Seq_strand_plus;
		else
			direction=Seq_strand_minus;
	}
	else{
		if (dmp->ms.oldPos.x>pt->x)
			direction=Seq_strand_minus;
		else
			direction=Seq_strand_plus;
	}
	/*'from' (first_bsp_coord) always less than 'to' (bsp_coord)*/
	if (first_bsp_coord>bsp_coord)
		swap(&first_bsp_coord,&bsp_coord);
	/*save the new position*/
	dmp->ms.oldPos=dmp->ms.newPos=*pt;
	dmp->ms.old_row=Line_num;
	dmp->ms.old_col=Disp_coord;
	dmp->ms.old_pgp=cur_pgp;


	/*now, we can send the Select message*/
	/* first_pgp, old_pgp & cur_pgp have the same sip... 
	they are on a same row; I can use one of them*/
	slp = SeqLocIntNew (first_bsp_coord, bsp_coord, direction, cur_pgp->sip);
	UDV_DecodeIdxFeat (dmp->MSA_d.pgp_l.entitiesTbl[Line_num-1], 
		&bsp_eID,&bsp_iID);
    if (bsp_eID!=0 && bsp_iID!=0 && slp){
		/*update InfoPanel*/
		dmwp=(DdvMainWinPtr)GetObjectExtra(dmp->hParent);
		if (dmwp){
			Char szAccess[21];
			
			DDV_GetSeqNameGivenRow(dmp->MSA_d.pgp_l.TableHead, Line_num,
				szAccess);
	        
			switch(dmp->MouseMode){
				case DDV_MOUSEMODE_SELECT:
					DDV_DispSelRangeInStatus(dmwp->InfoPanel,first_bsp_coord+1, 
						bsp_coord+1,Line_num,szAccess);/*+1 : switch to one-based value*/
					break;
				case DDV_MOUSEMODE_QUERY:
				case DDV_MOUSEMODE_EDIT:
					DDV_DispPositionInStatus(dmwp->InfoPanel, 
						bsp_pos+1,Line_num,szAccess);/*+1 : switch to one-based value*/
					break;
			}
		}
		/*send the message*/
        if (dmp->MouseMode==DDV_MOUSEMODE_SELECT){
			ObjMgrAlsoSelect (bsp_eID, bsp_iID, 
    	        OBJ_BIOSEQ,OM_REGION_SEQLOC, slp);
		}
    }
	else{
		if (slp) 
			SeqLocFree(slp);
	}
}


/*****************************************************************************

Function: DDV_ClickProc()

Purpose: manage Mouse-Click; this function is designed for the autonomous viewer
		to manage the InfoPanel and the Features List Dialog Box
		DON'T USE FOR EXTERNAL SOFTWARE...
		
Parameters:	p; panel handle (currently unused)
			pt; new mouse position

Return value: none

*****************************************************************************/
NLM_EXTERN void DDV_ClickProc(PaneL p, PoinT pt)
{
DdvMainPtr 	dmp;
RecT        rcP;
Int4        bsp_coord, SeqAlign_coord, Disp_coord, Line_num,bsp_start,bsp_stop,
            old_pos;
Uint1       uWhere;
ParaGPtr    cur_pgp=NULL;
SeqLocPtr   slp;
Uint2       bsp_eID,bsp_iID;
Boolean     bDeselectAll=TRUE,bIsAlreadyFullSel;

	/*get the panel data*/
	dmp = (DdvMainPtr) GetObjectExtra(p);
	if (dmp==NULL) return;
	
	ObjectRect(p,&rcP);

	if (ctrlKey || shftKey)
		bDeselectAll=FALSE;

	/*the mouse is located in a ParaG ?*/
	/*if YES : the user is selecting either sequence letter(s) or a feature*/
	if (pt.x>(rcP.left+dmp->GrData.udv_panel.cxName+
		dmp->GrData.udv_scale.cxLeftScale)){
		if (pt.y>(rcP.top+dmp->GrData.udv_panel.cyScale)){
			/*get the ParaG list and look for the ParaG where the mouse is located*/
			/*note : only implemented for DDV_DISP_HORZ type of display...*/
			if (dmp->MSA_d.pgp_l.DisplayType==DDV_DISP_HORZ){
				if (DDV_GetCoordsGivenAClick(dmp,&rcP,&pt,
					&bsp_coord,&SeqAlign_coord,&Disp_coord,&Line_num,&uWhere,&cur_pgp)){
					if (uWhere==PGP_REGION_SEQ){
						if (bDeselectAll && dmp->MouseMode==DDV_MOUSEMODE_SELECT)
							ObjMgrDeSelectAll ();
						if (shftKey && bsp_coord!=(Int4)-1){
							UDV_DecodeIdxFeat (dmp->MSA_d.pgp_l.entitiesTbl[Line_num-1], 
								&bsp_eID,&bsp_iID);
							if (bsp_eID!=0 && bsp_iID!=0){
								slp=DDV_GetClosetSeqLocGivenBspPos(cur_pgp->sip,bsp_eID, 
									bsp_iID, bsp_coord, &old_pos, TRUE);
								if (slp){
									if (dmp->MouseMode==DDV_MOUSEMODE_SELECT)
										ObjMgrAlsoSelect (bsp_eID, bsp_iID, 
											OBJ_BIOSEQ,OM_REGION_SEQLOC, slp);
									if (old_pos!=(Int4)-1)
										Disp_coord=old_pos;
								}
							}
							return;
						}
						/*save the positions */
						dmp->ms.oldPos=pt;
						dmp->ms.newPos=pt;
						dmp->ms.Action_type=MS_ACTION_SELECT_SEQ;
						dmp->ms.old_row=Line_num;
						dmp->ms.old_col=Disp_coord;
						dmp->ms.old_pgp=cur_pgp;
						dmp->ms.first_row=Line_num;
						dmp->ms.first_col=Disp_coord;
						dmp->ms.first_pgp=cur_pgp;
						/*update caret pos, if needed only*/
						if (dmp->MouseMode==DDV_MOUSEMODE_EDIT){
							dmp->dci.old_row=dmp->dci.new_row;
							dmp->dci.old_col=dmp->dci.new_col;
						}
						/*select a single letter*/
						DDV_SendBSPSelectMsg( dmp, bsp_coord, SeqAlign_coord,
							Disp_coord, Line_num, cur_pgp,&pt);
					}
				}
			}
		}
	}	
	/*the mouse is located within the Name list region (left of ParaG) ?*/
	else if (pt.x>=dmp->GrData.udv_panel.cxName-1 && 
			pt.x<=dmp->GrData.udv_panel.cxName+3){
				rcP.left=5*dmp->GrData.udv_font.cxChar; /*min value*/
				rcP.right=2*PANEL_NAME_WIDTH*dmp->GrData.udv_font.cxChar;/*max val*/
				pt.x=dmp->GrData.udv_panel.cxName;
				dmp->ms.oldPos=pt;
				dmp->ms.newPos=pt;
				dmp->ms.rcClip=rcP;
				dmp->ms.Action_type=MS_ACTION_RESIZE_WIN;
				InvertMode();	
				UDV_draw_double_cursor(dmp->ms.rcClip,
					dmp->ms.oldPos);
				CrossCursor();
	}
	/*the mouse is located on a sequence name ?*/
	else if (pt.x<dmp->GrData.udv_panel.cxName){
		if (dblClick){/*load the sequence viewer on that sequence*/
			Line_num=DDV_GetRowGivenMousePos(dmp, &rcP, &pt);
			if (Line_num!=(Int4)-1){
				UDV_DecodeIdxFeat (dmp->MSA_d.pgp_l.entitiesTbl[Line_num-1], 
					&bsp_eID,&bsp_iID);
				GatherProcLaunch (OMPROC_VIEW, FALSE, bsp_eID, bsp_iID, 
						OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
			}
		}
	}
	/*the mouse is located  just on the left of a sequence ? 
		If true, select the full sequence*/
	else if ((pt.x>=dmp->GrData.udv_panel.cxName+3) /*after the name of the sequence ?*/
			&&
			(pt.x<=dmp->GrData.udv_panel.cxName+3+dmp->GrData.udv_scale.cxLeftScale+
			dmp->GrData.udv_font.cxChar) /*before the sequence ?*/
			&& 
			(dmp->MouseMode==DDV_MOUSEMODE_SELECT)/*valid only if selection mode*/
			){
		if (pt.y>(rcP.top+dmp->GrData.udv_panel.cyScale)){
			/*get the ParaG list and look for the ParaG where the mouse is located*/
			/*note : only implemented for DDV_DISP_HORZ type of display...*/
			if (dmp->MSA_d.pgp_l.DisplayType==DDV_DISP_HORZ){
				Line_num=DDV_GetRowGivenMousePos(dmp, &rcP, &pt);
#ifdef _DDVEDIT
                if (Line_num!=(Int4)-1)
					dmp->Row = Line_num-1;
#endif /* _DDVEDIT */
				if (Line_num!=(Int4)-1){
					/*get min,max coord. of BSP*/
					UDV_GetBspRangeinPGPList(dmp->MSA_d.pgp_l.TableHead[Line_num-1],
						&bsp_start,&bsp_stop);
					if (bsp_start!=INT4_MAX && bsp_stop!=INT4_MIN){
						/*now, we can send the Select message*/
						/* first_pgp, old_pgp & cur_pgp have the same sip... 
						they are on a same row; I can use one of them*/
						cur_pgp=(ParaGPtr)(dmp->MSA_d.pgp_l.TableHead[Line_num-1]->data.ptrvalue);
						slp = SeqLocIntNew (bsp_start, bsp_stop, Seq_strand_plus, cur_pgp->sip);
						UDV_DecodeIdxFeat (dmp->MSA_d.pgp_l.entitiesTbl[Line_num-1], 
							&bsp_eID,&bsp_iID);
						if (bsp_eID!=0 && bsp_iID!=0 && slp){
							bIsAlreadyFullSel=DDV_IsFullBspAlreadySel(bsp_eID,bsp_iID,
								bsp_start,bsp_stop);
							if (!bIsAlreadyFullSel){
								if (!bDeselectAll){
									ObjMgrAlsoSelect (bsp_eID, bsp_iID, 
										OBJ_BIOSEQ,OM_REGION_SEQLOC, slp);
								}
								else{
									ObjMgrDeSelectAll ();
									ObjMgrSelect (bsp_eID, bsp_iID, 
										OBJ_BIOSEQ,OM_REGION_SEQLOC, slp);
								}
							}
							else{
								if (bDeselectAll)
									ObjMgrDeSelectAll ();
								ObjMgrDeSelect (bsp_eID, bsp_iID, 
									OBJ_BIOSEQ,OM_REGION_SEQLOC, slp);
							}
						}		
					}
					dmp->ms.Action_type=MS_ACTION_SELECT_FULL_SEQ;
				}
			}
		}
	}
	else{/*click outside "special" regions : deselect everything*/
		if (dmp->MouseMode==DDV_MOUSEMODE_SELECT)
			ObjMgrDeSelectAll ();
	}
}


/*****************************************************************************

Function: DDV_DragProc()

Purpose: manage Mouse-Drag

Parameters:	p; panel handle (currently unused)
			pt; new position

Return value: none

*****************************************************************************/
NLM_EXTERN void DDV_DragProc(PaneL p, PoinT pt)
{
DdvMainPtr 	dmp;
RecT        rcP;
ParaGPtr    cur_pgp=NULL;
Int4        bsp_coord,SeqAlign_coord, 
		    Disp_coord, Line_num;/*to get the coordinates*/
Uint1       uWhere; /*where the use clicked (seq, feat, name,...)*/
	        

	/*get the panel data*/
	dmp = (DdvMainPtr) GetObjectExtra(p);
	if (dmp==NULL) return;
	
	/*what is happening ?*/
	switch(dmp->ms.Action_type){
		case MS_ACTION_SELECT_SEQ:/*selection on a sequence*/
			if (shftKey) break;
			ObjectRect(p,&rcP);
			/*signle sequence selection: force the mouse to stay on the same row*/
			pt.y=dmp->ms.oldPos.y;
			/*get the ParaG list and look for the ParaG where the mouse is located*/
			/*note : only implemented for DDV_DISP_HORZ type of display...*/
			if (dmp->MSA_d.pgp_l.DisplayType==DDV_DISP_HORZ){
				if (DDV_GetCoordsGivenAClick(dmp,&rcP,&pt,
					&bsp_coord,&SeqAlign_coord,&Disp_coord,&Line_num,
					&uWhere,&cur_pgp)){
					/*if same position... bye,bye */
					if (dmp->ms.old_col==Disp_coord)
						break;
					if (uWhere==PGP_REGION_SEQ){/*we are on the sequence itself*/
						DDV_SendBSPSelectMsg( dmp, bsp_coord, SeqAlign_coord,
							 Disp_coord, Line_num, cur_pgp,&pt);
					}
				}
			}
			break;
		case MS_ACTION_RESIZE_WIN:/*the user moves the 3D line located between 
			name area and SeqAlign area*/
			InvertMode();	
			UDV_draw_double_cursor(dmp->ms.rcClip,
				dmp->ms.oldPos);
			dmp->ms.newPos=pt;
			UDV_draw_double_cursor(dmp->ms.rcClip,
				dmp->ms.newPos);
			dmp->ms.oldPos=dmp->ms.newPos;
			Update();
			break;
		case MS_ACTION_SELECT_FULL_SEQ:
			break;
	}
}

/*****************************************************************************

Function: DDV_ReleaseProc()

Purpose: manage Mouse-Release

Parameters:	p; panel handle (currently unused)
			pt; new position

Return value: none

*****************************************************************************/
NLM_EXTERN void DDV_ReleaseProc(PaneL p, PoinT pt)
{
WindoW 	      temport;
DdvMainPtr 	  dmp;
DdvMainWinPtr dmwp;
DDVUpdateMSGPtr dump;
RecT        rcP;
ParaGPtr    cur_pgp=NULL;
Int4        bsp_coord, SeqAlign_coord, 
		    Disp_coord, Line_num;/*to get the coordinates*/
Uint1       uWhere; /*where the use clicked (seq, feat, name,...)*/
Uint2       bsp_eID, bsp_iID;


	/*get the panel data*/
	dmp = (DdvMainPtr) GetObjectExtra(p);
	if (dmp==NULL) return;


    switch(dmp->ms.Action_type){
		case MS_ACTION_SELECT_FULL_SEQ:
            ObjectRect(p,&rcP);
			Line_num=DDV_GetRowGivenMousePos(dmp,&rcP,&pt);
			if (Line_num!=(Int4)-1){
                UDV_DecodeIdxFeat (dmp->MSA_d.pgp_l.entitiesTbl[Line_num-1], 
                    &bsp_eID,&bsp_iID);                    
                ObjMgrSendMsg(OM_MSG_MOUSEUP, bsp_eID, bsp_iID, OBJ_BIOSEQ);
			}
			break;
        case MS_ACTION_SELECT_SEQ:
            ObjectRect(p,&rcP);
            pt.y=dmp->ms.oldPos.y;
            if (dmp->MSA_d.pgp_l.DisplayType==DDV_DISP_HORZ){            
                if (DDV_GetCoordsGivenAClick(dmp,&rcP,&pt,
                        &bsp_coord,&SeqAlign_coord,&Disp_coord,&Line_num,&uWhere,&cur_pgp)){
					if (uWhere==PGP_REGION_SEQ){
                        UDV_DecodeIdxFeat (dmp->MSA_d.pgp_l.entitiesTbl[Line_num-1], 
                            &bsp_eID,&bsp_iID);                    
                        ObjMgrSendMsg(OM_MSG_MOUSEUP, bsp_eID, bsp_iID, OBJ_BIOSEQ);
						/*send a OM_MSG_UPDATE message to move the caret (if needed)*/
						if (dmp->MouseMode==DDV_MOUSEMODE_EDIT){
							dump=(DDVUpdateMSGPtr)MemNew(sizeof(DDVUpdateMSG));
							if (dump){
								dmp->dci.new_row=Line_num;
								dmp->dci.new_col=Disp_coord;
								dump->type=UPDATE_TYPE_CARETPOS;
								dump->data=NULL;
								ObjMgrSendProcMsg(OM_MSG_UPDATE, bsp_eID, bsp_iID, OBJ_BIOSEQ,
										0,0,(Pointer)dump);
								MemFree(dump);
							}
						}
                    }
                }
            }
			break;
		case MS_ACTION_RESIZE_WIN:
			temport=SavePort((WindoW)p);
			Select(p);
			InvertMode();	
			UDV_draw_double_cursor(dmp->ms.rcClip,
				dmp->ms.oldPos);
			/*redraw panel with new 'cxName' value*/
			if (PtInRect(dmp->ms.newPos,&dmp->ms.rcClip)){
				RecT rc;

				ObjectRect(p,&rc);
				rc.left=0;/*for obscure reasons, not == 0*/
				rc.top=0;
				dmp->GrData.udv_panel.cxName=dmp->ms.newPos.x;
				DDV_SetupWin (p,FALSE,NULL);
				InvalRect(&rc);
			}
			RestorePort(temport);
			break;
	}	

	MemSet(&(dmp->ms),0,sizeof(UDV_mouse_select));
	dmp->ms.Action_type=MS_ACTION_FEAT_NOTHING;
	/*update InfoPanel*/
	dmwp=(DdvMainWinPtr)GetObjectExtra(dmp->hParent);
	if (dmwp && dmp->MouseMode!=DDV_MOUSEMODE_EDIT)
		SetTitle(dmwp->InfoPanel,"Ready !");
	ArrowCursor();
}

/*****************************************************************************

Function: DDV_HoldProc()

Purpose: manage Mouse-Hold

Parameters:	p; panel handle (currently unused)
			pt; new position

Return value: none

*****************************************************************************/
NLM_EXTERN void DDV_HoldProc(PaneL p, PoinT pt)
{
DdvMainPtr 	dmp;
RecT        rcP;
ParaGPtr    cur_pgp=NULL;
Int4        bsp_coord, SeqAlign_coord, 
		    Disp_coord, Line_num;/*to get the coordinates*/
Uint1       uWhere; /*where the use clicked (seq, feat, name,...)*/

	/*get the panel data*/
	dmp = (DdvMainPtr) GetObjectExtra(p);
	if (dmp==NULL) return;
	/*what is happening ?*/
	switch(dmp->ms.Action_type){
		case MS_ACTION_SELECT_SEQ:
			if (shftKey) break;
			ObjectRect(p,&rcP);
			/*signle sequence selection: force the mouse to stay on the same row*/
			pt.y=dmp->ms.oldPos.y;
			/*get the ParaG list and look for the ParaG where the mouse is located*/
			/*note : only implemented for DDV_DISP_HORZ type of display...*/
			if (dmp->MSA_d.pgp_l.DisplayType==DDV_DISP_HORZ){
				if (DDV_GetCoordsGivenAClick(dmp,&rcP,&pt,
					&bsp_coord,&SeqAlign_coord,&Disp_coord,&Line_num,
					&uWhere,&cur_pgp)){
					/*if the mouse is on the far left/right; auto-scroll the panel*/
					DDV_AutoHScroll(dmp,Disp_coord);
					if (uWhere==PGP_REGION_SEQ){/*we are on the sequence itself*/
						/*if same position... bye,bye */
						if (dmp->ms.old_col==Disp_coord)
							break;
						DDV_SendBSPSelectMsg( dmp, bsp_coord, SeqAlign_coord,
							 Disp_coord, Line_num, cur_pgp,&pt);
					}
				}
			}
			break;
		case MS_ACTION_RESIZE_WIN:
		case MS_ACTION_SELECT_FULL_SEQ:
			break;
	}/*end switch()*/
}

/*****************************************************************************

Function: DDV_MoveCaretLR()

Purpose: move the caret to the right or to the left

Parameters:	

Return value: none

*****************************************************************************/
static void DDV_MoveCaretLR(DdvMainPtr dmp,Int4 decal_Hscroll,
	Boolean bMoveRight,BaR hsb,Int4 old_Hpos,Int4 from_col,Int4 to_col)
{
DDVUpdateMSGPtr dump;
Uint2           bsp_eID, bsp_iID;

	dump=(DDVUpdateMSGPtr)MemNew(sizeof(DDVUpdateMSG));
	if (dump){
		/*get bsp IDs*/
		UDV_DecodeIdxFeat(dmp->MSA_d.pgp_l.entitiesTbl[dmp->dci.new_row-1], 
            	&bsp_eID,&bsp_iID);
		/*save old position*/
		dmp->dci.old_col=dmp->dci.new_col;
		dmp->dci.old_row=dmp->dci.new_row;
		/*update the position*/
		if (bMoveRight)
			dmp->dci.new_col+=decal_Hscroll;
		else
			dmp->dci.new_col-=decal_Hscroll;
		if (dmp->dci.new_col<0){
			dmp->dci.new_col=0;
			Beep();
			return;
		}
		if (dmp->dci.new_col>dmp->MSA_d.pgp_l.LengthAli-1){
			dmp->dci.new_col=dmp->MSA_d.pgp_l.LengthAli-1;
			Beep();
			return;
		}
		/*auto-scroll the panel if needed*/
		if (dmp->dci.new_col<=from_col){
			SetValue(hsb,old_Hpos-(to_col-from_col)/2);
		}
		if (dmp->dci.new_col>=to_col){
			SetValue(hsb,old_Hpos+(to_col-from_col)/2);
		}
		/*if ok, move the caret*/
		dump->type=UPDATE_TYPE_CARETPOS;
		dump->data=NULL;
		ObjMgrSendProcMsg(OM_MSG_UPDATE, bsp_eID, bsp_iID, OBJ_BIOSEQ,
				0,0,(Pointer)dump);
		MemFree(dump);
	}
}

/*****************************************************************************

Function: DDV_MoveCaretLR()

Purpose: move the caret up or down

Parameters:	

Return value: none

*****************************************************************************/
static void DDV_MoveCaretUD(DdvMainPtr dmp,Int4 decal_Vscroll,
	Boolean bMoveDown,BaR vsb,Int4 old_Vpos,Int4 from_row,Int4 to_row)
{
DDVUpdateMSGPtr dump;
Uint2           bsp_eID, bsp_iID;

	dump=(DDVUpdateMSGPtr)MemNew(sizeof(DDVUpdateMSG));
	if (dump){
		/*get bsp IDs*/
		UDV_DecodeIdxFeat(dmp->MSA_d.pgp_l.entitiesTbl[dmp->dci.new_row-1], 
            	&bsp_eID,&bsp_iID);
		/*save old position*/
		dmp->dci.old_row=dmp->dci.new_row;
		dmp->dci.old_col=dmp->dci.new_col;
		/*update the position*/
		if (bMoveDown)
			dmp->dci.new_row+=decal_Vscroll;
		else
			dmp->dci.new_row-=decal_Vscroll;
		if (dmp->dci.new_row<1){
			dmp->dci.new_row=1;
			Beep();
			return;
		}
		if (dmp->dci.new_row>dmp->MSA_d.pgp_l.nBsp){
			dmp->dci.new_row=dmp->MSA_d.pgp_l.nBsp;
			Beep();
			return;
		}
		/*auto-scroll the panel if needed*/
		if (dmp->dci.new_row<=from_row){
			SetValue(vsb,old_Vpos-(to_row-from_row)/2);
		}
		if (dmp->dci.new_row>=to_row-1){
			SetValue(vsb,old_Vpos+(to_row-from_row)/2);
		}
		/*if ok, move the caret*/
		dump->type=UPDATE_TYPE_CARETPOS;
		dump->data=NULL;
		ObjMgrSendProcMsg(OM_MSG_UPDATE, bsp_eID, bsp_iID, OBJ_BIOSEQ,
				0,0,(Pointer)dump);
		MemFree(dump);
	}
}

/*****************************************************************************

Function: DDV_KeyboardProc()

Purpose: manage keyboard actions

Parameters:	s; panel handle (DDV panel)
			ch; key

Return value: none

*****************************************************************************/
NLM_EXTERN void DDV_KeyboardProc (SlatE s, Char ch)
{
DdvMainWinPtr dmwp; /*main window*/
DdvMainPtr dmp;     /*DDV panel associated data*/
BaR        hsb,vsb; /*DDV' scroll bars*/
Int4       old_Hpos,old_Vpos,/*scrolls positions*/
           decal_Hscroll,decal_Vscroll,/*used to scroll DDV's panel*/
           from_col,to_col,from_row,to_row,/*display size*/
           bsp_pos;
		   
	if ( (int) ch == 0 ) return;
	
	/*get the panel data*/
	dmp = (DdvMainPtr) GetObjectExtra(s);
	if (dmp==NULL) return;
	
	/*get scroll handles and data*/
	hsb = GetSlateHScrollBar ((SlatE) dmp->hWndDDV);
	old_Hpos=GetBarValue(hsb);
	vsb = GetSlateVScrollBar ((SlatE) dmp->hWndDDV);
	old_Vpos=GetBarValue(vsb);

	/*get the display size*/
	DDV_GetCurrentDispRange(dmp->hWndDDV,&(dmp->GrData),
		dmp->MSA_d.pgp_l.LengthAli,&from_col,&to_col,&from_row,&to_row);
	
	/*if pressed, scroll panel by one 'panel' size*/
	if (ctrlKey) {
		decal_Hscroll=to_col-from_col;
		decal_Vscroll=to_row-from_row;
	}
	else{
		decal_Hscroll=1;
		decal_Vscroll=1;
	}

	/*switch to one-bsed value because dmp->dci values 
	  are one-based values*/
	from_row++;to_row++;

	switch ((int) TO_UPPER(ch)){
		case NLM_LEFT:
			/*move the caret ?*/
			if (dmp->MouseMode==DDV_MOUSEMODE_EDIT){
				decal_Hscroll=1;/*doesn't accept ctrlKey motions*/
				DDV_MoveCaretLR(dmp,decal_Hscroll,FALSE,hsb,old_Hpos,
					from_col,to_col);
			}
			/*scroll the panel*/
			else{
				SetValue(hsb,old_Hpos-decal_Hscroll);
			}
			break;
		case NLM_RIGHT:
			/*move the caret ?*/
			if (dmp->MouseMode==DDV_MOUSEMODE_EDIT){
				decal_Hscroll=1;/*doesn't accept ctrlKey motions*/
				DDV_MoveCaretLR(dmp,decal_Hscroll,TRUE,hsb,old_Hpos,
					from_col,to_col);
			}
			/*scroll the panel*/
			else {
				SetValue(hsb,old_Hpos+decal_Hscroll);
			}
			break;
		case NLM_UP:
			/*move the caret ?*/
			if (dmp->MouseMode==DDV_MOUSEMODE_EDIT){
				decal_Vscroll=1;/*doesn't accept ctrlKey motions*/
				DDV_MoveCaretUD(dmp,decal_Vscroll,FALSE,vsb,old_Vpos,
					from_row,to_row);
			}
			/*scroll the panel*/
			else{
				SetValue(vsb,old_Vpos-decal_Vscroll);
			}
			break;
		case NLM_DOWN:
			/*move the caret ?*/
			if (dmp->MouseMode==DDV_MOUSEMODE_EDIT){
				decal_Vscroll=1;/*doesn't accept ctrlKey motions*/
				DDV_MoveCaretUD(dmp,decal_Vscroll,TRUE,vsb,old_Vpos,
					from_row,to_row);
			}
			/*scroll the panel*/
			else{
				SetValue(vsb,old_Vpos+decal_Vscroll);
			}
			break;
		default:
			Beep ();
			break;
	}

	/*update InfoPanel with position*/
	dmwp=(DdvMainWinPtr)GetObjectExtra(dmp->hParent);
	if (dmwp && dmp->MouseMode==DDV_MOUSEMODE_EDIT){
		Char szAccess[21];
		DDV_GetSeqNameGivenRow(dmp->MSA_d.pgp_l.TableHead, dmp->dci.new_row,
			szAccess);
		bsp_pos=DDV_GetBspCoordGivenPgpList(
			dmp->MSA_d.pgp_l.TableHead[dmp->dci.new_row-1],
			dmp->dci.new_col);
		if (bsp_pos!=(Int4)-1){
			DDV_DispPositionInStatus(dmwp->InfoPanel, 
				bsp_pos+1,dmp->dci.new_row,szAccess);
					/*+1 : switch to one-based value*/
		}
	}
}

