/*   udvgraph.c
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
* File Name:  udvgraph.c
*
* Author:  Patrick Durand
*
* Version Creation Date:   5/3/99
*
* $Revision: 6.14 $
*
* File Description: 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: udvgraph.c,v $
* Revision 6.14  1999/09/13 20:37:17  durand
* update UDV_Draw_scale for DDV
*
* Revision 6.13  1999/09/07 23:03:10  durand
* update UDV_Draw_scale to deal with discontinuous SEqAlign
*
* Revision 6.12  1999/07/30 20:08:57  durand
* updates for the new Entrez graphical viewer
*
* Revision 6.11  1999/07/19 20:35:35  durand
* switch ScalePositionfrom from Boolean to Uint1 in UDV_Draw_scale
*
* Revision 6.10  1999/06/29 14:59:52  durand
* update udv_draw_scale for DDV
*
* Revision 6.9  1999/06/16 13:07:00  durand
* update UDV functions to be used by DDV
*
* Revision 6.8  1999/06/08 13:52:35  durand
* update UDV data structures for the MSA editor
*
* Revision 6.7  1999/06/07 15:39:43  durand
* add LOG line to keep track of the history
*
*
*
* ==========================================================================
*/

#include <udviewer.h>
#include <salmedia.h>
#include <saledit.h>

static RecT  UDV_calc_RCdraw(Int4 StartLine,Int4 nLines, RecT rcP,
		Int2 decal_gauche,Int4 decal_haut,Int2 LineHeight);

/*****************************************************************************

Function: LogoFontCreate()

Purpose: create a simple set of fonts for the logo 

Return value: none

*****************************************************************************/
NLM_EXTERN void LogoFontCreate(FonT PNTR f1,FonT PNTR f2,FonT PNTR f3)
{

#ifdef WIN_MAC
		  *f1 = ParseFont ("Geneva,10");
		  *f2 = ParseFont ("Times,14,b");
		  *f3 = ParseFont ("Times,14,b,i");
#endif

#ifdef WIN_MSWIN
		  *f1 = ParseFont ("Arial,16");
		  *f2 = ParseFont ("Times New Roman,20,b");
		  *f3 = ParseFont ("Times New Roman,20,b,i");
#endif

#ifdef WIN_MOTIF
		  *f1 = ParseFont ("Helvetica,20,b");
		  *f2 = ParseFont ("Times,24,b");
		  *f3 = ParseFont ("Times,24,b,i");
#endif

}


/*****************************************************************************

Function: UDV_FontCreate()

Purpose: create a simple set of fonts for the viewer

Return value: handle to the font 

*****************************************************************************/
NLM_EXTERN FonT UDV_InitFont(UDVFontDataPtr udvfp)
{

#ifdef WIN_MAC
		udvfp->hFnt = ParseFont ("Monaco, 9");
#endif

#ifdef WIN_MSWIN
		udvfp->hFnt = ParseFont ("Courier, 7");
#endif

#ifdef WIN_MOTIF
		udvfp->hFnt = ParseFont ("fixed, 12");
#endif

	udvfp->UseDefaultColorLetter=TRUE;
	udvfp->LetterColor=(Uint4)FONT_DEF_COLOR;

	return(udvfp->hFnt);
}

/*****************************************************************************

Function: UDV_FontDim()

Purpose: computes cxChar and cyChar  

Return value: none

*****************************************************************************/
NLM_EXTERN void  UDV_FontDim(Int2Ptr cxChar,Int2Ptr cyChar)
{
	*cxChar=MaxCharWidth();
	*cyChar=LineHeight();
}


/*****************************************************************************

Function: UDV_ComputeLineHeight()

Purpose: compute the height of a single data line  

Return value: the height

*****************************************************************************/
NLM_EXTERN Int2 UDV_ComputeLineHeight(Int2 cyChar)
{
	/*return(3*cyChar/2);*/
	return(cyChar);
}

/*****************************************************************************

Function: UDV_ComputePanelSize()

Purpose: computes cxClient and cyClient  
         rc is the UDV panel
		 
Return value: none

*****************************************************************************/
NLM_EXTERN void  UDV_ComputePanelSize(RecT rc,Int2Ptr cxClient,
			Int2Ptr cyClient)
{
 
	InsetRect (&rc, (Int2)VIEWER_HORZ_MARGIN, (Int2)VIEWER_VERT_MARGIN);
	*cxClient=rc.right-rc.left+1;
	*cyClient=rc.bottom-rc.top+1;

}

/*****************************************************************************

Function: UDV_ComputeBlockByLine()

Purpose: computes the number of Ten-letter blocks as well as the number of 
		letters displayed on one line of the viewer  

Return value: none (results are in nCharByLine & nBlockByLine)

*****************************************************************************/
NLM_EXTERN void  UDV_ComputeBlockByLine(Int2 cxClient,Int2 cxName,
			Int2 cxLeftScale,Int2 cxChar,Int2Ptr nCharByLine,
			Int2Ptr nBlockByLine)
{
Int2 cx, cxWin,i,LB;

	LB=LETTER_BLOCK_WIDTH+1;/*+1: inter-block space*/
	cxWin=cxClient-cxName-cxLeftScale-2*VIEWER_HORZ_MARGIN;
	i=1;
	while(TRUE){/*count number of blocks within panel*/
		cx=i*LB*cxChar;
		if (cx>cxWin){
			i--;
			break;
		}
		i++;
	}
	
	*nCharByLine=i*LETTER_BLOCK_WIDTH;
	*nBlockByLine=i;
}

/*****************************************************************************

Function: UDV_Init_GraphData()

Purpose: Init Graphical values for the UnDViewer. Called one times at the
		start of the soft.

Return value: none (results are in GrDataPtr)

*****************************************************************************/
NLM_EXTERN void  UDV_Init_ScaleData(UnDViewerGraphDataPtr GrDataPtr)
{
	/*display options*/
	GrDataPtr->udv_panel.ShowFeatures=TRUE;
	GrDataPtr->udv_panel.ShowScale=TRUE;
	GrDataPtr->DisplayOptions=DDV_DISP_VERT;

	/*scale option*/
	GrDataPtr->udv_scale.ScalePosition=(Int2)SCALE_POS_BOTH;
	GrDataPtr->udv_scale.ShowMajorTick=TRUE;
	GrDataPtr->udv_scale.ShowMMinorTick=TRUE;
	GrDataPtr->udv_scale.MajorTickEvery=(Int2)SCALE_MAJOR_TICK_EVERY;
	GrDataPtr->udv_scale.MinorTickEvery=(Int2)SCALE_MINOR_TICK_EVERY;
	GrDataPtr->udv_scale.ScaleColor=(Uint4)SCALE_LETTER_COLOR;
	GrDataPtr->udv_scale.TickMajColor=(Uint4)SCALE_MAJTICK_COLOR;
	GrDataPtr->udv_scale.TickMinColor=(Uint4)SCALE_MINTICK_COLOR;
	GrDataPtr->udv_scale.cxLeftScale=(Int2)SCALE_WIDTH_IF_LEFT*
					GrDataPtr->udv_font.cxChar;


}

/*****************************************************************************

Function: UDV_Init_GraphData()

Purpose: Init Graphical values for the UnDViewer. Called one times at the
		start of the soft.

Return value: none (results are in GrDataPtr)

*****************************************************************************/
NLM_EXTERN void  UDV_Init_GraphData(PaneL p,UnDViewerGraphDataPtr GrDataPtr)
{
RecT rc;	
	
	UDV_Init_ScaleData(GrDataPtr);
	/*Panel Size & Lines to display*/
	if (p!=NULL){/*Vibrant viewer*/
		ObjectRect(p,&rc);
		UDV_ComputePanelSize(rc,&GrDataPtr->udv_panel.cxClient,
			&GrDataPtr->udv_panel.cyClient);
		GrDataPtr->udv_panel.cxName=PANEL_NAME_WIDTH*GrDataPtr->udv_font.cxChar;
		UDV_ComputeBlockByLine(GrDataPtr->udv_panel.cxClient,
			GrDataPtr->udv_panel.cxName,
			GrDataPtr->udv_scale.cxLeftScale,GrDataPtr->udv_font.cxChar,
			&GrDataPtr->udv_panel.nCharByLine,
			&GrDataPtr->udv_panel.nBlockByLine);
	}

}

/*****************************************************************************

Function: UDV_BuildFeatColorTable()

Purpose: build a simple colour table to draw Nuc/Prot features

Return value: a pointer to a colout table (each colour is encoded as an Uint4)

*****************************************************************************/
NLM_EXTERN Uint4Ptr  UDV_BuildFeatColorTable(void)
{
Uint4Ptr pClr;
Uint1 i;

	pClr=(Uint4Ptr)MemNew(FEATDEF_MAX*sizeof(Uint4));
	if (!pClr) return(NULL);
	for (i=0;i<FEATDEF_MAX;i++)	
		pClr[i]=GetColorRGB(0,0,0);

	/*Warning: values are one-based*/
	pClr[FEATDEF_GENE]=GetColorRGB(0,0,0);			/*black*/
	
	pClr[FEATDEF_precursor_RNA]=GetColorRGB(255,150,255);/*pink*/
	pClr[FEATDEF_misc_RNA]=GetColorRGB(255,150,255);
	pClr[FEATDEF_preRNA]=GetColorRGB(255,150,255);
	pClr[FEATDEF_mRNA]=GetColorRGB(255,150,255);
	pClr[FEATDEF_tRNA]=GetColorRGB(255,150,255);
	pClr[FEATDEF_rRNA]=GetColorRGB(255,150,255);
	pClr[FEATDEF_snRNA]=GetColorRGB(255,150,255);
	pClr[FEATDEF_scRNA]=GetColorRGB(255,150,255);
	pClr[FEATDEF_otherRNA]=GetColorRGB(255,150,255);
	pClr[FEATDEF_prim_transcript]=GetColorRGB(255,150,255);
	
	pClr[FEATDEF_CDS]=GetColorRGB(0,0,255);	/*blue*/
	pClr[FEATDEF_exon]=GetColorRGB(0,0,255);
	pClr[FEATDEF_intron]=GetColorRGB(0,0,255);

	pClr[FEATDEF_PROT]=GetColorRGB(0,128,0);		/*dk green*/
	pClr[FEATDEF_mat_peptide]=GetColorRGB(0,128,0);
	pClr[FEATDEF_sig_peptide]=GetColorRGB(0,128,0);
	pClr[FEATDEF_transit_peptide]=GetColorRGB(0,128,0);
	pClr[FEATDEF_preprotein]=GetColorRGB(0,128,0);
	pClr[FEATDEF_mat_peptide_aa]=GetColorRGB(0,128,0);
	pClr[FEATDEF_sig_peptide_aa]=GetColorRGB(0,128,0);
	pClr[FEATDEF_transit_peptide_aa]=GetColorRGB(0,128,0);


	pClr[FEATDEF_SITE]=GetColorRGB(255,0,0);		/*red*/

	pClr[FEATDEF_REGION]=GetColorRGB(210,154,14);		/*orange*/
	pClr[FEATDEF_mutation]=GetColorRGB(210,154,14);
	pClr[FEATDEF_variation]=GetColorRGB(210,154,14);

	pClr[FEATDEF_PSEC_STR]=GetColorRGB(104,201,220); /*cyan*/

	pClr[FEATDEF_HET]=GetColorRGB(128,128,0);	/*yellow*/

	pClr[FEATDEF_BOND]=GetColorRGB(255,92,255);	/*pink*/

	return(pClr);
}

/*****************************************************************************

Function: UDV_Build_NA_LayoutPalette()

Purpose: build a simple layout table to draw colour Nuc sequence

Return value: none (results store in GrData)

*****************************************************************************/
NLM_EXTERN void  UDV_Build_NA_LayoutPalette(UnDViewerGraphDataPtr GrData)
{

	GrData->NA_LayoutPal[NA_A_LAYOUT].LetClr=GetColorRGB(0,192,0);
	GrData->NA_LayoutPal[NA_A_LAYOUT].bkClr=0;
	GrData->NA_LayoutPal[NA_A_LAYOUT].AspLet=LAYOUT_LOWER_CASE;

	GrData->NA_LayoutPal[NA_G_LAYOUT].LetClr=GetColorRGB(192,192,0);
	GrData->NA_LayoutPal[NA_G_LAYOUT].bkClr=0;
	GrData->NA_LayoutPal[NA_G_LAYOUT].AspLet=LAYOUT_LOWER_CASE;

	GrData->NA_LayoutPal[NA_C_LAYOUT].LetClr=GetColorRGB(0,0,255);
	GrData->NA_LayoutPal[NA_C_LAYOUT].bkClr=0;
	GrData->NA_LayoutPal[NA_C_LAYOUT].AspLet=LAYOUT_LOWER_CASE;

	GrData->NA_LayoutPal[NA_T_LAYOUT].LetClr=GetColorRGB(255,0,0);
	GrData->NA_LayoutPal[NA_T_LAYOUT].bkClr=0;
	GrData->NA_LayoutPal[NA_T_LAYOUT].AspLet=LAYOUT_LOWER_CASE;

	GrData->NA_LayoutPal[NA_U_LAYOUT].LetClr=GetColorRGB(192,0,192);
	GrData->NA_LayoutPal[NA_U_LAYOUT].bkClr=0;
	GrData->NA_LayoutPal[NA_U_LAYOUT].AspLet=LAYOUT_LOWER_CASE;
}

/*****************************************************************************

Function: UDV_Build_AA_LayoutPalette()

Purpose: build a simple layout table to draw colour Prot sequence

Return value: none (results store in GrData)

*****************************************************************************/
NLM_EXTERN void  UDV_Build_AA_LayoutPalette(UnDViewerGraphDataPtr GrData)
{
RGBColoR LetClr[26]={{0,0,0},	/* A */
					{0,0,0},	/* B */
					{0,0,0},	/* C */
					{255,0,0},	/* D */
					{255,0,0},	/* E */
					{0,192,0},	/* F */
					{0,0,0},	/* G */
					{0,0,0},	/* H */
					{0,192,0},	/* I */
					{0,0,0},	/* J */
					{0,0,255},	/* K */
					{0,192,0},	/* L */
					{0,192,0},	/* M */
					{255,0,0},	/* N */
					{0,0,0},	/* O */
					{255,0,0},	/* P */
					{255,0,0},	/* Q */
					{0,0,255},	/* R */
					{0,0,0},	/* S */
					{0,0,0},	/* T */
					{0,0,0},	/* U */
					{0,192,0},	/* V */
					{0,192,0},	/* W */
					{0,0,0},	/* X */
					{0,192,0},	/* Y */
					{0,0,0},	/* Z */
					};

RGBColoR BkClr[26]={{255,255,255},	/* A */
					{255,255,255},	/* B */
					{255,255,255},	/* C */
					{255,255,255},	/* D */
					{255,255,255},	/* E */
					{255,255,255},	/* F */
					{255,255,255},	/* G */
					{255,255,255},	/* H */
					{255,255,255},	/* I */
					{255,255,255},	/* J */
					{255,255,255},	/* K */
					{255,255,255},	/* L */
					{255,255,255},	/* M */
					{255,255,255},	/* N */
					{255,255,255},	/* O */
					{255,255,255},	/* P */
					{255,255,255},	/* Q */
					{255,255,255},	/* R */
					{255,255,255},	/* S */
					{255,255,255},	/* T */
					{255,255,255},	/* U */
					{255,255,255},	/* V */
					{255,255,255},	/* W */
					{255,255,255},	/* X */
					{255,255,255},	/* Y */
					{255,255,255},	/* Z */
					};
Uint1 i;

	for(i=0;i<26;i++){
		GrData->AA_LayoutPal[i].LetClr=GetColorRGB(LetClr[i].red,
													LetClr[i].green,
													LetClr[i].blue);
		GrData->AA_LayoutPal[i].bkClr=GetColorRGB(BkClr[i].red,
													BkClr[i].green,
													BkClr[i].blue);
		GrData->AA_LayoutPal[i].AspLet=LAYOUT_UPPER_CASE;
	}	
}

/*****************************************************************************

Function: UDV_calc_pos()

Purpose: draw a vertical cursor to track position of a feature

Parameters:	GrData; graphical data object
			pt; mouse position
			rcClip; region where pgp is included
			pgp; ParaG data
			szPos; used to put pos as a string

Return value: pos (1-based value)

*****************************************************************************/
static Int4  UDV_calc_pos(UnDViewerGraphDataPtr GrData,PoinT pt,
		RecT rcClip,ParaGPtr pgp,CharPtr szPos)
{
Int2 xMargin,i;
Int4 pos,start,stop;


	xMargin=rcClip.left+GrData->udv_font.cxChar;

	pos=(Int4)((pt.x-xMargin+GrData->udv_font.cxChar/2)/
		GrData->udv_font.cxChar);

	for (i=0;i<GrData->udv_panel.nBlockByLine+1;i++){
		start=i*LETTER_BLOCK_WIDTH+i;
		stop=start+LETTER_BLOCK_WIDTH;
		if (pos>=start && pos<stop) {pos=pos-i; break;}
		if (pos==stop) {pos=pos-i-1; break;}
	}
	pos+=pgp->StartLetter;
	if (pos<pgp->StartLetter) pos=pgp->StartLetter;
	if (pos>pgp->StopLetter) pos=pgp->StopLetter;
	pos++;
	
	if (szPos!=NULL) sprintf(szPos,"%d",pos);
	
	return(pos);
}

/*****************************************************************************

Function: UDV_find_item()

Purpose: find if a feature is selected (generally in reponse to ObjMgr)

Parameters: ParaG; list to check
			bsp; identifier of the sequence
			entityID, 
			itemID; identifiers of the feature
			szFeatureName; name of the feature, if found (return value)
			pgpFeat; ParaG where the feature is located
			nLineDecal; line number within ParaG where the feature is located
			isp; identifier of the feature (eID, iID, iType)
			
Return value: TRUE if feature has been found in the ParaG list

*****************************************************************************/
static Boolean  UDV_find_item(ValNodePtr ParaG,BioseqPtr bsp,
		Uint2 entityID,Uint2 itemID,Uint2 itemType,CharPtr szFeatureName,
		ParaGPtr PNTR pgpFeat,Uint2Ptr nLineDecal,UDV_Item_SelectPtr isp)
{
ParaGPtr pgp;				/*data to check*/
Uint2 j;					/*little counter*/
ValNodePtr vnp,vnp2;		/*use to scan features*/
Uint2 iID,idx;				/*used to retrieve a desired Feature*/
SeqMgrFeatContext context;	/*used to retrieve feature data*/
Char szBuf[255]={""};		/*feature name */
Char szTmp[255]={""};		/*feature name */
Boolean bFound=FALSE;		/*TRUE if feature found*/
	
	if (ParaG){
		for(vnp=ParaG ; vnp!=NULL ; vnp=vnp->next){
			if (vnp->data.ptrvalue){
				pgp=(ParaGPtr)vnp->data.ptrvalue;
				/*look through features*/
				for(j=0,vnp2=pgp->pFeatList;j<pgp->nFeat;j++,vnp2=vnp2->next){
					if (vnp2 == NULL) break;

					UDV_BigDecodeIdxFeat ((Uint8)vnp2->data.bigintvalue, &iID,&idx,
							NULL,NULL);
					if (iID==itemID){
						if (!SeqMgrGetDesiredFeature(entityID,
							bsp,iID,idx,NULL,&context)) {
							break;
						}
						FeatDefLabel (context.sfp, szBuf, sizeof (szBuf) - 1, 
							OM_LABEL_BOTH);
						if (context.numivals>1){
							sprintf(szTmp,"%s [from %d to %d in %d segments]. ",
								szBuf,context.left+1,context.right+1,
								context.numivals);
						}
						else{
							sprintf(szTmp,"%s [from %d to %d]. ",
								szBuf,context.left+1,context.right+1);
						}
						StringCpy(szBuf,szTmp);
						sprintf(szTmp,"Length = %d. ",
							context.right-context.left+1);
						StringCat(szBuf,szTmp);

						if (context.sfp->comment) {
							sprintf(szTmp,"Note: %s. ",context.sfp->comment);
							StringCat(szBuf,szTmp);
						}
						/*store szBuf in szFeatureName*/
						if (szFeatureName) 
							StringCpy(szFeatureName,szBuf);
						*pgpFeat=pgp;
						*nLineDecal+=j;
						isp->eIDsel=entityID;
						isp->iIDsel=itemID;
						isp->iTypeSel=itemType;
						bFound=TRUE;
						break;
					}
				}
			}
			if (bFound) break;
		}		
	}
	return(bFound);
}

/*****************************************************************************

Function: UDV_click_item()

Purpose: find whether the user has clicked on a feature or not (retrieve
	values which will be sent to ObjMgr)

Parameters: GrData; graphical data object
			pt; mouse position
			rc; ParaG region
			pgp: ParaG data
			bsp_i; bsp data
			entityID, 
			itemID,
			itemType; identifiers of the features
			index; index as defined by the Feature Index for a specific feature

Return value: TRUE if  the user has clicked within a feat

*****************************************************************************/
static Boolean  UDV_click_item(UnDViewerGraphDataPtr GrData,PoinT pt,
		RecT rc,ParaGPtr pgp,BspInfo bsp_i,Uint2Ptr entityID,
		Uint2Ptr itemID,Uint2Ptr itemType,Uint2Ptr index)
{
Boolean InFeat=FALSE;		/*TRUE if the user has clicked within a feat*/
ValNodePtr vnp;				/*use to scan features*/
Int2 i,j,yDecal=0/*1*/;			/*use to find a feat region*/
Int4 pos,OccupyTo;			/*position (# of letters) of the mouse*/
Uint2 iID,idx;				/*used to retrieve a desired Feature*/
SeqMgrFeatContext context;	/*used to retrieve feature data*/
Boolean IsTransNeeded=FALSE;

	if (pgp->pFeatList==NULL) return(InFeat);

	/*go at the bottom of the sequence*/
	if (GrData->udv_scale.ScalePosition==SCALE_POS_TOP || 
					GrData->udv_scale.ScalePosition==SCALE_POS_BOTH) yDecal++;
	if (GrData->udv_scale.ShowMajorTick) yDecal++;
	
	yDecal=rc.top+yDecal*GrData->udv_font.LineHeight;
	OccupyTo=pgp->StopLetter;
	
	for(j=0,vnp=pgp->pFeatList;j<pgp->nFeat;j++,vnp=vnp->next){
		if (vnp == NULL) break;
		i=1;
		/*get desired feature given iID and idx*/
		UDV_BigDecodeIdxFeat ((Uint8)vnp->data.bigintvalue, &iID,&idx,NULL,NULL);
		if (!SeqMgrGetDesiredFeature(bsp_i.bsp_entityID,bsp_i.bsp,
            iID,idx,NULL,&context)) {
			break;
		}
		else{
			/*CDS occupies 1 or 2 lines*/
			if (context.sfp->data.choice==SEQFEAT_CDREGION){
				if (UDV_IsTranslationNeeded(&context,pgp)) i=2;
				else i=1;
			}
			else i=1;

			if (context.left<=OccupyTo)
				yDecal+=GrData->udv_font.LineHeight;

			if (pt.y>=yDecal && pt.y<=yDecal+i*GrData->udv_font.LineHeight){
				/*pos is a 1-based value*/
				if(pt.x<GrData->udv_panel.cxName-2)
					pos=context.left;
				else pos=UDV_calc_pos(GrData,pt,rc,pgp,NULL)-1;
				/*HET feature correction for the ends*/
				if(context.featdeftype== FEATDEF_HET){
					context.right=context.ivals[2*context.numivals-1];
				}

				if (pos>=context.left && pos<=context.right) {
					*entityID=bsp_i.bsp_entityID;
					*itemID=iID;
					*itemType=OBJ_SEQFEAT;
					*index=context.index;
					InFeat=TRUE;
					break;
				}
			}
			OccupyTo=context.right;
		}
	}
	return(InFeat);
}


/*****************************************************************************

Function: UDV_draw_vert_cursor()

Purpose: draw a vertical cursor to track position of a feature

Parameters:	rcClip; draw only when the mouse is within rcClip
			pos; new position to draw

Return value: none

*****************************************************************************/
static void UDV_draw_vert_cursor(RecT rcClip,PoinT pos)
{
	if (PtInRect(pos,&rcClip)){
		Dotted();
		MoveTo(pos.x,rcClip.top);
		LineTo(pos.x,rcClip.bottom);
		Solid();
	}
}

/*****************************************************************************

Function: UDV_draw_horz_cursor()

Purpose: draw a horizontal cursor to track position of a feature

Parameters:	GrData; graphical data object
			rcClip; draw only when the mouse is within rcClip
			pos; new position to draw

Return value: none

*****************************************************************************/
static void 	UDV_draw_horz_cursor(UnDViewerGraphDataPtr GrData,
			RecT rcClip,PoinT pos)
{
	
	if (PtInRect(pos,&rcClip)){
		PoinT arrow[3];
		Dotted();
		MoveTo((Int2)(rcClip.left-GrData->udv_scale.cxLeftScale),pos.y);
		LineTo(pos.x,pos.y);
		Solid();
		arrow[0].x=rcClip.left-GrData->udv_scale.cxLeftScale;
		arrow[0].y=pos.y;
		arrow[1].x=arrow[0].x+GrData->udv_font.cxChar/2;
		arrow[1].y=pos.y-GrData->udv_font.cxChar/2;
		arrow[2].x=arrow[1].x;
		arrow[2].y=pos.y+GrData->udv_font.cxChar/2;
		PaintPoly(3,arrow);
	}
}

/*****************************************************************************

Function: UDV_draw_double_cursor()

Purpose: draw a double vertical cursor to modify the cxName value (i.e. the
		width of the region where names are drawn)

Parameters:	rcClip; draw only when the mouse is within rcClip
			pos; new position to draw

Return value: none

*****************************************************************************/
static void UDV_draw_double_cursor(RecT rcClip,PoinT pos)
{
	if (PtInRect(pos,&rcClip)){
		Dotted();
		MoveTo(pos.x,rcClip.top);
		LineTo(pos.x,rcClip.bottom);
		MoveTo((Int2)(pos.x+3),rcClip.top);
		LineTo((Int2)(pos.x+3),rcClip.bottom);
		Solid();
	}
}

/*****************************************************************************

Function: UDV_select_feature()

Purpose: select a feature. This function is generally called in reponse to a
	message: " Hey Man, the user has clicked on a feature".

Parameters:	p; UDV panel
            vdp; main data struct
			entityID, itemID; the feature
			bRepos; move the undviewer panel if true. When the user clicks on
			this panel, this value is false. When the user clicks on the Feature
			Dialog Box, this value is TRUE in order to show to the nice user
			the start of the clicked feature.

Return value: none

*****************************************************************************/
NLM_EXTERN void UDV_select_feature(PaneL p,ViewerDialogDataPtr vdp,
			Uint2 entityID,Uint2 itemID,Boolean bRepos)
{
/*ViewerDialogDataPtr vdp;*/
ParaGPtr 			pgp_select;	/*selected pgp*/
RecT 				rcI;		/*InfoPanel Rect size*/
CharPtr 			szFeatName; /*used to modify the InfoPanel Text*/
WindoW 				temport;	/*to allow a safe draw*/
Uint2 				nLineDecal=0;
Int4 				decal;
Boolean				ShowTop=FALSE,
					ShowTick=FALSE;

	if (vdp == NULL) return;

	if (vdp->InfoPanel) szFeatName=(CharPtr) GetObjectExtra (vdp->InfoPanel);
	else szFeatName=NULL;

	/*if szFeatName is Null, not a problem: UDV_find_feature()
		does nothing with it*/
	vdp->Old_Item_select=vdp->Item_select;
	/*compute the position (number of lines) of the first feature within the ParaG
	by taking into account the scale*/
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
	nLineDecal=1;
	if (ShowTop) nLineDecal++;
	if (ShowTick) nLineDecal++;

	if (UDV_find_item(vdp->ParaG,vdp->bsp_i.bsp,
		entityID,itemID,OBJ_SEQFEAT,szFeatName,&pgp_select,&nLineDecal,
		&vdp->Item_select)){
			decal=pgp_select->StartLine+(Int4)nLineDecal;
			if ((decal<vdp->udv_graph.udv_vscrl.ScrollPos ||
				decal>(vdp->udv_graph.udv_vscrl.ScrollPage+
				vdp->udv_graph.udv_vscrl.ScrollPos)) && bRepos){
				Int4 CurPos;
				RecT rcP;
				
				CurPos=decal;
				UnDViewerVScrlUpdate(p,FALSE,CurPos);
				/*create new buffer*/
				UDV_create_buffer(&vdp->udv_graph,vdp->ParaG,&vdp->bsp_i,NULL);
				temport=SavePort((WindoW)p);
				Select(p);
				ObjectRect(p,&rcP);
				InvalRect(&rcP);
				RestorePort(temport);
			}
			else UDV_draw_viewer(p);

			vdp->Old_Item_select.eIDsel=(Uint2)-1;
			vdp->Old_Item_select.iIDsel=(Uint2)-1;
			vdp->Old_Item_select.iTypeSel=(Uint2)-1;
	}
	if (vdp->InfoPanel) {/*InfoPanel doesn't exist for slave viewer*/
		/*use only for the Autonomous version of the viewer*/
		temport=SavePort((WindoW)vdp->InfoPanel);
		Select(vdp->InfoPanel);
		ObjectRect(vdp->InfoPanel,&rcI);
		InvalRect(&rcI);
		RestorePort(temport);
	}
}

/*****************************************************************************

Function: UDV_ClickMouse()

Purpose: manage Mouse-Click

Parameters:	p; panel handle (currently unused)
			vdp; viewer data struct
			pt; new position

Note: this function MUST be called by external software using this viewer


Return value: none

*****************************************************************************/
NLM_EXTERN void UDV_ClickMouse(PaneL p, ViewerDialogDataPtr vdp,PoinT pt)
{
ValNodePtr 			vnp;
ParaGPtr 			pgp=NULL;
RecT 				rc /* ,rc2 */ ,rcP;
Int4 				decal_top;
Boolean 			bFind=FALSE;
Int2				decal_left;

	if (vdp == NULL) return;
	
	Select(p);
	
	/*Is mouse is in the ParaG or Name region ?*/
	if ((pt.x>(vdp->udv_graph.udv_panel.cxName+
			vdp->udv_graph.udv_scale.cxLeftScale+
			vdp->udv_graph.udv_font.cxChar)) /*|| 
			pt.x<vdp->udv_graph.udv_panel.cxName-2*/){

		ObjectRect(p,&rcP);
		decal_top=vdp->udv_graph.udv_vscrl.ScrollPos*
				vdp->udv_graph.udv_font.LineHeight-VIEWER_VERT_MARGIN;
		decal_left=vdp->udv_graph.udv_panel.cxName+rcP.left+
				vdp->udv_graph.udv_scale.cxLeftScale;

		/*look for the ParaG*/
		if (vdp->ParaG){
			for(vnp=vdp->ParaG ; vnp!=NULL ; vnp=vnp->next){
				if (vnp->data.ptrvalue){
					pgp=(ParaGPtr)vnp->data.ptrvalue;

					/*Compute Int2 rect*/
					if (((pgp->StartLine+pgp->nLines)*
							vdp->udv_graph.udv_font.LineHeight+rcP.top)
							>=decal_top){
						rc=UDV_calc_RCdraw(pgp->StartLine,pgp->nLines,rcP,
							decal_left,decal_top,
							vdp->udv_graph.udv_font.LineHeight);
						
						/*rc2.left=0;
						rc2.right=vdp->udv_graph.udv_panel.cxName-2;
						rc2.top=rc.top;
						rc2.bottom=rc.bottom;*/
						/*   ParaG           or  Name         region*/
						if (PtInRect(pt,&rc) /*|| PtInRect(pt,&rc2)*/){
							vdp->UDV_ms.rcClip=rc;
							bFind=TRUE;
							break;
						}
					}
				}
			}		
		}
		if (bFind){
			Boolean ClickFeat=FALSE;
			Uint2 entityID;
			Uint2 itemID;
			Uint2 itemType;
			Uint2 index;
			/*is the user clicks on a feature*/
			if (vdp->udv_graph.udv_panel.ShowFeatures){
				ClickFeat=UDV_click_item(&vdp->udv_graph,pt,rc,pgp,vdp->bsp_i,
					&entityID,&itemID,&itemType,&index);
				
			}
			if (ClickFeat){
				ObjMgrSendMsg(OM_MSG_SELECT,entityID,itemID,OBJ_SEQFEAT);
			}
			/*click outside a feature*/
			if (!ClickFeat){/*disable feat select*/
				/*prepare the deselection*/
				vdp->Old_Item_select=vdp->Item_select;
				vdp->Item_select.eIDsel=(Uint2)-1;
				vdp->Item_select.iIDsel=(Uint2)-1;
				vdp->Item_select.iTypeSel=(Uint2)-1;
				/*draw viewer: do the deselection*/
				UDV_draw_viewer(p);
				/*feat select is now invalidated*/
				vdp->Item_select=vdp->Old_Item_select;
			}
			if (!ClickFeat){/*show cursor*/
				WindoW temport;

				temport=SavePort((WindoW)p);
				Select(p);
				vdp->UDV_ms.oldPos=pt;
				vdp->UDV_ms.newPos=pt;
				vdp->UDV_ms.Action_type=MS_ACTION_FEAT_CURSOR;
				InvertMode();	
				UDV_draw_vert_cursor(vdp->UDV_ms.rcClip,
						vdp->UDV_ms.oldPos);
				UDV_draw_horz_cursor(&vdp->udv_graph,vdp->UDV_ms.rcClip,
						vdp->UDV_ms.oldPos);
				vdp->UDV_ms.pgp=pgp;
				MoveTo((Int2)(pt.x+2),(Int2)(pt.y-2));
				UDV_calc_pos(&vdp->udv_graph,pt,vdp->UDV_ms.rcClip,
					vdp->UDV_ms.pgp,vdp->UDV_ms.szPos);
#ifndef WIN_MSWIN
				if (PtInRect(vdp->UDV_ms.oldPos,&vdp->UDV_ms.rcClip)){
					SelectFont(vdp->udv_graph.udv_font.hFnt);
					PaintString(vdp->UDV_ms.szPos);
				}
#endif
				PlusCursor();
				RestorePort(temport);
				Update();
			}
		}
	}

	/*Is mouse within the Name list region (left of ParaG) ?*/
	if (pt.x>=vdp->udv_graph.udv_panel.cxName && 
		pt.x<=vdp->udv_graph.udv_panel.cxName+3){
			ObjectRect(p,&rc);
			rc.left=5*vdp->udv_graph.udv_font.cxChar; /*min value*/
			rc.right=2*PANEL_NAME_WIDTH*vdp->udv_graph.udv_font.cxChar;/*max val*/
			pt.x=vdp->udv_graph.udv_panel.cxName;
			vdp->UDV_ms.oldPos=pt;
			vdp->UDV_ms.newPos=pt;
			vdp->UDV_ms.rcClip=rc;
			vdp->UDV_ms.Action_type=MS_ACTION_RESIZE_WIN;
			InvertMode();	
			UDV_draw_double_cursor(vdp->UDV_ms.rcClip,
				vdp->UDV_ms.oldPos);
			CrossCursor();
	}

}

/*****************************************************************************

Function: UDV_ClickProc()

Purpose: manage Mouse-Click; this function is designed for the autonomous viewer
		to manage the InfoPanel and the Features List Dialog Box
		DON'T USE FOR EXTERNAL SOFTWARE...
		
Parameters:	p; panel handle (currently unused)
			pt; new mouse position

Return value: none

*****************************************************************************/
NLM_EXTERN void UDV_ClickProc(PaneL p, PoinT pt)
{
WindoW 				w;
ViewerDialogDataPtr vdp;
ValNodePtr 			vnp;
ParaGPtr 			pgp=NULL;
RecT 				rc /* ,rc2 */ ,rcP;
Int4 				decal_top;
Boolean 			bFind=FALSE;
ViewerMainPtr		vmp;
Int2				decal_left;

	/*get the parent Window of the Panel*/
	w=ParentWindow(p);
	
	if (w==NULL) return;

	/*get some usefull data...*/
	vmp = (ViewerMainPtr) GetObjectExtra (w);

	if (vmp==NULL) return;
	
	vdp=vmp->vdp;

	if (vdp == NULL) return;
	
	/*Select(p);*/
	
	/*Is mouse is in the ParaG or Name region ?*/
	if ((pt.x>(vdp->udv_graph.udv_panel.cxName+
			vdp->udv_graph.udv_scale.cxLeftScale+
			vdp->udv_graph.udv_font.cxChar)) /*|| 
			pt.x<vdp->udv_graph.udv_panel.cxName-2*/){

		ObjectRect(p,&rcP);
		decal_top=vdp->udv_graph.udv_vscrl.ScrollPos*
				vdp->udv_graph.udv_font.LineHeight-VIEWER_VERT_MARGIN;
		decal_left=vdp->udv_graph.udv_panel.cxName+rcP.left+
				vdp->udv_graph.udv_scale.cxLeftScale;

		/*look for the ParaG*/
		if (vdp->ParaG){
			for(vnp=vdp->ParaG ; vnp!=NULL ; vnp=vnp->next){
				if (vnp->data.ptrvalue){
					pgp=(ParaGPtr)vnp->data.ptrvalue;

					/*Compute Int2 rect*/
					if (((pgp->StartLine+pgp->nLines)*
							vdp->udv_graph.udv_font.LineHeight+rcP.top)
							>=decal_top){
						rc=UDV_calc_RCdraw(pgp->StartLine,pgp->nLines,rcP,
							decal_left,decal_top,
							vdp->udv_graph.udv_font.LineHeight);
						
						/*rc2.left=0;
						rc2.right=vdp->udv_graph.udv_panel.cxName-2;
						rc2.top=rc.top;
						rc2.bottom=rc.bottom;*/
						/*   ParaG           or  Name         region*/
						if (PtInRect(pt,&rc) /*|| PtInRect(pt,&rc2)*/){
							vdp->UDV_ms.rcClip=rc;
							bFind=TRUE;
							break;
						}
					}
				}
			}		
		}
		if (bFind){
			Boolean ClickFeat=FALSE;
			Uint2 entityID;
			Uint2 itemID;
			Uint2 itemType;
			Uint2 index;
			/*is the user clicks on a feature*/
			if (vdp->udv_graph.udv_panel.ShowFeatures){
				ClickFeat=UDV_click_item(&vdp->udv_graph,pt,rc,pgp,vdp->bsp_i,
					&entityID,&itemID,&itemType,&index);
				/*future implementation: if ClickFeat==TRUE, send a message 
				SELECT to ObjManager*/
				if (ClickFeat) {
					ObjMgrSendMsg(OM_MSG_SELECT,entityID,itemID,OBJ_SEQFEAT);
					if (vmp->hFeatDlg){/*update Features List Dlg if needed*/
						FLMDataPtr	pflm=
								(FLMDataPtr)GetObjectExtra(vmp->hFeatDlg);
						if (pflm){
							Int1 Choice;
							ScanFeatForSelect 	sffs;
							Boolean SeqFeatListAvail[SEQFEAT_MAX];

							Choice=pflm->SeqFeatClass[GetValue(pflm->pop)-1];
							sffs.eID=entityID;		
							sffs.iID=itemID;		
							sffs.index=0;
							sffs.compteur=0;
							MemSet(sffs.SeqFeatListAvail,0,
									sizeof(sffs.SeqFeatListAvail));
							if (Choice==0){/*0 means All*/
								MemSet(SeqFeatListAvail,1,
										sizeof(SeqFeatListAvail));
							}
							else{/*otherwise search only a particular feat.*/
								MemSet(SeqFeatListAvail,0,
										sizeof(SeqFeatListAvail));
								SeqFeatListAvail[Choice]=TRUE;
							}	
							SeqMgrExploreFeatures (vmp->vdp->bsp_i.bsp, 
								(Pointer) &sffs,UDV_FeaturesListBoxFind, 
								NULL, SeqFeatListAvail, NULL);
							if (sffs.index>0)
								OwnerDrawLbox_SelectItem(pflm->lbox,sffs.index);
						}
					}
					/*UDV_select_feature(vdp,entityID,itemID,FALSE);*/
				}
			}
			/*click outside a feature*/
			if (!ClickFeat){/*disable feat select*/
				WindoW temport;
				RecT rcI;
				CharPtr szFeatName;
		
				/*prepare the deselection*/
				vdp->Old_Item_select=vdp->Item_select;
				vdp->Item_select.eIDsel=(Uint2)-1;
				vdp->Item_select.iIDsel=(Uint2)-1;
				vdp->Item_select.iTypeSel=(Uint2)-1;
				/*draw viewer: do the deselection*/
				UDV_draw_viewer(p);
				/*erase the InfoPanel string properly*/
				if (vdp->InfoPanel){
					szFeatName=(CharPtr) GetObjectExtra (vdp->InfoPanel);
					MemSet(szFeatName,0,sizeof(szFeatName));
					temport=SavePort((WindoW)vdp->InfoPanel);
					Select(vdp->InfoPanel);
					ObjectRect(vdp->InfoPanel,&rcI);
					InvalRect(&rcI);
					RestorePort(temport);
				}
				/*feat select is now invalidated*/
				vdp->Item_select=vdp->Old_Item_select;
			}
			if (!ClickFeat){/*show cursor*/
				WindoW temport;

				temport=SavePort((WindoW)p);
				Select(p);
				vdp->UDV_ms.oldPos=pt;
				vdp->UDV_ms.newPos=pt;
				vdp->UDV_ms.Action_type=MS_ACTION_FEAT_CURSOR;
				InvertMode();	
				UDV_draw_vert_cursor(vdp->UDV_ms.rcClip,
						vdp->UDV_ms.oldPos);
				UDV_draw_horz_cursor(&vdp->udv_graph,vdp->UDV_ms.rcClip,
						vdp->UDV_ms.oldPos);
				vdp->UDV_ms.pgp=pgp;
				MoveTo((Int2)(pt.x+2),(Int2)(pt.y-2));
				UDV_calc_pos(&vdp->udv_graph,pt,vdp->UDV_ms.rcClip,
					vdp->UDV_ms.pgp,vdp->UDV_ms.szPos);
#ifndef WIN_MSWIN
				if (PtInRect(vdp->UDV_ms.oldPos,&vdp->UDV_ms.rcClip)){
					SelectFont(vdp->udv_graph.udv_font.hFnt);
					PaintString(vdp->UDV_ms.szPos);
				}
#endif
				PlusCursor();
				RestorePort(temport);
				Update();
			}
		}
	}

	/*Is mouse within the Name list region (left of ParaG) ?*/
	if (pt.x>=vdp->udv_graph.udv_panel.cxName && 
		pt.x<=vdp->udv_graph.udv_panel.cxName+3){
			ObjectRect(p,&rc);
			rc.left=5*vdp->udv_graph.udv_font.cxChar; /*min value*/
			rc.right=2*PANEL_NAME_WIDTH*vdp->udv_graph.udv_font.cxChar;/*max val*/
			pt.x=vdp->udv_graph.udv_panel.cxName;
			vdp->UDV_ms.oldPos=pt;
			vdp->UDV_ms.newPos=pt;
			vdp->UDV_ms.rcClip=rc;
			vdp->UDV_ms.Action_type=MS_ACTION_RESIZE_WIN;
			InvertMode();	
			UDV_draw_double_cursor(vdp->UDV_ms.rcClip,
				vdp->UDV_ms.oldPos);
			CrossCursor();
	}
}


/*****************************************************************************

Function: UDV_DragMouse()

Purpose: manage Mouse-Drag

Parameters:	p; panel handle (currently unused)
			vdp; viewer data struct
			pt; new position

Note: this function MUST be called by external software using this viewer

Return value: none

*****************************************************************************/
NLM_EXTERN void UDV_DragMouse(PaneL p, ViewerDialogDataPtr vdp,PoinT pt)
{
WindoW 				temport;

	if (vdp == NULL) return;

	if (vdp->UDV_ms.Action_type==MS_ACTION_FEAT_NOTHING) return;
	
	if (vdp->UDV_ms.Action_type==MS_ACTION_FEAT_CURSOR){
		InvertMode();	
		UDV_draw_vert_cursor(vdp->UDV_ms.rcClip,
			vdp->UDV_ms.oldPos);
		UDV_draw_horz_cursor(&vdp->udv_graph,vdp->UDV_ms.rcClip,
			vdp->UDV_ms.oldPos);
			MoveTo((Int2)(vdp->UDV_ms.oldPos.x+2),(Int2)(vdp->UDV_ms.oldPos.y-2));
#ifndef WIN_MSWIN
		if (PtInRect(vdp->UDV_ms.oldPos,&vdp->UDV_ms.rcClip)){
			SelectFont(vdp->udv_graph.udv_font.hFnt);
			PaintString(vdp->UDV_ms.szPos);
		}
#endif
		vdp->UDV_ms.newPos=pt;
		UDV_draw_vert_cursor(vdp->UDV_ms.rcClip,
			vdp->UDV_ms.newPos);
		UDV_draw_horz_cursor(&vdp->udv_graph,vdp->UDV_ms.rcClip,
			vdp->UDV_ms.newPos);
		MoveTo((Int2)(vdp->UDV_ms.newPos.x+2),(Int2)(vdp->UDV_ms.newPos.y-2));
		UDV_calc_pos(&vdp->udv_graph,pt,vdp->UDV_ms.rcClip,vdp->UDV_ms.pgp,
				vdp->UDV_ms.szPos);
#ifndef WIN_MSWIN
		if (PtInRect(vdp->UDV_ms.newPos,&vdp->UDV_ms.rcClip)){
			PaintString(vdp->UDV_ms.szPos);
		}
#endif
		vdp->UDV_ms.oldPos=vdp->UDV_ms.newPos;
	}
	
	if (vdp->UDV_ms.Action_type==MS_ACTION_RESIZE_WIN){
		InvertMode();	
		UDV_draw_double_cursor(vdp->UDV_ms.rcClip,
			vdp->UDV_ms.oldPos);
		vdp->UDV_ms.newPos=pt;
		UDV_draw_double_cursor(vdp->UDV_ms.rcClip,
			vdp->UDV_ms.newPos);
		vdp->UDV_ms.oldPos=vdp->UDV_ms.newPos;
	}
	Update();
}

/*****************************************************************************

Function: UDV_DragProc()

Purpose: manage Mouse-Drag

Parameters:	p; panel handle (currently unused)
			pt; new position

Return value: none

*****************************************************************************/
NLM_EXTERN void UDV_DragProc(PaneL p, PoinT pt)
{
WindoW 				w;
ViewerDialogDataPtr vdp;
ViewerMainPtr		vmp;

	/*get the parent Window of the Panel*/
	w=ParentWindow(p);
	
	if (w==NULL) return;

	/*get some usefull data...*/
	vmp = (ViewerMainPtr) GetObjectExtra (w);

	if (vmp==NULL) return;
	
	vdp=vmp->vdp;
	
	UDV_DragMouse(p,vdp,pt);
}

/*****************************************************************************

Function: UDV_ReleaseMouse()

Purpose: manage Mouse-Release

Parameters:	p; panel handle (currently unused)
			vdp; viewer data struct
			pt; new position

Note: this function MUST be called by external software using this viewer

Return value: none

*****************************************************************************/
NLM_EXTERN void UDV_ReleaseMouse(PaneL p,ViewerDialogDataPtr vdp, PoinT pt)
{
WindoW 	temport;


	if (vdp == NULL) return;
	
	if (vdp->UDV_ms.Action_type==MS_ACTION_FEAT_NOTHING) return;

	temport=SavePort((WindoW)p);
	Select(p);

	if (vdp->UDV_ms.Action_type==MS_ACTION_FEAT_CURSOR){
		InvertMode();	
		UDV_draw_vert_cursor(vdp->UDV_ms.rcClip,
			vdp->UDV_ms.oldPos);
		UDV_draw_horz_cursor(&vdp->udv_graph,vdp->UDV_ms.rcClip,
			vdp->UDV_ms.oldPos);
		MoveTo((Int2)(vdp->UDV_ms.oldPos.x+2),(Int2)(vdp->UDV_ms.oldPos.y-2));
		UDV_calc_pos(&vdp->udv_graph,pt,vdp->UDV_ms.rcClip,vdp->UDV_ms.pgp,
				vdp->UDV_ms.szPos);
#ifndef WIN_MSWIN
		if (PtInRect(vdp->UDV_ms.oldPos,&vdp->UDV_ms.rcClip)){
			SelectFont(vdp->udv_graph.udv_font.hFnt);
			PaintString(vdp->UDV_ms.szPos);
		}
#endif
	}

	if (vdp->UDV_ms.Action_type==MS_ACTION_RESIZE_WIN){
		InvertMode();	
		UDV_draw_double_cursor(vdp->UDV_ms.rcClip,
			vdp->UDV_ms.oldPos);
		/*redraw panel with new 'cxName' value*/
		if (PtInRect(vdp->UDV_ms.newPos,&vdp->UDV_ms.rcClip)){
			RecT rc;

			/*Select(p);*/
			ObjectRect(p,&rc);
			rc.left=0;/*for obscure reasons, not == 0*/
			rc.top=0;
			vdp->udv_graph.udv_panel.cxName=vdp->UDV_ms.newPos.x;
			UDV_resize_viewer( p, vdp);
			InvalRect(&rc);
		}

	}	

	vdp->UDV_ms.oldPos.x=0;
	vdp->UDV_ms.oldPos.y=0;
	vdp->UDV_ms.newPos.x=0;
	vdp->UDV_ms.newPos.y=0;
	LoadRect(&vdp->UDV_ms.rcClip,0,0,0,0);
	vdp->UDV_ms.Action_type=MS_ACTION_FEAT_NOTHING;
	vdp->UDV_ms.pgp=NULL;
	ArrowCursor();
	RestorePort(temport);
}

/*****************************************************************************

Function: UDV_ReleaseProc()

Purpose: manage Mouse-Release

Parameters:	p; panel handle (currently unused)
			pt; new position

Return value: none

*****************************************************************************/
NLM_EXTERN void UDV_ReleaseProc(PaneL p, PoinT pt)
{
WindoW 				w;
ViewerDialogDataPtr vdp;
ViewerMainPtr		vmp;

	/*get the parent Window of the Panel*/
	w=ParentWindow(p);
	
	if (w==NULL) return;

	/*get some usefull data...*/
	vmp = (ViewerMainPtr) GetObjectExtra (w);

	if (vmp==NULL) return;
	
	vdp=vmp->vdp;

	UDV_ReleaseMouse(p,vdp,pt);
}

/*******************************************************************************

  Function : UDV_draw_helixDNA()
  
  Purpose : draw a little peace of a DNA helix for the NCBI logo
  
  Parameters :	x,y; start position of the helix
  				cxChar,cyChar; size of the current font
				bEnd; if True; draw half an helix

  Return value : none

*******************************************************************************/
void static  UDV_draw_helixDNA(Int2 x,Int2 y,Int2 cxChar,Int2 cyChar,
		Boolean bEnd)
{

Int2 cxChar_2,cyChar_2;
PoinT box[4];

	cxChar_2=cxChar/2;
	cyChar_2=cyChar/2;
	SetColor(GetColorRGB(0,0,255));
	/*1*/
	box[0].x=x+cxChar_2;
	box[0].y=y+2*cyChar;
	box[1].x=x;
	box[1].y=box[0].y-cyChar;
	box[2].x=x+cxChar;
	box[2].y=box[1].y;
	box[3].x=box[2].x+cxChar_2;
	box[3].y=box[0].y;
	
	PaintPoly(4,box);

	SetColor(GetColorRGB(0,0,128));
	/*2*/
	box[0].x=box[3].x;
	box[0].y=box[3].y;
	box[1].x=box[2].x;
	box[1].y=box[2].y;
	box[2].x=box[1].x+2*cxChar+cxChar_2;
	box[2].y=box[0].y-5*cyChar;
	box[3].x=box[2].x+cxChar_2;
	box[3].y=box[2].y+cyChar;
	PaintPoly(4,box);
	
	if (!bEnd){
		SetColor(GetColorRGB(0,0,255));
		/*3*/
		box[0].x=box[3].x;
		box[0].y=box[3].y;
		box[1].x=box[2].x;
		box[1].y=box[2].y;
		box[2].x=box[1].x+cxChar;
		box[2].y=box[1].y;
		box[3].x=box[2].x+cxChar_2;
		box[3].y=box[0].y;
		PaintPoly(4,box);
	}

	SetColor(GetColorRGB(128,0,128));
	/*4*/
	box[0].x=x;
	box[0].y=y-2*cyChar;
	box[1].x=x+cxChar_2;
	box[1].y=box[0].y-cyChar;
	box[2].x=box[0].x+cxChar;
	box[2].y=box[0].y;
	PaintPoly(3,box);

	if (!bEnd){
		SetColor(GetColorRGB(192,92,192));
		/*5*/
		box[0].x=box[1].x;
		box[0].y=box[1].y;
		box[1].x=box[0].x+cxChar;
		box[1].y=box[0].y;
		box[2].x=box[1].x+3*cxChar;
		box[2].y=box[0].y+5*cyChar;
		box[3].x=box[2].x-cxChar;
		box[3].y=box[2].y;
		PaintPoly(4,box);

		SetColor(GetColorRGB(128,0,128));
		/*7*/
		box[0].x=box[2].x;
		box[0].y=box[2].y;
		box[1].x=box[0].x-cxChar_2;
		box[1].y=box[0].y-cyChar;
		box[2].x=box[1].x+2*cxChar;
		box[2].y=box[1].y-3*cyChar;
		box[3].x=box[2].x+cxChar;
		box[3].y=box[2].y;
		PaintPoly(4,box);

		SetColor(GetColorRGB(0,0,255));
		/*8*/
		box[0].x=x+4*cxChar;
		box[0].y=y-2*cyChar;
		box[1].x=box[0].x+cxChar;
		box[1].y=box[0].y;
		box[2].x=box[0].x+2*cxChar+cxChar_2+cxChar_2;
		box[2].y=box[0].y+3*cyChar;
		box[3].x=box[2].x-cxChar;
		box[3].y=box[2].y;
		PaintPoly(4,box);
	}
	Black();
}

/*******************************************************************************

  Function : draw_logo()
  
  Purpose : draw the NCBI logo for the main window
  
  Parameters :	rcP; where to put the logo
  				f1,f2,f3; the fonts used to display the logo text

  Return value : none

*******************************************************************************/
static void  draw_logo(RecT rcP,FonT f1,FonT f2,FonT f3,CharPtr szTxt0,
	CharPtr szTxt4)
{
/*Char szTxt0[]="UnD-Viewer";
Char szTxt4[]=", a sequence viewer for GenBank";*/
Char szTxt1[]="01101011";
Char szTxt2[]="National Center for";
Char szTxt3[]="Biotechnology Information";
Int2 pos,cxChar=7,cyChar=13,cxGlobal;
Int2 w1,w2,w3,y,x;

	/*some computation for positionning*/
	SelectFont(f1);
	w1=StringWidth(szTxt1);
	SelectFont(f2);
	w2=StringWidth(szTxt2);
	w3=StringWidth(szTxt3);
	cxGlobal=22*cxChar+w1+w2;
	
	x=(rcP.right-rcP.left)/2-(cxGlobal/2);
	y=(rcP.bottom-rcP.top)/2;
	
	/*double DNA helix*/
	UDV_draw_helixDNA(x,y,cxChar,cyChar,FALSE);
	UDV_draw_helixDNA((Int2)(x+6*cxChar),y,cxChar,cyChar,FALSE);
	UDV_draw_helixDNA((Int2)(x+12*cxChar),y,cxChar,cyChar,TRUE);

	/*Undviewer*/
	SelectFont(f3);
	Blue();
	MoveTo((Int2)(x+2*cxChar),(Int2)(y-6*cyChar));
	PaintString(szTxt0);

	Black();
	MoveTo((Int2)(x+StringWidth(szTxt0)+3*cxChar),(Int2)(y-6*cyChar));
	SelectFont(f1);
	PaintString(szTxt4);

	/*01101011*/
	pos=x+13*cxChar;
	SelectFont(f1);
	Black();
	MoveTo(pos,y);
	PaintString(szTxt1);
	pos+=w1;

	/*NCBI*/
	SelectFont(f2);
	pos+=(w3/2+2*cxChar);

	MoveTo((Int2)(pos-w2/2),(Int2)(y-cyChar));
	PaintString(szTxt2);
	MoveTo((Int2)(pos-w3/2),(Int2)(y+cyChar));
	PaintString(szTxt3);
}

/*******************************************************************************

  Function : UDV_draw_select()
  
  Purpose : draw a rect around a selected item
  
  Parameters :	rc; rect coordinates

  Return value : none

*******************************************************************************/
static void  UDV_draw_select(RecT rc)
{
	Dotted();
	FrameRect(&rc);
	Solid();
}

/*******************************************************************************

  Function : UDV_draw_ss_cont()
  
  Purpose : complete the end(s) of secondary structure to highlight the case
  			of a truncation at ParaG ends
  
  Parameters :	start; draw the arrow from...
				stop;... to.
				StartLetter; beginning of a ParaG (zero-based value)
				cxChar; width of a letter
				xMargin; draw in the ParaG on the right a this value
				y2; vertical position to draw the sequence
				ContinueRight; draw on the right only if TRUE
				ContinueLeft; draw on the left only if TRUE

  Return value : none

*******************************************************************************/
static void  UDV_draw_ss_cont(Int4 start,Int4 stop,Int4 StartLetter,
		Int2 cxChar,Int2 xMargin,Int2 y2,Boolean ContinueRight,
		Boolean ContinueLeft)
{
Int2 x;

	if (ContinueRight || ContinueLeft){
		Dotted();
		Red();
	}
	if (ContinueLeft){
		x=(((start-StartLetter)+(start-StartLetter)/
			LETTER_BLOCK_WIDTH)*cxChar)+xMargin-cxChar/2;
		
		MoveTo(x,y2);
		LineTo((Int2)(x-2*cxChar),y2);
	}
	if (ContinueRight){
		Int4 stop2;
		stop2=stop+((stop-start)/LETTER_BLOCK_WIDTH);
		x=(((stop2-StartLetter)+(stop2-StartLetter)/
			LETTER_BLOCK_WIDTH)*cxChar)+xMargin;
		
		MoveTo(x,y2);
		LineTo((Int2)(x+2*cxChar),y2);
	}
	
	if (ContinueRight || ContinueLeft){
		Solid();
		Black();
	}
}
/*******************************************************************************

  Function : UDV_draw_struc_strand()
  
  Purpose : draw the strand (prot struct only) as an arrow
  
  Parameters :	start; draw the arrow from...
				stop;... to.
				StartLetter; beginning of a ParaG (zero-based value)
				y; vertical position to draw the sequence
				cxChar; width of a letter
				cyChar; height of a letter
				LineHeight; height of a feature line
				xMargin; draw in the ParaG on the right a this value
				clr; color
				ContinueRight; draw on the right only if TRUE
				ContinueLeft; draw on the left only if TRUE

  Return value : none

*******************************************************************************/
static void  UDV_draw_struc_strand(Int4 start,Int4 stop,Int4 StartLetter,
			Int2 y,Int2 cxChar,Int2 cyChar,Int2 LineHeight,Int2 xMargin,
			Uint4 clr,Boolean ContinueLeft,Boolean ContinueRight)
{
PoinT box[7];
Int2 y2,x,cxChar_2,cyChar_2,cyChar_4;

	y2=y-LineHeight/2;
	cxChar_2=cxChar/2;
	cyChar_4=cyChar/4;
	cyChar_2=cyChar/2;
	
	/*ContinueRight || ContinueLef: draw thin line to highlight continuous
	secondary structure at ParaG ends*/
	if (ContinueLeft || ContinueRight)
		UDV_draw_ss_cont(start,stop,StartLetter,cxChar,xMargin,y2,ContinueRight,
		ContinueLeft);
	if (start!=stop){
		x=(((start-StartLetter)+(start-StartLetter)/
			LETTER_BLOCK_WIDTH)*cxChar)+xMargin+cxChar_2;

		box[0].x=x-cxChar_2;
		box[0].y=y2-cyChar_4;
		box[6].x=box[0].x;
		box[6].y=y2+cyChar_4;

		/*3d effect*/
		DkGray();
		MoveTo((Int2)(box[0].x-1),(Int2)(box[0].y+1));
		LineTo((Int2)(box[6].x-1),(Int2)(box[6].y+1));
		MoveTo((Int2)(box[0].x-2),(Int2)(box[0].y+2));
		LineTo((Int2)(box[6].x-2),(Int2)(box[6].y+2));

		/*stop+=((stop-start)/LETTER_BLOCK_WIDTH);*/
		x=(((stop-StartLetter)+(stop-StartLetter)/
			LETTER_BLOCK_WIDTH)*cxChar)+xMargin;
		x-=cxChar_2;

		box[1].x=x-cxChar_2;
		box[1].y=box[0].y;
		box[2].x=box[1].x;
		box[2].y=y2-cyChar_2;
		box[3].x=x+cxChar_2;
		box[3].y=y2;
		box[4].x=box[1].x;
		box[4].y=y2+cyChar_2;
		box[5].x=box[1].x;
		box[5].y=box[6].y;

		/*3d effect*/
		MoveTo((Int2)(box[6].x-1),(Int2)(box[6].y+1));
		LineTo((Int2)(box[5].x),(Int2)(box[5].y+1));
		LineTo((Int2)(box[4].x),(Int2)(box[4].y+1));
		MoveTo((Int2)(box[6].x-2),(Int2)(box[6].y+2));
		LineTo((Int2)(box[5].x-1),(Int2)(box[5].y+2));
		LineTo((Int2)(box[4].x-1),(Int2)(box[4].y+2));

		MoveTo(box[4].x,(Int2)(box[4].y+2));
		LineTo(box[3].x,(Int2)(box[3].y+2));

		if (clr!=(Uint4)-1) SetColor(clr);
		else Black();
		/*draw*/
		PaintPoly(7,box);
	}
	else{/*if start==stop, draw only the end of the arrow*/
		stop+=((stop-start)/LETTER_BLOCK_WIDTH);
		x=(((stop-StartLetter)+(stop-StartLetter)/
			LETTER_BLOCK_WIDTH)*cxChar)+xMargin/*-cxChar_2*/;

		box[0].x=x-cxChar_2;
		box[0].y=y2-cyChar_2;
		box[1].x=x+cxChar_2;
		box[1].y=y2;
		box[2].x=box[0].x;
		box[2].y=y2+cyChar_2;

		/*3d effect*/
		DkGray();
		MoveTo((Int2)(box[0].x-1),(Int2)(box[0].y+1));
		LineTo((Int2)(box[2].x-1),(Int2)(box[2].y+1));
		MoveTo((Int2)(box[0].x-2),(Int2)(box[0].y+2));
		LineTo((Int2)(box[2].x-2),(Int2)(box[2].y+2));
		MoveTo((Int2)(box[2].x),(Int2)(box[2].y+2));
		LineTo((Int2)(box[1].x),(Int2)(box[1].y+2));

		if (clr!=(Uint4)-1) SetColor(clr);
		else Black();
		/*draw*/
		PaintPoly(3,box);
	}
	Black();	
}


/*******************************************************************************

  Function : UDV_draw_struc_helix()
  
  Purpose : draw the helix (prot struct only) as an helix
  
  Parameters :	start; draw the arrow from...
				stop;... to.
				StartLetter; beginning of a ParaG (zero-based value)
				y; vertical position to draw the sequence
				cxChar; width of a letter
				cyChar; height of a letter
				LineHeight; height of a feature line
				xMargin; draw in the ParaG on the right a this value
				clr; color
				ContinueRight; draw on the right only if TRUE
				ContinueLeft; draw on the left only if TRUE
	   
  Return value : none

*******************************************************************************/
static void  UDV_draw_struc_helix(Int4 start,Int4 stop,Int4 StartLetter,
			Int2 y,Int2 cxChar,Int2 cyChar,Int2 LineHeight,Int2 xMargin,
			Uint4 clr,Boolean ContinueLeft,Boolean ContinueRight)
{
PoinT box[4],pt0,pt2;
Int4 i;
Int2 y2,x,xDecal,cxChar_2,cyChar_2,cxChar_4;
Int1 status=DRAW_HELIX_DOWN;
Boolean bContPrev=FALSE;

	cxChar_2=cxChar/2;
	cxChar_4=cxChar/4;
	cyChar_2=cyChar/2;

	y2=y-LineHeight/2;

	/*ContinueRight || ContinueLef: draw thin line to highlight continuous
	secondary structure at ParaG ends*/
	if (ContinueLeft || ContinueRight)
		UDV_draw_ss_cont(start,stop,StartLetter,cxChar,xMargin,y2,ContinueRight,
		ContinueLeft);

	/*stop+=((stop-start)/LETTER_BLOCK_WIDTH);*/
	stop+=((stop-StartLetter)/LETTER_BLOCK_WIDTH)-
		((start-StartLetter)/LETTER_BLOCK_WIDTH);

	xDecal=((start-StartLetter)+((start-StartLetter)/LETTER_BLOCK_WIDTH))
			*cxChar+xMargin;
	for(i=start;i<=stop;i++){
		x=(i-start)*cxChar+xDecal;
		switch(status){
			case DRAW_HELIX_DOWN:
				/*compute polygon*/
				box[0].x=x;
				box[0].y=y2+cyChar_2;
				box[1].x=x-cxChar_2;
				box[1].y=y2;
				box[2].x=x+cxChar_4;
				box[2].y=y2;
				box[3].x=x+cxChar_2+cxChar_4;
				box[3].y=box[0].y;
				/*set color*/
				if (clr!=(Uint4)-1) SetColor(clr);
				else LtGray();
				/*draw*/
				PaintPoly(4,box);
				if (bContPrev){
					box[1]=pt0;
					box[3]=box[2];
					box[2]=pt2;
					bContPrev=FALSE;
					PaintPoly(4,box);
				}
				/*after DRAW_HELIX_DOWN will trace a DRAW_HELIX_MIDDLE*/
				status=DRAW_HELIX_MIDDLE;				
				break;
			case DRAW_HELIX_MIDDLE:
				/*compute polygon*/
				box[0].x=x-cxChar_4;
				box[0].y=y2+cyChar_2;
				box[1].x=x-cxChar+cxChar_4;
				box[1].y=y2;
				box[2].x=x+cxChar_2-cxChar_4;
				box[2].y=y2-cyChar_2;
				box[3].x=x+cxChar-cxChar_4;
				box[3].y=y2;
				/*set color: always dark gray*/
				DkGray();
				/*draw*/
				PaintPoly(4,box);
				/*after DRAW_HELIX_MIDDLE will trace a DRAW_HELIX_UP*/
				status=DRAW_HELIX_UP;				
				break;
			case DRAW_HELIX_UP:
				/*compute polygon*/
				box[0].x=x-cxChar_4;
				box[0].y=y2;
				box[1].x=x-cxChar_2-cxChar_4;
				box[1].y=y2-cyChar_2;
				box[2].x=x;
				box[2].y=box[1].y;
				box[3].x=x+cxChar_2;
				box[3].y=y2;
				/*set color: always dark gray*/
				if (clr!=(Uint4)-1) SetColor(clr);
				else LtGray();
				/*draw*/
				PaintPoly(4,box);
				/*after DRAW_HELIX_UP will trace a DRAW_HELIX_DOWN*/
				status=DRAW_HELIX_DOWN;
				bContPrev=TRUE;
				pt0=box[0];
				pt2=box[2];
				break;
		}
	}
	Black();
}

/*******************************************************************************

  Function : UDV_draw_CDS_minus()
  
  Purpose : draw the translation of a CDS located on a MINUS strand
  
  Parameters : context; feature data (see explore.h for a def of the structure)
  				start; draw the sequence from...
				stop;... to.
				StartLetter; beginning of a ParaG (zero-based value)
				GrData; graphical data
				y; vertical position to draw the sequence
				i; order number of the exon (if CDS is exon-encoded)
				xMargin; draw in the ParaG on the right a this value
				UseDefClr; TRUE: draw the sequence using the DefClr value
					Otherwise use GrData->AA_LayoutPal[] values
  
  Return value : none

*******************************************************************************/
static void  UDV_draw_CDS_minus(SeqMgrFeatContextPtr context,Int4 start,
			Int4 stop,Int4 StartLetter,UnDViewerGraphDataPtr GrData,Int2 y,
			Int2 i,Int2 xMargin,Boolean UseDefClr,Uint4 DefClr)
{						
CharPtr      	str=NULL;
Int4			il=0,
				pos=0,
				n=0,
				stop_prot=0;
Int2 			nCompt=0,
				decal=0,
				x_prot=0,
				y_prot=0,
				numivals2=0,
				cxChar_2;
ByteStorePtr	bs=NULL;

	/*retrieve the protein sequence; need to be optimized in future release*/
	bs = ProteinFromCdRegion (context->sfp, TRUE);
	if (bs){
		str = (CharPtr) BSMerge (bs, NULL);
		BSFree (bs);
	}
	else return;

	if (!str) return;	

	/*count the lenght of the intron(s)*/
	numivals2=context->numivals*2;
	if (i>0){
		for(nCompt=i;nCompt>0;nCompt-=2){
			il+=(context->ivals[nCompt-2]-
				context->ivals[nCompt+1]-1);
		}
	}

	pos=context->ivals[1]-il-stop;
	/*place the current position within the correct translation frame*/
	if (pos % 3){
		if (!((pos-1) % 3)){
			pos=(pos-1)/3;
			decal=-1;
		}
		else if (!((pos+1) % 3)){
			pos=(pos+1)/3;
			decal=1;
		}						
	}else {
		pos=pos/3;
		decal=0;
	}

	n=pos;
	stop_prot=stop-decal-1; /*-1 == middle of the codon*/
	cxChar_2=GrData->udv_font.cxChar/2;
	
	/*vertical pos: just below the CDS colour box (pos: y)*/
	y_prot=y+GrData->udv_font.LineHeight-GrData->udv_font.cyChar/2;
	/*NOTICE : I draw from right to left; whereas normally I draw always 
		from left to right*/				
	if (UseDefClr) SetColor(DefClr);
	
	while(stop_prot>=start){
		x_prot=((stop_prot-StartLetter)+
				((stop_prot-StartLetter)/LETTER_BLOCK_WIDTH))*
				GrData->udv_font.cxChar+xMargin-cxChar_2;
		if (!UseDefClr){
			if (str[n]>='A' && str[n]<='Z'){
				SetColor(GrData->AA_LayoutPal[str[n]-'A'].LetClr);
			}
			else{
				Black();
			}
		}
		MoveTo(x_prot,y_prot);
		PaintChar(str[n]);
		stop_prot-=3;
		n++;
	}

	if (str){
		MemFree (str);
	}
	Black();
}


/*******************************************************************************

  Function : UDV_draw_CDS_plus()
  
  Purpose : draw the translation of a CDS located on a PLUS strand
  
  Parameters : context; feature data (see explore.h for a def of the structure)
  				start; draw the sequence from...
				stop;... to.
				StartLetter; beginning of a ParaG (zero-based value)
				GrData; graphical data
				y; vertical position to draw the sequence
				i; order number of the exon (if CDS is exon-encoded)
				xMargin; draw in the ParaG on the right a this value
				UseDefClr; TRUE: draw the sequence using the DefClr value
					Otherwise use GrData->AA_LayoutPal[] values
  
  Return value : none

*******************************************************************************/
static void  UDV_draw_CDS_plus(SeqMgrFeatContextPtr context,Int4 start,
			Int4 stop,Int4 StartLetter,UnDViewerGraphDataPtr GrData,Int2 y,
			Int2 i,Int2 xMargin,Boolean UseDefClr,Uint4 DefClr)
{						
CharPtr      	str=NULL;
Int4			il=0,
				pos=0,
				n=0,
				start_prot=0;
Int2 			n2=0,
				nCompt=0,
				decal=0,
				x_prot=0,
				y_prot=0,
				numivals2=0;
ByteStorePtr	bs=NULL;
Char szBuf[SZBUF_SIZE]={""};

	/*retrieve the protein sequence; need to be optimized in future release*/
	bs = ProteinFromCdRegion (context->sfp, TRUE);
	if (bs){
		str = (CharPtr) BSMerge (bs, NULL);
		BSFree (bs);
	}
	else return;

	if (!str) return;	
	MemSet(szBuf,' ',SZBUF_SIZE-1);	
	
	/*count the lenght of the intron(s)*/
	numivals2=context->numivals*2;
	if (i>0 && i<numivals2){
		for(nCompt=2;nCompt<i+2;nCompt+=2){
			il+=(context->ivals[nCompt]-
				context->ivals[nCompt-1]-1);
		}
	}

	pos=start-il-context->ivals[0];
						
	/*place the current position within the correct translation frame*/
	if (pos % 3){
		if (!((pos-1) % 3)){
			pos=(pos-1)/3;
			decal=-1;
		}
		else if (!((pos+1) % 3)){
			pos=(pos+1)/3;
			decal=1;
		}						
	}else {
		pos=pos/3;
		decal=0;
	}

	n=pos;
	start_prot=start+decal+1;

	/*vertical pos: just below the CDS colour box (pos: y)*/
	y_prot=y+GrData->udv_font.LineHeight-GrData->udv_font.cyChar/2;
	x_prot=((start_prot-StartLetter)+
			((start_prot-StartLetter)/LETTER_BLOCK_WIDTH))*
			GrData->udv_font.cxChar+xMargin-GrData->udv_font.cxChar/2;

	if (UseDefClr) SetColor(DefClr);
	nCompt=start_prot-StartLetter+(start_prot-StartLetter)/LETTER_BLOCK_WIDTH;	
	while(start_prot<=stop){
		n2=start_prot-StartLetter+(start_prot-StartLetter)/LETTER_BLOCK_WIDTH
			-nCompt;
		if (n2<SZBUF_SIZE-1)
			szBuf[n2]=str[n];
		start_prot+=3;
		/*nCompt+=3;*/
		n++;
	}

	MoveTo(x_prot,y_prot);
	PaintString(szBuf);
	if (str){
		MemFree (str);
	}
	Black();
}


/*******************************************************************************

  Function : UDV_draw_big_arrow_feat()
  
  Purpose : draw a single feature as an arrow
  
  Parameters :  start_x; draw from...
  				stop_x; ...to
				y: vertical position
				cxChar; width of a character
				cyChar; height of a character
				LineHeight; height of a line in a ParaG
				start_nat; start is a box, arrow or nothing ?
 				stop_nat; stop is a box, arrow or nothing ?
				clr; colour 

  Return value : none

*******************************************************************************/
static void  UDV_draw_big_arrow_feat(Int2 start_x,Int2 y,
			Int2 cxChar,Int2 cyChar,Uint4 clr)
{
PoinT arrow[3];

	if (clr!=(Uint4)-1) SetColor(clr);
	else Black();
	y-=3;
	arrow[0].x=start_x-cxChar/2;
	arrow[0].y=y;
	arrow[1].x=start_x+cxChar/2;
	arrow[1].y=y;
	arrow[2].x=start_x;
	arrow[2].y=y-(2*cyChar)/3;

	PaintPoly(3,arrow);
	Black();
}

/*******************************************************************************

  Function : UDV_draw_big_arrow_bond()
  
  Purpose : draw a bond feature as an arrow with a connected line
  
  Parameters :  start_x; draw from...
  				stop_x; ...to
				y: vertical position
				cxChar; width of a character
				cyChar; height of a character
				LineHeight; height of a line in a ParaG
				start_nat; start is a box, arrow or nothing ?
 				stop_nat; stop is a box, arrow or nothing ?
				clr; colour 

  Return value : none

*******************************************************************************/
static void  UDV_draw_big_arrow_bond(Int2 start_x,Int2 y,
			Int2 cxChar,Int2 cyChar,Uint4 clr,Int1 What)
{
PoinT arrow[3];

	if (clr!=(Uint4)-1) SetColor(clr);
	else Black();
	y-=3;
	
	arrow[0].x=start_x;
	arrow[0].y=y;
	if (What==BOND_LEFT) arrow[1].x=start_x+cxChar/2;
	else arrow[1].x=start_x-cxChar/2;
	arrow[1].y=y-cyChar/3;
	arrow[2].x=start_x;
	arrow[2].y=y-(2*cyChar)/3;
	PaintPoly(3,arrow);
	Black();
}


/*******************************************************************************

  Function : UDV_draw_big_line_feat()
  
  Purpose : draw a single feature as a box
  
  Parameters :  start_x; draw from...
  				stop_x; ...to
				y: vertical position
				cxChar; width of a character
				cyChar; height of a character
				LineHeight; height of a line in a ParaG
				start_nat; start is a box, arrow or nothing ?
 				stop_nat; stop is a box, arrow or nothing ?
				clr; colour 
 
  Return value : none

*******************************************************************************/
static void  UDV_draw_big_line_feat(Int2 start_x,Int2 stop_x,Int2 y,
			Int2 cxChar,Int2 LineHeight,Int2 start_nat,Int2 stop_nat,
			Uint4 clr)
{
RecT rcLine,rcBox;
Int2 y2,x2,l2,l3,l6;
PoinT arrow[3];

	y2=y-LineHeight/2;
	l2=LineHeight/2;
	l3=LineHeight/3;
	l6=LineHeight/6;
	
	if (clr!=(Uint4)-1) SetColor(clr);
	else Black();
	
	/*start end : box, arrow or nothing*/
	if (start_nat==FEATURE_START_BOX){
		x2=start_x;
		rcBox.left=x2-2;	
  		rcBox.top=y2-l3;	
		rcBox.right=x2+2;	
		rcBox.bottom=y2+l3;
		start_x=rcBox.right;	
		PaintRect(&rcBox);
	}
	else if (start_nat==FEATURE_START_ARROW){
		arrow[0].x=start_x;
		arrow[0].y=y2;
		arrow[1].x=start_x+cxChar;
		arrow[1].y=y2-l2;
		arrow[2].x=arrow[1].x;
		arrow[2].y=y2+l2;
		start_x=arrow[2].x;
		PaintPoly(3,arrow);
	}
	else if (start_nat==FEATURE_START_ARROW_END){
		arrow[0].x=start_x;
		arrow[0].y=y2;
		arrow[1].x=start_x+cxChar;
		arrow[1].y=y2-l2;
		arrow[2].x=arrow[1].x;
		arrow[2].y=y2+l2;
		start_x=arrow[2].x;
		PaintPoly(3,arrow);
		MoveTo(arrow[0].x,arrow[1].y);
		LineTo(arrow[0].x,arrow[2].y);
	}
	
	/*stop end : box, arrow or nothing*/
	if (stop_nat==FEATURE_START_BOX){
		x2=stop_x;
		rcBox.left=stop_x-2;	
		rcBox.top=y2-l3;	
		rcBox.right=x2+2;	
		rcBox.bottom=y2+l3;
		stop_x=rcBox.left;	
		PaintRect(&rcBox);
	}
	else if (stop_nat==FEATURE_START_ARROW){
		arrow[0].x=stop_x;
		arrow[0].y=y2;
		arrow[1].x=stop_x-cxChar;
		arrow[1].y=y2-l2;
		arrow[2].x=arrow[1].x;
		arrow[2].y=y2+l2;
		stop_x=arrow[2].x;
		PaintPoly(3,arrow);
	}
	else if (stop_nat==FEATURE_START_ARROW_END){
		arrow[0].x=stop_x;
		arrow[0].y=y2;
		arrow[1].x=stop_x-cxChar;
		arrow[1].y=y2-l2;
		arrow[2].x=arrow[1].x;
		arrow[2].y=y2+l2;
		stop_x=arrow[2].x;
		PaintPoly(3,arrow);
		MoveTo(arrow[0].x,arrow[1].y);
		LineTo(arrow[0].x,arrow[2].y);
	}

	if (start_x>=stop_x) return;
	rcLine.left=start_x;
	rcLine.top=y2-_max_(2,l6);
	rcLine.right=stop_x;
	rcLine.bottom=y2+_max_(2,l6);
	PaintRect(&rcLine);
	Black();	
}

/*******************************************************************************

  Function : UDV_draw_thin_line()
  
  Purpose : draw a thin black line (generally used to delineate the introns)
  
  Parameters : 	start_x; draw from...
  				stop_x; ...to
				y: vertical position
				LineHeight; height of a line in a ParaG
  
  Return value : none

*******************************************************************************/
static void  UDV_draw_thin_line(Int2 start_x,Int2 stop_x,Int2 y,
					Int2 LineHeight)
{
Int2 y2;

	y2=y-LineHeight/2;
	
	LtGray();
	
	MoveTo(start_x,y2);
	LineTo(stop_x,y2);

	Black();
}

/*******************************************************************************

  Function : UDV_Draw_features()
  
  Purpose : I thing the name of the function is meaningfull...
  
  Parameters : 	GrData; graphical data (font size, etc)
  				bsp_i ; general data of the Bioseq
				UseDefClr; TRUE: draw the sequence using the DefClr value
					Otherwise use GrData->AA_LayoutPal[] values
				pgp; data of the ParaG
				rc; rectangle containing this ParaG
  				pClr; colour table [FEATDEF_MAX size]
				is; item to select, if needed
				old_is; old item to deselect, if needed

  Note : don't try to understand this function alone... if you want to keep
         a good health !				

  Return value : none

*******************************************************************************/
static void  UDV_Draw_features(UnDViewerGraphDataPtr GrData,
					BspInfoPtr bsp_i,Boolean UseDefClr,Uint4 DefClr,
					ParaGPtr pgp,RecT PNTR rc,Uint4Ptr pClr,
					UDV_Item_Select is,UDV_Item_Select old_is)
{
Int2 y,ybase;				/*text position*/
Int2 xMargin;				/*initial margins*/
Int2 nTicks=0,nLet=0;		/*y decal*/
ValNodePtr vnp;				/*Features list of itemID and index values*/
Uint2 iID,idx,lineID;		/*used to retrieve a desired Feature*/
/*Uint8 index_g;			merged value containing itemID, index and lineID*/
SeqMgrFeatContextPtr context;	/*used to retrieve feature data*/
SeqMgrFeatContext context2;	/*used to retrieve feature data*/
Int4 start,stop/*,OccupyTo*/;	/*limit of the feature to draw*/
Int2 start_x,stop_x,		/*id. but in pixels*/
	start_nat,stop_nat;		/*ends of the feature : box, arrow*/
Int2 i,numivals2,i_decal,j;	/*counters*/
Uint4	clr;				/*color used to draw feature*/		
Boolean b_draw_line,bCDS,b_connect_left,b_connect_right,
		b_end_left,b_end_right;
Int2 idx1,idx2,idx3,idx4,idx5,idx6,nLines=0;
BioseqPtr parent;
SeqMgrSegmentContext contextPart;

	if (!pgp->pFeatList) return;
	
	/*compute position*/
	if (GrData->udv_scale.ShowMajorTick) nTicks++;
	if (GrData->udv_scale.ScalePosition==SCALE_POS_TOP || 
					GrData->udv_scale.ScalePosition==SCALE_POS_BOTH) nLet++;
	
	xMargin=rc->left+GrData->udv_font.cxChar;
	ybase=rc->top+(nTicks+nLet+1/*2*/)*GrData->udv_font.LineHeight;
	/*OccupyTo=pgp->StopLetter;*/

	/*the current bsp is just a segment ?*/
	parent=SeqMgrGetParentOfPart(bsp_i->bsp,&contextPart);
	context=NULL;

	/*draw : loop on all features in a ParaG*/
	for(j=0,vnp=pgp->pFeatList;j<pgp->nFeat;j++,vnp=vnp->next){
		if (vnp == NULL) break;
		/*index_g=(Uint8)vnp->data.bigintvalue;*/
		UDV_BigDecodeIdxFeat ((Uint8)vnp->data.bigintvalue, &iID,&idx,&lineID,NULL);
		/*get desired feature given iID and idx*/
		if (!SeqMgrGetDesiredFeature(bsp_i->bsp_entityID,
				(parent!=NULL ? parent : bsp_i->bsp),
                iID,idx,NULL,&context2)) continue;


		if (context) {
			MemFree(context->ivals);
			context=MemFree(context);
		}
		context=UDV_ConvertFeatContext(&context2,contextPart.cumOffset,
			bsp_i->bsp->length);
		if (!context) continue;
		/*depending on the strand, various things are possible :*/
			/*PLUS/MINUS : if region (start!=stop)-> draw big box, with arrow*/

			/*UNKNOWN : if region (start!=stop)-> draw big box, without arrow*/

			/*PLUS/MINUS/UNKNOWN : if single letter (start==stop)-> draw 
			vertical arrow*/
		
		/*HET feature correction for the ends*/
		if(context->featdeftype== FEATDEF_HET){
			context->right=context->ivals[2*context->numivals-1];
		}

		/*temporary situation; will be modified in the future*/
		if (context->strand>Seq_strand_minus ||
			context->strand==Seq_strand_unknown) context->strand=Seq_strand_plus;

		/*strand PLUS*/
		if (context->strand==Seq_strand_plus){
			numivals2=context->numivals*2;
			i=0;
			i_decal=2;
		}
		
		/*strand MINUS*/
		if (context->strand==Seq_strand_minus){
			numivals2=2*context->numivals-2;
			i=numivals2;
			i_decal=-2;
		}		

		bCDS=FALSE;
		/*if (_max_(context->ivals[i],pgp->StartLetter)<=OccupyTo){
			nLines++;
		}*/
		nLines=lineID;
		while (TRUE){
			b_connect_left=FALSE;
			b_connect_right=FALSE;
			b_end_left=FALSE;
			b_end_right=FALSE;
			b_draw_line=FALSE;

			/*if ivals.stop > end ParaG -> end of drawing*/
			if (context->ivals[i]>pgp->StopLetter) break;
			/*if ivals.stop<= start ParaG : not yet in the current ParaG*/
			
			if (context->ivals[i+1]<pgp->StartLetter) {
				if (context->strand==Seq_strand_plus || 
						context->strand==Seq_strand_unknown){
					if (numivals2>2 && i+2<numivals2){
					/*stop ParaG < start next ivals -> inter-region: fill 
					the ParaG with a thin line; this is the case
					for coding region: draw thin line to delineate the introns*/		
						if (context->ivals[i+2]>pgp->StopLetter){
							b_draw_line=TRUE;
						}
					}
				}
				if (context->strand==Seq_strand_minus){
					if (numivals2>2 && i-2>-1){
					/*stop ParaG < start next ivals -> inter-region: fill 
					the ParaG with a thin line; this is the case
					for coding region: draw thin line to delineate the introns*/		
						if (context->ivals[i-2]>pgp->StopLetter){
							b_draw_line=TRUE;
						}
					}
				}
				if (b_draw_line){
					start=0;
					/*nLines++;*/
					stop=pgp->StopLetter-pgp->StartLetter;
					start_x=start*GrData->udv_font.cxChar+xMargin;
					stop_x=(stop+(stop/LETTER_BLOCK_WIDTH))*
						GrData->udv_font.cxChar+xMargin;
					y=ybase+nLines*GrData->udv_font.LineHeight;
					UDV_draw_thin_line(start_x,stop_x,y,
							GrData->udv_font.LineHeight);
					/*OccupyTo=pgp->StopLetter;*/
				}
				if (context->strand==Seq_strand_plus || 
						context->strand==Seq_strand_unknown){
					i=i+i_decal;
					if (i>numivals2-2) break;
					else continue;
				}					
				if (context->strand==Seq_strand_minus){
					i=i+i_decal;
					if (i<0) break;
					else continue;
				}
			}
						
			/*compute the limits of the feature within ParaG*/
			start=_max_(context->ivals[i],pgp->StartLetter);
			stop=_min_(context->ivals[i+1],pgp->StopLetter);

			/*are there connections; exons/introns for example*/
				/*on the left*/
			if (context->strand==Seq_strand_plus || 
						context->strand==Seq_strand_unknown){
				if (numivals2>2 && i>0){
					if (start>pgp->StartLetter){
						idx1=i-1;
						idx2=1;
						idx3=-2;
						b_connect_left=TRUE;
						b_end_left=TRUE;
					}
					else if(start==pgp->StartLetter){
						b_end_left=TRUE;
					}
				}
			}
			if (context->strand==Seq_strand_minus){
				if (numivals2>2 && i<numivals2-1){
					if (start>pgp->StartLetter){
						idx1=i+3;
						idx2=3;
						idx3=0;
						b_connect_left=TRUE;
						b_end_left=TRUE;
					}
					else if(start==pgp->StartLetter){
						b_end_left=TRUE;
					}
				}
			}
				/*on the right*/
			if (context->strand==Seq_strand_plus || 
						context->strand==Seq_strand_unknown){
				if (numivals2>2 && i+2<numivals2){
					if (stop<pgp->StopLetter){
						idx4=i+2;
						idx5=1;
						idx6=-2;
						b_connect_right=TRUE;
						b_end_right=TRUE;
					}
					else if(stop==pgp->StopLetter){
						b_end_right=TRUE;
					}
				}
			}
			if (context->strand==Seq_strand_minus){
				if (numivals2>2 && i>0){
					if (stop<pgp->StopLetter){
						idx4=i-2;
						idx5=3;
						idx6=0;
						b_connect_right=TRUE;
						b_end_right=TRUE;
					}
					else if(stop==pgp->StopLetter){
						b_end_right=TRUE;
					}
				}
			}
			/*compute the 'nature' of start & stop: box, arrow or nothing*/
			/*I use only Seq_strand_minus or Seq_strand_plus !!!!*/
			if (context->strand==Seq_strand_plus || 
				context->strand==Seq_strand_minus){
				if (context->ivals[i]==start){
					if (context->strand==Seq_strand_plus){ 
						if (b_end_left) start_nat=FEATURE_START_NOTHING;
						else start_nat=FEATURE_START_BOX;
					}
					else {
						if (b_end_left) start_nat=FEATURE_START_ARROW;
						else start_nat=FEATURE_START_ARROW_END;
					}
				}
				else{
					start_nat=FEATURE_START_NOTHING;
				}

				if (context->ivals[i+1]==stop){
					if (context->strand==Seq_strand_minus){ 
						if (b_end_right) stop_nat=FEATURE_START_NOTHING;
						else stop_nat=FEATURE_START_BOX;
					}
					else{
						if (b_end_right) stop_nat=FEATURE_START_ARROW;
						else stop_nat=FEATURE_START_ARROW_END;
					}
				}
				else{
					stop_nat=FEATURE_START_NOTHING;
				}
			}
			else{
				start_nat=FEATURE_START_NOTHING;
				stop_nat=FEATURE_START_NOTHING;
			}

			/*compute limits of the drawing... zero based values from the left 
				of ParaG rc*/
			start_x=((start-pgp->StartLetter)+
				((start-pgp->StartLetter)/LETTER_BLOCK_WIDTH))*
				GrData->udv_font.cxChar+xMargin;
			stop_x=((stop-pgp->StartLetter)+((stop-pgp->StartLetter)/
				LETTER_BLOCK_WIDTH))*GrData->udv_font.cxChar+xMargin;

			/*colour*/
			if (pClr){
				clr=pClr[context->featdeftype];
			}
			else clr=(Uint4)-1;
			
			y=ybase+nLines*GrData->udv_font.LineHeight;
			/*draw feature*/
			if (start!=stop || (start==stop && context->left!=context->right) 
				/*|| context->strand!=Seq_strand_unknown*/){
				switch(context->featdeftype){
					case FEATDEF_HET:
						if (clr!=(Uint4)-1) SetColor(clr);
						MoveTo(start_x-GrData->udv_font.cxChar/2,y);
						PaintChar('H');
						/*UDV_draw_big_arrow_feat(start_x,y,GrData->udv_font.cxChar,
							GrData->udv_font.cyChar,clr);*/
						b_connect_left=FALSE;
						b_connect_right=FALSE;
						break;
					case FEATDEF_BOND:
						UDV_draw_thin_line(start_x,stop_x,y,
							GrData->udv_font.LineHeight);
						if (pgp->StartLetter<=context->left && 
							pgp->StopLetter>=context->left){
							UDV_draw_big_arrow_bond(start_x,y,
								GrData->udv_font.cxChar,
								GrData->udv_font.cyChar,
								clr,BOND_LEFT);
						}
						if (pgp->StartLetter<=context->right && 
							pgp->StopLetter>=context->right){
							UDV_draw_big_arrow_bond(stop_x,y,
								GrData->udv_font.cxChar,
								GrData->udv_font.cyChar,
								clr,BOND_RIGHT);
						}
						break;
					case FEATDEF_PSEC_STR:
						if (context->sfp){
							Boolean ContinueLeft=FALSE;
							Boolean ContinueRight=FALSE;
							
							if (context->sfp->data.value.intvalue==1){/*helix*/
								if (start!=context->left) ContinueLeft=TRUE;
								if (stop!=context->right)ContinueRight=TRUE;
								UDV_draw_struc_helix( start, stop,
									pgp->StartLetter,
									y,GrData->udv_font.cxChar,
									GrData->udv_font.cyChar,
									GrData->udv_font.LineHeight,xMargin,clr,
									ContinueLeft,ContinueRight);
							}
							if (context->sfp->data.value.intvalue==2){/*sheet*/
								if (start!=context->left) ContinueLeft=TRUE;
								if (stop!=context->right)ContinueRight=TRUE;
								UDV_draw_struc_strand( start, stop,
									pgp->StartLetter,
									y,GrData->udv_font.cxChar,
									GrData->udv_font.cyChar,
									GrData->udv_font.LineHeight,xMargin,clr,
									ContinueLeft,ContinueRight);
							}
						}
						break;
					case FEATDEF_CDS:
						bCDS=TRUE;
						if (context->strand==Seq_strand_plus)
							UDV_draw_CDS_plus(context,start,stop,
								pgp->StartLetter,
								GrData,y,i,xMargin, UseDefClr, DefClr);						
						if (context->strand==Seq_strand_minus)
							UDV_draw_CDS_minus(context,start,stop,
								pgp->StartLetter,
								GrData,y,i,xMargin, UseDefClr, DefClr);
						nLines++;
						/*WARNING: do not add anything between FEATDEF_CDS 
						and the following default case because after the 
						translation, I must draw a box*/						
					default:
						UDV_draw_big_line_feat(start_x,stop_x,y,
							GrData->udv_font.cxChar,
							GrData->udv_font.LineHeight,
							start_nat,stop_nat,clr);
						break;
				}
			}
			else{/*feature is one letter long*/
				UDV_draw_big_arrow_feat(start_x,y,GrData->udv_font.cxChar,
					GrData->udv_font.cyChar,clr);
			}
			/*OccupyTo=stop;*/
			/*select feature, if needed*/
			if (old_is.eIDsel!=-1){
				if (old_is.iIDsel==iID){
					RecT rcSel;

					rcSel.left=start_x-5;
					rcSel.top=y-GrData->udv_font.LineHeight+1;
					rcSel.right=stop_x+5;
					if (!bCDS) rcSel.bottom=y;
					else rcSel.bottom=y+GrData->udv_font.LineHeight;
					White();
					UDV_draw_select(rcSel);
				}
			}
			if (is.eIDsel!=-1){
				if (is.iIDsel==iID){
					RecT rcSel;

					rcSel.left=start_x-5;
					rcSel.top=y-GrData->udv_font.LineHeight+1;
					rcSel.right=stop_x+5;
					if (!bCDS) rcSel.bottom=y;
					else rcSel.bottom=y+GrData->udv_font.LineHeight;
					Red();
					UDV_draw_select(rcSel);
				}
			}
			Black();

			/*draw thin line if needed; used to connect exons/introns for 
			example*/
			if (b_connect_left){
				Int4 start2,stop2;
				stop2=start;
				start2=_max_(pgp->StartLetter,context->ivals[idx1]);
				start_x=(start2-pgp->StartLetter+
					((start2-pgp->StartLetter)/LETTER_BLOCK_WIDTH))*
					GrData->udv_font.cxChar+xMargin;
				stop_x=(stop2-pgp->StartLetter+
					((stop2-pgp->StartLetter)/
					LETTER_BLOCK_WIDTH))*GrData->udv_font.cxChar+xMargin;
				UDV_draw_thin_line((Int2)(start_x+idx2),(Int2)(stop_x+idx3),y,
						GrData->udv_font.LineHeight);
			}

			if (b_connect_right){
				Int4 start2,stop2;
				start2=stop;
				stop2=_min_(pgp->StopLetter,context->ivals[idx4]);
				start_x=(start2-pgp->StartLetter+
					((start2-pgp->StartLetter)/LETTER_BLOCK_WIDTH))*
					GrData->udv_font.cxChar+xMargin;
				stop_x=(stop2-pgp->StartLetter+
					((stop2-pgp->StartLetter)/
					LETTER_BLOCK_WIDTH))*GrData->udv_font.cxChar+xMargin;
				UDV_draw_thin_line((Int2)(start_x+idx5),(Int2)(stop_x+idx6),y,
					GrData->udv_font.LineHeight);
				/*OccupyTo=stop2;*/
			}

			/*if (context->sfp && context->sfp->data.choice==SEQFEAT_CDREGION)
					y+=GrData->udv_font.LineHeight;		*/

			if (context->strand==Seq_strand_plus || 
						context->strand==Seq_strand_unknown){
				i=i+i_decal;
				if (i>numivals2-2) break;
			}
			if (context->strand==Seq_strand_minus){
				i=i+i_decal;
				if (i<0) break;
			}
		}		
/*		y+=GrData->udv_font.LineHeight;
		if (bCDS) y+=GrData->udv_font.LineHeight;*/
	}
	Black();	
}

/*******************************************************************************

  Function : UDV_Draw_features_label()
  
  Purpose : draw the labels of each feature located in a ParaG
  
  Parameters :	GrData; graphical data (font size, etc)
  				bsp_i ; general data of the Bioseq
				pgp; data of the ParaG
				rc; rectangle containing this ParaG
  				pClr; colour table [FEATDEF_MAX size]
				is; item to select, if needed
				old_is; old item to select, if needed

  Return value : none

*******************************************************************************/
static void  UDV_Draw_features_label(UnDViewerGraphDataPtr GrData,
					BspInfoPtr bsp_i,ParaGPtr pgp,RecT PNTR rc,Uint4Ptr pClr,
					UDV_Item_Select is,UDV_Item_Select old_is)
{
Int2 x,y;					/*text position*/
Int2 xMargin;				/*initial margins*/
Int2 nTicks=0,nLet=0;		/*y decal*/
Char szBuf[51]={""};		/*name to draw*/
ValNodePtr vnp;				/*Features list of itemID and index values*/
Uint2 iID,idx;				/*used to retrieve a desired Feature*/
/*Uint4 index_g;			merged value containing itemID and index*/
SeqMgrFeatContext context;	/*used to retrieve feature data*/
Uint2 j;					/*counter*/

	if (!pgp->pFeatList) return;
	
	/*compute position*/
	if (GrData->udv_scale.ShowMajorTick) nTicks++;
	if (GrData->udv_scale.ScalePosition==SCALE_POS_TOP || 
					GrData->udv_scale.ScalePosition==SCALE_POS_BOTH) nLet++;
	
	xMargin=rc->left-2*GrData->udv_font.cxChar-GrData->udv_scale.cxLeftScale;
	y=rc->top+(nTicks+nLet+2)*GrData->udv_font.LineHeight;
	
	/*draw : loop on all features in a ParaG*/
	for(j=0,vnp=pgp->pFeatList;j<pgp->nFeat;j++,vnp=vnp->next){
		if (vnp == NULL) break;
		
		UDV_BigDecodeIdxFeat ((Uint8)vnp->data.bigintvalue, &iID,&idx,NULL,NULL);

		/*get desired feature given iID and idx*/
		if (!SeqMgrGetDesiredFeature(bsp_i->bsp_entityID,bsp_i->bsp,
                                iID,idx,NULL,&context)) continue;

		if (context.sfp){
			/*copy name*/
			if (FeatDefLabel (context.sfp, szBuf, sizeof (szBuf) - 1, 
						OM_LABEL_BOTH)){

				/*draw name*/
				x=xMargin-StringWidth(szBuf);
				MoveTo(x,y);
				if (pClr){
					SetColor(pClr[context.featdeftype]);
				}
				PaintString (szBuf);
				/*select if needed*/
				if (old_is.eIDsel!=-1){
					if (old_is.iIDsel==iID){
						RecT rcSel;
						
						rcSel.left=x-3;
						rcSel.top=y-GrData->udv_font.cyChar+1;
						rcSel.right=rcSel.left+StringWidth(szBuf)+4;
						rcSel.bottom=y+2;
						White();
						UDV_draw_select(rcSel);
					}
				}
				if (is.eIDsel!=-1){
					if (is.iIDsel==iID){
						RecT rcSel;
						
						rcSel.left=x-3;
						rcSel.top=y-GrData->udv_font.cyChar+1;
						rcSel.right=rcSel.left+StringWidth(szBuf)+4;
						rcSel.bottom=y+2;
						Red();
						UDV_draw_select(rcSel);
						
					}
				}
				Black();
				y+=GrData->udv_font.LineHeight;
				
			}
			if (context.sfp->data.choice==SEQFEAT_CDREGION){
				if (UDV_IsTranslationNeeded(&context,pgp)) 
					y+=GrData->udv_font.LineHeight;
			}
		}
	}
	Black();	
}


/*******************************************************************************

  Function : UDV_Draw_sequence_name()
  
  Purpose : draw the name of the sequence (on the left of it) of a ParaG
  
  Parameters : GrData; graphical data (font size, etc)
  				bsp_i ; general data of the Bioseq
				rc; rectangle containing this ParaG
  
  Return value : none

*******************************************************************************/
static void  UDV_Draw_sequence_name(UnDViewerGraphDataPtr GrData,
					BspInfoPtr bsp_i,RecT PNTR rc)
{
Int2 x,y,hdecal;			/*text position*/
Int2 nTicks=0,nLet=0;	/*y decal*/

	/*compute position*/
	if (GrData->udv_scale.ShowMajorTick) nTicks++;
	if (GrData->udv_scale.ScalePosition==SCALE_POS_TOP || 
					GrData->udv_scale.ScalePosition==SCALE_POS_BOTH) nLet++;
	
	hdecal=nTicks+nLet+1;
	
	/*draw name*/
	x=rc->left-2*GrData->udv_font.cxChar-GrData->udv_scale.cxLeftScale-
				StringWidth(bsp_i->bspAccNum);
	y=rc->top+hdecal*GrData->udv_font.LineHeight;
	MoveTo(x,y);
	PaintString (bsp_i->bspAccNum);
}

/*******************************************************************************

  Function : GiveClrFromIdx()
  
  Purpose : retrieve a colour
  
  Parameters : 	idx;position within the bsp (zero-based)
  				seq_color; colour table
				
  Return value : pointer to colour struct

*******************************************************************************/
static ResidueColorCellPtr GiveClrFromIdx(Int4 idx,ValNodePtr seq_color)
{
Int4  i;

	i=0;
	while(seq_color){
		if (i==idx) {
			return((ResidueColorCellPtr)seq_color->data.ptrvalue);	
		}
		seq_color=seq_color->next;
		i++;
	}

   return NULL;
}

/*******************************************************************************

  Function : UDV_Draw_sequence()
  
  Purpose : draw the Bioseq' sequence of a ParaG
  
  Parameters : 	GrData; graphical data (font size, etc)
				vnp_color; colours table
  				bsp_i ; general data of the Bioseq
				IsNuc;True if nuc sequence
				UseDefClr; TRUE: draw the sequence using the DefClr value
					Otherwise use GrData->AA_LayoutPal[] values
				pgp; data of the ParaG
				rc; rectangle containing this ParaG
			    start_decal; buffer start letter (0-based)
				StartLetter; ParaG start letter (0-based)
				szSequence; sequence buffer
				
  Return value : none

*******************************************************************************/
static void  UDV_Draw_sequence(UnDViewerGraphDataPtr GrData,ValNodePtr vnp_color,
					Boolean UseDefClr,Uint4 DefClr,
					ParaGPtr pgp,RecT PNTR rc,Int4 start_decal,Int4 StartLetter,
					CharPtr szSequence)
{
Int4 pos,taille,stop;	/*scale start at...(used to draw text)*/
Int2 x,y,ldecal,hdecal;	/*text position*/
Int2 xMargin;			/*initial margins*/
Int2 nTicks=0,nLet=0;	/*y decal*/
Uint1 residue;			/*letter and sequence length to retrieve*/
Uint2 nCompt=0,nCompt2=0;
CharPtr szBuf;
ResidueColorCellPtr rgb;
Uint4 blackColor = GetColorRGB(0,0,0),newColor,curColor;

	
	/*alloc a buffer*/
	taille=pgp->StopLetter-pgp->StartLetter+
		(pgp->StopLetter-pgp->StartLetter)/LETTER_BLOCK_WIDTH+5;
	szBuf=(CharPtr)MemNew(taille*sizeof(Char));
	if (!szBuf) return;
	MemSet(szBuf,0,taille*sizeof(Char));
	
	/*scale or not... to be or not to be ?*/
	if (GrData->udv_scale.ShowMajorTick) nTicks++;
	if (GrData->udv_scale.ScalePosition==SCALE_POS_TOP || 
					GrData->udv_scale.ScalePosition==SCALE_POS_BOTH) nLet++;
	
	pos=pgp->StartLetter+1;	/*remember : StartLetter is zero based*/
	xMargin=rc->left+GrData->udv_font.cxChar;
	y=rc->top;
	/*xDecal=xMargin;*/
	ldecal=GrData->udv_font.cxChar/2;
	hdecal=nTicks+nLet+1;

	if (UseDefClr || vnp_color==NULL) SetColor(DefClr);
	newColor=(Uint4)-1;
	stop=pgp->StopLetter+2;
	y=y+hdecal*GrData->udv_font.LineHeight;
	x=xMargin-ldecal;
	/*draw!*/
	if (vnp_color==NULL){
		while(pos<stop){
			residue=(Uint1)szSequence[StartLetter-start_decal+(nCompt++)];
			szBuf[nCompt2++]=residue;

			/*each LETTER_BLOCK_WIDTH, add a blank column*/
			if (!(pos % LETTER_BLOCK_WIDTH)) {
				szBuf[nCompt2]=' ';
				nCompt2++;
			}

			pos++;		
		}
		MoveTo(x,y);
		PaintString(szBuf);
	}
	else{ /*use vnp_color table of colours; complex colour scheme*/
		while(pos<stop){
			/*retrieve colour*/
			rgb=GiveClrFromIdx(pos-1,vnp_color);
			if(rgb != NULL ){
			   curColor = GetColorRGB (rgb->rgb[0], rgb->rgb[1],rgb->rgb[2]);
			}
			else curColor=blackColor;
			/*new colour?*/
			if (curColor!=newColor){
				newColor=curColor;
				if (szBuf[0]!='\0'){/*something to draw ?*/
					szBuf[nCompt2]='\0';/*CLOSE THE STRING*/
					MoveTo(x,y);
					SetColor(newColor);
					PaintString(szBuf);
				}
				x+=(nCompt2*GrData->udv_font.cxChar);
				nCompt2=0;
				szBuf[nCompt2]='\0';
			}
			residue=(Uint1)szSequence[StartLetter-start_decal+(nCompt++)];
			szBuf[nCompt2++]=residue;
			
			/*each LETTER_BLOCK_WIDTH, add a blank column*/
			if (!(pos % LETTER_BLOCK_WIDTH)) {
				szBuf[nCompt2]=' ';
				nCompt2++;
			}
			pos++;		
		}
		if (szBuf[0]!='\0'){/*something to draw ?*/
			szBuf[nCompt2]='\0';/*CLOSE THE STRING*/
			MoveTo(x,y);
			SetColor(newColor);
			PaintString(szBuf);
		}
	}
	Black();
	MemFree(szBuf);
}

/*******************************************************************************

  Function : UDV_Draw_scale()
  
  Purpose : draw the numerical scale of a ParaG
  
  Parameters : 	GrData; graphical data (font size, etc)
				ShowMajorTick;TRUE = show major ticks ( | )
				ShowMMinorTick;TRUE = show minor ticks ( . )
				ScalePosition; left, top, ...
				StartLetter; start scale at...
				StopLetter; ...and stop at
				rc; rectangle containing this ParaG
				LeftDecal;adjust position on the left, if needed
				ScaleMaxVal;  scale length
				UseBlockDisp; use the 10 by 10 letters block display
				
  Return value : none

*******************************************************************************/
NLM_EXTERN void UDV_Draw_scale(UnDViewerGraphDataPtr GrData,
		Boolean ShowMajorTick,Boolean ShowMMinorTick,Uint1 ScalePosition,
		Int4 StartLetter,Int4 StopLetter,RecT PNTR rc,Int2 LeftDecal,
		Int4 ScaleMaxVal,Int4 AlignPos,Boolean UseBlockDisp,Int2 ColWidth)
{
Int4 pos,AliPos;				/*scale start at...(used to draw text)*/
Int2 x,y,y2,y3,y4,xDecal;	/*text position*/
Int2 xMargin,yMargin;	/*initial margins*/
Char szBuf[15]={""};	/*scale value*/
Int2 nTicks=0,nLet=0;	/*y decal*/
Int2 cxChar_2,cyChar_4;

	if (ShowMajorTick) nTicks++;
	if (ScalePosition==SCALE_POS_TOP || 
					ScalePosition==SCALE_POS_BOTH) nLet++;
	
	pos=StartLetter+1;	/*remember : StartLetter is zero based*/
	xMargin=rc->left+LeftDecal;
	yMargin=rc->top;
	y=yMargin;
	xDecal=xMargin;
	cxChar_2=GrData->udv_font.cxChar/2;
	cyChar_4=GrData->udv_font.cyChar/4;
	y2=(Int2)(y+GrData->udv_font.LineHeight);
	AliPos=AlignPos+1;/*switch to 1-based value*/
	ScaleMaxVal++;
	/*scale on top ?*/
	if ((ScalePosition==SCALE_POS_TOP || 
					ScalePosition==SCALE_POS_BOTH)){
		while(pos<StopLetter+1){

			if (/*(pos==StartLetter+1) ||*/ !(AliPos % LETTER_BLOCK_WIDTH)
					&& AliPos!=ScaleMaxVal){
				sprintf(szBuf,"%d",AliPos);
				SetColor(GrData->udv_scale.ScaleColor);
				/*scale on top or both*/
				if(ShowMajorTick){/*center text on the tick*/
					x=xDecal-(StringWidth(szBuf))/2;
				}
				else{/*align text right on pos*/
					x=xDecal-StringWidth(szBuf)+cxChar_2;
				}
				if (ScalePosition==SCALE_POS_TOP || 
							pos!=StartLetter+1){
					MoveTo(x,y2);
					PaintString (szBuf);
				}
			}

			if (AliPos==ScaleMaxVal){/*end sequence only*/
				if (ScalePosition==SCALE_POS_TOP || 
						ScalePosition==SCALE_POS_BOTH){
					/*scale on top or both*/
					if(ShowMajorTick){
						/*center text on the tick*/
						sprintf(szBuf,"%d",AliPos);
						x=xDecal-(StringWidth(szBuf))/2;
						MoveTo(x,(Int2)(y2+GrData->udv_font.cyChar));
						SetColor(GrData->udv_scale.ScaleColor);
						PaintString (szBuf);
					}
					/*else {
						x=xDecal-StringWidth(szBuf)+GrData->udv_font.cxChar;
						MoveTo(x,y2);
					}*/
				}
			}

			xDecal+=ColWidth;

			/*each LETTER_BLOCK_WIDTH, add a blank column*/
			if (!(AliPos % LETTER_BLOCK_WIDTH)&&UseBlockDisp) xDecal+=ColWidth;
			pos++;AliPos++;
		}
	}

	/*scale on left only*/
	pos=StartLetter+1;	/*remember : StartLetter is zero based*/
	y=yMargin;
	xDecal=xMargin;
	y2=(Int2)(y+(nTicks+nLet+1)*GrData->udv_font.LineHeight);
	AliPos=AlignPos+1;/*switch to 1-based value*/
	if ((ScalePosition==SCALE_POS_LEFT || 
		ScalePosition==SCALE_POS_BOTH)){
		/*scale on the left; must be located on the BioSeq line*/
		SetColor(GrData->udv_scale.ScaleColor);
		sprintf(szBuf,"%d",AliPos);
		x=rc->left-GrData->udv_font.cxChar-StringWidth(szBuf);
		MoveTo(x,y2);
		PaintString (szBuf);
	}

	/*Major ticks*/
	pos=StartLetter+1;	/*remember : StartLetter is zero based*/
	y=yMargin;
	xDecal=xMargin;
	y2=(Int2)(y+((nLet+1)*GrData->udv_font.LineHeight)-cyChar_4);
	y3=(Int2)(y+(nLet*GrData->udv_font.LineHeight)+cyChar_4);
	y4=(Int2)(y+(nLet*GrData->udv_font.LineHeight)+GrData->udv_font.cyChar);
	if ((ScalePosition==SCALE_POS_TOP || 
		ScalePosition==SCALE_POS_BOTH) && 
		ShowMajorTick){
		
		while(pos<StopLetter+1){
			if (/*pos!=StartLetter+1 && */AliPos!=ScaleMaxVal && 
				!(AliPos % LETTER_BLOCK_WIDTH)){
				SetColor(GrData->udv_scale.TickMajColor);
				x=xDecal;
				MoveTo(x,y3);
				LineTo(x,y2);
			}
			/*else if(ScalePosition==SCALE_POS_TOP &&
					pos==StartLetter+1){
				SetColor(GrData->udv_scale.TickMajColor);
				x=xDecal;
				MoveTo(x,y3);
				LineTo(x,y2);
			}*/
			/*Minor ticks; only shown if Major present*/
			if (ShowMMinorTick && 
				!(AliPos % (LETTER_BLOCK_WIDTH/2)) && 
					(AliPos % LETTER_BLOCK_WIDTH)){
				SetColor(GrData->udv_scale.TickMinColor);
				x=xDecal;
				MoveTo(x,y4-2);
				LineTo(x,y2-2);
			}

			xDecal+=ColWidth;

			/*each LETTER_BLOCK_WIDTH, add a blank column*/
			if (!(AliPos % LETTER_BLOCK_WIDTH)&&UseBlockDisp) xDecal+=ColWidth;
			AliPos++;pos++;
		}
	}
	Black();
}

/*******************************************************************************

  Function : UDV_calc_decalRC()
  
  Purpose : compute the y decal for ParaG RecT correct positionning on
  		the viewer panel
  
  Parameters :	ParaG; list of pgp
  				ScrollPos; actual position of the Vscroll bar
				SrollMax; max value of the Vscroll bar
				ScrollPage; size of a scroll page by page
				nTotLines; total lines for all the pgp
				LineHeight; height of a single line (pixel unit) 
				vnp; first ParaG to draw in the panel
				stop; value to be used to stop the draw
				
  Return value : y decal value (pixel unit)

*******************************************************************************/
static Int4  UDV_calc_decalRC(ValNodePtr ParaG,Int4 ScrollPos,
		Int4 ScrollMax,Int4 ScrollPage,Int4 nTotLines,Int2 LineHeight,
		ValNodePtr PNTR vnp,Int4Ptr stop)
{
ValNodePtr vnp2;
Int4 decal=0;
ParaGPtr pgp;
		
	/*find the ParaG where ScrollPos is Located*/
	for(*vnp=ParaG ; *vnp != NULL ; *vnp=(*vnp)->next){
		if ((*vnp)->data.ptrvalue){
			pgp=(ParaGPtr)(*vnp)->data.ptrvalue;
			if ((ScrollPos>=pgp->StartLine) &&
				(ScrollPos<=(pgp->StartLine+pgp->nLines))){
					break;
			}
		}		
	}	

	if (ScrollMax){
		*stop=MIN(nTotLines,ScrollPos+ScrollPage);
		for(vnp2=(*vnp) ; vnp2 != NULL ; vnp2=vnp2->next){
			pgp=(ParaGPtr)vnp2->data.ptrvalue;
			if (pgp->StartLine+pgp->nLines>=(*stop)){
				break;
			}
		}
	}
	else{
		*stop=nTotLines;
	}

	decal=ScrollPos*LineHeight-VIEWER_VERT_MARGIN;

	return(decal);
}

/*******************************************************************************

  Function : UDV_calc_RCdraw()
  
  Purpose : compute the RC of a ParaG to be drawn  
  
  Parameters :	rc_pgp; absolute position of a ParaG
  				decal; y decal to place rc_pgp within the viewer panel 
					(this value is calculated by UDV_calc_decalRC())
  
  Return value : position of the ParaG within the viewer panel

*******************************************************************************/
static RecT  UDV_calc_RCdraw(Int4 StartLine,Int4 nLines, RecT rcP,
		Int2 decal_gauche,Int4 decal_haut,Int2 LineHeight)
{
RecT rc;

	rc.top=(Int2)(StartLine*LineHeight+rcP.top-decal_haut);
	rc.bottom=(Int2)(rc.top+nLines*LineHeight);
	rc.left=decal_gauche;
	rc.right=rcP.right-VIEWER_HORZ_MARGIN;

	return(rc);
}

/*******************************************************************************

  Function : UDV_draw_empty_panel()
  
  Purpose : draw a emplty window
  
  Parameters :	rectangle to fill 
  				GrDate; graphical data
  
  Return value : none

*******************************************************************************/
static void  UDV_draw_empty_panel(RecT rc,UnDViewerGraphDataPtr GrData,
	Boolean ShowText)
{
	White();
	PaintRect(&rc);

	if (!ShowText) return;
	Black();
	SelectFont(GrData->udv_font.hFnt);

	MoveTo((Int2)(rc.left+5),(Int2)(rc.top+GrData->udv_font.cyChar));
	PaintString("Window");
	MoveTo((Int2)(rc.left+5),(Int2)(rc.top+2*GrData->udv_font.cyChar));
	PaintString("too small");
	MoveTo((Int2)(rc.left+5),(Int2)(rc.top+3*GrData->udv_font.cyChar));
	PaintString("to draw.");
}

/*******************************************************************************

  Function : UDV_draw_panel()
  
  Purpose :  draw the panel of UnD viewer
  
  Parameters :	handle of the panel
  
  Note : this function can be used as a callback by external software using
  		UnD-Viewer
  
  Return value : none

*******************************************************************************/
static void UDV_draw_panel(PaneL p,ViewerDialogDataPtr vdp)	
{
WindoW 				temport;
Int4 				stop,decal_haut,from_row,to_row,nTotRow,nLinesDraw;
ValNodePtr 			vnp,vnp2;
ParaGPtr 			pgp;
RecT 				rc,rcP,rcP_old;
Int2				decal_gauche;
SeqEditViewProcsPtr svpp=NULL;
ValNodePtr          vnp_seqinfo=NULL;
ValNodePtr          vnp_color = NULL;
MediaInfoPtr        MIPtr=NULL;

	temport=SavePort(ParentWindow(p));
	Select (p);
	/*restrict panel drawing area: 'add' little margins*/
	ObjectRect(p,&rcP);
	rcP_old=rcP;
	InsetRect(&rcP,4,4);
	ClipRect(&rcP);


	if (vdp == NULL) {
		ResetClip();	
		return;
	}

	if (vdp->ParaG == NULL) {
		UDV_draw_empty_panel(rcP,&vdp->udv_graph,FALSE);
		ResetClip();	
		return;
	}

	if (vdp->udv_graph.udv_panel.nBlockByLine<1) {
		UDV_draw_empty_panel(rcP,&vdp->udv_graph,TRUE);
		ResetClip();	
		return;
	}

	/*decal_haut=UDV_calc_decalRC(vdp->ParaG,vdp->udv_graph.udv_vscrl.ScrollPos,
		vdp->udv_graph.udv_vscrl.ScrollMax,vdp->udv_graph.udv_vscrl.ScrollPage,
		vdp->udv_graph.udv_panel.nTotLines,vdp->udv_graph.udv_font.LineHeight,
		&vnp,&stop);*/
	
	decal_haut=vdp->udv_graph.udv_vscrl.ScrollPos*vdp->udv_graph.udv_font.LineHeight-
			VIEWER_VERT_MARGIN;	

/*printf("%d..%d\n",rcP.top,rcP.bottom);*/
	if (vdp->udv_graph.bFirst){
		from_row=0;
		to_row=(rcP.bottom-rcP.top)/vdp->udv_graph.udv_font.LineHeight+1;
		vdp->udv_graph.bFirst=FALSE;
	}
	else{
		from_row=(updateRect.top-rcP.top)/vdp->udv_graph.udv_font.LineHeight-1;
		to_row=(updateRect.bottom-rcP.top)/vdp->udv_graph.udv_font.LineHeight+1;
	}
	if (from_row<0 && to_row<0) return;		
	if (from_row<0) from_row=0;
	if (to_row<0) to_row=0;
	from_row+=vdp->udv_graph.udv_vscrl.ScrollPos;
	to_row+=vdp->udv_graph.udv_vscrl.ScrollPos;
	if (from_row>vdp->udv_graph.udv_panel.nTotLines) 
		from_row=vdp->udv_graph.udv_panel.nTotLines;
	if (to_row>vdp->udv_graph.udv_panel.nTotLines) 
		to_row=vdp->udv_graph.udv_panel.nTotLines;
	nTotRow=to_row-from_row;

	/*find the ParaG where ScrollPos is Located*/
	vnp=vdp->ParaG;
	while(vnp){
		if (vnp->data.ptrvalue){
			pgp=(ParaGPtr)vnp->data.ptrvalue;
			if ((pgp->StartLine<=to_row)&&((pgp->StartLine+pgp->nLines)>=from_row)){
					break;
			}
		}
		vnp=vnp->next;
	}	

	/*retrieve the colours table for the bsp*/
	svpp = (SeqEditViewProcsPtr) GetAppProperty ("SeqEditDisplayForm");
	if (svpp) {
		vnp_seqinfo = svpp->seqinfo;
		while(vnp_seqinfo){
			MIPtr=(MediaInfoPtr)vnp_seqinfo->data.ptrvalue;
			if (MIPtr->entityID==vdp->bsp_i.bsp_entityID && 
				MIPtr->itemID==vdp->bsp_i.bsp_itemID){
				vnp_color=MIPtr->seq_color;
				break;
		}
		vnp_seqinfo=vnp_seqinfo->next;
		}
	}
	
	/*draw ParaG until out of panel*/
	SelectFont(vdp->udv_graph.udv_font.hFnt);

	/*3D border to resize the cxName field*/
	LtGray();
	
	MoveTo((Int2)(vdp->udv_graph.udv_panel.cxName-1),rcP.top);
	LineTo((Int2)(vdp->udv_graph.udv_panel.cxName-1),rcP.bottom);

	LtGray();
	MoveTo((Int2)(vdp->udv_graph.udv_panel.cxName+1),rcP.top);
	LineTo((Int2)(vdp->udv_graph.udv_panel.cxName+1),rcP.bottom);
	LtGray();
	MoveTo((Int2)(vdp->udv_graph.udv_panel.cxName+2),rcP.top);
	LineTo((Int2)(vdp->udv_graph.udv_panel.cxName+2),rcP.bottom);

	Black();
	MoveTo((Int2)(vdp->udv_graph.udv_panel.cxName+3),rcP.top);
	LineTo((Int2)(vdp->udv_graph.udv_panel.cxName+3),rcP.bottom);

	decal_gauche=vdp->udv_graph.udv_panel.cxName+rcP_old.left+
			vdp->udv_graph.udv_scale.cxLeftScale;
	nLinesDraw=0;		
	for(vnp2=vnp ; vnp2 != NULL ; vnp2=vnp2->next){
		if (vnp2->data.ptrvalue){
			pgp=(ParaGPtr)vnp2->data.ptrvalue;
			
			rc=UDV_calc_RCdraw(pgp->StartLine,pgp->nLines,rcP_old,
					decal_gauche,decal_haut,vdp->udv_graph.udv_font.LineHeight);

		/*Numerical scale*/
			if (vdp->udv_graph.udv_panel.ShowScale)
				UDV_Draw_scale(&vdp->udv_graph,
					vdp->udv_graph.udv_scale.ShowMajorTick,
					vdp->udv_graph.udv_scale.ShowMMinorTick,
					vdp->udv_graph.udv_scale.ScalePosition,
					pgp->StartLetter,pgp->StopLetter,&rc,
					vdp->udv_graph.udv_font.cxChar,
					vdp->bsp_i.bspLength,pgp->StartLetter,
					TRUE,
					vdp->udv_graph.udv_font.cxChar);
		/*Sequence*/
			if (vdp->bsp_i.SeqBuf)
				UDV_Draw_sequence(&vdp->udv_graph,vnp_color,
					vdp->udv_graph.udv_font.UseDefaultColorLetter,
					vdp->udv_graph.udv_font.LetterColor,
					pgp,&rc,vdp->bsp_i.StartBuf,pgp->StartLetter,
					vdp->bsp_i.SeqBuf);
		/*Sequence name*/
			UDV_Draw_sequence_name(&vdp->udv_graph,&vdp->bsp_i,&rc);
		/*Features*/
			if (vdp->udv_graph.udv_panel.ShowFeatures){
				UDV_Draw_features(&vdp->udv_graph,&vdp->bsp_i,
					vdp->udv_graph.udv_font.UseDefaultColorLetter,
					vdp->udv_graph.udv_font.LetterColor,pgp,&rc,
					vdp->udv_graph.pClr,vdp->Item_select,vdp->Old_Item_select);
				/*Feature names
					UDV_Draw_features_label(&vdp->udv_graph,&vdp->bsp_i,pgp,
						&rc,
						vdp->udv_graph.pClr,vdp->Item_select,
						vdp->Old_Item_select);*/
			}
			nLinesDraw+=pgp->nLines;
			if (nLinesDraw>nTotRow) break;
			/*if (pgp->StartLine > stop) break;*/
		}		
	}	

	ResetClip();	
	RestorePort(temport);
}

/*******************************************************************************

  Function : UDV_draw_viewer()
  
  Purpose :  draw the panel of UnD viewer
  
  Parameters :	handle of the panel
  
  Note : this function is called when clickfeat; this is not a callback
  
  Return value : none

*******************************************************************************/
NLM_EXTERN void  UDV_draw_viewer (PaneL p)
{
WindoW 	w;
ViewerMainPtr		vmp;

	/*get the parent Window of the Panel*/
	w=ParentWindow(p);
	if (w==NULL) return;
	
	/*get some usefull data...*/
	vmp = (ViewerMainPtr) GetObjectExtra (w);

	if (vmp==NULL) return;

	UDV_draw_panel(p,vmp->vdp);
}


/*******************************************************************************

  Function : UDV_Logo_onDraw()
  
  Purpose : draw the software logo
  
  Parameters :	handle of the panel
  
  Return value : none

*******************************************************************************/
NLM_EXTERN void UDV_Logo_onDraw (PaneL p)
{
WindoW 				w,temport;
RecT 				/*rc,*/rcP;
/*ViewerMainPtr		vmp;*/
UDVLogoDataPtr		ldp;

	/*get the parent Window of the Panel*/
	w=ParentWindow(p);
	if (w==NULL) return;
	
	/*get some usefull data...*/
	/*vmp = (ViewerMainPtr) GetObjectExtra (w);*/
	ldp = (UDVLogoDataPtr) GetAppProperty ("UDVLogoData");

	if (ldp==NULL) return;
	
	temport=SavePort(ParentWindow(p));
	Select (p);
	/*restrict panel drawing area: 'add' little margins*/
	ObjectRect(p,&rcP);
	InsetRect(&rcP,4,4);
	ClipRect(&rcP);

	White();
	PaintRect(&rcP);
	if (ldp->f1 && ldp->f2 && ldp->f3) draw_logo(rcP,ldp->f1,ldp->f2,ldp->f3,
			ldp->szTitle,ldp->szDesc);
	ResetClip();	

	/*3D border*/
	Black();
	FrameRect(&rcP);
	InsetRect(&rcP,-1,-1);
	LtGray();
	FrameRect(&rcP);
	InsetRect(&rcP,-1,-1);
	FrameRect(&rcP);
	Black();
	InsetRect(&rcP,-2,-2);
	FrameRect(&rcP);
	RestorePort(temport);
}

/*******************************************************************************

  Function : UDV_InfoPanelDrawProc()
  
  Purpose : draw the content of the Info panel
  
  Parameters :	handle of the panel
  
  Return value : none

*******************************************************************************/
NLM_EXTERN void UDV_InfoPanelDrawProc(PaneL p)
{
RecT 				rc;
CharPtr 			szFeatName;
WindoW 				w,temport;
ViewerDialogDataPtr vdp;
ViewerMainPtr		vmp;

	w=ParentWindow(p);
	
	if (w==NULL) return;

	/*get some usefull data...*/
	vmp = (ViewerMainPtr) GetObjectExtra (w);

	if (vmp==NULL) return;
	temport=SavePort(w);

	Select(p);
	vdp=vmp->vdp;

	if (vdp) SelectFont(vdp->udv_graph.udv_font.hFnt);
	else return;
	
	szFeatName=(CharPtr) GetObjectExtra (p);

	ObjectRect(p,&rc);
	LtGray();
	PaintRect(&rc);
	Black();
	FrameRect(&rc);
	if (szFeatName){
		MoveTo((Int2)(rc.left+10),
				(Int2)(rc.bottom-vdp->udv_graph.udv_font.cyChar/4));
		Black();
		PaintString(szFeatName);
	}
	Black();
	RestorePort(temport);
}
