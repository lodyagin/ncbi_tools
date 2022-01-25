/*  $Id: pgppop.c,v 6.24 1999/09/16 18:52:26 durand Exp $
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
* File Name:  pgppop.c
*
* Author:  Patrick Durand
*
* Version Creation Date:   05/03/99
*
* $Revision: 6.24 $
*
* File Description: 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: pgppop.c,v $
* Revision 6.24  1999/09/16 18:52:26  durand
* redesign the PopSet viewer toolbar
*
* Revision 6.23  1999/09/07 13:41:55  durand
* update Entrez links for PopSet Viewer
*
* Revision 6.22  1999/09/02 19:08:56  chappey
* fixes in PrintSeqAlignCallback
*
* Revision 6.21  1999/09/01 21:04:20  durand
* call SeqAlignSetFree after PairSeqAlign2MultiSeqAlign
*
* Revision 6.20  1999/08/31 13:51:30  durand
* add PubMed link for PopSet Viewer
*
* Revision 6.19  1999/08/31 12:35:46  wheelan
* fixed bug in DDV_GetSeqAlign
*
* $Log: pgppop.c,v $
* Revision 6.24  1999/09/16 18:52:26  durand
* redesign the PopSet viewer toolbar
*
* Revision 6.23  1999/09/07 13:41:55  durand
* update Entrez links for PopSet Viewer
*
* Revision 6.22  1999/09/02 19:08:56  chappey
* fixes in PrintSeqAlignCallback
*
* Revision 6.21  1999/09/01 21:04:20  durand
* call SeqAlignSetFree after PairSeqAlign2MultiSeqAlign
*
* Revision 6.20  1999/08/31 13:51:30  durand
* add PubMed link for PopSet Viewer
*
* Revision 6.18  1999/08/30 21:24:23  durand
* display Unpublished study when a PopSet entry doesn't have a title
*
* Revision 6.17  1999/08/30 20:29:56  durand
* make PrintSeqAlignCallback estern function
*
* Revision 6.16  1999/08/30 18:20:11  durand
* update SeqAlignToBS
*
* Revision 6.15  1999/08/30 17:57:05  sirotkin
* changed mydata in PrintSeqAlignCallback
*
* Revision 6.14  1999/08/30 14:18:08  durand
* use ByteStore to format the SeqAlign output
*
* Revision 6.12  1999/08/27 14:54:04  durand
* fix memory leaks in PopSet Viewer
*
* Revision 6.11  1999/08/16 18:56:59  lewisg
* made DDV_GetSeqAlign extern and added prototype to header
*
* Revision 6.10  1999/08/08 15:54:35  chappey
* made DDV_GetSeqAlign static as there is no prototype
*
* Revision 6.9  1999/08/07 16:57:48  sicotte
* fixed compiler warnings
*
* Revision 6.8  1999/08/07 16:53:13  sicotte
* added includes sqnutils.h and alignval.h for new code
*
* Revision 6.7  1999/08/06 21:43:17  chappey
* SeqAlignToBS new function to save in ByteStore structure the text output of the SeqAlign(s) packaged in a SeqEntry
*
* Revision 6.6  1999/07/22 13:23:13  durand
* made DDV_SearchAli external function
*
* Revision 6.5  1999/07/21 21:52:12  durand
* add some functions to display a summary for a PopSet entry
*
* Revision 6.4  1999/07/19 21:16:01  durand
* add DDV_ResetParaGSeqAlignCoord to reset the seqalign coord in the display data structures of DDV
*
* Revision 6.3  1999/07/15 18:20:51  durand
* add display options to support BLAST outputs
*
* Revision 6.2  1999/07/13 20:46:56  durand
* comment out call to PairSeqAlign2MultiSeqAlign to avoid DDV compiling problems
*
* Revision 6.1  1999/07/09 13:59:58  durand
* move pgppop from desktop to api
*
* Revision 6.2  1999/07/07 19:12:26  durand
* fix a tiny bug in DDV_GetSequenceFromParaG
*
* Revision 6.1  1999/07/06 20:18:07  kans
* initial public checkin
*
* Revision 1.35  1999/07/06 18:54:27  durand
* add new features for the display of PopSet viewer
*
* Revision 1.34  1999/07/02 13:22:17  durand
* fix bugs for the display of minus strand sequences
*
* Revision 1.32  1999/06/29 16:48:10  shavirin
* Changed definition of function DDV_ShowSeqAlign()
*
* Revision 1.31  1999/06/28 22:07:20  durand
* add loader functions and clean the code with Lint and Purify
*
* Revision 1.30  1999/06/24 20:49:41  shavirin
* Added new function DDV_ShowSeqAlign().
*
* Revision 1.29  1999/06/23 18:11:17  durand
* fix a variable initialization problem under NT
*
* Revision 1.28  1999/06/23 17:24:23  durand
* use a binary encoding to manage the display styles
*
* Revision 1.27  1999/06/21 18:37:56  durand
* update DDV_DisplayDefaultAlign to produce full text output
*
* Revision 1.26  1999/06/19 18:36:13  durand
* new display procedure
*
* Revision 1.25  1999/06/16 13:10:48  durand
* update/add functions for Vibrant DDV
*
* Revision 1.24  1999/06/14 23:49:43  durand
* add function for Vibrant DDV
*
* Revision 1.23  1999/06/11 22:33:01  durand
* add new functions for Vibrant DDV
*
* Revision 1.22  1999/06/11 17:59:40  durand
* popset viewer uses more code from UDV
*
* Revision 1.21  1999/06/09 21:35:30  durand
* add constructors/destructors for BspInfo struct as well as read seq function
*
*
*
*
* ==========================================================================
*/
#include <stdio.h>
#include <ncbi.h>

/*
#include <biocol.h>
*/
#include <blocks.h>
#include <pgppop.h>
#include <tofasta.h>
#include <udvdef.h>
#include <udvseq.h>
#include <salstruc.h>

#include <sqnutils.h>
#include <alignval.h>

#define DDV_COLOR_MAX 33 
 static Uint1 DDV_PaletteRGB[DDV_COLOR_MAX][3] =
 {
     /* Red  Grn Blue     Color     ColorID             Old Hex  New Hex */
   255, 255, 255, /* default     0                 0xffffff - */ 
   255,  00, 153, /* hotpink     1                 0xff1493 0xff0099 */
   255,   0, 255, /* magenta     2                 0xff00ff - */
   153,  51, 255, /* purple      3                 0x9b30ff 0x9933ff */
     0,   0, 255, /* blue        4                 0x0000ff - */
    00, 153, 255, /* sky         5                 0x1e90ff 0x0099ff */
     0, 255, 255, /* cyan        6                 0x00ffff - */
     0, 255, 153, /* sea         7                 0x00ff8f 0x00ff99 */
     0, 255,   0, /* green       8                 0x00ff00 - */
   255, 255,   0, /* yellow      9                 0xffff00 - */
   255, 204,   0, /* gold       10                 0xffa500 0xffcc00 */
   255, 102,   0, /* orange     11                 0xff4500 0xff6600 */
   255,   0,   0, /* red        12                 0xff0000 - */
   255, 153, 153, /* pink       13                 0xff7256 0xff9999 */
   255, 204, 204, /* pinktint   14                 0xffaeb9 0xffcccc */
   255, 255, 255, /* white      15                0xffffff  -*/
     0,   0,   0, /* black      16                0x000000 - */
   204, 204, 255, /* bluetint   17                0xb0e2ff - 0xccccff*/
   153, 255, 153, /* greentint  18                0x9aff9a - 0x99ff99 */
   255, 255, 153, /* yellowtint 19                0xffec8b - 0xffff99 */
   102, 102, 102, /* gray       20                0x7d7d7d   0x666666 */
   153, 102,  51, /* brown      21                0x8b5742   01x996633 */
   255,   0,   0, /* perm colors 22: red          0xff0000 - */
     0, 255,   0, /* perm colors 23: green        0x00ff00 - */
   255,   0, 255, /* perm colors 24: magenta      0xff00ff - */
    00, 153, 255, /* perm colors 25: sky          0x1e90ff 0x0099ff */
   153,	 51, 255, /* perm colors 26: purple        0x9b30ff 0x9933ff */
     0, 255,   0, /* SS colors 27: helix, green   0x00ff00 - */
   255, 204,   0, /* SS colors 28: strand, gold   0xffa500 0xffcc00 */
   255, 102,   0, /* SS colors 29: turn, orange   0xff4500 0xff6600 */
     0, 255, 255, /* SS colors 30: coil, cyan     ox00ffff - */
   255, 255,   0, /* highlight colors 31: yellow  0xffff00 - */
     0,   0,   0 /* background colors 32: black  0x000000 - */
 };

static CharPtr DDV_PaletteRGB_NAV[DDV_COLOR_MAX] =
 {
     /* Red  Grn Blue     Color     ColorID             Old Hex  New Hex */
   "ffffff", /* default     0                 0xffffff - */ 
   "ff0099", /* hotpink     1                 0xff1493 0xff0099 */
   "ff00ff", /* magenta     2                 0xff00ff - */
   "9933ff", /* purple      3                 0x9b30ff 0x9933ff */
   "0000ff", /* blue        4                 0x0000ff - */
   "0099ff", /* sky         5                 0x1e90ff 0x0099ff */
   "00ffff", /* cyan        6                 0x00ffff - */
   "00ff99", /* sea         7                 0x00ff8f 0x00ff99 */
   "00ff00", /* green       8                 0x00ff00 - */
   "ffff00", /* yellow      9                 0xffff00 - */
   "ffcc00", /* gold       10                 0xffa500 0xffcc00 */
   "ff6600", /* orange     11                 0xff4500 0xff6600 */
   "ff0000", /* red        12                 0xff0000 - */
   "ff9999", /* pink       13                 0xff7256 0xff9999 */
   "ffcccc", /* pinktint   14                 0xffaeb9 0xffcccc */
   "ffffff", /* white      15                0xffffff  -*/
   "000000", /* black      16                0x000000 - */
   "ccccff", /* bluetint   17                0xb0e2ff - 0xccccff*/
   "99ff99", /* greentint  18                0x9aff9a - 0x99ff99 */
   "ffff99", /* yellowtint 19                0xffec8b - 0xffff99 */
   "666666", /* gray       20                0x7d7d7d   0x666666 */
   "996633", /* brown      21                0x8b5742   01x996633 */
   "ff0000", /* perm colors 22: red          0xff0000 - */
   "00ff00", /* perm colors 23: green        0x00ff00 - */
   "ff00ff", /* perm colors 24: magenta      0xff00ff - */
   "0099ff", /* perm colors 25: sky          0x1e90ff 0x0099ff */
   "9933ff", /* perm colors 26: purple        0x9b30ff 0x9933ff */
   "00ff00", /* SS colors 27: helix, green   0x00ff00 - */
   "ffcc00", /* SS colors 28: strand, gold   0xffa500 0xffcc00 */
   "ff6600", /* SS colors 29: turn, orange   0xff4500 0xff6600 */
   "00ffff", /* SS colors 30: coil, cyan     ox00ffff - */
   "ffff00", /* highlight colors 31: yellow  0xffff00 - */
   "000000" /* background colors 32: black  0x000000 - */
 };

#define DDVCOL_DEFAULT 0
#define DDVCOL_HOTPINK 1
#define DDVCOL_MAGENTA 2
#define DDVCOL_PURPLE 3
#define DDVCOL_BLUE 4
#define DDVCOL_SKY 5
#define DDVCOL_CYAN 6
#define DDVCOL_SEA 7
#define DDVCOL_GREEN 8
#define DDVCOL_YELLOW 9
#define DDVCOL_GOLD 10
#define DDVCOL_ORANGE 11
#define DDVCOL_RED 12
#define DDVCOL_PINK 13
#define DDVCOL_PINKTINT 14
#define DDVCOL_WHITE 15
#define DDVCOL_BLACK 16
#define DDVCOL_BLUETINT 17
#define DDVCOL_GREENTINT 18
#define DDVCOL_YELLOWTINT 19
#define DDVCOL_GRAY 20
#define DDVCOL_BROWN 21

#define DDVCOL_HIGHLIGHT 31
#define DDVCOL_BACKGROUND 32

typedef struct getbsppop{
	Int4 nCompt;
	FILE *fp;
	}GetBspPop, PNTR GetBspPopPtr;

/*Entrez and PopSet Viewer scripts*/
#define szEntrezScript "http://www.ncbi.nlm.nih.gov/entrez/utils/qmap.cgi?uid=%d&form=6&db=%s&Dopt=g"
#define szSeqNameStatus "<a href=\"%s\" onMouseOut=\"window.status=''\" \nonMouseOver=\"window.status='%s';return true\">"
#define WWWDDV_script "wwwddv.cgi?gi=%d&from=%d&to=%d&disp=%u"
#define WWWDDV_script2 "wwwddv.cgi?gis=%s&disp=%u"
#define WWWDDV_save_script1 "wwwddv.cgi"

static void DDV_ToolsScript(SeqEntryPtr sep, Uint4 disp_format,FILE *fp, 
	Boolean ForToolbar);

/*******************************************************************************

  Function : ConvertToLower()
  
  Purpose : convert an upper case string to a lower case one
  
  Parameters : header of the ValNodeList
  
  Note : this function is a copy of to_lower() (ncbi/tools). Because pgppop.c
         is in ncbi/api, I cannot use the code of salutil.
		 
  Return value : none 

*******************************************************************************/
static CharPtr ConvertToLower (CharPtr str)
{
  CharPtr tmp;
 
  if (str==NULL)
     return NULL;
  tmp = str;
  while (*tmp!='\n' && *tmp!='\0') {
     *tmp = TO_LOWER(*tmp);
     tmp++;
  }
  return str;
}

/*******************************************************************************

  Function : DDV_DeleteTxtList()
  
  Purpose : delete a list of sequence content in a ParaG
  
  Parameters : header of the ValNodeList
  
  Return value : none 

*******************************************************************************/
extern void DDV_DeleteTxtList(ValNodePtr PNTR vnp)
{
ValNodePtr 	vnp2;
MsaTxtDispPtr mtdp;

	for(vnp2=(*vnp) ; vnp2!=NULL ; vnp2=vnp2->next){
		mtdp=(MsaTxtDispPtr)vnp2->data.ptrvalue;
		if (mtdp) MemFree(mtdp);
	}
	ValNodeFree(*vnp);
	*vnp=NULL;
}

/*******************************************************************************

  Function : DDV_DeleteParaGList()
  
  Purpose : delete a list of ParaG
  
  Parameters : header of the ValNodeList
  
  Return value : none 

*******************************************************************************/
extern void DDV_DeleteParaGList(ValNodePtr PNTR vnp)
{
ValNodePtr 	vnp2;
ParaGPtr pgp;

	for(vnp2=(*vnp) ; vnp2!=NULL ; vnp2=vnp2->next){
		pgp=(ParaGPtr)vnp2->data.ptrvalue;
		if (pgp){
			if (pgp->ptxtList) DDV_DeleteTxtList(&pgp->ptxtList);
			MemFree(pgp);
		}
	}
	ValNodeFree(*vnp);
	*vnp=NULL;
}

/*******************************************************************************

  Function : DDV_DeleteDisplayList()
  
  Purpose : delete the complete lists used for the display
  
  Parameters : main data block used for the display
  
  Return value : none 

*******************************************************************************/
extern void DDV_DeleteDisplayList(MsaParaGPopListPtr mpplp)
{
Int4 i;

	/*ValNode with ParaG*/
	if (mpplp->TableHead){
		for(i=0;i<mpplp->nBsp;i++){
			if (mpplp->TableHead[i]){
				DDV_DeleteParaGList(&(mpplp->TableHead[i]));
			}
		}
		MemFree(mpplp->TableHead);
	}
}


/*******************************************************************************

  Function : DDV_ResetParaGSeqAlignCoord()
  
  Purpose : this function reset the SeqAlign coords in the display data, ParaG
  
  Parameters : mpplp;main data block used for the display
  				LineSize; max size of each ParaG. Must be the same value as
				the one passed in DDV_PopDisplay().
  
  Return value : none 

*******************************************************************************/
extern void DDV_ResetParaGSeqAlignCoord(MsaParaGPopListPtr mpplp,Int2 LineSize)
{
Int4 i,decal;
Int2 OldLineSize;
ValNodePtr vnp;
ParaGPtr pgp;

	for(i=0;i<mpplp->nBsp;i++){
		vnp=mpplp->TableHead[i];
		decal=0;
		while(vnp){
			pgp=(ParaGPtr)(vnp->data.ptrvalue);
			if (pgp){
				/*the last ParaG of a list may be smaller than LineSize*/
				/*So, I get the current ParaG size and reuse it two lines below*/
				OldLineSize=pgp->StopLetter-pgp->StartLetter+1;
				pgp->StartLetter=decal;
				pgp->StopLetter=pgp->StartLetter+_min_(LineSize,OldLineSize)-1;
			}
			decal+=(Int4)LineSize;
			vnp=vnp->next;
		}
	}
}

/*******************************************************************************

  Function : DDV_PopDisplay2()
  
  Purpose : (re)populate ParaGs from the indexed seqalign
  
  Parameters : 	sabp;header of the indexed seqalign
  				mpplp; where ParaGs list are stored
  				nBSP;sabp contains 'nBsp' sequences
				MasterScaleStyle;scale stype for the master sequence
				SlaveScaleStyle;scale stype for the slave sequence
				cbl; number or char/line for a each ParaG
				nBlockToPop; number of blocks to analyse
				ali_from; create ParaG list given a from-to in SEQALIGN COORD
				ali_to;create ParaG list given a from-to in SEQALIGN COORD
				
  Return value : TRUE if success

*******************************************************************************/
extern Boolean DDV_PopDisplay2(SABlockPtr sabp,MsaParaGPopListPtr mpplp,
			Int4 nBSP,Uint1 MasterScaleStyle,Boolean ShowTick,
			Int2 cbl,Int4 nBlockToPop,Int4 ali_from,Int4 ali_to)
{
SABlockPtr 	sabp2;	/*used to scan "SeqALign" Index blocks*/
SegmentPtr 	segp;	/*used to scan each Bsp segment in a block*/
Int4 		startPgpPos2;/*used to set the from value of each ParaG*/
Int4 		startcopy;/*used to copy the content of a segment in a ParaG*/
Int4 		stopcopy;/*idem*/
Int4 		segsize;/*idem*/
BoolPtr 	bFirst;/*idem*/
Int4 		nCharToCopy;/*idem*/
Int4 		countSegMax=0;/*used to count the number of sequence to display*/
Int4 		countSeg=0;/*idem*/
Int4		countBlock=0;
Int4 		AbsLengthAli=0;/*used to compute the size of the align to display*/
Int4 		i;/*just a little counter*/
ValNodePtr 	vnp,vnp2;/*used to delete additional ParaG nodes, if needed;
				used for REpopulation*/
ParaGPtr pgp1,pgp2;/*used to create new ParaGs*/
MsaTxtDispPtr mtdp;/*used to store the content of a line in a ParaG*/
ValNodePtr PNTR TableDesc=NULL;/*used to populate mpplp->TableHead*/
ValNodePtr PNTR TableText=NULL;/*used to populate ParaG with Text*/
Int4		gFrom,gTo,gDecal_L,gDecal_R;/*used to prepare a display from a range*/

/*
#ifdef DEBUG_DDV
	blk_PrintSABP(sabp);
#endif
*/
	/*some initialisations*/
	if (mpplp->TableHead==NULL){/*First time*/
		mpplp->TableHead=(ValNodePtr PNTR)MemNew(nBSP*sizeof(ValNodePtr));
		if (mpplp->TableHead==NULL) return(FALSE);
		MemSet(mpplp->TableHead,0,nBSP*sizeof(ValNodePtr));
	}
	else{/*reuse the structures*/
		/*same number of BSPs ?*/
		if (nBSP>mpplp->nBsp){/*increase the table*/
			ValNodePtr PNTR tmp;
			tmp=(ValNodePtr PNTR)MemNew(nBSP*sizeof(ValNodePtr));
			if (tmp==NULL) return(FALSE);
			MemSet(tmp,0,nBSP*sizeof(ValNodePtr));
			MemCopy(tmp,mpplp->TableHead,mpplp->nBsp*sizeof(ValNodePtr));
			MemFree(mpplp->TableHead);
			mpplp->TableHead=tmp;
		}
	}
	/*used to populate mpplp->TableHead*/
	TableDesc=(ValNodePtr PNTR)MemNew(nBSP*sizeof(ValNodePtr));
	if (TableDesc==NULL) return(FALSE);
	/*copy the ValNode heads in TableDesc*/
	for(i=0;i<nBSP;i++){
		TableDesc[i]=mpplp->TableHead[i];
	}

	bFirst=(BoolPtr)MemNew(nBSP*sizeof(BoolPtr));
	if (bFirst==NULL) return(FALSE);
	MemSet(bFirst,1,nBSP*sizeof(BoolPtr));

	/*used to populate ParaG with Text*/
	TableText=(ValNodePtr PNTR)MemNew(nBSP*sizeof(ValNodePtr));
	if (TableText==NULL) return(FALSE);
	MemSet(TableText,0,nBSP*sizeof(ValNodePtr));

	/*scan each block of the "SeqALign" Index*/
	for(sabp2=sabp ; sabp2!=NULL ; sabp2=sabp2->next){
		if (!sabp2->visible) continue;/*used only visible block*/
		if(ali_from!=(-1) && ali_to!=(-1)){/*use a range*/
			if (ali_from>sabp2->to) continue;/*not yet in range*/
			gFrom=_max_(ali_from,sabp2->from);
			gTo=_min_(ali_to,sabp2->to);
			gDecal_L=gFrom-sabp2->from;
			gDecal_R=sabp2->to-gTo;
		}
		else{/*display the full alignment*/
			gFrom=sabp2->from;
			gTo=sabp2->to;
			gDecal_L=gDecal_R=0;
		}
		/*used to compute the total align length to display*/
		AbsLengthAli+=(gTo-gFrom+1);
		/*used to compute the total number of bsp to display*/
		countSeg=0;
		countBlock++;
		/*scan each Bsp segment in a block*/
		for(segp=sabp2->segp_head ; segp!=NULL ; segp=segp->next){
			if (!segp->visible){/*hide sequence*/
				/*TableDesc[segp->BspID-1]=NULL;*/
				continue;
			}

			countSeg++;
			/*allocation ?*/
			if (TableDesc[segp->BspID-1]==NULL){
				mpplp->TableHead[segp->BspID-1]=ValNodeNew(NULL);
				TableDesc[segp->BspID-1]=
					mpplp->TableHead[segp->BspID-1];
			}
			/*is this segment a gap ?*/
			if (segp->gap==0) {/*no*/
				if(segp->strand!=Seq_strand_minus){
					startcopy=segp->from+gDecal_L;/*use Bsp coordinates*/
					stopcopy=segp->to+1-gDecal_R;/*use Bsp coordinates*/
				}
				else{
					startcopy=segp->to-gDecal_L;/*use Bsp coordinates*/
					stopcopy=segp->from-1+gDecal_R;/*use Bsp coordinates*/
				}
			}
			else {/*yes*/
				startcopy=gFrom;/*use Align coordinates*/
				stopcopy=gTo+1;/*use Align coordinates*/
			}
			startPgpPos2=gFrom;
			/*retrieve a ParaG*/
			pgp1=(ParaGPtr)TableDesc[segp->BspID-1]->data.ptrvalue;
			/*copy segment data to ParaG*/
			while(((segp->strand!=Seq_strand_minus || segp->gap!=0) ? 
					(startcopy<stopcopy) :
					(startcopy>stopcopy))){
				if (bFirst[segp->BspID-1]){
					bFirst[segp->BspID-1]=FALSE;
					/*is this ParaG already exist*/
					if (pgp1==NULL){/*no: create it*/
						pgp1=(ParaGPtr)MemNew(sizeof(ParaG));
						if (pgp1==NULL) return(FALSE);
						MemSet(pgp1,0,sizeof(ParaG));
						TableDesc[segp->BspID-1]->data.ptrvalue=(Pointer)pgp1;
					}

					/*is the Text fields exist ?*/
					if(pgp1->ptxtList==NULL){
						pgp1->ptxtList=ValNodeNew(NULL);
						if (!pgp1->ptxtList) return(FALSE);
					}
					TableText[segp->BspID-1]=pgp1->ptxtList;

					pgp1->StartLetter=pgp1->StopLetter=startPgpPos2;
					if (segp->BspID-1==0){/*master sequence*/
						pgp1->ScaleStyle=MasterScaleStyle;
						if (MasterScaleStyle&SCALE_POS_LEFT){ 
							pgp1->nLines=1;
						}
						else 
							if (MasterScaleStyle&SCALE_POS_TOP) {
								pgp1->nLines=2;
							}
						else 
							if (MasterScaleStyle&SCALE_POS_BOTH) {
								pgp1->nLines=2;
							}
						else 
							if (MasterScaleStyle&SCALE_POS_NONE) {
								pgp1->nLines=0;
							}
						if (ShowTick) pgp1->nLines++;
					}
					else{
						pgp1->nLines=1;
					}
					pgp1->sip=segp->sip;
				}
				else{
					if (!TableText[segp->BspID-1]->next){
						TableText[segp->BspID-1]->next=
							ValNodeAdd(&(TableText[segp->BspID-1]));
						if (!(TableText[segp->BspID-1]->next)) return(FALSE);
					}
					TableText[segp->BspID-1]=TableText[segp->BspID-1]->next;
				}
				if (segp->strand!=Seq_strand_minus || segp->gap!=0) 
					segsize=stopcopy-startcopy;
				else
					segsize=startcopy-stopcopy;
				
				nCharToCopy=MIN(segsize,cbl-(pgp1->StopLetter-pgp1->StartLetter));
				pgp1->StopLetter+=nCharToCopy;
				/*init a new block of info*/
				mtdp=(MsaTxtDispPtr)TableText[segp->BspID-1]->data.ptrvalue;
				if(!mtdp){
					mtdp=(MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
					if (mtdp==NULL) return(FALSE);
					TableText[segp->BspID-1]->data.ptrvalue=(Pointer)mtdp;
				}
				if (segp->strand!=Seq_strand_minus || segp->gap!=0){
					mtdp->from=startcopy;
					mtdp->to=startcopy+nCharToCopy-1;
				}
				else{
					mtdp->from=startcopy-nCharToCopy+1;
					mtdp->to=startcopy;
				}
				mtdp->SegID=sabp2->SegID;
				mtdp->BspID=segp->BspID;
				mtdp->IsGap=(segp->gap==0 ? FALSE : TRUE);
				mtdp->strand=segp->strand;
				if (mtdp->IsGap) mtdp->TextStyle=MSA_TXT_STYLE_GAP;
				else mtdp->TextStyle=MSA_TXT_STYLE_SEQ;
				
				/*add a next node if needed*/
				/*if (!TableText[segp->BspID-1]->next){
					TableText[segp->BspID-1]->next=
						ValNodeAdd(&(TableText[segp->BspID-1]));
					if (!(TableText[segp->BspID-1]->next)) return(FALSE);
				}
				TableText[segp->BspID-1]=TableText[segp->BspID-1]->next;*/
				pgp2=pgp1;
				vnp=TableDesc[segp->BspID-1];
				/*switch to a next ParaG*/
				if (pgp1->StopLetter-pgp1->StartLetter>=cbl){
					pgp1->StopLetter--;/*swith to zero-based*/
					if (!TableDesc[segp->BspID-1]->next){
						TableDesc[segp->BspID-1]->next=
							ValNodeAdd(&(TableDesc[segp->BspID-1]));
						if (!(TableDesc[segp->BspID-1]->next)) return(FALSE);
					}
					TableDesc[segp->BspID-1]=
						TableDesc[segp->BspID-1]->next;
					vnp2=TableDesc[segp->BspID-1];
					/*retrieve a ParaG*/
					pgp1=(ParaGPtr)TableDesc[segp->BspID-1]->data.ptrvalue;
					bFirst[segp->BspID-1]=TRUE;
					startPgpPos2+=nCharToCopy;
					/*delete additional ParaG->ptxtList if needed*/
					if (TableText[segp->BspID-1]->next){
						DDV_DeleteTxtList(&(TableText[segp->BspID-1]->next));
					}
				}
				else{
					vnp2=TableDesc[segp->BspID-1]->next;
				}
				/*used to compute the total number of bsp to display*/
				if (segp->strand!=Seq_strand_minus || segp->gap!=0)
					startcopy+=nCharToCopy;
				else
					startcopy-=nCharToCopy;
			}/*end while(segment)*/
			/*just in case: switch ParaG-> to zero based*/
			if (countBlock==nBlockToPop){
				if (pgp2 && ((pgp2->StopLetter-pgp2->StartLetter+1)>=AbsLengthAli)) {
					pgp2->StopLetter=pgp2->StartLetter+AbsLengthAli-1;
				}
				/*erase additional ParaG if needed*/
				if (vnp2){
					DDV_DeleteParaGList(&vnp2);
				}
				vnp->next=NULL;
				if (TableText[segp->BspID-1]->next){
					DDV_DeleteTxtList(&TableText[segp->BspID-1]->next);
				}
			}
		}/*end for (segment)*/
		if (countSeg>countSegMax) countSegMax=countSeg;
		if (ali_from!=(-1) && ali_to!=(-1)){/*out of range ?*/
			if (ali_to<=sabp2->to) break;
		}
	}/*end for(block)*/

	/*size of the align to display*/
	mpplp->LengthAli=AbsLengthAli;
	mpplp->nBsp=countSeg;
	MemFree(bFirst);	
	MemFree(TableDesc);
	MemFree(TableText);
	return(TRUE);
}
	/*small version of DDV_PopDisplay2*/
extern Boolean DDV_PopDisplay(SABlockPtr sabp,MsaParaGPopListPtr mpplp,
			Int4 nBSP,Uint1 MasterScaleStyle,Boolean ShowTick,
			Int2 cbl,Int4 nBlockToPop)
{
	return(DDV_PopDisplay2(sabp, mpplp,nBSP, MasterScaleStyle, ShowTick,
			 cbl, nBlockToPop,-1,-1));
}

/*******************************************************************************

  Function : GetClrFromLetter()
  
  Purpose : get the colour of a letter (na or aa)
  
  Parameters : 	residue; the letter
  				IsAA; True when prot
				szColor; HTML encoded color
				
  Return value : - (see szColor)

*******************************************************************************/
static void GetClrFromLetter(Char residue,Boolean IsAA,CharPtr szColor)
{

	/*AA colour code system*/
	if(IsAA){/*not yet implemented*/
		StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_BLACK]);
	}
	else{/*NA colour code*/
		switch(residue){
			case '-':
				StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_PINK]);
				break;
			case '.':
				StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_BACKGROUND]);
				break;
			case 'A':case 'a':
				StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_GREEN]);
				break;
			case 'G':case 'g':
				StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_BLACK]);
				break;
			case 'C':case 'c':
				StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_BLUE]);
				break;
			case 'T':case 't':
				StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_RED]);
				break;
			case 'N': case 'n':
				StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_HOTPINK]);
				break;
			default:
				StringCpy(szColor,DDV_PaletteRGB_NAV[DDVCOL_BLACK]);
				break;
		}		
	}
}

/*******************************************************************************

  Function : ConcatStr()
  
  Purpose : merge to CharPtr; realloc if needed
  
  Parameters : 	src; CharPtr source
  				newStr; CharPtr to add to src

  Return value : the merging result of src and newStr

*******************************************************************************/
static CharPtr ConcatStr(CharPtr src,CharPtr newStr)
{
CharPtr retScr;

	if(src==NULL){
		retScr=StringSave(newStr);
	}
	else{
		retScr=(CharPtr)MemNew((StringLen(src)+StringLen(newStr)+1)*sizeof(Char));
		if (!retScr) return(NULL);
		StringCpy(retScr,src);
		StringCat(retScr,newStr);
		MemFree(src);
	}
	
	return(retScr);
}


/*******************************************************************************

  Function : DDV_Print_Sequence()
  
  Purpose : display a sequence
  
  Parameters : 	szSeq; the sequence to display
  				bUseColour; print out a coloured sequence
  				IsAA; True is szSeq is a prot
				fp;file to put in the sequence

  Note: this function can be used to display both Ticks and Numbers for
        the top scale. See the 'what' parameter
				
  Return value : -

*******************************************************************************/
static CharPtr  DDV_Print_Sequence(CharPtr szSeq,Boolean bUseColour,Boolean IsAA)
{
CharPtr szBuf,szTmp,szSequence=NULL;
Int4    pos,taille,stop;			
Char    residue;					
Uint2   nCompt=0,nCompt2=0;
Char    curColor[10],newColor[10]={""};
	
	/*alloc a buffer*/
	taille=StringLen(szSeq)+1;
	szBuf=(CharPtr)MemNew(taille*sizeof(Char));
	if (!szBuf) {
		/*fprintf(fp,"%s",szSeq);*/
		szSequence=StringSave(szSeq);
		return(szSequence);
	}
	MemSet(szBuf,0,taille*sizeof(Char));
	
	/*prepare the colour display*/
	pos=1;	
	stop=taille+1;

	/*draw!*/
	if (!bUseColour){/*black-colour display*/
		szTmp=StringChr(szSeq,'\0');
		if (szTmp && *(szTmp-1)==' ') *(szTmp-1)='\0';
		/*fprintf(fp,"%s",szSeq);*/
		szSequence=StringSave(szSeq);
	}
	else{ /*use colour-coded alphabet fo NA & AA*/
		while(pos<stop){
			residue=szSeq[nCompt];
			/*retrieve colour*/
			GetClrFromLetter(residue,IsAA,curColor);
			/*new colour?*/
			if (StringCmp(curColor,newColor)){
				if (szBuf[0]!='\0'){/*something to draw ?*/
					szBuf[nCompt2]='\0';/*CLOSE THE STRING*/
					/*fprintf(fp,"<FONT COLOR=#%s>%s</FONT>",newColor,szBuf);*/
					szSequence=ConcatStr(szSequence,"<FONT COLOR=#");
					szSequence=ConcatStr(szSequence,newColor);
					szSequence=ConcatStr(szSequence,">");
					szSequence=ConcatStr(szSequence,szBuf);
					szSequence=ConcatStr(szSequence,"</FONT>");
				}
				StringCpy(newColor,curColor);
				nCompt2=0;
				szBuf[nCompt2]='\0';
			}
			szBuf[nCompt2++]=residue;
			
			/*each LETTER_BLOCK_WIDTH, add a blank column*/
			pos++;nCompt++;		
		}
		if (szBuf[0]!='\0'){/*something to draw ?*/
			szBuf[nCompt2]='\0';/*CLOSE THE STRING*/
			/*fprintf(fp,"<FONT COLOR=#%s>%s</FONT>",newColor,szBuf);*/
			szSequence=ConcatStr(szSequence,"<FONT COLOR=#");
			szSequence=ConcatStr(szSequence,newColor);
			szSequence=ConcatStr(szSequence,">");
			szSequence=ConcatStr(szSequence,szBuf);
			szSequence=ConcatStr(szSequence,"</FONT>");
		}
	}
	MemFree(szBuf);
	return(szSequence);
}

/*******************************************************************************

  Function : DDV_DisplayTopScale()
  
  Purpose : display a Top scale above the first sequence of the alignment
  
  Parameters : 	from; from 
  				to;  to ... I have no other explanation !
  				offset;used to locate correctly the ruler
				disp_format; display format options
				what;tick or number

  Note: this function can be used to display both Ticks and Numbers for
        the top scale. See the 'what' parameter
				
  Return value : -

*******************************************************************************/
#define DDV_AddBlank(_diff, _fp) \
do { Int4 _i;  for (_i=0; _i<(Int4)_diff; _i++) fprintf(_fp," "); } while(0)

#define DDV_AddBlank2(_diff, _str) \
do { Int4 _i;  for (_i=0; _i<(Int4)_diff; _i++) *_str=ConcatStr(*_str," "); } while(0)

static CharPtr DDV_DisplayTopScale(Int4 from,Int4 to,Int4 offset,Uint4 disp_format,
	Uint1 what)
{
CharPtr	szText,szScale=NULL;
Int4 	pos,i,stop,size,j,k;
Char 	szPos[15]={""};	
Char 	szBlank[TABLE_LEFT_MARGIN+5]={""};	

	size=to-from+1;
	pos=from+1;
	stop=to+2;
	i=0;
	if (disp_format&DISPE_SHOWBLOCK) size+=(size/LETTER_BLOCK_WIDTH);
	szText=(CharPtr)MemNew((size+1)*sizeof(Char));
	if (szText){
		MemSet(szText,' ',size*sizeof(Char));
		while(pos<stop){
			if (!(pos % LETTER_BLOCK_WIDTH)){
				switch(what){
					case SCALE_TICK:
						szText[i]='|';
						break;
					case SCALE_NUM:
						sprintf(szPos,"%i",pos);
						j=StringLen(szPos);
						for(k=0;k<j;k++){
							szText[i-j+k+1]=szPos[k];
						}
						break;
				}
				if (disp_format&DISPE_SHOWBLOCK) i++;
			}
			if (!(pos % LETTER_HALF_BLOCK_WIDTH) && 
				(pos % LETTER_BLOCK_WIDTH) && (what==SCALE_TICK)){
				szText[i]='.';
			}
			pos++;i++;
		}
		if (what==SCALE_NUM) /*fprintf(fp,"\n");*/
			szScale=ConcatStr(szScale,"\n");
		if(disp_format&DISP_FULL_HTML){
			/*fprintf(fp,"<font color=#700777>");*/
			szScale=ConcatStr(szScale,"<font color=#700777>");
			DDV_AddBlank2(offset,&szScale);
			szScale=ConcatStr(szScale,szText);
			/*fprintf(fp,"%s",szText);*/
			szScale=ConcatStr(szScale,"</font>\n");
			/*fprintf(fp,"</font>\n");*/
		}else if(disp_format&DISP_FULL_TXT){
			DDV_AddBlank2(offset,&szScale);
			/*fprintf(fp,"%s\n",szText);*/
			szScale=ConcatStr(szScale,szText);
			szScale=ConcatStr(szScale,"\n");
		}
		MemFree(szText);
	}
	return(szScale);
}

/*******************************************************************************

  Function : DDV_GetSequenceFromParaG()
  
  Purpose : retrieve a sequence from a ParaG 
  
  Parameters : pgp ; a ParaG
				szSequence; filled here
				bspLength; size of the sequence
				IsAA; true if the sequence is a protein
  				bsp_strand; strand of the ParaG
				bsp_start; from
				bsp_stop; to of the Parag (bsp coord)

  Note : bsp_strand,bsp_start,bsp_stop ; pointers can be NULL if you don't want
  		to retrieve that information.

  Return value : TRUE if success; see also bsp_strand,bsp_start,bsp_stop

*******************************************************************************/
extern Boolean DDV_GetSequenceFromParaG(ParaGPtr pgp,
	CharPtr PNTR szSequence,Int4 bspLength,Boolean IsAA,Uint1 PNTR bsp_strand,
	Int4 PNTR bsp_start,Int4 PNTR bsp_stop)
{
Int4          size,from,to,ret,seek;
CharPtr       SeqBuf;
ValNodePtr    vnp2;
MsaTxtDispPtr mtdp;
SeqLocPtr     slp=NULL;
SeqPortPtr    spp=NULL;
Boolean       bOk=FALSE;
SeqIdPtr      sip;
Uint1		  strand=Seq_strand_unknown;

	if (!szSequence || !pgp) 
		return(FALSE);

	/* Test validity of Blockfied structure - Debug only*/
/*
{{	
Int4 last_from=0,old_from=INT4_MAX;
	for(vnp2=pgp->ptxtList;vnp2!=NULL;vnp2=vnp2->next){
		mtdp=(MsaTxtDispPtr)(vnp2->data.ptrvalue);
		if (mtdp){
			if (mtdp->IsGap==0) {
				strand=mtdp->strand;
				from = mtdp->from;
				to = mtdp->to;
				if(from>to)
					fprintf(stderr,"Invalid block order of from,t\n");
				if (strand!=Seq_strand_minus) {
					if(from<last_from)
						fprintf(stderr,"+ non monotonic\n");
					last_from =from; 
				} else {
					if(from>old_from)
						fprintf(stderr,"- non-monotonic\n");
					old_from = from;
				}
			}
		}
	}
}}
*/
	/*analyse the sequence description; get the first block
	which describes a Bioseq, not a gap; then retrieves the from*/
	from=INT4_MAX;
	to=INT4_MIN;
	/*strand=(Uint1)-1;*/
	for(vnp2=pgp->ptxtList;vnp2!=NULL;vnp2=vnp2->next){
		mtdp=(MsaTxtDispPtr)(vnp2->data.ptrvalue);
		if (mtdp){
			strand=mtdp->strand;
			if (mtdp->IsGap==0) {
				from=MIN(from,mtdp->from);
				to=MAX(to,mtdp->to);
				if (!bOk) bOk=TRUE;
				/*break;*/
			}
		}
	}

	
	/*if bOk==FALSE, we have only a big gap to display*/
	if (bOk){
		/*open a seqport on the ParaG*/
		sip=SeqIdFindBest(pgp->sip,SEQID_GI);
		if (sip==NULL) sip=pgp->sip;
		slp=SeqLocIntNew(from,to,strand, sip);
		if (!slp) return(FALSE);
		spp = SeqPortNewByLoc (slp, (Uint1)(IsAA==TRUE ? Seq_code_iupacaa 
				: Seq_code_iupacna));
		if (!spp) return(FALSE);
	}

	/*build the text to display for the ParaG: AATAT---ATAT--AAA*/
	for(vnp2=pgp->ptxtList;vnp2!=NULL;vnp2=vnp2->next){
		mtdp=(MsaTxtDispPtr)vnp2->data.ptrvalue;
		if (mtdp){
			Char c;
			size=mtdp->to-mtdp->from+1;
			SeqBuf=(CharPtr)MemNew((size+1)*sizeof(Char));
			if (!SeqBuf) continue;
			if(mtdp->IsGap){
				switch(mtdp->TextStyle){
					case MSA_TXT_STYLE_1:
						c=' ';
						break;
					case MSA_TXT_STYLE_GAP:
						c='-';
						break;
				}
				MemSet(SeqBuf,c,size*sizeof(Char));
				SeqBuf[size]='\0';
			}
			else{
				if (strand == Seq_strand_minus)
					seek = to-mtdp->to;
				else
					seek = mtdp->from-from;
				SeqPortSeek(spp, seek, SEEK_SET);
				ret=SeqPortRead (spp,(Uint1Ptr)SeqBuf, (Int2)size);
				if (ret!=size){
					MemSet(SeqBuf,'?',size*sizeof(Char));
					SeqBuf[size]='\0';
				}
			}
			StringCat(*szSequence,SeqBuf);
			MemFree(SeqBuf);
			SeqBuf=NULL;
			
		}
	}

	if (spp) SeqPortFree (spp);
	if (slp) SeqLocFree (slp);

	if (bsp_strand) *bsp_strand=strand;
	if (bsp_start) *bsp_start=from;
	if (bsp_stop) *bsp_stop=to;

	return(TRUE);
}

/*******************************************************************************

  Function : DDV_DisplayParaG()
  
  Purpose : display the content of a ParaG
  
  Parameters : 	vnp; the node containing the ParaG
  				i; sequence number (from 1 to nSequence)
  				ShowScale;if true, show the top scale
				disp_format; display format options
				szSeqAcc; Accession number of the sequence
				szSeqName; sequence name
				IsAA; TRUE means protein
				bDispPhilouName;used for PHYLIP format
				szFormattedLine; see 'Return value'
				
  Return value : the next node in the ParaG ValNode List of a sequence;
                fill szFormattedLine with the info ready to display

*******************************************************************************/
static ValNodePtr DDV_DisplayParaG(ValNodePtr vnp,Int4 i,Boolean ShowScale,
		Uint4 disp_format,CharPtr szSeqAcc,CharPtr szSeqName,Boolean IsAA,
		Boolean bDispPhilouName,Int4 max_length,Int4 max_scale,
		CharPtr PNTR szFormattedLine)		
{
static CharPtr  szSeqMaster=NULL;
ParaGPtr 		pgp;	/*ParaG data*/
CharPtr			szSequence,szDisp,szFLine=NULL,szTmp;
BioseqPtr       bsp;
SeqIdPtr        sip;
Int4            bspLength,size,stop,nCompt2,pos,gi,diff,disp,bsp_start,bsp_stop;
Char 			szBuf4[WWW_SCRIPT_SIZE]={""};	/*Entrez query*/
Uint1			bsp_strand;
	
	if (vnp==NULL) return(NULL);

	pgp=(ParaGPtr)vnp->data.ptrvalue;
	if (!pgp) return(vnp->next);	

	szSequence=(CharPtr)MemNew((pgp->StopLetter-pgp->StartLetter+3)*sizeof(Char));

	if (szSequence) MemSet(szSequence,0,(pgp->StopLetter-pgp->StartLetter+2)
			*sizeof(Char));
	else return(NULL);

	if(i==1){/*to retrieve the sequence of the Master*/
		if (szSeqMaster) MemFree(szSeqMaster);
	}

	/*get some bsp info*/	
	sip=SeqIdFindBest(pgp->sip,0);
	if (sip==NULL) sip=pgp->sip;
	bsp=BioseqLockById(sip);
	if (!bsp) goto error;
	bspLength=BioseqGetLen(bsp);
	BioseqUnlockById(sip);
	if (!DDV_GetSequenceFromParaG(pgp,&szSequence,bspLength,IsAA,&bsp_strand,
			&bsp_start,&bsp_stop)) goto error;

	/*switch to lower case if needed*/
	if (disp_format&TEXT_LOWERCASE){
		szSequence=ConvertToLower(szSequence);
	}
	
	/*7 : StringLen("unknown"); see below, when I print szSeqAcc*/
	if (max_length<7) max_length=7;
	disp=0;

	if(disp_format&DISP_FULL_HTML || disp_format&DISP_FULL_TXT){/*output format: html*/
		/*compute the left part size of the line to display*/
		if (disp_format&DISP_ORDER_NUM) {disp+=7;max_length+=7;}
		if (*szSeqAcc) disp+=StringLen(szSeqAcc);
		else disp+=7;
		diff=max_length-disp;
		diff++;
		disp+=diff;
		if (disp_format&DISP_STRAND) disp+=2;
		if (disp_format&DISP_BSP_COORD)disp+=max_scale+2;
		/*display top scale if needed*/
		if ((disp_format&RULER_TOP) && ShowScale) {
			szTmp=DDV_DisplayTopScale(pgp->StartLetter,
				pgp->StopLetter,disp,disp_format,SCALE_NUM);
			if (szTmp){
				szFLine=ConcatStr(szFLine,szTmp);
				MemFree(szTmp);
			}
		}
		if ((disp_format&RULER_TICK) && ShowScale) {
			szTmp=DDV_DisplayTopScale(pgp->StartLetter,
				pgp->StopLetter,disp,disp_format,SCALE_TICK);
			if (szTmp){
				szFLine=ConcatStr(szFLine,szTmp);
				MemFree(szTmp);
			}
		}
		/*display order num if needed*/
		if (disp_format&DISP_ORDER_NUM){ /*fprintf(fp,"[%4i] ",i);*/
			sprintf(szBuf4,"[%4i] ",i);
			szFLine=ConcatStr(szFLine,szBuf4);
		}
		/*display sequence name*/
		if (*szSeqAcc) {
			if (disp_format&DISP_FULL_HTML){
				gi=GetGIForSeqId(pgp->sip);
				if(gi>0){
					sprintf(szBuf4,szEntrezScript,gi,(IsAA ? "p" : "n"));
				}
				else{
					sprintf(szBuf4,"javascript:void(0)");
				}
				if (szSeqName){
					szFLine=ConcatStr(szFLine,"<a href=\"");
					szFLine=ConcatStr(szFLine,szBuf4);
					szFLine=ConcatStr(szFLine,"\" onMouseOut=\"window.status=''\" \nonMouseOver=\"window.status='");
					szFLine=ConcatStr(szFLine,szSeqName);
					szFLine=ConcatStr(szFLine,"';return true\">");
					/*fprintf(fp,szSeqNameStatus,szBuf4,szSeqName);*/
				}
				else{
					szFLine=ConcatStr(szFLine,"<a href=\"");
					szFLine=ConcatStr(szFLine,szBuf4);
					szFLine=ConcatStr(szFLine,"\" onMouseOut=\"window.status=''\" \nonMouseOver=\"window.status='");
					szFLine=ConcatStr(szFLine,"unknown name");
					szFLine=ConcatStr(szFLine,"';return true\">");
					/*fprintf(fp,szSeqNameStatus,szBuf4,"unknown name");*/
				}
			}
			/*fprintf(fp,"%s",szSeqAcc);*/
			szFLine=ConcatStr(szFLine,szSeqAcc);
			if (disp_format&DISP_FULL_HTML) 
				/*fprintf(fp,"</a>");*/
				szFLine=ConcatStr(szFLine,"</a>");
		}
		else {
			/*fprintf(fp,"unknown");*/
			szFLine=ConcatStr(szFLine,"unknown");
		}
		DDV_AddBlank2(diff,&szFLine);
		/*display strand if needed*/
		if (disp_format&DISP_STRAND){
			/*fprintf(fp,"%-2s", (bsp_strand==Seq_strand_minus ? "<" : ">"));*/
			sprintf(szBuf4,"%-2s", (bsp_strand==Seq_strand_minus ? "<" : ">"));
			szFLine=ConcatStr(szFLine,szBuf4);
		}
		/*display bsp coord if needed*/
		if (disp_format&DISP_BSP_COORD){
			Char szBuf2[15]={""};
			Int4 iVal;
			
			if (bsp_strand==Seq_strand_minus) iVal=bsp_stop;
			else iVal=bsp_start;

			if (iVal!=INT4_MAX) {
				sprintf(szBuf2,"%i",++iVal);/*switch to base one*/
				/*fprintf(fp,"%i", iVal);*/
				szFLine=ConcatStr(szFLine,szBuf2);
				DDV_AddBlank2(MAX(2,max_scale-StringLen(szBuf2)+2),&szFLine);
			} 
			else DDV_AddBlank2(max_scale+2,&szFLine);
		}
	}else if (disp_format&DISP_PHYLIP_TXT){
		if (bDispPhilouName){/*Phylip format: seq. names for the first block only*/
			if (szSeqAcc) {/*name only*/
				/*fprintf(fp,"%-10s ",szSeqAcc);*/
				sprintf(szBuf4,"%-10s ",szSeqAcc);
				szFLine=ConcatStr(szFLine,szBuf4);
			}
			else {
				/*fprintf(fp,"%-10s ","unknown");*/
				szFLine=ConcatStr(szFLine,"unknown");
			}
		}
		else {
			/*fprintf(fp,"%-10s "," ");*/
			sprintf(szBuf4,"%-10s "," ");
			szFLine=ConcatStr(szFLine,szBuf4);
		}
	}
	if (i==1){
		szSeqMaster=StringSave(szSequence);
	}
	/*display sequence*/
	size=StringLen(szSequence);
	if (disp_format&DISPE_SHOWBLOCK){/*AFSDDGTFDG GDFDTEGDFG DFDGETDFGD*/
		stop=size+1;
		size+=(size/LETTER_BLOCK_WIDTH);
		szDisp=(CharPtr)MemNew((size+1)*sizeof(Char));
		if (szDisp){
			MemSet(szDisp,' ',size*sizeof(Char));
			pos=1;
			nCompt2=0;
			while(pos<stop){
				if(disp_format&VIEW_FULLSEQ){
					szDisp[nCompt2++]=szSequence[pos-1];
				}else if(disp_format&VIEW_VARIA){
					if(szSeqMaster && i!=1){
						if(szSequence[pos-1]=='-'){
							szSequence[pos-1]='-';
						}
						if(szSequence[pos-1]!=szSeqMaster[pos-1]){
							szDisp[nCompt2++]=szSequence[pos-1];
						}
						else{
							if(szSequence[pos-1]!='-') szDisp[nCompt2++]='.';
							else szDisp[nCompt2++]=szSequence[pos-1];
						}
					}
					else{
						szDisp[nCompt2++]=szSequence[pos-1];
					}
				}
				if (!(pos % LETTER_BLOCK_WIDTH)) {
					szDisp[nCompt2++]=' ';
				}
				pos++;
			}
			szTmp=DDV_Print_Sequence(szDisp,(Boolean)(disp_format&DISPE_COLOR),IsAA);
			if (szTmp){
				szFLine=ConcatStr(szFLine,szTmp);
				MemFree(szTmp);
			}
			MemFree(szDisp);
		}
		else{
			szTmp=DDV_Print_Sequence(szSequence,(Boolean)(disp_format&DISPE_COLOR),IsAA);
			if (szTmp){
				szFLine=ConcatStr(szFLine,szTmp);
				MemFree(szTmp);
			}
		}
	}
	else{/*AFSDDGTFDGGDFDTEGDFGDFDGETDFGD*/
		if(disp_format&VIEW_FULLSEQ){
			if(disp_format&DISP_FULL_HTML){/*output format*/
				szTmp=DDV_Print_Sequence(szSequence,(Boolean)(disp_format&DISPE_COLOR),IsAA);
				if (szTmp){
					szFLine=ConcatStr(szFLine,szTmp);
					MemFree(szTmp);
				}
			}else {
				/*fprintf(fp,"%s",szSequence);*/
				szFLine=ConcatStr(szFLine,szSequence);
			}
		}else if(disp_format&VIEW_VARIA){
			szDisp=(CharPtr)MemNew((size+1)*sizeof(Char));
			if (szDisp){
				MemSet(szDisp,' ',size*sizeof(Char));
				pos=0;
				stop=size+1;
				while(pos<stop){
					if(szSeqMaster && i!=1){
						if(szSequence[pos]!=szSeqMaster[pos]){
							szDisp[pos]=szSequence[pos];
						}
						else{
							if(pos==stop-1)
								szDisp[pos]=szSequence[pos];
							else
								szDisp[pos]='.';
						}
					}
					else{
						szDisp[pos]=szSequence[pos];
					}
					pos++;
				}
				if(disp_format&DISP_FULL_HTML){/*output format*/
					szTmp=DDV_Print_Sequence(szDisp,(Boolean)(disp_format&DISPE_COLOR),IsAA);
					if (szTmp){
						szFLine=ConcatStr(szFLine,szTmp);
						MemFree(szTmp);
					}
				}else{
					/*fprintf(fp,"%s",szDisp);*/
					szFLine=ConcatStr(szFLine,szDisp);
				}
				MemFree(szDisp);
			}
			else {
				if(disp_format&DISP_FULL_HTML){/*output format*/
					szTmp=DDV_Print_Sequence(szSequence,(Boolean)(disp_format&DISPE_COLOR),IsAA);
					if (szTmp){
						szFLine=ConcatStr(szFLine,szTmp);
						MemFree(szTmp);
					}
				}else{
					/*fprintf(fp,"%s",szSequence);*/
					szFLine=ConcatStr(szFLine,szSequence);
				}
			}
		}
	}
	/*display bsp coord if needed*/
	if (disp_format&DISP_BSP_COORD){
		Int4 iVal;

		if (bsp_strand==Seq_strand_minus) iVal=bsp_start;
		else iVal=bsp_stop;
		
		if (iVal!=INT4_MIN) {
			/*fprintf(fp,"  %i", iVal+1);*/
			sprintf(szBuf4,"  %i", iVal+1);
			szFLine=ConcatStr(szFLine,szBuf4);
		}
		else DDV_AddBlank2(max_scale+2,&szFLine);
	}
	
error:
	MemFree(szSequence);
	szFLine=ConcatStr(szFLine,"\n");
	/*fprintf(fp,"\n");*/
	*szFormattedLine=szFLine;
	
	return(vnp->next);
}

/*******************************************************************************

  Function : DDV_Nav_Arrows()
  
  Purpose : display the "navigator" links (to jump to next/previous blocks)
  
  Parameters : 	
				gi; gi number used by PopSet viewer to build the HTTP links
				
  				TotalAliLength; SeqAlign size
				numBlockAffich; block number displayed (when display by block is
						required; otherwise numBlockAffich=0).
				disp_format; full HTML, Phylip, ...
				LineSize; size of one text line (number of letter for the sequence)
				
  Return value : -

*******************************************************************************/
static CharPtr DDV_Nav_Arrows(Int4 gi,Int4 TotalAliLength,
		Int4 numBlockAffich,Uint4 disp_format,Int2 LineSize)
{
CharPtr	szWide,szFLine=NULL;
Int4 	k,wide,lGoLeft,lGoRight,fromL,toL,fromR,toR;
Char	szGoLeft[30]={""};
Char	szGoRight[30]={""};
Char 	szWWW[WWW_SCRIPT_SIZE]={""};

	k=TotalAliLength/LineSize;
	if (TotalAliLength%LineSize) k++;
	if (numBlockAffich>1){/*left arrow*/
		fromL=(numBlockAffich-2)*LineSize+1;
		toL=_min_((numBlockAffich-1)*LineSize,TotalAliLength);
		sprintf(szGoLeft,"[%d..%d]",fromL,toL);
		lGoLeft=StringLen(szGoLeft)+2;
	}
	else{
		StringCpy(szGoLeft," ");
		fromL=0;
		toL=0;
		lGoLeft=1;
	}
	if (numBlockAffich<k){/*right arrow*/
		fromR=(numBlockAffich)*LineSize+1;
		toR=_min_((numBlockAffich+1)*LineSize,TotalAliLength);
		sprintf(szGoRight,"[%d..%d]",fromR,toR);
		lGoRight=StringLen(szGoRight)+2;
	}
	else{
		StringCpy(szGoRight," ");
		fromR=0;
		toR=0;
		lGoRight=1;
	}
	wide=LineSize+28;
	if (disp_format&DISPE_SHOWBLOCK) wide+=(LineSize/LETTER_BLOCK_WIDTH);
	szWide=(CharPtr)MemNew(wide*sizeof(Char));
	MemSet(szWide,'.',(wide-lGoLeft-lGoRight)*sizeof(Char));
	szWide[wide-lGoLeft-lGoRight]='\0';
	/*fprintf(fp,"<pre>");*/
	szFLine=ConcatStr(szFLine,"<pre>");
	if(fromL!=0 && toL!=0){
		sprintf(szWWW,WWWDDV_script,gi,fromL,toL,disp_format);
		/*fprintf(fp,"&lt;&lt<a href=\"%s\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to view this region';return true\" >%s</a>",szWWW,szGoLeft);*/
		szFLine=ConcatStr(szFLine,"&lt;&lt<a href=\"");
		szFLine=ConcatStr(szFLine,szWWW);
		szFLine=ConcatStr(szFLine,"\" onMouseOut=\"window.status='';return true\" ");
		szFLine=ConcatStr(szFLine,"onMouseOver=\"window.status='Click to view this region';return true\" >");
		szFLine=ConcatStr(szFLine,szGoLeft);
		szFLine=ConcatStr(szFLine,"</a>");
	}
	else{
		/*fprintf(fp,"%s",szGoLeft);*/
		szFLine=ConcatStr(szFLine,szGoLeft);
	}
	/*fprintf(fp,"<FONT COLOR=#FFFFFF>%s</FONT>",szWide);*/
	szFLine=ConcatStr(szFLine,"<FONT COLOR=#FFFFFF>");
	szFLine=ConcatStr(szFLine,szWide);
	szFLine=ConcatStr(szFLine,"</FONT>");
	
	if(fromR!=0 && toR!=0){
		sprintf(szWWW,WWWDDV_script,gi,fromR,toR,disp_format);
		/*fprintf(fp,"<a href=\"%s\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to view this region';return true\" >%s</a>&gt;&gt",szWWW,szGoRight);*/
		szFLine=ConcatStr(szFLine,"<a href=\"");
		szFLine=ConcatStr(szFLine,szWWW);
		szFLine=ConcatStr(szFLine,"\" onMouseOut=\"window.status='';return true\" ");
		szFLine=ConcatStr(szFLine,"onMouseOver=\"window.status='Click to view this region';return true\" >");
		szFLine=ConcatStr(szFLine,szGoRight);
		szFLine=ConcatStr(szFLine,"</a>&gt;&gt");
	}
	else{
		/*fprintf(fp,"%s",szGoRight);*/
		szFLine=ConcatStr(szFLine,szGoRight);
	}
	/*fprintf(fp,"</PRE>");*/
	szFLine=ConcatStr(szFLine,"</pre>");

	MemFree(szWide);
	return(szFLine);
}

/*******************************************************************************

  Function : DDV_AffichageParaG()
  
  Purpose : function to create an HTML/TXT/PHYLIP display of an Indexed SeqAlign
  
  Parameters : 	mpplp; contain the data to create a display (sap, ...)
				gi; gi number used by PopSet viewer to build the HTTP links
				from - to; range of the SeqAlign to display
  				TotalAliLength; SeqAlign size
				numBlockAffich; block number displayed (when display by block is
						required; otherwise numBlockAffich=0).
				disp_format; full HTML, Phylip, ...
				LineSize; size of each line of sequence
				fp; file to put the ouptut
				bspp; filled here withthe formatted text
				 
  Return value : - (see bspp)

  Note : if fp==NULL, formatted text will be stored in bspp; Otherwise the
         text is directly sent to fp.
		 
*******************************************************************************/
extern void DDV_AffichageParaG(MsaParaGPopListPtr mpplp,Int4 gi,Int4 from,Int4 to,
		Int4 TotalAliLength,Int4 numBlockAffich,Uint4 disp_format,Int2 LineSize,
		FILE *fp,ByteStorePtr PNTR bspp)
{
Int4 			i,j=0,max_length=INT4_MIN,max_scale=INT4_MIN;
ValNodePtr PNTR TableDesc=NULL;
CharPtr 		szColor[]={"#E6EEEE","ccccff","FFFFFF"};
CharPtr 		*szSeqAcc;
CharPtr 		*szSeqName;
CharPtr 		szFline,szTmp;
Char 			szBuf[25],szBuf2[25];
ParaGPtr 		pgp;
BioseqPtr 		bsp;
BoolPtr 		IsAAp;
Boolean			bDispPhilouName=TRUE;
SeqIdPtr        sip;

	/*if bspp is not allocated*/
	if (fp==NULL && (*bspp)==NULL){
		(*bspp)=BSNew((Int4)1);
		if ((*bspp)==NULL) return;
	}
	
	/*prepare a list of ParaG heads, names, Acc Number and Mol_type*/
	TableDesc=(ValNodePtr PNTR)MemNew(mpplp->nBsp*sizeof(ValNodePtr));
	szSeqAcc=(CharPtr *)MemNew(mpplp->nBsp*sizeof(CharPtr));
	szSeqName=(CharPtr *)MemNew(mpplp->nBsp*sizeof(CharPtr));
	IsAAp=(BoolPtr)MemNew(mpplp->nBsp*sizeof(Boolean));

	if (!TableDesc || !szSeqAcc || !szSeqName || !IsAAp) goto error;
	MemSet(szSeqAcc,0,mpplp->nBsp*sizeof(CharPtr));	
	MemSet(szSeqName,0,mpplp->nBsp*sizeof(CharPtr));	
	MemSet(IsAAp,0,mpplp->nBsp*sizeof(Boolean));	
	if (numBlockAffich>0) j=2;
	else j=0;

	for(i=0;i<mpplp->nBsp;i++){
		TableDesc[i]=mpplp->TableHead[i];/*copy the ValNode heads in TableDesc*/
		/*get the access code*/
		pgp=(ParaGPtr)TableDesc[i]->data.ptrvalue;
		if (pgp){/*create the name table*/
			bsp=BioseqLockById(pgp->sip);
			if (bsp){
				Int2 j;
				if (disp_format&SEQNAME_BLAST_STD){ 
					if (i==0) 
						StringCpy(szBuf,"Query:");
					else
						StringCpy(szBuf,"Sbjct:");
				}
				else{
					sip = SeqIdFindBest(bsp->id, SEQID_GENBANK);
					if (!sip)
				   		sip = SeqIdFindBest(bsp->id, 0);
					SeqIdWrite(sip, szBuf,PRINTID_TEXTID_ACCESSION, 14);   
				}
				/*access code*/
				szSeqAcc[i]=StringSave(szBuf);
				max_length=MAX(max_length,(Int4)StringLen(szSeqAcc[i]));
				/*seq title*/
				szSeqName[i]=(CharPtr)MemNew(100*sizeof(Char));
				if (szSeqName[i]){
					FastaDefLine(bsp, szSeqName[i], 99, NULL, NULL, 0);
					j=0;
					while(szSeqName[i][j]){
						if (szSeqName[i][j]=='\'') szSeqName[i][j]=' ';
						j++;
					}
				}
				
				IsAAp[i]=ISA_aa(bsp->mol);
				/*size for the BSP coord display*/
				sprintf(szBuf2,"%i",BioseqGetLen(bsp));
				max_scale=MAX(max_scale,(Int4)StringLen(szBuf2));
				BioseqUnlock(bsp);
			}
		}
	}

	szFline=NULL;
	if(disp_format&DISP_FULL_HTML){/*output format*/
		if (numBlockAffich>0){/*nav arrow*/
			szTmp=DDV_Nav_Arrows(gi,TotalAliLength,numBlockAffich,disp_format,LineSize);
			szFline=ConcatStr(szFline,szTmp);
			MemFree(szTmp);
		}
		/*fprintf(fp,"<table>");*/
		/*fprintf(fp,"<tr><td bgcolor=%s><pre>",szColor[j]);*/
		szFline=ConcatStr(szFline,"<table><tr><td bgcolor=");
		szFline=ConcatStr(szFline,szColor[j]);
		szFline=ConcatStr(szFline,"><pre>");
	}else if(disp_format&DISP_PHYLIP_TXT){
		/*fprintf(fp,"%7d  %d\n",mpplp->nBsp,mpplp->LengthAli);*/
		sprintf(szBuf,"%7d  %d\n",mpplp->nBsp,mpplp->LengthAli);
		szFline=ConcatStr(NULL,szBuf);
	}

	if(szFline){
		if (fp)
			fprintf(fp,"%s",szFline);
		else
			BSWrite(*bspp, (VoidPtr)szFline, StringLen(szFline));
		MemFree(szFline);
	}

	i=1;
	while(TRUE){
		szFline=NULL;
		TableDesc[i-1]=DDV_DisplayParaG(TableDesc[i-1],i,(Boolean)(i==1 ? TRUE : FALSE),
			disp_format,szSeqAcc[i-1],szSeqName[i-1],IsAAp[i-1],bDispPhilouName,
			max_length,max_scale,&szFline);
		if (szFline){
			if (fp)
				fprintf(fp,"%s",szFline);
			else
				BSWrite(*bspp, (VoidPtr)szFline, StringLen(szFline));
			/*fprintf(fp,"%s",szFline);*/
			MemFree(szFline);
		}
		i++;
		if (i>mpplp->nBsp){
			i=1;
			if(disp_format&DISP_FULL_HTML){
				/*fprintf(fp,"</pre></td></tr>");*/
				szFline=ConcatStr(NULL,"</pre></td></tr>");
				if(szFline){
					if (fp)
						fprintf(fp,"%s",szFline);
					else
						BSWrite(*bspp, (VoidPtr)szFline, StringLen(szFline));
					/*fprintf(fp,"%s",szFline);*/
					MemFree(szFline);
				}
			}else if((disp_format&DISP_PHYLIP_TXT) || (disp_format&DISP_FULL_TXT)){
				/*fprintf(fp,"\n");*/
				szFline=ConcatStr(NULL,"\n");
				if(szFline){
					if (fp)
						fprintf(fp,"%s",szFline);
					else
						BSWrite(*bspp, (VoidPtr)szFline, StringLen(szFline));
					/*fprintf(fp,"%s",szFline);*/
					MemFree(szFline);
				}
			}
			if (numBlockAffich==0){
				if (j==1) j=0; else j=1;
			}
			if(disp_format&DISP_FULL_HTML){
				/*fprintf(fp,"<tr><td bgcolor=%s><pre>",szColor[j]);*/
				szFline=ConcatStr(NULL,"<tr><td bgcolor=");
				szFline=ConcatStr(szFline,szColor[j]);
				szFline=ConcatStr(szFline,"><pre>");
				if(szFline){
					if (fp)
						fprintf(fp,"%s",szFline);
					else
						BSWrite(*bspp, (VoidPtr)szFline, StringLen(szFline));
					/*fprintf(fp,"%s",szFline);*/
					MemFree(szFline);
				}
			}else if(disp_format&DISP_PHYLIP_TXT){
				bDispPhilouName=FALSE;
			}
			if (TableDesc[mpplp->nBsp-1]==NULL) break;
		} 
	}

	if(disp_format&DISP_FULL_HTML){
		/*fprintf(fp,"</pre></td></tr>");
		fprintf(fp,"</table>");*/
		szFline=ConcatStr(szFline,"</pre></td></tr></table>");
		if (numBlockAffich>0)/*nav arrow*/{
			szTmp=DDV_Nav_Arrows(gi,TotalAliLength,numBlockAffich,disp_format,LineSize);
			szFline=ConcatStr(szFline,szTmp);
			MemFree(szTmp);
		}
		if(szFline){
			if (fp)
				fprintf(fp,"%s",szFline);
			else
				BSWrite(*bspp, (VoidPtr)szFline, StringLen(szFline));
			/*fprintf(fp,"%s",szFline);*/
			MemFree(szFline);
		}
	}
	/*end of the display*/
error:
	if (TableDesc) MemFree(TableDesc);
	if (szSeqAcc){
		for(i=0;i<mpplp->nBsp;i++){
			if (szSeqAcc[i]){
				MemFree(szSeqAcc[i]);
			}
		}
		MemFree(szSeqAcc);
	}
	if (szSeqName){
		for(i=0;i<mpplp->nBsp;i++){
			if (szSeqName[i]){
				MemFree(szSeqName[i]);
			}
		}
		MemFree(szSeqName);
	}
	if (IsAAp) MemFree(IsAAp);
}

/*******************************************************************************

  Function : DDV_GetGappedSequence()
  
  Purpose : get a gapped sequence from the indexed seq align
  
  Parameters : 	segp; head of the nodes list of one bsp
  				szSequence; ALLOCATED buffer of size SeqAlign, where the sequence
							will be stored.
				IsAA; doesn't need an explanation, I think...
				
  Note : do not forget to allocate 'szSequence' before entering this function
  
  Return value : none (see szSequence)

*******************************************************************************/
static void DDV_GetGappedSequence(SegmentPtr segp,CharPtr szSequence,
		Boolean IsAA)
{
SegmentPtr sp;
CharPtr SeqBuf;/*I use a pointer because I cannot use a fixed array*/

	if (!szSequence || !segp) return;
	/*loop each segment of the bioseq*/
	for(sp=segp ; sp!=NULL ; sp=sp->bsp_next){
		if(sp->gap){/*gap segment*/
			SeqBuf=(CharPtr)MemNew((sp->gap+1)*sizeof(Char));
			if (SeqBuf){
				MemSet(SeqBuf,'-',(sp->gap)*sizeof(Char));
			}
		}
		else{/*sequence segment*/
			SeqBuf=UDV_Read_SequenceEx (sp->sip,sp->from,sp->to, 
				IsAA,(Int2)(sp->to-sp->from+1),(Uint1)
				(sp->strand == Seq_strand_minus ? Seq_strand_minus : Seq_strand_plus));
			if (!SeqBuf){
				SeqBuf=(CharPtr)MemNew((sp->to-sp->from+2)*sizeof(Char));
				if (SeqBuf){
					MemSet(SeqBuf,'?',(sp->to-sp->from+1)*sizeof(Char));
				}
			}
		}
		if (SeqBuf){
			StringCat(szSequence,SeqBuf);
			MemFree(SeqBuf);
			SeqBuf=NULL;
		}
	}
}

/*******************************************************************************

  Function : DDV_PrintFastaGappedSequence()
  
  Purpose : print a character string as a FASTA sequence format
  
  Parameters : 	szSequence; what do you think it is ?
  				fp; may be a real file or 'stdout'
  
  Return value : none

*******************************************************************************/
#define line_fas_size 80
static void DDV_PrintFastaGappedSequence(CharPtr szSequence,FILE *fp)
{
Char szBuf[line_fas_size+1]={""};
Int4 i=0,j=0,len=StringLen(szSequence);

	if (!szSequence || !fp) return;

	while(szSequence[i]){
		szBuf[j++]=szSequence[i++];
		if (i==len || j==line_fas_size){
			szBuf[j]='\0';
			fprintf(fp,"%s\n",szBuf);
			j=0;
		}
	}
}

/*******************************************************************************

  Function : DDV_GetFullFASTAforIdxAli()
  
  Purpose : create the FASTA file for a full Indexed SeqAlign
  
  Parameters :  sabp; indexed seqalign
  				FasType; see pgppop.h (display format )
				SAsize;size of the indexed seqalign
				fp; file to put the FASTA
  
  Return value : TRUE if success 

*******************************************************************************/
extern Boolean DDV_GetFullFASTAforIdxAli(SABlockPtr sabp,Uint4 disp_format,
	Int4 SAsize,FILE *fp)
{
SegmentPtr	segph;
BioseqPtr	bsp;	
Char DefLine[255];
CharPtr szBuf;

	if(disp_format&DISP_FASTA_NOGAP){
		/*scan the first block in SABP to retrieve each BSP*/
		for (segph=sabp->segp_head; segph!=NULL;segph=segph->next){
			bsp=BioseqLockById(segph->sip);
			if (bsp){
				BioseqToFasta(bsp,fp,(Boolean)ISA_na(bsp->mol));
				BioseqUnlock(bsp);
			}
		}
	}
	else{
		szBuf=(CharPtr)MemNew((SAsize+1)*sizeof(Char));
		if (!szBuf) return (FALSE);
		/*scan the first block in SABP to retrieve each BSP*/
		for (segph=sabp->segp_head; segph!=NULL;segph=segph->next){
			bsp=BioseqLockById(segph->sip);
			if (bsp){
				/*print the FASTA header for the bioseq*/
				FastaId(bsp,DefLine,254);
				fprintf(fp,">%s ",DefLine);
				CreateDefLine(NULL, bsp, DefLine, 254, 0,NULL, NULL);
				fprintf(fp,"%s\n",DefLine);
				/*get and print the sequence with gaps*/
				MemSet(szBuf,0,(SAsize+1)*sizeof(Char));
				DDV_GetGappedSequence(segph,szBuf,(Boolean)ISA_aa(bsp->mol));
				DDV_PrintFastaGappedSequence(szBuf,fp);
				BioseqUnlock(bsp);
			}
		}
	}
	return(TRUE);	
}

/*******************************************************************************

  Function : DDV_PrintStudyName()
  
  Purpose : print the name/author of the PopSet study
  
*******************************************************************************/
extern void DDV_PrintStudyName(CharPtr szPopSetName,CharPtr szPopSetAuth,
		Int4 pmid,FILE *fp)
{
	if (StringLen(szPopSetName)>0 && StringLen(szPopSetAuth)>0) {
		if (StringStr(szPopSetName,"* no title *") || 
			StringStr(szPopSetAuth,"* no title *")){
			fprintf(fp,"<h2>Unpublished study.</h2>");
		}
		else{
			if (pmid!=0){
				fprintf(fp,"<h2>%s (",szPopSetName);
				fprintf(fp,"<a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=%d&dopt=Abstract\">",pmid);
				fprintf(fp,"%s <i>et al.</i></a>)<BR></h2>",szPopSetAuth);
			}
			else{
				fprintf(fp,"<h2>%s (%s <i>et al.</i>, <i>unpublished</i>)<BR></h2>", 
					szPopSetName,szPopSetAuth);
			}
		}
	}
	else{
		fprintf(fp,"<h2>Unpublished study.</h2>");
	}
}

/*******************************************************************************

  Function : PrintHTML_Header()
  
  Purpose : print the header of the HTML file; used only for debugging
  
  Parameters : none
  
  Return value : none 

*******************************************************************************/
extern void PrintHTML_Header(Uint4 type,FILE *fp)
{
	if ((type&DISP_PHYLIP_TXT) || (type&DISP_FASTA_NOGAP) || 
			(type&DISP_FASTA_GAP) || (type&DISP_FULL_TXT)){/*FASTA/PHYLIP output*/
		fprintf(fp,"Content-type: SeqAlign\n\n");
		return;
	}
	fprintf(fp,"Content-type: text/html\n\n");
	fprintf(fp,"<HTML>\n");
	fprintf(fp,"<HEAD>\n");
	fprintf(fp,"<TITLE>WWW DeuxD-Viewer (beta)</TITLE>\n");
	fprintf(fp,"<link rel=\"stylesheet\" href=\"http://www.ncbi.nlm.nih.gov/corehtml/ncbi.css\">");
	fprintf(fp,"</HEAD>\n");
	fprintf(fp,"<BODY BGCOLOR=FFFFFF>\n");
}


/*******************************************************************************

  Function : PrintHTML_Tail()
  
  Purpose : print the tail of the HTML file; used only for debugging
  
  Parameters : none
  
  Return value : none 

*******************************************************************************/
extern void PrintHTML_Tail(Uint4 type,FILE *fp)
{
	if ((type&DISP_PHYLIP_TXT) || (type&DISP_FASTA_NOGAP) || 
			(type&DISP_FASTA_GAP) || (type&DISP_FULL_TXT)){
		return;
	}
	fprintf(fp,"\n</BODY>\n");
	fprintf(fp,"</HTML>\n");
}

static void printGIorGIs(URLDataPtr udp,FILE *fp){
	if (!(udp->disp_format&ALIGN_REBUILD)){
		fprintf(fp,"<INPUT TYPE=\"hidden\" NAME=\"gi\" VALUE=\"%d\"> \n",udp->gi);
	}
	else{
		fprintf(fp,"<INPUT TYPE=\"hidden\" NAME=\"gis\" VALUE=\"%s\"> \n",udp->szGIs);
	}
}
/******************************************************************************/
static void PrintSaveScript(URLDataPtr udp,FILE *fp)
{
Uint4 disp_format;


	/*form for text format*/
	disp_format=RULER_TOP;
	disp_format|=RULER_TICK;
	disp_format|=DISPE_SHOWBLOCK;
	if (udp->disp_format&VIEW_FULLSEQ) disp_format|=VIEW_FULLSEQ;
	else disp_format|=VIEW_VARIA;
	disp_format|=DISP_FULL_TXT;
	if (udp->disp_format&DISP_STRAND) disp_format|=DISP_STRAND;
	if (udp->disp_format&DISP_BSP_COORD) disp_format|=DISP_BSP_COORD;
	if (udp->disp_format&ALIGN_REBUILD) disp_format|=ALIGN_REBUILD;
	
	fprintf(fp,"<br><span class=\"GUTTER1\">Save Alignment:</span><br>\n");
	fprintf(fp,"<FORM NAME=\"alitxt\" ACTION=\"%s\" METHOD=\"GET\" \n",WWWDDV_save_script1);
	fprintf(fp,"ENCTYPE=\"application/x-www-form-urlencoded\"> \n");
	printGIorGIs(udp,fp);
	fprintf(fp,"<INPUT TYPE=\"hidden\" NAME=\"disp\" VALUE=\"%d\"> \n",disp_format);
	fprintf(fp,"<a href=\"javascript:document.alitxt.submit()\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to save the sequence alignment in a file (text format)';\
return true\" class=\"GUTTER2\">Text format</a>");
	fprintf(fp,"</FORM> \n");

	/*form for Phylip format*/
	disp_format=VIEW_FULLSEQ;
	disp_format|=DISP_PHYLIP_TXT;
	if (udp->disp_format&ALIGN_REBUILD) disp_format|=ALIGN_REBUILD;
	fprintf(fp,"<FORM NAME=\"philou\" ACTION=\"%s\" METHOD=\"GET\" \n",WWWDDV_save_script1);
	fprintf(fp,"ENCTYPE=\"application/x-www-form-urlencoded\"> \n");
	printGIorGIs(udp,fp);
	fprintf(fp,"<INPUT TYPE=\"hidden\" NAME=\"disp\" VALUE=\"%d\"> \n",disp_format);
	fprintf(fp,"<a href=\"javascript:document.philou.submit()\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to save the sequence alignment in a file (PHYLIP format)';\
return true\" class=\"GUTTER2\">PHYLIP format</a>");
	fprintf(fp,"</FORM> \n");

	/*form 1 for FASTA*/
	disp_format=DISP_FASTA_GAP;
	if (udp->disp_format&ALIGN_REBUILD) disp_format|=ALIGN_REBUILD;
	fprintf(fp,"<FORM NAME=\"fasta1\" ACTION=\"%s\" METHOD=\"GET\" \n",WWWDDV_save_script1);
	fprintf(fp,"ENCTYPE=\"application/x-www-form-urlencoded\"> \n");
	printGIorGIs(udp,fp);
	fprintf(fp,"<INPUT TYPE=\"hidden\" NAME=\"disp\" VALUE=\"%d\"> \n",disp_format);
	fprintf(fp,"<a href=\"javascript:document.fasta1.submit()\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to save the sequence alignment in a file (FASTA format with gaps)';\
return true\" class=\"GUTTER2\">FASTA format with gaps</a>");
	fprintf(fp,"</FORM> \n");

	/*form 2 for FASTA*/
	disp_format=DISP_FASTA_NOGAP;
	if (udp->disp_format&ALIGN_REBUILD) disp_format|=ALIGN_REBUILD;
	fprintf(fp,"<FORM NAME=\"fasta2\" ACTION=\"%s\" METHOD=\"GET\" \n",WWWDDV_save_script1);
	fprintf(fp,"ENCTYPE=\"application/x-www-form-urlencoded\"> \n");
	printGIorGIs(udp,fp);
	fprintf(fp,"<INPUT TYPE=\"hidden\" NAME=\"disp\" VALUE=\"%d\"> \n",disp_format);
	fprintf(fp,"<a href=\"javascript:document.fasta2.submit()\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to save the sequence alignment in a file (FASTA format without gaps)';\
return true\" class=\"GUTTER2\">FASTA format without gaps</a>");
	fprintf(fp,"</FORM>");
}

/******************************************************************************/
extern void PrintPopContHeader(URLDataPtr udp,FILE *fp)
{
	fprintf(fp,"<!--  the contents   --> \n");
	fprintf(fp,"<table border=0 width=615 cellspacing=0 cellpadding=0>\n");
	fprintf(fp,"  <tr valign=\"TOP\"> <!--  left column   --> \n");
	fprintf(fp,"    <td width=145> \n");
	fprintf(fp,"	<img src=\"spacer10.GIF\" width=145 height=1 border=0>\n");
	DDV_ToolsScript(udp->sep,udp->disp_format,fp,TRUE);
	fprintf(fp,"	</td><td>");
}

/******************************************************************************/
extern void PrintPopSizeAliTable(Int4 sapLength,Int4 sapNumBioseqs,FILE *fp)
{
	/*display the size of the SeqAlign*/
	fprintf(fp,"<TABLE border=1 cellpadding=6 cellspacing=0 bgcolor=#e6eeee>\n");
	fprintf(fp,"<tr valign=\"top\">\n");
    fprintf(fp,"<td class=\"SMALL1\" bgcolor=#CCCCFF>Alignment Length</td>\n");
	fprintf(fp,"<td class=\"SMALL1\">%d</td></tr>",sapLength);
	fprintf(fp,"<tr valign=\"top\">\n");
    fprintf(fp,"<td class=\"SMALL1\" bgcolor=#CCCCFF># of Sequences</td>\n");
	fprintf(fp,"<td class=\"SMALL1\">%d</td></tr>",sapNumBioseqs);
	fprintf(fp,"</TD></TR></TABLE>\n");
}

/******************************************************************************/
extern void PrintPopRightContent(URLDataPtr udp,Int4 sapLength,Int4 sapNumBioseqs,FILE *fp,
		Int4 numBlockAffich,Int2 LineSize)
{

Char szWWW[3*WWW_SCRIPT_SIZE]={""};
Int4 i,j,k,min,max;
Uint4 disp_format;
		
	fprintf(fp,"<!--  the contents   --> \n");
	fprintf(fp,"<table border=0 width=615 cellspacing=0 cellpadding=0>\n");
	fprintf(fp,"  <tr valign=\"TOP\"> <!--  left column   --> \n");
	fprintf(fp,"    <td width=125> \n");
	fprintf(fp,"<img src=\"spacer10.GIF\" width=125 height=1 border=0>\n");

	/*display the size of the SeqAlign*/
	PrintPopSizeAliTable(sapLength,sapNumBioseqs,fp);

	fprintf(fp,"<br><span class=\"GUTTER1\">Display Options:</span><br>\n");

	/*colour code*/
	if (udp->disp_format&DISPE_COLOR) disp_format=udp->disp_format^DISPE_COLOR;
	else disp_format=udp->disp_format|DISPE_COLOR;

	if (!(udp->disp_format&ALIGN_REBUILD))
		sprintf(szWWW,WWWDDV_script,udp->gi,udp->from+1,udp->to+1,disp_format);
	else
		sprintf(szWWW,WWWDDV_script2,udp->szGIs,disp_format);
	fprintf(fp,"<a href=\"%s\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to use sequence colour codes';return true\" \
class=\"GUTTER2\"> Colour Display</a>\n",szWWW);
	if (udp->disp_format&DISPE_COLOR) fprintf(fp,"<span class=\"GUTTER1\">  x</span>\n<BR>");
	else fprintf(fp,"<br>");
	/*fprintf(fp,"<br>");*/

	/*display the view style Full Seq*/
	if (udp->disp_format&VIEW_VARIA){
		disp_format=udp->disp_format^VIEW_VARIA;
		disp_format|=VIEW_FULLSEQ;
	}else if (udp->disp_format&VIEW_FULLSEQ){
		disp_format=udp->disp_format;
	}

	if (!(udp->disp_format&ALIGN_REBUILD))
		sprintf(szWWW,WWWDDV_script,udp->gi,udp->from+1,udp->to+1,disp_format);
	else
		sprintf(szWWW,WWWDDV_script2,udp->szGIs,disp_format);
		
	fprintf(fp,"<a href=\"%s\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to view the full sequences';return true\" \
class=\"GUTTER2\">\n",szWWW);
	fprintf(fp,"Full Sequences</a> \n");
	if (udp->disp_format&VIEW_FULLSEQ){
		fprintf(fp,"<span class=\"GUTTER1\">&gt;&gt;</span>\n");
	}
	fprintf(fp,"<br>");
	/*display the view style Varia*/
	if (udp->disp_format&VIEW_FULLSEQ){
		disp_format=udp->disp_format^VIEW_FULLSEQ;
		disp_format|=VIEW_VARIA;
	}else if (udp->disp_format&VIEW_VARIA){
		disp_format=udp->disp_format;
	}
	if (!(udp->disp_format&ALIGN_REBUILD))
		sprintf(szWWW,WWWDDV_script,udp->gi,udp->from+1,udp->to+1,disp_format);
	else
		sprintf(szWWW,WWWDDV_script2,udp->szGIs,disp_format);
	fprintf(fp,"<a href=\"%s\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to view only the variations';return true\" class=\"GUTTER2\">\n",szWWW);
	fprintf(fp,"Variations</a> \n");
	if (udp->disp_format&VIEW_VARIA){
		fprintf(fp,"<span class=\"GUTTER1\">&gt;&gt;</span>\n");
	}
	fprintf(fp,"<br><br>");

	/*display block style*/
	if (!(udp->disp_format&ALIGN_REBUILD)){
		
		fprintf(fp,"<TABLE border=0 cellpadding=0 cellspacing=0>\n");
		sprintf(szWWW,WWWDDV_script,udp->gi,1,sapLength,udp->disp_format);
		fprintf(fp,"<tr><td><a href=\"%s\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to view the full sequence alignment on a single page';\
return true\" class=\"GUTTER2\">\n",szWWW);
		fprintf(fp,"All Blocks</a></td>\n");
		if (numBlockAffich==0) fprintf(fp,"<td><span class=\"GUTTER1\">&gt;&gt;</span>\n");
		else fprintf(fp,"<td>");
		fprintf(fp,"</td></tr>");

		fprintf(fp,"<tr><td>");
		fprintf(fp,"<span class=\"GUTTER2\">By Block :</span></td><td></td></tr>\n");
		i=1;
		fprintf(fp,"<tr><td></td><td>");

		k=sapLength/LineSize;
		if (sapLength%LineSize) k++;

		for(j=0;j<k;j++,i++){
			min=j*LineSize+1;
			max=_min_((j+1)*LineSize,sapLength);
			sprintf(szWWW,WWWDDV_script,udp->gi,min,max,udp->disp_format);
			fprintf(fp,"<a href=\"%s\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to view this region';return true\" \
class=\"GUTTER2\">   [%d..%d]</a>\n",szWWW,min,max);
			if (numBlockAffich==i) fprintf(fp,"<span class=\"GUTTER1\">  &gt;&gt;</span>\n");
			if (j<(k-1)) fprintf(fp,"<br>");
		}
		fprintf(fp,"</td></tr></table>\n");
	}	
	
	/*tool bar*/
	fprintf(fp,"<BR>");
	if (!(udp->disp_format&ALIGN_REBUILD))
		DDV_ToolsScript(udp->sep,udp->disp_format,fp,TRUE);

	/*save script*/
	PrintSaveScript(udp,fp);
	
	
	fprintf(fp,"    <!-- right content column  -->\n"); 
	fprintf(fp,"    <!-- extra column to force things over the gif border -->\n"); 
	fprintf(fp,"    <td width=15><img src=spacer10.GIF width=15 height=1 border=0>\n");
	fprintf(fp,"    </td>\n");
	fprintf(fp,"    <td border=0> \n");

	return;
}

/******************************************************************************/
extern void PrintPopulationTail(Uint4 type,FILE *fp)
{
	if ((type&DISP_PHYLIP_TXT) || (type&DISP_FASTA_NOGAP) || 
			(type&DISP_FASTA_GAP) || (type&DISP_FULL_TXT)){
		return;
	}
	fprintf(fp,"<HR width=\"90%\" noshade>\n");
	fprintf(fp,"<a href=\"#top\">[Top]</a>");
	fprintf(fp,"<font size=-1><ADDRESS>Comments and suggestions to: [<A HREF=\"mailto:info@ncbi.nlm.nih.gov\">info@ncbi.nlm.nih.gov</A>]</ADDRESS></font></td></tr></table>\n");
	fprintf(fp,"<!--  end of content  -->\n"); 
	fprintf(fp,"</table>\n");
	fprintf(fp,"</body>\n");
	fprintf(fp,"</html>\n");
}

/******************************************************************************/
extern void PrintPopulationHeadG(Uint4 type,FILE *fp)
{
	if ((type&DISP_PHYLIP_TXT) || (type&DISP_FASTA_NOGAP) || 
			(type&DISP_FASTA_GAP) || (type&DISP_FULL_TXT)){
/*		fprintf(fp,"Content-type: application/octet-stream\n\n");*/
		fprintf(fp,"Content-type: SeqAlign\n\n");
		return;
	}
	fprintf(fp,"Content-type: text/html\n\n");
	fprintf(fp,"<html>\n");
	fprintf(fp,"<head>\n");
	fprintf(fp,"  <title>Population</title>\n");
	fprintf(fp,"  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n");
	fprintf(fp," <META NAME=\"keywords\" CONTENT=\"insert your keywords for the search engine\">\n");
	fprintf(fp," <META NAME=\"description\" CONTENT=\"insert the description to be displayed by the search engine.  Also searched by the search engine.\">\n");
	fprintf(fp,"  <link rel=\"stylesheet\" href=\"http://www.ncbi.nlm.nih.gov/corehtml/ncbi.css\">\n");
	fprintf(fp,"</head>\n");
}

/******************************************************************************/
extern void PrintPopulationBar(CharPtr db,FILE *fp)
{
	fprintf(fp,"<body bgcolor=\"#FFFFFF\" background=\"http://www.ncbi.nlm.nih.gov/corehtml/bkgd.gif\" text=\"#000000\" link=\"#000099\" vlink=\"#6666CC\">\n");
	fprintf(fp,"<a name=\"top\"></a>");
	fprintf(fp,"<!--  the header   -->\n"); 
	fprintf(fp,"<TABLE border=0><TR><TD>\n");
	fprintf(fp,"<table border=0 width=600 cellspacing=0 cellpadding=0>\n");
	fprintf(fp,"   <tr> <td width=140><a href=http://www.ncbi.nlm.nih.gov> <img src=\"http://www.ncbi.nlm.nih.gov/corehtml/left.GIF\" width=130 height=45 border=0></a></td>\n");
	/*if (StringCmp(db, "n") == 0) {
		fprintf(fp,"    <td width=360 class=\"head1\" valign=\"BOTTOM\"> <span class=\"H1\">Entrez Nucleotides</span></td>\n");
	} else if (StringCmp(db, "p") == 0) {
		fprintf(fp,"    <td width=360 class=\"head1\" valign=\"BOTTOM\"> <span class=\"H1\">Entrez Proteins</span></td>\n");
	}*/
	fprintf(fp,"    <td align=\"left\"> <IMG SRC=\"/entrez/query/static/gifs/entrez_popset.gif\"></td>\n");
	fprintf(fp,"  </tr>\n");
	fprintf(fp,"</table>\n");
	fprintf(fp,"<TD></TR><TR><TD>\n");
	fprintf(fp,"<!--  the quicklinks bar   -->\n"); 
	fprintf(fp,"<table CLASS=\"table-text\" border=0 width=600 cellspacing=0 cellpadding=3 bgcolor=\"#000000\">\n");
	fprintf(fp,"  <tr><td colspan=3>\n");
	fprintf(fp,"<table border=0 width=600 cellspacing=0 cellpadding=1>\n");
	fprintf(fp,"<tr bgcolor=\"#000000\">  \n"); 
	fprintf(fp,"    <td width=85><a href=http://www.ncbi.nlm.nih.gov/BLAST/ class=\"BAR\">BLAST</a></td>\n");
	fprintf(fp,"   <td width=85><a href=\"/entrez/query.fcgi?db=PubMed\" class=\"BAR\">PubMed</a></td>\n");
	fprintf(fp,"    <td width=85><a href=\"/entrez/query.fcgi?db=Nucleotide\" class=\"BAR\">Nucleotide</a></td>\n");
	fprintf(fp,"    <td width=85><a href=\"/entrez/query.fcgi?db=Protein\" class=\"BAR\">Protein</a></td>\n");
	fprintf(fp,"    <td width=85><a href=\"/entrez/query.fcgi?db=Genome\" class=\"BAR\">Genome</a></td>\n");
	fprintf(fp,"    <td width=85><a href=\"/entrez/query.fcgi?db=Structure\" class=\"BAR\">Structure</a></td>\n");
	fprintf(fp,"    <td width=85><a href=\"/entrez/query.fcgi?db=Popset\"  class=\"BAR\">PopSet</a></td>\n");
	fprintf(fp,"    <td width=85><a href=\"http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html\"  class=\"BAR\">Taxonomy</a></td>\n");
	fprintf(fp,"  </tr></table></td></tr>\n");
	fprintf(fp,"   </table></TD></TR></TABLE>\n");
}
/*******************************************************************************

  Function : DDV_GetFirstAuthor()
  
  Purpose : retrieve the first author of the PopSet
  
  Parameters : 	alp; valnode of authors
  				name; filled in this function
  
  Return value : - (see 'name')

*******************************************************************************/
static void DDV_GetFirstAuthor (AuthListPtr alp, CharPtr name)
{
ValNodePtr   names;
AuthorPtr    ap;
PersonIdPtr  pid;

	if (!alp) return;
	names = alp->names;/*this is a ValNodeList, but I only need to deal
	with the first node, i.e the first author*/
	if (alp->choice == 1) {
		ap = (AuthorPtr) names->data.ptrvalue;
		if (ap != NULL) {
			pid = ap->name;
			if (pid != NULL) {
				PersonIdLabel(pid, name, INFO_SIZE, PIDLABEL_GENBANK);
			}
		}
	}
	else if (alp->choice == 2 || alp->choice == 3) {
		StringNCpy(name,(CharPtr)names->data.ptrvalue,INFO_SIZE);
	}
}

/*******************************************************************************

  Function : DDV_AnalysePubDesc()
  
  Purpose : analyse a PubDesc pointer and get Author name and title
  
  Parameters : 	pdp; PubDesc pointer
  				szPopSetName; fills here with a title
				szPopSetAuth;fills here with the first author name
  
  Return value : - (see 'name')

*******************************************************************************/
static void	DDV_AnalysePubDesc(PubdescPtr pdp, CharPtr szPopSetName,
		CharPtr szPopSetAuth,Int4Ptr pmid)
{
ValNodePtr   vnp,vnp2;
CitArtPtr    cap;
CitGenPtr    cgp;
CitSubPtr    csp;
CitPatPtr    cpp;
CitBookPtr   cbp;
MedlineEntryPtr mep;

	/*tests for debug only*/
	/*PubLabel(pdp->pub, szPopSetAuth, INFO_SIZE, OM_LABEL_CONTENT);*/
	/*PubLabel (pdp->pub, szPopSetName, INFO_SIZE, OM_LABEL_SUMMARY);*/
	for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
		switch(vnp->choice){
			case PUB_Article:/*article*/
				cap = (CitArtPtr) vnp->data.ptrvalue;
				if (!cap) break;
				for (vnp2 = cap->title; vnp2 != NULL; vnp2 = vnp2->next) {
					if(vnp2->choice==Cit_title_name){
						StringNCpy(szPopSetName,(CharPtr)vnp2->data.ptrvalue,INFO_SIZE);
					}
				}
				DDV_GetFirstAuthor (cap->authors, szPopSetAuth);
				break;
			case PUB_Gen:/*gen*/
				cgp = (CitGenPtr) vnp->data.ptrvalue;
				if (!cgp) break;
				if (cgp->title && *(cgp->title))
					StringNCpy(szPopSetName,cgp->title,INFO_SIZE);
				DDV_GetFirstAuthor (cgp->authors, szPopSetAuth);
				break;
			case PUB_Sub:/*Sub*/
				csp = (CitSubPtr) vnp->data.ptrvalue;
				if (!csp) break;
				if (csp->descr && *(csp->descr))
					StringNCpy(szPopSetName,csp->descr,INFO_SIZE);
				DDV_GetFirstAuthor (csp->authors, szPopSetAuth);
				break;
			case PUB_Patent:/*Pat*/
				cpp = (CitPatPtr) vnp->data.ptrvalue;
				if (!cpp) break;
				if (cpp->title && *(cpp->title))
					StringNCpy(szPopSetName,cpp->title,INFO_SIZE);
				DDV_GetFirstAuthor (cpp->authors, szPopSetAuth);
				break;
			case PUB_Book:/*Book*/
				cbp = (CitBookPtr) vnp->data.ptrvalue;
				if (!cbp) break;
				for (vnp2 = cbp->title; vnp2 != NULL; vnp2 = vnp2->next) {
					if(vnp2->choice==Cit_title_name){
						StringNCpy(szPopSetName,(CharPtr)vnp2->data.ptrvalue,INFO_SIZE);
					}
				}
				DDV_GetFirstAuthor (cbp->authors, szPopSetAuth);
				break;
			case PUB_Medline:/*Medline*/
				mep = (MedlineEntryPtr) vnp->data.ptrvalue;
				if (!mep) break;
				cap = mep->cit;
				for (vnp2 = cap->title; vnp2 != NULL; vnp2 = vnp2->next) {
					if(vnp2->choice==Cit_title_name){
						StringNCpy(szPopSetName,(CharPtr)vnp2->data.ptrvalue,INFO_SIZE);
					}
				}
				DDV_GetFirstAuthor (cap->authors, szPopSetAuth);
				break;
			case PUB_PMid:/*PubMed Id*/
				*pmid=vnp->data.intvalue;
				break;
		}
	}
}

/*******************************************************************************

  Function : DDV_GetArticleInfo()
  
  Purpose : retrieve the article's title of the PopSet
  
  Parameters : 	sep; seqentry
  				szPopSetName; filled in this function
  				szPopSetAuth; filled in this function

  Note:
	A BioseqSetPtr contains the field 'descr'. This is a ValNodeList.
	Among the nodes, one is of type Seq_descr_pub, the PubMed data information.
	In this case, descr->data.ptrvalue if of type PubdescPtr. Let say
	pdp=(PubdescPtr)descr->data.ptrvalue. Now, pdp is a ValNode List of 
	PubMed fields where pdp->choice can be PUB_Article. In that case
	pdp->data.ptrvalue if of type CitArtPtr. This data block contains
	the article info: title and authors.

	This function is designed to retrieve the first title and the first
	author correponding to that title.

	I've enhanced this function to manage Gen, Sub, Pat, Book, Medline.
  
  Return value : - (see 'szPopSetName' & 'szPopSetAuth')

*******************************************************************************/
extern void	DDV_GetArticleInfo(SeqEntryPtr sep,CharPtr szPopSetName,
	CharPtr szPopSetAuth,Int4Ptr pmid)
{
BioseqSetPtr bssp;
PubdescPtr   pdp;
ValNodePtr   vnp3;

	if (!sep || !szPopSetName || !szPopSetAuth) return;
	
	*szPopSetName='\0';
	*szPopSetAuth='\0';
	*pmid=0;
	
	if (IS_Bioseq_set (sep)) {/*a Pop-set is a BSSP*/
		bssp = (BioseqSetPtr) sep->data.ptrvalue;		
		if (bssp){
			if(bssp->descr != NULL) {
				for (vnp3=bssp->descr ; vnp3 != NULL ; vnp3 = vnp3->next){
					if(vnp3->choice==Seq_descr_pub){
    					pdp = (PubdescPtr) vnp3->data.ptrvalue;
						if (pdp){
							DDV_AnalysePubDesc(pdp,szPopSetName,
								szPopSetAuth,pmid);
						}
					}
				}
			}
		}
	}
	if (*(szPopSetName)=='\0') StringCpy(szPopSetName,"* no title *");
	if (*(szPopSetAuth)=='\0') StringCpy(szPopSetAuth,"* no authors *");
}

/*******************************************************************************

  Function : DDV_AnalyseBiosource()
  
  Purpose : analyse a BioSource pointer and get TaxName and common name
  
  Parameters : 	biop; BioSource pointer
  				szTaxName; fills here with a TaxName
				szCommonName;fills here with the first common name
  
  Return value : - (see 'name')

*******************************************************************************/
static void	DDV_AnalyseBiosource(BioSourcePtr biop, CharPtr PNTR szTaxName,
		CharPtr PNTR szCommonName)
{

	if (!biop) return;

	if (biop->org && biop->org->taxname){
		*szTaxName=biop->org->taxname;
	}
	if (biop->org && biop->org->common){
		*szCommonName=biop->org->common;
	}
	
}

/*******************************************************************************

  Function : DDV_GetEntryBioSource()
  
  Purpose : retrieve the BioSource, given a sep
  
  Note : this function can be called several times to fill the valnode vnpp
  with all the Biosources of all SeqSet of an Entry. This is used, for example,
  by the PopSet Viewer to build the organim names list.
  
  Parameters : 	sep; seqentry
  				vnpp; list to fill. data.ptrvalue if of type PopSourcePtr (see pgppop.h)

*******************************************************************************/

extern void	DDV_GetEntryBioSource(SeqEntryPtr sep,ValNodePtr PNTR vnpp)
{
BioseqSetPtr bssp;
BioSourcePtr biop;
ValNodePtr   vnp;
PopSourcePtr psp;

	if (!sep || !vnpp) return;
	
	if (IS_Bioseq_set (sep)) {/*a Pop-set is a BSSP*/
		bssp = (BioseqSetPtr) sep->data.ptrvalue;		
		if (bssp){
			if(bssp->descr != NULL) {/*scan the 'desc' field of the bssp*/
				for (vnp=bssp->descr ; vnp != NULL ; vnp = vnp->next){
					if(vnp->choice==Seq_descr_source){/*BioSource ?*/
    					biop = (BioSourcePtr) vnp->data.ptrvalue;
						if (biop){
							psp=(PopSourcePtr)MemNew(sizeof(PopSource));
							if (psp){
								DDV_AnalyseBiosource(biop,&(psp->szTaxName),
									&(psp->szCommonName));
								if (!(*vnpp)) *vnpp=ValNodeAddPointer(NULL,0,(Pointer)psp);
								else ValNodeAddPointer(vnpp,0,(Pointer)psp);
							}
						}
					}
				}
			}
		}
	}
}

/*******************************************************************************

  Function : SearchBioSource()
  
  Purpose : caalback for SeqEntryExplore to get all the biosources of all set
  in a popset entry, fo example
  
  Parameters : see Toolbox doc !
  
  Return value : none 

*******************************************************************************/
static void SearchBioSource(SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)
{
ValNodePtr PNTR	vnpp;

	vnpp=(ValNodePtr PNTR)mydata;

	if (IS_Bioseq_set (sep)) {
		DDV_GetEntryBioSource(sep, vnpp);
	}
}


/*******************************************************************************

  Function : DDV_SearchBioseq_2()
  
  Purpose : store bsps of a SeqEntry ; release for DDV_ToolsScript (see below)
  
  Parameters : see Toolbox doc !
  
  Return value : none 

*******************************************************************************/
static Boolean LIBCALLBACK DDV_SearchBioseq_2 (BioseqPtr bsp, 
			SeqMgrBioseqContextPtr context)
{
ValNodePtr PNTR	BspTable;
Int4 		    uid=0;

	BspTable=(ValNodePtr PNTR)context->userdata;

	if (bsp) {
		uid=GetGIForSeqId(bsp->id);
		if (uid!=0){
			ValNodeAddInt(BspTable,1,uid);
		}
	}

	return(TRUE);
}


/*******************************************************************************

  Function : DDV_SearchBioseq_1()
  
  Purpose : store bsps of a SeqEntry ; release for DDV_GetBspList (see below)
  
  Parameters : see Toolbox doc !
  
  Return value : none 

*******************************************************************************/
static Boolean LIBCALLBACK DDV_SearchBioseq_1 (BioseqPtr bsp, 
			SeqMgrBioseqContextPtr context)
{
/*Int4Ptr pCompt;*/
Char szSeqAcc[15]={""};
Char szSeqName[100]={""};
Char szBuf4[WWW_SCRIPT_SIZE]={""};	/*Entrez query*/
Int4 gi;
SeqIdPtr sip;
GetBspPopPtr gbpp;	

	gbpp=(GetBspPopPtr)context->userdata;

	gbpp->nCompt++;
	
	/*prepare a HTML display*/
	fprintf(gbpp->fp,"[%4i] ",gbpp->nCompt);
	/*access code*/
	SeqIdWrite(bsp->id, szSeqAcc, PRINTID_TEXTID_ACCESSION, 14);
	/*seq title*/
	FastaDefLine(bsp, szSeqName, 99, NULL, NULL, 0);
	/*Gi number*/
	sip=SeqIdFindBest(bsp->id,0);
	if (sip==NULL) sip=bsp->id;
	if (sip){
		gi=GetGIForSeqId(sip);
		if(gi>0){
			sprintf(szBuf4,szEntrezScript,gi,(ISA_aa(bsp->mol) ? "p" : "n"));
		}
		else{
			sprintf(szBuf4,"javascript:void(0)");
		}
		if (StringLen(szSeqName)){
                Int2 j=0;
	     		while(szSeqName[j]){
        			if (szSeqName[j]=='\'') szSeqName[j]=' ';
				    j++;
		       	}
			fprintf(gbpp->fp,szSeqNameStatus,szBuf4,szSeqName);
		}
		else{
			fprintf(gbpp->fp,szSeqNameStatus,szBuf4,"unknown name");
		}
	}
	fprintf(gbpp->fp,"%10s ",szSeqAcc);
	if (sip) fprintf(gbpp->fp,"</a>  ");
	fprintf(gbpp->fp,": %s<br>\n",szSeqName);
	return(TRUE);
}

/*******************************************************************************

  Function : DDV_GetBspList()
  
  Purpose : retrieve (a list of) BSPs (if the SeqEntry doesn't contain any
  		SeqAligns).
  
  Parameters : 	sep; seqentry
  
  Return value : - 

*******************************************************************************/
extern void DDV_GetBspList(SeqEntryPtr sep,FILE *fp)
{

Int2 entityID;
GetBspPop gbp;

	gbp.nCompt=0;
	gbp.fp=fp;
		
	fprintf(fp,"<font color=\"#9933ff\"><h4>Please note that this population study ");
	fprintf(fp,"does not contain any sequence alignments. <br>");
	fprintf(fp,"Here is the list of the sequences reported in the entry:</h4></font><br>\n");
	/*scan for any BSPs*/
	entityID=SeqMgrIndexFeatures(0, (Pointer)sep);
	fprintf(fp,"<pre>\n");
	SeqMgrExploreBioseqs (entityID, NULL, (Pointer)&gbp, 
			DDV_SearchBioseq_1, TRUE, FALSE, TRUE);
	fprintf(fp,"<TABLE border=1 cellpadding=6 cellspacing=0 bgcolor=#e6eeee>\n");
	fprintf(fp,"<tr>");
    fprintf(fp,"<td>Total : %i sequence(s).</td>\n",gbp.nCompt);
	fprintf(fp,"</TD></TR></TABLE>\n");
	fprintf(fp,"</pre>\n");
}

/*******************************************************************************

  Function : DDV_ToolsScript()
  
  Purpose : print out a form to run a SeqAlign rebuild.
  
  Parameters : 	sep; seqentry to analyse
  				disp_format; format display
  				fp; file
				ForToolbar; TRUE only for the PopSet Viewer
				
  Return value : - 

*******************************************************************************/
static void DDV_ToolsScript(SeqEntryPtr sep, Uint4 disp_format,FILE *fp, 
	Boolean ForToolbar)
{
ValNodePtr uid_list=NULL,vnp;
Int2 entityID;
Uint4 disp_format2;

	entityID=SeqMgrIndexFeatures(0, (Pointer)sep);
	SeqMgrExploreBioseqs (entityID, NULL, (Pointer)&uid_list, 
			DDV_SearchBioseq_2, TRUE, FALSE, TRUE);
	if (ForToolbar)
		fprintf(fp,"<span class=\"GUTTER1\">Tools:</span><br>\n");
	
	fprintf(fp,"<FORM NAME=\"realign1\" ACTION=\"%s\" METHOD=\"GET\" \n",WWWDDV_save_script1);
	fprintf(fp,"ENCTYPE=\"application/x-www-form-urlencoded\"> \n");
	fprintf(fp,"<INPUT TYPE=\"hidden\" NAME=\"gis\" VALUE=\"");
	vnp=uid_list;
	while(vnp){
		fprintf(fp,"%lda",(Int4)vnp->data.intvalue);
		vnp=vnp->next;
	}
	fprintf(fp,"\">\n");
/*	disp_format2=disp_format;*/
	disp_format2=RULER_TOP;
	disp_format2|=RULER_TICK;
	disp_format2|=DISPE_SHOWBLOCK;
	disp_format2|=VIEW_VARIA;
	disp_format2|=DISP_FULL_HTML;
	disp_format2|=DISP_ORDER_NUM; 
	disp_format2|=DISP_STRAND;
	disp_format2|=DISP_BSP_COORD;
	disp_format2|=ALIGN_REBUILD;

	fprintf(fp,"<INPUT TYPE=\"hidden\" NAME=\"disp\" VALUE=\"%d\"> \n",disp_format2);
	fprintf(fp,"<a href=\"javascript:document.realign1.submit()\" onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to build a new sequence alignment';\
return true\" %s>%s</a>",(ForToolbar ? "class=\"GUTTER2\"" : " "),
	(ForToolbar ? "Align" : "- build a master-slave sequence alignment."));
	fprintf(fp,"</FORM> \n");
	
}

/*******************************************************************************

  Function : DDV_DisplayDefaultAlign()
  
  Purpose : build a default SeqAlign display.
  
  Parameters : 	sap; SeqAlign
				gi;gi number - used by PopSet viewer to build the HTTP links
				from - to; range of the Sap to display
				disp_format; display formats (see pgppop.h)
				disp_data; for future implementation
				fp; to send the output
  Return value : FALSE if failed 

*******************************************************************************/
extern Boolean DDV_DisplayDefaultAlign(SeqAlignPtr sap,Int4 gi,Int4 from,Int4 to,
		Uint4 disp_format,Pointer disp_data,FILE *fp)
{
SABlockPtr		sabp=NULL;/*indexed SeqAlign*/
Int4			sapLength, sapNumBioseqs, numBlocks;/*indexed SeqAlign size*/
MsaParaGPopList mppl;/*contains the data for the display*/
Int4			numBlockAffich=0;
Int2            LineSize;
DDVOptionsBlockPtr dobp;
ByteStorePtr   PNTR byteSpp;

	MemSet(&mppl,0,sizeof(MsaParaGPopList));

	/*build the indexed SeqAlign*/
	sabp=Blockify(sap);
	if (sabp==NULL) return(FALSE);

	/*retrieve the size of the alignment*/
	GetBlockInfo(sabp, &sapLength, &sapNumBioseqs,&numBlocks,NULL);
	
	if (from==-1) from=0;
	if (to==-1) to=sapLength-1;

	if (disp_format&DISP_FASTA_NOGAP || /*FASTA output*/
			disp_format&DISP_FASTA_GAP){/*don't need a ParaG structure*/
		DDV_GetFullFASTAforIdxAli(sabp,disp_format,sapLength,fp);
	}
	else{/*HTML and PHYLIP: need a ParaG structure*/
		/*build the ParaG list for the display*/
		if (disp_data){
			dobp=(DDVOptionsBlockPtr)disp_data;
			LineSize=dobp->LineSize;
			byteSpp=dobp->byteSpp;
		}
		else{
			LineSize=ParaG_Size;
			byteSpp=NULL;
		}
		if (!DDV_PopDisplay2(sabp,&mppl,sapNumBioseqs,SCALE_POS_TOP,
				TRUE,LineSize,numBlocks,from,to)) goto error;
				/*ParaG_Size is defined in pgppop.h*/
		/*show by block ?*/
		if ((to-from+1)<sapLength) {
			numBlockAffich=from/LineSize+1;
		}		

		/*Create the display */
		DDV_AffichageParaG(&mppl,gi,from,to,sapLength,numBlockAffich,disp_format,
			LineSize,fp,byteSpp);

		/*done... delete data*/
		DDV_DeleteDisplayList(&mppl);
	}
error:
	if (sabp) CleanupBlocks(sabp);
	
	return(TRUE);
}

/*******************************************************************************

  Function : DDV_ReadSeqGivenSpp()
  
  Purpose : read a sequence buffer.
  
  Parameters : 	spp; open seqport
  				from-to; range to read
				strand; one of the Seq_strand_... (see objloc.h)
				IsAA; TRUE if the sequence is a protein
				bError; set to TRUE only if an error occurs when reading the seq.
				
  Return value : a buffer with the sequence 

*******************************************************************************/
extern CharPtr DDV_ReadSeqGivenSpp(SeqPortPtr spp,Int4 from,Int4 to,Uint1 strand,
	Boolean IsAA,BoolPtr bError)
{
Int4 len,seek,i;
CharPtr szSeq=NULL;

Byte residue;

	if (!spp) return (NULL);
	if (from>to) return (NULL);

	len=to-from+1;
	szSeq=(CharPtr)MemNew((len+1)*sizeof(Char));		
	if (!szSeq) return (NULL);
	if(IsAA) {
		seek=from;
	} 
	else {
		if(strand == Seq_strand_minus)
			seek = to;
		else
			seek = from;
	}
	SeqPortSeek(spp, seek, SEEK_SET);
	/*ret=SeqPortRead (spp,(Uint1Ptr)szSeq, (Int2)len);*/
	for(i=0;i<len;i++){
		residue = SeqPortGetResidue(spp);
		if (residue!=SEQPORT_EOF && IS_residue(residue)) {
			szSeq[i] = (Char)residue;
		}
	}
/*	if (ret!=len) *bError=TRUE;
	else *bError=FALSE;*/

	return(szSeq);
}


/*******************************************************************************

  Function : DDV_OpenBspFullSeqPort()
  
  Purpose : open a seqport of a full sequence.
  
  Parameters : 	bsp; bioseq (be sure it's locked)
				strand; one of the Seq_strand_... (see objloc.h)

  Return value : an open seqport 

*******************************************************************************/
extern SeqPortPtr DDV_OpenBspFullSeqPort(BioseqPtr bsp,Uint1 strand)
{
SeqPortPtr spp;

	if (!bsp) return (NULL);
		
	if(bsp->mol!=Seq_mol_aa) {
    	spp= SeqPortNew(bsp, (Int4)0,(Int4) bsp->length-1, 
			strand, Seq_code_iupacna);
	}
	else {
    	spp= SeqPortNew(bsp, (Int4)0,(Int4) bsp->length-1, 
			0, Seq_code_ncbistdaa);
	}
	return(spp);
}

/*******************************************************************************

  Function : DDV_BspInfoNew()
  
  Purpose : BspInfo constructor.
  
  Parameters : 	InitInfo; TRUE means init the 'ascii' fields of BspInfo
  				InitSeqPort; TRUE means open a seqport
				sip; sequence identifier
				strand; one of the Seq_strand_... (see objloc.h)

  Return value : a pointer to an initialized BspInfo structure 

*******************************************************************************/
extern BspInfoPtr DDV_BspInfoNew(Boolean InitInfo,Boolean InitSeqPort,SeqIdPtr sip,
		Uint1 strand)
{
BspInfoPtr bip;


	if (!sip) return(NULL);

	/*alloc*/
	bip=(BspInfoPtr)MemNew(sizeof(BspInfo));
	if (!bip) return(NULL);
	
	/*lock the bsp in memory*/
	bip->bsp=BioseqLockById(sip);
	if (!bip->bsp) return(NULL);
	bip->sip=sip;
	/*init Bsp info fields*/
	if (InitInfo){
		UDV_ReadBspDataForViewer(bip);
	}
	
	/*open the seqport on the full sequence*/
	if (InitSeqPort){
		bip->spp=DDV_OpenBspFullSeqPort(bip->bsp,strand);
	}
	
	return(bip);
}


/*******************************************************************************

  Function : DDV_BspInfoDelete()
  
  Purpose : BspInfo destructor.
  
  Parameters : 	bip; a pointer to an initialized BspInfo structure

  Return value : NULL, always 

*******************************************************************************/
extern BspInfoPtr DDV_BspInfoDelete(BspInfoPtr bip)
{ 
	if(bip) {
		if (bip->spp) SeqPortFree(bip->spp);
		if (bip->bsp) BioseqUnlock(bip->bsp);
		if (bip->bspGeneticCode) MemFree(bip->bspGeneticCode);
		if (bip->SeqBuf) MemFree(bip->SeqBuf);
		MemFree(bip);
	}
	return(NULL);
}

/*******************************************************************************

  Function : DDV_BspInfoDeleteList()
  
  Purpose : BspInfo list destructor.
  
  Parameters : 	bip; a pointer to an initialized BspInfo structure

  Return value : NULL, always 

*******************************************************************************/
extern BspInfoPtr DDV_BspInfoDeleteList(BspInfoPtr bip)
{
BspInfoPtr next;

	while (bip != NULL){
		next = bip->next;
		DDV_BspInfoDelete(bip);
		bip = next;
	}
	return(NULL);
}

/*******************************************************************************

  Function : DDV_LocateParaG()
  
  Purpose : compute StartLine valies for PGP in a Full Horz display style 
  
  Parameters : -
  
  Return value : -

*******************************************************************************/
extern void DDV_LocateParaG(ValNodePtr PNTR TableHead,Int4 nBsp)
{
ParaGPtr pgp;
Int4 StartLine=0,StartLine2=0,i;
ValNodePtr vnp;
Boolean bFirst;
	
	if (TableHead==NULL || nBsp==0) return;
	
	for(i=0;i<nBsp;i++){
		StartLine=StartLine2;
		bFirst=TRUE;
		for(vnp=TableHead[i];vnp!=NULL;vnp=vnp->next){
			pgp=(ParaGPtr)(vnp->data.ptrvalue);
			if (pgp){
				pgp->StartLine=StartLine;
				if (bFirst) {
					StartLine2+=pgp->nLines;
					/*StartLine=StartLine2;*/
					bFirst=FALSE;
				}
			}
		}
	}	
}

/**** DEBUG only ****/
#ifdef DEBUG_DDV
static void display_ParaG_content(MsaParaGPopListPtr mpplp)
{
BioseqPtr bsp_debug=NULL;
Char 			szBuf[15];
ParaGPtr 		pgp;
MsaTxtDispPtr mtdp;
ValNodePtr    vnp,vnp2;
Uint4 i,j,k,size;
FILE *fp;
Boolean bClose=TRUE;

	fp=fopen("zpgplist.txt","w");
	if (fp==NULL) {fp=stdout;bClose=FALSE;}
	for(i=0;i<mpplp->nBsp;i++){
		/*for each line (i.e. bioseq) */
		vnp=mpplp->TableHead[i];
		if (vnp) pgp=(ParaGPtr)(vnp->data.ptrvalue);
		else pgp=NULL;
		if(pgp){
			if (pgp->sip) bsp_debug=BioseqLockById(pgp->sip);
			else bsp_debug=NULL;
			if (bsp_debug){
				/*get the access code*/
				SeqIdWrite(bsp_debug->id, szBuf, PRINTID_TEXTID_ACCESSION, 14);
				fprintf(fp,"[%5u] %s , size: %d, IsAA :%d \n",
					i,szBuf,BioseqGetLen(bsp_debug),
					ISA_aa(bsp_debug->mol));
				BioseqUnlock(bsp_debug);
			}
			else StringCpy(szBuf,"unknown");
		}
		k=1;size=1;
		while(vnp){
			pgp=(ParaGPtr)(vnp->data.ptrvalue);
			if (pgp){/*create the name table*/
				/*loop on the pgp*/
				fprintf(fp,"  ->ParaG[%d] , range (%7d..%7d):\n",k++,size,size+ParaG_Size-1);
				size+=ParaG_Size;
				vnp2=pgp->ptxtList;
				j=1;
				while(vnp2){
					mtdp=(MsaTxtDispPtr)vnp2->data.ptrvalue;
					if (mtdp){
					fprintf(fp,"    (%4u): range(%7d..%7d) , Gap: %2d, IDs (%4d,%4d) , Strand(%d)\n",
						j++,mtdp->from,mtdp->to,mtdp->IsGap,mtdp->SegID,mtdp->BspID,
						mtdp->strand);
					}
					vnp2=vnp2->next;
				}
			}
			vnp=vnp->next;
		}
		fprintf(fp,"\n");
	}
	if (bClose) fclose(fp);

}
#endif /*DEBUG_DDV*/

/*******************************************************************************

  Function : DDV_CreateDisplay()
  
  Purpose : create a display for DDV panel.
  
  Parameters : 	sap; seqalign
				mpplp; allocated structure filled here

  Return value : TRUE if success

*******************************************************************************/
extern Boolean DDV_CreateDisplay(SeqAlignPtr sap,MsaParaGPopListPtr mpplp)
{
Int4 numBlocks;
Boolean nRet;

	if (!sap) return(FALSE);

	/*build the indexed SeqAlign*/
	mpplp->sabp=Blockify(sap);
	if (!mpplp->sabp) return(FALSE);

	/*retrieve the size of the alignment*/
	GetBlockInfo(mpplp->sabp, &mpplp->LengthAli, &mpplp->nBsp,&numBlocks,NULL);
	
	/*build the ParaG list for the display*/
	nRet=DDV_PopDisplay2(mpplp->sabp,mpplp,mpplp->nBsp,SCALE_POS_NONE,
			TRUE,ParaG_Size,numBlocks,0,mpplp->LengthAli-1);

/*
#ifdef DEBUG_DDV
	display_ParaG_content(mpplp);
#endif	
*/
	DDV_LocateParaG(mpplp->TableHead,mpplp->nBsp);
	
	return(nRet);
}

Boolean DDV_ShowSeqAlign(SeqAlignPtr seqalign, Int4 gi, Int4 from, Int4 to,
                         Uint4 disp_format)
{
    
    SABlockPtr sabp;       /* indexed SeqAlign */
    Int4  sapLength, sapNumBioseqs, numBlocks;/* indexed SeqAlign size */
    Int4  restdiv;         /* used to modify the udp->from value */
    MsaParaGPopList mppl;  /* contains the data for the display */
    Char szPopSetName[INFO_SIZE+1]={""};/*article title of the Pop-Set*/
    /*article first author name of the Pop-Set*/    
    Char szPopSetAuth[INFO_SIZE+1]={""};
    Int4 numBlockAffich=0;
    
    /*seqalign = PairSeqAlign2MultiSeqAlign(seqalign);*/
    
    /*build the indexed SeqAlign*/
    if((sabp = Blockify(seqalign)) == NULL) {
        printf("Error No. 54. Could not show align");
        return FALSE;
    }
    
    /*retrieve the size of the alignment*/
    GetBlockInfo(sabp, &sapLength, &sapNumBioseqs,&numBlocks,NULL);
    
    if ((disp_format & DISP_FASTA_NOGAP) || /*FASTA output*/
        (disp_format & DISP_FASTA_GAP)){ /*don't need a ParaG structure*/
        DDV_GetFullFASTAforIdxAli(sabp, disp_format, sapLength, stdout);
    } else {   /* HTML or Phylip : need to create the ParaG structure */
        
        /* for a nice view, always start as a multiple of LETTER_BLOCK_WIDTH */
        restdiv = from % LETTER_BLOCK_WIDTH;
        from = from-restdiv;
        
        /*user didn't give from-to; from=0 and to=sapLength-1*/
        if (to < 1) to = sapLength-1;
        
        /*some security*/
        
        if (from >= to) 
            return FALSE;
        
        if (from < 0 || to < 0) 
            return FALSE;
        
        /* build the ParaG list for the display */
        DDV_PopDisplay2(sabp, &mppl, sapNumBioseqs, SCALE_POS_TOP,
                        TRUE, ParaG_Size, numBlocks, from, to);
        
        /*show by block ?*/
        if ((to - from + 1) < sapLength) {
            numBlockAffich = from/ParaG_Size+1;
        }	
        
        if (disp_format & DISP_FULL_HTML) {

            /* PrintPopulationBar(NULL, stdout); */ /* PubMed header */
            
            /*name of the output: title of the article*/
            /*DDV_GetArticleInfo(sep, szPopSetName, szPopSetAuth);*/
            /* PrintPopRightContent(szPopSetName,szPopSetAuth, udp, 
               sapLength,sapNumBioseqs, stdout, 
               numBlockAffich); */
        }
        
        /*Create the HTML display for the Population division*/
        DDV_AffichageParaG(&mppl, gi, from, to, sapLength,
                           numBlockAffich, disp_format, ParaG_Size,stdout,NULL);
        /*done... delete data*/
        DDV_DeleteDisplayList(&mppl);
    }
    
    if (sabp) CleanupBlocks(sabp);
    
    return TRUE;
}

/*******************************************************************************

  Function : DDV_DisplayBlastSAP()
  
  Purpose : test function to display a BLAST output.
  
  Parameters : 	sap; seqalign

  Return value : -

*******************************************************************************/
extern void DDV_DisplayBlastSAP(SeqAlignPtr sap)
{
Uint4 option;
FILE  *fp;
SeqAlignPtr sap2;
DDVOptionsBlock dob;


	if (!sap) return;

	fp=fopen("zblastout.txt","w");
	if (!fp) return;
  
  	/*init display format*/
	option =VIEW_FULLSEQ;
	option|=DISP_FULL_TXT;
	option|=DISP_BSP_COORD;
	option|=TEXT_LOWERCASE;
	option|=SEQNAME_BLAST_STD;
	dob.LineSize=(Int2)60;
	
	sap2=sap;
	while(sap2){
		DDV_DisplayDefaultAlign (sap2, 0, -1,-1, option, &dob, fp);
		sap2=sap2->next;
	}
	fclose(fp);

}

/*******************************************************************************

  Function : DDV_SearchAli()
  
  Purpose : callback of SeqEntryExplore from  DDV_GetSeqAlign()
  
  Parameters : see toolkit code...
  
  Return value : none 

*******************************************************************************/
extern void DDV_SearchAli(SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)
{
ValNodePtr PNTR	vnpp;
BioseqPtr     bsp;
BioseqSetPtr  bssp;
SeqAnnotPtr   sap,sap2;
SeqAlignPtr   salp;

	vnpp=(ValNodePtr PNTR)mydata;
	
	/*Bioseq or BioseqSet ?... the SeqAlign (if any) is in the seqannot part*/
	if (IS_Bioseq (sep)) {
		bsp = (BioseqPtr) sep->data.ptrvalue;
		sap = bsp->annot;
	} else if (IS_Bioseq_set (sep)) {
		bssp = (BioseqSetPtr) sep->data.ptrvalue;
		sap = bssp->annot;
	} else return;
	sap2=sap;
	while (sap2 != NULL) {
		if (sap2->type == 2) {/*seqalign*/
			salp=(SeqAlignPtr)sap2->data;
			while(salp){
				if (!(*vnpp)) *vnpp=ValNodeAddPointer(NULL,0,(Pointer)salp);
				else ValNodeAddPointer(vnpp,0,(Pointer)salp);
				salp=salp->next;
			}
		}
		sap2=sap2->next;
	}
}

/*******************************************************************************

  Function : DDV_GetSeqAlign()
  
  Purpose : given dataptr and datatype, look for any SeqAligns
  
  Parameters : dataptr; the data to analyse
  				datatype; OBJ_BIOSEQ,OBJ_SEQENTRY, etc.
				vnp; filled here. data.ptrvalue == SeqALignPtr
  
  Return value : none 

*******************************************************************************/
extern void DDV_GetSeqAlign(Pointer dataptr,Uint2 datatype,
		ValNodePtr PNTR vnp)
{
BioseqPtr     bsp;
BioseqSetPtr  bssp;
SeqEntryPtr   sep;
SeqAnnotPtr   sap;

	if(datatype==OBJ_SEQENTRY){
		SeqEntryExplore((SeqEntryPtr)dataptr,(Pointer) vnp,DDV_SearchAli);
	} 
	else if(datatype==OBJ_BIOSEQ || datatype==OBJ_BIOSEQSET){
		/*create a seqentry; then explore it*/
		sep = SeqEntryNew ();
		if (sep != NULL) {
			if (datatype == OBJ_BIOSEQ) {
				bsp = (BioseqPtr) dataptr;
				sep->choice = 1;
				sep->data.ptrvalue = bsp;
			} else if (datatype == OBJ_BIOSEQSET) {
				bssp = (BioseqSetPtr) dataptr;
				sep->choice = 2;
				sep->data.ptrvalue = bssp;
			}
			SeqEntryExplore(sep,(Pointer) vnp,DDV_SearchAli);
			SeqEntryFree (sep);
		}
	}
	else if(datatype==OBJ_SEQALIGN){
		ValNodeAddPointer(vnp,0,(SeqAlignPtr)dataptr);
	}
	else if(datatype==OBJ_SEQANNOT){
		sap=(SeqAnnotPtr)dataptr;
		while (sap != NULL) {
			if (sap->type == 2) {/*seqalign*/
				ValNodeAddPointer(vnp,0,sap->data);
			}
			sap=sap->next;
		}
	}
}

/*******************************************************************************

  Function : DDV_PrintPopSetSummary()
  
  Purpose : display a summary for a PopSet entry
  
*******************************************************************************/
extern void DDV_PrintPopSetSummary(SeqEntryPtr sep, Int4 gi, FILE *FileOut)
{
ValNodePtr  vnp_sap,vnp_biosrc,vnp_biosrc2,vnp,vnp2;
SeqAlignPtr sap;
SABlockPtr  sabp;
Int4		sapLength, sapNumBioseqs, numBlocks;
Uint2	    i;	
Boolean bPairwise=TRUE,bPrintTax,bPrintCom;
Uint4 disp_format;
PopSourcePtr psp,psp2;

	disp_format=RULER_TOP;
	disp_format|=RULER_TICK;
	disp_format|=DISPE_SHOWBLOCK;
	disp_format|=VIEW_FULLSEQ;
	disp_format|=DISP_FULL_HTML;
	disp_format|=DISP_ORDER_NUM; 
	disp_format|=DISP_STRAND;
	disp_format|=DISP_BSP_COORD;

	/*scan the SeqEntry and get all the SeqALigns*/
	vnp_sap=NULL;
	DDV_GetSeqAlign((Pointer) sep, OBJ_SEQENTRY, &vnp_sap);

	/*scan the SeqEntry and get all the biosources*/
	vnp_biosrc=NULL;
	SeqEntryExplore(sep,(Pointer) &vnp_biosrc,SearchBioSource);

	fprintf(FileOut,"<UL>\n");

	/*summary of biosource; if a name occurs several times, display only once*/
	if (vnp_biosrc){
		vnp_biosrc2=NULL;
		vnp=vnp_biosrc;
		fprintf(FileOut,"<B>Source of the sequences : </B><BR><UL>");
		while(vnp){
			psp=(PopSourcePtr)vnp->data.ptrvalue;
			if (psp->szTaxName) bPrintTax=TRUE;
			else bPrintTax=FALSE;
			if (psp->szCommonName) bPrintCom=TRUE;
			else bPrintCom=FALSE;
			if (bPrintTax || bPrintCom){
				if (vnp_biosrc2){
					vnp2=vnp_biosrc2;
					while(vnp2){
						psp2=(PopSourcePtr)vnp2->data.ptrvalue;
						if (psp->szTaxName && psp2->szTaxName){
							if (!StringCmp(psp->szTaxName,psp2->szTaxName)) bPrintTax=FALSE;
						}
						if (psp->szCommonName && psp2->szCommonName){
							if (!StringCmp(psp->szCommonName,psp2->szCommonName)) bPrintCom=FALSE;
						}
						vnp2=vnp2->next;
					}
					if (bPrintTax || bPrintCom){
						fprintf(FileOut,", \n");
						ValNodeAddPointer(&vnp_biosrc2,0,(Pointer)psp);
					}
				}
				else{
					vnp_biosrc2=ValNodeAddPointer(NULL,0,(Pointer)psp);
				}
				if (bPrintTax && bPrintCom){
					fprintf(FileOut,"%s (%s)",psp->szTaxName,psp->szCommonName);
				}
				else if (bPrintTax && !bPrintCom){
					fprintf(FileOut,"%s",psp->szTaxName);
				}
				else if (!bPrintTax && bPrintCom){
					fprintf(FileOut,"%s",psp->szCommonName);
				}
			}
			vnp=vnp->next;
		}
		fprintf(FileOut,".<BR></UL>\n");
		ValNodeFreeData(vnp_biosrc);
		if (vnp_biosrc2) ValNodeFreeData(vnp_biosrc2);
	}
	
	fprintf(FileOut,"<B>Content :</B><BR><BR><UL>\n");
	/*scan vnp_sap and display the summary for SeqAlign*/
	if(vnp_sap){
		i=0;
		vnp=vnp_sap;
		fprintf(FileOut,"<B>List of sequence alignment(s) : </B><BR><UL>");
		while(vnp){
			/*get a SeqAlign*/
			sap = (SeqAlignPtr) vnp->data.ptrvalue;			
			/*build the indexed SeqAlign*/
			sabp=Blockify(sap);
			if (sabp){
				/*retrieve the size of the alignment*/
				GetBlockInfo(sabp, &sapLength, &sapNumBioseqs,&numBlocks,NULL);
				if (sapNumBioseqs>2) 
					bPairwise=FALSE;
				fprintf(FileOut,"Sequence Alignment no. %i : %i sequences, %i letters ",
						i+1,sapNumBioseqs,sapLength);
				if (gi!=0) {
					fprintf(FileOut," (<a href=\"wwwddv.cgi?gi=%i\">Text view</a>).",gi);
					/*fprintf(FileOut,", graphical summary).");*/
				}
				fprintf(FileOut,"<BR>\n");
				/*if (*szTaxName && *szCommonName){
					fprintf(FileOut,"( %s - %s ).",szTaxName,szCommonName);
				}*/
				/*delete data*/
				CleanupBlocks(sabp);
			}
			vnp=vnp->next;i++;
		}
		if (i>1 && bPairwise){
			DDV_ToolsScript(sep, disp_format,FileOut, FALSE);
		}
		fprintf(FileOut,"<BR></UL>\n");
		ValNodeFree(vnp_sap);
	}
	else{/*summary if no SeqAlign*/
		
		fprintf(FileOut,"This entry contains no sequence alignment.<BR>");
		if (gi!=0) {
			fprintf(FileOut,"Yan can either :<BR><BR>");
			fprintf(FileOut,"<a href=\"wwwddv.cgi?gi=%i\"",gi);
			fprintf(FileOut," onMouseOut=\"window.status='';return true\" \
onMouseOver=\"window.status='Click to see the list of sequences'; return true\" >");
			fprintf(FileOut,"- view the list of sequences</a>, or <BR>\n");
			DDV_ToolsScript(sep, disp_format,FileOut, FALSE);
		}
	}
	fprintf(FileOut,"</UL></UL>\n");
}

/**************************************************************
***************************************************************
COLOMBE FOR KARL 
***************************************************************
***************************************************************/

extern void PrintSeqAlignCallback (SeqEntryPtr sep, Pointer mydata,
                                   Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqAnnotPtr        sap;
  SeqAlignPtr        salp,salp2;
  Uint4              option = 0;
  /*FILE               *fp;*/
  ByteStorePtr PNTR byteSpp;
  DDVOptionsBlockPtr dobp = (DDVOptionsBlockPtr) mydata;
  
  if (sep != NULL && sep->data.ptrvalue)
  {
     option = RULER_TOP;
     option|=DISPE_SHOWBLOCK;
     option|=VIEW_VARIA;
     option|=DISP_FULL_TXT;
     option|=DISP_STRAND;
     option|=DISP_BSP_COORD;

     if (IS_Bioseq(sep))
     {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           for (sap=bsp->annot;sap!=NULL;sap=sap->next)
           {
              if (sap->type == 2) {
                 salp = (SeqAlignPtr)sap->data;
                 salp2 = PairSeqAlign2MultiSeqAlign (salp);
                 if (salp2==NULL) {
                    salp2 = (SeqAlignPtr) AsnIoMemCopy (salp, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);   
                 }
                 DDV_DisplayDefaultAlign (salp2, 0, -1,-1, option, dobp, NULL);
                 SeqAlignSetFree(salp2);
                 BSWrite(*(dobp ->byteSpp),"\n//\n",4);
              }
           }
        }
     }
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           for (sap=bssp->annot;sap!=NULL;sap=sap->next)
           {
              if (sap->type == 2) {
                 salp = (SeqAlignPtr)sap->data;
                 salp2 = PairSeqAlign2MultiSeqAlign (salp);
                 if (salp2==NULL) {
                    salp2 = (SeqAlignPtr) AsnIoMemCopy (salp, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);   
                 }
                 DDV_DisplayDefaultAlign (salp2, 0, -1,-1, option, dobp, NULL);
                 SeqAlignSetFree(salp2);
                 BSWrite(*(dobp ->byteSpp),"\n//\n",4);
              }
           }
        }
     }
  }
}

extern ByteStorePtr SeqAlignToBS (Uint2 entityID)
{{
  SeqEntryPtr sep;
  ByteStorePtr bsp = NULL;
  FILE *fp;
  DDVOptionsBlockPtr dobp=MemNew(sizeof (DDVOptionsBlock));

  dobp->LineSize=50; 
  dobp->byteSpp = & bsp;
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep) {
     SeqEntryExplore (sep, (Pointer)dobp, PrintSeqAlignCallback);
  }
  dobp->byteSpp=NULL;
  MemFree(dobp);
  return bsp;
}}

/****/

