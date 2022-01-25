/*  $Id: ddvcreate.c,v 1.23 1999/12/29 22:55:03 lewisg Exp $
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
* File Name:  ddvcreate.c
*
* Author:  Sarah Wheelan
*
* Version Creation Date:   08/99
*
* $Revision: 1.23 $
*
* File Description: 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ddvcreate.c,v $
* Revision 1.23  1999/12/29 22:55:03  lewisg
* get rid of seqalign id
*
* Revision 1.22  1999/12/20 20:20:41  lewisg
* allow cn3d to do color and ddv to do case when both are running
*
* Revision 1.21  1999/12/20 14:45:37  durand
* add new functions for the new BLAST outputs
*
* Revision 1.20  1999/12/08 22:42:17  durand
* add the code to produce colored BLAST outputs
*
* Revision 1.19  1999/12/03 23:17:22  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 1.18  1999/12/03 18:24:08  durand
* gapped regions are displayed using lower cases; BLAST outputs
*
* Revision 1.17  1999/12/03 13:23:08  durand
* replace AlnMgrGetNumSeqs by AlnMgrGetNumRows
*
* Revision 1.16  1999/12/02 14:46:42  durand
* fix a bug in UABuildDescriptor; problem with unaligned region of size 0
*
* Revision 1.15  1999/11/24 21:26:29  vakatov
* Fixed for the C++ and/or MSVC DLL compilation
*
* Revision 1.14  1999/11/17 22:44:02  durand
* speed up the selection manager for large SeqAlign
*
* Revision 1.13  1999/11/09 22:14:01  shavirin
* Added parameter follower to the Blast score printing function
*
* Revision 1.12  1999/11/09 17:08:59  durand
* transfer some functions from ddvgraph to ddvcreate, so that ddvcreate remains Vibrant free and can be compiled with BLAST
*
* Revision 1.11  1999/11/03 21:29:48  durand
* add CTRL and SHFT keys for mouse selection. redesign the loader functions of DDV to properly register the message callbacks
*
* Revision 1.10  1999/10/29 14:15:40  durand
* add simple mouse selection functions
*
* Revision 1.9  1999/10/22 20:12:47  durand
* add Export command (text, HTML and Phylip formats)
*
* Revision 1.8  1999/10/22 14:57:25  durand
* use AlnMgrIndexSingleChildSeqAlign to build BLAST output
*
* Revision 1.7  1999/10/22 13:17:19  durand
* remove unreferenced variables to avoid warnings with Visual C
*
* Revision 1.6  1999/10/20 13:17:19  durand
* add display for disc. SeqAlign tails
*
* Revision 1.5  1999/10/15 21:57:35  durand
* add a UI for display options
*
* Revision 1.4  1999/10/12 15:01:29  lewisg
* resolve confict with internal/ddv
*
* Revision 1.3  1999/10/07 13:36:08  wheelan
* bug fixes -- AlnLengths were 1 too long
*
* Revision 1.2  1999/10/05 14:04:21  wheelan
* made DDV_CreateDisplayFromIndex_2, which creates the display for discontinuous alignments; minor bug fixes
*
* Revision 1.1  1999/09/30 14:10:25  durand
* add ddv to toolkit
*
* Revision 1.11  1999/09/30 13:38:09  durand
* DDV_CreateDisplayFromIndex takes ParaG_Size as an argument
*
*
*
* ==========================================================================
*/
#include <ddvcreate.h>
#include <pgppop.h>
#include <udvseq.h>
#include <tofasta.h>
#include <sqnutils.h>

typedef struct ddvdataforcolorfunc{
	SeqIdPtr sip;
	Int4     from;
	Int4     to;
	Uint1    strand;
	Boolean  IsUnAligned;
	Uint1    style;
	Uint1    rgb[3];
} DDVDataForColorFunc, PNTR DDVDataForColorFuncPtr;

#define MAX_NCBIstdaa 26
#define MAX_NCBI4na 17
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
#define DDV_COLOR_MAX 33 
 static Uint1 DDV_PaletteRGB[DDV_COLOR_MAX][3] =
 {
     /* Red  Grn Blue     Color     ColorID             Old Hex  New Hex */
   255, 255, 255, /* default     0                 0xffffff - */ 
   255,   0, 153, /* hotpink     1                 0xff1493 0xff0099 */
   255,   0, 255, /* magenta     2                 0xff00ff - */
   153,  51, 255, /* purple      3                 0x9b30ff 0x9933ff */
     0,   0, 255, /* blue        4                 0x0000ff - */
    00, 153, 255, /* sky         5                 0x1e90ff 0x0099ff */
     0, 255, 255, /* cyan        6                 0x00ffff - */
     0, 255, 153, /* sea         7                 0x00ff8f 0x00ff99 */
     0, 192,   0, /* green       8                 0x00ff00 - */
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
Uint1 DDV_STD_AAColor[MAX_NCBIstdaa] =
{DDVCOL_BROWN, /*-*/
DDVCOL_GRAY,   /*A*/
DDVCOL_ORANGE, /*B*/
DDVCOL_PINK, /*C*/
DDVCOL_RED,    /*D*/
DDVCOL_RED,    /*E*/
DDVCOL_PURPLE, /*F*/
DDVCOL_GRAY,   /*G*/
DDVCOL_BLUE,   /*H*/
DDVCOL_GRAY,   /*I*/
DDVCOL_BLUE,   /*K*/
DDVCOL_GRAY,   /*L*/
DDVCOL_PINK, /*M*/
DDVCOL_CYAN,   /*N*/
DDVCOL_MAGENTA,/*P*/
DDVCOL_CYAN,   /*Q*/
DDVCOL_BLUE,   /*R*/
DDVCOL_GREEN,  /*S*/
DDVCOL_GREEN,  /*T*/
DDVCOL_GRAY,   /*V*/
DDVCOL_PURPLE, /*W*/
DDVCOL_GRAY,   /*X*/
DDVCOL_PURPLE, /*Y*/
DDVCOL_ORANGE, /*Z*/
DDVCOL_GOLD,   /*U*/
DDVCOL_HOTPINK /***/};

Uint1 DDV_STD_NAColor[MAX_NCBI4na] = {
DDVCOL_PINK,   /* gap */
DDVCOL_GREEN,  /* A */
DDVCOL_BLUE,   /* C */
DDVCOL_CYAN,   /* A or C */
DDVCOL_BLACK,  /* G */
DDVCOL_CYAN,   /* G or A */
DDVCOL_CYAN,   /* G or C */
DDVCOL_PURPLE, /* G C or A */
DDVCOL_RED,    /* T */
DDVCOL_CYAN,   /* A or T */
DDVCOL_CYAN,   /* T or C */
DDVCOL_PURPLE, /* A C ot T */
DDVCOL_CYAN,   /* G ot T */
DDVCOL_PURPLE, /* G or A or T */
DDVCOL_PURPLE, /* G T or C */
DDVCOL_HOTPINK,/* G T A or C */
DDVCOL_GRAY /*?, but MAX_NCBI4na==17, see mmdbapi2.h*/
};


extern void display_ParaG_content(MsaParaGPopListPtr mpplp,Int2 LineSize)
{
BioseqPtr bsp_debug=NULL;
Char                    szBuf[15];
ParaGPtr                pgp;
MsaTxtDispPtr mtdp;
ValNodePtr    vnp,vnp2;
Uint4 j,k,size;
FILE *fp;
Boolean bClose=TRUE;
Int4 i;
	if (LineSize<=50)/*just for the integrity of the function...*/
		LineSize=ParaG_Size;
		
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
                                fprintf(fp,"  ->ParaG[%d] , range (%7d..%7d):\n",k++,size,size+LineSize-1);
                                size+=LineSize;
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

static Boolean DDV_CreateDisplayFromIndex_1(SeqAlignPtr sap, MsaParaGPopListPtr mpplp, Int2 LineSize)
{
   AlnMsgPtr      amp;
   Int4           i;
   Boolean        more;
   MsaTxtDispPtr  mtdp;
   Int4           n;
   ParaGPtr       pgp;
   SeqIdPtr       sip;
   ValNodePtr     vnp;
   ValNodePtr     vnp_head;
   ValNodePtr     PNTR vnp_list;  
   ValNodePtr     vnp_para;
   ValNodePtr     vnp_para_head;

	if (LineSize<=50)/*just for the integrity of the function...*/
		LineSize=ParaG_Size;
		
   amp = AlnMsgNew();
   vnp_list = (ValNodePtr PNTR)MemNew((mpplp->nBsp)*sizeof(ValNodePtr));
   for (n = 1; n<=(mpplp->nBsp); n++)
   {
      vnp_para = vnp_para_head = NULL;
      sip = AlnMgrGetNthSeqIdPtr(mpplp->sap, n);
      for (i = 0; i<mpplp->LengthAli; i+=LineSize)
      {
         amp = AlnMsgReNew(amp);
         pgp = (ParaGPtr)MemNew(sizeof(ParaG));
         pgp->StartLine = n - 1;
         pgp->nLines = 1;
         pgp->StartLetter = amp->from_m = i;
         if (i + LineSize >= mpplp->LengthAli)
            pgp->StopLetter = amp->to_m = mpplp->LengthAli - 1;
         else
            pgp->StopLetter = amp->to_m = i + LineSize -1;
         pgp->sip = sip;
         pgp->ScaleStyle = SCALE_POS_NONE;
         amp->which_bsq = NULL;
         amp->row_num = n;
         vnp = vnp_head = NULL;
         amp->place = 0;
         while (more = AlnMgrGetNextAlnBit(mpplp->sap, amp))
         {
            vnp = ValNodeAdd(&vnp);
            if (!vnp_head)
               vnp_head = vnp;
            mtdp = (MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
            mtdp->from = amp->from_b;
            mtdp->to = amp->to_b;
            mtdp->IsGap = (Boolean)(amp->gap);
            if (mtdp->IsGap)
               mtdp->TextStyle = MSA_TXT_STYLE_GAP;
            else
               mtdp->TextStyle = MSA_TXT_STYLE_SEQ;
            mtdp->strand = amp->strand;
			mtdp->IsUnAligned=FALSE;
            vnp->data.ptrvalue = mtdp;
         }
         pgp->ptxtList = vnp_head;
         vnp_para = ValNodeAdd(&vnp_para);
         vnp_para->data.ptrvalue = pgp;
         if (!vnp_para_head)
            vnp_para_head = vnp_para;
      }
      vnp_list[n-1]=vnp_para_head;
   }
   mpplp->TableHead = vnp_list;
   display_ParaG_content(mpplp,LineSize);
   DDV_LocateParaG(mpplp->TableHead,mpplp->nBsp);
   return TRUE;
}


static Boolean DDV_CreateDisplayFromIndex_2(SeqAlignPtr sap, MsaParaGPopListPtr mpplp, Int2 LineSize)
{
   AlnMsgPtr        amp;
   Int4             i;
   Boolean          more;
   MsaTxtDispPtr    mtdp;
   Int4             n;
   ParaGPtr         pgp;


   ValNodePtr       vnp;
   ValNodePtr       vnp_head;
   ValNodePtr       PNTR vnp_list;
   ValNodePtr       vnp_para;
   ValNodePtr       vnp_para_head;
   Int4             PopulateTo;
   Int4             ReportSpacer;
   Int4             PlaceSpacer;
   Int4             RealAlnLength;
   SeqIdPtr         sip;
   Boolean          last;


   amp = AlnMsgNew();
   vnp_list = (ValNodePtr PNTR)MemNew((mpplp->nBsp)*sizeof(ValNodePtr));
   RealAlnLength = AlnMgrGetAlnLength(mpplp->sap, FALSE);
   for (n = 1; n<=(mpplp->nBsp); n++)
   {
      vnp_para = vnp_para_head = NULL;
      sip = AlnMgrGetNthSeqIdPtr(mpplp->sap, n);
          ReportSpacer=0;
          amp->to_m=-1;
      amp = AlnMsgReNew(amp);
      last = FALSE;
      for (i = 0; i<mpplp->LengthAli; i+=LineSize)
      {
         amp->place = 0;
         pgp = (ParaGPtr)MemNew(sizeof(ParaG));
         pgp->StartLine = n - 1;
         pgp->nLines = 1;
         pgp->StartLetter = i;
         if (i + LineSize >= mpplp->LengthAli)
         {
            last = TRUE;
            pgp->StopLetter = mpplp->LengthAli;
         }
         else
            pgp->StopLetter = i + LineSize -1;
         pgp->sip = sip;
         pgp->ScaleStyle = SCALE_POS_NONE;
         amp->which_master = 0;
         if (amp->to_m != 0)
                 amp->from_m = amp->to_m + 1 ;
         else
                 amp->from_m = 0;
         amp->to_m = amp->from_m + LineSize -1;
         amp->row_num = n;
         if (amp->to_m > RealAlnLength)
            amp->to_m = RealAlnLength - 1;
         more = TRUE;
         vnp = vnp_head = NULL;
                 PopulateTo=0;
                 PlaceSpacer=0;
         while (amp->place == 0)
         {
            vnp = ValNodeAdd(&vnp);
            if (!vnp_head)
               vnp_head = vnp;
            if (amp->send_space)
            {
               if (!last)
                  amp->to_m=amp->to_m - PlaceSpacer;
            }
                        if (ReportSpacer){
                           amp->to_m -=ReportSpacer;
               mtdp = (MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
               mtdp->from = pgp->StartLetter;
               mtdp->to = mtdp->from+ReportSpacer-1;
                           PopulateTo+=(ReportSpacer+1);
                           ReportSpacer=0;
               mtdp->IsGap = TRUE;
               mtdp->TextStyle = MSA_TXT_STYLE_1;
               vnp->data.ptrvalue = mtdp;
               vnp = ValNodeAdd(&vnp);
                        }
            amp->send_space = FALSE;
            more = AlnMgrGetNextAlnBit(mpplp->sap, amp);
            mtdp = (MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
            mtdp->from = amp->from_b;
            mtdp->to = amp->to_b;
            PopulateTo+=(amp->to_b-amp->from_b+1);
                        mtdp->IsGap = (Boolean)(amp->gap);
            if (mtdp->IsGap)
               mtdp->TextStyle = MSA_TXT_STYLE_GAP;
            else
               mtdp->TextStyle = MSA_TXT_STYLE_SEQ;
            mtdp->strand = amp->strand;
            vnp->data.ptrvalue = mtdp;
            if (amp->send_space)
            {
               vnp = ValNodeAdd(&vnp);
               mtdp = (MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
               mtdp->from = pgp->StartLetter+PopulateTo;
               mtdp->to = mtdp->from+SPACER_TXT_BLANK-1;
                           if (mtdp->to>pgp->StopLetter){
                               ReportSpacer=mtdp->to-pgp->StopLetter;
                               mtdp->to=pgp->StopLetter;
               }
                           else{
                               ReportSpacer=0;
                           }
                           PlaceSpacer=SPACER_TXT_BLANK-ReportSpacer;
                           PopulateTo+=(mtdp->to-mtdp->from+1);
               mtdp->IsGap = TRUE;
               mtdp->TextStyle = MSA_TXT_STYLE_1;
				mtdp->IsUnAligned=FALSE;
               vnp->data.ptrvalue = mtdp;
            }
         }
         pgp->ptxtList = vnp_head;
         vnp_para = ValNodeAdd(&vnp_para);
         vnp_para->data.ptrvalue = pgp;
         if (!vnp_para_head)
            vnp_para_head = vnp_para;
      }
      vnp_list[n-1]=vnp_para_head;
   }
   mpplp->TableHead = vnp_list;
   display_ParaG_content(mpplp,LineSize);
   DDV_LocateParaG(mpplp->TableHead,mpplp->nBsp);
   return TRUE;
}

NLM_EXTERN Boolean DDV_CreateDisplayFromIndex(SeqAlignPtr sap, MsaParaGPopListPtr mpplp, 
		Int2 LineSize, DDV_Disp_OptPtr ddop)
{
   AMAlignIndexPtr  amaip;
   AlnMsgPtr        amp;
   Int4             i;
   Boolean          more;
   MsaTxtDispPtr    mtdp;
   Int4             n;
   ParaGPtr         pgp;


   ValNodePtr       vnp;
   ValNodePtr       vnp_head;
   ValNodePtr       PNTR vnp_list;
   ValNodePtr       vnp_para;
   ValNodePtr       vnp_para_head;
   Int4             PopulateTo;


   SeqIdPtr         sip;
   Boolean          last;

	if (LineSize<=50)/*just for the integrity of the function...*/
		LineSize=ParaG_Size;
		
   if (!sap)
      return FALSE;
   if (!AlnMgrIndexSeqAlign(sap))
      return FALSE;
   if (sap->saip->indextype == INDEX_PARENT)
   {
      amaip = (AMAlignIndexPtr)sap->saip;
      if (sap->type == SAT_PARTIAL || (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_SEGMENTED_MASTERSLAVE))
      {
         /*mpplp->sap = sap;
         if (sap->type == SAT_PARTIAL)
            mpplp->LengthAli = AlnMgrGetAlnLength(mpplp->sap, FALSE)+SPACER_TXT_BLANK*(amaip->alnsaps-1);
         else
            mpplp->LengthAli = AlnMgrGetAlnLength(mpplp->sap, FALSE) + SPACER_TXT_BLANK*(amaip->numseg-1);
         mpplp->nBsp = AlnMgrGetNumRows(mpplp->sap);
         mpplp->DisplayType = DDV_DISP_HORZ;*/
         if (!DDV_CreateDisplay_DiscAlign(sap,mpplp,LineSize,
				ddop))
         /*if (!DDV_CreateDisplayFromIndex_2(mpplp->sap, mpplp, LineSize))*/
            return FALSE;
         else
            return TRUE;
      } else if (sap->type == SAT_MASTERSLAVE && amaip->mstype == AM_MASTERSLAVE)
      {
         mpplp->sap = sap;
         mpplp->LengthAli = AlnMgrGetAlnLength(mpplp->sap, FALSE);
         mpplp->nBsp = AlnMgrGetNumRows(mpplp->sap);
         mpplp->DisplayType = DDV_DISP_HORZ;
         if (!DDV_CreateDisplayFromIndex_1(sap, mpplp,LineSize))
            return FALSE;
         else
            return TRUE;
      } else
         mpplp->sap = (SeqAlignPtr)sap->segs;
   } else 
      mpplp->sap = sap;
   mpplp->LengthAli = AlnMgrGetAlnLength(mpplp->sap, FALSE);
   mpplp->nBsp = AlnMgrGetNumRows(mpplp->sap);
   mpplp->DisplayType = DDV_DISP_HORZ;
   amp = AlnMsgNew();
   vnp_list = (ValNodePtr PNTR)MemNew((mpplp->nBsp)*sizeof(ValNodePtr));
   
   for (n = 1; n<=(mpplp->nBsp); n++)
   {
      vnp_para = vnp_para_head = NULL;
      sip = AlnMgrGetNthSeqIdPtr(mpplp->sap, n);
	  amp->to_m=-1;
      amp = AlnMsgReNew(amp);
      last = FALSE;
      for (i = 0; i<mpplp->LengthAli; i+=LineSize)
      {
         amp->place = 0;
         pgp = (ParaGPtr)MemNew(sizeof(ParaG));
         pgp->StartLine = n - 1;
         pgp->nLines = 1;
         pgp->StartLetter = i;
         if (i + LineSize >= mpplp->LengthAli)
         {
            last = TRUE;
            pgp->StopLetter = mpplp->LengthAli;
         }
         else
            pgp->StopLetter = i + LineSize -1;
         pgp->sip = sip;
         pgp->ScaleStyle = SCALE_POS_NONE;
         amp->which_master = 0;
         if (amp->to_m != 0)
		 amp->from_m = amp->to_m + 1 ;
         else
                 amp->from_m = 0;
	 amp->to_m = amp->from_m + LineSize -1;
         amp->row_num = n;
         if (amp->to_m > mpplp->LengthAli)
            amp->to_m = mpplp->LengthAli - 1;
         more = TRUE;
         vnp = vnp_head = NULL;
		 PopulateTo=0;
         while (more = AlnMgrGetNextAlnBit(mpplp->sap, amp))
         {
            vnp = ValNodeAdd(&vnp);
            if (!vnp_head)
               vnp_head = vnp;
            mtdp = (MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
            mtdp->from = amp->from_b;
            mtdp->to = amp->to_b;
			mtdp->IsGap = (Boolean)(amp->gap);
            if (mtdp->IsGap)
               mtdp->TextStyle = MSA_TXT_STYLE_GAP;
            else
               mtdp->TextStyle = MSA_TXT_STYLE_SEQ;
            mtdp->strand = amp->strand;
			mtdp->IsUnAligned=FALSE;
            vnp->data.ptrvalue = mtdp;
         }
         pgp->ptxtList = vnp_head;
         vnp_para = ValNodeAdd(&vnp_para);
         vnp_para->data.ptrvalue = pgp;
         if (!vnp_para_head)
            vnp_para_head = vnp_para;
      }
      vnp_list[n-1]=vnp_para_head;
   }
   mpplp->TableHead = vnp_list;
   display_ParaG_content(mpplp,LineSize);
   DDV_LocateParaG(mpplp->TableHead,mpplp->nBsp);
   return TRUE;
}


/*******************************************************************************

  Function : DDV_CreateDisplay_DiscAlign()
  
  Purpose : create a display for a Disc. SeqAlign

  Note :	Style : is used to determine how to display an unaligned region
			Justification : in case where we want to display the sequence
				  in an UA region, this parameter tells the function how to
				  justify the sequence.

  Return value : FALSE if failure

*******************************************************************************/
NLM_EXTERN Boolean DDV_CreateDisplay_DiscAlign(SeqAlignPtr sap, 
		MsaParaGPopListPtr mpplp, Int2 LineSize,DDV_Disp_OptPtr ddop)
{
Int4            TotLength,cumulpop,StartLetter,n,nPgp,from_bsp,to_bsp,MaxLength,
                BspLength,BspStart,BspStop,nBsp,AbsPos;
Boolean         more,bFirstMtdp,bFirstPgp,IsGap;
ValNodePtr      vnp,vnp2,vnp3,vnp_head,vnp_para,vnp_mtdp;
ValNodePtr PNTR vnp_list;
ParaGPtr        pgp;
SeqIdPtr        sip;
DescriDispPtr   ddp;
MsaTxtDispPtr   mtdp_pgp;
AlnMsgPtr       amp;
UAMsgPtr        uamp;
Uint1           strand,strand_tmp;

	/*descriptor*/
	nBsp = AlnMgrGetNumRows(sap);
	vnp_head=UABuildDescriptor(sap,nBsp,LineSize,ddop,&TotLength,
		ddop->ShowLeftTail,ddop->ShowRightTail);
	if (vnp_head==NULL)
		return(FALSE);

	/*now, use the previous descriptor to build the ParaG List*/
	mpplp->sap=sap;
	mpplp->nBsp = nBsp;
	mpplp->LengthAli = TotLength;
	mpplp->DisplayType = DDV_DISP_HORZ;
	amp = AlnMsgNew();
	vnp_list = (ValNodePtr PNTR)MemNew((mpplp->nBsp)*sizeof(ValNodePtr));
	if (!vnp_list || !amp){
		/*delete the descriptor*/
		ValNodeFreeData(vnp_head);
		return(FALSE);
	}

	cumulpop=0;
	
	/*loop on each sequence*/
	for (n = 1; n<=(mpplp->nBsp); n++){
		sip = AlnMgrGetNthSeqIdPtr(mpplp->sap, n);
		/*for each sequence, use the descriptors to build the ParaG list*/
		bFirstPgp=TRUE;
		bFirstMtdp=TRUE;
		StartLetter=0;AbsPos=0;/*disp coord*/
		nPgp=0;
		vnp=vnp_head;
		while(vnp){
			ddp=(DescriDispPtr)vnp->data.ptrvalue;
			if (ddp->TextStyle==MSA_TXT_STYLE_REG_ALIGN){/*aligned*/
				/*initialize the amp struct for the current aligned region*/
				amp = AlnMsgReNew(amp);
				amp->place = 0;
				amp->from_m = ddp->from;
				amp->to_m=ddp->to;
				amp->which_master=0;
				amp->row_num=n;
				amp->send_space = FALSE;
				more=TRUE;
				AlnMgrGetNextAlnBit(mpplp->sap, amp);
				while (more){
					/*list of TxT descriptors*/
					mtdp_pgp = (MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
					mtdp_pgp->IsGap = (Boolean)(amp->gap);
					if (mtdp_pgp->IsGap){/*disp coord*/
						mtdp_pgp->from=AbsPos;
						mtdp_pgp->to=mtdp_pgp->from+(amp->to_b-amp->from_b);
						mtdp_pgp->TextStyle = MSA_TXT_STYLE_GAP;
					}
					else{/*BSP coord*/
						mtdp_pgp->from = amp->from_b;
						mtdp_pgp->to = amp->to_b;
						mtdp_pgp->TextStyle = MSA_TXT_STYLE_SEQ;
					}
					mtdp_pgp->strand = amp->strand;
					mtdp_pgp->IsUnAligned=FALSE;
					strand=amp->strand;
					/*add the new node*/
					if (bFirstMtdp){
						vnp2=ValNodeAddPointer(NULL,0,(Pointer)mtdp_pgp);
						vnp_mtdp=vnp2;
						bFirstMtdp=FALSE;
					}
					else{
						vnp2=ValNodeAddPointer(&vnp2,0,(Pointer)mtdp_pgp);
					}
					more=AlnMgrGetNextAlnBit(mpplp->sap, amp);
					AbsPos+=(mtdp_pgp->to-mtdp_pgp->from+1);
				}
			}
			else{/*unaligned*/
				switch(ddop->DispDiscStyle){/*user's display choice*/
					case MSA_TXT_STYLE_1:/*just put a little spacer*/
						mtdp_pgp = (MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
						mtdp_pgp->from = StartLetter;/*put DISP coord*/
						mtdp_pgp->to = mtdp_pgp->from+(ddp->to-ddp->from);
						mtdp_pgp->IsGap = TRUE;
						mtdp_pgp->TextStyle = MSA_TXT_STYLE_1;
						mtdp_pgp->strand = strand;
						mtdp_pgp->IsUnAligned=TRUE;
						/*add the new node*/
						if (bFirstMtdp){
							vnp2=ValNodeAddPointer(NULL,0,(Pointer)mtdp_pgp);
							vnp_mtdp=vnp2;
							bFirstMtdp=FALSE;
						}
						else{
							vnp2=ValNodeAddPointer(&vnp2,0,(Pointer)mtdp_pgp);
						}
						break;
					case MSA_TXT_STYLE_2:/*display the unaligned sequence*/
						/*init the analysis process of an UA region*/
						switch(ddp->UAnum){
							case (Int4)-1:/*left tail*/
								uamp=UAMgrIntUAMsg(ddp->from, ddp->to, 
									DISP_JUST_RIGHT, strand);
								break;
							case (Int4)-2:/*right tail*/
								uamp=UAMgrIntUAMsg(ddp->from, ddp->to, 
									DISP_JUST_LEFT, strand);
								break;
							default:/*within SAP*/
								uamp=UAMgrIntUAMsg(ddp->from, ddp->to, 
									ddop->DiscJustification, strand);
								break;
						}
						
						MaxLength=ddp->UAMaxlength;
						
						if (ddp->UAnum>0){/*within SAP*/
							AlnMgrGetNthUnalignedForNthRow(sap, ddp->UAnum, n, &BspStart, 
								&BspStop);
						}
						else{/*Tails of the SAP*/
							AlnMgrGetNthRowTail(sap,n,
								(Uint1)(ddp->UAnum==(Int4)-1 ? LEFT_TAIL : RIGHT_TAIL),
								&BspStart,&BspStop,&strand_tmp);
						
						}
												
						if (BspStart!=(Int4)-1 && BspStop!=(Int4)-1){
							BspLength=BspStop-BspStart+1;
						}
						/*analyse the content of the UA region*/
						while(UAMgrGetNextUAbit(uamp,MaxLength,BspLength,BspStart,BspStop)){
							/*get data*/
							UAMgrUAMsgGetInfo(uamp,&from_bsp,&to_bsp,&IsGap);
							/*list of TxT descriptors*/
							mtdp_pgp = (MsaTxtDispPtr)MemNew(sizeof(MsaTxtDisp));
							mtdp_pgp->IsGap = IsGap;
							if (mtdp_pgp->IsGap){/*disp coord*/
								mtdp_pgp->TextStyle = MSA_TXT_STYLE_GAP;
								mtdp_pgp->from=AbsPos;
								mtdp_pgp->to=mtdp_pgp->from+(to_bsp-from_bsp);
							}
							else{/*BSP coord*/
								mtdp_pgp->from = from_bsp;
								mtdp_pgp->to = to_bsp;
								mtdp_pgp->TextStyle = MSA_TXT_STYLE_SEQ;
							}
							mtdp_pgp->strand = strand;
							mtdp_pgp->IsUnAligned=TRUE;
							/*add the new node*/
							if (bFirstMtdp){
								vnp2=ValNodeAddPointer(NULL,0,(Pointer)mtdp_pgp);
								vnp_mtdp=vnp2;
								bFirstMtdp=FALSE;
							}
							else{
								vnp2=ValNodeAddPointer(&vnp2,0,(Pointer)mtdp_pgp);
							}
							AbsPos+=(mtdp_pgp->to-mtdp_pgp->from+1);
						}
						uamp=UAMgrDelUAMsg(&uamp);
						break;
				}
			}
			cumulpop+=(ddp->to-ddp->from+1);
			StartLetter+=(ddp->to-ddp->from+1);
			if (cumulpop==LineSize || vnp->next==NULL){
				/*new ParaG*/
				pgp=(ParaGPtr)MemNew(sizeof(ParaG));
				pgp->StartLine=n-1;
				pgp->nLines=1;
				pgp->StartLetter=nPgp*LineSize;
				pgp->StopLetter=_min_(pgp->StartLetter+LineSize-1,mpplp->LengthAli-1);
				pgp->sip=sip;
				pgp->ScaleStyle=SCALE_POS_NONE;
				pgp->ptxtList=vnp_mtdp;
				if (bFirstPgp){
					vnp3=ValNodeAddPointer(NULL,0,(Pointer)pgp);
					vnp_para=vnp3;
					bFirstPgp=FALSE;
				}
				else{
					vnp3=ValNodeAddPointer(&vnp3,0,(Pointer)pgp);
				}
				bFirstMtdp=TRUE;
				cumulpop=0;
				nPgp++;
			}
			vnp=vnp->next;
		}/*while() on the descriptor list*/
		vnp_list[n-1]=vnp_para;
	}/*for() on the bsp */
	
	/*delete the descriptor*/
	ValNodeFreeData(vnp_head);

	mpplp->TableHead = vnp_list;
	/*just for debugging purpose*/
	display_ParaG_content(mpplp,LineSize);
	DDV_LocateParaG(mpplp->TableHead,mpplp->nBsp);
	MemFree(amp);

	return(TRUE);
}

/*******************************************************************************

  Function : UAMgrIntUAMsg()
  
  Purpose : constructor for a struct UAMsg
  
  Note : from and to are in UnAligned region coord. system. Such a coord system
         always starts at 0. Length = max length of the UA region
		 
  Return value : -

*******************************************************************************/
NLM_EXTERN UAMsgPtr UAMgrIntUAMsg(Int4 from, Int4 to, Uint1 DispType, Uint1 strand)
{
UAMsgPtr uamp=NULL;

	uamp=(UAMsgPtr)MemNew(sizeof(UAMsg));
	if (uamp){
		uamp->from_ua=from;
		uamp->to_ua=to;
		uamp->DispType=DispType;
		uamp->strand=strand;
		uamp->CurrentPos=(Int4)-1;
	}
	return(uamp);
}

/*******************************************************************************

  Function : UAMgrDelUAMsg()
  
  Purpose : destructor for UAMsg structure
    
  Return value : -

*******************************************************************************/
NLM_EXTERN UAMsgPtr UAMgrDelUAMsg(UAMsgPtr PNTR uamp)
{
	if (*uamp)
		MemFree(*uamp);
	return(NULL);
}

/*******************************************************************************

  Function : UAMgrUAMsgGetInfo()
  
  Purpose : query function for a UAMsg structure
  
  Note : from_bsp and to_bsp are in BSP coord system.
    
  Return value : -

*******************************************************************************/
NLM_EXTERN void UAMgrUAMsgGetInfo(UAMsgPtr uamp,Int4Ptr from_bsp, Int4Ptr to_bsp, 
		BoolPtr IsGap)
{
	*from_bsp=uamp->from_r;
	*to_bsp=uamp->to_r;
	*IsGap=uamp->IsGap;
}

/*******************************************************************************

  Function : UAMgrGetNextUAbit()
  
  Purpose : analytical function to retrieve information about an unaligned
            region. uamp must init. before entering this function.
    
  Return value : -

*******************************************************************************/
NLM_EXTERN Boolean UAMgrGetNextUAbit(UAMsgPtr uamp, Int4 MaxLength, Int4 BspLength,
		Int4 BspStart,Int4 BspStop)
{
Int4 start,stop,diff,diff2,halfBspLength,max,tmp;

	/*means we are outside the scanned region*/
	if (uamp->CurrentPos>uamp->to_ua) 
		return(FALSE);
	
	/*first time we enter that function*/
	if (uamp->CurrentPos==(Int4)-1){
		uamp->CurrentPos=uamp->from_ua;
	}

	diff2=0;
	
	if (BspStart!=(Int4)-1 && BspStop!=(Int4)-1){
		halfBspLength=BspLength/2;
		
		/*compute start and stop (UA coord system) of the sequence. It depends
		on the display option*/
		switch(uamp->DispType){
			case DISP_JUST_RIGHT:
				start=MaxLength-BspLength;
				stop=MaxLength-1;
				break;
			case DISP_JUST_CENTER:
				start=MaxLength/2-BspLength/2;
				stop=start+BspLength-1;
				break;
			case DISP_JUST_SPLIT:
				/*particular case; we have to deal with two BSP segments (one is
				left-justified, the other is right-justified). The Bioseq is splitted
				in two equal pieces : halfBspLength*/
				/*they are two exceptions here : don't split if BspLength==
				MaxLength and if BspLength==1; use DISP_JUST_LEFT instead*/
				if (BspLength!=MaxLength && BspLength!=1){
					/*here we compute the position of the GAP */
					start=halfBspLength;
					stop=MaxLength-halfBspLength-2;
					/*look if we are in one of the two BSP segments or in the gap*/
					if(uamp->from_ua<start){
						/*left aligned segment*/
						start=0;
						stop=halfBspLength-1;
						diff2=0;
					}
					else if (uamp->to_ua>stop){
						/*right aligned segment*/
						start=MaxLength-halfBspLength;
						stop=MaxLength-1;
						diff2=halfBspLength;
					}
					else{
						/*he, he, he : we are in the gap !*/
						start=(Int4)-1;
						stop=(Int4)-1;
					}
				}
				else{
					start=0;
					stop=BspLength-1;
				}
				break;
			case DISP_JUST_LEFT:
			default:
				start=0;
				stop=BspLength-1;
				break;
		}
	}
	else{
		/*he, he, he : we are in the gap !*/
		start=(Int4)-1;
		stop=(Int4)-1;
	}
	
	/*are we in the sequence*/
	if (start<=uamp->CurrentPos && uamp->CurrentPos<=stop){
		/*compute the BSP coordinates*/
		uamp->from_r=BspStart+(uamp->CurrentPos-start);
		if (uamp->DispType==DISP_JUST_SPLIT)
			uamp->from_r+=diff2;
		max=MIN(uamp->to_ua,stop);
		uamp->to_r=uamp->from_r+(max-uamp->CurrentPos);
		/*if minus strand (nuc. only); reverse the values*/
		if (uamp->strand==Seq_strand_minus){
			diff=uamp->to_r-uamp->from_r;
			uamp->from_r=uamp->from_r+(BspStop-uamp->from_r)-(uamp->from_r-BspStart);
			uamp->to_r=uamp->from_r-diff;
			tmp=uamp->to_r;
			uamp->to_r=uamp->from_r;
			uamp->from_r=tmp;
		}
		uamp->IsGap=FALSE;
		/*this will be used on the next function call*/
		uamp->CurrentPos=max+1;
	}
	/*he, he, he : we are in the gap !*/
	else{
		uamp->from_r=uamp->CurrentPos;
		/*be carefull : the right side of this segment could be in the
		middle of a gap...*/
		if (uamp->from_r<start){
			if (start<=uamp->to_ua){
				uamp->to_r=start-1;
				uamp->CurrentPos=start;
			}
			else{
				uamp->to_r=uamp->to_ua;
				uamp->CurrentPos=uamp->to_r+1;
			}
		}
		else{
			uamp->to_r=uamp->to_ua;
			uamp->CurrentPos=uamp->to_r+1;
		}
		uamp->IsGap=TRUE;
	}
	
	return(TRUE);
}

/*******************************************************************************

  Function : DDV_BuildBspEntitiesTbl()
  
  Purpose : build a table of eID, iID of each Bioseq of a SeqAlign.
    
  Return value : a table of entities

*******************************************************************************/
NLM_EXTERN Uint4Ptr DDV_BuildBspEntitiesTbl(ValNodePtr PNTR TableHead,Int4 nBsp)
{
ValNodePtr vnp;
ParaGPtr   pgp;
BioseqPtr  bsp;
Uint4Ptr    entitiesTbl;
Uint2      bsp_eID,bsp_iID;
Int4       i;

	entitiesTbl=(Uint4Ptr)MemNew(nBsp*sizeof(Uint4));
	if (!entitiesTbl) return(NULL);
	
	for (i=0;i<nBsp;i++){
		vnp=TableHead[i];
		if (vnp){
			pgp=(ParaGPtr)vnp->data.ptrvalue;
			
			bsp=BioseqLockById(pgp->sip);
			if (bsp){
				bsp_eID=ObjMgrGetEntityIDForPointer((Pointer)bsp);
				bsp_iID = GetItemIDGivenPointer (bsp_eID, 
						OBJ_BIOSEQ, (Pointer) bsp);	
				entitiesTbl[i]=UDV_EncodeIdxFeat(bsp_eID,bsp_iID);
				BioseqUnlockById(pgp->sip);
			}
		}
	}
	return(entitiesTbl);
}

static ValNodePtr DDV_BuildTailDescriptor(SeqAlignPtr sap, Int4 nBsp, Int2 LineSize,
	DDV_Disp_OptPtr ddop, Uint1 which_tail,Int4Ptr TotLength,
	Int4Ptr cumulpop)
{
Int4            startcopy,stopcopy,pop,i,MaxLength,start,stop;
Boolean         bError;
ValNodePtr      vnp,vnp_head;
DescriDispPtr   ddp;
Uint1           strand;

	vnp_head=NULL;
	bError=FALSE;
	MaxLength=0;
	/*get the size of the left UA region*/
	for(i=0;i<nBsp;i++){
		if (AlnMgrGetNthRowTail(sap,i+1,which_tail,&start,&stop,&strand)){
			if (start!=(Int4)-1 && stop!=(Int4)-1){
				MaxLength=MAX(MaxLength,ABS(stop-start)+1);
			}
		}
	}
	/*build a descriptor node for that region, if needed*/
	if (MaxLength>0){
		startcopy=0;
		stopcopy=MaxLength;
		*TotLength+=(stopcopy-startcopy);/*to be used for the display coord. system*/

		while(startcopy<stopcopy){
			/*store data*/
			pop=_min_(stopcopy-startcopy,LineSize-(*cumulpop));
			ddp=(DescriDispPtr)MemNew(sizeof(DescriDisp));
			if (ddp==NULL) {
				bError=TRUE;break;
			}
			ddp->from=startcopy;/*use relative coord, i.e. UA coord system*/
			ddp->to=ddp->from+pop-1;
			ddp->TextStyle=MSA_TXT_STYLE_REG_UNALIGN;
			if (which_tail==LEFT_TAIL)
				ddp->UAnum=(Int4)-1;/*UA region number; -1 means LEFT TAIL*/
			else
				ddp->UAnum=(Int4)-2;/*UA region number; -2 means RIGHT TAIL*/
			ddp->UAMaxlength=stopcopy;
			/*create a new node*/
			if (vnp_head==NULL) {
				vnp_head=ValNodeAddPointer(NULL,0,(Pointer)ddp);
				vnp=vnp_head;
			}
			else vnp=ValNodeAddPointer(&vnp,0,(Pointer)ddp);

			if (!vnp){
				bError=TRUE;
				break;
			}
			(*cumulpop)+=pop;
			if ((*cumulpop)==LineSize) (*cumulpop)=0;
			startcopy+=pop;
		}
	}

	if (bError){
		vnp_head=ValNodeFreeData(vnp_head);
	}
	
	return(vnp_head);
}

/*******************************************************************************

  Function : UABuildDescriptor()
  
  Purpose : build the descriptor of a SeqAlign. This function is ONLY designed 
            for discontinuous SeqAlign. The descriptor is a valnodelist. For
			each node, the field data.ptrvalue is of type MsaTxtDispPtr. It
			contains inportant info : start-stop of either an aligned or unaligned
			region. In the first case, SeqAlign coord system is used, otherwise
			UA coord system is used (UA=UnAligned). The descriptor should then
			be used as a guide to populate ParaG for the display. Note that this
			descriptor is a "pre"-layout of the ParaG structure : some descriptor's
			nodes put together = a ParaG in size.
    
  Return value : the descriptor (see also TotLength; total size of the SeqAlign,
	  display coordinates)

*******************************************************************************/
NLM_EXTERN ValNodePtr UABuildDescriptor(SeqAlignPtr sap, Int4 nBsp, Int2 LineSize,
	DDV_Disp_OptPtr ddop, Int4Ptr TotLength,Boolean AddLeftUAPart,
	Boolean AddRightUAPart)
{
Int4            length,r,cumulpop,StartLetter,startcopy,stopcopy,pop,UAnum;
Boolean         bUnAligned,bError;
ValNodePtr      vnp,vnp_head,vnp2;
DescriDispPtr   ddp;


	r=0;
	cumulpop=0;
	StartLetter=0;/*align coord*/
	bError=FALSE;
	vnp_head=NULL;
	*TotLength=0;
	/*MaxLength=0;*/


	/*build UnAligned left part, if needed*/
	if (AddLeftUAPart && ddop->DispDiscStyle!=MSA_TXT_STYLE_1){
		vnp_head=DDV_BuildTailDescriptor(sap,nBsp,LineSize,ddop,LEFT_TAIL,TotLength,
			&cumulpop);
		/*if (vnp_head==NULL) 
			goto erreur;*/
		/*get the last node; will be used in the SAP loop (see below)*/
		if (vnp_head){
			vnp=vnp_head;
			while(vnp->next){
				vnp=vnp->next;
			}
		}
	}

	/*loop through each aligned and unaligned regions to build
	a ValNode list of descriptors. These descriptors will be used to
	generate the ParaG list*/
	UAnum=0;
	while(AlnMgrGetNextLengthBit(sap,&length,&r)){

		startcopy=0;

		if (length<0){
			bUnAligned=TRUE;
			switch(ddop->DispDiscStyle){/*user's display choice*/
				case MSA_TXT_STYLE_1:/*put a spacer between 2 align blocks*/
					stopcopy=ddop->SpacerSize;
					break;
				case MSA_TXT_STYLE_2:/*put sequence between 2 align blocks*/
					stopcopy=ABS(length);
					break;
			}
			UAnum++;
		}
		else if (length>0) {
			bUnAligned=FALSE;
			stopcopy=length;
		} 
		else {
			bUnAligned=TRUE;
			stopcopy=length;
			UAnum++;
		}

		*TotLength+=(stopcopy-startcopy);/*to be used for the display coord. system*/
		
		while(startcopy<stopcopy){
			/*store data*/
			pop=_min_(stopcopy-startcopy,LineSize-cumulpop);
			ddp=(DescriDispPtr)MemNew(sizeof(DescriDisp));
			if (ddp==NULL) {
				bError=TRUE;break;
			}
			if (bUnAligned==FALSE){
				ddp->from=StartLetter;/*use align coord*/
				ddp->to=ddp->from+pop-1;
				ddp->TextStyle=MSA_TXT_STYLE_REG_ALIGN;
				StartLetter+=pop;
			}
			else{
				ddp->from=startcopy;/*use relative coord, i.e. UA coord system*/
				ddp->to=ddp->from+pop-1;
				ddp->TextStyle=MSA_TXT_STYLE_REG_UNALIGN;
				ddp->UAnum=UAnum;/*UA region number*/
				ddp->UAMaxlength=stopcopy;
			}			
			/*create a new node*/
			if (vnp_head==NULL) {
				vnp_head=ValNodeAddPointer(NULL,0,(Pointer)ddp);
				vnp=vnp_head;
			}
			else vnp=ValNodeAddPointer(&vnp,0,(Pointer)ddp);

			if (!vnp){
				bError=TRUE;
				break;
			}
			cumulpop+=pop;
			if (cumulpop==LineSize) cumulpop=0;
			startcopy+=pop;
		}
		if (bError) break;
	}

	if (bError){
		vnp_head=ValNodeFreeData(vnp_head);
		goto erreur;
	}

	/*build UnAligned right part, if needed*/
	if (AddRightUAPart && ddop->DispDiscStyle!=MSA_TXT_STYLE_1){
		vnp2=DDV_BuildTailDescriptor(sap,nBsp,LineSize,ddop,RIGHT_TAIL,TotLength,
			&cumulpop);
		/*add vnp2 at the end of vnp_head*/
		if (vnp2) {
			vnp=vnp_head;
			while(vnp->next){
				vnp=vnp->next;
			}
			vnp->next=vnp2;
		}
		/*else{
			vnp_head=ValNodeFreeData(vnp_head);
		}*/
	}

/*debug only : print the content of the descriptor*/

/*	vnp=vnp_head;
	while(vnp){
		ddp=(DescriDispPtr)vnp->data.ptrvalue;
		printf("[%4i..%4i],%2i (%4i) ; ",ddp->from,ddp->to,ddp->TextStyle,
			ddp->to-ddp->from+1);
		if (ddp->TextStyle==MSA_TXT_STYLE_REG_UNALIGN){
			printf("MaxLength = %i ; UAnum = %i\n",ddp->UAMaxlength,ddp->UAnum);
		}
		else{
			printf("\n");
		}
		vnp=vnp->next;
	}
*/
/*debug only - end*/

erreur :
	return(vnp_head);
}

/*******************************************************************************

  Function : DDV_DumpSAPInAFile()
  
  Purpose : write an Indexed SeqALign in a file (use the code of the Population
      Viewer).
  
  Return value : -

*******************************************************************************/
NLM_EXTERN void DDV_DumpSAPInAFile(MsaParaGPopListPtr mpplp,DDVOptionsBlockPtr dobp, 
		FILE *hFile,Uint4 option,DDV_ColorGlobal * gclr)
{
ByteStorePtr    PNTR byteSpp;
Int4Ptr PNTR    matrix  ;
Int2            LineSize;

	if (!mpplp || !hFile) return;
	
	/*get the optional data*/
	if (dobp){
		LineSize=dobp->LineSize;
		byteSpp=dobp->byteSpp;
		matrix=dobp->matrix;
	}
	else{
		LineSize=ParaG_Size;
		byteSpp=NULL;
	}

	/*do the display */
	DDV_AffichageParaG(mpplp,0,0,mpplp->LengthAli-1,mpplp->LengthAli,0,option,
		LineSize,hFile,byteSpp,matrix,gclr,NULL);

}

/*******************************************************************************

  Function : DDV_DumpSAPInFastaFile()
  
  Purpose : write an Indexed SeqALign in a file (use FASTA format only).
  
  Return value : -

*******************************************************************************/
NLM_EXTERN void DDV_DumpSAPInFastaFile(MsaParaGPopListPtr mpplp,DDVOptionsBlockPtr dobp, 
		FILE *fp,Boolean bPrintGap)
{
CharPtr    szBuf,szSeq;
ParaGPtr   pgp;
SeqIdPtr   sip;
ValNodePtr vnp;
BioseqPtr  bsp;
Int4       i,j,n,LineSize,bspLength;
Boolean    bFirst,IsAA;
Char       letter;
Char       DefLine[255];

	if (!mpplp || !fp) return;

	if (dobp){
		LineSize=dobp->LineSize;
	}
	else{
		LineSize=ParaG_Size;
	}
	
	szBuf=(CharPtr)MemNew((LineSize+1)*sizeof(Char));
	if (!szBuf) return;/*todo : error msg for the user*/
	
	for(i=0;i<mpplp->nBsp;i++){
		vnp=mpplp->TableHead[i];
		n=0;
		bFirst=TRUE;
		while(vnp){
			pgp=(ParaGPtr)vnp->data.ptrvalue;
			sip=SeqIdFindBest(pgp->sip,0);
			if (sip==NULL) sip=pgp->sip;
			bsp=BioseqLockById(sip);
			if (!bsp) return;/*todo : error msg for the user*/
			bspLength=BioseqGetLen(bsp);
			IsAA=ISA_aa(bsp->mol);
			if (bFirst){
				FastaId(bsp,DefLine,254);
				fprintf(fp,">%s ",DefLine);
				CreateDefLine(NULL, bsp, DefLine, 254, 0,NULL, NULL);
				fprintf(fp,"%s\n",DefLine);
				bFirst=FALSE;
			}
			BioseqUnlockById(sip);
			szSeq=(CharPtr)MemNew((pgp->StopLetter-pgp->StartLetter+2)*sizeof(Char));
			if (!szSeq) return;/*todo : error msg for the user*/
			DDV_GetSequenceFromParaG(pgp,&szSeq,bspLength,IsAA,NULL,NULL,NULL);
			j=0;
			while(szSeq[j]){
				letter=szSeq[j++];
				if (isalpha(letter)){
					szBuf[n++]=letter;
				}
				else{
					if (bPrintGap)
						szBuf[n++]=letter;
				}
				if (n==LineSize){
					szBuf[n]='\0';
					fprintf(fp,"%s\n",szBuf);
					n=0;
				}
			}
			MemFree(szSeq);
			vnp=vnp->next;
			if (vnp==NULL && n!=0){
				szBuf[n]='\0';
				fprintf(fp,"%s\n",szBuf);
			}
		}
	}
	MemFree(szBuf);
}

/*******************************************************************************

  Function : DDV_InitDefSAPdispStyles_BLAST()
  
  Purpose : set the default styles for SeqAligns.
  
  Return value : -

*******************************************************************************/
static void DDV_InitDefSAPdispStyles_BLAST(DDV_Disp_OptPtr ddop)
{
	/*use colors*/
	ddop->bUseColors=FALSE;
	/*disc SAP styles*/
	ddop->ShowLeftTail=FALSE;
	ddop->ShowRightTail=FALSE;
	ddop->DispDiscStyle=MSA_TXT_STYLE_2;
	ddop->SpacerSize=SPACER_TXT_BLANK;
    ddop->DiscJustification=DISP_JUST_LEFT;

	/*ruler style*/
	ddop->RulerStyle=RulerStyle_Continue_Start;
}


/*******************************************************************************

  Function : DDV_DisplayNewBLAST()
  
  Purpose : build a default SeqAlign display for BLAST.
  
  Parameters : 	sap; SeqAlign
				disp_format; display formats (see pgppop.h)
				disp_data; for future implementation
				fp; to send the output

  NOTE : this function uses the alignmgr.[ch] instead of blocks.[ch]
  
  Return value : FALSE if failed 

*******************************************************************************/
static Boolean DDV_DisplayNewBLAST(SeqAlignPtr sap, ValNodePtr mask,
		Int4Ptr PNTR matrix,Uint4 disp_format,Pointer disp_data,FILE *fp)
{
MsaParaGPopList mppl;/*contains the data for the display*/
Int2            LineSize;
DDVOptionsBlockPtr dobp;
ByteStorePtr   PNTR byteSpp;
DDV_ColorGlobal *layout;
DDV_Disp_Opt   ddo;
Boolean        bUseLayout;
Uint1          what,mask_style;

	MemSet(&mppl,0,sizeof(MsaParaGPopList));
	layout=NULL;

	/*get the optional data*/
	if (disp_data){
		dobp=(DDVOptionsBlockPtr)disp_data;
		LineSize=dobp->LineSize;
		byteSpp=dobp->byteSpp;
		bUseLayout=dobp->bUseLayout;
		what=dobp->LayoutType;
	}
	else{
		LineSize=ParaG_Size;
		byteSpp=NULL;
		bUseLayout=FALSE;
	}

	DDV_InitDefSAPdispStyles_BLAST(&ddo);
	
	/*build the indexed SeqAlign and the ParaG structure*/
	DDV_CreateDisplayFromIndex(sap, &mppl, LineSize,&ddo);
	
	/*create the layout*/
	if (bUseLayout){
		layout = DDV_CreateColorGlobal(TRUE);
		if (layout){
			if (dobp->LayoutType&DDV_USE_STDCLR){
				DDV_InitColour_When_Start(sap,&mppl,&layout, FALSE);
			}
			if (dobp->LayoutType&DDV_USE_UALAYOUT){
				DDV_LayoutUAregion(sap,&mppl,&layout,
					dobp->UAlayout.rgb,dobp->UAlayout.style);
			}
			if (dobp->LayoutType&DDV_USE_ISOLAYOUT){
				DDV_LayoutISOColors(layout,mppl.TableHead,mppl.nBsp,
					0,TRUE,matrix,dobp->ISOlayout.clr_ident,
					dobp->ISOlayout.clr_simil,dobp->ISOlayout.clr_other);
			}
			if (dobp->LayoutType&DDV_USE_GAPLAYOUT){
				DDV_LayoutGap(layout,mppl.TableHead[0],mppl.TableHead[1],
					dobp->GAPlayout.style,dobp->GAPlayout.rgb);
				DDV_LayoutGap(layout,mppl.TableHead[1],mppl.TableHead[0],
					dobp->GAPlayout.style,dobp->GAPlayout.rgb);
			}
			if (dobp->LayoutType&DDV_USE_MASKLAYOUT){
				DDV_LayoutMaskRegions(layout,mppl.TableHead[0],mask,
					dobp->MASKlayout.style,dobp->MASKlayout.rgb);
			}
		}
	}
	
	/*do the display */
	DDV_AffichageParaG(&mppl,0,0,mppl.LengthAli-1,mppl.LengthAli,0,disp_format,
		LineSize,fp,byteSpp,matrix,layout,mask);

	/*done... delete data*/
	DDV_DeleteDisplayList(&mppl);
	if (layout) DDV_DeleteColorGlobal(layout);

	
	return(TRUE);
}


/*******************************************************************************

  Function : DDV_DisplayBlastSAP()
  
  Purpose : test function to display a BLAST output.
  
  Parameters : 	sap; seqalign
			    mask; list of masked regions in the query
				fo; output file;
				is_na; TRUE means nuc sequence
				tx_option; some display options
				
  Return value : -

*******************************************************************************/
NLM_EXTERN Boolean DDV_DisplayBlastSAP(SeqAlignPtr sap, ValNodePtr mask,
	FILE *fp, Boolean is_na, Uint4 tx_option)
{
DDVOptionsBlock dob;
SeqAlignPtr     sap4;
SeqIdPtr        new_id = NULL, old_id = NULL;    
Int4Ptr PNTR    matrix;
Uint4           option,i;
Boolean         bRet,bUseLayout, follower= FALSE;

    if (!sap || !fp) 
        return(FALSE);
    
	/*get the matrix*/
    if (is_na == FALSE){
        matrix=load_default_matrix();
        if (!matrix)
            return(FALSE);
    } else {
        matrix=NULL;
    }
    
    /*init display format*/
    option =VIEW_FULLSEQ;
    option|=DISP_FULL_HTML;
    option|=DISP_BSP_COORD;
    option|=SEQNAME_BLAST_STD;
    option|=DISP_BLAST_STD;
	option|=DISP_BLAST_MIDLINE;
    option|=DISPE_COLOR;
    dob.LineSize=(Int2)60;
    dob.matrix=matrix;
	
    bRet=TRUE;
    sap4=sap;
    while(sap4) {
        /*build the Index*/
        if (sap4->segtype == SAS_DISC){
            if (!sap4 || !AlnMgrIndexSingleSeqAlign(sap4)){
                bRet=FALSE;
                break;
            }
            dob.bUseLayout=TRUE;
        }
        else{
            if (!sap4 || !AlnMgrIndexSingleChildSeqAlign(sap4)){
                bRet=FALSE;
                break;
            }
            dob.bUseLayout=TRUE;
        }
		
		/*set the layout type, if needed*/
		if (dob.bUseLayout){
			/*basic layout*/
			dob.LayoutType=DDV_USE_UALAYOUT;
			dob.LayoutType|=DDV_USE_GAPLAYOUT;
			dob.LayoutType|=DDV_USE_MASKLAYOUT;
			/*dob.LayoutType|=DDV_USE_ISOLAYOUT;*/
			/*ext. layout*/
			dob.UAlayout.style=DDV_TXTSTYLE_LOWERCASE;
			dob.UAlayout.style|=DDV_TXTSTYLE_ITALIC;
			dob.UAlayout.style|=DDV_TXTSTYLE_COLOR;
			dob.UAlayout.rgb[0]=255;
			dob.UAlayout.rgb[1]=0;
			dob.UAlayout.rgb[2]=0;
			dob.GAPlayout.style=DDV_TXTSTYLE_LOWERCASE;
			dob.GAPlayout.style|=DDV_TXTSTYLE_ITALIC;
			dob.GAPlayout.style|=DDV_TXTSTYLE_COLOR;
			dob.GAPlayout.rgb[0]=255;
			dob.GAPlayout.rgb[1]=0;
			dob.GAPlayout.rgb[2]=0;
			dob.MASKlayout.style=DDV_TXTSTYLE_BOLD;
			dob.MASKlayout.style|=DDV_TXTSTYLE_UNDERLINE;
			dob.ISOlayout.clr_ident=DDVCOL_ORANGE;
			dob.ISOlayout.clr_simil=DDVCOL_SKY;
			dob.ISOlayout.clr_other=DDVCOL_BLACK;
		}
					
        if (option&DISP_FULL_HTML){
			fprintf(fp,"<pre>\n");
		}
        /* Attempt to print score for the alignment */
        new_id = TxGetSubjectIdFromSeqAlign(sap4);
        if(old_id != NULL) {
            if(SeqIdMatch(new_id, old_id))
                follower = TRUE;
        }
		
		/*separator*/
		fprintf(fp,"<HR WIDTH=\"400\">");
        old_id = new_id;
        if(!FormatScoreFromSeqAlign(sap4, tx_option, fp, matrix, follower)){
            bRet=FALSE;
            break;
        }
        follower = FALSE;
        if (option&DISP_FULL_HTML){
			fprintf(fp,"</pre>\n");
		}

        /*display a SeqAlign*/
		if (!DDV_DisplayNewBLAST(sap4, mask, matrix, option, (Pointer) &dob, fp)){
            bRet=FALSE;
            break;
        }
        sap4 = sap4->next;
    }
    
    if (matrix){
        for(i = 0; i<TX_MATRIX_SIZE; ++i)
            MemFree(matrix[i]);
        MemFree(matrix);
    } 
    
    return(bRet);
}
/*******************************************************************************

  Function : DDV_InitDefSAPdispStyles()
  
  Purpose : set the default styles for SeqAligns.
  
  Return value : -

*******************************************************************************/
NLM_EXTERN void DDV_InitDefSAPdispStyles(DDV_Disp_OptPtr ddop)
{
	/*use colors*/
	ddop->bUseColors=TRUE;
	/*disc SAP styles*/
	ddop->ShowLeftTail=FALSE;
	ddop->ShowRightTail=FALSE;
	ddop->DispDiscStyle=MSA_TXT_STYLE_1;
	ddop->SpacerSize=SPACER_TXT_BLANK;
    ddop->DiscJustification=DISP_JUST_CENTER;

	/*ruler style*/
	ddop->RulerStyle=RulerStyle_Continue_Start;
}

/*******************************************************************************

  Function : DDV_InitCn3DSAPdispStyles()
  
  Purpose : set the styles for Cn3D SeqAligns.
  
  Return value : -

*******************************************************************************/
NLM_EXTERN void DDV_InitCn3DSAPdispStyles(DDV_Disp_OptPtr ddop)
{
	/*use colors*/
	ddop->bUseColors=TRUE;
	/*disc SAP styles*/
	ddop->ShowLeftTail=TRUE;
	ddop->ShowRightTail=TRUE;
	ddop->DispDiscStyle=MSA_TXT_STYLE_2;
	ddop->SpacerSize=0;
    ddop->DiscJustification=DISP_JUST_SPLIT;

	/*ruler style*/
	ddop->RulerStyle=RulerStyle_Continue_Start;
}

/*******************************************************************************

  Function : DDV_ReadSeqBin()
  
  Purpose : read a sequence using Seq_code_ncbistdaa or Seq_code_ncbi4na
    
  Return value : the sequence 

*******************************************************************************/
static Uint1Ptr DDV_ReadSeqBin (SeqIdPtr sip, Int4 from, Int4 to, 
		Boolean IsProt,Int2 len,Uint1 strand)
{
SeqLocPtr  slp;
SeqPortPtr spp;
Uint1Ptr    btr;
Uint1	   residue;
Uint2	   i=0;

	/*from always < than to*/
	slp = SeqLocIntNew (from, to, strand, sip);
	if (!slp) return(NULL);
	spp = SeqPortNewByLoc (slp, (Uint1)(IsProt==TRUE ? Seq_code_ncbistdaa 
			: Seq_code_ncbi4na));
	if (spp != NULL) {
		btr = (Uint1Ptr) MemNew (len * sizeof(Uint1));
		if (!btr) return(NULL);
		while ((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF) {
			if (IS_residue(residue)) {
				btr[i] = residue;
				i++;
			}
		}
		/*SeqPortRead(spp, btr, len);*/
		SeqPortFree (spp);
	}   
	SeqLocFree (slp);
	return btr;
}

/*******************************************************************************

  Function : DDV_SetStyle()
  
  Purpose : set the layout of a single letter 
  
  Return value : -.

*******************************************************************************/
static void DDV_SetStyle(DDV_ColorCell * dccp,Uint1 style, Uint1 *rgb)
{
	if (style&DDV_TXTSTYLE_ITALIC)
		dccp->UseItalic=TRUE;
    else dccp->UseItalic=FALSE;
	if (style&DDV_TXTSTYLE_BOLD)
		dccp->UseBold=TRUE;
    else dccp->UseBold=FALSE;
	if (style&DDV_TXTSTYLE_UNDERLINE)
		dccp->UseUnderline=TRUE;
    else dccp->UseUnderline=FALSE;
	if (style&DDV_TXTSTYLE_LOWERCASE)
		dccp->LowerCase=TRUE;
    else dccp->LowerCase=FALSE;
	if (style&DDV_TXTSTYLE_COLOR){
		dccp->rgb[0]=rgb[0];
		dccp->rgb[1]=rgb[1];
		dccp->rgb[2]=rgb[2];
	}
}

/*******************************************************************************

  Function : DDV_SetStd_AA_UA_ClrEx()
  
  Purpose : set the standard colours for a nuc or prot sequence; 
     or the same function can be used to layout UnAligned regions
    
  Return value : -

*******************************************************************************/
#define read_buf_size 3000
static void DDV_SetStd_AA_UA_ClrEx(DDV_ColorGlobal *pColorGlobal, void *pData,
                                DDV_Range *pRange,Uint1 what)
{
DDVDataForColorFuncPtr ddfcfp;
BioseqPtr              bsp;
Int4                   bsp_length,i,from,to,limit;
Boolean                IsAA;
Int4                   dp;
Uint1Ptr               SeqBuf;
Int2                   len;
Uint1                  idx;
DDV_ColorCell *        dccp;
Uint1                  rgb[3];

	ddfcfp=(DDVDataForColorFuncPtr)pData;
	
	if (!ddfcfp) return;
	
	/*get bioseq size and type */
	bsp=BioseqLockById(ddfcfp->sip);
	if (bsp){
		bsp_length=bsp->length;
		IsAA=ISA_aa(bsp->mol);
		BioseqUnlock(bsp);
	}
	
	/*compute the ends of the bioseq's region we want to analyse*/
	if (ddfcfp->from!=(Int4)-1)
		from=ddfcfp->from;
	else
		from=0;

	if (ddfcfp->to!=(Int4)-1)
		limit=ddfcfp->to;
	else
		limit=bsp_length-1;

	to=from+read_buf_size-1;
	if (to>limit) to=limit;
	
	/*read the sequence by chunk of 'read_buf_size' letters*/
	while(TRUE){
		len=(Int2)(to-from+1);
		/*get a sequence chunck ; from is always < than to*/
		SeqBuf=DDV_ReadSeqBin (ddfcfp->sip, from, to, IsAA,len,ddfcfp->strand);
		if (SeqBuf){/*scan each letter and get the appropriate colour*/
			for(i=0;i<len;i++){
				/*compute the correct position; depend on the strand*/
				if (ddfcfp->strand==Seq_strand_minus)
					dp=to-i;
				else
					dp=from+i;
				/*first, query Color Manager*/
                dccp=DDV_GetColor(pColorGlobal,ddfcfp->sip,-1,dp);
                if (!dccp) return;
                if(what & DDV_USE_STDCLR) {
                    if (IsAA) 
                        idx=DDV_STD_AAColor[*(SeqBuf+i)];/*AA*/
                    else
                        idx=DDV_STD_NAColor[*(SeqBuf+i)];/*NA*/
                    
                    if (ddfcfp->IsUnAligned)DDV_SetStyle(dccp,ddfcfp->style 
                        | DDV_TXTSTYLE_LOWERCASE,ddfcfp->rgb);
                    else DDV_SetStyle(dccp,ddfcfp->style,DDV_PaletteRGB[idx]);
                }
                else if(what & DDV_USE_UALAYOUT) {
                    if (ddfcfp->IsUnAligned)DDV_SetStyle(dccp,ddfcfp->style 
                        | DDV_TXTSTYLE_LOWERCASE,ddfcfp->rgb);
                    else DDV_SetStyle(dccp,ddfcfp->style,DDV_PaletteRGB[DDVCOL_BLACK]);
                }				
                /*set the colour*/
                DDV_SetColor(pColorGlobal, ddfcfp->sip, -1, dp, dccp);
            }
            MemFree(SeqBuf);
        }
		
		/*compute the next chunk*/		
		from+=read_buf_size;
		to+=read_buf_size;
		if (to>limit) to=limit;
		/*end ?*/
		if (from>limit) break;
	}

}


/*******************************************************************************

  Function : DDV_InitColourSAP()
  
  Purpose : init the colours for each Bsp of a SeqAlign
  
  Parameters : sap; a (indexed) SeqAlign
  		global colour data struct (see ddvcolor.[ch] in api)
        NoColor; don't color the bsp's, but do change the case.
		
  Note : usually, this function is called when DDV loads a SeqAlign
  
  Return value : FALSE is failure.

*******************************************************************************/
static Boolean DDV_InitColourSAP(SeqAlignPtr sap,MsaParaGPopListPtr mpplp,
	DDV_ColorGlobal * * gclr,Uint1 what,Uint1 * pRGB,Uint1 style)
{
DDVDataForColorFunc  ddfcf;
SeqIdPtr             sip;
ValNodePtr           vnp,vnp2;
ParaGPtr             pgp;
MsaTxtDispPtr        mtdp;
Boolean              cont,bAligned;
Int4                 n;
Int4                 start;
Int4                 stop;
Int4                 tmp;
Int4                 nRow;

    if (!sap)
    return (FALSE);

    if (!(*gclr)){
        *gclr = DDV_CreateColorGlobal(FALSE);
        if (!(*gclr))
            return(FALSE);
    }

    for(n=0;n<mpplp->nBsp;n++){
        vnp=mpplp->TableHead[n];
        while(vnp){
            pgp=(ParaGPtr)vnp->data.ptrvalue;
            if (pgp){
                vnp2=pgp->ptxtList;
                sip=pgp->sip;
                while(vnp2){
                    mtdp = (MsaTxtDispPtr)vnp2->data.ptrvalue;
                    if(mtdp && mtdp->TextStyle==MSA_TXT_STYLE_SEQ){
                        ddfcf.sip=sip;
                        ddfcf.from=mtdp->from;
                        ddfcf.to=mtdp->to;
                        ddfcf.strand=mtdp->strand;
                        ddfcf.IsUnAligned=mtdp->IsUnAligned;
                        ddfcf.style=style;
                        if (pRGB){
                            ddfcf.rgb[0]=pRGB[0];
                            ddfcf.rgb[1]=pRGB[1];
                            ddfcf.rgb[2]=pRGB[2];
                        }
                        else{
                            ddfcf.rgb[0]=0;
                            ddfcf.rgb[1]=0;
                            ddfcf.rgb[2]=0;
                        }
                        DDV_SetStd_AA_UA_ClrEx(*gclr, (Pointer)&ddfcf, NULL, what);
                    }
                    vnp2=vnp2->next;
                }
            }
            vnp=vnp->next;
        }
    }
    return(TRUE);
}

/*******************************************************************************

  Function : DDV_InitColour_When_Start()
  
  Purpose : standard init function for Vibrant DDV ; see DDV_InitColourSAP()
  
*******************************************************************************/
NLM_EXTERN Boolean DDV_InitColour_When_Start(SeqAlignPtr sap,MsaParaGPopListPtr mpplp,
	DDV_ColorGlobal * * gclr, Boolean NoColor)
{
    if(NoColor) return(DDV_InitColourSAP(sap,mpplp,gclr,DDV_USE_STDCLR,NULL,0));
    else return(DDV_InitColourSAP(sap,mpplp,gclr,DDV_USE_STDCLR,NULL,
        DDV_TXTSTYLE_COLOR));
}

/*******************************************************************************

  Function : DDV_LayoutUAregion()
  
  Purpose : set a layout for UnAligned regions; see DDV_InitColourSAP()
  
*******************************************************************************/
NLM_EXTERN Boolean DDV_LayoutUAregion(SeqAlignPtr sap,MsaParaGPopListPtr mpplp,
	DDV_ColorGlobal * * gclr,Uint1 * pRGB,Uint1 style)
{
	return(DDV_InitColourSAP(sap,mpplp,gclr,DDV_USE_UALAYOUT,pRGB,style));
}


/*******************************************************************************

  Function : UDV_BigEncodeIdxFeat4()
  
  Purpose : encode two Uint4 into one Uint8
  
*******************************************************************************/
static Uint8 UDV_BigEncodeIdxFeat4 (Uint4 val1,Uint4 val2)
{
Uint4 index_g[2];
	
	index_g[0]=val1;
	index_g[1]=val2;
	
	return *((Uint8Ptr) index_g);
	
}

/*******************************************************************************

  Function : UDV_BigEncodeIdxFeat4()
  
  Purpose : decode one Uint8 into two Uint4
  
*******************************************************************************/
static void  UDV_BigDecodeIdxFeat4 (Uint8 index_g, Uint4Ptr val1, Uint4Ptr val2)
{
Uint4Ptr  index_g2;

	index_g2 = (Uint4Ptr) (&index_g);
	if (val1) *val1 = (Uint4) index_g2 [0];
	if (val2) *val2 = (Uint4) index_g2 [1];
}

/*******************************************************************************

  Function : DDV_GetGapCoord()
  
  Purpose : get the list of gaps given a ParaG.
  
  Return value : -.

*******************************************************************************/
static void DDV_GetGapCoord(ParaGPtr pgp,ValNodePtr PNTR gap_list)
{
ValNodePtr    vnp;
MsaTxtDispPtr mtdp;

	if (!pgp) return;
	vnp=pgp->ptxtList;
	while(vnp){
		mtdp=(MsaTxtDispPtr)vnp->data.ptrvalue;
		if (mtdp){
			if(mtdp->IsGap && (mtdp->TextStyle==MSA_TXT_STYLE_GAP || 
				mtdp->TextStyle==MSA_TXT_STYLE_2)){
				ValNodeAddBigInt(gap_list,0,
					UDV_BigEncodeIdxFeat4 ((Uint4)mtdp->from,(Uint4)mtdp->to));
			}
		}
		vnp=vnp->next;
	}
}

/*******************************************************************************

  Function : DDV_AllGapInLowerCase()
  
  Purpose : analyse sequences by pair and switch to lower case all residues
            "aligned" with a gap.
  
  Return value : -.

*******************************************************************************/
NLM_EXTERN void DDV_LayoutGap(DDV_ColorGlobal *pcg,
	ValNodePtr vnp_seq1,ValNodePtr vnp_seq2,Uint1 style,Uint1 *rgb)
{
DDV_ColorCell * dccp;
ValNodePtr      vnp1,vnp2,vnp3,vnp4;
ParaGPtr        pgp1,pgp2;
Uint4           i,from,to;/*display coord.*/
Int4            bsp_coord;
DDV_ColorCell   dcc;

	if (pcg==NULL || vnp_seq1==NULL || vnp_seq2==NULL) return;
	/*scan the first sequence to modify the second*/
	vnp1=vnp_seq1;
	vnp2=vnp_seq2;
	vnp3=NULL;
	while(vnp1){
		/*is there a gap in the vnp1 ParaG*/
		pgp1=(ParaGPtr)vnp1->data.ptrvalue;
		pgp2=(ParaGPtr)vnp2->data.ptrvalue;
		DDV_GetGapCoord(pgp1,&vnp3);
		if (vnp3){
			vnp4=vnp3;
			while(vnp4){
				UDV_BigDecodeIdxFeat4((Uint8)vnp4->data.bigintvalue, 
					&from, &to);
				to++;
				for (i=from;i<to;i++){
					bsp_coord=(Int4)DDV_GetBspCoordGivenDispCoord(pgp2,(Int4)i);
					if (bsp_coord>=0){
						dccp=DDV_GetColor(pcg,pgp2->sip,-1,bsp_coord);
						if (dccp){/*if there is already a color, just
						  update the required parameter*/
							DDV_SetStyle(dccp,style,rgb);
							DDV_SetColor(pcg,pgp2->sip,-1,bsp_coord,dccp);
						}
						else{/*otherwise, set up a new layout for that letter*/
							memset(&dcc,0,sizeof(DDV_ColorCell));
							DDV_SetStyle(&dcc,style,rgb);
							DDV_SetColor(pcg, pgp2->sip, -1, bsp_coord, &dcc);
						}
					}
				}
				vnp4=vnp4->next;
			}
		}
		vnp1=vnp1->next;
		vnp2=vnp2->next;
		if (vnp3) 
			vnp3=ValNodeFree(vnp3);
	}
}

/*******************************************************************************

  Function : DDV_LayoutMaskRegions()
  
  Purpose : layout the masked regions of a sequence (usually the query sequence
       of BLAST).
  
  Return value : -.

*******************************************************************************/
NLM_EXTERN void DDV_LayoutMaskRegions(DDV_ColorGlobal *pcg,
	ValNodePtr vnp_query,ValNodePtr mask,Uint1 style,Uint1 *rgb)
{
DDV_ColorCell * dccp;
ValNodePtr      vnp;
ParaGPtr        pgp;
SeqIdPtr        sip;
SeqLocPtr       slp;
Int4            i,bsp_coord,bsp_start,bsp_stop,mask_start,mask_stop,scan_from,
                scan_to;
DDV_ColorCell   dcc;

	if (pcg==NULL || vnp_query==NULL || mask==NULL) return;

	/*get the sip*/
	pgp=(ParaGPtr)vnp_query->data.ptrvalue;
	if (!pgp) return;
	sip=pgp->sip;

	/*get query sequence range*/
	UDV_GetBspRangeinPGPList(vnp_query,&bsp_start,&bsp_stop);

	/*scan the mask list*/
	vnp=mask;
	while(vnp){
		slp=(SeqLocPtr)vnp->data.ptrvalue;
		if (slp){
			mask_start=SeqLocStart(slp);
			mask_stop=SeqLocStop(slp);
			/*check if the masked region is in the bsp range*/
			if (mask_stop>=bsp_start && mask_start<=bsp_stop){
				scan_from=_max_(mask_start,bsp_start);
				scan_to=_min_(mask_stop,bsp_stop)+1;
				for(i=scan_from;i<scan_to;i++){
					dccp=DDV_GetColor(pcg,sip,-1,i);
					if (dccp){/*if there is already a color, just
					  update the LowerCase parameter*/
						DDV_SetStyle(dccp,style,rgb);
						DDV_SetColor(pcg,sip,-1,i,dccp);
					}
					else{/*otherwies, set a new layout for that letter*/
						memset(&dcc,0,sizeof(DDV_ColorCell));
						DDV_SetStyle(&dcc,style,rgb);
						DDV_SetColor(pcg, sip, -1, i, &dcc);
					}
				}
			}
		}
		vnp=vnp->next;
	}
}

/*******************************************************************************

  Function : DDV_SetISOClr()
  
  Purpose : layout the colors for ident/simil/other (I/S/O).
  
*******************************************************************************/
static void DDV_SetISOClr(DDV_ColorGlobal * pcg, ParaGPtr pgp,Int4 disp_coord, Uint1 idx)
{					
DDV_ColorCell * dccp;
DDV_ColorCell   dcc;
Int4            bsp_coord;

	bsp_coord=DDV_GetBspCoordGivenDispCoord(pgp,disp_coord);
	if (bsp_coord>=0){
		dccp=DDV_GetColor(pcg,pgp->sip,-1,bsp_coord);
		if (dccp){/*if there is already a color, just
		  update the required parameter*/
			dccp->rgb[0]=DDV_PaletteRGB[idx][0];
			dccp->rgb[1]=DDV_PaletteRGB[idx][1];
			dccp->rgb[2]=DDV_PaletteRGB[idx][2];
			DDV_SetColor(pcg,pgp->sip,-1,bsp_coord,dccp);
		}
		else{/*otherwise, set up a new layout for that letter*/
			memset(&dcc,0,sizeof(DDV_ColorCell));
			dcc.rgb[0]=DDV_PaletteRGB[idx][0];
			dcc.rgb[1]=DDV_PaletteRGB[idx][1];
			dcc.rgb[2]=DDV_PaletteRGB[idx][2];
			DDV_SetColor(pcg, pgp->sip, -1, bsp_coord, &dcc);
		}
	}
}

/*******************************************************************************

  Function : DDV_LayoutIdentColors()
  
  Purpose : layout the colors for ident/simil/other (I/S/O).
  
  Return value : -.

*******************************************************************************/
NLM_EXTERN void DDV_LayoutISOColors(DDV_ColorGlobal *pcg,ValNodePtr * row_list,Int4 nRow,
	Int4 Master,Boolean bSetMaster,Int4Ptr * matrix,Uint1 IdentClr,Uint1 SimilClr,
	Uint1 OtherClr)
{
ValNodePtr * row;
ValNodePtr vnp;
ParaGPtr  pgpQuery,pgpSubject;
CharPtr   szQuery,szComp;
BioseqPtr bsp;
Boolean   IsAA;
Int4      i,j,len,bspLength;
Uint1     idx;

	row=(ValNodePtr *)MemNew(nRow*sizeof(ValNodePtr));
	if (!row) return;
	for (i=0;i<nRow;i++){
		row[i]=row_list[i];
	}

	vnp=row_list[Master];
	while(vnp){
		/*get the master sequence*/
		pgpQuery=(ParaGPtr)vnp->data.ptrvalue;
		szQuery=(CharPtr)MemNew((pgpQuery->StopLetter-pgpQuery->StartLetter+3)*sizeof(Char));
		if (!szQuery) goto error;
		bsp=BioseqLockById(pgpQuery->sip);
		if (!bsp) goto error;
		bspLength=BioseqGetLen(bsp);
		IsAA=ISA_aa(bsp->mol);
		BioseqUnlockById(pgpQuery->sip);
		if (!DDV_GetSequenceFromParaG(pgpQuery,&szQuery,bspLength,IsAA,NULL,
				NULL,NULL)) goto error;
		/*lopp on each row, get the sequence, set the colors*/
		for (i=0;i<nRow;i++){
			if (i==Master) continue;
			pgpSubject=(ParaGPtr)row[i]->data.ptrvalue;
			szComp=DDV_GetBLASTCompLine_3(szQuery,pgpSubject,matrix);
			if (szComp){
				len=StringLen(szComp)+1;
				for (j=0;j<len;j++){
					if (isalpha(szComp[j]) || szComp[j]=='|'){
						idx=IdentClr;
					}
					else if (szComp[j]=='+'){
						idx=SimilClr;
					}
					else{
						idx=OtherClr;
					}
					DDV_SetISOClr(pcg,pgpSubject,pgpSubject->StartLetter+j,idx);
					if (bSetMaster)
						DDV_SetISOClr(pcg,pgpQuery,pgpQuery->StartLetter+j,idx);
				}
				MemFree(szComp);
			}
			row[i]=row[i]->next;
		}
		szQuery=MemFree(szQuery);
		vnp=vnp->next;
	}
error:

	MemFree(row);
	if (szQuery) MemFree(szQuery);
}

