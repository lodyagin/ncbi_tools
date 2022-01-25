/*   dotvibrant.c
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
* File Name:  dotvibrant.c
*
* Author:  Fasika Aklilu
*
* Version Creation Date:   7/5/00
*
* $Revision: 6.10 $
*
* File Description: mouse management, graphic engine of the sequence viewer
*                   part of this code is also used for the WWW Entrez viewer
*                   (WWW_UDV define)
* Modifications:
* --------------------------------------------------------------------------
* $Log: dotvibrant.c,v $
* Revision 6.10  2000/08/07 16:34:51  kans
* added public domain notice
*
Revision 6.9  2000/08/07 13:46:59  sicotte
added revision

Revision 6.8  2000/08/07 13:46:34  sicotte
added revision

Revision 6.7  2000/08/07 13:46:05  sicotte
added Version logging

Revision 6.6  2000/08/07 13:45:22  sicotte
fixed (long) cast in sprintf and fprintf
*
*
* ==========================================================================
*/

/* dotvibrant.c */

#include <dotviewer.h>


 /****************************************************************************

     DEFINES                                                            
 ***************************************************************************/

#define MAXZOOMSCALEVAL 23
#define BLOCK_SIZE 10
#define VIS_LEN 100

#define dot_SEQVIEW 1
#define dot_FEATVIEW 2
 /****************************************************************************

     GLOBAL VARIABLES                                                             
 ***************************************************************************/

static CharPtr  na_names[] = {"A", "C", "G", "T", "-"};
static CharPtr  aa_names [] =  {"-", "A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W","X", "Y", "Z", "U" , "*"};
static Int4  zoomScaleVal [MAXZOOMSCALEVAL] = {
  1L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 20L,
  30L, 40L, 50L, 60L, 70L, 80L, 90L, 100L, 200L, 500L, 1000L
};
static PrimitivE prim_prev=NULL;

 /****************************************************************************

      FUNCTION DECLARATIONS                                                               
 ***************************************************************************/

static void DOT_DrawXAxis(SegmenT seg2, RecT  r, Int4 height, Int4 xstart,Int4 xstop, Int4 scale);
static void DOT_DrawYAxis(SegmenT seg2, RecT  r, Int4 height, Int4 ystart, Int4 ystop, Int4 scale, Int4 Fh);
static void DOT_StartDotPlotWithParams (DOTVibDataPtr vdp, BioseqPtr qbsp, BioseqPtr sbsp);
static Boolean DOT_StartDotPlot (BioseqPtr qbsp, BioseqPtr sbsp);
static Int4 DOT_GetValue (TexT t);
static void DOT_VCrossHairs (RecT rcP, DOTVibDataPtr vdp, Int4 VFrom,Int4 HFrom);
static void DOT_HCrossHairs (RecT rcP, DOTVibDataPtr vdp, Int4 VFrom,Int4 HFrom);
static void DOT_QViewerClickProc(VieweR v, SegmenT seg, PoinT pt);
static void DOT_SViewerClickProc(VieweR v, SegmenT seg, PoinT pt);
static void DOT_ExitAlign(DOTAlignInfoPtr alp);

/*_________________________________________(DOT_Compression)_____

  Purpose : Calculate sequence compression for window1

____________________________________________________________________*/

Int4 DOT_Compression (Int4  len, Int4 viewersize)
{
  Int4   cmpfactor;
  double f;

  

      f = (double)(len/viewersize);
      if (f <= 1)
        cmpfactor = 1;
      else
        cmpfactor = ceil(f);

      return cmpfactor;

}

/*_______________________________(DOT_SetupMenus)____________

  Purpose : Setup menus function for window1.

____________________________________________________________________*/
void DOT_SetupMenus ()

{
  MenU    m;  

#ifdef WIN_MAC
  m = AppleMenu (NULL);
  DeskAccGroup (m);
#endif

}



/*________________________________(DOT_DrawXGrids)____________________


  Purpose : Draw x-axis for DisplayHits, window1.

____________________________________________________________________*/
void DOT_DrawXGrids (RecT rupdate, RecT rcP, DOTVibDataPtr vdp, Int4 VFrom, Int4 HFrom, Int4 HTo, Int4 comp, Boolean GRID)
{
  Int4         x, y, y2, offset, start, end;
  Int4         scale_pos, pos, Hseq_pos;
  Char         scale_buf[15] = {""};	/*scale value*/
  Boolean      Decrement = FALSE;
  
  
  if (vdp->mip->qstrand == Seq_strand_minus)
    Decrement=TRUE;

  if (vdp->curr_qlen <50)
    return;
 
  pos = 100;
  offset=vdp->mip->q_start;
 /* select the font type */
  SelectFont(vdp->Fnt);

  Black(); /* X axis */
  MoveTo (rcP.left,rcP.top-10);
  LineTo (rcP.left+ HTo,rcP.top-10);


  /* write sequence length */
  if ( vdp->curr_qlen - HFrom == HTo)
    {

      Red();
      sprintf(scale_buf, "%ld", (long)vdp->mip->q_stop);
      MoveTo (rcP.left+ HTo +10 ,rcP.top-5);
      PaintString (scale_buf);
    }

  HTo += HFrom;
  
/*   if (!Decrement) */
 /*    { */
      for (scale_pos = HFrom; scale_pos <= HTo; scale_pos++)
        {
          
          /*  draw 0 on axis */
          if (scale_pos == 0) 
            { 
              Black();
              x =  rcP.left;
              y =  rcP.top;
              MoveTo (x, y-10);
              LineTo (x,y-20);
              Blue();
              MoveTo(x,y -25);
              sprintf(scale_buf, "%ld", (long)offset);
              PaintString (scale_buf);
              
              
            } 
          else  
            {
              if (!(scale_pos % pos))
                {
                  Hseq_pos = (!Decrement)?((scale_pos)*comp)+offset:offset - ((scale_pos)*comp);
                  
                  x = rcP.left + scale_pos-HFrom;
                  y = MAX(rupdate.top, rcP.top);
                  y2 = MIN(rupdate.bottom, rcP.top+vdp->curr_slen-VFrom);
                  
                  if ((!(HTo<rcP.left))&& (y<y2) && GRID)
                    {
                      LtGray();
                      Dotted();
                      MoveTo (x, y);
                      LineTo (x,y2);  
                    }
                  
                  /* add scale values */
                  Black();
                  Solid();
                  y = rcP.top-10;
                  MoveTo (x, y);
                  LineTo (x,y-10);
                  sprintf(scale_buf, "%ld", (long)Hseq_pos);
                  x = rcP.left + scale_pos -HFrom - (StringWidth(scale_buf)/2);
                  y = rcP.top -25;
                  if (x>rcP.left)
                    {
                      Blue();
                      MoveTo(x,y);
                      PaintString (scale_buf);
                    }
                }
              else if (!(scale_pos % (pos/2)))
                {
                  x = rcP.left + scale_pos -HFrom;
                  y = rcP.top-10;
                  MoveTo (x, y);
                  LineTo (x,y-7);
                  
                }
              else if (!(scale_pos % (pos/10)))
                {
                  Black();              
                  x = rcP.left + scale_pos -HFrom;
                  y = rcP.top-10;
                  MoveTo (x, y);
                  LineTo (x,y-5);
                  
                }
            }
        }


 end:
  Black();
    
    
}

/*________________________________(DOT_DrawYGrids)____________________

  Purpose : Draw y-axis for DisplayHits, window1.

____________________________________________________________________*/
static void DOT_DrawYGrids (RecT rupdate, RecT rcP, DOTVibDataPtr vdp, Int4 VFrom, Int4 HFrom, Int4 VTo, Int4 comp, Boolean GRID)
{
  Int4         x, y, x2, offset;
  Int4         scale_pos, pos, Vseq_pos;
  Char         scale_buf[15] = {""};	/*scale value*/
  Int4         fh;
  Boolean      Decrement = FALSE;
  

  if (vdp->mip->sstrand == Seq_strand_minus)
    Decrement = TRUE;

  if (vdp->curr_slen <50)
    return;

  pos = 100;
  offset=vdp->mip->s_start;

  SelectFont(vdp->Fnt);
  fh = FontHeight();

  Black();
  MoveTo (rcP.left-10,rcP.top);
  LineTo (rcP.left-10,rcP.top+ VTo);


  if (vdp->curr_slen - VFrom == VTo)
    {
  
      /* write the sequence length */
      Red();
      sprintf(scale_buf, "%ld", (long)vdp->mip->s_stop);
      MoveTo (rcP.left-10 -(StringWidth(scale_buf)/2),rcP.top+ VTo +10 + (fh/2)); 
      PaintString (scale_buf);
    }

  VTo += VFrom;

  for (scale_pos = VFrom; scale_pos <= VTo; scale_pos++)
    {
      if (scale_pos == 0) 
        {
          Black();
          x =  rcP.left-10;
          y =  rcP.top;
          MoveTo (x, y);
          LineTo (x-10,y);
          Blue();
          sprintf(scale_buf, "%ld", (long)offset);
          MoveTo(rcP.left-25-StringWidth(scale_buf),y+fh-4);
          PaintString (scale_buf);
        } 
      else if (!(scale_pos % pos))
        {
          Vseq_pos = (!Decrement)?((scale_pos)*comp)+offset: offset - ((scale_pos)*comp);
          
          x = MAX(rupdate.left, rcP.left);
          x2 = MIN(rupdate.right, rcP.left+vdp->curr_qlen-HFrom);
          y = rcP.top -VFrom + scale_pos;
         
          /* draw vertical grid */
          if ((scale_pos != 0) && (x<x2) && GRID)
            {
              LtGray();
              Dotted(); 
              MoveTo (x, y);
              LineTo (x2,y);
            }
          /* add scale values */
          Black(); 
          Solid();
          x = rcP.left -10;
          MoveTo (x, y);
          LineTo (x-10,y);
          sprintf(scale_buf, "%ld", (long)Vseq_pos);
          x = rcP.left -25 - StringWidth(scale_buf);
          y = rcP.top-VFrom + scale_pos+ (fh/2);
           
          if (y-fh > rcP.top)
            {
                Blue();
                MoveTo(x,y);
                PaintString (scale_buf);
            }
          
        }
      else if (!(scale_pos % (pos/2)))
        {
          Black(); 
          x = rcP.left - 10;
          y = rcP.top - VFrom + scale_pos;
          MoveTo (x, y);
          LineTo (x-7,y);
        }
      else if (!(scale_pos % (pos/10)))
        {
          Black();
          x = rcP.left - 10;
          y = rcP.top - VFrom + scale_pos;
          MoveTo (x, y);
          LineTo (x-5,y);
        }
    }
  Black();
}


/*_______________________________________________(DOT_UpdateLRBT)________________

  Purpose : Computes Left, Right, Bottom and Top values for DisplayHits.

____________________________________________________________________*/
void DOT_UpdateLRBT (RecT r, RecT rcP, Int4Ptr Left, Int4Ptr Right, Int4Ptr Bottom, Int4Ptr Top)
{
  
  if (r.left > rcP.left)
    *Left = r.left;
  else
    *Left = rcP.left;
  
  if (r.right < rcP.right)
    *Right = r.right;
  else
    *Right = rcP.right;
  
  if (r.bottom < rcP.bottom)
    *Bottom = r.bottom;
  else
    *Bottom = rcP.bottom;
  
  if (r.top > rcP.top)
    *Top = r.top;
  else
    *Top = rcP.top;
  
  return;
}
/*________________________________________(DOT_AddRectMargins)_____________

  Purpose : Add horizontal and vertical margins to rect, window1.

____________________________________________________________________*/
static void DOT_AddRectMargins (RectPtr r, DOTVibDataPtr vdp)
{
  r->left += vdp->HORZ_MARGIN;
  r->top += vdp->VERT_MARGIN;

  return;
}

/*____________________________________________(DOT_ChangeMainViewerCutoff)______

  Purpose : Change cutoff function for threshold-ramp, window1.

____________________________________________________________________*/
void DOT_ChangeMainViewerCutoff (BaR b, GraphiC g, Int2 new, Int2 old) 
{
  DOTVibDataPtr vdp;
  WindoW     w, temport;
  RecT       rcP;

  vdp = (DOTVibDataPtr) GetObjectExtra (b);

  Select(vdp->panel);
  ObjectRect(vdp->panel, &rcP);
  InsetRect(&rcP,4,4);

  w = (WindoW)ParentWindow(vdp->panel);
  temport = SavePort(w);

  vdp->sdp.TrampPos = new+20;
  
  DOT_AddRectMargins(&rcP, vdp);
  InvalRect (&rcP);
  RestorePort (temport);
  Update();

}

static PoinT curpnt;
static PoinT fstpnt;

/*________________________________________(DOT_RectOverlpRect)_____________

  Purpose : Find overlapping region between two rects.

____________________________________________________________________*/
static Boolean DOT_RectOverlpRect (RectPtr r1, RectPtr r2)
{
  if (r1->left <r2->left)
    r1->left = r2->left;
  if (r1->right >r2->right)
    r1->right = r2->right;
  if (r1->top < r2->top)
    r1->top = r2->top;
  if (r1->bottom > r2->bottom)
    r1->bottom = r2->bottom;

  return TRUE;
}

/*________________________________________(DOT_RectIntsRect)_____________

  Purpose : Checks if selected rect is in selectable rect, window1.

____________________________________________________________________*/
static Boolean DOT_RectIntsRect (RectPtr r1, RectPtr r2)
{
  
  if (r1->top > r2->bottom || r1->bottom < r2->top || r1->left > r2->right || r1->right < r2->left)
    return FALSE;
  else
    return TRUE;
}

/*________________________________________(DOT_LoopStop)_____________

  Purpose : calculate stop value in score_array for threshold-ramp.

____________________________________________________________________*/
static Int4 DOT_LoopStop (DOTVibDataPtr vdp)
{
  Int4   stop, pos;

  pos=ceil((double)(vdp->sdp.TrampPos*vdp->mip->unique)/100);

  if (pos>=vdp->mip->unique)
    return vdp->mip->index;
  else
    stop=vdp->mip->score_array[pos]+1;

  return stop;
}


/*_______________________________________________(DOT_DisplayHits)_________

Purpose : Draw function for window1.

____________________________________________________________________*/

void DOT_DisplayHits (PaneL p)
{

  Int4           i, x, y, x2, y2, index;
  RecT           rcP, rupd, rcS, rcs, dr, rcR, rcP_off;
  Int4           Xscore, Yscore, comp, stop;
  Uint1Ptr PNTR  scores;
  DOTVibDataPtr    vdp;
  DOTMainDataPtr   mip;
  DOTSelDataPtr      data;
  WindoW         w;
  BaR            Vsb, Hsb;
  Int4           VFrom, HFrom, VTo, HTo, cutoff;
  Int4           main_left, main_top, main_bottom, main_right;
  Int4           q_start, s_start, length;
  Int4           s_stop, q_stop, ycomp, xcomp;
  Int4           swdt, shgt;
  DOTDiagPtr     PNTR hitlist;
  DOTAlnPtr      PNTR Alnlist;
  Int4           Left, Right, Bottom, Top;
  Int2           dx, dy;
  DOTAlnPtr      PNTR alnL;
  DOTAlnPtr      aln;
  
  w = (WindoW)ParentWindow(p);
  rupd = updateRect;
  ObjectRect(p, &rcP);
  InsetRect(&rcP,4,4); 
  ClipRect(&rcP);
  
  if (!(vdp = (DOTVibDataPtr)GetObjectExtra(w))) return;
  
  VFrom = vdp->sdp.VFrom; 
  HFrom = vdp->sdp.HFrom;
  VTo = (MIN (vdp->sdp.PgLen+VFrom, vdp->curr_slen-VFrom));
  HTo = (MIN (vdp->sdp.PgWdth+HFrom, vdp->curr_qlen-HFrom));
  
  comp = vdp->comp;
  

  DOT_AddRectMargins (&rcP, vdp);
  DOT_UpdateLRBT (rupd, rcP, &Left, &Right, &Bottom, &Top);
  
  mip = vdp->mip;
  data = (DOTSelDataPtr)vdp->data;
  if (vdp->selectMode == dot_FEATVIEW)
    {
      if(data->selected)
        {
          DOT_VCrossHairs (rcP, vdp, VFrom, HFrom);
          DOT_HCrossHairs (rcP, vdp, VFrom, HFrom);
        }
    }
  else
    {
  if (data->selected) 
    {
      rcs=data->rcS;
      
      dx=HFrom-data->H_pos;
      dy=VFrom-data->V_pos;
      
      xcomp=vdp->originalcomp-comp;
      ycomp=xcomp;
      
      rcS.left=rcs.left-HFrom;
      rcS.right=rcs.right-HFrom;
      rcS.top=rcs.top-VFrom;
      rcS.bottom=rcs.bottom-VFrom;
      
      rcR.left = Left;
      rcR.right = Right;      
      rcR.top = Top;
      rcR.bottom = Bottom;
      
      if (DOT_RectIntsRect(&rcS, &rcR))
        {
          rcP_off = rcP;
          DOT_RectOverlpRect(&rcS, &rcR);
          SectRect (&rcS, &rcP_off, &dr);
          Yellow();
          PaintRect (&dr);
          Black();
        }
    }
    }
  
  
  DOT_DrawXGrids(rupd, rcP, vdp, VFrom, HFrom, HTo, comp, vdp->showGrid);
  DOT_DrawYGrids(rupd, rcP, vdp, VFrom, HFrom, VTo, comp, vdp->showGrid);
  
  if (vdp->showDotPlot)
    {
      hitlist = mip->hitlist;
      if (mip->unique<=1)
        stop = mip->index;
      else
        stop = DOT_LoopStop (vdp);
      
      for (i = 0; i<stop ; i++)
        {       
          length = hitlist[i]->length;
          q_start = ABS(mip->q_start-hitlist[i]->q_start);
          s_start = ABS(mip->s_start-hitlist[i]->s_start);
            
          q_start = (q_start/comp)-HFrom;
          s_start = (s_start/comp)-VFrom;
          length = length/comp;
          
          x = rcP.left + q_start;
          x2 = x+ length;

          y = rcP.top +  s_start;
          y2 = y + length;
         
          if (y > Bottom || y2 < Top || x > Right || x2 < Left)
            continue; /* outside of drawing Rgn */
          
          if (y < rcP.top) 
            {
              x = x+(rcP.top-y);
              y=rcP.top;
            }
          
          if (x< rcP.left)
            {
              y = y+(rcP.left-x); 
              x = rcP.left;
            }
          
          MoveTo(x, y);
          LineTo(x2, y2);
          
        }
    }

  if (vdp->showALIGN) /* overlay Blast hits */
    {
      Red();
      alnL = vdp->alp->Alnlist;
      stop = vdp->alp->index;
      main_left = rcP.left + mip->q_start;
      main_right = main_left + mip->qlen;
      main_top = rcP.top + mip->s_start;
      main_bottom = main_top + mip->slen;
      
      for (i = 0; i<stop; i++)
        {

          /* you should get aln->stops instead of using the lengths*/
          aln = alnL[i];
          length = aln->q_stop - aln->q_start;
          q_start = (aln->q_start/comp)-HFrom;
          s_start = (aln->s_start/comp)-VFrom;
          length = length/comp;
          
          x = rcP.left +  q_start;
          y = rcP.top +  s_start;
          
          x2 = x+ length-1;
          y2 = y+ length-1;
          
          if (y > main_bottom || y2 < main_top || x > main_right || x2 < main_left)
            continue; /* outside of drawing Rgn */
          
          if (y < rcP.top) 
            {
              x = x+(rcP.top-y);
              y=rcP.top;
            }
          
          if (x< rcP.left)
            {
              y = y+(rcP.left-x); 
              x = rcP.left;
            }
          
          if (x2>main_right)
            {
              y2=y2-(x2-main_right);
              x2=main_right;
            }
          
          if (y2>main_bottom)
            {
              x2 = x2-(y2-main_bottom);
              y2=main_bottom;
            }
          
          MoveTo(x, y);
          LineTo(x2, y2);
          
        }
      Black();
    }
  
  ResetClip();
  return;
}


/*________________________________________(DOT_SetCurrSeqlen)_____________


  Purpose : Update displayed seq-length after reduce/enlarge functions.

____________________________________________________________________*/
  /*set the current size of the sequence display */
static void DOT_SetCurrSeqlen (DOTVibDataPtr vdp)
{
  Int4 comp;
 
  comp = vdp->comp;

  vdp->curr_slen = (vdp->mip->slen)/comp;
  vdp->curr_qlen = (vdp->mip->qlen)/comp;
 
}


/*________________________________________(DOT_VScrlUpdate)_____________

  Purpose : Update function for vertical scroll proc, window1.

____________________________________________________________________*/

static void DOT_VScrlUpdate(DOTVibDataPtr vdp, BaR vsb, Int4 VCurPos)
{

  VCurPos = VCurPos*vdp->sdp.UnitY;

  /*set cursor position to new Units */
  VCurPos = VCurPos/vdp->sdp.UnitY;

  if (VCurPos<0) VCurPos=0;
  
  if (VCurPos >= vdp->sdp.YScrlMax)
    vdp->sdp.YScrlPos = vdp->sdp.YScrlMax;
  else 
    vdp->sdp.YScrlPos = VCurPos;
  
  vdp->sdp.VFrom=vdp->sdp.YScrlPos*vdp->sdp.UnitY;

  /*update scroll*/
  CorrectBarMax(vsb, vdp->sdp.YScrlMax);
  CorrectBarValue(vsb, vdp->sdp.YScrlPos);
  CorrectBarPage(vsb, vdp->sdp.YScrlPage, vdp->sdp.YScrlPage);
  
}


/*________________________________________(DOT_HScrlUpdate)_____________

  Purpose : Update function for horizontal scroll proc.

____________________________________________________________________*/

static void DOT_HScrlUpdate(DOTVibDataPtr vdp, BaR hsb, Int4 HCurPos)
{  
  HCurPos = HCurPos*vdp->sdp.UnitX;

  /*set cursor position to new Units */
  HCurPos = HCurPos/vdp->sdp.UnitX;

  if (HCurPos<0) HCurPos=0;
  
  if (HCurPos >= vdp->sdp.XScrlMax)
    vdp->sdp.XScrlPos = vdp->sdp.XScrlMax;
  else 
    vdp->sdp.XScrlPos = HCurPos;
  
  vdp->sdp.HFrom=vdp->sdp.XScrlPos*vdp->sdp.UnitX;

  /*update scroll*/
  CorrectBarMax(hsb, vdp->sdp.XScrlMax);
  CorrectBarValue(hsb, vdp->sdp.XScrlPos);
  CorrectBarPage(hsb, vdp->sdp.XScrlPage, vdp->sdp.XScrlPage);
  
}



/*________________________________________________(DOT_VscrlProc)_____________

  Purpose : Vertical Scroll proc, window1.

____________________________________________________________________*/
static void DOT_VscrlProc (BaR vsb, SlatE s, Int2 newval, Int2 oldval)
{
  WindoW 		temport, w;
  DOTVibDataPtr   vdp;
  RecT         rcP;
  Int2         dy, offset;
  Int2         visLines, vmargin, hmargin;
  PaneL        p;


  w = (WindoW)ParentWindow((PaneL)s);
  vdp = (DOTVibDataPtr)GetObjectExtra (w);
  p = vdp->panel;
  
  offset=vdp->sdp.HFrom;
  vmargin=vdp->VERT_MARGIN;
  hmargin=vdp->HORZ_MARGIN;

  temport = SavePort (w);

  Select(p);
  ObjectRect(p, &rcP);
  InsetRect(&rcP, 4, 4);
  ClipRect(&rcP);

  vdp->sdp.YScrlPos = newval;
  visLines = vdp->sdp.YScrlPage;
  vdp->sdp.VFrom=newval*vdp->sdp.UnitY;

  rcP.right = (MIN (rcP.left+hmargin+vdp->curr_qlen-offset, rcP.right));
  rcP.top += vmargin;

  dy = newval- oldval;
  if (ABS(dy) < vdp->sdp.YScrlPage)
    {
      ScrollRect(&rcP,  0, (Int2)((-dy)*vdp->sdp.UnitY));
    }
  else
    {
      InsetRect(&rcP, -1, -1);
      InvalRect(&rcP);
    }
  ResetClip();
  RestorePort(temport);
/*   Update(); */
  
}



/*________________________________________________(DOT_HscrlProc)_____________

  Purpose : Horizontal scroll proc, window1.

____________________________________________________________________*/
static void DOT_HscrlProc (BaR Hsb, SlatE s, Int2 newval, Int2 oldval)
{
  WindoW 		temport, w;
  DOTVibDataPtr   vdp;
  RecT         rcP;
  Int2         dx, visLines, offset, hmargin, vmargin;

  
  w = (WindoW)ParentWindow((PaneL)s);
  vdp = (DOTVibDataPtr)GetObjectExtra (w);

  temport = SavePort (w);
  Select(vdp->panel);
  ObjectRect(vdp->panel, &rcP);
  InsetRect(&rcP, 4, 4);
  ClipRect (&rcP);

  offset=vdp->sdp.VFrom;
  hmargin=vdp->HORZ_MARGIN;
  vmargin=vdp->VERT_MARGIN;

  rcP.bottom = (MIN (rcP.top+vmargin+vdp->curr_slen-offset, rcP.bottom));
  rcP.left += hmargin;

  vdp->sdp.XScrlPos = newval;
  visLines = vdp->sdp.XScrlPage;
  vdp->sdp.HFrom=newval*vdp->sdp.UnitX;

  dx = newval - oldval;

  if (ABS(dx) < vdp->sdp.XScrlPage)
    {
      ScrollRect(&rcP, (Int2)((-dx)*vdp->sdp.UnitX) , 0);
    }
  else
    {
      InsetRect(&rcP, -1, -1);
      InvalRect(&rcP);
    }
  
  ResetClip();
  RestorePort(temport);
/*   Update(); */
  
}


/*________________________________________(DOT_UpdateMainPanel)_____________

  Purpose : Update selection for window1.

____________________________________________________________________*/

static void DOT_UpdateMainPanel(DOTVibDataPtr vdp, Boolean update_all)
{
  WindoW     temport;
  RecT       rc;
  DOTSelDataPtr  data;

  data=(DOTSelDataPtr)vdp->data;

  temport = SavePort((WindoW)ParentWindow(vdp->panel));
  Select(vdp->panel);
  ObjectRect(vdp->panel, &rc);
  ClipRect(&rc);
  
  if (!update_all)
    {
      if (data->rm_lastselected)
        { 
          InsetRect(&data->old_rcS, -1, -1);
          InvalRect(&data->old_rcS);
          data->rm_lastselected=FALSE;
        }
      InsetRect(&data->rcS, -1, -1);
      InvalRect (&data->rcS);
    }
  else
    {
      InsetRect(&rc, -1, -1);
      InvalRect(&rc);
    }

  ResetClip();
  RestorePort(temport);
  Update();

}


/*______________________________________(DOT_VScroll)___________

  Purpose : Correct vertical scroll bar values, window1.

____________________________________________________________________*/
static void DOT_VScroll (DOTVibDataPtr vdp, BaR vsb)
{
 
  CorrectBarMax(vsb,vdp->sdp.YScrlMax);
  CorrectBarValue(vsb,vdp->sdp.YScrlPos);
  CorrectBarPage(vsb, vdp->sdp.YScrlPage, vdp->sdp.YScrlPage);

}
/*______________________________________(DOT_HScroll)___________

  Purpose : Correct horizontal scroll bar values, window1.

____________________________________________________________________*/
static void DOT_HScroll (DOTVibDataPtr vdp, BaR hsb)
{

  CorrectBarMax(hsb,vdp->sdp.XScrlMax);
  CorrectBarValue(hsb,vdp->sdp.XScrlPos);
  CorrectBarPage(hsb,vdp->sdp.XScrlPage, vdp->sdp.XScrlPage);

}
/*________________________________________(DOT_ComputePanelSize)_____________

  Purpose : Calculate panel size for scrolling functions, window1.

____________________________________________________________________*/
static void DOT_ComputePanelSize (RecT rcP, DOTVibDataPtr vdp, Int4Ptr PgWdth, Int4Ptr PgLen)
{

  InsetRect(&rcP, 4, 4);
  rcP.left += vdp->HORZ_MARGIN;
  rcP.top += vdp->VERT_MARGIN;

  *PgWdth =rcP.right-rcP.left;
  *PgLen  =rcP.bottom-rcP.top;

}

/*________________________________________(DOT_SetScrlVals)_____________

  Purpose : Set scroll values for window1.

____________________________________________________________________*/
static void DOT_SetScrlVals (DOTVibDataPtr vdp)
{
/*   Int4  scrollfctr = 40; */

/*   vdp->sdp.UnitY = vdp->sdp.PgLen/scrollfctr; */
  vdp->sdp.UnitY = 16; /*  constant value*/
  vdp->sdp.TotUnitsY = vdp->curr_slen/vdp->sdp.UnitY;
  vdp->sdp.YScrlPage = vdp->sdp.PgLen/vdp->sdp.UnitY;
  vdp->sdp.YScrlMax = vdp->sdp.TotUnitsY-(vdp->sdp.YScrlPage - vdp->VERT_MARGIN/vdp->sdp.UnitY);

/*   vdp->sdp.UnitX = vdp->sdp.PgWdth/scrollfctr; */
  vdp->sdp.UnitX = 16;/*  constant value*/
  vdp->sdp.TotUnitsX = vdp->curr_qlen/vdp->sdp.UnitY;
  vdp->sdp.XScrlPage = vdp->sdp.PgWdth/vdp->sdp.UnitX;
  vdp->sdp.XScrlMax = vdp->sdp.TotUnitsX-(vdp->sdp.XScrlPage - vdp->HORZ_MARGIN/vdp->sdp.UnitX);

  /* image is smaller than page size */
  if ((vdp->sdp.YScrlPage + vdp->HORZ_MARGIN/vdp->sdp.UnitX) > vdp->sdp.TotUnitsY)
    {
      vdp->sdp.YScrlMax = 0;
      vdp->sdp.YScrlPage = 0;
      vdp->sdp.YScrlPos = 0;
    }
  if (vdp->sdp.XScrlPage>vdp->sdp.TotUnitsX)
    {
      vdp->sdp.XScrlMax = 0;
      vdp->sdp.XScrlPage = 0;
      vdp->sdp.XScrlPos = 0;
    }
 
}

/*______________________________________(DOT_SetUpWin)___________

  Purpose : Scrolling info setup function for window1.

____________________________________________________________________*/
static void DOT_SetUpWin(WindoW w, PaneL p, DOTVibDataPtr vdp)
{

  Int4       Vcurpos, Hcurpos, height, width, viewersize, len;
  Int4       gap, lmargin, vsbWidth, hsbHeight;
  BaR        vsb, hsb;
  WindoW     temport;
  RecT       rcP, rcW, rcHsb, rcVsb;



  temport = SavePort (w);
  Select(p);
  ObjectRect(p, &rcP);

 /* Reset Panel Parameters */
  ObjectRect(w, &rcW);
  width = rcW.right-rcW.left;
  height = rcW.bottom-rcW.top;
  vsb = GetSlateVScrollBar ((SlatE) p);
  hsb = GetSlateHScrollBar ((SlatE) p);
  
  GetPosition(vsb,&rcVsb);
  GetPosition(hsb,&rcHsb);

  gap=2;
  lmargin=10;
  vsbWidth=rcVsb.right-rcVsb.left;
  hsbHeight=rcHsb.bottom-rcHsb.top;
    
  rcP.right = width - vsbWidth - gap;
  rcP.bottom = height - hsbHeight - gap;
  rcP.left=lmargin;

  SetPosition (vdp->panel, &rcP);
  AdjustPrnt (vdp->panel, &rcP, FALSE);

  viewersize=MIN(rcP.right-rcP.left,rcP.bottom-rcP.top)-vdp->HORZ_MARGIN;
  len=MAX(vdp->mip->qlen, vdp->mip->slen);
  vdp->comp=DOT_Compression(len, viewersize);
  vdp->originalcomp=vdp->comp;
 

  vdp->sdp.YScrlPos = 0;
  vdp->sdp.XScrlPos = 0;

  DOT_SetCurrSeqlen (vdp);
  DOT_ComputePanelSize(rcP, vdp, &(vdp->sdp.PgWdth), &(vdp->sdp.PgLen));
  DOT_SetScrlVals(vdp);

  DOT_VScroll (vdp, vsb);
  DOT_HScroll (vdp, hsb);

  RestorePort(temport);

}




/*________________________________________(DOT_CloseSequenceWindow)_____________

  Purpose : Close function for Sequence Window.

____________________________________________________________________*/
static void  DOT_CloseSequenceWindow (ButtoN b)
{

  DOTVibDataPtr vdp2, vdp;
  DOTSelDataPtr  data;

  Uint1Ptr   qseq, sseq;

  vdp2=(DOTVibDataPtr)GetObjectExtra(ParentWindow(b));

  data=(DOTSelDataPtr)vdp2->data;
  data->selected=FALSE;
  vdp=data->vdp;
  SetTitle(vdp->Infopanel, vdp->iInfo);

  DOT_UpdateMainPanel(vdp, FALSE);
 
  vdp2->sv->pict1=DeletePicture(vdp2->sv->pict1);
  vdp2->sv->pict2=DeletePicture(vdp2->sv->pict2);
  
  if (vdp2->sv->salp)
    MemFree(vdp2->sv->salp);
  if (vdp2->sv) MemFree(vdp2->sv);

  DOT_FreeMainInfo(vdp2->mip);
  if (vdp2->mip) MemFree(vdp2->mip);
  if (vdp2) MemFree(vdp2);
  vdp->ChildWin=Remove (vdp->ChildWin);
}



/*_____________________________________________________________________

  Purpose : Remove feature linked list

____________________________________________________________________*/
static void DOT_FreeFeatPointers(DOTRowPtr drp)
{
  DOTFeatPtr dfp_temp=NULL, dfp=NULL;

  dfp=drp->dfp;
  dfp_temp=dfp;
  dfp=dfp->next;
  while(dfp != NULL)
    {
      dfp_temp=MemFree(dfp_temp);
      dfp_temp = dfp;
      dfp = dfp->next;
    }
  if (dfp_temp!=NULL) MemFree(dfp_temp);

}

/*________________________________________(DOT_CloseFeatWindow)_____________

  Purpose : Close function for Feature Window.

____________________________________________________________________*/
static void  DOT_CloseFeatWindow (IteM i)
{
  DOTFeatListPtr flp;
  DOTSelDataPtr  data;


      flp=(DOTFeatListPtr)GetObjectExtra(ParentWindow(i));
      if (!flp) return;

      data=(DOTSelDataPtr)flp->data;
      data->selected=FALSE;
      DOT_UpdateMainPanel(data->vdp, TRUE);

      flp->segQuery=DeletePicture(flp->segQuery);
      flp->segSubject=DeletePicture(flp->segSubject);
      DOT_FreeFeatPointers(flp->query_drp);
      DOT_FreeFeatPointers(flp->subject_drp);
      if (flp->featindex) flp->featindex=MemFree(flp->featindex);
      if (flp->query_drp) flp->query_drp=MemFree(flp->query_drp);
      if (flp->subject_drp)flp->subject_drp=MemFree(flp->subject_drp);
      if (flp->FeatWin) flp->FeatWin=Remove (flp->FeatWin);
      data->vdp->ChildWin=NULL;
      if (flp) flp=MemFree(flp);
     
}

/*_____________________________________________________________________

  Purpose : Remove 'sequence' or 'feature' window

____________________________________________________________________*/

static WindoW DOT_ClearLastWindow(WindoW w, Boolean is_sequence)
{
  DOTVibDataPtr vdp;
  DOTFeatListPtr flp;
  DOTSelDataPtr  data;

  if (!w) return w;

  if (is_sequence)
    {
      vdp=(DOTVibDataPtr)GetObjectExtra(w);
      if (!vdp) return w;

      data=(DOTSelDataPtr)vdp->data;
      data->selected=FALSE;

      if(vdp->sv->v1) DeleteViewer(vdp->sv->v1);
      if(vdp->sv->v2) DeleteViewer(vdp->sv->v2);
      if(vdp->sv->pict1) DeletePicture(vdp->sv->pict1);
      if(vdp->sv->pict2) DeletePicture(vdp->sv->pict2);
      vdp->sv=MemFree(vdp->sv);

/*       if (vdp->mip)  */
/*         { */
/*           DOT_FreeMainInfo(vdp->mip); */
/*           MemFree(vdp->mip); */
/*         } */

    }
  else
    {
      flp=(DOTFeatListPtr)GetObjectExtra(w);
      if (!flp) return w;

      data=(DOTSelDataPtr)flp->data;
      data->selected=FALSE;

      DeletePicture(flp->segQuery);
      DeletePicture(flp->segSubject);
     /*  DOT_FreeFeatPointers(flp->query_drp); */
/*       DOT_FreeFeatPointers(flp->subject_drp); */
/*       MemFree(flp->query_drp); */
/*       MemFree(flp->subject_drp); */
      MemFree(flp);
    }

  DOT_UpdateMainPanel(data->vdp, TRUE);

  w=Remove(w);

  return NULL;
 
}
/*_____________________________________________________________________

  Purpose : New set of function to show features on dotviewer

____________________________________________________________________*/

void DOT_ModeProc(ChoicE i)
{
  WindoW      w, temport;
  RecT        rcP;
  DOTVibDataPtr  vdp;
  DOTSelDataPtr  data;
  

  w = (WindoW)ParentWindow(i);
  temport=SavePort(w);
  vdp = (DOTVibDataPtr)GetObjectExtra (w);
  if (!vdp) return;

  if (vdp->ChildWin !=NULL)
    {
      if (vdp->selectMode == 1)
        vdp->ChildWin=DOT_ClearLastWindow(vdp->ChildWin, TRUE);
      else
        vdp->ChildWin=DOT_ClearLastWindow(vdp->ChildWin, FALSE);
    }

  vdp->selectMode = GetValue(i);
  data=(DOTSelDataPtr)vdp->data;
  data->selected=FALSE;
  ObjectRect(vdp->panel, &rcP);
  Select(vdp->panel);
  InvalRect(&rcP);
  RestorePort(temport);
  Update();

}

static void DOT_VCrossHairs (RecT rcP, DOTVibDataPtr vdp, Int4 VFrom,Int4 HFrom)
{
  Int4   y, y2, x, x2, Vseq;
  Int4   cursor_size =15;
  Char   seq_pos[15];
/*   SelectPtr c_data; */

/*   c_data = cip->data; */

  /* vertical line */
  y = rcP.top-vdp->VERT_MARGIN+20;
  y2= rcP.bottom-1/* MIN(rcP.bottom - 1 , VTo) */;
  x2 = x = curpnt.x/* +HFrom */;
/*   x2 = x = rcP.left +c_h.x -cip->HORZ_MARGIN - cursor_size; */
  if (x>=rcP.left)
    {
      Magenta();
      MoveTo (x, y);
      LineTo (x2, y2);

    }
  Black();
  return;
  
}

static void DOT_HCrossHairs (RecT rcP, DOTVibDataPtr vdp, Int4 VFrom,Int4 HFrom)
{
  Int4  y, y2, x, x2;
  Int4  cursor_size = 15;
  Char  seq_pos[15];
/*   SelectPtr  c_data; */

/*   data = vdp->data; */

  /* horizontal line */
    y2 = y = curpnt.y /* + VFrom */;
/*   y2 = y = rcP.top + c_h.y  -cip->VERT_MARGIN - cursor_size; */
  x = rcP.left-vdp->HORZ_MARGIN+20;
  x2 = rcP.right-1/* MIN(rcP.right -1, HTo) */;
  if (y>=rcP.top)
    {
      Magenta ();
      MoveTo(x, y);
      LineTo (x2, y2);
    }
  Black();
  return;
  
}


static void DOT_MoveCrossHairs (RecT rcP, DOTVibDataPtr vdp)
{
  Int4   y, y2, x, x2;
  Int4   cursor_size = 15;
  Char   seq_pos[15];
  DOTSelDataPtr  data;

  SelectFont(vdp->Fnt);

  data = (DOTSelDataPtr) vdp->data;
  DOT_AddRectMargins(&rcP, vdp);

  /* vertical line */
  y = rcP.top-vdp->VERT_MARGIN+20;
  y2= rcP.bottom-1;
  x2 = x = curpnt.x;
  if (x>rcP.left)
    {
      MoveTo (x, y);
      LineTo (x2, y2);
    }

  /* horizontal line */
  y2 = y = curpnt.y;
  x = rcP.left-vdp->HORZ_MARGIN +20;
  x2 = rcP.right-1;
  if (y>rcP.top)
    {
      MoveTo(x, y);
      LineTo (x2, y2);

    }

  Black();


  return;
  
}


static void DOT_SelectLineProc (PaneL p)

{
  RecT        rcP;
  DOTSelDataPtr  data;
  DOTVibDataPtr  vdp;

  vdp = (DOTVibDataPtr) GetObjectExtra (ParentWindow(p));
  Dotted ();
  ObjectRect (p, &rcP);
  InsetRect (&rcP, 4, 4);
  DOT_MoveCrossHairs (rcP, vdp); 
}


/*________________________________________(DOT_SelectFrameProc)_____________

  Purpose : select frame for click and drag functions of window1.

____________________________________________________________________*/
static void DOT_SelectFrameProc (PaneL p)

{
  RecT  dr;
  RecT  or;
  RecT  r;

  Dotted ();
  ObjectRect (p, &or);
  InsetRect (&or, 2, 2);
  LoadRect (&r, fstpnt.x, fstpnt.y, curpnt.x, curpnt.y);
  SectRect (&r, &or, &dr);
  FrameRect (&dr);
  
}

/*________________________________________(DOT_DrawXAxis)_____________

  Purpose : Draw x-axis function for viewer1, window2.

____________________________________________________________________*/
void DOT_DrawXAxis(SegmenT seg2, RecT  r, Int4 height, Int4 xstart,Int4 xstop, Int4 scale)
{
  Int4         pos, xlen, x, y, scale_pos, i, j, bigtick, midtick, smalltick;
  Char         scale_buf[15] = {""};	/*scale value*/
  Boolean      Decrement=FALSE;


  if (xstart>xstop)
    Decrement=TRUE;

  if (scale==0)
    scale=1;

  pos=100*scale;
  bigtick=10*scale;
  midtick=7*scale;
  smalltick=5*scale;
 

  xlen=ABS(xstop-xstart);

  AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
  AddLine(seg2, r.left, height+10, r.left+xlen, height+10, FALSE,-1);
/*   AddLine(seg2, r.left, height-ylen, r.left+xlen, height-ylen, FALSE, -1); */
  sprintf(scale_buf, "%ld", (long)xstop);
  AddAttribute(seg2, COLOR_ATT, RED_COLOR, 0,0,0,0);
  AddLabel(seg2, r.left+xlen+10*scale, height+5, scale_buf, SMALL_TEXT, 0, MIDDLE_RIGHT, 0);
  
  AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
  if (!Decrement)
    {
      for (scale_pos = xstart, i=0; scale_pos <= xstop; scale_pos++, i++)
        {
          
          if (!(scale_pos % pos))
            {
              x = r.left + i;
              y = height+10;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x, y+bigtick, FALSE,-1);
              sprintf(scale_buf, "%ld", (long)scale_pos);
              AddAttribute(seg2, COLOR_ATT, BLUE_COLOR, 0,0,0,0);
              AddLabel(seg2, x, y+15*scale, scale_buf, SMALL_TEXT, 0, UPPER_CENTER, 0);
              
            }
          else if (!(scale_pos % (pos/2)))
            {
              x = r.left + i;
              y = height+10;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x, y+midtick, FALSE,-1);
              
            }
          else if (!(scale_pos % (pos/10)))
            {
              x = r.left + i;
              y = height+10;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x, y+smalltick, FALSE,-1);
            }
        }
    }
  else
    {
      for (scale_pos = xstop, i=0, j=xstart; scale_pos <= xstart; scale_pos++, i++, j--)
        {
          
          if (!(scale_pos % pos))
                {
                  x = r.left + i;
                  y = height+10;
                  AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
                  AddLine(seg2, x, y, x, y+bigtick, FALSE,-1);
                  sprintf(scale_buf, "%ld", (long)j);
                  AddAttribute(seg2, COLOR_ATT, BLUE_COLOR, 0,0,0,0);
                  AddLabel(seg2, x, y+15*scale, scale_buf, SMALL_TEXT, 0, UPPER_CENTER, 0);
                  
                }
              else if (!(scale_pos % (pos/2)))
                {
                  x = r.left + i;
                  y = height+10;
                  AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
                  AddLine(seg2, x, y, x, y+midtick, FALSE,-1);
                  
                }
              else if (!(scale_pos % (pos/10)))
                {
                  x = r.left + i;
                  y = height+10;
                  AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
                  AddLine(seg2, x, y, x, y+smalltick, FALSE,-1);
                }
        }
    }
  
}


/*________________________________________(DOT_DrawYAxis)_____________


  Purpose : Draw y-axis function for viewer1, window2.

____________________________________________________________________*/
void DOT_DrawYAxis(SegmenT seg2, RecT  r, Int4 height, Int4 ystart, Int4 ystop, Int4 scale, Int4 Fh)
{
  Int4         smalltick, midtick, bigtick;
  Int4         pos, ylen, x, y, scale_pos, i, j, Fh_2,Fh_4; 
  Char         scale_buf[15] = {""};	/*scale value*/
  Boolean      Decrement = FALSE;


  if (ystart>ystop)
    Decrement=TRUE;

  if (scale==0)
    scale=1;

  pos=100*scale;
  smalltick=5*scale;
  midtick=7*scale;
  bigtick=10*scale;
 
  Fh_2=Fh/2;
  Fh_4=Fh/4;
  

  ylen=ABS(ystop-ystart);

  AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
  AddLine(seg2, r.left-10, height, r.left-10, height-ylen, FALSE,-1);

  sprintf(scale_buf, "%ld", (long)ystop);
  AddAttribute(seg2, COLOR_ATT, RED_COLOR, 0,0,0,0);
  AddLabel(seg2, r.left-10, height-ylen-10*scale, scale_buf, SMALL_TEXT, 0, MIDDLE_CENTER, 0);

  
  if (!Decrement)
    {
      for (scale_pos = ystart, i=0; scale_pos <= ystop; scale_pos++, i++)
        {
          if (!(scale_pos % pos))
            {
              x = r.left-10;
              y = height-i/* +VFrom */;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x-bigtick, y, FALSE,-1);
              
              sprintf(scale_buf, "%ld", (long)scale_pos);
              y = y-Fh_2;
              AddAttribute(seg2, COLOR_ATT, BLUE_COLOR, 0,0,0,0);
              AddLabel(seg2, x-15*scale, y, scale_buf, SMALL_TEXT, 0, MIDDLE_LEFT, 0);
              
            }
          else if (!(scale_pos % (pos/2)))
            {
              x = r.left-10;
              y = height-i/* +VFrom */;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x-midtick, y, FALSE,-1);
              
            }
          else if (!(scale_pos % (pos/10)))
            {
              x = r.left-10;
              y = height-i/* +VFrom */;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x-smalltick, y, FALSE,-1);
            }
          
        }
    }
  else
    {

      for (scale_pos = ystop, i=0, j=ystart; scale_pos <= ystart; scale_pos++, i++, j--)
        {
          if (!(scale_pos % pos))
            {
              x = r.left-10;
              y = height-i/* +VFrom */;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x-bigtick, y, FALSE,-1);
              
              sprintf(scale_buf, "%ld", (long)j);
              y = y-Fh_2;
              AddAttribute(seg2, COLOR_ATT, BLUE_COLOR, 0,0,0,0);
              AddLabel(seg2, x-15*scale, y, scale_buf, SMALL_TEXT, 0, MIDDLE_LEFT, 0);
              
            }
          else if (!(scale_pos % (pos/2)))
            {
              x = r.left-10;
              y = height-i/* +VFrom */;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x-midtick, y, FALSE,-1);
              
            }
          else if (!(scale_pos % (pos/10)))
            {
              x = r.left-10;
              y = height-i/* +VFrom */;
              AddAttribute(seg2, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
              AddLine(seg2, x, y, x-smalltick, y, FALSE,-1);
            }
          
        }
    }


}


/*________________________________________(DOT_DisplayDiags)_____________


  Purpose : Draw function for viewer1, window2.

____________________________________________________________________*/

static Boolean LIBCALLBACK DOT_SVDisplayDiags(DOTVibDataPtr vdp2, DOTSelDataPtr data)
{
  RecT               rcP;
  DOTDiagPtr         PNTR hitlist;
  DOTVibDataPtr        vdp;
  VieweR             v;
  SegmenT            seg1, seg2, seg3;
  Char               buffer[50];
  Int4               stop, cutoff=0, q_start, s_start, length;
  Int4               x, y, x2, y2;
  Int4               i, s_stop, q_stop;
  Int4               p_VFrom, p_HFrom;
  Int4               width, height, primID;
  Int4               r_width, r_ht, temp, x_start, y_start, x_stop, y_stop;


  seg1=CreateSegment(vdp2->sv->pict1, 1, 0); /* diags */
  seg2=CreateSegment(vdp2->sv->pict1, 2, 0); /* axis */
  seg3=CreateSegment(vdp2->sv->pict1, 3, 0); /* diag coordinates */

  v=vdp2->sv->v1;
  vdp2->sv->seg1=seg1; /* diags */

  GetPosition(v, &rcP); 
  InsetRect(&rcP, 4, 4);
  vdp = data->vdp;

  width = rcP.right-rcP.left;
  height = rcP.bottom-rcP.top-2*(vdp2->VERT_MARGIN*vdp2->sv->scaleValue);

  p_VFrom=vdp->sdp.VFrom;
  p_HFrom=vdp->sdp.HFrom;

  /* Rect Parameters */
  rcP.left+=2*(vdp2->HORZ_MARGIN*vdp2->sv->scaleValue);

  hitlist = vdp2->mip->hitlist;
  if (vdp2->mip->unique<=1)
    stop=vdp2->mip->index;
  else
    stop=DOT_LoopStop (vdp2);
  
  DOT_DrawXAxis(seg2, rcP, height, data->q_start, data->q_stop, vdp2->sv->scaleValue);
  DOT_DrawYAxis(seg2, rcP, height, data->s_start, data->s_stop, vdp2->sv->scaleValue, vdp2->Fh);

  AddAttribute(seg1, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);


  for (i = 0; i<stop ; i++)
       {  
   
         length = hitlist[i]->length-1;
         x_start =  ABS(data->q_start- hitlist[i]->q_start);
         y_start =   ABS(data->s_start-  hitlist[i]->s_start);

         x = rcP.left + x_start;
         x2 = x+ length;
         y = height- y_start;
         y2 = y - length;

         AddLine(seg1, x, y, x2, y2, FALSE,(Uint2) i+1);
       }

  return TRUE;
}

/*________________________________________(DOT_WorldtoScreen)_____________

  Purpose : calculate sequence coordinates from screen coords.

____________________________________________________________________*/
static Int4 DOT_WorldtoScreen(Int4 wPos, Int4 chw_2)
{
  Int4  sPos;

  sPos=wPos*chw_2;

  return sPos;
}

/*________________________________________(DOT_DrawScale)_____________

  Purpose : Draw scale of aligned seqs, viewer2, window2.

____________________________________________________________________*/
void DOT_DrawScale (SegmenT sbSeg, DOTMainDataPtr mip, RecT  r, Int4 margin,Int4 res_cnt, Int4 bloc_cnt, Int4 Fh, Int4 q_pos, Int4 s_pos, Int4 chw_2, Int4 chw_4)
{
  Char      Buf[15]={""};
  Int4      x, x2,y1, y2, y3, y4, y5, y6, y7, y8;
  Int4      pos, col;


  col =DOT_WorldtoScreen(res_cnt+bloc_cnt, chw_2); 

  x = r.left+margin+col-chw_4;
  y1=7*Fh;
  y2=y1-Fh/2;
  y3=y2-Fh/2;
  AddAttribute(sbSeg, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
  AddLine(sbSeg, x, y2, x, y3, FALSE, 0);

  pos=(mip->qstrand==Seq_strand_plus)?res_cnt+q_pos:ABS(q_pos-res_cnt); 
  sprintf(Buf,"%d",(int)pos);
  x2=x-StringWidth(Buf)/2;

  AddAttribute(sbSeg, COLOR_ATT, BLUE_COLOR, 0,0,0,0);
  AddLabel(sbSeg, x2,y1, Buf, SMALL_TEXT, 0, UPPER_RIGHT, 0);

  y4=2*Fh;
  y5=y4+Fh+Fh/2;
  y6=y5+Fh/2;
  AddAttribute(sbSeg, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
  AddLine(sbSeg, x, y6, x, y5, FALSE, 0);

  pos=(mip->sstrand==Seq_strand_plus)?res_cnt+s_pos:ABS(s_pos-res_cnt);
  sprintf(Buf,"%d",(int)pos);
  x-=StringWidth(Buf)/2;

  AddAttribute(sbSeg, COLOR_ATT, BLUE_COLOR, 0,0,0,0);
  AddLabel(sbSeg, x2,y4, Buf, SMALL_TEXT, 0, UPPER_RIGHT, 0);

}

/*________________________________________(DOT_HlSeq)_____________

  Purpose : Draw line to show hit in alignment, viewer2, window2.

____________________________________________________________________*/
static void DOT_HlSeq(SegmenT sbSeg, RecT  rcP, Int4 margin, Int4 wPos, Int4 bloc_pos, Int4 hit_start, Int4 hit_stop, Int4 Fh, Int4 chw_2)
{
  Int4  xpos, len,  bloc_beg, begin, end;
  Int4  x1, x2, y1;


  bloc_beg=wPos-bloc_pos;
  begin=MAX(bloc_beg, hit_start);
  
  end=MIN(wPos,hit_stop);
  len=end-begin;
  xpos=DOT_WorldtoScreen(begin, chw_2);
 
  x1= rcP.left+margin+xpos;
  x2= x1+len*chw_2;
  y1 = 5*Fh;

  AddAttribute(sbSeg, COLOR_ATT, RED_COLOR, 0,0,0,0);
 AddLine(sbSeg, x1, y1, x2, y1, FALSE,0); 

}

/*________________________________________(DOT_HlMatch)_____________

  Purpose : Highlight matched residues, viewer2, window2.

____________________________________________________________________
*/
static void DOT_HlMatch(SegmenT sbSeg, RecT  rcP, Int4 margin, Int4 wPos, Int4 Fh, Int4 chw_2)
{
  Int4  xpos, x1, x2, y1, y2;

  xpos=DOT_WorldtoScreen(wPos, chw_2);

  x1= rcP.left+margin+xpos;
  x2= x1+chw_2;
  y1 = 6*Fh;
  y2 = 4*Fh;

  AddAttribute(sbSeg, COLOR_ATT, CYAN_COLOR, 0,0,0,0);
  AddRectangle (sbSeg, x1, y1, x2, y2, 0, TRUE, 0);
}


/*________________________________________(DOT_UpdatePt)_____________

  Purpose : corrects clicked points for click/drag/release procs of window1.

____________________________________________________________________*/

static void DOT_UpdatePt (void)
{
  Int4    temp;

  if (curpnt.x < fstpnt.x)
    {
      temp = curpnt.x;
      curpnt.x = fstpnt.x;
      fstpnt.x = temp;
    }

  if (curpnt.y < fstpnt.y)
    {
      temp = curpnt.y;
      curpnt.y = fstpnt.y;
      fstpnt.y = temp;
    }

  return;
  
}

/*________________________________________(DOT_FillNewSeqBufs)_____________

  Purpose : Fill seq buffers for window2.

____________________________________________________________________*/
static void DOT_FillNewSeqBufs (DOTVibDataPtr vdp, DOTVibDataPtr vdp2, Boolean is_zoom)
{

  Int4      qlen, slen, i;
  Uint1Ptr  qseq, sseq, q, s, qBuf, sBuf;
  DOTSelDataPtr data;
  DOTMainDataPtr mip1, mip2;
  

  data=(DOTSelDataPtr)vdp->data;
  mip1 =  vdp->mip;
  mip2 =  vdp2->mip;
  
  qlen=mip2->qlen;
  slen=mip2->slen;
  qseq = mip1->qseq;
  sseq = mip1->sseq;

  if (is_zoom)
    {
      mip2->qseq=MemFree(mip2->qseq);
      mip2->sseq=MemFree(mip2->sseq);
    }

  /* get sequence position relative to the vdp sequence buffer */

  q = qseq+ABS(data->q_start-mip1->q_start);
  s = sseq+ABS(data->s_start-mip1->s_start);


  if (!(qBuf = (Uint1Ptr) MemNew (sizeof(Uint1)*(qlen)))) return;
  if (!(sBuf = (Uint1Ptr) MemNew (sizeof(Uint1)*(slen)))) return;
  
   mip2->qseq=qBuf; 
   mip2->sseq=sBuf; 


  i=0;
  while(i< qlen)
    {  
      
      *qBuf=*q;

      i++;
      q++;
      qBuf++;
    }

  i=0;
  while(i< slen)
    {  
      
      *sBuf=*s;

      i++;
      s++;
      sBuf++;
    }
}

/*________________________________________(DOT_InitCInfo)_____________

  Purpose : Initialize vdp2 (DOTVibDataPtr for window2).

____________________________________________________________________*/

static void DOT_InitCInfo(DOTVibDataPtr vdp, DOTVibDataPtr vdp2, DOTSelDataPtr data)
{
  Char       colBuf[12]={""};
  Int4       i;
  DOTMainDataPtr mip1, mip2;

  /* set up second window parameters */
  mip1 =  vdp->mip;
  mip2=(DOTMainDataPtr) MemNew (sizeof(DOTMainData));

  vdp2->sv=(DOTSeqViewrPtr)MemNew(sizeof(DOTSeqViewr));
  vdp2->Fnt = vdp->Fnt;
  vdp2->HORZ_MARGIN=vdp->HORZ_MARGIN;
  vdp2->VERT_MARGIN=vdp->VERT_MARGIN;
  
  mip2=DOT_InitMainInfo (mip2, mip1->qbsp, mip1->sbsp, mip1->word_size, mip1->tree_limit, data->q_start, data->q_stop, data->s_start, data->s_stop);

  vdp2->sdp.TrampPos=75;
  mip2->qname=mip1->qname;
  mip2->sname=mip1->sname;
  mip2->matrix = mip1->matrix;

  /* sequence viewer initiatize */

  vdp2->sv->do_scale=TRUE;
  vdp2->sv->scaleValue=0;
  vdp2->sv->showLabels=FALSE;
  vdp2->sv->old_primID=-1;

  /* update with every new select */
  mip2->sstrand = mip1->sstrand;
  mip2->qstrand = mip1->qstrand;

  vdp2->curr_slen=data->slen;
  vdp2->curr_qlen=data->qlen;
  mip2->qlen=data->qlen;
  mip2->slen=data->slen;
  vdp2->data=data;
 
  /* hits info */

  vdp2->Fh=vdp->Fh;
  vdp2->charw=vdp->charw;
  
  vdp2->mip= mip2;

}

/*________________________________________(DOT_SVCleanupProc)_____________

  Purpose : Cleanup proc for window2.

____________________________________________________________________*/
/*
static void LIBCALLBACK DOT_SVCleanupProc(GraphiC g, VoidPtr data)
{
  DOTVibDataPtr vdp;

	vdp=(DOTVibDataPtr)data;
	
	if (vdp) DOT_FreeMainInfo(vdp->mip);
   MemFree(vdp->mip);
   MemFree(vdp);
}

*/

/*________________________________________(Init_bufs)_____________

  Purpose : Initialize sequence buffers to NULL.

____________________________________________________________________*/

static void Init_bufs(CharPtr qBuf, CharPtr sBuf, Int4 size)
{
  Int4   i;

  for (i = 0; i<size; i++)
    { 
      qBuf[i]= '\0'; 
      sBuf[i]= '\0'; 
    }
}


/*________________________________________(DOT_SVDisplaySequence)_____________

  Purpose : Draw function for viewer2 of second window.

____________________________________________________________________*/
static void  DOT_SVDisplaySequence(DOTVibDataPtr vdp2, DOTAlnPtr salp)
{
  RecT       rc;
  DOTVibDataPtr vdp;
  DOTSelDataPtr  data;
  DOTMainDataPtr mip1, mip2;
  SegmenT    pict2,nmSeg,hlSeg,clSeg, sbSeg;
  Boolean    match=FALSE;
  Uint1Ptr   q, s, end;
  Int4       a,b,c, pos;
  Boolean    ambig=FALSE;
  Boolean    get_barposition=TRUE;
  Int4       i, k, spaces, prot_threshold=0, bloc_cntr, res_cnt, bloc_size;
  Int4       x, y, x2, y2, q_ahead, s_ahead, qavail, savail;
  Int4       vis_2, buf_len,bufsize, sglen;
  Int4       xstart, xstop, ystart, ystop, hit_start, hit_stop, q_pos, s_pos;
  Int4       V2Height, Fh, margin;
  Int4       xdiff_left,xdiff_right,ydiff_left,ydiff_right,width,height;
  Int4       chw_2, chw_4;
  Int4       qbufstart, qbufstop, sbufstart, sbufstop;
  BaR        vsb;
  CharPtr    qBuf;
  CharPtr    sBuf;
  CharPtr PNTR residue_names;
  Int4Ptr PNTR matrix;
  Int4       num_cls, num_segs, seg_len;
  Int4       q_start, q_stop, s_start, s_stop;      

  data=(DOTSelDataPtr)vdp2->data;
  if (data==NULL) return;
  vdp=data->vdp;
  if (vdp==NULL)  return;
  mip1= vdp->mip;
  mip2= vdp2->mip;

  margin=MAX((StringWidth(mip1->qname)), (StringWidth(mip1->sname)))+10;

  GetPosition(vdp2->sv->v2, &rc);
  Fh=vdp2->Fh;
  vis_2=VIS_LEN/2;
  matrix=mip1->matrix;
  x = rc.left;
  y = 5*Fh; /* on top -cartesian */
  y2 = y -Fh;
  x2 = x + margin;

  chw_2=vdp2->charw/2;
  chw_4=vdp2->charw/4;

  pict2=vdp2->sv->pict2;

  if (mip1->is_na)
    residue_names=na_names;
  else
    residue_names=aa_names;
  
  width=ABS(data->q_stop-data->q_start);
  height=ABS(data->s_stop-data->s_start);

/*   if (salp->q_start<salp->stop) */
/*     { */
/*       qbufstart=ABS(salp->q_start-data->q_start); */
/*       sbufstart=ABS(salp->s_start-data->s_start); */
/*       sbufstop=ABS(salp->s_stop-data->s_start); */

  /* convert from sequence to buffer coordinates */
  q_start = ABS(data->q_start - salp->q_start);
  q_stop = ABS(data->q_start - salp->q_stop);
  s_start = ABS(data->s_start - salp->s_start);
  s_stop = ABS(data->s_start - salp->s_stop);

  xdiff_left=q_start;
  xdiff_right=width-q_stop;
  ydiff_left=s_start;
  ydiff_right=height-s_stop;

  if (xdiff_left<ydiff_left)
    {
      xstart=0;
      ystart=ydiff_left-xdiff_left;
    }
  else
    {
      ystart=0;
      xstart=xdiff_left-ydiff_left;
    }
  
  if (xdiff_right<ydiff_right)
    {
      xstop=q_stop+xdiff_right;
      ystop=s_stop+xdiff_right;
    }
  else
    {
      xstop=q_stop+ydiff_right;
      ystop=s_stop+ydiff_right;
    }

  if ((xstop-xstart)!=(ystop-ystart))/* these should be equal*/
    return;

  buf_len=xstop-xstart;
  hit_start=q_start-xstart;
  hit_start+=hit_start/10;
  hit_stop=q_stop-xstart;
  hit_stop+=hit_stop/10;

  
  q_pos=(mip1->qstrand==Seq_strand_plus)?data->q_start+xstart:data->q_start-xstart;
  s_pos=(mip1->qstrand==Seq_strand_plus)?data->s_start+ystart:data->s_start-ystart;

  s=mip2->sseq+ystart;
  q=mip2->qseq+xstart;
  end=s+buf_len;
  
  /* length of sequence attached to a seg */
  seg_len=(rc.right-rc.left)/vdp2->charw;

  /* total num of segs */
  num_segs=buf_len/seg_len;
  if (seg_len%buf_len)
    num_segs++;
  if (num_segs==0)
      num_segs=1; /* you need at least one */
      

  /* num of seg clusters */
  num_cls=num_segs/5;
  if (5%num_segs)
    num_cls++;
  if (num_cls==0)
      num_cls=1; /* you need at least one */

  bufsize=seg_len+(seg_len/10)+2;

  qBuf=(CharPtr)MemNew(sizeof(Char)*bufsize);
  sBuf=(CharPtr)MemNew(sizeof(Char)*bufsize);
  Init_bufs(qBuf,sBuf,seg_len);

  bufsize=seg_len+(seg_len/10)+1;
   

/*   Ambiguous=(is_na)?(*q > 3 || *s > 3):(*q >24 || *s >24 || *q<1 || *s<1 ); */

  i=0;
  a=0;b=0;c=0;
  spaces=0;res_cnt=0;bloc_cntr=0;bloc_size=1;
  k=1;
  
  for (a=0; a<num_cls; a++)
    {
      clSeg=CreateSegment(pict2,a+1,0);
      for(b=0; b<num_segs; b++)
        {
          sglen=seg_len;
          sbSeg=CreateSegment(clSeg,b+1,0);
          Init_bufs(qBuf,sBuf,seg_len);
   
          for(c=0; c<sglen;c++)
            {
              if (!(i<buf_len))
                {
                  a=num_cls;
                  b=num_segs;
                  sglen=c;
                  if (hit_start<=i && i<=(hit_stop+bloc_size))
                    {
                      DOT_HlSeq(sbSeg, rc, margin, i+1, bloc_size, hit_start, hit_stop, Fh,chw_2);
                    }
                  goto end;
                }

              if (mip1->is_na)
                {
                  if (*q > 3 || *s > 3)
                    ambig=TRUE;
                }
              else
                {
                  if (*q >24 || *s >24 || *q<1 || *s<1 )
                    ambig=TRUE;
                }

              if (ambig)
                {
                  if (s<end) /* not end */
                    {
                      qBuf[c]=*residue_names[(int)*q];
                      sBuf[c]=*residue_names[(int)*s];
                      ambig=FALSE;
                      goto skip; /* continue */
                    }
                  else
                    {
                      a=num_cls;
                      b=num_segs;
                      sglen=c;
                      if (hit_start<=i && i<=(hit_stop+bloc_size))
                        {
                          DOT_HlSeq(sbSeg, rc, margin, i+1, bloc_size, hit_start, hit_stop, Fh,chw_2);
                        }
                      goto end;
                    }
                }
                        
              if (*s == *q)
                match=TRUE;
            
            
              qBuf[c]=*residue_names[(int)*q];
              sBuf[c]=*residue_names[(int)*s];
              
              if (match)
                {
                  DOT_HlMatch(sbSeg, rc, margin, i, Fh, chw_2);
                  match=FALSE;
                }
              
            skip:

              if (!(k % BLOCK_SIZE))
                {
                  c++;
                  sglen++;
                  i++;
                  buf_len++;
                  bloc_size=0;

                  qBuf[c] = ' ';
                  sBuf[c] = ' ';
                  if(match)
                    spaces++;
                  
                  bloc_cntr++; 
                  res_cnt=i-bloc_cntr;
                  DOT_DrawScale(sbSeg, vdp2->mip, rc, margin, res_cnt, bloc_cntr, Fh, q_pos, s_pos, chw_2, chw_4);
                  if (hit_start<=i && i<(hit_stop+BLOCK_SIZE))
                    {
                      DOT_HlSeq(sbSeg, rc, margin, i, BLOCK_SIZE, hit_start, hit_stop, Fh,chw_2);
                      if (get_barposition)
                        {
                          pos=DOT_WorldtoScreen(i,chw_2);
                          vdp2->sv->barp=pos+rc.left+margin;
                          get_barposition=FALSE;
                        }
                    }
                }
              
              i++;
              k++;
              q++;
              s++;
              bloc_size++;
              
            }
        end:
          
          /* attach buffer */

          pos=DOT_WorldtoScreen(i-sglen, chw_2);
          AddAttribute(sbSeg, COLOR_ATT, RED_COLOR, 0,0,0,0);
          AddLabel(sbSeg, x2+pos, y, qBuf , SMALL_TEXT, 0, UPPER_RIGHT, 0);
          AddAttribute(sbSeg, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
          AddLabel(sbSeg, x2+pos, y2, sBuf , SMALL_TEXT, 0, UPPER_RIGHT, 0);
 

          
        }
    }

/*write seq names */

  nmSeg=CreateSegment(pict2, num_cls+1, 0);
  AddAttribute(nmSeg, COLOR_ATT, MAGENTA_COLOR, 0,0,0,0);
  AddLabel(nmSeg, x, y, mip1->qname, SMALL_TEXT, 0, UPPER_RIGHT, 0);
  AddLabel(nmSeg, x, y2, mip1->sname, SMALL_TEXT, 0, UPPER_RIGHT, 0);
    

  qBuf=(CharPtr)MemFree(qBuf);
  sBuf=(CharPtr)MemFree(sBuf);

  return;
    
}

/*________________________________________(DOT_SVPopulateSequenceViewer)______

  Purpose : Calls draw function for viewer2 of second window.

____________________________________________________________________*/
static Boolean DOT_SVPopulateSequenceViewer(DOTVibDataPtr vdp2)
{
  RecT    rc;
  Int4    x;
  BaR     hsb;
  Char    infoBuf[255];
  DOTSelDataPtr data;
  DOTAlnPtr salp;

  data=(DOTSelDataPtr)vdp2->data;
  if (data==NULL) return FALSE;

  ResetViewer(vdp2->sv->v2);
  vdp2->sv->pict2=DeletePicture(vdp2->sv->pict2);
  vdp2->sv->pict2=CreatePicture();
  GetPosition(vdp2->sv->v2, &rc);

  salp=vdp2->sv->salp;

  if (salp==NULL)
    {
      AddAttribute(vdp2->sv->pict2, COLOR_ATT, RED_COLOR, 0,0,0,0);
      AddLabel(vdp2->sv->pict2, rc.left, rc.bottom-rc.top, "- CLICK on a diagonal to view the aligned sequence -", SMALL_TEXT, 0, UPPER_RIGHT, 0);
       
      /*reset infopanel2 title*/

      SetTitle(vdp2->Infopanel, vdp2->iInfo);
    }
  else
    {
      DOT_SVDisplaySequence(vdp2, salp);

       /*reset infopanel2 title*/

      sprintf(infoBuf, "Selected ..  %s (horizontal) [%ld..%ld]   vs.\n             %s (vertical) [%ld..%ld]", vdp2->mip->qname, /* data->q_start+ */(long)salp->q_start, /* data->q_start+ */(long)salp->q_stop-1, vdp2->mip->sname, /* data->s_start+ */(long)salp->s_start, /* data->s_start+ */(long)salp->s_stop-1);
      SetTitle(vdp2->Infopanel,infoBuf);
    }
  
  
  AttachPicture(vdp2->sv->v2, vdp2->sv->pict2, vdp2->sv->barp, INT4_MIN, UPPER_CENTER, 1, 1, NULL);
  ArrowCursor();
  return TRUE;
 
}


/*________________________________________(DOT_SVFindHit)_____________

  Purpose : gets sequence coordinates of clicked diag (viewer2, window2)..

____________________________________________________________________*/
DOTAlnPtr DOT_SVFindHit(DOTVibDataPtr vdp2, Int4 primID)
{
  DOTAlnPtr  salp;
  Int4     index, s_start, q_start, q_stop, len;
  
  if (primID<1 || primID>vdp2->mip->index)
    return NULL;

  if (vdp2->sv->salp == NULL)
    salp=(DOTAlnPtr)MemNew(sizeof(DOTAln));
  else
    salp=vdp2->sv->salp;

  index=primID-1;
  q_start=vdp2->mip->hitlist[index]->q_start;
  s_start=vdp2->mip->hitlist[index]->s_start;
  len=vdp2->mip->hitlist[index]->length;

  salp->q_start=q_start;
  salp->q_stop=(vdp2->mip->qstrand == Seq_strand_plus)?q_start+len:q_start-len;
  salp->s_start=s_start;
  salp->s_stop=(vdp2->mip->sstrand==Seq_strand_plus)?s_start+len:s_start-len;
  salp->primID=(Uint2)primID;
  
  return salp;
}



/*________________________________________(DOT_SVGetDiag)_____________

  Purpose : Calls func to get sequence coordinates of clicked diag (viewer2, window2).

____________________________________________________________________*/

static Boolean DOT_SVGetDiag (Uint2 segID, Uint2 primID, VoidPtr userdata)

{
  DOTVibDataPtr  vdp2;
  DOTAlnPtr     salp;

  vdp2 = (DOTVibDataPtr) userdata;
  if (vdp2 == NULL || segID !=1) 
    return FALSE;

  salp=DOT_SVFindHit(vdp2, primID);
  userdata=(Pointer)salp;
  return TRUE;
}
/*________________________________________(DOT_SVClickProc)_____________

   Purpose : Click proc for viewer1 of second window.

____________________________________________________________________*/
static void DOT_SVClickProc(VieweR v, SegmenT seg, PoinT pt)
{
  DOTVibDataPtr  vdp2;
  SegmenT     new_seg;
  PrimitivE   old_prim, new_prim;
  Uint2       primCT, primID=0, segID;
  Int1        highlight;
  DOTAlnPtr     salp=NULL;

  vdp2=(DOTVibDataPtr)GetObjectExtra((WindoW)ParentWindow(v));
  new_seg=FindSegment(v, pt, &segID, &primID, &primCT);
  if (!(new_seg == vdp2->sv->seg1)) 
    new_seg=NULL;


  if (vdp2->sv->old_primID>=0) /* previous selection exists */
    {

      if (new_seg==NULL)
        {
          new_seg=FindSegment(v, vdp2->sv->old_pt, &segID, &primID, &primCT);
          old_prim = GetPrimitive (new_seg, vdp2->sv->old_primID);
          GetPrimDrawAttribute (old_prim, NULL, NULL, NULL, NULL, NULL, &highlight);
          if (highlight != PLAIN_PRIMITIVE) 
            HighlightPrimitive (v, seg, old_prim, PLAIN_PRIMITIVE);
          vdp2->sv->old_primID=-1;

          goto end;
        }
          old_prim = GetPrimitive (new_seg, vdp2->sv->old_primID);
          GetPrimDrawAttribute (old_prim, NULL, NULL, NULL, NULL, NULL, &highlight);
          if (highlight != PLAIN_PRIMITIVE) 
            HighlightPrimitive (v, seg, old_prim, PLAIN_PRIMITIVE);

    }
  
  if (new_seg!=NULL)
    {
      vdp2->sv->old_primID=primID;
      vdp2->sv->old_pt=pt;
      
      new_prim = GetPrimitive (new_seg, primID);
      GetPrimDrawAttribute (new_prim, NULL, NULL, NULL, NULL, NULL, &highlight);
      
      if (highlight != PLAIN_PRIMITIVE) 
        {
          HighlightPrimitive (v, new_seg, new_prim, PLAIN_PRIMITIVE);
          salp=NULL;
        } 
      else 
        {
          HighlightPrimitive (v, new_seg, new_prim, FRAME_PRIMITIVE);
          salp = DOT_SVFindHit(vdp2, primID);
        }
     }
 end:
  vdp2->sv->salp=salp;
  DOT_SVPopulateSequenceViewer (vdp2);

}


/*________________________________________(DOT_SVPopulateDiagViewer)_____________

  Purpose : Calls draw function for diags.

____________________________________________________________________*/
static Boolean DOT_SVPopulateDiagViewer(DOTVibDataPtr vdp2)
{
  Int4       scale, index;
  Char       str[16], infoBuf[255];
  DOTSelDataPtr  data;
  BaR        vsb;
  
  
  data=(DOTSelDataPtr)vdp2->data;
  
  ResetViewer(vdp2->sv->v1);
  vdp2->sv->pict1=DeletePicture(vdp2->sv->pict1);
  vdp2->sv->pict1=CreatePicture();
   
  if (DOT_SVDisplayDiags(vdp2, data)==FALSE)
    return FALSE;
  
  if (vdp2->sv->do_scale) 
    {
      for (index=1; index<MAXZOOMSCALEVAL; index++) 
        {
          sprintf (str, "%ld", (long) (zoomScaleVal [index]));
          PopupItem (vdp2->sv->scale, str);
        }
      SetValue (vdp2->sv->scale, vdp2->sv->scaleIndex);
      vdp2->sv->do_scale = FALSE;

    }

  SafeShow(vdp2->sv->scale);

 
  AttachPicture(vdp2->sv->v1, vdp2->sv->pict1, INT4_MIN, INT4_MAX, LOWER_RIGHT,  vdp2->sv->scaleValue , vdp2->sv->scaleValue , NULL);
  SetViewerProcs (vdp2->sv->v1, DOT_SVClickProc, NULL, NULL, NULL);
  ArrowCursor();
  return TRUE;
  
}


/*____________________________________________(DOT_ChangeSequenceWindowCutoff)______

  Purpose : Change threshold for diag using threshold ramp.

____________________________________________________________________*/
void DOT_ChangeSequenceViewerCutoff (BaR b, GraphiC g, Int2 new, Int2 old) 
{
  DOTVibDataPtr vdp2;
  WindoW     w;
  RecT       rcP;

  w=(WindoW)ParentWindow(b);
  vdp2 = (DOTVibDataPtr) GetObjectExtra (w);

  vdp2->sdp.TrampPos = new+20;
  DOT_SVPopulateDiagViewer(vdp2);

}

/*________________________________________(DOT_SVChangeScale)_____________

  Purpose : Change scale.

____________________________________________________________________*/
static void DOT_SVChangeScale (PopuP p)

{
  DOTVibDataPtr vdp2;
  Int4       index;

  vdp2 = (DOTVibDataPtr) GetObjectExtra (p);
  if (vdp2 != NULL) 
    {
      index = GetValue (vdp2->sv->scale);
      if (index <= MAXZOOMSCALEVAL && index > 0) 
        {
          vdp2->sv->scaleValue = zoomScaleVal [index];
        } 
      else 
        {
          vdp2->sv->scaleValue = 1;
        }

      DOT_SVPopulateDiagViewer (vdp2);
    }
}
/*________________________________________(DOT_SVChangeLabels)_____________

  Purpose : Show or Hide coordinate labels for diags.

____________________________________________________________________*/
static void DOT_SVChangeLabels (GrouP g)

{
  DOTVibDataPtr  vdp2;

  vdp2 = (DOTVibDataPtr) GetObjectExtra (g);
  if (vdp2 != NULL) 
  {
      vdp2->sv->showLabels=(Boolean)(GetValue(vdp2->sv->Labels)==1);
      DOT_SVPopulateDiagViewer (vdp2);
  }
}

/*________________________________________(DOT_SVCalculateScaling)_____________

  Purpose : Estimates size of picture.

____________________________________________________________________*/
static void DOT_SVCalculateScaling (DOTVibDataPtr vdp2)

{
  RecT   r, world;
  Int4   index, r_hgt, r_wdt, w_hgt, w_wdt, scale;
  double f1, f2;

  w_hgt=vdp2->mip->slen+(vdp2->mip->slen*0.15);
  w_wdt=vdp2->mip->qlen+(vdp2->mip->qlen*0.15);

  GetPosition(vdp2->sv->v1, &r);
  r_hgt=r.bottom-r.top;
  r_wdt=r.right-r.left;

  f1=(float)w_hgt/r_hgt;
  f2=(float)w_wdt/r_wdt;

  scale=MAX(ceil(f1), ceil(f2));

  for (index=1; index<MAXZOOMSCALEVAL; index++) 
    {
      if (zoomScaleVal [index]>= scale)
        {
          vdp2->sv->scaleValue=zoomScaleVal[index];
          vdp2->sv->scaleIndex=index;
          return;
        }
    }

  vdp2->sv->scaleValue=zoomScaleVal[MAXZOOMSCALEVAL-1];
  vdp2->sv->scaleIndex=MAXZOOMSCALEVAL-1;

}

/*________________________________________(DOT_ResizeFeatWindow)_____________

  Purpose : Resize function for window2.

____________________________________________________________________*/

static void DOT_ResizeFeatWindow(WindoW w)
{
  RecT     rcDlg,rcQry,rcSub, rcVsb,rcHsb, rcQi, rcSi;
  Int2     height,width,vsbWidth,in,gap,hsbHeight,QueryHeight,SubjectHeight, halfw;
  BaR      vsb,hsb;
  DOTFeatListPtr flp;
  WindoW   temport;


  flp=(DOTFeatListPtr)GetObjectExtra(w);
  temport=SavePort(w);

  ObjectRect(w,&rcDlg);
  width= rcDlg.right-rcDlg.left;
  halfw=width/2;
  height= rcDlg.bottom-rcDlg.top;
  
  SafeHide(flp->Query);
  SafeHide(flp->Subject);
  SafeHide(flp->QInfo);
  SafeHide(flp->SInfo);
  Update();
  
  vsb = GetSlateVScrollBar ((SlatE) flp->Query);
  hsb = GetSlateHScrollBar ((SlatE) flp->Subject);
  
  GetPosition(flp->Query,&rcQry);
  GetPosition(flp->Subject,&rcSub);
  GetPosition(flp->QInfo, &rcQi);
  GetPosition(flp->SInfo, &rcSi);
  GetPosition(vsb,&rcVsb);
  GetPosition(hsb,&rcHsb);

  in=2;
  gap=10;
  vsbWidth=rcVsb.right-rcVsb.left;
  hsbHeight=rcHsb.bottom-rcHsb.top;
  QueryHeight=rcQry.bottom-rcQry.top;
  SubjectHeight=rcSub.bottom-rcSub.top;
  
  /*new sizes for the viewers*/	
  rcSub.right=width-in-vsbWidth;
  rcSub.bottom=height-in-hsbHeight;
  rcSub.left=halfw+in;
  rcSi.left=rcSub.left;
  rcQry.right=halfw-in-vsbWidth;
  rcQry.bottom=height-in-hsbHeight;

 
  
  /*set the new sizes*/
  SetPosition(flp->Query,&rcQry);
  AdjustPrnt (flp->Query, &rcQry, FALSE);
  SetPosition(flp->Subject,&rcSub);
  AdjustPrnt (flp->Subject, &rcSub, FALSE);
  SetPosition(flp->SInfo,&rcSi);
  AdjustPrnt (flp->SInfo, &rcSi, FALSE);

  AttachPicture (flp->Query,flp->segQuery, INT4_MIN, flp->vert_Qpos, UPPER_LEFT,1 , 1, NULL); 
  SetViewerProcs (flp->Query, DOT_QViewerClickProc, NULL, NULL, NULL);

  AttachPicture (flp->Subject,flp->segSubject, INT4_MIN, flp->vert_Spos, UPPER_LEFT,1 , 1, NULL);  
  SetViewerProcs (flp->Subject, DOT_SViewerClickProc, NULL, NULL, NULL);

  SafeShow(flp->QInfo);
  SafeShow(flp->SInfo);
  SafeShow(flp->Query);
  SafeShow(flp->Subject);
  RestorePort(temport);
  Update();
}

/*________________________________________(DOT_ResizeSequenceWindow)_____________

  Purpose : Resize function for window2.

____________________________________________________________________*/

static void DOT_ResizeSequenceWindow(WindoW w)
{
  Int4     lmargin;
  RecT     rcDlg,rcV1,rcV2, rcVsb,rcHsb;
  Int2     height,width,gap,vsbWidth,in,hsbHeight,V1Height,V2Height;
  BaR      vsb,hsb;
  DOTVibDataPtr vdp2;
  Boolean  is_visible1, is_visible2;
  WindoW   temport;


  vdp2=(DOTVibDataPtr)GetObjectExtra(w);
  temport=SavePort(w);
  ObjectRect(w,&rcDlg);
  width= rcDlg.right-rcDlg.left;
  height= rcDlg.bottom-rcDlg.top;
  
  SafeHide(vdp2->sv->v1);
  SafeHide(vdp2->sv->v2);
  SafeHide(vdp2->Infopanel);
  Update();
  
  vsb = GetSlateVScrollBar ((SlatE) vdp2->sv->v1);
  hsb = GetSlateHScrollBar ((SlatE) vdp2->sv->v1);
  
  GetPosition(vdp2->sv->v1,&rcV1);
  GetPosition(vdp2->sv->v2,&rcV2);
  GetPosition(vsb,&rcVsb);
  GetPosition(hsb,&rcHsb);

  gap=2;
  in=vdp2->Fh;
  lmargin=10;
  vsbWidth=rcVsb.right-rcVsb.left;
  hsbHeight=rcHsb.bottom-rcHsb.top;
  V1Height=rcV1.bottom-rcV1.top;
  V2Height=9*vdp2->Fh;
  
  /*new sizes for the viewers*/	
  rcV1.left=lmargin;
  rcV1.right=width-gap-vsbWidth-rcV1.left;
  rcV2.bottom=height-gap-in-hsbHeight;
  rcV2.top=rcV2.bottom-V2Height;
  rcV1.bottom=rcV2.top-gap-hsbHeight;
  rcV2.left=lmargin;
  rcV2.right=width-gap-vsbWidth-rcV2.left;
  
  /*set the new sizes*/
  SetPosition(vdp2->sv->v1,&rcV1);
  AdjustPrnt (vdp2->sv->v1, &rcV1, FALSE);
  
  SetPosition(vdp2->sv->v2,&rcV2);
  AdjustPrnt (vdp2->sv->v2, &rcV2, FALSE);

  
  if ((rcV1.left<rcV1.right)&&(rcV1.top<rcV1.bottom))
    is_visible1=TRUE;
  if ((rcV2.left<rcV2.right)&&(rcV2.top<rcV2.bottom))
    is_visible2=TRUE;

  /*update viewers*/
  if (Visible (vdp2->sv->v1) && AllParentsVisible (vdp2->sv->v1))
    ViewerWasResized(vdp2->sv->v1);
  if (Visible (vdp2->sv->v2) && AllParentsVisible (vdp2->sv->v2))
    ViewerWasResized(vdp2->sv->v2);
  
  if (vdp2->sv->do_scale!=TRUE) /* window resized */
    {
      DOT_SVPopulateDiagViewer(vdp2);
      DOT_SVPopulateSequenceViewer(vdp2);
    }
                       
  SafeShow(vdp2->sv->v1);
  SafeShow(vdp2->sv->v2);
  SafeShow(vdp2->Infopanel);
  ArrowCursor();
  RestorePort(temport);
  Update();
}

 static Boolean DOT_DeleteCursorPrim (SegmenT seg, PrimitivE prim, Uint2 segID,
                           Uint2 primID, Uint2 primCt, VoidPtr userdata)
{
  DeletePrim(seg,prim);
  return TRUE;
}


static Int4 DOT_PointCursorOnFeature(SegmenT seg, DOTSelFeatPtr feat_list, Int2 fontHt, Boolean found)
{
  Int4  yBase, cursorHt;
  Boolean is_first =TRUE;
  Int4  vertbar_pos;

  cursorHt=fontHt/4;

  if (!found)
    {
      AddAttribute(seg, COLOR_ATT, BLUE_COLOR, 0, 0, 0, 0);
      yBase = (-1)*(feat_list->feat_num*(fontHt)-cursorHt);
      vertbar_pos=yBase+3*fontHt;
      AddRectangle(seg, 0,yBase, 2,yBase-2*cursorHt,RIGHT_ARROW,TRUE,0);
    }
  else
    {
      AddAttribute(seg, COLOR_ATT, RED_COLOR, 0, 0, 0, 0);      
      while (feat_list)
        {
          yBase = (-1)*((feat_list->feat_num-1)*fontHt+cursorHt);
          if (is_first)
            {
              vertbar_pos=yBase+3*fontHt;
              is_first=FALSE;
            }
          AddRectangle(seg, 0,yBase, 2,yBase-2*cursorHt,RIGHT_ARROW,TRUE,0);
          feat_list=feat_list->next;
        }
    }
  return (vertbar_pos);
}

static DOTFeatPtr   DOT_FindNextDfp(DOTFeatIndexPtr fdindex, DOTFeatPtr dfp, Boolean use_first)
{
  if (!use_first)
    dfp=dfp->next;
  while (dfp != NULL)
    {
      if (fdindex[dfp->type].show)
        return dfp;
      dfp=dfp->next;
    }
  return NULL;
}


static DOTSelFeatPtr DOT_FindFeatureInViewer(DOTFeatIndexPtr fdindex, DOTRowPtr drp, Int4 cursor_pos)
{
  DOTFeatPtr dfp_head;
  DOTSelFeatPtr feat_list=NULL;
  Int4       i;

  /* fix this so that it also works for minus strand */
  i=1;
  dfp_head = DOT_FindNextDfp(fdindex, drp->dfp, TRUE);
  while (dfp_head != NULL)
    {
      if (dfp_head->left<= cursor_pos && dfp_head->right >= cursor_pos)
        {
          if (!feat_list)
            {
              feat_list=(DOTSelFeatPtr)MemNew(sizeof(DOTSelFeat));
            }
          else
            {
              feat_list->next=(DOTSelFeatPtr)MemNew(sizeof(DOTSelFeat));
              feat_list=feat_list->next;
            }
          feat_list->feat_num = i;
        }
      if (dfp_head->right > cursor_pos)
        {
          if (feat_list) 
            return feat_list;
          else 
            return NULL;
        }
      i++;
      dfp_head=DOT_FindNextDfp(fdindex, dfp_head, FALSE);
    }
  return NULL;
}



static DOTSelFeatPtr DOT_FindBetweenFeats(DOTFeatIndexPtr fdindex, DOTRowPtr drp, Int4 cursor_pos)
{
  DOTFeatPtr dfp1, dfp2;
  DOTSelFeatPtr feat_list=NULL;
  Int4       i;

  /* fix this so that it also works for minus strand */
  i=1;
  dfp1 = DOT_FindNextDfp(fdindex, drp->dfp, TRUE);
  if (!dfp1) goto end;
/*   dfp2 = DOT_FindNextDfp(fdindex, dfp1, FALSE); */
/*   if (!dfp2) goto end; */

  while (dfp1 != NULL /* && dfp2 != NULL */)
    {
      if (dfp1->left > cursor_pos)
        {
          if (!feat_list)
            {
              feat_list=(DOTSelFeatPtr)MemNew(sizeof(DOTSelFeat));
              feat_list->feat_num = i-1;
              return (feat_list);
            }
        }
    
      i++;
      dfp1=DOT_FindNextDfp(fdindex, dfp1, FALSE);
/*       dfp2=DOT_FindNextDfp(fdindex, dfp1, FALSE); */
    }

 end:
  feat_list=(DOTSelFeatPtr)MemNew(sizeof(DOTSelFeat));
  feat_list->next=NULL;
  feat_list->feat_num = MAX(i-1, 1);
  return (feat_list);

}
/*******************************************************************************

  Function : DOT_AddFeatureToSegment()
  
  Purpose : analyse one feature and add it to the Feature Viewer

*******************************************************************************/
static void DOT_AddFeatureToSegment(SeqMgrFeatContextPtr context,
			DOTPopFeatPtr pfp) 
{
PrimitivE  prim;
Int4	   yBase, xMargin;
Char     str[50]; 


 sprintf(str, "%s (%ld - %ld)", pfp->dfp_cur->label, (long)pfp->dfp_cur->left, (long)pfp->dfp_cur->right); 
 AddAttribute(pfp->TopParentSeg, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);
 yBase = (-1)*pfp->nfeats*pfp->fontHeight;
 xMargin=5;

 if (context->numivals==1){

 prim = AddLabel(pfp->TopParentSeg, xMargin, yBase, str, SMALL_TEXT, 0, UPPER_RIGHT, 0);
 SetPrimitiveIDs(prim,context->entityID,context->itemID,OBJ_SEQFEAT, 0);
	}
 /*
	else{
		DOT_AddSegmentedFeature(pfp->TopParentSeg,context->numivals,context->ivals, context->strand,context->entityID,context->itemID,yBase);
	}
 */
}

/*******************************************************************************

  Function : DOT_AddFeaturesToViewer()
  
  Purpose : callback called by SeqMgr (see also DOT_PopFeatureViewers(), below)

*******************************************************************************/
static Boolean LIBCALLBACK DOT_AddFeaturesToViewer(SeqFeatPtr sfp, 
			SeqMgrFeatContextPtr context)
{
BioseqPtr            parent;
SeqMgrSegmentContext contextPart;
SeqMgrFeatContextPtr context2;
Int4                 segID;
DOTFeatPtr           dfp_new=NULL;
DOTRowPtr            drp=NULL;
DOTFeatIndexPtr      fdindex=NULL; 
DOTPopFeatPtr         pfp=NULL;

 
	pfp = (DOTPopFeatPtr) context->userdata;

	if (pfp==NULL){
		return(FALSE);
	}


   drp=pfp->drp;
   
   if (pfp->dfp_cur == NULL) 
     {
       drp->dfp=(DOTFeatPtr)MemNew(sizeof(DOTFeat));
       pfp->dfp_cur=drp->dfp;
       dfp_new = drp->dfp;
     }
   else
     {
       if(!(pfp->dfp_cur->next=(DOTFeatPtr)MemNew(sizeof(DOTFeat)))) return(FALSE);
       pfp->dfp_cur=pfp->dfp_cur->next;
       dfp_new = pfp->dfp_cur;
     }

   if (!dfp_new) return(FALSE);

       
   if (context->strand>Seq_strand_minus ||
       context->strand==Seq_strand_unknown) 
     {
       pfp->dfp_cur->strand=Seq_strand_plus;
     }

       dfp_new->label=(CharPtr)FeatDefTypeLabel(sfp);
       dfp_new->left=context->left;
       dfp_new->right=context->right;
       dfp_new->type=context->featdeftype;
  
       /*add a feature to the viewer*/
/*      DOT_AddFeatureToSegment(context,pfp);  */
       

       pfp->nfeats++;
       fdindex = pfp->featindex;
       
       if (!fdindex[context->featdeftype].present)
         {
           fdindex[context->featdeftype].present=TRUE;
           fdindex[context->featdeftype].show=TRUE;
           fdindex[context->featdeftype].label=dfp_new->label;
  
         }
       
	return(TRUE);
}


static Boolean  DOT_DeletePrims(SegmenT seg, PrimitivE prim, Uint2 segID,
                           Uint2 primID, Uint2 primCt, VoidPtr userdata)
{
  
  DeletePrim(seg, prim);
  return TRUE;
   
}

static void DOT_PlaceCursors(VieweR viewer, SegmenT seg, SegmenT segCursor, DOTFeatIndexPtr fdindex, DOTRowPtr drp, Int4 seqpos, Int2 fontHt, Int4Ptr ypos, Boolean from_click, Int4 fpos)
{
  Boolean found;
  DOTSelFeatPtr foundfeats;
  Int4    Y_pos;
  BaR     bar;
  WindoW  temport;
  RecT    rcP;

  if (!from_click)
    {
      foundfeats=DOT_FindFeatureInViewer(fdindex, drp, seqpos); 
      if (!foundfeats)
        {
          foundfeats=DOT_FindBetweenFeats(fdindex, drp, seqpos);
          found=FALSE;
        }
      else
        found=TRUE;
    }
  else
    {
      foundfeats=(DOTSelFeatPtr)MemNew(sizeof(DOTSelFeat));
      foundfeats->feat_num=fpos;
      foundfeats->next=NULL;
      found=TRUE;
    }
  ExploreSegment (segCursor, NULL, DOT_DeleteCursorPrim);
  Y_pos=DOT_PointCursorOnFeature(segCursor, foundfeats, fontHt, found);
  *ypos=Y_pos;

/*   if (from_main_viewer)  */
/*     { */
/*       bar=GetSlateHScrollBar((SlatE) viewer); */
/*       Y_pos=GetBarValue(bar); */
/*     } */
/*   else */
/*     SetBarValue(GetSlateHScrollBar((SlatE) viewer), Y_pos); */

  MemFree(foundfeats); /* fix this !!*/
  if (from_click)
    {
      temport = SavePort(ParentWindow(viewer));
      Select(viewer);
      GetPosition(viewer, &rcP);
      InsetRect(&rcP, -1, -1);
      InvalRect(&rcP);
      RestorePort(temport);
      Update();
    }
  else
    {
      AttachPicture (viewer, seg, INT4_MIN, Y_pos, UPPER_LEFT,1 , 1, NULL); 
    }
    
  
}


static void DOT_QViewerClickProc(VieweR v, SegmenT seg, PoinT pt)
{
  DOTFeatListPtr flp;
  DOTFeatPtr     dfp;
  DOTSelDataPtr  data;
  Int4           i, xmidPt, comp, j;
  Uint2          segID, primID, primCT;
  PrimitivE      prim;
  RecT           rc;
  Char           infoBuf[255];

  flp=(DOTFeatListPtr)GetObjectExtra((WindoW)ParentWindow(v));
  if (!flp) return;

  if (FindSegPrim(v, pt, NULL, NULL, &prim))
    {
      FindSegment(v, pt, &segID, &primID, &primCT);
      if (segID==1)
        dfp=flp->query_drp->dfp;
      else
        return;
      i=1;
      dfp=DOT_FindNextDfp(flp->featindex, dfp, TRUE);
      if (primID!=1)
        {
          for (; i<primID; i++)
            {
              dfp=DOT_FindNextDfp(flp->featindex, dfp, FALSE);
              if (dfp==NULL) return;
            }
        }


      if (dfp==NULL) return;

      data=(DOTSelDataPtr)flp->data;
      xmidPt=ABS(dfp->right+dfp->left)/2;
      data->q_start=xmidPt;
      GetPosition(data->vdp->panel, &rc);
      DOT_AddRectMargins(&rc, data->vdp);
      comp=data->vdp->comp;
      xmidPt=ABS(xmidPt-flp->mip->q_start);
      curpnt.x=MIN(rc.left+(xmidPt/comp), rc.right);
      
      DOT_PlaceCursors(v, seg, flp->segQCursor, flp->featindex, flp->query_drp, dfp->left, flp->fontHt, &(flp->vert_Qpos), TRUE, i);
      SetViewerProcs (v, DOT_QViewerClickProc, NULL, NULL, NULL);
      sprintf(infoBuf, "Hairs .. X-axis (%s) [%ld]  vs.  Y-axis (%s) [%ld]", flp->mip->qname, (long)data->q_start, flp->mip->sname, (long)data->s_start);
      SetTitle(data->vdp->Infopanel,infoBuf);
      DOT_UpdateMainPanel(data->vdp, TRUE);
    }
}

static void DOT_SViewerClickProc(VieweR v, SegmenT seg, PoinT pt)
{
  DOTFeatListPtr flp;
  DOTFeatPtr     dfp;
  DOTSelDataPtr  data;
  Char           infoBuf[255];
  RecT           rc;
  Int4           i, ymidPt, comp;
  Uint2          segID, primID, primCT;
  PrimitivE      prim;


  flp=(DOTFeatListPtr)GetObjectExtra((WindoW)ParentWindow(v));
  if (!flp) return;

  if (FindSegPrim(v, pt, NULL, NULL, &prim))
    {
      FindSegment(v, pt, &segID, &primID, &primCT);
      if (segID==1)
        dfp=flp->subject_drp->dfp;
      else
        return;
      i=1;
      dfp=DOT_FindNextDfp(flp->featindex, dfp, TRUE);
      if (primID!=1)
        {
          for (; i<primID; i++)
            {
              dfp=DOT_FindNextDfp(flp->featindex, dfp, FALSE);
              if (dfp==NULL) return;
            }
        }


      if (dfp==NULL) return;

      data=(DOTSelDataPtr)flp->data;
      ymidPt = ABS(dfp->right+dfp->left)/2;
      data->s_start=ymidPt;
      GetPosition(data->vdp->panel, &rc);
      DOT_AddRectMargins(&rc, data->vdp);
      ymidPt=ABS(ymidPt -flp->mip->s_start);
      comp=data->vdp->comp;
      curpnt.y=MIN(rc.top+(ymidPt/comp), rc.bottom);

      DOT_PlaceCursors(v, seg, flp->segSCursor, flp->featindex, flp->subject_drp, dfp->left, flp->fontHt, &(flp->vert_Spos), TRUE, i);
      SetViewerProcs (v, DOT_SViewerClickProc, NULL, NULL, NULL);
      sprintf(infoBuf, "Hairs .. X-axis (%s) [%ld]  vs.  Y-axis (%s) [%ld]", flp->mip->qname, (long)data->q_start, flp->mip->sname, (long)data->s_start);
      SetTitle(data->vdp->Infopanel,infoBuf);
      DOT_UpdateMainPanel(data->vdp, TRUE);
    }
}




static void DOT_PlaceFeat(DOTFeatListPtr flp, SegmenT seg, Int4 yBase, DOTFeatPtr dfp, Int4 primID)
{
  Char     str[50]; 
  Int4     xMargin;
  
  sprintf(str, "%s (%ld,%ld,%ld)", dfp->label, (long)dfp->left, (long)dfp->right, (long)ABS(dfp->right-dfp->left)); 
  AddAttribute(seg, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);
  xMargin=5;
    
  AddLabel(seg, xMargin, yBase, str, SMALL_TEXT, 0, UPPER_RIGHT, primID);
 
}



static void DOT_UpdateFeatViewer (DOTFeatListPtr flp, VieweR viewer, SegmenT seg, SegmenT segName, SegmenT segCursor, VwrClckProc click, DOTRowPtr drp, Int4 cur_pos, Int4Ptr vert_pos)
{
  Int4 i;
  DOTFeatPtr dfeatp;
  DOTFeatIndexPtr fdindex;

  ExploreSegment(seg, NULL, DOT_DeletePrims);
  dfeatp=drp->dfp;
  fdindex=flp->featindex;

  i = 1;
  
  while (dfeatp != NULL)
    {
      if (fdindex[dfeatp->type].show)
        {
          DOT_PlaceFeat(flp, segName, (-1)*i*flp->fontHt, dfeatp, i);
          i++;
        }
      dfeatp=dfeatp->next;
    }
  if (i>1)
    {
      DOT_PlaceCursors(viewer, seg, segCursor, flp->featindex, drp, cur_pos, flp->fontHt, vert_pos, FALSE, 0);
      SetViewerProcs (viewer, click, NULL, NULL, NULL);

    }
  else
    {
      AddAttribute(seg, COLOR_ATT, RED_COLOR, 0,0,0,0);
      AddLabel(seg, 0, 0, "- no features -", SMALL_TEXT, 0, UPPER_RIGHT, 0);
      AttachPicture (viewer, seg, INT4_MIN, INT4_MAX, UPPER_LEFT,1 , 1, NULL); 
    }
}
/*_________________________________________________________________

  Function : DOT_PopFeatureViewers()
  
  Purpose : populate the Feature Viewer with the features.

___________________________________________________________________*/
static void DOT_PopFeatureViewers(DOTFeatListPtr flp)
{
  BioseqPtr query_bsp, subject_bsp;
  DOTPopFeat dpf;
  DOTRowPtr  drp1, drp2;
  Int2       fontH;
  DOTMainDataPtr mip;

  mip=flp->mip;
  if (!mip) return;
  drp1 = (DOTRowPtr)MemNew(sizeof(DOTRow));
  flp->query_drp=drp1;
  drp2=(DOTRowPtr)MemNew(sizeof(DOTRow));
  flp->subject_drp=drp2;

  fontH = FontHeight();
  flp->fontHt = fontH;
  /* seqmgrexplorefeatures goes through features in order left->right*/
  query_bsp=mip->qbsp;
  memset(&dpf, 0, sizeof(dpf));
  dpf.TopParentView=flp->Query;
  dpf.TopParentSeg= CreateSegment(flp->segQuery, 1, 0);
  dpf.nfeats = 1;
  dpf.fontHeight=fontH;
  dpf.drp=flp->query_drp; 
  dpf.dfp_cur=NULL;
  dpf.featindex=flp->featindex;
  flp->qFeatscount=SeqMgrExploreFeatures (query_bsp, (Pointer) &dpf, DOT_AddFeaturesToViewer, NULL, NULL, NULL);

  subject_bsp=mip->sbsp;
  memset(&dpf, 0, sizeof(dpf));
  dpf.TopParentView=flp->Subject;
  dpf.TopParentSeg= CreateSegment(flp->segSubject, 1, 0) ;
  dpf.nfeats = 1;
  dpf.fontHeight=fontH;
  dpf.drp=flp->subject_drp; 
  dpf.drp->dfp=NULL; 
  dpf.dfp_cur=NULL;
  dpf.featindex=flp->featindex;
  flp->sFeatscount=SeqMgrExploreFeatures (subject_bsp, (Pointer)&dpf, 
                          DOT_AddFeaturesToViewer, NULL, NULL, NULL);
}

static void DOT_HideAcceptProc(ButtoN b)
{

}

static void DOT_HideFeatList(ButtoN b)
{
    WindoW hHideDlg;
    DOTFeatListPtr flp;
    Int4 i, numrows;
    Boolean status;
    DOTFeatIndexPtr fdindex;

	hHideDlg = (WindoW)ParentWindow(b);
	if(hHideDlg == NULL) return;	
	flp = (DOTFeatListPtr) GetObjectExtra(hHideDlg);
   if(flp == NULL) return;
   numrows = flp->numrows;
   fdindex = flp->featindex;

    for(i = 1; i <= numrows; i++)
      {
        if(!GetItemStatus(flp->featList, i))
          {
            fdindex[fdindex[i].deref].show=FALSE;
          }
        else
          fdindex[fdindex[i].deref].show=TRUE;
      }
    DOT_UpdateFeatViewer(flp, flp->Query, flp->segQuery, flp->segQName, flp->segQCursor, DOT_QViewerClickProc, flp->query_drp, ((DOTSelDataPtr)flp->data)->q_start, &(flp->vert_Qpos));
    DOT_UpdateFeatViewer(flp, flp->Subject, flp->segSubject, flp->segSName, flp->segSCursor, DOT_SViewerClickProc, flp->subject_drp, ((DOTSelDataPtr)flp->data)->s_start, &(flp->vert_Spos));
    
    Remove(hHideDlg);
}

static void DOT_SetHideList(DOTFeatListPtr flp)
{
  Int4 i;
  DOTFeatIndexPtr fdindex;
  
  fdindex=flp->featindex;

  for(i = 1; i <= FEATDEF_MAX; i++) 
    {
      if (fdindex[fdindex[i].deref].present)
        {
          if (fdindex[fdindex[i].deref].show)
            SetItemStatus(flp->featList, i, TRUE);
          else
            SetItemStatus(flp->featList, i, FALSE);
        }
    }
}

static void DOT_HideFeatDlg(DOTFeatListPtr flp)
{
    Char szName[181], szName2[161];
    Int4 i, numrows;
    GrouP g, hg;
    ButtoN b;
    DOTFeatIndexPtr fdindex;

    if (flp==NULL) return;
    flp->hFeatDlg=NULL;

    flp->hFeatDlg = MovableModalWindow(-50, -20, -10, -10, "Hide/Show Features", NULL);
    if (flp->hFeatDlg==NULL) return;

    SetObjectExtra (flp->hFeatDlg, (void *)flp, NULL);


    hg = HiddenGroup(flp->hFeatDlg, 1, 2, NULL);
    g = NormalGroup(hg, 1, 2, "Choose Features to show:", systemFont, NULL);

    fdindex = flp->featindex;
    flp->featList = MultiList(g,20, 6, NULL);
    numrows=1;
    for(i = 1; i < FEATDEF_MAX; i++) 
      {
        if (fdindex[i].present)
          {
            ListItem(flp->featList, fdindex[i].label);
            fdindex[numrows].deref=i;
            numrows++;
          }   
      }
    DOT_SetHideList(flp);

    flp->numrows=numrows-1;
    g = HiddenGroup(hg, 2, 1, NULL);
    SetGroupSpacing(g, 15, 15);
    b = DefaultButton(g, "OK", (BtnActnProc) DOT_HideFeatList);

    Show(flp->hFeatDlg);
    return;
}

static void DOT_HideFeatDlgItem(IteM i)
{
  DOTFeatListPtr flp;
  WindoW         FeatWin;

  FeatWin=(WindoW)ParentWindow(i);
  if (FeatWin==NULL)return;

  flp=(DOTFeatListPtr)GetObjectExtra(FeatWin);
  DOT_HideFeatDlg(flp);
}

/*________________________________________(DOT_BuildFeatGUI)_____________

  Purpose : Creates viewers for Features.

____________________________________________________________________*/

static WindoW DOT_BuildFeatGUI (DOTFeatListPtr flp)
{
  WindoW  FeatWin;
  DOTMainDataPtr mip;
  MenU    m1;
  VieweR  v1, v2;
  SegmenT seg1, seg2;
  GrouP   subg1, subg2, g;
  PrompT  pr1, pr2;
  Int2   Margins;
  Char   str[20];

	if (!flp) return(NULL);

	Margins=10*stdCharWidth;
	FeatWin = DocumentWindow(Margins,Margins ,-10, -10,"Features", NULL, DOT_ResizeFeatWindow);
   if (!FeatWin) return(NULL);
   SetObjectExtra (FeatWin, (Pointer) flp, NULL);

   flp->FeatWin=FeatWin;

   m1 = PulldownMenu (FeatWin, "Options");
   CommandItem(m1, "Hide ..", DOT_HideFeatDlgItem);
   CommandItem(m1, "Quit", DOT_CloseFeatWindow);

	mip= flp->mip;


	g=HiddenGroup(FeatWin,2,0,NULL);

   subg1=HiddenGroup(g, 0, 2, NULL);
   sprintf(str, "%s(X-axis)", mip->qname);
   flp->QInfo=StaticPrompt (subg1, str , 0, 0, systemFont, 'l');
	v1=CreateViewer(subg1, 150 , 200, TRUE, TRUE);	
	seg1=CreatePicture();
   
/*    AlignObjects(ALIGN_JUSTIFY, (HANDLE) pr1, (HANDLE) v1, NULL, NULL); */

	/*viewer for close up of features*/
   subg2=HiddenGroup(g, 0, 2, NULL);
   sprintf(str, "%s(Y-axis)", mip->sname);
   flp->SInfo= StaticPrompt (subg2, str , 0, 0, systemFont, 'l');
   v2 = CreateViewer(subg2, 150,200,TRUE,TRUE);
   seg2 = CreatePicture();
 
   flp->Query=v1;
   flp->segQuery=seg1;
	flp->Subject=v2; /* lines on the control panel*/
	flp->segSubject=seg2;
   flp->segQName =CreateSegment(flp->segQuery, 1, 0);
   flp->segSName=CreateSegment(flp->segSubject, 1, 0);
   flp->segQCursor=CreateSegment(flp->segQuery, 2, 0);
   flp->segSCursor=CreateSegment(flp->segSubject, 2, 0);

/*    AlignObjects(ALIGN_JUSTIFY, (HANDLE) pr2, (HANDLE) v2, NULL, NULL); */

   RealizeWindow(FeatWin);
   DOT_ResizeFeatWindow(FeatWin);

	return(FeatWin);	
}

/*________________________________________(DOT_SVBuildDiagViewer)_____________

  Purpose : Creates viewers for second window.

____________________________________________________________________*/

static WindoW DOT_SVBuildDiagViewer(DOTVibDataPtr vdp2)
{
  WindoW    wSequence;
  GrouP     g, s, s2, s3, s4;
  PrompT    pr1, pr2; 
  VieweR    v1,v2;
  RecT      rc;
  SegmenT   pict1,pict2;
  Int2      Margins, pixwidth; 
  MenU      menu; 
  BaR       hsb;
  Char      zoombuf[]={"Decrease scale to zoom in .."};


	if (!vdp2) return(NULL);   

	Margins=10*stdCharWidth;
	wSequence = DocumentWindow(Margins,Margins ,-10, -10, "Selected Region", NULL, DOT_ResizeSequenceWindow);
/*    menu = PulldownMenu (wSequence, "Display"); */
/*    CommandItem (menu, "Close", DOT_CloseSequenceWindow);  */ 
  
	if (!wSequence) return(NULL);
	
   GetPosition (wSequence,&rc);
   pixwidth=600; /* some approximate value */

   /* first top group */

   s = HiddenGroup (wSequence,1, 3, NULL);

   /*threshold bar*/

   s3 = HiddenGroup (s,5, 0, NULL);
   pr2=StaticPrompt (s3, "Threshold-Ramp:", 0, 10,vdp2->Fnt, 'l');
   pr1=StaticPrompt (s3, "  20%", 0, 0, vdp2->Fnt, 'l');
   vdp2->sdp.ScrollBar = ScrollBar (s3, 15, 5, DOT_ChangeSequenceViewerCutoff);
   pr1=StaticPrompt (s3, "100%", 0, 0, vdp2->Fnt, 'l');
   PushButton(s3, "Quit", DOT_CloseSequenceWindow);
   SetObjectExtra(vdp2->sdp.ScrollBar, vdp2, NULL);
   
   CorrectBarMax (vdp2->sdp.ScrollBar, 80); /* 100% */
   CorrectBarValue (vdp2->sdp.ScrollBar, 60);/* 100% */
   
   /* second top group */

   s2 = HiddenGroup (s,0 ,2, NULL);
/*    SetGroupSpacing(s2, 3,0); */
   s4=HiddenGroup (s2, 2, 0, NULL);
   pr1 = StaticPrompt (s4, zoombuf, StringWidth(zoombuf)+10 , popupMenuHeight, vdp2->Fnt, 'l');
#ifdef WIN_MAC
   vdp2->sv->scale = PopupList (s4, TRUE, DOT_SVChangeScale);
#endif

#ifndef WIN_MAC
   vdp2->sv->scale = PopupList (s4, FALSE, DOT_SVChangeScale);
#endif
   sprintf(vdp2->iInfo,"Diag not selected");
   vdp2->Infopanel= StaticPrompt (s2, vdp2->iInfo, pixwidth, vdp2->Fh, vdp2->Fnt, 'l');
   SetObjectExtra (vdp2->sv->scale, vdp2, NULL);
	v1=CreateViewer(s,600,500,TRUE,TRUE);
	pict1=CreatePicture();


   /* bottom group */

	g=HiddenGroup(wSequence,0,2,NULL);
	v2=CreateViewer(g,600,150,FALSE,TRUE);
	pict2=CreatePicture();

	vdp2->sv->w=wSequence;
	vdp2->sv->v1=v1;
	vdp2->sv->pict1=pict1;
	vdp2->sv->v2=v2;
	vdp2->sv->pict2=pict2;

	SetObjectExtra (vdp2->sv->w, (Pointer) vdp2, NULL);
	AlignObjects(ALIGN_JUSTIFY, (HANDLE) v2, (HANDLE)g, NULL, NULL);
	SetColorCell((GraphiC)vdp2->sv->v1, 0,0,255,192);

   RealizeWindow(wSequence);
   DOT_ResizeSequenceWindow(wSequence);

  /* calculate initial scale */
   DOT_SVCalculateScaling(vdp2);

	/*populate the viewer : hits*/
   if (DOT_SVPopulateDiagViewer(vdp2)==FALSE)
    goto end;

  if (DOT_SVPopulateSequenceViewer(vdp2)==FALSE)
    goto end;

  return(wSequence);	

 end:
   
   ErrPostEx (SEV_WARNING, 0, 0, "%s", "Display functions failed");
   return NULL;

}

/*________________________________________(DOT_UpdateDataRects)_____________

  Purpose : Updates main window rect params in data struct.

____________________________________________________________________
*/
static void DOT_UpdateDataRects (DOTSelDataPtr data, RecT rc, DOTVibDataPtr vdp, Boolean updateSelectedRect)
{
  Int4  width, height, VFrom=vdp->sdp.VFrom, HFrom=vdp->sdp.HFrom;
  Int4  comp, dx, dy, xstart, ystart;


  DOT_AddRectMargins(&rc, vdp);

  /* update the limits of selectable region of parent window */
  width=MIN (rc.right-rc.left, vdp->curr_qlen - HFrom );
  height=MIN (rc.bottom-rc.top, vdp->curr_slen - VFrom );
    
  data->rcP.left=rc.left;
  data->rcP.top=rc.top;
  data->rcP.right=rc.left+width;
  data->rcP.bottom=rc.top+height;

  /* update the size of selected rect */
  if (updateSelectedRect)
    {
      xstart=ABS(data->q_start-vdp->mip->q_start);
      ystart=ABS(data->s_start-vdp->mip->s_start);      
      comp=vdp->comp;
      dx=data->H_pos-HFrom;
      dy=data->V_pos-VFrom;
      width=ABS(data->q_stop-data->q_start)/comp;
      height=ABS(data->s_stop-data->s_start)/comp;

      data->rcS.left=rc.left+((xstart)/comp)/* -HFrom */;
      data->rcS.right=data->rcS.left+width;
      data->rcS.top=rc.top+((ystart)/comp)/* -VFrom */;
      data->rcS.bottom=data->rcS.top+height;
    }
}

static void DOT_InitFeatIndex(DOTFeatIndexPtr fdindex)
{
  Int4  i;

  for(i=1; i<FEATDEF_MAX; i++)
    {
      fdindex[i].show=TRUE;
    }
}
/*________________________________________(DOT_ClickProc)_____________

  Purpose : Click proc for main window - no action.

____________________________________________________________________*/
static void DOT_ClickProc (PaneL p, PoinT pt)
{
  DOTSelDataPtr   data;
  RecT        rc, prc;
  DOTVibDataPtr  vdp;


  ObjectRect(p, &prc);
/*   data = (DOTSelDataPtr) GetObjectExtra (p); */
  vdp = (DOTVibDataPtr)GetObjectExtra(ParentWindow(p));
  if(vdp==NULL) return;
  data=(DOTSelDataPtr)vdp->data;
  if (!data) return;

  /* specify clickable region */
  DOT_UpdateDataRects(data, prc, vdp, FALSE);

  
  rc = data->rcP;

  if (!PtInRect (pt, &rc)) 
    {
      if (pt.y < rc.top) pt.y = rc.top;
      if (pt.y > rc.bottom) pt.y = rc.bottom;
      if (pt.x < rc.left) pt.x = rc.left;
      if (pt.x > rc.right) pt.x = rc.right;
    }

  fstpnt = pt;
  curpnt = pt;
  
  InvertMode();
  if (vdp->selectMode==dot_SEQVIEW)
    DOT_SelectFrameProc (p);
  else
    DOT_SelectLineProc(p);

  if (data->selected)
    data->rm_lastselected=TRUE;
  else
    data->selected=TRUE;

  SetObjectExtra(p, data, NULL);

}

/*________________________________________(DOT_DragProc)_____________

  Purpose : Drag proc for main window - no action.

____________________________________________________________________*/
static void DOT_DragProc (PaneL p, PoinT pt)
{
  DOTSelDataPtr data;
  RecT      rc;
  DOTVibDataPtr vdp;

  vdp=(DOTVibDataPtr)GetObjectExtra(ParentWindow(p));
  if (vdp==NULL) return;
  data = (DOTSelDataPtr) GetObjectExtra (p);
  


  InvertMode();
  if (vdp->selectMode == dot_SEQVIEW)
    DOT_SelectFrameProc(p);
  else
    DOT_SelectLineProc(p);
  
  rc=data->rcP;


  if (!PtInRect (pt, &rc)) 
    {
      if (pt.y < rc.top) pt.y = rc.top;
      if (pt.y > rc.bottom) pt.y = rc.bottom;
      if (pt.x < rc.left) pt.x = rc.left;
      if (pt.x > rc.right) pt.x = rc.right;
    }

  curpnt = pt;
  
  if (vdp->selectMode==dot_SEQVIEW)
    DOT_SelectFrameProc(p);
  else
    DOT_SelectLineProc(p);
    
  
  if (data->selected)
    {
      data->rm_lastselected=TRUE;
    }
  else
    data->selected = TRUE;
  
  SetObjectExtra(p, data, NULL);

}


/*________________________________________(DOT_ReleaseProc)_____________

  Purpose : Release Proc for main window - calls up second window.

____________________________________________________________________*/
static void DOT_ReleaseProc(PaneL p, PoinT pt)
{
  DOTSelDataPtr   data;
  MenU        menu, gmenu;
  Int4        VFrom, HFrom, sFeatscount, qFeatscount;
  DOTVibDataPtr  vdp2=NULL, vdp=NULL;
  DOTMainDataPtr mip1;
  DOTFeatListPtr     flp;
  DOTDiagPtr  PNTR hitlist;
  RecT        rc, rcP;
  Int4        width, height, index, xscale, yscale, scale, Y_pos;
  Int2        dx, dy;
  GrouP       g;
  Char        infoBuf[255];
  


  vdp = (DOTVibDataPtr)GetObjectExtra(ParentWindow(p));
  if (!vdp) return;

  data = (DOTSelDataPtr) vdp->data;

  if (!data->selected) return;

  rc = data->rcP;

  if (!PtInRect (pt, &rc)) 
    {
      if (pt.y < rc.top) pt.y = rc.top;
      if (pt.y > rc.bottom) pt.y = rc.bottom;
      if (pt.x < rc.left) pt.x = rc.left;
      if (pt.x > rc.right) pt.x = rc.right;
    }
  curpnt = pt;
  mip1 = vdp->mip;

  VFrom  = vdp->sdp.VFrom; 
  HFrom = vdp->sdp.HFrom;


  
  if (vdp->selectMode == dot_SEQVIEW)
    {
      DOT_UpdatePt();
      if (vdp->ChildWin==NULL)
        vdp2=(DOTVibDataPtr) MemNew (sizeof(DOTVibData));
      else
        {
          vdp->ChildWin=DOT_ClearLastWindow(vdp->ChildWin, TRUE);
          data->selected=TRUE;
        }

      InvertMode();
      DOT_SelectFrameProc(p);

      dx=HFrom-data->H_pos;
      dy=VFrom-data->V_pos;
      
      /* previous rect coordinates */
      
      data->old_rcS.left=data->rcS.left-dx;
      data->old_rcS.right=data->rcS.right-dx;
      data->old_rcS.top=data->rcS.top-dy;
      data->old_rcS.bottom=data->rcS.bottom-dy;
      
      /* new rect coordinates on parent window */
      
      data->rcS.left = fstpnt.x;
      data->rcS.top = fstpnt.y;
      data->rcS.right = curpnt.x;
      data->rcS.bottom = curpnt.y;
      
      data->H_pos=HFrom;
      data->V_pos=VFrom;
      
  /* map selected region to sequence(world) coordinates 
     plus or minus one to account for errors when rounding off */

  if (mip1->qstrand == Seq_strand_plus)
    {
      data->q_start = MAX((fstpnt.x  - rc.left  + (HFrom-1))*vdp->comp, 0) + mip1->q_start;
      data->q_stop = MIN(((curpnt.x - rc.left +(HFrom+1))*vdp->comp)+mip1->q_start, mip1->q_stop);
    }
  else
    {
      data->q_start = mip1->q_start - MAX((fstpnt.x  - rc.left  + (HFrom-1))*vdp->comp, 0);
      data->q_stop = MAX(mip1->q_start-((curpnt.x - rc.left +(HFrom+1))*vdp->comp), mip1->q_stop);
    }

  if (mip1->sstrand == Seq_strand_plus)
    {
      data->s_start = MAX((fstpnt.y  - rc.top   + (VFrom-1))*vdp->comp, 0) + mip1->s_start;
      data->s_stop = MIN(((curpnt.y - rc.top +(VFrom+1))*vdp->comp)+mip1->s_start, mip1->s_stop);
    }
  else
    {
      data->s_start= mip1->s_start - MAX((fstpnt.y - rc.top + (VFrom-1))*vdp->comp, 0);
      data->s_stop=MAX(mip1->s_start-((curpnt.y - rc.top + (VFrom+1))*vdp->comp), mip1->s_stop);
    }

  data->qlen=ABS(data->q_stop-data->q_start)+1;
  data->slen=ABS(data->s_stop-data->s_start)+1;

  /* create new sequence buffers */
  DOT_InitCInfo(vdp, vdp2, data);
  DOT_FillNewSeqBufs(vdp, vdp2, FALSE);

  /* compute, store and sort hits*/
  if (DOT_BuildHitList(vdp2->mip, TRUE, TRUE)<0)
    {
      data->selected=FALSE;
      data->rm_lastselected =TRUE;
      SetTitle(vdp->Infopanel, vdp->iInfo);
      DOT_UpdateMainPanel(vdp, FALSE);
      Beep();
      return;/* no hits */
    }


    /*reset infopanel*/

  sprintf(infoBuf, "Selected ..   %s (horizontal) [%ld..%ld]   vs.   %s (vertical) [%ld..%ld]", mip1->qname, (long)data->q_start, (long)data->q_stop, mip1->sname, (long)data->s_start, (long)data->s_stop);
  SetTitle(vdp->Infopanel,infoBuf);

  DOT_UpdateMainPanel(vdp, FALSE);

 

      /* create second window */
  vdp->ChildWin=DOT_SVBuildDiagViewer(vdp2);

    }
  else
    {

      InvertMode();
      DOT_SelectLineProc(p);
 
      
      dx=HFrom-data->H_pos;
      dy=VFrom-data->V_pos;
      
      /* previous hair coordinates */
      
      data->old_rcS.left=data->rcS.left-dx;
      data->old_rcS.right=data->rcS.right-dx;
      data->old_rcS.top=data->rcS.top-dy;
      data->old_rcS.bottom=data->rcS.bottom-dy;
      
      /* new hair coordinates on parent window */
      
      data->rcS.left = curpnt.x;
      data->rcS.top = curpnt.y;
      data->rcS.right = curpnt.x;
      data->rcS.bottom = curpnt.y;
      
      data->H_pos=HFrom;
      data->V_pos=VFrom;

      if (mip1->qstrand == Seq_strand_plus)
        {
          data->q_start = MIN(MAX((curpnt.x - rc.left +(HFrom))*vdp->comp, 0)+mip1->q_start, mip1->q_stop);
        }
      else
        {
          data->q_start = MAX(mip1->q_start-((curpnt.x - rc.left +(HFrom))*vdp->comp), mip1->q_stop);
        }
      
      if (mip1->sstrand == Seq_strand_plus)
        {
          data->s_start = MIN(MAX((curpnt.y - rc.top +(VFrom))*vdp->comp, 0)+mip1->s_start, mip1->s_stop);
        }
      else
        {
          data->s_start = MAX(mip1->s_start-((curpnt.y - rc.top + (VFrom))*vdp->comp), mip1->s_stop);
        }

      DOT_UpdateMainPanel(vdp, TRUE);
      /* look for features in selected region */

      if (vdp->ChildWin !=NULL)
        {
          flp=(DOTFeatListPtr)GetObjectExtra(vdp->ChildWin);
        }
      else /* first pass */
        {
          flp=(DOTFeatListPtr)MemNew(sizeof(DOTFeatList));
          flp->data=vdp->data;
          flp->mip=vdp->mip;
          flp->featindex=(DOTFeatIndexPtr) MemNew(sizeof(DOTFeatIndex)*FEATDEF_MAX);
          DOT_InitFeatIndex(flp->featindex);
          vdp->ChildWin=DOT_BuildFeatGUI(flp);
          if (!vdp->ChildWin) return;
          DOT_PopFeatureViewers(flp);
          if (flp->qFeatscount==0 && flp->sFeatscount==0)
            {
              data->selected=FALSE;
              DOT_UpdateMainPanel(vdp, TRUE);
              MemFree(flp);
              ErrPostEx(SEV_WARNING, 0, 0, "no features on bioseqs");
              return;
            }
        }
      
      DOT_UpdateFeatViewer(flp, flp->Query, flp->segQuery, flp->segQName, flp->segQCursor, DOT_QViewerClickProc, flp->query_drp, ((DOTSelDataPtr)flp->data)->q_start, &(flp->vert_Qpos));
      DOT_UpdateFeatViewer(flp, flp->Subject, flp->segSubject, flp->segSName, flp->segSCursor, DOT_SViewerClickProc, flp->subject_drp, ((DOTSelDataPtr)flp->data)->s_start, &(flp->vert_Spos));
      sprintf(infoBuf, "Hairs .. X-axis (%s) [%ld]  vs.  Y-axis (%s) [%ld]", mip1->qname, (long)data->q_start, mip1->sname, (long)data->s_start);
      SetTitle(vdp->Infopanel,infoBuf);
      
    }

  Show(vdp->ChildWin);
}

/*________________________________________(DOT_InitDataStruct)_____________

  Purpose : Initialize data structure.

____________________________________________________________________*/
static void DOT_InitDataStruct (DOTVibDataPtr vdp)
{
  DOTSelDataPtr data;
  RecT       rc;
  
  
  ObjectRect(vdp->panel, &rc);
  InsetRect(&rc, 4, 4);

  data = (DOTSelDataPtr) MemNew (sizeof (DOTSelData));
  data->selected = FALSE;
  data->q_start = 0;
  data->q_stop = 0;
  data->s_start = 0;
  data->s_stop = 0;
  data->V_pos=0;
  data->H_pos=0;
  data->vdp = vdp;
  DOT_UpdateDataRects(data, rc, vdp, FALSE);
  /* initialize document left and top parameters */


  vdp->data = (VoidPtr)data;
  if (vdp->panel != NULL) 
    SetObjectExtra (vdp->panel, data, NULL);
  else
    return;
}
/*________________________________________(DOT_ReducesizeProc)_____________
  
Purpose : Increase compression of main window display.

____________________________________________________________________*/
static void DOT_ReduceSizeProc (IteM i) 
{
  WindoW      w, temport;
  RecT        rcP;
  DOTVibDataPtr  vdp;
  DOTSelDataPtr   data;
  Int4        xscale, yscale, VCurPos, HCurPos;
  BaR         vsb, hsb;

  

  w = (WindoW)ParentWindow(i);
  temport = SavePort(w);

  vdp = (DOTVibDataPtr)GetObjectExtra (w);
  if (vdp==NULL) return;

  data=(DOTSelDataPtr)vdp->data;

  Select(vdp->panel);
  ObjectRect(vdp->panel, &rcP);

  vdp->comp *=2;
 
  DOT_SetCurrSeqlen (vdp);
  DOT_SetScrlVals (vdp);

  vdp->sdp.XScrlPos = (vdp->sdp.XScrlPos*vdp->sdp.UnitX/2)/vdp->sdp.UnitX;
  vdp->sdp.YScrlPos = (vdp->sdp.YScrlPos*vdp->sdp.UnitY/2)/vdp->sdp.UnitY;

  /*current scroll status*/
  vsb = GetSlateVScrollBar ((SlatE) vdp->panel);
  VCurPos=GetBarValue(vsb);
  hsb = GetSlateHScrollBar ((SlatE) vdp->panel);
  HCurPos=GetBarValue(hsb);

  /* update scroll values*/
  DOT_VScrlUpdate(vdp, vsb, VCurPos);
  DOT_HScrlUpdate(vdp, hsb, HCurPos);

  DOT_UpdateDataRects(data, rcP, vdp, TRUE);

  InsetRect(&rcP, -1, -1);
  InvalRect (&rcP);
  RestorePort (temport);
  Update();
  
}

/*________________________________________(DOT_EnlargeSizeProc)_____________
  Purpose : Reduce compression of main window display.

____________________________________________________________________*/
static void DOT_EnlargeSizeProc (IteM i) 
{
  WindoW    w, temport;
  RecT      rcP;
  DOTVibDataPtr vdp;
  DOTSelDataPtr  data;
  Int4     xscale, yscale, VCurPos, HCurPos;
  BaR      vsb, hsb;

  w = (WindoW)ParentWindow(i);
  temport = SavePort(w);

  vdp = (DOTVibDataPtr)GetObjectExtra (w);
  if (vdp==NULL) return;

  data=(DOTSelDataPtr)vdp->data;

  Select(vdp->panel);
  ObjectRect(vdp->panel, &rcP);
  if (vdp->comp >= 2)
    {
      vdp->comp/= 2;
    }

  DOT_SetCurrSeqlen(vdp);
  DOT_SetScrlVals (vdp);
  DOT_UpdateDataRects(data, rcP, vdp, TRUE);

  vdp->sdp.XScrlPos = (vdp->sdp.XScrlPos*vdp->sdp.UnitX*2)/vdp->sdp.UnitX;
  vdp->sdp.YScrlPos = (vdp->sdp.YScrlPos*vdp->sdp.UnitY*2)/vdp->sdp.UnitY;

  /*current scroll status*/
  vsb = GetSlateVScrollBar ((SlatE) vdp->panel);
  VCurPos=GetBarValue(vsb);
  hsb = GetSlateHScrollBar ((SlatE) vdp->panel);
  HCurPos=GetBarValue(hsb);

  /* update scroll values*/
  DOT_VScrlUpdate(vdp, vsb, VCurPos);
  DOT_HScrlUpdate(vdp, hsb, HCurPos);

  InsetRect(&rcP, -1, -1);
  InvalRect (&rcP);
  RestorePort (temport);
  Update();
  
}

static void DOT_ImageSizeHasChanged(DOTVibDataPtr vdp)
{
  BaR      vsb, hsb;
  Int4     HCurPos, VCurPos;

  DOT_SetCurrSeqlen(vdp);
  DOT_SetScrlVals (vdp);

  /*current scroll status*/
  vsb = GetSlateVScrollBar ((SlatE) vdp->panel);
  VCurPos=GetBarValue(vsb);
  hsb = GetSlateHScrollBar ((SlatE) vdp->panel);
  HCurPos=GetBarValue(hsb);

  /* update scroll values*/
  DOT_VScrlUpdate(vdp, vsb, VCurPos);
  DOT_HScrlUpdate(vdp, hsb, HCurPos);

}

/*________________________________________(DOT_OriginalSizeProc)_____________

  Purpose : Resize the main window display to the original size.

____________________________________________________________________*/
static void DOT_OriginalSizeProc (IteM i) 
{
  WindoW    w, temport;
  RecT      rcP;
  DOTVibDataPtr vdp;
  DOTSelDataPtr  data;
  Int4     xscale, yscale;

  w = (WindoW)ParentWindow(i);
  temport = SavePort(w);

  vdp = (DOTVibDataPtr)GetObjectExtra (w);
  if (vdp==NULL) return;

  data=(DOTSelDataPtr)vdp->data;

  Select(vdp->panel);
  ObjectRect(vdp->panel, &rcP);
  
  vdp->comp = vdp->originalcomp;
  
  DOT_ImageSizeHasChanged(vdp);
  DOT_UpdateDataRects(data, rcP, vdp, TRUE);


  InsetRect(&rcP, -1, -1);
  InvalRect (&rcP);
  RestorePort (temport);
  Update();
  
}

/*_______________________________________________(DOT_GetInfoProc)___

  Purpose : 'About' information on software features.

____________________________________________________________________*/

static void DOT_GetInfoProc (IteM i)

{
  WindoW   w;

  w = FixedWindow(-50, -1, -1, -1, "Info", StdCloseWindowProc);
  StaticPrompt (w, "Dot plotter..", 0,  popupMenuHeight, programFont, 'l');

  RealizeWindow(w);
  Show (w);
}

/*_______________________________________________(DOT_FeatAnalysisProc)____
  

  Purpose : Bring up Feature analysis selection window.

____________________________________________________________________*/
static void DOT_FeatAnalysisProc(ButtoN b)
{
  WindoW  w;
  GrouP   g, g2;


  w = FixedWindow(-50, -1, -1, -1, "ACUTE", StdCloseWindowProc);

  g = HiddenGroup(w, 0, 2, NULL);
  g2 = HiddenGroup(g, 0, 2, NULL);
  RadioButton(g2, "Coils");
  RadioButton(g2, "SEG");

  RealizeWindow(w);
  Show (w);
}


/*_______________________________________________(DOT_GetValue)________

 Purpose: Get int value of input TexT

_______________________________________________________________________*/
Int4 DOT_GetValue (TexT t)
{
  Char str[20];
  Int4 val;

  GetTitle (t,  str,  sizeof(str));
  if (StringHasNoText(str))
    {
      ErrPostEx (SEV_WARNING, 0, 0, "%s", "missing parameter(s)");
      return -1;
    }

  val=atoi(str);

  return val;
}


/*_______________________________________________(DOT_DoParams)___

  Purpose : Recalculate dot plot with new parameters.

____________________________________________________________________*/

void DOT_DoParams(ButtoN b)
{
  WindoW w;
  DOTparamsinfoPtr pip;
  DOTVibDataPtr       vdp;
  DOTMainDataPtr      mip;
  RecT             rc;
  Char             wordsize[5], treelimit[20];

  WatchCursor();
  w=(WindoW)ParentWindow(b);
  if (!(pip=(DOTparamsinfoPtr)GetObjectExtra(b))) return;
  if (!(vdp=(DOTVibDataPtr)GetObjectExtra(w))) return;
  mip=vdp->mip;

  DOT_FreeHitsArray(mip->hitlist, mip->index);

  mip->word_size = DOT_GetValue(pip->word_size);
  mip->tree_limit = DOT_GetValue(pip->tree_limit);
  mip->first_pass=TRUE;
  mip->cutoff_score=0;  

  if (pip) MemFree(pip);

  if (DOT_BuildHitList(mip, TRUE, TRUE)<0) 
    {
      Message(MSG_ERROR, "DOT- No hits");
      return; /* no hits */
    }

  DOT_UpdateMainPanel(vdp, TRUE);
  ArrowCursor();
  Remove(w);
}

void DOT_CancelParams(ButtoN b)
{

  WindoW  w;
  DOTparamsinfoPtr pip;

  w=ParentWindow(b);
  pip=(DOTparamsinfoPtr)GetObjectExtra(b);
  if (pip) MemFree(pip);
  Remove(w);
}


void DOT_QuitProg(ButtoN b)
{

  WindoW  w;
  DOTparamsinfoPtr pip;
  DOTVibDataPtr       vdp2;

  w=ParentWindow(b);
  pip=(DOTparamsinfoPtr)GetObjectExtra(b);
  vdp2=(DOTVibDataPtr)GetObjectExtra(w);
  DOT_FreeMainInfo(vdp2->mip);
  if (vdp2->mip) MemFree(vdp2->mip);
  if (vdp2) MemFree(vdp2);
  if (pip) MemFree(pip);
  Remove(w);
  QuitProgram();
}

/*____________________________________________(DOT_StartDOTPLOT)_____________


  Purpose : Starts Dot Plot with user parameters.

____________________________________________________________________*/

void DOT_StartDOTPLOT(ButtoN b)
{
  WindoW w;
  DOTparamsinfoPtr pip;
  DOTVibDataPtr       vdp;
  DOTMainDataPtr      mip;
  RecT             rc;
  Char             wordsize[5], treelimit[20];
  BioseqPtr        qbsp, sbsp;


  WatchCursor();
  w=(WindoW)ParentWindow(b);
  pip=(DOTparamsinfoPtr)GetObjectExtra(b);
  vdp=(DOTVibDataPtr)GetObjectExtra(w);
  mip=(DOTMainDataPtr)MemNew(sizeof(DOTMainData));
/*   DOT_FreeHitsArray(mip->hitlist, mip->index); */

  mip->q_start= DOT_GetValue(pip->xstart);
  mip->q_stop=DOT_GetValue(pip->xstop);
  mip->s_start=DOT_GetValue(pip->ystart);
  mip->s_stop=DOT_GetValue(pip->ystop);
  mip->qlen=mip->q_stop-mip->q_start;
  mip->slen=mip->s_stop-mip->s_start;
  mip->word_size = DOT_GetValue(pip->word_size);
  mip->tree_limit = DOT_GetValue(pip->tree_limit);
  mip->first_pass=TRUE;
  mip->cutoff_score=0; 
  if (pip) MemFree(pip);

  qbsp=mip->qbsp;
  sbsp=mip->sbsp;
  Remove(w);

  DOT_CreateAndStore(mip, qbsp, sbsp, mip->q_start, mip->q_stop, mip->s_start, mip->s_stop, mip->word_size, mip->tree_limit, FALSE);

}



void DOT_SetupParamsWindow(DOTVibDataPtr vdp, Boolean is_startup, Boolean is_nuc)
{
  GrouP g, g1, g2,p1, p2;
  GrouP mainzoomg, zoomg1, zoomg2;
  ButtoN  b1, b2;
  DOTparamsinfoPtr pip;
  DOTMainDataPtr   mip;
  WindoW  w;
  Char    str1[20], str2[20], str3[20], str4[20];


  ArrowCursor();

  pip=(DOTparamsinfoPtr)MemNew(sizeof(DOTparamsinfo));
  mip=vdp->mip;

  w = FixedWindow(-50, -25, -1, -1, "Parameters", StdCloseWindowProc);
  SetObjectExtra(w, (Pointer)vdp, NULL);

  if (is_startup)
    {
      g=HiddenGroup(w, 1, 3, NULL);
      mainzoomg = NormalGroup(g, 1, 2, "Zoom Parameters", systemFont, NULL);
      sprintf(str1, "%ld", (long)mip->q_start);
      sprintf(str2, "%ld", (long)mip->s_start);
      sprintf(str3, "%ld", (long)mip->q_stop);
      sprintf(str4, "%ld", (long)mip->s_stop);
      zoomg1=HiddenGroup(mainzoomg, 3, 1, NULL);
      StaticPrompt(zoomg1, "x-axis:", 0, 0, systemFont, 'l');
      pip->xstart = DialogText(zoomg1, str1, 5, NULL);
      pip->xstop = DialogText(zoomg1, str3, 5, NULL);
      zoomg2=HiddenGroup(mainzoomg, 3, 1, NULL);
      StaticPrompt(zoomg2, "y-axis:", 0, 0, systemFont, 'l');
      pip->ystart = DialogText(zoomg2, str2, 5, NULL);
      pip->ystop = DialogText(zoomg2, str4, 5, NULL);
    }
  else
    g=HiddenGroup(w, 1, 2, NULL);

  g1 = NormalGroup(g, 1, 2, "Parameters", systemFont, NULL);
  sprintf(str1, "%ld", (long)mip->word_size);
  sprintf(str2, "%ld", (long)mip->tree_limit);
  if (is_nuc)
    sprintf(str3, "(4 - 11):");
  else
    sprintf(str3, "(1, 2 or 3):");

  p1=HiddenGroup(g1, 3, 1, NULL);
  StaticPrompt(p1, "Word size ", 0, 0, systemFont, 'l');
  StaticPrompt(p1, str3, 0, 0, systemFont, 'l');
  pip->word_size = DialogText(p1, str1, 5, NULL);

  p2=HiddenGroup(g1, 2, 1, NULL);
  StaticPrompt(p2, "Max number of hits:", 0, 0, systemFont, 'l');
  pip->tree_limit = DialogText(p2, str2, 5, NULL);

  if (!is_startup)
    {
      g2=HiddenGroup(g, 2, 1, NULL);
      b1 =PushButton(g2, "Accept", DOT_DoParams);
      b2 =PushButton(g2, "Cancel", DOT_CancelParams);
      SetObjectExtra(b1, (Pointer)pip, StdCleanupExtraProc);
      SetObjectExtra(b2, (Pointer)pip, StdCleanupExtraProc);
    }
  else
    {
      g2=HiddenGroup(g, 2, 1, NULL);
      b1 =PushButton(g2, "Accept", DOT_StartDOTPLOT);
      b2 =PushButton(g2, "Cancel", DOT_QuitProg);
      SetObjectExtra(b1, (Pointer)pip, StdCleanupExtraProc);
      SetObjectExtra(b2, (Pointer)pip, StdCleanupExtraProc);
    }
  
  Show(w);

}


static void DOT_ParametersProc (IteM i)
{
  DOTVibDataPtr          vdp;
  WindoW              temport, w;
  RecT                rc;
  
  w= (WindoW)ParentWindow(i);
  temport=SavePort(w);
  vdp=(DOTVibDataPtr)GetObjectExtra(w);

  if (vdp==NULL) return;

  DOT_SetupParamsWindow(vdp, FALSE, FALSE);

  Select(vdp->panel);
  ObjectRect(vdp->panel, &rc);
   
  InsetRect(&rc, -1, -1);
  InvalRect (&rc);
  RestorePort (temport);
  Update();

}


/*_______________________________________________(DOT_DoZoom)___

  Purpose : Zoom into a specific region of sequence.

____________________________________________________________________*/

void DOT_DoZoom(ButtoN b)
{
  WindoW w;
  DOTparamsinfoPtr zip;
  DOTVibDataPtr     vdp;
  DOTSelDataPtr  data;
  DOTMainDataPtr mip;
  Int4           viewersize, newLen;
  RecT           rc;
  Char   xstart[20], xstop[20], ystart[20], ystop[20];


  w=(WindoW)ParentWindow(b);
  if (!(zip=(DOTparamsinfoPtr)GetObjectExtra(b))) return;
  if (!(vdp=(DOTVibDataPtr)GetObjectExtra(w))) return;
  mip=vdp->mip;
  data=vdp->data;
  Remove(w);

  DOT_FreeHitsArray(mip->hitlist, mip->index);


  mip->q_start= DOT_GetValue(zip->xstart);
  mip->q_stop=DOT_GetValue(zip->xstop);
  mip->s_start=DOT_GetValue(zip->ystart);
  mip->s_stop=DOT_GetValue(zip->ystop);
  mip->qlen=mip->q_stop-mip->q_start;
  mip->slen=mip->s_stop-mip->s_start;
  mip->first_pass = TRUE;
  mip->cutoff_score=0;  

  DOT_GetSeqs(mip, TRUE);
  if (DOT_BuildHitList(mip, TRUE, TRUE)<0) 
    {
      Message(MSG_ERROR, "DOT- No hits");
      return; /* no hits */
    } 
 
  ObjectRect(vdp->panel, &rc);
  viewersize=MIN(rc.right-rc.left,rc.bottom-rc.top)-vdp->HORZ_MARGIN;
  newLen=MAX(mip->qlen, mip->slen);
  vdp->comp=DOT_Compression(newLen, viewersize);

  DOT_ImageSizeHasChanged(vdp);

  DOT_UpdateMainPanel(vdp, TRUE);

  if (zip) MemFree(zip);
  
}



static void DOT_CancelZoom (ButtoN b)
{
  WindoW  w;
  DOTparamsinfoPtr zip;

  w=ParentWindow(b);
  zip=(DOTparamsinfoPtr)GetObjectExtra(b);
  if (zip) MemFree(zip);
  Remove(w);
}


void DOT_SetupZoomWindow(DOTVibDataPtr vdp)
{
  GrouP          g, zoomg, bottomg, maingroup;
  ButtoN         b1, b2;
  DOTparamsinfoPtr zip;
  DOTMainDataPtr   mip;
  WindoW         w;
  Char           str1[20], str2[20], str3[20], str4[20];

  if (!(zip=(DOTparamsinfoPtr)MemNew(sizeof(DOTparamsinfo)))) return;
  mip=vdp->mip;


  w = FixedWindow(-50, -25, -1, -1, "Zoom", StdCloseWindowProc);
  SetObjectExtra(w, (Pointer)vdp, NULL);
  
  maingroup=HiddenGroup(w, -1, 2, NULL);

  g = NormalGroup(maingroup, 1, 2, "Zoom Parameters", systemFont, NULL);
  sprintf(str1, "%ld", (long)mip->q_start);
  sprintf(str2, "%ld", (long)mip->s_start);
  sprintf(str3, "%ld", (long)mip->q_stop);
  sprintf(str4, "%ld", (long)mip->s_stop);
  zoomg=HiddenGroup(g, 3, 2, NULL);
  StaticPrompt(zoomg, "x-axis:", 0, 0, systemFont, 'l');
  zip->xstart = DialogText(zoomg, str1, 5, NULL);
  zip->xstop = DialogText(zoomg, str3, 5, NULL);
  StaticPrompt(zoomg, "y-axis:", 0, 0, systemFont, 'l');
  zip->ystart = DialogText(zoomg, str2, 5, NULL);
  zip->ystop = DialogText(zoomg, str4, 5, NULL);

  bottomg=HiddenGroup(maingroup, 2, 1, NULL);
  b1 =PushButton (bottomg, "Accept", DOT_DoZoom);
  b2 =PushButton(bottomg, "Cancel", DOT_CancelZoom);
  SetObjectExtra(b1, (Pointer)zip, StdCleanupExtraProc);
  SetObjectExtra(b2, (Pointer)zip, StdCleanupExtraProc);

  Show(w);
}

static void DOT_ZoomProc (IteM i)
{
  DOTVibDataPtr          vdp;
  WindoW              temport, w;
  RecT                rc;
  
  w= (WindoW)ParentWindow(i);
  temport=SavePort(w);
  vdp=(DOTVibDataPtr)GetObjectExtra(w);
  if (vdp==NULL) return;

  DOT_SetupZoomWindow(vdp); 
  
}


/*_______________________________________________(DOT_DoDisplayOpts)___

  Purpose : Update Display

____________________________________________________________________*/


static void DOT_DisplayOptsProc(ChoicE i)

{
  DOTVibDataPtr          vdp;
  WindoW              w;
  Int4                value;
  

  w=(WindoW)ParentWindow(i);
  vdp=(DOTVibDataPtr)GetObjectExtra(w);
  
  value=GetValue(i);
  
  if (value==1)
    {
      vdp->showDotPlot=TRUE;
      vdp->showALIGN=FALSE;
    }
  else if (value==2)
    {
      vdp->showDotPlot=TRUE;
      vdp->showALIGN=TRUE;
    }
  else if (value==3)
    {
      vdp->showDotPlot=FALSE;
      vdp->showALIGN=TRUE;
    }

  DOT_UpdateMainPanel(vdp, TRUE);
  
}
/*____________________________________________(DOT_InitAlnIPtr)_____________


  Purpose : Initialize DOTAlignInfoPtr.

____________________________________________________________________*/
void DOT_InitAlnIPtr (DOTAlignInfoPtr alp)
{

  alp->VERT_MARGIN=50;
  alp->HORZ_MARGIN=80;
  alp->do_scale=TRUE;
  alp->showLabels=FALSE;

  alp->Fh=FontHeight();

}
/*___________________________(AlnMgrGetNextLine)_______


  Purpose : Get starts and stops for a seq_align.

____________________________________________________________________*/
static Boolean AlnMgrGetNextLine(SeqAlignPtr sap, AlnMsgPtr amp1, AlnMsgPtr amp2, Int4Ptr x1, Int4Ptr y1, Int4Ptr x2, Int4Ptr y2, Int4Ptr n) 
{ 
  Boolean  retval; 
  
  if (sap == NULL || amp1 == NULL || amp2 == NULL || x1 == NULL || y1 == NULL || x2 == NULL || y2 == NULL ||n == NULL) 
    return FALSE;
  amp1->row_num = 1; 
  retval = AlnMgrGetNextAlnBit(sap, amp1);
  if (retval == FALSE) 
    return FALSE; 
  if (amp1->gap == 0) 
    { 
      *x1 = amp1->from_b; 
      *x2 = amp1->to_b; 
    } else 
      { 
        if (*x2 == 0) 
          AlnMgrGetNthSeqRangeInSA(sap, 1, x2, NULL); 
        *x1 = *x2;
      } 
  amp2->row_num = 2;
  retval = AlnMgrGetNextAlnBit(sap, amp2);
  if (retval == FALSE) 
    return FALSE; 
  if (amp2->gap == 0) 
    { 
      *y1 = amp2->from_b;
      *y2 = amp2->to_b; 
    } else 
      { 
        if (*y2 == 0) 
          AlnMgrGetNthSeqRangeInSA(sap, 2, y2, NULL); 
        *y1 = *y2;
      } 
  *n++;
  return TRUE; 
} 

/*____________________________________________(DOT_GetSeqAln)___________


  Purpose : Get seq_align coordinates, store as DOTAlnList array.

____________________________________________________________________*/
static Boolean DOT_GetSeqAln (DOTAlignInfoPtr alp)
{ 
  AlnMsgPtr    amp1, amp2;
  OMProcControlPtr ompcp;
  BioseqPtr    bsp1, bsp2;
  SeqAlignPtr  sap, salp;
  SeqIdPtr     sip1, sip2;
  DOTAlnPtr      Aln;
  Boolean      more=FALSE, saved=FALSE;
  Int4         x1, y1, x2, y2, numlines, n, i;
  Int4         xlen, ylen;
  
  
  sap=alp->sap;


  numlines = AlnMgrGetNumSegments(sap); 
  /* query */
  sip1=AlnMgrGetNthSeqIdPtr(sap, 1); 
  sip2=AlnMgrGetNthSeqIdPtr(sap, 2); 
  bsp1=BioseqLockById(sip1); 
  bsp2=BioseqLockById(sip2); 
  alp->xlen=bsp1->length;  
  alp->ylen=bsp2->length; 

  amp1 = AlnMsgNew(); 
  amp2 = AlnMsgNew(); 
  amp1->to_m = amp2->to_m = -1; 
  
  alp->Alnlist=(DOTAlnPtr PNTR) MemNew(sizeof(DOTAlnPtr)*numlines);
  
  n = x1 = x2 = y1 = y2 = i =xlen=ylen=0; 
  if (sap->segtype==SAS_DISC)
    { 
      salp=(SeqAlignPtr)sap->segs; 
       
    } 
  else 
    {
      salp=sap; 
    } 

  while (salp)
    {
      while (more = AlnMgrGetNextLine(salp, amp1, amp2, &x1, &y1, &x2, &y2, &n)) 
        { 
          saved=TRUE;
          Aln=(DOTAlnPtr)MemNew(sizeof(DOTAln));
          Aln->q_start=x1;
          Aln->q_stop=x2;
          Aln->s_start=y1;
          Aln->s_stop=y2;
          Aln->primID=n;
          alp->Alnlist[i]=Aln;
          i++;

        }
      amp1 = AlnMsgReNew(amp1);
      amp2 = AlnMsgReNew(amp2);
      amp1->to_m = amp2->to_m = -1;
      salp=salp->next;
    }
  if (!saved)
    {
      ErrPostEx(SEV_WARNING, 0, 0, "Invalid SeqAlign format");
      return (FALSE);
    }
  alp->index=i;
  AlnMsgFree(amp1);
  AlnMsgFree(amp2);
  SeqAlignFree(sap);
  return (TRUE);
} 

/*____________________________________________(DOT_FillAlignInfoPointer)___________


  Purpose : Get seq_align coordinates, store as DOTAlnList array.

____________________________________________________________________*/
static Boolean DOT_FillAlignInfoPointer (DOTAlignInfoPtr alp)
{ 
  AlnMsgPtr    amp1, amp2;
  OMProcControlPtr ompcp;
  BioseqPtr    bsp1, bsp2;
  SeqAlignPtr  sap=NULL, salp=NULL, sap_tmp=NULL;
  SeqIdPtr     sip1, sip2;
  DOTAlnPtr      Aln;
  Boolean      more=FALSE, saved=FALSE;
  Int4         x1=0, y1=0, x2=0, y2=0, numlines=0, n=0, i;
  Int4         xlen=0, ylen=0, xmin=0, xmax=0, ymin=0, ymax=0;
  Char      text1[42];
  Char      text2[42];


  if (alp==NULL)
  {
    ArrowCursor();
    return(FALSE);
  }
  
  sap=alp->sap;
  if (sap->segtype==SAS_DISC)
    { 
      salp=(SeqAlignPtr)sap->segs; 
       
    } 
  else 
    {
      salp=sap; 
    } 

  alp->sip = AlnMgrGetUniqueSeqs(salp, &n);
   if (n == 0)
     {
       ArrowCursor();
       return(FALSE);
     }
   if (n > 2)
   {
      ErrPostEx(SEV_ERROR, 0, 0, "Can't display diags for more than 2 sequences");
      SeqIdSetFree(alp->sip);
      ArrowCursor();
      return FALSE;
   }
   if (alp->sip->next == NULL)
      sip2 = alp->sip;
   else
      sip2 = alp->sip->next;
   alp->title = (CharPtr)MemNew(50*sizeof(Char));
   SeqIdWrite(alp->sip, alp->title, PRINTID_FASTA_SHORT, 41);
   SeqIdWrite(sip2, text2, PRINTID_FASTA_SHORT, 41);
   StringCat(alp->title, " vs ");
   StringCat(alp->title, text2);
   

   sap_tmp=salp;
   while (sap_tmp){
     numlines +=AlnMgrGetNumSegments(sap_tmp);
     sap_tmp=sap_tmp->next;
   }

  amp1 = AlnMsgNew(); 
  amp2 = AlnMsgNew(); 
  amp1->to_m = amp2->to_m = -1; 
  
  alp->Alnlist=(DOTAlnPtr PNTR) MemNew(sizeof(DOTAlnPtr)*numlines);
  
  n = x1 = x2 = y1 = y2 = i =xlen=ylen=0; 

  while (salp)
    {
      while (more = AlnMgrGetNextLine(salp, amp1, amp2, &x1, &y1, &x2, &y2, &n)) 
        { 
          saved=TRUE;
          Aln=(DOTAlnPtr)MemNew(sizeof(DOTAln));
          
          Aln->q_start=x1;
          Aln->q_stop=x2;
          Aln->s_start=y1;
          Aln->s_stop=y2;
          Aln->primID=i+1;
          alp->Alnlist[i]=Aln;
          i++;

          if (xmin>x1)
            xmin=x1;
          
          if (xmax<x2)
            xmax=x2;
          
          if (ymin>y1)
            ymin=y1;
          
          if (ymax<y2)
            ymax=y2;
        }
      amp1 = AlnMsgReNew(amp1);
      amp2 = AlnMsgReNew(amp2);
      amp1->to_m = amp2->to_m = -1;
      salp=salp->next;
    }

  if (!saved)
    {
       ErrPostEx(SEV_WARNING, 0, 0, "Invalid SeqAlign format");
       return (FALSE);
  }
/*    AlnMgrGetNthSeqRangeInSA(salp, 1, &x1, &x2);  */
   alp->xlen = xmax - xmin + 1; 
   alp->xstart=xmin;
/*    AlnMgrGetNthSeqRangeInSA(salp, 2, &y1, &y2); */
   alp->ylen = ymax - ymin + 1; 
   alp->ystart=ymin;
   alp->index=i;
   AlnMsgFree(amp1);
   AlnMsgFree(amp2);
   SeqAlignFree(sap);
   ArrowCursor();
   return TRUE;
} 

/*_______________________________________________(DOT_DoBlast)___

  Purpose : Calls Blast2Seq.

____________________________________________________________________*/

static FloatHi DOT_get_eval(Int4 exp)
{
  FloatHi eval;
  Int4 i;

  eval = 1;
  for (i=1; i<=exp; i++)
  {
     eval = eval/10;
  }
  return eval;
}


static void DOT_DoBlast (ButtoN b)
{

  BLAST_OptionsBlkPtr options;
  SeqAlignPtr         sap, sap_final;
  DOTVibDataPtr       vdp;
  WindoW              temport, w;
  RecT                rc;
  DOTblastinfoPtr     bip;
  BioseqPtr           bsp1, bsp2;
  Int4                e;
  Int2                i;
  CharPtr             text;
  Boolean             is_local;
   Char                 eval[8];


  WatchCursor();
  w=(WindoW)ParentWindow(b);
  vdp=(DOTVibDataPtr)GetObjectExtra(w);
  if (vdp==NULL)
    goto end;

  Remove(w);

  bip = (DOTblastinfoPtr)GetObjectExtra(b);
  if (bip == NULL || bip->bsp1 == NULL || bip->bsp2 == NULL)
    goto end;
  bsp1 = bip->bsp1;
  bsp2 = bip->bsp2;

  i = GetValue(bip->localorglobal);
  if (i==1)
    is_local=TRUE;
  else
    is_local=FALSE;


  i = GetValue(bip->progname);
  if (i == 1)
    text = StringSave("blastn");
  else if (i == 2)
    text = StringSave("blastp");
  else if (i == 3)
    text = StringSave("blastx");
  else if (i == 4)
    text = StringSave("tblastn");
  else if (i == 5)
    text = StringSave("tblastx");
  else
    goto end;
  
  options = BLASTOptionNew(text, TRUE);
  
  i = GetValue(bip->gapped);
  if (i == 1)
    {
      options->gapped_calculation = TRUE;
    } else if (i == 2)
      {
        options = BLASTOptionNew(text, FALSE);
        options->gapped_calculation = FALSE;
      } else
        goto end;
  
  GetTitle(bip->eval, eval, 14);
  if (eval != NULL)
    {
      e = atoi(eval);
      options->expect_value = DOT_get_eval(e);
    }
  
  GetTitle(bip->wordsize, eval, 5);
  if (eval != NULL)
    options->wordsize = atoi(eval);
  if (GetStatus(bip->maskrep) == TRUE)
    {
      if (GetStatus(bip->masksimple) == TRUE)
        options->filter_string = StringSave("m L;R");
      else
        options->filter_string = StringSave("m R");
    } else if (GetStatus(bip->masksimple) == TRUE)
      options->filter_string = StringSave("m L");
 
   sap = BlastTwoSequences(bsp1, bsp2, text, options);
   MemFree(options);
   BLASTOptionDelete(options);
   ArrowCursor();

   if (vdp->alp){ /* previous alignment*/
     vdp->alp->pict=NULL;
     DOT_ExitAlign(vdp->alp);
     vdp->showALIGN=FALSE;
   }

   if (sap == NULL)
     {
       ErrPostEx(SEV_WARNING, 0, 0, "No BLAST hits found");
       goto end;
     }
   
   if (!AlnMgrIndexSeqAlign(sap)) goto end;
   if (!AlnMgrMakeMultipleByScore(sap)) goto end;
     
   AlnMgrDeleteHidden(sap, FALSE);

 
   sap_final = sap;

   bip=MemFree(bip);
   w=Remove(w);
   

   vdp->showALIGN=TRUE;
   SetValue(vdp->displayOpts2, 2);
   vdp->alp=(DOTAlignInfoPtr)MemNew(sizeof(DOTAlignInfo));
   vdp->alp->sap=sap_final;
   DOT_InitAlnIPtr(vdp->alp);
   if (!DOT_GetSeqAln (vdp->alp))
     {
       vdp->alp->pict=NULL;
       DOT_ExitAlign(vdp->alp);
       vdp->showALIGN=FALSE;
     }

   DOT_UpdateMainPanel(vdp, TRUE);

   end:
   if (text) MemFree(text);
   ArrowCursor();
   return;
}


static void DOT_CancelBlast (ButtoN b)
{
  WindoW  w;
  DOTblastinfoPtr bip;

  w=ParentWindow(b);
  bip=(DOTblastinfoPtr)GetObjectExtra(b);
  if (bip) MemFree(bip);
  Remove(w);
}


void DOT_SetupBlastWindow(DOTVibDataPtr vdp)
{
   DOTblastinfoPtr  bip;
   DOTMainDataPtr   mip;
   ButtoN            b;
   ButtoN            b1;
   GrouP             maingroup, topg, globalg, localg, gapsg, eANDwg;
   GrouP             submitg, maskg, blastg, bottomg, g1, g2, g3;
   CharPtr           title;
   WindoW            w;
   


   if (!(bip = (DOTblastinfoPtr)MemNew(sizeof(DOTblastinfo)))) return;
   mip=vdp->mip;

   bip->bsp1 = mip->qbsp;
   bip->bsp2 = mip->sbsp;
   title = StringCat(mip->qname, " vs ");
   title = StringCat(title, mip->sname);
   w = FixedWindow(-50, -25, -1, -1, title, StdCloseWindowProc);
   SetObjectExtra(w, vdp, NULL);

   maingroup = HiddenGroup(w, 1, 4, NULL);  
   StaticPrompt(maingroup, "Blast2Seqs Options ..", 0, popupMenuHeight, systemFont, 'l');

   topg = HiddenGroup (maingroup, -1, 2, NULL);
   blastg = NormalGroup(topg, 1,1, "Blast Program",  systemFont,NULL);
   bip->progname = HiddenGroup(blastg, 5, 0, NULL);
   RadioButton(bip->progname, "blastn");
   RadioButton(bip->progname, "blastp");
   RadioButton(bip->progname, "blastx");
   RadioButton(bip->progname, "tblastn");
   RadioButton(bip->progname, "tblastx");
   SetValue(bip->progname, 1);

   globalg = NormalGroup(topg,1, 1, "Alignment Type",  systemFont,NULL);
   bip->localorglobal = HiddenGroup(globalg, 2, 1, NULL);
   RadioButton(bip->localorglobal, "Local");
   RadioButton(bip->localorglobal, "Global");
   SetValue(bip->localorglobal, 1);

   bottomg=HiddenGroup(maingroup, 1, 2, NULL);
   localg = NormalGroup(bottomg, 1, 3, "Local Alignment Options", systemFont, NULL);
   g1 = NormalGroup(localg, 1,1, "",  systemFont,NULL);
   gapsg = HiddenGroup(g1, 3, 2, NULL);
   bip->gapped = HiddenGroup(gapsg, 0, 2, NULL);
   RadioButton(bip->gapped, "gapped");
   RadioButton(bip->gapped, "ungapped");
   SetValue(bip->gapped, 1);

   g2 = NormalGroup(localg, 1,1, "",  systemFont,NULL);
   maskg=HiddenGroup(g2, 0, 2, NULL);
   bip->maskrep = CheckBox(maskg, "Mask Repeats", NULL);
   SetStatus(bip->maskrep, FALSE);
   bip->masksimple = CheckBox(maskg, "Mask Simple Sequence", NULL);
   SetStatus(bip->masksimple, TRUE);

   g3 = NormalGroup(localg, 1,1, "",  systemFont,NULL);
   eANDwg = HiddenGroup(g3, 2, 2, NULL);
   StaticPrompt(eANDwg, "E-value:  e-", 0, 0, systemFont, 'l');
   bip->eval = DialogText(eANDwg, "1", 5, NULL);
   StaticPrompt(eANDwg, "wordsize:", 0, 0, systemFont, 'l');
   bip->wordsize = DialogText(eANDwg, "11", 5, NULL);


   submitg=HiddenGroup(bottomg, 2, 0, NULL);
   b = PushButton(submitg, "BLAST", DOT_DoBlast);
   b1 = PushButton(submitg, "Cancel", DOT_CancelBlast);
   SetObjectExtra(b1, (Pointer)bip, StdCleanupExtraProc);
   SetObjectExtra(b, (Pointer)bip, StdCleanupExtraProc);
   Show(w);

   return;


}


static void DOT_Blast2SeqProc (IteM i)

{
  DOTVibDataPtr          vdp;
  WindoW              w;
  
  w= (WindoW)ParentWindow(i);
  vdp=(DOTVibDataPtr)GetObjectExtra(w);

  if (vdp==NULL) return;

  Enable(vdp->displayOpts1);
  DOT_SetupBlastWindow(vdp); 
  
}

/*________________________________________(DOT_FreeAlnList)_____________


  Purpose : Free list of alignments in align list, then free ptr.

____________________________________________________________________*/
static void  DOT_FreeAlnList(DOTAlignInfoPtr alp)
{
  Int4 j, index;
  DOTAlnPtr PNTR alnL, aln;

  if (alp->Alnlist) 
    {
      alnL=alp->Alnlist;
      index=alp->index;
      
      for(j = 0; j < index; j++) 
        {
          aln=alnL[j];
          if (aln) MemFree(aln);
        }
      
      if (alnL) MemFree(alnL);
    }
}
/*____________________________________________(DOT_ExitAlign)___________


  Purpose : Clears alignment structure

____________________________________________________________________*/
static void DOT_ExitAlign(DOTAlignInfoPtr alp){


  if (alp->Alnlist)
    DOT_FreeAlnList(alp);
  
  if (alp->pict)
    DeletePicture(alp->pict);

  if (alp) MemFree(alp);

}
/*_______________________________________________(DOT_QuitMainWindow)___

  Purpose : Quit Main window.

____________________________________________________________________*/

static void DOT_QuitMainWindow (IteM i)

{
  WindoW  w;
  DOTVibDataPtr vdp;
  DOTSelDataPtr  data;

  w=(WindoW)ParentWindow(i);
  vdp=(DOTVibDataPtr)GetObjectExtra(w);  
  data=(DOTSelDataPtr)vdp->data;
  if (data) MemFree(data);

  if (vdp->mip->qname) MemFree(vdp->mip->qname);
  if (vdp->mip->sname) MemFree(vdp->mip->sname);

  DOT_FreeMainInfoPtrEx(vdp->mip);

  if (vdp->alp)
    {
      if (vdp->alp->Alnlist)
        DOT_FreeAlnList(vdp->alp);
      if (vdp->alp) MemFree(vdp->alp);
    }
  if (vdp) MemFree(vdp);
  
  QuitProgram ();
}
/*________________________________________(DOT_GridProc)_____________

  Purpose : Show or Hide the grid on Main window.

____________________________________________________________________*/

static void DOT_GridProc (ChoicE i) 
{
  WindoW      w, temport;
  RecT        rcP;
  DOTVibDataPtr  vdp;
  Int4        value;
  

  w = (WindoW)ParentWindow(i);
  vdp = (DOTVibDataPtr)GetObjectExtra (w);

  value = GetValue(i);
  if (value==1)
    vdp->showGrid = TRUE;
  else
    vdp->showGrid = FALSE;

  DOT_UpdateMainPanel(vdp, TRUE);
  
}

/*________________________________________(DOT_CreateWindowDetails)____

  Purpose : Create buttons on the main window.

*/
static void DOT_CreateWindowDetails (WindoW w, DOTVibDataPtr vdp)
{
  GrouP          g, g1, g2;
  Int4           index, maxzoom, wvals;
  Char           str[15];
  RecT           rc;
  ButtoN         blastButton;
  BaR 			  vsb;
  BaR 			  hsb;
  Int4           Xpixels, Ypixels, pixwidth;
  DOTSelDataPtr      data;
  DOTMainDataPtr     mip;
  PrompT         pr1, pr2;

  
  mip=vdp->mip;

  ObjectRect(w, &rc);
  pixwidth=rc.right-rc.left-10;
  
  g=HiddenGroup(w, 0, 3, NULL); 
  SetGroupSpacing(g, 0, 20);

  g1=HiddenGroup(g, 7, 0, NULL);

  g2=HiddenGroup (g,4, 0, NULL);  
  SetGroupSpacing(g2, 3, 0);
  
  pr2=StaticPrompt (g2, "Threshold-Ramp:", 0, 10, vdp->Fnt, 'l');
  pr1=StaticPrompt (g2, "  20%", 0, 0, vdp->Fnt, 'l');
  vdp->sdp.ScrollBar = ScrollBar (g2, 15, 5, DOT_ChangeMainViewerCutoff);
  pr1=StaticPrompt (g2, "100%", 0, 0, vdp->Fnt, 'l');
  
  sprintf(vdp->iInfo, "%s (horizontal) [%ld..%ld]  vs.   %s (vertical) [%ld..%ld]", mip->qname, (long)mip->q_start, (long)mip->q_stop, mip->sname, (long)mip->s_start, (long)mip->s_stop);
  
  vdp->Infopanel= StaticPrompt (g, vdp->iInfo, pixwidth, popupMenuHeight, vdp->Fnt, 'l');  
  SetObjectExtra(vdp->sdp.ScrollBar, vdp, NULL);

  CorrectBarMax (vdp->sdp.ScrollBar, 80); /* 100% */
  CorrectBarValue (vdp->sdp.ScrollBar, 60);/* 100% */

  /* Main Panel*/
 Xpixels = 600;
 Ypixels = 600;


 vdp->panel = AutonomousPanel (w, Xpixels, Ypixels, DOT_DisplayHits,  DOT_VscrlProc, DOT_HscrlProc , 0, NULL, NULL);
      
 DOT_SetUpWin (w, vdp->panel, vdp);
 DOT_InitDataStruct(vdp);
 SetPanelClick (vdp->panel, DOT_ClickProc, DOT_DragProc, NULL, DOT_ReleaseProc);


}

/*________________________________________(DOT_ResizeMainWindow)_____________

  Purpose : Resize function for main window.

____________________________________________________________________*/
static void DOT_ResizeMainWindow (WindoW w)
{
  RecT        rcP, rcW, rcHsb, rcVsb;
  DOTVibDataPtr  vdp;
  WindoW      temport;
  Int2        height, width, in , gap;
  Int4        VCurPos, HCurPos, vsbWidth, hsbHeight;
  PaneL       p;
  BaR         vsb, hsb;
  DOTSelDataPtr   data;


  vdp = (DOTVibDataPtr) GetObjectExtra (w);
  if (vdp == NULL) return;
  p = vdp->panel;
  
  if (p == NULL) return; 

  ObjectRect (w, &rcW);
  width = rcW.right - rcW.left;
  height = rcW.bottom - rcW.top;
  vsb = GetSlateVScrollBar ((SlatE) p);
  hsb = GetSlateHScrollBar ((SlatE) p);
  
  SafeHide(vdp->Infopanel);
  SafeHide(p);
  Update();

  GetPosition(vsb,&rcVsb);
  GetPosition(hsb,&rcHsb);
  GetPosition (p, &rcP);

  in=2;
  gap=10;
  vsbWidth=rcVsb.right-rcVsb.left;
  hsbHeight=rcHsb.bottom-rcHsb.top;

  
  rcP.right = width - vsbWidth - in;
  rcP.bottom = height - hsbHeight -in;

  SetPosition (p, &rcP);
  AdjustPrnt (p, &rcP, FALSE);
/*   Update(); */

  
  ObjectRect(p ,&rcP);
  temport=SavePort((WindoW)ParentWindow(p));
  Select(p);

  DOT_ComputePanelSize(rcP, vdp, &(vdp->sdp.PgWdth), &(vdp->sdp.PgLen));
  DOT_SetScrlVals (vdp);
  if (vdp->data)
    {
      data=(DOTSelDataPtr)vdp->data;
      DOT_UpdateDataRects(data, rcP, vdp, FALSE);
    }

  /*update scroll status*/
  VCurPos=GetBarValue(vsb);
  HCurPos=GetBarValue(hsb);
  DOT_VScrlUpdate(vdp, vsb, VCurPos);
  DOT_HScrlUpdate(vdp, hsb, HCurPos);
  
/*   InvalRect(&rcP); */
  RestorePort(temport);
  SafeShow(vdp->Infopanel);
  SafeShow(vdp->panel);
  Update();
  
}

/*_______________________________________________(DOT_NewAnalysis)___


  Purpose : Close old windows, begin new analysis.

____________________________________________________________________*/
static void DOT_NewAnalysis (ButtoN b)
{
  WindoW     w;
  DOTVibDataPtr vdp, vdp2;
  DOTSelDataPtr  data;
  
  w = (WindoW)ParentWindow (b);
  
  if (w!= NULL)
    {
      vdp = (DOTVibDataPtr)GetObjectExtra(w);
      data=(DOTSelDataPtr)vdp->data;
      if (data)
        MemFree(data);
      if (vdp)
        {
          if (vdp->mip)
            DOT_FreeMainInfoPtrEx(vdp->mip);
          MemFree(vdp);
        }
      Remove(w);
    }

  if (vdp->ChildWin!= NULL)
    {
      DOT_CloseSequenceWindow(b);
    }

}


/*____________________________________________(DOT_InitFont)_____________


  Purpose : Get font.

____________________________________________________________________*/

static void DOT_InitFont (DOTVibDataPtr vdp)
{
#ifdef WIN_MAC
  vdp->Fnt = ParseFont ("Monaco, 9");
#endif

#ifdef WIN_MSWIN
  vdp->Fnt = ParseFont ("Courier, 7");
#endif

#ifdef WIN_MOTIF
  vdp->Fnt = ParseFont ("fixed, 12");
#endif

	return;

}

/*____________________________________________(DOT_InitVibData)____________


  Purpose : Initialize DOTVibDataPtr for main window.

____________________________________________________________________*/

static void DOT_InitVibData(DOTVibDataPtr vdp)
{


  vdp->sdp.TrampPos=75;
  vdp->showGrid=TRUE;
  vdp->showDotPlot=TRUE;

  DOT_InitFont(vdp);
  vdp->charw=MaxCharWidth();
  vdp->Fh=FontHeight();
  vdp->VERT_MARGIN = 50;
  vdp->HORZ_MARGIN = 80;
  vdp->selectMode = 1;

}
/*________________________________________(DOT_OpenSeqAlignASN)_____________


  Purpose : Open a seqannot file with a seqalign and update the viewer.

____________________________________________________________________*/
static void DOT_OpenSeqAlignASN(IteM i)
{
  DOTVibDataPtr vdp;
  FILE*         afile;
  Char      path [PATH_MAX];
  Pointer         dataptr=NULL;
  Uint2           datatype;
  SeqAnnotPtr     sanp=NULL;
  SeqAlignPtr     sap=NULL;
  DOTAlignInfoPtr alp;


  vdp=(DOTVibDataPtr)GetObjectExtra(ParentWindow(i));
  
  if (GetInputFileName(path, sizeof(path), "", "TEXT"))
    {
      if (path != NULL)
        {
          if (!(afile = FileOpen(path, "r"))){     
            ErrPostEx(SEV_FATAL, 0, 0, "file not found");
            return;
          } 
          dataptr = ReadAsnFastaOrFlatFile (afile, &datatype, NULL, FALSE, FALSE, TRUE, FALSE);
          if (!dataptr){     
            ErrPostEx(SEV_FATAL, 0, 0, "no seqalign found");
              return;
          }
          sanp = (SeqAnnotPtr)(dataptr);
          sap = (SeqAlignPtr)(sanp->data);
          if (!sap){     
            ErrPostEx(SEV_FATAL, 0, 0, "no seqalign found");
              return;
          } 
          AlnMgrIndexSeqAlign(sap);
          Enable((HANDLE)vdp->displayOpts1);
          if (vdp->alp){ /* previous alignment*/
            vdp->alp->pict=NULL;
            DOT_ExitAlign(vdp->alp);
          }
          alp=(DOTAlignInfoPtr)MemNew(sizeof(DOTAlignInfo));
          alp->sap = sap;
          DOT_InitAlnIPtr(alp);
          if (!DOT_FillAlignInfoPointer(alp)){
            vdp->showALIGN=FALSE;
            goto end;
          }
          vdp->alp=alp;
          vdp->showALIGN=TRUE;
          SetValue(vdp->displayOpts2, 2);

        }

    end:
      DOT_UpdateMainPanel(vdp, TRUE);
      fclose(afile);
    }

  return;
}

/*________________________________________(DOT_MakeMainViewer)_____________


  Purpose : Create main window.

____________________________________________________________________*/

extern void DOT_MakeMainViewer (DOTMainDataPtr mip, SeqAlignPtr sap)
{
  WindoW         w;
  Int2           margins;
  MenU		     m1, m2, m3, m4, s1, s2;
  IteM           i;
  Int4           n;
  ChoicE         ch;
  DOTVibDataPtr  vdp;
  DOTAlignInfoPtr alp;
  Char           title[50];

  vdp=(DOTVibDataPtr)MemNew(sizeof(DOTVibData));

  /* cgeck this on Patrick's code */
  margins = 4*stdCharWidth;
  sprintf(title, "%s", mip->qname);
  StringCat(title, "  vs  ");
  StringCat(title, mip->sname);

#ifdef WIN_MAC
  DOT_SetupMenus ();
#endif
  w = DocumentWindow (margins, margins, 800, 800, title, StdCloseWindowProc,   DOT_ResizeMainWindow);

#ifndef WIN_MAC
  DOT_SetupMenus ();
#endif

  SetObjectExtra(w, (Pointer)vdp, NULL);
  m1 = PulldownMenu (w, "File"); 
/*   CommandItem (m1, "New analysis", DOT_NewAnalysis);  */
  i=CommandItem(m1, "Open ASN SeqAlign...", DOT_OpenSeqAlignASN);
  CommandItem (m1, "Quit",DOT_QuitMainWindow); 
  
  m2 = PulldownMenu (w, "Options"); 
  vdp->displayOpts1=SubMenu(m2,"Display");
  if (!sap) Disable((HANDLE)vdp->displayOpts1);
  vdp->displayOpts2=ChoiceGroup (vdp->displayOpts1, DOT_DisplayOptsProc);
  ChoiceItem (vdp->displayOpts2, "Dots ONLY");
  ChoiceItem (vdp->displayOpts2, "Dots & Aligns"); 
  ChoiceItem (vdp->displayOpts2, "Aligns ONLY"); 
  SetValue(vdp->displayOpts2,2); 
  SeparatorItem(m2);
  s1=SubMenu(m2,"Grid");
  ch=ChoiceGroup (s1, DOT_GridProc);
  ChoiceItem (ch, "show"); 
  ChoiceItem (ch, "hide");
  SeparatorItem(m2);
  s2=SubMenu(m2,"Select Mode");
  ch=ChoiceGroup (s2, DOT_ModeProc);
  ChoiceItem (ch, "Sequence"); 
  ChoiceItem (ch, "Features");
  SetValue(ch, 1);
  SeparatorItem(m2);
  i = CommandItem (m2, "Parameters ..", DOT_ParametersProc); 
  SetObjectExtra(i, vdp, NULL);
  m3 = PulldownMenu(w, "Resize");
  i = CommandItem (m3, "Reduce", DOT_ReduceSizeProc); 
  i = CommandItem (m3, "Enlarge", DOT_EnlargeSizeProc); 
  i = CommandItem (m3, "Original", DOT_OriginalSizeProc);
  i = CommandItem (m3, "Zoom into region..", DOT_ZoomProc); 
  SetObjectExtra(i, vdp, NULL);
  m4 = PulldownMenu (w, "Analysis"); 
  CommandItem(m4, "Blast2Seq ..", DOT_Blast2SeqProc);
  SetValue(ch, 1); 
  

  vdp->MainWin = w;
  vdp->mip= mip;
  if (sap){
    alp=(DOTAlignInfoPtr)MemNew(sizeof(DOTAlignInfo));
    alp->sap = sap;
    DOT_InitAlnIPtr(alp);
    if (!DOT_FillAlignInfoPointer(alp)){
      DOT_FreeMainInfoPtrEx(vdp->mip);
      QuitProgram();
    }
    vdp->alp=alp;
    vdp->showALIGN=TRUE;
    SetValue(vdp->displayOpts2, 2);

  }
  DOT_InitVibData(vdp);
  DOT_CreateWindowDetails(w, vdp);
  RealizeWindow(w);
  Show (w);
  ProcessEvents();
}



/*________________________________________(DOT_AlnDisplayDiags)_____________
 
 Purpose : Draw function for seq alignments.

____________________________________________________________________*/

Boolean LIBCALLBACK DOT_AlnDisplayDiags(DOTAlignInfoPtr alp)
{
  RecT               rcP;
  DOTAlnPtr    PNTR    alnL;
  DOTAlnPtr            aln;
  DOTVibDataPtr        vdp;
  VieweR             v;
  SegmenT            seg1, seg2, seg3;
  Char               buffer[50];
  Int4               index, x, y, x2, y2, i, length;
  Int4               q_start, q_stop, s_start, s_stop;
  Int4               width, height, primID;
  Int4               r_width, r_ht;


  seg1=CreateSegment(alp->pict, 1, 0); /* diags */
  seg2=CreateSegment(alp->pict, 2, 0); /* axis */
  seg3=CreateSegment(alp->pict, 3, 0); /* diag labels*/
  v=alp->v;  alp->seg1=seg1; /* diags */

  GetPosition(v, &rcP); 
  InsetRect(&rcP, 4, 4);

  width = rcP.right-rcP.left;
  height = rcP.bottom-rcP.top-2*(alp->HORZ_MARGIN*alp->scaleValue);

  /* Rect Parameters */
  rcP.left+=2*(alp->HORZ_MARGIN*alp->scaleValue);

  index = alp->index;
  alnL = alp->Alnlist;
  
  
  DOT_DrawXAxis(seg2, rcP, height, alp->xstart, alp->xlen+alp->xstart, alp->scaleValue);
  DOT_DrawYAxis(seg2, rcP, height, alp->ystart, alp->ylen+alp->ystart, alp->scaleValue, alp->Fh);

  AddAttribute(seg1, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);

  for (i = 0; i<index ; i++)
       {  

         aln=alnL[i];
         length= aln->q_stop - aln->q_start;
         q_start = rcP.left + aln->q_start;
         s_start = height - aln->s_start;
         
         q_stop = q_start + length;
         s_stop = s_start - length;

         if ( alp->showLabels )
           {
             x=aln->q_start;
             y=aln->s_start;
             sprintf(buffer, "(%ld,%ld,%ld)", (long)x, (long)y, (long)length);
             AddLabel(seg3, q_start+7, s_start, buffer, SMALL_TEXT, 0, UPPER_RIGHT, 0);
           }


         AddLine(seg1, q_start, s_start, q_stop, s_stop, FALSE,(Uint2) aln->primID);
       }

  return TRUE;
}


/*________________________________________(DOT_AlnCalculateScaling)________


  Purpose : Calculate approximate scale for alignment display.

____________________________________________________________________*/
static void DOT_AlnCalculateScaling (DOTAlignInfoPtr alp)

{
  RecT   r, world;
  Int4   index, r_hgt, r_wdt, w_hgt, w_wdt, scale;
  double f1, f2;

  

  w_hgt=alp->xlen+(alp->xlen*0.15);
  w_wdt=alp->ylen+(alp->ylen*0.15);

  GetPosition(alp->v, &r);
  r_hgt=r.bottom-r.top;
  r_wdt=r.right-r.left;

  f1=(float)w_hgt/r_hgt;
  f2=(float)w_wdt/r_wdt;

  scale=MAX(ceil(f1), ceil(f2));

  for (index=1; index<MAXZOOMSCALEVAL; index++) 
    {
      if (zoomScaleVal [index]>= scale)
        {
          alp->scaleValue=zoomScaleVal[index];
          alp->scaleIndex=index;
          return;
        }
    }

  alp->scaleValue=zoomScaleVal[MAXZOOMSCALEVAL-1];
  alp->scaleIndex=MAXZOOMSCALEVAL-1;

}


static Boolean DOT_DeselectPrim (SegmenT seg, PrimitivE prim, Uint2 segID,Uint2 primID, Uint2 primCt, VoidPtr userdata)

{

  Int1  highlight;
  VieweR  v;

  v = (VieweR)userdata;
  GetPrimDrawAttribute (prim, NULL, NULL, NULL, NULL, NULL, &highlight);
  if (highlight != PLAIN_PRIMITIVE) 
    HighlightPrimitive (v, seg, prim, PLAIN_PRIMITIVE);

  return TRUE;
}

static DOTAlnPtr DOT_FindAlignment(DOTAlignInfoPtr alp, Uint2 primID)
{
  DOTAlnPtr aln;
  DOTAlnPtr PNTR alnL;
  Int4      index, j;

  if (!alp->Alnlist) return NULL;

  alnL=alp->Alnlist;
  index=alp->index;
      
  for(j = 0; j < index; j++) 
    {
      aln=alnL[j];
      if (aln->primID==primID && primID==j+1) return(aln);
    }
      
  return NULL;
}



static void DOT_AlignClickProc(VieweR v, SegmenT seg, PoinT pt)
{
  DOTAlignInfoPtr  alp;
  PrimitivE  prim;
  Uint2       primID=0;
  Int1        highlight;
  DOTAlnPtr   Aln=NULL;
  Char        title[100], str[50];

  alp=(DOTAlignInfoPtr)GetObjectExtra((WindoW)ParentWindow(v));
/*   ExploreSegment (seg, (Pointer)&v, DOT_DeselectPrim);  */
  if (FindSegPrim(v, pt, NULL, NULL, &prim)){
     GetPrimitiveIDs(prim, NULL, NULL, NULL, &primID);
     if (prim_prev!=NULL){
       GetPrimDrawAttribute (prim_prev, NULL, NULL, NULL, NULL, NULL, &highlight);
       if (highlight != PLAIN_PRIMITIVE) 
         HighlightPrimitive (v, seg, prim_prev, PLAIN_PRIMITIVE);
     }
     if (primID!=0)
        {
          Aln=DOT_FindAlignment(alp, primID);
          if (!Aln)
            goto end;
          else{
            MemSet((Pointer)title, '\0', sizeof(title));
            StringCat(title, " X-axis[");
            sprintf(str, "%ld", (long)Aln->q_start);
            StringCat(title, str);
            StringCat(title, "-");
            sprintf(str, "%ld", (long)Aln->q_stop);
            StringCat(title, str);
            StringCat(title, "]     Y-axis[");
            sprintf(str, "%ld", (long)Aln->s_start);
            StringCat(title, str);
            StringCat(title, "-");
            sprintf(str, "%ld", (long)Aln->s_stop);
            StringCat(title, str);
            StringCat(title, "]");
            SetTitle(alp->Infopanel, title);
            GetPrimDrawAttribute (prim, NULL, NULL, NULL, NULL, NULL, &highlight);
            if (highlight == PLAIN_PRIMITIVE) 
              HighlightPrimitive (v, seg, prim, FRAME_PRIMITIVE);

            prim_prev=prim;
          }
        }
      else
        {
        end:
          SetTitle(alp->Infopanel, "");
          prim_prev=NULL;
        }
  }
      
}
/*________________________________________(DOT_PopulateAlnViewer)_____________


  Purpose : Calls draw function for alignment display.

____________________________________________________________________*/
static Boolean DOT_PopulateAlnViewer(DOTAlignInfoPtr alp)
{
  Int4       scale, index;
  Char       str[16];
  DOTSelDataPtr  data;
  BaR        vsb;

  
  ResetViewer(alp->v);
  alp->pict=DeletePicture(alp->pict);
  alp->pict=CreatePicture();
   
  if (DOT_AlnDisplayDiags(alp)==FALSE)
    return FALSE;

  /* initialize scale value*/
  
  if (alp->do_scale) 
    {

      for (index=1; index<MAXZOOMSCALEVAL; index++) 
        {
          sprintf (str, "%ld", (long) (zoomScaleVal [index]));
          PopupItem (alp->scale, str);
        }
      SetValue (alp->scale, alp->scaleIndex);
      alp->do_scale = FALSE;

    }

  SafeShow(alp->scale);

  AttachPicture(alp->v, alp->pict, INT4_MIN, INT4_MAX, LOWER_RIGHT,  alp->scaleValue , alp->scaleValue , NULL);
  SetViewerProcs (alp->v, DOT_AlignClickProc, NULL, NULL, NULL);
  ArrowCursor();
  return TRUE;
}

/*________________________________________(DOT_AlnChangeScale)_____________


  Purpose : Change scale function for alignment function.

____________________________________________________________________*/
static void DOT_AlnChangeScale (PopuP p)

{
  DOTAlignInfoPtr   alp;
  Int4       index;

  alp = (DOTAlignInfoPtr) GetObjectExtra (p);
  if (alp != NULL) 
    {
      index = GetValue (alp->scale);
      if (index <= MAXZOOMSCALEVAL && index > 0) 
        {
          alp->scaleValue = zoomScaleVal [index];
        } 
      else 
        {
          alp->scaleValue = 1;
        }

      DOT_PopulateAlnViewer (alp);
    }
}


/*________________________________________(DOT_CloseAlnWindow)_____________


  Purpose : Close function alignment window.

____________________________________________________________________*/
static void  DOT_QuitAlnWindow (IteM i)
{
  WindoW  alnw;
  DOTAlignInfoPtr alp;

  alnw=ParentWindow(i);

  alp=(DOTAlignInfoPtr)GetObjectExtra(alnw);
  DOT_ExitAlign(alp);
  alnw=Remove (alnw);
  QuitProgram();
}
/*________________________________________(DOT_ResizeAlnViewer)_____________


  Purpose : Resize function for alignment window.

____________________________________________________________________*/

static void DOT_ResizeAlnWindow(WindoW w)
{

  DOTAlignInfoPtr alp;
  RecT     rcDlg,rcV, rcVsb, rcHsb;
  Int2     height, width, gap, vsbWidth, hsbHeight, lmargin;
  BaR      vsb, hsb;
  WindoW   temport;
  
  

  alp=(DOTAlignInfoPtr)GetObjectExtra(w);
  temport=SavePort(w);
  ObjectRect(w,&rcDlg);
  width= rcDlg.right-rcDlg.left;
  height= rcDlg.bottom-rcDlg.top;
  
  SafeHide(alp->v);
  SafeHide(alp->Infopanel);
  Update();

  vsb = GetSlateVScrollBar ((SlatE) alp->v);
  hsb = GetSlateHScrollBar ((SlatE) alp->v);

  GetPosition(alp->v,&rcV);
  GetPosition(vsb,&rcVsb);
  GetPosition(hsb,&rcHsb);
  
  gap=2;
  lmargin=10;
  
  vsbWidth=rcVsb.right-rcVsb.left;
  hsbHeight=rcHsb.bottom-rcHsb.top;
  
  /*new sizes for the viewers*/	
  rcV.left=lmargin;
  rcV.right=width-gap-vsbWidth-rcV.left;
  rcV.bottom=height-gap-hsbHeight;
  
  
  /*set the new sizes*/
  SetPosition(alp->v,&rcV);
  AdjustPrnt (alp->v, &rcV, FALSE);

  if (Visible (alp->v) && AllParentsVisible (alp->v)) 
    {
      ViewerWasResized(alp->v);
    }

  if (alp->do_scale!=TRUE){
    DOT_PopulateAlnViewer (alp);
  }
  SafeShow(alp->v);
  SafeShow(alp->Infopanel);
  ArrowCursor();
  RestorePort(temport);
  Update();

}

/*________________________________________(DOT_BuildAlnViewer)_____________


  Purpose : Creates alignment window.

____________________________________________________________________*/

Boolean DOT_BuildAlnViewer(DOTAlignInfoPtr alp)
{
  WindoW    alnw;
  GrouP     s, s2;
  PrompT    pr1; 
  Int2      Margins; 
  MenU      menu; 
  BaR       hsb;
  Char      zoombuf[]={"Decrease scale to zoom in .."};
  Char      coordbuf[]={"                                       "};
  Char      title1[50], title2[50], text[100];



  if (!alp) return(FALSE);
  MemSet((Pointer)text, '\0', sizeof(text));
  StringCat(text, " X-axis:");
  SeqIdWrite(alp->sip, title1, PRINTID_FASTA_SHORT, 41);
  StringCat(text, title1);
  StringCat(text, "     Y-axis:");
  SeqIdWrite(alp->sip->next, title2, PRINTID_FASTA_SHORT, 41);
  StringCat(text, title2);
  
	Margins=10*stdCharWidth;
	alnw = DocumentWindow(Margins,Margins ,-10, -10, alp->title, StdCloseWindowProc, DOT_ResizeAlnWindow);
   menu = PulldownMenu (alnw, "Display");
   CommandItem (menu, "Quit", DOT_QuitAlnWindow); 

	if (!alnw) return(FALSE);
 
   s = HiddenGroup (alnw,0, 4, NULL);

   StaticPrompt (s, text,StringWidth(text)+10 , popupMenuHeight, programFont, 'l');
   alp->Infopanel= StaticPrompt (s, coordbuf,StringWidth(coordbuf)+10 , 0, programFont, 'l');

   s2 = HiddenGroup (s,3 ,0, NULL);
   SetGroupSpacing(s2, 3,0);
   pr1 = StaticPrompt (s2, zoombuf, StringWidth(zoombuf)+10 , popupMenuHeight, programFont, 'l');
#ifdef WIN_MAC
   alp->scale = PopupList (s2, TRUE, DOT_AlnChangeScale);
#endif

#ifndef WIN_MAC
   alp->scale = PopupList (s2, FALSE, DOT_AlnChangeScale);
#endif

   SetObjectExtra (alp->scale, alp, NULL);
 	/* viewer */
   alp->w=alnw;
	alp->v=CreateViewer(s,600,500,TRUE,TRUE);
	alp->pict=CreatePicture();
   
   SetObjectExtra (alnw, (Pointer) alp, NULL);
   RealizeWindow(alnw);
   DOT_ResizeAlnWindow(alnw);

   /* calculate scale values */
   DOT_AlnCalculateScaling(alp);

	/*populate the viewer : seq aligns*/
   if (DOT_PopulateAlnViewer(alp)==FALSE)
     goto end;

   return TRUE;	

 end:
   
   ErrPostEx (SEV_WARNING, 0, 0, "%s", "BuildAlnViewer failed");
   return FALSE;

}



/*____________________________________________(Start_AlignPlot)_____________


  Purpose : Start function for Alignment Display.

____________________________________________________________________*/
Boolean Start_AlignPlot (DOTAlignInfoPtr alp)
{
  
  DOT_InitAlnIPtr(alp);
  if (!DOT_GetSeqAln (alp))
     {
       alp->pict=NULL;
       DOT_ExitAlign(alp);
       return FALSE;
     }
  if (!DOT_BuildAlnViewer(alp))
    return FALSE;
  return TRUE;
}
/*____________________________________________(Run_DiagPlot)_____________


  Purpose : Run alignment display given seqalignpointer

____________________________________________________________________*/
extern void DOT_AlignPlotGivenSeqAlign(SeqAlignPtr sap)
{
  DOTAlignInfoPtr alp;

  alp=(DOTAlignInfoPtr)MemNew(sizeof(DOTAlignInfo));
  alp->sap = sap;
  DOT_InitAlnIPtr(alp);
  if (!DOT_FillAlignInfoPointer(alp)){
    MemFree(alp);
    QuitProgram();
  }
  if (!DOT_BuildAlnViewer(alp))
    return;
  Show(alp->w);
  ProcessEvents();
  return;
}

/*____________________________________________(Run_DiagPlot)_____________


  Purpose : Registered function for Alignment Display.

____________________________________________________________________*/
Int2 LIBCALLBACK Run_DiagPlot (Pointer data)
{
  DOTAlignInfoPtr          alp;
   OMProcControlPtr  ompcp;
   OMUserDataPtr     omudp;
   Uint2             userkey;
   SeqAlignPtr       sap;

   sap = NULL;
   ompcp = (OMProcControlPtr)data;
   if (ompcp == NULL || ompcp->proc == NULL)
   {
       ErrPostEx(SEV_ERROR, 0, 0, "Align- Run_DiagPlot with NULL ompcp");
       return OM_MSG_RET_ERROR;
   }
   if (ompcp->input_itemtype != OBJ_SEQALIGN)
   {
      ErrPostEx(SEV_ERROR, 0, 0, "Align- Must start Run_DiagPlot with a seqalign");
      return OM_MSG_RET_ERROR;
   }
   sap = ompcp->input_data;
   if (sap == NULL)
   {
      ErrPostEx(SEV_ERROR, 0, 0, "Align- NULL seqalign in Run_DiagPlot");
      return OM_MSG_RET_ERROR;
   }
   if (sap->saip == NULL)
   {
      if (!AlnMgrIndexSeqAlign(sap))
      {
         ErrPostEx(SEV_ERROR, 0, 0, "Align- Can't index input seqalign");
         return OM_MSG_RET_ERROR;
      }
   }
   userkey = OMGetNextUserKey();
   omudp = ObjMgrAddUserData(sap->idx.entityID, ompcp->proc->procid, OMPROC_VIEW, userkey);
   alp=(DOTAlignInfoPtr)MemNew(sizeof(DOTAlignInfo));
   omudp->userdata.ptrvalue = (Pointer)alp;
   alp->sap = sap;
   if (!Start_AlignPlot(alp))
      return OM_MSG_RET_ERROR;
   Show(alp->w);
   return OM_MSG_RET_OK;

}


/*____________________________________________(DOT_StartDotPlot)_____________


  Purpose : Start up Dot Plot, called by Registered function Run_DotPlot.

____________________________________________________________________*/
Boolean DOT_StartDotPlot (BioseqPtr qbsp, BioseqPtr sbsp)
{

  DOTMainDataPtr mip;
 


  if (!(mip = (DOTMainDataPtr) MemNew (sizeof(DOTMainData)))) return FALSE;

  DOT_CreateAndStore(mip, qbsp, sbsp, 0, qbsp->length, 0, sbsp->length, 0, 0, TRUE);

  DOT_MakeMainViewer(mip, NULL);

  return TRUE;

}


/*____________________________________________(DOT_StartDotPlotWithParams)_____________


  Purpose : Start up Dot Plot with a user input parameters window.

____________________________________________________________________*/



static void DOT_StartDotPlotWithParams (DOTVibDataPtr vdp, BioseqPtr qbsp, BioseqPtr sbsp)
{

  if (vdp==NULL) return;

  vdp->mip=DOT_InitMainInfo(vdp->mip, qbsp, sbsp, 8, 10000,0, qbsp->length-1, 0, sbsp->length-1);

  DOT_SetupParamsWindow(vdp, TRUE, vdp->mip->is_na);
/*   ProcessEvents(); */

}


/*____________________________________________(Run_DotPlot)_____________


  Purpose : Registered function for Dot Plot display.

____________________________________________________________________*/
Int2 LIBCALLBACK Run_DotPlot (Pointer data)
{
/*   DOTMainDataPtr       mip; */
  BioseqPtr         qbsp, sbsp;
  SeqEntryPtr       sep;
  SeqAlignPtr         sap;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  SelStructPtr      ssp;
  Uint2             userkey;


  ompcp = (OMProcControlPtr) data;

  /* checks */

  if (ompcp == NULL || ompcp->proc == NULL) 
    {
      ErrPostEx(SEV_ERROR, 0, 0, "DOT- Run_DotPlot with NULL ompcp");
      return OM_MSG_RET_ERROR;
    }
  
  if (ompcp->input_itemtype == OBJ_BIOSEQ)
    {
      
      /* get data */
      ssp = ObjMgrGetSelected ();

      if (ssp == NULL || ssp->itemtype != OBJ_BIOSEQ || ssp->next == NULL || ssp->next->itemtype != OBJ_BIOSEQ)
        {

          ErrPostEx(SEV_ERROR, 0, 0, "DOT- Must select two bioseqs");
          return OM_MSG_RET_ERROR;
        }

      userkey = OMGetNextUserKey(); 
      omudp = ObjMgrAddUserData(0, ompcp->proc->procid, OMPROC_EDIT, userkey);

  /* get the bsp from the passed in data */
      
      if (!(qbsp = (BioseqPtr)GetPointerForIDs (ssp->entityID, ssp->itemID, OBJ_BIOSEQ)) || (sbsp = (BioseqPtr)GetPointerForIDs (ssp->next->entityID, ssp->next->itemID, OBJ_BIOSEQ))) 
        {
          Message(MSG_ERROR, "DOT- Error in retrieving bioseqs");
          return OM_MSG_RET_ERROR;
        }
      
      DOT_StartDotPlot (qbsp, sbsp);
      /* you can also call DOT_StartDotPlotWithParams(DOTMainDataPtr mip, BioseqPtr qbsp, BioseqPtr sbsp); */
      
    }
  else
    {
      ErrPostEx(SEV_ERROR, 0, 0, "DOT- Must select item");
      return OM_MSG_RET_ERROR;
    }

  return OM_MSG_RET_OK;
}
