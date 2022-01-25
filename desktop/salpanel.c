/*   salpanel.c
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
* File Name:  salpanel.c
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.31 $
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
#include <saledit.h>
#include <salpanel.h>
#include <salutil.h>
#include <salfiles.h>
#include <salstruc.h>
#include <fstyle.h>

#define OBJ_VIRT 254

static Uint1 rectSym [] = {
  0xFE, 0x82, 0x82, 0x82, 0x82, 0x82, 0xFE, 0x00
};
static Uint1 diamondSym [] = {
  0x10, 0x28, 0x44, 0x82, 0x44, 0x28, 0x10, 0x00
};
static Uint1 ovalSym [] = {
  0x38, 0x44, 0x82, 0x82, 0x82, 0x44, 0x38, 0x00
};
static Uint1 leftTriSym [] = {
  0x06, 0x1A, 0x62, 0x82, 0x62, 0x1A, 0x06, 0x00
};
static Uint1 rightTriSym [] = {
  0xC0, 0xB0, 0x8C, 0x82, 0x8C, 0xB0, 0xC0, 0x00
};
static Uint1 upTriSym [] = {
  0x10, 0x28, 0x28, 0x44, 0x44, 0x82, 0xFE, 0x00
};
static Uint1 downTriSym [] = {
  0xFE, 0x82, 0x44, 0x44, 0x28, 0x28, 0x10, 0x00
};
static Uint1 rectFillSym [] = {
  0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0x00
};
static Uint1 diamondFillSym [] = {
  0x10, 0x38, 0x7C, 0xFE, 0x7C, 0x38, 0x10, 0x00
};
static Uint1 ovalFillSym [] = {
  0x38, 0x7C, 0xFE, 0xFE, 0xFE, 0x7C, 0x38, 0x00
};
static Uint1 leftTriFillSym [] = {
  0x06, 0x1E, 0x7E, 0xFE, 0x7E, 0x1E, 0x06, 0x00
};
static Uint1 rightTriFillSym [] = {
  0xC0, 0xF0, 0xFC, 0xFE, 0xFC, 0xF0, 0xC0, 0x00
};
static Uint1 upTriFillSym [] = {
  0x10, 0x38, 0x38, 0x7C, 0x7C, 0xFE, 0xFE, 0x00
};
static Uint1 downTriFillSym [] = {
  0xFE, 0xFE, 0x7C, 0x7C, 0x38, 0x38, 0x10, 0x00
};
static Uint1 rightOvalSym [] = {
  0x18, 0x14, 0x12, 0x12, 0x12, 0x14, 0x18, 0x00
};
static Uint1 leftOvalSym [] = {
  0x38, 0x44, 0x82, 0x82, 0x82, 0x44, 0x38, 0x00
};
static Uint1 rightOvalFillSym [] = {
  0x18, 0x1C, 0x1E, 0x1E, 0x1E, 0x1C, 0x18, 0x00
};
static Uint1 leftOvalFillSym [] = {
  0x30, 0x50, 0x90, 0x90, 0x90, 0x50, 0x30, 0x00
};

/*######################################################################
#
#       functions for setting up the color for different object
#
######################################################################*/
#define RGB_B(x) (Uint1)((x)&255);
#define RGB_G(x) (Uint1)(((x)>>8)&255);
#define RGB_R(x) (Uint1)(((x)>>16)&255);

static Boolean convert_color(Int4 val, Uint1Ptr color)
{

        if(val<0 || color == NULL)
                return FALSE;
        color[0] = RGB_R(val);
        color[1] = RGB_G(val);
        color[2] = RGB_B(val);
        return TRUE;
}

static Uint4 getcolor_fromstyles (Uint2 itemsubtype)
{
  Int4 c_val;
  Uint1 color [3];

  c_val = GetMuskCParam (itemsubtype, MSM_SEGMENT, MSM_COLOR);
  if ( convert_color (c_val, color) )
     return GetColorRGB (color[0], color[1], color[2]);
  return GetColorRGB (0, 0, 0);
}

static Uint4 getcolorforvirtfeat (Uint2 item)
{
  return GetColorRGB ((Uint1)((item % 3) * 40), 
                      (Uint1)(255-(item % 3) * 40), 
                      (Uint1)((item % 3) * 40)); 
}

/**********************************************
***   GetAlignEditData
**********************************************/
extern WindoW getwindow_frompanel (PaneL pnl)
{
/*
#ifdef WIN_MAC
  w = FrontWindow ();
#else
  w = ParentWindow (pnl);
#endif
*/
  return ParentWindow (pnl);
}

extern PaneL GetPanelFromWindow (WindoW w)
{
  SeqEditViewFormPtr wdp;

  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  if (wdp == NULL) { 
    return NULL; 
  }
  if ( wdp->pnl == NULL ) { 
    return NULL;
  }
  return wdp->pnl;
}

extern EditAlignDataPtr GetAlignEditData (WindoW w)
{
  PaneL            pnl;
  EditAlignDataPtr adp;

  pnl = GetPanelFromWindow (w);
  if (pnl == NULL) return NULL;
  GetPanelExtra (pnl, &adp);
  if ( adp == NULL )
     return NULL;
  if ( adp->firstssp == NULL )
     return NULL;
  return adp;
}

extern EditAlignDataPtr GetAlignDataPanel (PaneL pnl)
{
  EditAlignDataPtr adp;

  if ( pnl == NULL )
     return NULL;
  GetPanelExtra (pnl, &adp);
  if ( adp == NULL )
     return NULL;
  if ( adp->firstssp == NULL )
     return NULL;
  return adp;
}

/**********************************************
***   AlignDataSet_Restore
***     Inval Rect
***     RestorePort
**********************************************/
static void get_client_rect (PaneL p, RectPtr prc)
{
  ObjectRect (p, prc);
  InsetRect (prc, HRZ_BORDER_WIDTH, VER_BORDER_WIDTH);
}

static void get_client_rectxy (PaneL p, RectPtr prc, Int2 x, Int2 y)
{
  ObjectRect (p, prc);
  InsetRect (prc, x, y);
}

/**************************************************************
***
***   
**************************************************************/
static SelStructPtr get_firstline (SelEdStructPtr sesp1, SelStructPtr buffer)
{
  SelEdStructPtr sesp;
  SelStructPtr   buf;

  if (buffer == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in get_firstline [25]");
         return NULL;
  }
  if (sesp1 == NULL) {
         return buffer;
  }
  for (buf=buffer; buf!=NULL; buf=buf->next) 
  {
         sesp = (SelEdStructPtr) buf->region;
         if (is_sameId (sesp1->entityID, sesp1->itemID, sesp1->itemtype, 255, sesp->entityID, sesp->itemID, sesp->itemtype, 255) ) 
            break;
  }
  if (buf == NULL) {
         return buffer;
  }
  return buf;
}

/*********************************************************************
***  FindNextSegment
*********************************************************************/
static SelEdStructPtr FindNextSegment (SelEdStructPtr current)
{
  Uint2  ei, ii;
  if ( current->next == NULL ) return NULL;
  ei = current->entityID;
  ii = current->itemID;
  current = current->next;
  while ( current != NULL )
  {
     if (current->entityID == ei &&  current->itemID == ii) break;
     current = current->next;
  }
  return current;
}

static CharPtr get_substring (CharPtr str, Int4 drw_start, Int4 drw_width)
{
  Int4            width;
  Int4            stringlens;
  CharPtr         strp;

  if (str == NULL ) 
     return NULL; 
  if (str[0]=='\0')
     return NULL;
  stringlens = StringLen (str);
  if ( drw_start >= stringlens ) { 
     return NULL; } 
  strp = str + drw_start;
  stringlens = StringLen (strp);
  if (stringlens == 0) 
     return NULL; 
  width = MIN ((Int2) drw_width, (Int2) stringlens);
  if ( !not_empty_string (strp, width) ) 
     return NULL;
  return strp;
}

static SelStructPtr go_to_next_to_draw (EditAlignDataPtr adp, Boolean next, Int2 offset)
{
  TextAlignBufPtr  curtdp;
  SelStructPtr     curvnp;
  SelEdStructPtr   curssp = NULL;
  ValNodePtr       vnp = NULL;
  Int4             from_inbuf;  /* alignment coordinates in buffer */
  Int4             from_inseq;  /* alignment coordinates in buffer */
  Int2             drw_width;   /* length of drw_str */
  Uint2            itemsubtype;
  CharPtr          curstr = NULL;
  SeqAlignPtr      salp = (SeqAlignPtr) adp->sap_align->data;
  Int2             offsettmp;

  if (adp->voffset == 0)
     return adp->buffer;
  if (offset == 0)
     return adp->firstssp;
  from_inseq = adp->hoffset;
  from_inbuf = adp->hoffset - adp->bufferstart;
  drw_width = MIN ((Int2) adp->visibleWidth, (Int2) (adp->bufferlength - from_inbuf));
  if ( drw_width <= 0 ) {
         return NULL;
  }
  curvnp = adp->buffer;
  offsettmp = 0;
  while (offsettmp < offset && from_inseq < adp->length)
  {
         curssp = (SelEdStructPtr) curvnp->region;
         itemsubtype = curvnp->itemtype;
         if (itemsubtype == EDITDEF_SCA)      
                offsettmp++;
         else if (itemsubtype == EDITDEF_SCB) 
                offsettmp;
         else if (curssp->itemtype ==OBJ_BIOSEQ && itemsubtype == FEATDEF_BAD)
         {
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )   
                   offsettmp++;
         }
         if (offsettmp < offset) {
            if (next)
               curvnp = curvnp->next;
            else
               curvnp = curvnp->prev;
         }
         if (offsettmp < offset && curvnp == NULL) 
         {
            if (next) {
               curvnp = adp->buffer;
               adp->hoffset += adp->visibleWidth;
               from_inbuf += drw_width;
               from_inseq += drw_width;
            } else {
               curvnp = adp->buffertail;
               adp->hoffset -= adp->visibleWidth;
               from_inbuf -= drw_width;
               from_inseq -= drw_width;
            } 
            if (from_inseq >= adp->length) break;
            if (from_inseq <= 0 || from_inbuf <= 0 || adp->hoffset<=0) break;
            drw_width = MIN ((Int2) adp->visibleWidth, 
                (Int2)(adp->bufferlength -(adp->hoffset -adp->bufferstart)));
         }
  }
  if (from_inseq < 0 || from_inbuf < 0 || adp->hoffset<0)  {
     adp->hoffset = adp->bufferstart;
     return adp->buffer;
  }
  if (from_inseq >= adp->length) 
     return NULL;
  return curvnp;
}


static SelStructPtr next_to_draw (EditAlignDataPtr adp, Boolean next)
{
  TextAlignBufPtr  curtdp;
  SelStructPtr     curvnp;
  SelStructPtr     tmp;
  SelEdStructPtr   curssp = NULL;
  ValNodePtr       vnp = NULL;
  SeqLocPtr        curslp;
  SeqIdPtr         sip;
  CharPtr          curstr = NULL;
  Int4             start, stop;
  Int4             start2, stop2;
  Int4             from_inbuf;  /* alignment coordinates in buffer */
  Int4             from_inseq;  /* alignment coordinates in buffer */
  Int2             drw_width;   /* length of drw_str */
  Int2             chklocp;
  Uint2            itemsubtype;
  SeqAlignPtr      salp = (SeqAlignPtr) adp->sap_align->data;
  Boolean          empty_line;

  if (next) {
     if ( adp->firstssp->next == NULL)
     {
         adp->hoffset += adp->visibleWidth;
         return adp->buffer;
     }
     else tmp = adp->firstssp->next;
  }
  if (!next) {
     if ( adp->firstssp->prev == NULL)
     {
         if (adp->hoffset == 0) {
            adp->voffset = 0;
            return adp->buffer;
         }
         adp->hoffset -= adp->visibleWidth;
         tmp = adp->buffertail;
     }
     else tmp = adp->firstssp->prev;
  }
  from_inseq = adp->hoffset;
  from_inbuf = adp->hoffset - adp->bufferstart;
  drw_width = MIN ((Int2) adp->visibleWidth, 
                   (Int2) (adp->bufferlength - from_inbuf));
  if ( drw_width <= 0 ) {
         return NULL;
  }
  empty_line = FALSE;
  curvnp = tmp;
  while ( from_inseq < adp->length )
  {
         curssp = (SelEdStructPtr) curvnp->region;
         itemsubtype = curvnp->itemtype;
         if (itemsubtype == EDITDEF_SCA)      
                return curvnp;
         else if (itemsubtype == EDITDEF_SCB) 
                return curvnp;
         else if (curssp->itemtype ==OBJ_BIOSEQ && itemsubtype == FEATDEF_BAD)
         {
            if (curssp->data !=NULL) {
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )   
                   if (adp->draw_emptyline || (!adp->draw_emptyline && !stringhasnochar(curstr, 0, drw_width) )) {
                      empty_line = FALSE;
                      return curvnp;
                   } 
                   else 
                      empty_line = TRUE;
            }
         }
         else if (curssp->itemtype==OBJ_SEQFEAT && itemsubtype==FEATDEF_PROT)
         {
            if (curssp->data !=NULL && !empty_line) {
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )                   
                   return curvnp;      
            }
         }
         else if (itemsubtype == FEATDEF_TRSL) 
         {
            if (!empty_line ) {
             while (curssp!=NULL) 
             {
               if (curssp->region !=NULL) {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2=SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord(start2, sip, salp, 0, chklocp);
                  stop2=SeqLocStop(curslp);
                  chklocp =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord(stop2, sip, salp, 0, chklocp);
                  if (start<stop && start<from_inseq+drw_width && stop>from_inseq) 
                     return curvnp;      
               }
               curssp = FindNextSegment (curssp);
             }
            }
         }
         else if (curssp->itemtype ==OBJ_SEQFEAT)
         {
            if ( !empty_line ) {
             while (curssp!=NULL) 
             {
               if (curssp->region !=NULL) { 
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2=SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord(start2, sip, salp, 0, chklocp);
                  stop2=SeqLocStop(curslp);
                  chklocp =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord(stop2, sip, salp, 0, chklocp);
                  if (start<stop && start<from_inseq+drw_width && stop>=from_inseq) 
                     return curvnp;      
               }
               curssp = FindNextSegment (curssp);
             }
            }
         }
         else if (itemsubtype == EDITDEF_CPL) 
         {
            if (curssp->data !=NULL && !empty_line) { 
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring (curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )                        
                   return curvnp;      
            }
         }

         else if( itemsubtype==SEQFEAT_GENE || itemsubtype==SEQFEAT_RNA
              ||  itemsubtype==SEQFEAT_CDREGION)
         {
            if ( !empty_line ) {
             while (curssp!=NULL) 
             {
               if (curssp->region !=NULL) {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2=SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord(start2, sip, salp, 0, chklocp);
                  stop2=SeqLocStop(curslp);
                  chklocp =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord(stop2, sip, salp, 0, chklocp);
                  if (start<stop && start<from_inseq+drw_width && stop>from_inseq) 
                     return curvnp;
               }
               curssp = FindNextSegment (curssp);
             }
            }
         }
         else if (itemsubtype>=EDITDEF_RF1 && itemsubtype<=EDITDEF_RF6) 
         {
            if (curssp->region !=NULL && !empty_line) {
                curslp = (SeqLocPtr) curssp->region;
                if ( SeqLocStop (curslp) > from_inseq && curssp->data != NULL) 
                   return curvnp; 
            }
         }
         else {
            ErrPostEx (SEV_ERROR, 0, 0, "fail in next_to_draw [46] subtype %ld type %ld", itemsubtype, curssp->itemtype); 
            return NULL;
         }
         if (next)
            curvnp = curvnp->next;
         else
            curvnp = curvnp->prev;
         if (curvnp == NULL) 
         {
            if (next) {
               curvnp = adp->buffer;
               adp->hoffset += adp->visibleWidth;
               from_inbuf += drw_width;
               from_inseq += drw_width;
            } else {
               curvnp = adp->buffertail;
               adp->hoffset -= adp->visibleWidth;
               from_inbuf -= drw_width;
               from_inseq -= drw_width;
            } 
            if (from_inseq >= adp->length) break;
            drw_width = MIN ((Int2) adp->visibleWidth, 
                (Int2)(adp->bufferlength -(adp->hoffset -adp->bufferstart)));
         }
  }
  return NULL;
}

/********************************************************
***  SetupScrollBar   
***
***  
*********************************************************/
extern BaR SeqEdGetSlateScrollBar (PaneL pnl)
{
  EditAlignDataPtr   adp;

  if (pnl)
  {
     adp=GetAlignDataPanel (pnl);
     if(adp!=NULL)
     {
        if (adp->vscrollbar_mode) {
           return GetSlateVScrollBar ((SlatE)pnl);      
        }
        return GetSlateHScrollBar ((SlatE)pnl); 
     }
  }
  return NULL;
}

extern Int2 SeqEdGetValueScrollBar (PaneL pnl)
{
  return GetValue(SeqEdGetSlateScrollBar(pnl));
}

extern void SeqEdSetValueScrollBar (PaneL pnl, Int2 value)
{
  SetValue(SeqEdGetSlateScrollBar(pnl), value);
}

extern void SeqEdCorrectBarPage (PaneL pnl, Int2 page1, Int2 page2)
{
  BaR sb;

  sb = SeqEdGetSlateScrollBar (pnl);
  CorrectBarPage (sb, page1, page2);
}

extern void SeqEdCorrectBarValue (PaneL pnl, Int2 value)
{
  BaR sb;

  sb = SeqEdGetSlateScrollBar (pnl);
  CorrectBarValue (sb, value);
}

extern void SeqEdCorrectBarMax (PaneL pnl, Int2 value)
{
  BaR sb;

  sb = SeqEdGetSlateScrollBar (pnl);
  CorrectBarMax (sb, value);
}

extern void SeqEdSetCorrectBarMax (PaneL pnl, Int4 max, Int4 offset)
{
  BaR          sb;
  Int2         cbm = 0;

  if (max < 0) 
     cbm = 0;
  else if (max > 32767) 
     cbm = 32767;
  else 
     cbm = MAX ((Int2) 0, (Int2) (max -1));
  sb = SeqEdGetSlateScrollBar (pnl);
  CorrectBarMax (sb, cbm);
  SetValue (sb, (Int2)(offset));
}

static void count_feature_buf_line (ValNodePtr fnp_list, Int4 g_left, Int4 g_right, ValNodePtr PNTR feature_line)
{
        FeatNodePtr fnp;
        Int4 c_left, c_right;
        ValNodePtr vnp;
        Boolean found;

        if(fnp_list == NULL)
                return;
        
        while(fnp_list)
        {
           fnp = fnp_list->data.ptrvalue;
           c_left = fnp->extremes.left;
           c_right = fnp->extremes.right;
           if(!(c_left > g_right || c_right < g_left))
           {
                found = FALSE;
                for(vnp = *feature_line; vnp != NULL && !found; vnp = vnp->next)
                {
                        if(vnp->data.intvalue == (Int4)(fnp->itemID))
                                found = TRUE;
                }
                if(!found)
                        ValNodeAddInt(feature_line, 0, (Int4)(fnp->itemID));
           }
           fnp_list = fnp_list->next;
        }
}

static Int4 CountTextAlignNodeNum(AlignNodePtr anp, Int4 m_left, Int4 m_right)
{
        Int4 num_line = 0;
        Int4 g_left, g_right;

        AlignSegPtr asp;
        ValNodePtr feature_line, curr;  /*the number of lines for a feature*/

        g_left = anp->extremes.left;
        g_right = anp->extremes.right;
        if(m_left > g_right || m_right < g_left)
                return 0;

        num_line = 1;
        feature_line = NULL;

        /*process  the GAPs and the DIAGs segs*/
        for(asp = anp->segs; asp !=NULL; asp = asp->next)
        {
           g_left = asp->gr.left;
           g_right = asp->gr.right;
           if(!(g_left > m_right || g_right < m_left))
           {
              switch(asp->type)
              {  
                case GAP_SEG:
                   break;

                case REG_SEG:
                case DIAG_SEG:
                   g_left = MAX(m_left, g_left);
                   g_right = MIN(m_right, g_right);
                   count_feature_buf_line (asp->cnp, g_left, g_right, &feature_line);
                   break;
                default:
                   break;
              }
           }
           if(g_left > m_right)
                break;
        }
        if(feature_line != NULL)
        {
           for(curr = feature_line; curr != NULL; curr = curr->next)
                ++num_line;
           ValNodeFree(feature_line);
        }
 
        return num_line;
}
static Int4 addline_perblock (EditAlignDataPtr adp, Int4 diffs)
{
  Int4        line = 0;
  ValNodePtr  vnp;
  SeqParamPtr prm;
  Int1        j;

  if (adp->draw_scale) line += (Int4) diffs;
  if (adp->draw_bars)  line += (Int4) diffs;
  for (vnp = adp->params; vnp != NULL; vnp = vnp->next)
  {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if ( prm->complement ) line += (Int4) diffs;
     for (j=0; j<=6; j++) 
        if (prm->rf[j]) line += (Int4) diffs;
  }
  return line;
}  
 
static Int2 feat_linenum (Int4 slp_start, Int4 slp_stop, Int4 line_len, Int4 left,
Int4 right)
{
  Int4         modstart;
  Int4         modstop;
 
  slp_start = MAX (slp_start, left);
  modstart = slp_start % line_len;
  if ( modstart > 0) slp_start -= modstart;
 
  slp_stop = MIN (slp_stop, right);
  modstop = slp_stop % line_len;
  if ( modstop > 0) slp_stop += line_len;
 
  return (Int2)((slp_stop - slp_start) / line_len);
}
 
static Int4 CountFeatNum (ValNodePtr adpfeat, Int4 line_len, Int4 left, Int4 right){
  ValNodePtr   vnp;
  SelEdStructPtr sesp;
  SeqLocPtr    slp;
  Int4         line = 0;
 
  for (vnp = adpfeat; vnp != NULL; vnp = vnp->next)
  {
     sesp = (SelEdStructPtr) vnp->data.ptrvalue;
     for (; sesp != NULL; sesp = sesp->next) 
     {
        if (vnp->choice ==FEATDEF_CDS && sesp->regiontype ==OM_REGION_SEQLOC 
        && sesp->region !=NULL)
        {
           slp = (SeqLocPtr) sesp->region;;
           if (SeqLocStart(slp) > right || SeqLocStop(slp) < left) { }
           else {
              line+= feat_linenum (SeqLocStart(slp), SeqLocStop(slp), line_len, left, right);
           }
        }
     }
  }
  return line;
}
static Int4 count_nline (EditAlignDataPtr adp, ValNodePtr anp_list, Int4 line_len, Int4 left, Int4 right, Int4 voffset)
{
        AlignNodePtr anp;
        Int4 c_start, c_stop;
        Int4 line_num = 0;
        Int4 h_block = 0;
        ValNodePtr curr;

        if(anp_list == NULL)
                return h_block;
        if(voffset == 0)
                return h_block;
        anp = anp_list->data.ptrvalue;
        if(left == -1)
                left = anp->extremes.left;
        if(right == -1)
                right = anp->extremes.right;
        if(left > anp->extremes.right || right < anp->extremes.left)
                return h_block;
        left = MAX(left, anp->extremes.left);
        right = MIN(right, anp->extremes.right);
        if (left >= right)
                return h_block;
        c_start = left;
        while(line_num < voffset)
        {
                c_stop = MIN(right, (c_start+line_len-1));
                for(curr = anp_list; curr != NULL; curr = curr->next)
                {
                        anp = curr->data.ptrvalue;
                        line_num += CountTextAlignNodeNum(anp, c_start, c_stop);
                }
                line_num += (Int4) addline_perblock (adp, 1);
                line_num += (Int4) CountFeatNum (adp->feat, line_len, c_start, c_stop);
                line_num += (Int4) CountFeatNum (adp->seqfeat, line_len, c_start, c_stop);
                if (line_num > voffset) break;
                ++h_block;
                c_start = c_stop+1;
        }
  return c_start;
}

extern Int2 hoffset2voffset (EditAlignDataPtr adp, ValNodePtr anp_list, Int4 line_len, Int4 left, Int4 right, Int4 hoffset)
{
        AlignNodePtr anp;
        Int4 c_start, c_stop;
        Int4 line_num = 0;
        Int2 preline;
        Int2 h_block = 0;
        ValNodePtr curr;
 
        if(anp_list == NULL)
                return h_block;
        if(hoffset == 0)
                return h_block;
        anp = anp_list->data.ptrvalue;
        if(left == -1)
                left = anp->extremes.left;
        if(right == -1)
                right = anp->extremes.right;
        if(left > anp->extremes.right || right < anp->extremes.left)
                return h_block;
        left = MAX(left, anp->extremes.left);
        right = MIN(right, anp->extremes.right);
        if (left >= right)
                return h_block;
        c_start = left;
        while(c_start<hoffset)
        {
                preline = line_num;
                c_stop = MIN(right, (c_start+line_len-1));
                for(curr = anp_list; curr != NULL; curr = curr->next)
                {
                        anp = curr->data.ptrvalue;
                        line_num += CountTextAlignNodeNum(anp, c_start, c_stop);
                }
                line_num += (Int4) addline_perblock (adp, 1);
                line_num += (Int4) CountFeatNum (adp->feat, line_len, c_start, c_stop);
                line_num += (Int4) CountFeatNum (adp->seqfeat, line_len, c_start, c_stop);
                ++h_block;
                c_start = c_stop+1;
        }
  return preline;
}


/*******************************************
***   Scrolling functions  
***
*** function test whether in/out buffer: 
***          not used 
***
********************************************/
extern void VscrlProc (BaR sb, SlatE s, Int2 newval, Int2 oldval)
{
  EditAlignDataPtr adp;
  RecT         r;
  Int2         pixels;
  WindoW       tempPort;
  Int4         oldhoffset;
  Int2         oldvoffset;
  Int4         temp;
  Int2         x;

  if ( s == NULL ) { 
    return; 
  }
  if ( (adp = GetAlignDataPanel ((PaneL) s)) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  tempPort = SavePort ((PaneL) s);
  Select ((PaneL) s);
  ObjectRect ((PaneL) s, &r);
  if ((newval > oldval && newval - oldval <= adp->vPage ) 
   || (newval < oldval && oldval - newval <= adp->vPage )) 
  {
     InsetRect (&r, HRZ_BORDER_WIDTH, VER_BORDER_WIDTH);
     pixels = (oldval - newval) * adp->lineheight;
     r.bottom = r.top + adp->pnlLine * adp->lineheight +1;
     r.top = r.top + 1;
     ScrollRect (&r, 0, pixels);
  } 
  oldvoffset = adp->voffset;
  adp->voffset = GetValue (sb);
  if (abs(newval - oldval) == 1)
  {
     oldhoffset = adp->hoffset;
     temp = oldhoffset  + adp->visibleWidth + adp->visibleLength;
     temp = MIN (temp, adp->bufferstart + adp->bufferlength);
     if ((oldhoffset + adp->visibleWidth > adp->bufferstart 
     && temp < adp->bufferstart + adp->bufferlength) 
     || temp == adp->length)
     {
         adp->firstssp = next_to_draw (adp, (Boolean)((newval-oldval)>0));
     }
     else {
         data_collect_arrange (adp, TRUE);
         if (temp < adp->length)
            InvalRect (&r);
     }
     if (adp->hoffset == 0 && adp->firstssp->prev == NULL && adp->voffset != 0)
     {
         adp->voffset = 0;
         CorrectBarValue (sb, (Int2) 0);
     }
     if (adp->hoffset+adp->visibleWidth>adp->length && adp->voffset<adp->nlines)
     {
        adp->voffset = adp->nlines;
        CorrectBarValue (sb, (Int2) adp->nlines);
     }
  } 
  else {
     x = adp->seqnumber;
     if (adp->draw_scale) x++;
     if (adp->draw_bars) x++;
     adp->hoffset = count_nline (adp, adp->anp_list, adp->visibleWidth, 0, adp->length-1, (Int4)((FloatLo)adp->voffset * adp->vscrollunit));
     if (adp->hoffset > adp->bufferstart 
     && adp->hoffset+adp->visibleLength+ adp->visibleWidth < adp->bufferstart +adp->bufferlength) 
     {
         data_collect_arrange (adp, FALSE);
     } else {
         data_collect_arrange (adp, TRUE);
     }
     if (x == 0) {
        adp->firstssp=go_to_next_to_draw(adp, TRUE, (Int2)0);
     } else {
        adp->firstssp=go_to_next_to_draw(adp, TRUE, (Int2)(adp->voffset%x));
     }
     InvalRect (&r);
  }
  RestorePort(tempPort);
  Update ();
}

extern void HscrlProc (BaR sb, SlatE s, Int2 newval, Int2 oldval)
{
  EditAlignDataPtr adp;
  RecT         r;
  Int2         pixels;
  WindoW       tempPort;
  Int4         oldhoffset;

  if ( s == NULL ) {
    return;
  }
  if ( (adp = GetAlignDataPanel ((PaneL) s)) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  tempPort = SavePort ((PaneL) s);
  Select ((PaneL) s);
  ObjectRect ((PaneL) s, &r);
  if ((newval > oldval && newval - oldval <= adp->hPage )
   || (newval < oldval && oldval - newval <= adp->hPage ))
  {
     InsetRect (&r, HRZ_BORDER_WIDTH, VER_BORDER_WIDTH);
     pixels = (oldval - newval) * adp->charw +1;
     ScrollRect (&r, pixels, 0);
  } else {
     InvalRect (&r);
  }
  oldhoffset = adp->hoffset;
  adp->hoffset = GetValue (sb);
  RestorePort(tempPort);
  Update ();
}

/*******************************************************************
***  
***    do_resize
***  
*******************************************************************/
extern void do_resize_panel (PaneL pnl, EditAlignDataPtr adp, Int2 width, Int2 height, Boolean rearrange)
{
  Int4         old_voffset;
  Int4         new_buffer;
  Int4         j;
  Int4         lg;
  Int2         old_visibleWidth;
  Int2         x, y;

  x = (width - adp->margin.right) / adp->charw;
  if (x < 0) 
     x = 0;
  adp->pnlWidth = x;  
  if (adp->pnlWidth < adp->marginleft + 10) {
     adp->firstssp = NULL;
     return;
  }
  x = (height - adp->margin.bottom) / adp->lineheight;
  if (x < 0) 
     x = 0;
  adp->pnlLine = x; 
  if (adp->pnlLine < 3) {
     adp->firstssp = NULL;
     return;
  }
  y = 0; x = 0;
  if (adp->columnpcell > 0) {
     y = (Int2) (adp->pnlWidth -adp->marginleft) / (Int2) adp->columnpcell;
     x = (Int2) (adp->pnlWidth -adp->marginleft -y) % (Int2)(adp->columnpcell);
     if (x == 9) 
        x = -1;
  }
  old_visibleWidth = adp->visibleWidth;
  adp->visibleWidth = (Int2) (adp->pnlWidth -adp->marginleft -y -x);
  if (adp->visibleWidth < 10) {
     adp->firstssp = NULL;
     return;
  }
  if ( adp->seqnumber == 0 ) 
     return;
  if (old_visibleWidth != adp->visibleWidth) {
     old_voffset = adp->voffset;
     adp->voffset = (Int2)(((float)old_visibleWidth/(float)adp->visibleWidth) * (float)old_voffset);
  }
  new_buffer = adp->pnlLine * adp->visibleWidth;
  if (new_buffer * 3 > adp->minbufferlength)
  {
     adp->minbufferlength = new_buffer * 3;
     if ( adp->colonne != NULL ) adp->colonne = MemFree (adp->colonne);
     lg = adp->minbufferlength + adp->editbuffer + 4;
     adp->colonne = (Int4Ptr) MemNew ((size_t) (lg * sizeof(Int4)));
     for (j=0; j<adp->minbufferlength +adp->editbuffer; j++) adp->colonne[j] = -1;
     rearrange = TRUE;
  }
  adp->vPage = adp->pnlLine - 1;  
  adp->hPage = adp->visibleWidth - 1;
  data_collect_arrange (adp, rearrange);
  SeqEdSetCorrectBarMax (pnl, adp->nlines, adp->voffset);
  SeqEdCorrectBarPage (pnl, adp->vPage, adp->vPage);
  SeqEdCorrectBarValue (pnl, SeqEdGetValueScrollBar (pnl));
}

extern void do_resize_window (PaneL pnl, EditAlignDataPtr adp, Boolean rearrange)
{
  SeqEditViewFormPtr wdp;
  WindoW       w;
  RecT         rw;    /* window rect        */
  RecT         rp;    /* panel rect         */
  RecT         rb;    /* buttons rect       */
  RecT         rct;   /* new rect for panel */
  Int2         buttonwidth,
               buttonheight;
  Int2         x, y;

  w = getwindow_frompanel (pnl);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  if (wdp == NULL) 
     return;
  ObjectRect (w, &rw);
  get_client_rect (pnl, &rp);
  GetPosition (pnl, &rp);
  x = adp->xoff;
  y = rw.bottom - rw.top;
  LoadRect ( &rct, x, adp->yoff, (Int2)(rw.right - rw.left - adp->x), 
                                 (Int2)(y - adp->y));
  SetPosition (pnl, &rct );
  AdjustPrnt (pnl, &rct, FALSE);
  ObjectRect (pnl, &rp );

  GetPosition ((GrouP) wdp->btngp, &rb);
  buttonwidth  = rb.right  - rb.left;
  buttonheight = rb.bottom - rb.top;
  LoadRect (&rb, x, (Int2)(y -buttonheight - adp->ybutt - 1),
                     (Int2)(x+ buttonwidth), (Int2)(y-adp->ybutt-1));
  SetPosition(wdp->btngp, &rb);
  AdjustPrnt (wdp->btngp, &rb, FALSE);

  ResetClip ();
  Update ();
  InsetRect (&rp, 4, 4);
  do_resize_panel (pnl, adp, (Int2)(rp.right - rp.left), (Int2)(rp.bottom - rp.top), rearrange);
  return;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/***********************************************************************
***  Draw Scale 
***         draw_scale
***         draw_bars
***
***********************************************************************/
static void draw_scale (EditAlignDataPtr adp, Int4 hoffset, Int2 scalelength, PoinT *ptlh)
{
  Char   str[128];
  Int4   scal;
  Int2   ptx, pty;
  Int4   j;
  Int2   marqueediff;

  if ( !adp->draw_scale ) return;
  SetColor (adp->colorRefs[COLOR_SCALE]);
  ptx = ptlh->x + adp->margin.left + (Int2)(adp->charw * 0.5) - 2;
  pty = ptlh->y + adp->ascent;
  marqueediff = (Int2)( 2.00 / 6.00 * adp->lineheight);
  for ( j = hoffset; j < hoffset + scalelength; j++) 
  {
         if ( adp->colonne[j] > -1) 
         {
                scal = (Int4)(adp->gr.left + 1 + adp->colonne[j]);
                if (scal % 10 == 0) 
                {
                       sprintf (str, "%d", (int)scal);
                       MoveTo ((Int2)(ptx - StringWidth(str) + (adp->charw/2) +1), 
                               (Int2)(pty + marqueediff));
                       PaintString (str);
                }
         }
         if (adp->columnpcell > 0 && j > hoffset)
            if ((Int2) j % (Int2) adp->columnpcell == 0) ptx += adp->charw;
         ptx += adp->charw;
  }
  Black();
  return ;
}

static void draw_bars (EditAlignDataPtr adp, Int4 hoffset, Int2 scalelength, PoinT *ptlh)
{
  Int4   scal;
  Int2   ptx;
  Int2   y;
  Int4   j;
  Int2   marqueelong, marqueeshort, marqueediff;

  if ( !adp->draw_bars ) return;
  SetColor (adp->colorRefs[COLOR_SCALE]);
  ptx = ptlh->x + adp->margin.left + (Int2)(adp->charw * 0.5) - 1;
  y = ptlh->y + (Int2) (2.00 / 6.00 * adp->lineheight);
  marqueelong = (Int2) (4.00 / 6.00 * adp->lineheight);
  marqueeshort= (Int2) (2.00 / 6.00 * adp->lineheight);
  marqueediff = (Int2) (2.00 / 6.00 * adp->lineheight);
  for ( j = hoffset; j < hoffset + scalelength; j++) 
  {
         if ( adp->colonne[j] > -1) 
         {
                scal = (Int4)(adp->gr.left + 1 + adp->colonne[j]);
                if (scal % 10 == 0) 
                {
                       MoveTo (ptx , y );
                       LineTo (ptx , (Int2)(y + marqueelong));
                }
                else if (scal % 5 == 0) 
                {
                       MoveTo (ptx , (Int2)(y + marqueediff));
                       LineTo (ptx , (Int2)(y + marqueediff + marqueeshort));
                }
         }
         if (adp->columnpcell > 0 && j > hoffset)
            if ((Int2) j % (Int2) adp->columnpcell == 0) ptx += adp->charw;
         ptx += adp->charw;
  }
  Black();
  return ;
}
         
/***********************************************************************
***  draw_id  
************************************************************************/
static void draw_id (EditAlignDataPtr adp, PoinT *pt, Int2 index, CharPtr strid, Uint1 strand, Int4 pos, Uint2 itemtype, Boolean idselected, Boolean is_master, Int2 group)
{
  Char     str[128], str1[128];
  RecT     rct;
  CharPtr  tmp;
  Int2     stringlens;
  Int2     j;
  Int2     total = 0;
  Int2     posspace;

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!*/ idselected = FALSE;

  if (adp->length < 1000) { posspace = 4;
  } else if (adp->length < 10000) { posspace = 5;
  } else if (adp->length < 100000) { posspace = 6;
  } else if (adp->length < 1000000) { posspace = 7;
  } else { posspace = 8;
  }

  if (adp->marginwithindex) total += 3;
  if (adp->marginwithgroup) total += 4;
  if (adp->marginwithpos) total += posspace;
  total += 2;

  str[0] = '\0';
  tmp = str;
  if (adp->marginwithindex && index > 0) {
     sprintf (str1, "%3d ",(int) index);
  } else {
     sprintf (str1, " ");
  }
  tmp = StringMove (tmp, str1);
  tmp[0] = '\0';
  
  adp->marginwithgroup = FALSE;
  if (adp->marginwithgroup) {
     if (group >= 0)
        sprintf (str1, "%4d ",(int) group);
     else sprintf (str1, "    ");
     tmp = StringMove (tmp, str1);
     tmp[0] = '\0';
  }

  if (strid != NULL && strid[0]!='\0') {
     stringlens = (Int2)(StringLen (strid));
/******!!!!!!!!!!!!!!!
     if (stringlens > 4) {
         if (strid[0] == 'l' && strid[1] == 'c' && strid[2] == 'l') {
            strid += 4;
            stringlens -= 4;
         }
     }
!!!!!!!!!!!!!!!!!*********/
     for (j=0; j<stringlens; j++) 
        str1 [j] = strid [j]; 
     str1 [j] = '\0';
     if (stringlens < adp->size_labels) {
        for (j=stringlens; j<adp->size_labels; j++)
              str1 [j] = ' ';
        str1 [j] = '\0';
        stringlens = (Int2)(StringLen (str1));
     }
     stringlens = MIN (adp->marginleft - total, stringlens);
     str1 [stringlens] = '\0';
     tmp = StringMove (tmp, str1);         
  } 
  else if (adp->marginleft > 20) {
     sprintf (str1, "          ");
     tmp = StringMove (tmp, str1);
  }
/*
  symbol[0] = ' ';
  if (strand == Seq_strand_minus) 
     symbol[1] = '<';
  else 
     symbol[1] = '>';
  symbol[2] = '\0';
  tmp = StringMove (tmp, symbol);
*/
  if (adp->marginwithpos) {
     if (posspace <= 4) {
        if (pos > 0)
           sprintf (str1, "%4ld ", (long) pos);
        else sprintf (str1, "     ");
     } else if (posspace == 5) {
        if (pos > 0) 
           sprintf (str1, "%5ld ", (long) pos); 
        else sprintf (str1, "      ");
     } else if (posspace == 6) {
        if (pos > 0)  
           sprintf (str1, "%6ld ", (long) pos);  
        else sprintf (str1, "       "); 
     } else {
        if (pos > 0)  
           sprintf (str1, "%7ld ", (long) pos);
        else sprintf (str1, "        ");
     } 
     tmp = StringMove (tmp, str1);
  }
  *tmp = '\0';

  stringlens = (Int2)(StringLen(str));
  if (stringlens < adp->marginleft -1)
         for (j = stringlens; j < adp->marginleft; j++, tmp++) 
            *tmp = ' '; 
  str [adp->marginleft -1] = ' ';
  str [adp->marginleft] = '\0';
  if ( itemtype != OBJ_BIOSEQ && adp->displaytype )
         SelectColor (162, 163, 82);
  else if (is_master && !adp->all_sequences)
         SetColor (adp->colorRefs[COLOR_ID_MASTER]);
  else 
         SetColor (adp->colorRefs[COLOR_ID]);
  if ( !idselected ) 
  {
Black();
         MoveTo (pt->x, (Int2)(pt->y + adp->ascent));
         PaintString (str);
  }
  else  {
         InvertColors ();
         LoadRect(&rct, pt->x, pt->y, (Int2)(pt->x +stringlens *adp->charw -2), 
                         (Int2)(pt->y + adp->lineheight));
         EraseRect (&rct);
         MoveTo (pt->x, (Int2)(pt->y +adp->ascent));
         PaintString (str);
         InvertColors();
  } 
  Black();
  return;
}

/******************************************************************
***      paint_caret,  draw caret : draws caret
***
*** Now caret is a T between basis of letters
*** Before it was a bar between letters:

  pt1.x = pt.x + ( column + marginleft ) * charw -1;
  pt1.y = pt.y;
  pt2.x = pt1.x;
  pt2.y = pt1.y + lineheight -2;
***
*******************************************************************/
static void paint_caret (PoinT pt, Int2 column, Int2 charw, Int2 lineheight, Int2 marginleft)
{
  PoinT pt1, pt2;
  Int2  lenv = 3;
  Int2  lenh = 2;
  Int2  row;

  Red ();
  WidePen (2);
  row = pt.x + ( column + marginleft ) * charw -1;
  pt1.x = row; 
  pt1.y = pt.y  + lineheight -1;
  pt2.x = row;
  pt2.y = pt1.y  -lenv;
  DrawLine (pt1, pt2);
  pt2.x = row +lenh;
  pt2.y = pt1.y;
  DrawLine (pt1, pt2);
  pt2.x = row -lenh;
  pt2.y = pt1.y;
  DrawLine (pt1, pt2);
  WidePen (1);
  Black ();
} 

static void draw_caret (EditAlignDataPtr adp, SelEdStructPtr sesp, PoinT pt, Int2 line)
{
  SeqLocPtr    slp;
  Int2         column;  
  
  if (!adp->display_panel) 
  {
     if ( is_samess_ses (&(adp->caret), sesp))
     {  
        slp = (SeqLocPtr) adp->caret.region;
        if(SeqPosInLineColumn (SeqLocStart(slp), adp->alignline[line], &column, adp->hoffset, adp) )
        {
           paint_caret (pt, column, adp->charw, adp->lineheight, adp->marginleft);
        }
     }
  }
}

/*********************************************************************
***      PaintSubseq : paint visible sequence  
**********************************************************************/
static void PaintSubseq (Int4 k, Int4 lg, Int4 from, PoinT *pt, CharPtr vstr, Boolean invert, EditAlignDataPtr adp, Boolean draw_diff, CharPtr *masterstr, Boolean cmplt, EditCellPtr ecp) 
{
  CharPtr      strPtr;
  Char         str[512];
  Uint4        curColor;
  Int4         to = k + lg;
  Int2         caret_color = 0;
  CharPtr      masterstrptr = NULL;
  RecT         rct;
  Int4         lg2 = MIN (lg, adp->visibleWidth); 
  Uint4        blackColor = GetColorRGB(0,0,0),
               newColor;
  Boolean      pretty = (Boolean) (ecp!=NULL);

  dashedstring ((CharPtr) str, 512);
  strPtr = str;
  if (pretty) {
     curColor = GetColorRGB (ecp[from].r, ecp[from].g, ecp[from].b);
  }
  else {
     curColor = adp->colorRefs[(Uint1)(vstr[k] - '*')];
  } 
  if (masterstr != NULL) {
     if (*masterstr != NULL) 
        masterstrptr = *masterstr;
  }
  while ( k < to ) 
  {
     if (pretty) {
        newColor = GetColorRGB (ecp[k+from].r, ecp[k+from].g, ecp[k+from].b);
     }
     else {
        newColor = adp->colorRefs[(Uint1)(vstr[k] - '*')];
     }
     if ( adp->colonne [from +k -adp->bufferstart] < 0 ) 
     {
        *strPtr = 0;
        strPtr = str;
        if (cmplt) 
           strPtr = complement_string (strPtr);
        MoveTo ((Int2)(pt->x +adp->margin.left), (Int2)(pt->y + adp->ascent));
        PaintString (strPtr);
        pt->x += caret_color * adp->charw;
        caret_color = 1;
        strPtr = str;
        Black(); 
        if ((from +k -adp->bufferstart) > adp->length) 
           return;
        pt->x += adp->intersalpwidth * adp->charw;
        k += adp->intersalpwidth;
        SetColor (curColor);
     }
     else if (curColor != newColor)
     {
        *strPtr = 0;
        strPtr = str;
        if (cmplt) 
           strPtr = complement_string (strPtr);
        MoveTo ((Int2)(pt->x +adp->margin.left), (Int2)(pt->y +adp->ascent));
        if (invert)
        {
           SetColor (adp->colorRefs[COLOR_SELECT]);
           InvertColors ();
           LoadRect (&rct, (Int2)(pt->x+adp->margin.left), pt->y, (Int2)(pt->x +adp->margin.left+ caret_color*adp->charw), (Int2)(pt->y +adp->lineheight -1));
           EraseRect (&rct);
           InvertColors ();
           if (adp->colorRefs[COLOR_SELECT]==blackColor)
              White();
           else  {
              Black ();
           }
        } 
        else {
           SetColor (curColor);
        }
        PaintString (strPtr);
        pt->x += caret_color *adp->charw;
        caret_color = 1;
        strPtr = str;
        if ( draw_diff && masterstrptr != NULL ) 
        {
           if (vstr[k] !='-') {
              if ( *masterstrptr == vstr[k] || TO_LOWER(*masterstrptr) == TO_LOWER(vstr[k])) 
                 *strPtr = '.';
              else 
                 *strPtr = vstr[k];
           }
           else 
              *strPtr = vstr[k];
           masterstrptr++;
        } 
        else 
           *strPtr = vstr[k];
        strPtr++;
        if (adp->columnpcell > 0) 
        {
           if ((Int2) (k+1) % (Int2) adp->columnpcell == 0) {
              *strPtr = ' ';
              strPtr++;
              caret_color++;
           }
        }
        curColor = newColor;
        k++;
     } 
     else {
        if ( draw_diff && masterstrptr != NULL ) 
        {
           if (vstr[k] !='-') {
              if ( *masterstrptr == vstr[k] || TO_LOWER(*masterstrptr) == TO_LOWER(vstr[k]))
                 *strPtr = '.';
              else 
                 *strPtr = vstr[k];
           }
           else
              *strPtr = vstr[k];
           masterstrptr++;
        } 
        else 
           *strPtr = vstr[k];
        strPtr++;
        caret_color++;
        if (adp->columnpcell > 0) 
        {
           if ((Int2) (k+1) % (Int2) adp->columnpcell == 0) 
           {
              *strPtr = ' ';
              strPtr++;
              caret_color++;
           }
        }
        k++;
     }
  }
  *strPtr = 0;
  strPtr = str;
  if (cmplt) 
     strPtr = complement_string (strPtr);
  MoveTo ((Int2)(pt->x +adp->margin.left), (Int2)(pt->y +adp->ascent));
  if (invert)
  {
     SetColor (adp->colorRefs[COLOR_SELECT]); 
     InvertColors ();
     LoadRect (&rct, (Int2)(pt->x+adp->margin.left), pt->y, (Int2)(pt->x +adp->margin.left+ caret_color*adp->charw), (Int2)(pt->y +adp->lineheight -1));
     EraseRect (&rct);
     InvertColors ();
     if(adp->colorRefs[COLOR_SELECT]==blackColor)
        White();
     else {
        Black ();   
     }
  }
  else {
     SetColor (curColor);
  }
  PaintString (strPtr);
  Black ();
  pt->x += caret_color *adp->charw;
  *masterstr = masterstrptr;
  return;
} 

/*****************************************************************
***  draw_seq
******************************************************************/
static ValNodePtr deleteseqloc (ValNodePtr head,  ValNodePtr delp)
{
  ValNodePtr next = NULL,
             pre = NULL,
             tmp = NULL;
  
  tmp = head; 
  while (tmp!=NULL) 
  {
     next = tmp->next;
     if (tmp == delp) 
     {
        if (pre==NULL) {
           head = next;
        }
        else {
           pre->next = next;
        }
        tmp->next = NULL;
        tmp->data.ptrvalue = SeqLocFree (tmp->data.ptrvalue);
        ValNodeFree (tmp);
     }
     else 
        pre = tmp;
     tmp = next;
  }
  return head;
}


static ValNodePtr simplifySeqLocList (ValNodePtr valnode)
{
  ValNodePtr   vnpa = NULL,
               vnpb = NULL;
  SeqLocPtr    slpa,
               slpb;
  SeqIntPtr    sint;
  Boolean      check = TRUE,
               loopin = TRUE;

  while (check)
  {
     check = FALSE;
     loopin = TRUE;
     vnpa = valnode;
     while (vnpa != NULL && loopin)
     {
           for (vnpb = vnpa->next; vnpb != NULL; vnpb = vnpb->next)
           {
                 slpa = (SeqLocPtr)vnpa->data.ptrvalue;
                 slpb = (SeqLocPtr)vnpb->data.ptrvalue;
                 if (SeqLocCompare (slpa, slpb) == SLC_A_IN_B) {
                    valnode = deleteseqloc (valnode, vnpa);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_B_IN_A) {
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_A_EQ_B) {
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocCompare (slpa, slpb) == SLC_A_OVERLAP_B) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    if (SeqLocStart(slpa) < SeqLocStart(slpb)) {
                       sint->to = SeqLocStop(slpb);
                    }
                    else {
                       sint = (SeqIntPtr) slpa->data.ptrvalue;
                       sint->from = SeqLocStart(slpb);
                    }
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocStop(slpa) == SeqLocStart(slpb)-1) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    sint->to = SeqLocStop(slpb);
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
                 else if (SeqLocStart(slpa) == SeqLocStop(slpb)+1) {
                    sint = (SeqIntPtr) slpa->data.ptrvalue;
                    sint->from = SeqLocStart(slpb);
                    valnode = deleteseqloc (valnode, vnpb);
                    check = TRUE;
                    loopin = FALSE;
                    break;
                 }
           }
           if (loopin && vnpa != NULL)
              vnpa = vnpa->next;
     }
  }
  return valnode;
}


static Int4 getnextpos (ValNodePtr vnp, Int4 *stop)
{
  ValNodePtr tmp;
  SeqLocPtr  slp;
  SeqIntPtr  sint;
  Int4       valmin = (Int4)(INT4_MAX-2), 
             val = (Int4)-1,
             val2 = (Int4)-1;

  for (tmp=vnp; tmp!=NULL; tmp=tmp->next) 
  {
     slp = (SeqLocPtr) tmp->data.ptrvalue;
     if (SeqLocStart(slp) < valmin) {
        valmin = SeqLocStart(slp);
     }
  }
  if (valmin < (Int4)(INT4_MAX-2)) 
  {
     for (tmp=vnp; tmp!=NULL; tmp=tmp->next) 
     {
        slp = (SeqLocPtr) tmp->data.ptrvalue;
        if (SeqLocStart(slp) == valmin) {
           val = SeqLocStart(slp);
           val2 = SeqLocStop (slp);
           sint = (SeqIntPtr) slp->data.ptrvalue;
           sint->from = INT4_MAX;
           sint->to = INT4_MAX;
           break;
        }
     }
  }
  *stop = val2;
  return val;
}

static void draw_seq (EditAlignDataPtr adp, SelEdStructPtr sesp, Int4 from, Int4 drw_width, PoinT ptlh, CharPtr drw_str, Boolean draw_diff, CharPtr masterstr, Boolean cplmt, Uint2 itemtype, Boolean check_selection)
{
  ValNodePtr    vnp;
  SelEdStructPtr sesp2;
  EditCellPtr   ecp;
  SelStructPtr  ssptmp;
  SeqIdPtr      sip;
  Int4          start, stop;
  Int4          startb, stopb;
  Int4          left, right;
  Int4          seg_lg;
  Int4          k,
                k1;
  Int2          chklocp;
  Int2          ptx, pty;
  Int2          iCount = 0;
  SeqLocPtr     slp;
  ValNodePtr    startp = NULL;
  
  ptx = ptlh.x +adp->margin.left;
  pty = ptlh.y;
  adp->start_select = -1;
  start = stop = startb = stopb = -1;
  ssptmp = NULL;
  if (check_selection) 
  {
    ssptmp = ObjMgrGetSelected (); 
    for (; ssptmp != NULL; ssptmp = ssptmp->next) 
    {
     if ( checkssp_for_editor (ssptmp) && is_samess_ses (ssptmp, sesp) ) 
     {
         sip = SeqLocId ((SeqLocPtr) ssptmp->region);
         start = SeqLocStart ((SeqLocPtr) ssptmp->region);
         adp->start_select = start;
         chklocp =chkloc(sip, start, adp->sqloc_list, &start);
         start = SeqCoordToAlignCoord (start, sip, (SeqAlignPtr) adp->sap_align->data, 0, chklocp);
         stop  = SeqLocStop  ((SeqLocPtr) ssptmp->region);
         chklocp =chkloc(sip, stop, adp->sqloc_list, &stop);
         stop = SeqCoordToAlignCoord (stop, sip, (SeqAlignPtr) adp->sap_align->data, 0, chklocp);

         if (start<=stop && (start < from + drw_width && (stop + 1) > from)) 
         {
               slp = SeqLocIntNew (start, stop, Seq_strand_plus, sip); 
               ValNodeAddPointer (&startp, 0, (Pointer)slp);
               iCount ++;
         }
     }
    }
  }

  ecp = NULL;
  if (itemtype == OBJ_BIOSEQ && adp->copyalign!=NULL) {
     for (vnp=adp->copyalign; vnp!=NULL; vnp=vnp->next) {
        sesp2 = (SelEdStructPtr) vnp->data.ptrvalue;
        if (sesp2->entityID==sesp->entityID && sesp2->itemID==sesp->itemID) 
        {
           ecp=(EditCellPtr)(sesp2->data->data.ptrvalue);
           if (ecp!=NULL) {
              break;
           }
        }   
     }
  }
 
  if (iCount > 0) {
     k = 0;
     k1 = from;
     startp = simplifySeqLocList (startp);
     start = getnextpos (startp, &stop);
     while ( k < drw_width && start >= 0) {
        if ( k + from < start ) 
        {
           left = (Int4) k1;
           right= start;
           seg_lg = MIN (right - left, drw_width);
           PaintSubseq (k, seg_lg, from, &ptlh, drw_str, FALSE, adp, draw_diff, &masterstr, cplmt, ecp);
           k += seg_lg;
        } 
        else if ( k +from >= start && k +from < stop + 1) 
        {
           left = (Int4) MAX (k1, start);
           right= (Int4) MIN (from + drw_width, stop + 1);
           seg_lg = MIN (right - left, drw_width);
           PaintSubseq (k, seg_lg, from,  &ptlh, drw_str, TRUE, adp, draw_diff, &masterstr, cplmt, ecp);
           k += seg_lg;
           k1 = stop + 1;
           start = getnextpos (startp, &stop);
        }         
     }
     if ( k < drw_width) 
     {
        left = (Int4) k1;
        right= (Int4) from + drw_width;
        seg_lg = MIN (right - left, drw_width);
        PaintSubseq (k, seg_lg, from,  &ptlh, drw_str, FALSE, adp, draw_diff, &masterstr, cplmt, ecp);
     }
  }
  else {
     PaintSubseq (0, drw_width, from, &ptlh, drw_str, FALSE, adp, draw_diff, &masterstr, cplmt, ecp);
  }
  if (startp != NULL)
     ValNodeFree (startp);
  return;
}

/*****************************************************************
***  draw_line
******************************************************************/
static void draw_line (EditAlignDataPtr adp, Int4 start, Int4 stop, Uint1 strand, PoinT *pt, Int4 from, Int4 drw_width, Int2 line, Int2 alignline, Uint4 color, Uint2 wideline, Boolean selected, Boolean partialstart, Boolean partialstop, BoolPtr gapline)
{
  PoinT         pt1, pt2,
                pttmp;
  PoinT         oldpt1, oldpt2;
  FloatLo       hgt = (FloatLo)(2.5 / 4.0);
  RecT          leftrct, rightrct;
  Int4          left, right;
  Int2          yline;  
  Int2          column;  
  Int2          above, below;
  Int4          j, k;
  Int2          col1;
  Uint1         startin, stopin;
  BoolPtr       gap;
  
  left = (Int4) MAX ((Int4) from, (Int4) start);
  SeqPosInLineColumn (left, alignline, &column, adp->hoffset, adp);
  yline = pt->y + (Int2) (hgt * adp->ascent);
  pt1.x = pt->x + (column + adp->marginleft) * adp->charw;
  pt1.y = yline;

  if (adp->columnpcell > 0) {
     col1 = column - (column-1)/adp->columnpcell;
  }
  right = (Int4) MIN ((Int4)from+drw_width-1, (Int4)stop);
  SeqPosInLineColumn (right, alignline, &column, adp->hoffset, adp);
  pt2.x = pt->x + (column + adp->marginleft + 1) * adp->charw;
  pt2.y = yline;

  startin = stopin = 0;
  if ( start >= from ) {
     if (strand == Seq_strand_minus)  
        startin = 1;
     else
        startin = 2;
  } 
  if (stop < from+drw_width && stop != start) {
     if (strand == Seq_strand_minus) 
        stopin = 1;
     else 
        stopin = 2;
  }
  oldpt1.x = pt1.x;
  oldpt1.y = pt1.y;
  oldpt2.x = pt2.x;
  oldpt2.y = pt2.y;
  if (start<stop) {
   if (startin > 0)
     pt1.x += 5;
  
   SetColor (color);
   if (gapline == NULL) {
     if (stopin > 0)
        pt2.x -= 5;
     WidePen (wideline);
     DrawLine (pt1, pt2);
   } 
   else {
     pttmp.y = yline;
     j=k=left;
     for (j=0, gap=gapline; j<col1; j++) 
        gap++;
     while (j<right) {
        if (*gap) {
           while (*gap && k<right) {
              gap++;
              k++;
           }
           WidePen (1);
           Dashed ();
        }
        else {
           while (!(*gap) && k<right) {
              gap++;
              k++;
           }
           WidePen (wideline);
        }
        SeqPosInLineColumn (k, alignline, &column, adp->hoffset, adp);
        if (k==right && stopin == 0)  {
           column++;
        }
        pttmp.x = pt->x + (column + adp->marginleft) * adp->charw;
        DrawLine (pt1, pttmp);
        Solid ();
        pt1.x = pttmp.x;
        j = k;
     }
   }
  }
  if (selected) {
#ifdef WIN_MAC
         above = 2;
         below = 4;
#endif
#ifdef WIN_MSWIN
         above = 3;
         below = 3;
#endif
#ifdef WIN_MOTIF
         above = 3;
         below = 3;
#endif
         WidePen (1);
         Black ();
         pt1.y = pt2.y = yline + below;
         DrawLine (pt1, pt2);
         pt1.y = pt2.y = yline - above;
         DrawLine (pt1, pt2);
         White ();
         pt1.y = pt2.y = yline + below - 1;
         DrawLine (pt1, pt2);
         pt1.y = pt2.y = yline - above + 1;
         DrawLine (pt1, pt2);
  }
  WidePen (1);
  SetColor (color);
  pt1.x = oldpt1.x;
  pt1.y = oldpt1.y;
  pt2.x = oldpt2.x;
  pt2.y = oldpt2.y;
  if (startin == 1) {
#ifdef WIN_MAC
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4) );
     OffsetRect (&leftrct, 0, 1);
#endif
#ifdef WIN_MSWIN
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4) );
#endif
#ifdef WIN_MOTIF
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2) (pt1.y +4) );
#endif
     if (partialstart)
        CopyBits (&leftrct, leftTriSym);
     else 
        CopyBits (&leftrct, leftTriFillSym);
  } 
  else if (startin == 2) {
#ifdef WIN_MAC
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4));
     OffsetRect (&leftrct, 0, 1);
#endif
#ifdef WIN_MSWIN
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4));
#endif
#ifdef WIN_MOTIF
     LoadRect(&leftrct, pt1.x , (Int2)(pt1.y -3), (Int2)(pt1.x +7), (Int2)(pt1.y +4));
#endif
     if (partialstart)
        CopyBits (&leftrct, rectSym);
     else
        CopyBits (&leftrct, rectFillSym);
  }
  if (stopin == 1) {
#ifdef WIN_MAC
     LoadRect(&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
     OffsetRect (&rightrct, 0, 1);
#endif
#ifdef WIN_MSWIN
     LoadRect(&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
#endif
#ifdef WIN_MOTIF
     LoadRect(&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
#endif
     if (partialstop)
        CopyBits (&leftrct, rectSym);
     else 
        CopyBits (&rightrct, rectFillSym);
  }
  else if (stopin == 2) {
#ifdef WIN_MAC
     LoadRect (&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
     OffsetRect (&rightrct, 0, 1);
#endif
#ifdef WIN_MSWIN
     LoadRect (&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
#endif
#ifdef WIN_MOTIF
     LoadRect (&rightrct, (Int2)(pt2.x -7), (Int2)(pt2.y -3), pt2.x, (Int2)(pt2.y +4));
#endif
     if (partialstop)
        CopyBits (&rightrct,rightTriSym);
     else 
        CopyBits (&rightrct,rightTriFillSym);
  }
  Black ();
  WidePen (1);
}

/*********************************************************************
***  draw_cds
*********************************************************************/
static void draw_trans (EditAlignDataPtr adp, SelEdStructPtr cds, CharPtr trans, PoinT *pt, Int4 from, Int4 drw_width, Uint2 offset)
{
  SeqLocPtr     slp;
  CharPtr       transtr;
  CharPtr       strPtr, transPtr;
  Int4          left, right;
  Int2          k;  

  if ( trans == NULL || cds == NULL ) return;
  slp = (SeqLocPtr) cds->region;
  if (slp == NULL) return;
  left  = (Int4) MAX ((Int4) from, (Int4) SeqLocStart(slp));
  right = (Int4) MIN ((Int4) from+drw_width, (Int4) SeqLocStop(slp)+1);
  if (right < left) {
     ErrPostEx (SEV_ERROR, 0, 0, "fail in draw_trans: left %ld right %ld",left,right);
  }
  transtr = MemNew ((size_t)(512 * (sizeof (Char))));
  emptystring (transtr, 512);
  transtr[0] = ' ';
  transtr[511] = '\0';
  strPtr = transtr;
  strPtr += (Int2) (left - from);
  transPtr = trans;
  transPtr += cds->offset;   /*!!!!!!!!!!!!! += cds->codonstart ?? */
  if (SeqLocStart(slp) < from) 
     transPtr += (Int2)(from - SeqLocStart(slp));  
  if (adp->prot_mode == MPROTAA) {
     for (k =0; k < (right - left); k++, strPtr++, transPtr++) {
        if ( *transPtr=='M' || *transPtr=='*') *strPtr = *transPtr;
     }
  } 
  else {
     for (k =0; k < (right - left); k++, strPtr++, transPtr++) {
        *strPtr = *transPtr;
     }
  }
  strPtr++;
  *strPtr= '\0';
  draw_seq (adp, cds, from, drw_width, *pt, transtr, FALSE, NULL, TRUE, OBJ_SEQFEAT, FALSE);
  transtr = MemFree (transtr);
  return;
}

/*********************************************************************
***  draw_pept
*********************************************************************/
static void draw_pept (EditAlignDataPtr adp, SelEdStructPtr cds, PoinT *pt, Int4 from, Int4 drw_width, Int2 line, Uint2 offset)
{
  ValNodePtr  pept;
  
  if (cds!=NULL) {
     pept = cds->data;
     if (pept != NULL)  { 
        if (pept->data.ptrvalue != NULL)  { 
           draw_trans (adp, cds, (CharPtr) pept->data.ptrvalue, pt, from, drw_width, offset);
        }
     }
  }
}

/*********************************************************************
***  draw_feat
*********************************************************************/
static void draw_feat (EditAlignDataPtr adp, Uint2 eid, Uint2 iid, Uint2 it, Uint2 choice, Int4 start, Int4 stop, Uint1 strand, CharPtr label, PoinT *pt, Int4 from, Int4 drw_width, Int2 line, Int2 alignline, Uint4 color, Boolean partialstart, Boolean partialstop, BoolPtr gapline)
{
  Char      str[128];
  Uint2     wideline;
  Boolean   featselect;

  if ( start > stop || start >= from + drw_width || stop < from ) 
     return;
  if (it == OBJ_VIRT) 
  {
     wideline = 5;
     if (choice == SEQFEAT_CDREGION ) 
         sprintf (str, "CDS %d", (int) iid);
     else if (choice == SEQFEAT_GENE) 
         sprintf (str, "gene %d", (int) iid);
     else if (choice == SEQFEAT_RNA)
         sprintf (str, "mRNA %d", (int) iid);
     str[15] = '\0';
     draw_id (adp, pt, -1, str, strand, 0, 3, FALSE, FALSE, -1);
  }
  else {
     wideline = 3;
     draw_id (adp, pt, -1, label, strand, 0, 3, FALSE, FALSE, -1);
  }
  featselect = (Boolean) (is_selectedbyID (eid, iid, it) != NULL);
  draw_line (adp, start, stop, strand, pt, from, drw_width, line, alignline, color, wideline, featselect, partialstart, partialstop, gapline);
}

static void what_inline (EditAlignDataPtr adp, Int2 line, Uint2 ei, Uint2 ii, Uint2 it, Uint2 ist, Uint2 al)
{ 
  adp->seqEntity_id [line]  = ei;
  adp->item_id [line]  = ii;
  adp->itemtype[line] = it;
  adp->itemsubtype[line] = ist;
  adp->alignline [line]= al;
}

static Boolean id_is_invalid (EditAlignDataPtr adp, PoinT ptlh)
{
  RecT  rect; 
  LoadRect(&rect, ptlh.x, ptlh.y, 
                  (Int2)(ptlh.x + adp->marginleft * adp->charw), 
                  (Int2)(ptlh.y + adp->lineheight));
  return RectInRgn (&rect, updateRgn);
}

static Boolean row_is_invalid (EditAlignDataPtr adp, PoinT ptlh)
{
  RecT  rect; 
  LoadRect (&rect, (Int2)(ptlh.x + adp->margin.left - adp->charw), ptlh.y,  
                   (Int2)(ptlh.x + adp->pnlWidth * adp->charw + 1), 
                   (Int2)(ptlh.y + adp->lineheight));
  return RectInRgn (&rect, updateRgn);
}

static Int2 rang (Uint2 format, SelEdStructPtr ssp)
{
  if (format == OBJ_BIOSEQ)
     return ssp->entityID;
  if (format == OBJ_SEQALIGN)
     return ssp->itemID;
  return 0;
}

static Int2 getparam (ValNodePtr vnprm, Uint2 eID, Uint2 iID, Uint1 choice)
{
  ValNodePtr   vnp;
  SeqParamPtr  prm;

  for (vnp = vnprm; vnp != NULL; vnp = vnp->next) {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if (prm->entityID == eID && prm->itemID == iID)
     {
        if (choice == 1)
           return prm->group;
     }
  }
  return -1;
}

static CharPtr seqid_tolabel (SeqIdPtr sip, Uint2 choice)
{
  BioseqPtr bsp;
  SeqIdPtr  tmp = NULL;
  Char      str[120];
  CharPtr   strp;

  str[0] = '\0';
  if (choice > 0 && choice <= PRINTID_GIcc)
  {
     bsp = BioseqLockById(sip);
     if (bsp!=NULL) 
     {
        if (choice==PRINTID_GIcc)
           tmp = SeqIdFindBest (bsp->id, 0);
        else 
           tmp = bsp->id;
        if (tmp)
        {
           SeqIdWrite (tmp, str, choice, 120);
           BioseqUnlock(bsp);
           strp = StringSave (str);
           return strp;
        }
        BioseqUnlock(bsp);
     }
  }
  if (choice!=0 && sip) {
     SeqIdWrite (sip, str, choice, 120);
     strp = StringSave (str);
     return strp;
  } 
  strp = StringSave (str);
  return strp;
}

/*****************************************************************
***
***     Draw Alignment in Panel p 
***
******************************************************************/
static Char strcmpl[] = "complement";

extern void on_draw (PaneL p)
{
  PoinT            ptlh;        /* PoinT at the top of lineheight */
  RecT             rp;          /* Panel Rect */
  EditAlignDataPtr adp;
  TextAlignBufPtr  curtdp;
  SelStructPtr     curvnp;
  SelEdStructPtr   curssp = NULL;
  ValNodePtr       vnp = NULL;
  SeqLocPtr        curslp;
  SeqLocPtr        sip;
  SeqAlignPtr      salp;

  Int4          from_inbuf;  /* alignment coordinates in buffer */
  Int4          from_inseq;  /* alignment coordinates in buffer */
  Int4          stringlens;
  Int4          start,  stop;
  Int4          start2, stop2;

  Int2          line;        /* current line */
  Int2          curalgline;
  Int2          drw_width;   /* length of drw_str */
  Int2          group;
  Int2          chklocp,
                chklocp2;
  Uint2         itemsubtype;
  Uint2         offset;
 
  CharPtr       stridp;
  CharPtr       curstr = NULL;
  CharPtr       masterbuf =NULL, 
                masterstr =NULL;
  Int2          oldstyle;

  Boolean       invalid_id;
  Boolean       invalid_row;
  Boolean       is_master, draw_diff;
  Boolean       first = TRUE;
  Boolean       line_drawn;
  Boolean       empty_line;

  Int1          gapinline;
  BoolPtr       gapline,
                gaplinep;
  
  if ( (adp = GetAlignDataPanel (p)) == NULL ) 
         return;
  if ( adp->seqnumber == 0 ) return;
  salp = (SeqAlignPtr) adp->sap_align->data;

  if ( adp->firstssp == NULL ) {
     adp->firstssp = get_firstline (NULL, adp->buffer);
  }
  if ( adp->firstssp == NULL ) return;
  if ( adp->firstssp->region == NULL ) return;
  from_inseq = adp->hoffset;
  from_inbuf = adp->hoffset - adp->bufferstart;
  drw_width = MIN ((Int2) adp->visibleWidth, 
                   (Int2) (adp->bufferlength - from_inbuf));
  if ( drw_width <= 0 || adp->bufferstart < 0 ) {
         return;
  }

  oldstyle = GetMuskCurrentSt ();
  SetMuskCurrentSt (GetMuskStyleName (adp->styleNum));

  line = (drw_width +2) * sizeof(Boolean);
  gapline = (BoolPtr) MemNew ((size_t)(line));
  MemSet ((Pointer)gapline, (int)(FALSE), (size_t)(line)); 
  gaplinep = NULL;
  gapinline = LINE_NOGAP;
  empty_line = FALSE;
 
  curssp = (SelEdStructPtr) adp->firstssp->region;
  if (curssp->itemtype ==OBJ_SEQFEAT) {
     for (curvnp=adp->firstssp; curvnp!=NULL; curvnp=curvnp->prev) {
        curssp = (SelEdStructPtr) curvnp->region;
        if (curssp->itemtype==OBJ_BIOSEQ && curvnp->itemtype== FEATDEF_BAD)
        {
           if (curssp->data !=NULL) {
              vnp = (ValNodePtr) curssp->data;
              curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
              curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
              if ( curstr != NULL )
              {
                 gapinline=getgapsfromstring(curstr, 0, drw_width, &gapline);
                 empty_line =(Boolean)(gapinline == LINE_ONLYGAP);
                 if (gapinline == LINE_NOGAP)
                    gaplinep = NULL;
                 else
                    gaplinep = gapline;
              }
           }
           break; 
        }
     }
  }
  get_client_rect (p, &rp);
  SelectFont ((FonT)(adp->font));
  ptlh.x = rp.left;
  ptlh.y = rp.top;
  line = 0;
  curalgline = 0;
  adp->visibleLength = drw_width;
  masterbuf = get_master (adp->linebuff, adp->master.entityID, 
                          adp->master.itemID, adp->master.itemtype);
  curvnp = adp->firstssp;

  while ( line < adp->pnlLine )
  {
     curssp = (SelEdStructPtr) curvnp->region;
     if (curssp !=NULL) {
         adp->lastses = curssp;
         itemsubtype = curvnp->itemtype;
         invalid_row = row_is_invalid(adp, ptlh);
         invalid_id  = id_is_invalid (adp, ptlh);
         if (itemsubtype == EDITDEF_SCA) 
         {
                what_inline(adp, line, LINE0, LINE0, LINE0, LINE0, curalgline);
                if ( invalid_row ) {
                   draw_scale (adp, from_inbuf, adp->visibleWidth, &ptlh);
                }
                line++;
                ptlh.y += adp->lineheight;
         }
         else if (itemsubtype == EDITDEF_SCB) 
         {
                what_inline(adp, line, LINE0, LINE0, LINE0, LINE0, curalgline);
                if ( invalid_row ) {
                   draw_bars (adp, from_inbuf, adp->visibleWidth, &ptlh);
                }
                line++;
                ptlh.y += adp->lineheight;
         }
         else if (curssp->itemtype ==OBJ_BIOSEQ && itemsubtype == FEATDEF_BAD)
         {
            if (curssp->data !=NULL) {
               vnp = (ValNodePtr) curssp->data;
               curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
               curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
               if ( curstr != NULL )         
               {
                 gapinline=getgapsfromstring(curstr, 0, drw_width, &gapline);
                 empty_line =(Boolean)(gapinline == LINE_ONLYGAP);
                 if (gapinline == LINE_NOGAP) 
                    gaplinep = NULL;
                 else
                    gaplinep = gapline;
                 if (adp->draw_emptyline || (!adp->draw_emptyline && !empty_line)) {
                   what_inline (adp, line, curssp->entityID, curssp->itemID,  
                                curssp->itemtype, itemsubtype, curalgline);
                   is_master = is_samess_ses (&(adp->master), curssp);
                   draw_diff=!(!adp->charmode || (adp->charmode && is_master));
                   masterstr = NULL;
                   if ( draw_diff ) {
                      stringlens = StringLen (masterbuf);
                      if ( from_inbuf < stringlens )
                            masterstr = masterbuf + from_inbuf;
                   }
                   if ( invalid_id ) {
                     if (adp->displaytype) 
                     {
                        if (adp->marginwithgroup) 
                           group=getparam(adp->params,curssp->entityID, curssp->itemID, 1);
                        else group = -1;
                        sip = SeqLocId((SeqLocPtr) curssp->region);
                        start=AlignCoordToSeqCoord(from_inseq, sip,salp, adp->sqloc_list, 0) +1;
                        stridp = seqid_tolabel (sip, adp->printid);
                        draw_id (adp, &ptlh, rang(adp->input_format, curssp), stridp, curtdp->strand, start, curssp->itemtype, FALSE, is_master, group);
                     }
                   }
                   if ( invalid_row && adp->displaytype )
                   {
                      draw_seq (adp, curssp, from_inseq, drw_width, ptlh, curstr, draw_diff, masterstr, FALSE, curssp->itemtype, TRUE);
                      draw_caret (adp, curssp, ptlh, line);
                   }
                   line++;
                   ptlh.y += adp->lineheight;
                 }
               }
            }
         }
         else if (curssp->itemtype==OBJ_SEQFEAT && itemsubtype==FEATDEF_PROT)
         {
            if (curssp->data !=NULL && !empty_line) 
            {
               vnp = (ValNodePtr) curssp->data;
               curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
               curstr = get_substring(curtdp->buf, from_inbuf, drw_width);
               if ( curstr != NULL )                   
               {
                  what_inline (adp, line, curssp->entityID, curssp->itemID,  curssp->itemtype, itemsubtype, curalgline);
                  draw_id (adp, &ptlh, -1, curtdp->label, curtdp->strand, 0, curssp->itemtype, FALSE, FALSE, -1);
                  if ( invalid_row ) {
                     draw_seq (adp, curssp, from_inseq, drw_width, ptlh, curstr, draw_diff, masterstr, FALSE, curssp->itemtype, FALSE);
                  }
                  line++;
                  ptlh.y += adp->lineheight;
               }
            }
         }
         else if ( itemsubtype == FEATDEF_TRSL ) 
         {
            if (!empty_line) {
               first = TRUE;
               line_drawn = FALSE;
               while ( curssp != NULL && curssp->data != NULL) 
               {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2 = SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord (start2, sip, salp, 0, chklocp);
                  stop2 = SeqLocStop(curslp);
                  chklocp =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord (stop2, sip, salp, 0, chklocp);
                  if (start < from_inseq + drw_width && stop >= from_inseq)
                  {
                     if ( invalid_row ) {
                        draw_pept (adp, curssp, &ptlh, from_inseq, drw_width, line, 0);
                     }
                     if (first) {
                        what_inline (adp, line, curssp->entityID, curssp->itemID, curssp->itemtype, itemsubtype, curalgline);
                        line_drawn = TRUE;
                        first = FALSE;
                     }
                  }
                  curssp = FindNextSegment (curssp);
               }
               if (line_drawn) {
                  line++;
                  ptlh.y += adp->lineheight;
               }
            }
         }
         else if ( curssp->itemtype ==OBJ_SEQFEAT )
         {
            if (!empty_line) {
               first = TRUE;
               line_drawn = FALSE;
               while (curssp!=NULL) 
               {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2 = SeqLocStart(curslp);
                  chklocp = chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord (start2, sip, salp, 0, chklocp);
                  stop2 = SeqLocStop(curslp);
                  chklocp2 = chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord (stop2, sip, salp, 0, chklocp2);
                  if (start <= stop && start < from_inseq + drw_width && stop >= from_inseq)
                  {
                     if ( invalid_row ) 
                     {
                        draw_feat (adp, curssp->entityID, curssp->itemID, curssp->itemtype, itemsubtype, start, stop, SeqLocStrand(curslp), curssp->label, &ptlh, from_inseq, drw_width, line, curalgline, (Uint4) getcolor_fromstyles(itemsubtype), (Boolean)(chklocp==-1), (Boolean)(chklocp2==APPEND_RESIDUE), gaplinep);
                     }
                     if (first) {
                        what_inline(adp, line, curssp->entityID, curssp->itemID,
                                     curssp->itemtype, itemsubtype, curalgline);
                        line_drawn = TRUE;
                        first = FALSE;
                     }
                  }
                  curssp = FindNextSegment (curssp);
               }
               if (line_drawn) {
                  line++;
                  ptlh.y += adp->lineheight;
               }
            }
         }
         else if ( itemsubtype == EDITDEF_CPL ) 
         {
            if (curssp->data !=NULL && !empty_line) 
            {
                vnp = (ValNodePtr) curssp->data;
                curtdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                curstr = get_substring (curtdp->buf, from_inbuf, drw_width);
                if ( curstr != NULL )                         /* !!!!!!! */
                {
                   what_inline (adp, line, curssp->entityID, curssp->itemID,  
                            curssp->itemtype, itemsubtype, curalgline);
                   if ( invalid_id && (( !adp->displaytype 
                   && curssp->itemtype != OBJ_BIOSEQ) || adp->displaytype)) 
                      draw_id (adp, &ptlh, -1, strcmpl, 255, 0, curssp->itemtype, FALSE, FALSE, -1);
                   if ( invalid_row && (( !adp->displaytype 
                   && curssp->itemtype != OBJ_BIOSEQ) || adp->displaytype) ) 
                   {
                      draw_seq (adp, curssp, from_inseq, drw_width, ptlh, curstr, draw_diff, masterstr, TRUE, itemsubtype, FALSE);
                   }
                   line++;
                   ptlh.y += adp->lineheight;
                }
            }
         }
         else if( itemsubtype==SEQFEAT_GENE || itemsubtype==SEQFEAT_RNA
              ||  itemsubtype==SEQFEAT_CDREGION)
         {
            if (!empty_line) {
               first = TRUE;
               line_drawn = FALSE;
               while (curssp!=NULL) 
               {
                  curslp = (SeqLocPtr) curssp->region;
                  sip = SeqLocId (curslp);
                  start2 = SeqLocStart(curslp);
                  chklocp =chkloc(sip, start2, adp->sqloc_list, &start2);
                  start= SeqCoordToAlignCoord (start2, sip, salp, 0, chklocp);
                  stop2 = SeqLocStop(curslp);
                  chklocp2 =chkloc(sip, stop2, adp->sqloc_list, &stop2);
                  stop = SeqCoordToAlignCoord (stop2, sip, salp, 0, chklocp);
                  if (start <= from_inseq + drw_width && stop >= from_inseq)
                  {
                     if ( invalid_row ) 
                     {
                        draw_feat (adp, curssp->entityID, curssp->itemID, curssp->itemtype, itemsubtype, start, stop, SeqLocStrand(curslp), curssp->label, &ptlh, from_inseq, drw_width, line, curalgline, (Uint4) getcolorforvirtfeat(curssp->itemID), (Boolean)(chklocp==-1), (Boolean)(chklocp2==APPEND_RESIDUE), gaplinep);
                     }
                     if (first) {
                        what_inline (adp, line, curssp->entityID, curssp->itemID,
                                     curssp->itemtype, itemsubtype, curalgline);
                        line_drawn = TRUE;
                        first = FALSE;
                     }
                  }
                  curssp = FindNextSegment (curssp);
               }
               if (line_drawn) {
                  line++;
                  ptlh.y += adp->lineheight;
               }
            }
         }
         else if (itemsubtype>=EDITDEF_RF1 && itemsubtype<=EDITDEF_RF6) 
         {
            if (curssp->data !=NULL  && !empty_line) {
               curslp = (SeqLocPtr) curssp->region;
               vnp = (ValNodePtr) curssp->data;
               if ( SeqLocStop (curslp) > from_inseq && vnp != NULL) 
               {
                  what_inline (adp, line, curssp->entityID, curssp->itemID,
                                curssp->itemtype, itemsubtype, curalgline);
                  curstr = (CharPtr) vnp->data.ptrvalue;
                  offset = (Uint2)(itemsubtype - EDITDEF_RF1) % (Uint2)3;
                  if (itemsubtype>=EDITDEF_RF1 && itemsubtype<EDITDEF_RF4)
                     draw_id (adp, &ptlh, rang(adp->input_format, curssp), NULL,  Seq_strand_plus, 0, curssp->itemtype, FALSE, FALSE, -1);
                  else
                     draw_id (adp, &ptlh, rang(adp->input_format, curssp), NULL, Seq_strand_minus, 0, curssp->itemtype, FALSE, FALSE, -1);

                  if ( invalid_row ) 
                  {
                    draw_trans (adp, curssp, curstr, &ptlh, from_inseq, drw_width, offset);
                  }
               }
               line++;
               ptlh.y += adp->lineheight;
            }
         }
         else {
            ErrPostEx (SEV_ERROR, 0, 0, "fail in on_draw [46] subtype %ld type %ld", itemsubtype, curssp->itemtype); 
            return;
         }
         if (line == adp->pnlLine) 
            break;
         curvnp = curvnp->next;
         if (curvnp == NULL) 
         {
            if (adp->rowpcell > 0 
            && ((Int2)(curalgline+1) % (Int2)adp->rowpcell) == 0) 
            {
                what_inline(adp,line,LINE0,LINE0, LINE0, LINE0, curalgline);
                line++;
                ptlh.y += adp->lineheight;
            }
            curvnp = adp->buffer;
            curalgline++;
            from_inbuf += drw_width;
            from_inseq += drw_width;
            if (from_inseq >= adp->length) 
                break;
            drw_width = MIN ((Int2) adp->visibleWidth, (Int2)(adp->bufferlength 
                      - (adp->hoffset -adp->bufferstart) - adp->visibleLength));
            adp->visibleLength += drw_width;
         }
     }
  }
  MemFree (gapline);
  SetMuskCurrentSt (GetMuskStyleName (oldstyle));
  return;
}
