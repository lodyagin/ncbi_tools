/* $Id: seqpanel.c,v 6.11 2002/12/02 16:00:10 kans Exp $
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
* Author:  Vlad Lebedev
**
* ==========================================================================
*/

#include <bspview.h>
#include <seqpanel.h>
/*
#include <samutil.h>
*/
#include <objmgr.h>
#include <explore.h>

enum ESeqNum   { eNumNone=1, eNumSide=2, eNumTop=3 };
enum ELineType { eTypeTopSeqNumbers, eTypeTopScaleMarks, eTypeSequence, eTypeFeature };
enum EDrawGrid { eDrawGridOn=1, eDrawGridOff=2 };
enum EShowFeatures { eShowFeaturesOn=1, eShowFeaturesOff=1 };

#define GRID_LINE_SPACE       6
#define NO_GRID_LINE_SPACE    2
#define SEQ_GROUP_SIZE       10  /* Sequence group size           */
#define SEQ_START_POS       100  /* Draw Sequence from this x pos */
#define SEQ_X_OFFSET          4
#define SEQ_Y_OFFSET          4
#define PROT_PRODUCT_TYPE   255

static Uint1 FillRectangleSym [] = { 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0x00 };
static Uint1 FillLeftArrowSym [] = { 0x06, 0x1E, 0x7E, 0xFE, 0x7E, 0x1E, 0x06, 0x00 };
static Uint1 FillRightArrowSym[] = { 0xC0, 0xF0, 0xFC, 0xFE, 0xFC, 0xF0, 0xC0, 0x00 };

typedef struct seqParag { 
  ValNodePtr pFeatList; /* list of feature indexes */ 
} SeqParaG, PNTR SeqParaGPtr;


static BioseqViewPtr GetBioseqViewPtr(PaneL p)
{
  return &((BioseqViewFormPtr) GetObjectExtra (p))->bvd;
}

static void FreeSeqPanelLines(SeqPanLinePtr PNTR splp, BioseqViewPtr bvp)
{
  Int4 i;
  for (i=0; i<bvp->TotalLines; i++) MemFree(splp[i]);
  MemFree (splp);
}


static SeqPanLinePtr MakeSeqPanLine(Int2 type, Int4 line)
{
  SeqPanLinePtr plp = (SeqPanLinePtr)MemNew(sizeof(SeqPanLine));
  plp->lineType     = type;
  plp->bioSeqLine   = line;
  return plp;
}

static void ShowSeqView (BioseqViewPtr bvp, Boolean show)
{
  if (bvp == NULL) return;
  show ? SafeShow (bvp->seqView         ) : SafeHide (bvp->seqView         );
  show ? SafeShow (bvp->seqViewParentGrp) : SafeHide (bvp->seqViewParentGrp);
  show ? SafeShow (bvp->clickMe         ) : SafeHide (bvp->clickMe         );
  SafeHide (bvp->styleControlGrp);
  SafeHide (bvp->scaleControlGrp);
  SafeHide (bvp->findGeneGrp    );
}

static SeqPanLinePtr PNTR CreateSeqPanelLines(Int2 lineLength, BioseqViewPtr bvp)
{
  BioseqPtr          bsp;
  SeqParaGPtr   PNTR ref;
  SeqPanLinePtr PNTR splp;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp;
  ValNodePtr         vnp;
  Int2               fLines = GetValue(bvp->newNumControl)==eNumTop ? 3 : 1;  
  Int4               lCount = 0, i, pCount;
  
  bsp    = bvp->bsp;
  pCount = floor(bsp->length / lineLength) + 1; /* Total number of paragraphs */
  
  ref = (SeqParaGPtr*) MemNew( (size_t)(sizeof(SeqParaGPtr)*pCount) );
  for (i=0; i<pCount; i++) ref[i] = (SeqParaGPtr)MemNew(sizeof(SeqParaG));  

  if (GetValue(bvp->newFeatControl)==eShowFeaturesOn) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
    while (sfp != NULL) {
      if (fcontext.seqfeattype!=SEQFEAT_PUB  &&  fcontext.seqfeattype!=SEQFEAT_BIOSRC) {
        Boolean   coding   = fcontext.seqfeattype == SEQFEAT_CDREGION;
        Int4      paraFrom = floor (fcontext.left  / lineLength);      /* feature starting paragraph */
        Int4      paraTo   = ceil  (fcontext.right / lineLength) + 1;  /* feature ending paragraph   */
        /*BioseqPtr bsp_prot = BioseqFind (SeqLocId(sfp->product));*/
        
        for (i=paraFrom; i!=paraTo && i<pCount; i++) {
          ValNodeAddInt(&ref[i]->pFeatList, fcontext.seqfeattype, fcontext.itemID);
          if (coding) ValNodeAddInt(&ref[i]->pFeatList, PROT_PRODUCT_TYPE, fcontext.itemID); /* add space for prot product */
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);      
    }
  }  

  bvp->TotalLines = 0; /* go through all pararaphs and count total */
  for (i=0; i<pCount; i++) bvp->TotalLines += fLines+ValNodeLen (ref[i]->pFeatList); 

  splp = (SeqPanLinePtr*) MemNew( (size_t)(sizeof(SeqPanLinePtr)*bvp->TotalLines) );
  for (i=0; i<pCount; i++) {
    if (fLines==3) {
      splp[lCount++] = MakeSeqPanLine(eTypeTopSeqNumbers, i);
      splp[lCount++] = MakeSeqPanLine(eTypeTopScaleMarks, i);
    }

    splp[lCount++] = MakeSeqPanLine(eTypeSequence, i);         /* Add sequence line      */
    
   	for (vnp=ref[i]->pFeatList; vnp!=NULL; vnp=vnp->next) {
        SeqPanLinePtr plp = MakeSeqPanLine(eTypeFeature, i);   /* Add feature line       */
        plp->idx          = vnp->data.intvalue;
        plp->protProduct  = vnp->choice==PROT_PRODUCT_TYPE ? TRUE : FALSE;
        splp[lCount++]    = plp;
    }
  }

  for (i=0; i<pCount; i++) {
    ValNodeFree (ref[i]->pFeatList);
    MemFree (ref[i]);
  }
  MemFree (ref);
  return splp;
}

static void onSeqViewClick (PaneL p, PoinT pt)
{
  BioseqViewPtr bvp;

  bvp  = GetBioseqViewPtr (p);
  bvp->wasDoubleClick = dblClick;
  bvp->wasShiftKey = shftKey;
  bvp->old_rect_shown = FALSE;
  bvp->pnt_start = pt;
  bvp->pnt_stop = pt;
}


static void onSeqViewRelease (PaneL p, PoinT pt)
{
  RecT          r;
  BaR           sb;
  BioseqPtr     bsp;
  BioseqViewPtr bvp;
  SeqPanLinePtr splp;
  Int4          line;
  Uint2         entityID = 0;
  Uint4         itemID = 0;
  Uint2         itemtype = 0;
  SeqEntryPtr   sep;

  bvp  = GetBioseqViewPtr (p);
  bsp  = bvp->bsp;
  sb   = GetSlateVScrollBar ((SlatE)bvp->seqView);

  ObjectRect (bvp->seqView, &r);
  InsetRect (&r, 4, 4);

  line = (pt.y-r.top-SEQ_Y_OFFSET) / bvp->LineHeight + GetBarValue(sb);
  if(line>=bvp->TotalLines) return;
  
  splp = bvp->SeqPanLines[line];
  switch (splp->lineType)  {
    case eTypeSequence:
      /*
      printf("Sequence from %d to %d\n", splp->bioSeqLine*bvp->CharsAtLine+1, splp->bioSeqLine*bvp->CharsAtLine+bvp->CharsAtLine);
      */
      break;
    case eTypeFeature:
      /*
      if (splp->protProduct) printf("Prot. product: %d\n", splp->idx);
      else printf("Feature. Index: %d\n", splp->idx);
      */
      if (splp->protProduct) {
        /* just get CDS feature from protein product */
        entityID = ObjMgrGetEntityIDForPointer (bsp);
        itemID = splp->idx;
        itemtype = OBJ_SEQFEAT;
      } else {
        entityID = ObjMgrGetEntityIDForPointer (bsp);
        itemID = splp->idx;
        itemtype = OBJ_SEQFEAT;
      }
      break;
  }
  if (itemID < 1) {
    if (! bvp->wasShiftKey) {
      ObjMgrDeSelect (0, 0, 0, 0, NULL);
    }
    return;
  }

  if (bvp->wasDoubleClick) {
      sep = GetTopSeqEntryForEntityID (entityID);
      if (bvp->launchSubviewers) {
        WatchCursor ();
        Update ();
        LaunchNewBioseqViewer (bvp->bsp, entityID, itemID, itemtype);
        ArrowCursor ();
        Update ();
        return;
      } else if (LaunchViewerNotEditor (bvp, sep, entityID, itemID, itemtype)) {
        WatchCursor ();
        Update ();
        LaunchNewBioseqViewer (bvp->bsp, entityID, itemID, itemtype);
        ArrowCursor ();
        Update ();
        return;
      } else if (bvp->launchEditors) {
        WatchCursor ();
        Update ();
        GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, itemID,
                          itemtype, 0, 0, itemtype, 0);
        ArrowCursor ();
        Update ();
        return;
      } else {
        return;
      }
    if (! bvp->sendSelectMessages) return;
    if (bvp->wasShiftKey) {
      ObjMgrAlsoSelect (entityID, itemID, itemtype, 0, NULL);
    } else {
      ObjMgrSelect (entityID, itemID, itemtype, 0, NULL);
    }
  } else if (bvp->sendSelectMessages) {
      if (bvp->wasShiftKey) {
        ObjMgrAlsoSelect (entityID, itemID, itemtype, 0, NULL);
      } else {
        ObjMgrSelect (entityID, itemID, itemtype, 0, NULL);
      }
  }
}


static void ResizeSeqView (BioseqViewPtr bvp)
{
  RecT r;
  Int4 height, width;
  BaR  sb = GetSlateVScrollBar ((SlatE)bvp->seqView);

  ObjectRect (bvp->seqView, &r);
  InsetRect (&r, 4, 4);
  width  = r.right  - r.left;
  height = r.bottom - r.top - 2 * SEQ_Y_OFFSET;

  bvp->BlocksAtLine = (width - SEQ_X_OFFSET - SEQ_START_POS) / ((SEQ_GROUP_SIZE + 1) * bvp->CharWidth);
  if (bvp->BlocksAtLine == 0) bvp->BlocksAtLine = 1; /* Always at least 1 block */

  bvp->CharsAtLine = SEQ_GROUP_SIZE * bvp->BlocksAtLine;

  if (bvp->SeqPanLines) FreeSeqPanelLines (bvp->SeqPanLines, bvp);
  bvp->SeqPanLines = CreateSeqPanelLines (bvp->CharsAtLine, bvp);

  SetBarMax (sb, bvp->TotalLines - height / bvp->LineHeight);
  CorrectBarPage(sb, (height / bvp->LineHeight) - 1, (height / bvp->LineHeight) - 1);

  SetPanelClick(bvp->seqView, onSeqViewClick, NULL, NULL, onSeqViewRelease);
}


static void PopulateSeqView (BioseqViewPtr bvp)
{
  RecT r;

  SelectFont ((FonT)(bvp->displayFont));

  bvp->DrawGrid   = GetValue (bvp->newGridControl) == eDrawGridOn;
  bvp->LineSpace  = bvp->DrawGrid ? GRID_LINE_SPACE-1 : NO_GRID_LINE_SPACE;
  bvp->CharHeight = FontHeight ();
  bvp->CharWidth  = CharWidth ('A');
  bvp->LineHeight = bvp->CharHeight + bvp->LineSpace;

  CorrectBarValue (GetSlateVScrollBar ((SlatE)bvp->seqView), 0);
  ResizeSeqView (bvp);  

  Select (bvp->seqView);
  ObjectRect (bvp->seqView, &r);
  InsetRect (&r, 2, 2);
  InvalRect  (&r);

  SetPanelClick(bvp->seqView, onSeqViewClick, NULL, NULL, onSeqViewRelease);
}



static void onCloseSeqPanel (PaneL p)
{
  BioseqViewPtr bvp;
  bvp  = GetBioseqViewPtr (p);

  if (bvp->SeqPanLines) FreeSeqPanelLines (bvp->SeqPanLines, bvp);
  bvp->SeqPanLines = NULL;
}


static void DrawTopScaleMarks(Int2 x, Int2 y, Int4 line, BioseqViewPtr bvp)
{
  Int2 block, ctr=0;
  
  Magenta ();
  for (block=0;  block!=bvp->BlocksAtLine  &&  ctr>=0;  block++) {
    MoveTo(x+SEQ_X_OFFSET+SEQ_START_POS+(block+1)*SEQ_GROUP_SIZE*bvp->CharWidth+block*bvp->CharWidth - bvp->CharWidth/2 - 1, y);
    LineTo(x+SEQ_X_OFFSET+SEQ_START_POS+(block+1)*SEQ_GROUP_SIZE*bvp->CharWidth+block*bvp->CharWidth - bvp->CharWidth/2 - 1, y-bvp->CharHeight);

    MoveTo(x+SEQ_X_OFFSET+SEQ_START_POS+(block+1)*SEQ_GROUP_SIZE*bvp->CharWidth+block*bvp->CharWidth - SEQ_GROUP_SIZE*bvp->CharWidth/2 - bvp->CharWidth/2 - 1, y);
    LineTo(x+SEQ_X_OFFSET+SEQ_START_POS+(block+1)*SEQ_GROUP_SIZE*bvp->CharWidth+block*bvp->CharWidth - SEQ_GROUP_SIZE*bvp->CharWidth/2 - bvp->CharWidth/2 - 1, y-bvp->CharHeight/2);
  }
}

static void DrawTopSeqNums(Int2 x, Int2 y, Int4 line, BioseqViewPtr bvp)
{
  Int2 block, ctr=0;
  char buf[20];
  
  Magenta ();
  for (block=0;  block!=bvp->BlocksAtLine  &&  ctr>=0;  block++) {
    sprintf(buf, "%d", line * bvp->CharsAtLine + (block+1)*SEQ_GROUP_SIZE);
    PaintStringEx (buf, x+SEQ_X_OFFSET+SEQ_START_POS+(block+1)*SEQ_GROUP_SIZE*bvp->CharWidth+block*bvp->CharWidth - bvp->CharWidth*StrLen(buf), y);
  }
}


static void DrawSideLineNumbers(Int2 x, Int2 y, Int4 line, BioseqViewPtr bvp)
{
  char buf[20];
  sprintf(buf, "%d", line * bvp->CharsAtLine + 1);

  Magenta ();
  PaintStringEx (buf, x+ SEQ_X_OFFSET + (SEQ_START_POS-30-bvp->CharWidth*StrLen(buf)), y);
}

static void DrawSequence(Int2 x, Int2 y, Int4 line, SeqPortPtr spp, Uint1Ptr buf, BioseqViewPtr bvp)
{
  Int2 block, ctr=0, i;
  Char ch;
  Uint1Ptr ptr;

  Black ();
  for (block=0;  block!=bvp->BlocksAtLine  &&  ctr>=0;  block++) {
    SeqPortSeek (spp, line*bvp->CharsAtLine+block*SEQ_GROUP_SIZE, SEEK_SET);
    ctr = SeqPortRead (spp, buf, SEQ_GROUP_SIZE);
    if (ctr>0) {
      for (ptr = buf, i = 0; i < ctr; ptr++, i++) {
        ch = (Char) *ptr;
        ch = TO_LOWER (ch);
        *ptr = (Uint1) ch;
      }
      if (ctr<SEQ_GROUP_SIZE) MemSet(buf+ctr, '\0', 1);
      PaintStringEx ( (char*)buf, x+SEQ_X_OFFSET+SEQ_START_POS+block*SEQ_GROUP_SIZE*bvp->CharWidth+block*bvp->CharWidth, y);
    }
  }
}

static void DrawLineEx(Int2 x1, Int2 y1, Int2 x2, Int2 y2, Int2 width)
{
  WidePen (width);
  MoveTo  (x1, y1); LineTo  (x2, y2);
  WidePen (1);
}


static void DrawLtGrid(Int2 x1, Int2 y1, Int2 x2, Int2 y2) { LtGray(); MoveTo (x1, y1); LineTo (x2, y2); }
static void DrawDkGrid(Int2 x1, Int2 y1, Int2 x2, Int2 y2) { DkGray(); MoveTo (x1, y1); LineTo (x2, y2); }


static Int2 SeqPos2XCoord(Int2 x, Int4 seqXPos, BioseqViewPtr bvp)
{
  Int4 blocks = seqXPos / SEQ_GROUP_SIZE;
  Int2 xPos = x + SEQ_X_OFFSET + SEQ_START_POS + seqXPos * bvp->CharWidth + blocks*bvp->CharWidth;
  Int2 rPos = x + SEQ_X_OFFSET + SEQ_START_POS + bvp->CharsAtLine * bvp->CharWidth + (blocks-1) * bvp->CharWidth;
  return xPos > rPos ? rPos : xPos;
}

static Int4 GetFeatureX(Int2 x, Int4 feat_pos, Int4 bsStart, Int4 bsFinish, Boolean is_start)
{
  Int4 lineFeatPos = is_start ? (bsStart  < feat_pos ? feat_pos : bsStart ): 
                                (bsFinish < feat_pos ? bsFinish : feat_pos);
  return lineFeatPos;
}


static Boolean IsInRange(Int4 pos, Int4 min_pos, Int4 max_pos)
{
  return min_pos <= pos  &&  pos <= max_pos;
}


static void DrawFeature(Int2 x, Int2 y, Int4 line, Int2 itemID, Boolean protProduct, BioseqPtr bsp, BioseqViewPtr bvp)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  BioseqPtr         bsp_prot;
  CharPtr           str_prot;
  Int4              bsStart, bsFinish, ffStart, ffFinish, i;
  Int2              x1, x2;
  Char              featLabel[13];
    
  sfp      = SeqMgrGetDesiredFeature (0, bsp, itemID, 0, NULL, &fcontext);
  bsStart  = line    * bvp->CharsAtLine;
  bsFinish = bsStart + bvp->CharsAtLine;
  ffStart  = GetFeatureX(x, fcontext.left,  bsStart, bsFinish, TRUE ) - bsStart;
  ffFinish = GetFeatureX(x, fcontext.right, bsStart, bsFinish, FALSE) - bsStart;
  x1       = SeqPos2XCoord(x, ffStart,  bvp);
  x2       = SeqPos2XCoord(x, ffFinish, bvp);

  fcontext.seqfeattype == SEQFEAT_CDREGION ? Blue() : Black();
  StrNCpy( featLabel, fcontext.label, 12); featLabel[12]='\0';
  PaintStringEx (featLabel, x+10, y);

  if (protProduct) {
    bsp_prot = BioseqFind (SeqLocId(sfp->product));
    str_prot = GetSequenceByBsp (bsp_prot);
    if (str_prot==NULL) {
      PaintStringEx ("Protein sequence in not avaliable", SeqPos2XCoord(x, 0, bvp), y);
      return;
    }
    for (i=ffStart; i!=ffFinish; i++) PaintStringEx ("~", SeqPos2XCoord(x, i, bvp), y);
  }
  else DrawLineEx(x1, y-bvp->LineHeight/2+2, x2, y-bvp->LineHeight/2+2, 1);
  
  for (i=0; i!=fcontext.numivals; i++) {
    RecT rect;
    Int4 regStart  = fcontext.strand==Seq_strand_minus ? fcontext.ivals[i*2+1] : fcontext.ivals[i*2];
    Int4 regFinish = fcontext.strand==Seq_strand_minus ? fcontext.ivals[i*2] : fcontext.ivals[i*2+1];
    Int4 fStart    = GetFeatureX(x, regStart,  bsStart, bsFinish, TRUE );
    Int4 fFinish   = GetFeatureX(x, regFinish, bsStart, bsFinish, FALSE);
    
    if( !IsInRange(regStart, bsStart,  bsFinish ) && !IsInRange(regFinish, bsStart,  bsFinish) &&
        !IsInRange(bsStart,  regStart, regFinish) && !IsInRange(bsFinish,  regStart, regFinish) ) continue;       

    x1 = SeqPos2XCoord(x, fStart  - bsStart, bvp);
    x2 = SeqPos2XCoord(x, fFinish - bsStart, bvp);
    
    if (protProduct) { /* draw protein product */
      Int4 j, k, offPos;
      Int4 prodPos  = fStart - regStart;           /* offset in the array of protein product (letter) */
      Int4 prodPos2 = fStart - fcontext.left + 2;  /* offset in the screen position of the product    */
      Char tmp[2];
      
      if (prodPos > 0) prodPos++;
      for (k=0; k!=i; k++) { 
        prodPos += (fcontext.strand==Seq_strand_minus ? (fcontext.ivals[k*2]   - fcontext.ivals[k*2+1]) :
                                                        (fcontext.ivals[k*2+1] - fcontext.ivals[k*2])) + 1;
      }
      prodPos = prodPos / 3;
      
      if (prodPos2 % 3==0) offPos = 0;
      else offPos = (prodPos2 / 3 + 1) * 3 - prodPos2;

#ifdef WIN_MOTIF
      Gray();
      DrawLineEx(x1, y-1, x2, y-1, bvp->CharHeight);
#else
      White();
      DrawLineEx(x1, y-bvp->CharHeight, x2, y-bvp->CharHeight, bvp->CharHeight);
#endif 
      

      for (j=fStart; j<fFinish; j+=3) {
        Int4 xPos = j+offPos-bsStart;
        
        if (xPos>=bvp->CharsAtLine  ||  xPos>regFinish-bsStart) break;
        StrNCpy( tmp, &str_prot[prodPos++], 1); tmp[1]='\0';
        Blue();
        PaintStringEx ( tmp, SeqPos2XCoord(x, xPos, bvp), y);
      }
    } 
    else { /* draw feature line */

#ifdef WIN_MOTIF
      DrawLineEx(x1+1, y-bvp->LineHeight/2+2, x2, y-bvp->LineHeight/2+2, 3);
#else
      DrawLineEx(x1, y-bvp->LineHeight/2+1, x2, y-bvp->LineHeight/2+1, 3);
#endif      

      if (IsInRange(regStart, bsStart, bsFinish)) {
        LoadRect (&rect, x1, y-bvp->LineHeight/2-1, x1+7, y-bvp->LineHeight/2+6);
        CopyBits (&rect, fcontext.strand==Seq_strand_minus ? FillLeftArrowSym : FillRectangleSym );
      }
      if (IsInRange(regFinish, bsStart, bsFinish)) {
        LoadRect (&rect, x2, y-bvp->LineHeight/2-1, x2+7, y-bvp->LineHeight/2+6);
        CopyBits (&rect, fcontext.strand==Seq_strand_minus ? FillRectangleSym : FillRightArrowSym);
      }
    } /* feature or product */
  } /* for */
}


static void onDrawSeqPanel (PaneL p)
{
  BioseqViewPtr bvp;
  BioseqPtr     bsp;
  SeqPortPtr    spp;
  SeqPanLinePtr splp;
  Uint1Ptr      buf;
  BaR           sb;
  RecT          r;
  Int4          line;
  Int2          x, y;
 
  bvp = GetBioseqViewPtr (p);
  bsp = bvp->bsp;
  spp = SeqPortNew (bsp, 0, bsp->length-1, Seq_strand_plus, Seq_code_iupacna);
  buf = MemNew (SEQ_GROUP_SIZE+1);
  sb  = GetSlateVScrollBar ((SlatE)bvp->seqView);
  
  ObjectRect (p, &r);
  InsetRect (&r, 4, 4);
  x = r.left + 1;
  y = r.top  + bvp->CharHeight + SEQ_Y_OFFSET;
  
  SelectFont ((FonT)(bvp->displayFont));  
  for (line=GetBarValue(sb); line<bvp->TotalLines  &&  y <= r.bottom-2*SEQ_Y_OFFSET; line++) {
    if (IsInRange(y, updateRect.top,updateRect.bottom) || 
        IsInRange(y+bvp->LineHeight,updateRect.top,updateRect.bottom))
    {
      /* draw begin */
      splp = bvp->SeqPanLines[line];
      switch ( splp->lineType ) {
        case eTypeTopSeqNumbers:
  	      if (bvp->DrawGrid) DrawLtGrid(SEQ_START_POS+1, y+bvp->LineSpace/2, r.right, y+bvp->LineSpace/2);
          DrawTopSeqNums(x, y, splp->bioSeqLine, bvp);                               /* Draw top numbering   */
          break;
        case eTypeTopScaleMarks:
    	    if (bvp->DrawGrid) DrawLtGrid(SEQ_START_POS+1, y+bvp->LineSpace/2, r.right, y+bvp->LineSpace/2);
          DrawTopScaleMarks(x, y, splp->bioSeqLine, bvp);                            /* Draw top scale marks */
          break;
        case eTypeSequence:
          if (GetValue(bvp->newNumControl)!=eNumNone) DrawSideLineNumbers(x, y, splp->bioSeqLine, bvp); /* Draw line numbers    */
          DrawSequence(x, y, splp->bioSeqLine, spp, buf, bvp);                       /* Draw the sequence    */
          if (bvp->DrawGrid) DrawLtGrid(x, y+bvp->LineSpace/2, r.right, y+bvp->LineSpace/2);
          break;
        case eTypeFeature:
          if (bvp->DrawGrid) DrawLtGrid(x, y+bvp->LineSpace/2, r.right, y+bvp->LineSpace/2);
          DrawFeature(x, y, splp->bioSeqLine, splp->idx, splp->protProduct, bsp, bvp);    /* Draw Features        */
          break;
      }
      if(bvp->DrawGrid && line<bvp->TotalLines-1 && splp->bioSeqLine!=bvp->SeqPanLines[line+1]->bioSeqLine)
          DrawDkGrid(x, y+bvp->LineSpace/2, r.right, y+bvp->LineSpace/2);             /* Draw Horizontal Grid */
      /* draw end */
    }
    y += bvp->LineHeight;
  }

  if (bvp->DrawGrid) DrawDkGrid(SEQ_START_POS, r.top, SEQ_START_POS, r.bottom);       /* Draw Vertical Grid   */

  MemFree (buf);
  SeqPortFree (spp);
}

static void onVScrollBarSeqPanel (BaR sb, SlatE s, Int4 newval, Int4 oldval)
{
  BioseqViewPtr bvp;
  RecT r_scroll, r_redraw;
  Int2 pixels;
  Int2 lines;
  
  Select (s);
  bvp = GetBioseqViewPtr((PaneL)s);
  ObjectRect (bvp->seqView, &r_scroll);
  ObjectRect (bvp->seqView, &r_redraw);
  InsetRect (&r_scroll, 2, 2);
  InsetRect (&r_redraw, 2, 2);

  lines = (r_scroll.bottom - r_scroll.top - 2*SEQ_Y_OFFSET-3) / bvp->LineHeight;
  
  if (abs(oldval-newval)>lines) InvalRect  (&r_redraw);
  else {
    pixels = (oldval - newval) * bvp->LineHeight;

    if (pixels<0) {
      r_scroll.top   += SEQ_Y_OFFSET + 2;
      r_redraw.bottom = r_scroll.bottom = r_scroll.top + lines * bvp->LineHeight + 6;
      r_redraw.top    = r_scroll.bottom-abs(pixels)-bvp->LineHeight;
    }
    else {
      r_redraw.bottom = r_redraw.top + pixels + SEQ_Y_OFFSET;
      r_scroll.bottom = r_scroll.top + lines * bvp->LineHeight + 6;
    }
    ScrollRect (&r_scroll, 0, pixels);
    InvalRect (&r_redraw);
  }
  
  Update ();
}

/* extern functions */

PaneL CreateSeqViewPanel (GrouP g, Int2 w, Int2 h)
{
  PaneL pnl = AutonomousPanel4 (g, w, h, onDrawSeqPanel, onVScrollBarSeqPanel, NULL, sizeof (BioseqViewPtr), onCloseSeqPanel, NULL);
  SetPanelClick(pnl, onSeqViewClick, NULL, NULL, onSeqViewRelease);
  return pnl;
}

BioseqPageData seqpnlPageData = {
  "Sequence", TRUE, TRUE, TRUE, FALSE, -1,
  PopulateSeqView, ShowSeqView, NULL,
  NULL, NULL,
  NULL, NULL, ResizeSeqView, NULL
};

