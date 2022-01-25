/*   salparam.c
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
* File Name:  salparam.c
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.29 $
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
#include <salparam.h>
#include <saledit.h>
#include <salutil.h>
#include <salstruc.h>
#include <salsap.h>
#include <salpanel.h>
#include <salfiles.h>
#include <saldist.h>
#include <saled.h>
#include <dotmatrx.h>
#include <urkptpf.h>
#include <dlogutil.h>
#include <biosrc.h>

#define OBJ_VIRT 254

static void DotPlotBtn (ButtoN b);

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

/*********************************
***    Change Font Procedures
***********************************/
extern void SeqFontProc (IteM i)
{
  WindoW             temport;        
  SeqEditViewFormPtr wdp;
  EditAlignDataPtr   adp;
  FonT               f;
  Nlm_FontSpec       font;

  wdp = (SeqEditViewFormPtr) GetObjectExtra ((WindoW)ParentWindow (i));
  if ( ( adp = GetAlignDataPanel (wdp->pnl) ) == NULL ) return;
  GetFontSpec ((FonT)(adp->font), &font);
  if ( ChooseFont(&font, CFF_MONOSPACE, NULL)) 
  {
     f = CreateFont(&font);
     adp->font = (Handle)f; 
     SelectFont (f);
     adp->lineheight = LineHeight () + adp->interline;
     adp->leading = Leading ();
     adp->charw   = CharWidth ('0');
     adp->ascent  = Ascent ();
     adp->scaleheight  = 2 * adp->lineheight;
     adp->margin.left  = adp->marginleft * adp->charw;
     temport = SavePort((WindoW)ParentWindow (i));
     Select (wdp->pnl);
     do_resize_window (wdp->pnl, adp, FALSE);
     inval_panel (wdp->pnl, -1 ,-1);
     RestorePort (temport);
  } 
  Update ();
  return;
}

static void ChangeFontButton (ButtoN b)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  FonT             f;
  Nlm_FontSpec     font;

  wdialog = ParentWindow (b);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  GetFontSpec ((FonT)(adp->font), &font);
  if ( ChooseFont(&font, CFF_READ_FSP, NULL)) 
  {
         f = CreateFont(&font);
         adp->newstyle->font = (Handle)f;
  } 
  return;
}

static void ChangeInterlineButton (TexT text)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Char             str[16];

  wdialog = ParentWindow (text);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  GetTitle (text, str, 16);
  adp->newstyle->interline = (Int2) atoi (str);
  return;
}

static void ChangeMarginIndexButton (ButtoN bn)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  adp->newstyle->marginwithindex = GetStatus (bn);
  return;
}

static void ChangeMarginIdButton (ButtoN bn)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  adp->newstyle->marginwithIds = GetStatus (bn);
  return;
}

static void ChangeMarginPosButton (ButtoN bn)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  adp->newstyle->marginwithpos = GetStatus (bn);
  return;
}

static void DrawScaleButton (ButtoN bn)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  adp->newstyle->draw_scale = GetStatus (bn);
  return;
}

static void DrawBarScaleButton (ButtoN bn)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  adp->newstyle->draw_bars = GetStatus (bn);
  return;
}

static void RowCellButton (TexT text)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Int2             val;
  Char             str[16];

  wdialog = ParentWindow (text);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  GetTitle (text, str, 16);
  val = (Int2) atoi (str);
  if (val > 50) {
         val = 50;
         SetTitle (text, "40");
  }
  adp->newstyle->rowpcell = (Uint1) val;
  return;
}

static void ColumnCellButton (TexT text)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Int2             val;
  Char             str[16];

  wdialog = ParentWindow (text);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  GetTitle (text, str, 16);
  val = (Int2) atoi (str);
  if (val > 50) {
         val = 50;
         SetTitle (text, "40");
  }
  adp->newstyle->columnpcell = (Uint1)val;
  return;
}
static void ChangeMarginButton (TexT text)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Int2             margin;
  Char             str[16];

  wdialog = ParentWindow (text);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  GetTitle (text, str, 16);
  margin = (Int2) atoi (str);
  if (margin > 40) {
         margin = 40;
         SetTitle (text, "40");
  }
  adp->newstyle->marginleft = margin;
  return;
}

static void ChangeMarginPopup (PopuP p)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Int2             val;
  Char             str[8];

  wdialog = ParentWindow (p);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  switch ( (val = GetValue (p)) ) {
      case 2:
         adp->marginleft = SALSA_PHYLIP;
         break;
      case 3:
         adp->marginleft = SALSA_CLUSTALV;
         break;
      case 4:
         adp->marginleft = SALSA_MACAW;
         break;
      default:
         break;
  }
  sprintf (str, "%d", (int) adp->marginleft);
  SetTitle (dbdp->txt5, str);
  adp->margin.left = adp->marginleft * adp->charw;
  return;
}

/*********************************
***    SelectColorProc
***********************************/

extern void setUpLetterColor ( Uint4 *colorRefs, Uint1 letter,
                               Uint1 red, Uint1 green, Uint1 blue )
{
  colorRefs[letter] = GetColorRGB(red,green,blue);
}

static void SelectColorScaleButton (ButtoN bn)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Uint1            r, g, b;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  if ( ChooseColorDialog (&r, &g, &b, 0) ) {
     setUpLetterColor ( adp->newstyle->colorRefs, (Uint1)COLOR_SCALE, r, g, b);
  }
  return;
}

static void SelectColorIdButton (ButtoN bn)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Uint1            r, g, b;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  if ( ChooseColorDialog (&r, &g, &b, 0) ) {
      setUpLetterColor ( adp->newstyle->colorRefs, (Uint1)COLOR_ID, r, g, b );
  }
  return;
}

static void SelectColorIdMasterButton (ButtoN bn)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Uint1            r, g, b;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  if ( ChooseColorDialog (&r, &g, &b, 0) ) {
    setUpLetterColor(adp->newstyle->colorRefs, (Uint1)COLOR_ID_MASTER, r, g, b);
  }
  return;
}

static void SelectColorSelectionButton (ButtoN bn)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Uint1            r, g, b;

  wdialog = ParentWindow (bn);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  if ( ChooseColorDialog (&r, &g, &b, 0) ) {
     setUpLetterColor (adp->newstyle->colorRefs, (Uint1)COLOR_SELECT, r, g, b);
  }
  return;
}

static void SelectAaColorPopup (PopuP p)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr     adp;
  Int2             val;
  Uint1            r, g, b;

  wdialog = ParentWindow (p);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  if ( ChooseColorDialog (&r, &g, &b, 0) ) {
    switch ( ( val = GetValue(p) ) ) {
     case 1:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('A' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('C' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('D' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('E' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('F' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('G' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('H' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('I' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('K' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('L' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('M' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('N' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('P' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('Q' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('R' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('S' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('T' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('V' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('W' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('Y' - '*'), r, g, b);
        break;
     case 2:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('M' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('F' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('L' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('V' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('I' - '*'), r, g, b);
        break;
     case 3:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('R' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('H' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('E' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('K' - '*'), r, g, b);
        break;
     case 4:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('T' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('N' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('P' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('S' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('D' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('Q' - '*'), r, g, b);
        break;
     case 5:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('G' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('C' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('A' - '*'), r, g, b);
        break;
     case 6:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('*' - '*'), r, g, b);
        break;
     default:
        break;
    }
  }
  return;
}

static void SelectAnColorPopup (PopuP p)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr     adp;
  Int2             val;
  Uint1            r, g, b;

  wdialog = ParentWindow (p);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  if ( ChooseColorDialog (&r, &g, &b, 0) ) {
    switch ( ( val = GetValue(p) ) ) {
     case 1:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('a' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('c' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('t' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('g' - '*'), r, g, b);
        break;
     case 2:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('a' - '*'), r, g, b);
        break;
     case 3:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('c' - '*'), r, g, b);
        break;
     case 4:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('g' - '*'), r, g, b);
        break;
     case 5:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('t' - '*'), r, g, b);
        break;
     case 6:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('g' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('a' - '*'), r, g, b);
        break;
     case 7:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('t' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('c' - '*'), r, g, b);
        break;
     case 8:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('a' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('c' - '*'), r, g, b);
        break;
     case 9:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('a' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('t' - '*'), r, g, b);
        break;
     case 10:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('g' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('c' - '*'), r, g, b);
        break;
     case 11:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('g' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('t' - '*'), r, g, b);
        break;
     case 12:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('g' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('c' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('a' - '*'), r, g, b);
        break;
     case 13:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('a' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('c' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('t' - '*'), r, g, b);
        break;
     case 14:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('g' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('a' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('t' - '*'), r, g, b);
        break;
     case 15:
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('g' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('t' - '*'), r, g, b);
        setUpLetterColor(adp->newstyle->colorRefs,(Uint1)('c' - '*'), r, g, b);
        break;
     default:
        break;
    }
  }
  return;
}

static void ResetBnStyleProc (ButtoN b)
{
  WindoW             w, wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  Uint2              j;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  w = dbdp->w;
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;

  adp->newstyle->font = adp->font;
  adp->newstyle->interline = adp->interline;
  adp->newstyle->marginwithindex = adp->marginwithindex;
  adp->newstyle->marginwithpos = adp->marginwithpos;
  adp->newstyle->marginleft = adp->marginleft;
  adp->newstyle->draw_scale = adp->draw_scale;
  adp->newstyle->draw_bars = adp->draw_bars;
  adp->newstyle->rowpcell = adp->rowpcell;
  adp->newstyle->columnpcell = adp->columnpcell;
  for (j=0; j<256; j++)
         adp->newstyle->colorRefs[j] = adp->colorRefs[j];
  adp->newstyle = MemFree (adp->newstyle);
  Remove (wdialog);
  return;
} 

static void ApplyBnStyleProc (ButtoN b)
{
  WindoW             w, wdialog;        
  WindoW             temport;
  DialogBoxDataPtr   dbdp;
  SeqEditViewFormPtr wdp;
  EditAlignDataPtr   adp;
  Uint2              j;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  w = dbdp->w;
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);

  WatchCursor();
  if (adp->font != adp->newstyle->font)
  {
         adp->font = adp->newstyle->font;
         SelectFont ((FonT)(adp->font));
         adp->lineheight = LineHeight () + adp->interline;
         adp->leading = Leading ();
         adp->charw   = CharWidth ('0');
         adp->ascent  = Ascent ();
         adp->scaleheight  = 2 * adp->lineheight;
         adp->margin.left  = adp->marginleft * adp->charw;
  }
  if (adp->interline != adp->newstyle->interline) {
         adp->interline = adp->newstyle->interline;
         SelectFont ((FonT)(adp->font));
         adp->lineheight = LineHeight () + adp->interline;
  }
  if (adp->marginwithindex != adp->newstyle->marginwithindex)
         adp->marginwithindex = adp->newstyle->marginwithindex;
  if (adp->marginwithpos != adp->newstyle->marginwithpos)
         adp->marginwithpos = adp->newstyle->marginwithpos;

  if (adp->marginleft != adp->newstyle->marginleft) {
         adp->marginleft = adp->newstyle->marginleft;
         adp->margin.left = adp->marginleft * adp->charw;
  }
  if (adp->draw_scale != adp->newstyle->draw_scale)
         adp->draw_scale = adp->newstyle->draw_scale;
  if (adp->draw_bars != adp->newstyle->draw_bars)
         adp->draw_bars = adp->newstyle->draw_bars;
  adp->nscaleline = 0;
  if ( adp->draw_scale ) adp->nscaleline ++;
  if ( adp->draw_bars )  adp->nscaleline ++;
  if (adp->rowpcell != adp->newstyle->rowpcell) {
     adp->rowpcell = adp->newstyle->rowpcell;
  }
  if (adp->columnpcell != adp->newstyle->columnpcell) 
  {
     Int2 x, y;
     adp->columnpcell = adp->newstyle->columnpcell;

     y = 0; x = 0;
     if (adp->columnpcell > 0) {
        y = (Int2) (adp->pnlWidth -adp->marginleft) / (Int2) adp->columnpcell;
        x = (Int2)(adp->pnlWidth -adp->marginleft -y) % (Int2)(adp->columnpcell);
     }
     adp->visibleWidth = (Int2) (adp->pnlWidth -adp->marginleft -y -x);
     if (adp->visibleWidth < 10) {
        adp->firstssp = NULL;
     }
  }
  for (j=0; j<256; j++)
         adp->colorRefs[j] = adp->newstyle->colorRefs[j];
  adp->newstyle = MemFree (adp->newstyle);
  do_resize_window (wdp->pnl, adp, TRUE);
  temport = SavePort(w);
  Select (wdp->pnl);
  inval_panel (wdp->pnl, -1 ,-1); 
  RestorePort (temport);
  Update ();
  Remove (wdialog);
  Update ();
  ArrowCursor();
  return;
}
/***************************************************
*** Translation
****************************************************/
static void SetRfFunc (PaneL pnl, Uint2 rf)
{
  WindoW           temport;
  EditAlignDataPtr adp;
  SeqEditViewFormPtr wdp;
  SelStructPtr     ssp;
  ValNodePtr       vnp;
  SeqParamPtr      prm;
  Int2             j;

  wdp = (SeqEditViewFormPtr) GetObjectExtra ((WindoW)ParentWindow (pnl));
  if ( ( adp = GetAlignDataPanel (pnl) ) == NULL ) return;
  if (adp->seqnumber == 0 || ISA_aa(adp->mol_type))
    return;
  temport = SavePort(pnl);
  Select (pnl);
  inval_panel (pnl, -1, -1);
  if ( checkOMss_for_itemtype (OBJ_BIOSEQ) == 0 ) {
         ssp = &(adp->master);
  }
  else ssp = ObjMgrGetSelected();  
  for (; ssp != NULL; ssp = ssp->next) {
     if (checkssp_for_editor (ssp) && ssp->itemtype == OBJ_BIOSEQ) {
        for (vnp = adp->params; vnp != NULL; vnp = vnp->next) {
           prm = (SeqParamPtr) vnp->data.ptrvalue;
           if ( prm->entityID == ssp->entityID ) {
             if (rf < (Uint2) 6) {
                prm->rf[rf] = TRUE; 
             }
             switch ((int)rf) {
               case 0: 
                       prm->rf[1] = prm->rf[2] = FALSE; 
                       SetStatus (wdp->rfitem [1], FALSE);
                       SetStatus (wdp->rfitem [2], FALSE);
                       SetStatus (wdp->rfitem [6], FALSE);
                       SetStatus (wdp->rfitem [8], FALSE);
                       break;
               case 1: 
                       prm->rf[0] = prm->rf[2] = FALSE; 
                       SetStatus (wdp->rfitem [0], FALSE);
                       SetStatus (wdp->rfitem [2], FALSE);
                       SetStatus (wdp->rfitem [6], FALSE);
                       SetStatus (wdp->rfitem [8], FALSE);
                       break;
               case 2: 
                       prm->rf[0] = prm->rf[1] = FALSE; 
                       SetStatus (wdp->rfitem [0], FALSE);
                       SetStatus (wdp->rfitem [1], FALSE);
                       SetStatus (wdp->rfitem [6], FALSE);
                       SetStatus (wdp->rfitem [8], FALSE);
                       break;
               case 3: 
                       prm->rf[4] = prm->rf[5] = FALSE; 
                       SetStatus (wdp->rfitem [4], FALSE);
                       SetStatus (wdp->rfitem [5], FALSE);
                       SetStatus (wdp->rfitem [7], FALSE);
                       SetStatus (wdp->rfitem [8], FALSE);
                       break;
               case 4: 
                       prm->rf[3] = prm->rf[5] = FALSE; 
                       SetStatus (wdp->rfitem [3], FALSE);
                       SetStatus (wdp->rfitem [5], FALSE);
                       SetStatus (wdp->rfitem [7], FALSE);
                       SetStatus (wdp->rfitem [8], FALSE);
                       break;
               case 5: 
                       prm->rf[3] = prm->rf[4] = FALSE; 
                       SetStatus (wdp->rfitem [3], FALSE);
                       SetStatus (wdp->rfitem [4], FALSE);
                       SetStatus (wdp->rfitem [7], FALSE);
                       SetStatus (wdp->rfitem [8], FALSE);
                       break;
               case 6: if (GetStatus (wdp->rfitem [6]) == FALSE) {
                        prm->rf[0] =FALSE; prm->rf[1] =FALSE; prm->rf[2]=FALSE; 
                       } else {
                        prm->rf[0] =TRUE; prm->rf[1] =TRUE; prm->rf[2] =TRUE; 
                       }
                       SetStatus (wdp->rfitem [0], FALSE);
                       SetStatus (wdp->rfitem [1], FALSE);
                       SetStatus (wdp->rfitem [2], FALSE);
                       SetStatus (wdp->rfitem [8], FALSE);
                       break;
               case 7: if (GetStatus (wdp->rfitem [7]) == FALSE) {
                        prm->rf[3] =FALSE; prm->rf[4] =FALSE; prm->rf[5]=FALSE; 
                       } else {
                        prm->rf[3] = TRUE; prm->rf[4] = TRUE; prm->rf[5] =TRUE; 
                       }
                       SetStatus (wdp->rfitem [3], FALSE);
                       SetStatus (wdp->rfitem [4], FALSE);
                       SetStatus (wdp->rfitem [5], FALSE);
                       SetStatus (wdp->rfitem [8], FALSE);
                       break;
               case 8: 
                       if (GetStatus (wdp->rfitem [8]) == FALSE) {
                        prm->rf[0] =FALSE; prm->rf[1] =FALSE; prm->rf[2]=FALSE; 
                        prm->rf[3] =FALSE; prm->rf[4] =FALSE; prm->rf[5]=FALSE; 
                       } else {
                        prm->rf[0] =TRUE; prm->rf[1] =TRUE; prm->rf[2] =TRUE; 
                        prm->rf[3] =TRUE; prm->rf[4] =TRUE; prm->rf[5] =TRUE; 
                       }
                       for (j=0; j<8; j++) 
                          SetStatus (wdp->rfitem [0], FALSE);
                       break;
               case 9: 
                       prm->rf[0] =FALSE; prm->rf[1] =FALSE; prm->rf[2]=FALSE; 
                       prm->rf[3] =FALSE; prm->rf[4] =FALSE; prm->rf[5]=FALSE; 
                       for (j=0; j<9; j++) 
                          SetStatus (wdp->rfitem [j], FALSE);
                       break;
               default : break;
             } 
           } 
        } 
     }
  }
  data_collect_arrange (adp, TRUE);
  SeqEdSetCorrectBarMax (pnl, adp->nlines, adp->voffset);
  RestorePort (temport);
}

extern void rf1ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)0);
}
extern void rf2ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)1);
}
extern void rf3ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)2);
}
extern void rf4ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)3);
}
extern void rf5ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)4);
}
extern void rf6ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)5);
}
extern void rf7ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)6);
}
extern void rf8ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)7);
}
extern void rf9ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)8);
}
extern void rf10ItemProc (IteM i) {
   SetRfFunc ((PaneL)GetPanelFromWindow((WindoW)ParentWindow (i)), (Uint2)9);
}

/***************************************************
***   Dialog Proc
***     read number of sequences and length 
***       to allocate the rigth size       
***       SALSA_PHYLIP: read in the file
***       others: read from a dialog box            
****************************************************/
extern void SelectSequenceFormat (PopuP p)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;

  wdialog = ParentWindow (p);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  switch ( (Int2) GetValue (p) ) {
      case 1:
         dbdp->align_format = SALSA_FASTA;
         break;
      case 2:
         dbdp->align_format = SALSA_IDS;
         break;
      case 3:
         dbdp->align_format = SALSA_GBFF;
         break;
      default:
         break;
  }
  return;
}

extern void SelectAlignFormat (PopuP p)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;

  wdialog = ParentWindow (p);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  dbdp->align_format = (Int2) GetValue (p);
  return;
}

extern void SelectMolType (PopuP p)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;

  wdialog = ParentWindow (p);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  switch ( (Int2) GetValue (p) ) {
      case 1:
         dbdp->mol_type = Seq_mol_dna;
         break;
      case 2:
         dbdp->mol_type = Seq_mol_rna;
         break;
      case 3:
         dbdp->mol_type = Seq_mol_aa;
         break;
      case 4:
         dbdp->mol_type = Seq_mol_na;
         break;
      default:
         break;
  }
  return;
}

extern void FileInProc (ButtoN b)
{
  WindoW           wdialog;
  DialogBoxDataPtr dbdp;
  FILE             *fp;
  int              n;
  long             lens;
  Char             namep [PATH_MAX];
  Char             str [64];
 
  wdialog = ParentWindow (b);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if (!GetInputFileName (namep, PATH_MAX, "", "TEXT"))
    return;
  SetTitle (dbdp->txt1, namep);  
  if ( dbdp->align_format == SALSA_PHYLIP ) {
         if ( (fp = FileOpen (namep, "r")) == NULL) {
                return ;
         }
         fscanf(fp, "%d %d", &n, &lens);
         FileClose(fp);
         if (n < 2) {
                return ;
         }
         sprintf (str, "%d", (int) n);
         SetTitle (dbdp->txt2, str);
         sprintf (str, "%ld", (long) lens);
         SetTitle (dbdp->txt3, str);
  }
  return;
}

extern void DefinePanelDialog (IteM i)
{
  WindoW           w, wdialog;
  EditAlignDataPtr adp;
  DialogBoxDataPtr dbdp;
  GrouP            g;
  ButtoN           esBtn1, esBtn3, esBtn4;
  PopuP            popf;
  Char             str[16];
  Uint2            j;

  w = ParentWindow (i);
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;

  adp->newstyle = (AlignStylePtr) MemNew (sizeof (AlignStyle));
  adp->newstyle->font = adp->font;
  adp->newstyle->interline = adp->interline;
  adp->newstyle->marginwithindex = adp->marginwithindex;
  adp->newstyle->marginwithpos = adp->marginwithpos;
  adp->newstyle->marginleft = adp->marginleft;
  adp->newstyle->draw_scale = adp->draw_scale;
  adp->newstyle->draw_bars = adp->draw_bars;
  adp->newstyle->rowpcell = adp->rowpcell;
  adp->newstyle->columnpcell = adp->columnpcell;
  for (j=0; j< 256; j++) 
         adp->newstyle->colorRefs[j] = adp->colorRefs[j];

  wdialog = FixedWindow (-50, -33, -10, -10, "Sequence", NULL);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  g = HiddenGroup (wdialog, 1, 0, NULL);
  PushButton (g, "Font", ChangeFontButton );

  g = HiddenGroup (wdialog, 2, 0, NULL);
  esBtn4 = CheckBox (g, "draw scale", DrawScaleButton);  
  SetStatus (esBtn4, adp->draw_scale);
  esBtn4 = CheckBox (g, "draw bars", DrawBarScaleButton);  
  SetStatus (esBtn4, adp->draw_bars);

  g = HiddenGroup (wdialog, 4, 0, NULL);
  StaticPrompt (g, "column per cell",0,popupMenuHeight,systemFont, 'l');
  sprintf (str, "%d", (int) adp->columnpcell);
  dbdp->txt4 = DialogText (g, str, 5, ColumnCellButton);
  StaticPrompt (g, "Interline", 0, popupMenuHeight, systemFont, 'l');
  sprintf (str, "%d", (int) adp->interline);
  dbdp->txt3 = DialogText (g, str, 5, ChangeInterlineButton);

  g = HiddenGroup (wdialog, 4, 0, NULL);
  StaticPrompt (g, "Margin : ", 0, popupMenuHeight, systemFont, 'l');
  sprintf (str, "%d", (int) adp->marginleft);
  dbdp->txt5 = DialogText (g, str, 5, ChangeMarginButton);
  StaticPrompt (g, "columns", 0, popupMenuHeight,systemFont, 'l');
  popf = PopupList (g, TRUE, ChangeMarginPopup);
  PopupItem (popf, " ");
  PopupItem (popf, "PHYLIP");
  PopupItem (popf, "ClustalV");
  PopupItem (popf, "Macaw");
  SetValue  (popf, 1);

  g = HiddenGroup (wdialog, 2, 0, NULL);
  esBtn3 = CheckBox (g, "with positions", ChangeMarginPosButton);  
  SetStatus (esBtn3, adp->marginwithpos);
  esBtn1 = CheckBox (g, "with index  ", ChangeMarginIndexButton);  
  SetStatus (esBtn1, adp->marginwithindex);

if (ISA_na(adp->mol_type))
{
  g = HiddenGroup (wdialog, 7, 0, NULL);
  StaticPrompt (g, "Color:", 0, popupMenuHeight, systemFont, 'l');
  PushButton (g, "SeqId", SelectColorIdButton );
  PushButton (g, "master", SelectColorIdMasterButton );
  PushButton (g, "selection", SelectColorSelectionButton );
  PushButton (g, "scale", SelectColorScaleButton );
  popf = PopupList (g, TRUE, SelectAnColorPopup);
  PopupItem (popf, "ATGC");
  PopupItem (popf, "A");
  PopupItem (popf, "C");
  PopupItem (popf, "G");
  PopupItem (popf, "T");
  PopupItem (popf, "AG");
  PopupItem (popf, "TC");
  PopupItem (popf, "AC");
  PopupItem (popf, "AT");
  PopupItem (popf, "GC");
  PopupItem (popf, "GT");
  PopupItem (popf, "GCA");
  PopupItem (popf, "ACT");
  PopupItem (popf, "GAT");
  PopupItem (popf, "GTC");
  PopupItem (popf, "-");
  SetValue  (popf, 1);
  popf = PopupList (g, TRUE, SelectAaColorPopup);
  PopupItem (popf, "All");
  PopupItem (popf, "MFLVI");
  PopupItem (popf, "RHEK");
  PopupItem (popf, "TNPSDQ");
  PopupItem (popf, "GCA");
  PopupItem (popf, "*");
  SetValue  (popf, 1);
}
  g = HiddenGroup (wdialog, 2, 0, NULL);
  PushButton (g, "Cancel", ResetBnStyleProc);
  PushButton (g, "Accept", ApplyBnStyleProc);
  Show (wdialog); 
  return;
}

extern void inval_panel (PaneL pnl, Int2 x, Int2 y)
{
  RecT     r;

  get_client_rectxy (pnl, &r, x, y);
  InvalRect (&r);
  return;
}

extern void inval_window (WindoW w)
{
  PaneL pnl;
  pnl = GetPanelFromWindow (w);
  if (pnl != NULL) { 
     inval_panel (pnl, -1 ,-1);
  }
  return ;
}

extern void inval_rect (Int2 left, Int2 top, Int2 right, Int2 bottom)
{
  RecT         rm;
  LoadRect (&rm, left, top, right, bottom);
  InsetRect (&rm, -1, -1);
  InvalRect (&rm); 
}


extern void inval_allines (EditAlignDataPtr adp, RecT *rp)
{
  Int2         ptx, pty, Seqlens_pix, j;

  ptx = rp->left + adp->margin.left;
  pty = rp->top;
  j = adp->visibleWidth;
  if (adp->columnpcell > 0) 
     j += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;
  Seqlens_pix = j * adp->charw;
  for (j = 0; j < adp->pnlLine; j++) 
  {
     if ( adp->seqEntity_id[j] > 0 ) {
        inval_rect (ptx, pty, (Int2)(ptx + Seqlens_pix), (Int2)(pty + adp->lineheight));
     }
     pty += adp->lineheight;
  }
}

extern void inval_selstruct (EditAlignDataPtr adp, Uint2 ei, Uint2 ii, Uint2 it, Uint2 itemsubtype, RecT *rp, Int2 left, Int2 toright)
{
  Int2         ptx = rp->left + left, 
               pty = rp->top, 
               j;
  for (j = 0; j < adp->pnlLine; j++) 
  {
    if (is_sameId(adp->seqEntity_id[j], adp->item_id[j], adp->itemtype[j], 255, ei, ii, it, 255) )
    {
       if (itemsubtype == 255 || adp->itemsubtype[j] == itemsubtype) {
          inval_rect (ptx, pty, (Int2)(ptx + toright), (Int2)(pty + adp->lineheight));
       }
    }
    pty += adp->lineheight;
  }
}

extern void inval_selstructpos (EditAlignDataPtr adp, Uint2 ei, Uint2 ii, Uint2 it, RecT *rp, Int4 pos)
{
  Int2         ptx = rp->left, 
               pty;
  Int2         column1, line1;

  SeqPosToLineColumn (ii, ei, it, pos, &line1, &column1, adp->hoffset, adp);
  pty = rp->top + line1 * adp->lineheight;
  inval_rect ((Int2)(ptx + adp->margin.left + (column1-2) * adp->charw),  
              (Int2)(pty), (Int2)(ptx+adp->margin.left+(column1+2) *adp->charw),
              (Int2)(pty + adp->lineheight));
  return;
}

extern void inval_selstructpos_tobottom (EditAlignDataPtr adp, Uint2 ei, Uint2 ii, Uint2 it, RecT *rp, Int4 pos)
{
  Int2         ptx = rp->left, 
               pty;
  Int2         column1, line1;
  Int2         width;

  SeqPosToLineColumn (ii, ei, it, pos, &line1, &column1, adp->hoffset, adp);
  pty = rp->top + line1 * adp->lineheight;
  width = adp->visibleWidth;
  if (adp->columnpcell > 0)
       width += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;

  inval_rect ((Int2)(ptx + adp->margin.left + (column1-2) * adp->charw),  
              (Int2)(pty), (Int2)(ptx+adp->margin.left+ width *adp->charw),
              (Int2)(pty + adp->lineheight));
  pty += adp->lineheight;
  inval_rect ((Int2)(ptx),
              (Int2)(pty), (Int2)(ptx+adp->margin.left+ width *adp->charw),
              (Int2) rp->bottom);
  return;
}

extern void inval_selstructloc (EditAlignDataPtr adp, Uint2 ei, Uint2 ii, Uint2 it, Uint2 itemsubtype, RecT *rp, Int4 from, Int4 to)
{
  RecT         rid;
  Int2         ptx = rp->left, 
               pty;
  Int2         column1, column2, line1;
  Int2         width, j;
  Int2         aln,
               alnmax;
  Boolean      bool1, bool2;

  aln = -1;
  alnmax=-1;
  for (j=0; j<MAXLineWindow; j++) {
        if (ii == adp->item_id[j] && ei == adp->seqEntity_id[j] && it == adp->itemtype[j] ) {
           if (aln<0)
              aln++;  
           alnmax++;
        }
  }
  if (aln<0)
     return;
  if (from < adp->colonne[adp->hoffset - adp->bufferstart]) {     
     from = adp->colonne[adp->hoffset - adp->bufferstart + aln*adp->visibleWidth];
  }
  from += (aln * adp->visibleWidth);
  SeqPosToLineColumn (ii, ei, it, from, &line1, &column1, adp->hoffset, adp);
  if (line1<0)
     from = adp->colonne[adp->hoffset - adp->bufferstart + aln*adp->visibleWidth];

  while (line1 < 0 && aln < alnmax)
  {
     aln++;
     from += (Int4)(adp->visibleWidth);
     SeqPosToLineColumn (ii, ei, it, from, &line1, &column1, adp->hoffset, adp);
  }
  if (line1 < 0)
     return;
  width = adp->visibleWidth;
  if (adp->columnpcell > 0) 
       width += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;
  pty = rp->top + line1 * adp->lineheight;

  for (j = line1; j < adp->pnlLine; j++) 
  {
    if (is_sameId(adp->seqEntity_id[j], adp->item_id[j], adp->itemtype[j], 255, ei, ii, it, 255)) 
    {
     if (itemsubtype == 255 || adp->itemsubtype[j] == itemsubtype) 
     {
       bool1 = SeqPosInLineColumn (from, adp->alignline[j], &column1, adp->hoffset, adp);
       bool2 = SeqPosInLineColumn (to, adp->alignline[j], &column2, adp->hoffset, adp);
       if (bool1 && bool2) {
          LoadRect (&rid, (Int2)(ptx +adp->margin.left+(column1-1) *adp->charw),
                  pty, (Int2)(ptx + adp->margin.left + (column2+1) *adp->charw),
                       (Int2)(pty + adp->lineheight));
          InvalRect (&rid);
       }
       else if ( bool1 ) {
          LoadRect (&rid, (Int2)(ptx+adp->margin.left+(column1-1) *adp->charw), 
                  pty, (Int2)(ptx + adp->margin.left + width * adp->charw), 
                       (Int2)(pty + adp->lineheight));
          InvalRect (&rid);
       }
       else if ( bool2 ) {
          LoadRect (&rid, (Int2)(ptx + adp->margin.left), pty, 
                     (Int2)(ptx + adp->margin.left + (column2+1) * adp->charw), 
                     (Int2)(pty + adp->lineheight));
          InvalRect (&rid);
       }
       else if ( !bool1 && !bool2 && column1 == -1 && column2 == -2 ) { 
          LoadRect (&rid, (Int2)(ptx + adp->margin.left), pty, 
                          (Int2)(ptx + adp->margin.left + width * adp->charw), 
                          (Int2)(pty + adp->lineheight));
          InvalRect (&rid);
       }
     }
    }
    pty += adp->lineheight;
  }
}


extern void inval_selstructloc_forfeat (EditAlignDataPtr adp, Uint2 ei, Uint2 ii, Uint2 it, Uint2 itemsubtype, RecT *rp, Int4 from, Int4 to)
{
  RecT         rid;
  Int2         ptx = rp->left, 
               pty;
  Int2         column1, column2, line1;
  Int2         width, j;
  Uint2        aln;
  Boolean      bool1, bool2;

  if (from < adp->colonne[adp->hoffset - adp->bufferstart]) {     
     aln = 0;
  for (j=0; j<MAXLineWindow; j++) {
        if (ii == adp->item_id[j] && ei == adp->seqEntity_id[j] && it == adp->itemtype[j] ) {
           aln = adp->alignline [j]; 
           break;
        }
     }
     from = adp->colonne[adp->hoffset - adp->bufferstart + aln*adp->visibleWidth];
  }
  SeqPosToLineColumn (ii, ei, it, from, &line1, &column1, adp->hoffset, adp);
  if (line1 < 0) {
     return;
  }
  width = adp->visibleWidth;
  if (adp->columnpcell > 0) 
       width += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;
  pty = rp->top + line1 * adp->lineheight;
  
  for (j = line1; j < adp->pnlLine; j++) 
  {

    if (is_sameId(adp->seqEntity_id[j], adp->item_id[j], adp->itemtype[j], 255, ei, ii, it, 255)) 
    {
     if (itemsubtype == 255 || adp->itemsubtype[j] == itemsubtype) 
     {
       bool1 = SeqPosInLineColumn (from, adp->alignline[j], &column1, adp->hoffset, adp);
       bool2 = SeqPosInLineColumn (to, adp->alignline[j], &column2, adp->hoffset, adp);
       if (bool1 && bool2) {
          LoadRect (&rid, (Int2)(ptx +adp->margin.left+(column1-1) *adp->charw),
                  pty, (Int2)(ptx + adp->margin.left + (column2+1) *adp->charw),
                       (Int2)(pty + adp->lineheight));
          InvalRect (&rid);
       }
       else if ( bool1 ) {
          LoadRect (&rid, (Int2)(ptx+adp->margin.left+(column1-1) *adp->charw), 
                  pty, (Int2)(ptx + adp->margin.left + width * adp->charw), 
                       (Int2)(pty + adp->lineheight));
          InvalRect (&rid);
       }
       else if ( bool2 ) {
          LoadRect (&rid, (Int2)(ptx + adp->margin.left), pty, 
                     (Int2)(ptx + adp->margin.left + (column2+1) * adp->charw), 
                     (Int2)(pty + adp->lineheight));
          InvalRect (&rid);
       }
       else if ( !bool1 && !bool2 && column1 == -1 && column2 == -2 ) { 
          LoadRect (&rid, (Int2)(ptx + adp->margin.left), pty, 
                          (Int2)(ptx + adp->margin.left + width * adp->charw), 
                          (Int2)(pty + adp->lineheight));
          InvalRect (&rid);
       }
     }
    }
    pty += adp->lineheight;
  }
}
extern void inval_rectidselected (EditAlignDataPtr adp, RecT *rp)
{
  SelStructPtr ssp = NULL;
  SelStruct    sp;
  Int2         ptx, pty;
  Int2         width, j;

  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) return;
  ptx = rp->left;
  pty = rp->top;
  for (j = 0; j < adp->pnlLine; j++) 
  {
        sp.entityID = adp->seqEntity_id[j];
        sp.itemID = adp->item_id[j];
        sp.itemtype = adp->itemtype[j];
        sp.regiontype = 0;
        if( (ssp = is_selectedbyID (sp.entityID, sp.itemID, sp.itemtype))!=NULL)
        {
           width = adp->visibleWidth;
           if (adp->columnpcell > 0) 
                 width += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;
           inval_rect (ptx, pty, (Int2)(ptx +adp->margin.left +width *adp->charw), 
                            (Int2)(pty + adp->lineheight));
        }
        pty += adp->lineheight;
  }
}

extern void inval_pep (EditAlignDataPtr adp, RecT *rp)
{
  Int2         ptx, pty;
  Int2         width, j;

  ptx = rp->left;
  pty = rp->top;
  for (j = 0; j < adp->pnlLine; j++) {
        if((adp->itemsubtype[j] >=EDITDEF_RF1 
         && adp->itemsubtype[j] <=EDITDEF_RF6)
        || adp->itemsubtype[j] == FEATDEF_TRSL ) {
           width = adp->visibleWidth;
           if (adp->columnpcell > 0) 
                 width +=(Int2) adp->visibleWidth / (Int2) adp->columnpcell;
           inval_rect ((Int2)(ptx +adp->margin.left), pty,
                       (Int2)(ptx +adp->margin.left + width *adp->charw), 
                       (Int2)(pty + adp->lineheight));
        }
        pty += adp->lineheight;
  }
}

/******************************************************************/
extern void inval_all (EditAlignDataPtr adp, RecT *rp, Uint2 itemtype1, Uint2 itemtype2, Uint2 itemtype3, Int2 width)
{
  SelStructPtr ssp = NULL;
  ssp = ObjMgrGetSelected ();
  if (ssp == NULL) return;
  while (ssp != NULL) 
  {
         if ( checkssp_for_editor(ssp) && (ssp->itemtype == itemtype1 
         || ssp->itemtype == itemtype2 || ssp->itemtype == itemtype3)) {
            if (ssp->itemtype == OBJ_SEQFEAT || ssp->itemtype == OBJ_VIRT) 
               inval_selstruct (adp, ssp->entityID, ssp->itemID, ssp->itemtype, 255, rp, adp->margin.left, (Int2)(width * adp->charw));
            else
               inval_selstructloc (adp, ssp->entityID, ssp->itemID, ssp->itemtype, 255, rp, SeqLocStart ((SeqLocPtr) ssp->region), SeqLocStop((SeqLocPtr) ssp->region)); 
         }
         ssp = ssp->next;
  }
}

/******************************************************************/
/******************************************************************/
static void refreshlst (LisT lst, EditAlignDataPtr adp)
{
  AlignNodePtr     anp;
  SeqParamPtr      prm;
  ValNodePtr       vnp, vnp2;
  Char             str [68];
  Char             str2 [24];
  CharPtr          tmp;

  Reset (lst);
  vnp2 = adp->params;
  for (vnp = adp->anp_list; vnp!=NULL && vnp2!=NULL; vnp = vnp->next)
  {
     str [0] = '\0';
     tmp = str;
     prm = (SeqParamPtr) vnp2->data.ptrvalue;
     sprintf (str2, "%d ", prm->group);
     tmp = StringMove (tmp, str2);
     anp = (AlignNodePtr) vnp->data.ptrvalue;
     SeqIdWrite (anp->sip, str2, PRINTID_REPORT, sizeof (str2));
     tmp = StringMove (tmp, str2);
     ListItem (lst, str);
     vnp2 = vnp2->next;
  }
}

static void grouplst (ButtoN b)
{
  DialogBoxDataPtr dbdp;
  LisT             lst;
  EditAlignDataPtr adp;
  ValNodePtr       vnp;
  SeqParamPtr      prm;
  Int2             n;
  Int2             j, jmax;
  Boolean          val,
                   select = FALSE;

  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (b));
  if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) {
     lst = dbdp->lst1;
     n = adp->ngroup + 1;
     vnp = adp->params;
     jmax = CountItems (lst);
     for (j=1; j <= jmax && vnp != NULL; j++, vnp = vnp->next)
     {   
        val = GetItemStatus (lst, j);
        if (val) {
           prm = (SeqParamPtr) vnp->data.ptrvalue;
           prm->group = n;
           select = TRUE;
        }
     }
     if (select) {
        adp->ngroup ++;
        refreshlst (lst, adp);
     }
  }
}

static void ungrouplst (ButtoN b)
{
  DialogBoxDataPtr dbdp;
  LisT             lst;
  EditAlignDataPtr adp;
  ValNodePtr       vnp;
  SeqParamPtr      prm;
  Int2             j, jmax;
  Boolean          val;

  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (b));
  if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) {
     lst = dbdp->lst1;
     vnp = adp->params;
     jmax = CountItems (lst);
     for (j=1; j <= jmax && vnp != NULL; j++, vnp = vnp->next)
     {   
        val = GetItemStatus (lst, j);
        if (val) {
           prm = (SeqParamPtr) vnp->data.ptrvalue;
           prm->group = 0;
        }
     }
     adp->ngroup = 0;
     refreshlst (lst, adp);
  }
}

static void mergelst (ButtoN b)
{
  DialogBoxDataPtr dbdp;
  LisT             lst;
  EditAlignDataPtr adp;
  ValNodePtr       vnp=NULL;
  SeqParamPtr      prm;
  Int2             j, jmax;
  Boolean          val;

  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (b));
  if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) {
     lst = dbdp->lst1;
     jmax = CountItems (lst);
     for (j=1; j <= jmax && vnp != NULL; j++, vnp = vnp->next)
     {   
        val = GetItemStatus (lst, j);
        if (val) {
/** ???????? **/
        }
     }
  }
}

static void unmergelst (ButtoN b)
{
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  ValNodePtr       vnp;
  SeqParamPtr      prm;

  dbdp = (DialogBoxDataPtr) GetObjectExtra (ParentWindow (b));
  if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) {
/** ??????? **/
  }
}
static void DismissGroupButton (ButtoN b)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  ValNodePtr         vnp;
  SeqParamPtr        prm;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  {
     for (vnp = adp->params; vnp != NULL; vnp = vnp->next)
     {   
        prm = (SeqParamPtr) vnp->data.ptrvalue;
        prm->group = 0;
     }
     adp->ngroup = 0;
  }  
  Remove (wdialog);
  return;
}

static void AcceptGroupButton (ButtoN b)
{
  WindoW             w, wdialog;        
  WindoW             temport;
  DialogBoxDataPtr   dbdp;
  SeqEditViewFormPtr wdp;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  w = dbdp->w;
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  temport = SavePort(w);
  Select (wdp->pnl);
  inval_panel (wdp->pnl, -1 ,-1);
  RestorePort (temport);
  Update ();
  Remove (wdialog);
  Update ();
  return;
}
/******************************************************************/
extern void FormulaDialog (IteM i)
{
  WindoW           w, wdialog;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  GrouP            g1, g2, g3;
  ButtoN           b1, b2, b3, b4, b5;
  PrompT           p;
  Int2             numberlines = 15;

  w = ParentWindow (i);
  wdialog = FixedWindow (-50, -33, -10, -10, " ", StdCloseWindowProc);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;
 
  g1 = HiddenGroup (wdialog, 2, 0, NULL);

  g2 = HiddenGroup (g1, 0, 2, NULL);
  p = StaticPrompt (g2, "Select sequences", 0, dialogTextHeight, systemFont, 'c');         
  dbdp->lst1 = MultiList (g2, numberlines, 10, NULL);
  refreshlst (dbdp->lst1, adp);
  AlignObjects (ALIGN_CENTER, (HANDLE) dbdp->lst1, (HANDLE) p, NULL);

  g3 = HiddenGroup (g1, 0, 6, NULL);
  b1 = PushButton (g3, "Group", grouplst);
  b2 = PushButton (g3, "Ungroup", ungrouplst);
  b3 = PushButton (g3, "Merge", mergelst);
  b3 = PushButton (g3, "All visible", unmergelst);
  b4 = PushButton (g3, "Accept", AcceptGroupButton);
  b5 = PushButton (g3, "Dismiss", DismissGroupButton);
  AlignObjects (ALIGN_CENTER, (HANDLE) b1, (HANDLE) b2, (HANDLE) b3, (HANDLE) b4, NULL);

  RealizeWindow (wdialog);
  Show (wdialog);
  return;
}
/******************************************************************/
/******************************************************************/
static void BoolButton (GrouP g)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;

  wdialog = ParentWindow (g);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  dbdp->bool = (Boolean)((GetValue (g)) == 1);
  return;
}

/******************************************************************
static void SelectMatrixPopup (PopuP p)
{
  WindoW           wdialog;        
  DialogBoxDataPtr dbdp;

  wdialog = ParentWindow (p);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  dbdp->matrix = GetValue (p);
  return;
}

static void AlignToFasta2Button (ButtoN b)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  SeqEntryPtr        sep;
  ValNodePtr         sqloc_out;
  SeqAlignPtr        salp;
  Int2               nseq;
  Boolean            new;

  WatchCursor ();
  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  sep = FastaRead (NULL, adp->mol_type);
  if (sep != NULL) {
     sqloc_out = SeqEntryToSeqLoc (sep, &nseq, adp->mol_type);
     if (sqloc_out != NULL) {
        salp = AlignToFunc(adp, sqloc_out, dbdp->bool, 2);
        LaunchAlignEditor(salp);
     }
  }
  Remove (wdialog);
  Update ();
  ArrowCursor ();
  return;
}

static void AlignTogi2Button (ButtoN b)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  ValNodePtr         sqloc_out;
  SeqAlignPtr        salp;
  Boolean            new;

  WatchCursor ();
  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  sqloc_out = IdRead (NULL);
  if (sqloc_out != NULL) {
     salp = AlignToFunc (adp, sqloc_out, dbdp->bool, 2);
     LaunchAlignEditor(salp);
  }
  Remove (wdialog);
  Update ();
  ArrowCursor ();
  return;
}

static void AlignToAsn2Button (ButtoN b)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  SeqEntryPtr        sep;
  ValNodePtr         sqloc_out;
  SeqAlignPtr        salp;
  Int2               nseq;
  Boolean            new;

  WatchCursor ();
  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  sep = seqentry_read (NULL);
  if (sep != NULL) { 
     sqloc_out = SeqEntryToSeqLoc (sep, &nseq, adp->mol_type);
     if (sqloc_out != NULL) {
        salp = AlignToFunc(adp, sqloc_out, dbdp->bool, 2);
        LaunchAlignEditor(salp);
     }
  }
  Remove (wdialog);
  Update ();
  ArrowCursor ();
  return;
}

static void AlignToFasta3Button (ButtoN b)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  SeqEntryPtr        sep;
  ValNodePtr         sqloc_out;
  SeqAlignPtr        salp;
  Int2               nseq;
  Boolean            new;

  WatchCursor ();
  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  sep = FastaRead (NULL, adp->mol_type);
  if (sep != NULL) {
     sqloc_out = SeqEntryToSeqLoc (sep, &nseq, adp->mol_type);
     if (sqloc_out != NULL) {
        salp = AlignToFunc(adp, sqloc_out, dbdp->bool, 3);
        LaunchAlignEditor(salp);
     }
  }
  Remove (wdialog);
  Update ();
  ArrowCursor ();
  return;
}

static void AlignTogi3Button (ButtoN b)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  ValNodePtr         sqloc_out;
  SeqAlignPtr        salp;
  Boolean            new;

  WatchCursor ();
  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  sqloc_out = IdRead (NULL);
  if (sqloc_out != NULL) {
     salp = AlignToFunc (adp, sqloc_out, dbdp->bool, 3);
     LaunchAlignEditor(salp);
  }
  Remove (wdialog);
  Update ();
  ArrowCursor ();
  return;
}

static void AlignToAsn3Button (ButtoN b)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  SeqEntryPtr        sep;
  ValNodePtr         sqloc_out;
  SeqAlignPtr        salp;
  Int2               nseq;
  Boolean            new;

  WatchCursor ();
  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  sep = seqentry_read (NULL);
  if (sep != NULL) { 
     sqloc_out = SeqEntryToSeqLoc (sep, &nseq, adp->mol_type);
     if (sqloc_out != NULL) {
        salp = AlignToFunc(adp, sqloc_out, dbdp->bool, 3);
        LaunchAlignEditor(salp);
     }
  }
  Remove (wdialog);
  Update ();
  ArrowCursor ();
  return;
}

static void BandAlignButton (ButtoN b)
{
  WindoW             wdialog;        
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  ValNodePtr         sqloc_out=NULL;
  SeqAlignPtr        salp;
  Boolean            new;

  WatchCursor ();
  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) != NULL ) { 
     salp = AlignToFunc (adp, sqloc_out, dbdp->bool, 5);
     if (new) 
        LaunchAlignEditor(salp);
     else
        repopulate_panel (dbdp->w, adp, salp);
  }
  Remove (wdialog);
  Update ();
  ArrowCursor ();
  return;
}

extern void AlignDialog (IteM i)
{
  WindoW           w, wdialog;
  SeqEditViewFormPtr wdp;
  DialogBoxDataPtr dbdp;
  GrouP            g1, g2, g3, g4, g5;
  GrouP            sg1, sg2;
  PrompT           p;
  PopuP            popf;
#ifdef SALSA_DEBUG
  ButtoN           b;
#endif

  w = ParentWindow (i);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  wdialog = FixedWindow (-50, -33, -10, -10, "Alignment", StdCloseWindowProc);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  dbdp->bool = TRUE;
  dbdp->separator = 2;
  dbdp->matrix = 4;

  g3 = HiddenGroup (wdialog, 1, 0, NULL);
  StaticPrompt (g3, "Recompute alignment of selected sequences", 0, dialogTextHeight, systemFont, 'l');

  g5 = HiddenGroup (wdialog, 2, 0, NULL);
  StaticPrompt (g5, "Similar nucleotide seq.", 0, dialogTextHeight, systemFont, 'l');  
  PushButton (g5, "Sim3", Sim3Button);
  StaticPrompt(g5,"Divergent nucleotide seq.", 0, dialogTextHeight, systemFont, 'l');
  PushButton (g5, "Sim2", Sim2Button);
  StaticPrompt (g5, "Stringent method ", 0, dialogTextHeight, systemFont, 'l');  
  PushButton (g5, "Sim", SimButton);
  StaticPrompt (g5, "Align translation ", 0, dialogTextHeight, systemFont, 'l'); 
  PushButton (g5, "Sim", SimThruProtButton);
#ifdef SALSA_DEBUG
  if (wdp->extended_align_menu) {
     b = PushButton (g5, "ClustalW", ClustalwButton);
     b = PushButton (g5, "ClustalW translation", ClustalwProtButton);
  }
#endif

  g3 = HiddenGroup (wdialog, 1, 0, NULL);
  StaticPrompt (g3, "Import and recompute alignment", 0, dialogTextHeight, systemFont, 'l');

  g5 = HiddenGroup (wdialog, 3, 0, NULL);
  StaticPrompt (g5, "Import FASTA", 0, dialogTextHeight, systemFont, 'l'); 
  PushButton (g5, "Sim2", AlignToFasta2Button);
  PushButton (g5, "Sim3", AlignToFasta3Button);
  StaticPrompt (g5, "Import NCBI SeqID", 0, dialogTextHeight, systemFont, 'l'); 
  PushButton (g5, "Sim2", AlignTogi2Button);
  PushButton (g5, "Sim3", AlignTogi3Button);
  StaticPrompt (g5, "Import ASN.1", 0, dialogTextHeight, systemFont, 'l'); 
  PushButton (g5, "Sim2", AlignToAsn2Button);
  PushButton (g5, "Sim3", AlignToAsn3Button);
#ifdef SALSA_DEBUG
  PushButton (g5, "BandAlign", BandAlignButton);
#endif
  g4 = HiddenGroup (wdialog, 2, 0, NULL);
  sg1 = HiddenGroup (g4, 1, 0, NULL);
  p = StaticPrompt (sg1, "Add extremities", 0, dialogTextHeight, systemFont, 'c');  
  sg2 = HiddenGroup (g4, 2, 0, BoolButton);
  RadioButton (sg2, "yes");
  RadioButton (sg2, "no");
  SetValue (sg2, 1);

  g1 = HiddenGroup (wdialog, 1, 0, NULL);
  StaticPrompt (g1, "Protein alignment", 0, dialogTextHeight, systemFont, 'l');

  g1 = HiddenGroup (wdialog, 2, 0, NULL);
  p = StaticPrompt (g1, "Matrix", 0, dialogTextHeight, systemFont, 'c');
  popf = PopupList (g1, TRUE, SelectMatrixPopup);
  PopupItem (popf, "identity");
  PopupItem (popf, "Blosum30");
  PopupItem (popf, "Blosum40");
  PopupItem (popf, "Blosum62");
  PopupItem (popf, "Blosum80");
  PopupItem (popf, "pam60");
  PopupItem (popf, "pam120");
  PopupItem (popf, "pam250");
  PopupItem (popf, "pam350");
  PopupItem (popf, "gonnet80");
  PopupItem (popf, "gonnet250");
  PopupItem (popf, "gonnet350");
  SetValue  (popf, dbdp->matrix);

  g2 = HiddenGroup (wdialog, 2, 0, NULL);
  sg1 = HiddenGroup (g2, 0, 1, NULL);
  p = StaticPrompt (sg1, "Separator", 0, dialogTextHeight, systemFont, 'c');  
  sg2 = HiddenGroup (g2, 0, 3, NULL);
  RadioButton (sg2, "space");
  RadioButton (sg2, "comma");
  RadioButton (sg2, "tab");
  SetValue (sg2, dbdp->separator);
  g3 = HiddenGroup (wdialog, 1, 0, NULL);
  PushButton (g3, "Dismiss", StdCancelButtonProc);

  RealizeWindow (wdialog);
  Show (wdialog);
  return;
}
******************************************************************/

/*******************************************************
***
***   FIND functions: John Kuzio function
***
*******************************************************/
#define FIND_MAXRANGE 10000

static ValNodePtr JK_NTPattern (CharPtr pattern, ValNodePtr sqloc_list, Boolean flagInvert, Uint1 mol_type)
{
  ValNodePtr      patlist = NULL;
  ValNodePtr      vnp;
  SeqLocPtr       slp, slptmp,
                  next;
  SeqIdPtr        sip;
  SeqIntPtr       sit;
  Uint1Ptr        sequence;
  ComPatPtr       cpp;
  Int4            from, to;
  Int4            new_lens;
  Int2            jk_moltype;
  SeqAlignPtr     salp;

  if (pattern == NULL)
     return NULL;
  CleanPattern (pattern);
  if (StringLen (pattern) == 0)
     return NULL;
  if (ISA_na(mol_type)) jk_moltype = 0;
  else jk_moltype = 1;
  cpp = CompilePattern (pattern, jk_moltype);
  if (cpp != NULL && flagInvert)
     cpp = InvertPattern (cpp);
  if (cpp == NULL)
     return NULL;
     
  for (vnp = sqloc_list; vnp != NULL; vnp = vnp->next)
  {
     slp = (SeqLocPtr)vnp->data.ptrvalue; 
     sip = SeqLocId(slp);
     from = SeqLocStart(slp);
     to = SeqLocStop(slp);
     sequence = (Uint1Ptr)load_seq_data (sip, from, to, (Boolean)jk_moltype, &new_lens);
     if (sequence != NULL) {
        salp = PatternMatch (sequence, 0, Seq_strand_plus, sip, cpp, 0, 0, TRUE);
        if (salp != NULL) {
           slp = MatchSa2Sl (&salp);
           slptmp=slp;
           while (slptmp!=NULL) {
              next = slptmp->next;
              slptmp->next = NULL;
              if (from > 0)
              {
                 sit = (SeqIntPtr)slptmp->data.ptrvalue;
                 sit->from += from;
                 sit->to += from;
              }
              ValNodeAddPointer (&patlist, 0, (Pointer)slptmp);
              slptmp = next;
           }
        }
        MemFree (sequence);
     }
   }
  cpp = ComPatFree (cpp);
  return patlist;
}

typedef struct ccid {
  SeqIdPtr    sip;
  SeqEntryPtr sep;
} CcId, PNTR CcIdPtr;
 
static void FindSeqEntryForSeqIdCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  CcIdPtr            cip;
 
  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     cip = (CcIdPtr)mydata;
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL && ISA_na (bsp->mol)) {
           sip = SeqIdFindBest(bsp->id, 0);
           if (SeqIdForSameBioseq(cip->sip, sip))
              cip->sep = sep;
        }
     }   
  }
}
 
static Int2 CC_SeqEntryToGeneticCode (Uint2 entityID, SeqIdPtr sip)
{
  SeqEntryPtr sep_head,
              sep;
  CcId        ci;
  Int2        genCode = 0;

  sep_head  = GetTopSeqEntryForEntityID (entityID);
  ci.sip = SeqIdDup (sip);
  ci.sep = NULL;
  SeqEntryExplore(sep_head,(Pointer)&ci, FindSeqEntryForSeqIdCallback);
  sep = ci.sep;
  SeqIdFree (ci.sip);
  if (sep!=NULL)
     genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
  return genCode;
}

static ValNodePtr JK_NTPattern2 (CharPtr pattern, ValNodePtr sqloc_list, Boolean flagInvert, Uint2 entityID)
{
  ByteStorePtr    bsp;
  ValNodePtr      patlist = NULL;
  ValNodePtr      vnp;
  SeqAlignPtr     salp;
  SeqLocPtr       slp, slptmp,
                  next,
                  slp2,
                  slp3;
  SeqIdPtr        sip;
  Uint1Ptr        sequence;
  ComPatPtr       cpp;
  Int4            from, 
                  to;
  Int2            jk_moltype;
  Int2            genCode = Seq_code_ncbieaa;
  Uint1           frame;

  if (pattern == NULL)
     return NULL;
  CleanPattern (pattern);
  if (StringLen (pattern) == 0)
     return NULL;
  jk_moltype = 1;
  cpp = CompilePattern (pattern, jk_moltype);
  if (cpp != NULL && flagInvert)
     cpp = InvertPattern (cpp);
  if (cpp == NULL)
     return NULL;
     
  for (vnp = sqloc_list; vnp != NULL; vnp = vnp->next)
  {
     slp = (SeqLocPtr)vnp->data.ptrvalue; 
     sip = SeqLocId(slp);
     genCode = CC_SeqEntryToGeneticCode (entityID, sip);
     if (genCode == 0)
        genCode = Seq_code_ncbieaa;
     from = SeqLocStart(slp);
     to = SeqLocStop(slp);
     slp2 = fuzz_loc (from, to, Seq_strand_plus, sip, TRUE, TRUE);
     if (slp2!=NULL) {
      for (frame=1; frame<4; frame++) 
      {
        bsp = cds_to_pept (slp2, frame, genCode, TRUE);
        if (bsp!=NULL) {
         sequence = (Uint1Ptr) BSMerge (bsp, NULL);
         BSFree (bsp);
         if (sequence != NULL) {
           salp = PatternMatch (sequence, 0, Seq_strand_plus, sip, cpp, 0, 0, TRUE);
           if (salp!=NULL) {
              slp = MatchSa2Sl (&salp);
              slptmp=slp;
              while (slptmp!=NULL) {
                 next = slptmp->next;
                 slptmp->next = NULL;
                 slp3 = SeqLocIntNew ((Int4)(from + 3*SeqLocStart(slptmp)+frame-1), (Int4)(from + 3*SeqLocStop(slptmp)+frame+1), SeqLocStrand(slptmp), SeqLocId(slptmp));
                 ValNodeFree (slptmp);
                 ValNodeAddPointer (&patlist, 0, (Pointer)slp3);
                 slptmp = next;
              }
           }
           MemFree (sequence);
         }
        }
      }
      ValNodeFree (slp2);
     }
  }
  cpp = ComPatFree (cpp);
  return patlist;
}

/******************************************************************/
static CharPtr get_string_fromtext (TexT txt)
{ 
  Char               str [128];
  CharPtr            strp;
  Int4               lens;
  Int4               j;

  GetTitle (txt, str, sizeof(str)-1);
  if (str == NULL)
     return NULL;
  lens = StringLen(str);
  if (lens < 1)
     return NULL;
  for (j=0, strp = str; j<lens; j++, strp++) {
     *strp = TO_UPPER(*strp);
  }
  strp = (CharPtr) MemNew ((StringLen(str)+4) * sizeof (Char));
  StringCpy (strp, str);
  return strp;
}

static void FindPatternButton (ButtoN b)
{
  WindoW             wdialog;
  WindoW             temport;
  DialogBoxDataPtr   dbdp;
  SeqEditViewFormPtr wdp;
  EditAlignDataPtr   adp;
  ValNodePtr         vnp;
  SelStructPtr       ssp;
  Char               str[128];
  CharPtr            strp;
  SeqParamPtr        prm;
  Uint1              frame;
  Int2               count = 0;
  Boolean            invert= FALSE,
                     translate= FALSE;
 
  wdialog = ParentWindow (b);
/****************
  Hide (wdialog);
*****************/
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  strp = get_string_fromtext (dbdp->txt1);
  if (strp == NULL) {
     SetTitle (dbdp->prompt, "No item found");
     return;
  }
  wdp = (SeqEditViewFormPtr) GetObjectExtra (dbdp->w);
  if ((adp = GetAlignEditData (dbdp->w)) == NULL )
     return;
  if (adp->current_pattern != NULL) {
     if ((StringStr(adp->current_pattern, strp)) && (StringLen(adp->current_pattern) == StringLen (strp)))
     {
        MemFree (strp);
        adp->cur_pat = ShowNextPattern (adp->match_pat, adp->cur_pat, adp->edit_pos);
        return;
     }
     else {
        MemFree (adp->current_pattern);
     }
  }
  adp->current_pattern = strp;
  if (adp->match_pat != NULL)
     adp->match_pat = SelStructDelList (adp->match_pat);
  if (dbdp->bt!=NULL)
     invert = (Boolean) GetStatus (dbdp->bt);
  if (dbdp->bt2!=NULL)
     translate = (Boolean) GetStatus (dbdp->bt2);
  if (translate && !invert) {
     vnp = JK_NTPattern2 (strp, adp->sqloc_list, invert, adp->input_entityID);
  }
  else 
     vnp = JK_NTPattern (strp, adp->sqloc_list, invert, adp->mol_type);
  if (vnp != NULL) {
     adp->match_pat = SetupMatchPatList (vnp, adp->anp_list);
     for (ssp=adp->match_pat; ssp!=NULL; ssp=ssp->next)
        count++;
     if (count ==0)
        SetTitle (dbdp->prompt, "No item found");
     else if (count ==1) {
        SetTitle (dbdp->prompt, "1 item found");
     } else {
        sprintf (str, "%d items found   -   next=CTRL->, previous=CTRL<- ", (int) count);
        SetTitle (dbdp->prompt, str);
     }
     if (translate) {
        for (frame=0; frame<10; frame++) {
           SetStatus (wdp->rfitem [frame], FALSE);
           for (vnp = adp->params; vnp != NULL; vnp = vnp->next) {
              prm = (SeqParamPtr) vnp->data.ptrvalue;
              prm->rf[frame] = FALSE;
           }  
        }
        for (ssp=adp->match_pat; ssp!=NULL; ssp=ssp->next) {
           frame = SeqLocStart((SeqLocPtr)ssp->region)%3;
           if (frame <8) {
              SetStatus (wdp->rfitem [frame], TRUE);
              for (vnp = adp->params; vnp != NULL; vnp = vnp->next) {
                 prm = (SeqParamPtr) vnp->data.ptrvalue;
                 prm->rf[frame] = TRUE;
              }
           }
        }   
        do_resize_window (wdp->pnl, adp, TRUE);
        temport = SavePort(dbdp->w);
        Select (wdp->pnl);
        inval_panel (wdp->pnl, -1 ,-1);
        RestorePort (temport);
        Update ();
     }
     adp->cur_pat = ShowNextPattern (adp->match_pat, NULL, adp->edit_pos);
  }
  else {
     SetTitle (dbdp->prompt, "No item found");
  }
/*****************
  Remove (wdialog);
******************/
  Update ();
  return;
}

extern void FindPatternDialog (IteM i)
{
  WindoW             w, wdialog;
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  GrouP              g1, g2, g3;

  w = ParentWindow (i);
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;
  wdialog=FixedWindow (-50, -33, -10, -10, "Find", StdCloseWindowProc);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  dbdp->bt = NULL;
  dbdp->bt2 = NULL;

  g1 = HiddenGroup (wdialog, 1, 0, NULL);
  StaticPrompt (g1, "Find pattern:", 0, popupMenuHeight, programFont, 'l');

  g1 = HiddenGroup (wdialog, 1, 0, NULL);
  if (adp->current_pattern == NULL) {
     dbdp->txt1 = ScrollText (g1, 25, 8, programFont, TRUE, NULL);
  } else {
     dbdp->txt1 = ScrollText (g1, 25, 8, programFont, TRUE, NULL);
  } 
  g3 = HiddenGroup (wdialog, 1, 0, NULL);
  dbdp->prompt = StaticPrompt (g3, "", (Int2)(25*stdCharWidth), stdLineHeight, programFont, 'l'); 
  if (ISA_na(adp->mol_type))
  {
     g2 = HiddenGroup (wdialog, 4, 0, NULL);
     dbdp->bt = CheckBox (g2, "reverse complement", NULL);  
     SetStatus (dbdp->bt, FALSE);
     dbdp->bt2 = CheckBox (g2, "translate sequence", NULL);
     SetStatus (dbdp->bt2, FALSE); 
     PushButton (g2, "Find Next", FindPatternButton);
     PushButton (g2, "Dismiss", StdCancelButtonProc);
  }
  else {
     g2 = HiddenGroup (wdialog, 2, 0, NULL);
     PushButton (g2, "Find Next", FindPatternButton);
     PushButton (g2, "Close", StdCancelButtonProc);
  }
  RealizeWindow (wdialog);
  Show (wdialog);
  return;
}

/******************************************************************/
static void GetInt4Value (TexT text)
{
  WindoW           wdialog;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Char             str[16];
 
  wdialog = ParentWindow (text);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ((adp = GetAlignEditData (dbdp->w)) == NULL ) return;
  GetTitle (text, str, 16);
  adp->int4value = (Int2) atoi (str);
  return;
}

static void GetInt4Value2 (TexT text)
{
  WindoW           wdialog;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Char             str[16];
 
  wdialog = ParentWindow (text);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ((adp = GetAlignEditData (dbdp->w)) == NULL ) return;
  GetTitle (text, str, 16);
  adp->int4value2 = (Int2) atoi (str);
  return;
}
/******************************************************************
***
***   DOT MATRIX
***
******************************************************************/
static int LIBCALLBACK callback(SeqAlignPtr seqalign)
{
  if (seqalign != NULL)
     LaunchAlignEditor (seqalign);
  return 1;
}
 
static void DotPlotBtn (ButtoN b)
{
  WindoW             wdialog;
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  SeqEntryPtr        sep;
  BioseqPtr          bsp;
  ValNodePtr         vnp;
  SeqLocPtr          slp = NULL;
  SeqIdPtr           sipnew = NULL;
  Int2               status;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  if (ISA_na(adp->mol_type) && adp->seqnumber == 1) {
     if (adp->align_format == SALSA_FASTA)
        sep = FastaRead (NULL, Seq_mol_na);
     else if (adp->align_format == SALSA_ASN1) {
        sep = AsnReadForSalsa (NULL);
        sep = getfirst_sep (sep, adp->mol_type);
     }   
     if (sep != NULL && IS_Bioseq (sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        sipnew = SeqIdFindBest(bsp->id, 0);   
        vnp = adp->sqloc_list;
        if (vnp!=NULL) {
           slp = (SeqLocPtr) vnp->data.ptrvalue;
           WatchCursor ();
           status = DotMatrixSearch(SeqLocId(slp), sipnew, callback);
           ArrowCursor ();
        }
     }
  } 
  else
     ErrPostEx (SEV_ERROR, 0, 0, "There should be 1 sequence in the editor");
  Remove (wdialog);
  Update ();
}

extern void DotPlotDialog (IteM i)
{
  WindoW           w, wdialog;
  DialogBoxDataPtr dbdp;
  GrouP            g1, g2;
  ButtoN           bt1, bt2;

  w = ParentWindow (i);
  wdialog = FixedWindow (-50, -33, -10, -10, "Dot Matrix", StdCloseWindowProc);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  g1 = HiddenGroup (wdialog, 3, 0, NULL);
  StaticPrompt (g1, "Input format:", 0, popupMenuHeight, systemFont, 'l');
  bt1 = RadioButton (g1, "Fasta");
  bt2 = RadioButton (g1, "ASN.1");
  Disable (bt2);
  SetValue (g1, 2);

  g2 = HiddenGroup (wdialog, 4, 0, NULL);
  PushButton (g2, "Process", DotPlotBtn);
  PushButton (g2, "Dismiss", StdCancelButtonProc);
  RealizeWindow (wdialog);
  Show (wdialog);
  return;
}

/******************************************************************/
extern void DotPlotItem (IteM i)
{
  EditAlignDataPtr   adp;
  SelStructPtr       ssp;
  AlignNodePtr       anp;
  SeqLocPtr          slp;
  SeqIdPtr           sip1 = NULL, 
                     sip2 = NULL;
  Int2               status;
  Int2               j;

  adp = GetAlignEditData ((WindoW) ParentWindow (i));
  if (adp->seqnumber == 0 ) return;
  if (GetNumberObjMgrSelect() != 2) {
     Message (MSG_OK, "Dot Matrix: Please select 2 sequences");
     return;
  }
  for (j=0, ssp = ObjMgrGetSelected(); j<2 && ssp!= NULL; ssp = ssp->next)
  {
     if ( checkssp_for_editor (ssp) ) {
        anp = (AlignNodePtr) AlignNodeFind (adp->anp_list, ssp->entityID, ssp->itemID, ssp->itemtype);
        if ( anp != NULL )
        {
           slp = CollectSeqLocFromAlignNode (anp);
           if (j==0) sip1 = SeqLocId (slp);
           else sip2 = SeqLocId (slp);
        }
     }
     j++;
  }
  if (sip1!=NULL && sip2!=NULL) {
     WatchCursor ();
     status = DotMatrixSearch(sip1, sip2, callback);
     ArrowCursor ();
  }
}

/******************************************
***  EXPORT FUNCTION
*******************************************/
static Boolean FindBspFromItem (GatherContextPtr gcp)

{
  BioseqPtr  PNTR bspp;

  bspp = (BioseqPtr PNTR) gcp->userdata;
  if (bspp != NULL && gcp->thistype == OBJ_BIOSEQ) {
    *bspp = (BioseqPtr) gcp->thisitem;
  }
  return TRUE;
}

extern void ExportTextFunc (ButtoN b)
{
  WindoW           wdialog;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  Char             namep [PATH_MAX];
  SelStructPtr     ssp,
                   ssphead = NULL;
  BioseqPtr        bsp;
  SeqAnnotPtr      sap;
  SeqAlignPtr      salp;
  SeqLocPtr        slp;
  AlignNodePtr     anp;
  ValNodePtr       vnp;
  Int4             start1, stop1;
  Int4             start=-1, 
                   stop=-1;
  FILE             *fout;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  if (!GetOutputFileName (namep, PATH_MAX, NULL))
     return;
  WatchCursor ();
  ssphead = NULL;
  if ( checkOMss_for_itemtype (OBJ_BIOSEQ) != 0 ) {
     ssphead = (Pointer)ObjMgrGetSelected ();
  }
  start1 = adp->int4value - 1;
  stop1 = adp->int4value2 - 1;   
  fout = FileOpen (namep, "w");
  if (fout != NULL) {
        if (adp->align_format == SALSA_SHWTXT) {
           ShowAlignmentText (fout, adp, ssphead, adp->marginleft, start1, stop1, FALSE);
        }
        else if (adp->align_format == SALSA_FASTA) {
           salp = (SeqAlignPtr)adp->sap_align->data;
           if (ssphead==NULL) {
              for (vnp = adp->anp_list; vnp != NULL; vnp = vnp->next) {
                 if ( (anp = (AlignNodePtr) vnp->data.ptrvalue) != NULL)
                 {
                    GatherItem (anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ, (Pointer) (&bsp), FindBspFromItem);
                    if (bsp!=NULL) {
                       start = AlignCoordToSeqCoord (start1, bsp->id, salp, adp->sqloc_list, 0);
                       stop = AlignCoordToSeqCoord (stop1, bsp->id, salp, adp->sqloc_list, 0);
                       EditBioseqToFasta (bsp, fout, (Boolean)ISA_na(adp->mol_type), start, stop);
                    }
                 }
              }
           } else for (ssp=ssphead;ssp!=NULL; ssp=ssp->next) {
              slp = (SeqLocPtr)ssp->region;
              if (slp!=NULL) {
                 if (GetAlignCoordFromSeqLoc(slp, salp, &start, &stop)) {
                    GatherItem (ssp->entityID, ssp->itemID, ssp->itemtype, (Pointer) (&bsp), FindBspFromItem);
                    if (bsp != NULL) {
                       start = AlignCoordToSeqCoord (start1, bsp->id, salp, adp->sqloc_list, 0);
                       stop = AlignCoordToSeqCoord (stop1, bsp->id, salp, adp->sqloc_list, 0);
                       EditBioseqToFasta (bsp, fout, (Boolean)ISA_na(adp->mol_type), start, stop);
                    }
                 }
              }
           }
        }
        else if (adp->align_format == SALSA_FASTGAP) {
           salp = (SeqAlignPtr)adp->sap_align->data;
           showfastagap_fromalign (salp, 60, fout);
        }
        else if (adp->align_format == SALSA_PHYLIP) {
           ShowAlignmentText (fout, adp, ssphead, 11, start1, stop1, FALSE);
        }
        else if (adp->align_format == SALSA_ASN1) {
           sap = SeqAnnotBoolSegToDenseSeg (adp->sap_align);
           if (sap != NULL) {
              seqannot_write (sap, namep);
              SeqAnnotFree (sap);
           }
        }
        FileClose (fout);
  }
  ArrowCursor ();
  Remove (wdialog);
  Update ();
}
 
extern void SalsaExportDialog (IteM i)
{
  WindoW           w, wdialog;
  DialogBoxDataPtr dbdp;
  GrouP            g;
  EditAlignDataPtr adp;
  Char             str [255], str2[10];
  CharPtr          tmp;
  SelStructPtr     ssp = NULL,
                   ssphead = NULL;
  ValNodePtr       vnp;
  AlignNodePtr     anp;
  Int4             j;
  Int2             selected_seq = 0;
  Boolean          first = TRUE;
 
  w = ParentWindow (i);
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;
  if ( checkOMss_for_itemtype (OBJ_BIOSEQ) == 0 ) {
     for (vnp = adp->anp_list; vnp != NULL; vnp = vnp->next) {
         if ( (anp = (AlignNodePtr) vnp->data.ptrvalue) != NULL)
         {
            ssp=SelStructNew(anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ, 0, 0, anp->sip, 0, FALSE);
            ssphead = (Pointer)SelStructAdd(ssphead, ssp);
         }
     }
  } else {
     ssphead = (Pointer)ObjMgrGetSelected ();
     selected_seq = checkOMss_for_itemtype (OBJ_BIOSEQ);
  }
  if (ssphead ==NULL) {
     return;
  }
  wdialog = FixedWindow (-50, -33, -10, -10, "Export", StdCloseWindowProc);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  g = HiddenGroup (wdialog, 2, 0, NULL);
  StaticPrompt (g, "Export:  ", 0, popupMenuHeight, systemFont, 'l');
  ssp=ssphead;
  while (first && ssp!=NULL) {
     tmp = str;
     for (j=0; j<10 && ssp!=NULL; j++, ssp=ssp->next) {
        SeqIdWrite (SeqLocId((SeqLocPtr)ssp->region), str2, PRINTID_REPORT, sizeof (str2));
        tmp = StringMove (tmp, str2);
        tmp = StringMove (tmp, " ");
     }
     if (ssp==NULL && !first)
        tmp = StringMove (tmp, " ...");
     if (!first)
        g = HiddenGroup (wdialog, 1, 0, NULL);
     else first = FALSE;
     StaticPrompt (g, str, 0, popupMenuHeight, programFont, 'l');
  }
  if (selected_seq == 1) {
     adp->int4value = SeqLocStart((SeqLocPtr)ssphead->region) +1;
     adp->int4value2 = SeqLocStop((SeqLocPtr)ssphead->region) +1;
  }
  else {
     adp->int4value = 1;
     adp->int4value2= adp->length;
  }
  if (selected_seq == 0)
     ssphead = SelStructDelList (ssphead);
  g = HiddenGroup (wdialog, 5, 0, NULL);
  StaticPrompt (g, "Positions:  ", 0, popupMenuHeight, systemFont, 'l');
  StaticPrompt (g, "from", 0, popupMenuHeight, programFont, 'l');
  sprintf (str2, "%ld", (long) adp->int4value);
  DialogText (g, str2, (Int2)6, GetInt4Value);
  StaticPrompt (g, "to", 0, popupMenuHeight, programFont, 'l');
  sprintf (str2, "%ld", (long) adp->int4value2);
  DialogText (g, str2, (Int2)6, GetInt4Value2);
  g = HiddenGroup (wdialog, 2, 0, NULL);
  PushButton (g, "Proceed", ExportTextFunc); 
  PushButton (g, "Cancel", StdCancelButtonProc);
  RealizeWindow (wdialog);
  Show (wdialog);
  return;
}

static void SelectVisibleSeq (CharPtr str, Int4 start, Int4 stop, Boolean show, EditAlignDataPtr adp)
{
  SeqParamPtr      prm;
  ValNodePtr       vnp,
                   curr;
  AlignNodePtr     anp;
  CharPtr          tmp,
                   tmp2;
  Int4             j=0;
/******* ???????????
  curr = adp->anp_list;
  for (vnp = adp->params; vnp != NULL; vnp = vnp = vnp->next) {
     anp = (AlignNodePtr) curr->data.ptrvalue;
     tmp = load_seq_data (anp->sip, start-1, stop+1, (Boolean)ISA_aa(adp->mol_type), &j);
     tmp2 = StringStr(tmp, str);
     if (show && tmp2==NULL) {
        prm = (SeqParamPtr) vnp->data.ptrvalue;
     }
     else if (!show && tmp2!=NULL) {
        prm = (SeqParamPtr) vnp->data.ptrvalue;
     }
     curr = curr->next;
  }
******/
  return; 
}

static void SelectSeqFunc (ButtoN b)
{
  WindoW           wdialog;
  WindoW           temport;
  DialogBoxDataPtr dbdp;
  SeqEditViewFormPtr wdp;
  EditAlignDataPtr adp;
  SelStructPtr     ssp;
  SeqLocPtr        slp;
  CharPtr          str;
  Int4             j;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (dbdp->w);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  ssp = (SelStructPtr) adp->extra_data;
  if (ssp==NULL)
     return;
  slp = (SeqLocPtr) ssp->region;
  if (slp==NULL) {
     adp->extra_data = NULL;
     return;
  }
  str = load_seq_data (SeqLocId(slp), SeqLocStart(slp), SeqLocStop(slp), (Boolean)ISA_aa(adp->mol_type), &j);
  adp->extra_data = NULL;
  SelectVisibleSeq(str,  SeqLocStart(slp), SeqLocStop(slp), dbdp->bool, adp);
  data_collect_arrange (adp, TRUE);
  temport = SavePort (wdp->pnl);
  Select (wdp->pnl);
  inval_panel (wdp->pnl, -1, -1);
  RestorePort (temport);
  Update ();
  Remove (wdialog);
  Update ();
  ArrowCursor();
  return;
}

static void HighlightSelectionFunc (ButtoN b)
{
  WindoW           wdialog;
  DialogBoxDataPtr dbdp;
  SeqEditViewFormPtr wdp;
  EditAlignDataPtr adp;
  SelStructPtr     ssp;
  Int4             j;
  Uint1Ptr         rgb;

  ValNodePtr       vnp;
  SeqLocPtr        slp;
  AlignNodePtr     anp;
  CharPtr          str,
                   tmp,
                   tmp2;
  SeqAlignPtr      salp;
  Int4             start, stop;
  Int4             start2, stop2;
  Int4             starta, stopa;
  Boolean          show;

  rgb=(Uint1Ptr)MemNew((size_t)(3*sizeof(Uint1)));
  rgb[0]=0;
  rgb[1]=255;
  rgb[2]=0;
  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (dbdp->w);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  ssp = (SelStructPtr) adp->extra_data;
  if (ssp==NULL)
     return;
  slp = (SeqLocPtr) ssp->region;
  if (slp==NULL) {
     adp->extra_data = NULL;
     return;
  }
  start = SeqLocStart(slp);
  stop = SeqLocStop(slp);
  str = load_seq_data (SeqLocId(slp), start, stop, (Boolean)ISA_aa(adp->mol_type), &j);
  adp->extra_data = NULL;
  show = dbdp->bool;
/*
  ObjMgrDeferUpdate ();
*/
  salp=(SeqAlignPtr) adp->sap_align->data;
  starta= SeqCoordToAlignCoord (start, SeqLocId(slp), salp, 0, 0);
  stopa = SeqCoordToAlignCoord (stop, SeqLocId(slp), salp, 0, 0);
  for (vnp = adp->anp_list; vnp != NULL; vnp = vnp = vnp->next) {
     anp = (AlignNodePtr) vnp->data.ptrvalue;
     start2 = AlignCoordToSeqCoord (starta, anp->sip, salp, adp->sqloc_list, 0);
     stop2 = AlignCoordToSeqCoord (stopa, anp->sip, salp, adp->sqloc_list, 0); 
     tmp = load_seq_data (anp->sip, start2-1, stop2+1, (Boolean)ISA_aa(adp->mol_type), &j);
     tmp2 = StringStr(tmp, str);
     if (show && tmp2!=NULL) {
        slp = SeqLocIntNew (start2, stop2, Seq_strand_plus, anp->sip);
        ObjMgrSetColor (anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ, OM_REGION_SEQLOC, (Pointer)slp, rgb);  
        SeqLocFree (slp);
     }
     else if (!show && tmp2==NULL) {
        slp = SeqLocIntNew (start2, stop2, Seq_strand_plus, anp->sip);
        ObjMgrSetColor (anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ, OM_REGION_SEQLOC, (Pointer)slp, rgb);  
        SeqLocFree (slp);
     }
  }
/*
  EnableUpdates ();
*/
  return;
}

extern void SelectSeqDialog (IteM i)
{
  WindoW           w, wdialog;
  DialogBoxDataPtr dbdp;
  GrouP            g,
                   g2, g3, g4, g5;
  EditAlignDataPtr adp;
  Char             str [64]; 
  CharPtr          tmp;
  SelStructPtr     ssphead = NULL;
  SeqLocPtr        slp;
  Int4             j;
 
  ssphead = (Pointer)ObjMgrGetSelected ();
  if (ssphead == NULL)
     return;
  if (ssphead->region == NULL)
     return;
  w = ParentWindow (i);
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;
  str[0] = '\0';
  tmp = str;
  adp->int4value = -1;
  adp->int4value2 = -1;
  slp = (SeqLocPtr)ssphead->region;
  tmp = load_seq_data (SeqLocId(slp), SeqLocStart(slp), SeqLocStop(slp), (Boolean)ISA_aa(adp->mol_type), &j);
  adp->int4value = SeqLocStart(slp) +1;
  adp->int4value2 = SeqLocStop(slp) +1;
  adp->extra_data = (Pointer) ssphead;

  wdialog = FixedWindow (-50, -33, -10, -10, "Export", StdCloseWindowProc);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  dbdp->bool = TRUE;

  g = HiddenGroup (wdialog, 3, 0, NULL);
  StaticPrompt(g, "Select the sequences ", 0, popupMenuHeight, systemFont,'l');

  g2 = HiddenGroup (g, 2, 0, BoolButton);
  RadioButton (g2, "showing");
  RadioButton (g2, "not showing");
  SetValue (g2, 1);

  g3 = HiddenGroup (wdialog, 2, 0, NULL);
  sprintf (str, "the pattern  %s", tmp);
  StaticPrompt (g3, str, 0, popupMenuHeight, systemFont, 'l');

  g4 = HiddenGroup (wdialog, 5, 0, NULL);
  sprintf(str, "from position  %d  to  %d", adp->int4value, adp->int4value2);
  StaticPrompt (g4, str, 0, popupMenuHeight, systemFont, 'l');

  g5 = HiddenGroup (wdialog, 3, 0, NULL);
  PushButton (g5, "Show", HighlightSelectionFunc); 
  PushButton (g5, "Select", SelectSeqFunc); 
  PushButton (g5, "Cancel", StdCancelButtonProc);
  RealizeWindow (wdialog);
  Show (wdialog);
  return;
}
static void SelectRegionFunc (ButtoN b)
{
  WindoW           w;      /* window from where the region is selected */
  WindoW           wdialog;/* window of the dialog box */
  DialogBoxDataPtr dbdp;   /*        extra data */
  EditAlignDataPtr adp;
  SeqAlignPtr      salp;
  Char             str[16];
  Int4             left, right;

  wdialog = (WindoW) ParentWindow (b);
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  w = dbdp->w;
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;
  GetTitle (dbdp->txt2, str, 16);
  left = (Int4) atol (str) - 1;
  GetTitle (dbdp->txt3, str, 16);
  right = (Int4) atol (str);
  if (left > 0 && right < adp->length)
  {
     Hide (wdialog);
     Update();
     if ( adp->seqnumber == 0 ) return;
     salp = SeqAlignTrunc ((SeqAlignPtr) adp->sap_align->data, left, right);
     if (salp != NULL) {
        LaunchAlignEditor (salp);
     }
     Remove (wdialog);
  }
  return;    
}


extern void SelectRegionDialog (IteM i)
{
  WindoW           w, wdialog;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  GrouP            g;
  Char             str[16];
 
  w = (WindoW) ParentWindow (i);
  wdialog = FixedWindow (-50, -33, -10, -10, "Sequence", NULL);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;

  g = HiddenGroup (wdialog, 1, 0, NULL);
  StaticPrompt (g, "Select positions in the alignment:", 0, dialogTextHeight, systemFont, 'l');

  g = HiddenGroup (wdialog, 4, 0, NULL);
  StaticPrompt (g, "from ", 0, dialogTextHeight, systemFont, 'l');
  sprintf (str, "%ld", (long) ( adp->gr.left + 1 ));
  dbdp->txt2 = DialogText (g, str, 8, NULL);
  StaticPrompt (g, "to ", 0, dialogTextHeight, systemFont, 'l');
  sprintf (str, "%ld", (long) adp->gr.right);
  dbdp->txt3 = DialogText (g, str, 8, NULL);

  g = HiddenGroup (wdialog, 2, 0, NULL);
  PushButton (g, "OK", SelectRegionFunc);
  PushButton (g, "Cancel", StdCancelButtonProc);
  Show (wdialog);
  return;
}



static void ColorIdentity (ButtoN b)
{
  WindoW           w;      /* window from where the region is selected */
  WindoW           wdialog;/* window of the dialog box */
  WindoW           temport;
  DialogBoxDataPtr dbdp;   /*        extra data */
  SeqEditViewFormPtr wdp;
  EditAlignDataPtr adp;
  CharPtr PNTR     buf;
  ValNodePtr       vnp;
  SelEdStructPtr   sesp;
  EditCellPtr      ecp;
  Int2             seqnumber;
  Int2             quorump[256];
  Int2             quorum;
  Int4             pos;
  Int2             k, kmax, max;

  wdialog = (WindoW) ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  w = dbdp->w;
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  if ( ( adp = GetAlignDataPanel (wdp->pnl) ) == NULL ) return;
  if ( adp->seqnumber == 0 || adp->edit_mode != PRETTY_EDIT) 
     return;
  for (k=0; k<10; k++) {
     adp->popcolor[k] = GetValue (dbdp->color[k]);
     if (adp->popcolor[k] == 1) {
        adp->col[k].r = adp->col[k].g = adp->col[k].b = 0;
     } else if (adp->popcolor[k] == 2) {
        adp->col[k].r = adp->col[k].g = adp->col[k].b = 0;
     } else if (adp->popcolor[k] == 3) {
        adp->col[k].r = 183; adp->col[k].g = 199; adp->col[k].b = 242;
     } else if (adp->popcolor[k] == 4) {
        adp->col[k].r = 0; adp->col[k].g = adp->col[k].b = 255;
     } else if (adp->popcolor[k] == 5) {
        adp->col[k].r = 255; adp->col[k].g = adp->col[k].b = 0;
     } else if (adp->popcolor[k] == 6) {
        adp->col[k].r = adp->col[k].g = 255; adp->col[k].b = 0;
     } else if (adp->popcolor[k] == 7) {
        adp->col[k].r = 0; adp->col[k].g = 255; adp->col[k].b = 0;
     }
  }
  buf = buf2array (adp->linebuff, 0, adp->seqnumber);
  if (buf!=NULL) { 
     for (pos=0; pos<adp->bufferlength; pos++) {
        for (k=0; k<256; k++) 
           quorump [k]=0;
        for (k=0; k<adp->seqnumber; k++) {
           quorump [(int)buf[k][pos]] ++;
        }
        max=0;
        for (k=0; k<256; k++) {
           if (k != '-' && quorump[k]>max) {
              max = quorump[k];
              kmax=k;
           }
        }
        seqnumber = adp->seqnumber - quorump[(int)('-')]; 
        if (seqnumber > 1 && max > 0) 
        {
           vnp=adp->copyalign;
           for (k=0; k<adp->seqnumber; k++) 
           {
              if (buf[k][pos] == (char)kmax) {
                 sesp = (SelEdStructPtr) vnp->data.ptrvalue;
                 ecp=(EditCellPtr)(sesp->data->data.ptrvalue); 
                 if (ecp!=NULL) { 
                  quorum = 1000*max/seqnumber;
                  if(quorum > 900) {
                    ecp[pos].r = adp->col[9].r; 
                    ecp[pos].g = adp->col[9].g; ecp[pos].b = adp->col[9].b;
                  }
                  else if(quorum > 800) {
                    ecp[pos].r = adp->col[8].r; 
                    ecp[pos].g = adp->col[8].g; ecp[pos].b = adp->col[8].b;
                  }
                  else if(quorum > 700) {
                    ecp[pos].r = adp->col[7].r; 
                    ecp[pos].g = adp->col[7].g; ecp[pos].b = adp->col[7].b;
                  }
                  else if(quorum > 600) {
                    ecp[pos].r = adp->col[6].r; 
                    ecp[pos].g = adp->col[6].g; ecp[pos].b = adp->col[6].b;
                  }
                  else if(quorum > 500) {
                    ecp[pos].r = adp->col[5].r; 
                    ecp[pos].g = adp->col[5].g; ecp[pos].b = adp->col[5].b;
                  }
                  else if(quorum > 400) {
                    ecp[pos].r = adp->col[4].r; 
                    ecp[pos].g = adp->col[4].g; ecp[pos].b = adp->col[4].b;
                  }
                  else if(quorum > 300) {
                    ecp[pos].r = adp->col[3].r; 
                    ecp[pos].g = adp->col[3].g; ecp[pos].b = adp->col[3].b;
                  }
                  else if(quorum > 200) {
                    ecp[pos].r = adp->col[2].r; 
                    ecp[pos].g = adp->col[2].g; ecp[pos].b = adp->col[2].b;
                  }
                  else if(quorum > 100) {
                    ecp[pos].r = adp->col[1].r; 
                    ecp[pos].g = adp->col[1].g; ecp[pos].b = adp->col[1].b;
                  }
                  else { 
                    ecp[pos].r = adp->col[0].r; 
                    ecp[pos].g =  adp->col[0].g; ecp[pos].b =  adp->col[0].b;
                  }
                 }  
              }
              vnp=vnp->next;
           }
        }
     }
  }
  temport = SavePort(w);
  Select (wdp->pnl);
  inval_panel (wdp->pnl, -1 ,-1);
  RestorePort (temport);
  Update ();
  Remove (wdialog);
  Update ();
}

extern void ColorIdentityDialog (WindoW w)
{
  WindoW           wdialog;
  DialogBoxDataPtr dbdp;
  EditAlignDataPtr adp;
  GrouP            g,
                   g2;
  Char             str[16];
  Int2             j, k, x;
 
  wdialog = FixedWindow (-50, -33, -10, -10, "Color Identity", NULL);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  if ( ( adp = GetAlignEditData (w) ) == NULL ) return;
  g = HiddenGroup (wdialog, 1,0,NULL);
  StaticPrompt (g, "Select a color for range of conservation", 0, dialogTextHeight, programFont, 'l');
  x = 0;
  g = HiddenGroup (wdialog, 6, 5, NULL);
  for (j=0; j<100; j += 20) {
     for (k=0; k<2; k++) {
        if (k==0)
           sprintf (str, "%d - %d%c", j, j+10, '%');
        else 
           sprintf (str, "  %d - %d%c", j+10, j+20, '%');
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');

        StaticPrompt (g, "[color]", 0, dialogTextHeight, programFont, 'l'); 
        dbdp->color[x] = PopupList (g, TRUE, NULL);
        PopupItem (dbdp->color[x], "none");
        PopupItem (dbdp->color[x], "new");
        PopupItem (dbdp->color[x], "marine");
        PopupItem (dbdp->color[x], "blue");
        PopupItem (dbdp->color[x], "red");
        PopupItem (dbdp->color[x], "yellow");
        PopupItem (dbdp->color[x], "green");
        SetValue  (dbdp->color[x], adp->popcolor[x]);
        x++;
     } 
  }
  g2 = HiddenGroup (wdialog, 2, 0, NULL);
  PushButton (g2, "Accept", ColorIdentity);
  PushButton (g2, "Dismiss", StdCancelButtonProc);
  Show (wdialog);
  return;
}

extern void ColorIdentityDialogItem (IteM i)
{
  ColorIdentityDialog ((WindoW) ParentWindow (i));
}

