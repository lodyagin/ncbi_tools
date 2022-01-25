/*   ncbidraw.h
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
* File Name:  ncbidraw.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/1/91
*
* $Revision: 6.5 $
*
* File Description: 
*       Vibrant drawing procedure definitions
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: ncbidraw.h,v $
* Revision 6.5  1998/07/14 16:44:26  vakatov
* Added VibrantIsGUI() and <internal> Nlm_VibrantSetGUI()
*
* Revision 6.4  1998/07/01 18:27:42  vakatov
* Use "const" qualifier somewhere
*
* Revision 6.3  1998/06/10 22:05:26  lewisg
* added PROTO() around *Quadrant()
*
* Revision 6.2  1998/06/01 17:20:13  vakatov
* Added code to draw/fill 90-degree arc/pie (quadrants)
*
* Revision 6.1  1997/09/16 20:24:53  vakatov
* Added Nlm_FitStringWidth()
*
* Revision 5.8  1997/07/23 19:42:22  vakatov
* Added Nlm_PaintStringEx function(text + position)
*
* Revision 5.7  1997/07/18 15:18:45  vakatov
* [!WIN_GIF] +Dummy Nlm_CreateGIF() and Nlm_SaveGIF();  /return FALSE/
*
* Revision 5.6  1997/06/18 22:10:10  vakatov
* [WIN_GIF] Moved Create/SaveGIF() decl. from "ncbigif.h" to "ncbidraw.h"
*
* Revision 5.2  1997/06/04 00:05:20  kans
* support for Japanese by Tomo Koike of DDBJ
*
* Revision 5.1  1996/09/30 19:56:00  vakatov
* + Nlm_SetPenDash PROTO and synopsis
*
* Revision 4.1  1996/05/10  21:14:21  vakatov
* Added UpdateColorTable() function to allow the user to read his color
* preferencies from an external config.-file
*
* Revision 2.12  1995/07/14  17:48:26  kans
* new CopyPixmap (AS)
*
* ==========================================================================
*/

#ifndef _NCBIDRAW_
#define _NCBIDRAW_

#ifdef __cplusplus
extern "C" {
#endif

/***  PORTABLE GRAPHIC PRIMITIVE OBJECT TYPES  ***/

/*
*  Platform independent point and rectangle data structures.
*/

typedef  struct  Nlm_point {
  Nlm_Int2  x;
  Nlm_Int2  y;
} Nlm_PoinT, PNTR Nlm_PointPtr;

typedef  struct  Nlm_rect {
  Nlm_Int2  left;
  Nlm_Int2  top;
  Nlm_Int2  right;
  Nlm_Int2  bottom;
} Nlm_RecT, PNTR Nlm_RectPtr;

typedef  Nlm_Handle  Nlm_RegioN;

typedef  struct  Nlm_font {
  Nlm_VoidPtr  dummy;
} HNDL Nlm_FonT;

/***  GLOBAL VARIABLES  ***/

/*
*  The update region contains the update region on drawing requests,
*  and contains at least the visible region during other callbacks.
*  Clipping is set to the update region on drawing events.  The
*  update rectangle has the minimum rectangular boundary of the
*  update region.  These variables are kept up to date by the Vibrant
*  event loop.
*/

extern  Nlm_RegioN  Nlm_updateRgn;
extern  Nlm_RecT    Nlm_updateRect;

/*
*  The standard systemFont and programFont variables are used to
*  specify fonts for dialog objects and scrolling text objects.
*/

extern  Nlm_FonT  Nlm_systemFont;
extern  Nlm_FonT  Nlm_programFont;

/*
*  Miscellaneous constants for pixel sizes of the standard font.
*/

extern  Nlm_Int2  Nlm_stdAscent;
extern  Nlm_Int2  Nlm_stdDescent;
extern  Nlm_Int2  Nlm_stdLeading;
extern  Nlm_Int2  Nlm_stdFontHeight;
extern  Nlm_Int2  Nlm_stdLineHeight;
extern  Nlm_Int2  Nlm_stdCharWidth;

/***  DRAWING PROCEDURES  ***/

void         Nlm_SetUpDrawingTools PROTO((void));
void         Nlm_CleanUpDrawingTools PROTO((void));

/*
*  It is not necessary to create a new font when switching colors.  The
*  family for GetFont can be Roman, Swiss, Modern, Script, or Decorative.
*/

/*
*  ScrollRect will erase and invalidate all invalid regions.
*/

void         Nlm_ResetDrawingTools PROTO((void));

void         Nlm_CopyMode PROTO((void));
void         Nlm_MergeMode PROTO((void));
void         Nlm_InvertMode PROTO((void));
void         Nlm_EraseMode PROTO((void));

void         Nlm_Black PROTO((void));
void         Nlm_Red PROTO((void));
void         Nlm_Green PROTO((void));
void         Nlm_Blue PROTO((void));
void         Nlm_Cyan PROTO((void));
void         Nlm_Magenta PROTO((void));
void         Nlm_Yellow PROTO((void));
void         Nlm_White PROTO((void));
void         Nlm_Gray PROTO((void));
void         Nlm_LtGray PROTO((void));
void         Nlm_DkGray PROTO((void));
void         Nlm_SelectColor PROTO((Nlm_Uint1 red, Nlm_Uint1 green, Nlm_Uint1 blue));
Nlm_Uint4    Nlm_GetColorRGB PROTO((Nlm_Uint1 red, Nlm_Uint1 green, Nlm_Uint1 blue));
Nlm_Uint4    Nlm_GetColor PROTO((void));
void         Nlm_SetColor PROTO((Nlm_Uint4 color));
void         Nlm_InvertColors PROTO((void));
void         Nlm_DecodeColor PROTO((Nlm_Uint4 color, Nlm_Uint1Ptr red, Nlm_Uint1Ptr green, Nlm_Uint1Ptr blue));

void         Nlm_Solid PROTO((void));
void         Nlm_Dark PROTO((void));
void         Nlm_Medium PROTO((void));
void         Nlm_Light PROTO((void));
void         Nlm_Empty PROTO((void));
void         Nlm_SetPenPattern PROTO((Nlm_VoidPtr pattern));
void         Nlm_Dotted PROTO((void));
void         Nlm_Dashed PROTO((void));
void         Nlm_WidePen PROTO((Nlm_Int2 width));

/* Under X11       -- full functionality(pen offset, dash length and gap length);
   Under Win-NT    -- parameter "offset" ignored(always zero)
   Other platforms -- exactly analogous to "Nlm_Dotted()"
   */
void Nlm_SetPenDash PROTO((Nlm_Uint1 offset, Nlm_Uint1 dash, Nlm_Uint1 gap));

/* esl++ */
/***  FONT HANDLING PROCEDURES  ***/

/* font styles */
#define STYLE_REGULAR       0
#define STYLE_BOLD          1
#define STYLE_ITALIC        2
#define STYLE_BOLD_ITALIC   3
#define STYLE_UNDERLINE     4
/* (other bits are used for platform-specific styles) */

/* font charset */
#define CHARSET_NULL        0
#define CHARSET_ANSI        1
#define CHARSET_SYMBOL      2
#define CHARSET_KANJI		3	/* Japanese Kanji */
#define CHASET_HANGUL		4	/* Korean character set */

/* font pitch */
#define PITCH_NULL          0
#define PITCH_FIXED         1
#define PITCH_VARIABLE      2

/* font family */
#define FAMILY_NULL         0
#define FAMILY_ROMAN        1
#define FAMILY_SWISS        2
#define FAMILY_MODERN       3
#define FAMILY_SCRIPT       4
#define FAMILY_DECORATIVE   5
#define FAMILY_GOTHIC		6	/* Japanese Kanji or Korean Hangle */
#define FAMILY_MINCHOU		7	/* Japanese Kanji */

#define FONT_NAME_SIZE      32

typedef struct Nlm_fontspec {
  Nlm_Char name [FONT_NAME_SIZE];
  Nlm_Int2 size;
  Nlm_Uint1 style;
  Nlm_Uint1 charset;
  Nlm_Uint1 pitch;
  Nlm_Uint1 family;
} Nlm_FontSpec, PNTR Nlm_FontSpecPtr;

Nlm_FonT     Nlm_CreateFont PROTO((Nlm_FontSpecPtr fsp));
Nlm_FonT     Nlm_GetResidentFont PROTO((Nlm_FontSpecPtr fsp));
Nlm_FonT     Nlm_CopyFont PROTO((Nlm_FonT font));
Nlm_FonT     Nlm_DeleteFont PROTO((Nlm_FonT font));
Nlm_FonT     Nlm_FindNextResidentFont PROTO((Nlm_FonT font));
Nlm_Boolean  Nlm_GetFontSpec PROTO((Nlm_FonT font, Nlm_FontSpecPtr fsp));
Nlm_Boolean  Nlm_EqualFontSpec PROTO((Nlm_FontSpecPtr fsp1, Nlm_FontSpecPtr fsp2));

/*
*  The functions below return the specifications for common fonts.
*  The return value points to the static buffer that is owerwritten
*  by next call to any of these functions.
*  Example: FonT f = CreateFont (Times (24, STYLE_BOLD_ITALIC));
*/
extern Nlm_FontSpecPtr Nlm_Helvetica PROTO((Nlm_Int2 size, Nlm_Uint1 style));
extern Nlm_FontSpecPtr Nlm_Times PROTO((Nlm_Int2 size, Nlm_Uint1 style));
extern Nlm_FontSpecPtr Nlm_Courier PROTO((Nlm_Int2 size, Nlm_Uint1 style));
extern Nlm_FontSpecPtr Nlm_Symbol PROTO((Nlm_Int2 size, Nlm_Uint1 style));
extern Nlm_FontSpecPtr Nlm_Gothic PROTO((Nlm_Int2 size, Nlm_Uint1 style));
extern Nlm_FontSpecPtr Nlm_Minchou PROTO((Nlm_Int2 size, Nlm_Uint1 style));
extern Nlm_FontSpecPtr Nlm_GothicFixed PROTO((Nlm_Int2 size, Nlm_Uint1 style));
extern Nlm_FontSpecPtr Nlm_MinchouFixed PROTO((Nlm_Int2 size, Nlm_Uint1 style));

/* esl++ end */

Nlm_FonT     Nlm_GetFont PROTO((Nlm_CharPtr name, Nlm_Int2 size, Nlm_Boolean bld, Nlm_Boolean itlc, Nlm_Boolean undrln, Nlm_CharPtr fmly));
Nlm_FonT     Nlm_GetFontEx PROTO((Nlm_CharPtr name, Nlm_Int2 size, Nlm_Boolean bld, Nlm_Boolean itlc, Nlm_Boolean undrln, Nlm_CharPtr fmly, Nlm_CharPtr chset, Nlm_Boolean fixed));
Nlm_FonT     Nlm_ParseFont PROTO((Nlm_CharPtr spec));
Nlm_FonT     Nlm_ParseFontEx PROTO((Nlm_CharPtr scrSpec, Nlm_CharPtr prtSpec));
void         Nlm_SelectFont PROTO((Nlm_FonT f));
void         Nlm_AssignPrinterFont PROTO((Nlm_FonT scrnFont, Nlm_FonT prtrFont));

Nlm_Int2     Nlm_CharWidth PROTO((Nlm_Char ch));
Nlm_Int2     Nlm_StringWidth PROTO((const Nlm_Char* text));
Nlm_Int2     Nlm_TextWidth PROTO((const Nlm_Char* text, size_t len));
Nlm_Int2     Nlm_Ascent PROTO((void));
Nlm_Int2     Nlm_Descent PROTO((void));
Nlm_Int2     Nlm_Leading PROTO((void));
Nlm_Int2     Nlm_FontHeight PROTO((void));
Nlm_Int2     Nlm_LineHeight PROTO((void));
Nlm_Int2     Nlm_MaxCharWidth PROTO((void));

/* Calculate(based on the presently active font) and return maximum
 * number of characters from the string "str" which can be fit into
 * "*max_width" pixels.
 * Return 0 if the "str" is NULL or empty or if the curr.font is NULL.
 */
size_t Nlm_FitStringWidth PROTO((const Nlm_Char PNTR str, Nlm_Int4 max_width));

void         Nlm_SetPen PROTO((Nlm_PoinT pt));
void         Nlm_GetPen PROTO((Nlm_PointPtr pt));

void         Nlm_PaintChar PROTO((Nlm_Char ch));
void         Nlm_PaintString PROTO((Nlm_CharPtr text));
void         Nlm_PaintStringEx PROTO((Nlm_CharPtr text, Nlm_Int2 x, Nlm_Int2 y));
void CDECL   Nlm_PaintText VPROTO((char *format, ...));

void         Nlm_DrawString PROTO((Nlm_RectPtr r, Nlm_CharPtr text, Nlm_Char jst, Nlm_Boolean gray));
void         Nlm_DrawText PROTO((Nlm_RectPtr r, Nlm_CharPtr text, size_t len, Nlm_Char jst, Nlm_Boolean gray));

void         Nlm_MoveTo PROTO((Nlm_Int2 x, Nlm_Int2 y));
void         Nlm_LineTo PROTO((Nlm_Int2 x, Nlm_Int2 y));
void         Nlm_DrawLine PROTO((Nlm_PoinT pt1, Nlm_PoinT pt2));

void         Nlm_LoadPt PROTO((Nlm_PointPtr pt, Nlm_Int2 x, Nlm_Int2 y));
void         Nlm_AddPt PROTO((Nlm_PoinT src, Nlm_PointPtr dst));
void         Nlm_SubPt PROTO((Nlm_PoinT src, Nlm_PointPtr dst));
Nlm_Boolean  Nlm_EqualPt PROTO((Nlm_PoinT p1, Nlm_PoinT p2));
Nlm_Boolean  Nlm_PtInRect PROTO((Nlm_PoinT pt, Nlm_RectPtr r));
Nlm_Boolean  Nlm_PtInRgn PROTO((Nlm_PoinT pt, Nlm_RegioN rgn));

void         Nlm_LoadRect PROTO((Nlm_RectPtr r, Nlm_Int2 lf, Nlm_Int2 tp, Nlm_Int2 rt, Nlm_Int2 bt));
void         Nlm_UpsetRect PROTO((Nlm_RectPtr r, Nlm_Int2 lf, Nlm_Int2 tp, Nlm_Int2 rt, Nlm_Int2 bt));
void         Nlm_OffsetRect PROTO((Nlm_RectPtr r, Nlm_Int2 dx, Nlm_Int2 dy));
void         Nlm_InsetRect PROTO((Nlm_RectPtr r, Nlm_Int2 dx, Nlm_Int2 dy));
Nlm_Boolean  Nlm_SectRect PROTO((Nlm_RectPtr src1, Nlm_RectPtr src2, Nlm_RectPtr dst));
Nlm_Boolean  Nlm_UnionRect PROTO((Nlm_RectPtr src1, Nlm_RectPtr src2, Nlm_RectPtr dst));
Nlm_Boolean  Nlm_EqualRect PROTO((Nlm_RectPtr r1, Nlm_RectPtr r2));
Nlm_Boolean  Nlm_EmptyRect PROTO((Nlm_RectPtr r));
Nlm_Boolean  Nlm_RectInRect PROTO((Nlm_RectPtr r1, Nlm_RectPtr r2));
Nlm_Boolean  Nlm_RectInRgn PROTO((Nlm_RectPtr r, Nlm_RegioN rgn));

void         Nlm_EraseRect PROTO((Nlm_RectPtr r));
void         Nlm_FrameRect PROTO((Nlm_RectPtr r));
void         Nlm_PaintRect PROTO((Nlm_RectPtr r));
void         Nlm_InvertRect PROTO((Nlm_RectPtr r));
void         Nlm_ScrollRect PROTO((Nlm_RectPtr r, Nlm_Int2 dx, Nlm_Int2 dy));

void         Nlm_EraseOval PROTO((Nlm_RectPtr r));
void         Nlm_FrameOval PROTO((Nlm_RectPtr r));
void         Nlm_PaintOval PROTO((Nlm_RectPtr r));
void         Nlm_InvertOval PROTO((Nlm_RectPtr r));

void         Nlm_EraseRoundRect PROTO((Nlm_RectPtr r, Nlm_Int2 ovlWid, Nlm_Int2 ovlHgt));
void         Nlm_FrameRoundRect PROTO((Nlm_RectPtr r, Nlm_Int2 ovlWid, Nlm_Int2 ovlHgt));
void         Nlm_PaintRoundRect PROTO((Nlm_RectPtr r, Nlm_Int2 ovlWid, Nlm_Int2 ovlHgt));
void         Nlm_InvertRoundRect PROTO((Nlm_RectPtr r, Nlm_Int2 ovlWid, Nlm_Int2 ovlHgt));

void         Nlm_EraseArc PROTO((Nlm_RectPtr r, Nlm_PoinT start, Nlm_PoinT end));
void         Nlm_FrameArc PROTO((Nlm_RectPtr r, Nlm_PoinT start, Nlm_PoinT end));
void         Nlm_PaintArc PROTO((Nlm_RectPtr r, Nlm_PoinT start, Nlm_PoinT end));
void         Nlm_InvertArc PROTO((Nlm_RectPtr r, Nlm_PoinT start, Nlm_PoinT end));

/* Special case of an arc(90 grad) */
typedef enum {
  eQ_RightTop    = 1,
  eQ_LeftTop     = 2,
  eQ_LeftBottom  = 3,
  eQ_RightBottom = 4
} Nlm_EQuadrant;
#define EQuadrant Nlm_EQuadrant

void         Nlm_EraseQuadrant  PROTO((Nlm_RectPtr r, Nlm_EQuadrant quadrant));
void         Nlm_FrameQuadrant  PROTO((Nlm_RectPtr r, Nlm_EQuadrant quadrant));
void         Nlm_PaintQuadrant  PROTO((Nlm_RectPtr r, Nlm_EQuadrant quadrant));
void         Nlm_InvertQuadrant PROTO((Nlm_RectPtr r, Nlm_EQuadrant quadrant));


void         Nlm_ErasePoly PROTO((Nlm_Int2 num, Nlm_PointPtr pts));
void         Nlm_FramePoly PROTO((Nlm_Int2 num, Nlm_PointPtr pts));
void         Nlm_PaintPoly PROTO((Nlm_Int2 num, Nlm_PointPtr pts));
void         Nlm_InvertPoly PROTO((Nlm_Int2 num, Nlm_PointPtr pts));

Nlm_RegioN   Nlm_CreateRgn PROTO((void));
Nlm_RegioN   Nlm_DestroyRgn PROTO((Nlm_RegioN rgn));
void         Nlm_ClearRgn PROTO((Nlm_RegioN rgn));
void         Nlm_LoadRectRgn PROTO((Nlm_RegioN rgn, Nlm_Int2 lf, Nlm_Int2 tp, Nlm_Int2 rt, Nlm_Int2 bt));
void         Nlm_OffsetRgn PROTO((Nlm_RegioN rgn, Nlm_Int2 dx, Nlm_Int2 dy));
void         Nlm_SectRgn PROTO((Nlm_RegioN src1, Nlm_RegioN src2, Nlm_RegioN dst));
void         Nlm_UnionRgn PROTO((Nlm_RegioN src1, Nlm_RegioN src2, Nlm_RegioN dst));
void         Nlm_DiffRgn PROTO((Nlm_RegioN src1, Nlm_RegioN src2, Nlm_RegioN dst));
void         Nlm_XorRgn PROTO((Nlm_RegioN src1, Nlm_RegioN src2, Nlm_RegioN dst));
Nlm_Boolean  Nlm_EqualRgn PROTO((Nlm_RegioN rgn1, Nlm_RegioN rgn2));
Nlm_Boolean  Nlm_EmptyRgn PROTO((Nlm_RegioN rgn));

void         Nlm_EraseRgn PROTO((Nlm_RegioN rgn));
void         Nlm_FrameRgn PROTO((Nlm_RegioN rgn));
void         Nlm_PaintRgn PROTO((Nlm_RegioN rgn));
void         Nlm_InvertRgn PROTO((Nlm_RegioN rgn));

void         Nlm_ClipRect PROTO((Nlm_RectPtr r));
void         Nlm_ClipRgn PROTO((Nlm_RegioN rgn));
void         Nlm_ResetClip PROTO((void));

void         Nlm_ValidRect PROTO((Nlm_RectPtr r));
void         Nlm_InvalRect PROTO((Nlm_RectPtr r));
void         Nlm_ValidRgn PROTO((Nlm_RegioN rgn));
void         Nlm_InvalRgn PROTO((Nlm_RegioN rgn));

void         Nlm_CopyBits PROTO((Nlm_RectPtr r, Nlm_VoidPtr source));

typedef struct {
  Nlm_Uint1 red;
  Nlm_Uint1 green;
  Nlm_Uint1 blue;
} Nlm_RGBColoR, PNTR Nlm_RGBColoRPtr;

void         Nlm_CopyPixmap PROTO((Nlm_RectPtr r, Nlm_Int1Ptr source, 
                                   Nlm_Uint1 totalColors, 
                                   Nlm_RGBColoRPtr colorTable));

/*
 * Try to read alternate color set from the user-specified file;
 * <table_index> <red> <green> <blue> [<name>].
 * Return number of updated colors
 */ 
size_t Nlm_UpdateColorTable PROTO((Nlm_RGBColoR table[], size_t table_size,
                                   const Nlm_Char PNTR filename));
#define UpdateColorTable Nlm_UpdateColorTable


Nlm_Boolean Nlm_CreateGIF PROTO((Nlm_Int2 width, Nlm_Int2 height,
                                 Nlm_Boolean transparent));
#define CreateGIF     Nlm_CreateGIF

Nlm_Boolean Nlm_SaveGIF PROTO((FILE PNTR out));
#define SaveGIF       Nlm_SaveGIF

/* If the application is using GUI or just drawing(picture) functionality */
Nlm_Boolean Nlm_VibrantIsGUI PROTO((void));
#define VibrantIsGUI Nlm_VibrantIsGUI

#ifdef __cplusplus
}
#endif

#endif

