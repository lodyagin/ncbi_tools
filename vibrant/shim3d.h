/*   shim3d.h
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
* File Name:  shim3d.h
*
* Author:  Lewis Geer
*
* Version Creation Date:   1/29/99
*
* $Revision: 6.2 $
*
* File Description: 
*  header file for shims to replace Viewer3D with OpenGL
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: shim3d.h,v $
* Revision 6.2  1999/04/06 17:17:17  lewisg
* fix typo
*
* Revision 6.1  1999/04/06 14:23:29  lewisg
* add opengl replacement for viewer3d
*

* ==========================================================================
*/

#ifndef _SHIM3D_
#define _SHIM3D_

#ifdef _OPENGL

#include <ncbilcl.h>
#include <ncbistd.h>
#include <vibmouse.h>

#ifdef __cplusplus
extern "C" {
#endif

#define OGL_DEFAULT_SIZE 5.0

#define OGLTEXT3D_MIDDLE 0x20
#define OGLTEXT3D_CENTER 0x2

#ifndef X_ROTATE_SBAR
#define X_ROTATE_SBAR 0x1
#endif

#ifndef Y_ROTATE_SBAR
#define Y_ROTATE_SBAR 0x2
#endif

#ifndef Z_ROTATE_SBAR
#define Z_ROTATE_SBAR 0x4
#endif

#define OGLMAXLAYERS 128
#define OGLFONTBASE 1000
#define OGLSELECTBUFFER 1000

typedef enum
{
    MouseOGL_DoNothing = 0,          
    MouseOGL_RotateYX,
    MouseOGL_RotateZX,
    MouseOGL_RotateYZ,
    MouseOGL_Move,
    MouseOGL_Zoom,        
    MouseOGL_NumStd
}  Nlm_enumStdMAOGL;

typedef Boolean (*Nlm_MAInitOGLFunc)(MAPtr ma, Nlm_VoidPtr data);


typedef struct _OGL_BoundBox 
/* bounds a volume */
{
	Nlm_FloatLo x[2];
	Nlm_FloatLo y[2];
	Nlm_FloatLo z[2];
	Nlm_Boolean set;
} TOGL_BoundBox;

typedef struct _OGL_ColorCell
/* contains a color.  duplicates ResidueColorCell, but we don't need the dependency */
{
    Nlm_Uint1 rgb[3];
} TOGL_ColorCell;

typedef struct _OGL_PaletteIndex
/* points into the palette created for OpenGL running in Index color mode */
{
    TOGL_ColorCell ColorCell;  /* the color */
    Nlm_Int4 Begin;  /* the beginning index into the palette */
    Nlm_Int4 End;  /* the end of the pallette */
} TOGL_PaletteIndex;

typedef struct _TOGL_Layers TOGL_Layers;

typedef struct _OGL_Data
/* general runtime information */
{
	Nlm_FloatLo Translate[3];			/* current translation */
    Nlm_FloatLo ViewPoint[3];           /* the current viewpoint. used for zoom and move */
    Nlm_FloatLo MaxSize;                /* biggest side of the bound box */
    TOGL_BoundBox BoundBox;         /* the containing box of the molecule */
    Nlm_WindoW ParentWindow;        
    Nlm_PaneL Panel;                /* needed to set palette */
    ValNodePtr PaletteExpanded;             /* the palette itself */
    ValNodePtr PaletteIndex;        /* palette index. type is TOGL_PaletteIndex */
    Nlm_Boolean NewPalette;  /* is the palette new? */
    Nlm_Int2 Width;   /* width and height of image */
    Nlm_Int2 Height;
    Nlm_BaR Z_rotate;  /* z rotation scroll bar */
    MAPtr       ma;
    MA_GroupPtr ma_std_group[MouseOGL_NumStd];
    Nlm_VoidPtr ModelMatrix;  /* temporary copy of modelview matrix */
    Nlm_Boolean IsPlaying;  /* is the animation running? */
    TOGL_Layers * Layers;  /* the layers and their state */
    Nlm_Int4 SpaceWidth;  /* width of space character */
    Nlm_Int4 SpaceHeight;  /* height of space character */
    Nlm_VoidPtr Space;  /* bitmap containing a space */
    Nlm_Boolean IndexMode;  /* number of bits per pixel.  If < 16, used color index mode */
    TOGL_ColorCell Background;  /* background color */
    /* note that the highlight color is the responsibility of the application */
    Nlm_Boolean SelectMode;  /* are we doing a selection? */
    Nlm_VoidPtr SelectBuffer;  /* buffer where the selection are dumped */
    Nlm_PoinT SelectPoint;  /* the point on the screen that was clicked for selection */
    Nlm_Uint4 SelectHits;  /* the number of hits */
} TOGL_Data;



extern Nlm_FloatHi * OGL_CrossProduct( Nlm_FloatHi * v1, Nlm_FloatHi * v2);
extern void OGL_Normalize( Nlm_FloatHi * v);
extern Nlm_FloatHi * OGL_MakeNormal (Nlm_FloatHi * origin, Nlm_FloatHi * v1, Nlm_FloatHi * v2);
extern void OGL_Reset (TOGL_Data * OGL_Data);
extern void OGL_AddQuad3D (TOGL_Data * OGL_Data, TOGL_ColorCell * color, Nlm_FloatHi * v1, Nlm_FloatHi * v2,
                                Nlm_FloatHi * v3, Nlm_FloatHi * v4);
extern void OGL_ClearBoundBox(TOGL_BoundBox *);
extern void OGL_ClearOGL_Data(TOGL_Data *);
extern void OGL_DrawViewer3D(TOGL_Data *);
extern void OGL_AddCylinder3D(TOGL_Data * OGL_Data, TOGL_ColorCell * color,
                                    Nlm_FloatHi x1, Nlm_FloatHi y1, Nlm_FloatHi z1, Nlm_FloatHi x2,
                                    Nlm_FloatHi y2, Nlm_FloatHi z2, Nlm_FloatHi radius);
extern void OGL_AddLine3D(TOGL_Data * OGL_Data, TOGL_ColorCell * color,
                                Nlm_FloatHi x1, Nlm_FloatHi y1, Nlm_FloatHi z1,
                                Nlm_FloatHi x2, Nlm_FloatHi y2, Nlm_FloatHi z2);
extern void OGL_AddSphere3D(TOGL_Data * OGL_Data, TOGL_ColorCell * color,
                                  Nlm_FloatHi x, Nlm_FloatHi y, Nlm_FloatHi z, Nlm_FloatHi radius);
extern void OGL_AddText3D(TOGL_Data * OGL_Data, TOGL_ColorCell * color, Nlm_CharPtr string,
                                Nlm_FloatHi x, Nlm_FloatHi y, Nlm_FloatHi z, Nlm_Int2 flags);
extern Nlm_Boolean OGL_SetPosition3D(TOGL_Data * OGL_Data, Nlm_RectPtr rect);
extern void OGL_DataFromBound(TOGL_Data * OGL_Data);
extern void OGL_Start();
extern void OGL_End();
extern void OGL_SetLayers(TOGL_Data * OGL_Data, Nlm_Boolean Status);
extern void OGL_SetLayer(TOGL_Data * OGL_Data, Nlm_Int4 i, Nlm_Boolean Status);
extern Nlm_Boolean OGL_GetLayer(TOGL_Data * OGL_Data, Nlm_Int4 i);
extern void OGL_SetLayerTop3D(TOGL_Data * OGL_Data, Nlm_Int4 TopLayer);
extern void OGL_AllLayerOnProc(TOGL_Data * OGL_Data);
extern void OGL_RewindLayerProc(TOGL_Data * OGL_Data);
extern void OGL_PrevLayerProc(TOGL_Data * OGL_Data);
extern void OGL_NextLayerProc(TOGL_Data * OGL_Data);
extern void OGL_Play(TOGL_Data * OGL_Data);
extern ValNodePtr OGL_SearchPaletteIndex(ValNodePtr PaletteIndex, TOGL_ColorCell * ColorCell);
extern TOGL_Data * OGL_CreateViewer(Nlm_GrouP prnt, Nlm_Uint2Ptr width, Nlm_Uint2 height,
                                   Nlm_Int4 flags, Nlm_MenU ma_group_menu, Nlm_MenU ma_action_menu,
                                   Nlm_MAInitOGLFunc ma_init_func, Nlm_VoidPtr ma_init_data);
extern void OGL_LoadName(Nlm_VoidPtr PtrValue);
extern Nlm_VoidPtr OGL_Hit(TOGL_Data * OGL_Data);
extern void OGL_Select(TOGL_Data * OGL_Data, Nlm_Boolean SelectMode);
extern void OGL_SetSelectPoint(TOGL_Data * OGL_Data, Nlm_PoinT Point);
extern void OGL_Redraw(TOGL_Data * OGL_Data);


#ifdef __cplusplus
}
#endif

#endif /* _OPENGL */

#endif

