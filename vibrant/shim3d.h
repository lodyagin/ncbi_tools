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
* $Revision: 6.21 $
*
* File Description:
*  header file for shims to replace Viewer3D with OpenGL
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: shim3d.h,v $
* Revision 6.21  2000/01/14 21:40:39  lewisg
* add translucent spheres, ion labels, new cpk, fix misc bugs
*
* Revision 6.20  2000/01/07 00:22:46  thiessen
* fixes for LessTif and OpenGL X visual selection
*
* Revision 6.19  2000/01/06 17:23:36  thiessen
* various OpenGL improvements
*
* Revision 6.18  2000/01/06 00:04:41  lewisg
* selection bug fixes, update message outbound, animation APIs moved to vibrant
*
* Revision 6.17  1999/12/17 20:25:01  thiessen
* put in preliminary PNG output (for Cn3D)
*
* Revision 6.16  1999/12/15 19:18:49  thiessen
* bug fix for previous revision
*
* Revision 6.14  1999/12/08 22:58:00  thiessen
* added quality settings for OpenGL rendering
*
* Revision 6.13  1999/11/23 19:24:16  thiessen
* better solution to OpenGL render rect setting on Mac
*
* Revision 6.12  1999/11/19 18:07:24  thiessen
* added label capability for OpenGL version on Mac and Motif
*
* Revision 6.11  1999/11/08 16:43:21  thiessen
* major rearrangement of OpenGL color/material/lighting; also added 3-d (thick) brick
*
* Revision 6.10  1999/10/31 22:39:34  thiessen
* added wifreframe worm capability to viewer3d
*
* Revision 6.9  1999/10/15 17:37:44  thiessen
* put in splined 'worm' model for virtual BB
*
* Revision 6.8  1999/09/22 14:59:21  thiessen
* another minor fix of forward decl. for SGI
*
* Revision 6.7  1999/09/22 14:22:49  lewisg
* fixed forward declaration of TOGL_Layers to compile on SGI
*
* Revision 6.6  1999/09/21 14:11:22  thiessen
* added #include <GL/gl.h>
*
* Revision 6.5  1999/09/21 13:45:31  thiessen
* port of Lewis's OpenGL code to X/Motif
*
* Revision 6.4  1999/09/20 20:12:56  lewisg
* change typedefs for a colorcell, add triangle generator, fix incorrect return values
*
* Revision 6.3  1999/06/14 23:15:11  lewisg
* moved useful helper functions out of the ifdef
*
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

#include <ncbilcl.h>
#include <ncbistd.h>
#include <vibmouse.h>
#include <ddvcolor.h>

#ifndef _OPENGL
#include <pictur3d.h>
#endif

/*

  Do not include gl.h here.  On Windows, gl.h is dependent on windows.h.
  However, some functions in vibrant have the same name as some of the
  functions in windows.h, so files that include both can have conflicts.

*/


#ifdef __cplusplus
extern "C" {
#endif

void  OGL_CreateCTransform(Nlm_FloatHi x1, Nlm_FloatHi y1, Nlm_FloatHi z1,
                           Nlm_FloatHi x2, Nlm_FloatHi y2, Nlm_FloatHi z2,
                           Nlm_FloatHi * Rotate, Nlm_FloatHi * Translate,
                           Nlm_FloatHi * length);

extern Nlm_FloatHi *OGL_CrossProduct(Nlm_FloatHi * v1, Nlm_FloatHi * v2);

#define OGL_SQR(oglx) ((oglx)*(oglx))


/* vibrant version of this function */
#ifndef _OPENGL
extern void Nlm_AddHalfWorm3D(Nlm_Picture3D pic,
                              Nlm_Segment3D segment, BigScalar userData,
                              Uint1 layer, Uint1 color,
                              Nlm_FloatHi x0, Nlm_FloatHi y0, Nlm_FloatHi z0,
                              Nlm_FloatHi x1, Nlm_FloatHi y1, Nlm_FloatHi z1,
                              Nlm_FloatHi x2, Nlm_FloatHi y2, Nlm_FloatHi z2,
                              Nlm_FloatHi x3, Nlm_FloatHi y3, Nlm_FloatHi z3,
                              Nlm_FloatHi radius, Nlm_Int4 segments);
#endif




#ifdef _OPENGL

/*
*  the following 2 includes are a subset of vibrant.h, since vibdefns.h
*  and vibforms.h confict with windows.h
*/

#include <vibtypes.h>
#include <vibprocs.h>


/* defines */

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

typedef enum {
    MouseOGL_DoNothing = 0,
    MouseOGL_RotateYX,
    MouseOGL_RotateZX,
    MouseOGL_RotateYZ,
    MouseOGL_Move,
    MouseOGL_Zoom,
    MouseOGL_NumStd
} Nlm_enumStdMAOGL;

typedef Boolean(*Nlm_MAInitOGLFunc) (MAPtr ma, Nlm_VoidPtr data);


typedef struct _OGL_BoundBox
/* bounds a volume */
{
    Nlm_FloatLo x[2];
    Nlm_FloatLo y[2];
    Nlm_FloatLo z[2];
    Nlm_Boolean set;
} TOGL_BoundBox;


typedef struct _OGL_PaletteIndex
/* points into the palette created for OpenGL running in Index color mode */
{
    DDV_ColorCell ColorCell;    /* the color */
    Nlm_Int4 Begin;             /* the beginning index into the palette */
    Nlm_Int4 End;               /* the end of the pallette */
} TOGL_PaletteIndex;

struct _TOGL_Layers;

typedef struct _TOGL_Layers TOGL_Layers_;


typedef struct _OGL_Data
/* general runtime information */
{
    /* the current viewpoint. used for zoom and move */
    Nlm_FloatLo CameraDistance; /* distance (on Z-axis) from origin */
    Nlm_FloatLo CameraAngle;    /* in radians */
    Nlm_FloatLo CameraDirection[2]; /* point in Z=0 plane camera points at */

    Nlm_FloatLo MaxSize;        /* biggest side of the bound box */
    TOGL_BoundBox BoundBox;     /* the containing box of the molecule */
    Nlm_WindoW ParentWindow;
    Nlm_PaneL Panel;            /* needed to set palette */
    ValNodePtr PaletteExpanded; /* the palette itself */
    ValNodePtr PaletteIndex;    /* palette index. type is TOGL_PaletteIndex */
    Nlm_Boolean NewPalette;     /* is the palette new? */
    Nlm_BaR Z_rotate;           /* z rotation scroll bar */
    MAPtr ma;
    MA_GroupPtr ma_std_group[MouseOGL_NumStd];
    Nlm_VoidPtr ModelMatrix;    /* temporary copy of modelview matrix */
    Nlm_Boolean IsPlaying;      /* is the animation running? */
    TOGL_Layers_ *Layers;       /* the layers and their state */
    Nlm_Int4 SpaceWidth;        /* width of space character */
    Nlm_Int4 SpaceHeight;       /* height of space character */
    Nlm_Boolean IndexMode;      /* number of bits per pixel.  If < 16, used color index mode */
    DDV_ColorCell Background;   /* background color */
    /* note that the highlight color is the responsibility of the application */
    Nlm_Boolean SelectMode;     /* are we doing a selection? */
    Nlm_VoidPtr SelectBuffer;   /* buffer where the selection are dumped */
    Nlm_PoinT SelectPoint;      /* the point on the screen that was clicked for selection */
    Nlm_Uint4 SelectHits;       /* the number of hits */

#ifdef WIN_MOTIF
    /* various info on X stuff related to OpenGL rendering context */
    /* these are void pointers, because including the X headers above this
       causes all sorts of name conflicts in various modules. Big pain! */
    void *display;              /* is actually a Display*      */
    void *visinfo;              /* is actually an XVisualInfo* */
#endif

} TOGL_Data;


extern void Nlm_AddHalfWorm3D(TOGL_Data * OGL_Data, DDV_ColorCell * color,
                              Nlm_FloatHi x0, Nlm_FloatHi y0, Nlm_FloatHi z0,
                              Nlm_FloatHi x1, Nlm_FloatHi y1, Nlm_FloatHi z1,
                              Nlm_FloatHi x2, Nlm_FloatHi y2, Nlm_FloatHi z2,
                              Nlm_FloatHi x3, Nlm_FloatHi y3, Nlm_FloatHi z3,
                              Nlm_FloatHi radius, 
                              Nlm_Int4 segments, Nlm_Int4 sides);

extern void OGL_Normalize(Nlm_FloatHi * v);
extern Nlm_FloatHi *OGL_MakeNormal(Nlm_FloatHi * origin, Nlm_FloatHi * v1,
                                   Nlm_FloatHi * v2);
extern void OGL_Reset(TOGL_Data * OGL_Data);
extern void OGL_AddQuad3D(TOGL_Data * OGL_Data, DDV_ColorCell * color,
                          Nlm_FloatHi * v1, Nlm_FloatHi * v2,
                          Nlm_FloatHi * v3, Nlm_FloatHi * v4);
extern void OGL_AddBrick3D(TOGL_Data * OGL_Data, DDV_ColorCell * color,
                          Nlm_FloatHi * v1, Nlm_FloatHi * v2,
                          Nlm_FloatHi * v3, Nlm_FloatHi * v4,
                          Nlm_FloatHi thickness);
extern void OGL_AddTri3D(TOGL_Data * OGL_Data, DDV_ColorCell * color,
                         Nlm_FloatHi * v1, Nlm_FloatHi * v2,
                         Nlm_FloatHi * v3, Nlm_FloatHi * Normal);

extern void OGL_ClearBoundBox(TOGL_BoundBox *);
extern void OGL_ClearOGL_Data(TOGL_Data *);
extern void OGL_DrawViewer3D(TOGL_Data *);
extern void OGL_AddCylinder3D(TOGL_Data * OGL_Data, DDV_ColorCell * color,
                              Nlm_FloatHi x1, Nlm_FloatHi y1,
                              Nlm_FloatHi z1, Nlm_FloatHi x2,
                              Nlm_FloatHi y2, Nlm_FloatHi z2,
                              Nlm_FloatHi radius, Nlm_Int4 sides);
extern void OGL_AddLine3D(TOGL_Data * OGL_Data, DDV_ColorCell * color,
                          Nlm_FloatHi x1, Nlm_FloatHi y1, Nlm_FloatHi z1,
                          Nlm_FloatHi x2, Nlm_FloatHi y2, Nlm_FloatHi z2);
extern void OGL_AddSphere3D(TOGL_Data * OGL_Data, DDV_ColorCell * color,
                            Nlm_FloatHi x, Nlm_FloatHi y, Nlm_FloatHi z,
                            Nlm_FloatHi radius, Nlm_Int4 slices, 
                            Nlm_Int4 stacks, Nlm_FloatHi alpha);
extern void OGL_AddText3D(TOGL_Data * OGL_Data, DDV_ColorCell * color,
                          Nlm_CharPtr string, Nlm_FloatHi x, Nlm_FloatHi y,
                          Nlm_FloatHi z, Nlm_Int2 flags);
extern Nlm_Boolean OGL_SetPosition3D(TOGL_Data * OGL_Data,
                                     Nlm_RectPtr rect);
extern void OGL_DataFromBound(TOGL_Data * OGL_Data);
extern void OGL_Start(TOGL_Data * OGL_Data, Nlm_Int1 List);
extern void OGL_End();
extern void OGL_SetLayers(TOGL_Data * OGL_Data, Nlm_Boolean status);
extern void OGL_SetLayer(TOGL_Data * OGL_Data, Nlm_Int4 i,
                         Nlm_Boolean status);
extern Nlm_Boolean OGL_GetLayer(TOGL_Data * OGL_Data, Nlm_Int4 i);
extern void OGL_SetLayerTop3D(TOGL_Data * OGL_Data, Nlm_Int4 TopLayer);
extern void OGL_AllLayerOnProc(TOGL_Data * OGL_Data);
extern void OGL_RewindLayerProc(TOGL_Data * OGL_Data);
extern void OGL_PrevLayerProc(TOGL_Data * OGL_Data);
extern void OGL_NextLayerProc(TOGL_Data * OGL_Data);
extern void OGL_Play(TOGL_Data * OGL_Data);
extern ValNodePtr OGL_SearchPaletteIndex(ValNodePtr PaletteIndex,
                                         DDV_ColorCell * ColorCell);
extern TOGL_Data *OGL_CreateViewer(Nlm_GrouP prnt, Nlm_Uint2Ptr width,
                                   Nlm_Uint2 height, Nlm_Int4 flags,
                                   Nlm_MenU ma_group_menu,
                                   Nlm_MenU ma_action_menu,
                                   Nlm_MAInitOGLFunc ma_init_func,
                                   Nlm_VoidPtr ma_init_data);
extern void OGL_LoadName(Nlm_VoidPtr PtrValue);
extern Nlm_VoidPtr OGL_Hit(TOGL_Data * OGL_Data);
extern void OGL_Select(TOGL_Data * OGL_Data, Nlm_Boolean SelectMode);
extern void OGL_SetSelectPoint(TOGL_Data * OGL_Data, Nlm_PoinT Point);
NLM_EXTERN Nlm_Boolean OGL_IsPlaying(TOGL_Data *pOGL_Data);
NLM_EXTERN void OGL_StopPlaying(TOGL_Data *pOGL_Data);
NLM_EXTERN void OGL_StartPlaying(TOGL_Data *pOGL_Data);

#ifdef _PNG
void Nlm_SaveImagePNG(Nlm_Char *defname); /* for PNG export */
#endif /* _PNG */

#endif /* _OPENGL */

#ifdef __cplusplus
}
#endif

#endif