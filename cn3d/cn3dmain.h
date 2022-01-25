/*   cn3dmain.h
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
* File Name:  cn3dmain.h
*
* Author:  Christopher Hogue
*
* Version Creation Date:   1/31/96
*
* $Revision: 6.11 $
*
* File Description: Main entry point for Cn3d 
*                   
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: cn3dmain.h,v $
* Revision 6.11  1999/04/06 16:06:11  coremake
* Fixes of typos
*
* Revision 6.10  1999/04/06 14:23:30  lewisg
* add opengl replacement for viewer3d
*
* Revision 6.9  1999/02/11 18:48:15  lewisg
* delete color index functions
*
* Revision 6.8  1999/02/10 23:49:43  lewisg
* use RGB values instead of indexed palette
*
* Revision 6.7  1999/01/14 19:07:17  kans
* network availability is configurable
*
* Revision 6.6  1998/08/26 18:28:37  kans
* fixed -v -fd warnings
*
* Revision 6.5  1998/04/28 15:14:27  lewisg
* moved OpenMimeFileWithDeletion to cn3dopen
*
* Revision 6.4  1998/04/27 23:23:04  lewisg
* added ability to open mime files
*
* Revision 6.3  1998/04/20 18:36:10  lewisg
* moved extern for Viewer3d to cn3dmain.h
*
* Revision 6.2  1998/03/06 23:19:16  lewisg
* codewarrior fixes
*
* Revision 6.1  1998/03/06 01:19:47  lewisg
* merge
*
* Revision 6.0  1997/08/25 18:13:33  madden
* Revision changed to 6.0
*
* Revision 5.4  1997/07/29 21:17:09  vakatov
* [WIN32,DLL]  Made Cn3D's stubbed functions be DLL-exportable
*
* Revision 5.3  1997/03/20 19:03:44  vakatov
* Non-standalone(Entrez-specific) protos has been moved to "cn3dentr.h":
* Cn3D_SetQueryCallback() and Cn3DWin(), and the latter has been renamed
* to Cn3DWin_Entrez().
* The remainig protos concern to the "cn3dwin.c" code.
*
 * Revision 5.2  1996/07/29  21:13:02  epstein
 * add prototype for function Cn3D_SetQueryCallback()
 *
 * Revision 5.1  1996/06/12  14:31:50  hogue
 * Added Cn3DWin function for integration into Entrez
 *
 * Revision 5.0  1996/05/28  14:05:44  ostell
 * Set to revision 5.0
 *
 * Revision 1.4  1996/04/18  16:56:37  hogue
 * Altered color palette for multi-structure display
 *
 * Revision 1.3  1996/03/30  23:40:45  hogue
 * Redraw now saves camera
 *
 * Revision 1.2  1996/03/29  20:00:45  hogue
 * Integrated 3d viewing, menus & controls for algorithmic rendering
 *
 * Revision 1.1  1996/02/01  18:47:38  kans
 * Initial revision
 *
*
* ==========================================================================
*/

#ifndef _CN3DMAIN_
#define _CN3DMAIN_

/* number of colors used for secondary structure */
#define CN3D_COLOR_SS 4
/* colors used for secondary structure */
#define CN3D_COLOR_HELIX 0
#define CN3D_COLOR_STRAND 1
#define CN3D_COLOR_TURN 2
#define CN3D_COLOR_COIL 3
/* total number of colors in fixed palette */
#define CN3D_COLOR_MAX 33
/* total number of colors allowed in indexed palettes. keep room for background */
#define CN3D_MAX_PALETTE 100

#ifdef _OPENGL
#include <vibrant.h>
#include <shim3d.h>
#else
#include <viewer3d.h>
#endif
#include <math.h>
#include <objalign.h>
#include <mmdbapi.h>
#include <salmedia.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

    /* standard color entry */
typedef struct _Cn3D_Color {
    ResidueColorCell ColorCell;  /* pointer to the color cell */
    Char * Name;  /* the name of the color */
    Char * Paints;  /* what the color paints */
    Int4 Index; /* an index assigned to the colors after they are allocated */
} TCn3D_Color;

typedef struct _Cn3D_ColorData {
    TCn3D_Color SSColors[CN3D_COLOR_SS];  /* secondary structure palette */
    TCn3D_Color Highlight;  /* this is the highlight color */
    ValNodePtr Palette;  /* valnode list of Cn3D_Colors.  type is TCn3D_Color */
#ifdef _OPENGL
    TOGL_Data * OGL_Data; /* pointer to OGL data */
#endif
} TCn3D_ColorData;

extern ResidueColorCell Cn3d_PaletteRGB[];   /* yanli */

extern TCn3D_ColorData Cn3D_ColorData;  /* where all dynamic color info is kept */
#ifndef _OPENGL
    extern Viewer3D Cn3D_v3d;  /* the 3d view pane */
#endif

extern void LIBCALL Cn3D_EnableFileOps(void);
extern void LIBCALL Cn3D_DisableFileOps(void);
extern void LIBCALL Cn3D_DisableMenus(void);
extern void LIBCALL Cn3D_EnableMenus(void);
extern Boolean LIBCALL readErrors (void);
#ifndef _OPENGL
extern void LIBCALL Cn3D_SaveActiveCam(void);
#endif
NLM_EXTERN void LIBCALL Cn3D_Redraw(Boolean New); 
NLM_EXTERN void LIBCALL Cn3D_ResetActiveStrucProc(void);
extern WindoW LIBCALL Cn3DWin(WndActnProc on_close, MenU *file_menu, ItmActnProc netconfig, Boolean usingEntrez);

/*extern void LaunchAlignViewer (SeqAlignPtr salp);*/
extern Boolean LIBCALL VASTInit (void);
extern BiostrucAnnotSetPtr LIBCALL VASTBsAnnotSetGet (Int4 uid);
extern void Cn3DResizeProc (WindoW w);

/* for salsa */
/*extern Int2 LIBCALLBACK AlgViewFunc (Pointer data); 
extern Int2 LIBCALLBACK AlgEditFunc (Pointer data);
extern Int2 LIBCALLBACK SeqEditFunc PROTO((Pointer data)); */

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
