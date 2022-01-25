/*   cn3dshim.h
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
* File Name:  cn3dshim.h
*
* Author: Lewis Geer
*
* Version Creation Date:   1/26/99
*
* File Description: Files use in Viewer3D -> OpenGL shim
*
* Modifications:
* $Log: cn3dshim.h,v $
* Revision 6.2  1999/04/06 14:23:30  lewisg
* add opengl replacement for viewer3d
*
* Revision 6.1  1999/02/10 23:49:44  lewisg
* use RGB values instead of indexed palette
*
*
* ===========================================================================  */

#ifndef _CN3DSHIM_
#define _CN3DSHIM_

#ifdef __cplusplus
extern "C" {
#endif

#include <shim3d.h>

/* constructor for TCn3D_Color structure */
extern void Cn3D_ConstructColor(TCn3D_Color * Color, Nlm_Uint1 Red, Nlm_Uint1 Green, Nlm_Uint1 Blue);
/* destructor for Color structure */
extern void Cn3D_DestructColor(TCn3D_Color * Color);

/* constructor for TCn3D_Color structure */
extern void Cn3D_ConstructColorData(TCn3D_ColorData * ColorData
#ifdef _OPENGL
                                    , TOGL_Data * OGL_Data
#endif
                                    );
/* destructor for Color structure */
extern void Cn3D_DestructColorData(TCn3D_ColorData * ColorData);

extern void Cn3D_FreePalette(ValNodePtr * Palette);
extern void Cn3D_SetColorChoice(ValNodePtr Palette);
extern ValNodePtr Cn3D_SearchColor(ValNodePtr Palette, ResidueColorCell * ColorCell);
extern void Cn3D_RequestColor(TCn3D_ColorData * ColorData, ResidueColorCell * ColorCell);
extern Nlm_Int4 Cn3D_ColorIndex(TCn3D_ColorData * ColorData, ResidueColorCell * ColorCell);
extern void Cn3D_CopyColorCell(ResidueColorCell * Destination,  ResidueColorCell * Source);
extern void Cn3D_SetColorCell(ResidueColorCell * Destination, Nlm_Uint1 Red, Nlm_Uint1 Green, Nlm_Uint1 Blue);
extern void Cn3D_DefaultSSColor(TCn3D_ColorData * ColorData);


#ifdef _OPENGL

extern void LIBCALL Cn3D_Size (TOGL_BoundBox * BoundBox, PDNMS pdnmsThis);
extern Nlm_GrouP LIBCALL  OGL_Controls ( Nlm_GrouP prnt);
extern Nlm_Boolean OGL_IsPlaying(void);
extern void OGL_StopPlaying(void);


#endif /* _OPENGL */

#ifdef __cplusplus
}
#endif


#endif /* _CN3DSHIM_ */
