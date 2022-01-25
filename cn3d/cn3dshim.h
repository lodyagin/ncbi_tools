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
* Revision 6.14  2000/01/06 00:04:43  lewisg
* selection bug fixes, update message outbound, animation APIs moved to vibrant
*
* Revision 6.13  2000/01/04 15:55:51  lewisg
* don't hang on disconnected network and fix memory leak/hang at exit
*
* Revision 6.12  1999/12/27 23:14:12  lewisg
* add colormgr show/hide in Cn3D
*
* Revision 6.11  1999/12/23 21:40:33  lewisg
* move animation controls to dialog
*
* Revision 6.10  1999/12/08 22:58:01  thiessen
* added quality settings for OpenGL rendering
*
* Revision 6.9  1999/12/03 23:17:23  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 6.8  1999/12/01 16:15:54  lewisg
* interim checkin to fix blocking memory leak
*
* Revision 6.7  1999/11/22 14:46:44  thiessen
* moved _OPENGL code blocks to only vibrant and ncbicn3d libraries
*
* Revision 6.6  1999/10/29 14:15:31  thiessen
* ran all Cn3D source through GNU Indent to prettify
*
* Revision 6.5  1999/10/18 15:32:50  lewisg
* move ClearSequences() to cn3dshim.c
*
* Revision 6.4  1999/10/05 23:18:25  lewisg
* add ddv and udv to cn3d with memory management
*
* Revision 6.3  1999/09/21 18:09:15  lewisg
* binary search added to color manager, various bug fixes, etc.
*
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

#ifndef _OPENGL
#include <viewer3d.h>
extern Viewer3D Cn3D_v3d;   /* the 3d view pane */
extern void /*LIBCALL*/ Cn3D_SaveActiveCam(void);
#endif

typedef struct _Cn3D_AnimateDlg {
    WindoW Cn3D_wAnimate;
} Cn3D_AnimateDlg;


typedef struct _Cn3D_ColorData {
    DDV_ColorGlobal *pDDVColorGlobal;
    Boolean IsUserData;     /* is the DDV_ColorGlobal userdata? */
    SeqAnnot *sap;          /* the current seqalign */
    ValNode *pvnsep;        /* the current seqentry */
    Boolean StandAlone;     /* is Cn3D running standalone? */
    Uint2 sapprocid, sepprocid, userkey;
#ifdef _OPENGL
    TOGL_Data *OGL_Data;    /* pointer to OGL data */
#endif
    WindoW Cn3D_w;
    Cn3D_AnimateDlg AnimateDlg;
    Boolean UseEntrez;  /* turn on entrez use */
    Boolean EntrezOn;  /* is entrez on? */
} TCn3D_ColorData;

extern TCn3D_ColorData Cn3D_ColorData; /* where all dynamic color info is kept */


/*****************************************************************************

Function: Cn3D_UseNetwork()

Purpose:  Determines if Cn3D should use the network
  
Returns:  TRUE if yes

*****************************************************************************/

NLM_EXTERN Boolean Cn3D_UseNetwork();


/*****************************************************************************

Function: Cn3D_SetVisible()

Purpose: Sets the visible bit for a chain in the biostruc AND the color
        manager
  
Parameters: pColorGlobal
            pmmdThis: the chain
            fVisible: true if the chain is to be visible

*****************************************************************************/

NLM_EXTERN void Cn3D_SetVisible(DDV_ColorGlobal *pColorGlobal, PMMD pmmdThis,
                                Boolean fVisible);


/*****************************************************************************

Function: Cn3D_IsVisible()

Purpose: Gets the visible bit for a chain in the biostruc.
  
Parameters: pColorGlobal
            pmmdThis: the chain

Returns:  TRUE if visible

*****************************************************************************/

NLM_EXTERN Boolean Cn3D_IsVisible(DDV_ColorGlobal *pColorGlobal,
                                   PMMD pmmdThis);


/*****************************************************************************

Function: ClearSequences()

Purpose: Deletes the Cn3D messagefunc from userdata on the SeqAnnots and
         SeqEntries presently displayed.
  
*****************************************************************************/
extern void ClearSequences(void);


extern void Cn3D_ConstructColorData(TCn3D_ColorData * ColorData
#ifdef _OPENGL
                                    , TOGL_Data * OGL_Data
#endif
                                    , Boolean StandAlone);
/* destructor for Color structure */
extern void Cn3D_DestructColorData(TCn3D_ColorData * ColorData);

#ifdef _OPENGL
extern void LIBCALL Cn3D_Size(TOGL_BoundBox * BoundBox,
                              PDNMS pdnmsThis);
NLM_EXTERN void Cn3D_Animate(IteM i);
extern Nlm_GrouP LIBCALL OGL_Quality(Nlm_GrouP prnt);
#endif                          /* _OPENGL */

extern void fnCHLresidue(PDNMG pdnmgThis,
#ifndef _OPENGL
                         Viewer3D vvv,
#endif
                         Boolean highlight);

#ifdef __cplusplus
}
#endif
#endif                          /* _CN3DSHIM_ */