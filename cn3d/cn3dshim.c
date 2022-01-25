/*   cn3dshim.c
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
* File Name:  cn3dshim.c
*
* Authors:  Lewis Geer
*
* Version Creation Date:
*
* File Description: Shim functions necessary to convert to OpenGL
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cn3dshim.c,v $
* Revision 6.22  2000/01/18 22:49:16  lewisg
* send OM_MSG_FLUSH to ddv/udv, tweak CPK coloration, misc bugs
*
* Revision 6.21  2000/01/08 00:47:53  lewisg
* fixes to selection, update, color
*
* Revision 6.20  2000/01/06 00:04:42  lewisg
* selection bug fixes, update message outbound, animation APIs moved to vibrant
*
* Revision 6.19  2000/01/04 15:55:51  lewisg
* don't hang on disconnected network and fix memory leak/hang at exit
*
* Revision 6.18  1999/12/29 22:55:03  lewisg
* get rid of seqalign id
*
* Revision 6.17  1999/12/28 15:08:44  lewisg
* remove remaining mediainfo code
*
* Revision 6.16  1999/12/27 23:14:12  lewisg
* add colormgr show/hide in Cn3D
*
* Revision 6.15  1999/12/23 21:40:33  lewisg
* move animation controls to dialog
*
* Revision 6.14  1999/12/02 20:31:59  lewisg
* put seqentries into bioseqset and fix calling convention in alignmgr.c
*
* Revision 6.13  1999/12/01 16:15:54  lewisg
* interim checkin to fix blocking memory leak
*
* Revision 6.12  1999/11/30 22:46:37  vakatov
* Cast callback arg in Cn3D_SizeCB() lest to cast the callback func.type
*
* Revision 6.11  1999/11/12 16:06:34  lewisg
* fix sequentry to valnode conversion
*
* Revision 6.10  1999/11/10 23:19:42  lewisg
* rewrite of selection code for ddv
*
* Revision 6.9  1999/11/02 23:06:08  lewisg
* fix cn3d to launch correctly if there is no seqentry associated with bioseq
*
* Revision 6.8  1999/10/29 14:15:30  thiessen
* ran all Cn3D source through GNU Indent to prettify
*
* Revision 6.7  1999/10/18 15:32:50  lewisg
* move ClearSequences() to cn3dshim.c
*
* Revision 6.6  1999/10/15 20:56:40  lewisg
* append DDV_ColorGlobal as userdata.  free memory when cn3d terminates.
*
* Revision 6.5  1999/10/05 23:18:24  lewisg
* add ddv and udv to cn3d with memory management
*
* Revision 6.4  1999/09/22 20:07:39  thiessen
* minor fix for mac compiler
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
*/


#include <cn3dmain.h>
#include <cn3dmsg.h>
#include <ncbi.h>
#include <cn3dshim.h>
#include <seqmgr.h>
#include <salpacc.h>

/*****************************************************************************

Function: Cn3D_UseNetwork()

Purpose:  Determines if Cn3D should use the network
  
Returns:  TRUE if yes

*****************************************************************************/

NLM_EXTERN Boolean Cn3D_UseNetwork()
{
    Char str[64];

    if (GetAppParam
        ("CN3D", "SETTINGS", "NETWORKAVAILABLE", NULL, str, sizeof(str))) {
        if (StringICmp(str, "TRUE") == 0) return TRUE;
    }
    return FALSE;
}

/*****************************************************************************

Function: Cn3D_SetVisible()

Purpose: Sets the visible bit for a chain in the biostruc AND the color
        manager
  
Parameters: pColorGlobal
            pmmdThis: the chain
            fVisible: TRUE if the chain is to be visible

*****************************************************************************/

NLM_EXTERN void Cn3D_SetVisible(DDV_ColorGlobal *pColorGlobal, PMMD pmmdThis,
                                Boolean fVisible)
{
    if(pmmdThis == NULL) return;
    if(pmmdThis->pSeqId == NULL) return;
    DDV_SetVisible(pColorGlobal, pmmdThis->pSeqId, -1, fVisible);
    pmmdThis->bVisible = fVisible;
}


/*****************************************************************************

Function: Cn3D_IsVisible()

Purpose: Gets the visible bit for a chain in the biostruc.
  
Parameters: pColorGlobal
            pmmdThis: the chain

Returns:  TRUE if visible

*****************************************************************************/

NLM_EXTERN Boolean Cn3D_IsVisible(DDV_ColorGlobal *pColorGlobal,
                                   PMMD pmmdThis)
{
    if(pmmdThis == NULL) return FALSE;
    return pmmdThis->bVisible;
}

/*****************************************************************************

Function: ClearSequences()

Purpose: Deletes the Cn3D messagefunc from userdata on the SeqAnnots and
         SeqEntries presently displayed.
  
*****************************************************************************/

void ClearSequences(void)
{
    Uint2 entityID;
    SeqAnnot *sap;
    SeqAlign *salp;
    
    if (Cn3D_ColorData.sap != NULL) {
        for (sap = Cn3D_ColorData.sap; sap != NULL; sap = sap->next) {        
            if (sap->data == NULL) continue;
            salp = sap->data;
            entityID = ObjMgrGetEntityIDForPointer(salp);
            ObjMgrFreeUserData(entityID, Cn3D_ColorData.sapprocid, OMPROC_EDIT,
                Cn3D_ColorData.userkey);               
            ObjMgrSendMsg(OM_MSG_FLUSH, entityID, 0, 0);
        }
    }
        
    entityID =
        ObjMgrGetEntityIDForPointer((void *) Cn3D_ColorData.pvnsep);
    ObjMgrFreeUserData(entityID, Cn3D_ColorData.sepprocid, OMPROC_VIEW,
        Cn3D_ColorData.userkey);
    ObjMgrSendMsg(OM_MSG_FLUSH, entityID, 0, 0);
    
    if(Cn3D_ColorData.IsUserData == FALSE)
        DDV_DeleteColorGlobal(Cn3D_ColorData.pDDVColorGlobal);
    Cn3D_ColorData.pDDVColorGlobal = NULL;
    
    return;
}


/* constructor for DDV_ColorEntry structure */
void Cn3D_ConstructColorData(TCn3D_ColorData * ColorData
#ifdef _OPENGL
                             , TOGL_Data * OGL_Data
#endif
                             , Boolean StandAlone)
{

    if (ColorData == NULL)
        return;
#ifdef _OPENGL
    if (OGL_Data == NULL)
        return;
#endif
    ColorData->sap = NULL;
    ColorData->pvnsep = NULL;
    ColorData->StandAlone = StandAlone;
    ColorData->pDDVColorGlobal = NULL;
    ColorData->IsUserData = FALSE;
#ifdef _OPENGL
    ColorData->OGL_Data = OGL_Data;
#endif
    ColorData->AnimateDlg.Cn3D_wAnimate = NULL;
    ColorData->Cn3D_w = NULL;
}

/* destructor for Color structure */
void Cn3D_DestructColorData(TCn3D_ColorData * ColorData)
{
    if (ColorData == NULL)
        return;
    DDV_DeleteColorGlobal(ColorData->pDDVColorGlobal);
}


/************************* OpenGL specific functions *****************************/

#ifdef _OPENGL


static void LIBCALLBACK Cn3D_SizeCB(PFB pfbThis, Nlm_Int4 iModel,
                                    Nlm_Int4 iIndex,
                                    Pointer ptr)
/* callback used to find the bounding box of the atoms of a structure */
{
    TOGL_BoundBox* BoundBox = (TOGL_BoundBox*) ptr;
    PMAD pmadThis = NULL;
    PALD paldThis = NULL;
    Nlm_FloatLoPtr pflvDataThis = NULL;
    Nlm_Int4 i;

    if (pfbThis == NULL || BoundBox == NULL)
        return;
    if (pfbThis->bMe == (Byte) AM_MAD) { /* is this molecular atom data? *//* iIndex isn't used */
        pmadThis = (PMAD) pfbThis;
        paldThis = GetAtomLocs(pmadThis, iModel);
        while (paldThis) {
            pflvDataThis = paldThis->pflvData;
            if (BoundBox->set == FALSE) { /* if the bounding box has not been used, set it */
                for (i = 0; i < 6; i++) {
                    BoundBox->x[i] = pflvDataThis[i / (int) 2];
                    BoundBox->set = TRUE;
                }
            } else {
                for (i = 0; i < 6; i += 2) { /* check to see if the atom lies outside the bound box */
                    if (pflvDataThis[i / (Nlm_Int4) 2] < BoundBox->x[i])
                        BoundBox->x[i] = pflvDataThis[i / (Nlm_Int4) 2];
                    if (pflvDataThis[i / (Nlm_Int4) 2] >
                        BoundBox->x[i + 1]) BoundBox->x[i + 1] =
                            pflvDataThis[i / (Nlm_Int4) 2];
                }
            }
            paldThis = paldThis->next;
        }
    }
}




void LIBCALL Cn3D_Size(TOGL_BoundBox * BoundBox, PDNMS pdnmsThis)
/* find coordinates of the bounding box of a structure */
{
    if (BoundBox == NULL || pdnmsThis == NULL)
        return;
    OGL_ClearBoundBox(BoundBox);
    BoundBox->set = FALSE;
    TraverseModels(pdnmsThis, TRAVERSE_ATOM, 0, BoundBox, Cn3D_SizeCB);
}


/* Buttons on the OGL Viewer UI */

static ButtoN OGL_allButton;
static ButtoN OGL_rewindButton;
static ButtoN OGL_prevButton;
static ButtoN OGL_nextButton;
#ifndef WIN_MAC
static ButtoN OGL_playButton;
#endif                          /* ndef WIN_MAC */


/* the OpenGL UI */

static void AllLayerOnProc(Nlm_ButtoN b)
{
    OGL_StopPlaying(Cn3D_ColorData.OGL_Data);
    OGL_AllLayerOnProc(Cn3D_ColorData.OGL_Data);
    OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
}

static void RewindLayerProc(Nlm_ButtoN b)
{
    OGL_StopPlaying(Cn3D_ColorData.OGL_Data);
    OGL_RewindLayerProc(Cn3D_ColorData.OGL_Data);
    OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
}

static void PrevLayerProc(Nlm_ButtoN b)
{
    OGL_StopPlaying(Cn3D_ColorData.OGL_Data);
    OGL_PrevLayerProc(Cn3D_ColorData.OGL_Data);
    OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
}

static void NextLayerProc(Nlm_ButtoN b)
{
    OGL_StopPlaying(Cn3D_ColorData.OGL_Data);
    OGL_NextLayerProc(Cn3D_ColorData.OGL_Data);
    OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
}

static void PlayOGL(Nlm_ButtoN b)
{
#ifndef WIN_MAC
    if( OGL_IsPlaying(Cn3D_ColorData.OGL_Data)) {
        OGL_StopPlaying(Cn3D_ColorData.OGL_Data);
        return;
    }
    OGL_StartPlaying(Cn3D_ColorData.OGL_Data);
#endif
}

static void Cn3D_AnimateCloseProc(ButtoN b)
{
    WindoW hAnimateDlg;

	hAnimateDlg= (WindoW)ParentWindow(b);
	if(hAnimateDlg == NULL) return;	

    OGL_StopPlaying(Cn3D_ColorData.OGL_Data);
    Remove(Cn3D_ColorData.AnimateDlg.Cn3D_wAnimate);
    Cn3D_ColorData.AnimateDlg.Cn3D_wAnimate = NULL;
    return;
}


NLM_EXTERN void Cn3D_Animate(IteM i)
{
    GrouP h;
   
    if(Cn3D_ColorData.AnimateDlg.Cn3D_wAnimate != NULL) return;
    Cn3D_ColorData.AnimateDlg.Cn3D_wAnimate =
        FixedWindow(-30, -20, -10, -10, "Model and Animation Controls",
                    NULL);
    OGL_StopPlaying(Cn3D_ColorData.OGL_Data);
    h = HiddenGroup(Cn3D_ColorData.AnimateDlg.Cn3D_wAnimate, 
#ifndef WIN_MAC
        6,
#else
        7,
#endif
        1, NULL);


    OGL_rewindButton = Nlm_PushButton(h, "|<<", RewindLayerProc);
    OGL_prevButton = Nlm_PushButton(h, "|< ", PrevLayerProc);
#ifndef WIN_MAC
    OGL_playButton = Nlm_PushButton(h, " > ", PlayOGL);
#endif
    OGL_nextButton = Nlm_PushButton(h, ">| ", NextLayerProc);
    OGL_allButton = Nlm_PushButton(h, "All", AllLayerOnProc);
    PushButton(h, "Close", (BtnActnProc) Cn3D_AnimateCloseProc);

    Show(Cn3D_ColorData.AnimateDlg.Cn3D_wAnimate);
    return;
}


#endif                          /* _OPENGL */
