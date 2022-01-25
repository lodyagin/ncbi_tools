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

/* constructor for TCn3D_Color structure */
void Cn3D_ConstructColor(TCn3D_Color * Color, Nlm_Uint1 Red, Nlm_Uint1 Green, Nlm_Uint1 Blue)
{
    if(Color == NULL) return;
    Color->Name = NULL;
    Color->Paints = NULL;
    Color->Index = 0;
    Cn3D_SetColorCell(&(Color->ColorCell), Red, Green, Blue);
}


void Cn3D_CopyColorCell(ResidueColorCell * Destination,  ResidueColorCell * Source)
{
    if(Destination == NULL || Source == NULL) return;
    MemCpy(Destination, Source, sizeof(ResidueColorCell));
}


void Cn3D_SetColorCell(ResidueColorCell * ColorCell, Nlm_Uint1 Red, Nlm_Uint1 Green, Nlm_Uint1 Blue)
{
    if(ColorCell == NULL) return;
    ColorCell->rgb[0] = Red;
    ColorCell->rgb[1] = Green;
    ColorCell->rgb[2] = Blue;
}


/* destructor for Color structure */
void Cn3D_DestructColor(TCn3D_Color * Color)
{
    if (Color == NULL) return;
    MemFree(Color->Name);
    MemFree(Color->Paints);
}


void Cn3D_DefaultSSColor(TCn3D_ColorData * ColorData)
{
    if(ColorData == NULL) return;
    Cn3D_SetColorCell(&(ColorData->SSColors[CN3D_COLOR_HELIX].ColorCell), 0, 255, 0);
    Cn3D_SetColorCell(&(ColorData->SSColors[CN3D_COLOR_STRAND].ColorCell), 255, 165, 0);
    Cn3D_SetColorCell(&(ColorData->SSColors[CN3D_COLOR_TURN].ColorCell), 255, 69, 0);
    Cn3D_SetColorCell(&(ColorData->SSColors[CN3D_COLOR_COIL].ColorCell), 0, 255, 255);
}


/* constructor for TCn3D_Color structure */
void Cn3D_ConstructColorData(TCn3D_ColorData * ColorData
#ifdef _OPENGL
                             , TOGL_Data * OGL_Data
#endif
                             )
{
    if (ColorData == NULL) return;
#ifdef _OPENGL
    if(OGL_Data == NULL) return;
#endif
    Cn3D_ConstructColor(&(ColorData->SSColors[CN3D_COLOR_HELIX]), 0, 255, 0);
    ColorData->SSColors[CN3D_COLOR_HELIX].Name = StringSave("Helix");
    Cn3D_ConstructColor(&(ColorData->SSColors[CN3D_COLOR_STRAND]), 255, 165, 0);
    ColorData->SSColors[CN3D_COLOR_STRAND].Name = StringSave("Strand");
    Cn3D_ConstructColor(&(ColorData->SSColors[CN3D_COLOR_TURN]), 255, 69, 0);
    ColorData->SSColors[CN3D_COLOR_TURN].Name = StringSave("Turn");
    Cn3D_ConstructColor(&(ColorData->SSColors[CN3D_COLOR_COIL]), 0, 255, 255);
    ColorData->SSColors[CN3D_COLOR_COIL].Name = StringSave("Coil");
    Cn3D_ConstructColor(&(ColorData->Highlight), 255, 255, 0);
    ColorData->SSColors[CN3D_COLOR_COIL].Name = StringSave("Highlight");
    
    ColorData->Palette = NULL;
#ifdef _OPENGL
    ColorData->OGL_Data = OGL_Data;
#endif
}

/* destructor for Color structure */
void Cn3D_DestructColorData(TCn3D_ColorData * ColorData)
{
    Nlm_Int4 i;
    
    if(ColorData == NULL) return;
    for (i=0;i<CN3D_COLOR_SS; i++) 
        Cn3D_DestructColor(&(ColorData->SSColors[i]));
    Cn3D_FreePalette(&(ColorData->Palette));
}


void Cn3D_FreePalette(ValNodePtr * Palette)
/* free a palette */
{
    ValNodePtr Index;
    
    if(Palette == NULL) return;
    Index = *Palette;
    while (Index) {
        Cn3D_DestructColor((TCn3D_Color *)(Index->data.ptrvalue));
        Index = Index->next;
    }
    ValNodeFree(*Palette);
    *Palette = NULL;
}


void Cn3D_RequestColor(TCn3D_ColorData * ColorData, ResidueColorCell * ColorCell)
/* ask for a new color to be put on the palette */
{
    TCn3D_Color * NewColor;
    
    if(ColorData == NULL || ColorCell == NULL) return;
    if ( Cn3D_SearchColor(ColorData->Palette, ColorCell) == NULL ) {
        NewColor = MemNew(sizeof(TCn3D_Color));
        if (!NewColor) return;
        Cn3D_ConstructColor(NewColor, ColorCell->rgb[0], ColorCell->rgb[1], ColorCell->rgb[2]);
        ValNodeAddPointer(&(ColorData->Palette), 0, NewColor);
    }
    return;
}


ValNodePtr Cn3D_SearchColor(ValNodePtr Palette, ResidueColorCell * ColorCell)
/* search for a color in the palette */
{
    if(ColorCell == NULL) return NULL;
    while(Palette) {
        if(((TCn3D_Color *)(Palette->data.ptrvalue))->ColorCell.rgb[0] == ColorCell->rgb[0] &&
            ((TCn3D_Color *)(Palette->data.ptrvalue))->ColorCell.rgb[1] == ColorCell->rgb[1]  &&
            ((TCn3D_Color *)(Palette->data.ptrvalue))->ColorCell.rgb[2] == ColorCell->rgb[2] )
            return Palette;
        Palette = Palette->next;
    }
    return NULL;
}


void Cn3D_SetColorChoice(ValNodePtr Palette)
/* sets the valnode choice value according to position in palette */
{
    Nlm_Int4 i = 0;
    
    while(Palette) {
        ((TCn3D_Color *)(Palette->data.ptrvalue))->Index = i;
        i++;
        Palette = Palette->next;
    }
}


Nlm_Int4 Cn3D_ColorIndex(TCn3D_ColorData * ColorData, ResidueColorCell * ColorCell)
/* returns the index into the palette for a given color */
{
    ValNodePtr Color;
    
    if(ColorData == NULL || ColorCell == NULL) return 0;
    Color = Cn3D_SearchColor(ColorData->Palette, ColorCell);
    if ( Color == NULL ) return 0;
    else return ((TCn3D_Color *)(Color->data.ptrvalue))->Index;
}


/************************* OpenGL specific functions *****************************/

#ifdef _OPENGL


static void LIBCALLBACK Cn3D_SizeCB(PFB pfbThis, Nlm_Int4 iModel, Nlm_Int4 iIndex, TOGL_BoundBox * BoundBox)
/* callback used to find the bounding box of the atoms of a structure */
{
    PMAD pmadThis = NULL;
    PALD paldThis = NULL;
    Nlm_FloatLoPtr pflvDataThis = NULL;
    Nlm_Int4 i;
    
    if(pfbThis == NULL || BoundBox == NULL) return;
    if (pfbThis->bMe == (Byte) AM_MAD)  /* is this molecular atom data? */
    {  /* iIndex isn't used */
        pmadThis = (PMAD) pfbThis;
        paldThis = GetAtomLocs(pmadThis, iModel);
        while (paldThis)
        {  
            pflvDataThis = paldThis->pflvData;
            if(BoundBox->set == FALSE ) {  /* if the bounding box has not been used, set it */
                for (i=0; i<6; i++ ) {
                    BoundBox->x[i] = pflvDataThis[i/(int)2];
                    BoundBox->set = TRUE;
                }
            }
            else {
                for (i=0; i<6; i += 2 ) {  /* check to see if the atom lies outside the bound box */
                    if( pflvDataThis[i/(Nlm_Int4)2] < BoundBox->x[i] ) BoundBox->x[i] = pflvDataThis[i/(Nlm_Int4)2];
                    if( pflvDataThis[i/(Nlm_Int4)2] > BoundBox->x[i+1] ) BoundBox->x[i+1] = pflvDataThis[i/(Nlm_Int4)2];
                }
            }
            paldThis = paldThis->next;
        }
    }
}




void LIBCALL Cn3D_Size(TOGL_BoundBox * BoundBox, PDNMS pdnmsThis)
/* find coordinates of the bounding box of a structure */
{
    if(BoundBox == NULL || pdnmsThis == NULL) return;
    OGL_ClearBoundBox(BoundBox);
    BoundBox->set = FALSE;
    TraverseModels(pdnmsThis, TRAVERSE_ATOM, 0, BoundBox,  Cn3D_SizeCB);
}


/* Buttons on the OGL Viewer UI */

static ButtoN  OGL_allButton;
static ButtoN  OGL_rewindButton;
static ButtoN  OGL_prevButton;
static ButtoN  OGL_nextButton;
#ifndef WIN_MAC
static ButtoN  OGL_playButton;
#endif /* ndef WIN_MAC */


/* the OpenGL UI */

static void AllLayerOnProc(Nlm_ButtoN b)
{   
    OGL_AllLayerOnProc((TOGL_Data *)Nlm_GetObjectExtra( b ));
    OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
}

static void RewindLayerProc(Nlm_ButtoN b)
{   
    OGL_RewindLayerProc((TOGL_Data *)Nlm_GetObjectExtra( b ));
    OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
}

static void PrevLayerProc(Nlm_ButtoN b)
{   
    OGL_PrevLayerProc((TOGL_Data *)Nlm_GetObjectExtra( b ));
    OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
}

static void NextLayerProc(Nlm_ButtoN b)
{   
    OGL_NextLayerProc((TOGL_Data *)Nlm_GetObjectExtra( b ));
    OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
}


extern Nlm_Boolean OGL_IsPlaying(void)
{
#ifdef WIN_MAC
    return FALSE;
#else
    return (Nlm_Boolean)(OGL_playButton  &&  Nlm_GetStatus(OGL_playButton));
#endif
}

extern void OGL_StopPlaying(void) 
{
#ifndef WIN_MAC
    if ( OGL_playButton ) Nlm_SetStatus(OGL_playButton, FALSE);
#endif
}

static void PlayOGL(Nlm_ButtoN b)
{
#ifndef WIN_MAC
    TOGL_Data * OGL_Data;
    
    OGL_Data = (TOGL_Data *)Nlm_GetObjectExtra( b );
    
    while (!Nlm_QuittingProgram()  &&  OGL_IsPlaying())
    {
        OGL_Play(OGL_Data);
        OGL_DrawViewer3D(Cn3D_ColorData.OGL_Data);
        
        while ( Nlm_EventAvail() )
            Nlm_ProcessAnEvent();
    }
#endif
}


Nlm_GrouP LIBCALL  OGL_Controls ( Nlm_GrouP prnt)
{
    Nlm_GrouP g, h;
    
    g = Nlm_HiddenGroup ( prnt, -1, 0, NULL );
    if (!g) return NULL;
    
    h = Nlm_HiddenGroup (g, 1, 5, NULL);
    
    /* set up the procedures to run */
    
    OGL_allButton    = Nlm_PushButton(h, "all",  AllLayerOnProc);
    OGL_rewindButton = Nlm_PushButton(h, "<<", RewindLayerProc);
    OGL_prevButton   = Nlm_PushButton(h, "<",  PrevLayerProc);
    OGL_nextButton   = Nlm_PushButton(h, ">",  NextLayerProc);
#ifndef WIN_MAC
    OGL_playButton =  Nlm_CheckBox(h, "Animate", PlayOGL);
#endif
    
    /* add pointers to global info */
    
    Nlm_SetObjectExtra(OGL_allButton,    (Nlm_VoidPtr) (Cn3D_ColorData.OGL_Data), NULL);
    Nlm_SetObjectExtra(OGL_rewindButton, (Nlm_VoidPtr) (Cn3D_ColorData.OGL_Data), NULL);
    Nlm_SetObjectExtra(OGL_prevButton,   (Nlm_VoidPtr) (Cn3D_ColorData.OGL_Data), NULL);
    Nlm_SetObjectExtra(OGL_nextButton,   (Nlm_VoidPtr) (Cn3D_ColorData.OGL_Data), NULL);
#ifndef WIN_MAC
    Nlm_SetObjectExtra(OGL_playButton,   (Nlm_VoidPtr) (Cn3D_ColorData.OGL_Data), NULL);
#endif    
    
    return g;
}


#endif /* _OPENGL */
