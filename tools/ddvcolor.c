/*   $Id: ddvcolor.c,v 6.4 1999/07/20 14:38:11 durand Exp $
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
* File Name:  $Id: ddvcolor.c,v 6.4 1999/07/20 14:38:11 durand Exp $
*
* Author:  Lewis Geer
*
* Version Creation Date:   6/4/99
*
* $Revision: 6.4 $
*
* File Description: Shared color information for viewers
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ddvcolor.c,v $
* Revision 6.4  1999/07/20 14:38:11  durand
* update DDV_InRange to include the ends in the test
*
* Revision 6.3  1999/07/19 19:12:21  lewisg
* fix bug adding new pColor in old pMediaInfo.  Also added ability to detect overlaps when allocating new pColor
*
* Revision 6.2  1999/07/19 14:37:03  lewisg
* fix bug when allocating mediainfo with from != 0
*
* Revision 6.1  1999/07/16 18:46:46  lewisg
* moved ddvcolor from api to tools
*
* Revision 1.4  1999/07/13 23:24:48  lewisg
* added DDV_SipList
*
* Revision 1.3  1999/07/13 14:38:39  lewisg
* separated out networking code
*
* Revision 1.4  1999/06/24 17:48:28  lewisg
* added SpecialColors to global
*
* Revision 1.3  1999/06/16 13:51:50  lewisg
* added palette management functions from cn3d
*
* Revision 1.2  1999/06/14 17:20:47  lewisg
* minor bug fix for leaks
*
* Revision 1.1  1999/06/11 23:37:09  lewisg
* color management functions
*
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <sequtil.h>
#include <ddvcolor.h>
#include <salutil.h>



/*****************************************************************************
*
*   Set default values for DDV_Position.  
*
*****************************************************************************/

Int4 DDV_DefaultPosition(DDV_Position *pPosition, Int4 lPosition, 
                         Boolean fAlign)
{
    if(pPosition == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DefaultPosition");
        return 0;
    }

    pPosition->lPosition = lPosition;
    pPosition->fAlign = fAlign;
    return 1;
}

/*****************************************************************************
*
*   Set default values for DDV_Range.  
*
*****************************************************************************/

Int4 DDV_DefaultRange(DDV_Range *pRange, Int4 lFrom, Int4 lTo, 
                         Boolean fAlign)
{
    if(pRange == NULL || lFrom > lTo ) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DefaultRange");
        return 0;
    }
    
    pRange->lFrom = lFrom;
    pRange->lTo = lTo;
    pRange->fAlign = fAlign;
    return 1;
}

/*****************************************************************************
*
*   Constructor for DDV_ColorEntry structure.  DDV_ColorEntry is a
*   (Red, Green, Blue) color triplet and associated szName.
*
*****************************************************************************/

DDV_ColorEntry * DDV_CreateColorEntry(Char *szName, Nlm_Uint1 Red,
                          Nlm_Uint1 Green, Nlm_Uint1 Blue)
{
    DDV_ColorEntry *pColor;
    
    pColor = MemNew(sizeof(DDV_ColorEntry));
    if (pColor == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid malloc in DDV_CreateColorEntry");
        return NULL;
    }

    pColor->Name = StringSave(szName);
    pColor->Paints = NULL;
    pColor->Index = 0;
    DDV_SetColorCell(&pColor->ColorCell, Red, Green, Blue);

    return pColor;
}


/*****************************************************************************
*
*   Destructor for DDV_ColorEntry Color.
*
*****************************************************************************/

/* destructor for Color structure */
void DDV_DeleteColorEntry(DDV_ColorEntry *Color)
{
    MemFree(Color->Name);
    MemFree(Color->Paints);
    MemFree(Color);
}


/*****************************************************************************
*
*   Frees a Palette.  Palette is a ValNode list of DDV_ColorEntries
*   Note that it frees the entire Valnode list pointed to by Palette and
*   sets Palette to NULL.
*
*****************************************************************************/

void DDV_FreePalette(ValNode **Palette)
{
    ValNode *Index;
    
    if(Palette == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_FreePalette");
        return;
    }

    Index = *Palette;
    while (Index) {
        DDV_DeleteColorEntry((DDV_ColorEntry *)(Index->data.ptrvalue));
        Index = Index->next;
    }
    ValNodeFree(*Palette);
    *Palette = NULL;
}

/*****************************************************************************
*
*   Puts ColorCell on the Palette list if it isn't there already.  
*
*****************************************************************************/

void DDV_RequestColor(ValNode **Palette, DDV_ColorCell *ColorCell)
{
    DDV_ColorEntry *NewColor;
    
    if(Palette == NULL || ColorCell == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_RequestColor");
        return;
    }

    if ( DDV_SearchColor(*Palette, ColorCell) == NULL ) {
        NewColor = DDV_CreateColorEntry("", ColorCell->rgb[0], ColorCell->rgb[1],
            ColorCell->rgb[2]);
        if (!NewColor) return;
        ValNodeAddPointer(Palette, 0, NewColor);
    }
    return;
}

/*****************************************************************************
*
*   Puts ColorEntry on the Palette list if the sZname isn't there already,
*   Otherwise replaces the entry.  Does NOT make a copy of the DDV_ColorEntry
*
*****************************************************************************/

void DDV_RequestColorbyName(ValNode **Palette, DDV_ColorEntry *pColorIn)
{
    DDV_ColorEntry *pColorEntry;
    ValNode *pvn;
    
    if(Palette == NULL || pColorIn == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_RequestColorbyName");
        return;
    }

    pColorEntry = DDV_SearchColorbyName(*Palette, pColorIn->Name);
    if (pColorEntry == NULL) {
        ValNodeAddPointer(Palette, 0, pColorIn);
        return;
    }
    else {
        for(pvn = *Palette; pvn->data.ptrvalue != pColorEntry;
            pvn = pvn->next);
        DDV_DeleteColorEntry((DDV_ColorEntry*)pvn->data.ptrvalue);
        pvn->data.ptrvalue = pColorIn;
    }
}


/*****************************************************************************
*
*   Looks for ColorCell on the Palette.  Returns * to the 
*   DDV_ColorEntry.
*
*****************************************************************************/

DDV_ColorEntry * DDV_SearchColor(ValNode *Palette, DDV_ColorCell *ColorCell)
{
    ValNode *pvnPalette = Palette;

    if(ColorCell == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SearchColor");
        return NULL;
    }
 
    while(pvnPalette) {
        if(((DDV_ColorEntry *)(pvnPalette->data.ptrvalue))->ColorCell.rgb[0]
            == ColorCell->rgb[0] &&
            ((DDV_ColorEntry *)(pvnPalette->data.ptrvalue))->ColorCell.rgb[1]
            == ColorCell->rgb[1]  &&
            ((DDV_ColorEntry *)(pvnPalette->data.ptrvalue))->ColorCell.rgb[2]
            == ColorCell->rgb[2] )
            return (DDV_ColorEntry *)(pvnPalette->data.ptrvalue);
        pvnPalette = pvnPalette->next;
    }
    return NULL;
}

/*****************************************************************************
*
*   Looks for a DDV_ColorEntry on the Palette by szName.  Returns *
*   to the DDV_ColorEntry.
*
*****************************************************************************/

DDV_ColorEntry * DDV_SearchColorbyName(ValNode *Palette, Char *szName)
{
    ValNode *pvnPalette = Palette;

    while(pvnPalette) {
        if(StrCmp(((DDV_ColorEntry *)(pvnPalette->data.ptrvalue))->Name,
            szName) == 0 )
            return (DDV_ColorEntry *)(pvnPalette->data.ptrvalue);
        pvnPalette = pvnPalette->next;
    }
    return NULL;
}

/*****************************************************************************
*
*   Looks for a DDV_ColorEntry on the Palette by szName.  Returns *
*   to ColorCell in the DDV_ColorEntry.
*
*****************************************************************************/

DDV_ColorCell * DDV_SearchColorCellbyName(ValNode *Palette, Char *szName)
{
    DDV_ColorEntry *pColorEntry;

    pColorEntry = DDV_SearchColorbyName(Palette, szName);
    if(pColorEntry == NULL) return NULL;

    return &pColorEntry->ColorCell;
}

/*****************************************************************************
*
*   Sets the Index in a DDV_ColorEntry according to position in Palette
*
*****************************************************************************/

void DDV_SetColorChoice(ValNode *Palette)
/* sets the valnode choice value according to position in palette */
{
    Int4 i = 0;
    ValNode *pvnPalette = Palette;

    if(Palette == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SetColorChoice");
        return;
    }

    while(pvnPalette) {
        ((DDV_ColorEntry *)(pvnPalette->data.ptrvalue))->Index = i;
        i++;
        pvnPalette = pvnPalette->next;
    }
}

/*****************************************************************************
*
*   Returns the index into the Palette for a given ColorCell.
*   Return index on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_ColorIndex(ValNode *Palette, DDV_ColorCell *ColorCell)
/* returns the index into the palette for a given color */
{
    DDV_ColorEntry *Color;
    
    if(Palette == NULL || ColorCell == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_ColorIndex");
        return 0;
    }

    Color = DDV_SearchColor(Palette, ColorCell);
    if ( Color == NULL ) return 0;
    else return Color->Index;
}


/*****************************************************************************
*
*   Copy a DDV_ColorCell.
*   Return 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_CopyColorCell(DDV_ColorCell *pDestination,  DDV_ColorCell *pSource)
{
    if(pDestination == NULL || pSource == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_CopyColorCell");
        return 0;
    }

    MemCpy(pDestination, pSource, sizeof(DDV_ColorCell));
    return 1;
}

/*****************************************************************************
*
*   Set a color cell using Uint1 values
*   Return 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_SetColorCell(DDV_ColorCell *pColorCell, Nlm_Uint1 Red,
                      Nlm_Uint1 Green, Nlm_Uint1 Blue)
{
    if(pColorCell == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SetColorCell");
        return 0;
    }

    pColorCell->rgb[0] = Red;
    pColorCell->rgb[1] = Green;
    pColorCell->rgb[2] = Blue;
    return 1;
}


/*****************************************************************************
*
*   Create and maintain the MediaInfo2 structure which contain coloring
*   information, visible/invisible, aligned information for a given object.
*
*
*****************************************************************************/

/*****************************************************************************
*
*   Checks to see if a pPosition is inside pRange.
*   If it is, return 1
*   If it isn't, return 0
*   If pPosition and pRange don't have the same coordinate system, return -1
*
*****************************************************************************/

Int4 DDV_InRange(DDV_Position *pPosition, DDV_Range *pRange)
{
    if(pPosition == NULL || pRange == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_InRange");
        return -1;
    }

    if(pPosition->fAlign != pRange->fAlign) return -1;
    if(pPosition->lPosition >= pRange->lFrom && pPosition->lPosition <= pRange->lTo)
        return 1;
    return 0;
}

/*****************************************************************************
*
*   Checks to see if pRange1 overlaps pRange2.
*   If it does completely, return DDV_TOTALLAP
*   If it doesn't, return DDV_NOLAP
*   If pRange1 overlaps the front of pRange2, return DDV_FRONTLAP
*   If pRange1 overlaps the rear of pRange2, return DDV_BACKLAP
*   If pRange1 and pRange2 don't have the same coordinate system, return -1
*
*****************************************************************************/

Int4 DDV_RangeOverlap(DDV_Range *pRange1, DDV_Range *pRange2)
{
    if(pRange1 == NULL || pRange2 == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_RangeOverlap");
        return -1;
    }

    if(pRange1->fAlign != pRange2->fAlign) return -1;

    if(pRange1->lFrom >= pRange2->lFrom) {
        if(pRange1->lTo <= pRange2->lTo) return DDV_TOTALLAP;
        else return DDV_BACKLAP;
    }
    else if( pRange1->lTo >= pRange2->lFrom) return DDV_FRONTLAP;

    return DDV_NOLAP;
}

/*****************************************************************************
*
*   Returns a DDV_ColorCell for a given position if it is inside of pColor.
*   Returns NULL otherwise.
*
*****************************************************************************/

DDV_ColorCell * DDV_GetCellByPosition(DDV_Color *pColor,
                                      DDV_Position *pPosition)
{
    if(pColor == NULL || pPosition == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_GetCellByPosition");
        return NULL;
    }
    
    if(DDV_InRange(pPosition, &pColor->Range) != 1) return NULL;
    return &(pColor->pColorCell[pPosition->lPosition - pColor->Range.lFrom]);
}

/*****************************************************************************
*
*   Sets a DDV_ColorCell for a given position if it is inside of pColor.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetCellByPosition(DDV_Color * pColor, DDV_Position *pPosition,
                 DDV_ColorCell *pColorCell)
{
    if(pColor == NULL || pColorCell == NULL || pPosition == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SetCellByPosition");
        return 0;
    }
    
    if(DDV_InRange(pPosition, &pColor->Range) != 1) return 0;
    
    MemCpy(&pColor->pColorCell[pPosition->lPosition - pColor->Range.lFrom],
        pColorCell, sizeof(DDV_ColorCell));

    return 1;
}

/*****************************************************************************
*
*   Creates a DDV_Color that spans the entire object given by sip with
*   coordinates lFrom to lTo and using fAlign coordinates.
*   The default color is black.
*   Returns NULL on error.
*
*****************************************************************************/

DDV_Color * DDV_CreateDefaultColor(DDV_ColorGlobal *pColorGlobal,
                                        SeqId *sip, Int4 lFrom, Int4 lTo,
                                        Boolean fAlign)
{
    DDV_Color *pColor;
    Int4 i;
    DDV_ColorCell *pColorCell;

    if (sip == NULL || pColorGlobal == NULL || lFrom > lTo) {
        ErrPostEx(SEV_ERROR, 0, 0, 
            "Invalid call on DDV_CreateDefaultColor");
        return NULL;
    }

    pColor = MemNew(sizeof(DDV_Color));
    if (pColor == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, 
            "Unable to malloc DDV_Color in DDV_CreateDefaultColor");
        return NULL;
    }

    pColor->Range.lFrom = lFrom; 
    pColor->Range.lTo = lTo;
    pColor->Range.fAlign = fAlign;
    pColor->pColorCell = MemNew(sizeof(DDV_ColorCell) * (lTo - lFrom + 1));
    if( pColor->pColorCell == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, 
            "Unable to malloc DDV_ColorCells in DDV_CreateDefaultColor");
        return NULL;
    }

    pColorCell = DDV_SearchColorCellbyName(pColorGlobal->pvnSpecialColors,
            "Default");
    if(pColorCell == NULL) return NULL;

    for(i = 0; i <= lTo - lFrom; i++) DDV_CopyColorCell(
        &pColor->pColorCell[i], pColorCell);
   
    return pColor;
}

/*****************************************************************************
*
*   Deletes a DDV_Color.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_DeleteColor(DDV_Color * pColor)
{
    if(pColor == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DeleteColor");
        return 0;
    }

    MemFree(pColor->pColorCell);
    MemFree(pColor);
    return 1;
}

/* todo: delete mediainfo for a given sip */
/* todo: lock count for mediainfo */
/* todo: lock count for colorglobal */


/*****************************************************************************
*
*   Creates a new DDV_MediaInfo structure.  Creates a default color structure
*   for the given object specified by sip.
*   Inserts the given pColor and sets fVisible.  Duplicates the sip.
*
*   Returns NULL on error.
*
*****************************************************************************/

DDV_MediaInfo * DDV_NewMediaInfo(SeqId *sip, DDV_Color *pColor,
                                  Boolean fVisible)
{
    DDV_MediaInfo *pMediaInfo;

    if (sip == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_NewMediaInfo");
        return NULL;
    }


    pMediaInfo = MemNew(sizeof(DDV_MediaInfo));
    if (pMediaInfo == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid Malloc in DDV_NewMediaInfo");
        return NULL;
    }

    pMediaInfo->sip = SeqIdDupList(sip);
    pMediaInfo->pvnColor = NULL;
    if(pColor != NULL) ValNodeAddPointer(&pMediaInfo->pvnColor, 0, pColor);
    pMediaInfo->fVisible = fVisible;
    
    return pMediaInfo;
}

/*****************************************************************************
*
*   Returns a ValNode list, each element of which points to the SeqId of a
*   sequence that has a corresponding DDV_MediaInfo structure in pColorGlobal.
*
*   The Valnode list should be freed by the calling routine.  However, the 
*   SeqIds pointed to by the list should *not* be freed.
*
*   Returns NULL on error.
*
*****************************************************************************/

ValNode * DDV_SipList(DDV_ColorGlobal *pColorGlobal)
{
    ValNode *pvn, *pvnReturn = NULL;
    DDV_MediaInfo *pMediaInfo;

    if(pColorGlobal == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SipList");
        return NULL;
    }

    for(pvn = pColorGlobal->pvnMediaInfo ; pvn != NULL; pvn = pvn->next) {
        pMediaInfo = (DDV_MediaInfo *)pvn->data.ptrvalue;
        ValNodeAddPointer(&pvnReturn, 0, pMediaInfo->sip);
    }

    return pvnReturn;
}

/*****************************************************************************
*
*   Creates a default MediaInfo structure and adds it to pColorGlobal
*   for the given object specified by sip. Duplicates the sip.
*   The sequence is visible by default.
*   Coordinates lFrom to lTo and using fAlign coordinates.
*
*   Returns 0 on error.
*
*****************************************************************************/

Int4 DDV_DefaultMediaInfoByLen(DDV_ColorGlobal * pColorGlobal, SeqId *sip,
                               Int4 lFrom, Int4 lTo, Boolean fAlign)
{
    DDV_MediaInfo *pMediaInfo;
    DDV_Color *pColor;
    ValNode *pvn;
    Int4 lOverlap;
    DDV_Range Range;

    if(pColorGlobal == NULL || sip == NULL || lFrom > lTo) {
        ErrPostEx(SEV_ERROR, 0, 0, 
            "Invalid call on DDV_DefaultMediaInfoByLen");
        return 0;
    }

    Range.fAlign = fAlign;
    Range.lFrom = lFrom;
    Range.lTo = lTo;

    /* check to see if it already exits */
    pMediaInfo = DDV_GetMediaInfo(pColorGlobal, sip);
    if (pMediaInfo == NULL) {
        pMediaInfo = DDV_NewMediaInfo(sip, NULL, TRUE);
        if( pMediaInfo == NULL) {
            ErrPostEx(SEV_ERROR, 0, 0, 
                "Unable to create pMediaInfo in DDV_DefaultMediaInfoByLen");
            return 0;
        }
        DDV_AddMediaInfo(pColorGlobal, pMediaInfo);
    }
    
    for(pvn = pMediaInfo->pvnColor; pvn != NULL; pvn = pvn->next) {
        pColor = (DDV_Color *)pvn->data.ptrvalue;
        lOverlap = DDV_RangeOverlap(&Range, &pColor->Range);
        if(lOverlap == DDV_NOLAP) continue;
        else if(lOverlap == DDV_TOTALLAP ) return 1;
        else if(lOverlap == DDV_FRONTLAP) {
            lTo = pColor->Range.lFrom - 1;
            break;
        }
        else if(lOverlap == DDV_BACKLAP) {
            lFrom = pColor->Range.lTo + 1;
            break;
        }
        else return 0;
    }
    
    pColor = DDV_CreateDefaultColor(pColorGlobal, sip, lFrom, lTo, fAlign);
    if(pColor == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, 
            "Unable to create pColor in DDV_DefaultMediaInfoByLen");
        return 0;
    }
    ValNodeAddPointer(&pMediaInfo->pvnColor, 0, pColor);

    return 1;
}

/* simplified version of above that retrieves the length and sets the coordinate
   system to sequence */

Int4 DDV_DefaultMediaInfo(DDV_ColorGlobal * pColorGlobal, SeqId *sip)
{
    Int4 lLength;
    Bioseq *bsp;

    if(pColorGlobal == NULL || sip == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DefaultMediaInfo");
        return 0;
    }

    bsp = BioseqLockById(sip);
    if(bsp == NULL) {
       ErrPostEx (SEV_ERROR, 0, 0, 
           "BioseqLockById failed in DDV_DefaultMediaInfo");
       return 0;
    }

    lLength = BioseqGetLen(bsp) - 1;
    BioseqUnlock (bsp);

    return DDV_DefaultMediaInfoByLen(pColorGlobal, sip, 0, lLength, FALSE);
}

/*****************************************************************************
*
*   Deletes a DDV_MediaInfo structure.
*
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_DeleteMediaInfo(DDV_MediaInfo *pMediaInfo)
{
    ValNode *pvn;

    if(pMediaInfo == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DeleteMediaInfo");
        return 0;
    }

    SeqIdFree(pMediaInfo->sip);
    
    for(pvn = pMediaInfo->pvnColor; pvn != NULL; pvn = pvn->next) 
        DDV_DeleteColor((DDV_Color *)pvn->data.ptrvalue);
    ValNodeFree(pMediaInfo->pvnColor);

    MemFree(pMediaInfo);
    return 1;
}




/*****************************************************************************
*
*   The color function queue.  The queue is used to 
*   register color function by szName (like "Color by Conservation") and to 
*   give each function an lPriority.  The algorithm to determine what function
*   to call is:
*   - if requested, a default coloration is performed.
*   - look for an fOverride algorithm with the correct szName. If found, run it,
*     then stop.
*   - look for all function that match the correct szName OR have fCallMe set.
*     Run each algorithm in priority order.
*   
*   Viewers must take care to set the visible bits beforehand as this
*   influences coloration.
*
*   pfnDDV_ColorFunc takes as arguments pData for user data and a range
*   lFrom to lTo.  If lTo < lFrom, the range should be ignored.
*
*****************************************************************************/


/*****************************************************************************
*
*   Returns a new color queue entry with name szName, priority lPriority.
*   fOverride is explained above.
*
*****************************************************************************/

DDV_ColorQueue * DDV_NewColorQueue(DDV_ColorFunc pfnColorFunc, Char * szName,
                              Int4 lPriority, Boolean fOverride)
{
    DDV_ColorQueue *pColorQueue;

    if(pfnColorFunc == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_NewColorFn");
        return NULL;
    }

    pColorQueue = MemNew(sizeof(DDV_ColorQueue));
    if(pColorQueue == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Unable to Malloc in DDV_NewColorFn");
        return NULL;
    }

    pColorQueue->fOverride = fOverride;
    pColorQueue->lPriority = lPriority;
    pColorQueue->pfnColorFunc = pfnColorFunc;
    pColorQueue->szName = StringSave(szName);
    pColorQueue->next = NULL;

    return pColorQueue;
}

/*****************************************************************************
*
*   Deletes a DDV_ColorQueue.  
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_DeleteColorQueue(DDV_ColorQueue *pColorQueue)
{
    if(pColorQueue == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DeleteColorFn");
        return 0;
    }
    MemFree(pColorQueue->szName);
    MemFree(pColorQueue);
    return 1;
}



/*****************************************************************************
*
*   Sets up default secondary structure colors in 
*   pColorGlobal->pvnSpecialColors palette.
*   Returns 0 on failure.
*
*****************************************************************************/

Int4 DDV_DefaultSSColor(DDV_ColorGlobal * pColorGlobal)
{
    DDV_ColorEntry *pColorEntry;

    if(pColorGlobal == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid Call on DDV_DefaultSSColor");
        return 0;
    }

    pColorEntry = DDV_CreateColorEntry("Helix", 0, 255, 0);
    ValNodeAddPointer(&pColorGlobal->pvnSpecialColors, 0, pColorEntry);
    pColorEntry = DDV_CreateColorEntry("Strand", 255, 165, 0);
    ValNodeAddPointer(&pColorGlobal->pvnSpecialColors, 0, pColorEntry);
    pColorEntry = DDV_CreateColorEntry("Turn", 255, 69, 0);
    ValNodeAddPointer(&pColorGlobal->pvnSpecialColors, 0, pColorEntry);
    pColorEntry = DDV_CreateColorEntry("Coil", 0, 255, 255);
    ValNodeAddPointer(&pColorGlobal->pvnSpecialColors, 0, pColorEntry);

    return 1;
}

/*****************************************************************************
*
*   Returns a new DDV_ColorGlobal structure.
*   Returns NULL on failure
*
*****************************************************************************/

DDV_ColorGlobal * DDV_CreateColorGlobal(Boolean fDefaultColor)
{
    DDV_ColorGlobal *pColorGlobal;
    DDV_ColorEntry *pColorEntry;

    pColorGlobal = MemNew(sizeof(DDV_ColorGlobal));
    if(pColorGlobal == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid Malloc in DDV_CreateColorGlobal");
        return NULL;
    }

    pColorGlobal->fColorByMaster = FALSE;
    pColorGlobal->fDefaultColor = fDefaultColor;
    pColorGlobal->MasterSeqId = NULL;
    pColorGlobal->pvnColorQueue = NULL;
    pColorGlobal->pvnAllColorFns = NULL;
    pColorGlobal->pvnMediaInfo = NULL;
    pColorGlobal->Palette = NULL;
    pColorGlobal->pvnSpecialColors = NULL;

    DDV_DefaultSSColor(pColorGlobal);
    pColorEntry = DDV_CreateColorEntry("Highlight", 255, 255, 0);
    ValNodeAddPointer(&pColorGlobal->pvnSpecialColors, 0, pColorEntry);
    pColorEntry = DDV_CreateColorEntry("Default", 0, 0, 0);
    ValNodeAddPointer(&pColorGlobal->pvnSpecialColors, 0, pColorEntry);

    return pColorGlobal;
}

/*****************************************************************************
*
*   Deletes a  DDV_ColorGlobal structure.
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_DeleteColorGlobal(DDV_ColorGlobal *pColorGlobal)
{
    ValNode *pvn;

    if(pColorGlobal == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DeleteColorGlobal");
        return 0;
    }

    for(pvn = pColorGlobal->pvnColorQueue; pvn != NULL; pvn = pvn->next) 
        DDV_DeleteColorQueue((DDV_ColorQueue *)pvn->data.ptrvalue);
    ValNodeFree(pColorGlobal->pvnColorQueue);

    for(pvn = pColorGlobal->pvnAllColorFns; pvn != NULL; pvn = pvn->next) 
        DDV_DeleteColorQueue((DDV_ColorQueue *)pvn->data.ptrvalue);
    ValNodeFree(pColorGlobal->pvnAllColorFns);

    for(pvn = pColorGlobal->pvnMediaInfo ; pvn != NULL; pvn = pvn->next) 
        DDV_DeleteMediaInfo((DDV_MediaInfo *)pvn->data.ptrvalue);
    ValNodeFree(pColorGlobal->pvnMediaInfo);

    DDV_FreePalette(&pColorGlobal->Palette);
    DDV_FreePalette(&pColorGlobal->pvnSpecialColors);

    MemFree(pColorGlobal);
    return 1;
}

/*****************************************************************************
*
*   Delete all information for the given SeqId.
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_RemoveMediaInfo(DDV_ColorGlobal *pColorGlobal, SeqId *sip)
{
    ValNode *pvn, *pvnPrevious = NULL;
    DDV_MediaInfo *pMediaInfo;

    if(pColorGlobal == NULL || sip == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_RemoveMediaInfo");
        return 0;
    }

    for(pvn = pColorGlobal->pvnMediaInfo ; pvn != NULL; pvn = pvn->next) {
        pMediaInfo = (DDV_MediaInfo *)pvn->data.ptrvalue;
        if(DDV_SeqIdCompare(sip, pMediaInfo->sip)) {
            DDV_DeleteMediaInfo(pMediaInfo);
            if(pvnPrevious != NULL) {  /* splice out node */
                pvnPrevious->next = pvn->next;
            }
            else {  /* if first node */
                pColorGlobal->pvnMediaInfo = pvn->next;
            }
            MemFree(pvn);
            break;
        }
        pvnPrevious = pvn;
    }

    return 1;
}

/*****************************************************************************
*
*   Add MediaInfo to DDV_ColorGlobal
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_AddMediaInfo(DDV_ColorGlobal *pColorGlobal,
                              DDV_MediaInfo *pMediaInfo)
{
    if(pColorGlobal == NULL || pMediaInfo == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_AddMediaInfo");
        return 0;
    }

    ValNodeAddPointer(&pColorGlobal->pvnMediaInfo, 0, pMediaInfo);
    return 1;
}

/*****************************************************************************
*
*   Compare two SeqId to make sure all valnodes compare exactly
*
*****************************************************************************/

Boolean DDV_SeqIdCompare(SeqId *sip1, SeqId *sip2)
{
    SeqId *sip;
    Boolean retval = TRUE;

    if(sip1 == NULL || sip2 == NULL) return FALSE;
    if(ValNodeLen(sip1) != ValNodeLen(sip2)) return FALSE;

    for(sip = sip1; sip != NULL; sip = sip->next)
        if(!SeqIdIn(sip, sip2)) retval = FALSE;

    return retval;
}

/*****************************************************************************
*
*   Given an SeqId, return the associated DDV_MediaInfo.
*
*****************************************************************************/

DDV_MediaInfo * DDV_GetMediaInfo(DDV_ColorGlobal *pColorGlobal, SeqId * sip)
{
    ValNode *pvn;
    DDV_MediaInfo *pMediaInfo;

    if(pColorGlobal == NULL || sip == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_GetMediaInfo");
        return NULL;
    }
    
    for(pvn = pColorGlobal->pvnMediaInfo ; pvn != NULL; pvn = pvn->next) {
        pMediaInfo = (DDV_MediaInfo *)pvn->data.ptrvalue;
        if(DDV_SeqIdCompare(sip, pMediaInfo->sip)) break;
    }

    if(pvn == NULL) return NULL;
    else return pMediaInfo;
}

/*****************************************************************************
*
*   Given an SeqId, return the visible state of the SeqId.
*
*****************************************************************************/

Boolean DDV_IsVisible(DDV_ColorGlobal *pColorGlobal, SeqId *sip)
{
    DDV_MediaInfo *pMediaInfo;

    if(pColorGlobal == NULL || sip == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_IsVisible");
        return TRUE;
    }

    pMediaInfo = DDV_GetMediaInfo(pColorGlobal, sip);
    if(pMediaInfo == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "No MediaInfo in DDV_IsVisible");
        return TRUE;
    }

    return pMediaInfo->fVisible;
}

/*****************************************************************************
*
*   Set the visible state (fVisible) of a SeqId.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetVisible(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                    Boolean fVisible)
{
    DDV_MediaInfo *pMediaInfo;

    if(pColorGlobal == NULL || sip == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SetVisible");
        return 0;
    }

    pMediaInfo = DDV_GetMediaInfo(pColorGlobal, sip);
    if(pMediaInfo == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "No MediaInfo in DDV_SetVisible");
        return TRUE;
    }

    pMediaInfo->fVisible = fVisible;

    return 1;
}

/*****************************************************************************
*
*   Retrieve a color for an SeqId at position pPosition or at lPosition in
*   the fAlign coordinate system.
*   Returns a DDV_ColorCell containing the color.  NULL on failure.
*
*   The DDV_ColorCell returned is read only.
*
*****************************************************************************/

DDV_ColorCell * DDV_GetColorSimple(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                             Int4 lPosition, Boolean fAlign)
{
    DDV_Position Position;

    Position.lPosition = lPosition;
    Position.fAlign = fAlign;
    return DDV_GetColor(pColorGlobal, sip, &Position);
}

DDV_ColorCell * DDV_GetColor(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                             DDV_Position *pPosition)
{
    DDV_MediaInfo *pMediaInfo;
    ValNode *pvn;
    DDV_Color *pColor;
    DDV_ColorCell *pColorCell = NULL;

    if(pColorGlobal == NULL || sip == NULL || pPosition == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_GetColor");
        return NULL;
    }

    pMediaInfo = DDV_GetMediaInfo(pColorGlobal, sip);
    if(pMediaInfo == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "No MediaInfo in DDV_GetColor");
        return NULL;
    }

    for(pvn = pMediaInfo->pvnColor; pvn != NULL; pvn = pvn->next) {
        pColor = (DDV_Color *)pvn->data.ptrvalue;
        pColorCell = DDV_GetCellByPosition(pColor, pPosition);
        if(pColorCell != NULL) break;
    }

    if(pColorCell == NULL) {
        pColorCell = DDV_SearchColorCellbyName(pColorGlobal->pvnSpecialColors,
            "Default");
    }
    return pColorCell;
    
}

/*****************************************************************************
*
*   Set a color for in a DDV_MediaInfo at position pPosition.  The color set is
*   pColorCell.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetColorInMediaInfo(DDV_ColorGlobal *pColorGlobal,
                             DDV_MediaInfo *pMediaInfo, DDV_Position *pPosition,
                             DDV_ColorCell *pColorCell)
{
    ValNode *pvn;
    DDV_Color *pColor;
    Int4 retval = 0;

    if(pColorGlobal == NULL || pColorCell == NULL || pMediaInfo == NULL ||
        pPosition == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SetColorInMediaInfo");
        return 0;
    }

    for(pvn = pMediaInfo->pvnColor; pvn != NULL; pvn = pvn->next) {
        pColor = (DDV_Color *)pvn->data.ptrvalue;
        retval = DDV_SetCellByPosition(pColor, pPosition, pColorCell);
        if(retval) break;
    }
    if(retval) DDV_RequestColor(&pColorGlobal->Palette, pColorCell);

    return retval;
}

/*****************************************************************************
*
*   Set a color for an SeqId at position pPositionor at lPosition in
*   the fAlign coordinate system.  The color set is
*   pColorCell.  Makes a copy of pColorCell.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetColorSimple(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                  Int4 lPosition, Boolean fAlign, DDV_ColorCell *pColorCell)
{
    DDV_Position Position;

    Position.lPosition = lPosition;
    Position.fAlign = fAlign;
    return DDV_SetColor(pColorGlobal, sip, &Position, pColorCell);
}

Int4 DDV_SetColor(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                  DDV_Position *pPosition, DDV_ColorCell *pColorCell)
{
    DDV_MediaInfo *pMediaInfo;

    if(pColorGlobal == NULL || pColorCell == NULL || sip == NULL ||
        pPosition == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SetColor");
        return 0;
    }

    pMediaInfo = DDV_GetMediaInfo(pColorGlobal, sip);
    if(pMediaInfo == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "No MediaInfo in DDV_SetColor");
        return 0;
    }

    return DDV_SetColorInMediaInfo(pColorGlobal, pMediaInfo, pPosition,
                                   pColorCell);

}

/*****************************************************************************
*
*   Set a color for an SeqId at lPosition in
*   the fAlign coordinate system.  Uses Uint1's for input.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetColorbyColor(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                         Int4 lPosition, Boolean fAlign, Uint1 red,
                         Uint1 green, Uint1 blue)
{
    DDV_MediaInfo *pMediaInfo;
    DDV_ColorCell ColorCell;
    DDV_Position Position;

    if(pColorGlobal == NULL || sip == NULL ) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_SetColor");
        return 0;
    }

    Position.fAlign = fAlign;
    Position.lPosition = lPosition;

    pMediaInfo = DDV_GetMediaInfo(pColorGlobal, sip);
    if(pMediaInfo == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "No MediaInfo in DDV_SetColor");
        return 0;
    }

    ColorCell.rgb[0] = red;
    ColorCell.rgb[1] = green;
    ColorCell.rgb[2] = blue;
    return DDV_SetColorInMediaInfo(pColorGlobal, pMediaInfo, &Position,
                                   &ColorCell);

}


/*****************************************************************************
*
*   Calls the color functions.  pData is for user data.
*
*   Returns 1 on success, 0 on failure.
*
*   Note that this only updates the colors.  It doesn't tell the applications
*   to paint the colors.  This is done through OM_MSG_SETCOLOR with a new
*   regiontype TBD.
*
*   fColorByMaster in DDV_ColorGlobal (i.e. color all sequence like the master
*   sequence) is NOT handled by this function.  This
*   is because this code has, in general, no idea if it is using alignment
*   coordinates or alignment coordinates.  To color by master, the user must
*   register a specific function to do this.
*
*   lFrom and lTo are the range coloring is done over.
*
*****************************************************************************/

Int4 DDV_ColorExecute(DDV_ColorGlobal *pColorGlobal, void * pData,
                      DDV_Range *pRange)
{
    DDV_ColorQueue  *pcq, *pcqPrev;
    ValNode *pvn, *pvnPrev, *pvnPrevPrev, *pvnTemp;
    Int4 lCount = 0, i;

    if(pColorGlobal == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_ColorExecute");
        return 0;
    }

    DDV_FreePalette(&pColorGlobal->Palette);

    /* color by default? */
    if(pColorGlobal->fDefaultColor) DDV_DefaultColor(pColorGlobal,
        NULL, NULL);

    /* look for an override */
    for(pvn = pColorGlobal->pvnColorQueue; pvn != NULL; pvn = pvn->next) {
        pcq = (DDV_ColorQueue *)pvn->data.ptrvalue;
        if(pcq->fOverride) {  /* this function is king */
            (*pcq->pfnColorFunc)(pColorGlobal, pData, pRange);
            return 1;
        }
        lCount++;
    }
    
    /* sort the color functions by priority */

    for(i = 0; i < lCount - 1; i++) {
        for(pvnPrevPrev = NULL,
            pvnPrev = pColorGlobal->pvnColorQueue, 
            pvn = pColorGlobal->pvnColorQueue->next;
            pvn != NULL;
            pvnPrevPrev = pvnPrev,
            pvnPrev = pvn,
            pvn = pvn->next) {

            pcq = (DDV_ColorQueue *)pvn->data.ptrvalue;
            pcqPrev = (DDV_ColorQueue *)pvnPrev->data.ptrvalue;

            if(pcq->lPriority < pcqPrev->lPriority) {
                pvnTemp = pvn->next;
                pvn->next = pvnPrev;
                pvnPrev->next = pvnTemp;
                if(pvnPrevPrev != NULL) pvnPrevPrev->next = pvn;
                else pColorGlobal->pvnColorQueue = pvn;
                pvnTemp = pvnPrev;
                pvnPrev = pvn;
                pvn = pvnTemp;
            }
        }
    }

    /* call the color functions in priority order */
    for(pvn = pColorGlobal->pvnColorQueue; pvn != NULL; pvn = pvn->next) {
        pcq = (DDV_ColorQueue *)pvn->data.ptrvalue;
        (*pcq->pfnColorFunc)(pColorGlobal, pData, pRange);
    }

    return 1;
}

/*****************************************************************************
*
*   Callback for coloring all ColorCells by the default color.
*
*
*****************************************************************************/

void DDV_DefaultColor(DDV_ColorGlobal *pColorGlobal, void *pData,
                                DDV_Range *pRange)
{
    ValNode *pvnMediaInfo, *pvnColor;
    DDV_MediaInfo *pMediaInfo;
    DDV_Color *pColor;
    DDV_ColorCell *pccDefault;
    Int4 i;

    if(pColorGlobal == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DefaultColor");
        return;
    }

    pccDefault = DDV_SearchColorCellbyName(pColorGlobal->pvnSpecialColors,
        "Default");

    if(pColorGlobal == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "DDV_DefaultColor: no default color");
        return;
    }

    for(pvnMediaInfo = pColorGlobal->pvnMediaInfo; pvnMediaInfo != NULL;
        pvnMediaInfo = pvnMediaInfo->next) {
        pMediaInfo = (DDV_MediaInfo *)pvnMediaInfo->data.ptrvalue;
        for(pvnColor = pMediaInfo->pvnColor; pvnColor != NULL; 
            pvnColor = pvnColor->next) {
            pColor = (DDV_Color *)pvnColor->data.ptrvalue;
                for(i = 0; i <= pColor->Range.lTo - pColor->Range.lFrom; i++)
                    DDV_CopyColorCell(&pColor->pColorCell[i], pccDefault);
        }
    }
        /* put the default color on palette */
    DDV_RequestColor(&pColorGlobal->Palette, pccDefault);
    return;
}

/*****************************************************************************
*
*   Adds a DDV_ColorQueue to pColorGlobal
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_AddColorFunc(DDV_ColorGlobal *pColorGlobal,
                       DDV_ColorQueue *pColorQueue)
{
    if(pColorGlobal == NULL || pColorQueue == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_AddColorQueue");
        return 0;
    }

    ValNodeAddPointer(&pColorGlobal->pvnColorQueue, 0, pColorQueue);
    return 1;
}

/*****************************************************************************
*
*   Deletes the given DDV_ColorFunc out of the ColorQueue
*   Returns number deleted on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_DeleteColorFunc(DDV_ColorGlobal *pColorGlobal,
                       DDV_ColorFunc  pfnColorFunc)
{
    ValNode *pvn, *pvnPrevious = NULL;
    DDV_ColorQueue *pColorQueue;
    Int4 retval = 0;

    if(pColorGlobal == NULL || pfnColorFunc == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid call on DDV_DeleteColorQueue");
        return 0;
    }

    for(pvn = pColorGlobal->pvnColorQueue; pvn != NULL; pvn = pvn->next) {
        pColorQueue = (DDV_ColorQueue *)pvn->data.ptrvalue;
        if(pColorQueue->pfnColorFunc == pfnColorFunc) {
            DDV_DeleteColorQueue(pColorQueue);
            retval++;
            if(pvnPrevious != NULL) {  /* splice out node */
                pvnPrevious->next = pvn->next;
            }
            else {  /* if first node */
                pColorGlobal->pvnColorQueue = pvn->next;
            }
            MemFree(pvn);
            break;
        }
        pvnPrevious = pvn;
    }

    return retval;
}

