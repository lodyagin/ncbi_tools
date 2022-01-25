/*   $Id: ddvcolor.h,v 6.2 1999/07/19 19:12:21 lewisg Exp $
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
* File Name:  $Id: ddvcolor.h,v 6.2 1999/07/19 19:12:21 lewisg Exp $
*
* Author:  Lewis Geer
*
* Version Creation Date:   6/2/99
*
* $Revision: 6.2 $
*
* File Description: Shared color information for viewers
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ddvcolor.h,v $
* Revision 6.2  1999/07/19 19:12:21  lewisg
* fix bug adding new pColor in old pMediaInfo.  Also added ability to detect overlaps when allocating new pColor
*
* Revision 6.1  1999/07/16 18:46:46  lewisg
* moved ddvcolor from api to tools
*
* Revision 1.4  1999/07/13 23:24:48  lewisg
* added DDV_SipList
*
* Revision 1.3  1999/07/13 17:02:00  lewisg
* fix forward declaration error in sgi compiler
*
* Revision 1.2  1999/07/13 14:38:39  lewisg
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
* Revision 1.1  1999/06/11 23:37:10  lewisg
* color management functions
*
*
*
* ==========================================================================
*/


#ifndef DDVCOLOR_H
#define DDVCOLOR_H


#ifdef __cplusplus
extern "C" {
#endif

#include <ncbilcl.h>
#include <ncbistd.h>
#include <objall.h>


/* a single coordinate on a sequence */
typedef struct _DDV_Position {
    Int4 lPosition;
    Boolean fAlign;  /* is this an alignment coord?  otherwise bioseq coord */
} DDV_Position;

/* a range on a sequence */
typedef struct _DDV_Range {
    Int4 lFrom;
    Int4 lTo;
    Boolean fAlign;  /* is this an alignment coord?  otherwise bioseq coord */
} DDV_Range;

typedef struct _DDV_ColorCell {
  Uint1 rgb[3];  /* standard rgb color cell */
} DDV_ColorCell;

typedef struct _DDV_Color {
    DDV_Range Range;  /* the coordinates of the entity covered by the colorcell */
    DDV_ColorCell *pColorCell;  /* the array of colors */
} DDV_Color;

    /* standard color entry */
typedef struct _DDV_ColorEntry {
    DDV_ColorCell ColorCell;  /* the color cell */
    Char * Name;  /* the name of the color */
    Char * Paints;  /* what the color paints */
    Int4 Index; /* an index assigned to the colors after they are allocated */
} DDV_ColorEntry;

typedef struct _DDV_MediaInfo {
SeqId *sip;  /* this should be a duplicate */
ValNode *pvnColor;  /* ValNode list of DDV_Color associated with the object */
Boolean fVisible;  /* is the sequence visible? */
} DDV_MediaInfo;

/*****************************************************************************
*
*   This is the global structure that contains all of the color info
*
*****************************************************************************/

typedef struct _DDV_ColorGlobal {
ValNode *pvnMediaInfo;  /* ValNode list of DDV_MediaInfo */
SeqId *MasterSeqId;  /* this refers to the master sequence */
Boolean fColorByMaster;  /* use only the master sequence to color */
Boolean fDefaultColor;   /* fDefaultColor means run the default coloration */
ValNode *pvnColorQueue;  /* Valnode list of DDV_ColorQueue */
ValNode *pvnAllColorFns; /* Valnode list of all color functions that can
                            be exposed externally to the viewers.  to be used
                            by SAM */
ValNode *Palette;         /* ValNode list of DDV_ColorEntries */
ValNode *pvnSpecialColors;  /* ValNode list of DDV_ColorEntries.  used for
                            highlight, secondary structure, etc. */
} DDV_ColorGlobal;

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
*   DDV_ColorFunc takes as arguments pData for user data and a range
*   pRange.  If To < From, the range should be ignored.
*
*****************************************************************************/

/* standard lPriority's */
#define DDVPRIHI 20
#define DDVPRIMED 0
#define DDVPRILO -20

typedef void (*DDV_ColorFunc)(DDV_ColorGlobal *pColorGlobal, void *pData,
                                DDV_Range *pRange);

typedef struct _DDV_ColorQueue {
Char * szName;      /* the menu name.  must be unique if used by SAM to add
                        and delete functions */
DDV_ColorFunc  pfnColorFunc;
Int4 lPriority;     /* order that the functions are called */
Boolean fOverride;  /* ignore priority and only call this function */
struct _DDV_ColorQueue * next; /* these are used to sort the color functions */
struct _DDV_ColorQueue * prev;  /* *not* for storing them together. */
} DDV_ColorQueue;

/*****************************************************************************
*
*   Set default values for DDV_Position.  
*
*****************************************************************************/

Int4 DDV_DefaultPosition(DDV_Position *pPosition, Int4 lPosition, 
                         Boolean fAlign);

/*****************************************************************************
*
*   Set default values for DDV_Range.  
*
*****************************************************************************/

Int4 DDV_DefaultRange(DDV_Range *pRange, Int4 lFrom, Int4 lTo, 
                         Boolean fAlign);

/*****************************************************************************
*
*   Constructor for DDV_ColorEntry structure.  DDV_ColorEntry is a
*   (Red, Green, Blue) color triplet and associated szName.
*
*****************************************************************************/

DDV_ColorEntry * DDV_CreateColorEntry(Char *szName, Nlm_Uint1 Red,
                          Nlm_Uint1 Green, Nlm_Uint1 Blue);

/*****************************************************************************
*
*   Destructor for DDV_ColorEntry Color.
*
*****************************************************************************/

void DDV_DeleteColorEntry(DDV_ColorEntry *Color);

/*****************************************************************************
*
*   Frees a Palette.
*   Note that it frees the entire Valnode list pointed to by Palette and
*   sets Palette to NULL.
*
*****************************************************************************/

void DDV_FreePalette(ValNode **Palette);

/*****************************************************************************
*
*   Sets the Index in a DDV_ColorEntry according to position in Palette
*
*****************************************************************************/

void DDV_SetColorChoice(ValNode *Palette);

/*****************************************************************************
*
*   Looks for ColorCell on the Palette.  Returns * to the DDV_ColorEntry
*
*****************************************************************************/

DDV_ColorEntry * DDV_SearchColor(ValNode *Palette, DDV_ColorCell *ColorCell);

/*****************************************************************************
*
*   Looks for a DDV_ColorEntry on the Palette by szName.  Returns *
*   to the DDV_ColorEntry.
*
*****************************************************************************/

DDV_ColorEntry * DDV_SearchColorbyName(ValNode *Palette, Char *szName);

/*****************************************************************************
*
*   Looks for a DDV_ColorEntry on the Palette by szName.  Returns *
*   to ColorCell in the DDV_ColorEntry.
*
*****************************************************************************/

DDV_ColorCell * DDV_SearchColorCellbyName(ValNode *Palette, Char *szName);

/*****************************************************************************
*
*   Puts ColorCell on the Palette list if it isn't there already.
*
*****************************************************************************/

void DDV_RequestColor(ValNode **Palette, DDV_ColorCell *ColorCell);

/*****************************************************************************
*
*   Puts ColorEntry on the Palette list if the sZname isn't there already,
*   Otherwise replaces the entry.  Does NOT make a copy of the DDV_ColorEntry
*
*****************************************************************************/

void DDV_RequestColorbyName(ValNode **Palette, DDV_ColorEntry *pColorIn);

/*****************************************************************************
*
*   Returns the index into the Palette for a given ColorCell.
*   Return index on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_ColorIndex(ValNode *Palette, DDV_ColorCell *ColorCell);


/*****************************************************************************
*
*   Copy a DDV_ColorCell.
*   Return 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_CopyColorCell(DDV_ColorCell * pDestination,  DDV_ColorCell * pSource);

/*****************************************************************************
*
*   Set a color cell using Uint1 values
*   Return 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_SetColorCell(DDV_ColorCell * pColorCell, Nlm_Uint1 Red,
                      Nlm_Uint1 Green, Nlm_Uint1 Blue);

/*****************************************************************************
*
*   Checks to see if a pPosition is inside pRange.
*   If it is, return 1
*   If it isn't, return 0
*   If pPosition and pRange don't have the same coordinate system, return -1
*
*****************************************************************************/

Int4 DDV_InRange(DDV_Position *pPosition, DDV_Range *pRange);

/*****************************************************************************
*
*   Checks to see if pRange1 overlaps pRange2.
*   If it does completely, return DDV_TOTALLAP
*   If it doesn't, return DDV_NOLAP
*   If pRange1 overlaps the front of pRange2, return DDV_FRONTLAP
*   If pRange1 overlaps the rear of pRange2, return DDV_BACKLAP
*   If pPosition and pRange don't have the same coordinate system, return -1
*
*****************************************************************************/

Int4 DDV_RangeOverlap(DDV_Range *pRange1, DDV_Range *pRange2);

#define DDV_NOLAP 0
#define DDV_TOTALLAP 1
#define DDV_FRONTLAP 2
#define DDV_BACKLAP 3

/*****************************************************************************
*
*   Returns a DDV_ColorCell for a given position if it is inside of pColor.
*   Returns NULL otherwise.
*
*****************************************************************************/

DDV_ColorCell * DDV_GetCellByPosition(DDV_Color *pColor, DDV_Position *pPosition);

/*****************************************************************************
*
*   Sets a DDV_ColorCell for a given position if it is inside of pColor.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetCellByPosition(DDV_Color * pColor, DDV_Position *pPosition,
                 DDV_ColorCell *pColorCell);

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
                                        Boolean fAlign);

/*****************************************************************************
*
*   Deletes a DDV_Color.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_DeleteColor(DDV_Color *pColor);

/*****************************************************************************
*
*   Creates a new DDV_MediaInfo structure.  Creates a default color structure
*   for the given object specified by sip.  Duplicates the sip.
*   Inserts the given pColor as valnode entry.  Sets fVisible.
*
*   Returns NULL on error.
*
*****************************************************************************/

DDV_MediaInfo * DDV_NewMediaInfo(SeqId *sip, DDV_Color *pColor,
                                  Boolean fVisible);

/*****************************************************************************
*
*   Creates a default MediaInfo structure and adds it to pColorGlobal
*   for the given object specified by sip. Duplicates the sip.  Sequence is
*   visible by default.  Coordinates lFrom to lTo and using fAlign coordinates.
*
*   Returns 0 on error.
*
*****************************************************************************/

Int4 DDV_DefaultMediaInfoByLen(DDV_ColorGlobal * pColorGlobal, SeqId *sip,
                               Int4 lFrom, Int4 lTo, Boolean fAlign);

/* simplified version of above that retrieves the length and sets the coordinate
   system to sequence */

Int4 DDV_DefaultMediaInfo(DDV_ColorGlobal * pColorGlobal, SeqId *sip);

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

ValNode * DDV_SipList(DDV_ColorGlobal *pColorGlobal);

/*****************************************************************************
*
*   Deletes a DDV_MediaInfo structure.
*
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_DeleteMediaInfo(DDV_MediaInfo *pMediaInfo);


/*****************************************************************************
*
*   Returns a new color queue entry with name szName, priority lPriority.
*   fOverride is explained above.
*
*****************************************************************************/

DDV_ColorQueue * DDV_NewColorQueue(DDV_ColorFunc pfnColorFunc, Char * szName,
                              Int4 lPriority, Boolean fOverride);

/*****************************************************************************
*
*   Deletes a DDV_ColorQueue.  
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_DeleteColorQueue(DDV_ColorQueue *pColorQueue);

/*****************************************************************************
*
*   Sets up default secondary structure colors in 
*   pColorGlobal->pvnSpecialColors palette.
*   Returns 0 on failure.
*
*****************************************************************************/

Int4 DDV_DefaultSSColor(DDV_ColorGlobal * pColorGlobal);

/*****************************************************************************
*
*   Returns a new DDV_ColorGlobal structure.
*   Returns NULL on failure
*
*****************************************************************************/

DDV_ColorGlobal * DDV_CreateColorGlobal(Boolean fDefaultColor);

/*****************************************************************************
*
*   Deletes a  DDV_ColorGlobal structure.
*   Returns 1 on success, 0 on failture
*
*****************************************************************************/

Int4 DDV_DeleteColorGlobal(DDV_ColorGlobal *pColorGlobal);

/*****************************************************************************
*
*   Delete all information for the given SeqId.
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_RemoveMediaInfo(DDV_ColorGlobal *pColorGlobal, SeqId *sip);

/*****************************************************************************
*
*   Add MediaInfo to DDV_ColorGlobal
*   Returns 1 on success, 0 on failure
*
*****************************************************************************/

Int4 DDV_AddMediaInfo(DDV_ColorGlobal *pColorGlobal,
                              DDV_MediaInfo *pMediaInfo);

/*****************************************************************************
*
*   Compare two SeqId to make sure all valnodes compare exactly
*
*****************************************************************************/

Boolean DDV_SeqIdCompare(SeqId *sip1, SeqId *sip2);

/*****************************************************************************
*
*   Given an SeqId, return the associated DDV_MediaInfo.
*
*****************************************************************************/

DDV_MediaInfo * DDV_GetMediaInfo(DDV_ColorGlobal *pColorGlobal, SeqId *sip);

/*****************************************************************************
*
*   Given an SeqId, return the visible state of the SeqId.
*
*****************************************************************************/

Boolean DDV_IsVisible(DDV_ColorGlobal *pColorGlobal, SeqId *sip);

/*****************************************************************************
*
*   Set the visible state (fVisible) of a SeqId.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetVisible(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                    Boolean fVisible);

/*****************************************************************************
*
*   Retrieve a color for an SeqId at position pPositionor at lPosition in
*   the fAlign coordinate system.
*   Returns a DDV_ColorCell containing the color.  NULL on failure
*
*   The DDV_ColorCell returned is read only.
*
*****************************************************************************/

DDV_ColorCell * DDV_GetColorSimple(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                             Int4 lPosition, Boolean fAlign);

DDV_ColorCell * DDV_GetColor(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                             DDV_Position *pPosition);

/*****************************************************************************
*
*   Set a color for in a DDV_MediaInfo at position pPosition.  The color set is
*   pColorCell.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetColorInMediaInfo(DDV_ColorGlobal *pColorGlobal,
                             DDV_MediaInfo *pMediaInfo, DDV_Position *pPosition,
                             DDV_ColorCell *pColorCell);

/*****************************************************************************
*
*   Set a color for an SeqId at position pPosition or at lPosition in
*   the fAlign coordinate system.  The color set is
*   pColorCell.  Makes a copy of pColorCell.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetColorSimple(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                  Int4 lPosition, Boolean fAlign, DDV_ColorCell *pColorCell);

Int4 DDV_SetColor(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                  DDV_Position *pPosition, DDV_ColorCell *pColorCell);

/*****************************************************************************
*
*   Set a color for an SeqId at lPosition in
*   the fAlign coordinate system.  Uses Uint1's for input.
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_SetColorbyColor(DDV_ColorGlobal *pColorGlobal, SeqId *sip,
                         Int4 lPosition, Boolean fAlign, Uint1 red,
                         Uint1 green, Uint1 blue);

/*****************************************************************************
*
*   Calls the color functions.  pData is for user data.
*
*   Returns 1 on success, 0 on failure.
*
*   fColorByMaster in DDV_ColorGlobal (i.e. color all sequence like the master
*   sequence) is NOT handled by this function.  This
*   is because this code has, in general, no idea if it is using alignment
*   coordinates or alignment coordinates.  To color by master, the user must
*   register a specific function to do this.
*
*   pRange the range coloring is done over.
*
*****************************************************************************/

Int4 DDV_ColorExecute(DDV_ColorGlobal *pColorGlobal, void * pData,
                      DDV_Range *pRange);

/*****************************************************************************
*
*   Callback for coloring all ColorCells by the default color.
*
*
*****************************************************************************/

void DDV_DefaultColor(DDV_ColorGlobal *pColorGlobal, void *pData,
                                DDV_Range *pRange);

/*****************************************************************************
*
*   Adds a DDV_ColorQueue to pColorGlobal
*   Returns 1 on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_AddColorFunc(DDV_ColorGlobal *pColorGlobal,
                       DDV_ColorQueue *pColorQueue);

/*****************************************************************************
*
*   Deletes the given DDV_ColorFunc out of the ColorQueue
*   Returns number deleted on success, 0 on failure.
*
*****************************************************************************/

Int4 DDV_DeleteColorFunc(DDV_ColorGlobal *pColorGlobal,
                       DDV_ColorFunc  pfnColorFunc);


#ifdef __cplusplus
}
#endif

#endif /* DDVCOLOR_H */
