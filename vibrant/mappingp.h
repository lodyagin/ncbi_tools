/*   mappingP.h
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
* File Name:  mappingP.h
*
* Author:  Jonathan Kans, Alex Smirnov, Jill Shermer
*
* Version Creation Date:   1/19/93
*
* $Revision: 6.0 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* $Log: mappingp.h,v $
* Revision 6.0  1997/08/25 18:55:52  madden
* Revision changed to 6.0
*
* Revision 5.0  1996/05/28 13:45:08  ostell
* Set to revision 5.0
*
 * Revision 4.1  1995/12/14  17:29:30  kans
 * moved SearchSegment from viewer.c to picture.c, made extern
 *
 * Revision 4.0  1995/07/26  13:51:04  ostell
 * force revision to 4.0
 *
 * Revision 1.10  1995/05/17  15:15:14  kans
 * added Log line
 *
*
* ==========================================================================
*/

#ifndef _MAPPINGP_
#define _MAPPINGP_

#ifndef _VIBRANT_
#include <vibrant.h>
#endif

#ifndef _PICTURE_
#include <picture.h>
#endif

#ifndef _PICTUREP_
#include <picturep.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
*
*   STRUCTURE TYPEDEFS
*
*****************************************************************************/

typedef struct Nlm_scaleinfo {
  Nlm_BoxInfo  worldWindow;
  Nlm_BoxInfo  worldWindow16;
  Nlm_Int4     scaleX;
  Nlm_Int4     scaleY;
  Nlm_Int4     offsetX;
  Nlm_Int4     offsetY;
} Nlm_ScaleInfo, PNTR Nlm_ScalePtr;

typedef struct Nlm_drawinfo {
  Nlm_AttPPtr   primattrib;
  Nlm_ScaleInfo scale;
  Nlm_AttPData  curattrib;
  Nlm_Boolean   checked;
  Nlm_Int1      highlight;
} Nlm_DrawInfo, PNTR Nlm_DrawInfoPtr;

typedef struct Nlm_vscaleinfo {
  Nlm_BoxInfo  world;
  Nlm_RecT     view;
  Nlm_BoxInfo  port;
  Nlm_Int4     scaleX;
  Nlm_Int4     scaleY;
  Nlm_Int2     scrollX;
  Nlm_Int2     scrollY;
  Nlm_Boolean  force;
} Nlm_VScaleInfo, PNTR Nlm_VScalePtr;

extern void Nlm_LoadBox PROTO((BoxPtr box, Int4 left, Int4 top, Int4 right, 
                               Int4 bottom));
extern void Nlm_OutsetBox PROTO((BoxPtr box, Int4 dX, Int4 dY ));
extern void Nlm_MapWorldPointToPixel PROTO((Nlm_PointPtr pt, Nlm_PntPtr pnt, 
                               Nlm_ScalePtr scale));
extern void Nlm_MapPixelPointToWorld PROTO((Nlm_PntPtr pnt, Nlm_PointPtr pt, 
                               Nlm_ScalePtr scale));
extern void Nlm_MapWorldBoxToRect PROTO((Nlm_RectPtr r, Nlm_BoxPtr box, 
                               Nlm_ScalePtr scale));
extern void Nlm_MapRectToWorldBox PROTO(( Nlm_BoxPtr box, Nlm_RectPtr r,
                               Nlm_ScalePtr scale));
extern Nlm_Boolean Nlm_BoxInViewport PROTO((Nlm_RectPtr rct, Nlm_BoxPtr box, 
                               Nlm_VScalePtr scale));
extern Nlm_Boolean Nlm_LineIntoVPort PROTO(( Nlm_Int4Ptr x1, Nlm_Int4Ptr y1, 
                        Nlm_Int4Ptr x2, Int4Ptr y2, Nlm_BoxPtr worldWindow ));
extern Nlm_Boolean Nlm_IsLineInVPort PROTO(( Nlm_Int4 x1, Nlm_Int4 y1, 
                        Nlm_Int4 x2, Nlm_Int4 y2, Nlm_BoxPtr worldWindow ));
extern Nlm_SegmenT Nlm_SearchSegment PROTO((Nlm_SegmenT segment,
                        Nlm_ScalePtr scalePtr, Nlm_PrimitivE PNTR prPtr ));

#define ScaleInfo Nlm_ScaleInfo
#define ScalePtr Nlm_ScalePtr
#define DrawInfo Nlm_DrawInfo
#define DrawInfoPtr Nlm_DrawInfoPtr 
#define VScaleInfo Nlm_VScaleInfo 
#define VScalePtr Nlm_VScalePtr
#define LoadBox Nlm_LoadBox
#define OutsetBox Nlm_OutsetBox
#define MapWorldPointToPixel Nlm_MapWorldPointToPixel
#define MapPixelPointToWorld Nlm_MapPixelPointToWorld
#define MapWorldBoxToRect Nlm_MapWorldBoxToRect
#define MapRectToWorldBox Nlm_MapRectToWorldBox
#define BoxInViewport Nlm_BoxInViewport
#define LineIntoVPort Nlm_LineIntoVPort
#define IsLineInVPort Nlm_IsLineInVPort
#define SearchSegment Nlm_SearchSegment

#ifdef __cplusplus
}
#endif

#endif
