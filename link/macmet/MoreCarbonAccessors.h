/*   carbonmacros.h
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* File Name:  carbonmacros.h
*
* Author:  Joshua Juran
*
* Version Creation Date:   2001-03-21
*
* $Revision: 1.3 $
*
* File Description:
*		header for supporting Carbon calls in pre-Carbon targets
*		Mac only
*
* ==========================================================================
*/
#ifndef _MORECARBONACCESSORS_
#define _MORECARBONACCESSORS_

#include <ConditionalMacros.h>

// Carbonization
// -------------

// The problem:
// The Carbon API requires that we treat Toolbox data structures as opaque -- 
// that we not access their members directly, but use accessor functions.
// Since we'd like to use the same code across API targets, we'll use these accessors
// for pre-Carbon Mac OS as well as Carbon.  Some of the Carbon accessor functions 
// are defined as inline functions in the Universal Interfaces for just this purpose, 
// but some are not.  For these, we define our own macros.

// This file must be included AFTER ConditionalMacros.h.  To make life simple,
// we include it.

#if !ACCESSOR_CALLS_ARE_FUNCTIONS
# if UNIVERSAL_INTERFACES_VERSION < 0x0335  // Universal Interfaces 3.4 beta
// GetWindowFromPort is defined in UI 3.4
#  define GetWindowFromPort(grafPtr)                (grafPtr)
# endif
# define GetPortVisibleRegion(grafPtr, rgnHandle)  CopyRgn((grafPtr)->visRgn, (rgnHandle))
# define GetRegionBounds(rgnH, rectPtr)            (*(rectPtr) = (*(RgnHandle)(rgnH))->rgnBBox)
# define GetPortBounds(grafPtr, rectPtr)           (*(rectPtr) = (grafPtr)->portRect)
# define GetQDGlobalsArrow(cursorPtr)              (*(cursorPtr) = qd.arrow)
# define GetQDGlobalsScreenBits(bitmapPtr)         (*(bitmapPtr) = qd.screenBits)
# define GetPortBitMapForCopyBits(grafPtr)         (&(grafPtr)->portBits)
# define GetPortTextFont(grafPtr)                  ((grafPtr)->txFont)
# define GetPortTextFace(grafPtr)                  ((grafPtr)->txFace)
# define GetPortTextSize(grafPtr)                  ((grafPtr)->txSize)
#endif  // ACCESSOR_CALLS_ARE_FUNCTIONS

#if !TARGET_API_MAC_CARBON
# define EnableMenuItem(theMenu, item)             EnableItem(theMenu, item)
# define DisableMenuItem(theMenu, item)            DisableItem(theMenu, item)
#endif  // TARGET_API_MAC_CARBON

#endif
