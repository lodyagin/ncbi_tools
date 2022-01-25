/*  $Id: udvdef.h,v 6.3 1999/07/27 13:57:41 durand Exp $
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
* File Name:  udvdef.h
*
* Author:  Patrick Durand
*
* Version Creation Date:   07/09/99
*
* $Revision: 6.3 $
*
* File Description: this file is the companion of udviewer.h
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: udvdef.h,v $
* Revision 6.3  1999/07/27 13:57:41  durand
* modify display type defines
*
* Revision 6.2  1999/07/27 13:08:53  durand
* add display type defines
*
* Revision 6.1  1999/07/09 13:54:55  durand
* this is a companion file of udviewer.h
*
*
*
* ==========================================================================
*/

#ifndef _UDVDEF_
#define _UDVDEF_

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************

  INCLUDE SECTION 

*******************************************************************************/
#include <ncbi.h>

/*******************************************************************************

  DEFINE SECTION 

*******************************************************************************/

	/*number of letter by block: default is ten-letter blocks*/
#define LETTER_BLOCK_WIDTH 10

	/*margins around the UnDViewer*/
#define VIEWER_HORZ_MARGIN 10
#define VIEWER_VERT_MARGIN 10

	/*Numerical scale */
#define SCALE_POS_LEFT 1	/*scale position*/
#define SCALE_POS_TOP  2
#define SCALE_POS_BOTH  3
#define SCALE_POS_NONE  4
#define SCALE_MAJOR_TICK_EVERY 10	/*position of the ticks*/
#define SCALE_MINOR_TICK_EVERY 5
#define SCALE_LETTER_COLOR (GetColorRGB(0,0,255))	/*scale colours*/
#define SCALE_MAJTICK_COLOR (GetColorRGB(0,0,0))
#define SCALE_MINTICK_COLOR (GetColorRGB(128,128,128))
#define SCALE_WIDTH_IF_LEFT 9 /*nb. of cxChar (see Font struct); used to
                                   the scale on the left*/

	/*Font*/	
#define FONT_DEF_COLOR (GetColorRGB(0,0,0))/*default font colour*/

	/*Panel*/	
#define PANEL_NAME_WIDTH 10 /*nb. of cxChar (see Font struct) MAX: 50 !!;
	                            used to put names pn the left on the UDV window*/

	/*Info Panel*/	
#define BUFFER_SIZE  255 /*nb. of cxChar for Info Buffer*/

	/*NA and AA color layout*/	
#define NA_A_LAYOUT	0
#define NA_G_LAYOUT	1
#define NA_C_LAYOUT	2
#define NA_T_LAYOUT	3
#define NA_U_LAYOUT	4
#define LAYOUT_UPPER_CASE	1
#define LAYOUT_LOWER_CASE	2
	
	/*Mouse actions*/
#define MS_ACTION_FEAT_NOTHING	1 /*no action*/
#define MS_ACTION_FEAT_CURSOR	2 /* double cursor for features*/
#define MS_ACTION_RESIZE_WIN	3 /*resize cxName region*/	

	/*************************************************************************

	  I defined _min_ & _max_ because I needed the test '>=' in both 
	  _max_ and _min_

	*************************************************************************/
#define _min_(a,b)        ((a)>=(b)?(b):(a))
#define _max_(a,b)        ((a)>=(b)?(a):(b))

	/*used to draw features*/
		/*     |===>......===>.....===>|   : a big arrow feature*/
#define FEATURE_START_BOX		1	/*used to draw big arrows*/
#define FEATURE_START_ARROW		2
#define FEATURE_START_ARROW_END	3
#define FEATURE_START_NOTHING	4
		/* /\/\/\/\   : a helix feature (2D structures)*/
#define DRAW_HELIX_DOWN		1  /*used to draw helix faeture*/
#define DRAW_HELIX_MIDDLE	2
#define DRAW_HELIX_UP		3
		/*  >------<  : a bond feature*/
#define BOND_RIGHT 1	/*used to draw a bond feature*/
#define BOND_LEFT  2
#define SZBUF_SIZE 250
	
	/*sequence buffer management ; depending on VScroll up/down*/
#define BUFFER_REPOP_VCRL_LUP	1   /*line up*/
#define BUFFER_REPOP_VCRL_LDN	2   /*line down*/
#define BUFFER_REPOP_VCRL_PUP	3   /*page up*/
#define BUFFER_REPOP_VCRL_PDN	4   /*page down*/

	/*UDV/DDV panel region types*/
#define UDV_REGION_NAMELIST      ((Uint1)1)
#define UDV_REGION_SEPARATOR     ((Uint1)2)
#define UDV_REGION_PARAGLIST     ((Uint1)3)
#define UDV_REGION_SEQALIGNRULER ((Uint1)4) /*use by DDV only*/

#define DDV_DISP_HORZ ((Uint4)1)
#define DDV_DISP_VERT ((Uint4)2)
	
#ifdef __cplusplus
}
#endif

#endif /* _UDVDEF_ */

