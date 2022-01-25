/*   algorend.c
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
* File Name:  algorend.c
*
* Authors:  Christopher Hogue, Lewis Geer, Yanli Wang
*
* Version Creation Date:
*
* File Description: Algorithmic rendering routines and Vibrant controls for Cn3d
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: algorend.c,v $
* Revision 6.68  1999/05/11 17:54:17  ywang
* move Cn3d_AnyPrim out of OpenGL ifdef block in RenderAnAtom
*
* Revision 6.67  1999/05/05 21:12:57  ywang
* fix typing error bug in RenderGraph
*
* Revision 6.66  1999/04/14 16:40:00  ywang
* color objects by secondary structure for strucseqs
*
* Revision 6.65  1999/04/06 22:27:38  lewisg
* fix endless loop
*
* Revision 6.64  1999/04/06 16:06:10  coremake
* Fixes of typos
*
* Revision 6.63  1999/04/06 14:23:29  lewisg
* add opengl replacement for viewer3d
*
* Revision 6.62  1999/03/30 22:36:20  ywang
* add functions to color salsa for NcbiMimeAsn1_strucseqs & code reorganization
*
* Revision 6.61  1999/03/22 22:41:53  ywang
* redesign feature page, fix bugs
*
* Revision 6.60  1999/03/09 18:12:56  lewisg
* turn on termini labels
*
* Revision 6.59  1999/03/08 21:16:22  ywang
* initialize variables
*
* Revision 6.58  1999/03/01 20:21:31  ywang
* put in options for residue on/off control
*
* Revision 6.57  1999/02/24 23:06:37  ywang
* minor name change for MGD feature
*
* Revision 6.56  1999/02/10 23:49:41  lewisg
* use RGB values instead of indexed palette
*
* Revision 6.55  1999/01/27 21:51:51  ywang
* add label to 'Model' menu
*
* Revision 6.54  1999/01/26 17:14:35  ywang
* redesign Display menu and add 'display highlight residues only' function
*
* Revision 6.53  1999/01/20 22:57:24  ywang
* customize color for secondary structure & rearrange Option menu
*
 * Revision 6.52  1998/12/23  19:40:56  ywang
 * fix layer bug while drawing user-defined feature
 *
 * Revision 6.51  1998/12/16  19:44:00  kans
 * improvements to labels - combine name, letter, pdb controls into single popup
 *
* Revision 6.50  1998/12/16 01:32:08  kans
* first pass at improving label controls
*
* Revision 6.49  1998/12/14 18:39:55  kans
* enumerated alists to limit render and color choices
*
* Revision 6.48  1998/12/10 15:51:07  kans
* improved layout of render controls, safe set value/status to avoid flicker
*
* Revision 6.47  1998/12/10 01:32:06  kans
* first pass at improved render dialog
*
* Revision 6.46  1998/11/09 22:03:32  ywang
* fix bugs for modeling
*
 * Revision 6.45  1998/11/06  23:01:06  ywang
 * fix bugs for modeling
 *
 * Revision 6.44  1998/11/04  00:06:21  ywang
 * add function for modeling: change render/color for special residue(s)
 *
 * Revision 6.43  1998/10/28  19:02:46  kans
 * included salparam.h
 *
* Revision 6.42  1998/10/27 15:56:55  ywang
* add functions for testing color by sequence conservation
*
 * Revision 6.41  1998/10/21  21:18:00  ywang
 * color salsa by cycle_molecule  and by residue
 *
 * Revision 6.40  1998/10/21  15:51:28  ywang
 * reset residue color for salsa before cn3d redraws so that residues shown in salsa will become black if they are not shown in cn3d window
 *
 * Revision 6.39  1998/10/01  21:55:40  ywang
 * put function for object display control
 *
 * Revision 6.38  1998/09/30  22:10:47  ywang
 * control display on three levels: structure, chain, domain
 *
 * Revision 6.37  1998/09/24  15:52:49  kans
 * included cn3dmodl.h
 *
* Revision 6.36  1998/09/23 22:04:04  ywang
* synchronize show/hide between cn3d and salsa when display complexity is changed
*
 * Revision 6.35  1998/09/23  18:38:51  ywang
 * add functions to control display on domain level
 *
 * Revision 6.34  1998/09/22  17:53:46  ywang
 * add display control on MM level
 *
 * Revision 6.33  1998/09/17  20:53:51  ywang
 * set default color display as by conservation for struc+align case and replace the 'magenta' color with 'red'
 *
 * Revision 6.32  1998/08/26  18:28:35  kans
 * fixed -v -fd warnings
 *
* Revision 6.31  1998/07/29 22:17:36  lewisg
* change color of 1st aligned structure
*
* Revision 6.30  1998/07/28 22:17:39  lewisg
* turned on messaging when apply button pressed in render pane
*
* Revision 6.29  1998/07/28 21:14:17  lewisg
* renamed some items in the color by dropdown
*
* Revision 6.28  1998/07/09 21:50:23  ywang
* get rid of printf
*
 * Revision 6.26  1998/07/08  16:16:10  ywang
 * fix rendering bug for P backbone only model
 *
 * Revision 6.25  1998/07/07  19:31:59  ywang
 * fix rendering problem for model with Alpha Carbon onlyalgorend.c
 *
 * Revision 6.24  1998/07/02  15:05:17  lewisg
 * fix bug with resetting presentation on null structure
 *
* Revision 6.23  1998/06/16 18:00:23  lewisg
* moved rendering menus and created a reset presentation menu item
*
* Revision 6.22  1998/05/29 19:27:09  ywang
* fix potential bug
*
 * Revision 6.19  1998/05/26  21:51:00  ywang
 * limit color message
 *
 * Revision 6.18  1998/05/26  21:35:14  lewisg
 * added defaults to render menu, got rid of mouse 3D actions menu item
 *
* Revision 6.17  1998/05/22 22:23:21  lewisg
* always show heterogens
*
* Revision 6.16  1998/05/12 21:46:59  lewisg
* stricter conservation coloring
*
* Revision 6.15  1998/05/08 17:30:15  lewisg
* new conservation coloring
*
* Revision 6.14  1998/05/05 20:05:15  ywang
* add color by Protein to salsa
*
 * Revision 6.12  1998/04/30  15:23:40  ywang
 * call ColorSalsa from RenderGraph
 *
 * Revision 6.11  1998/04/28  22:47:16  lewisg
 * master/slave color in sync
 *
* Revision 6.10  1998/04/27 17:50:01  lewisg
* added color by conservation
*
* Revision 6.9  1998/04/21 23:00:51  lewisg
* added show aligned/unaligned
*
* Revision 6.8  1998/04/20 22:08:57  lewisg
* got rid of dead code
*
* Revision 6.7  1998/04/20 18:36:06  lewisg
* moved extern for Viewer3d to cn3dmain.h
*
* Revision 6.6  1998/04/18 00:33:46  lewisg
* added ability to turn slaves on/off
*
* Revision 6.5  1998/04/15 00:51:29  lewisg
* bug fixes for multiple alignment mode and alignment pane
*
* Revision 6.4  1998/04/06 04:25:20  lewisg
* added color by alignment
*
* Revision 6.3  1998/04/01 23:26:08  lewisg
* added new startup mode + fixed slave rendering
*
* Revision 6.2  1998/03/06 23:19:17  lewisg
* codewarrior fixes
*
* Revision 6.1  1998/03/06 01:16:43  lewisg
* merge
*
* Revision 6.0  1997/08/25 18:13:21  madden
* Revision changed to 6.0
*
* Revision 5.5  1996/08/19 21:04:29  vakatov
* Made functions of type "pNodeFunc" to be of that type(sharp)
* (and thus fixed fatal bug under Borland/WIN32); removed all
* castings to (pNodeFunc) -- to let compiler to check for such bugs
*
 * Revision 5.4  1996/08/06  16:04:06  hogue
 * Fixed a mistaken ==
 *
 * Revision 5.3  1996/07/31  18:35:50  hogue
 * Selection and segment-grouping to atom nodes added.
 *
 * Revision 5.2  1996/07/22  00:24:10  hogue
 * Added an origin 3D item for no-primitives condition and general use.
 *
 * Revision 5.1  1996/06/03  19:51:20  vakatov
 * VIEWSCALE decreased (caused "long" overfloat in the 3D drawing functions)
 *
 * Revision 5.0  1996/05/28  14:05:44  ostell
 * Set to revision 5.0
 *
 * Revision 1.13  1996/05/23  14:32:34  hogue
 * Made RenderLabel static for Codewarrior
 *
 * Revision 1.12  1996/05/22  22:10:34  kans
 * changed NULL to '\0'
 *
 * Revision 1.11  1996/05/22  21:46:25  hogue
 * Added White button to label controls.
 *
 * Revision 1.10  1996/05/22  20:48:53  hogue
 * Debugged the label controls, added some label features.
 *
 * Revision 1.9  1996/05/22  15:55:20  hogue
 * Changed label control panel, many other changes to support label drawing.
 * Also changed some FloatHi's to FloatLo's
 *
 * Revision 1.8  1996/05/14  15:19:41  hogue
 * Added label contols - but are not hooked up yet.
 *
 * Revision 1.7  1996/05/09  15:40:40  hogue
 * Domain rendering enabled.
 *
 * Revision 1.6  1996/04/26  21:44:48  vakatov
 * just casting...
 *
 * Revision 1.5  1996/04/26  18:42:08  vakatov
 * CN3D sources ported to MS-Windows;
 * the portability errors and warnings fixed, etc.
 *
*
* ==========================================================================
*/
#ifdef _OPENGL
#include <shim3d.h>
#else
#include <viewer3d.h>
#endif
#include <cn3dmain.h>
#include <math.h>
#include <mmdbapi.h>
#include <algorend.h>
#include <cn3dpane.h>
#include <cn3dmsg.h>
#include <cn3dmodl.h>
#include <saledit.h>
#include <salparam.h>
#include <cn3dshim.h>

static ENUM_ALIST(empty_alist)
  {"              ",     1},
END_ENUM_ALIST

static ENUM_ALIST(all_render_alist)
  {"Wireframe",          1},
  {"Tubes",              2},
  {"Ball & Stick",       3},
  {"Fat Tubes",          4},
  {"Space Fill   ",      5},
END_ENUM_ALIST

static ENUM_ALIST(ions_render_alist)
  {"Fat Tubes",          4},
  {"Space Fill   ",      5},
END_ENUM_ALIST

static ENUM_ALIST(connections_render_alist)
  {"Wireframe",          1},
  {"Tubes",              2},
  {"Ball & Stick",       3},
END_ENUM_ALIST

static ENUM_ALIST(solvents_render_alist)
  {"Ball & Stick",       3},
  {"Fat Tubes",          4},
  {"Space Fill   ",      5},
END_ENUM_ALIST

static ENUM_ALIST(all_color_alist)
  {"Molecule Cycle",     1},
  {"Sec. Structure",     2},
  {"Domain",             3},
  {"Residue",            4},
  {"Hydrophobicity",     5},
  {"CPK",                6},
  {"Temperature",        7},
  {"Object Cycle",       8},
  {"Structure",          9},
  {"Neighbor",          10},
END_ENUM_ALIST

static ENUM_ALIST(heterogen_color_alist)
  {"Molecule Cycle",     1},
  {"Residue",            4},
  {"Hydrophobicity",     5},
  {"CPK",                6},
  {"Temperature",        7},
  {"Object Cycle",       8},
  {"Structure",          9},
END_ENUM_ALIST

static ENUM_ALIST(ions_color_alist)
  {"Molecule Cycle",     1},
  {"Residue",            4},
  {"Hydrophobicity",     5},
  {"CPK",                6},
  {"Temperature",        7},
  {"Structure",          9},
END_ENUM_ALIST

static ENUM_ALIST(solvents_color_alist)
  {"Molecule Cycle",     1},
  {"Residue",            4},
  {"Hydrophobicity",     5},
  {"CPK",                6},
  {"Temperature",        7},
  {"Object Cycle",       8},
END_ENUM_ALIST

static ENUM_ALIST(objects_color_alist)
  {"Molecule Cycle",     1},
  {"Sec. Structure",     2},
  {"Domain",             3},
  {"Object Cycle",       8},
  {"Structure",          9},
END_ENUM_ALIST



/*define DEBUG_N 1 */

#define VIEWSCALE 1000000.0

#ifndef _OPENGL
    static Picture3D pic = NULL;
#endif

static Int2 Cn3d_ParsNow = 0;

static Int2 Cn3d_RenderNow = 0;
static Uint1 Cn3d_LabelNow = 1;
static Int2 Cn3d_ColorNow = 0;
static Int1 Cn3d_LayerNow = 0;
static Boolean Cn3d_DoHydrogens = FALSE;
static Boolean Cn3d_CopyToNode = FALSE;
static Boolean Cn3d_ColorPass = FALSE;  /* gathering unique colors for palette */
static Boolean Cn3d_AnyPrim = FALSE;
static int Cn3d_lSlaveNum = 0;  /* which slave being iterated over 0=master */


static PopuP    Cn3D_pupPBB;   		/* protein backbone options and on/off */
static PopuP	Cn3D_pupNABB;  		/* nucl. acid bb options and on/off */
static ButtoN 	Cn3D_lOnOffItem [10]; /* pieces parts to draw - first 2 empty */

static PopuP    Cn3D_pupRenderStyle [10];
static PopuP    Cn3D_pupColorStyle [10];
static EnumFieldAssocPtr renderAlist [10] =
  {all_render_alist, all_render_alist, all_render_alist, all_render_alist,
   all_render_alist, ions_render_alist, connections_render_alist,
   solvents_render_alist, empty_alist, empty_alist};
static EnumFieldAssocPtr colorAlist [10] =
  {all_color_alist, all_color_alist, all_color_alist, all_color_alist,
   heterogen_color_alist, ions_color_alist, empty_alist, solvents_color_alist,
   objects_color_alist, empty_alist};
static Int2 renderDefaults [10] = {1, 1, 1, 1, 1, 5, 1, 5, 1, 1};
static Int2 colorDefaults [10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


static ButtoN   Cn3D_bRedraw;  		/* the button to redraw */
static ButtoN   Cn3D_immediateRedraw;

static ButtoN   Cn3D_lOnOffLabel [4];
static PopuP    Cn3D_pupLabelAA;
static PopuP    Cn3D_pupLabelNT;
static PopuP    Cn3D_bLName [4];
static ButtoN   Cn3D_bLNum [4];
static ButtoN   Cn3D_bLTop [4];
static ButtoN   Cn3D_bLWhite [4];
static PopuP    Cn3D_pupLabelSize [4];

static ButtoN   Cn3D_bLRedraw;
static ButtoN   Cn3D_LimmediateRedraw;

extern CharPtr NCBIstdaaUC;
extern CharPtr NCBI4naUC;
extern CharPtr NCBI4naLC;
extern CharPtr NCBIstdaaLC;

extern Int1 KinAAColor[];
extern Int2 KinNAColor[];
extern Int1 KinAtoms[];
extern Int1 ColorNumKinBB[];
extern Int1 ColorNumKinSC[];
extern Int1 ColorNumKinAC[];
extern CharPtr KineColors[];
extern Int1 ElementKinColors[];
extern Int1 ThermKine[];
extern Int4 TempsKine[];
extern SpecialFeaturePtr sfp;
SpecialFeaturePtr sfpThis = NULL;
SpecialFeatureInfoPtr sfipThis = NULL;


Int1 bColorAlignments[NUM_SLAVES] = {C_sea, C_sky, C_brown, C_cyan, C_blue, C_gray, C_purple, C_green};  /* alignment colors */
Int1 bColorConservation[NUM_SLAVES] = {C_blue, C_sea, C_cyan, C_greentint, C_gold, C_yellow, C_yellowtint, C_white }; /* conservation colors */

Int1 PhobeAAColor[MAX_NCBIstdaa] = {21,20,14,
	9,12,12,20,21,4,20,4,20,20,14,14,14,4,14,14,20,20,21,20,14,10,1 };
/* PINK RED BLUE YELLOW GOLD selcys BROWN GREY 26 */


/***********************************/

ResidueColorCell * LIBCALL ColorFromSS(PMGD pmgdThis)
{
  Int2 rsult = 0;

  switch ((int) (pmgdThis->bNCBISecStru)) {
     case SS_HELIX:
        return &(Cn3D_ColorData.SSColors[CN3D_COLOR_HELIX].ColorCell);
     case SS_SHEET:
     case SS_STRAND:
        return &(Cn3D_ColorData.SSColors[CN3D_COLOR_STRAND].ColorCell);
     case SS_TURN:
        return &(Cn3D_ColorData.SSColors[CN3D_COLOR_TURN].ColorCell);
     default:
        return &(Cn3D_ColorData.SSColors[CN3D_COLOR_COIL].ColorCell);
  }
  return NULL;

}

void LIBCALL SetDefaultAlgorRender(PARS pars)
{
  if (!pars) return;

  pars->HydrogensOn = FALSE;
  pars->DisorderOn = FALSE;
  pars->AnimateOn = TRUE;

  pars->PVirtualBBOn = TRUE;
  pars->PRealBBOn = FALSE;
  pars->PExtraBBOn = FALSE;
  pars->PBBRender = R_WIRE;
  pars->PBBColor = C_BYSSTRU;

  pars->PBBLabelOn = 1;
  pars->PBBLabelJust = (Uint1) (LA_CENTER);
  pars->PBBLabelStyle = (Uint1) L_NAME | L_NUM | L_3LETR;
  pars->PBBLabelScale = 2;

  pars->PTermLabelOn = FALSE;
  pars->PTermLabelJust = (Uint1) (LA_CENTER);
  pars->PTermLabelStyle = (Uint1) 0;
  pars->PTermLabelScale = 4;

  pars->PResiduesOn = FALSE;
  pars->PResRender = R_WIRE;
  pars->PResColor = C_BYSSTRU;

  pars->NTVirtualBBOn = FALSE;
  pars->NTRealBBOn = FALSE;
  pars->NTExtraBBOn = TRUE;
  pars->NTBBRender = R_THICKWIRE;
  pars->NTBBColor = C_CPK;


  pars->NTBBLabelOn = 1;
  pars->NTBBLabelJust = (Uint1) (LA_CENTER);
  pars->NTBBLabelStyle = (Uint1) L_NAME | L_NUM | L_1LETR;
  pars->NTBBLabelScale = 2;

  pars->NTResiduesOn = TRUE;
  pars->NTResRender = R_WIRE;
  pars->NTResColor = C_CPK;

  pars->NTTermLabelOn = FALSE;
  pars->NTTermLabelJust = (Uint1) (LA_CENTER);
  pars->NTTermLabelStyle = (Uint1) 0;
  pars->NTTermLabelScale = 4;

  pars->HeterogensOn = TRUE;
  pars->HetRender = R_THICKWIRE;
  pars->HetColor = C_CPK;

  pars->IonsOn = TRUE;
  pars->IonRender = R_SPACE;
  pars->IonColor = C_CPK;

  pars->ConnectOn = FALSE;
  pars->ConnectRender = R_WIRE;
  pars->ConnectColor = 9; /* all same color?? */

  pars->SolventOn = FALSE;
  pars->SolventRender = R_BALLNSTICK;
  pars->SolventColor = C_CPK;

  pars->ObjectOn = TRUE;
  pars->ObjectRender = R_DEFAULT;
  pars->ObjectColor = C_BYSSTRU;

   return;
}

void LIBCALL SetAlignAlgorRender(PARS pars)
{
  if (!pars) return;

  pars->HydrogensOn = FALSE;
  pars->DisorderOn = FALSE;
  pars->AnimateOn = TRUE;

  pars->PVirtualBBOn = TRUE;
  pars->PRealBBOn = FALSE;
  pars->PExtraBBOn = FALSE;
  pars->PBBRender = R_THICKWIRE;
/*pars->PBBColor = C_BYALIGN;  */
  pars->PBBColor = C_BYCONS;    /*  yanli */

  pars->PBBLabelOn = 1;
  pars->PBBLabelJust = (Uint1) (LA_CENTER);
  pars->PBBLabelStyle = (Uint1) L_NAME | L_NUM | L_3LETR;
  pars->PBBLabelScale = 2;

  pars->PTermLabelOn = FALSE;
  pars->PTermLabelJust = (Uint1) (LA_CENTER);
  pars->PTermLabelStyle = (Uint1) 0;
  pars->PTermLabelScale = 4;

  pars->PResiduesOn = FALSE;
  pars->PResRender = R_WIRE;
/*pars->PResColor = C_BYALIGN;  */   /* yanli */
  pars->PResColor = C_BYCONS;

  pars->NTVirtualBBOn = FALSE;
  pars->NTRealBBOn = FALSE;
  pars->NTExtraBBOn = TRUE;
  pars->NTBBRender = R_THICKWIRE;
  pars->NTBBColor = C_CPK;


  pars->NTBBLabelOn = 1;
  pars->NTBBLabelJust = (Uint1) (LA_CENTER);
  pars->NTBBLabelStyle = (Uint1) L_NAME | L_NUM | L_1LETR;
  pars->NTBBLabelScale = 2;

  pars->NTResiduesOn = TRUE;
  pars->NTResRender = R_WIRE;
  pars->NTResColor = C_CPK;

  pars->NTTermLabelOn = FALSE;
  pars->NTTermLabelJust = (Uint1) (LA_CENTER);
  pars->NTTermLabelStyle = (Uint1) 0;
  pars->NTTermLabelScale = 4;

  pars->HeterogensOn = TRUE;
  pars->HetRender = R_THICKWIRE;
  pars->HetColor = C_CPK;

  pars->IonsOn = TRUE;
  pars->IonRender = R_SPACE;
  pars->IonColor = C_CPK;

  pars->ConnectOn = FALSE;
  pars->ConnectRender = R_WIRE;
  pars->ConnectColor = 9; /* all same color?? */

  pars->SolventOn = FALSE;
  pars->SolventRender = R_BALLNSTICK;
  pars->SolventColor = C_CPK;

  pars->ObjectOn = FALSE;
  pars->ObjectRender = R_DEFAULT;
  pars->ObjectColor = C_BYALIGN;

   return;
}

void LIBCALL SetStrucSeqsAlgorRender(PARS pars)
{

  if (!pars) return;

  pars->HydrogensOn = FALSE;
  pars->DisorderOn = FALSE;
  pars->AnimateOn = TRUE;

  pars->PVirtualBBOn = TRUE;
  pars->PRealBBOn = FALSE;
  pars->PExtraBBOn = FALSE;
  pars->PBBRender = R_THICKWIRE;
/*pars->PBBColor = C_BYALIGN;  */
  pars->PBBColor = C_BYCONS;    /*  yanli */

  pars->PBBLabelOn = 1;
  pars->PBBLabelJust = (Uint1) (LA_CENTER);
  pars->PBBLabelStyle = (Uint1) L_NAME | L_NUM | L_3LETR;
  pars->PBBLabelScale = 2;

  pars->PTermLabelOn = FALSE;
  pars->PTermLabelJust = (Uint1) (LA_CENTER);
  pars->PTermLabelStyle = (Uint1) 0;
  pars->PTermLabelScale = 4;

  pars->PResiduesOn = FALSE;
  pars->PResRender = R_WIRE;
  pars->PResColor = C_BYCONS;

  pars->NTVirtualBBOn = FALSE;
  pars->NTRealBBOn = FALSE;
  pars->NTExtraBBOn = TRUE;
  pars->NTBBRender = R_THICKWIRE;
  pars->NTBBColor = C_CPK;

  pars->NTBBLabelOn = 1;
  pars->NTBBLabelJust = (Uint1) (LA_CENTER);
  pars->NTBBLabelStyle = (Uint1) L_NAME | L_NUM | L_1LETR;
  pars->NTBBLabelScale = 2;

  pars->NTResiduesOn = TRUE;
  pars->NTResRender = R_WIRE;
  pars->NTResColor = C_CPK;

  pars->NTTermLabelOn = FALSE;
  pars->NTTermLabelJust = (Uint1) (LA_CENTER);
  pars->NTTermLabelStyle = (Uint1) 0;
  pars->NTTermLabelScale = 4;

  pars->HeterogensOn = TRUE;
  pars->HetRender = R_THICKWIRE;
  pars->HetColor = C_CPK;

  pars->IonsOn = TRUE;
  pars->IonRender = R_SPACE;
  pars->IonColor = C_CPK;

  pars->ConnectOn = FALSE;
  pars->ConnectRender = R_WIRE;
  pars->ConnectColor = 9; /* all same color?? */

  pars->SolventOn = FALSE;
  pars->SolventRender = R_BALLNSTICK;
  pars->SolventColor = C_CPK;

  pars->ObjectOn = TRUE;
  pars->ObjectRender = R_DEFAULT;
  pars->ObjectColor = C_BYSSTRU;

   return;
}


PARS LIBCALL NewAlgorRenderSet(void)
{
  PARS par = NULL;
  par = (PARS) MemNew((size_t)(sizeof(ARS)));
  if (!par) return NULL;
  SetDefaultAlgorRender(par);
  return par;
}

PARS LIBCALL NewAlignRenderSet(void)
{
  PARS par = NULL;
  par = (PARS) MemNew((size_t)(sizeof(ARS)));
  if (!par) return NULL;
  SetAlignAlgorRender(par);
  return par;
}

PARS LIBCALL NewStrucSeqsRenderSet(void)
{
  PARS par = NULL;
  par = (PARS) MemNew((size_t)(sizeof(ARS)));
  if (!par) return NULL;
  SetStrucSeqsAlgorRender(par);
  return par;
}

void LIBCALL FreeAlgorRenderSet(PARS pars)
{
  MemFree(pars);
}


static PARS GetThisSpecialAlgorRenderSet(Int2 choice)
{
   SpecialFeaturePtr sfpThis = NULL;
   SpecialFeatureInfoPtr sfipThis = NULL;
   PARS pars = NULL;

   sfpThis = sfp;
   while(sfpThis){
      if(sfpThis->choice == choice){
         sfipThis = sfpThis->data.ptrvalue;
         pars = sfipThis->parsSpecial;
         return(pars);
      }
      sfpThis = sfpThis->next;
   }

   return(pars);
}

PARS LIBCALL GetAlgorRenderSet(PDNMS pdnmsThis)
/* returns rendering information stucture.  lyg */
{
   PARS pars = NULL;
   PMSD pmsdThis;
   if (!pdnmsThis) return NULL;
   pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

   if (pmsdThis->pExtra == NULL)
     {
       pmsdThis->pExtra = NewAlgorRenderSet();
       if (pmsdThis->pExtra == NULL) return NULL;
     }
   return (PARS) pmsdThis->pExtra;
}


static void fnCn3D_RedrawWrapper(ButtoN b);

static void MasterChangeRenderProc (void)
{
  Int2 i, j, k;
  PARS pars = NULL;
  PMSD pmsdThis = NULL;
  PDNMS pdnmsThis = NULL;
  UIEnum val;

 /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;

  pars->PResiduesOn = GetStatus (Cn3D_lOnOffItem [2]);
  pars->NTResiduesOn = GetStatus (Cn3D_lOnOffItem [3]);
  pars->HeterogensOn = GetStatus (Cn3D_lOnOffItem [4]);
  pars->IonsOn = GetStatus (Cn3D_lOnOffItem [5]);
  pars->ConnectOn = GetStatus (Cn3D_lOnOffItem [6]);
  pars->SolventOn = GetStatus (Cn3D_lOnOffItem [7]);
  pars->ObjectOn = GetStatus (Cn3D_lOnOffItem [8]);
  pars->HydrogensOn = GetStatus (Cn3D_lOnOffItem [9]);

  pars->PVirtualBBOn = FALSE;
  pars->PRealBBOn = FALSE;
  pars->PExtraBBOn = FALSE;

  i = GetValue (Cn3D_pupPBB);
  switch (i)
   {
   	case 1:
    	   pars->PVirtualBBOn = TRUE; /* alpha c trace */
  	   break;
        case 2:
           pars->PRealBBOn = TRUE; /* partial atoms */
           break;
        case 3:
           pars->PExtraBBOn = TRUE; /* all atoms */
           break;
        default: ; /* none */
   }

  pars->NTVirtualBBOn = FALSE;
  pars->NTRealBBOn = FALSE;
  pars->NTExtraBBOn = FALSE;

  j = GetValue (Cn3D_pupNABB);
  switch (j)
   {
   	case 1:
    	   pars->NTVirtualBBOn = TRUE;
  	   break;
        case 2:
           pars->NTRealBBOn = TRUE;
           break;
        case 3:
           pars->NTExtraBBOn = TRUE;
           break;
        default: ;
   }

  for (i = 0; i < 9; i++) {
  if (GetEnumPopup (Cn3D_pupColorStyle [i], colorAlist [i], &val)) {
    j = (Int2) val;
  } else {
    j = 1;
  }
  /* j = GetValue (Cn3D_pupColorStyle [i]); */
   switch (j)
	  {
            case 1:
	      k = C_BYCHAIN;
	      break;
            case 2:
	      k = C_BYSSTRU;
	      break;
            case 3:
	      k = C_BYDOMAIN;
	      break;
	    case 4:
	       k = C_BYRES;
	      break;
	    case 5:
	      k = C_BYHYDRO;
	      break;
	    case 6:
	       k = C_CPK;
	      break;
	    case 7:
	       k = C_BYTEMP;
	      break;
	    case 8:
	      k = C_BYOBJECT;
	      break;
		case 9:
	      k = C_BYALIGN;
	      break;
		case 10:
	      k = C_BYCONS;
	      break;
	    default:
	      k = 0;
	  }
  switch (i)
    {
     case 0: /* prot bb */
         pars->PBBColor = k;
         break;
     case 1: /* na bb */
         pars->NTBBColor = k;
		 break;
     case 2: /* prot sc */
          pars->PResColor = k;
          break;
     case 3: /* na sc */
         pars->NTResColor = k;
         break;
     case 4: /* hets */
         pars->HetColor = k;
         break;
     case 5: /* ions */
         pars->IonColor = k;
         break;
     case 6: /* connections */
         /* pars->ConnectColor = k; */
         pars->ConnectColor = C_yellow;
         break;
     case 7: /* solvent */
         pars->SolventColor = k;
         break;
     case 8: /* object */
         pars->ObjectColor = k;
     default:
     ;
   }
   }

  for (i = 0; i < 9; i++) {
  if (GetEnumPopup (Cn3D_pupRenderStyle [i], renderAlist [i], &val)) {
    j = (Int2) val;
  } else {
    j = 1;
  }
  /* j = GetValue(Cn3D_pupRenderStyle [i]); */
  switch (j)
	 {
	   case 1:
	     k =  R_WIRE;
	     break;
	   case 2:
	     k = R_STICK;
	     break;
	   case 3:
	     k = R_BALLNSTICK;
	     break;
	   case 4:
	     k = R_THICKWIRE;
	     break;
	   case 5:
	     k = R_SPACE;
	    break;
     default:
     ;
	  }
  switch (i)
    {
     case 0: /* prot bb */
         pars->PBBRender = k;
         break;
     case 1: /* na bb */
         pars->NTBBRender = k;
         break;
     case 2: /* prot sc */
          pars->PResRender = k;
          break;
     case 3: /* na sc */
         pars->NTResRender = k;
         break;
     case 4: /* hets */
         pars->HetRender = k;
         break;
     case 5: /* ions */
         pars->IonRender = k;
         break;
     case 6: /* connections */
         pars->ConnectRender = k;
         break;
     case 7: /* solvent */
         pars->SolventRender = k;
         break;
     case 8: /* object */
         /* pars->ObjectRender = k; */
         pars->ObjectRender = R_DEFAULT;
     default:
	   ;
   }
   }

  if (GetStatus (Cn3D_immediateRedraw)) {
    fnCn3D_RedrawWrapper (Cn3D_bRedraw);
  }
}

static void ButtonOnOffProc (ButtoN b)

{
  MasterChangeRenderProc ();
}

static void PopupOnOffProc (PopuP p)

{
  MasterChangeRenderProc ();
}

void Cn3D_RenStruc(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
	pars->HydrogensOn = FALSE;
	pars->ObjectOn = TRUE;
	pars->SolventOn = FALSE;
	pars->PRealBBOn = FALSE;
	pars->PResiduesOn = FALSE;
	pars->PVirtualBBOn = TRUE;
	pars->PExtraBBOn = FALSE;
	pars->NTVirtualBBOn = FALSE;
	pars->NTRealBBOn = FALSE;
	pars->NTExtraBBOn = TRUE;
	pars->NTResiduesOn = FALSE;
	pars->HeterogensOn = TRUE;
	pars->IonsOn = TRUE;
	pars->ConnectOn = FALSE;
	pars->PBBRender = R_WIRE;
	pars->PResRender = R_WIRE;
	pars->NTBBRender = R_WIRE;
	pars->NTResRender = R_WIRE;
	pars->HetRender = R_WIRE;
	pars->IonRender = R_SPACE;
	pars->ConnectRender = R_WIRE;
	pars->SolventRender = R_BALLNSTICK;
	pars->ObjectRender = R_DEFAULT;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected(); 
  return;
}

void Cn3D_RenWire(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;
  PMSD pmsdThis = NULL;
  PDNML pdnmlThis = NULL;
  PMLD  pmldThis = NULL;


  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
 	pars->HydrogensOn = FALSE;
	pars->ObjectOn = FALSE;
	pars->SolventOn = FALSE;
	pars->PRealBBOn = FALSE;
	pars->PVirtualBBOn = FALSE;
	pars->PExtraBBOn = TRUE;
	pars->PResiduesOn = TRUE;
	pars->NTVirtualBBOn = FALSE;
	pars->NTRealBBOn = FALSE;
        pars->NTResiduesOn = TRUE;
	pars->NTExtraBBOn = TRUE;
	pars->HeterogensOn = TRUE;
	pars->IonsOn = TRUE;
	pars->ConnectOn = TRUE;
	pars->PBBRender = R_WIRE;
	pars->PResRender = R_WIRE;
	pars->NTBBRender = R_WIRE;
	pars->NTResRender = R_WIRE;
	pars->HetRender = R_WIRE;
	pars->IonRender = R_SPACE;
	pars->ConnectRender = R_WIRE;
	pars->SolventRender = R_BALLNSTICK;
	pars->ObjectRender = R_DEFAULT;

  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  return;
}

void Cn3D_RenTube(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
      	pars->HydrogensOn = FALSE;
	pars->ObjectOn = FALSE;
	pars->SolventOn = FALSE;
	pars->PRealBBOn = FALSE;
	pars->PVirtualBBOn = FALSE;
	pars->PExtraBBOn = TRUE;
	pars->PResiduesOn = TRUE;
	pars->NTVirtualBBOn = FALSE;
	pars->NTRealBBOn = FALSE;
        pars->NTResiduesOn = TRUE;
	pars->NTExtraBBOn = TRUE;
	pars->HeterogensOn = TRUE;
	pars->IonsOn = TRUE;
	pars->ConnectOn = TRUE;
	pars->PBBRender = R_THICKWIRE;
  	pars->PResRender = R_STICK;
 	pars->NTBBRender = R_THICKWIRE;
  	pars->NTResRender = R_STICK;
  	pars->HetRender = R_THICKWIRE;
  	pars->IonRender = R_SPACE;
 	pars->ConnectRender = R_STICK;
  	pars->SolventRender = R_BALLNSTICK;
 	pars->ObjectRender = R_DEFAULT;
   ResetRenderCtrls();
   Cn3D_Redraw(FALSE);
   Cn3dObjMgrGetSelected();
  return;
}


void Cn3D_RenHier(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
	pars->HydrogensOn = FALSE;
	pars->ObjectOn = FALSE;
	pars->SolventOn = FALSE;
	pars->PRealBBOn = FALSE;
	pars->PVirtualBBOn = FALSE;
	pars->PExtraBBOn = TRUE;
	pars->PResiduesOn = TRUE;
	pars->NTVirtualBBOn = FALSE;
	pars->NTRealBBOn = FALSE;
        pars->NTResiduesOn = TRUE;
	pars->NTExtraBBOn = TRUE;
	pars->HeterogensOn = TRUE;
	pars->IonsOn = TRUE;
	pars->ConnectOn = TRUE;
 	pars->PBBRender = R_THICKWIRE;
 	pars->PResRender = R_WIRE;
   	pars->NTBBRender = R_THICKWIRE;
  	pars->NTResRender = R_WIRE;
   	pars->HetRender = R_STICK;
   	pars->IonRender = R_SPACE;
   	pars->ConnectRender = R_STICK;
  	pars->SolventRender = R_BALLNSTICK;
  	pars->ObjectRender = R_DEFAULT;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  return;
}

void Cn3D_RenSpace(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
  	pars->HydrogensOn = TRUE;
	pars->ObjectOn = FALSE;
	pars->SolventOn = FALSE;
	pars->PRealBBOn = FALSE;
	pars->PVirtualBBOn = FALSE;
	pars->PExtraBBOn = TRUE;
	pars->PResiduesOn = TRUE;
	pars->NTVirtualBBOn = FALSE;
	pars->NTRealBBOn = FALSE;
        pars->NTResiduesOn = TRUE;
	pars->NTExtraBBOn = TRUE;
	pars->HeterogensOn = TRUE;
	pars->IonsOn = TRUE;
	pars->ConnectOn = TRUE;
	pars->PBBRender = R_SPACE;
	pars->PResRender = R_SPACE;
 	pars->NTBBRender = R_SPACE;
 	pars->NTResRender = R_SPACE;
 	pars->HetRender = R_SPACE;
 	pars->IonRender = R_SPACE;
 	pars->ConnectRender = R_WIRE;
 	pars->SolventRender = R_BALLNSTICK;
 	pars->ObjectRender = R_DEFAULT;
   ResetRenderCtrls();
   Cn3D_Redraw(FALSE);
   Cn3dObjMgrGetSelected();
  return;
}

void Cn3D_RenBS(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
	pars->HydrogensOn = FALSE;
	pars->ObjectOn = FALSE;
	pars->SolventOn = FALSE;
	pars->PRealBBOn = FALSE;
	pars->PVirtualBBOn = FALSE;
	pars->PExtraBBOn = TRUE;
	pars->PResiduesOn = TRUE;
	pars->NTVirtualBBOn = FALSE;
	pars->NTRealBBOn = FALSE;
        pars->NTResiduesOn = TRUE;
	pars->NTExtraBBOn = TRUE;
	pars->HeterogensOn = TRUE;
	pars->IonsOn = TRUE;
	pars->ConnectOn = TRUE;
	pars->PBBRender = R_BALLNSTICK;
 	pars->PResRender = R_BALLNSTICK;
 	pars->NTBBRender = R_BALLNSTICK;
 	pars->NTResRender = R_BALLNSTICK;
 	pars->HetRender = R_BALLNSTICK;
 	pars->IonRender = R_SPACE;
 	pars->ConnectRender = R_STICK;
 	pars->SolventRender = R_BALLNSTICK;
 	pars->ObjectRender = R_DEFAULT;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  return;
}

/**
Render with alignment settings.
*/

void Cn3D_RenAlign(IteM i)   
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
	pars->HydrogensOn = FALSE;
	pars->ObjectOn = FALSE;
	pars->SolventOn = FALSE;
	pars->PRealBBOn = FALSE;
	pars->PVirtualBBOn = TRUE;
	pars->PExtraBBOn = FALSE;
	pars->PResiduesOn = FALSE;
	pars->NTVirtualBBOn = FALSE;
	pars->NTRealBBOn = FALSE;
        pars->NTResiduesOn = TRUE;
	pars->NTExtraBBOn = TRUE;
	pars->HeterogensOn = TRUE;
	pars->IonsOn = TRUE;
	pars->ConnectOn = FALSE;
	pars->PBBRender = R_THICKWIRE;
 	pars->PResRender = R_WIRE;
 	pars->NTBBRender = R_THICKWIRE;
 	pars->NTResRender = R_WIRE;
 	pars->HetRender = R_THICKWIRE;
 	pars->IonRender = R_SPACE;
 	pars->ConnectRender = R_WIRE;
 	pars->SolventRender = R_BALLNSTICK;
 	pars->ObjectRender = R_WIRE;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  return;
}

void Cn3D_RenDefault(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;
  PMSD pmsdThis = NULL;


  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();
 
  if (!pdnmsThis) return;
  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
  if (pmsdThis->pdnmsSlaves) SetAlignAlgorRender(pars);
  else SetDefaultAlgorRender(pars);
  ResetRenderCtrls();
  ResetDisplayCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  return;
}



void Cn3D_ColCPK(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  Cn3D_ReColor = TRUE;

  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
        pars->PBBColor = C_CPK;
	pars->PResColor = C_CPK;
	pars->NTBBColor = C_CPK;
	pars->NTResColor = C_CPK;
	pars->HetColor = C_CPK;
	pars->IonColor = C_CPK;
	pars->SolventColor = C_CPK;
  ResetRenderCtrls();
  ResetSalsaColor(); 
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  Cn3D_ReColor = FALSE;
  return;
}


void Cn3D_ColDomain(IteM i)
{
   PARS pars = NULL;
   PDNMS pdnmsThis = NULL;

   Cn3D_ReColor = TRUE;

   /* fetch the active structure */
   pdnmsThis = GetSelectedModelstruc();

   if (!pdnmsThis) return;
   pars = GetAlgorRenderSet(pdnmsThis);
   if (!pars) return;
         pars->PBBColor = C_BYDOMAIN;
 	pars->PResColor = C_BYDOMAIN;
 	pars->ObjectColor = C_BYDOMAIN;
   ResetRenderCtrls();
   Cn3D_Redraw(FALSE);
   Cn3dObjMgrGetSelected();
   Cn3D_ReColor = FALSE;
   return;
}


void Cn3D_ColCy(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;

  Cn3D_ReColor = TRUE;

  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
        pars->PBBColor = C_BYCHAIN;
	pars->PResColor = C_BYCHAIN;
	pars->NTBBColor = C_BYCHAIN;
	pars->NTResColor = C_BYCHAIN;
	pars->HetColor = C_BYCHAIN;
	pars->ObjectColor = C_BYCHAIN;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
   Cn3D_ReColor = FALSE;
  return;
}

void Cn3D_ColStru(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;

  Cn3D_ReColor = TRUE;

  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
        pars->PBBColor = C_BYSSTRU;
	pars->PResColor = C_BYSSTRU;
	pars->ObjectColor = C_BYSSTRU;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  Cn3D_ReColor = FALSE;
  return;
}

void Cn3D_ColAlign(IteM i)  /* menu function for setting up color by alignement */
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  Cn3D_ReColor = TRUE;

  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
        pars->PBBColor = C_BYALIGN;
	pars->PResColor = C_BYALIGN;
	pars->ObjectColor = C_BYALIGN;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  Cn3D_ReColor = FALSE;
  return;
}

void Cn3D_ColCons(IteM i)  /* menu function for setting up color by conservation */
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  Cn3D_ReColor = TRUE;

  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
        pars->PBBColor = C_BYCONS;
	pars->PResColor = C_BYCONS;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  Cn3D_ReColor = FALSE;
  return;
}

void Cn3D_ColSeqCons(IteM i)  /* menu function for setting up color by conservation */
{
  
  SeqEditViewProcsPtr svpp;
  IteM  conscolor;
  IteM  iByAlign = NULL;

  Cn3D_ColAlign(iByAlign);

  svpp = (SeqEditViewProcsPtr) GetAppProperty ("SeqEditDisplayForm");
  conscolor = svpp->conscolor;
  ColorIdentityDialogItem ( conscolor );
  
  return;
}

void Cn3D_DoColSeqCons(void)
{  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;
 
  Cn3D_ReColor = TRUE;

  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
        pars->PBBColor = C_BYSEQCONS;
    pars->PResColor = C_BYSEQCONS;
  ResetRenderCtrls();

  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  Cn3D_ReColor = FALSE;
  return;
}

void Cn3D_ColRes(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  Cn3D_ReColor = TRUE;

  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
	pars->PBBColor = C_BYRES;
	pars->PResColor = C_BYRES;
	pars->NTBBColor = C_BYRES;
	pars->NTResColor = C_BYRES;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  Cn3D_ReColor = FALSE;
  return;
}

void Cn3D_ColHydro(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;

 
  Cn3D_ReColor = TRUE;
  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;

        pars->PBBColor = C_BYHYDRO;
	pars->PResColor = C_BYHYDRO;
  ResetRenderCtrls();
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  Cn3D_ReColor = FALSE;
  return;
}


void Cn3D_ColTemp(IteM i)
{
  PARS pars = NULL;
  PDNMS pdnmsThis = NULL;


  Cn3D_ReColor = TRUE;

  /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;
  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;
	pars->PBBColor = C_BYTEMP;
	pars->PResColor = C_BYTEMP;
	pars->NTBBColor = C_BYTEMP;
	pars->NTResColor = C_BYTEMP;
	pars->HetColor = C_BYTEMP;
	pars->IonColor = C_BYTEMP;
	pars->SolventColor = C_BYTEMP;
  ResetRenderCtrls();
  ResetSalsaColor(); 
  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
  Cn3D_ReColor = FALSE;

  return;
}



/* turns the alignment on and off using hack global.  lyg */
/*static void Cn3D_SetAlign(ButtoN b)
{
	PMSD pmsdThis = NULL;
	PDNMS pdnmsThis = NULL;
	
	/* fetch the active structure 
	pdnmsThis = GetSelectedModelstruc();

	if (!pdnmsThis) return;
	pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
	fShowAlignment = GetStatus(Cn3D_bAlign);
	if(fShowAlignment)  xu(pmsdThis);
	return;
}
*/

void LIBCALL ResetRenderCtrls(void)
{
  Int2 i, j, k;
  PARS pars = NULL;
  PMSD pmsdThis = NULL;
  PDNMS pdnmsThis = NULL;
  PDNML pdnmlThis = NULL;
  PMLD  pmldThis = NULL;
  

 /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) goto setout;

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars)
   {
    setout:
      for (k = 2; k < 10; k++) {
        SafeSetStatus (Cn3D_lOnOffItem [k], FALSE);
      }

      SafeSetValue(Cn3D_pupPBB, 4);
      SafeSetValue(Cn3D_pupNABB, 4);
      for (k = 0; k < 9; k++) {
        SetEnumPopup (Cn3D_pupRenderStyle [k], renderAlist [k], (UIEnum) renderDefaults [k]);
        SetEnumPopup (Cn3D_pupColorStyle [k], colorAlist [k], (UIEnum) colorDefaults [k]);
        /*
        SafeSetValue(Cn3D_pupRenderStyle [k], 1);
        SafeSetValue(Cn3D_pupColorStyle [k], 1);
        */
      }

      return;
   }

  SafeSetStatus (Cn3D_lOnOffItem [2], pars->PResiduesOn);
  SafeSetStatus (Cn3D_lOnOffItem [3], pars->NTResiduesOn);
  SafeSetStatus (Cn3D_lOnOffItem [4], pars->HeterogensOn);
  SafeSetStatus (Cn3D_lOnOffItem [5], pars->IonsOn);
  SafeSetStatus (Cn3D_lOnOffItem [6], pars->ConnectOn);
  SafeSetStatus (Cn3D_lOnOffItem [7], pars->SolventOn);
  SafeSetStatus (Cn3D_lOnOffItem [8], pars->ObjectOn);
  SafeSetStatus (Cn3D_lOnOffItem [9], pars->HydrogensOn);

/* yanli added */
/* alpha C trace only control */
  pdnmlThis = pmsdThis->pdnmlModels;
  pmldThis = (PMLD)pdnmlThis->data.ptrvalue;
  if(pmldThis->iType == Model_type_ncbi_backbone) {
     pars->PVirtualBBOn = TRUE;    
     pars->NTVirtualBBOn = TRUE;
     SafeSetValue(Cn3D_pupPBB, 1);
     SafeSetValue(Cn3D_pupNABB, 1);
  }
    

/*BACKBONE controls */
  if (pars->PVirtualBBOn)  SafeSetValue(Cn3D_pupPBB, 1);
  else
  if (pars->PRealBBOn) SafeSetValue(Cn3D_pupPBB, 2);
  else
  if (pars->PExtraBBOn) SafeSetValue(Cn3D_pupPBB, 3);
  else SafeSetValue(Cn3D_pupPBB, 4);

/* set status of BB items from PARS */
  if (pars->NTVirtualBBOn)  SafeSetValue(Cn3D_pupNABB, 1);
  else
  if (pars->NTRealBBOn) SafeSetValue(Cn3D_pupNABB, 2);
  else
  if (pars->NTExtraBBOn) SafeSetValue(Cn3D_pupNABB, 3);
  else SafeSetValue(Cn3D_pupNABB, 4);


  for (i = 0; i < 9; i++) {
  j = k = 0;
  switch (i)
    {
     case 0: /* prot bb */
         j = pars->PBBRender;
	   k = pars->PBBColor;
       break;
     case 1: /* na bb */
         j = pars->NTBBRender;
         k = pars->NTBBColor;
       break;
     case 2: /* prot sc */
         j = pars->PResRender;
         k = pars->PResColor;
       break;
     case 3: /* na sc */
         j = pars->NTResRender;
         k = pars->NTResColor;
       break;
     case 4: /* hets */
         j = pars->HetRender;
         k = pars->HetColor;
       break;
     case 5: /* ions */
         j = pars->IonRender;
         k = pars->IonColor;
       break;
     case 6: /* connections */
         j = pars->ConnectRender;
         k = pars->ConnectColor;
       break;
     case 7: /* solvent */
         j = pars->SolventRender;
         k = pars->SolventColor;
        break;
     case 8: /* object */
         j = pars->ObjectRender;
         k = pars->ObjectColor;
       break;
   }

  switch (j)
	 {
	   case R_SPACE:
	     SetEnumPopup (Cn3D_pupRenderStyle [i], renderAlist [i], (UIEnum) 5);
	    break;
	   case R_STICK:
	     SetEnumPopup (Cn3D_pupRenderStyle [i], renderAlist [i], (UIEnum) 2);
	    break;
	   case R_BALLNSTICK:
	     SetEnumPopup (Cn3D_pupRenderStyle [i], renderAlist [i], (UIEnum) 3);
	    break;
	   case R_THICKWIRE:
	     SetEnumPopup (Cn3D_pupRenderStyle [i], renderAlist [i], (UIEnum) 4);
	    break;
	   case R_DEFAULT:
	   case R_WIRE:
	     SetEnumPopup (Cn3D_pupRenderStyle [i], renderAlist [i], (UIEnum) 1);
	  }

   switch (k)
	  {
	    case C_BYCHAIN:
	      SetEnumPopup (Cn3D_pupColorStyle [i], colorAlist [i], (UIEnum) 1);
	      break;
	    case C_BYSSTRU:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 2);
	      break;
	    case C_BYDOMAIN:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 3);
	      break;
	    case C_BYRES:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 4);
	      break;
	    case C_BYHYDRO:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 5);
	      break;
	    case C_CPK:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 6);
	      break;
	    case C_BYTEMP:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 7);
	      break;
	    case C_BYOBJECT:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 8);
	      break;
		case C_BYALIGN:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 9);
	      break;
		case C_BYCONS:
	      SetEnumPopup (Cn3D_pupColorStyle [i],  colorAlist [i], (UIEnum) 10);
	      break;
	    default:
	      SetEnumPopup (Cn3D_pupColorStyle [i], colorAlist [i], (UIEnum) 1);
	  }
  }

}


static void MasterChangeLabelsProc (void)
{
  Int2 i, k, nameval;
  Uint1 codeval;
  Boolean NumOn, NameOn, PDBOn, WhiteOn, TopOn;
  PARS pars = NULL;
  PMSD pmsdThis = NULL;
  PDNMS pdnmsThis = NULL;


 /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) return;

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars) return;

  pars->PBBLabelOn = (Uint1) GetValue(Cn3D_pupLabelAA);
  pars->NTBBLabelOn = (Uint1) GetValue(Cn3D_pupLabelNT);

  /* set termini labels */
  pars->PTermLabelOn = (Uint1) GetStatus(Cn3D_lOnOffLabel[2]);
  pars->NTTermLabelOn = (Uint1) GetStatus(Cn3D_lOnOffLabel[3]);

  for (i = 0; i < 4; i++) {
  nameval = GetValue(Cn3D_bLName [i]);
  if (nameval == 3) codeval = (Uint1) L_1LETR;
  else codeval = (Uint1) L_3LETR;
  switch (i)
    {
     case 0: /* aa */
          /* clear bits */
         pars->PBBLabelStyle =  pars->PBBLabelStyle & ~((Uint1)( L_3LETR | L_1LETR)) ;
          /* set bit */
         pars->PBBLabelStyle = pars->PBBLabelStyle | (Uint1) codeval;
        break;
     case 1:  /* na   */
         pars->NTBBLabelStyle =  pars->NTBBLabelStyle & ~((Uint1)( L_3LETR | L_1LETR)) ;
         pars->NTBBLabelStyle = pars->NTBBLabelStyle | (Uint1) codeval;
        break;
     case 2: /* term prot*/
         pars->PTermLabelStyle =  pars->PTermLabelStyle & ~((Uint1)( L_3LETR | L_1LETR)) ;
         pars->PTermLabelStyle = pars->PTermLabelStyle  | (Uint1) codeval;
        break;
     case 3:  /* term nt   */
         pars->NTTermLabelStyle =  pars->NTTermLabelStyle& ~((Uint1)( L_3LETR | L_1LETR)) ;
         pars->NTTermLabelStyle = pars->NTTermLabelStyle  | (Uint1) codeval;
        break;
     default: ;
    }

  NumOn = GetStatus(Cn3D_bLNum [i]);
  switch (i)
    {
     case 0: /* aa */
         if (NumOn) pars->PBBLabelStyle =  pars->PBBLabelStyle | (Uint1) L_NUM;
         else pars->PBBLabelStyle = pars->PBBLabelStyle & (Uint1) ~L_NUM;
        break;
     case 1:  /* na   */
         if (NumOn) pars->NTBBLabelStyle =  pars->NTBBLabelStyle | (Uint1) L_NUM;
         else pars->NTBBLabelStyle = pars->NTBBLabelStyle & (Uint1) ~L_NUM;
        break;
     case 2: /* term prot*/
         if (NumOn) pars->PTermLabelStyle =  pars->PTermLabelStyle | (Uint1) L_NUM;
         else pars->PTermLabelStyle = pars->PTermLabelStyle & (Uint1) ~L_NUM;
        break;
     case 3:  /* term nt   */
          if (NumOn) pars->NTTermLabelStyle =  pars->NTTermLabelStyle | (Uint1) L_NUM;
         else pars->NTTermLabelStyle = pars->NTTermLabelStyle & (Uint1) ~L_NUM;
        break;
     default: ;
    }
  NameOn = (Boolean) (nameval > 1);
  switch (i)
    {
     case 0: /* aa */
         if (NameOn) pars->PBBLabelStyle =  pars->PBBLabelStyle | (Uint1) L_NAME;
         else pars->PBBLabelStyle = pars->PBBLabelStyle & (Uint1) ~L_NAME;
        break;
     case 1:  /* na   */
         if (NameOn) pars->NTBBLabelStyle =  pars->NTBBLabelStyle | (Uint1) L_NAME;
         else pars->NTBBLabelStyle = pars->NTBBLabelStyle & (Uint1) ~L_NAME;
        break;
     case 2: /* term prot*/
         if (NameOn) pars->PTermLabelStyle =  pars->PTermLabelStyle | (Uint1) L_NAME;
         else pars->PTermLabelStyle = pars->PTermLabelStyle & (Uint1) ~L_NAME;
        break;
     case 3:  /* term nt   */
          if (NameOn) pars->NTTermLabelStyle =  pars->NTTermLabelStyle | (Uint1) L_NAME;
         else pars->NTTermLabelStyle = pars->NTTermLabelStyle & (Uint1) ~L_NAME;
        break;
     default: ;
    }
  PDBOn = (Boolean) (nameval == 4);
  switch (i)
    {
     case 0: /* aa */
         if (PDBOn) pars->PBBLabelStyle =  pars->PBBLabelStyle | (Uint1) L_PDB;
         else pars->PBBLabelStyle = pars->PBBLabelStyle & (Uint1) ~L_PDB;
        break;
     case 1:  /* na   */
         if (PDBOn) pars->NTBBLabelStyle =  pars->NTBBLabelStyle | (Uint1) L_PDB;
         else pars->NTBBLabelStyle = pars->NTBBLabelStyle & (Uint1) ~L_PDB;
        break;
     case 2: /* term prot*/
         if (PDBOn) pars->PTermLabelStyle =  pars->PTermLabelStyle | (Uint1) L_PDB;
         else pars->PTermLabelStyle = pars->PTermLabelStyle & (Uint1) ~L_PDB;
        break;
     case 3:  /* term nt   */
          if (PDBOn) pars->NTTermLabelStyle =  pars->NTTermLabelStyle | (Uint1) L_PDB;
         else pars->NTTermLabelStyle = pars->NTTermLabelStyle & (Uint1) ~L_PDB;
        break;
     default: ;
    }
  WhiteOn = GetStatus(Cn3D_bLWhite [i]);
  switch (i)
    {
     case 0: /* aa */
         if (WhiteOn) pars->PBBLabelStyle =  pars->PBBLabelStyle | (Uint1) L_WHITE;
         else pars->PBBLabelStyle = pars->PBBLabelStyle & (Uint1) ~L_WHITE;
        break;
     case 1:  /* na   */
         if (WhiteOn) pars->NTBBLabelStyle =  pars->NTBBLabelStyle | (Uint1) L_WHITE;
         else pars->NTBBLabelStyle = pars->NTBBLabelStyle & (Uint1) ~L_WHITE;
        break;
     case 2: /* term prot*/
         if (WhiteOn) pars->PTermLabelStyle =  pars->PTermLabelStyle | (Uint1) L_WHITE;
         else pars->PTermLabelStyle = pars->PTermLabelStyle & (Uint1) ~L_WHITE;
        break;
     case 3:  /* term nt   */
          if (WhiteOn) pars->NTTermLabelStyle =  pars->NTTermLabelStyle | (Uint1) L_WHITE;
         else pars->NTTermLabelStyle = pars->NTTermLabelStyle & (Uint1) ~L_WHITE;
        break;
     default: ;
    }
  TopOn = GetStatus(Cn3D_bLTop [i]);
  switch (i)
    {
     case 0: /* aa */
         if (TopOn) pars->PBBLabelJust =  pars->PBBLabelJust | (Uint1) LA_FRONT;
         else pars->PBBLabelJust = pars->PBBLabelJust & (Uint1) ~LA_FRONT;
        break;
     case 1:  /* na   */
         if (TopOn) pars->NTBBLabelJust =  pars->NTBBLabelJust | (Uint1) LA_FRONT;
         else pars->NTBBLabelJust = pars->NTBBLabelJust & (Uint1) ~LA_FRONT;
        break;
     case 2: /* term prot*/
         if (TopOn) pars->PTermLabelJust =  pars->PTermLabelJust | (Uint1) LA_FRONT;
         else pars->PTermLabelJust = pars->PTermLabelJust & (Uint1) ~LA_FRONT;
        break;
     case 3:  /* term nt   */
          if (TopOn) pars->NTTermLabelJust =  pars->NTTermLabelJust | (Uint1) LA_FRONT;
         else pars->NTTermLabelJust = pars->NTTermLabelJust & (Uint1) ~LA_FRONT;
        break;
     default: ;
    }
  k = (Int2) GetValue(Cn3D_pupLabelSize [i]);
   switch (i)
    {
     case 0: /* prot bb */
         pars->PBBLabelScale = k;
         break;
     case 1: /* na bb */
         pars->NTBBLabelScale = k;
         break;
     case 2: /* prot term */
         pars->PTermLabelScale = k;
         break;
     case 3: /* na term */
         pars->NTTermLabelScale = k;
         break;
     default:
     ;
   }
  }


  if (GetStatus (Cn3D_LimmediateRedraw)) {
    Cn3D_RedrawProc (Cn3D_bLRedraw);
  }
}

static void LabelPopupOnOffProc (PopuP p)

{
  MasterChangeLabelsProc ();
}

static void LabelButtonOnOffProc (ButtoN b)

{
  MasterChangeLabelsProc ();
}

static void LabelGroupOnOffProc (GrouP g)

{
  MasterChangeLabelsProc ();
}

void LIBCALL ResetLabelCtrls(void)
{
  Int2 i, val, style, size;
  PARS pars = NULL;
  PMSD pmsdThis = NULL;
  PDNMS pdnmsThis = NULL;

 /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) goto setout;

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars)
   {
    setout:

      SafeSetValue(Cn3D_pupLabelAA, 1);
      SafeSetValue(Cn3D_pupLabelNT, 1);

      SafeSetStatus(Cn3D_lOnOffLabel [2], FALSE);
      SafeSetStatus(Cn3D_lOnOffLabel [3], FALSE);

      for (i = 0; i < 4; i++) {
        SafeSetValue(Cn3D_bLName [i], 1);
        SafeSetStatus(Cn3D_bLNum [i], FALSE);
        SafeSetStatus(Cn3D_bLTop [i], FALSE);
        SafeSetStatus(Cn3D_bLWhite [i], FALSE);
        SafeSetValue(Cn3D_pupLabelSize [i], 1);
      }

      return;
   }

  SafeSetValue(Cn3D_pupLabelAA, pars->PBBLabelOn);
  SafeSetValue(Cn3D_pupLabelNT, pars->NTBBLabelOn);

  SafeSetStatus(Cn3D_lOnOffLabel [2], pars->PTermLabelOn);
  SafeSetStatus(Cn3D_lOnOffLabel [3], pars->NTTermLabelOn);

  for (i = 0; i < 4; i++) {
    val = 1;
    style = 1;
    size = 1;
    switch (i) {
      case 0 :
        style = pars->PBBLabelStyle;
        size = pars->PBBLabelScale;
        break;
      case 1 :
        style = pars->NTBBLabelStyle;
        size = pars->NTBBLabelScale;
        break;
      case 2 :
        style = pars->PTermLabelStyle;
        size = pars->PTermLabelScale;
        break;
      case 3 :
        style = pars->NTTermLabelStyle;
        size = pars->NTTermLabelScale;
        break;
      default :
        break;
    }
    if (! (Boolean) (style & (Uint1) (L_NAME))) {
      val = 1;
    } else if ((Boolean) (style & (Uint1) (L_3LETR))) {
      val = 2;
    } else if ((Boolean) (style & (Uint1) (L_1LETR))) {
      val = 3;
    } else if ((Boolean) (style & (Uint1) L_PDB)) {
      val = 4;
    }
    SafeSetValue(Cn3D_bLName [i], val);
    SafeSetStatus(Cn3D_bLNum [i], (Boolean) (style & (Uint1) L_NUM));
    SafeSetStatus(Cn3D_bLTop [i], (Boolean) (style & (Uint1) LA_FRONT));
    SafeSetStatus(Cn3D_bLWhite [i], (Boolean) (style & (Uint1) ~L_WHITE));
    SafeSetValue(Cn3D_pupLabelSize [i], size);
  }
}


void Cn3D_RedrawProc(ButtoN b)
{
  Int2 k;
  PARS pars = NULL;
  PMSD pmsdThis = NULL;
  PDNMS pdnmsThis = NULL;
  PDNML pdnmlThis = NULL;
  PMLD  pmldThis = NULL;

/* yanli added */
   /* fetch the active structure */
  pdnmsThis = GetSelectedModelstruc();

  if (!pdnmsThis) goto setout;

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

  pars = GetAlgorRenderSet(pdnmsThis);
  if (!pars)
   {
    setout:
      for (k = 2; k < 10; k++) {
        SafeSetStatus (Cn3D_lOnOffItem [k], FALSE);
      }

      SafeSetValue(Cn3D_pupPBB, 4);
      SafeSetValue(Cn3D_pupNABB, 4);
      for (k = 0; k < 9; k++) {
        SetEnumPopup (Cn3D_pupRenderStyle [k], renderAlist [k], (UIEnum) renderDefaults [k]);
        SetEnumPopup (Cn3D_pupColorStyle [k], colorAlist [k], (UIEnum) colorDefaults [k]);
        /*
        SafeSetValue(Cn3D_pupRenderStyle [k], 1);
        SafeSetValue(Cn3D_pupColorStyle [k], 1);
        */
      }

      return;
   }

/* alpha C trace only control */
  pdnmlThis = pmsdThis->pdnmlModels;
  pmldThis = (PMLD)pdnmlThis->data.ptrvalue;
  if(pmldThis->iType == Model_type_ncbi_backbone) {
     pars->PVirtualBBOn = TRUE;
     pars->PRealBBOn = FALSE;
     pars->PExtraBBOn = FALSE;

     pars->NTVirtualBBOn = TRUE;
     pars->NTRealBBOn = FALSE;
     pars->NTExtraBBOn = FALSE;

     SafeSetValue(Cn3D_pupPBB, 1);
     SafeSetValue(Cn3D_pupNABB, 1);

     sfpThis = sfp;
     while(sfpThis){
        sfipThis = sfpThis->data.ptrvalue;
        pars = sfipThis->parsSpecial;

        pars->PVirtualBBOn = TRUE;
        pars->PRealBBOn = FALSE;
        pars->PExtraBBOn = FALSE;

        pars->NTVirtualBBOn = TRUE;
        pars->NTRealBBOn = FALSE;
        pars->NTExtraBBOn = FALSE;

/*      SafeSetValue(Cn3D_pupModelPBB, 1);
        SafeSetValue(Cn3D_pupModelNABB, 1);  */

        sfpThis = sfpThis->next;
     }
  }

  Cn3D_Redraw(FALSE);
  Cn3dObjMgrGetSelected();
}


static void AlignLabelButtons (Handle master, Handle slave)

{
  RecT pptPos, btnPos;
  Int2 delta;

  GetPosition (master, &pptPos);
  GetPosition (slave, &btnPos);
  delta = (pptPos.bottom + pptPos.top) / 2 - (btnPos.bottom + btnPos.top) / 2;
  if (delta != 0) {
    OffsetRect (&btnPos, 0, delta);
    SetPosition (slave, &btnPos);
    AdjustPrnt ((Nlm_GraphiC) slave, &btnPos, FALSE);
  }
}

GrouP LIBCALL LabelControls ( Nlm_GrouP prnt)
{
  GrouP g, h;
  Int2  k;
  PrompT ppt [4];
  PrompT ppt1, ppt2, ppt3, ppt4, ppt7, ppt8, ppt9;

  g = HiddenGroup ( prnt, -1, 0, NULL );
  StaticPrompt(g, "", 0, stdLineHeight, systemFont, 'l');
  StaticPrompt(g, "", 0, stdLineHeight, systemFont, 'l');
  if (!g) return NULL;
  SetGroupSpacing ( g, 3, 9 );
/*
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif
*/

  h = HiddenGroup (g, 0, -5, NULL);
  SetGroupSpacing (h, 5, 3);

  Cn3D_lOnOffLabel [0] = NULL; /* Cn3D_pupLabelAA includes off setting */
  Cn3D_lOnOffLabel [1] = NULL; /* Cn3D_pupLabelNT includes off setting */

  ppt1 = StaticPrompt (h, "Category", 0, stdLineHeight + 5, systemFont, 'c');
  ppt [0] = StaticPrompt (h, "Amino Acid", 0, popupMenuHeight, programFont, 'l');
  ppt [1] = StaticPrompt (h, "Nucleic Acid", 0, popupMenuHeight, programFont, 'l');
  ppt [2] = StaticPrompt (h, "N & C Termini", 0, popupMenuHeight, programFont, 'l');
  ppt [3] = StaticPrompt (h, "5' & 3' Termini", 0, popupMenuHeight, programFont, 'l');
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt1, (HANDLE) ppt [3], NULL);

  ppt2 = StaticPrompt (h, "On/Off", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_pupLabelAA = PopupList(h, FALSE, LabelPopupOnOffProc);
  PopupItem(Cn3D_pupLabelAA, "none");
  PopupItem(Cn3D_pupLabelAA, "every AA");
  PopupItem(Cn3D_pupLabelAA, "every 5");
  PopupItem(Cn3D_pupLabelAA, "every 10");
  PopupItem(Cn3D_pupLabelAA, "every 20");
  PopupItem(Cn3D_pupLabelAA, "every 50");
  Cn3D_pupLabelNT= PopupList(h, FALSE, LabelPopupOnOffProc);
  PopupItem(Cn3D_pupLabelNT, "none");
  PopupItem(Cn3D_pupLabelNT, "every NA");
  PopupItem(Cn3D_pupLabelNT, "every 5");
  PopupItem(Cn3D_pupLabelNT, "every 10");
  PopupItem(Cn3D_pupLabelNT, "every 20");
  PopupItem(Cn3D_pupLabelNT, "every 50");
  Cn3D_lOnOffLabel [2] = CheckBox (h ,"", LabelButtonOnOffProc);
  Cn3D_lOnOffLabel [3] = CheckBox (h ,"", LabelButtonOnOffProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt2, (HANDLE) Cn3D_pupLabelAA, (HANDLE) Cn3D_pupLabelNT,
                (HANDLE) Cn3D_lOnOffLabel [2], (HANDLE) Cn3D_lOnOffLabel [3], NULL);

  ppt3 = StaticPrompt (h, "Name", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 4; k++) {
    Cn3D_bLName [k] = PopupList (h, FALSE, LabelPopupOnOffProc);
    PopupItem (Cn3D_bLName [k], "none");
    PopupItem (Cn3D_bLName [k], "3 letter");
    PopupItem (Cn3D_bLName [k], "1 letter");
    PopupItem (Cn3D_bLName [k], "PDB");
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt3, (HANDLE) Cn3D_bLName [0], (HANDLE) Cn3D_bLName [1],
                (HANDLE) Cn3D_bLName [2], (HANDLE) Cn3D_bLName [3], NULL);

  ppt4 = StaticPrompt (h, "Num", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 4; k++) {
    Cn3D_bLNum [k] = CheckBox(h, "", LabelButtonOnOffProc );
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt4, (HANDLE) Cn3D_bLNum [0], (HANDLE) Cn3D_bLNum [1],
                (HANDLE) Cn3D_bLNum [2], (HANDLE) Cn3D_bLNum [3], NULL);

  ppt7 = StaticPrompt (h, "On Top", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 4; k++) {
    Cn3D_bLTop [k] = CheckBox(h, "", LabelButtonOnOffProc );
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt7, (HANDLE) Cn3D_bLTop [0], (HANDLE) Cn3D_bLTop [1],
                (HANDLE) Cn3D_bLTop [2], (HANDLE) Cn3D_bLTop [3], NULL);

  ppt8 = StaticPrompt (h, "White", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 4; k++) {
    Cn3D_bLWhite [k] = CheckBox(h, "",  LabelButtonOnOffProc  );
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt8, (HANDLE) Cn3D_bLWhite [0], (HANDLE) Cn3D_bLWhite [1],
                (HANDLE) Cn3D_bLWhite [2], (HANDLE) Cn3D_bLWhite [3], NULL);

  ppt9 = StaticPrompt (h, "Size", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 4; k++) {
    Cn3D_pupLabelSize [k] = PopupList(h, FALSE, LabelPopupOnOffProc);
    PopupItem(Cn3D_pupLabelSize [k], "1");
    PopupItem(Cn3D_pupLabelSize [k], "2");
    PopupItem(Cn3D_pupLabelSize [k], "3");
    PopupItem(Cn3D_pupLabelSize [k], "4");
    PopupItem(Cn3D_pupLabelSize [k], "5");
    PopupItem(Cn3D_pupLabelSize [k], "6");
    PopupItem(Cn3D_pupLabelSize [k], "7");
    PopupItem(Cn3D_pupLabelSize [k], "8");
    PopupItem(Cn3D_pupLabelSize [k], "9");
    PopupItem(Cn3D_pupLabelSize [k], "10");
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt9, (HANDLE) Cn3D_pupLabelSize [0], (HANDLE) Cn3D_pupLabelSize [1],
                (HANDLE) Cn3D_pupLabelSize [2], (HANDLE) Cn3D_pupLabelSize [3], NULL);


  Cn3D_LimmediateRedraw = CheckBox (g, "Immediate Redraw on Every Change", NULL);

  Cn3D_bLRedraw = PushButton(g, "Apply!", Cn3D_RedrawProc);

  AlignLabelButtons (Cn3D_pupLabelSize [2], Cn3D_lOnOffLabel [2]);
  AlignLabelButtons (Cn3D_pupLabelSize [3], Cn3D_lOnOffLabel [3]);
  for (k = 0; k < 4; k++) {
    AlignLabelButtons (Cn3D_pupLabelSize [k], Cn3D_bLName [k]);
    AlignLabelButtons (Cn3D_pupLabelSize [k], Cn3D_bLNum [k]);
    AlignLabelButtons (Cn3D_pupLabelSize [k], Cn3D_bLTop [k]);
    AlignLabelButtons (Cn3D_pupLabelSize [k], Cn3D_bLWhite [k]);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) Cn3D_bLRedraw,
                (HANDLE) Cn3D_LimmediateRedraw, (HANDLE) h, NULL);

  ResetLabelCtrls();
  return g;
}

static void fnCn3D_RedrawWrapper(ButtoN b)
{
  Cn3D_ReColor = TRUE;
  Cn3D_RedrawProc(b);
  Cn3D_ReColor = FALSE;
}

static CharPtr categorynames [] = {
  "", "", "Protein Sidechains", "Nucleotide Sidechains", "Heterogen",
  "Ions", "Connections", "Solvents", "3D Objects", "Hydrogens"
};

static PopuP RenderStyle (GrouP h, Int2 i)

{
  PopuP  p;

  if (renderAlist [i] != empty_alist) {
    p = PopupList (h, FALSE, PopupOnOffProc);
  } else {
    p = PopupList (h, FALSE, NULL);
  }
  InitEnumPopup (p, renderAlist [i], NULL);
  /*
  PopupItem (p, "Wireframe");
  PopupItem (p, "Tubes");
  PopupItem (p, "Ball&Stick");
  PopupItem (p, "Fat Tubes");
  PopupItem (p, "SpaceFill");
  */
  return p;
}

static PopuP ColorStyle (GrouP h, Int2 i)

{
  PopuP  p;

  if (colorAlist [i] != empty_alist) {
    p = PopupList (h, FALSE, PopupOnOffProc);
  } else {
    p = PopupList (h, FALSE, NULL);
  }
  InitEnumPopup (p, colorAlist [i], NULL);
  /*
  PopupItem (p, "Molecule Cycle");
  PopupItem (p, "Sec. Structure");
  PopupItem (p, "Domain");
  PopupItem (p, "Residue");
  PopupItem (p, "Hydrophobicity");
  PopupItem (p, "CPK");
  PopupItem (p, "Temperature");
  PopupItem (p, "Object Cycle");
  PopupItem (p, "Structure");
  PopupItem (p, "Neighbor");
  */
  return p;
}

GrouP LIBCALL  RenderControls ( Nlm_GrouP prnt)
{
  GrouP g, h;
  Int2  k;
  PrompT ppt [10];
  PrompT ppt1, ppt2, ppt3, ppt4;
  RecT pptPos, btnPos;
  Int2 delta;

  g = HiddenGroup ( prnt, -1, 0, NULL );
  if (!g) return NULL;
  SetGroupSpacing ( g, 3, 9 );
/*
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif
*/

  h = HiddenGroup (g, 0, -11, NULL);
  SetGroupSpacing (h, 5, 3);

  Cn3D_lOnOffItem [0] = NULL; /* Cn3D_pupPBB includes off setting */
  Cn3D_lOnOffItem [1] = NULL; /* Cn3D_pupNABB includes off setting */

  ppt1 = StaticPrompt (h, "Category", 0, stdLineHeight + 5, systemFont, 'c');
  ppt [0] = StaticPrompt (h, "Protein Backbone", 0, popupMenuHeight, programFont, 'l');
  ppt [1] = StaticPrompt (h, "Nucleotide Backbone", 0, popupMenuHeight, programFont, 'l');
  for (k = 2; k < 10; k++) {
    ppt [k] = StaticPrompt (h, categorynames [k], 0, popupMenuHeight, programFont, 'l');
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt1, (HANDLE) ppt [1], NULL);

  ppt2 = StaticPrompt (h, "On/Off", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_pupPBB = PopupList(h, FALSE, PopupOnOffProc);
  PopupItem(Cn3D_pupPBB, "alpha C trace");
  PopupItem(Cn3D_pupPBB, "partial atom");
  PopupItem(Cn3D_pupPBB, "all atoms");
  PopupItem(Cn3D_pupPBB, "none");
  Cn3D_pupNABB = PopupList(h, FALSE, PopupOnOffProc);
  PopupItem(Cn3D_pupNABB, "P trace");
  PopupItem(Cn3D_pupNABB, "partial atom");
  PopupItem(Cn3D_pupNABB, "all atoms");
  PopupItem(Cn3D_pupNABB, "none");
  for (k = 2; k < 10; k++) {
    Cn3D_lOnOffItem [k] = CheckBox (h, "", ButtonOnOffProc);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt2, (HANDLE) Cn3D_pupPBB, NULL);

  ppt3 = StaticPrompt (h, "Render", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_pupRenderStyle [0] = RenderStyle (h, 0);
  Cn3D_pupRenderStyle [1] = RenderStyle (h, 1);
  for (k = 2; k < 10; k++) {
    if (k != 8 && k != 9) {
      Cn3D_pupRenderStyle [k] = RenderStyle (h, k);
    } else {
      StaticPrompt (h, " ", 0, popupMenuHeight, programFont, 'l');
      Cn3D_pupRenderStyle [k] = NULL;
    }
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt3, (HANDLE) Cn3D_pupRenderStyle [0], NULL);

  ppt4 = StaticPrompt (h, "Color", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_pupColorStyle [0] = ColorStyle (h, 0);
  Cn3D_pupColorStyle [1] = ColorStyle (h, 1);
  for (k = 2; k < 10; k++) {
    if (k != 6 && k != 9) {
      Cn3D_pupColorStyle [k] = ColorStyle (h, k);
    } else {
      StaticPrompt (h, " ", 0, popupMenuHeight, programFont, 'l');
      Cn3D_pupColorStyle [k] = NULL;
    }
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt4, (HANDLE) Cn3D_pupColorStyle [0], NULL);

  Cn3D_pupRenderStyle [9] = NULL; /* not for hydrogens */
  Cn3D_pupColorStyle [9] = NULL; /* not for hydrogens */

  Cn3D_immediateRedraw = CheckBox (g, "Immediate Redraw on Every Change", NULL);

  Cn3D_bRedraw = PushButton (g, "Apply!", fnCn3D_RedrawWrapper);

  for (k = 2; k < 10; k++) {
    GetPosition (Cn3D_lOnOffItem [k], &btnPos);
    GetPosition (ppt [k], &pptPos);
    delta = (pptPos.bottom + pptPos.top) / 2 - (btnPos.bottom + btnPos.top) / 2;
    if (delta != 0) {
      OffsetRect (&btnPos, 0, delta);
      SetPosition (Cn3D_lOnOffItem [k], &btnPos);
      AdjustPrnt ((Nlm_GraphiC) Cn3D_lOnOffItem [k], &btnPos, FALSE);
    }
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) Cn3D_pupPBB, (HANDLE) Cn3D_pupNABB,
                (HANDLE) Cn3D_lOnOffItem [2], (HANDLE) Cn3D_lOnOffItem [3],
                (HANDLE) Cn3D_lOnOffItem [4], (HANDLE) Cn3D_lOnOffItem [5],
                (HANDLE) Cn3D_lOnOffItem [6], (HANDLE) Cn3D_lOnOffItem [7],
                (HANDLE) Cn3D_lOnOffItem [8], (HANDLE) Cn3D_lOnOffItem [9],
                NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) Cn3D_bRedraw,
                (HANDLE) Cn3D_immediateRedraw, (HANDLE) h, NULL);

  ResetRenderCtrls();
  return g;
}


/****************************************************************************/




PRK LIBCALL NewRenderKeep(void)
{
  PRK prkThis = NULL;
  prkThis = (PRK) MemNew((size_t)(sizeof(RK)));
  if (!prkThis) return NULL;
  prkThis->NodeWhat = (Byte) CONVERT_ALL;
  prkThis->NodeType = (Byte) CONVERT_ALL;
  /*prkThis->Color = (Int2) 0; /* a fixed color */
  prkThis->Bond = (Byte)  NO_BOND;  /*  use define */
  prkThis->Atom = (Byte)  ATOM_NONE;
  prkThis->BondWidth = (float)0;
  prkThis->LJust =  0;
  prkThis->LStyle = 0;
  prkThis->LScale = 1;
  return prkThis;
}

PRK LIBCALL CopyRenderKeep(PRK prkThis)
{
    PRK prkNew = NULL;
    prkNew = NewRenderKeep();
    if (!prkNew) return NULL;
    prkNew->NodeWhat = prkThis->NodeWhat;
    prkNew->NodeType = prkThis->NodeType;
    Cn3D_CopyColorCell(&(prkNew->Color), &(prkThis->Color)); /* a fixed color */
    prkNew->Bond = prkThis->Bond;  /*  use define */
    prkNew->Atom = prkThis->Atom;
    prkNew->BondWidth = prkThis->BondWidth;
    prkNew->LJust = prkThis->LJust;
    prkNew->LStyle = prkThis->LStyle;
    prkNew->LScale = prkThis->LScale;
    return prkNew;
}

void LIBCALL FreeRenderKeep(PRK prkThis)
{
  MemFree(prkThis);
}

static void RotTransScale(PMSD pmsdThis,   FloatLo fX,
					   FloatLo fY,
					   FloatLo fZ,
					   Int4Ptr piX,
					   Int4Ptr piY,
					   Int4Ptr piZ)
{
    FloatHi fXTemp,  fYTemp,  fZTemp;

    if (!pmsdThis) return;

    *piX = 0;
    *piY = 0;
    *piZ = 0;
    if (pmsdThis->pflTranslate)
	/* do translation */
	{
	  fX = fX - pmsdThis->pflTranslate[0];
	  fY = fY - pmsdThis->pflTranslate[1];
	  fZ = fZ - pmsdThis->pflTranslate[2];
        }
    if (pmsdThis->ppflRotate)
	/* do rotation */
	{
	  fXTemp = (FloatHi) (fX * pmsdThis->ppflRotate[0][0] +
			      fY * pmsdThis->ppflRotate[1][0] +
			      fZ * pmsdThis->ppflRotate[2][0]);
	  fYTemp= (FloatHi) ( fX * pmsdThis->ppflRotate[0][1] +
			      fY * pmsdThis->ppflRotate[1][1] +
			      fZ * pmsdThis->ppflRotate[2][1]);
	  fZTemp = (FloatHi)( fX * pmsdThis->ppflRotate[0][2] +
			      fY * pmsdThis->ppflRotate[1][2] +
			      fZ * pmsdThis->ppflRotate[2][2]);
          fX = (FloatLo)fXTemp;
          fY = (FloatLo)fYTemp;
          fZ = (FloatLo)fZTemp;
	}
	     /* scale */

     *piX = (Int4) (fX  *VIEWSCALE);
     *piY = (Int4) (fY  *VIEWSCALE);
     *piZ = (Int4) (fZ  *VIEWSCALE);
   return;
}


/***********************************************************/
/* This is the lowest-level call to the Viewer3D library   */
/* at this point all the decisions have been made          */
/* this performs any necessary model-space transformations */
/* if there are any rotation-translation matrices present  */
/* It may be necessary later to instantiate a linked-list  */
/* of such transformations for moving structures relative  */
/* to each other */


void LIBCALL RenderObject(PVNMO pvnmoThis)
{
    PMOD pmodThis = NULL;
    PMSD pmsdThis = NULL;
    PMGD pmgdThis = NULL;
    PMMD pmmdThis = NULL;
    PFB pfbThis = NULL;
    PFB pfbParent = NULL;
    PARS pars = NULL;
    Int4 i = 0;
    FloatLo fXFrom,  fYFrom,  fZFrom;
    FloatLo fXTo,  fYTo,  fZTo;
    Int4  iXFrom, iYFrom, iZFrom;
    Int4  iXTo, iYTo, iZTo;
    Int4 iCylRadius;
#ifndef _OPENGL
    Nlm_Prim3D poly;
#endif
    ResidueColorCell iColor;

#ifdef DEBUG_N
	  printf("RenderObject\n" );
#endif
    if (!pvnmoThis) return;
    pmodThis = (PMOD) pvnmoThis->data.ptrvalue;
    if (!pmodThis) return;

    if (!pmodThis->ppflObject) return;
    if(!pmodThis->bVisible) return;

    pmsdThis = ToMSDParent((PFB)pmodThis);
    if (!pmsdThis) return;
    
    pars = GetAlgorRenderSet(pmsdThis->pdnmsLink);
    if (!pars) return;  /* cannot render */
    
    Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[pars->ObjectColor]));
    if ( Cn3d_ColorNow == C_BYOBJECT)
    {
      Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ColorNumKinBB[(pmodThis->pmldCoordSet->iNoCoordSet % KIN_COLOR_NUM)]]));
    }
    if ( Cn3d_ColorNow == C_BYCHAIN)
    {
      if (pmodThis->pvnContains == NULL)
        Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ColorNumKinBB[0]]));
      else
      {
        pfbThis = (PFB) pmodThis->pvnContains->data.ptrvalue;
        if (!pfbThis) return;
        if (!IsMoleculeNode(pfbThis))
          pfbParent = (PFB) GetParentMol(pfbThis);
        else
          pfbParent = pfbThis;
        if (pfbParent)
        {
          pmmdThis = (PMMD) pfbParent;
          Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ColorNumKinBB[(pmmdThis->pdnmmLink->choice % KIN_COLOR_NUM)]]));
        }
        else
          Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ColorNumKinBB[0]]));
      }
    }
    if (Cn3d_ColorNow == C_BYSSTRU   || Cn3d_ColorNow == C_BYDOMAIN)
      if (pmodThis->pvnContains == NULL)
        Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ColorNumKinBB[0]]));
      else
      {
        pfbThis = (PFB) pmodThis->pvnContains->data.ptrvalue;
        if (!pfbThis) return;
        if (!IsGraphNode(pfbThis))
          pfbParent = (PFB) GetParentGraph(pfbThis);
        else
          pfbParent = pfbThis;
        if (pfbParent)
          if (IsGraphAminoAcid(pfbParent))
          {
            pmgdThis = (PMGD) pfbParent;
            if (Cn3d_ColorNow == C_BYSSTRU)
               Cn3D_CopyColorCell(&iColor, ColorFromSS(pmgdThis));
            if (Cn3d_ColorNow == C_BYDOMAIN)
              Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ColorNumKinBB[(pmgdThis->iDomain % KIN_COLOR_NUM)]]));
          }
          else
            Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ColorNumKinBB[0]]));
          else
            Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ColorNumKinBB[0]]));
      }
      if (Cn3d_ColorNow == C_BYALIGN)
      {
        Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[bColorAlignments[(Cn3d_lSlaveNum % NUM_SLAVES)]]));
      }
      if (Cn3d_ColorPass)  { Cn3D_RequestColor(&Cn3D_ColorData, &iColor); return; }
      

    switch (pmodThis->bWhat)
      {
	   case OBJ_CYLINDER:
#ifdef _OPENGL
	      fXFrom = (FloatLo) (pmodThis->ppflObject[0][0]);
	      fYFrom = (FloatLo) (pmodThis->ppflObject[0][1]);
	      fZFrom = (FloatLo) (pmodThis->ppflObject[0][2]);
	      fXTo = (FloatLo) (pmodThis->ppflObject[1][0]);
	      fYTo = (FloatLo) (pmodThis->ppflObject[1][1]);
	      fZTo = (FloatLo) (pmodThis->ppflObject[1][2]);
	      RotTransScale( pmsdThis,     fXFrom,
					   fYFrom,
					   fZFrom,
					   &iXFrom,
					   &iYFrom,
					   &iZFrom);
          fXFrom = iXFrom/VIEWSCALE;
          fYFrom = iYFrom/VIEWSCALE;
          fZFrom = iZFrom/VIEWSCALE;

	      RotTransScale( pmsdThis,     fXTo,
					   fYTo,
					   fZTo,
					   &iXTo,
					   &iYTo,
					   &iZTo);
          fXTo = iXTo/VIEWSCALE;
          fYTo = iYTo/VIEWSCALE;
          fZTo = iZTo/VIEWSCALE;

OGL_AddCylinder3D(Cn3D_ColorData.OGL_Data, &iColor, (FloatHi)fXFrom, (FloatHi)fYFrom, (FloatHi)fZFrom,
        (FloatHi)fXTo,(FloatHi)fYTo,(FloatHi)fZTo, (FloatHi)(pmodThis->flRadius));
#else /* _OPENGL */

	      fXFrom = (FloatLo) (pmodThis->ppflObject[0][0]);
	      fYFrom = (FloatLo) (pmodThis->ppflObject[0][1]);
	      fZFrom = (FloatLo) (pmodThis->ppflObject[0][2]);
	      fXTo = (FloatLo) (pmodThis->ppflObject[1][0]);
	      fYTo = (FloatLo) (pmodThis->ppflObject[1][1]);
	      fZTo = (FloatLo) (pmodThis->ppflObject[1][2]);
	      RotTransScale( pmsdThis,     fXFrom,
					   fYFrom,
					   fZFrom,
					   &iXFrom,
					   &iYFrom,
					   &iZFrom);
	      RotTransScale( pmsdThis,     fXTo,
					   fYTo,
					   fZTo,
					   &iXTo,
					   &iYTo,
					   &iZTo);
	      iCylRadius = (Int4) (pmodThis->flRadius*VIEWSCALE);
	      AddCylinder3D (pic, NULL, (BigScalar) pmodThis,
				      Cn3d_LayerNow ,
				      Cn3D_ColorIndex(&Cn3D_ColorData, &iColor),
				      iXFrom,iYFrom,iZFrom,
                                      iXTo,iYTo,iZTo,
				      iCylRadius  );
#endif /* else _OPENGL */
	     return;
	    case  OBJ_BRICK:
#ifdef _OPENGL
      {
          double v1[3], v2[3], v3[3], v4[3];
          Int4 xv, yv, zv;
	     RotTransScale( pmsdThis,
                       pmodThis->ppflObject[0][0],
					   pmodThis->ppflObject[0][1],
					   pmodThis->ppflObject[0][2],
					   &xv,
					   &yv,
					   &zv);
          v1[0] = xv/1000000.0;
          v1[1] = yv/1000000.0;
          v1[2] = zv/1000000.0;
	     RotTransScale( pmsdThis,
                       pmodThis->ppflObject[1][0],
					   pmodThis->ppflObject[1][1],
					   pmodThis->ppflObject[1][2],
					   &xv,
					   &yv,
					   &zv);
          v2[0] = xv/1000000.0;
          v2[1] = yv/1000000.0;
          v2[2] = zv/1000000.0;
	     RotTransScale( pmsdThis,
                       pmodThis->ppflObject[5][0],
					   pmodThis->ppflObject[5][1],
					   pmodThis->ppflObject[5][2],
					   &xv,
					   &yv,
					   &zv);
          v3[0] = xv/1000000.0;
          v3[1] = yv/1000000.0;
          v3[2] = zv/1000000.0;
	     RotTransScale( pmsdThis,
                       pmodThis->ppflObject[4][0],
					   pmodThis->ppflObject[4][1],
					   pmodThis->ppflObject[4][2],
					   &xv,
					   &yv,
					   &zv);
          v4[0] = xv/1000000.0;
          v4[1] = yv/1000000.0;
          v4[2] = zv/1000000.0;

    OGL_AddQuad3D(Cn3D_ColorData.OGL_Data, &iColor, v1, v2, v3, v4);

      }
#else /* _OPENGL */
	     fXFrom = pmodThis->ppflObject[0][0];
	     fYFrom = pmodThis->ppflObject[0][1];
	     fZFrom = pmodThis->ppflObject[0][2];
	     fXTo = pmodThis->ppflObject[1][0];
	     fYTo = pmodThis->ppflObject[1][1];
	     fZTo = pmodThis->ppflObject[1][2];
	     RotTransScale( pmsdThis,      fXFrom,
					   fYFrom,
					   fZFrom,
					   &iXFrom,
					   &iYFrom,
					   &iZFrom);
	     RotTransScale( pmsdThis,      fXTo,
					   fYTo,
					   fZTo,
					   &iXTo,
					   &iYTo,
					   &iZTo);
             Cn3d_AnyPrim = TRUE;
	     poly = AddPoly3D (pic, NULL, (BigScalar) pmodThis,
	               Cn3d_LayerNow, Cn3D_ColorIndex(&Cn3D_ColorData, &iColor),
				      iXFrom,iYFrom,iZFrom,
                                      iXTo,iYTo,iZTo);
#ifdef DEBUG_N
	  printf("%ld %ld %ld\n",  (long) iXFrom,  (long)iYFrom, (long) iZFrom);
	  printf("%ld %ld %ld\n",  (long) iXTo,  (long)iYTo, (long) iZTo);
#endif
      for (i=2; i<8; i++)
	      {
		fXTo =  pmodThis->ppflObject[i][0];
		fYTo =  pmodThis->ppflObject[i][1];
		fZTo =  pmodThis->ppflObject[i][2];
		RotTransScale( pmsdThis,   fXTo,
					   fYTo,
					   fZTo,
					   &iXTo,
					   &iYTo,
					   &iZTo);
                Cn3d_AnyPrim = TRUE;
		AddVertPoly3D (pic, poly, iXTo,iYTo,iZTo);
#ifdef DEBUG_N
	  printf("%ld %ld %ld\n",  (long) iXTo, (long) iYTo, (long) iZTo);
#endif
	      }
#endif /* else _OPENGL */

#ifdef DEBUG_N
	  printf("%ld %ld %ld\n",  (long) iXFrom, (long) iYFrom, (long) iZFrom);
#endif
	     return;
	  case  OBJ_SPHERE:
	  case OBJ_CONE:
	  case OBJ_TMESH:
          case OBJ_TRIANGLES:
	  default:
#ifdef DEBUG_N
	  printf("Other  \n");
#endif
          return;
       }
}



  /* Mode determines how to combine PALD's */
#define RL_CENTER      1  /* put at paldCenter */
#define RL_CENTERPLUSY 2  /* put at paldCenter plus 1x Scale in y direction */
#define RL_EXTRAPOL    3  /* paldTo - paldFrom + paldTo */
#define RL_BETACARB    4  /* paldCenter + ((paldFrom - paldCenter)/2)
                                  + (((paldTo -paldCenter)/2) */

static void RenderLabel(PDNMS pdnmsThis, CharPtr pclabel,
			  PALD paldCenter, PALD paldFrom,  PALD paldTo,
			  ResidueColorCell * iColor, Uint1 Just, Int2 Scale, Int2 Mode)
{
   FloatLo fXFrom,  fYFrom,  fZFrom;
   Int4  iX, iY, iZ;
   FloatLo fXTo,  fYTo,  fZTo;
   FloatLo fXCen, fYCen, fZCen;
   FloatLo fXt,  fYt,  fZt;
   PMSD pmsdThis = NULL;
   PALD paldPrim = NULL;
   Int2 flags;


  if (!pdnmsThis) return;
  if (!pclabel) return;
  if (((Mode == RL_CENTER) || (Mode == RL_BETACARB)) && (paldCenter == NULL)) return;
  if ((Mode == RL_CENTERPLUSY) && (paldCenter == NULL)) return;
  if (((Mode == RL_EXTRAPOL) || (Mode == RL_BETACARB)) && ((paldFrom == NULL) || (paldTo == NULL))) return;

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

   if (paldFrom)
     {
      fXFrom = (FloatLo) (paldFrom->pflvData[0]);
      fYFrom = (FloatLo) (paldFrom->pflvData[1]);
      fZFrom = (FloatLo) (paldFrom->pflvData[2]);
     }
   if (paldTo)
     {
      fXTo = (FloatLo) (paldTo->pflvData[0]);
      fYTo = (FloatLo) (paldTo->pflvData[1]);
      fZTo = (FloatLo) (paldTo->pflvData[2]);
      paldPrim = paldTo;

     }
   if (paldCenter)
     {
      fXCen = (FloatLo) (paldCenter->pflvData[0]);
      fYCen = (FloatLo) (paldCenter->pflvData[1]);
      fZCen = (FloatLo) (paldCenter->pflvData[2]);
      paldPrim = paldCenter;

     }


  switch (Mode)
   {
    case RL_CENTER:
      break;
    case RL_CENTERPLUSY:
         fYCen = fYCen + Scale;
      break;
    case RL_EXTRAPOL:
         fXCen = fXTo - fXFrom + fXTo;
         fYCen = fYTo - fYFrom + fYTo;
         fZCen = fZTo - fZFrom + fZTo;
      break;
    case RL_BETACARB:
         fXt = fXCen - ((fXFrom - fXCen)/2) - ((fXTo - fXCen)/2);
         fYt = fYCen - ((fYFrom - fYCen)/2) - ((fYTo - fYCen)/2);
         fZt = fZCen - ((fZFrom - fZCen)/2) - ((fZTo - fZCen)/2);
	 fXCen = fXt;
	 fYCen = fYt;
	 fZCen = fZt;
      break;

    }

   RotTransScale(pmsdThis,  fXCen,
			    fYCen,
			    fZCen,
			     &iX,
			     &iY,
			     &iZ);

#ifdef _OPENGL
   flags = 0;
   if (Just & LA_CENTER) flags =  OGLTEXT3D_MIDDLE | OGLTEXT3D_CENTER;
   /* if (Just & LA_FRONT) flags = flags | OGLTEXT3D_FRONT;  hard to do in OpenGL
   todo: get rid of the ui for this
   */

   OGL_AddText3D(Cn3D_ColorData.OGL_Data, iColor, pclabel, (FloatHi)fXCen, (FloatHi)fYCen, (FloatHi)fZCen, flags);
#else /* _OPENGL */
   flags = 0;
   if (Just & LA_LEFT) flags = flags | TEXT3D_LEFT;
   if (Just & LA_RIGHT) flags = flags | TEXT3D_RIGTH;
   if (Just & LA_UPPER) flags = flags | TEXT3D_UPPER;
   if (Just & LA_LOWER) flags = flags | TEXT3D_LOWER;
   /* only LA_CENTER and LA_FRONT are used */
   if (Just & LA_CENTER) flags =  TEXT3D_MIDDLE | TEXT3D_CENTER;
   if (Just & LA_FRONT) flags = flags | TEXT3D_FRONT;
   Cn3d_AnyPrim = TRUE;
   if (paldPrim->pGraphic == NULL)
       paldPrim->pGraphic = (Pointer) AddSegment3D(pic, NULL, (BigScalar) paldPrim, Cn3d_LayerNow);

   AddText3D (pic,  (Segment3D) paldPrim->pGraphic, (BigScalar) 0, Cn3d_LayerNow, Cn3D_ColorIndex(&Cn3D_ColorData, iColor),
                 pclabel,  iX, iY , iZ , (Uint4) (Scale*VIEWSCALE/2) , 0, flags);
#endif /* else _OPENGL */
   return;

}


void LIBCALL RenderBond( PALD paldFrom,  PALD paldTo, ResidueColorCell * iColor,
			FloatLo fCylRadius)
{

  /* paldFrom owns the bond */

   Int4  iXFrom, iYFrom, iZFrom;
   FloatLo fXFrom,  fYFrom,  fZFrom;
   FloatLo fXTempF,  fYTempF,  fZTempF;
   FloatLo fXTempT,  fYTempT,  fZTempT;
   Int4  iXTo, iYTo, iZTo;
   FloatLo fXTo,  fYTo,  fZTo;
   Int4 iCylRadius;
   PMSD pmsdParent = NULL;
   PMAD pmadFrom = NULL;
   PMAD pmadTo = NULL;

#ifdef DEBUG_N
/*printf("%f\n", (float) fCylRadius);*/
#endif

  pmadFrom = (PMAD) paldFrom->pfbParent;
   pmadTo = (PMAD) paldTo->pfbParent;
   if ((!pmadFrom) || (!pmadTo)) return;
   if ((AtomicNumber(pmadFrom) == Atom_element_h) &&
      (Cn3d_DoHydrogens == FALSE)) return;
   if ((AtomicNumber(pmadTo) == Atom_element_h) &&
      (Cn3d_DoHydrogens == FALSE)) return;
   pmsdParent = ToMSDParent((PFB)pmadFrom);
   fXFrom = (FloatLo) (paldFrom->pflvData[0]);
   fYFrom = (FloatLo) (paldFrom->pflvData[1]);
   fZFrom = (FloatLo) (paldFrom->pflvData[2]);
   fXTo = (FloatLo) (paldTo->pflvData[0]);
   fYTo = (FloatLo) (paldTo->pflvData[1]);
   fZTo = (FloatLo) (paldTo->pflvData[2]);

   if (pmsdParent->pflTranslate)
     /* do translation */
     {
	fXFrom = fXFrom - pmsdParent->pflTranslate[0];
	fYFrom = fYFrom - pmsdParent->pflTranslate[1];
	fZFrom = fZFrom - pmsdParent->pflTranslate[2];
	fXTo = fXTo - pmsdParent->pflTranslate[0];
	fYTo = fYTo - pmsdParent->pflTranslate[1];
	fZTo = fZTo - pmsdParent->pflTranslate[2];
     }
  if (pmsdParent->ppflRotate)
     /* do rotation */
     {
	fXTempF = (FloatLo) (fXFrom * pmsdParent->ppflRotate[0][0] +
	         fYFrom * pmsdParent->ppflRotate[1][0] +
		 fZFrom * pmsdParent->ppflRotate[2][0]);
	fYTempF= (FloatLo) ( fXFrom * pmsdParent->ppflRotate[0][1] +
	         fYFrom * pmsdParent->ppflRotate[1][1] +
		 fZFrom * pmsdParent->ppflRotate[2][1]);
	fZTempF = (FloatLo) ( fXFrom * pmsdParent->ppflRotate[0][2] +
	         fYFrom * pmsdParent->ppflRotate[1][2] +
		 fZFrom * pmsdParent->ppflRotate[2][2]);
        fXTempT = (FloatLo) ( fXTo * pmsdParent->ppflRotate[0][0] +
	         fYTo * pmsdParent->ppflRotate[1][0] +
		 fZTo * pmsdParent->ppflRotate[2][0]);
	fYTempT = (FloatLo) ( fXTo * pmsdParent->ppflRotate[0][1] +
	         fYTo * pmsdParent->ppflRotate[1][1] +
		 fZTo * pmsdParent->ppflRotate[2][1]);
	fZTempT = (FloatLo) ( fXTo * pmsdParent->ppflRotate[0][2] +
	         fYTo * pmsdParent->ppflRotate[1][2] +
		 fZTo * pmsdParent->ppflRotate[2][2]);
        fXFrom = fXTempF;
	fYFrom = fYTempF;
	fZFrom = fZTempF;
	fXTo   = fXTempT;
	fYTo   = fYTempT;
	fZTo   = fZTempT;
     }
      /* scale */
    iXFrom = (Int4) (fXFrom*VIEWSCALE);
    iYFrom = (Int4) (fYFrom*VIEWSCALE);
    iZFrom = (Int4) (fZFrom*VIEWSCALE);
    iXTo = (Int4) (fXTo*VIEWSCALE);
    iYTo = (Int4) (fYTo*VIEWSCALE);
    iZTo = (Int4) (fZTo*VIEWSCALE);
    iCylRadius = (Int4) (fCylRadius*VIEWSCALE);
    Cn3d_AnyPrim = TRUE;
#ifndef _OPENGL
    if (paldFrom->pGraphic == NULL)
       paldFrom->pGraphic = (Pointer) AddSegment3D(pic, NULL,
			(BigScalar) paldFrom, Cn3d_LayerNow);
#endif /* ndef _OPENGL */
    if (fCylRadius  <  CYL_THRESHOLD)
#if _OPENGL
         OGL_AddLine3D(Cn3D_ColorData.OGL_Data, iColor, (FloatHi)fXFrom, (FloatHi)fYFrom, (FloatHi)fZFrom,
        (FloatHi)fXTo, (FloatHi)fYTo, (FloatHi)fZTo);

#else /* _OPENGL */
       AddLine3D (pic, (Segment3D) paldFrom->pGraphic , (BigScalar) paldFrom, Cn3d_LayerNow, Cn3D_ColorIndex(&Cn3D_ColorData, iColor), iXFrom,iYFrom,iZFrom,
                                           iXTo,iYTo,iZTo);
#endif /* else _OPENGL */
    else
#if _OPENGL
        OGL_AddCylinder3D(Cn3D_ColorData.OGL_Data, iColor, (FloatHi)fXFrom, (FloatHi)fYFrom, (FloatHi)fZFrom,
        (FloatHi)fXTo, (FloatHi)fYTo, (FloatHi)fZTo, (FloatHi)fCylRadius);

#else /* _OPENGL */
       AddCylinder3D (pic, (Segment3D) paldFrom->pGraphic, (BigScalar) paldFrom, Cn3d_LayerNow , Cn3D_ColorIndex(&Cn3D_ColorData, iColor), iXFrom,iYFrom,iZFrom,
                                           iXTo,iYTo,iZTo, iCylRadius);
#endif /* else _OPENGL */
   return;
}




void  LIBCALL  RenderAnAtom(PALD paldAtom, ResidueColorCell * iColor,
			FloatLo fRadius)
{
    Int4  iXAtom, iYAtom, iZAtom;
   FloatLo fXAtom,  fYAtom,  fZAtom;
   FloatLo fXTempF,  fYTempF,  fZTempF;
   Int4 iRadius;
   PMSD pmsdParent = NULL;
   PMAD pmadFrom = NULL;
   PMAD pmadAtom = NULL;
   if (!paldAtom) return;
   pmadAtom = (PMAD) paldAtom->pfbParent;
   if (!pmadAtom) return;
   if ((AtomicNumber(pmadAtom) == Atom_element_h) &&
      (Cn3d_DoHydrogens == FALSE)) return;
   pmsdParent = ToMSDParent((PFB)pmadAtom);
   fXAtom = (FloatLo) (paldAtom->pflvData[0]);
   fYAtom = (FloatLo) (paldAtom->pflvData[1]);
   fZAtom = (FloatLo) (paldAtom->pflvData[2]);

   if (pmsdParent->pflTranslate)
     /* do translation */
     {
	fXAtom = fXAtom - pmsdParent->pflTranslate[0];
	fYAtom = fYAtom - pmsdParent->pflTranslate[1];
	fZAtom = fZAtom - pmsdParent->pflTranslate[2];

     }
  if (pmsdParent->ppflRotate)
     /* do rotation */
     {
	fXTempF = (FloatLo) (fXAtom * pmsdParent->ppflRotate[0][0] +
	         fYAtom * pmsdParent->ppflRotate[1][0] +
		 fZAtom * pmsdParent->ppflRotate[2][0]);
	fYTempF= (FloatLo) ( fXAtom * pmsdParent->ppflRotate[0][1] +
	         fYAtom * pmsdParent->ppflRotate[1][1] +
		 fZAtom * pmsdParent->ppflRotate[2][1]);
	fZTempF = (FloatLo) ( fXAtom * pmsdParent->ppflRotate[0][2] +
	         fYAtom * pmsdParent->ppflRotate[1][2] +
		 fZAtom * pmsdParent->ppflRotate[2][2]);
        fXAtom = fXTempF;
	fYAtom = fYTempF;
	fZAtom = fZTempF;
     }

#if _OPENGL
  
    OGL_AddSphere3D(Cn3D_ColorData.OGL_Data, (TOGL_ColorCell *)iColor, (FloatHi)fXAtom, (FloatHi)fYAtom, (FloatHi)fZAtom,
        (FloatHi)fRadius);
#else /* _OPENGL */

    Cn3d_AnyPrim = TRUE;  /* this probably can be deleted */
                    /* yanli moved it out of OpenGL ifdef block */  
      /* scale */
    iXAtom = (Int4) (fXAtom*VIEWSCALE);
    iYAtom = (Int4) (fYAtom*VIEWSCALE);
    iZAtom = (Int4) (fZAtom*VIEWSCALE);
    iRadius = (Int4) (fRadius*VIEWSCALE);
    if (paldAtom->pGraphic == NULL)
       paldAtom->pGraphic = (Pointer) AddSegment3D(pic, NULL, (BigScalar) paldAtom, Cn3d_LayerNow);
    AddSphere3D (pic, (Segment3D) paldAtom->pGraphic, (BigScalar) paldAtom, Cn3d_LayerNow,  Cn3D_ColorIndex(&Cn3D_ColorData, iColor),
                iXAtom,iYAtom,iZAtom, iRadius);
#endif /* else _OPENGL */

   return;
}



/*******************************************************/

static void RenderAllAtom(PFB pfbThis, Int4 iModel,  Int4 iIndex, Pointer ptr)
{
  PRK prKeep;
  PVNMA pvnmaThis = NULL;
  PMAD pmadThis = NULL;
  PALD paldThis = NULL;
  PVNMB pvnmbThis = NULL;
  PMBD pmbdThis = NULL;
  PMMD pmmdTo = NULL;
  PMMD pmmdThis;
  FloatLo flAtemp = (FloatLo)0;
  Int4 iTemp = 0;
  Int2 iBin = 0;
  FloatLo fRadius = (FloatLo)0;
  Boolean bDraw = TRUE;
  ValNodePtr pvnB;
  PALD paldMid = NULL;
  PALD paldDrawTo = NULL;
  PMAD pmadDrawTo = NULL;
  ResidueColorCell iColor;
/*  PMGD pmgdParent; */

  prKeep = (PRK) ptr;
  if (!prKeep) return;
  if (!IsAtomNode(pfbThis)) return;
  if (!(pfbThis->bMe & prKeep->NodeType)) return;
  if (prKeep->NodeWhat != CONVERT_ALL)
    if (!(pfbThis->bWhat & prKeep->NodeWhat)) return;
  if ((int) iIndex == RESIDUES)
      {
        if (IsAtomBackBone(pfbThis)) return;
      }
  pmadThis = (PMAD) pfbThis;
  if ((int) iIndex == REALBB)
    {  /* throw out the carbonyls */
	if (pmadThis->bWhat & AM_OCARBNYL) return;
	if (pmadThis->bWhat & AM_C1RIBOSE) return;
    }
  if ((int) iIndex == REALXTRABB)
    {
        if (pmadThis->bWhat & AM_C1RIBOSE) return;
    }
   if ((int) iIndex == IONSON)
    {
      if (!(pmadThis->bUpdate & AM_ION)) return;
    }
  if ((int) iIndex == HETSON)
     {
	 if ((pmadThis->bUpdate & AM_ION)) return;
     }
  paldThis = GetAtomLocs(pmadThis, iModel);
  if (!paldThis) return;  /* no location corresponding to this model */


  /* make the atomic-level color decisions here */
  Cn3D_CopyColorCell(&iColor, &(prKeep->Color));
  if (Cn3d_ColorNow == C_CPK)
    {
      Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ElementKinColors[(Int1)pmadThis->pvnmaLink->choice]]));
      if (Cn3d_ColorPass) {Cn3D_RequestColor(&Cn3D_ColorData, &iColor); return; }
    }

  if (Cn3d_ColorNow == C_BYTEMP)
    {
	/* calculate temp factor */
	if (paldThis->iFloatNo == 4) /* istotropic */
            {
	     iTemp = (Int4) paldThis->pflvData[4] * 100;
	    }
	if (paldThis->iFloatNo == 9) /* anisotropic */
	         {  /* calculate the isotropic temp factor */
		     flAtemp = (FloatLo) (((paldThis->pflvData[4] +
					    paldThis->pflvData[5] +
					    paldThis->pflvData[6]) / 3));
		     iTemp = (Int4) flAtemp * 100;
		 }
#ifdef DEBUG_Q
  printf("temp=%d\n", (int) iTemp);
#endif
	if (iTemp < TempsKine[0]) return;  /* too low */
	if (iTemp > TempsKine[15]) return; /* too high */
	/* find the bin */
	for (iBin = 1; iBin < KIN_COLOR_THERM; iBin++)
	  {
	      if ((iTemp < TempsKine[iBin+1]) && (iTemp >= TempsKine[iBin]))
	        break;
	  }
	Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ThermKine[iBin]]));
        if (Cn3d_ColorPass)  {Cn3D_RequestColor(&Cn3D_ColorData, &iColor); return; }
    }

  /* otherwise the color should have been set by the calling routine */

   bDraw = TRUE;
  /* set up the atom drawing size */
   switch (prKeep->Atom)
     {

       case ATOM_SPACE:
           fRadius = (FloatLo) ((ElementSize((Int1)pmadThis->pvnmaLink->choice)));
         break;
       case ATOM_2XBOND:
          fRadius = prKeep->BondWidth * 2;
         break;
       case ATOM_ISBOND:
          fRadius = prKeep->BondWidth * (float)0.9;
         break;
       case ATOM_NONE:
        default:
	  bDraw = FALSE;
     }
#ifdef DEBUG_N
/*printf("atom radius %f ", (float) fRadius);*/
#endif
   if ((Cn3d_DoHydrogens == FALSE) && (pmadThis->pvnmaLink->choice == 1)) bDraw = FALSE;
   if (Cn3d_ColorPass)  { Cn3D_RequestColor(&Cn3D_ColorData, &iColor); return; }
   if (bDraw) RenderAnAtom(paldThis, &iColor, fRadius);

   /* the rest of this routine handles half-bond drawing where the bonds */
   /* are owned by the parent atom */
   bDraw = TRUE;
   if (prKeep->Bond)
     {
         pvnB = pmadThis->pvnBonds; /* local bond list */
 	 while (pvnB)
	  {
	      pmbdThis = (PMBD) pvnB->data.ptrvalue;
	      if (pmbdThis->pmadTo == pmadThis) pmadDrawTo = pmbdThis->pmadFrom;
	      else pmadDrawTo = pmbdThis->pmadTo;
	      if (!pmadDrawTo) goto nextbond;
	      /*  bond handler */
	      bDraw = TRUE;
	      if ((int) iIndex == VIRTUALBB)
	        {   /* don't draw virtual bonds */
		   if (!IsBondVirtual(pmbdThis)) bDraw = FALSE;
		}
              if ((int) iIndex == RESIDUES)
               {
                 if (IsBondVirtual(pmbdThis)) bDraw = FALSE;
	       }
	      if ((int) iIndex == REALXTRABB)
	        {  /* don't draw to other non backbone atoms */
		    if (!(IsAtomBackBone(pmadDrawTo))) bDraw = FALSE;
 		    if (IsBondVirtual(pmbdThis)) bDraw = FALSE;
		}
	      if ((int) iIndex == REALBB)
	        {  /* don't draw to non-backbone or carbonyls */
		    if (!(IsAtomBackBone(pmadDrawTo))) bDraw = FALSE;
		    if (IsAtomOCarbonyl(pmadDrawTo)) bDraw = FALSE;
		    if (IsBondVirtual(pmbdThis)) bDraw = FALSE;
		}
	      if ((Cn3d_DoHydrogens == FALSE) && (pmadDrawTo->pvnmaLink->choice == 1)) bDraw = FALSE;
	      if (bDraw)
	        {
		    paldDrawTo = GetAtomLocs(pmadDrawTo,  iModel);
		    if (!paldDrawTo) goto nextbond;
	 	    if (!paldDrawTo->pflvData) goto nextbond;
		    pmmdTo = GetParentMol((PFB)paldDrawTo);
		    pmmdThis = GetParentMol((PFB)paldThis);
		    if ((pmmdTo != pmmdThis) && ((int) iIndex != CONNECTON))
		       {
		         goto nextbond; /* don't connect inter-mol bonds */
		       }
 		    if (((int) iIndex == RESIDUES) && (IsAtomBackBone(pmadDrawTo)))
		      {
 			/* draw a whole bond to bridge */
		        RenderBond(paldThis, paldDrawTo, &iColor, prKeep->BondWidth);
		      }
		    else
		      {
		    	paldMid = NewALD(); /* phoney location */
		    	if (!paldMid) goto nextbond;
		    	paldMid->pfbParent = (PFB) pmadThis; /* link to this  atom */
		    	paldMid->iFloatNo = (Int1) 3;
		    	paldMid->cAltConf = ' ';
		    	paldMid->pflvData = FLVector(0,3);
		    	if (!paldMid->pflvData) {FreeALD(paldMid); return;}
		    	/* midpoint between the two atoms */
		    	paldMid->pflvData[0] = (paldThis->pflvData[0] + paldDrawTo->pflvData[0]) / (float)2.0;
		    	paldMid->pflvData[1] = (paldThis->pflvData[1] + paldDrawTo->pflvData[1]) / (float)2.0;
		    	paldMid->pflvData[2] = (paldThis->pflvData[2] + paldDrawTo->pflvData[2]) / (float)2.0;
		    	if (paldMid && paldThis)
			  {
					
					RenderBond(paldThis, paldMid, &iColor,  prKeep->BondWidth);
			  }
                    	FreeALD(paldMid);
		      }
		}
	    nextbond:
            pvnB = pvnB->next;
	  }
     }
   return;
}


Int2 LIBCALL  GetGraphNCBIstdaa(PMGD pmgdThis)
{

    Int2 i;

    if (!(IsGraphAminoAcid(pmgdThis))) return 21;
    for (i=0; i<MAX_NCBIstdaa; i++)
      if (pmgdThis->pcIUPAC[0] == NCBIstdaaUC[i])
        {
          return i;
        }
    return 21;  /* Xxx */
}

Int2 LIBCALL GetGraphNCBI4na(PMGD pmgdThis)
{
    Int2 i;
    if (!(IsGraphNABase(pmgdThis))) return 16;
    for (i=0; i<MAX_NCBIstdaa; i++)
      if (pmgdThis->pcIUPAC[0] == NCBI4naUC[i]) return i;
    return 16;  /* N (any) */
}




static void LIBCALLBACK RenderGraph(PFB pfbThis, Int4 iModel,  Int4 iIndex, Pointer ptr)
{	/* traverser callback */
    /* Sets color scheme for rendering a graph */

  PRK prKeep;
  PRK prkNew = NULL;
  Int2 iResId;
  Int2 iEnd;
  Int2 iSkip;
  PDNMG pdnmgThis = NULL;
  PMSD pmsdThis = NULL;
  PMMD pmmdThis = NULL;
  PMGD pmgdThis = NULL;
  PVNMA pvnmaThis = NULL;
  PMAD pmadThis = NULL;
  PVNMB pvnmbThis = NULL;
  PMBD pmbdThis = NULL;
  CharPtr pcLabel = NULL;
  CharPtr pcLNum = NULL;
  CharPtr pcDash = NULL;
  CharPtr pcTemp = NULL;
  CharPtr pcTemp2 = NULL;
  CharPtr pcL = NULL;
  PDNMG pdnmgFrom = NULL;
  PDNMG pdnmgTo = NULL;
  PMAD pmadFrom = NULL;
  PMAD pmadTo = NULL;
  PALD paldThis = NULL;
  PALD paldFrom = NULL;
  PALD paldTo = NULL;
  Byte bReservedThis = 0;
  /* Uint1Ptr rgb; */

  SpecialFeatureInfoPtr sfipThis = NULL;

  prKeep = (PRK) ptr;
  if (!prKeep) return;
  prkNew = CopyRenderKeep(prKeep);
  if (!prkNew) return;
  if (!IsGraphNode(pfbThis)) goto cyalater;
  if (!(pfbThis->bMe & prKeep->NodeType)) goto cyalater;
  if (!(pfbThis->bWhat & prKeep->NodeWhat)) goto cyalater;
  pdnmgThis =  DNFromPFB(pfbThis);
  pmgdThis = (PMGD) pfbThis;

  if(Cn3D_DisplayHighlight) {
     if(pmgdThis && pmgdThis->bHighlighted !=1) goto cyalater;
  }

  if(Cn3d_ParsNow == 0 && pmgdThis->FeatureOn > 0) goto cyalater;

  if(Cn3d_ParsNow > 0) {
     if( pmgdThis->iUserDefinedFeature == 0 || pmgdThis->FeatureOn == 0 || pmgdThis->iUserDefinedFeature != Cn3d_ParsNow) goto cyalater; 
  }

  if(pmgdThis->bVisible != 1) goto cyalater;
  if(pmgdThis->bTurnedOff == 1) goto cyalater;
     /* control display show/off on MG level-- Yanli */

  pmsdThis = ToMSDParent(pfbThis);
  ASSERT(pmsdThis != NULL);


  /* check to see if we want to show aligned/unaligned regions */
  if(IsGraphAminoAcid(pmgdThis))
  {
    if(pmsdThis->bMaster) 
    {
      if (pmgdThis->bReserved && (pmgdThis->bReserved == pmsdThis->bAligned) && !Cn3D_fAlignOn && pmgdThis->iUserDefinedFeature != Cn3d_ParsNow) goto cyalater;
      if ((!(pmgdThis->bReserved) || (pmgdThis->bReserved != pmsdThis->bAligned)) && !Cn3D_fUnalignOn&& pmgdThis->iUserDefinedFeature != Cn3d_ParsNow) goto cyalater;
    }
    else
    {
      if (pmgdThis->bReserved)
      {
        if ( (*(pmgdThis->pbMasterReserved) == *(pmsdThis->pbAligned)) && !Cn3D_fAlignOn && pmgdThis->iUserDefinedFeature != Cn3d_ParsNow) goto cyalater;
        if ( (*(pmgdThis->pbMasterReserved) != *(pmsdThis->pbAligned)) && !Cn3D_fUnalignOn && pmgdThis->iUserDefinedFeature != Cn3d_ParsNow) goto cyalater;
      }
      if (!(pmgdThis->bReserved) && !Cn3D_fUnalignOn && pmgdThis->iUserDefinedFeature != Cn3d_ParsNow) goto cyalater;
    }
  }
  
  /* else if (pmgdThis->bReserved)
  {
    if ((*(pmgdThis->pbMasterReserved) == *(pmsdThis->pbAligned)) && !Cn3D_fAlignOn && pmgdThis->iUserDefinedFeature != Cn3d_ParsNow) goto cyalater;
  }
  else if ((!(pmgdThis->bReserved) || (*(pmgdThis->pbMasterReserved) != *(pmsdThis->pbAligned))) && !Cn3D_fUnalignOn && pmgdThis->iUserDefinedFeature != Cn3d_ParsNow) goto cyalater;
*/
  /* set the color */
#ifdef _OPENGL
  if(pmgdThis->bHighlighted) {
      if (IsGraphAminoAcid(pmgdThis)) {
          Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3D_ColorData.Highlight.ColorCell));
          if (Cn3d_ColorPass) { Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
          if( Cn3D_ReColor ){
              if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
                  if(iIndex != NTBBLABELSON && iIndex != PBBLABELSON){
                      /*rgb = (Uint1Ptr) GetRGB((Int2) prkNew->Color); lyg*/
                      ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);
                  }
              }
          }
      }
  }
  else 
#endif /* _OPENGL */
  if(Cn3d_ParsNow != 0 && pmgdThis->iUserDefinedFeature == Cn3d_ParsNow && pmgdThis->FeatureOn > 0){
     if (IsGraphAminoAcid(pmgdThis)) {
        Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[Cn3d_ColorNow]));
        if (Cn3d_ColorPass) { Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
        if( Cn3D_ReColor ){
           if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
              if(iIndex != NTBBLABELSON && iIndex != PBBLABELSON){
                 /*rgb = (Uint1Ptr) GetRGB((Int2) prkNew->Color); lyg*/
                 ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);
              }
           }
        }
     }
  }
  else if (Cn3d_ColorNow == C_BYRES)
  {
    if IsGraphAminoAcid(pmgdThis)
    {
      iResId = GetGraphNCBIstdaa(pmgdThis);
      Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[KinAAColor[iResId]]));
      if (Cn3d_ColorPass) { Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
      if( Cn3D_ReColor ){
         if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
            /*rgb = (Uint1Ptr) GetRGB((Int2) prkNew->Color); lyg */
            ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);
         }                    /* yanli */
      }
      if(Salsa_BioseqUpdate) Salsa_BioseqUpdate = FALSE;
    }
    else
      if IsGraphNABase(pmgdThis)
      {
        iResId = GetGraphNCBI4na(pmgdThis);
        Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[KinNAColor[iResId]]));
        if (Cn3d_ColorPass) {  Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
      }
  }
  /* if (Cn3d_ColorNow == C_BYLEN)
  {
  pmmdThis = GetParentMol((PFB)pmgdThis);
  iColorBin = (pmmdThis->iResCount-1) / KIN_COLOR_THERM;
  prkNew->Color = (Int2)ThermKine[
  (IndexFromNode((PFB) pmgdThis)-1 / iColorBin)];
  if (Cn3d_ColorPass) {  parsColor->IndexRGB[prkNew->Color]= 1 ; goto cyalater; }
  }
  */
  else if (Cn3d_ColorNow == C_BYSSTRU)
  {
    if (!IsGraphAminoAcid(pmgdThis)) goto cyalater;
    Cn3D_CopyColorCell(&(prkNew->Color), ColorFromSS(pmgdThis));
    if (Cn3d_ColorPass) {  Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
    if( Cn3D_ReColor ){
       if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
          /* rgb = (Uint1Ptr) GetRGB((Int2) prkNew->Color); lyg */
          ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);
       }                    /* yanli */
    }
    if(Salsa_BioseqUpdate) Salsa_BioseqUpdate = FALSE;
  }
  /*set the color for alignment.  skip nucleic acids. */
  else if (Cn3d_ColorNow == C_BYALIGN)
  {
    if (!IsGraphAminoAcid(pmgdThis)) goto cyalater;
    Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[bColorAlignments[(Cn3d_lSlaveNum % NUM_SLAVES)]]));
    if (Cn3d_ColorPass) {  Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
    if( Cn3D_ReColor ){
       if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
          /* rgb = (Uint1Ptr) GetRGB((Int2) prkNew->Color); lyg */
          ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);
       }                    /* yanli */
    }
    if(Salsa_BioseqUpdate) Salsa_BioseqUpdate = FALSE;
  }

  else if (Cn3d_ColorNow == C_BYCONS)
  {
    if (!IsGraphAminoAcid(pmgdThis)) goto cyalater;
/*    if (pmgdThis->pbMasterReserved) bReservedThis = *(pmgdThis->pbMasterReserved);
    else bReservedThis = pmgdThis->bReserved;
    prkNew->Color = (Int2) bColorConservation[(bReservedThis % NUM_SLAVES)];
*/
    if(pmsdThis->bMaster) 
    {
      if (pmgdThis->bReserved && (pmgdThis->bReserved == pmsdThis->bAligned)) 
      {
        Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[C_red])); /* change conserved color to red */
      }
      else
      {
        Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[bColorAlignments[(Cn3d_lSlaveNum % NUM_SLAVES)]]));
      }
    }
    else
    {
      if (pmgdThis->bReserved) 
      {
        if (*(pmgdThis->pbMasterReserved) == *(pmsdThis->pbAligned)) 
        {
          Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[C_red])); /* change conserved color to red */
        }
        else
        {
          Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[bColorAlignments[(Cn3d_lSlaveNum % NUM_SLAVES)]]));
        }
      }
      else
      {
        Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[bColorAlignments[(Cn3d_lSlaveNum % NUM_SLAVES)]]));
      }
      
    }
    
    if (Cn3d_ColorPass) {  Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
    if( Cn3D_ReColor ){
       if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
          /* rgb = (Uint1Ptr) GetRGB((Int2) prkNew->Color); lyg */
          ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);    
       }                    /* yanli */
    }
    if(Salsa_BioseqUpdate) Salsa_BioseqUpdate = FALSE;
  }
  
  else if (Cn3d_ColorNow == C_BYSEQCONS)
  {
    if(pmgdThis->bReserved) {
       Cn3D_CopyColorCell(&(prkNew->Color), GetColorIndexForMG(pmgdThis));
    }
    else
    {
       Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[bColorAlignments[(Cn3d_lSlaveNum % NUM_SLAVES)]]));
    }
    if (Cn3d_ColorPass) {  Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
    if(Salsa_BioseqUpdate) Salsa_BioseqUpdate = FALSE;
  }

  else if (Cn3d_ColorNow == C_BYDOMAIN)
  {
    if (!IsGraphAminoAcid(pmgdThis)) goto cyalater;
    Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[ColorNumKinBB[(pmgdThis->iDomain % KIN_COLOR_NUM)]])) ;
    if (Cn3d_ColorPass) {  Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
    if( Cn3D_ReColor ){
       if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
          if(prkNew->Color.rgb[0] == 255 && prkNew->Color.rgb[1] == 255 
              && prkNew->Color.rgb[2] == 255) Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[C_black]));
          /* white is not visible in salsa */
          ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);
       }                    /* yanli */
    }
    if(Salsa_BioseqUpdate) Salsa_BioseqUpdate = FALSE;
  }
  else if (Cn3d_ColorNow == C_BYHYDRO)
  {
    if IsGraphAminoAcid(pmgdThis)
    {
      iResId = GetGraphNCBIstdaa(pmgdThis);
      Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[PhobeAAColor[iResId]]));
      if (Cn3d_ColorPass) { Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
    if( Cn3D_ReColor ){
       if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
          /* rgb = (Uint1Ptr) GetRGB((Int2) prkNew->Color); */
          ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);
       }                    /* yanli */
    }
    if(Salsa_BioseqUpdate) Salsa_BioseqUpdate = FALSE;
    }
    else
      if IsGraphNABase(pmgdThis)
      {
        iResId = GetGraphNCBI4na(pmgdThis);
        Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[KinNAColor[iResId]]));
        if (Cn3d_ColorPass) { Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
      }
  }
  else if (Cn3d_ColorNow == C_BYCHAIN)
  {
    if( Cn3D_ReColor ){
       if(!Cn3d_ColorPass && Cn3D_ObjMgrOpen && !Salsa_BioseqUpdate){
          /* rgb = (Uint1Ptr) GetRGB((Int2) prkNew->Color);*/
          ColorSalsa_BYMG(pmgdThis, prkNew->Color.rgb);
       }                    /* yanli */
    }
    if(Salsa_BioseqUpdate) Salsa_BioseqUpdate = FALSE;
  }

   /* do the labeling */

  if ((iIndex == NTBBLABELSON) || (iIndex == PBBLABELSON))
    {

       if (Cn3d_ColorPass)
          {
	     /*if (prkNew->Color < C_top)*/
               { Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto cyalater; }
          }

       switch ((int) Cn3d_LabelNow)
         {
           case 2: iSkip = 1;
		break;
           case 3: iSkip = 5;
		break;
           case 4: iSkip = 10;
		break;
           case 5: iSkip = 20;
		break;
           case 6: iSkip = 50;
		break;
           default: iSkip = 1;
	 }
       if ( ((int) pdnmgThis->choice % (Int4) iSkip) != 0) goto cyalater;

       if (IsGraphAminoAcid(pmgdThis) && (prkNew->LStyle & L_NAME))
         {

	    iResId = GetGraphNCBIstdaa(pmgdThis);

	    if  (prkNew->LStyle & L_PDB)
               {
	          pcLabel = StringSave(pmgdThis->pcGraphName);
	 	  if (!pcLabel)
		    pcLabel = StringSave("UNK");
		  if (StringLen(pcLabel) > 3) pcLabel[3] = '\0'; /* truncate SER COOH dict. */
               }
            else
            if  (prkNew->LStyle & L_1LETR)
	         pcLabel = StringSave(pmgdThis->pcIUPAC);
            else
	        pcLabel = StringSave(AminoAcidNameFromIdx(iResId, USE_MIXCASE));
	 }


       if (IsGraphNABase(pmgdThis) && (prkNew->LStyle & L_NAME))
         {
            if ((prkNew->LStyle & L_1LETR) && (!(prkNew->LStyle & L_PDB)))
	        pcLabel = StringSave(pmgdThis->pcIUPAC);
	    else
               {
		  pcLabel = StringSave(pmgdThis->pcGraphName);
		  if (!pcLabel)
		    pcLabel = StringSave("UNK");
		  pcTemp = pcLabel;
		  if ((pmgdThis->bWhat & (Byte) DICT_LOCAL) &&
		    ((pcLabel[0]=='D') || (pcLabel[0]=='R')) &&
		    (pcLabel[1]=='N') &&
		    (pcLabel[2]=='A'))
		   {
		      pcTemp = (CharPtr)&pcLabel[3];
		   }
		  if (StringLen(pcTemp) > 3) pcTemp[3] = '\0'; /* truncate SER COOH dict. */
                  pcTemp2 = pcLabel;
                  pcLabel = StringSave( pcTemp );
                  MemFree(pcTemp2);
	       }
	 }

       if (prkNew->LStyle & L_NUM)
         {
	   if (prkNew->LStyle & L_PDB)
	     {
	       pcLNum = StringSave(pmgdThis->pcGraphNum);
	       if (!pcLNum) pcLNum = StringSave("?");
	        /* remove leading spaces */
               pcTemp = pcLNum;
               while ((pcTemp[0] == ' ') && (pcTemp[1] != '\0'))
                 {
                     pcTemp = (CharPtr) &pcTemp[1];
                 }
               pcTemp2 = pcLNum;
               pcLNum = StringSave(pcTemp);
               MemFree(pcTemp2);
	     }

	   else
	     {
		 pcLNum = (CharPtr) MemNew((size_t)(INTSTRLEN));
                 if (!pcLNum) goto cyalater;
		 sprintf(pcLNum, "%ld", (long) pdnmgThis->choice);
	     }
	 }

       if ((pcLabel == NULL) && (pcLNum == NULL))
         pcL = StringSave(" ");
       else
       if (pcLabel == NULL)
	  {
	         pcL = StringSave(pcLNum);
	         MemFree(pcLNum);
	  }
       else
         if (pcLNum == NULL)
	   {
		   pcL = StringSave(pcLabel);
		   MemFree(pcLabel);
	   }
         else
           {
	     pcDash = StringSave("-");
             if (pcLNum == NULL) pcLNum = StringSave(" ");
             pcL = (CharPtr) MemNew((size_t) (StringLen(pcLabel) + StringLen(pcDash) + StringLen(pcLNum) + 4));
             if (!pcL) goto cyalater;
             StringCpy(pcL,  pcLabel);
             StringCat(pcL,  pcDash);
             StringCat(pcL,  pcLNum);
             MemFree(pcLabel);
             MemFree(pcDash);
             MemFree(pcLNum);
	   }
          /* pcL is label ready to go */

       pmadThis = GetMainAtom(pdnmgThis);
       paldThis = GetAtomLocs(pmadThis,  iModel);
       if (!paldThis) goto nolabel;
       iEnd = RL_BETACARB;

       pdnmgTo = GetPrevGraph(pdnmgThis);
       pmadTo = GetMainAtom(pdnmgTo);
       paldTo = GetAtomLocs(pmadTo,  iModel);
       if (!paldTo) iEnd = RL_CENTERPLUSY;

       pdnmgFrom = GetNextGraph(pdnmgThis);
       pmadFrom = GetMainAtom(pdnmgFrom);
       paldFrom = GetAtomLocs(pmadFrom,  iModel);
       if (!paldFrom) iEnd = RL_CENTERPLUSY;

       RenderLabel(pmsdThis->pdnmsLink, pcL,
			  paldThis, paldFrom,  paldTo,
			  &(prkNew->Color), prkNew->LJust, prkNew->LScale, iEnd);
       MemFree(pcL);

    }

 nolabel:
#ifdef _OPENGL
		 OGL_LoadName((Nlm_VoidPtr)pmgdThis);
#endif

   /* do the rendering */
   /* deal with setting the NodeType/NodeWhat for virtual/backbone/realbb/residues */

   if ((prkNew->Bond == HALF_BOND) || (prkNew->Atom != ATOM_NONE))
       {  /* draw atoms or half-bonds */
         prkNew->NodeType = (Byte) AM_MAD;
	 switch ((int) iIndex)
	   {
	    case REALXTRABB:
	    case REALBB:
	        prkNew->NodeWhat = (Byte) (AM_BACKBONE);
	      break;
	    case VIRTUALBB:
		prkNew->NodeWhat = (Byte) (AM_CALPHA | AM_PALPHA);  /* virtual */
	      break;
	    default:
	        prkNew->NodeWhat = (Byte) CONVERT_ALL;  /* do all the atoms */
	    }
         pvnmaThis = pmgdThis->pvnmaAHead;
	 while (pvnmaThis)
          {

		 /*looks like this is where the color decision has to be made! */

		 pmadThis = (PMAD) pvnmaThis->data.ptrvalue;
         	RenderAllAtom((PFB)pmadThis, iModel, iIndex, prkNew);
	    pvnmaThis = pvnmaThis->next;
	  }
       }

   if (!(
             ((int) iIndex == VIRTUALBB)
          || ((int) iIndex == REALBB)
          || ((int) iIndex == REALXTRABB)
          )
          && (prkNew->Bond))
      {  /*  Traverse the InterResidueBonds & draw them as half bonds */
	 /* connecting to this residue in its residue color */
	 prkNew->NodeType = (Byte) AM_MAD;
	 prkNew->NodeWhat = (Byte) CONVERT_ALL;
         pmmdThis = GetParentMol((PFB)pmgdThis);
   	 pvnmbThis = pmmdThis->pvnmbIRBHead;
  	 while (pvnmbThis)  /* walk the inter-res bond list by hand */
          {
	     pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
             if (GetParentGraph((PFB)pmbdThis->pmadFrom) == pmgdThis)
               {  /* this inter-res bond is in this graph */
	         RenderAllAtom((PFB)pmbdThis->pmadFrom, iModel, iIndex, prkNew);
	       }
	     if (GetParentGraph((PFB)pmbdThis->pmadTo) == pmgdThis)
	       {
	         RenderAllAtom((PFB)pmbdThis->pmadTo,  iModel,  iIndex, prkNew);
	       }
	     pvnmbThis = pvnmbThis->next;
	  }
      }

     cyalater:
     FreeRenderKeep(prkNew);
     return;

}



static void LIBCALLBACK RenderMolecule(PFB pfbThis, Int4 iModel,  Int4 iIndex, Pointer ptr)
{	/* traverser callback */
  PRK prKeep;
  PRK prkNew = NULL;
  PMSD pmsdThis = NULL;
  PDNMM pdnmmThis = NULL;
  PMMD pmmdThis = NULL;
  PVNMB pvnmbThis = NULL;
  PMBD pmbdThis = NULL;
  Int2 iTest;
  Byte bHold;
  PDNMG pdnmgFrom;
  PMAD pmadFrom;
  PALD paldFrom;
  PDNMG pdnmgTo;
  PMAD pmadTo;
  PALD paldTo;
  CharPtr pcN,  pcC;
  CharPtr pcChain = NULL;
  CharPtr pcTemp = NULL;


  prKeep = (PRK) ptr;
  if (!prKeep) return;
  prkNew = CopyRenderKeep(prKeep);
  if (!prkNew) return;
  if (!IsMoleculeNode(pfbThis)) goto byelater;
  if (!(pfbThis->bMe & prKeep->NodeType)) goto byelater;  /* same test for molecule */
  pmsdThis = ToMSDParent(pfbThis);
  ASSERT(pmsdThis != NULL);

  bHold = prKeep->NodeWhat;

  /* special case for anded DNA/RNA */
  if  ((prKeep->NodeWhat == (Byte) (AM_DNA | AM_RNA)) &&
      ((pfbThis->bWhat & AM_DNA) || (pfbThis->bWhat & AM_RNA)))
       bHold = pfbThis->bWhat;

  /* special case for anded WAT/SOL */
  else if  ((prKeep->NodeWhat == (Byte) (AM_WAT | AM_SOL)) &&
      ((pfbThis->bWhat & AM_WAT) || (pfbThis->bWhat & AM_SOL)))
       bHold = pfbThis->bWhat;

  /* special case for HET/ION some ions are also hets - e.g. heme(Fe) */
  else if ((prKeep->NodeWhat == (Byte) (AM_HET)) &&
        ((pfbThis->bWhat == (Byte) AM_ION) || pfbThis->bWhat == (Byte) AM_HET) ||
         (pfbThis->bWhat == (Byte) AM_POLY))
	bHold = pfbThis->bWhat;

  if ((pfbThis->bWhat != bHold)) goto byelater;  /* branch bound Prot & NT from others */

  if ((iIndex == HETSON) || (iIndex == IONSON))
     if (!((pfbThis->bWhat == AM_HET) ||
           (pfbThis->bWhat == AM_POLY) || (pfbThis->bWhat == AM_ION)))
        goto byelater;

  pdnmmThis =  DNFromPFB(pfbThis);
  pmmdThis = (PMMD) pfbThis;
 
  if(pmmdThis->bVisible != 1) goto byelater;     
        /* to control Cn3D display on MM level -- Yanli */

/*if(pmmdThis->bVisible == 0) goto byelater; */
        /* to control show/off on MM level -- Yanli */

  if (Cn3d_ColorNow == C_BYCHAIN)
    {  /* homogenous color scheme */
      if ((iIndex == RESIDUES))
        Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[ColorNumKinSC[(pdnmmThis->choice % KIN_COLOR_NUM)]]));
      else
        Cn3D_CopyColorCell(&(prkNew->Color), &(Cn3d_PaletteRGB[ColorNumKinBB[(pdnmmThis->choice % KIN_COLOR_NUM)]]));
      if (Cn3d_ColorPass) {  Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color)); goto byelater; }
    }

  if ((iIndex == PTERMLABON))
    {
	pcN = "N";
	pcC =  "C";
    }

  if  (iIndex == NTTERMLABON)
    {
	pcN  = "5'";
	pcC = "3'";
    }


  if (iIndex == NTTERMLABON || iIndex == PTERMLABON )
    {  /* take care of   TERMINI here */
             if (Cn3d_ColorPass)
              {
	         /*if (prkNew->Color < C_top) */
		   {
                     Cn3D_RequestColor(&Cn3D_ColorData,&(prkNew->Color));
		     goto byelater;
		   }
              }
 	    else
              {
                if (prkNew->LStyle & L_NAME)
                   pcChain = pmmdThis->pcMolName;
		if (!pcChain) pcChain = " ";
                pcTemp = (CharPtr) MemNew((size_t) (StringLen(pcChain) + StringLen(pcN) + 4));
                if (!pcTemp) goto byelater;
                StringCpy(pcTemp, pcChain);
                StringCat(pcTemp, " ");
                StringCat(pcTemp, pcN);
	        pdnmgTo = GetFirstGraph(pmmdThis);
	    	if (!pdnmgTo) goto byelater;
		pmadTo = GetMainAtom(pdnmgTo);
		paldTo = GetAtomLocs(pmadTo,  iModel);
    		while (paldTo == NULL)
		  {  /* walk in */
		      pdnmgTo = GetNextGraph(pdnmgTo);
		      if (!pdnmgTo) goto byelater; /* bail if no next graph */
		      pmadTo = GetMainAtom(pdnmgTo);
		      paldTo = GetAtomLocs(pmadTo,  iModel);
	          }
		pdnmgFrom = GetNextGraph(pdnmgTo);
		if (!pdnmgFrom) goto byelater;
	    	pmadFrom = GetMainAtom(pdnmgFrom);
	    	paldFrom = GetAtomLocs(pmadFrom,  iModel);
	    	RenderLabel(pmsdThis->pdnmsLink, pcTemp,
			  NULL, paldFrom,  paldTo,
			  &(prkNew->Color), prkNew->LJust, prkNew->LScale, RL_EXTRAPOL);
                MemFree(pcTemp);
	    	pdnmgTo = NULL;
	    	pdnmgFrom = NULL;
	    	pmadFrom = NULL;
	    	pmadTo = NULL;
	    	paldFrom = NULL;
	    	paldTo = NULL;
                pcTemp = NULL;
		pcTemp = (CharPtr) MemNew((size_t) (StringLen(pcChain) + StringLen(pcC) + 4));
                if (!pcTemp) goto byelater;
                StringCpy(pcTemp, pcChain);
                StringCat(pcTemp, " ");
                StringCat(pcTemp, pcC);
	    	pdnmgTo = GetLastGraph(pmmdThis);
		if (!pdnmgTo) goto byelater;
		pmadTo = GetMainAtom(pdnmgTo);
	    	paldTo = GetAtomLocs(pmadTo,  iModel);
		while (paldTo == NULL)
		  {
		      pdnmgTo = GetPrevGraph(pdnmgTo);
		      if (!pdnmgTo) goto byelater;
		      pmadTo = GetMainAtom(pdnmgTo);
		      paldTo = GetAtomLocs(pmadTo,  iModel);
		  }
	    	pdnmgFrom = GetPrevGraph(pdnmgTo);
		if (!pdnmgFrom) goto byelater;
	    	pmadFrom = GetMainAtom(pdnmgFrom);
	    	paldFrom = GetAtomLocs(pmadFrom,  iModel);
	    	RenderLabel(pmsdThis->pdnmsLink, pcTemp,
			  NULL, paldFrom, paldTo,
			  &(prkNew->Color), prkNew->LJust, prkNew->LScale,
			  RL_EXTRAPOL);
		MemFree(pcTemp);
	      }
    }



    prkNew->NodeType = (Byte) AM_MGD;
    prkNew->NodeWhat = (Byte) CONVERT_ALL;
    iTest = TraverseOneModel(pmmdThis->pdnmgHead, TRAVERSE_GRAPH, iModel,
                  iIndex,  prkNew,  RenderGraph);



  byelater:
     FreeRenderKeep(prkNew);
     return;
}




static void SwitchRender(PDNMS pdnmsThis, Int2 iModel,
			   Byte bMolecule, Int2 iBackType)
{
    PRK prKeep = NULL;
    PVNMB pvnmbThis;
    PMBD pmbdThis;
    PMSD pmsdThis;
    PARS pars = NULL;

    Int2 iTest;
    if (!pdnmsThis) return;
    
    if(Cn3d_ParsNow == 0){
       pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
       pars = GetAlgorRenderSet(pdnmsThis);
       if (!pars) return;
    }
    else {
        pars = GetThisSpecialAlgorRenderSet(Cn3d_ParsNow);
        if(!pars) return;
    }

    prKeep = NewRenderKeep();
    if (!prKeep) return;
    if (iBackType == CONNECTON)
      {  /* do inter-molecule bonds */
	  Cn3D_CopyColorCell(&(prKeep->Color), &Cn3d_PaletteRGB[Cn3d_ColorNow]); /* yellow */
      /* this really needs to be fixed.  Cn3d_ColorNow should only hold the *type* of coloring */
	  prKeep->Bond = (Byte) HALF_BOND;
	  prKeep->Atom = (Byte) ATOM_NONE;
	  switch (Cn3d_RenderNow)
	    {
	         case R_STICK:
		    prKeep->BondWidth = HET_BOND_WIDTH;
		   break;
		 case R_BALLNSTICK:
		 case R_THICKWIRE:
		    prKeep->BondWidth = VIRT_BOND_WIDTH;
	    	   break;
	         case R_DEFAULT:
		 case R_WIRE:
		 default:
		    prKeep->BondWidth = (FloatLo)0;
	    }
	   pvnmbThis = pmsdThis->pvnmbIMBHead;
  	   while (pvnmbThis)  /* walk the inter-mol bond list by hand */
            {
	     pmbdThis = (PMBD) pvnmbThis->data.ptrvalue;
             RenderAllAtom((PFB)pmbdThis->pmadFrom, iModel, iBackType, prKeep);
	     RenderAllAtom((PFB)pmbdThis->pmadTo,  iModel,  iBackType, prKeep);
	     pvnmbThis = pvnmbThis->next;
	    }
      }
    else
    switch (Cn3d_ColorNow)
		{
		    case C_CPK:
		    case C_BYTEMP:
		       Cn3D_CopyColorCell(&(prKeep->Color), &Cn3d_PaletteRGB[0]);
		       switch (Cn3d_RenderNow)
		         {
			     case R_SPACE:
			   	prKeep->Bond = (Byte) NO_BOND;
				prKeep->BondWidth = (FloatLo)0;
				prKeep->Atom = (Byte) ATOM_SPACE;
				break;
			     case R_STICK:
			   	prKeep->Bond = (Byte) HALF_BOND;
				prKeep->BondWidth = HET_BOND_WIDTH;
				prKeep->Atom = (Byte) ATOM_ISBOND;
				break;
			     case R_BALLNSTICK:
			   	prKeep->Bond = (Byte) HALF_BOND;
				prKeep->BondWidth = HET_BOND_WIDTH;
				prKeep->Atom = (Byte) ATOM_2XBOND;
				break;
			     case R_THICKWIRE:
				prKeep->Bond = (Byte) HALF_BOND;
				prKeep->BondWidth = VIRT_BOND_WIDTH;
				prKeep->Atom = (Byte) ATOM_ISBOND;
				break;
			     case R_DEFAULT:
			     case R_WIRE:
			     default:
			   	prKeep->Bond = (Byte) HALF_BOND;
				prKeep->BondWidth = (FloatLo)0;
				prKeep->Atom = (Byte) ATOM_NONE;

			 }
			prKeep->NodeType = (Byte) AM_MMD;
			prKeep->NodeWhat = (Byte) bMolecule;
		       	iTest = TraverseOneModel(pdnmsThis, TRAVERSE_MOLECULE,
					iModel, (Int4) iBackType , prKeep,
                    RenderMolecule);

		       break;
		    case C_BYRES:
			case C_BYALIGN:
      case C_BYCONS:
		    case C_BYLEN:
		    case C_BYDOMAIN:
	            case C_BYSSTRU:
		    case C_BYHYDRO:
		    case C_BYCHAIN:
		    default:
		       Cn3D_CopyColorCell(&(prKeep->Color), &Cn3d_PaletteRGB[Cn3d_ColorNow]);
                       /* if (Cn3d_ColorNow < C_top) 
                          prKeep->Color = Cn3d_ColorNow; */
		       switch (Cn3d_RenderNow)
		         {

			     case R_NAME:
			     case R_NUMBER:
			     case R_PDBNUMBER:
				prKeep->Bond = (Byte) NO_BOND;
			        prKeep->BondWidth = (FloatLo) 0;
			        prKeep->Atom = (Byte) ATOM_NONE;
                                switch (iBackType)
				 {
				    case PBBLABELSON:
					prKeep->LJust = pars->PBBLabelJust;
				        prKeep->LStyle = pars->PBBLabelStyle;
					prKeep->LScale = pars->PBBLabelScale;
					break;
				    case NTBBLABELSON:
					prKeep->LJust = pars->NTBBLabelJust;
				        prKeep->LStyle = pars->NTBBLabelStyle;
					prKeep->LScale = pars->NTBBLabelScale;
					break;
			            case NTTERMLABON:
				        prKeep->LJust = pars->NTTermLabelJust;
				        prKeep->LStyle = pars->NTTermLabelStyle;
					prKeep->LScale = pars->NTTermLabelScale;
					break;
				    case PTERMLABON:
				        prKeep->LJust = pars->PTermLabelJust;
				        prKeep->LStyle = pars->PTermLabelStyle;
					prKeep->LScale = pars->PTermLabelScale;
				 }
				break;
			     case R_SPACE:
			   	prKeep->Bond = (Byte) NO_BOND;
				prKeep->BondWidth = (FloatLo)0;
				prKeep->Atom = (Byte) ATOM_SPACE;
				break;
			     case R_STICK:
			   	prKeep->Bond = (Byte) HALF_BOND;
				prKeep->BondWidth = HET_BOND_WIDTH;
				prKeep->Atom = (Byte) ATOM_ISBOND;
				break;
			     case R_BALLNSTICK:
			   	prKeep->Bond = (Byte) HALF_BOND;
				prKeep->BondWidth = HET_BOND_WIDTH;
				prKeep->Atom = (Byte) ATOM_2XBOND;
				break;
			     case R_THICKWIRE:
			 	prKeep->Bond = (Byte) HALF_BOND;
				prKeep->BondWidth = VIRT_BOND_WIDTH;
				prKeep->Atom = (Byte) ATOM_ISBOND;
				break;
			     case R_DEFAULT:
			     case R_WIRE:
			     default:
			   	prKeep->Bond = (Byte) HALF_BOND;
				prKeep->BondWidth = (FloatLo)0;
				prKeep->Atom = (Byte) ATOM_NONE;

			 }
		       prKeep->NodeType = (Byte) AM_MMD;
		       prKeep->NodeWhat = (Byte) bMolecule;
		       iTest = TraverseOneModel(pdnmsThis, TRAVERSE_MOLECULE,
					iModel, (Int4) iBackType , prKeep,
                    RenderMolecule);
	        }  /* switch Cn3d_ColorNow */

  FreeRenderKeep(prKeep);
  return;
}



void LIBCALL MakePaletteModel(PDNMS pdnmsThis, Int2 iModel )
{
  
  PMSD pmsdThis = NULL;
  PVNMO pvnmoThis = NULL;
  PARS pars = NULL;
    ResidueColorCell iColor;
  /* fetch the active structure */
  /*pdnmsThis = GetSelectedModelstruc();*/

  if (!pdnmsThis) return;

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
  if(Cn3d_ParsNow == 0){
     pars = GetAlgorRenderSet(pdnmsThis);
     if (!pars) return;  /* cannot render */
  }
  else {
     pars = GetThisSpecialAlgorRenderSet(Cn3d_ParsNow);
     if(!pars) return;
  }

/*parsColor = pars;   */   /* this has to be changed since we have multiple */
                          /* pars now, instead parsColor will be accocated */
                          /* at the beginning of the redraw process, and */
                          /* freed at the end of the redraw process */

#ifdef DEBUG_N
printf("traversing to make palette\n");
#endif
	    if (pars->ObjectOn)
            {
	  	Cn3d_RenderNow = pars->ObjectRender;
	  	Cn3d_ColorNow = pars->ObjectColor;
          	pvnmoThis = pmsdThis->pvnmoHead;
	  	while (pvnmoThis)
  	    	{
		  if ((Int2)pvnmoThis->choice == iModel)
		    {
		     	RenderObject(pvnmoThis);
		    }
		  pvnmoThis = pvnmoThis->next;
		}
	    }

	  Cn3d_DoHydrogens = pars->HydrogensOn;

	  if ((pars->PVirtualBBOn) || (pars->PRealBBOn) || (pars->PExtraBBOn ))
	    {
	      Cn3d_RenderNow = pars->PBBRender;
      	      Cn3d_ColorNow = pars->PBBColor;
	      if (Cn3d_ColorNow == C_CPK)
	        {
		   /* assume C O N H are protein backbone colors used */
              if (Cn3d_DoHydrogens) {
                  Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ElementKinColors[1]])); /* H */
                  Cn3D_RequestColor(&Cn3D_ColorData, &iColor);
              }
		   Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ElementKinColors[6]])); /* C */
           Cn3D_RequestColor(&Cn3D_ColorData, &iColor);
		   if (pars->PRealBBOn || pars->PExtraBBOn)
		     {
		       Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ElementKinColors[7]])); /* N */
               Cn3D_RequestColor(&Cn3D_ColorData, &iColor);
		       Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[ElementKinColors[8]])); /* O */
               Cn3D_RequestColor(&Cn3D_ColorData, &iColor);
		     }
		}
	    }
	  if (pars->PVirtualBBOn)
	    {
	     /* chop each one out into a separate procedure */
	     /* pass pointers, Model Id's etc, but use LayerNow, RenderNow, ColorNow */
	     if (!(Cn3d_ColorNow == C_CPK))
		SwitchRender(pdnmsThis,iModel, AM_PROT, VIRTUALBB);
	    }
	  if (pars->PRealBBOn || pars->PExtraBBOn)
	    {
	     if (!((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)))
	       SwitchRender(pdnmsThis, iModel, AM_PROT,
			(Byte)(pars->PExtraBBOn ? REALXTRABB : REALBB) );
	    }
	  if (pars->PResiduesOn)
	    {
	       Cn3d_RenderNow = pars->PResRender;
	       Cn3d_ColorNow = pars->PResColor;

	       SwitchRender(pdnmsThis, iModel, AM_PROT,  RESIDUES);
	    }
	  if ((pars->NTVirtualBBOn) || (pars->NTRealBBOn) ||  (pars->NTExtraBBOn))
	    {
	      Cn3d_RenderNow = pars->NTBBRender;
 	      Cn3d_ColorNow = pars->NTBBColor;
	    }
	  if (pars->NTVirtualBBOn)
	    {
		SwitchRender(pdnmsThis, iModel, (Byte) (AM_RNA | AM_DNA), VIRTUALBB);
	    }
	  if (pars->NTRealBBOn || pars->NTExtraBBOn)
	    {
		SwitchRender(pdnmsThis, iModel, (Byte)(AM_RNA | AM_DNA),
			(Byte)(pars->NTExtraBBOn ? REALXTRABB : REALBB) );
	    }
	  if (pars->NTResiduesOn)
	    {
	      Cn3d_RenderNow = pars->NTResRender;
	      Cn3d_ColorNow = pars->NTResColor;
   	      SwitchRender(pdnmsThis, iModel, (Byte)(AM_RNA | AM_DNA), RESIDUES);
	    }
	  if (pars->HeterogensOn)
	    {
	  	Cn3d_RenderNow = pars->HetRender;
	  	Cn3d_ColorNow = pars->HetColor;
		SwitchRender(pdnmsThis, iModel, (Byte) (AM_HET),  HETSON);
                SwitchRender(pdnmsThis, iModel, (Byte) (AM_POLY), HETSON);
	    }
	  if (pars->IonsOn)
	    {
		Cn3d_RenderNow = pars->IonRender;
	  	Cn3d_ColorNow = pars->IonColor;
	    	SwitchRender(pdnmsThis, iModel, (Byte)(AM_ION),  IONSON);
	    }
	  if (pars->ConnectOn)
 	   {
	        Cn3d_RenderNow = pars->ConnectRender;
	        Cn3d_ColorNow = pars->ConnectColor;
	 	SwitchRender(pdnmsThis, iModel,  0,  CONNECTON);
	    }
	  if (pars->SolventOn)
	    {
	        Cn3d_RenderNow =  pars->SolventRender;
	        Cn3d_ColorNow =  pars->SolventColor;
		SwitchRender(pdnmsThis, iModel, (Byte)(AM_SOL | AM_WAT), SOLVENTON);
	    }
	  if (pars->PBBLabelOn > 1)
            {
                Cn3d_RenderNow =  R_NAME;
	        Cn3d_ColorNow =   pars->PResColor;
	        Cn3d_LabelNow =   pars->PBBLabelOn;
          	if (((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)
                  || (pars->PBBLabelStyle & L_WHITE)))
                   Cn3d_ColorNow = C_white;
		SwitchRender(pdnmsThis, iModel, (Byte)(AM_PROT), PBBLABELSON);
            }
          if (pars->NTBBLabelOn > 1)
            {
                Cn3d_RenderNow =  R_NAME;
	        Cn3d_ColorNow =   pars->NTResColor;
		Cn3d_LabelNow =   pars->NTBBLabelOn;
		if (((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)
                  || (pars->NTBBLabelStyle & L_WHITE)))
                   Cn3d_ColorNow = C_white;
		SwitchRender(pdnmsThis, iModel, (Byte)(AM_DNA | AM_RNA), NTBBLABELSON);
            }
	  if (pars->PTermLabelOn)
	    {
		Cn3d_RenderNow =  R_NAME;
	        Cn3d_ColorNow =   pars->PBBColor;
		if (((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)
                  || (pars->PTermLabelStyle & L_WHITE)))
                   Cn3d_ColorNow = C_white;
		SwitchRender(pdnmsThis, iModel, (Byte)(AM_PROT), PTERMLABON);
	    }
	  if (pars->NTTermLabelOn)
	    {
		Cn3d_RenderNow =  R_NAME;
	        Cn3d_ColorNow =   pars->NTBBColor;
		if (((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)
                  || (pars->NTTermLabelStyle & L_WHITE)))
                   Cn3d_ColorNow = C_white;
		SwitchRender(pdnmsThis, iModel, (Byte)(AM_DNA | AM_RNA), NTTERMLABON);
	    }


return;
}


static void LIBCALLBACK DoGraphicNull(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PMAD pmadAtom;
  PALD paldLoc;

     if (IsAtomNode(pfbThis))
      {
         pmadAtom = (PMAD) pfbThis;
         paldLoc = GetAtomLocs(pmadAtom, iModel);
          while (paldLoc)
            {
             paldLoc->pGraphic = NULL;
             paldLoc = paldLoc->next; /* get next location */
            }
      }
}

void LIBCALL fnMSPLoop (PDNMS pdnmsThis)
{
    PMSD pmsdThis = NULL;
    PDNML pdnmlThis = NULL;
    PMLD pmldThis = NULL;
    PARS pars = NULL;
    
    
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
    
    
    pdnmlThis = pmsdThis->pdnmlModels;
    
    /* set up for doing one model or animation */
    while (pdnmlThis)
    {
        pmldThis = (PMLD) pdnmlThis->data.ptrvalue;
        if (pmldThis->bSelected & (Byte) 0x01)
        {
            Cn3d_ParsNow = 0;
            pars = GetAlgorRenderSet(pdnmsThis);
            if (!pars) return; 
            MakePaletteModel( pdnmsThis, pdnmlThis->choice );
            
            sfpThis = sfp;
            while(sfpThis){
                sfipThis = sfpThis->data.ptrvalue;
                if(sfipThis->On) {
                    Cn3d_ParsNow = sfpThis->choice;
                    
                    pars = GetThisSpecialAlgorRenderSet(Cn3d_ParsNow);
                    if(!pars) return;
                    MakePaletteModel( pdnmsThis, pdnmlThis->choice );
                }
                sfpThis = sfpThis->next;
            }
            
        }
        pdnmlThis = pdnmlThis->next;
    }
    return;
}


void LIBCALL MakeStrucPalette(void)
{

    PARS pars = NULL;
    PMSD pmsdThis = NULL;
    PDNMS pdnmsThis = NULL;
    PDNMS pdnmsThisSlave = NULL;
    Int2 i;

#ifndef _OPENGL
    pic = NULL;  /* this should not acesss pic ! - just a color traversing pass */
#endif
    pdnmsThis = GetSelectedModelstruc();
    if (!pdnmsThis) return;
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
    
    pars = GetAlgorRenderSet(pdnmsThis);
    if (!pars) return;  /* cannot render */
    
    /* clear out the static list of colors used */
    Cn3D_FreePalette(&(Cn3D_ColorData.Palette));
    Cn3d_ColorPass = TRUE;
    Cn3d_lSlaveNum = 0;  /* go over the master */
    
    fnMSPLoop (pdnmsThis);
    
    if(AreNeighborsOn())
    {
        pdnmsThisSlave = pmsdThis->pdnmsSlaves; /* go through the slave structures */
        while(pdnmsThisSlave) 
        {
            Cn3d_lSlaveNum++; 
            if(((PMSD)(pdnmsThisSlave->data.ptrvalue))->bVisible) fnMSPLoop(pdnmsThisSlave);
            pdnmsThisSlave = pdnmsThisSlave->next;
        }
    }
    
    Cn3d_ColorPass = FALSE;
    return;
}


#ifndef _OPENGL
Picture3D LIBCALL Do3DOrigin(Picture3D p3d)
{
    FloatLo fFrom[3];
    FloatLo fTo[3];
    Int4  iXFrom, iYFrom, iZFrom;
    Int4  iXTo, iYTo, iZTo;
    CharPtr pclabel[] = {"X", "Y",  "Z"};
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    PARS parsColor = NULL;
    Int4 i;
    ResidueColorCell iColor;

    /* appends origin onto p3d */

    if (p3d == NULL)
      {  /* origin forced because no other primitives */

	 Cn3d_LayerNow = 0;
	 SetLayerTop3D(0);
         Cn3d_ColorPass = FALSE;
      }

    pdnmsThis = GetSelectedModelstruc();
    if (!pdnmsThis) return NULL;
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
    Cn3D_CopyColorCell(&iColor, &(Cn3d_PaletteRGB[C_white]));
    if (Cn3d_ColorPass)
          {
            Cn3D_RequestColor(&Cn3D_ColorData, &iColor); return NULL;
          }


        for (i = 0; i<3; i++)
            {
               fFrom[i] = 0.0;
               fTo[i] = 0.0;
            }
         for (i = 0; i<3; i++)
            {
               fFrom[i] = 5.0;
               fTo[i] = -5.0;
               RotTransScale( pmsdThis,    fFrom[0],
					   fFrom[1],
					   fFrom[2],
					   &iXFrom,
					   &iYFrom,
					   &iZFrom);
               RotTransScale( pmsdThis,    fTo[0],
					   fTo[1],
					   fTo[2],
					   &iXTo,
					   &iYTo,
					   &iZTo);
	        AddLine3D (p3d, NULL, (BigScalar) 0, Cn3d_LayerNow,
		Cn3D_ColorIndex(&Cn3D_ColorData, &iColor), iXFrom,iYFrom,iZFrom,
                                        iXTo,iYTo,iZTo);
		AddText3D (p3d, NULL, (BigScalar) 0, Cn3d_LayerNow, Cn3D_ColorIndex(&Cn3D_ColorData, &iColor),
                 pclabel[i],  iXFrom, iYFrom , iZFrom , (Uint4) (4*VIEWSCALE/2) , 0, 0);
               fFrom[i] = 0.0;
               fTo[i] = 0.0;
	    }

   return p3d;
}
#endif /* ndef _OPENGL */


static void fnARModelLoop (PDNMS pdnmsThis, PDNML pdnmlThis, PARS pars)
{
    PMSD pmsdThis = NULL;
    PMLD pmldThis = NULL;
    PVNMO pvnmoThis = NULL;
    Int2  iModel = 0;
    Int2  iTest = 0;
    Int2 i = 0;
    
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
    pmldThis = (PMLD) pdnmlThis->data.ptrvalue;
    
    iModel = (Int2) pdnmlThis->choice;
    SetActiveModel(((PFB) pmsdThis), iModel);
    /* do the object traversal */
    if (pars->ObjectOn)
    {
        Cn3d_RenderNow = pars->ObjectRender;
        Cn3d_ColorNow = pars->ObjectColor;
        pvnmoThis = pmsdThis->pvnmoHead;
        while (pvnmoThis)
        {
            if ((Int2)pvnmoThis->choice == iModel)
            {
                RenderObject(pvnmoThis);
            }
            pvnmoThis = pvnmoThis->next;
        }
    }
    
    Cn3d_DoHydrogens = pars->HydrogensOn;
    
    if ((pars->PVirtualBBOn) || (pars->PRealBBOn) || (pars->PExtraBBOn) )
    {
        Cn3d_RenderNow = pars->PBBRender;
        Cn3d_ColorNow = pars->PBBColor;
    }
    if (pars->PVirtualBBOn)
    {
        /* chop each one out into a separate procedure */
        /* pass pointers, Model Id's etc, but use LayerNow, RenderNow, ColorNow */
        SwitchRender(pdnmsThis,iModel, AM_PROT, VIRTUALBB);
    }
    if (pars->PRealBBOn || pars->PExtraBBOn)
    {
        SwitchRender(pdnmsThis, iModel, AM_PROT,
            (Byte)(pars->PExtraBBOn ? REALXTRABB : REALBB) );
    }
    if (pars->PResiduesOn)
    {
        Cn3d_RenderNow = pars->PResRender;
        Cn3d_ColorNow = pars->PResColor;
        SwitchRender(pdnmsThis, iModel, AM_PROT,  RESIDUES);
    }
    if ((pars->NTVirtualBBOn) || (pars->NTRealBBOn) ||  (pars->NTExtraBBOn))
    {
        Cn3d_RenderNow = pars->NTBBRender;
        Cn3d_ColorNow = pars->NTBBColor;
    }
    if (pars->NTVirtualBBOn)
    {
        SwitchRender(pdnmsThis, iModel, (Byte) (AM_RNA | AM_DNA), VIRTUALBB);
    }
    if (pars->NTRealBBOn || pars->NTExtraBBOn)
    {
        SwitchRender(pdnmsThis, iModel, (Byte)(AM_RNA | AM_DNA),
            (Byte)(pars->NTExtraBBOn ? REALXTRABB : REALBB) );
    }
    if (pars->NTResiduesOn)
    {
        Cn3d_RenderNow = pars->NTResRender;
        Cn3d_ColorNow = pars->NTResColor;
        SwitchRender(pdnmsThis, iModel, (Byte)(AM_RNA | AM_DNA), RESIDUES);
    }
    if (pars->HeterogensOn)
    {
        Cn3d_RenderNow = pars->HetRender;
        Cn3d_ColorNow = pars->HetColor;
        SwitchRender(pdnmsThis, iModel, (Byte) (AM_HET), HETSON);
        SwitchRender(pdnmsThis, iModel, (Byte) (AM_POLY), HETSON);
    }
    if (pars->IonsOn)
    {
        Cn3d_RenderNow = pars->IonRender;
        Cn3d_ColorNow = pars->IonColor;
        SwitchRender(pdnmsThis, iModel, (Byte)(AM_ION),  IONSON);
    }
    if (pars->ConnectOn)
    {
        Cn3d_RenderNow = pars->ConnectRender;
        Cn3d_ColorNow = pars->ConnectColor;
        SwitchRender(pdnmsThis, iModel,  0,  CONNECTON);
    }
    if (pars->SolventOn)
    {
        Cn3d_RenderNow =  pars->SolventRender;
        Cn3d_ColorNow =  pars->SolventColor;
        SwitchRender(pdnmsThis, iModel, (Byte)(AM_SOL | AM_WAT), SOLVENTON);
    }
    if (pars->PBBLabelOn > 1)
    {
        Cn3d_RenderNow =  R_NAME;
        Cn3d_ColorNow =   pars->PResColor;
        Cn3d_LabelNow =   pars->PBBLabelOn;
        if (((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)
            || (pars->PBBLabelStyle & L_WHITE)))
            Cn3d_ColorNow = C_white;
        SwitchRender(pdnmsThis, iModel, (Byte)(AM_PROT), PBBLABELSON);
    }
    if (pars->NTBBLabelOn > 1)
    {
        Cn3d_RenderNow =  R_NAME;
        Cn3d_ColorNow =   pars->NTResColor;
        Cn3d_LabelNow =   pars->NTBBLabelOn;
        if (((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)
            || (pars->NTBBLabelStyle & L_WHITE)))
            Cn3d_ColorNow = C_white;
        SwitchRender(pdnmsThis, iModel, (Byte)(AM_DNA | AM_RNA), NTBBLABELSON);
    }
    if (pars->PTermLabelOn)
    {
        Cn3d_RenderNow =  R_NAME;
        Cn3d_ColorNow =   pars->PBBColor;
        if (((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)
            || (pars->PTermLabelStyle & L_WHITE)))
            Cn3d_ColorNow = C_white;
        SwitchRender(pdnmsThis, iModel, (Byte)(AM_PROT), PTERMLABON);
    }
    if (pars->NTTermLabelOn)
    {
        Cn3d_RenderNow =  R_NAME;
        Cn3d_ColorNow =   pars->NTBBColor;
        if (((Cn3d_ColorNow == C_CPK) || (Cn3d_ColorNow == C_BYTEMP)
            || (pars->NTTermLabelStyle & L_WHITE)))
            Cn3d_ColorNow = C_white;
        SwitchRender(pdnmsThis, iModel, (Byte)(AM_DNA | AM_RNA), NTTERMLABON);
    }
}


void LIBCALL fnARLoop (PDNMS pdnmsThis )
{  
    /* the rendering engine  */
    PDNML pdnmlThis = NULL;
    PARS pars = NULL;
    Int2 iTest;
    PMSD pmsdThis = NULL;
    PMLD pmldThis = NULL;
    
    
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
    
    iTest = TraverseModels(pdnmsThis, TRAVERSE_ATOM, 0, NULL, DoGraphicNull);
    
    pdnmlThis = pmsdThis->pdnmlModels;
    while (pdnmlThis)
    {
        pmldThis = (PMLD) pdnmlThis->data.ptrvalue;
        if (pmldThis->bSelected & (Byte) 0x01)
        {  /* draw it */
#ifdef _OPENGL
            OGL_Start(Cn3D_ColorData.OGL_Data, Cn3d_LayerNow);
#endif /* _OPENGL */
            
            Cn3d_ParsNow = 0;
            pars = GetAlgorRenderSet(pdnmsThis);
            if (!pars) return; 
            fnARModelLoop (pdnmsThis, pdnmlThis, pars);
            
            sfpThis = sfp;
            while(sfpThis){
                sfipThis = sfpThis->data.ptrvalue;
                if(sfipThis->On) {
                    Cn3d_ParsNow = sfpThis->choice;
                    
                    pars = GetThisSpecialAlgorRenderSet(Cn3d_ParsNow);
                    if(!pars) return;
                    fnARModelLoop (pdnmsThis, pdnmlThis, pars);
                }
                sfpThis = sfpThis->next;
            }
            
#ifdef _OPENGL
            OGL_End();
#endif /* _OPENGL */
            Cn3d_LayerNow++;
        }  /* if (pmldThis->bSelected & (Byte) 0x01) */
        pdnmlThis = pdnmlThis->next;
    } /* while(pdnmlThis) */
}


#ifdef _OPENGL
void LIBCALL AlgorithmicRendering(void)
#else
Picture3D LIBCALL AlgorithmicRendering(Picture3D p3d)
#endif
{ 
    PDNMS pdnmsThis = NULL, pdnmsThisSlave = NULL;
    PARS pars = NULL;
    PMSD pmsdThis = NULL;
    Int2 Cn3d_LayerNow_sofar = 0;
    
    Int2 i = 0;

#ifndef _OPENGL
    pic = p3d;  /* pic is the global 3d picture for algorend.c */
#endif
    Cn3d_LayerNow = 0;
    Cn3d_AnyPrim = FALSE;  /* has any primitive been drawn? */
    Cn3d_ColorPass = FALSE;
    Cn3d_lSlaveNum = 0;  /* go over the master */
    /* fetch the active structure */
    pdnmsThis = GetSelectedModelstruc();
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
    
    if (!pdnmsThis) return NULL;
    
    pars = GetAlgorRenderSet(pdnmsThis);
    if (!pars) return NULL;  /* cannot render */
    
#ifdef _OPENGL
      OGL_SetLayers(Cn3D_ColorData.OGL_Data, FALSE);  /* null all the layers out */
#endif /* _OPENGL */

    if(pmsdThis->bVisible) fnARLoop(pdnmsThis);
    
    if(AreNeighborsOn())
    {
        pdnmsThisSlave = pmsdThis->pdnmsSlaves; /* go through the slave structures */
        
        while(pdnmsThisSlave) 
        {
            Cn3d_lSlaveNum++;
            if(((PMSD)(pdnmsThisSlave->data.ptrvalue))->bVisible) {
                fnARLoop(pdnmsThisSlave);
            }
        pdnmsThisSlave = pdnmsThisSlave->next;
        }
    }
    
    
#ifdef _OPENGL
    OGL_SetLayerTop3D(Cn3D_ColorData.OGL_Data, Cn3d_LayerNow - 1);
    return;
#else /* _OPENGL */
    SetLayerTop3D(Cn3d_LayerNow - 1);
    if (Cn3d_AnyPrim) return p3d; /* has a primitive been drawn? */
    else {
        if (p3d != NULL)
            DeletePicture3D ( p3d );
        return NULL;
    }
#endif /* else _OPENGL */

}  


