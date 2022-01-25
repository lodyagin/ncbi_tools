/*   algorend.h
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
* File Name:  algorend.h
*
* Author:  Christopher Hogue
*
* Version Creation Date:   1/31/96
*
* File Description: algorithmic rendering structures
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: algorend.h,v $
* Revision 6.12  1998/11/04 00:06:23  ywang
* add function for modeling: change render/color for special residue(s)
*
 * Revision 6.11  1998/10/28  19:29:03  ywang
 * add C_BYSEQCONS macro
 *
 * Revision 6.10  1998/10/28  19:02:07  kans
 * added two prototypes
 *
* Revision 6.8  1998/06/16 18:00:28  lewisg
* moved rendering menus and created a reset presentation menu item
*
* Revision 6.7  1998/05/26 21:35:19  lewisg
* added defaults to render menu, got rid of mouse 3D actions menu item
*
* Revision 6.6  1998/04/27 17:50:06  lewisg
* added color by conservation
*
* Revision 6.5  1998/04/20 22:09:02  lewisg
* got rid of dead code
*
* Revision 6.4  1998/04/15 03:06:14  lewisg
* get rid of dos line breaks
*
* Revision 6.3  1998/04/01 23:26:13  lewisg
* added new startup mode + fixed slave rendering
*
* Revision 6.2  1998/03/06 23:19:22  lewisg
* codewarrior fixes
*
* Revision 6.1  1998/03/06 01:16:58  lewisg
* merge
*
* Revision 6.0  1997/08/25 18:13:26  madden
* Revision changed to 6.0
*
* Revision 5.2  1996/07/22 00:24:10  hogue
* Added an origin 3D item for no-primitives condition and general use.
*
 * Revision 5.1  1996/06/03  21:21:26  hogue
 * Made tubes bigger so they are less likely to look ball-n-stick like
 * than before.
 *
 * Revision 5.0  1996/05/28  14:05:44  ostell
 * Set to revision 5.0
 *
 * Revision 1.11  1996/05/22  21:46:48  hogue
 * Added white button to label controls.
 *
 * Revision 1.10  1996/05/22  20:47:01  hogue
 * Removed HetLabel variables
 *
 * Revision 1.9  1996/05/22  15:56:55  hogue
 * Altered the label structures to make them more useful.
 *
 * Revision 1.8  1996/05/14  15:19:14  hogue
 * Added LabelControls
 *
 * Revision 1.7  1996/05/09  18:33:28  vakatov
 * included <viewer3d.h> to get know the CAMERA_SIZE_I4 actual value
 *
 * Revision 1.6  1996/05/09  15:40:40  hogue
 * Domain rendering enabled.
 *
 * Revision 1.5  1996/04/26  18:41:47  vakatov
 * CN3D sources ported to MS-Windows;
 * the portability errors and warnings fixed, etc.
 *
* ==========================================================================
*/
  
#ifndef _ALGOREND_
#define _ALGOREND_ 1

#include <viewer3d.h>
#include <cn3dmain.h>


#ifdef __cplusplus
extern "C" {
#endif  


/***ASN.1 & ANNMM compatible values for Rendering*****/

#define R_DEFAULT 0
#define R_WIRE  1
#define R_SPACE 2
#define R_STICK 3
#define R_BALLNSTICK 4
#define R_THICKWIRE 5
#define R_NAME 10
#define R_NUMBER 11
#define R_PDBNUMBER 12

#define C_hotpink  1
#define C_magenta  2
#define C_purple   3
#define C_blue     4
#define C_sky      5
#define C_cyan     6
#define C_sea      7
#define C_green    8
#define C_yellow   9
#define C_gold    10
#define C_orange  11
#define C_red     12
#define C_pink    13
#define C_pinktint 14
#define C_white    15
#define C_black    16
#define C_bluetint 17
#define C_greentint 18
#define C_yellowtint 19
#define C_gray    20
#define C_brown   21
#define C_top     22
#define C_CPK 225
#define C_BYCHAIN 226
#define C_BYTEMP 227
#define C_BYRES 228
#define C_BYLEN 229
#define C_BYSSTRU 230
#define C_BYHYDRO 231
#define C_BYOBJECT 246
#define C_BYDOMAIN 247
#define C_BYALIGN 248	/* color by alignment */
#define C_BYCONS 249 /* color by structure conservation */
#define C_BYSEQCONS 250 /* color by sequence conservation */

/* the number of colors available to color aligned structures */
#define NUM_SLAVES 8

/* these set bond draw styles */
#define NO_BOND 0
#define HALF_BOND 1

/* these set atom widths */
#define ATOM_NONE   0
#define ATOM_SPACE  1
#define ATOM_2XBOND 2
#define ATOM_ISBOND 3

#define HET_BOND_WIDTH   (float)0.4
#define VIRT_BOND_WIDTH  (float)0.6
#define SUPER_BOND_WIDTH (float)1.0
#define CYL_THRESHOLD    (float)0.1
#define EXPAND_ATOM      (float)1.8

#define CONNECTON  0
#define VIRTUALBB  1
#define REALBB     2
#define REALXTRABB 3
#define RESIDUES   4
#define IONSON     5
#define HETSON     6
#define SOLVENTON  7
#define PBBLABELSON  8
#define NTBBLABELSON 9
#define HETLABELSON 10
#define PTERMLABON  11
#define NTTERMLABON 12


/* from the ASN.1 spec ...
        default         (0),  -- Default view
        wire            (1),  -- use wireframe 
        space           (2),  -- use spacefill
        stick           (3),  -- use stick model (thin cylinders)
        ballNStick      (4),  -- use ball & stick model
        thickWire       (5),  -- thicker wireframe
        bbOnly          (6),  -- show BB atoms only
        bbVirtual       (7),  -- show BB as virtual c-alpha bonds
        hide            (9),  -- don't show this
        name            (10), -- display its name next to it
        number          (11), -- display its number next to it 
        pdbNumber       (12), -- display its PDB number next to it
        hBond           (13), -- show any hBonds this forms
        objWireFrame    (150), -- display MMDB surface object as wireframe
        objPolygons     (151), -- display MMDB surface object as polygons
        bbWire          (200), -- backbone in wireframe
        bbSpace         (201), -- backbone in spacefill
        bbStick         (202), -- backbone in stick (thin cylinders)
        bbCurve         (203), -- backbone interpreted as a curve
        bbWorm          (204), -- backbone interpreted like a worm
        bbRibbon        (205), -- backbone as a ribbon (general)
        bbLineRibbon    (206), -- backbone as a vectorized ribbon
        bbCylRibbon     (207), -- bacbbone ribbon cyl-x-section
        bbRectRibbon    (208), -- backbone ribbon rect-x-section
        bbOvalRibbon    (209), -- backbone ribbon oval-x-section
        bbBallNStick    (210), -- backbone as ball and stick
        bbThickWire     (211), -- backbone in thick wireframe
        ssCylinder      (215), -- show helices as cylinders
        ssRectArrow     (216), -- show sheets with rect-x-section arrows
        ssSquareArrow   (217), -- show sheets with square-x-section arrows
        ssOvalArrow     (218), -- show sheets with oval-x-section arrows
        colorsetCPK     (225), -- color atoms like CPK models
        colorsetbyChain (226), -- color each chain different
        colorsetbyTemp  (227), -- color using isotropic Temp factors 
        colorsetbyRes   (228), -- color using residue properties
        colorsetbyLen   (229), -- color changes along chain length
        colorsetbySStru (230), -- color by secondary structure
        colorsetbyHydro (231), -- color by hydrophobicity
        colorsetTurn    (232), -- color Turn
        colorsetLoop    (233), -- color Loop
        colorsetHelix   (234), -- color Helix
        colorsetSheet   (235), -- color Sheet
        colorsetSite    (236), -- color Site
        colorsetHet     (237), -- color Het
        colorsetIon     (238), -- color metal Ion
        colorsetDNA     (239), -- color DNA
        colorsetRNA     (240), -- color RNA
        colorsetProtein (241), -- color Protein
        colorsetSugar   (242), -- color Sugar
        colorsetBasic   (243), -- color Basic pH
        colorsetAcidic  (244), -- color Acidic pH
        colorsetNeutr   (245), -- color Neutral pH
        other           (255) 

*/


/* this AlgorRenderSet gets hung off of an MSD structure for the first go-round */
/* will later feature-ize it */

#define   L_NAME      0x01
#define   L_NUM       0x02
#define   L_PDB       0x04
#define   L_WHITE     0x20
#define   L_3LETR     0x40
#define   L_1LETR     0x80


/* 3LETR is both bits set */

#define LA_LEFT   0x01
#define LA_RIGHT  0x02
#define LA_UPPER  0x04
#define LA_LOWER  0x08
#define LA_CENTER 0x20
#define LA_FRONT  0x40
 
typedef struct AlgorRenderSet {

/* this embeds a camera structure inside */

  Int4 dummy[CAMERA_SIZE_I4];
   
/* global settings */

  Boolean 	   HydrogensOn;  
  Boolean          DisorderOn;
  Boolean          AnimateOn;

/* an On Boolean signifies a separate pass through the structure data with a callback */
/* Protein Renderings */
  Boolean          PVirtualBBOn;
  Boolean          PRealBBOn;
  Boolean          PExtraBBOn;  /* carbonyl oxygens */
  Boolean          PResiduesOn;
  Int2		   PBBRender;  /* valid are {0,1-5} */
  Int2             PBBColor;  /* valid are {0,225-231} */
  Int2 		   PResRender;  /* valid are {0, 1-5} */
  Int2		   PResColor;  /* valid are {0,225-231} */


  Uint1            PBBLabelOn;  /* uses PResColor */
  Uint1            PBBLabelJust;   
  Uint1            PBBLabelStyle;
  Int2             PBBLabelScale;
  
  Boolean          PTermLabelOn; /* uses PBBColor */
  Uint1		   PTermLabelJust;
  Uint1            PTermLabelStyle;
  Int2             PTermLabelScale;
  
 
/* DNA/RNA  Renderings */
  Boolean          NTVirtualBBOn;
  Boolean          NTRealBBOn;
  Boolean          NTExtraBBOn;  /* phosponyl oxygens */
  Boolean          NTResiduesOn;
  Int2             NTBBRender; /* valid are {0,1-5} */
  Int2             NTBBColor;  /* valid are {0,225-231} */
  Int2		   NTResRender;/* valid are {0,1-5} */
  Int2 		   NTResColor; /* valid are {0,225-231} */ 
  Int2				DomainColor;
 
  Uint1            NTBBLabelOn;  /* uses NTResColor */
  Uint1            NTBBLabelJust;
  Uint1            NTBBLabelStyle;
  Int2             NTBBLabelScale;

  Boolean          NTTermLabelOn;  /* Uses NTBBColor */
  Uint1            NTTermLabelJust;
  Uint1            NTTermLabelStyle;
  Int2             NTTermLabelScale;
 
  Boolean          HeterogensOn;  
  Int2             HetRender; /* valid are {0,1-5} */
  Int2             HetColor; /* valid are {0,225,226,227} */ 

  Boolean          IonsOn;   
  Int2             IonRender;  /* valid are {0,2,4} */
  Int2             IonColor;  /* valid are {0,225,226,227} */

  Boolean          ConnectOn;    
  Int2             ConnectRender;  /* valid are {0, 1,3,4,5} */
  Int2             ConnectColor;   /* valid are {0} fixed yellow color */
  
  Boolean          SolventOn;   
  Int2             SolventRender;  /* {} from hets, 0 inherits from HetRender*/
  Int2             SolventColor;  /* {} from hets,  0 inherits fro HetColor */ 

  Boolean	   ObjectOn;
  Int2		   ObjectRender;
  Int2 	           ObjectColor;  /* valid are BYCHAIN, BYSSTRU */
  
  Uint1 NumColors;
  Uint1 IndexRGB[CN3D_COLOR_MAX];

} ARS, PNTR PARS;



/* this structure keeps data for the rendering callbacks */
typedef struct RenderKeep {
   Byte NodeWhat;
   Byte NodeType;
   Int2 Color; /* a fixed color */
   Byte Bond;  /*  use define */
   Byte Atom;
   FloatLo BondWidth; 
   Uint1 LJust;
   Uint1 LStyle;
   Int2  LScale;
} RK,  PNTR PRK;


/************function prototypes***********/

extern void  LIBCALL SetDefaultAlgorRender PROTO((PARS pars));
extern void  LIBCALL SetAlignAlgorRender PROTO((PARS pars));
extern PARS  LIBCALL NewAlgorRenderSet     PROTO((void));
extern PARS LIBCALL NewAlignRenderSet PROTO((void));
extern void  LIBCALL FreeAlgorRenderSet    PROTO((PARS pars));
extern PARS  LIBCALL GetAlgorRenderSet     PROTO((PDNMS pdnmsThis));
extern void  LIBCALL ResetRenderCtrls      PROTO((void));
extern GrouP LIBCALL RenderControls        PROTO((Nlm_GrouP prnt));
extern void  LIBCALL ResetLabelCtrls       PROTO((void));
extern GrouP LIBCALL LabelControls         PROTO((Nlm_GrouP prnt));
extern PRK   LIBCALL NewRenderKeep         PROTO((void));
extern PRK   LIBCALL CopyRenderKeep        PROTO((PRK prkThis));
extern void  LIBCALL FreeRenderKeep        PROTO((PRK prkThis));
extern void  LIBCALL RenderObject          PROTO((PVNMO pvnmoThis));
extern void  LIBCALL RenderBond   PROTO((PALD paldFrom, PALD paldTo,
                                         Int2 iColor, FloatLo fCylRadius));
extern void  LIBCALL RenderAnAtom PROTO((PALD paldAtom,
                                         Int2 iColor, FloatLo fRadius));
extern Int2  LIBCALL GetGraphNCBIstdaa     PROTO((PMGD pmgdThis));
extern Int2  LIBCALL GetGraphNCBI4na       PROTO((PMGD pmgdThis));
extern void  LIBCALL MakePaletteModel      PROTO((PDNMS pdnmsThis, Int2 iModel ));
extern void  LIBCALL MakeStrucPalette      PROTO((void)); 
extern Picture3D LIBCALL Do3DOrigin PROTO((Picture3D p3d));
extern Picture3D LIBCALL AlgorithmicRendering PROTO((Picture3D p3d));
extern void Cn3D_RedrawProc PROTO((ButtoN b));
extern void LIBCALL fnMSPLoop PROTO((PDNMS pdnmsThis));
extern void LIBCALL fnARLoop PROTO((PDNMS pdnmsThis ));

extern void Cn3D_RenStruc  PROTO ((IteM i));
extern void Cn3D_RenWire   PROTO ((IteM i));
extern void Cn3D_RenTube   PROTO ((IteM i));
extern void Cn3D_RenHier   PROTO ((IteM i));
extern void Cn3D_RenSpace  PROTO ((IteM i));
extern void Cn3D_RenBS     PROTO ((IteM i));
extern void Cn3D_RenDefault PROTO((IteM i));
extern void Cn3D_RenAlign  PROTO((IteM i));


extern void Cn3D_ColCPK    PROTO ((IteM i));
extern void Cn3D_ColDomain PROTO ((IteM i));
extern void Cn3D_ColCy     PROTO ((IteM i));
extern void Cn3D_ColStru   PROTO ((IteM i));
extern void Cn3D_ColRes    PROTO ((IteM i));
extern void Cn3D_ColHydro  PROTO ((IteM i));
extern void Cn3D_ColTemp   PROTO ((IteM i));
extern void Cn3D_ColAlign   PROTO ((IteM i));
extern void Cn3D_ColCons   PROTO ((IteM i));
extern void Cn3D_ColSeqCons PROTO ((IteM i));
extern void Cn3D_DoColSeqCons PROTO ((void));


#ifdef __cplusplus
}
#endif

#endif
