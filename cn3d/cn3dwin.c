/*  $Id: cn3dwin.c,v 6.56 1998/09/23 18:38:50 ywang Exp $
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
* Author:  Christopher Hogue
*
* Version Creation Date:   1/31/96
*
* File Description:  Cn3D GUI API
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cn3dwin.c,v $
* Revision 6.56  1998/09/23 18:38:50  ywang
* add functions to control display on domain level
*
 * Revision 6.55  1998/09/22  17:53:13  ywang
 * add menu for display control
 *
 * Revision 6.54  1998/09/17  20:51:17  ywang
 * add functions for edit background and highlight color
 *
 * Revision 6.53  1998/09/04  19:50:53  ywang
 * reorganize codes for highlighting
 *
 * Revision 6.52  1998/09/01  18:19:15  ywang
 * initialize IteM i = NULL in LaunchSequenceWindow
 *
 * Revision 6.51  1998/08/26  18:28:38  kans
 * fixed -v -fd warnings
 *
* Revision 6.50  1998/08/17 18:45:48  lewisg
* change version to 2.01
*
* Revision 6.49  1998/08/05 19:12:36  ywang
* comment out MK_Shift/MA_DClick, MK_Ctrl/MA_DClick
*
 * Revision 6.48  1998/07/23  14:34:30  chappey
 * resize of salsa window
 *
* Revision 6.47  1998/07/22 20:18:34  lewisg
* resize of salsa window
*
* Revision 6.46  1998/07/16 18:18:39  chappey
* SeqIsIn -> SeqIdForSameBioseq
*
* Revision 6.45  1998/07/13 23:19:57  ywang
* check hidden sequence whenever launching salsa
*
 * Revision 6.44  1998/07/09  21:47:37  ywang
 * reduce color messages for salsa
 *
 * Revision 6.43  1998/06/30  23:29:22  ywang
 * fix bugs regarding to read in more structures
 *
 * Revision 6.42  1998/06/30  20:49:02  ywang
 * improve performance and prepare for close salsa
 *
 * Revision 6.41  1998/06/22  18:58:07  chappey
 * moved GetAppProperty SeqEditDisplayForm to cn3dwin.c
 *
* Revision 6.40  1998/06/17 17:42:25  lewisg
* moved a menu item
*
* Revision 6.39  1998/06/16 18:00:30  lewisg
* moved rendering menus and created a reset presentation menu item
*
* Revision 6.38  1998/06/10 22:04:32  ywang
* remove obsolete code
*
 * Revision 6.37  1998/06/10  22:01:08  ywang
 * remove obsolete code
 *
 * Revision 6.36  1998/06/04  16:33:26  ywang
 * autamatially launch salsa window
 *
 * Revision 6.35  1998/05/28  22:06:00  ywang
 * maintain highlight upon rendering switch when salsa is on
 *
 * Revision 6.33  1998/05/26  22:06:13  ywang
 * salsa get cn3d color when it is launched
 *
 * Revision 6.32  1998/05/26  21:35:21  lewisg
 * added defaults to render menu, got rid of mouse 3D actions menu item
 *
* Revision 6.31  1998/05/18 22:09:14  ywang
* move codes around
*
 * Revision 6.30  1998/05/18  16:45:52  ywang
 * allocate memory for mediadata
 *
 * Revision 6.29  1998/05/14  14:53:42  ywang
 * fix bugs
 *
 * Revision 6.27  1998/05/06  23:50:27  lewisg
 * fixed launching problem with sequin
 *
* Revision 6.26  1998/05/05 20:06:36  ywang
* set yellow as highlight color in cn3d
*
 * Revision 6.25  1998/04/30  15:22:29  ywang
 * start to store Num_ActiveSlave
 *
 * Revision 6.24  1998/04/29  18:03:06  lewisg
 * new menus
 *
* Revision 6.23  1998/04/28 22:47:24  lewisg
* master/slave color in sync
*
* Revision 6.22  1998/04/28 18:51:07  ywang
* slight modification
*
 * Revision 6.21  1998/04/27  17:49:57  lewisg
 * added color by conservation
 *
* Revision 6.20  1998/04/20 16:03:27  ywang
* launch sequence viewer correctly in one biostruc/multiple chain case
*
 * Revision 6.19  1998/04/18  00:33:51  lewisg
 * added ability to turn slaves on/off
 *
* Revision 6.18  1998/04/17 01:08:13  ywang
* call LaunchAlignEditor in one pair alignment case while call LaunchAnnotAlignEditor in multiple alignment case
*
 * Revision 6.16  1998/04/15  19:22:51  chappey
 * Add REGISTER_NEW_SEQANNOT_EDIT to be able to launch multiple alignment from the SeqAnnot
 *
* Revision 6.15  1998/04/14 21:16:56  ywang
* try to pass the head of alignment data link list to salsa
*
 * Revision 6.14  1998/04/09  17:07:23  lewisg
 * added a version # and got rid of File/Vast Alignments
 *
* Revision 6.13  1998/04/03 18:06:24  ywang
* show NA sequence also
*
 * Revision 6.12  1998/04/02  23:01:51  kans
 * code warrior distinguishes 0 from NULL
 *
* Revision 6.11  1998/04/02 22:27:18  ywang
* multiple sequence viewer for multiple chain protein
*
 * Revision 6.9  1998/04/01  01:26:19  ywang
 * sequence view and get blast alignment-Colombe
 *
* Revision 6.8  1998/03/30 23:32:28  ywang
* Set Hightlight Color as Red instead of Blue; Associate simple double click action with one residue highlight
*
 * Revision 6.7  1998/03/30  16:02:11  kans
 * changed printf to ErrPostEx
 *
* Revision 6.6  1998/03/27 23:17:38  ywang
* Get psaAlignment from pmsdMaster
*
 * Revision 6.5  1998/03/26  22:27:42  kans
 * fixed CodeWarrior complaints and missing prototypes
 *
* Revision 6.4  1998/03/26 20:44:20  ywang
* start cn3d messager
*
 * Revision 6.3  1998/03/07  20:43:53  kans
 * moved Cn3D_fEntrezOn to cn3dwin.c
 *
* Revision 6.2  1998/03/06 19:03:00  kans
* needed to add two includes
*
* Revision 6.1  1998/03/06 01:22:59  lewisg
* merge
*
* Revision 6.0  1997/08/25 18:13:48  madden
* Revision changed to 6.0
*
* Revision 1.4  1997/07/29 21:17:11  vakatov
* [WIN32,DLL]  Made Cn3D's stubbed functions be DLL-exportable
*
* Revision 1.3  1997/03/31 16:53:14  vakatov
* Use Z-rotation scrollbar in the 3D-viewer controls.
* Changed the CN3D release from 1.0 to 1.1.
*
 * Revision 1.2  1997/03/20  21:21:12  vakatov
 * [WIN_MAC] Cn3DResizeProc():  take into account the window menubar height
 *
 * Revision 1.1  1997/03/20  16:24:36  vakatov
 * Initial revision
 *
 * Revision 5.18  1996/08/19  21:05:26  vakatov
 * Made functions of type "pNodeFunc" to be of that type(sharp)
 * (and thus fixed fatal bug under Borland/WIN32); removed all
 * castings to (pNodeFunc) -- to let compiler to check for such bugs
 *
 * Revision 5.17  1996/07/31  18:36:32  hogue
 * Segment highlighting code added, MIME-type code from Network Entrez
 * adapted for STANDALONE Cn3D.
 *
 * Revision 5.16  1996/07/29  21:12:44  epstein
 * add logic to hide query window when cn3d MIME viewer starts up
 *
 * Revision 5.15  1996/07/26  18:57:41  kans
 * hide and update, then show and update
 *
 * Revision 5.14  1996/07/22  00:26:37  hogue
 * Default 3D viewer made smaller, Help menu refers to WWW-site,
 * Quit Entrez menu item, and About Structure menu item which
 * displays Structure Summary + PDB remarks.  Also trapped
 * no-primitives condition.
 *
 * Revision 5.13  1996/06/14  14:54:19  vakatov
 * [WIN_MOTIF]  RestrictMotifColorsTo( 32 ) call added before creating
 * 3D-Viewer shell window -- to provide proper Motif/3D-V color sharing
 *
 * Revision 5.12  1996/06/13  21:26:45  kans
 * fixed mac resize
 *
 * Revision 5.11  1996/06/13  21:05:31  kans
 * fixed mac-specific typos
 *
 * Revision 5.10  1996/06/13  20:49:10  hogue
 * Removed Beta designation .
 *
 * Revision 5.9  1996/06/13  20:44:17  hogue
 * Hides window now when starting up in Entrez, fix for Mac window
 * resize problem.
 *
 * Revision 5.8  1996/06/13  16:36:42  kans
 * menus now in the window for all platforms, including Mac
 *
 * Revision 5.7  1996/06/12  14:31:33  hogue
 * Added Cn3DWin function for integration into Entrez
 *
 * Revision 5.6  1996/06/05  16:28:45  hogue
 * Fixed a gif file name save bug.
 *
 * Revision 5.5  1996/06/03  22:36:08  hogue
 * Rearragned top menus for ease of use & cosmetic reasons.
 *
 * Revision 5.4  1996/06/03  21:20:09  hogue
 * Changed FixedWindow to DocumentWindwo for Win/Mac resize-ability,
 * Added call to MyQuit to stop playing layers if was on.  Also
 * added an initial call to the resizer to clean up sizing problems
 * on startup.
 *
 * Revision 5.3  1996/06/03  20:00:57  kans
 * command item was missing a closing double quo;te
 *
 * Revision 5.2  1996/06/03  19:43:52  vakatov
 * Quit the program more careful
 *
 * Revision 5.1  1996/05/29  19:16:02  vakatov
 * Cn3DResizeProc():  LinkAnimCtrls3D() call added to relink 3D-controls
 * to the new 3D-viewer
 *
 * Revision 5.0  1996/05/28  14:05:44  ostell
 * Set to revision 5.0
 *
 * Revision 1.17  1996/05/23  14:33:49  hogue
 * Added Image menu with Zoom, Move, All, Save fn;
 * Added MS-Win Edit menu (Copy only) and Print fn on Image menu
 * Added call to stop animation playing for redraw.
 *
 * Revision 1.16  1996/05/22  21:23:54  hogue
 * Added Watch/Arrow cursor calls, ProgMon.
 *
 * Revision 1.15  1996/05/21  22:50:10  hogue
 * Capture camera of current structure before resize.
 *
 * Revision 1.14  1996/05/21  22:18:49  hogue
 * Added call to ResetLabelCtrls
 *
 * Revision 1.13  1996/05/21  22:12:28  vakatov
 * "Cn3DResizeProc()" rewritten and it is now able to count the menu-bar height
 * "Viewer3DGroups()" modified to better control over the groups positioning
 *
 * Revision 1.12  1996/05/14  15:45:16  hogue
 * Added UpdateColorTable call & changed palette structure accordingly,
 * Added Cn3dResizeProc and set up the Label,
 * Render & Viewer Controls and the 3D window to move around as planned-
 *  which uncovers new bugs in Viewer3d code...
 *
 * Revision 1.11  1996/05/09  15:41:10  hogue
 * Domain rendering enabled.
 *
 * Revision 1.10  1996/05/07  18:30:58  vakatov
 * Viewer3DGroups() -- slightly changed;  + casting...
 *
 * Revision 1.9  1996/04/26  21:44:21  vakatov
 * Tune the 3D-viewer size to fit the screen
 *
 * Revision 1.8  1996/04/26  18:42:24  vakatov
 * CN3D sources ported to MS-Windows;
 * the portability errors and warnings fixed, etc.
 *
 * Revision 1.7  1996/04/18  16:57:01  hogue
 * Altered color palette for multi-structure display, preparing for neighbors...
 *
 * Revision 1.6  1996/04/04  21:05:21  hogue
 * rearranged menus, fixed camera calls, added NCBI logo
 *
 * Revision 1.5  1996/03/30  23:40:19  hogue
 * Redraw now saves camera
 *
 * Revision 1.4  1996/03/29  20:00:06  hogue
 * Integrated 3d viewing, menus & controls for algorithmic rendering
 *
 * Revision 1.2  1996/02/02  19:39:32  hogue
 * Initial Revision
 *
*
* ==========================================================================
*/

#include <ncbi.h>
#include <accentr.h>
#include <cn3dmain.h>
#include <cn3dopen.h>
#include <cn3dxprt.h>
#include <cn3dwipe.h>
#include <cn3dslct.h>
#include <cn3dsave.h>
#include <algorend.h>
#include <naybor.h>
#include <objalign.h>
#include <objseq.h>
#include <objmgr.h>
#include <saledit.h>
#include <cn3dpane.h>
#include <lsqfetch.h>
#include <salutil.h>
#include <cn3dmsg.h>
#include <cn3dmodl.h>

static Uint2    Cn3D_Vy, Cn3D_Rx;
static MenU	Cn3D_sOpen;
static IteM	Cn3D_iSelStruc;
static IteM	Cn3D_iClearStruc;
static IteM     Cn3D_iRendCtrl;
static IteM     Cn3D_iLabelCtrl;
static IteM Cn3D_iAlignCtrl;
static IteM     Cn3D_iDisplayCtrl;     /* For display control, Yanli */
static IteM     Cn3D_iViewCtrl;
static MenU	Cn3D_sExport;
static WindoW 	Cn3D_w = NULL;
static MenU  	Cn3D_ma_group_menu;
static MenU  	Cn3D_ma_action_menu;
static MenU	Cn3D_sSave;
/*static MenU     Cn3D_sNaybor;*/
static MenU     Cn3D_mRender;
static MenU     Cn3D_mColor;
static MenU     Cn3D_mControls;

static Uint1   errNum;
static Int1    errType;
static CharPtr errMsg;

Viewer3D  Cn3D_v3d = NULL;
static Picture3D Cn3D_pMain;

static GrouP    Cn3D_gWinGP;
static GrouP    Cn3D_gViewer;
static GrouP    Cn3D_gRendCtrl;
static GrouP    Cn3D_gLabelCtrl;
static GrouP Cn3D_gAlignCtrl;
static GrouP    Cn3D_gDisplayCtrl;    /* For display control, Yanli  */
static GrouP    Cn3D_gViewCtrl;

static Nlm_Controls3D Cn3D_left;

Uint1 Cn3d_IndexRGB[CN3D_COLOR_MAX];
ValNodePtr Cn3d_ColorNames = NULL; /* choice holds table number */


Nlm_RGBColoR Cn3d_PaletteRGB[CN3D_COLOR_MAX] =
{
  255, 255, 255, /* default     0 */
  255,  20, 147, /* hotpink     1 */
  255,   0, 255, /* magenta     2 */
  155,	48, 255, /* purple      3 */
    0,   0, 255, /* blue        4 */
   30, 144, 255, /* sky         5 */
    0, 255, 255, /* cyan        6 */
    0, 255, 127, /* sea         7 */
    0, 255,   0, /* green       8 */
  255, 255,   0, /* yellow      9 */
  255, 165,   0, /* gold       10 */
  255,  69,   0, /* orange     11 */
  255,   0,   0, /* red        12 */
  255, 114,  86, /* pink       13 */
  255, 174, 185, /* pinktint   14 */
  255, 255, 255, /* white      15 */
    0,   0,   0, /* black      16 */
  176, 226, 255, /* bluetint   17 */
  154, 255, 154, /* greentint  18 */
  255, 236, 139, /* yellowtint 19 */
  125, 125, 125, /* gray       20 */
  139,  87,  66, /* brown      21 */
  255, 255, 255, /* user colors 22 */
  255, 255, 255, /* user colors 23 */
  255, 255, 255, /* user colors 24 */
  255, 255, 255, /* user colors 25 */
  255, 255, 255, /* user colors 26 */
  255, 255, 255, /* user colors 27 */
  255, 255, 255, /* user colors 28 */
  255, 255, 255, /* user colors 29 */
  255, 255, 255, /* user colors 30 */
  255, 255, 255, /* user colors 31 */
  255, 255, 255, /* user colors 32 */
  255, 255, 255, /* user colors 33 */
  255, 255, 255, /* user colors 34 */
  255, 255, 255, /* user colors 35 */
  255, 255, 255, /* user colors 36 */
  255, 255, 255, /* user colors 37 */
  255, 255, 255, /* user colors 38 */
  255, 255, 255, /* user colors 39 */
  255, 255, 255, /* user colors 40 */
  255, 255, 255, /* user colors 41 */
  255, 255, 255, /* user colors 42 */
  255, 255, 255, /* user colors 43 */
  255, 255, 255, /* user colors 44 */
  255, 255, 255, /* user colors 45 */
  255, 255, 255, /* user colors 46 */
  255, 255, 255, /* user colors 47 */
  255, 255, 255, /* user colors 48 */
  255, 255, 255, /* user colors 49 */
  255, 255, 255, /* user colors 50 */
  255, 255, 255, /* user colors 51 */
  255, 255, 255, /* user colors 52 */
  255, 255, 255, /* user colors 53 */
  255, 255, 255, /* user colors 54 */
  255, 255, 255, /* user colors 55 */
  255, 255, 255, /* user colors 56 */
  255, 255, 255, /* user colors 57 */
  255, 255, 255, /* user colors 58 */
  255, 255, 255, /* user colors 59 */
  255, 255, 255, /* user colors 60 */
  255, 255, 255, /* user colors 61 */
  255, 255, 255, /* user colors 62 */
  255, 255, 255  /* user colors 63 */
};


/*
Int2 LIBCALL Add3DColor(Uint1 Red, Uint1 Green, Uint1 Blue, ChatPtr name)
{

Finds next free color  slot in table.  Adds color.
Adds name to valnode list of names.

}


Boolean LIBCALL Remove3DColor(Int2 index)
{

Removes color number only if > 21

}

Boolean LIBCALL Replace3DColor(Int2 index)
{

Replaces RGB values of any named color

}

*/


extern Boolean LIBCALL readErrors(void)
{
  Uint1 i;

  errNum = DiagGetRecordCount();
  for (i = 0;  i < errNum;  i++)
    {
      errType = DiagGetRecordType( i );
      errMsg  = DiagGetRecordStr ( i );
      Message(MSG_OK, "%d %d %s", errNum, errType, errMsg);
    }

  return DiagHasErrorRec();
}


void LIBCALL Cn3D_EnableFileOps(void)
{
  if (GetFirstModelstruc() == NULL)
    { /* nothing in memory - disable stuff */
      Cn3D_DisableFileOps();
      /* just leave the "opener enabled" */
      Enable( Cn3D_sOpen );
    }
  else
    {
      Enable(Cn3D_sOpen);
      Enable(Cn3D_sSave);
      Enable(Cn3D_sExport);
      Enable(Cn3D_iSelStruc);
      Enable(Cn3D_iClearStruc);
      /* Enable(Cn3d_sNaybor); */
    }

  /*
  if (AreNeighborsOn())
    {
      DisableFileOps();
      Enable(Cn3D_sNaybor);
    }
  */

  return;
}

void LIBCALL Cn3D_DisableFileOps(void)
{
  Disable(Cn3D_sOpen);
  Disable(Cn3D_sSave);
  Disable(Cn3D_sExport);
  Disable(Cn3D_iSelStruc);
  Disable(Cn3D_iClearStruc);
  /* Disable(Cn3D_sNaybor); */
  return;
}


void LIBCALL Cn3D_DisableMenus(void)
{
  Cn3D_DisableFileOps();
  Disable(Cn3D_mRender);
  Disable(Cn3D_mColor);
  Disable(Cn3D_mControls);
}


void LIBCALL Cn3D_EnableMenus(void)
{
  Cn3D_EnableFileOps();
  Enable(Cn3D_mRender);
  Enable(Cn3D_mColor);
  Enable(Cn3D_mControls);
}


static void Cn3D_AboutProc(IteM i)
{
  MsgAlert(KEY_OK, SEV_INFO, "About Cn3D",
"Cn3D\n\nA 3-D Viewer for MMDB\nthe Molecular Modelling Database\nVersion 2.01.0000\n\nThe National Center for Biotechnology Information\ninfo@ncbi.nlm.nih.gov");
}


static void Cn3D_HelpProc(IteM i)
{
  MsgAlert(KEY_OK, SEV_INFO, "Cn3D Online Manual",
"An WWW-Based Manual for Cn3D can be viewed at:\n\n\
http://www.ncbi.nlm.nih.gov/Structure/cn3d.html\n\n\
 Send questions to:\n\
info@ncbi.nlm.nih.gov");
}


static void Cn3D_AllCB(IteM i)
{
  ZoomAll3D( Cn3D_v3d );
}

static void Cn3D_ZoomCB(IteM i)
{
  MsgAlert(KEY_OK, SEV_INFO, "Cn3D Online Manual",
           "Use the [Control] modifier key to zoom with mouse\n\n\
(see also Image/Mouse3D_Groups)");
}

static void Cn3D_MoveCB(IteM i)
{
  MsgAlert(KEY_OK, SEV_INFO, "Cn3D Online Manual",
           "Use the [Shift] modifier key to move with mouse\n\n\
(see also Image/Mouse3D_Groups)");
}

static void Cn3D_BgColor(IteM i)
{

    if(Cn3D_v3d != NULL) BgColorDlg3D(Cn3D_v3d);

}

static void Cn3D_HLColor(IteM i)
{

    Uint1 colorR, colorG, colorB;

    ChooseColorDialog(&colorR, &colorG, &colorB,0);
    SetHLColor3D(Cn3D_v3d, colorR, colorG, colorB);
}


/* Compose the starting(NCBI-Logo) 3D-Picture
 */
static void LogoProc(Nlm_ButtoN b)
{
  Nlm_Prim3D p[4];
  Int4       x1, x2;
  double     ang1, ang2;
  double     r;

  if (Cn3D_pMain != NULL)
    DeletePicture3D( Cn3D_pMain );
  Cn3D_pMain = CreatePicture3D();
  if ( readErrors() )
    return;

  ResetPicture3D( Cn3D_pMain );
  if ( readErrors() )
    return;

  AllocPalette3D(Cn3D_pMain, 4);
  if ( readErrors() )
    return;

  SetColor3D(Cn3D_pMain, 0, 255,   0, 255);
  SetColor3D(Cn3D_pMain, 1,   0,   0, 255);
  SetColor3D(Cn3D_pMain, 2, 255, 255,   0);
  SetColor3D(Cn3D_pMain, 3,   0,   0, 255);

  x1 = 0;
  x2 = 4000000;
  ang1 = 0.2;
  p[0] = AddPoly3D(Cn3D_pMain, NULL, 0, 0, 0,
                   x2, (Int4)(15000000*cos(ang1)), (Int4)(15000000*sin(ang1)),
                   x1, (Int4)(15000000*cos(ang1)), (Int4)(15000000*sin(ang1)));
  ang2  = 3.2;
  p[1] = AddPoly3D(Cn3D_pMain, NULL, 0, 0, 1,
                   x2, (Int4)(15000000*cos(ang2)), (Int4)(15000000*sin(ang2)),
                   x1, (Int4)(15000000*cos(ang2)), (Int4)(15000000*sin(ang2)));

  ang1 = 0.0;
  for (r = 1;  r <= 63;  r++, ang1 += 0.2)
    {
      x1 += 700000;
      x2 += 700000;
      AddVertPoly3D(Cn3D_pMain, p[0],
                    x2, (Int4)(15000000*cos(ang1)),(Int4)(15000000*sin(ang1)));
      AddVertPoly3D(Cn3D_pMain, p[0],
                    x1, (Int4)(15000000*cos(ang1)),(Int4)(15000000*sin(ang1)));
    }

  x1 = 0;
  x2 = 4000000;
  ang1 = 3.0;
  for (r = 1;  r <= 79;  r++, ang1 += 0.2)
    {
      x1 += 700000;
      x2 += 700000;
      AddVertPoly3D(Cn3D_pMain, p[1],
                    x2, (Int4)(15000000*cos(ang1)),(Int4)(15000000*sin(ang1)));
      AddVertPoly3D(Cn3D_pMain, p[1],
                    x1, (Int4)(15000000*cos(ang1)),(Int4)(15000000*sin(ang1)));
    }

  AddSphere3D(Cn3D_pMain, NULL, (BigScalar)0, 0 , 2, 45000000, 0, 0, 5000000);
  AddCylinder3D(Cn3D_pMain, NULL, (BigScalar)0, 0 , 2,
                52000000,  5000000, 0,
                52000000, -5000000, 0, 1000000);
  AddCylinder3D(Cn3D_pMain, NULL, (BigScalar)0, 0, 2,
                58000000, 5000000,0,
                58000000,-5000000,0, 1000000);
  AddSphere3D(Cn3D_pMain, NULL, (BigScalar)0, 0, 2, 65000000, 0, 0, 5000000);
  AddCylinder3D(Cn3D_pMain, NULL, (BigScalar)0, 0, 2,
                72000000,  5000000, 0,
                72000000, -5000000, 0, 1000000);
  AddSphere3D(Cn3D_pMain, NULL, (BigScalar)0, 0, 2, 80000000, 0, 0, 5000000);
  AddCylinder3D (Cn3D_pMain, NULL, (BigScalar)0, 0 , 2, 87000000,5000000,0,
                                         87000000, -5000000, 0, 1000000);

  AddCylinder3D(Cn3D_pMain, NULL, (BigScalar)0, 0 , 2,
                93000000, 5000000, 0,
                93000000,-5000000,0, 1000000);


  if ( readErrors() )
    return;
  SetHLColor3D(Cn3D_v3d, 0, 0, 0);

  AttachPicture3D(Cn3D_v3d, Cn3D_pMain, NULL);

  readErrors();
  SetLayerTop3D( 0 );
}



NLM_EXTERN void LIBCALL Cn3D_ResetActiveStrucProc(void)
{
  Cn3D_CountDomainProc();

  ResetRenderCtrls();
  ResetLabelCtrls();
  ResetAlignCtrls();
  ResetDisplayCtrls();
  Cn3D_Redraw( TRUE ); /* always a new structure */
  Cn3dObjMgrGetSelected();
}


void LIBCALL Cn3D_SaveActiveCam(void)
{
  PDNMS pdnmsThis = GetSelectedModelstruc();
  PARS  pars      = GetAlgorRenderSet( pdnmsThis );
    
  if ( pars )
    GetViewerInfo3D(Cn3D_v3d, NULL, (Nlm_Camera3DPtr)pars, NULL);
}


/* Reset global palette index
 */
static void Cn3d_ResetPalette(void)
{
  Int2 i;
  for (i = 0;  i < CN3D_COLOR_MAX;  i++)
    Cn3d_IndexRGB[i] = 0;
}


/* Merge the active structure palette with the global palette
 */
static void Cn3d_MergePalette(void)
{
  PARS pars;
  Int2 i;

  /* fetch the active structure */
  PDNMS pdnmsThis = GetSelectedModelstruc();
  if ( !pdnmsThis )
    return;

  pars = GetAlgorRenderSet( pdnmsThis );
  if ( !pars )
    return;

  for (i = 0;  i < CN3D_COLOR_MAX;  i++)
    {
      if ( pars->IndexRGB[i] )
        Cn3d_IndexRGB[i] = 1;
    }
}


static void Cn3d_Lock3DPalette(Picture3D ppic)
{
  /* last step before traversing to draw */
  Int2 i, j;
  Int2 iColorCount = 0;

  for (i = 0;  i < CN3D_COLOR_MAX;  i++)
    {
      if ( Cn3d_IndexRGB[i] )
        iColorCount++;
    }

  if (iColorCount == 0)
    {  /* allocate the default color */
      AllocPalette3D(ppic, 1);
      SetColor3D(ppic, 0, 255, 255, 255);
      return;
    }

  AllocPalette3D(ppic, (Uint1)iColorCount);
  if ( readErrors() )
    return ;

  for (i = 0, j = 0;  j < CN3D_COLOR_MAX  &&  i < iColorCount;  j++)
    if ( Cn3d_IndexRGB[j] )
      {
        SetColor3D(ppic,  (Uint1)i,
                   Cn3d_PaletteRGB[j].red,
                   Cn3d_PaletteRGB[j].green,
                   Cn3d_PaletteRGB[j].blue);
        Cn3d_IndexRGB[j] = (Uint1)i;
        i++;

        if ( readErrors() )
          return;
      }
  ASSERT ( i == iColorCount );
}


NLM_EXTERN void LIBCALL Cn3D_Redraw(Boolean  New)
{
  /* fetch the active structure */
  PARS pars = NULL;
  PDNMS pdnmsThis = GetSelectedModelstruc();
  if (pdnmsThis == NULL  ||  (pars = GetAlgorRenderSet(pdnmsThis)) == NULL)
    {
      LogoProc( NULL );
      return;
    }

  if ( !New )
    GetViewerInfo3D(Cn3D_v3d, NULL,  (Camera3DPtr)pars, NULL);

  WatchCursor();
  ProgMon( "Removing 3D image ..." );
  if ( IsPlaying3D() )
    StopPlaying3D();
  Cn3D_DisableFileOps();

  if ( Cn3D_pMain )
    DeletePicture3D( Cn3D_pMain );
  Cn3D_pMain = CreatePicture3D();
  if ( readErrors() )
    return;

  ResetPicture3D( Cn3D_pMain );
  if ( readErrors() )
    return;

  Cn3d_ResetPalette();  /* clear global palette */
  MakeStrucPalette();  /* make palette  for the active structure */
  Cn3d_MergePalette(); /* merge it to global palette */
  Cn3d_Lock3DPalette( Cn3D_pMain ); /* Allocates the Picture3D palette  */

  ProgMon( "Rendering Structure..." );
  Cn3D_pMain = AlgorithmicRendering( Cn3D_pMain );
  if ( readErrors() )
    return;
  if (Cn3D_pMain == NULL)
    { /* must have something; do xyz origin... reset palette with white */
      Cn3D_pMain = CreatePicture3D ();
      if ( readErrors() )
        return;
      ResetPicture3D( Cn3D_pMain );
      if ( readErrors() )
        return;
      Cn3d_ResetPalette();  /* clear global palette */
      Cn3d_IndexRGB[C_white] = 1;  /* use white */
      Cn3d_Lock3DPalette( Cn3D_pMain ); /* Allocates the Picture3D palette  */
      Cn3D_pMain = Do3DOrigin( Cn3D_pMain ); /* make the origin */
    }
  SetHLColor3D(Cn3D_v3d, 255, 255, 0);   /* use yellow for highlight */

  ProgMon( "Redrawing 3D image ..." );
  AttachPicture3D(Cn3D_v3d, Cn3D_pMain, (Camera3DPtr)pars);
  readErrors();
  Cn3D_EnableFileOps();
  ArrowCursor();
}


/* Mouse Action Callbacks
 */

static void LIBCALLBACK DoHighlightSeg(PFB pfbThis,
                                       Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PMAD     pmadAtom = (PMAD)pfbThis;
  Viewer3D vvv      = (Viewer3D)ptr;
  PALD paldLoc;

  if ( !IsAtomNode(pfbThis) )
    return;

  paldLoc = GetAtomLocs(pmadAtom, (Int2)iModel);
  while ( paldLoc )
    {
      Segment3D seg = (Segment3D)paldLoc->pGraphic;
      if (seg != NULL)
        HighlightSeg3D(vvv, seg,
                       (Boolean)(!IsSeg3DHlighted(vvv, seg)));

      paldLoc = paldLoc->next; /* get next location */
    }
}


static Boolean GetPointData(MAPtr ma, PoinT point,
                            Prim3D    PNTR prim,
                            Segment3D PNTR seg,
                            BigScalar PNTR data)
{
  Prim3D    x_prim = NULL;
  Segment3D x_seg  = NULL;
  BigScalar x_data = 0;

  do /* TRY */
    {{
      Viewer3D  vvv = Nlm_MAToViewer3D( ma );
      Picture3D ppp = NULL;

      Nlm_GetViewerInfo3D(vvv, &ppp, NULL, NULL);
      if ( !ppp )
        break;

      {{
        Uint2 n_prim = FindPrim3D(vvv, point);
        if ( !n_prim )
          break;
      }}
      
      x_prim = GetFoundPrim3D(vvv, 0);
      if ( !x_prim )
        break;

      if (data || seg)
        GetPrimInfo3D(ppp, x_prim, &x_data, NULL, NULL, &x_seg, NULL);
    }}  while ( 0 );

  if ( prim )
    *prim = x_prim;
  if ( seg )
    *seg = x_seg;
  if ( data )
    *data = x_data;

  return (Boolean)(x_data != 0);
}


static void HLatom_MA(MAPtr ma,
                      MA_TracePtr trace, PoinT point, VoidPtr extra)
{
  Viewer3D  vvv  = MAToViewer3D( ma );
  Prim3D    prim = NULL;
  Segment3D seg  = NULL;
  BigScalar data = 0;

  if (GetPointData(ma, point, &prim, &seg, &data)  &&
      IsAtomLocNode( (PFB)data ))
    HighlightSeg3D(vvv, seg, (Boolean)(!IsSeg3DHlighted(vvv, seg)));

  else if (data  &&  IsObjectNode( (PFB)data ))
    HighlightPrim3D(vvv, prim, (Boolean)(!IsPrim3DHlighted(vvv, prim)));

  else
    {
      if ( !prim )
        BgColorDlg3D( vvv ); 
      return;
    }

  RedrawViewer3D( vvv );
}
/*--------------------- yanli --------------------*/
void fnCHLresidueRedraw(void)
{
   RedrawViewer3D(Cn3D_v3d);
}
/*--------------------- yanli --------------------*/
void DoCHighlightSeg(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr, Boolean highlight)
{  PMAD     pmadAtom = (PMAD)pfbThis;
  Viewer3D vvv      = (Viewer3D)ptr;
  PALD paldLoc;

  if ( !IsAtomNode(pfbThis) )
    return;

  paldLoc = GetAtomLocs(pmadAtom, (Int2)iModel);
  while ( paldLoc )
    {
      Segment3D seg = (Segment3D)paldLoc->pGraphic;
      if (seg != NULL){
          if(highlight && !IsSeg3DHlighted(vvv, seg)) {
          HighlightSeg3D(vvv, seg, (Boolean)(!IsSeg3DHlighted(vvv, seg)));
          }
          else if(!highlight && IsSeg3DHlighted(vvv, seg)){
          HighlightSeg3D(vvv, seg, (Boolean)(!IsSeg3DHlighted(vvv, seg)));
          }
      }

      paldLoc = paldLoc->next; /* get next location */
    }
}
/*--------------------- yanli --------------------*/
void fnCHLresidue(PDNMG pdnmgThis, Viewer3D  vvv, Boolean highlight)
{
                 /* highlight residues corresponding to those in Salsa Window */

  PMGD pmgdThis = NULL;
  PVNMA pvnmaThis = NULL;
  PMAD pmadThis = NULL;
  PVNAL pvnalThis = NULL;
  PALD paldThis = NULL;
  
  pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
  pvnmaThis = pmgdThis->pvnmaAHead;

  while ( pvnmaThis )
  {
         pmadThis = (PMAD)pvnmaThis->data.ptrvalue;
         if(pmadThis == NULL) goto setout;
         pvnalThis = pmadThis->pvnalLocate;
         if(pvnalThis == NULL) goto setout;
         paldThis = pvnalThis->data.ptrvalue;
         if(paldThis == NULL) goto setout;
         DoCHighlightSeg((PFB)pmadThis, (Int4)paldThis->pvnalLink->choice, 0, vvv, highlight);
         setout:
         pvnmaThis = pvnmaThis->next;
  }
  return;
}
/*--------------------- yanli --------------------*/
void fnPreCHLresidue(PDNMG pdnmgThis, Boolean highlight)
{

   Segment3D seg;   /* dead code */

/* seg = (Segment3D)paldThis->pGraphic;  */   /* dead code */
      fnCHLresidue(pdnmgThis,  Cn3D_v3d, highlight);
}
/*-----------------------------------------------*/
static void fnHLresidue (PDNMG pdnmgThis, PALD paldThis, Viewer3D  vvv)
/* helper function for HLresidue_MA.  highlights a single residue. */
{
  
  PMGD pmgdThis = (PMGD) pdnmgThis->data.ptrvalue; 
  PVNMA pvnmaThis = pmgdThis->pvnmaAHead;
  while ( pvnmaThis )
  {
    PMAD pmadThis = (PMAD)pvnmaThis->data.ptrvalue;
    DoHighlightSeg((PFB)pmadThis,
      (Int4)paldThis->pvnalLink->choice, 0, vvv);
    pvnmaThis = pvnmaThis->next;
  }
  return;
}            
/*------------- yanli ----------------*/
static void HLresidue_MA(MAPtr ma,
                         MA_TracePtr trace, PoinT point, VoidPtr extra)
/* highlight a residue */
/* this is the old version of Chris-highlight one residue only each double click */
{
  Viewer3D  vvv  = MAToViewer3D( ma );
  Prim3D    prim = NULL;
  BigScalar data = 0;

  Boolean  highlight = FALSE; 

  if (GetPointData(ma, point, &prim, NULL, &data)  &&
      IsAtomLocNode( (PFB)data ) )
    {
      PALD   paldThis = (PALD)data;
      PMGD   pmgdThis = GetParentGraph( (PFB)data );
      PDNMG  pdnmgThis = pmgdThis->pdnmgLink;
      Segment3D seq = paldThis->pGraphic;
      if(!IsSeg3DHlighted(vvv, seq)) highlight = TRUE;
      else highlight = FALSE;

      if(Cn3D_ObjMgrOpen) MediaObjSelect(pdnmgThis, paldThis, highlight);
      else {
/*        fnPreCHLresidue(pdnmgThis, paldThis, highlight);   */
          fnPreCHLresidue(pdnmgThis, highlight); 
               /* go to fnPreCHLresidue directly */
      }

    }    

  else if (data  &&  IsObjectNode( (PFB)data ))
    HighlightPrim3D(vvv, prim, (Boolean)(!IsPrim3DHlighted(vvv, prim)));

  else
    return;

  RedrawViewer3D( vvv );
}
/*------------- yanli ----------------*/
static void HLresidue_MA_future(MAPtr ma,
                         MA_TracePtr trace, PoinT point, VoidPtr extra)
/* highlight a residue */  
/* this is the lewis's version for features in the future */
/* HLresidue_MA in this cn3dwin.c is that from Chris's old version */
  /*-- highlight one residue each time */
{
  Viewer3D  vvv  = MAToViewer3D( ma );
  Prim3D    prim = NULL;
  BigScalar data = 0;
  Byte bToggle;

  if (GetPointData(ma, point, &prim, NULL, &data)  &&
    IsAtomLocNode( (PFB)data ) )
  {
    PALD paldThis = (PALD)data;
    PMGD pmgdThis = GetParentGraph( (PFB)data );
    PDNMG pdnmgCounter;
    PMMD pmmdThis = GetParentMol( (PFB) data);
    
    if (pmmdThis->pdnmgStartSelect &&  pmmdThis->pdnmgEndSelect)
      /* if already selected, clear the selection */
    {
      pdnmgCounter = pmmdThis->pdnmgHead;
      bToggle = 0;
      while(pdnmgCounter) 
      {
        if (pmmdThis->pdnmgStartSelect == pdnmgCounter  ||
          pmmdThis->pdnmgEndSelect == pdnmgCounter ) bToggle++;
        if (bToggle) 
        {
          fnHLresidue(pdnmgCounter, paldThis, vvv);
          if (bToggle == 2) bToggle = 0;
        }
        pdnmgCounter = pdnmgCounter->next;
      }
      pmmdThis->pdnmgStartSelect = NULL;
      pmmdThis->pdnmgEndSelect = NULL;
    }
    
    if (pmmdThis->pdnmgStartSelect == NULL) 
      /* if no selection yet */
    {
      pmmdThis->pdnmgStartSelect = pmgdThis->pdnmgLink;
      fnHLresidue (pmmdThis->pdnmgStartSelect, paldThis, vvv);
    }
    else
      /* if start selection exists, but new end selection */
    {
      pmmdThis->pdnmgEndSelect = pmgdThis->pdnmgLink;
      pdnmgCounter = pmmdThis->pdnmgHead;
      bToggle = 0;
      while(pdnmgCounter) 
      {
        if (pmmdThis->pdnmgStartSelect == pdnmgCounter  ||
          pmmdThis->pdnmgEndSelect == pdnmgCounter ) bToggle++;
        if (bToggle  && pmmdThis->pdnmgStartSelect != pdnmgCounter) 
          fnHLresidue(pdnmgCounter, paldThis, vvv);
        if (bToggle == 2) bToggle = 0;
        pdnmgCounter = pdnmgCounter->next;
      }
      
    }
  }
 
  else if (data  &&  IsObjectNode( (PFB)data ))
    HighlightPrim3D(vvv, prim, (Boolean)(!IsPrim3DHlighted(vvv, prim)));
  
  else
    return;
  
  RedrawViewer3D( vvv );
}


static void HLmolecule_MA(MAPtr ma,
                          MA_TracePtr trace, PoinT point, VoidPtr extra)
{
  Viewer3D  vvv  = Nlm_MAToViewer3D( ma );
  Prim3D    prim = NULL;
  BigScalar data = 0;

  if (GetPointData(ma, point, &prim, NULL, &data)  &&
      IsAtomLocNode( (PFB)data ) )
    {
      PALD paldThis = (PALD)data;
      PMMD pmmdThis = GetParentMol( (PFB)data );
      if ( pmmdThis )
        TraverseAtoms(pmmdThis->pdnmgHead,
                      (Int4)paldThis->pvnalLink->choice,
                      0, vvv, DoHighlightSeg);
    }

  else if (data  &&  IsObjectNode( (PFB)data ))
    HighlightPrim3D(vvv, prim, (Boolean)(!IsPrim3DHlighted(vvv, prim)));

  else
    return;

  RedrawViewer3D( vvv );
}


static Boolean Cn3D_InitMA(MAPtr ma, VoidPtr data)
{
/*MActionPtr hl_atom     = MA_AddAction(ma, MK_Shift, MA_DClick,
                                        HLatom_MA, data, "Highlight-Atom"); */

  MActionPtr hl_residue  = MA_AddAction(ma, MK_Normal,  MA_DClick,
                                        HLresidue_MA, data, "Highlight-Res");

/*MActionPtr hl_molecule = MA_AddAction(ma, MK_Ctrl,   MA_DClick,
                                        HLmolecule_MA, data, "Highlight-Mol"); */

  return (Boolean)(MA_SetAction(hl_residue,  TRUE));
/*return (Boolean)(MA_SetAction(hl_atom,     TRUE)  &&
                   MA_SetAction(hl_residue,  TRUE)  &&
                   MA_SetAction(hl_molecule, TRUE)); */
}



static void Cn3D_RenderCtrlProc(IteM i);
static void Cn3D_LabelCtrlProc(IteM i);
static void Cn3D_AlignCtrlProc(IteM i);
static void Cn3D_DisplayCtrlProc(IteM i);
static void Cn3D_ViewerCtrlProc(IteM i);

static GrouP Viewer3DGroups(WindoW w, Uint2Ptr width, Uint2 height)
{
  Int2  groups = 0;
  PoinT pnt;
  RecT  Cn3D_rRC, Cn3D_rVC;

  GrouP g = HiddenGroup(w, 0, 0, NULL);

#ifdef WIN_MOTIF
  SetGroupMargins(g, 0, 0);
  SetGroupSpacing(g, 8, 1);
#else
  SetGroupMargins (g, 1, 1);
  SetGroupSpacing (g, 0, 0);
#endif

  Cn3D_gViewCtrl = NormalGroup(g, 0, 0, "Viewer", systemFont, NULL);

  Cn3D_left = CreateControls3D(Cn3D_gViewCtrl, FALSE, TRUE, NULL);

  GetPosition(Cn3D_gViewCtrl, &Cn3D_rVC);
  Cn3D_Vy = (Uint2)(Cn3D_rVC.bottom - Cn3D_rVC.top + 10);
  Break( g );

  Cn3D_gRendCtrl = RenderControls( g );
  GetPosition(Cn3D_gRendCtrl, &Cn3D_rRC);
  Cn3D_Rx = (Uint2)(Cn3D_rRC.right - Cn3D_rRC.left + 10);

  pnt.x = Cn3D_rRC.left;
  pnt.y = Cn3D_rRC.top;
  SetNextPosition(g, pnt);
  Cn3D_gLabelCtrl = LabelControls( g );

  GetPosition(Cn3D_gLabelCtrl, &Cn3D_rRC);
  pnt.x = Cn3D_rRC.left;
  pnt.y = Cn3D_rRC.top;
  SetNextPosition(g, pnt);
  Cn3D_gAlignCtrl = AlignControls( g );

  GetPosition(Cn3D_gLabelCtrl, &Cn3D_rRC);
  pnt.x = Cn3D_rRC.left;
  pnt.y = Cn3D_rRC.top;
  SetNextPosition(g, pnt);
  Cn3D_gDisplayCtrl = DisplayControls( g );

  pnt.x = 10;
  pnt.y = 15;
  SetNextPosition(g, pnt);
  Cn3D_gViewer = HiddenGroup(g, 0, 0, NULL); /* the viewer */

  Cn3D_v3d = CreateViewer3D(Cn3D_gViewer, width, height,
                            X_ROTATE_SBAR | Y_ROTATE_SBAR | Z_ROTATE_SBAR,
                            Cn3D_ma_group_menu, Cn3D_ma_action_menu,
                            Cn3D_InitMA, NULL);
  if (Cn3D_v3d == NULL) {
    Message ( MSG_OK, "Cn3D Viewer - Insufficient Memory For Structures" );
    return NULL;
  }

  if (GetStatus(Cn3D_iViewCtrl) == TRUE)
    LinkControls3D(Cn3D_left, Cn3D_v3d);

  return HiddenGroup(w, 1, 0, NULL);
}


void Cn3DResizeProc(WindoW w)
{
  RecT r;
  ObjectRect(Cn3D_w, &r);
  OffsetRect(&r, (Int2)(-r.left), (Int2)(-r.top));

  InsetRect(&r, 5, 5);

  if ((GetStatus(Cn3D_iRendCtrl ) == TRUE)  ||
      (GetStatus(Cn3D_iLabelCtrl) == TRUE)  ||
      (GetStatus(Cn3D_iAlignCtrl) == TRUE)  ||
      (GetStatus(Cn3D_iDisplayCtrl) == TRUE ))
    r.left += Cn3D_Rx;

#ifdef WIN_MAC
  {{
    extern Handle Nlm_GetWindowMenuBar(WindoW w);

    Handle hdl = Nlm_GetWindowMenuBar( Cn3D_w );
    if ( hdl )
      {
        RecT mbr;
        ObjectRect(hdl, &mbr);
        r.top += (mbr.bottom - mbr.top);
      }
  }}
#endif

  if (GetStatus(Cn3D_iViewCtrl) == TRUE)
    r.top  += Cn3D_Vy;

  SetPosition3D(Cn3D_v3d, &r);
}


static void Cn3D_GifSaveProc(IteM i)
{
  Char fname[PATH_MAX];
  Char defname[32];
  PDNMS pdnmsThis = GetSelectedModelstruc();

  fname[0] = '\0';
  defname[0] = '\0';

  if (pdnmsThis == NULL)
    StringNCpy_0(defname,
                 pdnmsThis ? GetStrucStrings(pdnmsThis, PDB_ACC) : "cn3d",
                 sizeof(defname) - 4);
  StringCat(defname, ".gif");

  if ( GetOutputFileName(fname, sizeof(fname), defname) )
    {
      SaveImageGIF(Nlm_GetViewerImage3D(Cn3D_v3d), fname);
    }
}


#ifdef WIN_MSWIN
static void Cn3D_ImageCopyProc(IteM I)
{
 CopyViewer3D( Cn3D_v3d );
}

/*
static void Cn3D_ImagePrintProc(IteM I)
{
 PrintViewer3D( Cn3D_v3d );
}
*/
#endif



#include <document.h>
static ParData Cn3D_ParFmt = { FALSE, FALSE, FALSE, FALSE, TRUE, 0, 0 };
static ColData Cn3D_ColFmt = { 0, 0, 40, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static void Cn3D_AboutQuit(WindoW w)
{
  Remove( w );
  return;
}


static void Cn3D_AboutSize(WindoW w)
{
  Int4 height, width;
  RecT r;

  DoC  d = (DoC)GetObjectExtra( w );
  if ( !d )
    return;

  WatchCursor();

  ObjectRect(w, &r);
  width  = r.right  - r.left;
  height = r.bottom - r.top;

  GetPosition(d, &r);
  r.right  = (Int2)(width  - r.left);
  r.bottom = (Int2)(height - r.top);
  SetPosition(d, &r);

  AdjustPrnt(d, &r, FALSE);
  ObjectRect(d, &r);
  InsetRect(&r, 4, 4);
  Cn3D_ColFmt.pixWidth = (Int2)(screenRect.right - screenRect.left);
  Cn3D_ColFmt.pixInset = 8;
  if (Visible(d)  &&  AllParentsVisible(d))
    UpdateDocument(d, 0, 0);

  ArrowCursor();
  Update();
}

Boolean Cn3D_fEntrezOn; /* global */

static void Cn3D_AlignEdit(IteM i)
/* launches the sequin editor, salsa */
{
  PDNMS pdnmsMaster = NULL;
  PDNMS pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;
  SeqAnnotPtr psaAlign = NULL;
  SeqAlignPtr salp = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  Int2 iCount = 0;

  Int2  handled;

  Uint2 entityID, itemID, itemtype;    
/*********************/
  SeqEditViewProcsPtr  svpp;

  if(!Cn3D_ObjMgrOpen) {
    return;
  }

  svpp = (SeqEditViewProcsPtr) GetAppProperty ("SeqEditDisplayForm");
  if (svpp != NULL) {
     svpp->minPixelWidth = 650;
     svpp->minPixelHeight = 120;
     svpp->viewer_mode = TRUE;
  }
/************************/
  pdnmsMaster = GetSelectedModelstruc();
  pmsdMaster = (PMSD) pdnmsMaster->data.ptrvalue;

  pdnmmHead = pmsdMaster->pdnmmHead;
  pmmdThis = pdnmmHead->data.ptrvalue;

  if(pmsdMaster->pdnmsSlaves != NULL) 
  {
     pmsdSlave = pmsdMaster -> pdnmsSlaves -> data.ptrvalue;
  }
  if(pmsdMaster->psaAlignment != NULL) 
     psaAlign = pmsdMaster->psaAlignment;
  else if( pmsdSlave!=NULL) {
     if (pmsdSlave->psaAlignment != NULL) 
        psaAlign = pmsdSlave->psaAlignment;
  }

  if (psaAlign == NULL)
  { 
     for(iCount = 0; iCount < Num_Bioseq; iCount++){
        handled = GatherProcLaunch (OMPROC_EDIT, FALSE, mediadata[iCount]->entityID, mediadata[iCount]->itemID, OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
     /* Cn3D_Redraw(FALSE); */  /* cause problem? */
     /* Cn3dObjMgrGetSelected();  */
     }
     Cn3D_ReColor = TRUE;
     Cn3D_Redraw(FALSE);     /* cause problem? */
     Cn3D_ReColor = FALSE;
     Cn3dObjMgrGetSelected();
     Cn3D_SalsaOpen = TRUE;
  }
 
  if (psaAlign!=NULL && psaAlign->data!=NULL) 
  {
    salp=(SeqAlignPtr)psaAlign->data;
    LaunchSalsa(salp);
  
    pdnmsSlave = pmsdMaster->pdnmsSlaves;
    iCount = 1;
    while(pdnmsSlave){
        pmsdSlave = (PMSD) pdnmsSlave->data.ptrvalue;
        pdnmmHead = pmsdSlave->pdnmmHead;
        while(pdnmmHead){
           pmmdThis = pdnmmHead->data.ptrvalue;
           if(SeqIdForSameBioseq(pmmdThis->pSeqId, mediadata[iCount]->sip)) {

              Salsa_BioseqUpdate = FALSE;

              if(!pmsdSlave->bVisible) {
                 entityID = BioseqFindEntity(pmmdThis->pSeqId, &itemID);
                 itemtype = OBJ_BIOSEQ; 
                 ObjMgrSendMsg(OM_MSG_HIDE, entityID, itemID, itemtype);
                 Salsa_BioseqUpdate = TRUE;
              }
           }
           pdnmmHead = pdnmmHead->next;
       }

       iCount++;
       pdnmsSlave = pdnmsSlave->next;
    }

    Cn3D_ReColor = TRUE;
    Cn3D_Redraw(FALSE);     /* cause problem?  */
    Cn3D_ReColor = FALSE;
    Cn3dObjMgrGetSelected();
  }    

  return;
}

void LaunchSequenceWindow(void)
{
  IteM i = NULL;
  Cn3D_AlignEdit(i);
}

static void Cn3D_AboutStruc(IteM i)
{
  PMSD   pmsdThis;
  DoC    d;
  Char   path[PATH_MAX];
  FILE   *fp;
  WindoW Cn3D_wAbout;

  PDNMS pdnmsThis = GetSelectedModelstruc();
  if ( !pdnmsThis )
    return;

  pmsdThis = (PMSD)pdnmsThis->data.ptrvalue;

  Cn3D_wAbout = DocumentWindow(-66, -50,  -10, -10,
                               GetStrucStrings(pdnmsThis, PDB_ACC),
                               Cn3D_AboutQuit, Cn3D_AboutSize);

  d = DocumentPanel(Cn3D_wAbout,
                    (Int2)(35 * stdCharWidth + 17),
                    (Int2)(14 * stdLineHeight));
  SetObjectExtra(Cn3D_wAbout, (Pointer)d, NULL);
  Cn3D_ColFmt.pixWidth = (Int2)(60 * stdCharWidth);
  Cn3D_ColFmt.font = programFont;

  path[0] = '\0';
  TmpNam( path );
  fp = FileOpen(path,"w");
  if (fp != NULL)
    {
      WriteStructSummary(pdnmsThis, fp);
      fprintf(fp, "\n\nPDB Remarks (non-REFERENCE)\n\n");
      WritePDBRemarks(pdnmsThis, fp);
      fprintf(fp, "\n\n\n\n\n\n\n\n\n\n");
      fflush( fp );
      FileClose( fp );
     }

  DisplayFancy(d, path, &Cn3D_ParFmt, &Cn3D_ColFmt, programFont, 4);
  FileRemove( path );
  Show( Cn3D_wAbout );
}


/* Create the generic CN3D menu system;
 * return handler to the "File" top menu
 */
static MenU Cn3D_Menus(WindoW w)
{
  MenU file_menu, menu;

  /* FILE top menu
   */

  file_menu = menu = PulldownMenu(w, "File/F");

  Cn3D_sOpen = Cn3D_OpenSub( menu );  /* Open submenu(see cn3dopen.c) */
  /* Import menu item would go here */
  Cn3D_sSave = Cn3D_SaveSub( menu );

  SeparatorItem( menu );
  Cn3D_sExport = Cn3D_ExportSub( menu );
  CommandItem(menu, "Save GIF/S", Cn3D_GifSaveProc);

  SeparatorItem( menu );
  Cn3D_iSelStruc =
    CommandItem(menu, "Active /A", Cn3D_SelectDlg);  /* see cn3dslct.c */
  Cn3D_iClearStruc =
    CommandItem(menu, "Clear /C", Cn3D_ClearSelProc); /* see cn3dwipe.c */



  /* EDIT top menu
   */

  menu = PulldownMenu(w, "Edit/E");
#ifdef WIN_MSWIN
  CommandItem(menu, "Copy/C", Cn3D_ImageCopyProc);
  SeparatorItem( menu );
#endif
  Cn3D_iRendCtrl = StatusItem(menu,
                              "Rendering Settings/R", Cn3D_RenderCtrlProc);
  Cn3D_iLabelCtrl = StatusItem(menu,
                               "Label Settings/L", Cn3D_LabelCtrlProc);

  Cn3D_iAlignCtrl = StatusItem(menu, "Neighbor Controls/N", Cn3D_AlignCtrlProc);

  Cn3D_iDisplayCtrl = StatusItem(menu, "Display Controls/N", Cn3D_DisplayCtrlProc);
                      /* Yanli */

  Cn3D_ma_group_menu  = SubMenu(menu, "Mouse Settings/M" );
  Cn3D_ma_action_menu = NULL /* SubMenu(menu, "Mouse3D Actions")*/;



  /* View top menu
   */

  menu = PulldownMenu(w, "View/V");

  CommandItem(menu, "Sequence Window/S", Cn3D_AlignEdit);
  Cn3D_iViewCtrl = StatusItem(menu,
                              "Animation Controls/A", Cn3D_ViewerCtrlProc);
  CommandItem(menu, "Structure Info/I", Cn3D_AboutStruc);

  SeparatorItem( menu );

  CommandItem(menu, "Reset Perspective/P",  Cn3D_AllCB );
  CommandItem(menu, "Reset Presentation/N",      Cn3D_RenDefault);



  /* RENDER top menu
   */

  menu = PulldownMenu(w, "Structure/S");

  CommandItem(menu, "Sec. Structure/S", Cn3D_RenStruc);
  CommandItem(menu, "WireFrame/W",      Cn3D_RenWire);
  CommandItem(menu, "Tubular/T",        Cn3D_RenTube);
  CommandItem(menu, "Hierarchy/H",      Cn3D_RenHier);
  CommandItem(menu, "Spacefill/P",      Cn3D_RenSpace);
  CommandItem(menu, "Ball and Stick/B", Cn3D_RenBS);
  CommandItem(menu, "Neighbor/N", Cn3D_RenAlign);

  /* COLOR top menu
   */

  menu = PulldownMenu(w, "Color/C");

  CommandItem(menu, "Sec. Structure/S", Cn3D_ColStru);
  CommandItem(menu, "Domain/D",         Cn3D_ColDomain);
  CommandItem(menu, "Cycle Molecule/C", Cn3D_ColCy);
  CommandItem(menu, "Residue/R",        Cn3D_ColRes);
  CommandItem(menu, "Hydrophobicity/H", Cn3D_ColHydro);
  CommandItem(menu, "CPK/K",            Cn3D_ColCPK);
  CommandItem(menu, "Temperature/T",    Cn3D_ColTemp);
  CommandItem(menu, "Neighbor/O",    Cn3D_ColCons);
  CommandItem(menu, "Structure/U",    Cn3D_ColAlign);

 /* OPTIONS top menu 
  */

  menu = PulldownMenu(w, "Option/O");
 
  CommandItem(menu, "Background Color/B", Cn3D_BgColor);
  SeparatorItem(menu);
  CommandItem(menu, "Highlight Color/H", Cn3D_HLColor);
  
  /* CONTROLS top menu
   */

   menu = PulldownMenu(w, "Help/H");

   CommandItem(menu,"Help/H", Cn3D_HelpProc);
   CommandItem(menu, "About Cn3D/B", Cn3D_AboutProc);

  SetStatus(Cn3D_iRendCtrl, FALSE);
  SetStatus(Cn3D_iViewCtrl, FALSE);
  SetStatus(Cn3D_iAlignCtrl, FALSE);


  return file_menu;
}


static void Cn3D_RenderCtrlProc(IteM i)
{
  /* use hide/show */
  if (i == NULL)
    {
      SetStatus(Cn3D_iRendCtrl, FALSE);
      Hide( Cn3D_gRendCtrl );
      Update();
      return;
    }

  if (GetStatus(Cn3D_iRendCtrl) == TRUE)
    {
     Hide( Cn3D_gLabelCtrl );
     Hide( Cn3D_gAlignCtrl );
     Hide( Cn3D_gDisplayCtrl );
     Show( Cn3D_gRendCtrl );
     Update();
     if ( GetStatus(Cn3D_iLabelCtrl) == TRUE || GetStatus(Cn3D_iAlignCtrl) == TRUE || GetStatus(Cn3D_iDisplayCtrl) == TRUE)
       {
         SetStatus(Cn3D_iLabelCtrl, FALSE);
         SetStatus(Cn3D_iAlignCtrl, FALSE);
         SetStatus(Cn3D_iDisplayCtrl, FALSE);
         return;
       }
    }
  else
    {
     Hide( Cn3D_gRendCtrl );
     Update();
    }

  Cn3DResizeProc( Cn3D_w );
}



static void Cn3D_LabelCtrlProc(IteM i)
{
  /* use hide/show */
  if (i == NULL)
    {
      SetStatus(Cn3D_iLabelCtrl, FALSE);
      Hide( Cn3D_gLabelCtrl );
      Update();
      return;
    }

  if (GetStatus(Cn3D_iLabelCtrl) == TRUE)
    {
     Hide( Cn3D_gRendCtrl );
     Hide( Cn3D_gAlignCtrl );
     Hide( Cn3D_gDisplayCtrl );
     Show( Cn3D_gLabelCtrl );
     Update();
     if (GetStatus(Cn3D_iRendCtrl) == TRUE || GetStatus(Cn3D_iAlignCtrl) == TRUE || GetStatus(Cn3D_iDisplayCtrl) == TRUE)
        {
          SetStatus(Cn3D_iRendCtrl, FALSE);
          SetStatus(Cn3D_iAlignCtrl, FALSE);
          SetStatus(Cn3D_iDisplayCtrl, FALSE);
          return;
        }
    }
  else
    {
     Hide( Cn3D_gLabelCtrl );
     Update();
    }

  Cn3DResizeProc( Cn3D_w );
}

static void Cn3D_AlignCtrlProc(IteM i)
{
  /* use hide/show */
  if (i == NULL)
    {
      SetStatus(Cn3D_iAlignCtrl, FALSE);
      Hide( Cn3D_gAlignCtrl );
      Update();
      return;
    }

  if (GetStatus(Cn3D_iAlignCtrl) == TRUE)
    {
     Hide( Cn3D_gRendCtrl );
     Hide(Cn3D_gLabelCtrl);
     Hide(Cn3D_gDisplayCtrl);
     Show( Cn3D_gAlignCtrl );
     Update();
     if (GetStatus(Cn3D_iRendCtrl) == TRUE || GetStatus(Cn3D_iLabelCtrl) == TRUE || GetStatus(Cn3D_iDisplayCtrl) == TRUE)
        {
          SetStatus(Cn3D_iRendCtrl, FALSE);
          SetStatus(Cn3D_iLabelCtrl, FALSE);
          SetStatus(Cn3D_iDisplayCtrl, FALSE);
          return;
        }
     
    }
  else
    {
     Hide( Cn3D_gAlignCtrl );
     Update();
    }

  Cn3DResizeProc( Cn3D_w );
}

static void Cn3D_DisplayCtrlProc(IteM i)
{
  /* use hide/show */
  if (i == NULL)
    {
      SetStatus(Cn3D_iDisplayCtrl, FALSE);
      Hide( Cn3D_gDisplayCtrl );
      Update();
      return;
    }

  if (GetStatus(Cn3D_iDisplayCtrl) == TRUE)
    {
     Hide( Cn3D_gRendCtrl );
     Hide( Cn3D_gLabelCtrl );
     Hide( Cn3D_gAlignCtrl );
     Show( Cn3D_gDisplayCtrl );
     Update();
     if (GetStatus(Cn3D_iRendCtrl) == TRUE || GetStatus(Cn3D_iLabelCtrl) == TRUE || GetStatus(Cn3D_iAlignCtrl) == TRUE)
        {
          SetStatus(Cn3D_iRendCtrl, FALSE);
          SetStatus(Cn3D_iLabelCtrl, FALSE);
          SetStatus(Cn3D_iAlignCtrl, FALSE);
          return;
        }

    }
  else
    {
     Hide( Cn3D_gDisplayCtrl );
     Update();
    }

  Cn3DResizeProc( Cn3D_w );
}


static void Cn3D_ViewerCtrlProc(IteM i)
{
  if (i == NULL)
    {
      SetStatus(Cn3D_iViewCtrl, FALSE);
      Hide( Cn3D_gViewCtrl );
      Update();
      return;
    }

  if (GetStatus(Cn3D_iViewCtrl) == TRUE)
    {
      Show( Cn3D_gViewCtrl );
      Update();
      LinkControls3D(Cn3D_left, Cn3D_v3d);
    }
  else
    {
     Hide( Cn3D_gViewCtrl );
     Update();
    }

  Cn3DResizeProc( Cn3D_w );
}


static void ControlsShowing(void)
{
  Cn3D_RenderCtrlProc( NULL );
  Cn3D_LabelCtrlProc ( NULL );
  Cn3D_AlignCtrlProc(NULL);
  Cn3D_DisplayCtrlProc(NULL);
  Cn3D_ViewerCtrlProc( NULL );
}



/*  Default Quits
 */
static void Cn3D_QuitProc(IteM i)
{
  QuitProgram();
}

static void Cn3D_Quit(WindoW w)
{
  QuitProgram();
}


/* Create a complete CN3D window and GUI environment
 */
extern WindoW LIBCALL Cn3DWin(WndActnProc on_close, MenU *file_menu)
{
  static Boolean Cn3D_Window_Alive = FALSE;

  fnMMDBCn3Dmode();  /* make mmdbapi run in cn3d mode */
  Cn3D_fAlignOn = TRUE;
  Cn3D_fUnalignOn = TRUE;

  UpdateColorTable(Cn3d_PaletteRGB, sizeof(Cn3d_PaletteRGB), "ncbi_rgb.txt");

  if ( Cn3D_Window_Alive )
    return (Handle)Cn3D_w;

  /* to ensure that Motif will allocate only the colors having
   * index less than 32;  all other colorcells will be used (and
   * are to be redefined) by 3D-viewer
   */
  RestrictMotifColorsTo( 32 );





  /* CN3d window and menus
   */

  {{
    MenU menu;

    /* CN3d window and menubar
     */
    Cn3D_w = DocumentWindow(-33, -10, -10, -10, "Cn3D 2.01",
                            (on_close ? on_close : Cn3D_Quit), NULL);

    /* CN3D general menu set
     */
    menu = Cn3D_Menus( Cn3D_w );

    if ( file_menu )
      *file_menu = menu;
    else
      { /* standard quit */
        SeparatorItem( menu );
        CommandItem(menu, "Quit/Q", Cn3D_QuitProc);
      }
  }}


  {{
    Uint2 Cn3D_uSize;
    Int2  Cn3D_size = (Int2)MIN(screenRect.right, screenRect.bottom);
    Cn3D_size -= 128;
    if (Cn3D_size < 200)
      Cn3D_size = 200;
    else if (Cn3D_size > 400)
      Cn3D_size = 400;

    Cn3D_uSize = (Uint2)Cn3D_size;
    Cn3D_gWinGP = Viewer3DGroups(Cn3D_w, &Cn3D_uSize, Cn3D_uSize);
  }}

  ProcessUpdatesFirst( FALSE );
  if (Cn3D_gWinGP == NULL)
    return (Handle)NULL;

  Cn3D_EnableFileOps();
  ControlsShowing();
  RealizeWindow( Cn3D_w );
  Cn3D_Redraw( TRUE );
  Cn3dObjMgrGetSelected();
  SetResize(Cn3D_w, Cn3DResizeProc);
  Cn3DResizeProc( Cn3D_w );
  Cn3D_Window_Alive = TRUE;

  return Cn3D_w;
}

