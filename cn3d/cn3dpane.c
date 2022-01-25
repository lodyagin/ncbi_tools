/*   cn3dpane.c
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
Description: side panel for aligned structures.  lyg.

* $Log: cn3dpane.c,v $
* Revision 6.28  1999/01/20 18:21:20  ywang
* include salmedia.h due to the move around of MediaInfo from cn3dmsg.h to the new created salmedia.h
*
 * Revision 6.27  1998/10/07  23:10:45  ywang
 * merge align control with general display control
 *
 * Revision 6.26  1998/09/23  22:10:39  ywang
 * to record checkin logs
 *
*/

#include <viewer3d.h>
#include <cn3dmain.h>
#include <math.h>
#include <mmdbapi.h>
#include <cn3dpane.h>
#include <algorend.h>
#include <cn3dmsg.h>
#include <salmedia.h>
#include <salutil.h>
#include <cn3dopen.h>

static ButtoN Cn3D_bShowAlign; 
static ButtoN   Cn3D_bLRedraw;
static LisT Cn3D_lSlaves;

static ButtoN Cn3D_bAlignOn_bak, Cn3D_bUnalignOn_bak;

Boolean Cn3D_fAlignOn, Cn3D_fUnalignOn; /* globals  for above buttons */
 
static void fnNeighborMode(ButtoN b)
{

if(GetStatus(Cn3D_bShowAlign)) SetNeighborOn();
else SetNeighborOff();
SetStatus(Cn3D_bShowAlign, AreNeighborsOn());
Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */
Cn3D_Redraw(FALSE);
Cn3dObjMgrGetSelected();
}

static void fnAlignOn(ButtoN b)
{
Cn3D_fAlignOn = GetStatus(Cn3D_bAlignOn_bak);
Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */
Cn3D_Redraw(FALSE);
Cn3dObjMgrGetSelected();
}

static void fnUnalignOn(ButtoN b)
{
Cn3D_fUnalignOn = GetStatus(Cn3D_bUnalignOn_bak);
Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */
Cn3D_Redraw(FALSE);
Cn3dObjMgrGetSelected();
}


GrouP LIBCALL AlignControls_bak ( Nlm_GrouP prnt)
{
  GrouP g;

  g = HiddenGroup ( prnt, 2, 0, NULL );
  if (!g) return NULL;
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif


 
/*Cn3D_bShowAlign = CheckBox(g, "Alignment mode", fnNeighborMode );    

  Cn3D_bAlignOn = CheckBox(g, "Aligned", fnAlignOn);
  Break(g);
  Cn3D_bUnalignOn_bak = CheckBox(g, "Unaligned", fnUnalignOn);
*/

  Cn3D_bAlignOn_bak = CheckBox(g, "Aligned", NULL);
  Break(g);
  Cn3D_bUnalignOn_bak = CheckBox(g, "Unaligned", NULL);
 
/* set the values for the above */
/*ResetAlignCtrls();  */

  return g;
}


void LIBCALL ResetAlignCtrls(void)
{

/*SetStatus(Cn3D_bShowAlign, AreNeighborsOn()); */
  SetStatus(Cn3D_bAlignOn_bak, Cn3D_fAlignOn);
  SetStatus(Cn3D_bUnalignOn_bak, Cn3D_fUnalignOn);
  
}

