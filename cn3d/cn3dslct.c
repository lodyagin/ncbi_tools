/*   cn3dslct.c
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
* File Name:  cn3dslct.c
*
* Author:  Christopher Hogue
*
* Version Creation Date:   1/31/96
*
* File Description: Cn3d structure selection dialog 
*                   
* Modifications:  
* --------------------------------------------------------------------------
* $Log: cn3dslct.c,v $
* Revision 6.8  1999/01/20 18:21:20  ywang
* include salmedia.h due to the move around of MediaInfo from cn3dmsg.h to the new created salmedia.h
*
 * Revision 6.7  1998/06/30  19:30:39  ywang
 * fix bugs for launching salsa when active structure changed
 *
 * Revision 6.6  1998/06/19  15:26:42  kans
 * needed to include cn3dmsg.h
 *
* Revision 6.5  1998/06/18 16:59:45  ywang
* add Cn3dObjRegister() to Cn3D_SelectStrucProc
*
 * Revision 6.4  1998/04/20  18:36:12  lewisg
 * moved extern for Viewer3d to cn3dmain.h
 *
* Revision 6.3  1998/04/16 00:32:26  lewisg
* corrected neighbor mode bugs
*
* Revision 6.2  1998/04/04 00:55:48  lewisg
* fixed active and clear dialog boxes to work on the master, not the slave
*
* Revision 6.1  1998/03/06 01:21:00  lewisg
* merge
*
* Revision 6.0  1997/08/25 18:13:45  madden
* Revision changed to 6.0
*
* Revision 5.0  1996/05/28 14:05:44  ostell
* Set to revision 5.0
*
 * Revision 1.5  1996/05/09  15:42:23  hogue
 * Fixed SGI compiler warnings.
 *
 * Revision 1.4  1996/04/26  18:42:43  vakatov
 * CN3D sources ported to MS-Windows;
 * the portability errors and warnings fixed, etc.
 *
 * Revision 1.3  1996/03/30  23:41:23  hogue
 * Redraw now saves camera
 *
 * Revision 1.2  1996/03/29  20:02:04  hogue
 * Added calls to get camera and reset active structure.
 *
 * Revision 1.1  1996/02/01  18:47:38  kans
 * Initial revision
 *
*
* ==========================================================================
*/

#include <viewer3d.h>
#include <mmdbapi.h>
#include <cn3dslct.h>
#include <cn3dmain.h>
#include <cn3dmsg.h>
#include <salmedia.h>
#include <algorend.h>


static WindoW	Cn3D_wStruSelect;
static LisT	Cn3D_lSelect;
static Boolean  Cn3D_Select_InUse = FALSE;
 
static void  Cn3D_SelectStrucProc(ButtoN B)
{
    Int4 iCount = 1;
    PDNMS pdnmsThis;
     
     
/* save the current camera for the active model before picking a new one */

    Cn3D_SaveActiveCam();
 
    pdnmsThis = GetFirstModelstruc();
    while (pdnmsThis)
     {
	 if (iCount == GetValue(Cn3D_lSelect))
	   {
	     if(AreNeighborsOn()) SetMasterModelstruc(pdnmsThis);
       /*else*/SetSelectedModelstruc(pdnmsThis);
       /* change the viewing mode, a bit repetitious, but what the heck */
       if(((PMSD)(pdnmsThis->data.ptrvalue))->pdnmsSlaves != NULL)
       {
         SetMasterModelstruc(pdnmsThis);
         SetNeighborOn();
       }
	     Remove(Cn3D_wStruSelect); 
	     Cn3D_EnableFileOps();
	     Cn3D_Select_InUse = FALSE;

         Cn3dObjRegister();  /* yanli */ 
         LaunchSequenceWindow();

 	     Cn3D_ResetActiveStrucProc(); 

	     return;
	   }
	 iCount++;
	 pdnmsThis = GetNextModelstruc();
     }
  Remove(Cn3D_wStruSelect);
  Cn3D_EnableFileOps();
  Cn3D_Select_InUse = FALSE;

  Cn3dObjRegister();  /* yanli */ 
  LaunchSequenceWindow();

  Cn3D_ResetActiveStrucProc(); 


  return;  
}


static void Cn3D_CancelSelectProc(ButtoN B)
{
  Remove(Cn3D_wStruSelect);
  Cn3D_EnableFileOps();
  Cn3D_Select_InUse = FALSE;
  return;
}


static void  Cn3D_AboutStrucProc(ButtoN B)
{
    Int4 iCount = 1;
    PDNMS pdnmsThis;
    PMSD pmsdThis;
    CharPtr pcContents;
    CharPtr pcClass;
    CharPtr pcSource;
    Char pcMMDB[INTSTRLEN];
     
    pdnmsThis = GetFirstModelstruc();
    while (pdnmsThis)
     {   /* find the one that was selected */
	 if (iCount == GetValue(Cn3D_lSelect))
	   {
	      break;
	   }
	 iCount++;
	 pdnmsThis = GetNextModelstruc();
     }
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;
    pcContents = StringSave( GetStrucStrings(pdnmsThis, LONG_NAME));
    pcClass = StringSave( GetStrucStrings(pdnmsThis, PDB_CLASS));
    pcSource = StringSave( GetStrucStrings(pdnmsThis, PDB_SOURCE));
    sprintf(pcMMDB,"%ld", (long) pmsdThis->iMMDBid);

    
    /* use a fancy displayer from Vibrant to arrange this information */
    
    MsgAlert(KEY_OK,SEV_INFO, "Structure is...", 
		   pcContents);
   
    if (pcContents) MemFree(pcContents);
    if (pcClass) MemFree(pcClass);
    if (pcSource) MemFree(pcSource);
    
    return;
}


void Cn3D_SelectDlg(IteM i)
{ /* a generic dialog for selecting a structure from in-memory ones */
  
    ValNodePtr pvnStruStrings = NULL;
    ValNodePtr pvnTemp = NULL;
    PDNMS pdnmsThis = NULL;
    Int2 iCount = 0;
    Int2 iSelected = 1;
    GrouP g;
    ButtoN b;
    PARS pars = NULL;
  
 
    if (Cn3D_Select_InUse) return;

    Cn3D_Select_InUse = TRUE;

    pdnmsThis = GetFirstModelstruc();
    if (!pdnmsThis)   
      {
         MsgAlert(KEY_NONE,SEV_ERROR, "No Structures", 
		  "No Structures in Memory");
 	 Cn3D_Select_InUse = FALSE;
	 return;
      }
    while (pdnmsThis)
	{
	    iCount++;
	     ValNodeCopyStr(&pvnStruStrings, iCount, 
				GetStrucStrings(pdnmsThis, PDB_ACC));
	    if (pdnmsThis == GetSelectedModelstruc())
		iSelected = iCount;	    
	 /*   printf("%d\n", iCount);  */
 	    pdnmsThis = GetNextModelstruc();
	}
    
    /* now we have a linked-list of the structure names */
      
    Cn3D_wStruSelect = ModalWindow(-20, -13,  -10,  -10, NULL);
    
    /* set up a group encolosing structures - models selection lists and - "info strings" */
    
    g = NormalGroup(Cn3D_wStruSelect, 0, 4, "Loaded Structures:",  systemFont, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 10, 5);
    Cn3D_lSelect = SingleList(g,  10, 5,NULL);
    b = DefaultButton(g, "OK", Cn3D_SelectStrucProc);
    b = PushButton(g, "Cancel", Cn3D_CancelSelectProc);
    b = PushButton(g, "More...", Cn3D_AboutStrucProc);
    pvnTemp = pvnStruStrings;
    while(pvnTemp)
      {
 	  ListItem(Cn3D_lSelect, (CharPtr) pvnTemp->data.ptrvalue);  
	  pvnTemp = pvnTemp->next;
      }
    if (pvnStruStrings) ValNodeFreeData(pvnStruStrings);
    SetValue(Cn3D_lSelect, iSelected);
    Select(Cn3D_lSelect);
    Cn3D_DisableFileOps();

    Show(Cn3D_wStruSelect);
    return;
} 

