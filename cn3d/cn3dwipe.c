/*   cn3dwipe.c
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
* File Name:  cn3dwipe.c
*
* Author:  Christopher Hogue
*
* Version Creation Date:   1/31/96
*
* File Description: Cn3d clear structures dialog 
*                   
* Modifications:  
* --------------------------------------------------------------------------
* $Log: cn3dwipe.c,v $
* Revision 6.4  1998/05/06 14:20:12  lewisg
* get rid of NULL structure bugs
*
* Revision 6.3  1998/04/16 00:32:27  lewisg
* corrected neighbor mode bugs
*
* Revision 6.2  1998/04/04 00:55:47  lewisg
* fixed active and clear dialog boxes to work on the master, not the slave
*
* Revision 6.1  1998/03/06 01:23:23  lewisg
* merge
*
* Revision 6.0  1997/08/25 18:13:50  madden
* Revision changed to 6.0
*
* Revision 5.1  1996/07/22 00:27:42  hogue
* Now upon clear it finds last structure loaded & refreshes 3D view.
*
 * Revision 5.0  1996/05/28  14:05:44  ostell
 * Set to revision 5.0
 *
 * Revision 1.3  1996/05/09  15:42:42  hogue
 * Fixed SGI compiler warnings
 *
 * Revision 1.2  1996/04/26  18:42:56  vakatov
 * CN3D sources ported to MS-Windows;
 * the portability errors and warnings fixed, etc.
 *
 * Revision 1.1  1996/02/01  18:47:38  kans
 * Initial revision
 *
*
* ==========================================================================
*/

#include <vibrant.h>
#include <mmdbapi.h>
#include <cn3dwipe.h>
#include <cn3dmain.h>

 
static Boolean  Cn3D_Clear_InUse = FALSE;
static WindoW	Cn3D_wStruClear;
static LisT	Cn3D_lClearStruc;




static void Cn3D_ClearAllProc(ButtoN b)
{
  ClearStructures();
  Remove(Cn3D_wStruClear);
  Cn3D_EnableFileOps();
  Cn3D_Clear_InUse = FALSE;
  Cn3D_ResetActiveStrucProc();       
  return;
}


static void Cn3D_CancelClearProc(ButtoN b)
{
  Remove(Cn3D_wStruClear);
  Cn3D_EnableFileOps();
  Cn3D_Clear_InUse = FALSE;
  return;
}

static void Cn3D_ClearStrucProc(ButtoN b)
{   
    Int2 iCount = 1;
    PDNMS pdnmsThis = NULL;
    PDNMS pdnmsNext = NULL;
     
     
    pdnmsThis = GetFirstModelstruc();
    while (pdnmsThis)
     {
	 if (GetItemStatus(Cn3D_lClearStruc, iCount) == TRUE)
	   {
	     FreeAModelstruc(pdnmsThis);
	   }
	 iCount++;
	 ASSERT ( iCount > 0 );
	 pdnmsThis = GetNextModelstruc();
     }
 
  pdnmsThis = GetSelectedModelstruc();
  if (pdnmsThis == NULL) /* selected structure was cleared */
     { /* go to the last Modelstruc on the list */
        pdnmsThis = GetFirstModelstruc();
        if (pdnmsThis != NULL)
          {
             pdnmsNext = GetNextModelstruc();
             while (pdnmsNext)
               {
                  pdnmsThis = pdnmsNext;
                  pdnmsNext = GetNextModelstruc();
               }
          }
        if(AreNeighborsOn()) SetMasterModelstruc(pdnmsThis);
        /*else*/ SetSelectedModelstruc(pdnmsThis);
       /* change the viewing mode, a bit repetitious, but what the heck */
        if(pdnmsThis)
        {
          if(((PMSD)(pdnmsThis->data.ptrvalue))->pdnmsSlaves != NULL)
          {
            SetMasterModelstruc(pdnmsThis);
            SetNeighborOn();
          }
        }

        Remove(Cn3D_wStruClear);   
        Cn3D_EnableFileOps();
        Cn3D_Clear_InUse = FALSE; 
        Cn3D_ResetActiveStrucProc();       
        return;
     }
  /* otherwise same selected structure remains */
  Remove(Cn3D_wStruClear);   
  Cn3D_EnableFileOps();
  Cn3D_Clear_InUse = FALSE; 
  return;
}




void Cn3D_ClearSelProc(IteM i)
{
   
    ValNodePtr pvnStruStrings = NULL;
    ValNodePtr pvnTemp = NULL;
    PDNMS pdnmsThis = NULL;
    Int2 iCount = 0;
    Int2 iSelected = 0;
    GrouP g;
    ButtoN b;
    
    if (Cn3D_Clear_InUse) return;
    else Cn3D_Clear_InUse = TRUE;
   
    pdnmsThis = GetFirstModelstruc();
    if (!pdnmsThis)   
      {
         MsgAlert(KEY_NONE,SEV_ERROR, "No Structures", 
		  "No Structures in Memory");
	 Cn3D_Clear_InUse = FALSE;
	 return;
      }
    while (pdnmsThis)
	{
	    iCount++;
	     ValNodeCopyStr(&pvnStruStrings, iCount, 
				GetStrucStrings(pdnmsThis, PDB_ACC));
	    if (pdnmsThis == GetSelectedModelstruc())
		iSelected = iCount;	    
 	    pdnmsThis = GetNextModelstruc();
	}
    
    /* now we have a linked-list of the structure names */
      
    Cn3D_wStruClear = ModalWindow(-20, -13,  -10,  -10, NULL);
    
    /* set up a group encolosing structures - models selection lists and - "info strings" */
    
    g = NormalGroup(Cn3D_wStruClear, 0, 4, "Clear Structures:",  systemFont, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 10, 5);
    Cn3D_lClearStruc = MultiList(g,  10, 5,NULL);
     b = DefaultButton(g, "OK", Cn3D_ClearStrucProc);
     b = PushButton(g, "Clear ALL", Cn3D_ClearAllProc);
     b = PushButton(g, "Cancel", Cn3D_CancelClearProc);
     
    pvnTemp = pvnStruStrings;
    while(pvnTemp)
      {
 	  ListItem(Cn3D_lClearStruc, (CharPtr) pvnTemp->data.ptrvalue);  
	  pvnTemp = pvnTemp->next;
      }
    if (pvnStruStrings) ValNodeFreeData(pvnStruStrings);
    if (iSelected) SetItemStatus(Cn3D_lClearStruc, iSelected, TRUE);
    Select(Cn3D_lClearStruc);
    
    /* disable appropriate stuff here */
    Cn3D_DisableFileOps();
    Show(Cn3D_wStruClear);
    return;  
    
}
