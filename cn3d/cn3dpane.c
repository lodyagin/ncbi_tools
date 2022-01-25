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
* Revision 6.26  1998/09/23 22:10:39  ywang
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
#include <salutil.h>
#include <cn3dopen.h>

static ButtoN Cn3D_bShowAlign; 
static ButtoN   Cn3D_bLRedraw;
static LisT Cn3D_lSlaves;
static ButtoN  Cn3D_bAlignOn, Cn3D_bUnalignOn; 

Boolean Cn3D_fAlignOn, Cn3D_fUnalignOn; /* globals  for above buttons */
 

static void fnAlignList(LisT l)
/* set the values of the alignment pane */
{
  Int4 iCount;
  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;
  BiostrucFeaturePtr pbsfThis;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  Byte bVisible;
  SeqIdPtr sipThis;
  Uint2 entityID, itemID, itemtype;         /* yanli */

  Cn3D_SaveActiveCam();

  pdnmsMaster = GetSelectedModelstruc();
  if (pdnmsMaster != NULL)
  {

    pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
    pmsdMaster->bAligned = 0;
    pdnmsSlave = pmsdMaster->pdnmsSlaves;
    TraverseGraphs( pdnmsMaster, 0, 0, NULL, fnClearMarkedResidues);
    pbsfThis =  pmsdMaster->pdnsfFeatures->data.ptrvalue;
    iCount = 1;
    while(pdnmsSlave)
    {
      pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;

      bVisible = pmsdSlave->bVisible;   /* yanli */

      TraverseGraphs( pdnmsSlave, 0, 0, NULL, fnClearMarkedResidues);

      pmsdSlave->bVisible = GetItemStatus(Cn3D_lSlaves, iCount);

      if (pmsdSlave->bVisible) {
	fnMarkAlignedResidues(pdnmsMaster, pdnmsSlave, pbsfThis);
	pmsdMaster->bAligned++;
      }
      pbsfThis = pbsfThis->next;

      if(bVisible != pmsdSlave->bVisible){
         pdnmmHead = pmsdSlave->pdnmmHead;
         while(pdnmmHead){
            pmmdThis = pdnmmHead->data.ptrvalue;
            if(pmmdThis){
               if(Num_Bioseq == 0) return;  /* important */

               sipThis = pmmdThis->pSeqId;

               if(SeqIdForSameBioseq(sipThis, mediadata[iCount]->sip)) {
                  sipThis = mediadata[iCount]->sip;
                  entityID = BioseqFindEntity(sipThis, &itemID);
                  itemtype = OBJ_BIOSEQ;

                  Salsa_BioseqUpdate = FALSE;

                  if(pmsdSlave->bVisible) {
/*                   ObjMgrAlsoSelect(entityID, itemID, itemtype, 0, NULL); */
                     ObjMgrSendMsg(OM_MSG_SHOW, entityID, itemID, itemtype);
                     Salsa_BioseqUpdate = TRUE;
                  }
                  else if(!pmsdSlave->bVisible) {

                     ObjMgrSendMsg(OM_MSG_HIDE, entityID, itemID, itemtype);
                     Salsa_BioseqUpdate = TRUE;

/*                   if(!is_selectedbyID (entityID, itemID, OBJ_BIOSEQ)){ 
                        ObjMgrAlsoSelect(entityID, itemID, itemtype, 0, NULL);
                        ObjMgrDeSelect(entityID, itemID, itemtype, 0, NULL);   
                        ObjMgrSendMsg(OM_MSG_HIDE, entityID, itemID, itemtype);
                        Salsa_BioseqUpdate = TRUE;
                     }
                     else {
                        ObjMgrDeSelect(entityID, itemID, itemtype, 0, NULL);    
                        ObjMgrSendMsg(OM_MSG_HIDE, entityID, itemID, itemtype);
                        Salsa_BioseqUpdate = TRUE;
                     }   */
                  }
                  break;
               }
            }

            pdnmmHead = pdnmmHead->next;

         }
      }                        /* yanli */

      iCount++;
      pdnmsSlave = pdnmsSlave->next;
    }
  }
  Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */

  Salsa_BioseqUpdate = FALSE;    /* Yanli Jul. 17 */
          /* this will be cleaned up when salsa decides whether to keep color */

  Cn3D_ReColor = TRUE;   /* for Colombe's notice */
  Cn3D_Redraw(FALSE);
  Cn3D_ReColor = FALSE;  
  Cn3dObjMgrGetSelected();
}



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
Cn3D_fAlignOn = GetStatus(Cn3D_bAlignOn);
Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */
Cn3D_Redraw(FALSE);
Cn3dObjMgrGetSelected();
}

static void fnUnalignOn(ButtoN b)
{
Cn3D_fUnalignOn = GetStatus(Cn3D_bUnalignOn);
Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */
Cn3D_Redraw(FALSE);
Cn3dObjMgrGetSelected();
}


GrouP LIBCALL AlignControls ( Nlm_GrouP prnt)
{
  GrouP g;

  g = NormalGroup ( prnt, 1, 0, "Aligned Structures", systemFont, NULL );
  if (!g) return NULL;
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif


 
  Cn3D_bShowAlign = CheckBox(g, "Alignment mode", fnNeighborMode );

  Cn3D_bAlignOn = CheckBox(g, "Show Aligned", fnAlignOn);

  Cn3D_bUnalignOn = CheckBox(g, "Show Unaligned", fnUnalignOn);
 

  Cn3D_lSlaves = MultiList(g, 10, 15,  fnAlignList);
  

/* set the values for the above */
  ResetAlignCtrls();

  return g;
}


void LIBCALL ResetAlignCtrls(void)
/* set the values of the alignment pane */
{
  Int4 iCount;
  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;

  SetStatus(Cn3D_bShowAlign, AreNeighborsOn());
  SetStatus(Cn3D_bAlignOn, Cn3D_fAlignOn);
  SetStatus(Cn3D_bUnalignOn, Cn3D_fUnalignOn);
  
  Reset(Cn3D_lSlaves);
  pdnmsMaster = GetSelectedModelstruc();
  if (pdnmsMaster != NULL)
  {
    pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
    pdnmsSlave = pmsdMaster->pdnmsSlaves;
    
    iCount = 1;
    while(pdnmsSlave)
    {
      pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;
      ListItem(Cn3D_lSlaves, (CharPtr) pmsdSlave->pcPDBName);  
      if(pmsdSlave->bVisible) SetItemStatus(Cn3D_lSlaves, iCount, TRUE);
      iCount++;
      pdnmsSlave = pdnmsSlave->next;
    }
  }
}

