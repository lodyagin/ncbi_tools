/*   cn3dmsg.c
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
* File Name:  cn3dmsg.c
*
* Author: Yanli Wang 
*
* Version Creation Date:   3/26/98
*
* File Description: Main functions for building up cn3d/salsa communication
*
* Modifications:
* $Log: cn3dmsg.c,v $
* Revision 6.32  1998/09/23 22:08:42  ywang
* to record checkin log
*
* ===========================================================================  */

#include <vibrant.h>
#include <document.h>
#include <vsm.h>
#include <sequtil.h>   /* for sequence load funcs */
#include <objsub.h>
#include <string.h>
#include <saledit.h>
#include <objmgr.h>
#include <cn3dmsg.h>
#include <mmdbapi.h>
#include <accentr.h>
#include <lsqfetch.h>
#include <salsap.h>


Boolean Cn3D_ObjMgrOpen;
Boolean Cn3D_SalsaOpen;
Boolean Salsa_BioseqUpdate;
Boolean Cn3D_ReColor;
Int4 Num_Bioseq, Num_Biostruc;
Int4 Num_ActiveSlave;

/* MediaInfo mediadata[MaxObj];  */
MediaInfo **mediadata;
Uint2 sap_entityID, sap_itemID;

extern Boolean Cn3D_fEntrezOn;

/*-----------------------------------------*/
void LaunchMediaViewer(BioseqPtr bsp)
{
  Uint2            entityID = 0;
  Int2             handled = 0;
  Uint2            options = 0;

  Uint2            itemID = 0;

  if(bsp != NULL) {
      entityID = BioseqFindEntity(bsp->id, &itemID);
      if(entityID == 0) entityID = ObjMgrRegister(OBJ_BIOSEQ, (Pointer) bsp);    
      options = ObjMgrGetOptions(entityID);
      options |= OM_OPT_FREE_IF_NO_VIEW;
      ObjMgrSetOptions(options, entityID);
      handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1, OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
      BioseqUnlock( bsp );
   }
}
/*-----------------------------------------*/
extern void MediaLaunch(void)
{
  BioseqPtr bsp = NULL;

  if(Num_Bioseq == 0) return;
  bsp = BioseqLockById(mediadata[0]->sip);
  if(bsp == NULL) ErrPostEx (SEV_ERROR, 0, 0, " BioseqLockById failed in MediaLaunch!\n");
  else {
     LaunchMediaViewer ((BioseqPtr) bsp);
  }

}
/*-----------------------------------------*/
static SeqIdPtr make_sip(Int4 id) {

  SeqIdPtr     sip=NULL;

  sip = ValNodeNew (NULL);
  if (sip != NULL) {
     sip->choice = SEQID_GI;
     sip->data.intvalue = id;
  }
  return sip;
}
/*-----------------------------------------*/
extern void MediaDataLoad2(PDNMM pdnmmHead)
{
  PDNMM pdnmmHead_tmp = NULL;
  PMMD  pmmdThis = NULL;
  BioseqPtr bsp;
  SeqIdPtr sip;
  Int4 thisGi;
  Int2 iCount = 0;
  Uint2 entityID, itemID;

                                          /* allocate memory for mediadata */     
  pdnmmHead_tmp = pdnmmHead;
  while(pdnmmHead_tmp){
     pmmdThis = pdnmmHead_tmp->data.ptrvalue;
     if(pmmdThis == NULL) goto errot_tmp;
     if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto errot_tmp;
     iCount++;
     errot_tmp:
     pdnmmHead_tmp = pdnmmHead_tmp->next;
  }

  Num_Bioseq = iCount;

  if(Num_Bioseq == 0) return;

  mediadata = (Pointer) MemNew((size_t) (Num_Bioseq * sizeof(Pointer)));
  for(iCount = 0; iCount < Num_Bioseq; iCount++){
     mediadata[iCount] = (MediaInfoPtr) MemNew(sizeof (MediaInfo));
  }

  iCount = 0;
  while(pdnmmHead){
     pmmdThis = pdnmmHead->data.ptrvalue;
     if(pmmdThis == NULL) goto errot;
     if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto errot;
     else thisGi = pmmdThis->iGi;     
               /* pmmdThis->iGi is either 0 or what it should be from mmdbapi */
               /* in one biostruc-sequence case from MMDB, it should NOT be zero for NA and Protein */
/*   sip = make_sip (thisGi);  */
     sip = pmmdThis->pSeqId;
     mediadata[iCount]->sip = sip;
     mediadata[iCount]->Gi = thisGi;
     bsp = BioseqLockById (sip); /* Remember always lock first! */
                                 /* and by the way to get bsp->length */
     if (bsp!=NULL) {
        mediadata[iCount]->length = bsp->length;
        BioseqUnlock (bsp);
     }
     entityID = BioseqFindEntity(sip, &itemID);
     mediadata[iCount]->entityID = entityID;
     mediadata[iCount]->itemID = itemID;
     mediadata[iCount]->bVisible = 1;
     iCount++;
     errot:
     pdnmmHead = pdnmmHead->next;
  }

  Num_ActiveSlave = iCount - 1;
}
/*-----------------------------------------*/
extern void MediaDataLoad(SeqAlignPtr salp)
{
  SeqIdPtr      sip = NULL, sip2 = NULL;
  Int4 iCount = 0;
  Boolean FirstNode = TRUE;
  Uint2 entityID, itemID;
  BioseqPtr bsp;
  SeqAlignPtr salp_tmp;

                                    /* allocate memory for mediadata */
  salp_tmp = salp;
  while(salp_tmp){
     sip2 = SeqAlignId (salp_tmp, 0);
     if (sip2 != NULL) {
       if(!FirstNode) sip2 = sip2->next;
       for(sip = sip2; sip != NULL; sip = sip->next){
          iCount++;
          FirstNode = FALSE;
       }
     }
     salp_tmp = salp_tmp->next;
  }
  Num_Bioseq = iCount;

  if(Num_Bioseq == 0) return;  /* important */

  mediadata = (Pointer) MemNew((size_t)(Num_Bioseq * sizeof(Pointer)));
  for(iCount = 0; iCount < Num_Bioseq; iCount++){
     mediadata[iCount] = (MediaInfoPtr) MemNew(sizeof(MediaInfo));
  }

  iCount = 0;
  FirstNode = TRUE;

  while(salp){
     sip2 = SeqAlignId (salp, 0);
     if (sip2 != NULL) {
       if(!FirstNode) sip2 = sip2->next;
       for(sip = sip2; sip != NULL; sip = sip->next){
          mediadata[iCount]->sip = sip;
          mediadata[iCount]->Gi = GetGIForSeqId(sip);
                    /* GetGIForSeqId(sip) is 0 if no Gi found */
          mediadata[iCount]->bVisible = TRUE; 
          iCount++;
          FirstNode = FALSE;
       }
     }
     salp = salp->next;
  }

  for(iCount = 0; iCount < Num_Bioseq; iCount++){
     bsp = BioseqLockById(mediadata[iCount]->sip);  /* lock first, and by the way to get bsp->length */
     if(bsp != NULL) {
        mediadata[iCount]->length = bsp->length;
        mediadata[iCount]->sip = SeqIdDup(SeqIdFindBest(bsp->id,0));
              /* to use the complete set of SeqId */
        BioseqUnlock( bsp );
     }
     entityID = BioseqFindEntity(mediadata[iCount]->sip, &itemID);
     mediadata[iCount]->entityID = entityID;
     mediadata[iCount]->itemID = itemID;
     mediadata[iCount]->bVisible = 1;
  }  
  
}
/*-----------------------------------------*/
extern void MediaRegister(void)
{
  REGISTER_BIOSEQ_BIOSTRUC_MEDIA;

}
/*-----------------------------------------*/
extern void Cn3dObjRegister(void)
{
  PDNMS pdnmsMaster = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;
  SeqAnnotPtr psaAlign = NULL;
  SeqAlignPtr salp = NULL;
  PDNMM pdnmmHead = NULL, pdnmmHead_tmp = NULL;
  PMMD  pmmdThis = NULL;

  ObjMgrPtr omp;

               /* turn off/on dispatcher here */
  if(Cn3D_fEntrezOn){
     EntrezBioseqFetchDisable();        

            /* hard lesson here: once FetchBS is used, cn3d can not get bsp by BioseqLockById though "EntrezBioseqFetchEnable" is called before "FetchBS" , as a result I have to do "EntrezBioseqFetchDisable", then "EntrezBioseqFetchEnable" again to be able to get bsp */ 

   }
/* if(!EntrezBioseqFetchEnable("Cn3D", FALSE)) ErrPostEx (SEV_ERROR, 0, 0, " EntrezBioseqFetchEnable failed!\n");   */
   EntrezBioseqFetchEnable("Cn3D", FALSE);
  
/* if(!EntrezInit("Cn3D", FALSE, NULL)) ErrPostEx (SEV_ERROR, 0, 0, " EntrezInit failed!\n");     
   if(!BioseqFetchInit(TRUE)) ErrPostEx (SEV_ERROR, 0, 0, " BioseqFetchInit failed!\n");     */      /* not neccessary */

  if(!Cn3D_fEntrezOn){
                            /* to add the ability to open salsa */
     MediaRegister();                      
     SalsaRegister();

     Cn3D_fEntrezOn = TRUE;

     omp = ObjMgrGet();
     omp->maxtemp = 60;
     ObjMgrSetHold();     /* magic line? */
  }

/*pdnmsMaster = GetMasterModelstruc();    
  if(pdnmsMaster == NULL) pdnmsMaster = (PDNMS) GetSelectedModelstruc(); */

  Cn3D_ObjMgrOpen = TRUE;   

  pdnmsMaster = (PDNMS) GetSelectedModelstruc();
          /* always use this--lewis's opinion */

  if(pdnmsMaster == NULL) {
/*   ErrPostEx (SEV_ERROR, 0, 0, " GetSelectedModelstruc is NULL: return! ");  */
     return;
  }

  pmsdMaster = (PMSD) pdnmsMaster->data.ptrvalue;

  pdnmmHead = pmsdMaster->pdnmmHead;

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
     MediaDataLoad2 (pdnmmHead);
     MediaLaunch();
     return;
  }

  if (psaAlign!=NULL && psaAlign->data!=NULL)
  {
    salp=(SeqAlignPtr)psaAlign->data;

    MediaDataLoad (salp);
    MediaLaunch();

  }  

  return;
}
/*-----------------------------------------*/
static Int2 LIBCALLBACK SeqStrucMediaMsgFunc(OMMsgStructPtr ommsp)
{
  OMUserDataPtr      omudp = NULL;
  MediaViewPtr       mvp = NULL;
 
  SelStructPtr sel = NULL;

  Boolean highlight = FALSE;

  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  mvp = (MediaViewPtr) omudp->userdata.ptrvalue;
  if(mvp == NULL){
      return OM_MSG_RET_ERROR;
  }

  switch (ommsp->message)
  {
      case OM_MSG_UPDATE:
           break;

      case OM_MSG_SELECT:

           if(ommsp->itemtype == OBJ_SEQALIGN) return OM_MSG_RET_OK;

           sel = ObjMgrGetSelected();
           while(sel != NULL) {
               highlight = TRUE;
               MediaHL(sel, highlight);  /* highlight in cn3d */
               sel = sel->next;
           }
           break;

      case OM_MSG_DESELECT:

           if(ommsp->itemtype == OBJ_SEQALIGN) return OM_MSG_RET_OK;

           sel = (SelStructPtr) MemNew(sizeof(SelStruct));
           sel->entityID = ommsp->entityID;
           sel->itemtype = ommsp->itemtype;
           sel->itemID = ommsp->itemID;
           sel->region = ommsp->region;
           sel->regiontype = ommsp->regiontype;
           if(sel) {
               highlight = FALSE;
               MediaHL(sel, highlight);  /* highlight in cn3d */
           }
           break;

      case OM_MSG_SETCOLOR:
           break;

  }

  return OM_MSG_RET_OK;

}
/*-----------------------------------------*/
extern Int2 LIBCALLBACK SeqStrucMediaFunc(Pointer data)
{
  MediaViewPtr        mvp = NULL;
  OMProcControlPtr    ompcp = NULL;
  OMUserDataPtr       omudp = NULL;
  SelStruct           ss;
  BioseqPtr           bsp = NULL;

  Uint2 slave_entityID, slave_itemID;
  Int4 iCount;
 

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) {
      Message (MSG_ERROR, "Data NULL [1]");
      return OM_MSG_RET_ERROR;
  }

  switch (ompcp->input_itemtype) {
      case OBJ_BIOSEQ :
           bsp = (BioseqPtr) ompcp->input_data;
           break;
      case 0 :
           return OM_MSG_RET_ERROR;
      default :
           return OM_MSG_RET_ERROR;
  }

  if (bsp == NULL) {
      Message (MSG_ERROR, "Data NULL [2]");
      return OM_MSG_RET_ERROR;
  }

  mvp = (MediaViewPtr) MemNew(sizeof (MediaView));

  ss.entityID = ompcp->input_entityID;
  ss.itemID = ompcp->input_itemID;
  ss.itemtype = ompcp->input_itemtype;
  ss.regiontype =0;
  ss.region = NULL;

  if(mvp != NULL){
      mvp->input_entityID = ompcp->input_entityID;
      mvp->input_itemID = ompcp->input_itemID;
      mvp->input_itemtype = ompcp->input_itemtype;
      mvp->this_itemtype = OBJ_SEQANNOT;
      mvp->this_subtype = 0;
      mvp->procid = ompcp->proc->procid;
      mvp->proctype = ompcp->proc->proctype;
      mvp->userkey = OMGetNextUserKey ();
      omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid, OMPROC_EDIT, mvp->userkey);
      if (omudp != NULL) {
          omudp->userdata.ptrvalue = (Pointer) mvp;
          omudp->messagefunc =  SeqStrucMediaMsgFunc;
      }
  }

  for(iCount = 1; iCount < Num_Bioseq; iCount++){  /* skip first one */
      slave_entityID = BioseqFindEntity(mediadata[iCount]->sip, &slave_itemID);
      omudp = ObjMgrAddUserData(slave_entityID, ompcp->proc->procid, OMPROC_EDIT, mvp->userkey);
      if(omudp != NULL){
          omudp->userdata.ptrvalue = (Pointer) mvp;
          omudp->messagefunc =  SeqStrucMediaMsgFunc;
      }
  }

  return OM_MSG_RET_DONE;

}
