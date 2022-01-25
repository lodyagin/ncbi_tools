/*   cn3dmatn.c
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
* File Name:  cn3dmatn.c
*
* Author: Yanli Wang
*
* Version Creation Date:   3/26/98
*
* File Description: Functions for cn3d/salsa communication
*
* Modifications:
* $Log: cn3dmatn.c,v $
* Revision 6.50  1998/09/23 22:07:32  ywang
* fix file name error
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
#include <mmdbapi.h>   /* the MMDB-API header */
#include <mmdbapi1.h>
#include <mmdbapi2.h>
#include <mmdbapi3.h>
#include <mmdbapi4.h>
#include <mmdbdata.h>
#include <cn3dmsg.h>
#include <sqnutils.h>

typedef struct {
  Uint1 rgb[3];
} MyColor;

#define CN3D_COLOR_MAX 64
MyColor ObjmgrColorPalette[CN3D_COLOR_MAX] = 
{
  255, 255, 255, /* default     0 */
  255,  20, 147, /* hotpink     1 */
  255,   0, 255, /* magenta     2 */
  155,  48, 255, /* purple      3 */
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

#define Default_Color 16          /* black */

/*----------------------------------------------*/
void SalsaRegister(void)
{

    REGISTER_NEW_BIOSEQ_EDIT;
    REGISTER_NEW_SEQALIGN_EDIT;
    REGISTER_NEW_SEQALIGN_VIEW;
    REGISTER_NEW_SEQANNOT_EDIT;

}
/*----------------------------------------------*/
void Cn3DLaunchAnnotAlignEditor (SeqAnnotPtr sap)
{
  Uint2            entityID,
                   itemID;
  Int2             handled;
  Uint2            options;

  if (sap != NULL) {
     entityID = ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
     options = ObjMgrGetOptions(entityID);
     options |= OM_OPT_FREE_IF_NO_VIEW;
     ObjMgrSetOptions(options, entityID);
     itemID = GetItemIDGivenPointer (entityID, OBJ_SEQANNOT, (Pointer) sap);
     handled = GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, itemID, OBJ_SEQANNOT, 0, 0, OBJ_SEQANNOT, 0);

     sap_entityID = entityID;
     sap_itemID = itemID;
     Cn3D_SalsaOpen = TRUE;
  }
}
/*----------------------------------------------*/
void LaunchSalsa(SeqAlignPtr salp)
{
  SeqAnnotPtr   sap, sap2;

    sap = SeqAnnotNew ();
    sap->type = 2;
    sap->data = (Pointer) salp;
    sap2 = (SeqAnnotPtr) AsnIoMemCopy ((Pointer) sap, (AsnReadFunc) SeqAnnotAsnRead, (AsnWriteFunc) SeqAnnotAsnWrite);
    sap->data = NULL;
    SeqAnnotFree (sap);
    Cn3DLaunchAnnotAlignEditor (sap2);
}
/*----------------------------------------------*/
PMMD FindMM(SeqIdPtr sip, Int4 iCount){

  PDNMS pdnmsThis = NULL, pdnmsThisSlave = NULL;
  PMSD  pmsdThis = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNML pdnmlThis = NULL;

  Boolean MM_found = FALSE;
  Int4 tmp_GI;
  Int4 slaveCount = 0;

  SeqIdPtr tmp_sip;

  Boolean bSingleMS = FALSE, bMaster = FALSE;

  pdnmsThis = GetSelectedModelstruc();
  if(pdnmsThis == NULL) {
/*   ErrPostEx (SEV_ERROR, 0, 0, " GetSelectedModelstruc is NULL: return! "); */
     return(NULL);
  }

  pdnmsThisSlave = ((PMSD)(pdnmsThis->data.ptrvalue))->pdnmsSlaves;
  if(pdnmsThisSlave == NULL) bSingleMS = TRUE;
  else if(iCount == -1){
     bMaster = TRUE;
  }
  else pdnmsThis = pdnmsThisSlave;
 

/*if(iCount == -1) {
      pdnmsThis = GetMasterModelstruc();
      if(pdnmsThis == NULL) {
         bSingleMS = TRUE;
         pdnmsThis = GetSelectedModelstruc();
      }
      else bMaster = TRUE;
      if(pdnmsThis == NULL) {
          return NULL;
      } 
  }
  else
  {
    if(pdnmsThis = GetMasterModelstruc())
    {
      pdnmsThis = ((PMSD)(pdnmsThis->data.ptrvalue))->pdnmsSlaves;
    }
    else {
       bSingleMS = TRUE;
       pdnmsThis = GetSelectedModelstruc();
    }
  }  */

  while (pdnmsThis && !MM_found){
      if(!bSingleMS && !bMaster){
         if(slaveCount != iCount) goto errot;
      }
      pmsdThis = pdnmsThis->data.ptrvalue;
      if(!pmsdThis->bVisible) goto errot;
      pdnmmHead = pmsdThis->pdnmmHead;
      while(pdnmmHead){
          pmmdThis = pdnmmHead->data.ptrvalue;
          if(pmmdThis){
              tmp_sip = pmmdThis->pSeqId;
              if(SeqIdForSameBioseq(tmp_sip, sip)){
                  MM_found = TRUE;
                  return(pmmdThis);
              }
           }

           pdnmmHead = pdnmmHead->next;
       }

      errot:
      pdnmsThis = pdnmsThis->next;    
      if(!bSingleMS) slaveCount++;
  }

  return(NULL);

}
/*----------------------------------------------*/
void DoMediaHL(PMMD  pmmdThis, Int4 from, Int4 to, Boolean highlight)
{
  PDNMG pdnmgThis = NULL;

  PMGD  pmgdThis = NULL;
  PVNMA pvnmaThis = NULL;
  PMAD  pmadThis = NULL;
  PVNAL pvnalThis = NULL;
  PALD  paldThis = NULL;

  Byte MainAtom = 0;


  pdnmgThis = pmmdThis->pdnmgHead;
  if(pdnmgThis == NULL) return;

  pmgdThis = pdnmgThis->data.ptrvalue;
  if ((pmgdThis->bWhat & (Byte) RES_RNA) || (pmgdThis->bWhat & (Byte) RES_DNA)) MainAtom = AM_PALPHA;
  if (pmgdThis->bWhat & (Byte) RES_AA) MainAtom = AM_CALPHA;

  if(from <= 1) from = 1;
  while(pdnmgThis){
      if(pdnmgThis->choice <= to && pdnmgThis->choice >= from){
          pmgdThis = pdnmgThis->data.ptrvalue;
          if(pmgdThis == NULL) goto setout;

                    /* following to deal with the simplest CA-only model */
                    /* similar check should do for DNA/RNA in future */
                    /* GetMainAtom could be used here */

/*        while(pvnmaThis){  */
/*           pmadThis = pvnmaThis->data.ptrvalue;  */
/*           if(StringCmp(pmadThis->pcAName, " CA ") == 0) break;  */
/*           if (pmadThis->bWhat & MainAtom) break;  */
/*           pvnmaThis = pvnmaThis->next;  */
/*        }   */

          fnPreCHLresidue(pdnmgThis, highlight);
      }
      setout:
      pdnmgThis = pdnmgThis->next;
  }

  fnCHLresidueRedraw();

}
/*----------------------------------------------*/
void MediaHL(SelStructPtr sel, Boolean highlight)
{
                                     /* media action -- highlight... */
  SeqLocPtr  slp;
  SeqIntPtr  sintp;
  SeqIdPtr   sip;
  Int4 obj_GI;
  Int4 from, to; 
  Int4 length, iCount = 0;
  Boolean MM_found = FALSE;

Char str[100];

  PMMD  pmmdThis = NULL;

  from = SeqLocStart(sel->region);
  to   = SeqLocStop(sel->region);

  from = from + 1; to = to + 1;
    /* residue starts with number 1 in structure, but with number 0 in sequence */
  slp = sel->region;
  if(slp == NULL) return;
  sintp = slp->data.ptrvalue;
  sip = sintp->id;

  for(iCount = 0; iCount < Num_Bioseq; iCount++){

/** for debugging porpuse **/
SeqIdWrite(sip, str, PRINTID_REPORT, sizeof (str));
SeqIdWrite(mediadata[iCount]->sip, str, PRINTID_REPORT, sizeof (str));

     if (SeqIdForSameBioseq(sip, mediadata[iCount]->sip) )
     {
       if(to == -1){
          length = mediadata[iCount]->length;
          to = length;   /* correct end position, see objmgr APPEND_RESIDUE */
       }
       pmmdThis = FindMM(mediadata[iCount]->sip, iCount - 1);
       if(pmmdThis != NULL){
          DoMediaHL(pmmdThis, from, to, highlight);
       }
     }
  }

}
/*----------------------------------------------*/
void MediaObjSelect(PDNMG pdnmgThis, PALD paldThis, Boolean highlight)
{
  BioseqPtr bspThis;
  PMGD pmgdThis = NULL;
  PMMD pmmdThis = NULL;
  SeqIdPtr  sip;
  SeqLocPtr  slp;
  SelStructPtr sel;
  SeqLocPtr  sel_slp;
  SeqIntPtr  sel_sintp;

  PMSD  pmsdThis = NULL;
  PDNMS pdnmsMaster = NULL;
  PMSD  pmsdMaster = NULL;
  PDNMS pdnmsSlaves = NULL;

  Int4 Gi = 0, iCount = 0;
  Uint2 entityID, itemID, itemtype;
  Int4 from, to;
  Boolean bMaster = FALSE, bSlave = FALSE, RegisterThis = FALSE, bSingleMS = FALSE;
  
  Boolean select_success = FALSE;
  SeqIdPtr sipThis;

  if(!Cn3D_ObjMgrOpen){
       fnPreCHLresidue(pdnmgThis, highlight);
       return;
  }

  if(Num_Bioseq == 0) return;   /* important */

  pmgdThis = pdnmgThis->data.ptrvalue;
  pmmdThis = GetParentMol((PFB)pmgdThis);

  if(pmmdThis == NULL) {
      return;
  }

  sip = pmmdThis->pSeqId;

  pmsdThis = ToMSDParent((PFB) pmgdThis);

  pdnmsMaster = GetSelectedModelstruc();
           /* always use GetSelectedMOdelstruc */
  if(pdnmsMaster){
     if(pmsdThis == pdnmsMaster->data.ptrvalue) {
        pmsdMaster = pdnmsMaster->data.ptrvalue;
        if(pmsdMaster->pdnmsSlaves != NULL) bMaster = TRUE;
        else bSingleMS = TRUE;
     }
     else {
        pmsdMaster = pdnmsMaster->data.ptrvalue;
        pdnmsSlaves = pmsdMaster->pdnmsSlaves;

        iCount = 0;

        while(pdnmsSlaves) {
           if(pmsdThis == pdnmsSlaves->data.ptrvalue) {
              bSlave = TRUE;
              break;
           }
           iCount++;
           pdnmsSlaves = pdnmsSlaves->next;
        }
      }
  }

  if(!bSingleMS){
     if(!bMaster && !bSlave) {
        fnPreCHLresidue(pdnmgThis, highlight);
        return;
     }
     if(bMaster) {
        if(SeqIdForSameBioseq(sip, mediadata[0]->sip)){
           sipThis = mediadata[0]->sip;        
           RegisterThis = TRUE;
        }
     }
     else{
        iCount = iCount + 1;
        if(SeqIdForSameBioseq(sip, mediadata[iCount]->sip)) {      
           sipThis = mediadata[iCount]->sip;
           RegisterThis = TRUE;
                  /* get rid of Gi, use SeqId to match Seqeuence and Structure */
        }        
    }
  }
  else {
     for(iCount = 0; iCount < Num_Bioseq; iCount++){
        if(SeqIdForSameBioseq(sip, mediadata[iCount]->sip)) {
           sipThis = mediadata[iCount]->sip;
           RegisterThis = TRUE;
           break;
                     /* could get rid of Gi later on, but since for */
                     /* one struc-seq from MMDB Gi should be there for */
                     /* protein/NA, let Gi be here for the moment */
        }
     }    
  }
 
  if(!RegisterThis){
      fnPreCHLresidue(pdnmgThis, highlight);
      return;
  }

  entityID = BioseqFindEntity(sipThis, &itemID);
  itemtype = OBJ_BIOSEQ;
  from = pdnmgThis->choice; to = pdnmgThis->choice;
  from = from -1; to = to -1;   
  /* residue starts with number 1 in structure, but with number 0 in sequence */
  sel = ObjMgrGetSelected();
  if(sel == NULL) {
      slp = SeqLocIntNew(from, to, Seq_strand_unknown, SeqIdDup(sipThis));
      select_success = ObjMgrSelect(entityID, itemID, itemtype, OM_REGION_SEQLOC, slp);
  }
  else {
      while(sel != NULL) {
         if(sel->entityID == entityID && sel->itemID == itemID && sel->itemtype == itemtype){
            sel_slp = (SeqLocPtr)sel->region;
            if(sel_slp != NULL && sel_slp->data.ptrvalue !=NULL) {
               if(SeqLocStop(sel_slp) == -2){
                  if(from >= SeqLocStart(sel_slp)){
                     select_success = ObjMgrDeSelect(entityID, itemID, itemtype, OM_REGION_SEQLOC, sel->region);
                     return;
                  }
               }
               else if(to <= SeqLocStop(sel_slp)  && from >= SeqLocStart(sel_slp))
               {
                  select_success = ObjMgrDeSelect(entityID, itemID, itemtype, OM_REGION_SEQLOC, sel->region);
                  return;
               }
            }
         }
         sel = sel->next;
      }
      slp = SeqLocIntNew(from, to, Seq_strand_unknown, SeqIdDup(sipThis));
      select_success = ObjMgrAlsoSelect(entityID, itemID, itemtype, OM_REGION_SEQLOC, slp);
  }

}
/*-----------------------------------------------*/
Uint1Ptr GetRGB(Int2 iColor)
{
  Uint1Ptr rgb;

  rgb = (ObjmgrColorPalette + iColor)->rgb;
 
  return(rgb);

}
/*-----------------------------------------------*/
void ColorSalsa(Uint2 entityID, Uint2 itemID, SeqIdPtr sip, Int4 from, Int4 to, Uint1Ptr rgb)
{
  SelStructPtr  sel;
  SeqIdPtr      sip_dup;

  if(Num_Bioseq == 0) return;    /* important */

  sel = (SelStructPtr)MemNew((size_t)sizeof(SelStruct));
  if(sel != NULL) {
     sel->entityID = entityID;
     sel->itemtype = OBJ_BIOSEQ;
     sel->itemID = itemID;
     sel->regiontype = OM_REGION_SEQLOC;
     sip_dup = SeqIdDup(sip);
     sel->region = (SeqLocPtr)SeqLocIntNew(from, to, Seq_strand_unknown, sip_dup);

     ObjMgrSetColor(sel->entityID, sel->itemID, sel->itemtype, sel->regiontype, sel->region,  rgb);

     sel->next = NULL;
     if(sel->region != NULL) SeqLocFree ((SeqLocPtr) sel->region);
     sel = MemFree(sel);

     sip_dup = SeqIdFree(sip_dup);
  }
}
/*-----------------------------------------------*/
void PrepareColorMsg(Int4 iCount, Int4 from, Int4 to, Uint1Ptr rgb)
{

  SeqIdPtr  sip;
  Uint2 entityID, itemID;

  if(Num_Bioseq == 0) return; /* important */

  sip = mediadata[iCount]->sip;
  entityID = mediadata[iCount]->entityID;
  itemID = mediadata[iCount]->itemID;

  ColorSalsa(entityID, itemID, sip, from, to, rgb);

}
/*-----------------------------------------------*/
void ColorSalsa_BYMG(PMGD pmgdThis, Uint1Ptr rgb)
{
  SelStructPtr  sel;
  Int4  iCount = 0, Gi = 0;
  Int4  from = 0, to = 0;

  Boolean bMaster = FALSE, bSingleMS = FALSE;

  PMMD pmmdThis = NULL;
  PDNMG pdnmgThis = NULL;
  PDNMS pdnmsMaster = NULL;
  PMSD  pmsdMaster = NULL;
  PMSD  pmsdThis = NULL;

  PDNMS pdnmsSlaves = NULL;
 
  SeqIdPtr sip;

  if(!Cn3D_ObjMgrOpen) return;

  if(Num_Bioseq == 0) return; /* important */

  pmsdThis = ToMSDParent((PFB) pmgdThis);

  pdnmsMaster = GetSelectedModelstruc();
              /* always use GetSelectedModelstruc */

  if(pdnmsMaster){
     if(pmsdThis == pdnmsMaster->data.ptrvalue) {
        pmsdMaster = pdnmsMaster->data.ptrvalue;
        if(pmsdMaster->pdnmsSlaves != NULL) bMaster = TRUE;
        else bSingleMS = TRUE;
     }
     else {
        pmsdMaster = pdnmsMaster->data.ptrvalue;
        pdnmsSlaves = pmsdMaster->pdnmsSlaves;

        iCount = 0;

        while(pdnmsSlaves) {
           if(pmsdThis == pdnmsSlaves->data.ptrvalue) {
              break;
           }
           iCount++;
           pdnmsSlaves = pdnmsSlaves->next;
        }
     }
  }

  pmmdThis = GetParentMol((PFB)pmgdThis);
  if(pmmdThis == NULL) {
     return;
  }
  sip = pmmdThis->pSeqId;
  if(sip == NULL) {
/*   ErrPostEx (SEV_ERROR, 0, 0, " SeqId is NULL: return! ");  */
     return;
  }

  pdnmgThis = pmgdThis->pdnmgLink;
  from = pdnmgThis->choice - 1;
  to   = pdnmgThis->choice - 1;
 
  if(!bSingleMS){
     if(bMaster) {
         if(SeqIdForSameBioseq(sip, mediadata[0]->sip)){
            PrepareColorMsg(0, from, to, rgb);
                     /* here the order and SeqId matters */
         }
     }
     else{
        iCount = iCount + 1;
        if(SeqIdForSameBioseq(sip, mediadata[iCount]->sip)) {
           PrepareColorMsg(iCount, from, to, rgb);
                    /* here the order and SeqId matters */
         }
     }
  }
  else{
     for(iCount = 0; iCount < Num_Bioseq; iCount++){
        if(SeqIdForSameBioseq(sip, mediadata[iCount]->sip)) {
           PrepareColorMsg(iCount, from, to, rgb);
           break;      /* could get rid of Gi later on */
                       /* but since for one biostruc-seq (from MMDB) case, */
                       /* Gi should be there for Rna/Protein, let it be there */
        }
     }
  }
         
}
/*-----------------------------------------------*/
void ResetSalsaColor(void)
{
  Int4 from = 0, to = 0;
  Int4 length = 0, iCount = 0;
  Uint2 entityID, itemID;
  SeqIdPtr sip;
  Uint1Ptr rgb;

  for(iCount = 0; iCount < Num_Bioseq; iCount++){
     if(!mediadata[iCount]->bVisible) continue;
     sip = mediadata[iCount]->sip;
     length = (Int4) mediadata[iCount]->length;
     from = 0; to = length - 1;
     entityID = mediadata[iCount]->entityID;
     itemID = mediadata[iCount]->itemID;
     rgb = (Uint1Ptr) GetRGB( Default_Color );  /* black--default color */
     
     ColorSalsa(entityID, itemID, sip, from, to, (Uint1Ptr)rgb); 
  }

}
/*-----------------------------------------------*/
void Cn3dObjMgrGetSelected(void)
{
  SelStructPtr sel = NULL;
  Boolean highlight = FALSE;

  sel = ObjMgrGetSelected();
  while(sel != NULL) {
     highlight = TRUE;
     MediaHL(sel, highlight);  /* highlight in cn3d */
     sel = sel->next;
  }

}
