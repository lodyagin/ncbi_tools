/* ==========================================================================
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
* =========================================================================== 
*
* File Name: cn3dmodl.c
*
* Author: Yanli Wang
*
* Version Creation Date: 18/9/1998
*
* $Log: cn3dmodl.c,v $
* Revision 6.20  1999/01/20 18:21:19  ywang
* include salmedia.h due to the move around of MediaInfo from cn3dmsg.h to the new created salmedia.h
*
 * Revision 6.19  1998/12/16  22:49:37  ywang
 * fix compiling warnings on Win32
 *
 * Revision 6.17  1998/11/09  22:03:31  ywang
 * fix bugs for modeling
 *
 * Revision 6.16  1998/11/06  23:01:05  ywang
 * fix bugs for modeling
 *
 * Revision 6.15  1998/11/04  13:03:30  kans
 * fixed variable name
 *
* Revision 6.14  1998/11/04 00:06:19  ywang
* add function for modeling: change render/color for special residue(s)
*
 * Revision 6.13  1998/10/21  15:48:23  ywang
 * rearrange View Control menu
 *
 * Revision 6.12  1998/10/20  17:50:37  ywang
 * fix MMDB/VAST domain number inconsistency problem
 *
 * Revision 6.11  1998/10/20  16:54:39  ywang
 * fix bug in AssignDomainAlignedStatus
 *
 * Revision 6.10  1998/10/13  14:44:17  ywang
 * fix bugs for single struc case in AssignDomainAlignedStatus
 *
 * Revision 6.9  1998/10/08  21:39:15  ywang
 * fix bug for assign domain aligned status
 *
 * Revision 6.8  1998/10/08  00:07:19  ywang
 * fix domain aligned status bug
 *
 * Revision 6.6  1998/10/02  18:25:15  kans
 * some functions did not have void as a return, which the Mac compiler detected
 *
* Revision 6.5  1998/10/01 21:55:39  ywang
* put function for object display control
*
 * Revision 6.4  1998/09/30  22:10:46  ywang
 * control display on three levels: structure, chain, domain
 *
 * Revision 6.3  1998/09/23  22:04:02  ywang
 * synchronize show/hide between cn3d and salsa when display complexity is changed
 *
 * Revision 6.2  1998/09/23  18:38:48  ywang
 * add functions to control display on domain level
 *
 * Revision 6.1  1998/09/22  18:02:54  ywang
 * panels and functions for display control
 *
*/

#include <ncbi.h>
#include <viewer3d.h>
#include <cn3dmain.h>
#include <math.h>
#include <mmdbapi.h>
#include <mmdbapi1.h>
#include <mmdbapi2.h>
#include <mmdbapi3.h>
#include <mmdbapi4.h>
#include <cn3dpane.h>
#include <algorend.h>
#include <cn3dmsg.h>
#include <salmedia.h>
#include <salutil.h>
#include <cn3dopen.h>
#include <cn3dmodl.h>
#include <cn3dpane.h>

static LisT     Cn3D_lModelOnOff;        /* pieces parts to draw */
static PopuP    Cn3D_pupModelPBB;        /* protein backbone options */
static PopuP    Cn3D_pupModelNABB;       /* nucl. acid bb options */
static PopuP    Cn3D_pupModelStyleItem;
static PopuP    Cn3D_pupModelRenderStyle;
static PopuP    Cn3D_pupModelColorStyle;
static TexT     FeatureTitle; 
static ButtoN   bResetFeatureTitle; 

static ButtoN Cn3D_bDisplayApply;
static ButtoN Cn3D_bModelApply;

static ButtoN Cn3D_bDisplayByStruc;
static LisT Cn3D_lStruc;

static GrouP Cn3D_gDisplayAlign;
static ButtoN  Cn3D_bAlignOn, Cn3D_bUnalignOn;
static ButtoN Cn3D_bDisplayAlignedDomain;

static ButtoN Cn3D_bDisplayDomain;
static ButtoN Cn3D_bDisplayOthers;
static LisT Cn3D_lFeature;

Int4 iDomainCount = 0;
DomainInfo **domaindata;
extern Boolean Cn3D_ReColor;

SpecialFeaturePtr sfp = NULL;
PARS parsSpecial = NULL;
Byte iAddedFeature = 0;

/*---------------------------------------------------------*/
void GetAlignStatus(void)
{
  Cn3D_fAlignOn = GetStatus(Cn3D_bAlignOn);
  Cn3D_fUnalignOn = GetStatus(Cn3D_bUnalignOn);
}
/*---------------------------------------------------------*/
static void Cn3D_ListOthersProc(ButtoN b)
{

  SpecialFeaturePtr sfpThis = NULL;
  SpecialFeatureInfoPtr sfipThis = NULL;

/*SetStatus(Cn3D_bDisplayOthers, FALSE);
  Message(MSG_OK, "This is assumed to display a list of user specified features! ");  */

  SetStatus(Cn3D_bDisplayDomain, FALSE);
  Reset(Cn3D_lFeature);
  if(!GetStatus(Cn3D_bDisplayOthers)) return;      

  sfpThis = sfp;
  while(sfpThis){
     sfipThis = sfpThis->data.ptrvalue;
     ListItem(Cn3D_lFeature, (CharPtr) sfipThis->title);
     SetItemStatus(Cn3D_lFeature, sfpThis->choice, sfipThis->On);
    
     sfpThis = sfpThis->next;
  }

}
/*---------------------------------------------------------*/
static void Cn3D_ListDomainProc(ButtoN b)
{
  Char  MolName[20];
  Int4 iCount = 0, iCountActive = 0;
  Int4 i;

  SetStatus(Cn3D_bDisplayOthers, FALSE);
  Reset(Cn3D_lFeature);

  if(!GetStatus(Cn3D_bDisplayDomain)) return;

  for(iCount = 0; iCount <iDomainCount; iCount++){
        domaindata[iCount]->bVisibleParent = GetItemStatus(Cn3D_lStruc, domaindata[iCount]->iStrucIndex + 1);
  }

  for(iCount = 0; iCount <iDomainCount; iCount++){
     if(domaindata[iCount]->iDomain == 0) sprintf(MolName, "%4s  %2s", domaindata[iCount]->pcPDBName, domaindata[iCount]->pcMolName);
     else sprintf(MolName, "%4s  %2s  %2d", domaindata[iCount]->pcPDBName, domaindata[iCount]->pcMolName, domaindata[iCount]->iDomain);
     if(domaindata[iCount]->bVisibleParent) {
        ListItem(Cn3D_lFeature, (CharPtr) MolName);
        iCountActive++;
        if(!GetStatus(Cn3D_bDisplayAlignedDomain)) SetItemStatus(Cn3D_lFeature, iCountActive, domaindata[iCount]->bVisible);      
        else {
           if(domaindata[iCount]->bAligned) {
              SetItemStatus(Cn3D_lFeature, iCountActive, TRUE);
              domaindata[iCount]->bVisible = TRUE;
           }
           else if(!domaindata[iCount]->bAligned) {
              SetItemStatus(Cn3D_lFeature, iCountActive, FALSE);
              domaindata[iCount]->bVisible = FALSE;
           }
        }
     }

  }

  return;
}
/*---------------------------------------------------------*/
static void Cn3D_ListAlignedDomainProc(ButtoN b)
{
  ButtoN b1;
  Int4 iCount;

  Cn3D_ListDomainProc(b1);
}
/*---------------------------------------------------------*/
static void FillFeatureListProc(LisT l)
{
  ButtoN b;

  if(GetStatus(Cn3D_bDisplayDomain)) Cn3D_ListDomainProc(b);
  else if(GetStatus(Cn3D_bDisplayOthers)) Cn3D_ListOthersProc(b);
}
/*---------------------------------------------------------*/
static void CheckSalsaDisplay(PMMD pmmdThis, Int4 iCountStruc, Byte bVisible)
{
   SeqIdPtr sipThis = NULL;
   Uint2 entityID, itemID, itemtype;

   sipThis = pmmdThis->pSeqId;

   if(sipThis){
      if(SeqIdForSameBioseq(sipThis, mediadata[iCountStruc]->sip)){
         entityID = BioseqFindEntity(sipThis, &itemID);
         itemtype = OBJ_BIOSEQ;
         if(bVisible) {
            ObjMgrSendMsg(OM_MSG_SHOW, entityID, itemID, itemtype);
         }
         else if(!bVisible) ObjMgrSendMsg(OM_MSG_HIDE, entityID, itemID, itemtype);
         mediadata[iCountStruc]->bVisible = bVisible;
      }
   }

   return;
}
/*---------------------------------------------------------*/
static void SetDomainParentStatus(Int4 Index, Byte bVisible)
{
  Int4 iCount = 0;

  for(iCount = 0; iCount < iDomainCount; iCount++){
     if(domaindata[iCount]->iStrucIndex == Index)domaindata[iCount]->bVisibleParent = bVisible;
  }

}
/*---------------------------------------------------------*/
static void Cn3D_ReSetModVisibleStatus( PVNMO pvnmoThis, Byte bVisible)
{
  PMOD  pmodThis = NULL;

  while(pvnmoThis){
     pmodThis = (PMOD) pvnmoThis->data.ptrvalue;
     pmodThis->bVisible = bVisible; 

     pvnmoThis = pvnmoThis->next;
  }

}
/*---------------------------------------------------------*/
static void Cn3D_ReSetMmVisibleStatus(PDNMM pdnmmHead, Byte bVisible)
{
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg;
  PMGD  pmgd;

  while(pdnmmHead){
      pmmdThis = pdnmmHead->data.ptrvalue;
      if(pmmdThis){
         if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout1;
         pmmdThis->bVisible = 1;

         for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {              pmgd = (PMGD) pdnmg->data.ptrvalue;
             pmgd->bVisible = 1;
         }
      }
      setout1:
      pdnmmHead = pdnmmHead->next;
  }

}
/*---------------------------------------------------------*/
static void Cn3D_ReSetVisibleStatus(void)
{
  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg;
  PMGD  pmgd;

  PVNMO pvnmoThis = NULL;

  pdnmsMaster = GetSelectedModelstruc();
  if (pdnmsMaster != NULL)
  {
    pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;

    pmsdMaster->bVisible = 1;

    pvnmoThis = pmsdMaster->pvnmoHead;
    Cn3D_ReSetModVisibleStatus(pvnmoThis, (Byte) pmsdMaster->bVisible);

    pdnmmHead = pmsdMaster->pdnmmHead;
    Cn3D_ReSetMmVisibleStatus(pdnmmHead, pmsdMaster->bVisible);

    pdnmsSlave = pmsdMaster->pdnmsSlaves;
    while(pdnmsSlave)
    {
      pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;

      pmsdSlave->bVisible = 1;

      pvnmoThis = pmsdSlave->pvnmoHead;
      Cn3D_ReSetModVisibleStatus(pvnmoThis, pmsdSlave->bVisible = 1);

      pdnmmHead = pmsdSlave->pdnmmHead;
      Cn3D_ReSetMmVisibleStatus(pdnmmHead, pmsdSlave->bVisible = 1);

      pdnmsSlave = pdnmsSlave->next;
    }
  }

}
/*---------------------------------------------------------*/
static void ReSetSalsaDisplay(void)
{
   SeqIdPtr sipThis = NULL;
   Uint2 entityID, itemID, itemtype;

   Int4 iCount = 0;

   for(iCount = 0; iCount < Num_Bioseq; iCount++){
      if(!mediadata[iCount]->bVisible){
         itemtype = OBJ_BIOSEQ;
         ObjMgrSendMsg(OM_MSG_SHOW, mediadata[iCount]->entityID, mediadata[iCount]->itemID, itemtype); 
         mediadata[iCount]->bVisible = 1;
       }
   }

}
/*---------------------------------------------------------*/
static void Cn3D_ListStrucProc(ButtoN b)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster, pmsdSlave;

  Int4 iCount;

  Reset(Cn3D_lStruc);

  Cn3D_ReSetVisibleStatus();

  ReSetSalsaDisplay();

  iCount = 0;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){
     pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
     ListItem(Cn3D_lStruc, pmsdMaster->pcPDBName);
     SetItemStatus(Cn3D_lStruc, iCount + 1, pmsdMaster->bVisible);  

     iCount++;

     pdnmsSlave = pmsdMaster->pdnmsSlaves;
     while(pdnmsSlave) {
        pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;
        ListItem(Cn3D_lStruc, pmsdSlave->pcPDBName);
        SetItemStatus(Cn3D_lStruc, iCount + 1, pmsdSlave->bVisible);     

        iCount++;

        pdnmsSlave = pdnmsSlave->next;
     }
  }
     
  Salsa_BioseqUpdate = FALSE;
  Cn3D_ReColor = TRUE;
  Cn3D_RedrawProc(FALSE);
  Cn3D_ReColor = FALSE;

  return;
}
/*---------------------------------------------------------*/
void AssignDomainAlignedStatus(void)
{
  Int4 iCount = 0, iCountStruc = 0;
  PDNMS pdnmsMaster = NULL;
  PMSD  pmsdMaster = NULL;
  BiostrucAnnotSetPtr pbsaThis = NULL;
  BiostrucFeatureSetPtr pbsfsThis = NULL;
  BiostrucFeaturePtr pbsfThis = NULL;

  Int2 iDomainMaster = 0, iDomainSlave = 0;

  Char iDomain[2];

  pdnmsMaster = GetSelectedModelstruc();
  if(!pdnmsMaster) return;
  else pmsdMaster = pdnmsMaster->data.ptrvalue;

  pbsaThis = pmsdMaster->psaStrucAlignment;
  if(pbsaThis == NULL) return;
  pbsfsThis = pbsaThis->features; 
  iDomainMaster = (Int2) (pbsfsThis->id % 100);

  for(iCount = 0; iCount < iDomainCount; iCount++){
     pbsfThis =  pbsfsThis->features;
     iCountStruc = 0;
     while(pbsfThis){
        iDomainSlave = (Int2) (((pbsfThis->id - pbsfThis->id % 10) / 10 ) % 100 );
        iCountStruc++;

        if(domaindata[iCount]->iStrucIndex == 0){
           if((StringNCmp(pbsfThis->name, domaindata[iCount]->pcPDBName, 4)
== 0) && (pbsfThis->name[4] == domaindata[iCount]->pcMolName[0])) {
/*            if((Char) pbsfThis->name[5] == '0'){   */
              if(iDomainMaster == 0){
                 domaindata[iCount]->bAligned = TRUE;
              }
              else if (iDomainMaster == domaindata[iCount]->iDomain) {
                 domaindata[iCount]->bAligned = TRUE;
              }
           }
        }
        else if(domaindata[iCount]->iStrucIndex == iCountStruc){
           if((StringNCmp(pbsfThis->name + 7, domaindata[iCount]->pcPDBName, 4) == 0) && (pbsfThis->name[11] == domaindata[iCount]->pcMolName[0])){
              if(iDomainSlave == 0){
                 domaindata[iCount]->bAligned = TRUE;
              }
              else if ( iDomainSlave == domaindata[iCount]->iDomain) {     
              domaindata[iCount]->bAligned = TRUE;
              }
           }
        }
      
        pbsfThis = pbsfThis->next;
     }
  }

}
/*---------------------------------------------------------*/
void LIBCALL Cn3D_CountDomainProc(void)
{
  Int4 iCount = 0;
  Int4 idom = 0;
  Int4 iCountStruc = 0, iCountChain = 0;

  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg;
  PMGD  pmgd;

  Boolean EncounterZeroDomain = FALSE;

  pdnmsMaster = GetSelectedModelstruc();
  if (pdnmsMaster != NULL)
  {
    pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
    pdnmmHead = pmsdMaster->pdnmmHead;

    while(pdnmmHead){
       pmmdThis = pdnmmHead->data.ptrvalue;
       if(pmmdThis){
          if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto pre_setout1;

          idom = 0;

          EncounterZeroDomain = FALSE;
          for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
              pmgd = (PMGD) pdnmg->data.ptrvalue;
              if(!EncounterZeroDomain){
                 if(pmgd->iDomain == 0) {
                    iCount++;
                    EncounterZeroDomain = TRUE;     
                 }
              }
              if(pmgd->iDomain > idom){
                    idom = (Int4) pmgd->iDomain;
                    iCount++;
              }
          }

       }
       pre_setout1:
       pdnmmHead = pdnmmHead->next;
    }
 
    pdnmsSlave = pmsdMaster->pdnmsSlaves;
    while(pdnmsSlave)
    {
      pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;
      pdnmmHead = pmsdSlave->pdnmmHead;

      while(pdnmmHead){
         pmmdThis = pdnmmHead->data.ptrvalue;
         if(pmmdThis){
            if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto pre_setout2;

            idom = 0;

            EncounterZeroDomain = FALSE;

            for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
              pmgd = (PMGD) pdnmg->data.ptrvalue;
              if(!EncounterZeroDomain){
                 if(pmgd->iDomain == 0) {
                    iCount++;
                    EncounterZeroDomain = TRUE;
                 }
              }
               pmgd = (PMGD) pdnmg->data.ptrvalue;
               if(pmgd->iDomain > idom){
                  idom = (Int4) pmgd->iDomain;
                  iCount++;
               }
            }

         }
         pre_setout2:
         pdnmmHead = pdnmmHead->next;
      }

      pdnmsSlave = pdnmsSlave->next;
    }
  }

  iDomainCount = iCount;

  /* to allocate memory for domaindata */

  domaindata = (Pointer) MemNew((size_t) (iDomainCount * sizeof(Pointer)));
  for(iCount = 0; iCount < iDomainCount; iCount++){
     domaindata[iCount] = (DomainInfoPtr) MemNew(sizeof (DomainInfo));
  }

  /* to load pdbname, chain id, domain id into domaindata */
  
  iCount = 0; iCountStruc = 0; iCountChain = 0;

  pdnmsMaster = GetSelectedModelstruc();
  if (pdnmsMaster != NULL)
  {
    pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
    pdnmmHead = pmsdMaster->pdnmmHead;

    while(pdnmmHead){
       pmmdThis = pdnmmHead->data.ptrvalue;
       if(pmmdThis){
          if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout1;
          pmmdThis->bVisible = 1;

          idom = 0;

          EncounterZeroDomain = FALSE;

          for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
              pmgd = (PMGD) pdnmg->data.ptrvalue;
              if(!EncounterZeroDomain){
                 if(pmgd->iDomain == 0){
                    StringCpy(domaindata[iCount]->pcPDBName, pmsdMaster->pcPDBName);
                    StringCpy(domaindata[iCount]->pcMolName, pmmdThis->pcMolName);
                    domaindata[iCount]->iDomain = 0;
                    domaindata[iCount]->iStrucIndex = iCountStruc;
                    domaindata[iCount]->iChainIndex = iCountChain;
                    domaindata[iCount]->bVisible = TRUE;
                    domaindata[iCount]->bVisibleParent = TRUE;
                    iCount++;
                    EncounterZeroDomain = TRUE;
                 }
              }
              if(pmgd->iDomain > idom){
                 idom = (Int4) pmgd->iDomain;
                 StringCpy(domaindata[iCount]->pcPDBName, pmsdMaster->pcPDBName);
                 StringCpy(domaindata[iCount]->pcMolName, pmmdThis->pcMolName);
                 domaindata[iCount]->iDomain = pmgd->iDomain;
                 domaindata[iCount]->iStrucIndex = iCountStruc;
                 domaindata[iCount]->iChainIndex = iCountChain;
                 domaindata[iCount]->bVisible = TRUE;
                 domaindata[iCount]->bVisibleParent = TRUE;
                    
                 iCount++;
              }
          }
          iCountChain++;
       }
       setout1:
       pdnmmHead = pdnmmHead->next;
    }

    iCountStruc++;
    pdnmsSlave = pmsdMaster->pdnmsSlaves;
    while(pdnmsSlave)
    {
      pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;
      pdnmmHead = pmsdSlave->pdnmmHead;

      while(pdnmmHead){
         pmmdThis = pdnmmHead->data.ptrvalue;
         if(pmmdThis){
            if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout2;

            idom = 0;

            EncounterZeroDomain = FALSE;

            for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
               pmgd = (PMGD) pdnmg->data.ptrvalue;
               if(!EncounterZeroDomain){
                  if(pmgd->iDomain == 0){
                     StringCpy(domaindata[iCount]->pcPDBName, pmsdSlave->pcPDBName);
                     StringCpy(domaindata[iCount]->pcMolName, pmmdThis->pcMolName);
                     domaindata[iCount]->iDomain = 0;
                     domaindata[iCount]->iStrucIndex = iCountStruc;
                     domaindata[iCount]->iChainIndex = iCountChain;
                     domaindata[iCount]->bVisible = TRUE;
                     domaindata[iCount]->bVisibleParent = TRUE;
                     iCount++;
                     EncounterZeroDomain = TRUE;
                   }
               }
               if(pmgd->iDomain > idom){
                  idom = (Int4) pmgd->iDomain;
                  StringCpy(domaindata[iCount]->pcPDBName, pmsdSlave->pcPDBName);
                  StringCpy(domaindata[iCount]->pcMolName, pmmdThis->pcMolName);
                  domaindata[iCount]->iDomain = pmgd->iDomain;
                  domaindata[iCount]->iStrucIndex = iCountStruc;
                  domaindata[iCount]->iChainIndex = iCountChain;
                  domaindata[iCount]->bVisible = TRUE;
                  domaindata[iCount]->bVisibleParent = TRUE;

                  iCount++;
               }
/*             if(pmgd->bReserved && !domaindata[iCount - 1]->bAligned) domaindata[iCount - 1]->bAligned = TRUE;  */
            }
            iCountChain++;
         }
         setout2:
         pdnmmHead = pdnmmHead->next;
      }

      iCountStruc++;
      pdnmsSlave = pdnmsSlave->next;
    }
  }

  AssignDomainAlignedStatus();

}
/*---------------------------------------------------------*/
void SynchronizeModVisibleStatus(PFB pfbThis)
{
   PMSD pmsdThis = NULL;
   PMMD pmmdThis = NULL;
   PDNMG pdnmg;
   PMGD pmgdThis = NULL;

   PVNMO pvnmoThis = NULL;
   PMOD pmodThis = NULL;
   
   ValNodePtr pvnThis = NULL;

   pmsdThis = ToMSDParent(pfbThis); 
   pvnmoThis = pmsdThis->pvnmoHead;

   if(pfbThis->bMe == AM_MMD) {
      pmmdThis = (PMMD) pfbThis;
      while(pvnmoThis){
         pmodThis = pvnmoThis->data.ptrvalue;
         if(IsGraphNode((PFB) pmodThis->pvnContains->data.ptrvalue)){
            if(pmmdThis == GetParentMol((PFB) pmodThis->pvnContains->data.ptrvalue)){
               pmodThis->bVisible = pmmdThis->bVisible;
            }
         }

         pvnmoThis = pvnmoThis->next;
      }
   } 
   else if(pfbThis->bMe == AM_MGD){
       while(pvnmoThis){
          pmodThis = pvnmoThis->data.ptrvalue;
          pmgdThis = pmodThis->pvnContains->data.ptrvalue;
          pmodThis->bVisible = pmgdThis->bVisible;

          pvnmoThis = pvnmoThis->next;
       }
   }
}
/*---------------------------------------------------------*/
static void SetDomainDataItemStatus(LisT l)
{
  Int4 iCount = 0, iCountActive = 0;

  if(GetStatus(Cn3D_bDisplayDomain)){
     for(iCount = 0; iCount < iDomainCount; iCount++){
        if(domaindata[iCount]->bVisibleParent) {
           domaindata[iCount]->bVisible = GetItemStatus(Cn3D_lFeature, iCountActive + 1);
           iCountActive++;
        }
     }
  }

  return; 
}
/*---------------------------------------------------------*/
Int4 iCountItemGet(PMSD pmsdThis, PMMD pmmdThis, PMGD pmgdThis)
{
  Int4 iCount = 0, iCountLive = 0;
  Int2 iDomain = 0;

  if(pmgdThis == NULL) iDomain = 0;
  else iDomain = pmgdThis->iDomain;

  for(iCount = 0; iCount < iDomainCount; iCount++){
     if(!domaindata[iCount]->bVisibleParent) continue;
     if(StringCmp(pmsdThis->pcPDBName, domaindata[iCount]->pcPDBName) == 0 && StringCmp(pmmdThis->pcMolName, domaindata[iCount]->pcMolName) ==0 && iDomain ==
domaindata[iCount]->iDomain) return (iCountLive + 1) ;
     iCountLive++;
    
  }

  return 0;
}
/*---------------------------------------------------------*/
static void fnAlignList(LisT l)
/* set the values of the alignment pane */
{
  Int4 iCount = 0;
  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;
  BiostrucFeaturePtr pbsfThis;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  Byte bVisible;
  SeqIdPtr sipThis;
  Uint2 entityID, itemID, itemtype;         

  ButtoN b;

  Cn3D_SaveActiveCam();

  iCount = 0;
  pdnmsMaster = GetSelectedModelstruc();
  if (pdnmsMaster != NULL)
  {

    pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;

    bVisible = pmsdMaster->bVisible;
    pmsdMaster->bVisible = GetItemStatus(Cn3D_lStruc, iCount + 1);
    SetDomainParentStatus(iCount, pmsdMaster->bVisible);

    pdnmmHead = pmsdMaster->pdnmmHead;
    Cn3D_ReSetMmVisibleStatus(pdnmmHead, pmsdMaster->bVisible);

    while(pdnmmHead){
       pmmdThis = pdnmmHead->data.ptrvalue;
       if(pmmdThis){
          if(pmmdThis->bWhat == AM_PROT || pmmdThis->bWhat == (Byte) AM_RNA || pmmdThis->bWhat == AM_DNA){
            if(Num_Bioseq != 0 && bVisible != pmsdMaster->bVisible) CheckSalsaDisplay(pmmdThis, iCount, pmsdMaster->bVisible);
          }
       }
       pdnmmHead = pdnmmHead->next;
    }

    TraverseGraphs( pdnmsMaster, 0, 0, NULL, fnClearMarkedResidues);
    pmsdMaster->bAligned = 0;
    pdnmsSlave = pmsdMaster->pdnmsSlaves;

    if(pmsdMaster->psaStrucAlignment != NULL) { 
       pbsfThis =  pmsdMaster->psaStrucAlignment->features->features;
       iCount++;
       while(pdnmsSlave)
       {
         pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;

         bVisible = pmsdSlave->bVisible;   

         TraverseGraphs( pdnmsSlave, 0, 0, NULL, fnClearMarkedResidues);

         pmsdSlave->bVisible = GetItemStatus(Cn3D_lStruc, iCount + 1);
         SetDomainParentStatus(iCount, pmsdSlave->bVisible);

         if (pmsdMaster->bVisible && pmsdSlave->bVisible) {
            fnMarkAlignedResidues(pdnmsMaster, pdnmsSlave, pbsfThis);
            pmsdMaster->bAligned++;
         }
         pbsfThis = pbsfThis->next;

         if(bVisible != pmsdSlave->bVisible){
            pdnmmHead = pmsdSlave->pdnmmHead;
            Cn3D_ReSetMmVisibleStatus(pdnmmHead, pmsdSlave->bVisible);
            while(pdnmmHead){
               pmmdThis = pdnmmHead->data.ptrvalue;
               if(pmmdThis){
                  if(pmmdThis->bWhat == AM_PROT || pmmdThis->bWhat == (Byte) AM_RNA || pmmdThis->bWhat == AM_DNA){
                  if(Num_Bioseq != 0 && bVisible != pmsdSlave->bVisible) CheckSalsaDisplay(pmmdThis, iCount, pmsdSlave->bVisible);
                  }
               }

               pdnmmHead = pdnmmHead->next;

            }
         }                                      

         iCount++;
         pdnmsSlave = pdnmsSlave->next;
       }
    }
  }

  Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */

}
/*---------------------------------------------------------*/
static void Cn3D_DisplayProc(ButtoN b )
{
  GrouP g;
  Int2 val;

  Int4 iCountItem = 0, iCountStruc = 0;
  Int4 idom = 0;

  SpecialFeaturePtr sfpThis;
  SpecialFeatureInfoPtr sfipThis = NULL;

  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg;
  PMGD  pmgd;
  PVNMO pvnmoThis;

  SeqIdPtr sipThis;
  Uint2 entityID, itemID, itemtype;

  Boolean bDomainVisible = FALSE;
  Boolean FirstMG = FALSE;
  Boolean AllDomainShow = FALSE;

  LisT l;

  Byte bVisible;

  GetAlignStatus();

     fnAlignList(l);

     if(GetStatus(Cn3D_bDisplayDomain) == TRUE) {     

        SetDomainDataItemStatus(l);        
        pdnmsMaster = GetSelectedModelstruc();
        if (pdnmsMaster != NULL) {
           pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
           pdnmmHead = pmsdMaster->pdnmmHead;

           while(pdnmmHead){
              pmmdThis = pdnmmHead->data.ptrvalue;
              if(pmmdThis){
                 if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte)
AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setoutB1;
                 pmmdThis->bVisible = 1;
                
                 for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next){
                    pmgd = (PMGD) pdnmg->data.ptrvalue;
                    iCountItem = iCountItemGet(pmsdMaster, pmmdThis, pmgd);
                    pmgd->bVisible = GetItemStatus(Cn3D_lFeature, iCountItem);
                    if (pmgd->bReserved && (pmgd->bReserved == pmsdMaster->bAligned) && !Cn3D_fAlignOn) pmgd->bVisible = 0;
                    if ((!(pmgd->bReserved) || (pmgd->bReserved != pmsdMaster->bAligned)) && !Cn3D_fUnalignOn) pmgd->bVisible = 0;
                 }
              }

              SynchronizeModVisibleStatus((PFB) pmgd);

              setoutB1:
              pdnmmHead = pdnmmHead->next;
           }

           pdnmsSlave = pmsdMaster->pdnmsSlaves;
           while(pdnmsSlave) {
              pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;

              pdnmmHead = pmsdSlave->pdnmmHead;
              while(pdnmmHead){
                 pmmdThis = pdnmmHead->data.ptrvalue;
                 if(pmmdThis){
                    if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setoutB2;

                    pmmdThis->bVisible = 1;

                    for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next){
                       pmgd = (PMGD) pdnmg->data.ptrvalue;
                       iCountItem = iCountItemGet(pmsdSlave, pmmdThis, pmgd);
                       pmgd->bVisible = GetItemStatus(Cn3D_lFeature, iCountItem);
                       if (pmgd->bReserved){
                          if ( (*(pmgd->pbMasterReserved) == *(pmsdSlave->pbAligned)) && !Cn3D_fAlignOn) pmgd->bVisible = 0;
                          if ( (*(pmgd->pbMasterReserved) != *(pmsdSlave->pbAligned)) && !Cn3D_fUnalignOn) pmgd->bVisible = 0;
                       }
                       if (!(pmgd->bReserved) && !Cn3D_fUnalignOn) pmgd->bVisible = 0;
                    }
                 }

                 SynchronizeModVisibleStatus((PFB) pmgd);

                 setoutB2:
                 pdnmmHead = pdnmmHead->next;
              }

              pdnmsSlave = pdnmsSlave->next;
           }
        }
     }      
     else if(GetStatus(Cn3D_bDisplayOthers) == TRUE) {
        sfpThis = sfp;
        while(sfpThis){
           sfipThis = sfpThis->data.ptrvalue;
           sfipThis->On = GetItemStatus(Cn3D_lFeature, sfpThis->choice);
           sfpThis = sfpThis->next;
        }

        pdnmsMaster = GetSelectedModelstruc();
        if (pdnmsMaster != NULL) {
           pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
           pdnmmHead = pmsdMaster->pdnmmHead;

           while(pdnmmHead){
              pmmdThis = pdnmmHead->data.ptrvalue;
              if(pmmdThis){
                 if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte)AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setoutB3;

                 for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next){
                    pmgd = (PMGD) pdnmg->data.ptrvalue;
                    pmgd->FeatureOn = GetItemStatus(Cn3D_lFeature, pmgd->iFeature);
                 }
              }
              setoutB3:
              pdnmmHead = pdnmmHead->next;
           }

           pdnmsSlave = pmsdMaster->pdnmsSlaves;
           while(pdnmsSlave) {
              pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;

              pdnmmHead = pmsdSlave->pdnmmHead;
              while(pdnmmHead){
                 pmmdThis = pdnmmHead->data.ptrvalue;
                 if(pmmdThis){
                    if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setoutB4;

                    for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next){
                       pmgd = (PMGD) pdnmg->data.ptrvalue;
                       pmgd->FeatureOn = GetItemStatus(Cn3D_lFeature, pmgd->iFeature);
                    }
                 }

                 setoutB4:
                 pdnmmHead = pdnmmHead->next;
              }
              pdnmsSlave = pdnmsSlave->next;
           }
        }
     }

  Salsa_BioseqUpdate = FALSE;
  Cn3D_ReColor = TRUE;
  Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */
  Cn3D_RedrawProc(FALSE);
  Cn3D_ReColor = FALSE;
  Cn3dObjMgrGetSelected();

  return;

} 
/*---------------------------------------------------------*/
GrouP LIBCALL AlignControls( Nlm_GrouP prnt)
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

  Cn3D_bAlignOn = CheckBox(g, "Aligned", NULL);
  Break(g);
  Cn3D_bUnalignOn = CheckBox(g, "Unaligned", NULL);
  Break(g);
  Cn3D_bDisplayAlignedDomain = CheckBox(g, "Aligned Domain Only", Cn3D_ListAlignedDomainProc);

  return g;
}     
/*---------------------------------------------------------*/
GrouP LIBCALL DisplayControls ( Nlm_GrouP prnt)
{
  GrouP g;
  RecT r1, r2;

  g = NormalGroup ( prnt, 1, 0, "View", systemFont, NULL );
  if (!g) return NULL;
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif

  Cn3D_bDisplayApply = PushButton(g, "Apply", Cn3D_DisplayProc);

  StaticPrompt(g, "Molecule", 0,0, Nlm_systemFont, 'l');
  Cn3D_lStruc = MultiList(g, 10, 4, FillFeatureListProc);    

  StaticPrompt(g, "Alignment", 0,0, Nlm_systemFont, 'l'); 
/*Cn3D_gDisplayAlign = AlignControls(g);   */
  Cn3D_bAlignOn = CheckBox(g, "Aligned", NULL);
  Break(g);
  Cn3D_bUnalignOn = CheckBox(g, "Unaligned", NULL);
  Break(g);
  Cn3D_bDisplayAlignedDomain = CheckBox(g, "Aligned Domain Only", Cn3D_ListAlignedDomainProc);
  SetStatus(Cn3D_bAlignOn, Cn3D_fAlignOn);
  SetStatus(Cn3D_bUnalignOn, Cn3D_fUnalignOn);

  StaticPrompt(g, "Feature", 0,0, Nlm_systemFont, 'l'); 
  Cn3D_bDisplayDomain = CheckBox(g, "Domain", Cn3D_ListDomainProc);      
  Cn3D_bDisplayOthers = CheckBox(g, "Others", Cn3D_ListOthersProc);      
  Cn3D_lFeature = MultiList(g, 10, 8, NULL);
  SetStatus(Cn3D_bDisplayDomain, TRUE);

  ResetDisplayCtrls();     

  return g;
}
/*---------------------------------------------------------*/
void LIBCALL ResetDisplayCtrls(void)
/* set the values of the alignment pane */
{
  ButtoN b;

  Reset(Cn3D_lStruc);      

  Cn3D_ListStrucProc(b);

  if(GetStatus(Cn3D_bDisplayDomain)) {
     Cn3D_ListDomainProc(b);
  }
  else if(GetStatus(Cn3D_bDisplayOthers)) Cn3D_ListOthersProc(b);

  SetStatus(Cn3D_bAlignOn, Cn3D_fAlignOn);
  SetStatus(Cn3D_bUnalignOn, Cn3D_fUnalignOn);
  
  return;
}
/*---------------------------------------------------------*/
void LIBCALLBACK DoLinkSpecialFeatureWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  pmgdThis = (PMGD) pfbThis;
  if(pmgdThis && pmgdThis->bHighlighted == 1){
     pmgdThis->iFeature = iAddedFeature;
     pmgdThis->FeatureOn = 1;
  }

}
/*---------------------------------------------------------*/
void  LinkSpecialFeatureWithMGD(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){

     pmsdMaster = pdnmsMaster->data.ptrvalue;
     TraverseGraphs(pdnmsMaster, 0, 0, NULL, DoLinkSpecialFeatureWithMGD);

      pdnmsSlave = pmsdMaster->pdnmsSlaves;
      while(pdnmsSlave) {
         TraverseGraphs(pdnmsSlave, 0, 0, NULL, DoLinkSpecialFeatureWithMGD);
         pdnmsSlave = pdnmsSlave->next;
      }
  }
  
}
/*---------------------------------------------------------*/
SpecialFeatureInfoPtr LIBCALL SpecialFeatureInfoNew()
{
  SpecialFeatureInfoPtr sfipNew = NULL;
  PRK prkNew = NULL;

  sfipNew = (SpecialFeatureInfoPtr) MemNew(sizeof(SpecialFeatureInfo)); 
  
  if(!sfipNew) return NULL;

  prkNew =  NewRenderKeep();
  sfipNew->prkSpecial = prkNew;

  return sfipNew;
}
/*---------------------------------------------------------*/
void fnCn3D_ModelRedrawWrapper(ButtoN b)
{
  SpecialFeatureInfoPtr sfip = NULL;
  PRK prKeep;
  Char str[1024];
  MsgAnswer ans;

  Int2 j, k;

  if(!parsSpecial) return;

  sfip = SpecialFeatureInfoNew();
  sfip->On = TRUE;

  GetTitle(FeatureTitle, str, sizeof(str));
  sfip->title = StringSave(str);
  
  sfip->parsSpecial = parsSpecial;

  iAddedFeature++;
  ValNodeAddPointer(&sfp, iAddedFeature, (SpecialFeatureInfoPtr) sfip);

  LinkSpecialFeatureWithMGD();
  
  Cn3D_ReColor = TRUE;
  Cn3D_RedrawProc(FALSE);
  Cn3D_ReColor = FALSE;

  parsSpecial = NULL;

  SetTitle(FeatureTitle, "name?");

  ResetModelCtrls();
}
/*---------------------------------------------------------*/
PARS GetSpecialAlgorRenderSet()
{
  if(parsSpecial == NULL){
     parsSpecial = NewAlgorRenderSet();

     parsSpecial->PVirtualBBOn = TRUE;
     parsSpecial->PRealBBOn = FALSE;
     parsSpecial->PExtraBBOn = FALSE;

     parsSpecial->PBBColor = C_red;
/*   parsSpecial->PResiduesOn = TRUE;  */
     parsSpecial->PResiduesOn = FALSE;
     parsSpecial->PResColor = C_red;

     parsSpecial->NTBBColor = C_red;
     parsSpecial->NTResColor = C_red;

     parsSpecial->IonsOn = FALSE;
     parsSpecial->ConnectOn = FALSE;
     parsSpecial->ObjectOn = FALSE; 

     ResetModelCtrls();
  }

  return parsSpecial;

}
/*---------------------------------------------------------*/
void FeatureTitleProc(TexT b)
{
  ResetModelCtrls;
}
/*---------------------------------------------------------*/
void ResetFeatureTitleProc(ButtoN b)
{
  SetTitle(FeatureTitle, "");
}
/*---------------------------------------------------------*/
static void ModelPBBOnOffProc(PopuP p)
{
  Int2 i;
  PARS pars = NULL;

  pars = GetSpecialAlgorRenderSet();
  if(!pars) return;

  pars->PVirtualBBOn = FALSE;
  pars->PRealBBOn = FALSE;
  pars->PExtraBBOn = FALSE;
  i = GetValue(Cn3D_pupModelPBB);

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
  return;
}
/*---------------------------------------------------------*/
static void ModelNTBBOnOffProc(PopuP p)
{
  Int2 i;
  PARS pars = NULL;

  pars = GetSpecialAlgorRenderSet();
  if (!pars) return;


  pars->NTVirtualBBOn = FALSE;
  pars->NTRealBBOn = FALSE;
  pars->NTExtraBBOn = FALSE;
  i = GetValue(Cn3D_pupModelNABB);

  switch (i)
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
  return;
}
/*---------------------------------------------------------*/
static void ModelOnOffProc(LisT l)
{
  Int2 i;
  PARS pars = NULL;

  pars = GetSpecialAlgorRenderSet();
  if (!pars) return;

  pars->PResiduesOn = GetItemStatus(Cn3D_lModelOnOff, 1);
  pars->NTResiduesOn = GetItemStatus(Cn3D_lModelOnOff, 2);

  return;
}
/*---------------------------------------------------------*/
static void Cn3D_SetModelStyle(PopuP p)
{
  Int2 i, j, k;
  PARS pars = NULL;

  pars = GetSpecialAlgorRenderSet();
  if (!pars) return;

  i = GetValue(Cn3D_pupModelStyleItem);
  j = k = 0;
  switch (i)
    {
     case 1: /* prot bb */
         j = pars->PBBRender;
     k = pars->PBBColor;
       break;
     case 2: /* prot sc */
         j = pars->PResRender;
     k = pars->PResColor;
       break;
     case 3: /* na bb */
         j = pars->NTBBRender;
     k = pars->NTBBColor;
       break;
     case 4: /* na sc */
         j = pars->NTResRender;
     k = pars->NTResColor;
       break;
   }
  switch (j)
     {
       case R_SPACE:
         SetValue(Cn3D_pupModelRenderStyle, 5);
        break;
       case R_STICK:
         SetValue(Cn3D_pupModelRenderStyle, 2);
        break;
       case R_BALLNSTICK:
         SetValue(Cn3D_pupModelRenderStyle, 3);
        break;
       case R_THICKWIRE:
         SetValue(Cn3D_pupModelRenderStyle, 4);
        break;
       case R_DEFAULT:
       case R_WIRE:
         SetValue(Cn3D_pupModelRenderStyle, 1);
      }
  switch(k)
     {
       case C_red:
          SetValue(Cn3D_pupModelColorStyle, 1);
          break;
       case C_green:
          SetValue(Cn3D_pupModelColorStyle, 2);
          break;
       case C_magenta:
          SetValue(Cn3D_pupModelColorStyle, 3);
          break;
       case C_sky:
          SetValue(Cn3D_pupModelColorStyle, 4);
          break;
       case C_purple:
          SetValue(Cn3D_pupModelColorStyle, 5);
          break;
     }
     
}
/*---------------------------------------------------------*/
static void Cn3D_SetModelRenderStyle(PopuP p)
{
  Int2 i,  j, k;
  PARS pars = NULL;

  pars = GetSpecialAlgorRenderSet();
  if(!pars) return;

  i = GetValue(Cn3D_pupModelStyleItem);
  j = GetValue(Cn3D_pupModelRenderStyle);
  switch (j)
     {
       case 5:
         k = R_SPACE;
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
       case 1:
         k =  R_WIRE;
      }


  switch (i)
    {
     case 1: /* prot bb */
         pars->PBBRender = k;
         break;
     case 2: /* prot sc */
          pars->PResRender = k;
          break;
     case 3: /* na bb */
         pars->NTBBRender = k;
     break;
     case 4: /* na sc */
         pars->NTResRender = k;
     break;
     case 5: /* hets */
         pars->HetRender = k;
     break;
     case 6: /* ions */
         pars->IonRender = k;
     break;
     case 7: /* connections */
         pars->ConnectRender = k;
     break;
     case 8: /* solvent */
         pars->SolventRender = k;
         break;
     case 9: /* object */
         pars->ObjectRender = k;
     default:
       ;
   }
  return;
}
/*---------------------------------------------------------*/
static void Cn3D_SetModelColorStyle(PopuP p)
{
  Int2 i,  j, k;
  PARS pars = NULL;

  pars = GetSpecialAlgorRenderSet();
  if(!pars) return;
  i = GetValue(Cn3D_pupModelStyleItem);
  j = GetValue(Cn3D_pupModelColorStyle);
   switch (j)
      {
            case 1:
          k = C_red;
          break;
            case 2:
          k = C_green;
          break;
            case 3:
          k = C_magenta;
          break;
        case 4:
           k = C_sky;
          break;
        case 5:
          k = C_purple;
          break;
        default:
          k = 0;
      }
  switch (i)
    {
     case 1: /* prot bb */
         pars->PBBColor = k;
         break;
     case 2: /* prot sc */
          pars->PResColor = k;
          break;
     case 3: /* na bb */
         pars->NTBBColor = k;
     break;
     case 4: /* na sc */
         pars->NTResColor = k;
     break;
   }
  return;
}
/*---------------------------------------------------------*/
GrouP LIBCALL ModelControls ( Nlm_GrouP prnt)
{
  GrouP g;


  g = NormalGroup ( prnt, 1, 0, "Add Feature", systemFont, NULL );
  if (!g) return NULL;
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif

  Cn3D_bModelApply = PushButton(g, "Apply", fnCn3D_ModelRedrawWrapper);

  StaticPrompt(g, "Title ",0,0,Nlm_systemFont,'l');
  FeatureTitle = DialogText(g, "", 10, FeatureTitleProc);
  SetTitle(FeatureTitle, "name?");
/*Advance(g);
  bResetFeatureTitle = DefaultButton(g, "Reset", ResetFeatureTitleProc);     */

/*BACKBONE controls */
  StaticPrompt (g, "Prot. Backbone",0,0,Nlm_systemFont,'l');
  Cn3D_pupModelPBB = PopupList(g, FALSE, ModelPBBOnOffProc);
  PopupItem(Cn3D_pupModelPBB, "alpha C trace");
  PopupItem(Cn3D_pupModelPBB, "partial atom");
  PopupItem(Cn3D_pupModelPBB, "all atoms");
  PopupItem(Cn3D_pupModelPBB, "none");

  StaticPrompt (g, "Nucl. Acid BBone",0,0,Nlm_systemFont,'l');
  Cn3D_pupModelNABB = PopupList(g, FALSE, ModelNTBBOnOffProc);
  PopupItem(Cn3D_pupModelNABB, "P trace");
  PopupItem(Cn3D_pupModelNABB, "partial atom");
  PopupItem(Cn3D_pupModelNABB, "all atoms");
  PopupItem(Cn3D_pupModelNABB, "none");

  StaticPrompt (g, "On/Off", 0,0, Nlm_systemFont, 'l');
  Cn3D_lModelOnOff = MultiList(g ,10,1, ModelOnOffProc);
  ListItem(Cn3D_lModelOnOff, "Prot. Sidechains");
  ListItem(Cn3D_lModelOnOff, "Nucl. Acid Bases");

  StaticPrompt(g, "-------",0,0,Nlm_systemFont,'c');

  StaticPrompt(g, "Style Detail:",0,0,Nlm_systemFont,'l');

  Cn3D_pupModelStyleItem = PopupList(g, FALSE,  Cn3D_SetModelStyle );
  PopupItem(Cn3D_pupModelStyleItem, "Prot. Backbone,");
  PopupItem(Cn3D_pupModelStyleItem, "Prot. Sidechains,");
  PopupItem(Cn3D_pupModelStyleItem, "Nucl. Acid BBone,");
  PopupItem(Cn3D_pupModelStyleItem, "Nucl. Acid Bases,");

  StaticPrompt(g, "Drawn with",0,0,Nlm_systemFont,'l');
  Cn3D_pupModelRenderStyle = PopupList(g , FALSE,  Cn3D_SetModelRenderStyle);
  PopupItem(Cn3D_pupModelRenderStyle, "Wireframe,");
  PopupItem(Cn3D_pupModelRenderStyle, "Tubes,");
  PopupItem(Cn3D_pupModelRenderStyle, "Ball&Stick,");
  PopupItem(Cn3D_pupModelRenderStyle, "Fat Tubes,");
  PopupItem(Cn3D_pupModelRenderStyle, "SpaceFill,");

  StaticPrompt(g, "           ",0,0,Nlm_systemFont,'l');
  StaticPrompt(g, "Colored by ",0,0,Nlm_systemFont,'l');
  Cn3D_pupModelColorStyle = PopupList(g, FALSE, Cn3D_SetModelColorStyle);
  PopupItem(Cn3D_pupModelColorStyle, "Red");
  PopupItem(Cn3D_pupModelColorStyle, "green");
  PopupItem(Cn3D_pupModelColorStyle,  "magenta");
  PopupItem(Cn3D_pupModelColorStyle, "sky");
  PopupItem(Cn3D_pupModelColorStyle, "purple");
  StaticPrompt(g, "           ",0,0,Nlm_systemFont,'l');

  ResetModelCtrls();
  return g;

}
/*---------------------------------------------------------*/
void LIBCALL ResetModelCtrls(void)
{

  PARS pars = NULL;;

  SetValue(Cn3D_pupModelRenderStyle, 1);
  SetValue(Cn3D_pupModelColorStyle, 1);

  pars = GetSpecialAlgorRenderSet();
  if (!pars)
   {
    setout:
      SetItemStatus(Cn3D_lModelOnOff, 1, FALSE);
      SetItemStatus(Cn3D_lModelOnOff, 2, FALSE);

      SetValue(Cn3D_pupModelPBB, 4);
      SetValue(Cn3D_pupModelNABB, 4);
      SetValue(Cn3D_pupModelStyleItem, 1);
      SetValue(Cn3D_pupModelRenderStyle, 1);
      SetValue(Cn3D_pupModelColorStyle, 1);

      return;
   }
  SetItemStatus(Cn3D_lModelOnOff, 1, pars->PResiduesOn);
  SetItemStatus(Cn3D_lModelOnOff, 2, pars->NTResiduesOn);

/* alpha C trace only control */
/*  pdnmlThis = pmsdThis->pdnmlModels;
  pmldThis = (PMLD)pdnmlThis->data.ptrvalue;
  if(pmldThis->iType == Model_type_ncbi_backbone) {
     pars->PVirtualBBOn = TRUE;
     pars->NTVirtualBBOn = TRUE;
     SetValue(Cn3D_pupModelPBB, 1);
     SetValue(Cn3D_pupModelNABB, 1);
  }  */

/*BACKBONE controls */
  if (pars->PVirtualBBOn)  SetValue(Cn3D_pupModelPBB, 1);
  else
  if (pars->PRealBBOn) SetValue(Cn3D_pupModelPBB, 2);
  else
  if (pars->PExtraBBOn) SetValue(Cn3D_pupModelPBB, 3);
  else SetValue(Cn3D_pupModelPBB, 4);

/* set status of BB items from PARS */
  if (pars->NTVirtualBBOn)  SetValue(Cn3D_pupModelNABB, 1);
  else
  if (pars->NTRealBBOn) SetValue(Cn3D_pupModelNABB, 2);
  else
  if (pars->NTExtraBBOn) SetValue(Cn3D_pupModelNABB, 3);
  else SetValue(Cn3D_pupModelNABB, 4);


  SetValue(Cn3D_pupModelStyleItem, 1);
  Cn3D_SetModelStyle(NULL);


 return;
}
