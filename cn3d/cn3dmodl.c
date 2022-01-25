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
* =========================================================================== 
*
* File Name: cn3dmodl.c
*
* Author: Yanli Wang
*
* Version Creation Date: 18/9/1998
*
* $Log: cn3dmodl.c,v $
* Revision 6.3  1998/09/23 22:04:02  ywang
* synchronize show/hide between cn3d and salsa when display complexity is changed
*
 * Revision 6.2  1998/09/23  18:38:48  ywang
 * add functions to control display on domain level
 *
 * Revision 6.1  1998/09/22  18:02:54  ywang
 * panels and functions for display control
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
#include <cn3dmodl.h>

static LisT Cn3D_lStruc;
static GrouP Cn3D_bShow;
static ButtoN Cn3D_DisplayApply;

Int4 iDomainCount = 0;
DomainInfo **domaindata;
extern Boolean Cn3D_ReColor;

void LIBCALL Cn3D_CountDomainProc(void)
{
  Int4 iCount = 0;
  Int4 idom = 0;

  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg;
  PMGD  pmgd;

  pdnmsMaster = GetSelectedModelstruc();
  if (pdnmsMaster != NULL)
  {
    pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
    pdnmmHead = pmsdMaster->pdnmmHead;

    iCount = 0;

    while(pdnmmHead){
       pmmdThis = pdnmmHead->data.ptrvalue;
       if(pmmdThis){
          if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto pre_setout1;

          iCount++;
          idom = 0;

          for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
              pmgd = (PMGD) pdnmg->data.ptrvalue;
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

            iCount++;
            idom = 0;

            for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
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
  pdnmsMaster = GetSelectedModelstruc();
  if (pdnmsMaster != NULL)
  {
    pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
    pdnmmHead = pmsdMaster->pdnmmHead;

    iCount = 0;

    while(pdnmmHead){
       pmmdThis = pdnmmHead->data.ptrvalue;
       if(pmmdThis){
          if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout1;
          StringCpy(domaindata[iCount]->pcPDBName, pmsdMaster->pcPDBName);
          StringCpy(domaindata[iCount]->pcMolName, pmmdThis->pcMolName);
          domaindata[iCount]->iDomain = 0;
          pmmdThis->bVisible = 1;

          iCount++;
          idom = 0;

          for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
              pmgd = (PMGD) pdnmg->data.ptrvalue;
              if(pmgd->iDomain > idom){
                 idom = (Int4) pmgd->iDomain;
                 StringCpy(domaindata[iCount]->pcPDBName, pmsdMaster->pcPDBName);
                 StringCpy(domaindata[iCount]->pcMolName, pmmdThis->pcMolName);
                 domaindata[iCount]->iDomain = pmgd->iDomain;
                    
                 iCount++;
              }
          }
       }
       setout1:
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
            if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout2;
            StringCpy(domaindata[iCount]->pcPDBName, pmsdSlave->pcPDBName);
            StringCpy(domaindata[iCount]->pcMolName, pmmdThis->pcMolName);
            domaindata[iCount]->iDomain = 0;

            iCount++;
            idom = 0;

            for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
               pmgd = (PMGD) pdnmg->data.ptrvalue;
               if(pmgd->iDomain > idom){
                  idom = (Int4) pmgd->iDomain;
                  StringCpy(domaindata[iCount]->pcPDBName, pmsdSlave->pcPDBName);
                  StringCpy(domaindata[iCount]->pcMolName, pmmdThis->pcMolName);
                  domaindata[iCount]->iDomain = pmgd->iDomain;

                  iCount++;
               }
            }
         }
         setout2:
         pdnmmHead = pdnmmHead->next;
      }

      pdnmsSlave = pdnmsSlave->next;
    }
  }

}

Int4 iCountItemGet(PMSD pmsdThis, PMMD pmmdThis, PMGD pmgdThis)
{
  Int4 iCount = 0;
  Int2 iDomain = 0;

  if(pmgdThis == NULL) iDomain = 0;
  else iDomain = pmgdThis->iDomain;

  for(iCount = 0; iCount < iDomainCount; iCount++){
     if(StringCmp(pmsdThis->pcPDBName, domaindata[iCount]->pcPDBName) == 0 && StringCmp(pmmdThis->pcMolName, domaindata[iCount]->pcMolName) ==0 && iDomain == domaindata[iCount]->iDomain) return (iCount + 1) ;
     
  }

  return 0;
}

static void SynchronizeItemStatus(PMSD pmsdThis, PMMD pmmdThis)
{
  Int4 iCount = 0;

  if(pmsdThis == NULL || pmmdThis == NULL) return;

  for(iCount = 0; iCount < iDomainCount; iCount++){
     if(StringCmp(pmsdThis->pcPDBName, domaindata[iCount]->pcPDBName) == 0 && StringCmp(pmmdThis->pcMolName, domaindata[iCount]->pcMolName) ==0) {
        SetItemStatus(Cn3D_lStruc, iCount + 1, TRUE);
     }
  }

  return;
}
static void Cn3D_DisplayProc(ButtoN b )
{
  GrouP g;
  Int2 val;

  Int4 iCount = 0;
  Int4 idom = 0;

  Int4 iCountMedia = 0;
  Byte bVisibleBefore;

  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg;
  PMGD  pmgd;

  SeqIdPtr sipThis;
  Uint2 entityID, itemID, itemtype;

  Boolean bDomainVisible = FALSE;
  Boolean FirstMG = FALSE;
  Boolean AllDomainShow = FALSE;

  g = GetObjectExtra(Cn3D_DisplayApply);

  val = GetValue(g);
  if(val == 0 || val == 2) ResetDisplayCtrls();
  else{
     pdnmsMaster = GetSelectedModelstruc();
     if (pdnmsMaster != NULL) {
        pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
        pdnmmHead = pmsdMaster->pdnmmHead;

        iCountMedia = 0;

        while(pdnmmHead){
           pmmdThis = pdnmmHead->data.ptrvalue;
           bDomainVisible = FALSE;
           AllDomainShow = TRUE;

           if(pmmdThis){
              if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout1;
    
              bVisibleBefore = pmmdThis->bVisible;
 
              iCount = iCountItemGet(pmsdMaster, pmmdThis, NULL);
              pmmdThis->bVisible = GetItemStatus(Cn3D_lStruc, iCount); 

              FirstMG = TRUE;
              for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next){
                 pmgd = (PMGD) pdnmg->data.ptrvalue;
                 if(pmgd->iDomain == 0 || pmmdThis->bVisible == 1) {
                    pmgd->bVisible = pmmdThis->bVisible;
                    if(pmgd->iDomain !=0 && pmgd->bVisible == 1 && FirstMG){
                       SynchronizeItemStatus(pmsdMaster, pmmdThis);
                       FirstMG = FALSE;
                    }
                 }
                 else {
                    iCount = iCountItemGet(pmsdMaster, pmmdThis, pmgd);
                    pmgd->bVisible = GetItemStatus(Cn3D_lStruc, iCount); 
                 }
                 if(pmgd->bVisible == 1) bDomainVisible = TRUE;
                 if(pmgd->bVisible != 1) AllDomainShow = FALSE;
              }

              if(bDomainVisible) pmmdThis->bVisible = 1;
              if(AllDomainShow) {
                 iCount = iCountItemGet(pmsdMaster, pmmdThis, NULL);
                 SetItemStatus(Cn3D_lStruc, iCount, TRUE);
               }

              if(Num_Bioseq != 0 && bVisibleBefore != pmmdThis->bVisible){
                 sipThis = pmmdThis->pSeqId;

                 if(SeqIdForSameBioseq(sipThis, mediadata[iCountMedia]->sip)) {
                     entityID = BioseqFindEntity(sipThis, &itemID);
                     itemtype = OBJ_BIOSEQ;
                     if(pmmdThis->bVisible) ObjMgrSendMsg(OM_MSG_SHOW, entityID, itemID, itemtype);
                     else if(!pmmdThis->bVisible) ObjMgrSendMsg(OM_MSG_HIDE, entityID, itemID, itemtype);
                     mediadata[iCountMedia]->bVisible = pmmdThis->bVisible;
                 }
              }

           }
           setout1:
           pdnmmHead = pdnmmHead->next;
        }

        pdnmsSlave = pmsdMaster->pdnmsSlaves;
        while(pdnmsSlave) {
           pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;
           pdnmmHead = pmsdSlave->pdnmmHead;
 
           iCountMedia++;
           while(pdnmmHead){
              pmmdThis = pdnmmHead->data.ptrvalue;

              bDomainVisible = FALSE;
              AllDomainShow = TRUE;

              if(pmmdThis){
                 if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout2;

                 bVisibleBefore = pmmdThis->bVisible;

                 iCount = iCountItemGet(pmsdSlave, pmmdThis, NULL);
                 pmmdThis->bVisible = GetItemStatus(Cn3D_lStruc, iCount);
                 FirstMG = TRUE;
                 for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next){   
                   pmgd = (PMGD) pdnmg->data.ptrvalue;
                   if(pmgd->iDomain == 0 || pmmdThis->bVisible == 1) {
                      pmgd->bVisible = pmmdThis->bVisible;
                      if(pmgd->iDomain !=0 && pmgd->bVisible == 1 && FirstMG){
                         SynchronizeItemStatus(pmsdSlave, pmmdThis);
                         FirstMG = FALSE;
                      }
                   }
                   else {
                      iCount = iCountItemGet(pmsdSlave, pmmdThis, pmgd);
                      pmgd->bVisible = GetItemStatus(Cn3D_lStruc, iCount);
                   }
                   if(pmgd->bVisible == 1) bDomainVisible = TRUE;
                   if(pmgd->bVisible != 1) AllDomainShow = FALSE;
                 }
                
                 if(bDomainVisible) pmmdThis->bVisible = 1;
                 if(AllDomainShow){
                    iCount = iCountItemGet(pmsdSlave, pmmdThis, NULL);
                    SetItemStatus(Cn3D_lStruc, iCount, TRUE);
                 }

                 if(Num_Bioseq != 0 && bVisibleBefore != pmmdThis->bVisible){
                    sipThis = pmmdThis->pSeqId;

                    if(SeqIdForSameBioseq(sipThis, mediadata[iCountMedia]->sip)) {
                        entityID = BioseqFindEntity(sipThis, &itemID);
                        itemtype = OBJ_BIOSEQ;
                        if(pmmdThis->bVisible) ObjMgrSendMsg(OM_MSG_SHOW, entityID, itemID, itemtype);
                        else if(!pmmdThis->bVisible) ObjMgrSendMsg(OM_MSG_HIDE, entityID, itemID, itemtype);
                        mediadata[iCountMedia]->bVisible = pmmdThis->bVisible;
                    }
                 }

              }
              setout2:
              pdnmmHead = pdnmmHead->next;

           }

           pdnmsSlave = pdnmsSlave->next;
 
        }
     }
  } 


  Salsa_BioseqUpdate = FALSE;
  Cn3D_ReColor = TRUE;
  ResetSalsaColor();
  Cn3D_RedrawProc(FALSE);
  Cn3D_ReColor = FALSE;

  return;

} 

GrouP LIBCALL DisplayControls ( Nlm_GrouP prnt)
{
  GrouP g;

  g = NormalGroup ( prnt, 1, 0, "Display", systemFont, NULL );
  if (!g) return NULL;
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif

  Cn3D_DisplayApply = PushButton(g, "Apply", Cn3D_DisplayProc);
  RadioButton(g, "All / On");
  RadioButton(g, "Only / On");
  SetValue(g, 2);
  SetObjectExtra(Cn3D_DisplayApply, g, NULL);
 
  Cn3D_lStruc = MultiList(g, 10, 15,  NULL);
  

/* set the values for the above */

  ResetDisplayCtrls();     

  return g;
}


void LIBCALL ResetDisplayCtrls(void)
/* set the values of the alignment pane */
{
  Int4 iCount = 0;
  Int4 idom = 0;

  PDNMS pdnmsMaster, pdnmsSlave;
  PMSD pmsdMaster, pmsdSlave;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg;
  PMGD  pmgd;

  Char  MolName[20];
  
  Reset(Cn3D_lStruc);      

  for(iCount = 0; iCount <iDomainCount; iCount++){
     if(domaindata[iCount]->iDomain == 0) sprintf(MolName, "%4s  %2s", domaindata[iCount]->pcPDBName, domaindata[iCount]->pcMolName); 
     else sprintf(MolName, "%4s  %2s  %2d", domaindata[iCount]->pcPDBName, domaindata[iCount]->pcMolName, domaindata[iCount]->iDomain); 
     ListItem(Cn3D_lStruc, (CharPtr) MolName);
     SetItemStatus(Cn3D_lStruc, iCount + 1, TRUE);
  }
     
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

          for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
              pmgd = (PMGD) pdnmg->data.ptrvalue;
              pmgd->bVisible = 1;
          }
       }
       setout1:
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
            if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout2;
            pmmdThis->bVisible = 1;

            for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next) {
               pmgd = (PMGD) pdnmg->data.ptrvalue;
               pmgd->bVisible = 1;

            }
         }
         setout2:
         pdnmmHead = pdnmmHead->next;
      }

      pdnmsSlave = pdnmsSlave->next;
    }
  }

}
