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
* Revision 6.63  1999/09/10 13:49:38  lewisg
* eliminate spaces before preprocessor directives
*
* Revision 6.62  1999/07/27 14:19:09  ywang
* fix previous in another way so that not to break any existant code
*
* Revision 6.61  1999/07/26 21:42:47  vakatov
* Fixed a multiple ValNode list deallocation in
* ClearSpecialFeature/SpecialFeatureFree
*
* Revision 6.60  1999/07/15 15:59:41  ywang
* improve feature control, call Cn3D_DeHighlightAll after ObjMgrDeSelectAll
*
* Revision 6.59  1999/07/12 14:15:37  ywang
* initialize ButtoN b for internal call Cn3D_ListDomainProc & Cn3D_ListStrucProc
*
* Revision 6.58  1999/07/07 20:45:35  ywang
* clear domaindata, mediadata, special feature before reading in new data in cn3d
*
* Revision 6.57  1999/07/01 16:01:08  ywang
* improve display control panel
*
* Revision 6.56  1999/06/16 15:13:29  ywang
* for feature editing: direct return for 'Cancel' 'Delete' 'Edit' 'Save' if no residues or no feature name selected
*
* Revision 6.55  1999/06/07 21:25:14  ywang
* add iUserDefinedFeatureOld, FeatureOnOld to protect saved features being overwritten by temporary feature
*
* Revision 6.54  1999/06/07 20:26:11  ywang
* change prompt for featue editing
*
* Revision 6.53  1999/06/07 19:39:50  ywang
* remove obsolete user defined features, an obsolete feature means a feature whose region has been completely overwirtten by a later defined feature
*
* Revision 6.52  1999/06/03 20:29:24  ywang
* fix bug in aligned domain display when they are from same structure
*
* Revision 6.51  1999/05/27 16:11:22  ywang
* initilize all local variables at defined
*
* Revision 6.50  1999/05/03 16:54:23  ywang
* fix bug for recording and matching mol-id for user defined features
*
* Revision 6.49  1999/04/22 21:10:44  ywang
* check NULL pointer on domain counting
*
* Revision 6.48  1999/04/21 21:15:09  ywang
* replace Cn3D_DeHighlightAll with ObjMgrDeSelectAll
*
* Revision 6.46  1999/04/16 22:21:23  ywang
* update residue aligned status on mediadata upon sequence SHOW/HIDE
*
* Revision 6.45  1999/04/15 21:37:35  ywang
* add AssignDomainAlignedStatusForStrucSeqs
*
* Revision 6.44  1999/04/14 17:43:26  ywang
* synchronize object show/hide upon action 'Show Selected Only' action
*
* Revision 6.43  1999/04/14 17:19:15  ywang
* synchronize object show/hide upon action 'Show Selected Only' action
*
* Revision 6.42  1999/04/13 23:06:14  ywang
* fix bug for coloring strucseqs data
*
* Revision 6.41  1999/04/06 20:11:02  lewisg
* more opengl
*
* Revision 6.40  1999/04/02 20:51:57  ywang
* fix warnings
*
* Revision 6.39  1999/03/22 22:41:52  ywang
* redesign feature page, fix bugs
*
* Revision 6.38  1999/03/18 23:13:03  ywang
* modify feature window
*
* Revision 6.37  1999/03/18 22:28:56  ywang
* add functions for saveout+readin+index user defined features
*
* Revision 6.36  1999/03/08 21:16:21  ywang
* initialize variables
*
* Revision 6.35  1999/03/08 19:34:35  ywang
* redesign Feature page
*
* Revision 6.34  1999/03/03 23:17:21  lewisg
* one master struct at a time, list slaves in structure info, bug fixes
*
* Revision 6.33  1999/03/01 20:21:30  ywang
* put in options for residue on/off control
*
* Revision 6.32  1999/02/24 23:00:02  ywang
* minor name change for user defined feature at MGD nodecn3dmodl.c
*
* Revision 6.31  1999/02/12 15:42:56  ywang
* change void static to static void for ClearSpecialFeaturewithMGD(void)
*
* Revision 6.30  1999/02/11 22:38:51  ywang
* fix bug on display highlight residues only--if no res are highlighted, cn3d sets that button status as FALSE and draw whole structurescn3dwin.c
*
* Revision 6.29  1999/02/08 18:42:25  ywang
* improve feature editing
*
* Revision 6.28  1999/02/04 16:14:44  ywang
* support delete added features
*
* Revision 6.27  1999/02/03 23:09:24  ywang
* add functions to allow editing 'added feature' in Model menu
*
* Revision 6.26  1999/02/02 22:25:27  ywang
* improve feature edit function
*
* Revision 6.25  1999/02/01 21:10:09  ywang
* correct mistake
*
* Revision 6.24  1999/02/01 20:43:25  ywang
* improve 'Model' menu
*
* Revision 6.23  1999/01/27 21:51:51  ywang
* add label to 'Model' menu
*
* Revision 6.22  1999/01/27 16:21:55  ywang
* redesign 'Model' menu
*
* Revision 6.21  1999/01/26 17:14:35  ywang
* redesign Display menu and add 'display highlight residues only' function
*
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
#ifdef _OPENGL
#include <shim3d.h>
#else
#include <viewer3d.h>
#endif
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
#include <objmime.h>

static ENUM_ALIST(empty_alist)
  {"              ",     1},
END_ENUM_ALIST

static ENUM_ALIST(all_render_alist)
  {"Wireframe",          1},
  {"Tubes",              2},
  {"Ball & Stick",       3},
  {"Fat Tubes",          4},
  {"Space Fill   ",      5},
END_ENUM_ALIST

static ENUM_ALIST(all_color_alist)
  {"red",           1},
  {"green",         2},
  {"magenta",       3},
  {"sky",           4},
  {"purple",        5},
END_ENUM_ALIST

static EnumFieldAssocPtr renderAlist [4] =
  {all_render_alist, all_render_alist, all_render_alist, all_render_alist};

static EnumFieldAssocPtr colorAlist [4] =
  {all_color_alist, all_color_alist, all_color_alist, all_color_alist};

static ButtoN   Cn3D_lModelOnOffItem[10];        /* pieces parts to draw */
static PopuP    Cn3D_pupModelPBB;        /* protein backbone options */
static PopuP    Cn3D_pupModelNABB;       /* nucl. acid bb options */
static PopuP    Cn3D_pupModelStyleItem;
static PopuP    Cn3D_pupModelRenderStyle[10];
static PopuP    Cn3D_pupModelColorStyle[10];
static TexT     FeatureTitle; 
static TexT     FeatureDescription; 
static ButtoN   bResetFeatureTitle; 

static ButtoN   Cn3D_lOnOffLabel [4];
static PopuP    Cn3D_pupLabelAA;
static PopuP    Cn3D_pupLabelNT;
static PopuP    Cn3D_bLName [4];
static ButtoN   Cn3D_bLNum [4];
static ButtoN   Cn3D_bLTop [4];
static ButtoN   Cn3D_bLWhite [4];
static PopuP    Cn3D_pupLabelSize [4];

static ButtoN Cn3D_bDisplayApply;
static ButtoN Cn3D_bModelApply;

static ButtoN Cn3D_bDisplayHighlight;
static ButtoN Cn3D_bHideHighlight;
static ButtoN Cn3D_bResidueDisplayReset;

static ButtoN Cn3D_bDisplayByStruc;
static LisT Cn3D_lStruc;

static GrouP Cn3D_gDisplayAlign;
static ButtoN  Cn3D_bAlignOn, Cn3D_bUnalignOn;
static ButtoN Cn3D_bDisplayAlignedDomain;

static ButtoN Cn3D_bDisplaySS;
static ButtoN Cn3D_bDisplayOthers;
static LisT Cn3D_lFeature;

static ButtoN bFeatureAdd;
static ButtoN bFeatureCancel;
static ButtoN bFeatureDelete;
static ButtoN bFeatureEdit;
static LisT Cn3D_lFeature2;
static LisT Cn3D_lFeature3;

Int4 iDomainCount = 0;
DomainInfo **domaindata;
extern Boolean Cn3D_ReColor;

SpecialFeaturePtr sfp = NULL;
SpecialFeaturePtr sfpThis_head = NULL, sfpThis_tail = NULL;
PARS parsSpecial = NULL;
Byte iAddedFeature = 0, iEditedFeature = 0;

Boolean Cn3D_DisplayHighlight = FALSE;
Boolean Cn3D_NoSingleHL = FALSE;
Boolean HLStatus = FALSE;
Boolean JustHLStatus = FALSE;
Boolean FeatureOverwritten = FALSE;
/*-----------------------------------------------*/
void LIBCALLBACK Cn3DCheckNoSingleHighlight(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PMGD pmgdThis = NULL;
  PDNMG  pdnmgThis = NULL;

  pmgdThis = (PMGD) pfbThis;
  if(pmgdThis->bHighlighted == 1) Cn3D_NoSingleHL = FALSE;
}
/*-----------------------------------------------*/
void Cn3DCheckHighlighted(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

       /* replace old function which depends on ObjMgr */
       /* by doing so, highlight for non ObjMgr registered residues will */
       /* also be picked up when do Cn3D_Redraw */

  Cn3D_NoSingleHL = TRUE;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){
      pmsdMaster = pdnmsMaster->data.ptrvalue;

      TraverseGraphs(pdnmsMaster, 0, 0, NULL, Cn3DCheckNoSingleHighlight);
      pdnmsSlave = pmsdMaster->pdnmsSlaves;
      while(pdnmsSlave) {
         TraverseGraphs(pdnmsSlave, 0, 0, NULL, Cn3DCheckNoSingleHighlight);
         pdnmsSlave = pdnmsSlave->next;
      }
   }
}
/*---------------------------------------------------------*/
void LIBCALLBACK Cn3DCheckMGTurnOnOffStatus(PFB pfbThis,Int4 iModel, Int4 iIndex,
Pointer ptr)
{
   PMGD pmgdThis = NULL;
   PDNMG  pdnmgThis = NULL;

   Int4 iCount = 0;

   pmgdThis = (PMGD) pfbThis;
  
   if(GetStatus(Cn3D_bHideHighlight)){
      if(pmgdThis->bHighlighted == 1) {
         pmgdThis->bVisible = 0;
         pmgdThis->bTurnedOff = 1;
      }
/*    else if(pmgdThis->bTurnedOff == 1){
         pmgdThis->bVisible = 0;
      } */ /* control in RenderGraph */
   }
   else if(GetStatus(Cn3D_bResidueDisplayReset)){
      for(iCount = 0; iCount < iDomainCount; iCount++){
         if(pmgdThis->iDomain == domaindata[iCount]->iDomain) {
            pmgdThis->bVisible = domaindata[iCount]->bVisible;
            pmgdThis->bTurnedOff = 0;
            break;
         }
      }
   }
}
/*---------------------------------------------------------*/
void Cn3DCheckTurnOnOffStatus(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){
      pmsdMaster = pdnmsMaster->data.ptrvalue;

      TraverseGraphs(pdnmsMaster, 0, 0, NULL, Cn3DCheckMGTurnOnOffStatus);
      pdnmsSlave = pmsdMaster->pdnmsSlaves;
      while(pdnmsSlave) {
         TraverseGraphs(pdnmsSlave, 0, 0, NULL, Cn3DCheckMGTurnOnOffStatus);
         pdnmsSlave = pdnmsSlave->next;
      }
   }
}
/*---------------------------------------------------------*/
void GetAlignStatus(void)
{
  Cn3D_fAlignOn = GetStatus(Cn3D_bAlignOn);
  Cn3D_fUnalignOn = GetStatus(Cn3D_bUnalignOn);
}
/*---------------------------------------------------------*/
static void Cn3D_ListSSProc(ButtoN b)
{
  SetStatus(Cn3D_bDisplayOthers, FALSE);

  ErrPostEx(SEV_ERROR, 0, 0, "Not function yet!  Worth?  Please comment, Thank you!"); 
}
/*---------------------------------------------------------*/
static void Cn3D_ListOthersProc(void)
{

  SpecialFeaturePtr sfpThis = NULL;
  SpecialFeatureInfoPtr sfipThis = NULL;
  Int4 iCount = 1;

/*SetStatus(Cn3D_bDisplayOthers, FALSE);
  Message(MSG_OK, "This is assumed to display a list of user specified features! ");  */

  Reset(Cn3D_lFeature3);

  sfpThis = sfp;
  while(sfpThis){
     sfipThis = sfpThis->data.ptrvalue;
     ListItem(Cn3D_lFeature3, (CharPtr) sfipThis->title);
     SetItemStatus(Cn3D_lFeature3, iCount, sfipThis->On);
     iCount++;
    
     sfpThis = sfpThis->next;
  }

}
/*---------------------------------------------------------*/
static void Cn3D_SelectDomainProc(LisT l)
{
  SetStatus(Cn3D_bDisplayAlignedDomain, FALSE);
}
/*---------------------------------------------------------*/
static void Cn3D_ListDomainProc(ButtoN b)
{
  Char  MolName[20];
  Int4 iCount = 0, iCountActive = 0;

  SetStatus(Cn3D_bDisplayOthers, FALSE);
  SetStatus(Cn3D_bDisplaySS, FALSE);
  Reset(Cn3D_lFeature);

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
  ButtoN b1 = NULL;

  Cn3D_ListDomainProc(b1);
}
/*---------------------------------------------------------*/
static void FillFeatureListProc(LisT l)
{
  ButtoN b = NULL;

  Cn3D_ListDomainProc(b);
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
  PDNMG pdnmg = NULL;
  PMGD  pmgd = NULL;

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
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg = NULL;
  PMGD  pmgd = NULL;

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
   PDNMS pdnmsMaster = NULL;

   Int4 iCount = 0;

   pdnmsMaster = GetSelectedModelstruc();
   if (pdnmsMaster != NULL)
       {

	   for(iCount = 0; iCount < Num_Bioseq; iCount++){
	       if(!mediadata[iCount]->bVisible){
		   itemtype = OBJ_BIOSEQ;
		   ObjMgrSendMsg(OM_MSG_SHOW, mediadata[iCount]->entityID,
				 mediadata[iCount]->itemID, itemtype); 
		   mediadata[iCount]->bVisible = 1;
	       }
	   }
       }
}
/*---------------------------------------------------------*/
static void Cn3D_ListStrucProc(ButtoN b)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

  Int4 iCount = 0;

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
void AssignDomainAlignedStatusForStrucSeqs(void)
{
  PDNMS pdnmsMaster = NULL;
  PMSD pmsdMaster = NULL;
  PDNMM pdnmmHead = NULL, pdnmmHead_temp = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg = NULL;
  PMGD  pmgd = NULL;

  Int4 iCount = 0, iCountChain = 0;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster == NULL) return;

  pmsdMaster = (PMSD)pdnmsMaster->data.ptrvalue;
  pdnmmHead = pmsdMaster->pdnmmHead;  

  for(iCount = 0; iCount < iDomainCount; iCount++){
     domaindata[iCount]->bAligned = FALSE;
     pdnmmHead_temp = pdnmmHead;
    
     iCountChain = 0;

     while(pdnmmHead_temp){
        pmmdThis = pdnmmHead_temp->data.ptrvalue;
        if(pmmdThis == NULL) break;
        if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout1;

        if(domaindata[iCount]->iChainIndex == iCountChain){
           for(pdnmg = pmmdThis->pdnmgHead; pdnmg != NULL; pdnmg = pdnmg->next){
              pmgd = (PMGD) pdnmg->data.ptrvalue;
              if(pmgd->iDomain == domaindata[iCount]->iDomain && pmgd->bReserved){
                 domaindata[iCount]->bAligned = TRUE;
                 goto setout2;
              }
           }
        }

        setout1:
        pdnmmHead_temp = pdnmmHead_temp->next;
        iCountChain++;
    }
    setout2: 
    continue;
  }
     
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

  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg = NULL;
  PMGD  pmgd = NULL;

  Boolean EncounterZeroDomain = FALSE;

  pdnmsMaster = GetSelectedModelstruc();

  if(pdnmsMaster == NULL) return;

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

  if(pmsdMaster->iMimeType == NcbiMimeAsn1_strucseqs) {
     AssignDomainAlignedStatusForStrucSeqs();
  }
  else AssignDomainAlignedStatus();
}
/*---------------------------------------------------------*/
void SynchronizeModVisibleStatus(PFB pfbThis)
{
   PMSD pmsdThis = NULL;
   PMMD pmmdThis = NULL;
   PDNMG pdnmg = NULL;
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
          pmodThis->bVisible = 0;
          pvnThis = pmodThis->pvnContains;
          while(pvnThis){
             pmgdThis = pvnThis->data.ptrvalue;
             if(Cn3D_DisplayHighlight) pmodThis->bVisible = pmgdThis->bHighlighted;
             else pmodThis->bVisible = pmgdThis->bVisible;
             if (pmodThis->bVisible == 1) break;
             pvnThis = pvnThis->next;
          }

          pvnmoThis = pvnmoThis->next;
       }
   }
}
/*---------------------------------------------------------*/
static void SetDomainDataItemStatus(LisT l)
{
  Int4 iCount = 0, iCountActive = 0;

     for(iCount = 0; iCount < iDomainCount; iCount++){
        if(domaindata[iCount]->bVisibleParent) {
           domaindata[iCount]->bVisible = GetItemStatus(Cn3D_lFeature, iCountActive + 1);
           iCountActive++;
        }
     }

  return; 
}
/*---------------------------------------------------------*/
Int4 iCountItemGet(PMSD pmsdThis, PMMD pmmdThis, PMGD pmgdThis, Int4 iCountStruc)
{
  Int4 iCount = 0, iCountLive = 0;
  Int2 iDomain = 0;

  if(pmgdThis == NULL) iDomain = 0;
  else iDomain = pmgdThis->iDomain;

  for(iCount = 0; iCount < iDomainCount; iCount++){
     if(!domaindata[iCount]->bVisibleParent) continue;
     if(StringCmp(pmsdThis->pcPDBName, domaindata[iCount]->pcPDBName) == 0 && StringCmp(pmmdThis->pcMolName, domaindata[iCount]->pcMolName) ==0 && iDomain ==
domaindata[iCount]->iDomain && iCountStruc == domaindata[iCount]->iStrucIndex) return (iCountLive + 1) ;
     iCountLive++;
    
  }

  return 0;
}
/*---------------------------------------------------------*/
static void fnAlignList(LisT l)
/* set the values of the alignment pane */
{
  Int4 iCount = 0;
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;
  BiostrucFeaturePtr pbsfThis = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  Byte bVisible = 0;
  SeqIdPtr sipThis = NULL;
  Uint2 entityID = 0, itemID = 0, itemtype = 0;         

#ifndef _OPENGL
  Cn3D_SaveActiveCam();
#endif

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

    if(pmsdMaster->iMimeType != NcbiMimeAsn1_strucseqs){
       TraverseGraphs( pdnmsMaster, 0, 0, NULL, fnClearMarkedResidues);
       pmsdMaster->bAligned = 0;
    }
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

  Cn3DCheckAlignmentStatusForStrucSeqsForMasterSeq();
  Cn3DCheckAlignmentStatusForStrucSeqs();
#ifndef _OPENGL
  Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */
#endif
}
/*---------------------------------------------------------*/
void Cn3D_bDisplayHighlightStatusSet(Boolean Yes)
{
  SetStatus(Cn3D_bDisplayHighlight, Yes);
}
/*---------------------------------------------------------*/
static void Cn3D_DisplayHighlightProc(ButtoN b )
{
  Cn3D_DisplayHighlight = GetStatus(Cn3D_bDisplayHighlight);
  if(GetStatus(Cn3D_bDisplayHighlight)){
     SetStatus(Cn3D_bHideHighlight, FALSE);
     SetStatus(Cn3D_bResidueDisplayReset, FALSE);
  }
}
/*---------------------------------------------------------*/
static void Cn3D_HideHighlightProc(ButtoN b)
{
  if(GetStatus(Cn3D_bHideHighlight)){
     SetStatus(Cn3D_bDisplayHighlight, FALSE);
     SetStatus(Cn3D_bResidueDisplayReset, FALSE);
  }
}
/*---------------------------------------------------------*/
static void Cn3D_ResidueDisplayResetProc(ButtoN b)
{
  if(GetStatus(Cn3D_bResidueDisplayReset)){
    SetStatus(Cn3D_bDisplayHighlight, FALSE);
    SetStatus(Cn3D_bHideHighlight, FALSE);
  }

}
/*---------------------------------------------------------*/
static void Cn3D_DisplayProc(ButtoN b )
{
  GrouP g;
  Int2 val;

  Int4 iCountItem = 0, iCountStruc = 0;
  Int4 idom = 0;

  SpecialFeaturePtr sfpThis = NULL;
  SpecialFeatureInfoPtr sfipThis = NULL;

  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg = NULL;
  PMGD  pmgd = NULL;
  PVNMO pvnmoThis = NULL;

  SeqIdPtr sipThis = NULL;
  Uint2 entityID = 0, itemID = 0, itemtype = 0;

  Boolean bDomainVisible = FALSE;
  Boolean FirstMG = FALSE;
  Boolean AllDomainShow = FALSE;

  LisT l = NULL;

  Byte bVisible = 0;

  GetAlignStatus();

  fnAlignList(l);

  if(Cn3D_DisplayHighlight) Cn3DCheckHighlighted();
  if(Cn3D_NoSingleHL) {
     Cn3D_DisplayHighlight = FALSE;
     Cn3D_bDisplayHighlightStatusSet(FALSE);
  }

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
                    iCountItem = iCountItemGet(pmsdMaster, pmmdThis, pmgd, iCountStruc);
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
           iCountStruc++;
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
                       iCountItem = iCountItemGet(pmsdSlave, pmmdThis, pmgd, iCountStruc);
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
              iCountStruc++;
           }
        }

  if(GetStatus(Cn3D_bHideHighlight) || GetStatus(Cn3D_bResidueDisplayReset)) Cn3DCheckTurnOnOffStatus();

  Salsa_BioseqUpdate = FALSE;
  Cn3D_ReColor = TRUE;
#ifndef _OPENGL
  Cn3D_v3d->is_zoomed = TRUE;  /* keep the proteins from moving */
#endif
  Cn3D_RedrawProc(FALSE);
  Cn3D_ReColor = FALSE;
  Cn3dObjMgrGetSelected();

  return;

} 
/*---------------------------------------------------------*/
static void Cn3D_bAlignOnProc(ButtoN b)
{
  SetStatus(Cn3D_bResidueDisplayReset, FALSE);
}
/*---------------------------------------------------------*/
static void Cn3D_bUnalignOnProc(ButtoN b)
{
  SetStatus(Cn3D_bResidueDisplayReset, FALSE);
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

  Cn3D_bAlignOn = CheckBox(g, "Aligned", Cn3D_bAlignOnProc);
  Break(g);
  Cn3D_bUnalignOn = CheckBox(g, "Unaligned", Cn3D_bUnalignOnProc);
  Break(g);
  Cn3D_bDisplayAlignedDomain = CheckBox(g, "Aligned Domain Only", Cn3D_ListAlignedDomainProc);

  return g;
}     
/*---------------------------------------------------------*/
GrouP LIBCALL DisplayControls ( Nlm_GrouP prnt)
{
  GrouP g, h;
  RecT r1, r2;
  PrompT ppt1, ppt2, ppt3, ppt4;

  g = HiddenGroup ( prnt, -1, 0, NULL );
  if (!g) return NULL;
  SetGroupSpacing ( g, 3, 9 );
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif

  StaticPrompt(g, "", 0, stdLineHeight, systemFont, 'l');    
  StaticPrompt(g, "", 0, stdLineHeight, systemFont, 'l');       

  h = HiddenGroup (g, 0, -4, NULL);
  SetGroupSpacing (h, 30, 10);

/*Cn3D_bDisplayApply = PushButton(g, "Apply", Cn3D_DisplayProc); */

  ppt1 = StaticPrompt(h, "Molecule", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_lStruc = MultiList(h, 6, 10, FillFeatureListProc);    
  StaticPrompt(h, "", 0, stdLineHeight, systemFont, 'l');    
  StaticPrompt(h, "", 0, stdLineHeight, systemFont, 'l');       
  AlignObjects (ALIGN_LEFT, (HANDLE) ppt1, (HANDLE) Cn3D_lStruc, NULL);      
  
  ppt3 = StaticPrompt(h, "Domain", 0, stdLineHeight + 5, systemFont, 'c');
/*Cn3D_bDisplaySS = CheckBox(h, "Secondary Structure", Cn3D_ListSSProc);
  Cn3D_bDisplayOthers = CheckBox(h, "User Defined", Cn3D_ListOthersProc); */
  Cn3D_bDisplayAlignedDomain = CheckBox(h, "Aligned Domain Only", Cn3D_ListAlignedDomainProc);    
  Cn3D_lFeature = MultiList(h, 6, 8, Cn3D_SelectDomainProc);
  StaticPrompt(h, "", 0, stdLineHeight, systemFont, 'l');       
  AlignObjects (ALIGN_LEFT, (HANDLE) ppt3, (HANDLE) Cn3D_lFeature, NULL);

  ppt4 = StaticPrompt(h, "Residue", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_bAlignOn = CheckBox(h, "Aligned", Cn3D_bAlignOnProc);   
  Cn3D_bUnalignOn = CheckBox(h, "Unaligned", Cn3D_bUnalignOnProc);   
  Cn3D_bDisplayHighlight =  CheckBox(h, "Show Selected Only", Cn3D_DisplayHighlightProc);
/*Cn3D_bHideHighlight = CheckBox(h, "Hide Selected Only", Cn3D_HideHighlightProc);
  Cn3D_bResidueDisplayReset = CheckBox(h, "Reset", Cn3D_ResidueDisplayResetProc); 
  AlignObjects (ALIGN_LEFT, (HANDLE) ppt4, (HANDLE) Cn3D_bDisplayHighlight, (HANDLE) Cn3D_bHideHighlight, (HANDLE) Cn3D_bResidueDisplayReset, NULL);  */
  AlignObjects (ALIGN_LEFT, (HANDLE) ppt4, (HANDLE) Cn3D_bDisplayHighlight, (HANDLE) Cn3D_bHideHighlight, (HANDLE) Cn3D_bResidueDisplayReset, NULL);
  SetStatus(Cn3D_bDisplayHighlight, FALSE); 
  SetStatus(Cn3D_bAlignOn, Cn3D_fAlignOn);
  SetStatus(Cn3D_bUnalignOn, Cn3D_fUnalignOn);
/*SetStatus(Cn3D_bHideHighlight, FALSE); 
  SetStatus(Cn3D_bResidueDisplayReset, FALSE);  */

/*ppt2 = StaticPrompt(h, "Alignment", 0, stdLineHeight + 5, systemFont, 'c');   
  Cn3D_bAlignOn = CheckBox(h, "Aligned", Cn3D_bAlignOnProc);   
  Cn3D_bUnalignOn = CheckBox(h, "Unaligned", Cn3D_bUnalignOnProc);   
  Cn3D_bDisplayAlignedDomain = CheckBox(h, "Aligned Domain Only", Cn3D_ListAlignedDomainProc);    
  StaticPrompt(h, "", 0, stdLineHeight, systemFont, 'l');    
  AlignObjects (ALIGN_LEFT, (HANDLE) ppt2, (HANDLE) Cn3D_bAlignOn, (HANDLE) Cn3D_bUnalignOn, (HANDLE) Cn3D_bDisplayAlignedDomain, NULL);  */
  SetStatus(Cn3D_bAlignOn, Cn3D_fAlignOn);
  SetStatus(Cn3D_bUnalignOn, Cn3D_fUnalignOn);

  StaticPrompt(g, "", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_bDisplayApply = PushButton(g, "Apply!", Cn3D_DisplayProc);
  AlignObjects (ALIGN_CENTER, (HANDLE)Cn3D_bDisplayApply, (HANDLE) h, NULL);

  ResetDisplayCtrls();     

  return g;
}
/*---------------------------------------------------------*/
void LIBCALL ResetDisplayCtrls(void)
/* set the values of the alignment pane */
{
  ButtoN b = NULL;

  Reset(Cn3D_lStruc);      

  Cn3D_ListStrucProc(b);
  Cn3D_ListDomainProc(b);

  SetStatus(Cn3D_bAlignOn, Cn3D_fAlignOn);
  SetStatus(Cn3D_bUnalignOn, Cn3D_fUnalignOn);
  
  return;
}
/*---------------------------------------------------------*/
static Int2 GetRenderItemIndex(Int2 k)
{
  Int2 i = 1;

  switch(k){
    case R_WIRE:
       i = 1; 
       break;
    case R_STICK:
       i = 2;
       break;
    case R_BALLNSTICK:
       i = 3;
       break;
    case R_THICKWIRE:
       i = 4;
       break;
    case R_SPACE:
       i = 5;
       break;
    default:
       ;
  }

  return i;

}
/*---------------------------------------------------------*/
static Int2 GetColorItemIndex(Int2 k)
{
  Int2 i = 1;

  switch(k){
     case C_red:
        i = 1;
        break;
     case  C_green:
        i = 2;
        break;
     case C_magenta:
        i = 3;
        break;
     case C_sky:
        i = 4;
        break;
     case C_purple:
        i = 5;
        break;
     default:
        ;
  }

  return i;

}
/*---------------------------------------------------------*/
void UpdateFeatureList2()
{
  SpecialFeaturePtr sfpThis = NULL;
  SpecialFeatureInfoPtr sfipThis = NULL;

  Reset(Cn3D_lFeature2);
  sfpThis = sfp;
  while(sfpThis){
     sfipThis = sfpThis->data.ptrvalue;
     ListItem(Cn3D_lFeature2, (CharPtr) sfipThis->title);   
     sfpThis = sfpThis->next;
  }
  SetValue(Cn3D_lFeature2, 0);

}
/*---------------------------------------------------------*/
void  LIBCALLBACK DoCleanSpecialFeatureWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  pmgdThis = (PMGD) pfbThis;
  if(pmgdThis && pmgdThis->iUserDefinedFeature == iAddedFeature){
     pmgdThis->iUserDefinedFeature = pmgdThis->iUserDefinedFeatureOld;
     pmgdThis->FeatureOn = pmgdThis->iUserDefinedFeatureOld;

     pmgdThis->iUserDefinedFeatureOld = 0;
     pmgdThis->iUserDefinedFeatureOld = 0;
  }

}
/*---------------------------------------------------------*/
static void ClearSpecialFeaturewithMGD(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){

     pmsdMaster = pdnmsMaster->data.ptrvalue;
     TraverseGraphs(pdnmsMaster, 0, 0, NULL, DoCleanSpecialFeatureWithMGD);

      pdnmsSlave = pmsdMaster->pdnmsSlaves;
      while(pdnmsSlave) {
         TraverseGraphs(pdnmsSlave, 0, 0, NULL, DoCleanSpecialFeatureWithMGD);
         pdnmsSlave = pdnmsSlave->next;
      }
  }

}
/*---------------------------------------------------------*/
void Cn3D_CheckHLStatus(void)
{
  PDNMS pdnmsThis = NULL, pdnmsSlave = NULL;
  PMSD  pmsdMaster = NULL, pmsdSlave = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  HLStatus = FALSE;

  pdnmsThis = GetSelectedModelstruc();
  if(pdnmsThis) pmsdMaster = (PMSD) pdnmsThis->data.ptrvalue;
  else return;

  pdnmmHead = pmsdMaster->pdnmmHead;
  while(pdnmmHead){
     pmmdThis = pdnmmHead->data.ptrvalue;
     if(pmmdThis){
        pdnmgThis = pmmdThis->pdnmgHead;
        while(pdnmgThis){
           pmgdThis = pdnmgThis->data.ptrvalue;
           if(pmgdThis && pmgdThis->bHighlighted == 1){
              HLStatus = TRUE;
              return;  
           }
           pdnmgThis = pdnmgThis->next;
        }
     }
     pdnmmHead = pdnmmHead->next;
  }

  pdnmsSlave = pmsdMaster->pdnmsSlaves;
  while(pdnmsSlave){
     pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;
     pdnmmHead = pmsdSlave->pdnmmHead;
      while(pdnmmHead){
         pmmdThis = pdnmmHead->data.ptrvalue;
         if(pmmdThis){
            pdnmgThis = pmmdThis->pdnmgHead;
            while(pdnmgThis){
               pmgdThis = pdnmgThis->data.ptrvalue;
               if(pmgdThis && pmgdThis->bHighlighted == 1){
                  HLStatus = TRUE;
                  return;
               }
               pdnmgThis = pdnmgThis->next;
            }
         }
         pdnmmHead = pdnmmHead->next;
     }
     pdnmsSlave = pdnmsSlave->next;
  }      

}
/*---------------------------------------------------------*/
void Cn3D_CheckJustHLStatus(void)
{
  PDNMS pdnmsThis = NULL, pdnmsSlave = NULL;
  PMSD  pmsdMaster = NULL, pmsdSlave = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  JustHLStatus = FALSE;

  pdnmsThis = GetSelectedModelstruc();
  if(pdnmsThis) pmsdMaster = (PMSD) pdnmsThis->data.ptrvalue;
  else return;

  pdnmmHead = pmsdMaster->pdnmmHead;
  while(pdnmmHead){
     pmmdThis = pdnmmHead->data.ptrvalue;
     if(pmmdThis){
        pdnmgThis = pmmdThis->pdnmgHead;
        while(pdnmgThis){
           pmgdThis = pdnmgThis->data.ptrvalue;
           if(pmgdThis && pmgdThis->bJustHighlighted == 1){
              JustHLStatus = TRUE;
              return;
           }
           pdnmgThis = pdnmgThis->next;
        }
     }
     pdnmmHead = pdnmmHead->next;
  }

  pdnmsSlave = pmsdMaster->pdnmsSlaves;
  while(pdnmsSlave){
     pmsdSlave = (PMSD)pdnmsSlave->data.ptrvalue;
     pdnmmHead = pmsdSlave->pdnmmHead;
      while(pdnmmHead){
         pmmdThis = pdnmmHead->data.ptrvalue;
         if(pmmdThis){
            pdnmgThis = pmmdThis->pdnmgHead;
            while(pdnmgThis){
               pmgdThis = pdnmgThis->data.ptrvalue;
               if(pmgdThis && pmgdThis->bJustHighlighted == 1){
                  JustHLStatus = TRUE;
                  return;
               }
               pdnmgThis = pdnmgThis->next;
            }
         }
         pdnmmHead = pdnmmHead->next;
     }
     pdnmsSlave = pdnmsSlave->next;
  }     

}
/*---------------------------------------------------------*/
void LIBCALLBACK DoDeHighlightWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  pmgdThis = (PMGD) pfbThis;
  if(pmgdThis) {
     pdnmgThis = pmgdThis->pdnmgLink;
     if(pmgdThis->bHighlighted == 1){
        MediaObjSelect(pdnmgThis, FALSE);
        pmgdThis->bHighlighted = 0;
     }
  }   

}
/*---------------------------------------------------------*/
void Cn3D_DeHighlightAll(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){

     pmsdMaster = pdnmsMaster->data.ptrvalue;
     TraverseGraphs(pdnmsMaster, 0, 0, NULL, DoDeHighlightWithMGD);

      pdnmsSlave = pmsdMaster->pdnmsSlaves;
      while(pdnmsSlave) {
         TraverseGraphs(pdnmsSlave, 0, 0, NULL, DoDeHighlightWithMGD);
         pdnmsSlave = pdnmsSlave->next;
      }
  }
}
/*---------------------------------------------------------*/
void LIBCALLBACK DoCleanJustHLStatusWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  pmgdThis = (PMGD) pfbThis;
  if(pmgdThis) pmgdThis->bJustHighlighted = 0;

}
/*---------------------------------------------------------*/
void Cn3D_CleanJustHLStatus(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){

     pmsdMaster = pdnmsMaster->data.ptrvalue;
     TraverseGraphs(pdnmsMaster, 0, 0, NULL, DoCleanJustHLStatusWithMGD);

      pdnmsSlave = pmsdMaster->pdnmsSlaves;
      while(pdnmsSlave) {
         TraverseGraphs(pdnmsSlave, 0, 0, NULL, DoCleanJustHLStatusWithMGD);
         pdnmsSlave = pdnmsSlave->next;
      }
  }
}
/*---------------------------------------------------------*/
void LIBCALLBACK DoTurnOnSpecialFeatureWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  pmgdThis = (PMGD) pfbThis;
  if((pmgdThis) && pmgdThis->iUserDefinedFeature == iEditedFeature) pmgdThis->FeatureOn = 1;

}
/*---------------------------------------------------------*/
void  TurnOnSpecialFeatureWithMGD(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){

     pmsdMaster = pdnmsMaster->data.ptrvalue;
     TraverseGraphs(pdnmsMaster, 0, 0, NULL, DoTurnOnSpecialFeatureWithMGD);

      pdnmsSlave = pmsdMaster->pdnmsSlaves;
      while(pdnmsSlave) {
         TraverseGraphs(pdnmsSlave, 0, 0, NULL, DoTurnOnSpecialFeatureWithMGD);
         pdnmsSlave = pdnmsSlave->next;
      }
  }

}
/*---------------------------------------------------------*/
void LIBCALLBACK DoUnLinkSpecialFeatureWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  pmgdThis = (PMGD) pfbThis;
  if(pmgdThis && pmgdThis->iUserDefinedFeature == iEditedFeature){
     pmgdThis->iUserDefinedFeature = 0;
     pmgdThis->FeatureOn = 0;
  }

}
/*---------------------------------------------------------*/
void  UnLinkSpecialFeatureWithMGD(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;

  pdnmsMaster = GetSelectedModelstruc();
  if(pdnmsMaster != NULL){

     pmsdMaster = pdnmsMaster->data.ptrvalue;
     TraverseGraphs(pdnmsMaster, 0, 0, NULL, DoUnLinkSpecialFeatureWithMGD);

      pdnmsSlave = pmsdMaster->pdnmsSlaves;
      while(pdnmsSlave) {
         TraverseGraphs(pdnmsSlave, 0, 0, NULL, DoUnLinkSpecialFeatureWithMGD);
         pdnmsSlave = pdnmsSlave->next;
      }
  }

}
/*---------------------------------------*/
SpecialFeatureInfoPtr GetThisSpecialAlgorRenderSet(Int2 id)
{
 SpecialFeaturePtr sfpThis = NULL;
 SpecialFeatureInfoPtr sfipThis = NULL;

 sfpThis = sfp;
 while(sfpThis){
    if(sfpThis->choice == id){
       return sfpThis->data.ptrvalue;
    }
    sfpThis = sfpThis->next;
 }

 return NULL;

}
/*---------------------------------------------------------*/
void LIBCALLBACK DoLinkSpecialFeatureWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr)
{
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  SpecialFeatureInfoPtr sfipThis = NULL;

  pmgdThis = (PMGD) pfbThis;
  if(pmgdThis && pmgdThis->bHighlighted == 1){
     if(pmgdThis->iUserDefinedFeature != 0) {
        if(GetStatus(bFeatureAdd)){
           sfipThis = GetThisSpecialAlgorRenderSet(pmgdThis->iUserDefinedFeature);
           if(sfipThis != NULL) sfipThis->iRes--;
           FeatureOverwritten = TRUE;
           pmgdThis->iUserDefinedFeatureOld = 0;
           pmgdThis->FeatureOnOld = 0;
        }
        else {
           pmgdThis->iUserDefinedFeatureOld = pmgdThis->iUserDefinedFeature; 
           pmgdThis->FeatureOnOld = pmgdThis->FeatureOn;
        }
     }
     pmgdThis->iUserDefinedFeature = iAddedFeature;
     pmgdThis->FeatureOn = 1;
     pmgdThis->bJustHighlighted = 1;
     sfipThis = GetThisSpecialAlgorRenderSet(pmgdThis->iUserDefinedFeature);
     if(sfipThis != NULL) sfipThis->iRes++;
  }
  else if(pmgdThis && pmgdThis->bJustHighlighted == 1){
     if(pmgdThis->iUserDefinedFeature != 0) {
        if(GetStatus(bFeatureAdd)){
           sfipThis = GetThisSpecialAlgorRenderSet(pmgdThis->iUserDefinedFeature);
           if(sfipThis != NULL) sfipThis->iRes--;
           FeatureOverwritten = TRUE;
           pmgdThis->iUserDefinedFeatureOld = 0;
           pmgdThis->FeatureOnOld = 0;
        }
        else {
           pmgdThis->iUserDefinedFeatureOld = pmgdThis->iUserDefinedFeature; 
           pmgdThis->FeatureOnOld = pmgdThis->FeatureOn;
        }
     }
     pmgdThis->iUserDefinedFeature = iAddedFeature;
     pmgdThis->FeatureOn = 1;
     sfipThis = GetThisSpecialAlgorRenderSet(pmgdThis->iUserDefinedFeature);
     if(sfipThis != NULL) sfipThis->iRes++;
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
SpecialFeatureInfoPtr SpecialFeatureInfoFree( SpecialFeatureInfoPtr sfipThis)
{
  PRK prkThis = NULL;
  PARS parsThis = NULL;

  if(sfipThis->title) MemFree(sfipThis->title);
  if(sfipThis->description) MemFree(sfipThis->description);    
  if(sfipThis->parsSpecial) FreeAlgorRenderSet(sfipThis->parsSpecial);
/*if(sfipThis->prkSpecial)  FreeRenderKeep(sfipThis->prkSpecial);    */
  sfipThis = MemFree(sfipThis);
  
  return(sfipThis);
}
/*---------------------------------------------------------*/
SpecialFeaturePtr LIBCALL SpecialFeatureFree(SpecialFeaturePtr sfpThis)
{
  SpecialFeatureInfoPtr sfipThis = NULL;

  if(sfpThis == NULL) return NULL;
  sfipThis = sfpThis->data.ptrvalue;
  if(sfipThis)sfipThis = SpecialFeatureInfoFree(sfipThis);
 
  sfpThis->data.ptrvalue = NULL;
  sfpThis = ValNodeFree(sfpThis);
       /* or should be sfpThis = MemFree(sfpThis)-here sfpThis->next is */
       /* already set to know or is NULL by design */
       /* done this way so no redundant code needed */

  return(sfpThis);

}
/*---------------------------------------------------------*/
SpecialFeatureInfoPtr LIBCALL SpecialFeatureInfoNew()
{
  SpecialFeatureInfoPtr sfipNew = NULL;
  PRK prkNew = NULL;

  sfipNew = (SpecialFeatureInfoPtr) MemNew(sizeof(SpecialFeatureInfo)); 
  
  if(!sfipNew) return NULL;

/*prkNew =  NewRenderKeep();
  sfipNew->prkSpecial = prkNew;     */

  return sfipNew;
}
/*--------------------------------------------------------*/
SpecialFeatureInfoPtr GetEditedSpecialFeatureInfo()
{
  Int4 iFeature = 0, iCount = 0;
  SpecialFeatureInfoPtr sfipThis = NULL;
  SpecialFeaturePtr sfpThis = NULL;

  iFeature = GetValue(Cn3D_lFeature2);
  iCount = 0;

  sfpThis = sfp;
  while(sfpThis){
     iCount++;
     if(iFeature == iCount){
        sfipThis = sfpThis->data.ptrvalue;
        iEditedFeature = sfpThis->choice;
        return sfipThis;
     }
     sfpThis = sfpThis->next;
  }

  return NULL;

}
/*---------------------------------------------------------*/
PARS GetEditedSpecialFeatureInfoPAR(void) 
{
  SpecialFeatureInfoPtr sfipThis = NULL;
  PARS parsThis = NULL;

  sfipThis = GetEditedSpecialFeatureInfo();
  if(sfipThis) {
     parsThis = sfipThis->parsSpecial;
     if(parsThis) return parsThis;
     else {
        return NULL; 
     }
  }

  return NULL;
}
/*---------------------------------------------------------*/
PARS GetDefaultSpecialAlgorRenderSet(void)
{
  PARS parsThis = NULL;

  parsThis = NewAlgorRenderSet();

  parsThis->PVirtualBBOn = TRUE;
  parsThis->PRealBBOn = FALSE;
  parsThis->PExtraBBOn = FALSE;

  parsThis->PBBColor = C_red;
  parsThis->PResiduesOn = FALSE;
  parsThis->PResColor = C_red;

  parsThis->NTBBColor = C_red;
  parsThis->NTResColor = C_red;

  parsThis->IonsOn = FALSE;
  parsThis->ConnectOn = FALSE;
  parsThis->ObjectOn = FALSE;   /* these init setting is obsolete now */

  return parsThis;

}
/*---------------------------------------------------------*/
PARS LIBCALL GetSpecialAlgorRenderSet(void)
{
  SpecialFeaturePtr sfpThis = NULL;
  PARS parsThis = NULL;
  if(GetStatus(bFeatureEdit)){
     sfpThis = sfp;
     if(sfpThis == NULL) return NULL;
     parsThis = GetEditedSpecialFeatureInfoPAR();
     return parsThis;
  }
  else if(parsSpecial == NULL){
     parsSpecial = GetDefaultSpecialAlgorRenderSet();

     ResetModelCtrls();

     return parsSpecial;
  }
  else {
     return parsSpecial;
  }

}
/*---------------------------------------------------------*/
void FeatureTitleProc(TexT b)
{
  ResetModelCtrls();
}
/*---------------------------------------------------------*/
void ResetFeatureTitleProc(ButtoN b)
{
  SetTitle(FeatureTitle, "");
}
/*---------------------------------------------------------*/
static void ChangeSpecialRenderProc (void)
{
  Int2 i = 0, j = 0, k = 0;
  PARS pars = NULL;
  UIEnum val;

  pars = GetSpecialAlgorRenderSet();
  if (!pars) return;

  pars->PResiduesOn = GetStatus (Cn3D_lModelOnOffItem [2]);
  pars->NTResiduesOn = GetStatus (Cn3D_lModelOnOffItem [3]);

  pars->PVirtualBBOn = FALSE;
  pars->PRealBBOn = FALSE;
  pars->PExtraBBOn = FALSE;

  i = GetValue (Cn3D_pupModelPBB);
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

  pars->NTVirtualBBOn = FALSE;
  pars->NTRealBBOn = FALSE;
  pars->NTExtraBBOn = FALSE;

  j = GetValue (Cn3D_pupModelNABB);
  switch (j) {
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

  for (i = 0; i < 4; i++) {
     if (GetEnumPopup (Cn3D_pupModelColorStyle [i], colorAlist [i], &val)) {
        j = (Int2) val;
     } 
     else {
        j = 1;
     }

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
         case 0: /* prot bb */
            pars->PBBColor = k;
            break;
         case 1: /* prot sc */
            pars->NTBBColor = k;
            break;
         case 2: /* na bb */
            pars->PResColor = k;
            break;
         case 3: /* na sc */
            pars->NTResColor = k;
            break;
       }
   }

  for (i = 0; i < 4; i++) {
      if (GetEnumPopup (Cn3D_pupModelRenderStyle [i], renderAlist [i], &val)) {
         j = (Int2) val; } 
      else {
         j = 1;
      }
      switch (j)
        {
       case 1:
         k =  R_WIRE;
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
       case 5:
         k = R_SPACE;
         break;
       default:
      ;
     }
    switch (i)
      {
     case 0: /* prot bb */
         pars->PBBRender = k;
         break;
     case 1: /* na bb */
         pars->NTBBRender = k;
         break;
     case 2: /* prot sc */
         pars->PResRender = k;
         break;
     case 3: /* na sc */
         pars->NTResRender = k;
         break;
     default:
       ;
    }
  }

}
/*---------------------------------------------------------*/
static void ModelPopupOnOffProc (PopuP p)
{
   ChangeSpecialRenderProc();
}
/*---------------------------------------------------------*/
static void ModelButtonOnOffProc (ButtoN b)
{
   ChangeSpecialRenderProc();
}
/*---------------------------------------------------------*/
static PopuP ModelRenderStyle (GrouP h, Int2 i)
{
  PopuP  p;

  if (renderAlist [i] != empty_alist) {
    p = PopupList (h, FALSE, ModelPopupOnOffProc);
  } else {
    p = PopupList (h, FALSE, NULL);
  }
  InitEnumPopup (p, renderAlist [i], NULL);
  /*
  PopupItem (p, "Wireframe");
  PopupItem (p, "Tubes");
  PopupItem (p, "Ball&Stick");
  PopupItem (p, "Fat Tubes");
  PopupItem (p, "SpaceFill");
  */
  return p;
}
/*---------------------------------------------------------*/
static PopuP ModelColorStyle (GrouP h, Int2 i)
{
  PopuP  p;

  if (colorAlist [i] != empty_alist) {
    p = PopupList (h, FALSE, ModelPopupOnOffProc);
  } else {
    p = PopupList (h, FALSE, NULL);
  }
  InitEnumPopup (p, colorAlist [i], NULL);
  /*
  PopupItem (p, "red");
  PopupItem (p, "green");
  ......
  */
  return p;

}
/*---------------------------------------------------------*/
static void ChangeSpecialLabelsProc (void)
{
  Int2 i = 0, k = 0, nameval = 0;
  Uint1 codeval;
  Boolean NumOn = FALSE, NameOn = FALSE, PDBOn = FALSE, WhiteOn = FALSE, TopOn = FALSE;
  PARS pars = NULL;

  pars = GetSpecialAlgorRenderSet();
  if (!pars) return;

  pars->PBBLabelOn = (Uint1) GetValue(Cn3D_pupLabelAA);
  pars->NTBBLabelOn = (Uint1) GetValue(Cn3D_pupLabelNT);

  for (i = 0; i < 4; i++) {
  nameval = GetValue(Cn3D_bLName [i]);
  if (nameval == 3) codeval = (Uint1) L_1LETR;
  else codeval = (Uint1) L_3LETR;
  switch (i)
    {
     case 0: /* aa */
          /* clear bits */
         pars->PBBLabelStyle =  pars->PBBLabelStyle & ~((Uint1)( L_3LETR | L_1LETR)) ;
          /* set bit */
         pars->PBBLabelStyle = pars->PBBLabelStyle | (Uint1) codeval;
        break;
     case 1:  /* na   */
         pars->NTBBLabelStyle =  pars->NTBBLabelStyle & ~((Uint1)( L_3LETR |
L_1LETR)) ;
         pars->NTBBLabelStyle = pars->NTBBLabelStyle | (Uint1) codeval;
        break;
      default:;
     }

  NumOn = GetStatus(Cn3D_bLNum [i]);
  switch (i)
    {
     case 0: /* aa */
         if (NumOn) pars->PBBLabelStyle =  pars->PBBLabelStyle | (Uint1) L_NUM;
         else pars->PBBLabelStyle = pars->PBBLabelStyle & (Uint1) ~L_NUM;
        break;
     case 1:  /* na   */
         if (NumOn) pars->NTBBLabelStyle =  pars->NTBBLabelStyle | (Uint1) L_NUM;
         else pars->NTBBLabelStyle = pars->NTBBLabelStyle & (Uint1) ~L_NUM;
        break;
      default:
      ;
     }
  NameOn = (Boolean) (nameval > 1);
  switch (i)
    {
     case 0: /* aa */
         if (NameOn) pars->PBBLabelStyle =  pars->PBBLabelStyle | (Uint1) L_NAME;
         else pars->PBBLabelStyle = pars->PBBLabelStyle & (Uint1) ~L_NAME;
        break;
     case 1:  /* na   */
         if (NameOn) pars->NTBBLabelStyle =  pars->NTBBLabelStyle | (Uint1) L_NAME;
         else pars->NTBBLabelStyle = pars->NTBBLabelStyle & (Uint1) ~L_NAME;
        break;
      default:
      ;
    }
  PDBOn = (Boolean) (nameval == 4);
  switch (i)
    {
     case 0: /* aa */
         if (PDBOn) pars->PBBLabelStyle =  pars->PBBLabelStyle | (Uint1) L_PDB;
         else pars->PBBLabelStyle = pars->PBBLabelStyle & (Uint1) ~L_PDB;
        break;
     case 1:  /* na   */
         if (PDBOn) pars->NTBBLabelStyle =  pars->NTBBLabelStyle | (Uint1) L_PDB;
         else pars->NTBBLabelStyle = pars->NTBBLabelStyle & (Uint1) ~L_PDB;
        break;
      default:
      ;
     }
  WhiteOn = GetStatus(Cn3D_bLWhite [i]);
  switch (i)
    {
     case 0: /* aa */
         if (WhiteOn) pars->PBBLabelStyle =  pars->PBBLabelStyle | (Uint1) L_WHITE;
         else pars->PBBLabelStyle = pars->PBBLabelStyle & (Uint1) ~L_WHITE;
        break;
     case 1:  /* na   */
         if (WhiteOn) pars->NTBBLabelStyle =  pars->NTBBLabelStyle | (Uint1)
L_WHITE;
         else pars->NTBBLabelStyle = pars->NTBBLabelStyle & (Uint1) ~L_WHITE;        break;
   default:
   ;
  }
  TopOn = GetStatus(Cn3D_bLTop [i]);
  switch (i)
    {
     case 0: /* aa */
         if (TopOn) pars->PBBLabelJust =  pars->PBBLabelJust | (Uint1) LA_FRONT;
         else pars->PBBLabelJust = pars->PBBLabelJust & (Uint1) ~LA_FRONT;
        break;
     case 1:  /* na   */
         if (TopOn) pars->NTBBLabelJust =  pars->NTBBLabelJust | (Uint1) LA_FRONT;
         else pars->NTBBLabelJust = pars->NTBBLabelJust & (Uint1) ~LA_FRONT;
        break;
      default:
      ;
     }
  k = (Int2) GetValue(Cn3D_pupLabelSize [i]);
   switch (i)
    {
     case 0: /* prot bb */
         pars->PBBLabelScale = k;
         break;
     case 1: /* na bb */
         pars->NTBBLabelScale = k;
         break;
      default:
      ;
     }
  }

}
/*---------------------------------------------------------*/
static void ModelLabelPopupOnOffProc (PopuP p)
{
  ChangeSpecialLabelsProc ();
}
/*---------------------------------------------------------*/
static void ModelLabelButtonOnOffProc (ButtoN b)
{
  ChangeSpecialLabelsProc ();
}
/*---------------------------------------------------------*/
static void ModelSetRenderCtrls(PARS pars)
{
  Int2 i, j, k;
  PMSD pmsdThis = NULL;
  PDNMS pdnmsThis = NULL;
  PDNML pdnmlThis = NULL;
  PMLD  pmldThis = NULL;

  if(!pars) return;

  i = 4;
  if(pars->PVirtualBBOn) i = 1;
  else if(pars->PRealBBOn) i = 2;
  else if(pars->PExtraBBOn) i = 3;
  SafeSetValue(Cn3D_pupModelPBB, i);

  i = 4;
  if(pars->NTVirtualBBOn) i = 1;
  else if(pars->NTRealBBOn) i = 2;
  else if(pars->NTExtraBBOn) i = 3;
  SafeSetValue(Cn3D_pupModelNABB, i);
  
  SafeSetStatus(Cn3D_lModelOnOffItem[2], pars->PResiduesOn);
  SafeSetStatus(Cn3D_lModelOnOffItem[3], pars->NTResiduesOn);

  i = GetRenderItemIndex(pars->PBBRender);
  SafeSetValue(Cn3D_pupModelRenderStyle[0], i);
  i = GetColorItemIndex(pars->PBBColor);
  SafeSetValue(Cn3D_pupModelColorStyle[0], i);

  i = GetRenderItemIndex(pars->NTBBRender);
  SafeSetValue(Cn3D_pupModelRenderStyle[1], i);
  i = GetColorItemIndex(pars->NTBBColor);
  SafeSetValue(Cn3D_pupModelColorStyle[1], i); 

  i = GetRenderItemIndex(pars->PResRender);
  SafeSetValue(Cn3D_pupModelRenderStyle[2], i);
  i = GetColorItemIndex(pars->PResColor);
  SafeSetValue(Cn3D_pupModelColorStyle[2], i); 

  i = GetRenderItemIndex(pars->NTResRender);
  SafeSetValue(Cn3D_pupModelRenderStyle[3], i);
  i = GetColorItemIndex(pars->NTResColor);
  SafeSetValue(Cn3D_pupModelColorStyle[3], i); 

}
/*---------------------------------------------------------*/
static void  ModelSetLabelCtrls(PARS pars)
{
  Int2 i = 0, val = 0, style = 0, size = 0;
  PMSD pmsdThis = NULL;
  PDNMS pdnmsThis = NULL;

  if (!pars) return;

  SafeSetValue(Cn3D_pupLabelAA, pars->PBBLabelOn);
  SafeSetValue(Cn3D_pupLabelNT, pars->NTBBLabelOn);

  SafeSetStatus(Cn3D_lOnOffLabel [2], pars->PTermLabelOn);
  SafeSetStatus(Cn3D_lOnOffLabel [3], pars->NTTermLabelOn);

  for (i = 0; i < 2; i++) {
    val = 1;
    style = 1;
    size = 1;
    switch (i) {
      case 0 :
        style = pars->PBBLabelStyle;
        size = pars->PBBLabelScale;
        break;
      case 1 :
        style = pars->NTBBLabelStyle;
        size = pars->NTBBLabelScale;
        break;
      default :
        break;
    }
    if (! (Boolean) (style & (Uint1) (L_NAME))) {
      val = 1;
    } else if ((Boolean) (style & (Uint1) (L_3LETR))) {
      val = 2;
    } else if ((Boolean) (style & (Uint1) (L_1LETR))) {
      val = 3;
    } else if ((Boolean) (style & (Uint1) L_PDB)) {
      val = 4;
    }
    SafeSetValue(Cn3D_bLName [i], val);
    SafeSetStatus(Cn3D_bLNum [i], (Boolean) (style & (Uint1) L_NUM));
    SafeSetStatus(Cn3D_bLTop [i], (Boolean) (style & (Uint1) LA_FRONT));
    SafeSetStatus(Cn3D_bLWhite [i], (Boolean) (style & (Uint1) ~L_WHITE));
    SafeSetValue(Cn3D_pupLabelSize [i], size);
  }
}
/*---------------------------------------------------------*/
static void  ModelResetLabelCtrls(void)
{
  PARS pars = NULL;

  pars = GetSpecialAlgorRenderSet();
  if (!pars) return;

  ModelSetLabelCtrls(pars);
  
}
/*---------------------------------------------------------*/
void bFeatureAddAction(ButtoN b)
{
  if(GetStatus(bFeatureAdd)){
     SetStatus(bFeatureEdit, FALSE);
     SetStatus(bFeatureDelete, FALSE);
     SetStatus(bFeatureCancel, FALSE);
  }

}
/*---------------------------------------------------------*/
void bFeatureDeleteAction(ButtoN b)
{
/*ErrPostEx(SEV_ERROR, 0, 0, "Not function yet, but coming soon!  ");
  SetStatus(bFeatureDelete, FALSE);     */

  if(GetStatus(bFeatureDelete)){
     SetStatus(bFeatureEdit, FALSE);
     SetStatus(bFeatureAdd, FALSE);
     SetStatus(bFeatureCancel, FALSE);
  }

}
/*---------------------------------------------------------*/
void lFeature2Action(LisT l){
  SpecialFeatureInfoPtr sfip = NULL;
  PARS parsThis = NULL;
 
  Int4 val = 0;

  if(GetStatus(bFeatureEdit) || GetStatus(bFeatureDelete)){
     sfip = GetEditedSpecialFeatureInfo();
     if(sfip){
        SetTitle(FeatureTitle, sfip->title);
        SetTitle(FeatureDescription, sfip->description);
     }
     parsThis = GetEditedSpecialFeatureInfoPAR();
     ModelSetRenderCtrls(parsThis);
     ModelSetLabelCtrls(parsThis);
  }

}
/*---------------------------------------------------------*/
static void bFeatureCancelAction(ButtoN b)
{
  if(GetStatus(bFeatureCancel)){
     SetStatus(bFeatureAdd, FALSE);
     SetStatus(bFeatureDelete, FALSE);
     SetStatus(bFeatureEdit, FALSE);
  }
  SetValue(Cn3D_lFeature2, 0);

}
/*---------------------------------------------------------*/
void bFeatureEditAction(ButtoN b)
{
  Int4 iFeature = 0, iCount = 0;
  SpecialFeaturePtr sfpThis = NULL;
  SpecialFeatureInfoPtr sfip = NULL;
  PARS parsThis = NULL;

  if(GetStatus(bFeatureEdit)){
     SetStatus(bFeatureAdd, FALSE);
     SetStatus(bFeatureDelete, FALSE); 
  }

  if(GetStatus(bFeatureEdit)){
     sfip = GetEditedSpecialFeatureInfo();
     if(sfip){
        SetTitle(FeatureTitle, sfip->title);
        SetTitle(FeatureDescription, sfip->description);
     }
     parsThis = GetEditedSpecialFeatureInfoPAR();
     ModelSetRenderCtrls(parsThis);
     ModelSetLabelCtrls(parsThis);
  }

}
/*---------------------------------------------------------*/
void  Cn3D_GetlFeature3InfoProc(void)
{
  PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
  PMSD pmsdMaster = NULL, pmsdSlave = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmg = NULL;
  PMGD  pmgd = NULL;

  SeqIdPtr sipThis = NULL;
  SpecialFeaturePtr sfpThis = NULL;
  SpecialFeatureInfoPtr sfipThis = NULL;

  Int4 iCount = 1;

        sfpThis = sfp;
        while(sfpThis){
           sfipThis = sfpThis->data.ptrvalue;
           sfipThis->On = GetItemStatus(Cn3D_lFeature3, iCount);
           sfpThis = sfpThis->next;
           iCount++;
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
                    if(pmgd->iUserDefinedFeature > 0){
                       iCount = 1;
                       sfpThis = sfp;
                       while(sfpThis){
                          if(pmgd->iUserDefinedFeature == sfpThis->choice) {
                             pmgd->FeatureOn = GetItemStatus(Cn3D_lFeature3, iCount);
                             break;
                          }
                          sfpThis = sfpThis->next;
                          iCount++;
                       }
                    }
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
                       if(pmgd->iUserDefinedFeature > 0){
                          iCount = 1;
                          sfpThis = sfp;
                          while(sfpThis){
                             if(pmgd->iUserDefinedFeature == sfpThis->choice) {
                                pmgd->FeatureOn = GetItemStatus(Cn3D_lFeature3, iCount);
                                break;
                             }
                             sfpThis = sfpThis->next;
                             iCount++;
                          }
                       }
                    }
                 }

                 setoutB4:
                 pdnmmHead = pdnmmHead->next;
              }
              pdnmsSlave = pdnmsSlave->next;
           }
        }

}
/*---------------------------------------------------------*/
void RemoveInvalideUserDefinedFeatures(void)
{
  SpecialFeatureInfoPtr sfipThis = NULL;
  SpecialFeaturePtr sfpThis = NULL, sfp_curr = NULL;

  if(sfp->next == NULL){
     sfipThis = sfp->data.ptrvalue;
     if(sfipThis->iRes == 0) sfp = SpecialFeatureFree(sfp);
  }
  else {
     sfpThis = sfp;
     while(sfpThis){
        sfp_curr = sfpThis;
        sfpThis = sfpThis->next;

        sfipThis = sfp_curr->data.ptrvalue;
        if(sfipThis->iRes == 0) {
           sfp_curr = ValNodeExtract(&sfp, sfp_curr->choice);
           sfp_curr = SpecialFeatureFree(sfp_curr);
        }
     }
  }
}
/*---------------------------------------------------------*/
void Cn3D_ModelRedrawWrapper(ButtoN b)
{
  SpecialFeatureInfoPtr sfip = NULL;
  SpecialFeaturePtr sfpThis = NULL;
  PARS pars = NULL;
  PRK prKeep;
  Char str[1024];
  MsgAnswer ans;

  Int2 j, k, choice = 0;
  Int4 iFeature = 0, iCount = 0;

  if(GetStatus(bFeatureCancel)){

     SetValue(Cn3D_lFeature2, 0);
     SetStatus(bFeatureCancel, FALSE);    

     Cn3D_CheckHLStatus();
     Cn3D_CheckJustHLStatus();
     if(!HLStatus && !JustHLStatus) return;
     else {
        Cn3D_CleanJustHLStatus();
        ObjMgrDeSelectAll();
        Cn3D_DeHighlightAll();      

        Cn3D_ReColor = TRUE;
        Cn3D_RedrawProc(FALSE);
        Cn3D_ReColor = FALSE;

        return;
     }
  }
  
  if(GetStatus(bFeatureAdd)){
     Cn3D_CheckHLStatus();
     Cn3D_CheckJustHLStatus();
     if(!JustHLStatus && !HLStatus) {
        SetStatus(bFeatureAdd, FALSE);
        return;
     }
  } 

  Cn3D_GetlFeature3InfoProc();     

  pars = GetSpecialAlgorRenderSet();
  if (!pars) {
     if(GetStatus(bFeatureEdit)) {
        SetStatus(bFeatureEdit, FALSE);
     }
     return;
  }

  ChangeSpecialRenderProc();
  ChangeSpecialLabelsProc();

  if(GetStatus(bFeatureDelete)) {
     SetStatus(bFeatureDelete, FALSE);
     if(sfp == NULL) {
        return;
     }
     else if(sfp->next == NULL) {
        sfip = GetEditedSpecialFeatureInfo();
        if(sfip == NULL) return;

        iEditedFeature = sfp->choice;
        sfp = SpecialFeatureFree(sfp);
     }
     else {
        sfip = GetEditedSpecialFeatureInfo();
        if(sfip == NULL) return;
      
        sfpThis = ValNodeExtract(&sfp, iEditedFeature);
        sfpThis = SpecialFeatureFree(sfpThis);
     }
     UnLinkSpecialFeatureWithMGD();

     UpdateFeatureList2();
     Cn3D_ListOthersProc();

     Cn3D_ReColor = TRUE;
     Cn3D_RedrawProc(FALSE);
     Cn3D_ReColor = FALSE;

     return;
  }

  if(GetStatus(bFeatureEdit)){
     sfip = GetEditedSpecialFeatureInfo();
     if(!sfip->On){
        sfip->On = TRUE;
        TurnOnSpecialFeatureWithMGD();
     }
  }
  else {
     sfip = SpecialFeatureInfoNew();
     sfip->On = TRUE;
     sfip->parsSpecial = parsSpecial;  /* or pars */
     iAddedFeature++;
     ValNodeAddPointer(&sfp, iAddedFeature, (SpecialFeatureInfoPtr) sfip);

     Cn3D_CheckHLStatus();
     if(HLStatus) Cn3D_CleanJustHLStatus();

     LinkSpecialFeatureWithMGD();
  }
 
  Cn3D_ReColor = TRUE;
  Cn3D_RedrawProc(FALSE);
  Cn3D_ReColor = FALSE;

  if(GetStatus(bFeatureAdd) || GetStatus(bFeatureEdit)){
     GetTitle(FeatureTitle, str, sizeof(str));
     sfip->title = StringSave(str);
     GetTitle(FeatureDescription, str, sizeof(str));
     sfip->description = StringSave(str);    

     if(GetStatus(bFeatureAdd)) {
        if(HLStatus) {
           ObjMgrDeSelectAll();
           Cn3D_DeHighlightAll();
        }
          /* to prevent mistake on second apply push on user's side */
        Cn3D_CleanJustHLStatus();
        SetStatus(bFeatureAdd, FALSE);
     }
     else if (GetStatus(bFeatureEdit)) SetStatus(bFeatureEdit, FALSE);

     parsSpecial = NULL;
     SetTitle(FeatureTitle, "name?");
     SetTitle(FeatureDescription, "");
     UpdateFeatureList2();
     ResetModelCtrls();
  }
  else {

      if(sfp->next == NULL) sfp = SpecialFeatureFree(sfp);
      else { 
         sfpThis = sfp;
         while(sfpThis){
            choice = sfpThis->choice;
            sfpThis = sfpThis->next;
         }
         sfpThis = ValNodeExtract(&sfp, choice);
         sfpThis = SpecialFeatureFree(sfpThis);
      }
      ClearSpecialFeaturewithMGD();
      parsSpecial = NULL;
      iAddedFeature--;
  }    

  if(FeatureOverwritten) RemoveInvalideUserDefinedFeatures();
  FeatureOverwritten = FALSE;

  UpdateFeatureList2();
  Cn3D_ListOthersProc();
  
}
/*---------------------------------------------------------*/
GrouP LIBCALL ModelControls ( Nlm_GrouP prnt)
{
  GrouP g, h1, h2, h3, h4;
  GrouP h5, h6, h7;
  Int2  k = 0;
  PrompT ppt [10], ppt0[10];
  PrompT ppt1, ppt2, ppt3, ppt4, ppt5, ppt6, ppt7, ppt8, ppt9, ppt10, ppt11;
  PrompT ppt12, ppt13, ppt14;
  RecT pptPos, btnPos;
  Int2 delta = 0;

  g = HiddenGroup ( prnt, -1, 0, NULL );
  if (!g) return NULL;
  SetGroupSpacing ( g, 3, 9 );
#ifdef WIN_MOTIF
      SetGroupMargins ( g, 4, 1 );
      SetGroupSpacing ( g, 2, 1 );
#else
      SetGroupMargins ( g, 1, 1 );
      SetGroupSpacing ( g, 0, 0 );
#endif

  h5 = HiddenGroup (g, 0, -3, NULL);
  SetGroupSpacing (h5, 20, 3);
  ppt12 = StaticPrompt (h5, "Define Features: Use Mouse to Select Residues First", 0, stdLineHeight + 5, systemFont, 'c');
  h2 = HiddenGroup (h5, 0, -1, NULL);
  SetGroupSpacing (h2, 5, 3);
  StaticPrompt(h2, "    Name             ",0, stdLineHeight + 5, systemFont, 'c');
  FeatureTitle = DialogText(h2, "", 10, FeatureTitleProc);
  SetTitle(FeatureTitle, "name?");

  h3 = HiddenGroup (h5, 0, -1, NULL);
  SetGroupSpacing (h3, 5, 3);
  StaticPrompt(h3, "    Description ",0, stdLineHeight + 5, systemFont, 'c');
  FeatureDescription = DialogText(h3, "", 20, NULL) ;
  AlignObjects (ALIGN_LEFT, (HANDLE) ppt12, (HANDLE)h2, (HANDLE)h3, NULL);
  AlignObjects (ALIGN_LEFT, (HANDLE)h2, (HANDLE) h3, NULL);

/*StaticPrompt(g, "         ",0, stdLineHeight + 5, systemFont, 'l'); */
  h1 = HiddenGroup (g, 0, -5, NULL);
  SetGroupSpacing (h1, 5, 3);

  Cn3D_lModelOnOffItem [0] = NULL; /* Cn3D_pupModelPBB includes off setting */
  Cn3D_lModelOnOffItem [1] = NULL; /* Cn3D_pupModelNABB includes off setting */

  ppt1 = StaticPrompt (h1, "    Set Render", 0, stdLineHeight + 5, systemFont, 'c');
  ppt [0] = StaticPrompt (h1, "  Protein Backbone", 0, popupMenuHeight, programFont, 'l');
  ppt [1] = StaticPrompt (h1, "  Nucleotide Backbone", 0, popupMenuHeight, programFont, 'l');
  ppt [2] = StaticPrompt (h1, "  Protein Sidechains", 0, popupMenuHeight, programFont, 'l');
  ppt [3] = StaticPrompt (h1, "  Nucleotide Sidechains", 0, popupMenuHeight, programFont, 'l');
/*AlignObjects (ALIGN_CENTER, (HANDLE) ppt1, (HANDLE) ppt [1], NULL); */

  ppt2 = StaticPrompt (h1, "On/Off", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_pupModelPBB = PopupList(h1, FALSE, ModelPopupOnOffProc);
  PopupItem(Cn3D_pupModelPBB, "alpha C trace");
  PopupItem(Cn3D_pupModelPBB, "partial atom");
  PopupItem(Cn3D_pupModelPBB, "all atoms");
  PopupItem(Cn3D_pupModelPBB, "none");
  SafeSetValue(Cn3D_pupModelPBB, 1);
  Cn3D_pupModelNABB = PopupList(h1, FALSE, ModelPopupOnOffProc);
  PopupItem(Cn3D_pupModelNABB, "P trace");
  PopupItem(Cn3D_pupModelNABB, "partial atom");
  PopupItem(Cn3D_pupModelNABB, "all atoms");
  PopupItem(Cn3D_pupModelNABB, "none");
  SafeSetValue(Cn3D_pupModelNABB, 3);

  for (k = 2; k < 4; k++) {
    Cn3D_lModelOnOffItem [k] = CheckBox (h1, "", ModelButtonOnOffProc);
  }
  SafeSetStatus(Cn3D_lModelOnOffItem [3], TRUE);
  for (k = 2; k < 4; k++) {
     GetPosition (Cn3D_lModelOnOffItem [k], &btnPos);
     GetPosition (ppt [k], &pptPos);
     delta = (pptPos.bottom + pptPos.top) / 2 - (btnPos.bottom + btnPos.top) / 2;
     if(delta != 0) {
        OffsetRect (&btnPos, 0, delta);
        SetPosition (Cn3D_lModelOnOffItem [k], &btnPos);
        AdjustPrnt ((Nlm_GraphiC)Cn3D_lModelOnOffItem [k], &btnPos, FALSE);
     }
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt2, (HANDLE) Cn3D_pupModelPBB, (HANDLE)Cn3D_lModelOnOffItem [2], (HANDLE)Cn3D_lModelOnOffItem [3], NULL);

  ppt3 = StaticPrompt (h1, "Render", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_pupModelRenderStyle [0] = ModelRenderStyle (h1, 0);
  Cn3D_pupModelRenderStyle [1] = ModelRenderStyle (h1, 1);
  for (k = 2; k < 4; k++) {
      Cn3D_pupModelRenderStyle [k] = ModelRenderStyle (h1, k);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt3, (HANDLE) Cn3D_pupModelRenderStyle [0], NULL);

  ppt4 = StaticPrompt (h1, "Color", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_pupModelColorStyle [0] = ModelColorStyle (h1, 0);
  Cn3D_pupModelColorStyle [1] = ModelColorStyle (h1, 1);
  for (k = 2; k < 4; k++) {
     Cn3D_pupModelColorStyle [k] = ModelColorStyle (h1, k);
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt4, (HANDLE) Cn3D_pupModelColorStyle [0], NULL);

/*StaticPrompt(g, "",0,0,Nlm_systemFont,'l'); */

  h4 = HiddenGroup (g, 0, -3, NULL);
  SetGroupSpacing (h4, 5, 3);
  
  ppt5 = StaticPrompt (h4, "    Set Label", 0, stdLineHeight + 5, systemFont, 'c');
  ppt0 [0] = StaticPrompt (h4, "  Amino Acid", 0, popupMenuHeight, programFont, 'l');
  ppt0 [1] = StaticPrompt (h4, "  Nucleic Acid", 0, popupMenuHeight, programFont, 'l');
/*AlignObjects (ALIGN_CENTER, (HANDLE) ppt5, (HANDLE) ppt0 [1], NULL);  */

  ppt6 = StaticPrompt (h4, "On/Off", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_pupLabelAA = PopupList(h4, FALSE,  ModelLabelPopupOnOffProc);
  PopupItem(Cn3D_pupLabelAA, "none");
  PopupItem(Cn3D_pupLabelAA, "every AA");
  PopupItem(Cn3D_pupLabelAA, "every 5");
  PopupItem(Cn3D_pupLabelAA, "every 10");
  PopupItem(Cn3D_pupLabelAA, "every 20");
  PopupItem(Cn3D_pupLabelAA, "every 50");
  Cn3D_pupLabelNT= PopupList(h4, FALSE, ModelLabelPopupOnOffProc);
  PopupItem(Cn3D_pupLabelNT, "none");
  PopupItem(Cn3D_pupLabelNT, "every NA");
  PopupItem(Cn3D_pupLabelNT, "every 5");
  PopupItem(Cn3D_pupLabelNT, "every 10");
  PopupItem(Cn3D_pupLabelNT, "every 20");
  PopupItem(Cn3D_pupLabelNT, "every 50");
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt6, (HANDLE) Cn3D_pupLabelAA, (HANDLE) Cn3D_pupLabelNT, NULL);

  ppt7 = StaticPrompt (h4, "Name", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 2; k++) {
    Cn3D_bLName [k] = PopupList (h4, FALSE, ModelLabelPopupOnOffProc);
    PopupItem (Cn3D_bLName [k], "none");
    PopupItem (Cn3D_bLName [k], "3 letter");
    PopupItem (Cn3D_bLName [k], "1 letter");
    PopupItem (Cn3D_bLName [k], "PDB");
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt7, (HANDLE) Cn3D_bLName [0], (HANDLE) Cn3D_bLName [1], NULL);

  ppt8 = StaticPrompt (h4, "Num", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 2; k++) {
    Cn3D_bLNum [k] = CheckBox(h4, "", ModelLabelButtonOnOffProc ); 
  }
  for (k = 0; k < 2; k++) {
    GetPosition (Cn3D_bLNum [k], &btnPos);
    GetPosition (ppt0 [k], &pptPos);
    delta = (pptPos.bottom + pptPos.top) / 2 - (btnPos.bottom + btnPos.top) / 2;
    if (delta != 0) {
      OffsetRect (&btnPos, 0, delta);
      SetPosition (Cn3D_bLNum [k], &btnPos);
      AdjustPrnt ((Nlm_GraphiC) Cn3D_bLNum [k], &btnPos, FALSE);
    }
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt8, (HANDLE) Cn3D_bLNum [0], (HANDLE) Cn3D_bLNum [1], NULL);

  ppt9 = StaticPrompt (h4, "On Top", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 2; k++) {
    Cn3D_bLTop [k] = CheckBox(h4, "", ModelLabelButtonOnOffProc );
  }
  for (k = 0; k < 2; k++) {
    GetPosition (Cn3D_bLTop [k], &btnPos);
    GetPosition (ppt0 [k], &pptPos);
    delta = (pptPos.bottom + pptPos.top) / 2 - (btnPos.bottom + btnPos.top) / 2;
    if (delta != 0) {
      OffsetRect (&btnPos, 0, delta);
      SetPosition (Cn3D_bLTop [k], &btnPos);
      AdjustPrnt ((Nlm_GraphiC) Cn3D_bLTop [k], &btnPos, FALSE);
    }
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt9, (HANDLE) Cn3D_bLTop [0], (HANDLE) Cn3D_bLTop [1], NULL);

  ppt10 = StaticPrompt (h4, "White", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 2; k++) {
    Cn3D_bLWhite [k] = CheckBox(h4, "",  ModelLabelButtonOnOffProc); 
  }
  for (k = 0; k < 2; k++) {
    GetPosition (Cn3D_bLWhite [k], &btnPos);
    GetPosition (ppt0 [k], &pptPos);
    delta = (pptPos.bottom + pptPos.top) / 2 - (btnPos.bottom + btnPos.top) / 2;
    if (delta != 0) {
      OffsetRect (&btnPos, 0, delta);
      SetPosition (Cn3D_bLWhite [k], &btnPos);
      AdjustPrnt ((Nlm_GraphiC) Cn3D_bLWhite [k], &btnPos, FALSE);
    }
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt10, (HANDLE) Cn3D_bLWhite [0], (HANDLE) Cn3D_bLWhite [1], NULL);

  ppt11 = StaticPrompt (h4, "Size", 0, stdLineHeight + 5, systemFont, 'c');
  for (k = 0; k < 2; k++) {
    Cn3D_pupLabelSize [k] = PopupList(h4, FALSE, ModelLabelPopupOnOffProc);
    PopupItem(Cn3D_pupLabelSize [k], "1");
    PopupItem(Cn3D_pupLabelSize [k], "2");
    PopupItem(Cn3D_pupLabelSize [k], "3");
    PopupItem(Cn3D_pupLabelSize [k], "4");
    PopupItem(Cn3D_pupLabelSize [k], "5");
    PopupItem(Cn3D_pupLabelSize [k], "6");
    PopupItem(Cn3D_pupLabelSize [k], "7");
    PopupItem(Cn3D_pupLabelSize [k], "8");
    PopupItem(Cn3D_pupLabelSize [k], "9");
    PopupItem(Cn3D_pupLabelSize [k], "10");
  }
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt11, (HANDLE) Cn3D_pupLabelSize [0],
(HANDLE) Cn3D_pupLabelSize [1], NULL);

  StaticPrompt(g, "",0,0,Nlm_systemFont,'l');    

  h7 = HiddenGroup (g, 0, -2, NULL);
  ppt13 = StaticPrompt (h7, "Edit Features:", 0, stdLineHeight + 5, systemFont, 'c');

  h6 = HiddenGroup (h7, 0, -2, NULL);
  StaticPrompt(h6, "   ",0,0,Nlm_systemFont,'l');
  StaticPrompt(h6, "   ",0,0,Nlm_systemFont,'l');

  Cn3D_lFeature2 = SingleList(h6, 5, 3, lFeature2Action);
  StaticPrompt(h6, "",0,0,Nlm_systemFont,'l');

  bFeatureAdd = CheckBox(h6, "Save", bFeatureAddAction);
  bFeatureCancel = CheckBox(h6, "Cancel", bFeatureCancelAction);

  bFeatureEdit = CheckBox(h6, "Edit", bFeatureEditAction);
  bFeatureDelete = CheckBox(h6, "Delete", bFeatureDeleteAction);

/*AlignObjects (ALIGN_LEFT, (HANDLE)bFeatureAdd, (HANDLE)bFeatureDelete, (HANDLE)bFeatureEdit, NULL); */

  ppt14 = StaticPrompt (h7, " Display Features:", 0, stdLineHeight + 5, systemFont, 'c');
  Cn3D_lFeature3 = MultiList(h7, 5, 3, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE)g, (HANDLE)ppt14, (HANDLE)Cn3D_lFeature3, NULL);

  Cn3D_bModelApply = PushButton(g, "Apply!", Cn3D_ModelRedrawWrapper);

  AlignObjects (ALIGN_CENTER, (HANDLE) Cn3D_bModelApply, (HANDLE) h1,  NULL);

  ModelResetLabelCtrls();

  ResetModelCtrls();         

  return g;

}
/*---------------------------------------------------------*/
void LIBCALL ResetModelCtrls(void)
{
                                /* this function may not be neccessary */
  Int2 k = 0;
  PARS pars = NULL;;

  pars = GetSpecialAlgorRenderSet();
  if (!pars)
   {
      SafeSetValue(Cn3D_pupModelPBB, 1);
      SafeSetValue(Cn3D_pupModelNABB, 1);

      for(k = 2; k < 4; k++){
         SetItemStatus(Cn3D_lModelOnOffItem[k], 1, FALSE);
      }

      SafeSetValue(Cn3D_pupModelPBB, 4);
      SafeSetValue(Cn3D_pupModelNABB, 4);

      for(k = 0; k < 4; k++){
         SafeSetValue(Cn3D_pupModelRenderStyle[k], 1);
         SafeSetValue(Cn3D_pupModelColorStyle[k], 1);
      }

      return;
   }

  return;
}
/*-------------------------------------------*/
BiostrucFeatureSetDescrPtr GetDescrForUserDefinedFeature(SpecialFeaturePtr sfpThis)
{
  ValNodePtr descr = NULL, descr_head = NULL, descr_tail = NULL;
  SpecialFeatureInfoPtr sfipThis = NULL;

  if(!sfpThis) return NULL;
  sfipThis = sfpThis->data.ptrvalue;
  if(!sfipThis) return NULL;

  descr = ValNodeNew(NULL);
  descr->choice = 1;
  descr->data.ptrvalue = StringSave("User Defined Features");

  descr_head = descr;
  descr_tail = descr;

  if(sfipThis->description) {
     descr = ValNodeNew(NULL);
     descr->choice = 3;
     descr->data.ptrvalue = StringSave(sfipThis->description);
     descr_tail->next = descr;
     descr_tail = descr;
  }

  if(sfipThis->On){
     descr = ValNodeNew(NULL);
     descr->choice = 3;
     descr->data.ptrvalue = StringSave("On");
     descr_tail->next = descr;
     descr_tail = descr;
  }

  return descr_head;
}
/*-------------------------------------------*/
ValNodePtr GetBiostrucFeatureLocationForUserDefinedFeature(PDNMS pdnmsThis, Int2 iFeature)
{

  PMSD pmsdThis = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  ResidueIntervalPntrPtr rsip = NULL, rsip_head = NULL, rsip_tail = NULL;
  ValNodePtr rsp = NULL;  /* ResiduePntrsPtr */
  ChemGraphPntrsPtr pcgpThis = NULL;
  ValNodePtr llp = NULL;

  Int2 iMol = 0, iMolCurrent = 0;
  Int4 from = 0, to = 0;

  Boolean Started = FALSE, SingleInterval, UserDefinedFeatureFound = FALSE;
  Boolean UserDefinedFeatureFoundInThisMolecule = FALSE;

  pmsdThis = pdnmsThis->data.ptrvalue;
  pdnmmHead = pmsdThis->pdnmmHead;
  iMol = 0;
  while(pdnmmHead){
     pmmdThis = pdnmmHead->data.ptrvalue;
     if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout1;
     Started = FALSE;
     UserDefinedFeatureFoundInThisMolecule = FALSE;
     SingleInterval = TRUE;

     for(pdnmgThis = pmmdThis->pdnmgHead; pdnmgThis != NULL; pdnmgThis = pdnmgThis->next) {
        pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
        if(pmgdThis->iUserDefinedFeature == iFeature){
           UserDefinedFeatureFound = TRUE;
           iMolCurrent = iMol + 1;
           if(!Started) {
              from = pdnmgThis->choice;
              to   = pdnmgThis->choice;
              Started = TRUE;
              UserDefinedFeatureFoundInThisMolecule = TRUE;
           }
           else {
              if((pdnmgThis->choice - to) == 1) {
                 to = pdnmgThis->choice;
              }
              else {
                 SingleInterval = FALSE;
                 rsip = ResidueIntervalPntrNew();
/*               rsip->molecule_id = iMolCurrent; */
                                                  /* 5/3/99/  */
                 rsip->molecule_id = pmmdThis->iChainId;
                 rsip->from = from;
                 rsip->to = to;
                 if(!rsip_head){
                    rsip_head = rsip;
                    rsip_tail = rsip;
                 }
                 else {
                    rsip_tail->next = rsip;
                    rsip_tail = rsip;
                 }
 
                 from = pdnmgThis->choice;
                 to   = pdnmgThis->choice;
                 iMolCurrent = iMol + 1;
              }
           }
        }
     }
     if(UserDefinedFeatureFoundInThisMolecule){
        rsip = ResidueIntervalPntrNew();
/*      rsip->molecule_id = iMolCurrent;  */
        rsip->molecule_id = pmmdThis->iChainId;
                               /* 5/3/99 */
        rsip->from = from;
        rsip->to = to;
        if(!rsip_head){
           rsip_head = rsip;
           rsip_tail = rsip;
        }
        else {
           rsip_tail->next = rsip;
           rsip_tail = rsip;
        }                          /* for either last region or single region */
     }
        
     iMol++;
     setout1:
     pdnmmHead = pdnmmHead->next;
  }

  if(!UserDefinedFeatureFound) return NULL;

  rsp = ValNodeNew(NULL);
  rsp->choice = ResiduePntrs_interval;
  rsp->data.ptrvalue = rsip_head;

  pcgpThis = ValNodeNew(NULL);
  pcgpThis->choice = ChemGraphPntrs_residues;
  pcgpThis->data.ptrvalue = rsp;

  llp = ValNodeNew(NULL);
  llp->choice = Location_location_subgraph;
  llp->data.ptrvalue = pcgpThis;

  return llp;

}
/*-------------------------------------------*/
UserFieldPtr FillUserObjectField(UserFieldPtr ufp, ObjectIdPtr oip)
{
  ufp->choice = 1;
  ufp->label = oip;

  return ufp;

}
/*-------------------------------------------*/
UserObjectPtr GetUserObjectForUserDefinedFeature(PARS pars)
{
  UserObjectPtr ubp = NULL;
  UserFieldPtr head_ufp = NULL, curr_ufp = NULL, tail_ufp = NULL;
  ObjectIdPtr oip = NULL;
  Char str[20];
  Uint1Ptr rgb = NULL;

  head_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Backbone Representative");
  FillUserObjectField(head_ufp, oip);

  if(pars->PVirtualBBOn) head_ufp->data.ptrvalue = StringSave("Alpha C Trace");
  else if(pars->PRealBBOn) head_ufp->data.ptrvalue = StringSave("partial atom");
  else if(pars->PExtraBBOn) head_ufp->data.ptrvalue = StringSave("all atoms");
  else head_ufp->data.ptrvalue = StringSave("none");

  curr_ufp = head_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Backbone Render");
  FillUserObjectField(tail_ufp, oip);
  if(pars->PBBRender == R_WIRE) tail_ufp->data.ptrvalue = StringSave("Wire Frame");
  else if(pars->PBBRender == R_STICK) tail_ufp->data.ptrvalue = StringSave("Tubes");
  else if(pars->PBBRender == R_BALLNSTICK) tail_ufp->data.ptrvalue = StringSave("Ball & Stick");
  else if(pars->PBBRender == R_THICKWIRE) tail_ufp->data.ptrvalue = StringSave("Fat Tubes");
  else if(pars->PBBRender == R_SPACE) tail_ufp->data.ptrvalue = StringSave("Space Fill");

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Backbone Color");
  FillUserObjectField(tail_ufp, oip);
  rgb = Cn3d_PaletteRGB[pars->PBBColor].rgb;
  sprintf(str, "%d %d %d", *rgb,  *(rgb +1), *(rgb + 2));
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Side Chain Status");
  FillUserObjectField(tail_ufp, oip);
  if(pars->PResiduesOn) tail_ufp->data.ptrvalue = StringSave("On");
  else tail_ufp->data.ptrvalue = StringSave("Off");

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Side Chain Render");
  FillUserObjectField(tail_ufp, oip);
  if(pars->PResRender == R_WIRE) tail_ufp->data.ptrvalue = StringSave("Wire Frame");
  else if(pars->PResRender == R_STICK) tail_ufp->data.ptrvalue = StringSave("Tubes");
  else if(pars->PResRender == R_BALLNSTICK) tail_ufp->data.ptrvalue = StringSave("Ball & Stick");
  else if(pars->PResRender == R_THICKWIRE) tail_ufp->data.ptrvalue = StringSave("Fat Tubes");
  else if(pars->PResRender == R_SPACE) tail_ufp->data.ptrvalue = StringSave("Space Fill");

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Side Chain Color");
  FillUserObjectField(tail_ufp, oip);
  rgb = Cn3d_PaletteRGB[pars->PResColor].rgb;
  sprintf(str, "%d %d %d", *rgb,  *(rgb +1), *(rgb + 2));
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Label Status");
  FillUserObjectField(tail_ufp, oip);
  sprintf(str, "%d", (Uint1) pars->PBBLabelOn);
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Label Just");
  FillUserObjectField(tail_ufp, oip);
  sprintf(str, "%d", (Uint1) pars->PBBLabelJust);
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Label Scale");
  FillUserObjectField(tail_ufp, oip);
  sprintf(str, "%d", (Uint1) pars->PBBLabelScale);
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Protein Label Style");
  FillUserObjectField(tail_ufp, oip);
  sprintf(str, "%d", (Uint1) pars->PBBLabelStyle);
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Backbone Representative");
  FillUserObjectField(tail_ufp, oip);

  if(pars->NTVirtualBBOn) tail_ufp->data.ptrvalue = StringSave("Alpha C Trace");
  else if(pars->NTRealBBOn) tail_ufp->data.ptrvalue = StringSave("partial atom");
  else if(pars->NTExtraBBOn) tail_ufp->data.ptrvalue = StringSave("all atoms");
  else tail_ufp->data.ptrvalue = StringSave("none");

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Backbone Render");
  FillUserObjectField(tail_ufp, oip);
  if(pars->NTBBRender == R_WIRE) tail_ufp->data.ptrvalue = StringSave("Wire Frame");
  else if(pars->NTBBRender == R_STICK) tail_ufp->data.ptrvalue = StringSave("Tubes");
  else if(pars->NTBBRender == R_BALLNSTICK) tail_ufp->data.ptrvalue = StringSave("Ball & Stick");
  else if(pars->NTBBRender == R_THICKWIRE) tail_ufp->data.ptrvalue = StringSave("Fat Tubes");
  else if(pars->NTBBRender == R_SPACE) tail_ufp->data.ptrvalue = StringSave("Space Fill");

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Backbone Color");
  FillUserObjectField(tail_ufp, oip);
  rgb = Cn3d_PaletteRGB[pars->NTBBColor].rgb;
  sprintf(str, "%d %d %d", *rgb,  *(rgb +1), *(rgb + 2));
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Side Chain Status");
  FillUserObjectField(tail_ufp, oip);
  if(pars->NTResiduesOn) tail_ufp->data.ptrvalue = StringSave("On");
  else tail_ufp->data.ptrvalue = StringSave("Off");

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Side Chain Render");
  FillUserObjectField(tail_ufp, oip);
  if(pars->NTResRender == R_WIRE) tail_ufp->data.ptrvalue = StringSave("Wire Frame");
  else if(pars->NTResRender == R_STICK) tail_ufp->data.ptrvalue = StringSave("Tubes");
  else if(pars->NTResRender == R_BALLNSTICK) tail_ufp->data.ptrvalue = StringSave("Ball & Stick");
  else if(pars->NTResRender == R_THICKWIRE) tail_ufp->data.ptrvalue = StringSave("Fat Tubes");
  else if(pars->NTResRender == R_SPACE) tail_ufp->data.ptrvalue = StringSave("Space Fill");

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Side Chain Color");
  FillUserObjectField(tail_ufp, oip);
  rgb = Cn3d_PaletteRGB[pars->NTResColor].rgb;
  sprintf(str, "%d %d %d", *rgb,  *(rgb +1), *(rgb + 2));
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Label Status");
  FillUserObjectField(tail_ufp, oip);
  sprintf(str, "%d", (Uint1) pars->NTBBLabelOn);
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Label Just");
  FillUserObjectField(tail_ufp, oip);
  sprintf(str, "%d", (Uint1) pars->NTBBLabelJust);
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Label Scale");
  FillUserObjectField(tail_ufp, oip);
  sprintf(str, "%d", (Uint1) pars->NTBBLabelScale);
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;

  tail_ufp = UserFieldNew();
  oip = ObjectIdNew();
  oip->str = StringSave("Nucleic Acid Label Style");
  FillUserObjectField(tail_ufp, oip);
  sprintf(str, "%d", (Uint1) pars->NTBBLabelStyle);
  tail_ufp->data.ptrvalue = StringSave(str);

  curr_ufp->next = tail_ufp;
  curr_ufp = tail_ufp;
  oip = ObjectIdNew();
  oip->str = StringSave("Cn3D Rendering");

  ubp = UserObjectNew();
  ubp->type = oip;
  ubp->data = head_ufp;
  ubp->_class = StringSave ("Property for User Defined Feature");

  return (ubp);
}
/*-------------------------------------------*/
ValNodePtr GetBiostrucFeaturePropertyForUserDefinedFeature(SpecialFeaturePtr sfpThis)
{
  ValNodePtr ppp = NULL;
  SpecialFeatureInfoPtr sfipThis = NULL;
  PARS parsThis = NULL;

  sfipThis = sfpThis->data.ptrvalue;
  parsThis = sfipThis->parsSpecial;

  ppp = ValNodeNew(NULL);
  ppp->choice = 6;
  ppp->data.ptrvalue = GetUserObjectForUserDefinedFeature(parsThis);

  return ppp;
} 
/*-------------------------------------------*/
BiostrucFeatureSetPtr GetBiostrucFeatureSetForUserDefinedFeature(PDNMS pdnmsThis, BiostrucPtr bsp, SpecialFeaturePtr sfpThis, Int2 iCount)
{
  BiostrucFeatureSetPtr bsfsp = NULL, bsfsp_head = NULL, bsfsp_tail = NULL;
  BiostrucFeaturePtr bsfp = NULL;

  SpecialFeatureInfoPtr sfipThis = NULL;

  sfipThis = sfpThis->data.ptrvalue;

  bsfp = BiostrucFeatureNew();
  bsfp->name = StringSave(sfipThis->title);
  bsfp->type = Feature_type_other;
  bsfp->Location_location = GetBiostrucFeatureLocationForUserDefinedFeature(pdnmsThis, sfpThis->choice);    

  if(!bsfp->Location_location){
     bsfp = BiostrucFeatureFree(bsfp);
     return NULL;
  }

  bsfp->Property_property = GetBiostrucFeaturePropertyForUserDefinedFeature(sfpThis);

  bsfsp = BiostrucFeatureSetNew();
  bsfsp->id = iCount;
  bsfsp->descr = GetDescrForUserDefinedFeature(sfpThis); 
  bsfsp->features = bsfp;

  return bsfsp;

}
/*-------------------------------------------*/
PDNMS Cn3DAddUserDefinedFeature(PDNMS pdnmsThis)
{
  BiostrucPtr bsp = NULL;
  PMSD pmsdThis = NULL;
  BiostrucFeatureSetPtr bsfsp = NULL, bsfsp_head = NULL, bsfsp_tail = NULL;
  SpecialFeaturePtr sfpThis = NULL;
  Int2 iCount;

  pmsdThis = pdnmsThis->data.ptrvalue;
  bsp = pmsdThis->pbsBS;

  sfpThis = sfp;
  iCount = 1;
  while(sfpThis){
     bsfsp = GetBiostrucFeatureSetForUserDefinedFeature(pdnmsThis, bsp, sfpThis, iCount);
     if(bsfsp){
        if(!bsfsp_head){
           bsfsp_head = bsfsp;
           bsfsp_tail = bsfsp;
        }
        else {
           bsfsp_tail->next = bsfsp;
           bsfsp_tail = bsfsp;
        }
     }
     
     sfpThis = sfpThis->next;
     iCount++;
  }

  bsfsp = bsp->features;
  while(bsfsp->next){
     bsfsp = bsfsp->next;
  }

  bsfsp->next = bsfsp_head;

  return pdnmsThis;

}
/*---------------------------------------*/
SpecialFeatureInfoPtr AdditionalUserDefinedFeature(Int2 id)
{
 SpecialFeaturePtr sfpThis = NULL;
 SpecialFeatureInfoPtr sfipThis = NULL;

 sfpThis = sfpThis_head;
 while(sfpThis){
    if(sfpThis->choice == id){
       return sfpThis->data.ptrvalue;
    }
    sfpThis = sfpThis->next;
 }

 return NULL;

}
/*-------------------------------------*/
void Cn3DIndexUserDefinedFeatureForMGD(PMSD pmsdThis, ValNodePtr llp, Int2 iFeature, Boolean FeatureOn)
{
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;
  PDNMG pdnmgThis = NULL;
  PMGD  pmgdThis = NULL;

  ChemGraphPntrsPtr pcgpThis = NULL;
  ValNodePtr rsp = NULL;  /* ResiduePntrsPtr */
  ResidueIntervalPntrPtr rsip = NULL;

  SpecialFeatureInfoPtr sfipThis = NULL;

  sfipThis = AdditionalUserDefinedFeature(iFeature);

  pcgpThis = llp->data.ptrvalue;
  rsp = pcgpThis->data.ptrvalue;
  rsip = rsp->data.ptrvalue;
  while(rsip){
     pdnmmHead = pmsdThis->pdnmmHead;
     while(pdnmmHead){
        pmmdThis = pdnmmHead->data.ptrvalue;
        if(pmmdThis->bWhat != (Byte) AM_PROT && pmmdThis->bWhat != (Byte) AM_RNA && pmmdThis->bWhat != (Byte) AM_DNA) goto setout1; 
        if(rsip->molecule_id == pmmdThis->iChainId){ 
           for(pdnmgThis = pmmdThis->pdnmgHead; pdnmgThis != NULL; pdnmgThis = pdnmgThis->next) {
              if(pdnmgThis->choice >= rsip->from && pdnmgThis->choice <= rsip->to){ 
                 pmgdThis = (PMGD) pdnmgThis->data.ptrvalue;
                 pmgdThis->iUserDefinedFeature = iFeature;
                 pmgdThis->FeatureOn = FeatureOn;
                 if(sfipThis) sfipThis->iRes++;
              }
           }
        }
        setout1:
        pdnmmHead = pdnmmHead->next;

     }

     rsip = rsip->next;
  }


}
/*----------------------------------------------*/
Int2 GetCn3d_PaletteRGBIndex(Uint1Ptr rgb)
{
  ResidueColorCellPtr prgb = NULL;
  Int2 iCount = 0, iColor = 0;

  prgb = Cn3d_PaletteRGB;
  for(iCount = 0; iCount < CN3D_COLOR_MAX; iCount++){
     if(*rgb == prgb->rgb[0] &&  *(rgb + 1) == prgb->rgb[1] && *(rgb + 2) == prgb->rgb[2]) {
        iColor = iCount;
        break;
     }

     prgb++;
  }

  return(iColor);

}
/*-------------------------------------*/
PARS GetSpecialParsForUserDefinedFeature(UserObjectPtr ubp)
{

  PARS pars = NULL;
  UserFieldPtr ufp = NULL;
  ObjectIdPtr oip = NULL;
  
  Uint1 rgb[3];
  short rgb_short[3], num;

/*pars = (PARS) MemNew((size_t)(sizeof(ARS))); */
  pars = GetDefaultSpecialAlgorRenderSet();

  ufp = ubp->data;
  while(ufp){
     oip = ufp->label;
     if(StringCmp(oip->str, "Protein Backbone Representative") == 0){
        if(StringCmp(ufp->data.ptrvalue, "Alpha C Trace") == 0) {
           pars->PVirtualBBOn = TRUE;
           pars->PRealBBOn = FALSE;
           pars->PExtraBBOn = FALSE;
        }
        else if (StringCmp(ufp->data.ptrvalue, "all atoms") == 0) {
           pars->PRealBBOn = FALSE;
           pars->PVirtualBBOn = FALSE;
           pars->PExtraBBOn = TRUE;
        }
        else if (StringCmp(ufp->data.ptrvalue, "partial atom") == 0) {
           pars->PExtraBBOn = FALSE;
           pars->PVirtualBBOn = FALSE;
           pars->PRealBBOn = TRUE;
        }
     }
     else if(StringCmp(oip->str, "Protein Backbone Render") == 0){
        if(StringCmp(ufp->data.ptrvalue, "Wire Frame") == 0) pars->PBBRender = R_WIRE;
        else if(StringCmp(ufp->data.ptrvalue, "Tubes") == 0) pars->PBBRender = R_STICK;
        else if(StringCmp(ufp->data.ptrvalue, "Ball & Stick") == 0) pars->PBBRender = R_BALLNSTICK;
        else if(StringCmp(ufp->data.ptrvalue, "Fat Tubes") == 0) pars->PBBRender = R_THICKWIRE;
        else if(StringCmp(ufp->data.ptrvalue, "Space Fill") == 0) pars->PBBRender = R_SPACE;
     }
     else if(StringCmp(oip->str, "Protein Backbone Color") == 0){
        sscanf(ufp->data.ptrvalue, "%hd %hd %hd", &rgb_short[0], &rgb_short[1], &rgb_short[2]);        
        rgb[0] = (Uint1) rgb_short[0]; rgb[1] = (Uint1) rgb_short[1]; rgb[2] = (Uint1) rgb_short[2];    
        pars->PBBColor = GetCn3d_PaletteRGBIndex(rgb); 
    
     }
     else if(StringCmp(oip->str, "Protein Side Chain Status") == 0){
        if(StringCmp(ufp->data.ptrvalue, "On") == 0) pars->PResiduesOn = TRUE;
        else pars->PResiduesOn = FALSE;
     }
     else if(StringCmp(oip->str, "Protein Side Chain Render") == 0){
        if(StringCmp(ufp->data.ptrvalue, "Wire Frame") == 0) pars->PResRender = R_WIRE;
        else if(StringCmp(ufp->data.ptrvalue, "Tubes") == 0) pars->PResRender =
R_STICK;
        else if(StringCmp(ufp->data.ptrvalue, "Ball & Stick") == 0) pars->PResRender = R_BALLNSTICK;
        else if(StringCmp(ufp->data.ptrvalue, "Fat Tubes") == 0) pars->PResRender = R_THICKWIRE;
        else if(StringCmp(ufp->data.ptrvalue, "Space Fill") == 0) pars->PResRender = R_SPACE;
     }
     else if(StringCmp(oip->str, "Protein Side Chain Color") == 0){
        sscanf(ufp->data.ptrvalue, "%hd %hd %hd", &rgb_short[0], &rgb_short[1], &rgb_short[2]);        
        rgb[0] = (Uint1) rgb_short[0]; rgb[1] = (Uint1) rgb_short[1]; rgb[2] = (Uint1) rgb_short[2];    
        pars->PResColor = GetCn3d_PaletteRGBIndex(rgb);
     }
     else if(StringCmp(oip->str, "Protein Label Status") == 0){
        sscanf(ufp->data.ptrvalue, "%hd", &num);
        pars->PBBLabelOn = (Uint1) num;
     }
     else if(StringCmp(oip->str, "Protein Label Just") == 0){
        sscanf(ufp->data.ptrvalue, "%hd", &num);
        pars->PBBLabelJust = (Uint1) num;
     }
     else if(StringCmp(oip->str, "Protein Label Scale") == 0){
        sscanf(ufp->data.ptrvalue, "%hd", &num);
        pars->PBBLabelScale = (Uint1) num;
     }
     else if(StringCmp(oip->str, "Protein Label Style") == 0){
        sscanf(ufp->data.ptrvalue, "%hd", &num);
        pars->PBBLabelStyle = (Uint1) num;
     } 
    if(StringCmp(oip->str, "Nucleic Acid Backbone Representative") == 0){
        if(StringCmp(ufp->data.ptrvalue, "Alpha C Trace") == 0) {
           pars->NTVirtualBBOn = TRUE;
           pars->NTRealBBOn = FALSE;
           pars->NTExtraBBOn = FALSE;
        }
        else if (StringCmp(ufp->data.ptrvalue, "all atoms") == 0) {
           pars->NTRealBBOn = FALSE;
           pars->NTVirtualBBOn = FALSE;
           pars->NTExtraBBOn = TRUE;
        }
        else if (StringCmp(ufp->data.ptrvalue, "partial atom") == 0) {
           pars->NTExtraBBOn = FALSE;
           pars->NTVirtualBBOn = FALSE;
           pars->NTRealBBOn = TRUE;
        }
     }
     else if(StringCmp(oip->str, "Nucleic Acid Backbone Render") == 0){
        if(StringCmp(ufp->data.ptrvalue, "Wire Frame") == 0) pars->NTBBRender = R_WIRE;
        else if(StringCmp(ufp->data.ptrvalue, "Tubes") == 0) pars->NTBBRender = R_STICK;
        else if(StringCmp(ufp->data.ptrvalue, "Ball & Stick") == 0) pars->NTBBRender = R_BALLNSTICK;
        else if(StringCmp(ufp->data.ptrvalue, "Fat Tubes") == 0) pars->NTBBRender = R_THICKWIRE; 
        else if(StringCmp(ufp->data.ptrvalue, "Space Fill") == 0) pars->NTBBRender = R_SPACE;
     }
     else if(StringCmp(oip->str, "Nucleic Acid Backbone Color") == 0){
        sscanf(ufp->data.ptrvalue, "%hd %hd %hd", &rgb_short[0], &rgb_short[1], &rgb_short[2]);

        rgb[0] = (Uint1) rgb_short[0]; rgb[1] = (Uint1) rgb_short[1]; rgb[2] = (Uint1) rgb_short[2];
        pars->NTBBColor = GetCn3d_PaletteRGBIndex(rgb);

     }
     else if(StringCmp(oip->str, "Nucleic Acid Side Chain Status") == 0){
        if(StringCmp(ufp->data.ptrvalue, "On") == 0) pars->NTResiduesOn = TRUE;
        else pars->NTResiduesOn = FALSE;
     }
     else if(StringCmp(oip->str, "Nucleic Side Chain Render") == 0){
        if(StringCmp(ufp->data.ptrvalue, "Wire Frame") == 0) pars->NTResRender = R_WIRE;
        else if(StringCmp(ufp->data.ptrvalue, "Tubes") == 0) pars->NTResRender =
R_STICK;
        else if(StringCmp(ufp->data.ptrvalue, "Ball & Stick") == 0) pars->NTResRender = R_BALLNSTICK;
        else if(StringCmp(ufp->data.ptrvalue, "Fat Tubes") == 0) pars->NTResRender = R_THICKWIRE;
        else if(StringCmp(ufp->data.ptrvalue, "Space Fill") == 0) pars->NTResRender = R_SPACE;
     }
     else if(StringCmp(oip->str, "Nucleic Acid Side Chain Color") == 0){
        sscanf(ufp->data.ptrvalue, "%hd %hd %hd", &rgb_short[0], &rgb_short[1], &rgb_short[2]);

        rgb[0] = (Uint1) rgb_short[0]; rgb[1] = (Uint1) rgb_short[1]; rgb[2] = (Uint1) rgb_short[2];
        pars->NTResColor = GetCn3d_PaletteRGBIndex(rgb);
     }
     else if(StringCmp(oip->str, "Nucleic Acid Label Status") == 0){
        sscanf(ufp->data.ptrvalue, "%hd", &num);
        pars->NTBBLabelOn = (Uint1) num;
     }
     else if(StringCmp(oip->str, "Nucleic Acid Label Just") == 0){
        sscanf(ufp->data.ptrvalue, "%hd", &num);
        pars->NTBBLabelJust = (Uint1) num;
     }
     else if(StringCmp(oip->str, "Nucleic Acid Label Scale") == 0){
        sscanf(ufp->data.ptrvalue, "%hd", &num);
        pars->NTBBLabelScale = (Uint1) num;
     }
     else if(StringCmp(oip->str, "Nucleic Acid Label Style") == 0){
        sscanf(ufp->data.ptrvalue, "%hd", &num);
        pars->NTBBLabelStyle = (Uint1) num;
     }

     ufp = ufp->next;
  }

  return pars;
}
/*-------------------------------------*/
void Cn3DIndexUserDefinedFeatureForMSD(PDNMS pdnmsThis)
{
  PMSD pmsdThis = NULL;

  BiostrucPtr bsp = NULL;
  BiostrucFeatureSetPtr bsfsp = NULL;
  BiostrucFeatureSetPtr bsfsp_prev = NULL, bsfsp_curr = NULL;
  BiostrucFeaturePtr bsfp = NULL;
  ValNodePtr descr = NULL;  
  SpecialFeaturePtr sfpThis = NULL; 
  SpecialFeatureInfoPtr sfipThis = NULL;

  UserObjectPtr ubp = NULL;
  ValNodePtr llp = NULL, ppp = NULL;

  Boolean  UserDefinedFeatureFound = FALSE;
  Boolean NewUserDefinedFeature = FALSE;

  if(pdnmsThis == NULL) return;


  pmsdThis = pdnmsThis->data.ptrvalue;
  bsp = pmsdThis->pbsBS;
  bsfsp = bsp->features;
  while(bsfsp){
     descr = bsfsp->descr;
     UserDefinedFeatureFound = FALSE;
     NewUserDefinedFeature = FALSE;
     while(descr){
        if(descr->choice == 1){
           if(StringCmp(descr->data.ptrvalue, "User Defined Features") == 0){
              UserDefinedFeatureFound = TRUE;
              sfipThis = AdditionalUserDefinedFeature(bsfsp->id);
              if(sfipThis == NULL){
                 sfipThis = SpecialFeatureInfoNew();
                 sfipThis->On = FALSE;
                 NewUserDefinedFeature = TRUE;
              }
           }
        }
        else if(descr->choice == 3){
           if(NewUserDefinedFeature){
              if(StringCmp(descr->data.ptrvalue, "On") == 0){
                 sfipThis->On = TRUE;
              }
              else if(descr->data.ptrvalue != NULL) {
                 sfipThis->description = StringSave(descr->data.ptrvalue);
              }
           }
        }

        descr = descr->next;
     }
     if(UserDefinedFeatureFound){
        bsfp = bsfsp->features;
        llp = bsfp->Location_location;

        if(NewUserDefinedFeature){
           sfipThis->title = StringSave(bsfp->name);
           ppp = bsfp->Property_property;
           ubp = ppp->data.ptrvalue;
           sfipThis->parsSpecial = GetSpecialParsForUserDefinedFeature(ubp); 

           sfpThis = ValNodeNew(NULL);
           sfpThis->choice = bsfsp->id;
           sfpThis->data.ptrvalue = sfipThis;
           if(!sfpThis_head){
              sfpThis_head = sfpThis;
              sfpThis_tail = sfpThis;
           }
           else {
              sfpThis_tail->next = sfpThis;
              sfpThis_tail = sfpThis;
           }        
        }

        Cn3DIndexUserDefinedFeatureForMGD(pmsdThis, llp, bsfsp->id, sfipThis->On);
     }
    
     bsfsp = bsfsp->next;
  }
              
  bsfsp = bsp->features;
  while(bsfsp){
     descr = bsfsp->descr;
     if(descr->choice == 1){
        if(StringCmp(descr->data.ptrvalue, "User Defined Features") == 0){
           if(bsfsp_prev) bsfsp_prev->next = NULL;
           bsfsp_curr = bsfsp;
           break;
         }
      }
      bsfsp_prev = bsfsp;
      bsfsp = bsfsp->next;
  }

  bsfsp_curr = BiostrucFeatureSetFree(bsfsp_curr);

}
/*-------------------------------------*/
void Cn3DIndexUserDefinedFeature()
{
  PDNMS pdnmsThis = NULL, pdnmsThisSlave = NULL;
  PMSD  pmsdThis = NULL;
  PDNMM pdnmmHead = NULL;
  PMMD  pmmdThis = NULL;

  SpecialFeaturePtr sfpThis = NULL;     

  pdnmsThis = GetSelectedModelstruc();
  if(!pdnmsThis) return;

  Cn3DIndexUserDefinedFeatureForMSD(pdnmsThis);
  
  pdnmsThisSlave = ((PMSD)(pdnmsThis->data.ptrvalue))->pdnmsSlaves;
  while(pdnmsThisSlave){
     Cn3DIndexUserDefinedFeatureForMSD(pdnmsThisSlave);
     pdnmsThisSlave = pdnmsThisSlave->next;
  } 

  sfp = sfpThis_head;

  sfpThis = sfp;
  iAddedFeature = 0;
  while(sfpThis){
     if(sfpThis->choice > iAddedFeature) iAddedFeature = sfpThis->choice;
     sfpThis = sfpThis->next;
  }     

  Cn3D_ListOthersProc();
  UpdateFeatureList2();

}    
/*-------------------------------------------*/
void ClearSpecialFeature(void)
{
  SpecialFeaturePtr sfpThis = NULL, next = NULL;

  sfpThis = sfp;
  while(sfpThis){
     next = sfpThis->next;
     sfpThis->next = NULL;
     sfpThis = SpecialFeatureFree(sfpThis);
     sfpThis = next;
  }

/*ValNodeFree(sfp); */

  sfp = NULL; sfpThis_head = NULL; sfpThis_tail = NULL;
  parsSpecial = NULL;
  iAddedFeature = 0; 
  iEditedFeature = 0;

}
/*-------------------------------------------*/
void ClearDomainData(void)
{
  Int4 iCount = 0;

  for(iCount = 0; iCount < iDomainCount; iCount++){
     MemFree(domaindata[iCount]);
  }

  MemFree(domaindata);

  domaindata = NULL;
  iDomainCount = 0;

  return;
}
/*-------------------------------------------*/
void ClearRest(void)
{

  if(sfp != NULL) ClearSpecialFeature();

  if(iDomainCount > 0) ClearDomainData();

  Cn3D_DisplayHighlight = FALSE;
  Cn3D_NoSingleHL = FALSE;
  HLStatus = FALSE;
  JustHLStatus = FALSE;
  FeatureOverwritten = FALSE;

}
