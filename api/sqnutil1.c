/*   sqnutil1.c
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
* File Name:  sqnutil1.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   9/2/97
*
* $Revision: 6.39 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <sqnutils.h>
#include <gather.h>
#include <subutil.h>
#include <objfdef.h>
#include <seqport.h>
#include <objproj.h>
/* #include <objmmdb1.h> */
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <tofasta.h>
#include <parsegb.h>
#include <utilpars.h>
#include <validatr.h>
#include <explore.h>

NLM_EXTERN DatePtr DateAdvance (DatePtr dp, Uint1 monthsToAdd)

{
  if (dp == NULL) {
    dp = DateCurr ();
  }
  if (dp != NULL && dp->data [0] == 1 && dp->data [1] > 0) {
    while (monthsToAdd > 12) {
      monthsToAdd--;
      (dp->data [1])++;
    }
    if (dp->data [2] < 13 - monthsToAdd) {
      (dp->data [2]) += monthsToAdd;
    } else {
      (dp->data [1])++;
      (dp->data [2]) -= (12 - monthsToAdd);
    }
    if (dp->data [2] == 0) {
      dp->data [2] = 1;
    }
    if (dp->data [3] == 0) {
      switch (dp->data [2]) {
        case 4 :
        case 6 :
        case 9 :
        case 11 :
          dp->data [3] = 30;
          break;
        case 2 :
          dp->data [3] = 28;
          break;
        default :
          dp->data [3] = 31;
          break;
      }
    }
  }
  switch (dp->data [2]) {
    case 4 :
    case 6 :
    case 9 :
    case 11 :
      if (dp->data [3] > 30) {
        dp->data [3] = 30;
      }
      break;
    case 2 :
      if (dp->data [3] > 28) {
        dp->data [3] = 28;
      }
      break;
    default :
      if (dp->data [3] > 31) {
        dp->data [3] = 31;
      }
      break;
  }
  return dp;
}

typedef struct orgscan {
  ObjMgrPtr     omp;
  Int2          nuclCode;
  Int2          mitoCode;
  Boolean       mito;
  Char          taxname [64];
  BioSourcePtr  biop;
} OrgScan, PNTR OrgScanPtr;

static Boolean OrgScanGatherFunc (GatherContextPtr gcp)

{
  BioSourcePtr   biop;
  Boolean        doCodes = FALSE;
  Boolean        doMito = FALSE;
  Boolean        doTaxname = FALSE;
  Boolean        mito = FALSE;
  Int2           mitoCode = 0;
  Int2           nuclCode = 0;
  ObjMgrTypePtr  omtp;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  OrgScanPtr     osp;
  ValNodePtr     sdp;
  SeqFeatPtr     sfp;
  Uint2          subtype;
  CharPtr        taxname = NULL;
  Int2           val;
  ValNodePtr     vnp;

  if (gcp == NULL || gcp->thisitem == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT  && gcp->thistype != OBJ_SEQDESC) return TRUE;

  osp = (OrgScanPtr) gcp->userdata;
  if (osp == NULL) return TRUE;

  subtype = 0;
  if (gcp->thistype == OBJ_SEQFEAT  || gcp->thistype == OBJ_SEQDESC) {
    omtp = ObjMgrTypeFind (osp->omp, gcp->thistype, NULL, NULL);
    if (omtp == NULL) {
      return TRUE;
    }
    if (omtp->subtypefunc != NULL) {
      subtype = (*(omtp->subtypefunc)) (gcp->thisitem);
    }
  }

  orp = NULL;
  biop = NULL;
  switch (gcp->thistype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      switch (subtype) {
        case FEATDEF_ORG :
          orp = (OrgRefPtr) sfp->data.value.ptrvalue;
          break;
        case FEATDEF_BIOSRC :
          biop = (BioSourcePtr) sfp->data.value.ptrvalue;
          break;
        default :
          break;
      }
      break;
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) gcp->thisitem;
      switch (subtype) {
        case Seq_descr_modif :
          vnp = (ValNodePtr) sdp->data.ptrvalue;
          while (vnp != NULL) {
            val = (Int2) vnp->data.intvalue;
            if (val == MODIF_mitochondrial || val == MODIF_kinetoplast) {
              mito = TRUE;
              doMito = TRUE;
              /* osp->mito = TRUE; */
            }
            vnp = vnp->next;
          }
          break;
        case Seq_descr_org :
          orp = (OrgRefPtr) sdp->data.ptrvalue;
          break;
        case Seq_descr_source :
          biop = (BioSourcePtr) sdp->data.ptrvalue;
          break;
        default :
          break;
      }
      break;
    default :
      break;
  }

  if (orp == NULL && biop != NULL) {
    orp = biop->org;
    mito = (Boolean) (biop->genome == 4 || biop->genome == 5);
    doMito = TRUE;
    /* osp->mito = (Boolean) (biop->genome == 4 || biop->genome == 5); */
  }
  if (orp != NULL) {
    taxname = orp->taxname;
    doTaxname = TRUE;
    /* StringNCpy_0 (osp->taxname, orp->taxname, sizeof (osp->taxname)); */
    onp = orp->orgname;
    if (onp != NULL) {
      nuclCode = onp->gcode;
      mitoCode = onp->mgcode;
      doCodes = TRUE;
      /* osp->nuclCode = onp->gcode;
      osp->mitoCode = onp->mgcode; */
    }
  }
  if (biop != NULL) {
    if (osp->biop == NULL || biop->is_focus) {
      osp->biop = biop;
      if (doMito) {
        osp->mito = mito;
      }
      if (doCodes) {
        osp->nuclCode = nuclCode;
        osp->mitoCode = mitoCode;
      }
      if (doTaxname) {
        StringNCpy_0 (osp->taxname, taxname, sizeof (osp->taxname));
      }
    }
  }

  return TRUE;
}

static Int2 SeqEntryOrEntityIDToGeneticCode (SeqEntryPtr sep, Uint2 entityID, BoolPtr mito,
                                             CharPtr taxname, size_t maxsize,
                                             BioSourcePtr PNTR biopp)

{
  GatherScope  gs;
  OrgScan      osp;

  if (mito != NULL) {
    *mito = FALSE;
  }
  if (taxname != NULL && maxsize > 0) {
    *taxname = '\0';
  }
  osp.mito = FALSE;
  osp.nuclCode = 0;
  osp.mitoCode = 0;
  osp.omp = ObjMgrGet ();
  osp.taxname [0] = '\0';
  osp.biop = NULL;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  if (sep != NULL) {
    gs.scope = sep;
    GatherSeqEntry (sep, (Pointer) &osp, OrgScanGatherFunc, &gs);
  } else if (entityID > 0) {
    GatherEntity (entityID, (Pointer) &osp, OrgScanGatherFunc, &gs);
  }
  if (mito != NULL) {
    *mito = osp.mito;
  }
  if (taxname != NULL && maxsize > 0) {
    StringNCpy_0 (taxname, osp.taxname, maxsize);
  }
  if (biopp != NULL) {
    *biopp = osp.biop;
  }
  if (osp.mito) {
    return osp.mitoCode;
  } else {
    return osp.nuclCode;
  }
}

NLM_EXTERN Int2 EntityIDToGeneticCode (Uint2 entityID, BoolPtr mito, CharPtr taxname, size_t maxsize)

{
  return SeqEntryOrEntityIDToGeneticCode (NULL, entityID, mito, taxname, maxsize, NULL);
}

NLM_EXTERN Int2 SeqEntryToGeneticCode (SeqEntryPtr sep, BoolPtr mito, CharPtr taxname, size_t maxsize)

{
  return SeqEntryOrEntityIDToGeneticCode (sep, 0, mito, taxname, maxsize, NULL);
}

NLM_EXTERN Int2 SeqEntryToBioSource (SeqEntryPtr sep, BoolPtr mito, CharPtr taxname, size_t maxsize, BioSourcePtr PNTR biopp)

{
  return SeqEntryOrEntityIDToGeneticCode (sep, 0, mito, taxname, maxsize, biopp);
}

static Boolean FindBspItem (GatherContextPtr gcp)

{
  BioseqPtr  PNTR bspp;

  bspp = (BioseqPtr PNTR) gcp->userdata;
  if (bspp != NULL && gcp->thistype == OBJ_BIOSEQ) {
    *bspp = (BioseqPtr) gcp->thisitem;
  }
  return TRUE;
}

NLM_EXTERN BioseqPtr GetBioseqGivenIDs (Uint2 entityID, Uint2 itemID, Uint2 itemtype)

{
  BioseqPtr  bsp;

  bsp = NULL;
  if (entityID > 0 && itemID > 0 && itemtype == OBJ_BIOSEQ) {
    GatherItem (entityID, itemID, itemtype, (Pointer) (&bsp), FindBspItem);
  }
  return bsp;
}

NLM_EXTERN BioseqPtr GetBioseqGivenSeqLoc (SeqLocPtr slp, Uint2 entityID)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;

  if (slp == NULL) return NULL;
  bsp = NULL;
  sip = SeqLocId (slp);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
  } else if (entityID > 0) {
    slp = SeqLocFindNext (slp, NULL);
    if (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          sep = GetBestTopParentForData (entityID, bsp);
          if (sep != NULL) {
            sep = FindNucSeqEntry (sep);
            if (sep != NULL && sep->choice == 1) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
            }
          }
        }
      }
    }
  }
  return bsp;
}

typedef struct tripletdata {
	Uint2      entityID;
	Uint2      itemID;
	Uint2      itemtype;
	Pointer    lookfor;
} TripletData, PNTR TripletDataPtr;

static Boolean FindIDsFromPointer (GatherContextPtr gcp)

{
  TripletDataPtr  tdp;

  tdp = (TripletDataPtr) gcp->userdata;
  if (tdp != NULL && gcp->thisitem == tdp->lookfor) {
    tdp->entityID = gcp->entityID;
    tdp->itemID = gcp->itemID;
    tdp->itemtype = gcp->thistype;
  }
  return TRUE;
}

NLM_EXTERN Uint2 GetItemIDGivenPointer (Uint2 entityID, Uint2 itemtype, Pointer lookfor)

{
  GatherScope  gs;
  TripletData  td;

  if (entityID > 0 && itemtype > 0 && itemtype < OBJ_MAX && lookfor != NULL) {
    td.entityID = 0;
    td.itemID = 0;
    td.itemtype = 0;
    td.lookfor = lookfor;
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    gs.get_feats_location = FALSE;
    MemSet ((Pointer)(gs.ignore), (int)(FALSE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    /* gs.ignore[itemtype] = FALSE; */
    GatherEntity (entityID, (Pointer) (&td), FindIDsFromPointer, &gs);
    if (td.entityID == entityID && td.itemID > 0 && td.itemtype == itemtype) {
      return td.itemID;
    }
  }
  return 0;
}

static void AddNucPart (BioseqPtr segseq, BioseqSetPtr parts, SeqEntryPtr addme)

{
  BioseqPtr    bsp;
  SeqLocPtr    slp;
  SeqEntryPtr  tmp;

  if (segseq == NULL || addme == NULL) return;
  if (addme->choice != 1 || addme->data.ptrvalue == NULL) return;
  bsp = (BioseqPtr) addme->data.ptrvalue;

  slp = ValNodeNew ((ValNodePtr) segseq->seq_ext);
  if (slp == NULL) return;
  if (segseq->seq_ext == NULL) {
    segseq->seq_ext = (Pointer) slp;
  }
  if (bsp->length >= 0) {
    segseq->length += bsp->length;
    slp->choice = SEQLOC_WHOLE;
    slp->data.ptrvalue = (Pointer) SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
  } else {
    slp->choice = SEQLOC_NULL;
    addme = SeqEntryFree (addme);
    return;
  }

  if (parts == NULL) {
    addme = SeqEntryFree (addme);
    return;
  }
  if (parts->seq_set != NULL) {
    tmp = parts->seq_set;
    while (tmp->next != NULL) {
      tmp = tmp->next;
    }
    tmp->next = addme;
  } else {
    parts->seq_set = addme;
  }
}

NLM_EXTERN void GetSeqEntryParent (SeqEntryPtr target, Pointer PNTR parentptr, Uint2Ptr parenttype)

{
  ObjMgrPtr      omp;
  ObjMgrDataPtr  omdp;

  if (parentptr == NULL || parenttype == NULL) return;
  *parenttype = 0;
  *parentptr = NULL;
  if (target == NULL || target->data.ptrvalue == NULL) return;
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  omdp = ObjMgrFindByData (omp, target->data.ptrvalue);
  if (omdp == NULL) return;
  *parenttype = omdp->parenttype;
  *parentptr = omdp->parentptr;
}

NLM_EXTERN void SaveSeqEntryObjMgrData (SeqEntryPtr target, ObjMgrDataPtr PNTR omdptopptr, ObjMgrData PNTR omdataptr)

{
  ObjMgrPtr	     omp;
  ObjMgrDataPtr  omdp, omdptop = NULL;

  if (target == NULL || omdptopptr == NULL || omdataptr == NULL) return;
  *omdptopptr = NULL;
  MemSet ((Pointer) omdataptr, 0, sizeof (ObjMgrData));
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  omdp = ObjMgrFindByData (omp, target->data.ptrvalue);
  if (omdp == NULL) return;
  omdptop = ObjMgrFindTop (omp, omdp);
  if (omdptop == NULL) return;
  if (omdptop->EntityID == 0) return;
  *omdptopptr = omdptop;
  MemCopy ((Pointer) omdataptr, omdptop, sizeof (ObjMgrData));
  omdptop->userdata = NULL;
}

extern void ObjMgrRemoveEntityIDFromRecycle (Uint2 entityID, ObjMgrPtr omp);
NLM_EXTERN void RestoreSeqEntryObjMgrData (SeqEntryPtr target, ObjMgrDataPtr omdptop, ObjMgrData PNTR omdataptr)

{
  ObjMgrPtr	omp;
  ObjMgrDataPtr omdp, omdpnew = NULL;

  if (target == NULL || omdptop == NULL || omdataptr == NULL) return;
  if (omdataptr->EntityID == 0) return;
  omp = ObjMgrGet ();
  if (omp == NULL) return;
  omdp = ObjMgrFindByData (omp, target->data.ptrvalue);
  if (omdp == NULL) return;
  omdpnew = ObjMgrFindTop (omp, omdp);
  if (omdpnew == NULL) return;
  if (omdpnew != omdptop) {
    omdpnew->EntityID = omdataptr->EntityID;
    omdptop->EntityID = 0;
    omdpnew->lockcnt = omdataptr->lockcnt;
    omdpnew->tempload = omdataptr->tempload;
    omdpnew->clipboard = omdataptr->clipboard;
    omdpnew->dirty = omdataptr->dirty;
    omdpnew->being_freed = omdataptr->being_freed;
    omdpnew->free = omdataptr->free;
    omdpnew->options = omdataptr->options;
    ObjMgrRemoveEntityIDFromRecycle (omdpnew->EntityID, omp);
  }
  omdpnew->userdata = omdataptr->userdata;
}

NLM_EXTERN void AddSeqEntryToSeqEntry (SeqEntryPtr target, SeqEntryPtr insert, Boolean relink)

{
  SeqEntryPtr    first;
  BioseqPtr      insertbsp;
  BioseqSetPtr   nuc_prot;
  Uint2          parenttype;
  Pointer        parentptr;
  BioseqSetPtr   parts;
  BioseqPtr      seg;
  BioseqSetPtr   segs;
  BioseqPtr      targetbsp;
  BioseqSetPtr   targetbssp;
  SeqEntryPtr    the_nuc;
  SeqEntryPtr    tmp;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;

  if (target == NULL || insert == NULL) return;
  if (target->data.ptrvalue == NULL || insert->data.ptrvalue == NULL) return;

  if (relink) {
    SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
    GetSeqEntryParent (target, &parentptr, &parenttype);
  }

  if (IS_Bioseq (target) && IS_Bioseq (insert)) {
    targetbsp = (BioseqPtr) target->data.ptrvalue;
    insertbsp = (BioseqPtr) insert->data.ptrvalue;
    if (ISA_na (targetbsp->mol)) {
      if (ISA_na (insertbsp->mol)) {

        seg = BioseqNew ();
        if (seg == NULL) return;
        seg->mol = targetbsp->mol;
        seg->repr = Seq_repr_seg;
        seg->seq_ext_type = 1;
        seg->length = 0;
        /* seg->id = MakeSeqID ("SEG_dna"); */
        /* seg->id = MakeNewProteinSeqId (NULL, NULL); */
        seg->id = MakeUniqueSeqID ("segseq_");
        SeqMgrAddToBioseqIndex (seg);

        the_nuc = SeqEntryNew ();
        if (the_nuc == NULL) return;
        the_nuc->choice = 1;
        the_nuc->data.ptrvalue = (Pointer) seg;

        segs = BioseqSetNew ();
        if (segs == NULL) return;
        segs->_class = 2;
        segs->seq_set = the_nuc;

        parts = BioseqSetNew ();
        if (parts == NULL) return;
        parts->_class = 4;

        tmp = SeqEntryNew ();
        if (tmp == NULL) return;
        tmp->choice = 2;
        tmp->data.ptrvalue = (Pointer) parts;
        the_nuc->next = tmp;

        first = SeqEntryNew ();
        if (first == NULL) return;
        first->choice = 1;
        first->data.ptrvalue = (Pointer) targetbsp;
        target->choice = 2;
        target->data.ptrvalue = (Pointer) segs;

        AddNucPart (seg, parts, first);
        AddNucPart (seg, parts, insert);

      } else if (ISA_aa (insertbsp->mol)) {

        nuc_prot = BioseqSetNew ();
        if (nuc_prot == NULL) return;
        nuc_prot->_class = 1;

        the_nuc = SeqEntryNew ();
        if (the_nuc == NULL) return;
        the_nuc->choice = 1;
        the_nuc->data.ptrvalue = (Pointer) targetbsp;
        target->choice = 2;
        target->data.ptrvalue = (Pointer) nuc_prot;
        nuc_prot->seq_set = the_nuc;

        the_nuc->next = insert;

      }
    }
  } else if (IS_Bioseq_set (target)) {
    targetbssp = (BioseqSetPtr) target->data.ptrvalue;
    if (targetbssp->_class == 1 && IS_Bioseq (insert)) {
     insertbsp = (BioseqPtr) insert->data.ptrvalue;
     if (ISA_aa (insertbsp->mol)) {

        nuc_prot = targetbssp;
        if (nuc_prot->seq_set != NULL) {
          tmp = nuc_prot->seq_set;
          while (tmp->next != NULL) {
            tmp = tmp->next;
          }
          tmp->next = insert;
        } else {
          nuc_prot->seq_set = insert;
        }

      }
    } else if (targetbssp->_class == 2 && IS_Bioseq (insert)) {
      insertbsp = (BioseqPtr) insert->data.ptrvalue;
      if (ISA_na (insertbsp->mol)) {

        the_nuc = FindNucSeqEntry (target);
        if (the_nuc != NULL && the_nuc->next != NULL) {
          tmp = the_nuc->next;
          if (tmp->choice == 2 && tmp->data.ptrvalue != NULL) {
            parts = (BioseqSetPtr) tmp->data.ptrvalue;
            if (parts->_class == 4 && the_nuc->choice == 1) {
              seg = (BioseqPtr) the_nuc->data.ptrvalue;
              AddNucPart (seg, parts, insert);
            }
          }
        }

      } else if (ISA_aa (insertbsp->mol)) {

        nuc_prot = BioseqSetNew ();
        if (nuc_prot == NULL) return;
        nuc_prot->_class = 1;

        first = SeqEntryNew ();
        if (first == NULL) return;
        first->choice = 2;
        first->data.ptrvalue = (Pointer) targetbssp;
        target->choice = 2;
        target->data.ptrvalue = (Pointer) nuc_prot;
        nuc_prot->seq_set = first;

        first->next = insert;

      }
    } else if (targetbssp->_class == 7) {

      if (targetbssp->seq_set != NULL) {
        tmp = targetbssp->seq_set;
        while (tmp->next != NULL) {
          tmp = tmp->next;
        }
        tmp->next = insert;
      } else {
        targetbssp->seq_set = insert;
      }
    } else if (targetbssp->_class >= 13 && targetbssp->_class <= 15) {

      if (targetbssp->seq_set != NULL) {
        tmp = targetbssp->seq_set;
        while (tmp->next != NULL) {
          tmp = tmp->next;
        }
        tmp->next = insert;
      } else {
        targetbssp->seq_set = insert;
      }

    } else if (targetbssp->_class == BioseqseqSet_class_gen_prod_set) {

      if (targetbssp->seq_set != NULL) {
        tmp = targetbssp->seq_set;
        while (tmp->next != NULL) {
          tmp = tmp->next;
        }
        tmp->next = insert;
      } else {
        targetbssp->seq_set = insert;
      }

    }
  }

  if (relink) {
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }
}

NLM_EXTERN void ReplaceSeqEntryWithSeqEntry (SeqEntryPtr target, SeqEntryPtr replaceWith, Boolean relink)

{
  Uint2          parenttype;
  Pointer        parentptr;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;

  if (target == NULL || replaceWith == NULL) return;

  if (relink) {
    SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
    GetSeqEntryParent (target, &parentptr, &parenttype);
  }

  if (target->choice == 1) {
    BioseqFree ((BioseqPtr) target->data.ptrvalue);
  } else if (target->choice == 2) {
    BioseqSetFree ((BioseqSetPtr) target->data.ptrvalue);
  }
  target->choice = replaceWith->choice;
  target->data.ptrvalue = replaceWith->data.ptrvalue;
  MemFree (replaceWith);

  if (relink) {
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }
}

static void SeqEntryRemoveLoop (SeqEntryPtr sep, SeqEntryPtr del, SeqEntryPtr PNTR prev)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   next;

  while (sep != NULL) {
    next = sep->next;
    if (sep == del) {
      *prev = sep->next;
      sep->next = NULL;
      SeqEntryFree (sep);
    } else {
      prev = (SeqEntryPtr PNTR) &(sep->next);
      if (IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL) {
          SeqEntryRemoveLoop (bssp->seq_set, del, &(bssp->seq_set));
        }
      }
    }
    sep = next;
  }
}

NLM_EXTERN void RemoveSeqEntryFromSeqEntry (SeqEntryPtr top, SeqEntryPtr del, Boolean relink)

{
  SeqEntryPtr    dummy;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;

  if (top == NULL || del == NULL) return;
  if (top->data.ptrvalue == NULL || del->data.ptrvalue == NULL) return;

  if (relink) {
    SaveSeqEntryObjMgrData (top, &omdptop, &omdata);
    GetSeqEntryParent (top, &parentptr, &parenttype);
  }

  dummy = NULL;
  SeqEntryRemoveLoop (top, del, &dummy);

  if (relink) {
    SeqMgrLinkSeqEntry (top, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (top, omdptop, &omdata);
  }
}

NLM_EXTERN void DeleteMultipleTitles (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Boolean       hastitle;
  ValNodePtr    nextsdp;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  hastitle = FALSE;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_title) {
      if (hastitle) {
        *(prevsdp) = sdp->next;
        sdp->next = NULL;
        SeqDescFree (sdp);
      } else {
        hastitle = TRUE;
        prevsdp = (Pointer PNTR) &(sdp->next);
      }
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

NLM_EXTERN void RenormalizeNucProtSets (SeqEntryPtr sep, Boolean relink)

{
  SeqAnnotPtr    annot;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     descr;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqAnnotPtr    sap;
  SeqEntryPtr    seqentry;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        RenormalizeNucProtSets (sep, relink);
      }
      return;
    }
    if (bssp != NULL && bssp->_class == 1) {
      seqentry = bssp->seq_set;
      if (seqentry != NULL && seqentry->next == NULL) {

        if (relink) {
          SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
          GetSeqEntryParent (sep, &parentptr, &parenttype);
        }

        descr = bssp->descr;
        bssp->descr = NULL;
        annot = bssp->annot;
        bssp->annot = NULL;

        sep->choice = seqentry->choice;
        sep->data.ptrvalue = seqentry->data.ptrvalue;
        seqentry->data.ptrvalue = NULL;
        bssp->seq_set = NULL;
        SeqEntryFree (seqentry);
        BioseqSetFree (bssp);

        sap = NULL;
        if (IS_Bioseq (sep)) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          ValNodeLink (&(bsp->descr), descr);
          if (bsp->annot == NULL) {
            bsp->annot = annot;
            annot = NULL;
          } else {
            sap = bsp->annot;
          }
        } else if (IS_Bioseq_set (sep)) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          ValNodeLink (&(bssp->descr), descr);
          if (bssp->annot == NULL) {
            bssp->annot = annot;
            annot = NULL;
          } else {
            sap = bssp->annot;
          }
        }
        if (sap != NULL) {
          while (sap->next != NULL) {
            sap = sap->next;
          }
          sap->next = annot;
        }
        DeleteMultipleTitles (sep, NULL, 0, 0);

        if (relink) {
          SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
          RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
        }
      }
    }
  }
}

NLM_EXTERN ValNodePtr ExtractBioSourceAndPubs (SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    descr;
  ValNodePtr    last;
  ValNodePtr    nextsdp;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return NULL;
  descr = NULL;
  last = NULL;
  sdp = NULL;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return NULL;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_pub || sdp->choice == Seq_descr_source) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      if (descr == NULL) {
        descr = sdp;
        last = descr;
      } else if (last != NULL) {
        last->next = sdp;
        last = last->next;
      }
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
  return descr;
}

NLM_EXTERN void ReplaceBioSourceAndPubs (SeqEntryPtr sep, ValNodePtr descr)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    last;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;

  if (sep == NULL || descr == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  last = descr;
  while (last->next != NULL) {
    last = last->next;
  }
  last->next = sdp;
  *(prevsdp) = descr;
}

typedef struct targetdata {
  BioseqPtr    bsp;
  SeqEntryPtr  nps;
  Boolean      skipGenProdSet;
} TargetData, PNTR TargetDataPtr;

static Boolean ReturnStackToItem (GatherContextPtr gcp)

{
  BioseqSetPtr   bssp;
  Int2           i;
  Uint2          itemtype;
  TargetDataPtr  tdp;

  if (gcp == NULL) return TRUE;
  tdp = (TargetDataPtr) gcp->userdata;
  if (tdp == NULL) return TRUE;
  if (gcp->gatherstack != NULL && gcp->numstack > 0) {
    for (i = 0; i < gcp->numstack; i++) {
      itemtype = gcp->gatherstack [i].itemtype;
      if (itemtype == OBJ_BIOSEQ || itemtype == OBJ_BIOSEQSET) {
        tdp->nps = SeqMgrGetSeqEntryForData (gcp->gatherstack [i].thisitem);
        if (gcp->gatherstack [i].itemtype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) gcp->gatherstack [i].thisitem;
          if (bssp->_class != 7 && bssp->_class != 13 &&
              bssp->_class != 14 && bssp->_class != 15 &&
              (bssp->_class != BioseqseqSet_class_gen_prod_set || (! tdp->skipGenProdSet))) {
            return FALSE;
          }
        } else if (gcp->gatherstack [i].itemtype == OBJ_BIOSEQ) {
          return FALSE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean GetStackToTarget (GatherContextPtr gcp)

{
  TargetDataPtr  tdp;

  if (gcp == NULL) return TRUE;
  tdp = (TargetDataPtr) gcp->userdata;
  if (tdp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ) {
    if (tdp->bsp == (BioseqPtr) gcp->thisitem) {
      return ReturnStackToItem (gcp);
    }
  }
  return TRUE;
}

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForDataEx (Uint2 entityID, BioseqPtr bsp, Boolean skipGenProdSet)

{
  GatherScope  gs;
  TargetData   td;

  td.bsp = bsp;
  td.nps = NULL;
  td.skipGenProdSet = skipGenProdSet;
  if (entityID > 0 && bsp != NULL) {
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    MemSet ((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    GatherEntity (entityID, (Pointer) &td, GetStackToTarget, &gs);
  }
  return td.nps;
}

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForData (Uint2 entityID, BioseqPtr bsp)

{
  return GetBestTopParentForDataEx (entityID, bsp, FALSE);
}

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForItemIDEx (Uint2 entityID, Uint2 itemID, Uint2 itemtype, Boolean skipGenProdSet)

{
  TargetData  td;

  td.bsp = NULL;
  td.nps = NULL;
  td.skipGenProdSet = skipGenProdSet;
  if (entityID > 0 && itemID > 0 && itemtype > 0) {
    GatherItem (entityID, itemID, itemtype, (Pointer) &td, ReturnStackToItem);
  }
  return td.nps;
}

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForItemID (Uint2 entityID, Uint2 itemID, Uint2 itemtype)

{
  return GetBestTopParentForItemIDEx (entityID, itemID, itemtype, FALSE);
}

NLM_EXTERN SeqEntryPtr LIBCALL GetTopSeqEntryForEntityID (Uint2 entityID)

{
  ObjMgrDataPtr  omdp;
  SeqSubmitPtr   ssp;

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL) {
    switch (omdp->datatype) {
      case OBJ_SEQSUB :
        ssp = (SeqSubmitPtr) omdp->dataptr;
        if (ssp != NULL && ssp->datatype == 1) {
          return (SeqEntryPtr) ssp->data;
        }
        break;
      case OBJ_BIOSEQ :
        return (SeqEntryPtr) omdp->choice;
      case OBJ_BIOSEQSET :
        return (SeqEntryPtr) omdp->choice;
      default :
        break;
    }
  }
  return NULL;
}

NLM_EXTERN Boolean CheckSeqLocForPartial (SeqLocPtr location, BoolPtr p5ptr, BoolPtr p3ptr)

{
  SeqLocPtr   firstSlp;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  Boolean     partial5;
  Boolean     partial3;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  SeqPntPtr   spp;

  partial5 = FALSE;
  partial3 = FALSE;
  if (location != NULL) {
    firstSlp = NULL;
    lastSlp = NULL;
    slp = SeqLocFindNext (location, NULL);
    while (slp != NULL) {
      if (firstSlp == NULL) {
        firstSlp = slp;
      }
      lastSlp = slp;
      slp = SeqLocFindNext (location, slp);
    }
    if (firstSlp != NULL) {
      if (firstSlp->choice == SEQLOC_INT && firstSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) firstSlp->data.ptrvalue;
        if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
          ifp = sip->if_to;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial5 = TRUE;
          }
        } else {
          ifp = sip->if_from;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial5 = TRUE;
          }
        }
      } else if (firstSlp->choice == SEQLOC_PNT && firstSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) firstSlp->data.ptrvalue;
        if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial5 = TRUE;
          }
        } else {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial5 = TRUE;
          }
        }
      }
    }
    if (lastSlp != NULL) {
      if (lastSlp->choice == SEQLOC_INT && lastSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) lastSlp->data.ptrvalue;
        if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
          ifp = sip->if_from;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial3 = TRUE;
          }
        } else {
          ifp = sip->if_to;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial3 = TRUE;
          }
        }
      } else if (lastSlp->choice == SEQLOC_PNT && lastSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) lastSlp->data.ptrvalue;
        if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 2) {
            partial3 = TRUE;
          }
        } else {
          ifp = spp->fuzz;
          if (ifp != NULL && ifp->choice == 4 && ifp->a == 1) {
            partial3 = TRUE;
          }
        }
      }
    }
  }
  if (p5ptr != NULL) {
    *p5ptr = partial5;
  }
  if (p3ptr != NULL) {
    *p3ptr = partial3;
  }
  return (Boolean) (partial5 || partial3);
}

NLM_EXTERN void SetSeqLocPartial (SeqLocPtr location, Boolean partial5, Boolean partial3)

{
  SeqLocPtr   firstSlp;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  SeqPntPtr   spp;

  if (location != NULL) {
    firstSlp = NULL;
    lastSlp = NULL;
    slp = SeqLocFindNext (location, NULL);
    while (slp != NULL) {
      if (firstSlp == NULL) {
        firstSlp = slp;
      }
      lastSlp = slp;
      slp = SeqLocFindNext (location, slp);
    }
    if (firstSlp != NULL) {
      if (firstSlp->choice == SEQLOC_INT && firstSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) firstSlp->data.ptrvalue;
        if (partial5) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
              sip->if_to = ifp;
              ifp->a = 1;
            } else {
              sip->if_from = ifp;
              ifp->a = 2;
            }
          }
        } else {
          if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
            sip->if_to = IntFuzzFree (sip->if_to);
          } else {
            sip->if_from = IntFuzzFree (sip->if_from);
          }
        }
      } else if (firstSlp->choice == SEQLOC_PNT && firstSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) firstSlp->data.ptrvalue;
        if (partial5) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
              spp->fuzz = ifp;
              ifp->a = 1;
            } else {
              spp->fuzz = ifp;
              ifp->a = 2;
            }
          }
        } else {
          if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
            spp->fuzz = IntFuzzFree (spp->fuzz);
          } else {
            spp->fuzz = IntFuzzFree (spp->fuzz);
          }
        }
      }
    }
    if (lastSlp != NULL) {
      if (lastSlp->choice == SEQLOC_INT && lastSlp->data.ptrvalue != NULL) {
        sip = (SeqIntPtr) lastSlp->data.ptrvalue;
        if (partial3) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
              sip->if_from = ifp;
              ifp->a = 2;
            } else {
              sip->if_to = ifp;
              ifp->a = 1;
            }
          }
        } else {
          if (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev) {
            sip->if_from = IntFuzzFree (sip->if_from);
          } else {
            sip->if_to = IntFuzzFree (sip->if_to);
          }
        }
      } else if (lastSlp->choice == SEQLOC_PNT && lastSlp->data.ptrvalue != NULL) {
        spp = (SeqPntPtr) lastSlp->data.ptrvalue;
        if (partial3) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
              spp->fuzz = ifp;
              ifp->a = 2;
            } else {
              spp->fuzz = ifp;
              ifp->a = 1;
            }
          }
        } else {
          if (spp->strand == Seq_strand_minus || spp->strand == Seq_strand_both_rev) {
            spp->fuzz = IntFuzzFree (spp->fuzz);
          } else {
            spp->fuzz = IntFuzzFree (spp->fuzz);
          }
        }
      }
    }
  }
}

/* begin PromoteXrefs section */

typedef struct geneextendlist {
  GeneRefPtr  grp;
  SeqLocPtr   slp;
  ObjMgrPtr   omp;
  Boolean     rsult;
  Char        label [41];
} GeneExtendList, PNTR GeneExtendPtr;

static Boolean GeneExtendFunc (GatherContextPtr gcp)

{
  BioseqPtr      bsp;
  GeneExtendPtr  gep;
  GeneRefPtr     grp;
  Boolean        hasNulls;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  SeqLocPtr      slp;
  Char           thislabel [41];

  if (gcp == NULL) return TRUE;

  gep = (GeneExtendPtr) gcp->userdata;
  if (gep == NULL ) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      omtp = ObjMgrTypeFind (gep->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        if (StringICmp (thislabel, gep->label) == 0) {
          if (SeqLocCompare (gep->slp, sfp->location) != SLC_NO_MATCH) {
            bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
            if (bsp != NULL) {
              slp = SeqLocMerge (bsp, sfp->location, gep->slp, TRUE, FALSE, FALSE);
              if (slp != NULL) {
                sfp->location = SeqLocFree (sfp->location);
                sfp->location = slp;
                if (bsp->repr == Seq_repr_seg) {
                  slp = SegLocToPartsEx (bsp, sfp->location, TRUE);
                  sfp->location = SeqLocFree (sfp->location);
                  sfp->location = slp;
                  hasNulls = LocationHasNullsBetween (sfp->location);
                  sfp->partial = (sfp->partial || hasNulls);
                }
                FreeAllFuzz (slp);
                gep->rsult = TRUE;
              }
            }
          }
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

static Boolean ExtendGene (GeneRefPtr grp, SeqEntryPtr nsep, SeqLocPtr slp)

{
  GeneExtendList  gel;
  GatherScope     gs;
  ObjMgrTypePtr   omtp;
  SeqFeatPtr      sfp;

  if (grp == NULL || nsep == NULL || slp == NULL) return FALSE;
  gel.grp = grp;
  gel.slp = slp;
  gel.omp = ObjMgrGet ();
  gel.label [0] = '\0';
  gel.rsult = FALSE;
  omtp = ObjMgrTypeFind (gel.omp, OBJ_SEQFEAT, NULL, NULL);
  if (omtp != NULL && omtp->labelfunc != NULL) {
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sfp->data.choice = SEQFEAT_GENE;
      sfp->data.value.ptrvalue = (Pointer) grp;
      (*(omtp->labelfunc)) ((Pointer) sfp, gel.label, 40, OM_LABEL_CONTENT);
      sfp->data.value.ptrvalue = NULL;
      SeqFeatFree (sfp);
    }
  }
  MemSet ((Pointer)(&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherSeqEntry (nsep, (Pointer) &gel, GeneExtendFunc, &gs);
  return gel.rsult;
}

NLM_EXTERN void SetEmptyGeneticCodes (SeqAnnotPtr sap, Int2 genCode)

{
  CdRegionPtr     crp;
  GeneticCodePtr  gc;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp;

  if (sap == NULL || sap->type != 1) return;
  for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_CDREGION) {
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        gc = crp->genetic_code;
        if (gc != NULL) {
          vnp = gc->data.ptrvalue;
          if (vnp != NULL && vnp->choice == 2) {
            vnp->data.intvalue = (Int4) genCode;
            /*
            if (vnp->data.intvalue == 0) {
              vnp->data.intvalue = (Int4) genCode;
            }
            */
          }
        }
      }
    }
  }
}

NLM_EXTERN void PromoteXrefs (SeqFeatPtr sfp, BioseqPtr bsp, Uint2 entityID)

{
  ByteStorePtr         bs;
  Char                 ch;
  CdRegionPtr          crp;
  Int2                 ctr = 1;
  ValNodePtr           descr;
  GBQualPtr            gbq;
  SeqFeatPtr           gene;
  GeneRefPtr           grp;
  Int2                 i;
  Char                 id [64];
  SeqEntryPtr          last = NULL;
  MolInfoPtr           mip;
  SeqFeatXrefPtr       next;
  GBQualPtr            nextqual;
  SeqEntryPtr          old;
  ObjMgrDataPtr        omdptop;
  ObjMgrData           omdata;
  Uint2                parenttype;
  Pointer              parentptr;
  Boolean              partial5;
  Boolean              partial3;
  BioseqPtr            pbsp;
  SeqFeatXrefPtr PNTR  prev;
  GBQualPtr PNTR       prevqual;
  SeqFeatPtr           prot;
  CharPtr              protseq;
  ProtRefPtr           prp;
  SeqEntryPtr          psep;
  CharPtr              ptr;
  SeqEntryPtr          sep;
  SeqIdPtr             sip = NULL;
  SeqEntryPtr          target = NULL;
  Uint4                version = 0;
  long int             val;
  ValNodePtr           vnp;
  SeqFeatXrefPtr       xref;

  if (sfp == NULL || bsp == NULL) return;
  while (sfp != NULL) {
    prev = &(sfp->xref);
    xref = sfp->xref;
    while (xref != NULL) {
      next = xref->next;
      if (xref->data.choice == SEQFEAT_GENE && sfp->data.choice != SEQFEAT_GENE) {
        grp = (GeneRefPtr) xref->data.value.ptrvalue;
        if (grp != NULL && SeqMgrGeneIsSuppressed (grp)) {
        } else {
          xref->data.value.ptrvalue = NULL;
          if (grp != NULL) {
            sep = SeqMgrGetSeqEntryForData (bsp);
            if (ExtendGene (grp, sep, sfp->location)) {
              GeneRefFree (grp);
            } else {
              gene = CreateNewFeature (sep, NULL, SEQFEAT_GENE, NULL);
              if (gene != NULL) {
                gene->data.value.ptrvalue = (Pointer) grp;
                gene->location = SeqLocFree (gene->location);
                gene->location = AsnIoMemCopy (sfp->location,
                                               (AsnReadFunc) SeqLocAsnRead,
                                               (AsnWriteFunc) SeqLocAsnWrite);
              }
            }
          }
          *(prev) = next;
          xref->next = NULL;
          xref->data.choice = 0;
          SeqFeatXrefFree (xref);
        }
      } else if (xref->data.choice == SEQFEAT_PROT &&
                 sfp->data.choice == SEQFEAT_CDREGION &&
                 sfp->product == NULL) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
        xref->data.value.ptrvalue = NULL;
        if (prp != NULL) {
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
/**
            crp->frame = 0;
**/
            bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
            if (bs != NULL) {
              protseq = BSMerge (bs, NULL);
              bs = BSFree (bs);
              if (protseq != NULL) {
                ptr = protseq;
                ch = *ptr;
                while (ch != '\0') {
                  *ptr = TO_UPPER (ch);
                  ptr++;
                  ch = *ptr;
                }
                i = (Int2) StringLen (protseq);
                if (i > 0 && protseq [i - 1] == '*') {
                  protseq [i - 1] = '\0';
                }
                bs = BSNew (1000);
                if (bs != NULL) {
                  ptr = protseq;
                  /*
                  if (protseq [0] == '-') {
                    ptr++;
                  }
                  */
                  BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
                }
                protseq = MemFree (protseq);
              }
              pbsp = BioseqNew ();
              if (pbsp != NULL) {
                pbsp->repr = Seq_repr_raw;
                pbsp->mol = Seq_mol_aa;
                pbsp->seq_data_type = Seq_code_ncbieaa;
                pbsp->seq_data = bs;
                pbsp->length = BSLen (bs);
                bs = NULL;
                sep = GetBestTopParentForData (entityID, bsp);
                old = SeqEntrySetScope (sep);
                gbq = sfp->qual;
                prevqual = (GBQualPtr PNTR) &(sfp->qual);
                id [0] = '\0';
                while (gbq != NULL) {
                  nextqual = gbq->next;
                  if (StringICmp (gbq->qual, "protein_id") == 0) {
                    *(prevqual) = gbq->next;
                    gbq->next = NULL;
                    StringNCpy_0 (id, gbq->val, sizeof (id));
                    GBQualFree (gbq);
                  } else {
                    prevqual = (GBQualPtr PNTR) &(gbq->next);
                  }
                  gbq = nextqual;
                }
                if (! StringHasNoText (id)) {
                  if (StringChr (id, '|') != NULL) {
                    sip = SeqIdParse (id);
                  } else {
                    ptr = StringChr (id, '.');
                    if (ptr != NULL) {
                      *ptr = '\0';
                      ptr++;
                      if (sscanf (ptr, "%ld", &val) == 1) {
                        version = (Uint4) val;
                      }
                    }
                    sip = SeqIdFromAccession (id, version, NULL);
                  }
                }
                if (sip != NULL) {
                  pbsp->id = sip;
                } else {
                  pbsp->id = MakeNewProteinSeqIdEx (sfp->location, NULL, &ctr);
                }
                CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
                SeqMgrAddToBioseqIndex (pbsp);
                SeqEntrySetScope (old);
                psep = SeqEntryNew ();
                if (psep != NULL) {
                  psep->choice = 1;
                  psep->data.ptrvalue = (Pointer) pbsp;
                  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) pbsp, psep);
                  mip = MolInfoNew ();
                  if (mip != NULL) {
                    mip->biomol = 8;
                    mip->tech = 8;
                    if (partial5 && partial3) {
                      mip->completeness = 5;
                    } else if (partial5) {
                      mip->completeness = 3;
                    } else if (partial3) {
                      mip->completeness = 4;
                    }
                    vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
                    if (vnp != NULL) {
                      vnp->data.ptrvalue = (Pointer) mip;
                    }
                  }
                  if (last == NULL) {

                    /* the first protein may change the set/seq structure,
                    so goes through AddSeqEntryToSeqEntry */

                    descr = ExtractBioSourceAndPubs (sep);
                    AddSeqEntryToSeqEntry (sep, psep, TRUE);
                    ReplaceBioSourceAndPubs (sep, descr);
                    target = sep;
                    SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
                    GetSeqEntryParent (target, &parentptr, &parenttype);
                  } else {
                    last->next = psep;
                    last = psep;
                  }
                  SetSeqFeatProduct (sfp, pbsp);
                  psep = SeqMgrGetSeqEntryForData (pbsp);
                  if (psep != NULL) {
                    last = psep;
                    prot = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
                    if (prot != NULL) {
                      prot->data.value.ptrvalue = (Pointer) prp;
                    }
                  }
                }
              }
            }
          }
        }
        *(prev) = next;
        xref->next = NULL;
        xref->data.choice = 0;
        SeqFeatXrefFree (xref);
      } else {
        prev = &(xref->next);
      }
      xref = next;
    }
    sfp = sfp->next;
  }
  if (target != NULL) {
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }
}

/* begin BasicSeqEntryCleanup section */

static Boolean HasNoText (CharPtr str)

{
  Uchar  ch;	/* to use 8bit characters in multibyte languages */

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static CharPtr TrimSpacesAndSemicolons (CharPtr str)

{
  CharPtr  amp;
  Uchar    ch;	/* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && ch <= ' ') {
      while (ch != '\0' && ch <= ' ') {
        ptr++;
        ch = *ptr;
      }
      while (ch != '\0') {
        *dst = ch;
        dst++;
        ptr++;
        ch = *ptr;
      }
      *dst = '\0';
    }
    amp = NULL;
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '&') {
        amp = ptr;
        dst = NULL;
      } else if (ch == ' ') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
      } else if (ch == ';') {
        if (dst == NULL && amp == NULL) {
          dst = ptr;
        }
      } else {
        dst = NULL;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static void CleanVisString (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesAndSemicolons (*strp);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void CleanDoubleQuote (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    if (ch == '"') {
      *str = '\'';
    }
    str++;
    ch = *str;
  }
}

static Boolean AlreadyInVnpList (ValNodePtr head, ValNodePtr curr)

{
  if (head == NULL || curr == NULL) return FALSE;
  /* since we cannot sort these lists, must check against all previous entries */
  while (head != curr && head != NULL) {
    if (StringICmp (head->data.ptrvalue, curr->data.ptrvalue) == 0) return TRUE;
    head = head->next;
  }
  return FALSE;
}

static void CleanVisStringList (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimSpacesAndSemicolons (vnp->data.ptrvalue);
    if (HasNoText (vnp->data.ptrvalue) || AlreadyInVnpList (*vnpp, vnp)) {
      *prev = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      prev = &(vnp->next);
    }
    vnp = next;
  }
}

static void CleanDoubleQuoteList (ValNodePtr vnp)

{
  while (vnp != NULL) {
    CleanDoubleQuote ((CharPtr) vnp->data.ptrvalue);
    vnp = vnp->next;
  }
}

static Boolean HandledGBQualOnGene (SeqFeatPtr sfp, GBQualPtr gbq)

{
  Int2        choice = 0;
  GeneRefPtr  grp;

  if (StringICmp (gbq->qual, "pseudo") == 0) {
    choice = 1;
  } else if (StringICmp (gbq->qual, "map") == 0) {
    choice = 2;
  } else if (StringICmp (gbq->qual, "allele") == 0) {
    choice = 3;
  }
  if (choice > 0) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp == NULL) return FALSE;
    switch (choice) {
      case 1 :
        grp->pseudo = TRUE;
        break;
      case 2 :
        if (grp->maploc != NULL) return FALSE;
        if (StringHasNoText (gbq->val)) return FALSE;
        grp->maploc = StringSave (gbq->val);
        break;
      case 3 :
        if (grp->allele != NULL) return FALSE;
        if (StringHasNoText (gbq->val)) return FALSE;
        grp->allele = StringSave (gbq->val);
        break;
      default :
        break;
    }
    return TRUE;
  }
  return FALSE;
}

/* code break parser functions from the flatfile parser */

static Uint1 GetQualValueAa (CharPtr qval)

{
   CharPtr  str, eptr, ptr;
   Uint1    aa;

	str = StringStr(qval, "aa:");
	if (str != NULL) {
   		str += 3;
   		while (*str == ' ')
       		++str;
   		for (eptr = str; *eptr != ')' && *eptr != ' ' && *eptr != '\0';  eptr++) continue;
	}

    ptr = TextSave(str, eptr-str);
    aa = ValidAminoAcid(ptr);
    MemFree(ptr);  

    return (aa);
}

static CharPtr SimpleValuePos (CharPtr qval)

{
   CharPtr bptr, eptr;

   if ((bptr = StringStr(qval, "(pos:")) == NULL) {
   		return NULL;
   }
    
   bptr += 5;
   while (*bptr == ' ')
       ++bptr;
   for (eptr = bptr; *eptr != ',' && *eptr != '\0'; eptr++) continue;

   return (TextSave(bptr, eptr-bptr));
}

extern Boolean ParseCodeBreak (SeqFeatPtr sfp, CharPtr val);
extern Boolean ParseCodeBreak (SeqFeatPtr sfp, CharPtr val)

{
  CodeBreakPtr  cbp;
  Choice        cp;
  CdRegionPtr   crp;
  Int4          diff;
  CodeBreakPtr  lastcbp;
  Boolean       locmap;
  int           num_errs;
  CharPtr       pos;
  Boolean       pos_range = FALSE;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  Boolean       sitesmap;
  SeqLocPtr     slp;
  Uint1         strand;
  Int4          temp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;
  if (StringHasNoText (val)) return FALSE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return FALSE;

  /* find SeqId to use */
  sip = SeqLocId (sfp->location);
  if (sip == NULL) {
    slp = SeqLocFindNext (sfp->location, NULL);
    if (slp != NULL) {
      sip = SeqLocId (slp);
    }
  }
  if (sip == NULL) return FALSE;

  cp.choice = 1; /* ncbieaa */
  cp.value.intvalue = (Int4) GetQualValueAa (val);
  cbp = CodeBreakNew ();
  if (cbp == NULL) return FALSE;
  cbp->aa = cp;

  /* parse location */
  pos = SimpleValuePos (val);
  cbp->loc = Nlm_gbparseint (pos, &locmap, &sitesmap, &num_errs, sip);
  MemFree (pos);
  if (cbp->loc == NULL) {
    CodeBreakFree (cbp);
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "transl_except parsing failed, %s, drop the transl_except", pos);
    return FALSE;
  }
  diff = SeqLocStop(cbp->loc) - SeqLocStart(cbp->loc);
  sintp = cbp->loc->data.ptrvalue;
  if (sintp == NULL) return FALSE;
  if (sintp->from > sintp->to) {
    temp = sintp->from;
    sintp->from = sintp->to;
    sintp->to = temp;
  }
  sintp->strand = SeqLocStrand (sfp->location);
  strand = sintp->strand;
  if ((diff != 2 && (strand != Seq_strand_minus)) ||
      (diff != -2 && (strand == Seq_strand_minus))) {
    pos_range = TRUE;
  }
  if (num_errs > 0 || pos_range) {
    CodeBreakFree (cbp);
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "transl_except range is wrong, %s, drop the transl_except", pos);
    return FALSE;
  }
  if (SeqLocCompare (sfp->location, cbp->loc) != SLC_B_IN_A) {
    CodeBreakFree (cbp);
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "/transl_except not in CDS: %s", val);
    return FALSE;
  }

  /* add to code break list */
  lastcbp = crp->code_break;
  if (lastcbp == NULL) {
    crp->code_break = cbp;
  } else {
     while (lastcbp->next != NULL) {
      lastcbp = lastcbp->next;
    }
    lastcbp->next = cbp;
  }
  return TRUE;
}

static Boolean HandledGBQualOnCDS (SeqFeatPtr sfp, GBQualPtr gbq)

{
  BioseqPtr       bsp;
  Int2            choice = 0;
  Boolean         pos_range = FALSE;
  SeqFeatPtr      prot = NULL;
  ProtRefPtr      prp = NULL;
  SeqAnnotPtr     sap;
  SeqFeatPtr      tmp;
  ValNode         vn;
  SeqFeatXrefPtr  xref;

  if (StringICmp (gbq->qual, "product") == 0) {
    choice = 1;
  } else if (StringICmp (gbq->qual, "function") == 0) {
    choice = 2;
  } else if (StringICmp (gbq->qual, "EC_number") == 0) {
    choice = 3;
  }
  if (choice > 0) {
    /* first see if protein product exists */
    bsp = BioseqFindFromSeqLoc (sfp->product);
    if (bsp != NULL && bsp->repr == Seq_repr_raw) {
      vn.choice = SEQLOC_WHOLE;
      vn.data.ptrvalue = (Pointer) SeqIdFindBest (bsp->id, 0);
      vn.next = NULL;
      for (sap = bsp->annot; sap != NULL && prot == NULL; sap = sap->next) {
        if (sap->type == 1) {
          for (tmp = (SeqFeatPtr) sap->data; tmp != NULL && prot == NULL; tmp = tmp->next) {
            if (tmp->data.choice == SEQFEAT_PROT) {
              if (SeqLocCompare (tmp->location, &vn)) {
                /* find first protein feature packaged on and located on bioseq */
                prot = tmp;
              }
            }
          }
        }
      }
    }
    if (prot != NULL) {
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
    }
    if (prp == NULL) {
      /* otherwise make cross reference */
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
        xref = xref->next;
      }
      if (xref == NULL) {
        prp = ProtRefNew ();
        if (prp == NULL) return FALSE;
        xref = SeqFeatXrefNew ();
        if (xref == NULL) return FALSE;
        xref->data.choice = SEQFEAT_PROT;
        xref->data.value.ptrvalue = (Pointer) prp;
        xref->next = sfp->xref;
        sfp->xref = xref;
      }
      if (xref != NULL) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
      }
    }
    if (prp == NULL) return FALSE;
    switch (choice) {
      case 1 :
        ValNodeCopyStr (&(prp->name), 0, gbq->val);
        break;
      case 2 :
        ValNodeCopyStr (&(prp->activity), 0, gbq->val);
        break;
      case 3 :
        ValNodeCopyStr (&(prp->ec), 0, gbq->val);
        break;
      default :
        break;
    }
    return TRUE;
  }

  if (StringICmp (gbq->qual, "transl_except") == 0) {
    return ParseCodeBreak (sfp, gbq->val);
  }
  return FALSE;
}

static Boolean HandledGBQualOnRNA (SeqFeatPtr sfp, GBQualPtr gbq)

{
  Uint1      aa;
  Int2       j;
  CharPtr    name;
  RnaRefPtr  rrp;
  tRNAPtr    trp;

  if (StringICmp (gbq->qual, "product") == 0 ||
      StringICmp (gbq->qual, "standard_name") == 0) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp == NULL) return FALSE;
    if (rrp->type == 0) {
      rrp->type = 255;
    }
    if (rrp->type == 3 && rrp->ext.choice == 1) {
      name = (CharPtr) rrp->ext.value.ptrvalue;
      aa = ParseTRnaString (name, NULL);
      if (aa != 0) {
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        trp = (tRNAPtr) MemNew (sizeof (tRNA));
        if (trp != NULL) {
          trp->aatype = 2;
          for (j = 0; j < 6; j++) {
            trp->codon [j] = 255;
          }
          trp->aa = aa;
          rrp->ext.choice = 2;
          rrp->ext.value.ptrvalue = (Pointer) trp;
        }
      }
    }
    if (rrp->type == 3 && rrp->ext.choice == 0) {
      AddQualifierToFeature (sfp, "product", gbq->val);
      return TRUE;
    }
    if (rrp->type == 3 && rrp->ext.choice == 2) {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL && trp->aatype == 2) {
        if (trp->aa == ParseTRnaString (gbq->val, NULL)) {
          return TRUE;
        }
      }
    }
    if (rrp->ext.choice != 0 && rrp->ext.choice != 1) return FALSE;
    if (! HasNoText ((CharPtr) rrp->ext.value.ptrvalue)) {
      if (StringICmp ((CharPtr) rrp->ext.value.ptrvalue, gbq->val) == 0) {
        return TRUE;
      }
      return FALSE;
    }
    if (rrp->ext.choice == 1 && rrp->ext.value.ptrvalue != NULL) {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    }
    if (rrp->ext.choice == 0 || rrp->ext.choice == 1) {
      rrp->ext.choice = 1;
      rrp->ext.value.ptrvalue = StringSave (gbq->val);
      return TRUE;
    }
  }
  return FALSE;
}

static Boolean HandledGBQualOnProt (SeqFeatPtr sfp, GBQualPtr gbq)

{
  Int2        choice = 0;
  ProtRefPtr  prp;
  ValNodePtr  vnp;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return FALSE;
  if (StringICmp (gbq->qual, "product") == 0) {
    choice = 1;
  } else if (StringICmp (gbq->qual, "function") == 0) {
    choice = 2;
  } else if (StringICmp (gbq->qual, "EC_number") == 0) {
    choice = 3;
  }
  if (choice == 1) {
    vnp = prp->name;
    if (vnp != NULL && (! HasNoText (vnp->data.ptrvalue))) return FALSE;
    ValNodeCopyStr (&(prp->name), 0, gbq->val);
    vnp = prp->name;
    if (vnp != NULL && prp->desc != NULL) {
      if (StringICmp (vnp->data.ptrvalue, prp->desc) == 0) {
        prp->desc = MemFree (prp->desc);
      }
    }
    return TRUE;
  } else if (choice == 2) {
    ValNodeCopyStr (&(prp->activity), 0, gbq->val);
    return TRUE;
  } else if (choice == 3) {
    ValNodeCopyStr (&(prp->ec), 0, gbq->val);
    return TRUE;
  }
  return TRUE; /* gbquals not appropriate on protein features */
}

static Boolean HandledGBQualOnImp (SeqFeatPtr sfp, GBQualPtr gbq)

{
  Char        ch;
  ImpFeatPtr  ifp;
  Int4        len;
  CharPtr     ptr;

  if (StringICmp (gbq->qual, "rpt_unit") == 0) {
    if (HasNoText (gbq->val)) return FALSE;
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (ifp == NULL) return FALSE;
    if (StringICmp (ifp->key, "repeat_region") != 0) return FALSE;
    len = SeqLocLen (sfp->location);
    if (len != (Int4) StringLen (gbq->val)) return FALSE;
    ptr = gbq->val;
    ch = *ptr;
    while (ch != '\0') {
      if (StringChr ("ACGTNacgtn", ch) == NULL) return FALSE;
      ptr++;
      ch = *ptr;
    }
    return TRUE;
  }
  return FALSE;
}

static CharPtr TrimParenthesesAndCommasAroundString (CharPtr str)

{
  Uchar    ch;	/* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && (ch < ' ' || ch == '(' || ch == ',')) {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ')' && ch != ',') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static CharPtr CombineSplitQual (CharPtr origval, CharPtr newval)

{
  size_t   len;
  CharPtr  str = NULL;

  len = StringLen (origval) + StringLen (newval) + 5;
  str = MemNew (sizeof (Char) * len);
  if (str == NULL) return origval;
  TrimParenthesesAndCommasAroundString (origval);
  TrimParenthesesAndCommasAroundString (newval);
  StringCpy (str, "(");
  StringCat (str, origval);
  StringCat (str, ",");
  StringCat (str, newval);
  StringCat (str, ")");
  /* free original string, knowing return value will replace it */
  MemFree (origval);
  return str;
}

static Boolean StringIsJustQuotes (CharPtr str)

{
  Nlm_Uchar  ch;	/* to use 8bit characters in multibyte languages */

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ' && ch != '"' && ch != '\'') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static void CleanupFeatureGBQuals (SeqFeatPtr sfp)

{
  DbtagPtr        db;
  GBQualPtr       gbq;
  GeneRefPtr      grp;
  size_t          len;
  GBQualPtr       nextqual;
  ObjectIdPtr     oip;
  GBQualPtr PNTR  prevqual;
  CharPtr         ptr;
  GBQualPtr       rpt_type = NULL;
  GBQualPtr       rpt_unit = NULL;
  CharPtr         str;
  CharPtr         tag;
  Boolean         unlink;
  GBQualPtr       usedin = NULL;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    CleanVisString (&(gbq->qual));
    CleanVisString (&(gbq->val));
    if (gbq->qual == NULL) {
      gbq->qual = StringSave ("");
    }
    if (StringIsJustQuotes (gbq->val)) {
      gbq->val = MemFree (gbq->val);
    }
    if (gbq->val == NULL) {
      gbq->val = StringSave ("");
    }
    nextqual = gbq->next;
    unlink = TRUE;
    if (StringICmp (gbq->qual, "partial") == 0) {
      sfp->partial = TRUE;
    } else if (StringICmp (gbq->qual, "evidence") == 0) {
      if (StringICmp (gbq->val, "experimental") == 0) {
        if (sfp->exp_ev != 2) {
          sfp->exp_ev = 1;
        }
      } else if (StringICmp (gbq->val, "not_experimental") == 0) {
        sfp->exp_ev = 2;
      }
    } else if (StringICmp (gbq->qual, "exception") == 0) {
      sfp->excpt = TRUE;
      if (! HasNoText (gbq->val)) {
        if (StringICmp (gbq->val, "TRUE") != 0) {
          if (sfp->except_text == NULL) {
            sfp->except_text = StringSaveNoNull (gbq->val);
          }
        }
      }
    } else if (StringICmp (gbq->qual, "note") == 0) {
      if (sfp->comment == NULL) {
        sfp->comment = gbq->val;
        gbq->val = NULL;
      } else {
        len = StringLen (sfp->comment) + StringLen (gbq->val) + 5;
        str = MemNew (sizeof (Char) * len);
        StringCpy (str, sfp->comment);
        StringCat (str, "; ");
        StringCat (str, gbq->val);
        sfp->comment = MemFree (sfp->comment);
        sfp->comment = str;
      }
    } else if (StringICmp (gbq->qual, "db_xref") == 0) {
      vnp = ValNodeNew (NULL);
      db = DbtagNew ();
      vnp->data.ptrvalue = db;
      tag = gbq->val;
      ptr = StringChr (tag, ':');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
        db->db = StringSave (tag);
        oip = ObjectIdNew ();
        oip->str = StringSave (ptr);
        db->tag = oip;
      } else {
        db->db = StringSave ("?");
        oip = ObjectIdNew ();
        oip->str = StringSave (tag);
        db->tag = oip;
      }
      vnp->next = sfp->dbxref;
      sfp->dbxref = vnp;
    } else if (StringICmp (gbq->qual, "gdb_xref") == 0) {
      vnp = ValNodeNew (NULL);
      db = DbtagNew ();
      vnp->data.ptrvalue = db;
      db->db = StringSave ("GDB");
      oip = ObjectIdNew ();
      oip->str = StringSave (gbq->val);
      db->tag = oip;
    } else if (sfp->data.choice == SEQFEAT_GENE && HandledGBQualOnGene (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_CDREGION && HandledGBQualOnCDS (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_RNA && HandledGBQualOnRNA (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_PROT && HandledGBQualOnProt (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_IMP && HandledGBQualOnImp (sfp, gbq)) {
    } else if (StringICmp (gbq->qual, "rpt_type") == 0) {
      if (rpt_type == NULL) {
        rpt_type = gbq;
        unlink = FALSE;
      } else {
        rpt_type->val = CombineSplitQual (rpt_type->val, gbq->val);
      }
    } else if (StringICmp (gbq->qual, "rpt_unit") == 0) {
      if (rpt_unit == NULL) {
        rpt_unit = gbq;
        unlink = FALSE;
      } else {
        rpt_unit->val = CombineSplitQual (rpt_unit->val, gbq->val);
      }
    } else if (StringICmp (gbq->qual, "usedin") == 0) {
      if (usedin == NULL) {
        usedin = gbq;
        unlink = FALSE;
      } else {
        usedin->val = CombineSplitQual (usedin->val, gbq->val);
      }
    } else if (StringICmp (gbq->qual, "pseudo") == 0) {
      sfp->pseudo = TRUE;
      unlink = FALSE; /* field, but do not remove from gbqual list until editors fixed */
    } else if (StringICmp (gbq->qual, "gene") == 0 && (! StringHasNoText (gbq->val))) {
      grp = GeneRefNew ();
      grp->locus = StringSave (gbq->val);
      xref = SeqFeatXrefNew ();
      xref->data.choice = SEQFEAT_GENE;
      xref->data.value.ptrvalue = (Pointer) grp;
      xref->next = sfp->xref;
      sfp->xref = xref;
    } else {
      unlink = FALSE;
    }
    if (unlink) {
      *(prevqual) = gbq->next;
      gbq->next = NULL;
      gbq->qual = MemFree (gbq->qual);
      gbq->val = MemFree (gbq->val);
      GBQualFree (gbq);
    } else {
      prevqual = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = nextqual;
  }
}

static int LIBCALLBACK SortByGBQualKey (VoidPtr ptr1, VoidPtr ptr2)

{
  GBQualPtr  gbq1;
  GBQualPtr  gbq2;
  CharPtr    str1;
  CharPtr    str2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  gbq1 = *((GBQualPtr PNTR) ptr1);
  gbq2 = *((GBQualPtr PNTR) ptr2);
  if (gbq1 == NULL || gbq2 == NULL) return 0;
  str1 = (CharPtr) gbq1->qual;
  str2 = (CharPtr) gbq2->qual;
  if (str1 == NULL || str2 == NULL) return 0;
  return StringICmp (str1, str2);
}

static Boolean GBQualsAlreadyInOrder (GBQualPtr list)

{
  GBQualPtr  curr;
  GBQualPtr  next;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    if (StringICmp (curr->qual, next->qual) > 0) return FALSE;
    curr = next;
    next = curr->next;
  }
  return TRUE;
}

static GBQualPtr SortFeatureGBQuals (GBQualPtr list)

{
  size_t     count, i;
  GBQualPtr  gbq, PNTR head;

  if (list == NULL) return NULL;
  if (GBQualsAlreadyInOrder (list)) return list;

  for (gbq = list, count = 0; gbq != NULL; gbq = gbq->next, count++) continue;
  head = MemNew (sizeof (GBQualPtr) * (count + 1));

  for (gbq = list, i = 0; gbq != NULL && i < count; i++) {
    head [i] = gbq;
    gbq = gbq->next;
  }

  HeapSort (head, count, sizeof (GBQualPtr), SortByGBQualKey);

  for (i = 0; i < count; i++) {
    gbq = head [i];
    gbq->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}

/* this identifies gbquals that should have been placed into special fields */

#define NUM_ILLEGAL_QUALS 15

/* StringICmp use of TO_UPPER means translation should go before transl_XXX */

static CharPtr illegalGbqualList [NUM_ILLEGAL_QUALS] = {
  "allele",
  "anticodon",
  "citation",
  "codon",
  "codon_start",
  "cons_splice",
  "db_xref",
  "frequency",
  "gene",
  "note",
  "PCR_conditions",
  "protein_id",
  "translation"
  "transl_except",
  "transl_table",
};

static Int2 QualifierIsIllegal (CharPtr qualname)

{
  Int2  L, R, mid;

  if (qualname == NULL || *qualname == '\0') return FALSE;

  L = 0;
  R = NUM_ILLEGAL_QUALS - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (illegalGbqualList [mid], qualname) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (illegalGbqualList [R], qualname) == 0) {
    return TRUE;
  }

  return FALSE;
}

static void GbqualLink (GBQualPtr PNTR head, GBQualPtr qual)

{
  GBQualPtr  gbq;

  if (head == NULL || qual == NULL) return;
  gbq = *head;
  if (gbq != NULL) {
    while (gbq->next != NULL) {
      gbq = gbq->next;
    }
    gbq->next = qual;
  } else {
    *head = qual;
  }
}

static GBQualPtr SortIllegalGBQuals (GBQualPtr list)

{
  GBQualPtr  gbq, next, legal = NULL, illegal = NULL;

  gbq = list;
  while (gbq != NULL) {
    next = gbq->next;
    gbq->next = NULL;
    if (QualifierIsIllegal (gbq->qual)) {
      GbqualLink (&illegal, gbq);
    } else {
      GbqualLink (&legal, gbq);
    }
    gbq = next;
  }
  GbqualLink (&legal, illegal);
  return legal;
}

static void CleanOrgModList (OrgModPtr PNTR ompp)

{
  OrgModPtr       next;
  OrgModPtr       omp;
  OrgModPtr PNTR  prev;

  if (ompp == NULL) return;
  prev = ompp;
  omp = *ompp;
  while (omp != NULL) {
    next = omp->next;
    CleanVisString (&(omp->subname));
    CleanVisString (&(omp->attrib));
    if (HasNoText (omp->subname)) {
      *prev = omp->next;
      omp->next = NULL;
      OrgModFree (omp);
    } else {
      prev = &(omp->next);
    }
    omp = next;
  }
}

static int LIBCALLBACK SortByOrgModSubtype (VoidPtr ptr1, VoidPtr ptr2)

{
  OrgModPtr  omp1;
  OrgModPtr  omp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  omp1 = *((OrgModPtr PNTR) ptr1);
  omp2 = *((OrgModPtr PNTR) ptr2);
  if (omp1 == NULL || omp2 == NULL) return 0;
  if (omp1->subtype > omp2->subtype) {
    return 1;
  } else if (omp1->subtype < omp2->subtype) {
    return -1;
  }
  return 0;
}

static Boolean OrgModsAlreadyInOrder (OrgModPtr list)

{
  OrgModPtr  curr;
  OrgModPtr  next;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    if (curr->subtype > next->subtype) return FALSE;
    curr = next;
    next = curr->next;
  }
  return TRUE;
}

static OrgModPtr SortOrgModList (OrgModPtr list)

{
  size_t     count, i;
  OrgModPtr  omp, PNTR head;

  if (list == NULL) return NULL;
  if (OrgModsAlreadyInOrder (list)) return list;

  for (omp = list, count = 0; omp != NULL; omp = omp->next, count++) continue;
  head = MemNew (sizeof (OrgModPtr) * (count + 1));

  for (omp = list, i = 0; omp != NULL && i < count; i++) {
    head [i] = omp;
    omp = omp->next;
  }

  HeapSort (head, count, sizeof (OrgModPtr), SortByOrgModSubtype);

  for (i = 0; i < count; i++) {
    omp = head [i];
    omp->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}

static void CleanSubSourceList (SubSourcePtr PNTR sspp)

{
  SubSourcePtr       next;
  SubSourcePtr PNTR  prev;
  SubSourcePtr       ssp;

  if (sspp == NULL) return;
  prev = sspp;
  ssp = *sspp;
  while (ssp != NULL) {
    next = ssp->next;
    if (ssp->subtype != 14 && ssp->subtype != 15) {
      CleanVisString (&(ssp->name));
    }
    CleanVisString (&(ssp->attrib));
    if (HasNoText (ssp->name) && ssp->subtype != 14 && ssp->subtype != 15) {
      *prev = ssp->next;
      ssp->next = NULL;
      SubSourceFree (ssp);
    } else {
      prev = &(ssp->next);
    }
    ssp = next;
  }
}

static int LIBCALLBACK SortBySubSourceSubtype (VoidPtr ptr1, VoidPtr ptr2)

{
  SubSourcePtr  ssp1;
  SubSourcePtr  ssp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  ssp1 = *((SubSourcePtr PNTR) ptr1);
  ssp2 = *((SubSourcePtr PNTR) ptr2);
  if (ssp1 == NULL || ssp2 == NULL) return 0;
  if (ssp1->subtype > ssp2->subtype) {
    return 1;
  } else if (ssp1->subtype < ssp2->subtype) {
    return -1;
  }
  return 0;
}

static Boolean SubSourceAlreadyInOrder (SubSourcePtr list)

{
  SubSourcePtr  curr;
  SubSourcePtr  next;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    if (curr->subtype > next->subtype) return FALSE;
    curr = next;
    next = curr->next;
  }
  return TRUE;
}

static SubSourcePtr SortSubSourceList (SubSourcePtr list)

{
  size_t        count, i;
  SubSourcePtr  ssp, PNTR head;

  if (list == NULL) return NULL;
  if (SubSourceAlreadyInOrder (list)) return list;

  for (ssp = list, count = 0; ssp != NULL; ssp = ssp->next, count++) continue;
  head = MemNew (sizeof (SubSourcePtr) * (count + 1));

  for (ssp = list, i = 0; ssp != NULL && i < count; i++) {
    head [i] = ssp;
    ssp = ssp->next;
  }

  HeapSort (head, count, sizeof (SubSourcePtr), SortBySubSourceSubtype);

  for (i = 0; i < count; i++) {
    ssp = head [i];
    ssp->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}

static void OrpModToSubSource (ValNodePtr PNTR vnpp, SubSourcePtr PNTR sspp)

{
  SubSourcePtr  ssp;
  ValNodePtr    vnp;

  if (vnpp == NULL || sspp == NULL) return;
  for (vnp = *vnpp; vnp != NULL; vnp = vnp->next) {
    ssp = SubSourceNew ();
    if (ssp != NULL) {
      ssp->subtype = 255;
      ssp->name = (CharPtr) vnp->data.ptrvalue;
      vnp->data.ptrvalue = NULL;
      if (*sspp != NULL) {
        ssp->next = *sspp;
      }
      *sspp = ssp;
    }
  }
  vnp = *vnpp;
  *vnpp = NULL;
  ValNodeFree (vnp);
}

static int LIBCALLBACK SortDbxref (VoidPtr ptr1, VoidPtr ptr2)

{
  int          compare;
  DbtagPtr     dbt1;
  DbtagPtr     dbt2;
  ObjectIdPtr  oip1;
  ObjectIdPtr  oip2;
  CharPtr      str1;
  CharPtr      str2;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  dbt1 = (DbtagPtr) vnp1->data.ptrvalue;
  dbt2 = (DbtagPtr) vnp2->data.ptrvalue;
  if (dbt1 == NULL || dbt2 == NULL) return 0;
  str1 = (CharPtr) dbt1->db;
  str2 = (CharPtr) dbt2->db;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  if (compare != 0) return compare;
  oip1 = dbt1->tag;
  oip2 = dbt2->tag;
  if (oip1 == NULL || oip2 == NULL) return 0;
  str1 = oip1->str;
  str2 = oip2->str;
  if (str1 != NULL && str2 != NULL) {
    return StringICmp (str1, str2);
  } else if (str1 == NULL && str2 == NULL) {
    if (oip1->id > oip2->id) {
      return 1;
    } else if (oip1->id < oip2->id) {
      return -1;
    }
  }
  return 0;
}

static void FixNumericDbxrefs (ValNodePtr vnp)

{
  Char         ch;
  DbtagPtr     dbt;
  Boolean      isNum;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  long         val;

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      oip = dbt->tag;
      if (oip != NULL) {
        ptr = oip->str;
        if (ptr != NULL) {
          isNum = TRUE;
          ch = *ptr;
          while (ch != '\0') {
            if ((! IS_DIGIT (ch)) && (! IS_WHITESP (ch))) {
              isNum = FALSE;
            }
            ptr++;
            ch = *ptr;
          }
          if (isNum) {
            if (sscanf (oip->str, "%ld", &val) == 1) {
              oip->id = (Int4) val;
              oip->str = MemFree (oip->str);
            }
          }
        }
      }
    }
    vnp = vnp->next;
  }
}

static void CleanupDuplicateDbxrefs (ValNodePtr PNTR prevvnp)

{
  DbtagPtr     dbt;
  DbtagPtr     last = NULL;
  ValNodePtr   nextvnp;
  ObjectIdPtr  oip1;
  ObjectIdPtr  oip2;
  CharPtr      str1;
  CharPtr      str2;
  Boolean      unlink;
  ValNodePtr   vnp;

  if (prevvnp == NULL) return;
  vnp = *prevvnp;
  FixNumericDbxrefs (vnp);
  while (vnp != NULL) {
    nextvnp = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      unlink = FALSE;
      if (last != NULL) {
        str1 = (CharPtr) dbt->db;
        str2 = (CharPtr) last->db;
        if (str1 != NULL && str2 != NULL && StringICmp (str1, str2) == 0) {
          oip1 = dbt->tag;
          oip2 = last->tag;
          if (oip1 != NULL && oip2 != NULL) {
            str1 = oip1->str;
            str2 = oip2->str;
            if (str1 != NULL && str2 != NULL) {
              if (StringICmp (str1, str2) == 0) {
                unlink = TRUE;
              }
            } else if (str1 == NULL && str2 == NULL) {
              if (oip1->id == oip2->id) {
                unlink = TRUE;
              }
            }
          }
        }
      } else {
        last = dbt;
      }
      if (unlink) {
        *prevvnp = vnp->next;
        vnp->next = NULL;
        DbtagFree (dbt);
        ValNodeFree (vnp);
      } else {
        last = dbt;
        prevvnp = (ValNodePtr PNTR) &(vnp->next);
      }
    }
    vnp = nextvnp;
  }
}

static void NormalizeAuthors (AuthListPtr alp)

{
  AuthorPtr    ap;
  Char         ch;
  CharPtr      initials;
  Int2         j;
  Int2         k;
  size_t       len;
  ValNodePtr   names;
  NameStdPtr   nsp;
  CharPtr      periods;
  PersonIdPtr  pid;

  if (alp == NULL || alp->choice != 1) return;
  for (names = alp->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (ap != NULL) {
      pid = ap->name;
      if (pid != NULL && pid->choice == 2) {
        nsp = pid->data;
        if (nsp != NULL && nsp->names [4] != NULL) {
          initials = nsp->names [4];
          len = MAX ((size_t) (StringLen (initials) * 2 + 4), (size_t) 64);
          periods = MemNew (len);
          if (periods == NULL) return;
          periods [0] = '\0';
          j = 0;
          k = 0;
          ch = initials [j];
          while (ch != '\0') {
            if (ch == '-') {
              periods [k] = ch;
              k++;
              j++;
              ch = initials [j];
            } else if (ch == '.') {
              j++;
              ch = initials [j];
            } else if (ch == ' ') {
              j++;
              ch = initials [j];
            } else {
              periods [k] = ch;
              k++;
              j++;
              ch = initials [j];
              periods [k] = '.';
              k++;
            }
          }
          periods [k] = '\0';
          nsp->names [4] = MemFree (nsp->names [4]);
          nsp->names [4] = StringSave (periods);
          MemFree (periods);
          CleanVisString (&(nsp->names [0]));
          CleanVisString (&(nsp->names [1]));
          CleanVisString (&(nsp->names [2]));
          CleanVisString (&(nsp->names [3]));
          CleanVisString (&(nsp->names [4]));
          CleanVisString (&(nsp->names [5]));
          CleanVisString (&(nsp->names [6]));
        }
      }
    }
  }
}

/* from utilpub.c */
static Boolean empty_citgen(CitGenPtr  cit)
{
	if (cit == NULL)
		return TRUE;
	if (cit->cit)
		return FALSE;
	if (cit->authors)
		return FALSE;
	if (cit->muid > 0)
		return FALSE;
	if (cit->journal)
		return FALSE;
	if (cit->volume)
		return FALSE;
	if (cit->issue)
		return FALSE;
	if (cit->pages)
		return FALSE;
	if (cit->date)
		return FALSE;
	if (cit->serial_number > 0)
		return FALSE;
	if (cit->title)
		return FALSE;
	if (cit->pmid > 0)
		return FALSE;
	return TRUE;
}

static void NormalizeAPub (ValNodePtr vnp, Boolean stripSerial)

{
  AuthListPtr  alp;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitPatPtr    cpp;
  CitSubPtr    csp;
  ImprintPtr   imp;

  if (vnp == NULL) return;
  if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) return;
  if (vnp->data.ptrvalue == NULL) return;
  switch (vnp->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cgp->authors);
      if (stripSerial) {
        cgp->serial_number = -1; /* but does not remove if empty */
      }
      if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
        cgp->date = DateFree (cgp->date); /* remove date if unpublished */
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) vnp->data.ptrvalue;
      NormalizeAuthors (csp->authors);
      alp = csp->authors;
      imp = csp->imp;
      if (alp != NULL && alp->affil == NULL && imp != NULL && imp->pub != NULL) {
        alp->affil = imp->pub;
        imp->pub = NULL;
      }
      if (csp->date == NULL && imp != NULL && imp->date != NULL) {
        csp->date = imp->date;
        imp->date = NULL;
      }
      if (imp != NULL && imp->pub == NULL) {
        csp->imp = ImprintFree (csp->imp);
      }
      break;
    case PUB_Article :
      cap = (CitArtPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cap->authors);
      break;
    case PUB_Book :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cbp->authors);
      break;
    case PUB_Man :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      if (cbp->othertype == 2 && cbp->let_type == 3) {
        NormalizeAuthors (cbp->authors);
      }
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cpp->authors);
      NormalizeAuthors (cpp->applicants);
      NormalizeAuthors (cpp->assignees);
      break;
    default :
      break;
  }
}

static void NormalizePubdesc (PubdescPtr pdp, Boolean stripSerial)

{
  CitGenPtr        cgp;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (pdp == NULL) return;
  prev = &(pdp->pub);
  vnp = pdp->pub;
  if (vnp != NULL && vnp->next == NULL && vnp->choice == PUB_Gen) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    NormalizeAuthors (cgp->authors);
    if (stripSerial) {
      cgp->serial_number = -1;
    }
    if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
      cgp->date = DateFree (cgp->date); /* remove date if unpublished */
    }
    return; /* but does not remove if empty and only element of Pub */
  }
  while (vnp != NULL) {
    next = vnp->next;
    NormalizeAPub (vnp, stripSerial);
    if (vnp->choice == PUB_Gen && empty_citgen ((CitGenPtr) vnp->data.ptrvalue)) {
      *prev = vnp->next;
      vnp->next = NULL;
      PubFree (vnp);
    } else {
      prev = &(vnp->next);
    }
    vnp = next;
  }
}

extern Uint1 FindTrnaAA3 (CharPtr str);

static void CleanupTrna (SeqFeatPtr sfp, tRNAPtr trp)

{
  Uint1           aa;
  Uint1           from = 0;
  Boolean         justTrnaText;
  SeqMapTablePtr  smtp;

  /* look for tRNA-OTHER with actual amino acid in comment */

  if (sfp == NULL || sfp->comment == NULL || trp == NULL) return;
  if (trp->aatype == 2) {
    aa = trp->aa;
  } else {
    switch (trp->aatype) {
      case 0 :
        from = 0;
        break;
      case 1 :
        from = Seq_code_iupacaa;
        break;
      case 2 :
        from = Seq_code_ncbieaa;
        break;
      case 3 :
        from = Seq_code_ncbi8aa;
        break;
      case 4 :
        from = Seq_code_ncbistdaa;
        break;
      default:
        break;
    }
    smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
    if (smtp != NULL) {
      aa = SeqMapTableConvert (smtp, trp->aa);
    }
  }
  if (aa != 'X') return;
  aa = ParseTRnaString (sfp->comment, &justTrnaText);
  if (aa == 0) return;
  trp->aa = aa;
  if (justTrnaText) {
    sfp->comment = MemFree (sfp->comment);
  }
}

static void CleanupFeatureStrings (SeqFeatPtr sfp, Boolean stripSerial)

{
  Uint1         aa;
  BioSourcePtr  biop;
  GeneRefPtr    grp;
  ImpFeatPtr    ifp;
  Int2          j;
  Boolean       justTrnaText;
  size_t        len;
  CharPtr       name;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  PubdescPtr    pdp;
  ProtRefPtr    prp;
  RnaRefPtr     rrp;
  CharPtr       str;
  tRNAPtr       trp;

  if (sfp == NULL) return;
  CleanVisString (&(sfp->comment));
  CleanVisString (&(sfp->title));
  CleanVisString (&(sfp->except_text));
  CleanDoubleQuote (sfp->comment);
  switch (sfp->data.choice) {
    case SEQFEAT_BOND :
    case SEQFEAT_SITE :
    case SEQFEAT_PSEC_STR :
    case SEQFEAT_COMMENT:
      return;
    default :
      break;
  }
  if (sfp->data.value.ptrvalue == NULL) return;
  orp = NULL;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      CleanVisString (&(grp->locus));
      CleanVisString (&(grp->allele));
      CleanVisString (&(grp->desc));
      CleanVisString (&(grp->maploc));
      CleanVisStringList (&(grp->syn));
      CleanDoubleQuote (grp->locus);
      CleanDoubleQuote (grp->allele);
      CleanDoubleQuote (grp->desc);
      CleanDoubleQuote (grp->maploc);
      CleanDoubleQuoteList (grp->syn);
      grp->db = ValNodeSort (grp->db, SortDbxref);
      CleanupDuplicateDbxrefs (&(grp->db));
      break;
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      break;
    case SEQFEAT_CDREGION :
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      CleanVisString (&(prp->desc));
      CleanVisStringList (&(prp->name));
      CleanVisStringList (&(prp->ec));
      CleanVisStringList (&(prp->activity));
      CleanDoubleQuote (prp->desc);
      CleanDoubleQuoteList (prp->name);
      CleanDoubleQuoteList (prp->ec);
      CleanDoubleQuoteList (prp->activity);
      prp->db = ValNodeSort (prp->db, SortDbxref);
      CleanupDuplicateDbxrefs (&(prp->db));
      if (prp->processed != 3 && prp->processed != 4 &&
          prp->name == NULL && sfp->comment != NULL) {
        ValNodeAddStr (&(prp->name), 0, sfp->comment);
        sfp->comment = NULL;
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->ext.choice == 1) {
        CleanVisString ((CharPtr PNTR) &(rrp->ext.value.ptrvalue));
        CleanDoubleQuote ((CharPtr) rrp->ext.value.ptrvalue);
        if (rrp->ext.value.ptrvalue == NULL) {
          rrp->ext.choice = 0;
        } else if (rrp->type == 4) {
          name = (CharPtr) rrp->ext.value.ptrvalue;
          len = StringLen (name);
          if (len > 5) {
            if (StringNICmp (name + len - 5, " rRNA", 5) == 0) {
              str = MemNew (len + 10);
              if (str != NULL) {
                StringNCpy (str, name, len - 5);
                StringCat (str, " ribosomal RNA");
                rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                rrp->ext.value.ptrvalue = (Pointer) str;
              }
            }
          }
        } else if (rrp->type == 3) {
          name = (CharPtr) rrp->ext.value.ptrvalue;
          aa = ParseTRnaString (name, &justTrnaText);
          if (aa != 0) {
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
            trp = (tRNAPtr) MemNew (sizeof (tRNA));
            if (trp != NULL) {
              trp->aatype = 2;
              for (j = 0; j < 6; j++) {
                trp->codon [j] = 255;
              }
              trp->aa = aa;
              rrp->ext.choice = 2;
              rrp->ext.value.ptrvalue = (Pointer) trp;
            }
          }
        }
      } else if (rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        CleanupTrna (sfp, trp);
      } else if (rrp->type == 3 && (! StringHasNoText (sfp->comment))) {
        aa = ParseTRnaString (sfp->comment, &justTrnaText);
        if (aa != 0) {
          trp = (tRNAPtr) MemNew (sizeof (tRNA));
          if (trp != NULL) {
            trp->aatype = 2;
            for (j = 0; j < 6; j++) {
              trp->codon [j] = 255;
            }
            trp->aa = aa;
            rrp->ext.choice = 2;
            rrp->ext.value.ptrvalue = (Pointer) trp;
            if (justTrnaText) {
              sfp->comment = MemFree (sfp->comment);
            }
          }
        }
      }
      if (rrp->ext.choice == 0 && sfp->comment != NULL && rrp->type != 3) {
        len = StringLen (sfp->comment);
        if (len > 5 && len < 20) {
          if (StringNICmp (sfp->comment + len - 3, "RNA", 3) == 0) {
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = sfp->comment;
            sfp->comment = NULL;
          }
        }
      }
      if (rrp->ext.choice == 1 && rrp->ext.value.ptrvalue != NULL) {
        if (StringICmp ((CharPtr) rrp->ext.value.ptrvalue, sfp->comment) == 0) {
          sfp->comment = MemFree (sfp->comment);
        }
      }
      break;
    case SEQFEAT_PUB :
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      CleanVisString (&(pdp->comment));
      CleanDoubleQuote (pdp->comment);
      NormalizePubdesc (pdp, stripSerial);
      break;
    case SEQFEAT_SEQ :
      break;
    case SEQFEAT_IMP :
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      CleanVisString (&(ifp->key));
      CleanVisString (&(ifp->loc));
      CleanVisString (&(ifp->descr));
      break;
    case SEQFEAT_REGION :
      CleanVisString ((CharPtr PNTR) &(sfp->data.value.ptrvalue));
      CleanDoubleQuote ((CharPtr) sfp->data.value.ptrvalue);
      if (sfp->data.value.ptrvalue == NULL) {
        sfp->data.choice = SEQFEAT_COMMENT;
      }
      break;
    case SEQFEAT_COMMENT :
      break;
    case SEQFEAT_BOND :
      break;
    case SEQFEAT_SITE :
      break;
    case SEQFEAT_RSITE :
      break;
    case SEQFEAT_USER :
      break;
    case SEQFEAT_TXINIT :
      break;
    case SEQFEAT_NUM :
      break;
    case SEQFEAT_PSEC_STR :
      break;
    case SEQFEAT_NON_STD_RESIDUE :
      break;
    case SEQFEAT_HET :
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      orp = biop->org;
      if (orp != NULL) {
        CleanVisStringList (&(orp->mod));
        OrpModToSubSource (&(orp->mod), &(biop->subtype));
      }
      CleanSubSourceList (&(biop->subtype));
      biop->subtype = SortSubSourceList (biop->subtype);
      break;
    default :
      break;
  }
  if (orp != NULL) {
    CleanVisString (&(orp->taxname));
    CleanVisString (&(orp->common));
    CleanVisStringList (&(orp->mod));
    CleanVisStringList (&(orp->syn));
    orp->db = ValNodeSort (orp->db, SortDbxref);
    CleanupDuplicateDbxrefs (&(orp->db));
    onp = orp->orgname;
    while (onp != NULL) {
      CleanVisString (&(onp->attrib));
      CleanVisString (&(onp->lineage));
      CleanVisString (&(onp->div));
      CleanOrgModList (&(onp->mod));
      onp->mod = SortOrgModList (onp->mod);
      onp = onp->next;
    }
  }
}

static void CleanupDescriptorStrings (ValNodePtr sdp, Boolean stripSerial)

{
  BioSourcePtr  biop;
  GBBlockPtr    gbp;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  PubdescPtr    pdp;

  if (sdp == NULL) return;
  switch (sdp->choice) {
    case Seq_descr_mol_type :
    case Seq_descr_method :
      return;
    default :
      break;
  }
  if (sdp->data.ptrvalue == NULL) return;
  orp = NULL;
  switch (sdp->choice) {
    case Seq_descr_mol_type :
      break;
    case Seq_descr_modif :
      break;
    case Seq_descr_method :
      break;
    case Seq_descr_name :
      break;
    case Seq_descr_title :
      break;
    case Seq_descr_org :
      orp = (OrgRefPtr) sdp->data.ptrvalue;
      break;
    case Seq_descr_comment :
      break;
    case Seq_descr_num :
      break;
    case Seq_descr_maploc :
      break;
    case Seq_descr_pir :
      break;
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      CleanVisStringList (&(gbp->extra_accessions));
      CleanVisStringList (&(gbp->keywords));
      CleanVisString (&(gbp->source));
      CleanVisString (&(gbp->origin));
      CleanVisString (&(gbp->date));
      CleanVisString (&(gbp->div));
      CleanVisString (&(gbp->taxonomy));
      break;
    case Seq_descr_pub :
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      CleanVisString (&(pdp->comment));
      NormalizePubdesc (pdp, stripSerial);
      break;
    case Seq_descr_region :
      break;
    case Seq_descr_user :
      break;
    case Seq_descr_sp :
      break;
    case Seq_descr_dbxref :
      break;
    case Seq_descr_embl :
      break;
    case Seq_descr_create_date :
      break;
    case Seq_descr_update_date :
      break;
    case Seq_descr_prf :
      break;
    case Seq_descr_pdb :
      break;
    case Seq_descr_het :
      break;
    case Seq_descr_source :
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      orp = biop->org;
      if (orp != NULL) {
        CleanVisStringList (&(orp->mod));
        OrpModToSubSource (&(orp->mod), &(biop->subtype));
      }
      CleanSubSourceList (&(biop->subtype));
      biop->subtype = SortSubSourceList (biop->subtype);
      break;
    case Seq_descr_molinfo :
      break;
    default :
      break;
  }
  if (orp != NULL) {
    CleanVisString (&(orp->taxname));
    CleanVisString (&(orp->common));
    CleanVisStringList (&(orp->mod));
    CleanVisStringList (&(orp->syn));
    onp = orp->orgname;
    while (onp != NULL) {
      CleanVisString (&(onp->attrib));
      CleanVisString (&(onp->lineage));
      CleanVisString (&(onp->div));
      CleanOrgModList (&(onp->mod));
      onp->mod = SortOrgModList (onp->mod);
      onp = onp->next;
    }
  }
}

static void CheckForSwissProtID (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  SeqIdPtr   sip;
  BoolPtr    stripSerial;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    stripSerial = (BoolPtr) mydata;
    if (stripSerial == NULL) return;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_SWISSPROT) {
        *stripSerial = FALSE;
      }
    }
  }
}

NLM_EXTERN void BasicSeqEntryCleanup (SeqEntryPtr sep)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Char          div [10];
  GBBlockPtr    gbp;
  ImpFeatPtr    ifp;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  SeqAnnotPtr   sap = NULL;
  ValNodePtr    sdp = NULL;
  SeqFeatPtr    sfp;
  Boolean       stripSerial = TRUE;
  SeqEntryPtr   tmp;

  if (sep == NULL) return;
  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtID);
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    sap = bsp->annot;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      BasicSeqEntryCleanup (tmp);
    }
    sap = bssp->annot;
    sdp = bssp->descr;
  } else return;
  biop = NULL;
  orp = NULL;
  gbp = NULL;
  div [0] = '\0';
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_IMP) {
          ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
          if (ifp != NULL && StringCmp (ifp->key, "CDS") == 0) {
            sfp->data.value.ptrvalue = ImpFeatFree (ifp);
            sfp->data.choice = SEQFEAT_CDREGION;
            sfp->data.value.ptrvalue = CdRegionNew ();
          }
        }
        CleanupFeatureGBQuals (sfp);
        sfp->qual = SortFeatureGBQuals (sfp->qual);
        sfp->qual = SortIllegalGBQuals (sfp->qual);
        CleanupFeatureStrings (sfp, stripSerial);
        sfp->dbxref = ValNodeSort (sfp->dbxref, SortDbxref);
        CleanupDuplicateDbxrefs (&(sfp->dbxref));
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    switch (sdp->choice) {
      case Seq_descr_org :
        orp = (OrgRefPtr) sdp->data.ptrvalue;
        break;
      case Seq_descr_genbank :
        gbp = (GBBlockPtr) sdp->data.ptrvalue;
        break;
      case Seq_descr_source :
        biop = (BioSourcePtr) sdp->data.ptrvalue;
        orp = biop->org;
        break;
      default :
        break;
    }
    CleanupDescriptorStrings (sdp, stripSerial);
    sdp = sdp->next;
  }

  /* copy genbank block division into biosource, if necessary */

  if (orp != NULL && gbp != NULL) {
    StringNCpy_0 (div, gbp->div, sizeof (div));
    if (StringHasNoText (div)) return;
    onp = orp->orgname;
    while (onp != NULL) {
      if (StringHasNoText (onp->div)) {
        onp->div = MemFree (onp->div);
        onp->div = StringSaveNoNull (div);
      }
      onp = onp->next;
    }
  }
}

/* end BasicSeqEntryCleanup section */

/* SeqIdStripLocus removes the SeqId.name field if accession is set */

NLM_EXTERN SeqIdPtr SeqIdStripLocus (SeqIdPtr sip)

{
  TextSeqIdPtr  tip;

  if (sip != NULL) {
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_OTHER :
        tip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tip != NULL) {
          if (! HasNoText (tip->accession)) {
            tip->name = MemFree (tip->name);
          }
        }
        break;
      default :
        break;
    }
  }
  return sip;
}

NLM_EXTERN SeqLocPtr StripLocusFromSeqLoc (SeqLocPtr location)

{
  SeqLocPtr      loc;
  SeqLocPtr      next;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  SeqPntPtr      spp;

  if (location == NULL) return NULL;
  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    next = SeqLocFindNext (location, slp);
    switch (slp->choice) {
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        SeqIdStripLocus (sip);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          SeqIdStripLocus (sinp->id);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          sip = SeqLocId (loc);
          SeqIdStripLocus (sip);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = sbp->a;
          if (spp != NULL) {
            SeqIdStripLocus (spp->id);
          }
          spp = sbp->b;
          if (spp != NULL) {
            SeqIdStripLocus (spp->id);
          }
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          SeqIdStripLocus (spp->id);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          SeqIdStripLocus (psp->id);
        }
        break;
      default :
        break;
    }
    slp = next;
  }
  return location;
}

static void GetRidOfLocusCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  sap = NULL;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1 && sap->data != NULL) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        StripLocusFromSeqLoc (sfp->location);
        StripLocusFromSeqLoc (sfp->product);
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

NLM_EXTERN void GetRidOfLocusInSeqIds (Uint2 entityID, SeqEntryPtr sep)

{
  if (entityID < 1 && sep == NULL) return;
  if (entityID > 0 && sep == NULL) {
    sep = GetTopSeqEntryForEntityID (entityID);
  }
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, GetRidOfLocusCallback);
}

/* Mac can now use static parse tables by using
   Make Strings Read-Only and Store Static Data in TOC
#ifdef OS_MAC
#define ASNLOAD_NEEDED 1
#endif
*/
#if defined(OS_DOS) || defined(WIN16)
#define ASNLOAD_NEEDED 1
#endif

static Boolean FileExists (CharPtr dirname, CharPtr subname, CharPtr filename)

{
  Char  path [PATH_MAX];

  StringNCpy_0 (path, dirname, sizeof (path));
  FileBuildPath (path, subname, NULL);
  FileBuildPath (path, NULL, filename);
  return (Boolean) (FileLength (path) > 0);
}

static Boolean CheckAsnloadPath (CharPtr dirname, CharPtr subdir)

{
#ifdef ASNLOAD_NEEDED
  Char  fname [16];
  int   i;

  for (i = 60; i <= 99; ++i) {
    sprintf (fname, "asnmedli.l%02d", (int) i);
    if (FileExists (dirname, subdir, fname)) {
      return TRUE;
    }
  }
  return FALSE;
#else
  return TRUE;
#endif
}

static Boolean CheckDataPath (CharPtr dirname, CharPtr subdir)

{
  return (Boolean) (FileExists (dirname, subdir, "seqcode.val"));
}

static Boolean CheckErrMsgPath (CharPtr dirname, CharPtr subdir)

{
  return (Boolean) (FileExists (dirname, subdir, "valid.msg"));
}

static void SetTransientPath (CharPtr dirname, CharPtr subname, CharPtr file,
                              CharPtr section, CharPtr type)

{
  Char  path [PATH_MAX];

  StringNCpy_0 (path, dirname, sizeof (path));
  FileBuildPath (path, subname, NULL);
  TransientSetAppParam (file, section, type, path);
}

NLM_EXTERN Boolean UseLocalAsnloadDataAndErrMsg (void)

{
  Boolean  asnFound;
  Boolean  dataFound;
  Char     path [PATH_MAX];
  CharPtr  ptr;

  ProgramPath (path, sizeof (path));
  ptr = StringRChr (path, DIRDELIMCHR);
  if (ptr != NULL) {
    ptr++;
    *ptr = '\0';
  }
  asnFound = CheckAsnloadPath (path, "asnload");
  dataFound = CheckDataPath (path, "data");
  if (! (asnFound && dataFound)) {
    if (ptr != NULL) {
      ptr--;
      *ptr = '\0';
      ptr = StringRChr (path, DIRDELIMCHR);
      if (ptr != NULL) {
        ptr++;
        *ptr = '\0';
      }
      asnFound = CheckAsnloadPath (path, "asnload");
      dataFound = CheckDataPath (path, "data");
    }
  }
  if (asnFound && dataFound) {
    SetTransientPath (path, "asnload", "NCBI", "NCBI", "ASNLOAD");
    SetTransientPath (path, "data", "NCBI", "NCBI", "DATA");
    if (CheckErrMsgPath (path, "errmsg")) {
      SetTransientPath (path, "errmsg", "NCBI", "ErrorProcessing", "MsgPath");
      TransientSetAppParam ("NCBI", "ErrorProcessing", "EO_BEEP", "No");
    }
    return TRUE;
  }
  return FALSE;
}

typedef struct miscdata {
  SeqEntryPtr  sep;
  Int2         count;
  Int2         desired;
} MiscData, PNTR MiscDataPtr;

static void FindNthSeqEntryCallback (SeqEntryPtr sep, Pointer mydata,
                                     Int4 index, Int2 indent)

{
  MiscDataPtr  mdp;

  if (sep != NULL && mydata != NULL) {
    mdp = (MiscDataPtr) mydata;
    (mdp->count)++;
    if (mdp->count == mdp->desired) {
      mdp->sep = sep;
    }
  }
}

NLM_EXTERN SeqEntryPtr LIBCALL FindNthSeqEntry (SeqEntryPtr sep, Int2 seq)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = seq;
  if (sep != NULL) {
    SeqEntryExplore (sep, (Pointer) (&md), FindNthSeqEntryCallback);
  }
  return md.sep;
}

NLM_EXTERN SeqEntryPtr LIBCALL FindNthBioseq (SeqEntryPtr sep, Int2 seq)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = seq;
  if (sep != NULL) {
    BioseqExplore (sep, (Pointer) (&md), FindNthSeqEntryCallback);
  }
  return md.sep;
}

NLM_EXTERN SeqEntryPtr LIBCALL FindNthSequinEntry (SeqEntryPtr sep, Int2 seq)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = seq;
  if (sep != NULL) {
    SequinEntryExplore (sep, (Pointer) (&md), FindNthSeqEntryCallback);
  }
  return md.sep;
}

static void FindNucSeqEntryCallback (SeqEntryPtr sep, Pointer mydata,
                                     Int4 index, Int2 indent)

{
  BioseqPtr    bsp;
  MiscDataPtr  mdp;

  if (sep != NULL && sep->choice == 1 && mydata != NULL) {
    mdp = (MiscDataPtr) mydata;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL && ISA_na (bsp->mol)) {
      if (mdp->sep == NULL) {
        mdp->sep = sep;
      }
    }
  }
}

NLM_EXTERN SeqEntryPtr LIBCALL FindNucSeqEntry (SeqEntryPtr sep)

{
  MiscData  md;

  md.sep = NULL;
  md.count = 0;
  md.desired = 0;
  if (sep != NULL) {
    BioseqExplore (sep, (Pointer) (&md), FindNucSeqEntryCallback);
  }
  return md.sep;
}

typedef struct kinddata {
  Boolean  hasNuc;
  Boolean  hasProt;
} KindData, PNTR KindPtr;

static void HasNucOrProtCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  KindPtr    kptr;

  if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL && mydata != NULL) {
    kptr = (KindPtr) mydata;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_na (bsp->mol)) {
      kptr->hasNuc = TRUE;
    } else if (ISA_aa (bsp->mol)) {
      kptr->hasProt = TRUE;
    }
  }
}

NLM_EXTERN Boolean LIBCALL SeqEntryHasNucs (SeqEntryPtr sep)

{
  KindData  kd;

  kd.hasNuc = FALSE;
  kd.hasProt = FALSE;
  if (sep != NULL) {
    BioseqExplore (sep, (Pointer) (&kd), HasNucOrProtCallback);
  }
  return kd.hasNuc;
}

NLM_EXTERN Boolean LIBCALL SeqEntryHasProts (SeqEntryPtr sep)

{
  KindData  kd;

  kd.hasNuc = FALSE;
  kd.hasProt = FALSE;
  if (sep != NULL) {
    BioseqExplore (sep, (Pointer) (&kd), HasNucOrProtCallback);
  }
  return kd.hasProt;
}

static Boolean CheckForAlignments (GatherContextPtr gcp)

{
  BoolPtr  boolptr;

  if (gcp == NULL) return TRUE;

  boolptr = (BoolPtr) gcp->userdata;
  if (boolptr == NULL ) return TRUE;

  switch (gcp->thistype) {
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      *boolptr = TRUE;
      return TRUE;
    default :
      break;
  }
  return TRUE;
}

NLM_EXTERN Boolean LIBCALL SeqEntryHasAligns (Uint2 entityID, SeqEntryPtr sep)

{
  GatherScope  gs;
  Boolean      rsult;

  rsult = FALSE;
  if (entityID == 0 || sep == NULL) return FALSE;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQALIGN] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQHIST] = FALSE;
  gs.ignore[OBJ_SEQHIST_ALIGN] = FALSE;
  gs.scope = sep;
  GatherEntity (entityID, (Pointer) (&rsult), CheckForAlignments, &gs);
  return rsult;
}

static void FindPowerBLASTAsnCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  AnnotDescrPtr  desc;
  ObjectIdPtr    oip;
  SeqAnnotPtr    sap;
  BoolPtr        rsult;

  if (sep == NULL || sep->data.ptrvalue == NULL || mydata == NULL) return;
  rsult = (BoolPtr) mydata;
  sap = (IS_Bioseq (sep)) ?
         ((BioseqPtr) sep->data.ptrvalue)->annot :
         ((BioseqSetPtr) sep->data.ptrvalue)->annot;
  while (sap != NULL) {
    if (sap->type == 2) {
      desc = NULL;
      while ((desc = ValNodeFindNext (sap->desc, desc, Annot_descr_user)) != NULL) {
        if (desc->data.ptrvalue != NULL) {
          oip = ((UserObjectPtr) desc->data.ptrvalue)->type;
          if (oip != NULL && StringCmp (oip->str, "Hist Seqalign") == 0) {
            *rsult = TRUE;
          }
        }
      }
    }
    sap = sap->next;
  }
}

NLM_EXTERN Boolean LIBCALL PowerBLASTASN1Detected (SeqEntryPtr sep)

{
  Boolean  rsult;

  rsult = FALSE;
  SeqEntryExplore (sep, (Pointer) &rsult, FindPowerBLASTAsnCallback);
  return rsult;
}

NLM_EXTERN SeqLocPtr CreateWholeInterval (SeqEntryPtr sep)

{
  BioseqPtr  bsp;
  SeqIntPtr  sip;
  SeqLocPtr  slp;

  slp = NULL;
  if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    slp = ValNodeNew (NULL);
    if (slp != NULL) {
      sip = SeqIntNew ();
      if (sip != NULL) {
        slp->choice = SEQLOC_INT;
        slp->data.ptrvalue = (Pointer) sip;
        sip->from = 0;
        sip->to = bsp->length - 1;
        sip->strand = Seq_strand_plus;
        sip->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
      }
    }
  }
  return slp;
}

NLM_EXTERN void FreeAllFuzz (SeqLocPtr location)

{
  SeqIntPtr  sip;
  SeqLocPtr  slp;

  if (location == NULL) return;
  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_INT) {
      sip = (SeqIntPtr) slp->data.ptrvalue;
      if (sip != NULL) {
        sip->if_to = IntFuzzFree (sip->if_to);
        sip->if_from = IntFuzzFree (sip->if_from);
      }
    }
    slp = SeqLocFindNext (location, slp);
  }
}

NLM_EXTERN Boolean LocationHasNullsBetween (SeqLocPtr location)

{
  SeqLocPtr  slp;

  if (location == NULL) return FALSE;
  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_NULL) return TRUE;
    slp = SeqLocFindNext (location, slp);
  }
  return FALSE;
}

NLM_EXTERN Uint2 FindFeatFromFeatDefType (Uint2 subtype)

{
  switch (subtype) {
    case FEATDEF_GENE :
      return SEQFEAT_GENE;
    case FEATDEF_ORG :
      return SEQFEAT_ORG;
    case FEATDEF_CDS :
      return SEQFEAT_CDREGION;
    case FEATDEF_PROT :
      return SEQFEAT_PROT;
    case FEATDEF_PUB :
      return SEQFEAT_PUB;
    case FEATDEF_SEQ :
      return SEQFEAT_SEQ;
    case FEATDEF_REGION :
      return SEQFEAT_REGION;
    case FEATDEF_COMMENT :
      return SEQFEAT_COMMENT;
    case FEATDEF_BOND :
      return SEQFEAT_BOND;
    case FEATDEF_SITE :
      return SEQFEAT_SITE;
    case FEATDEF_RSITE :
      return SEQFEAT_RSITE;
    case FEATDEF_USER :
      return SEQFEAT_USER;
    case FEATDEF_TXINIT :
      return SEQFEAT_TXINIT;
    case FEATDEF_NUM :
      return SEQFEAT_NUM;
    case FEATDEF_PSEC_STR :
      return SEQFEAT_PSEC_STR;
    case FEATDEF_NON_STD_RESIDUE :
      return SEQFEAT_NON_STD_RESIDUE;
    case FEATDEF_HET :
      return SEQFEAT_HET;
    case FEATDEF_BIOSRC :
      return SEQFEAT_BIOSRC;
    default :
      if (subtype >= FEATDEF_preRNA && subtype <= FEATDEF_otherRNA) {
        return SEQFEAT_RNA;
      }
      if (subtype >= FEATDEF_IMP && subtype <= FEATDEF_site_ref) {
        return SEQFEAT_IMP;
      }
      if (subtype >= FEATDEF_preprotein && subtype <= FEATDEF_transit_peptide_aa) {
        return SEQFEAT_PROT;
      }
  }
  return 0;
}

NLM_EXTERN SeqIdPtr MakeSeqID (CharPtr str)

{
  long int     num;
  ObjectIdPtr  oid;
  SeqIdPtr     sip;
  CharPtr      tmp;

  sip = NULL;
  if (str != NULL && *str != '\0') {
    if (StringChr (str, '|') != NULL) {
      sip = SeqIdParse (str);
    } else {
      oid = ObjectIdNew ();
      if (oid != NULL) {
        for (tmp = str; *tmp != '\0'; tmp++) {
          if (! IS_DIGIT (*tmp)) {
            oid->str = StringSave (str);
            break;
          }
        }
        if (oid->str == NULL) {
          sscanf (str, "%ld", &num);
          oid->id = (Int4) num;
        }
        sip = ValNodeNew (NULL);
        if (sip != NULL) {
          sip->choice = SEQID_LOCAL;
          sip->data.ptrvalue = (Pointer) oid;
        } else {
          ObjectIdFree (oid);
        }
      }
    }
  }
  return sip;
}

NLM_EXTERN SeqIdPtr MakeUniqueSeqID (CharPtr prefix)

{
	Char buf[40];
	CharPtr tmp;
	Int2 ctr;
	ValNodePtr newid;
	ObjectIdPtr oid;
	ValNode vn;
	TextSeqId tsi;
	ValNodePtr altid;
	size_t len;

	altid = &vn;
	vn.choice = SEQID_GENBANK;
	vn.next = NULL;
	vn.data.ptrvalue = &tsi;
	tsi.name = NULL;
	tsi.accession = NULL;
	tsi.release = NULL;
	tsi.version = INT2_MIN;

	len = StringLen (prefix);
	if (len > 0 && len < 32) {
		tmp = StringMove(buf, prefix);
	} else {
		tmp = StringMove(buf, "tmpseq_");
	}

	newid = ValNodeNew(NULL);
	oid = ObjectIdNew();
	oid->str = buf;   /* allocate this later */
	newid->choice = SEQID_LOCAL;
	newid->data.ptrvalue = oid;

	tsi.name = buf;   /* check for alternative form */

	for (ctr = 1; ctr < 32000; ctr++)
	{
		sprintf(tmp, "%d", (int)ctr);
		if ((BioseqFindCore(newid) == NULL) && (BioseqFindCore(altid) == NULL))
		{
			oid->str = StringSave(buf);
			return newid;
		}
	}

	return NULL;
}

#define NUM_ORDER 16

NLM_EXTERN SeqIdPtr SeqIdFindWorst (SeqIdPtr sip)

{
  Uint1  order [NUM_ORDER];

  SeqIdBestRank (order, NUM_ORDER);
  order [SEQID_LOCAL] = 10;
  order [SEQID_GENBANK] = 5;
  order [SEQID_EMBL] = 5;
  order [SEQID_PIR] = 5;
  order [SEQID_SWISSPROT] = 5;
  order [SEQID_DDBJ] = 5;
  order [SEQID_PRF] = 5;
  order [SEQID_PDB] = 5;
  order [SEQID_PATENT] = 10;
  order [SEQID_OTHER] = 8;
  order [SEQID_GENERAL] = 15;
  order [SEQID_GIBBSQ] = 15;
  order [SEQID_GIBBMT] = 15;
  order [SEQID_GIIM] = 20;
  order [SEQID_GI] = 20;
  return SeqIdSelect (sip, order, NUM_ORDER);
}

NLM_EXTERN SeqFeatPtr CreateNewFeature (SeqEntryPtr sep, SeqEntryPtr placeHere,
                             Uint1 choice, SeqFeatPtr useThis)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqFeatPtr    prev;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->choice != 1) return NULL;
  sfp = NULL;
  bsp = NULL;
  bssp = NULL;
  if (placeHere == NULL) {
    placeHere = sep;
  }
  if (placeHere != NULL && placeHere->data.ptrvalue != NULL) {
    if (placeHere->choice == 1) {
      bsp = (BioseqPtr) placeHere->data.ptrvalue;
      sap = bsp->annot;
      while (sap != NULL && (sap->name != NULL || sap->desc != NULL || sap->type != 1)) {
        sap = sap->next;
      }
      if (sap == NULL) {
        sap = SeqAnnotNew ();
        if (sap != NULL) {
          sap->type = 1;
          sap->next = bsp->annot;
          bsp->annot = sap;
        }
        sap = bsp->annot;
      }
    } else if (placeHere->choice == 2) {
      bssp = (BioseqSetPtr) placeHere->data.ptrvalue;
      sap = bssp->annot;
      while (sap != NULL && (sap->name != NULL || sap->desc != NULL || sap->type != 1)) {
        sap = sap->next;
      }
      if (sap == NULL) {
        sap = SeqAnnotNew ();
        if (sap != NULL) {
          sap->type = 1;
          sap->next = bssp->annot;
          bssp->annot = sap;
        }
        sap = bssp->annot;
      }
    } else {
      return NULL;
    }
    if (sap != NULL) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (useThis != NULL) {
        sfp = useThis;
      } else {
        sfp = SeqFeatNew ();
      }
      if (sap->data != NULL) {
        prev = sap->data;
        while (prev->next != NULL) {
          prev = prev->next;
        }
        prev->next = sfp;
      } else {
        sap->data = (Pointer) sfp;
      }
      if (sfp != NULL) {
        sfp->data.choice = choice;
        if (useThis == NULL) {
          sfp->location = CreateWholeInterval (sep);
        }
      }
    }
  }
  return sfp;
}

NLM_EXTERN SeqFeatPtr CreateNewFeatureOnBioseq (BioseqPtr bsp, Uint1 choice, SeqLocPtr slp)

{
  SeqEntryPtr  sep;
  SeqFeatPtr   sfp;

  if (bsp == NULL) return NULL;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return NULL;
  sfp = CreateNewFeature (sep, NULL, choice, NULL);
  if (sfp == NULL) return NULL;
  if (slp != NULL) {
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = AsnIoMemCopy (slp, (AsnReadFunc) SeqLocAsnRead,
                                  (AsnWriteFunc) SeqLocAsnWrite);
  }
  return sfp;
}

NLM_EXTERN ValNodePtr CreateNewDescriptor (SeqEntryPtr sep, Uint1 choice)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Uint1         _class;
  ValNodePtr    descr;
  SeqEntryPtr   seqentry;
  ValNodePtr    vnp;

  vnp = NULL;
  if (sep != NULL) {
    descr = NULL;
    vnp = NULL;
    bsp = NULL;
    bssp = NULL;
    seqentry = sep;
    while (seqentry != NULL) {
      if (seqentry->choice == 1) {
        bsp = (BioseqPtr) seqentry->data.ptrvalue;
        if (bsp != NULL) {
          descr = bsp->descr;
          vnp = ValNodeNew (descr);
          if (descr == NULL) {
            bsp->descr = vnp;
          }
        }
        seqentry = NULL;
      } else if (seqentry->choice == 2) {
        bssp = (BioseqSetPtr) seqentry->data.ptrvalue;
        if (bssp != NULL) {
          _class = bssp->_class;
          if (_class == 7) {
            descr = bssp->descr;
            vnp = ValNodeNew (descr);
            if (descr == NULL) {
              bssp->descr = vnp;
            }
            seqentry = NULL;
          } else if ((_class >= 5 && _class <= 8) || _class == 11 /* || _class == BioseqseqSet_class_gen_prod_set */) {
            seqentry = bssp->seq_set;
          } else {
            descr = bssp->descr;
            vnp = ValNodeNew (descr);
            if (descr == NULL) {
              bssp->descr = vnp;
            }
            seqentry = NULL;
          }
        } else {
          seqentry = NULL;
        }
      } else {
        seqentry = NULL;
      }
    }
    if (vnp != NULL) {
      vnp->choice = choice;
    }
  }
  return vnp;
}

NLM_EXTERN ValNodePtr CreateNewDescriptorOnBioseq (BioseqPtr bsp, Uint1 choice)

{
  SeqEntryPtr  sep;

  if (bsp == NULL) return NULL;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return NULL;
  return CreateNewDescriptor (sep, choice);
}

