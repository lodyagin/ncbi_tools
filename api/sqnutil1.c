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
* $Revision: 6.2 $
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
  Uchar    ch;	/* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
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
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ' && ch != ';') {
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
    if (HasNoText (vnp->data.ptrvalue)) {
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
  GeneRefPtr  grp;

  if (StringICmp (gbq->qual, "pseudo") == 0) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp == NULL) return FALSE;
    grp->pseudo = TRUE;
    return TRUE;
  }
  return FALSE;
}

static Boolean HandledGBQualOnCDS (SeqFeatPtr sfp, GBQualPtr gbq)

{
  ProtRefPtr      prp;
  SeqFeatXrefPtr  xref;

  if (StringICmp (gbq->qual, "product") == 0) {
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
    if (prp == NULL) return FALSE;
    ValNodeCopyStr (&(prp->name), 0, gbq->val);
    return TRUE;
  }
  return FALSE;
}

static Boolean HandledGBQualOnRNA (SeqFeatPtr sfp, GBQualPtr gbq)

{
  RnaRefPtr  rrp;

  if (StringICmp (gbq->qual, "product") == 0 ||
      StringICmp (gbq->qual, "standard_name") == 0) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp == NULL) return FALSE;
    if (rrp->type == 0) {
      rrp->type = 255;
    }
    if (rrp->ext.choice != 0 && rrp->ext.choice != 1) return FALSE;
    if (! HasNoText (rrp->ext.value.ptrvalue)) return FALSE;
    if (rrp->ext.choice == 1 && rrp->ext.value.ptrvalue != NULL) {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    }
    if (rrp->ext.choice == 0 || rrp->ext.choice == 1) {
      rrp->ext.choice = 1;
      rrp->ext.value.ptrvalue = StringSave (gbq->val);
    }
  }
  return FALSE;
}

static Boolean HandledGBQualOnProt (SeqFeatPtr sfp, GBQualPtr gbq)

{
  ProtRefPtr  prp;
  ValNodePtr  vnp;

  if (StringICmp (gbq->qual, "product") == 0) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp == NULL) return FALSE;
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
  }
  return FALSE;
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

static void CleanupFeatureGBQuals (SeqFeatPtr sfp)

{
  DbtagPtr        db;
  GBQualPtr       gbq;
  size_t          len;
  GBQualPtr       nextqual;
  ObjectIdPtr     oip;
  GBQualPtr PNTR  prevqual;
  CharPtr         ptr;
  CharPtr         str;
  CharPtr         tag;
  Boolean         unlink;
  ValNodePtr      vnp;

  if (sfp == NULL) return;
  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    CleanVisString (&(gbq->qual));
    CleanVisString (&(gbq->val));
    if (gbq->qual == NULL) {
      gbq->qual = StringSave ("");
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
        sfp->exp_ev = 1;
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
    } else if (sfp->data.choice == SEQFEAT_GENE && HandledGBQualOnGene (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_CDREGION && HandledGBQualOnCDS (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_RNA && HandledGBQualOnRNA (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_PROT && HandledGBQualOnProt (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_IMP && HandledGBQualOnImp (sfp, gbq)) {
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

static void NormalizeAPub (ValNodePtr vnp)

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
      cgp->serial_number = -1; /* but does not remove if empty */
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

static void NormalizePubdesc (PubdescPtr pdp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (pdp == NULL) return;
  prev = &(pdp->pub);
  vnp = pdp->pub;
  while (vnp != NULL) {
    next = vnp->next;
    NormalizeAPub (vnp);
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

static void CleanupFeatureStrings (SeqFeatPtr sfp)

{
  BioSourcePtr  biop;
  GeneRefPtr    grp;
  ImpFeatPtr    ifp;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  PubdescPtr    pdp;
  ProtRefPtr    prp;
  RnaRefPtr     rrp;

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
        }
      }
      break;
    case SEQFEAT_PUB :
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      CleanVisString (&(pdp->comment));
      CleanDoubleQuote (pdp->comment);
      NormalizePubdesc (pdp);
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
      CleanSubSourceList (&(biop->subtype));
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
      onp = onp->next;
    }
  }
}

static void CleanupDescriptorStrings (ValNodePtr sdp)

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
      NormalizePubdesc (pdp);
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
      CleanSubSourceList (&(biop->subtype));
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
      onp = onp->next;
    }
  }
}

NLM_EXTERN void BasicSeqEntryCleanup (SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ImpFeatPtr    ifp;
  SeqAnnotPtr   sap = NULL;
  ValNodePtr    sdp = NULL;
  SeqFeatPtr    sfp;
  SeqEntryPtr   tmp;

  if (sep == NULL) return;
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
        CleanupFeatureStrings (sfp);
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    CleanupDescriptorStrings (sdp);
    sdp = sdp->next;
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
    if (slp->choice == SEQLOC_NULL)  return TRUE;
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

