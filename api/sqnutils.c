/*   sqnutils.c
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
* File Name:  sqnutils.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   9/2/97
*
* $Revision: 6.69 $
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

static Boolean SqnStringHasNoText (CharPtr str)

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
          if (! SqnStringHasNoText (tip->accession)) {
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

#ifdef OS_MAC
#define ASNLOAD_NEEDED 1
#endif
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
  if (! (asnFound || dataFound)) {
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
  order [SEQID_OTHER] = 10;
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

static CharPtr SqnTrimSpacesAroundString (CharPtr str)

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
      if (ch != ' ') {
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

static CharPtr SqnStringSave (CharPtr from)

{
  size_t  len;
  CharPtr to;

  to = NULL;
  len = StringLen (from);
  if (len > 0) {
    to = (CharPtr) MemGet (len + 1, FALSE);
    if (to != NULL) {
      MemCpy (to, from, len + 1);
      SqnTrimSpacesAroundString (to);
    }
  }
  return to;
}

NLM_EXTERN void UpdateLocalId (BioseqPtr bsp, CharPtr localId)

{
  Char         ch;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  SeqIdPtr     sip;
  long         val;

  if (bsp != NULL) {
    if (localId != NULL) {
      sip = bsp->id;
      while (sip != NULL && sip->choice != SEQID_LOCAL) {
        sip = sip->next;
      }
      oip = NULL;
      if (sip != NULL) {
        oip = (ObjectIdPtr) sip->data.ptrvalue;
      } else {
        sip = ValNodeNew (bsp->id);
        if (bsp->id == NULL) {
          bsp->id = sip;
        }
        if (sip != NULL) {
          oip = ObjectIdNew ();
          sip->choice = SEQID_LOCAL;
          sip->data.ptrvalue = (Pointer) oip;
        }
      }
      if (oip != NULL) {
        oip->str = MemFree (oip->str);
        if (sscanf (localId, "%ld", &val) == 1) {
          oip->id = (Int4) val;
        } else {
          oip->str = SqnStringSave (localId);
          ptr = oip->str;
          ch = *ptr;
          while (ch != '\0') {
            if (ch == '|') {
              *ptr = '~';
            }
            ptr++;
            ch = *ptr;
          }
        }
      }
      SeqMgrReplaceInBioseqIndex (bsp);
    }
  }
}

NLM_EXTERN void UpdateTitle (BioseqPtr bsp, CharPtr title)

{
  ValNodePtr  vnp;

  if (bsp != NULL) {
    if (title != NULL) {
      vnp = NULL;
      if (bsp->descr != NULL) {
        vnp = ValNodeFindNext (bsp->descr, NULL, Seq_descr_title);
      }
      if (vnp == NULL) {
        vnp = ValNodeNew (bsp->descr);
        if (vnp != NULL) {
          vnp->choice = Seq_descr_title;
        }
        if (bsp->descr == NULL) {
          bsp->descr = vnp;
        }
      }
      if (vnp != NULL) {
        vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        vnp->data.ptrvalue = SqnStringSave (title);
      }
    }
  }
}

NLM_EXTERN GeneRefPtr CreateNewGeneRef (CharPtr locus, CharPtr allele,
                             CharPtr desc, Boolean pseudo)

{
  GeneRefPtr  geneRef;

  geneRef = GeneRefNew ();
  if (geneRef != NULL) {
    geneRef->locus = SqnStringSave (locus);
    geneRef->allele = SqnStringSave (allele);
    geneRef->desc = SqnStringSave (desc);
    geneRef->pseudo = pseudo;
    if (geneRef->locus == NULL && geneRef->allele == NULL && geneRef->desc == NULL) {
      geneRef = GeneRefFree (geneRef);
    }
  }
  return geneRef;
}

NLM_EXTERN ProtRefPtr CreateNewProtRef (CharPtr name, CharPtr desc,
                             CharPtr ec, CharPtr activity)

{
  ProtRefPtr  protRef;
  ValNodePtr  vnp;

  protRef = ProtRefNew ();
  if (protRef != NULL) {
    if (name != NULL && *name != '\0') {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->data.ptrvalue = SqnStringSave (name);
        protRef->name = vnp;
      }
    }
    protRef->desc = SqnStringSave (desc);
    if (ec != NULL && *ec != '\0') {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->data.ptrvalue = SqnStringSave (ec);
        protRef->ec = vnp;
      }
    }
    if (activity != NULL && *activity != '\0') {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->data.ptrvalue = SqnStringSave (activity);
        protRef->activity = vnp;
      }
    }
    if (protRef->name == NULL && protRef->desc == NULL &&
        protRef->ec == NULL && protRef->activity == NULL) {
      protRef = ProtRefFree (protRef);
    }
  }
  return protRef;
}

NLM_EXTERN CdRegionPtr CreateNewCdRgn (Int2 frame, Boolean orf, Int2 genCode)

{
  CdRegionPtr  cdRgn;
  ValNodePtr   code;
  ValNodePtr   vnp;

  cdRgn = CdRegionNew ();
  if (cdRgn != NULL) {
    cdRgn->orf = orf;
    cdRgn->conflict = FALSE;
    cdRgn->frame = frame;
    cdRgn->gaps = 0;
    cdRgn->mismatch = 0;
    cdRgn->stops = 0;
    code = ValNodeNew (NULL);
    if (code != NULL) {
      code->choice = 254;
      vnp = ValNodeNew (NULL);
      code->data.ptrvalue = vnp;
      if (vnp != NULL) {
        vnp->choice = 2;
        vnp->data.intvalue = (Int4) genCode;
      }
    }
    cdRgn->genetic_code = code;
    cdRgn->code_break = NULL;
  }
  return cdRgn;
}

NLM_EXTERN void SetSeqFeatData (SeqFeatPtr sfp, Pointer data)

{
  if (sfp != NULL) {
    sfp->data.value.ptrvalue = (Pointer) data;
  }
}

NLM_EXTERN void SetSeqFeatProduct (SeqFeatPtr sfp, BioseqPtr bsp)

{
  ValNodePtr  slp;

  if (sfp != NULL) {
    sfp->product = SeqLocFree (sfp->product);
    if (bsp != NULL && bsp->id != NULL) {
      slp = ValNodeNew (NULL);
      if (slp != NULL) {
        slp->choice = 3;
        slp->data.ptrvalue = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
      }
      sfp->product = slp;
    }
  }
}

NLM_EXTERN void ResetSeqFeatInterval (SeqFeatPtr sfp)

{
  if (sfp != NULL) {
    sfp->location = SeqLocFree (sfp->location);
  }
}

NLM_EXTERN void AddSeqFeatInterval (SeqFeatPtr sfp, BioseqPtr bsp, Int4 from,
                         Int4 to, Boolean partial5, Boolean partial3)

{
  Int2  fuzz_from;
  Int2  fuzz_to;
  Int2  strand;
  Int4  tmp;

  if (sfp != NULL && bsp != NULL) {
    strand = Seq_strand_plus;
    if (from > to) {
      tmp = from;
      from = to;
      to = tmp;
      strand = Seq_strand_minus;
    }
    fuzz_from = -1;
    fuzz_to = -1;
    if (partial5) {
      fuzz_from = 2;
    }
    if (partial3) {
      fuzz_to = 1;
    }
    AddIntToSeqFeat (sfp, from - 1, to - 1, bsp, fuzz_from, fuzz_to, strand);
  }
}

NLM_EXTERN void AddSeqFeatPoint (SeqFeatPtr sfp, BioseqPtr bsp, Int4 location,
                      Boolean fuzz_before, Boolean fuzz_after, Int2 strand)

{
	SeqLocPtr slp, tmp, tmp2;
	SeqPntPtr spp;
	IntFuzzPtr ifp;
	Int2 fuzz;

	spp = SeqPntNew();
	spp->point = location - 1;
	spp->id = SeqIdDup(SeqIdFindBest(bsp->id, 0));
	spp->strand = (Uint1)strand;

	fuzz = -1;
    if (fuzz_before) {
      fuzz = 4;        /* tl */
    } else if (fuzz_after) {
      fuzz = 3;        /* tr */
    }
	if (fuzz >= 0)
	{
		ifp = IntFuzzNew();
		ifp->choice = 4;   /* lim */
		ifp->a = (Int4)fuzz;
		spp->fuzz = ifp;
	}

	slp = ValNodeNew(NULL);
	slp->choice = SEQLOC_PNT;
	slp->data.ptrvalue = (Pointer)spp;

	if (sfp->location == NULL)
	{
		sfp->location = slp;
		return;
	}

	tmp = sfp->location;
	if (tmp->choice == SEQLOC_MIX)   /* second one already */
	{
		tmp2 = (ValNodePtr)(tmp->data.ptrvalue);
		while (tmp2->next != NULL)
			tmp2 = tmp2->next;
		tmp2->next = slp;
	}
	else                             /* create a chain */
	{
		tmp2 = ValNodeNew(NULL);
		tmp2->choice = SEQLOC_MIX;
		tmp2->data.ptrvalue = (Pointer)tmp;
		tmp->next = slp;
		sfp->location = tmp2;
	}
}

typedef struct seqlocrange {
  Int4		left;
  Int4		right;
  Uint1		strand;
  struct seqlocrange PNTR next;
 } SeqLocRange, PNTR SeqLocRangePtr;
 
static SeqLocRangePtr SeqLocRangeFree (SeqLocRangePtr slrp)

{
  SeqLocRangePtr  next;

  while (slrp != NULL) {
    next = slrp->next;
    MemFree (slrp);
    slrp = next;
  }
  return NULL;
}

static SeqLocRangePtr CollectRanges (BioseqPtr target, SeqLocPtr slp)

{
  SeqLocPtr       curr;
  SeqLocRangePtr  head;
  SeqLocRangePtr  last;
  SeqLocRangePtr  slrp;

  head = NULL;
  last = NULL;
  curr = SeqLocFindNext (slp, NULL);
  while (curr != NULL) {
    if (curr->choice != SEQLOC_NULL) {
      slrp = MemNew (sizeof (SeqLocRange));
      if (slrp != NULL) {
        slrp->left = GetOffsetInBioseq (curr, target, SEQLOC_LEFT_END);
        slrp->right = GetOffsetInBioseq (curr, target, SEQLOC_RIGHT_END);
        slrp->strand = SeqLocStrand (curr);
        if (head == NULL) {
          head = slrp;
        } else if (last != NULL) {
          last->next = slrp;
        } else {
          ErrPostEx (SEV_ERROR, 0, 0, "SeqLocMerge list problem");
          SeqLocRangeFree (head);
          return NULL;
        }
        last = slrp;
      }
    }
    curr = SeqLocFindNext (slp, curr);
  }
  return head;
}

static int LIBCALLBACK CompareRanges (VoidPtr ptr1, VoidPtr ptr2)

{
  SeqLocRangePtr   slrp1;
  SeqLocRangePtr   slrp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    slrp1 = *((SeqLocRangePtr PNTR) ptr1);
    slrp2 = *((SeqLocRangePtr PNTR) ptr2);
    if (slrp1 != NULL && slrp2 != NULL) {
      if (slrp1->left > slrp2->left) {
        return 1;
      } else if (slrp1->left < slrp2->left) {
        return -1;
      } else if (slrp1->right > slrp2->right) {
        return 1;
      } else if (slrp1->right < slrp2->right) {
        return -1;
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static int LIBCALLBACK CompareReverseRanges (VoidPtr ptr1, VoidPtr ptr2)

{
  return (0 - CompareRanges (ptr1, ptr2));
}

static SeqLocRangePtr SortRanges (SeqLocRangePtr list, Boolean reverse)

{
  size_t          count;
  SeqLocRangePtr  PNTR head;
  size_t          i;
  SeqLocRangePtr  tmp;

  if (list != NULL) {
    count = 0;
    tmp = list;
    while (tmp != NULL) {
      count++;
      tmp = tmp->next;
    }
    if (count > 0) {
      head = MemNew ((count + 1) * sizeof (SeqLocRangePtr));
      if (head != NULL) {
        tmp = list;
        i = 0;
        while (tmp != NULL && i < count) {
          head [i] = tmp;
          tmp = tmp->next;
          i++;
        }
        if (reverse) {
          HeapSort (head, count, sizeof (SeqLocRangePtr), CompareReverseRanges);
        } else {
          HeapSort (head, count, sizeof (SeqLocRangePtr), CompareRanges);
        }
        for (i = 0; i < count; i++) {
          tmp = head [i];
          tmp->next = head [i + 1];
        }
        list = head [0];
        MemFree (head);
      }
    }
  }
  return list;
}

static SeqLocRangePtr MergeOverlaps (SeqLocRangePtr list, Boolean fuse_joints)

{
  SeqLocRangePtr  last;
  SeqLocRangePtr  next;
  SeqLocRangePtr  this;

  if (list != NULL) {
    this = list->next;
    last = list;
    while (this != NULL) {
      next = this->next;
      if (this->left < last->right) {
        last->right = MAX (this->right, last->right);
        MemFree (this);
        last->next = next;
      } else if (fuse_joints && this->left == last->right) {
        last->right = MAX (this->right, last->right);
        MemFree (this);
        last->next = next;
      } else {
        last = this;
      }
      this = next;
    }
  }
  return list;
}

static SeqLocPtr SeqLocFromRange (SeqLocRangePtr head, BioseqPtr target,
                                  Boolean partial5, Boolean partial3,
                                  Boolean add_null)

{
  SeqLocPtr   firstSlp;
  Int4        from;
  Int2        fuzz_from;
  Int2        fuzz_to;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  Boolean     notFirst;
  SeqFeatPtr  sfp;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  Int2        strand;
  Int4        tmp;
  SeqLocPtr   tmploc1;
  SeqLocPtr   tmploc2;
  Int4        to;

  if (head == NULL) return NULL;
  slp = NULL;
  sfp = SeqFeatNew ();
  if (sfp != NULL) {
    notFirst = FALSE;
    while (head != NULL) {
      fuzz_from = -1;
      fuzz_to = -1;
      from = head->left;
      to = head->right;
      strand = head->strand;
      if (from > to) {
        tmp = from;
        from = to;
        to = tmp;
      }
      if (add_null  && notFirst) {
        slp = ValNodeNew (NULL);
        if (slp != NULL) {
          slp->choice = SEQLOC_NULL;
          tmploc1 = sfp->location;
          if (tmploc1 != NULL) {
            if (tmploc1->choice == SEQLOC_MIX) {
              tmploc2 = (ValNodePtr) (tmploc1->data.ptrvalue);
              if (tmploc2 != NULL) {
                while (tmploc2->next != NULL) {
                  tmploc2 = tmploc2->next;
                }
                tmploc2->next = slp;
              }
            } else {
              tmploc2 = ValNodeNew (NULL);
              if (tmploc2 != NULL) {
                tmploc2->choice = SEQLOC_MIX;
                tmploc2->data.ptrvalue = (Pointer) tmploc1;
                tmploc1->next = slp;
                sfp->location = tmploc2;
              }
            }
          }
        }
      }
      AddIntToSeqFeat (sfp, from, to, target,
                       fuzz_from, fuzz_to, strand);
      notFirst = TRUE;
      head = head->next;
    }
    firstSlp = NULL;
    lastSlp = NULL;
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      if (firstSlp == NULL) {
        firstSlp = slp;
      }
      lastSlp = slp;
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (firstSlp != NULL && firstSlp->choice == SEQLOC_INT &&
        firstSlp->data.ptrvalue != NULL && partial5) {
      sip = (SeqIntPtr) firstSlp->data.ptrvalue;
      ifp = IntFuzzNew ();
      if (ifp != NULL) {
        ifp->choice = 4;
        if (sip->strand == Seq_strand_minus ||
            sip->strand == Seq_strand_both_rev) {
          sip->if_to = ifp;
          ifp->a = 1;
        } else {
          sip->if_from = ifp;
          ifp->a = 2;
        }
      }
    }
    if (lastSlp != NULL && lastSlp->choice == SEQLOC_INT &&
        lastSlp->data.ptrvalue != NULL && partial3) {
      sip = (SeqIntPtr) lastSlp->data.ptrvalue;
      ifp = IntFuzzNew ();
      if (ifp != NULL) {
        ifp->choice = 4;
        if (sip->strand == Seq_strand_minus ||
            sip->strand == Seq_strand_both_rev) {
          sip->if_from = ifp;
          ifp->a = 2;
        } else {
          sip->if_to = ifp;
          ifp->a = 1;
        }
      }
    }
    slp = sfp->location;
    sfp->location = NULL;
    SeqFeatFree (sfp);
  }
  return slp;
}

NLM_EXTERN SeqLocPtr SeqLocMerge (BioseqPtr target, SeqLocPtr to, SeqLocPtr from,
                       Boolean single_interval, Boolean fuse_joints,
                       Boolean add_null)

{
  SeqLocRangePtr  head;
  SeqLocRangePtr  last;
  Boolean         partial5;
  Boolean         partial3;
  SeqLocPtr       slp;
  Int2            strand;

  if (target == NULL) return NULL;
  if (to == NULL && from == NULL) return NULL;

  slp = NULL;
  partial5 = FALSE;
  partial3 = FALSE;
  head = CollectRanges (target, to);
  if (head == NULL) {
    head = CollectRanges (target, from);
  } else {
    last = head;
    while (last->next != NULL) {
      last = last->next;
    }
    last->next = CollectRanges (target, from);
  }
  if (head != NULL) {
    strand = head->strand;
    head = SortRanges (head, FALSE);
    head = MergeOverlaps (head, fuse_joints);
    if (single_interval) {
      last = head;
      while (last->next != NULL) {
        last = last->next;
      }
      head->left = MIN (head->left, last->left);
      head->right = MAX (head->right, last->right);
      MemFree (head->next);
      head->next = NULL;
    }
    last = head;
    while (last != NULL) {
      last->strand = strand;
      last = last->next;
    }
    if (strand == Seq_strand_minus) {
      head = SortRanges (head, TRUE);
    }
    slp = SeqLocFromRange (head, target, partial5, partial3, add_null);
    head = SeqLocRangeFree (head);
  }
  return slp;
}

NLM_EXTERN Boolean SeqLocBadSortOrder (BioseqPtr bsp, SeqLocPtr slp)

{
  SeqLocRangePtr  curr;
  SeqLocRangePtr  head;
  SeqLocRangePtr  last;

  if (bsp == NULL || slp == NULL) return FALSE;
  if (SeqLocCheck (slp) == SEQLOCCHECK_WARNING) return FALSE;
  /*
  if (SeqLocId (slp) == NULL) return FALSE;
  */
  head = CollectRanges (bsp, slp);
  if (head == NULL || head->next == NULL) return FALSE;
  last = head;
  curr = head->next;
  while (curr != NULL) {
    if (curr->strand == Seq_strand_minus) {
      if (last->right < curr->right) return TRUE;
    } else {
      if (last->left > curr->left) return TRUE;
    }
    last = curr;
    curr = curr->next;
  }
  return FALSE;
}

static void ConvertToFeatsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent, Boolean toProts)

{
  BioseqPtr   bsp;
  SeqFeatPtr  sfp;
  ValNodePtr  vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (ISA_aa (bsp->mol) && (! toProts)) return;
  vnp = (ValNodePtr) mydata;
  if (vnp == NULL) return;
  while (vnp != NULL) {
    switch (vnp->choice) {
      case Seq_descr_pub :
        sfp = CreateNewFeature (sep, NULL, SEQFEAT_PUB, NULL);
        if (sfp != NULL) {
          sfp->data.value.ptrvalue = AsnIoMemCopy ((Pointer) vnp->data.ptrvalue,
                                                   (AsnReadFunc) PubdescAsnRead,
                                                   (AsnWriteFunc) PubdescAsnWrite);
        }
        break;
      case Seq_descr_source :
        sfp = CreateNewFeature (sep, NULL, SEQFEAT_BIOSRC, NULL);
        if (sfp != NULL) {
          sfp->data.value.ptrvalue = AsnIoMemCopy ((Pointer) vnp->data.ptrvalue,
                                                   (AsnReadFunc) BioSourceAsnRead,
                                                   (AsnWriteFunc) BioSourceAsnWrite);
        }
        break;
      case Seq_descr_comment :
        sfp = CreateNewFeature (sep, NULL, SEQFEAT_COMMENT, NULL);
        if (sfp != NULL) {
          sfp->comment = StringSave ((CharPtr) vnp->data.ptrvalue);
        }
        break;
      default :
        break;
    }
    vnp = vnp->next;
  }
}

static void ConvertToFeatsOnNucsAndProts (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  ConvertToFeatsCallback (sep, mydata, index, indent, TRUE);
}

static void ConvertToFeatsOnNucsOnly (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  ConvertToFeatsCallback (sep, mydata, index, indent, FALSE);
}

NLM_EXTERN Boolean ConvertPubSrcComDescsToFeats (SeqEntryPtr sep, Boolean pub, Boolean src, Boolean com, Boolean toProts)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    head;
  ValNodePtr    nextsdp;
  Pointer PNTR  prevsdp;
  Boolean       rsult;
  ValNodePtr    sdp;

  rsult = FALSE;
  if (! (pub || src || com)) return FALSE;
  if (sep == NULL || sep->data.ptrvalue == NULL) return FALSE;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        if (ConvertPubSrcComDescsToFeats (sep, pub, src, com, toProts)) {
          rsult = TRUE;
        }
      }
      return rsult;
    }
  }
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return FALSE;
  head = NULL;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if ((sdp->choice == Seq_descr_pub && pub) ||
        (sdp->choice == Seq_descr_source && src) ||
        (sdp->choice == Seq_descr_comment && com)) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      ValNodeLink (&head, sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
  rsult = (head != NULL);
  if (toProts) {
    BioseqExplore (sep, head, ConvertToFeatsOnNucsAndProts);
  } else {
    BioseqExplore (sep, head, ConvertToFeatsOnNucsOnly);
  }
  for (sdp = head; sdp != NULL; sdp = sdp->next) {
    switch (sdp->choice) {
      case Seq_descr_pub :
        PubdescFree ((PubdescPtr) sdp->data.ptrvalue);
        break;
      case Seq_descr_source :
        BioSourceFree ((BioSourcePtr) sdp->data.ptrvalue);
        break;
      case Seq_descr_comment :
        MemFree (sdp->data.ptrvalue);
        break;
      default :
        break;
    }
  }
  ValNodeFree (head);
  return rsult;
}

/* from Colombe */

static CharPtr sqn_string_complement (CharPtr str)
{
  CharPtr strp;

  for (strp = str; *strp != '\0'; strp++) {
         if (*strp == 'A') *strp = 'T';
         else if (*strp == 'T') *strp = 'A';
         else if (*strp == 'C') *strp = 'G';
         else if (*strp == 'G') *strp = 'C';
  }
  *strp = '\0';
  return str;
}

static CharPtr sqn_string_reverse (CharPtr str)
{
  Char    car;
  Int4    j;
  Int4    k;

  j = 0;
  k = StringLen (str) - 1;
  while (j < k) {
    car = str[j]; str[j] = str[k]; str[k] = car;
    j++;
    k--;
  }
  return str;
}

static Int4 sqn_ReadBufferFromSep (SeqPortPtr spp, CharPtr buffer, Int4 from, Int4 to, Int4 buffsegstart)
{
  Uint1    residue;
  Int4     k;
  Int4     pos;

  SeqPortSeek (spp, from, SEEK_SET);
  k = buffsegstart;
  pos = from;
  residue = SeqPortGetResidue(spp);
  while (pos < to && residue != SEQPORT_EOF)
  {
    if ( ! IS_residue(residue)) {
      /*
      switch (residue)
      {  
           case SEQPORT_VIRT:
              Message(MSG_OK,"SEQPORT_VIRT [%d=%c] at %ld\n", (int)residue, (char)residue, (long)pos);
              break;
           case SEQPORT_EOS:
              Message(MSG_OK,"[EOS]\n");
              break;
           default:
              Message(MSG_OK,"unknown char\n");
              break;
      }  
      pos++;
      */
    } else {
      buffer[k] = (Char) residue;
      k++;  
      pos++;
    }
    residue = SeqPortGetResidue(spp);
  }
  buffer[k] = '\0';
  return k;
}
 
static CharPtr sqn_load_seq_data (SeqIdPtr sip, Int4 from, Int4 to, Boolean is_prot, Int4 *lenp)
{
  BioseqPtr        bsp;
  SeqLocPtr        slp;
  SeqPortPtr       spp;
  CharPtr          str = NULL;
  Int4             lens;

  if (from > -1 && to > -1 && from >= to)
     return NULL;
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
     if (from < 0 || from > bsp->length -1)
        from = 0;
     if (to < 0 || to > bsp->length -1)
        to = bsp->length -1;
     BioseqUnlock (bsp);
     slp = SeqLocIntNew (from, to, Seq_strand_plus, sip);
     if (is_prot)
        spp = SeqPortNewByLoc (slp, Seq_code_ncbistdaa);
     else
        spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
     if (spp != NULL) {
        str = MemNew ((to-from+4) * sizeof(Char));
        lens = sqn_ReadBufferFromSep (spp, str, 0, to -from +1, 0);
        SeqPortFree (spp);
        if (lenp != NULL)
           *lenp = lens;
     }   
     SeqLocFree (slp);
  }
  return str;
}

static Int4 getlengthforid (SeqIdPtr sip)
{
  BioseqPtr        bsp;
  Int4             lens=0;

  if (sip==NULL)
     return 0;
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
     lens = bsp->length;
     BioseqUnlock (bsp);
  }
  return lens;
}

NLM_EXTERN SeqLocPtr StringSearchInBioseq (SeqIdPtr sip, CharPtr sub)
{
  SeqLocPtr slp=NULL;
  CharPtr   strdb,
            strtmp;
  Int4      lenbsp,
            fromp, top;
  Int4      lensub,
            lens,
            offset,
            shiftRange,
            maxRange;
  Boolean   firstpass=TRUE;

  lensub = StringLen (sub);
  maxRange = (Int4) MAX ((Int4)1000, lensub);
  lenbsp = getlengthforid (sip);
  while (slp == NULL) 
  {
   fromp = 0;
   top = MIN ((Int4)(fromp+maxRange), lenbsp) -1;
   while (fromp <= lenbsp && slp == NULL)
   {
     strdb = sqn_load_seq_data (sip, fromp, top, FALSE, &lens); 
     if (strdb != NULL)
     {
        offset = 0;
        strtmp = StringISearch (strdb, sub);
        if (strtmp != NULL) {
           offset =(Int4)abs(abs((long)strdb)-abs((long)strtmp));
           offset += fromp;
           if (offset > 0) {
              if (firstpass)
                 slp = SeqLocIntNew (offset, offset+lensub-1, Seq_strand_plus, sip);
              else 
                 slp = SeqLocIntNew (offset, offset+lensub-1, Seq_strand_minus, sip);
           }
        }
        MemFree (strdb);
     }
     shiftRange = maxRange - lensub;
     fromp = fromp + shiftRange;
     top = MIN ((Int4)(fromp+maxRange), lenbsp);
   }
   if (!firstpass) {
      sub = sqn_string_complement (sub);
      sub = sqn_string_reverse (sub);
      break;
   }
   firstpass=FALSE;
   sub = sqn_string_complement (sub);
   sub = sqn_string_reverse (sub);
  }
  return slp;
}

/*****************************************************************************
*
*   SequinEntryList (sep, mydata, mycallback, index, indent)
*       traverses all Seq-entry nodes beginning with sep
*       calls mycallback () at each node
*       Does enter BioseqSets of _class "parts", but ignores the
*       parts set itself
*
*****************************************************************************/

NLM_EXTERN Int4 SequinEntryList (SeqEntryPtr sep, Pointer mydata,
                                 SeqEntryFunc mycallback,
                                 Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return index;
  if (IS_Bioseq (sep)) {
    if (mycallback != NULL)
      (*mycallback) (sep, mydata, index, indent);
    return index + 1;
  }
  /*
  if (Bioseq_set_class (sep) == 4) return index;
  index++;
  */
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  sep = bssp->seq_set;
  indent++;
  while (sep != NULL) {
    index = SequinEntryList (sep, mydata, mycallback, index, indent);
    sep = sep->next;
  }
  return index;
}

/* PHRAP file reading functions */

static CharPtr ReadALine (CharPtr str, size_t size, FILE *fp)

{
  Char     ch;
  CharPtr  ptr;
  CharPtr  rsult;

  if (str == NULL || size < 1 || fp == NULL) return NULL;
  *str = '\0';
  rsult = fgets (str, size, fp);
  if (rsult != NULL) {
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch != '\n' && ch != '\r') {
      ptr++;
      ch = *ptr;
    }
    *ptr = '\0';
  }
  return rsult;
}

static Boolean HasNoText (CharPtr str)

{
  Char  ch;

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

static SeqEntryPtr ReadPhrapDNA (FILE *fp, CharPtr id)

{
  ByteStorePtr  bs = NULL;
  BioseqPtr     bsp = NULL;
  Char          buf [256];
  Char          ch;
  Boolean       goOn = TRUE;
  CharPtr       p;
  CharPtr       q;
  SeqEntryPtr   sep = NULL;
  CharPtr       str;

  if (fp == NULL || HasNoText (id)) return NULL;
  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;
  bsp = BioseqNew ();
  if (bsp == NULL) return NULL;
  bs = BSNew (1000);
  if (bs == NULL) return NULL;

  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

  bsp->mol = Seq_mol_na;
  bsp->seq_data_type = Seq_code_iupacna;
  bsp->repr = Seq_repr_raw;
  bsp->length = 0;
  bsp->id = MakeSeqID (id);
  SeqMgrAddToBioseqIndex (bsp);

  goOn = TRUE;
  while (goOn) {
    str = ReadALine (buf, sizeof (buf), fp);
    if (HasNoText (str)) {
      goOn = FALSE;
    } else {
      p = str;
      q = str;
      ch = *p;
      while (ch != '\0') {
        if (! (IS_ALPHA (ch))) {
          p++;
        } else {
          ch = TO_UPPER (ch);
          if (ch == 'X') {
            ch = 'N';
          }
          *q = ch;
          p++;
          q++;
        }
        ch = *p;
      }
      *q = '\0';
      BSWrite (bs, (VoidPtr) str, (Int4) StringLen (str));
    }
  }

  bsp->seq_data = bs;
  bsp->length = BSLen (bs);

  BioseqPack (bsp);
  return sep;
}

static SeqGraphPtr ReadPhrapQuality (FILE *fp, BioseqPtr bsp)

{
  ByteStorePtr  bs = NULL;
  Char          buf [256];
  Uint1         bytes [128];
  Char          ch;
  Boolean       goOn = TRUE;
  Int2          i;
  Int2          max = INT2_MIN;
  Int2          min = INT2_MAX;
  CharPtr       p;
  CharPtr       q;
  SeqGraphPtr   sgp = NULL;
  SeqIntPtr     sintp;
  CharPtr       str;
  int           val;

  if (fp == NULL || bsp == NULL) return NULL;
  sgp = SeqGraphNew ();
  if (sgp == NULL) return NULL;
  bs = BSNew (1000);
  if (bs == NULL) return NULL;

  goOn = TRUE;
  while (goOn) {
    str = ReadALine (buf, sizeof (buf), fp);
    if (HasNoText (str)) {
      goOn = FALSE;
    } else {
      i = 0;
      p = str;
      ch = *p;
      while (ch != '\0') {
        while (IS_WHITESP (ch)) {
          p++;
          ch = *p;
        }
        q = p;
        ch = *q;
        while (IS_DIGIT (ch)) {
          q++;
          ch = *q;
        }
        *q = '\0';
        q++;

        if (*p != '\0') {
          if (sscanf (p, "%d", &val) == 1) {
            if (val < 0 || val > 255) {
              /* error */
              val = 0;
            }
            bytes [i] = (Uint1) val;
            i++;
            max = MAX (max, (Int2) val);
            min = MIN (min, (Int2) val);
          }
        }
        p = q;
        ch = *p;
      }
      if (i > 0) {
        BSWrite (bs, (Pointer) bytes, (Int4) i);
      }
    }
  }

  sgp->numval = BSLen (bs);
  BSPutByte (bs, EOF);
  sgp->title = StringSave ("Phrap Quality");
  if (bsp->length != sgp->numval) {
    sgp->flags [0] = 1;
    sgp->compl = (bsp->length) / sgp->numval;
  } else {
    sgp->flags [0] = 0;
    sgp->compl = 1;
  }
  sgp->flags [1] = 0;
  sgp->flags [2] = 3;
  sgp->axis.intvalue = 0;
  sgp->min.intvalue = min;
  sgp->max.intvalue = max;
  sgp->a = 1.0;
  sgp->b = 0;
  sgp->values = (Pointer) bs;

  sintp = SeqIntNew ();
  sintp->from = 0;
  sintp->to = bsp->length - 1;
  sintp->id = SeqIdDup (bsp->id);
  ValNodeAddPointer (&(sgp->loc), SEQLOC_INT, (Pointer) sintp);

  return sgp;
}

static Boolean PhrapSequenceHasClipping (FILE *fp)

{
  Char     buf [256];
  Boolean  goOn = TRUE;
  Boolean  rsult = FALSE;
  CharPtr  str;

  if (fp == NULL) return FALSE;
  goOn = TRUE;
  while (goOn) {
    str = ReadALine (buf, sizeof (buf), fp);
    if (HasNoText (str)) {
      goOn = FALSE;
    } else {
      if (StringNCmp (str, "Clipping", 8) == 0) {
        rsult = TRUE;
      }
    }
  }
  return rsult;
}

static CharPtr BioseqGetLocalIdStr (BioseqPtr bsp)

{
  ObjectIdPtr  oip;
  SeqIdPtr     sip;

  if (bsp == NULL) return NULL;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) sip->data.ptrvalue;
      if (oip != NULL && oip->str != NULL) {
        return oip->str;
      }
    }
  }
  return NULL;
}

static SeqAnnotPtr NewGraphSeqAnnot (CharPtr name, SeqGraphPtr sgp)

{
  SeqAnnotPtr  sap = NULL;

  if (sgp == NULL) return NULL;
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;

  if (! HasNoText (name)) {
    ValNodeAddPointer (&(sap->desc), Annot_descr_name, StringSave (name));
  }
  sap->type = 3;
  sap->data = (Pointer) sgp;

  return sap;
}

static CharPtr taglist [] = {
  "", "DNA", "CO", "BaseQuality", "BQ", "Sequence", NULL
};

/* Phrap reading function based on sample code supplied by C. Magness */
NLM_EXTERN SeqEntryPtr ReadPhrapFile (FILE *fp)

{
  BioseqPtr    bsp;
  Char         buf [256];
  Boolean      goOn = TRUE;
  SeqEntryPtr  head = NULL;
  Int2         i;
  SeqEntryPtr  lastsep;
  SeqGraphPtr  lastsgp;
  CharPtr      p;
  CharPtr      q;
  SeqAnnotPtr  sap;
  SeqEntryPtr  sep = NULL;
  SeqGraphPtr  sgp;
  CharPtr      str;
  Int2         tag;

  if (fp == NULL) return NULL;
  goOn = TRUE;
  while (goOn) {
    str = ReadALine (buf, sizeof (buf), fp);
    if (str == NULL) {
      goOn = FALSE;
    } else if (! HasNoText (str)) {
      p = StringChr (str, ' ');
      if (p != NULL) {
        *p = '\0';
        p++;
      }
      tag = 0;
      for (i = 0; taglist [i] != NULL; i++) {
        if (StringCmp (str, taglist [i]) == 0) {
          tag = i;
        }
      }
      if (tag != 0) {
        if (p != NULL) {
          q = StringChr (p, ' ');
          if (q != NULL) {
            *q = '\0';
          }
        }
        switch (tag) {
          case 1 :
          case 2 :
            if (p != NULL) {
              sep = ReadPhrapDNA (fp, p);
              ValNodeLink (&head, sep);
            }
            /* for new format, sep points to current sequence */
            break;
          case 3 :
            if (p != NULL) {
              sep = head;
              while (sep != NULL && StringCmp (p, BioseqGetLocalIdStr ((BioseqPtr) sep->data.ptrvalue)) != 0) {
                sep = sep->next;
              }
            }
            /* and flow through to case 4 */
          case 4 :
            if (sep != NULL) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              sgp = ReadPhrapQuality (fp, bsp);
              if (sgp != NULL) {
                for (sap = bsp->annot; sap != NULL; sap = sap->next) {
                  if (sap->type == 3) {
                    for (lastsgp = sap->data; lastsgp->next != NULL; lastsgp = lastsgp->next) {
                      continue;
                    }
                    lastsgp->next = sgp;
                    break;
                  }
                }
                if (sap == NULL) {
                  if (bsp->annot != NULL) {
                    for (sap = bsp->annot; sap->next != NULL; sap = sap->next) {
                      continue;
                    }
                    sap->next = NewGraphSeqAnnot ("Graphs", sgp);
                  } else {
                    bsp->annot = NewGraphSeqAnnot ("Graphs", sgp);
                  }
                }
              }
            }
            break;
          case 5 :
            /* unlinkes and removes sep if Clipping line present */
            if (p != NULL) {
              if (PhrapSequenceHasClipping (fp)) {
                sep = head;
                lastsep = NULL;
                while (sep != NULL && StringCmp (p, BioseqGetLocalIdStr ((BioseqPtr) sep->data.ptrvalue)) != 0) {
                  lastsep = sep;
                  sep = sep->next;
                }
                if (sep != NULL) {
                  if (lastsep != NULL) {
                    lastsep->next = sep->next;
                    sep->next = NULL;
                    SeqEntryFree (sep);
                  } else {
                    head = sep->next;
                    sep->next = NULL;
                    SeqEntryFree (sep);
                  }
                }
              }
            }
            break;
          default :
            break;
        }
      }
    }
  }
  return head;
}

static ValNodePtr ParseContigOrFeatureTableString (CharPtr contigs, Boolean tabDelimited)

{
  Char        ch;
  Int2        i, j, k;
  CharPtr     str;
  Char        tmp [1024];
  ValNodePtr  vnp;

  vnp = NULL;
  i = 0;
  while (StringLen (contigs + i) > 0) {
    str = contigs + i;
    k = 0;
    ch = str [k];
    while (ch == ' ') {
      k++;
      ch = str [k];
    }
    j = 0;
    if (tabDelimited) {
      while (ch != '\0' && ch != '\t') {
        j++;
        ch = str [j + k];
      }
    } else {
      while (ch != '\0' && ch != ',' && (! (IS_WHITESP (ch)))) {
        j++;
        ch = str [j + k];
      }
    }
    if (ch == '\0') {
      i += j + k;
    } else {
      str [j + k] = '\0';
      i += j + k + 1;
    }
    StringNCpy_0 (tmp, str + k, sizeof (tmp));
    SqnTrimSpacesAroundString (tmp);
    if (HasNoText (tmp)) {
      ValNodeAdd (&vnp);
    } else {
      ValNodeCopyStr (&vnp, 0, tmp);
    }
  }
  if (vnp != NULL) {
    vnp->choice = (Uint1) ValNodeLen (vnp);
  }
  return vnp;
}

NLM_EXTERN SeqEntryPtr SetPhrapContigOrder (SeqEntryPtr head, CharPtr contigs)

{
  SeqEntryPtr  sep, lastsep, nextsep, newhead;
  ValNodePtr   vnphead, vnp;

  if (head == NULL || contigs == NULL) return head;
  vnphead = ParseContigOrFeatureTableString (contigs, FALSE);
  if (vnphead == NULL) return head;
  newhead = NULL;
  for (vnp = vnphead; vnp != NULL; vnp = vnp->next) {
    sep = head;
    lastsep = NULL;
    while (sep != NULL && StringCmp ((CharPtr) vnp->data.ptrvalue,
           BioseqGetLocalIdStr ((BioseqPtr) sep->data.ptrvalue)) != 0) {
      lastsep = sep;
      sep = sep->next;
    }
    if (sep != NULL) {
      if (lastsep != NULL) {
        lastsep->next = sep->next;
        sep->next = NULL;
        ValNodeLink (&newhead, sep);
      } else {
        head = sep->next;
        sep->next = NULL;
        ValNodeLink (&newhead, sep);
      }
    }
  }
  for (sep = head; sep != NULL; sep = nextsep) {
    nextsep = sep->next;
    sep->next = NULL;
    SeqEntryFree (sep);
    sep = nextsep;
  }
  ValNodeFreeData (vnphead);
  return newhead;
}

/* Build contig section */
static void  AddNucToContig (CharPtr accnString, Int4 from, Int4 to,
                             Int4 size, Int2 strand, BioseqPtr segseq,
                             BoolPtr hasgaps, Boolean isgap)

{
  Boolean       allDigits;
  Char          ch;
  DbtagPtr      dp;
  CharPtr       ptr;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  TextSeqIdPtr  tsip;
  long int      val;

  slp = ValNodeNew ((ValNodePtr) segseq->seq_ext);
  if (slp == NULL) return;
  if (segseq->seq_ext == NULL) {
    segseq->seq_ext = (Pointer) slp;
  }

  sintp = SeqIntNew ();
  sintp->from = from;
  sintp->to = to;
  sintp->strand = (Uint1) strand;

  slp->choice = SEQLOC_INT;
  slp->data.ptrvalue = (Pointer) sintp;

  sip = ValNodeNew (NULL);
  if (isgap) {
    /* sip = MakeUniqueSeqID ("gap_"); */
    dp = DbtagNew ();
    dp->db = StringSave ("SeqLit");
    dp->tag = ObjectIdNew ();
    dp->tag->id = 0;
    dp->tag->str = NULL;
    sip->choice = SEQID_GENERAL;
    sip->data.ptrvalue = dp;
    if (hasgaps != NULL) {
      *hasgaps = TRUE;
    }
  } else {
    allDigits = TRUE;
    ptr = accnString;
    ch = *ptr;
    while (ch != '\0' && allDigits) {
      if (! IS_DIGIT (ch)) {
        allDigits = FALSE;
      }
      ptr++;
      ch = *ptr;
    }
    if (allDigits && sscanf (accnString, "%ld", &val) == 1) {
      sip->choice = (Uint1) SEQID_GI;
      sip->data.intvalue = val;
    } else {
      tsip = TextSeqIdNew ();
      tsip->accession = StringSave (accnString);
      sip->choice = (Uint1) SEQID_GENBANK;
      sip->data.ptrvalue = tsip;
    }
  }

  sintp->id = sip;

  segseq->length += size;
}

#define accnString field [0]
#define startString field [1]
#define stopString field [2]
#define sizeString field [3]
#define strandString field [4]

static void AdjustContigValues (ValNodePtr line)

{
  Int2        i;
  ValNodePtr  nextAccn;
  ValNodePtr  nextStart;
  long int    num;
  ValNodePtr  thisStop;
  Char        tmp [32];
  Int4        val;

  if (line == NULL) return;
  for (i = 0, thisStop = line->data.ptrvalue; i < 2 && thisStop != NULL; i++, thisStop = thisStop->next) {
    continue;
  }
  line = line->next;
  while (line != NULL && line->data.ptrvalue == NULL) {
    line = line->next;
  }
  if (line == NULL) {
    if (thisStop != NULL) {
      if (sscanf ((CharPtr) thisStop->data.ptrvalue, "%ld", &num) == 1) {
        val = (Int4) num;
        val++;
        sprintf (tmp, "%ld", (long) val);
        thisStop->data.ptrvalue = MemFree (thisStop->data.ptrvalue);
        thisStop->data.ptrvalue = StringSave (tmp);
      }
    }
    return;
  }
  nextAccn = line->data.ptrvalue;
  if (nextAccn != NULL && StringICmp (nextAccn->data.ptrvalue, "gap") == 0) return;
  for (i = 0, nextStart = line->data.ptrvalue; i < 1 && nextStart != NULL; i++, nextStart = nextStart->next) {
    continue;
  }
  if (thisStop != NULL && nextStart != NULL) {
    thisStop->data.ptrvalue = MemFree (thisStop->data.ptrvalue);
    thisStop->data.ptrvalue = StringSave ((CharPtr) nextStart->data.ptrvalue);
  }
}

static void ProcessOneContigLine (ValNodePtr line, BioseqPtr segseq, Int4 lineNum,
                                  BoolPtr hasgaps, Boolean coordsOnMaster)

{
  Boolean     badNumber;
  CharPtr     field [5];
  Int2        i;
  Boolean     isgap;
  long int    num;
  Int4        size;
  Int4        start;
  Int4        stop;
  Int2        strand = Seq_strand_unknown;
  Int4        tmp;
  ValNodePtr  vnp;

  if (line == NULL || segseq == NULL) return;
  vnp = line->data.ptrvalue;
  if (vnp != NULL) {
    for (i = 0; i < 5; i++) {
      field [i] = NULL;
    }
    start = -1;
    stop = -1;
    size = -1;
    for (i = 0, vnp = line->data.ptrvalue; i < 5 && vnp != NULL; i++, vnp = vnp->next) {
      if (field [i] == NULL && (! HasNoText ((CharPtr) vnp->data.ptrvalue))) {
        field [i] = (CharPtr) vnp->data.ptrvalue;
      }
    }
  }

  if (HasNoText (accnString)) return;

  badNumber = FALSE;
  if (sizeString != NULL && sscanf (sizeString, "%ld", &num) == 1) {
    size = num;
  } else {
    size = -1;
  }
  if (startString != NULL && sscanf (startString, "%ld", &num) == 1) {
    start = num;
  } else {
    start = -1;
    badNumber = TRUE;
  }
  if (stopString != NULL && sscanf (stopString, "%ld", &num) == 1) {
    stop = num;
  } else {
    stop = -1;
    badNumber = TRUE;
  }
  if (start < 1 || stop < 1) {
    badNumber = TRUE;
  }
  isgap = FALSE;
  if (StringICmp (accnString, "gap") == 0) {
    if (size >= 0) {
      isgap = TRUE;
      badNumber = FALSE;
      start = 1;
      stop = size;
    }
  }

  if (badNumber) {
    if (startString == NULL) startString = "";
    if (stopString == NULL) stopString = "";
    if (start < 1 && stop < 1) {
      Message (MSG_POST, "Bad number in line %ld - start '%s', stop '%s'",
               (long) lineNum, startString, stopString);
    } else if (start < 1) {
      Message (MSG_POST, "Bad number in line %ld - start '%s'", (long) lineNum, startString);
    } else if (stop < 1) {
      Message (MSG_POST, "Bad number in line %ld - stop '%s'", (long) lineNum, stopString);
    } else {
      Message (MSG_POST, "Bad number in line %ld", (long) lineNum);
    }
    return;
  }

  if (isgap) {
    start = 0;
    stop = size - 1;
  } else {
    if (coordsOnMaster && start == stop) {
      Message (MSG_POST, "Ignoring accession %s", accnString);
      return;
    }

    start--;
    stop--;

    strand = Seq_strand_plus;
    if (strandString != NULL) {
      if (StringStr (strandString, "minus") ||
          StringChr (strandString, '-') ||
          StringStr (strandString, "complement")) {
        strand = Seq_strand_minus;
      }
    }
    if (start > stop) {
      tmp = start;
      start = stop;
      stop = tmp;
      strand = Seq_strand_minus;
    }
    if (strandString != NULL) {
      if (StringStr (strandString, "plus") || StringChr (strandString, '+')) {
        strand = Seq_strand_plus;
      }
    }

    if (coordsOnMaster) {
      stop -= (start + 1);
      start = 0;
    }

    size = ABS (stop - start) + 1;
  }

  AddNucToContig (accnString, start, stop, size, strand, segseq, hasgaps, isgap);
}

static void FreeFeatureTable (ValNodePtr head)

{
  ValNodePtr  vnp;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    vnp->data.ptrvalue = ValNodeFreeData (vnp->data.ptrvalue);
  }
  ValNodeFreeData (head);
}

NLM_EXTERN SeqEntryPtr ReadContigList (FILE *fp, Boolean coordinatesOnMaster)

{
  BioseqPtr    bsp;
  DeltaSeqPtr  dsp;
  Boolean      hasgaps;
  ValNodePtr   head = NULL;
  Char         line [1023];
  Int4         lineNum;
  Int4         pos;
  SeqEntryPtr  sep;
  CharPtr      str;
  ValNodePtr   vnp;

  if (fp == NULL) return NULL;

  pos = ftell (fp);
  str = ReadALine (line, sizeof (line), fp);
  if (str != NULL && StringNICmp (line, ">Contig", 7) == 0) {
    str = ReadALine (line, sizeof (line), fp);
  }
  while (str != NULL) {
    if (! HasNoText (line)) {
      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringNCmp (line, "//", 2) == 0 ||
          StringStr (line, "::=") != NULL) {
        fseek (fp, pos, SEEK_SET);
        break;
      }
      vnp = ParseContigOrFeatureTableString (line, TRUE);
      if (vnp != NULL) {
        ValNodeAddPointer (&head, 0, (Pointer) vnp);
      }
    }
    pos = ftell (fp);
    str = ReadALine (line, sizeof (line), fp);
  }
  if (head == NULL) return NULL;

  bsp = BioseqNew ();
  if (bsp == NULL) {
    FreeFeatureTable (head);
    return NULL;
  }
  bsp->mol = Seq_mol_dna;
  bsp->repr = Seq_repr_seg;
  bsp->seq_ext_type = 1;
  bsp->length = 0;
  bsp->id = MakeUniqueSeqID ("contig_");

  if (coordinatesOnMaster) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      if (vnp->data.ptrvalue != NULL) {
        AdjustContigValues (vnp);
      }
    }
  }

  lineNum = 0;
  hasgaps = FALSE;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    lineNum++;
    if (vnp->data.ptrvalue != NULL) {
      ProcessOneContigLine (vnp, bsp, lineNum, &hasgaps, coordinatesOnMaster);
    }
  }

  FreeFeatureTable (head);

  if (bsp->seq_ext == NULL) return NULL;

  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;
  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

  if (hasgaps) {
    dsp = GappedSeqLocsToDeltaSeqs (bsp->seq_ext);
    if (dsp != NULL) {
      bsp->seq_ext = SeqLocSetFree ((ValNodePtr) bsp->seq_ext);
      bsp->repr = Seq_repr_delta;
      bsp->seq_ext_type = 4;
      bsp->seq_ext = (Pointer) dsp;
    }
  }

  return sep;
}

/* ReadAsnFastaOrFlatFile section */

/* GetSeqId skips past LOCUS or ID, or starts past >, skips any white space, then
takes the next token as the seqID.  The return value points to the remaining
copied text, which for feature tables may contain a desired Seq-annot name. */

static CharPtr GetSeqId (CharPtr seqid, CharPtr str, size_t max, Boolean skiptag)

{
  Char     ch;
  CharPtr  ptr;

  if (seqid != NULL) {
    *seqid = '\0';
  }
  if (str == NULL || seqid == NULL) return FALSE;
  if (skiptag) {
    ch = *str;
    while (ch != '\0' && (! IS_WHITESP (ch))) {
      str++;
      ch = *str;
    }
  }
  ch = *str;
  while (IS_WHITESP (ch)) {
    str++;
    ch = *str;
  }
  StringNCpy_0 (seqid, str, max);
  str = seqid;
  while (ch != '\0' && (! IS_WHITESP (ch))) {
    str++;
    ch = *str;
  }
  if (ch == '\0') return NULL;
  *str = '\0';
  str++;
  ch = *str;
  while (ch != '\0' && (IS_WHITESP (ch))) {
    str++;
    ch = *str;
  }
  ptr = str;
  while (ch != '\0' && (! IS_WHITESP (ch))) {
    ptr++;
    ch = *ptr;
  }
  *ptr = '\0';
  return str;
}

/* PreCheckSeqForProteinType saves the current file position, then reads lines of
sequence, checking each character for letters that appear only in proteins.  It then
restores the file position, and returns true if it thinks it found a protein. */

static Boolean PreCheckSeqForProteinType (FILE *fp)

{
  Char     ch;
  Boolean  isProt = FALSE;
  Char     line [1023];
  CharPtr  p;
  Int4     pos;
  CharPtr  str;

  if (fp == NULL) return FALSE;

  pos = ftell (fp);
  str = ReadALine (line, sizeof (line), fp);
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringNCmp (line, "//", 2) == 0 ||
          StringStr (line, "::=") != NULL) {
        fseek (fp, pos, SEEK_SET);
        return isProt;
      }

      p = line;
      ch = *p;
      while (ch != '\0') {
        if (! (IS_ALPHA (ch))) {
          p++;
        } else {
          ch = TO_UPPER (ch);
          if (StringChr ("EFILPQZ", ch) != NULL) {
            isProt = TRUE;
          }
          p++;
        }
        ch = *p;
      }

    }

    str = ReadALine (line, sizeof (line), fp);
  }

  fseek (fp, pos, SEEK_SET);
  return isProt;
}

/* ReadFlatFileDNA reads lines of sequence into a byte store.  Unless it is forced to be
treated as a nucleotide or a protein, it first calls PreCheckSeqForProteinType to look at
the sequence in advance, checking for protein-specific letters. */

static ByteStorePtr ReadFlatFileDNA (FILE *fp, BoolPtr protPtr, Boolean forceNuc, Boolean forceProt)

{
  Char           ch;
  ByteStorePtr   bs = NULL;
  Boolean        isProt = FALSE;
  Char           line [1023];
  CharPtr        p;
  CharPtr        q;
  Int4           pos;
  CharPtr        str;

  if (fp == NULL) return NULL;
  bs = BSNew (1000);
  if (bs == NULL) return NULL;

  if (forceNuc) {
    isProt = FALSE;
  } else if (forceProt) {
    isProt = TRUE;
  } else {
    isProt = PreCheckSeqForProteinType (fp);
  }
  if (protPtr != NULL) {
    *protPtr = isProt;
  }

  pos = ftell (fp);
  str = ReadALine (line, sizeof (line), fp);
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        fseek (fp, pos, SEEK_SET);
        return bs;
      } else if (StringNCmp (line, "//", 2) == 0) {
        return bs;
      }

      p = line;
      q = line;
      ch = *p;
      while (ch != '\0') {
        if (! (IS_ALPHA (ch))) {
          if (isProt && (ch == '*' || ch == '-')) {
            *q = ch;
            q++;
          }
          p++;
        } else {
          ch = TO_UPPER (ch);
          if (! isProt) {
            if (ch == 'U') {
              ch = 'T';
            }
          }
          *q = ch;
          p++;
          q++;
        }
        ch = *p;
      }
      *q = '\0';
      BSWrite (bs, (VoidPtr) line, (Int4) StringLen (line));

    }

    pos = ftell (fp);
    str = ReadALine (line, sizeof (line), fp);
  }

  return bs;
}

/* ReadFeatureTable reads lines of feature intervals and qualifiers into a Seq-annot. */

#define NUM_FTABLE_COLUMNS  6

#define START_TAG           0
#define STOP_TAG            1
#define FEAT_TYPE_TAG       2
#define QUAL_TYPE_TAG       3
#define QUAL_VAL_TAG        4
#define STRAND_TAG          5

#define startStr   field [START_TAG]
#define stopStr    field [STOP_TAG]
#define featType   field [FEAT_TYPE_TAG]
#define qualType   field [QUAL_TYPE_TAG]
#define qualVal    field [QUAL_VAL_TAG]
#define strandStr  field [STRAND_TAG]

static Boolean ParseFeatTableLine (CharPtr line, Int4Ptr startP, Int4Ptr stopP,
                                   CharPtr PNTR featP, CharPtr PNTR qualP, CharPtr PNTR valP)

{
  Boolean     badNumber;
  CharPtr     field [NUM_FTABLE_COLUMNS];
  Int2        i;
  ValNodePtr  parsed;
  Int4        start;
  Int4        stop;
  Int4        tmp;
  long int    val;
  ValNodePtr  vnp;

  if (line == NULL || HasNoText (line)) return FALSE;
  parsed = ParseContigOrFeatureTableString (line, TRUE);
  if (parsed == NULL) return FALSE;

  for (i = 0; i < NUM_FTABLE_COLUMNS; i++) {
    field [i] = NULL;
  }
  start = -1;
  stop = -1;
  vnp = parsed;
  for (i = 0; i < NUM_FTABLE_COLUMNS && vnp != NULL; i++) {
    if (field [i] == NULL) {
      if (! HasNoText ((CharPtr) vnp->data.ptrvalue)) {
        field [i] = (CharPtr) vnp->data.ptrvalue;
      }
    }
    vnp = vnp->next;
  }

  badNumber = FALSE;
  if (startStr != NULL && sscanf (startStr, "%ld", &val) == 1) {
    start = val;
  } else {
    start = -1;
    badNumber = TRUE;
  }
  if (stopStr != NULL && sscanf (stopStr, "%ld", &val) == 1) {
    stop = val;
  } else {
    stop = -1;
    badNumber = TRUE;
  }

  if (badNumber) {
    start = -1;
    stop = -1;
  } else {
    start--;
    stop--;
    if (strandStr != NULL) {
      if (StringStr (strandStr, "minus") ||
          StringChr (strandStr, '-') ||
          StringStr (strandStr, "complement")) {
        if (start < stop) {
          tmp = start;
          start = stop;
          stop = tmp;
        }
      }
    }
  }

  *startP = start;
  *stopP = stop;
  *featP = featType;
  *qualP = qualType;
  *valP = qualVal;

  ValNodeFree (parsed);
  return TRUE;
}

static CharPtr aaList [] = {
  "-", "Gap", "Gap",        /* cannot be recognized because we split tRNA-xxx */
  "A", "Ala", "Alanine",
  "B", "Asx", "Asp or Asn",
  "C", "Cys", "Cysteine",
  "D", "Asp", "Aspartic Acid",
  "E", "Glu", "Glutamic Acid",
  "F", "Phe", "Phenylalanine",
  "G", "Gly", "Glycine",
  "H", "His", "Histidine",
  "I", "Ile", "Isoleucine",
  "K", "Lys", "Lysine",
  "L", "Leu", "Leucine",
  "M", "Met", "Methionine",
  "N", "Asn", "Asparagine",
  "P", "Pro", "Proline",
  "Q", "Gln", "Glutamine",
  "R", "Arg", "Arginine",
  "S", "Ser", "Serine",
  "T", "Thr", "Threonine",
  "V", "Val", "Valine",
  "W", "Trp", "Tryptophan",
  "X", "Xxx", "Undetermined or atypical",
  "Y", "Tyr", "Tyrosine",
  "Z", "Glx", "Glu or Gln",
  "U", "Sec", "Selenocysteine",
  "*", "Ter", "Termination",
  NULL, NULL, NULL
};

static Uint1 FindTrnaAA (CharPtr str)

{
  Uint1    aa;
  Int2     i;
  Int2     j;
  CharPtr  ptr;
  Char     tmp [128];

  if (HasNoText (str)) return 0;
  StringNCpy_0 (tmp, str, sizeof (tmp));
  SqnTrimSpacesAroundString (tmp);
  for (i = 0; aaList [i] != NULL; i += 3) {
    for (j = 0; j < 3; j++) {
      if (StringICmp (aaList [i + j], tmp) == 0) {
        ptr = aaList [i];
        aa = (Uint1) ptr [0];
        return aa;
      }
    }
  }
  return 0;
}

static ValNodePtr TokenizeTRnaString (CharPtr strx)

{
  Char        ch;
  ValNodePtr  head;
  Int2        i, j, k;
  CharPtr     ptr;
  Char        str [256];
  CharPtr     strs;
  Char        tmp [128];

  if (HasNoText (strx)) return NULL;
  strs = StringSave (strx);
  head = NULL;
  ptr = strs;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '*') {  /* keep possible terminator tRNA symbol */
    } else if (! IS_ALPHA (ch)) {
     *ptr = ' ';
    }
    ptr++;
    ch = *ptr;
  }
  i = 0;
  while (StringLen (strs + i) > 0) {
    StringNCpy_0 (str, strs + i, sizeof (str));
    k = 0;
    ch = str [k];
    while (ch == ' ') {
      k++;
      ch = str [k];
    }
    j = 0;
    while (ch != '\0' && ch != ' ') {
      j++;
      ch = str [j + k];
    }
    if (ch == ' ') {
      str [j + k] = '\0';
      i += j + k + 1;
    } else {
      i += j + k;
    }
    StringNCpy_0 (tmp, str + k, sizeof (tmp));
    SqnTrimSpacesAroundString (tmp);
    if (! HasNoText (tmp)) {
      ValNodeCopyStr (&head, 0, tmp);
    }
  }
  MemFree (strs);
  return head;
}

static CharPtr bondList [] = {
  "", "disulfide", "thiolester", "xlink", "thioether", NULL
};

static CharPtr siteList [] = {
  "", "active", "binding", "cleavage", "inhibit", "modified",
  "glycosylation", "myristoylation", "mutagenized", "metal binding",
  "phosphorylation", "acetylation", "amidation", "methylation",
  "hydroxylation", "sulfatation", "oxidative deamination",
  "pyrrolidone carboxylic acid", "gamma carboxyglutamic acid",
  "blocked", "lipid binding", "np binding", "DNA binding",
  "signal peptide", "transit peptide", "transmembrane region", NULL
};

static void StripHyphens (CharPtr str)

{
  Char     ch;
  CharPtr  ptr;

  if (str == NULL) return;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '-') {
      *ptr = ' ';
    }
    ptr++;
    ch = *ptr;
  }
}

static void AddQualifierToFeature (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  Uint1           aa;
  Uint1           curraa;
  GBQualPtr       gbq;
  GeneRefPtr      grp;
  ValNodePtr      head;
  Boolean         isGeneSyn = FALSE;
  Int2            j;
  GBQualPtr       last;
  ProtRefPtr      prp = NULL;
  Int2            qnum;
  RnaRefPtr       rrp;
  tRNAPtr         trna;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || HasNoText (qual) || HasNoText (val)) return;
  qnum = GBQualNameValid (qual);
  if (qnum <= -1) {
    if (StringNCmp (qual, "gene_syn", 8) == 0) {
      qnum = GBQUAL_gene;
      isGeneSyn = TRUE;
    }
  }
  if (qnum <= -1) {
    if (StringCmp (qual, "region_name") == 0 && sfp->data.choice == SEQFEAT_REGION) {
      sfp->data.value.ptrvalue = MemFree (sfp->data.value.ptrvalue);
      sfp->data.value.ptrvalue = StringSave (val);
    } else if (StringCmp (qual, "bond_type") == 0 && sfp->data.choice == SEQFEAT_BOND) {
      StripHyphens (val);
      sfp->data.value.intvalue = 255;
      for (j = 0; bondList [j] != NULL; j++) {
        if (StringNICmp (val, bondList [j], StringLen (bondList [j])) == 0) {
          sfp->data.value.intvalue = j;
        }
      }
    } else if (StringCmp (qual, "site_type") == 0 && sfp->data.choice == SEQFEAT_SITE) {
      StripHyphens (val);
      sfp->data.value.intvalue = 255;
      for (j = 0; siteList [j] != NULL; j++) {
        if (StringNICmp (val, siteList [j], StringLen (siteList [j])) == 0) {
          sfp->data.value.intvalue = j;
        }
      }
    }
    return;
  }
  if (qnum == GBQUAL_note) {
    sfp->comment = StringSave (val);
    return;
  } else if (qnum == GBQUAL_gene && sfp->data.choice != SEQFEAT_GENE) {
    xref = sfp->xref;
    while (xref != NULL && xref->data.choice != SEQFEAT_GENE) {
      xref = xref->next;
    }
    if (xref == NULL) {
      grp = GeneRefNew ();
      xref = SeqFeatXrefNew ();
      if (xref != NULL) {
        xref->data.choice = SEQFEAT_GENE;
        xref->data.value.ptrvalue = (Pointer) grp;
        xref->next = sfp->xref;
        sfp->xref = xref;
      }
    }
    if (xref != NULL) {
      grp = (GeneRefPtr) xref->data.value.ptrvalue;
      if (grp == NULL) return;
      if (isGeneSyn) {
        ValNodeCopyStr (&(grp->syn), 0, val);
      } else if (grp->locus == NULL) {
        grp->locus = StringSave (val);
      } else {
        ValNodeCopyStr (&(grp->syn), 0, val);
      }
    }
    return;
  }

  if (sfp->data.choice == SEQFEAT_GENE) {
    if (qnum == GBQUAL_gene) {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        if (isGeneSyn) {
          ValNodeCopyStr (&(grp->syn), 0, val);
        } else if (grp->locus == NULL) {
          grp->locus = StringSave (val);
        } else {
          ValNodeCopyStr (&(grp->syn), 0, val);
        }
      }
    }
    return;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    if (qnum == GBQUAL_function || qnum == GBQUAL_EC_number || qnum == GBQUAL_product) {
      xref = sfp->xref;
      while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
        xref = xref->next;
      }
      if (xref == NULL) {
        prp = ProtRefNew ();
        xref = SeqFeatXrefNew ();
        if (xref != NULL) {
          xref->data.choice = SEQFEAT_PROT;
          xref->data.value.ptrvalue = (Pointer) prp;
          xref->next = sfp->xref;
          sfp->xref = xref;
        }
      }
      if (xref != NULL) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
      }
      if (prp == NULL) return;
      if (qnum == GBQUAL_function) {
        ValNodeCopyStr (&(prp->activity), 0, val);
      } else if (qnum == GBQUAL_EC_number) {
        ValNodeCopyStr (&(prp->ec), 0, val);
      } else if (qnum == GBQUAL_product) {
        ValNodeCopyStr (&(prp->name), 0, val);
      }
      return;
    }
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    if (qnum == GBQUAL_product) {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp == NULL) return;
      if (rrp->type == 3) {
        rrp->ext.choice = 2;
        trna = (tRNAPtr) MemNew (sizeof (tRNA));
        rrp->ext.value.ptrvalue = (Pointer) trna;
        if (trna != NULL) {
          trna->aatype = 2;
          head = TokenizeTRnaString (val);
          aa = 0;
          for (vnp = head; (aa == 0 || aa == 'A') && vnp != NULL; vnp = vnp->next) {
            curraa = FindTrnaAA (vnp->data.ptrvalue);
            if (curraa != 0) {
              aa = curraa;
            }
          }
          if (aa == 0) {
            aa = 'X';
          }
          trna->aa = aa;
          for (j = 0; j < 6; j++) {
            trna->codon [j] = 255;
          }
        }
      } else {
        if (rrp->ext.choice == 1) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        }
        rrp->ext.choice = 1;
        rrp->ext.value.ptrvalue = StringSave (val);
      }
      return;
    }
  }

  gbq = GBQualNew ();
  if (gbq == NULL) return;
  gbq->qual = StringSave (qual);
  gbq->val = StringSave (val);
  if (sfp->qual == NULL) {
    sfp->qual = gbq;
  } else {
    last = sfp->qual;
    while (last->next != NULL) {
      last = last->next;
    }
    last->next = gbq;
  }
}

static SeqLocPtr AddIntervalToLocation (SeqLocPtr loc, SeqIdPtr sip, Int4 start, Int4 stop)

{
  Int4       flip;
  SeqIntPtr  sintp;
  SeqLocPtr  slp;
  Int2       strand;
  SeqLocPtr  tmp;

  if (sip == NULL) return NULL;

  sintp = SeqIntNew ();
  strand = Seq_strand_plus;
  if (start > stop) {
    flip = start;
    start = stop;
    stop = flip;
    strand = Seq_strand_minus;
  }
  sintp->from = start;
  sintp->to = stop;
  sintp->strand = strand;
  sintp->id = SeqIdDup (sip);
  slp = ValNodeAddPointer (NULL, SEQLOC_INT, (Pointer) sintp);
  if (loc == NULL) {
    return slp;
  }

  if (loc->choice == SEQLOC_MIX) {
    tmp = (ValNodePtr) (loc->data.ptrvalue);
    while (tmp->next != NULL) {
      tmp = tmp->next;
    }
    tmp->next = slp;
    return loc;
  } else {
    tmp = ValNodeNew (NULL);
    tmp->choice = SEQLOC_MIX;
    tmp->data.ptrvalue = (Pointer) loc;
    loc->next = slp;
    return tmp;
  }
}

static SeqAnnotPtr ReadFeatureTable (FILE *fp, CharPtr seqid, CharPtr annotname)

{
  CdRegionPtr    crp;
  AnnotDescrPtr  desc;
  CharPtr        feat;
  GeneRefPtr     grp;
  ImpFeatPtr     ifp;
  Char           line [1023];
  Int4           pos;
  SeqFeatPtr     prev;
  CharPtr        qual;
  Uint1          rnatype;
  RnaRefPtr      rrp;
  SeqAnnotPtr    sap = NULL;
  SeqFeatPtr     sfp = NULL;
  SeqIdPtr       sip;
  Int4           start;
  Int4           stop;
  CharPtr        str;
  CharPtr        val;

  if (fp == NULL || seqid == NULL) return NULL;
  sip = SeqIdFindBest (MakeSeqID (seqid), 0);
  if (sip == NULL) return NULL;

  pos = ftell (fp);
  str = ReadALine (line, sizeof (line), fp);
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        fseek (fp, pos, SEEK_SET);
        SeqIdFree (sip);
        return sap;
      } else if (StringNCmp (line, "//", 2) == 0) {
        SeqIdFree (sip);
        return sap;
      }

      if (ParseFeatTableLine (line, &start, &stop, &feat, &qual, &val)) {
        if (feat != NULL && start >= 0 && stop >= 0) {

          if (sap == NULL) {
            sap = SeqAnnotNew ();
            if (sap != NULL) {
              sap->type = 1;
              if (! HasNoText (annotname)) {
                desc = ValNodeNew (NULL);
                if (desc != NULL) {
                  desc->choice = Annot_descr_name;
                  desc->data.ptrvalue = StringSave (annotname);
                  sap->desc = desc;
                }
              }
            }
          }

          sfp = SeqFeatNew ();
          if (sfp != NULL && sap != NULL) {
            if (sap->data != NULL) {
              prev = sap->data;
              while (prev->next != NULL) {
                prev = prev->next;
              }
              prev->next = sfp;
            } else {
              sap->data = (Pointer) sfp;
            }

            if (StringCmp (feat, "gene") == 0) {

              sfp->data.choice = SEQFEAT_GENE;
              grp = GeneRefNew ();
              if (grp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) grp;
              }

            } else if (StringCmp (feat, "CDS") == 0) {

              sfp->data.choice = SEQFEAT_CDREGION;
              crp = CreateNewCdRgn (1, FALSE, 0);
              if (crp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) crp;
              }

            } else if (StringStr (feat, "RNA") != NULL) {

              sfp->data.choice = SEQFEAT_RNA;
              rrp = RnaRefNew ();
              if (rrp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) rrp;
                rnatype = 255;
                if (StringCmp (feat, "precursor_RNA") == 0) {
                  rnatype = 1;
                } else if (StringCmp (feat, "mRNA") == 0) {
                  rnatype = 2;
                } else if (StringCmp (feat, "tRNA") == 0) {
                  rnatype = 3;
                } else if (StringCmp (feat, "rRNA") == 0) {
                  rnatype = 4;
                } else if (StringCmp (feat, "misc_RNA") == 0) {
                  rnatype = 255;
                }
                rrp->type = rnatype;
              }

            } else if (StringCmp (feat, "Region") == 0) {

              sfp->data.choice = SEQFEAT_REGION;
              sfp->data.value.ptrvalue = StringSave ("?");

            } else if (StringCmp (feat, "Bond") == 0) {

              sfp->data.choice = SEQFEAT_BOND;
              sfp->data.value.intvalue = 255;

            } else if (StringCmp (feat, "Site") == 0) {

              sfp->data.choice = SEQFEAT_SITE;
              sfp->data.value.intvalue = 255;

            } else {
              sfp->data.choice = SEQFEAT_IMP;
              ifp = ImpFeatNew ();
              if (ifp != NULL) {
                ifp->key = StringSave (feat);
                sfp->data.value.ptrvalue = (Pointer) ifp;
              }
            }

            sfp->location = AddIntervalToLocation (NULL, sip, start, stop);
          }

        } else if (start >= 0 && stop >= 0 && feat == NULL && qual == NULL && val == NULL) {

          sfp->location = AddIntervalToLocation (sfp->location, sip, start, stop);

        } else if (qual != NULL && val != NULL) {

          AddQualifierToFeature (sfp, qual, val);

        }
      }

/* Uncomment next lines if Purify reports memory leaks here. */

      /*
      MemFree (feat);
      MemFree (qual);
      MemFree (val);
      */

    }

    pos = ftell (fp);
    str = ReadALine (line, sizeof (line), fp);
  }

  SeqIdFree (sip);
  return sap;
}

/* ReadRestrictionSiteTable reads lines of restriction enzyme names or cut sites into a Seq-annot. */

static SeqLocPtr AddPointToLocation (SeqLocPtr loc, SeqIdPtr sip, Int4 pt)

{
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;

  if (sip == NULL) return NULL;

  if (loc == NULL) {
    pspp = PackSeqPntNew ();
    pspp->id = SeqIdDup (sip);
    slp = ValNodeNew (NULL);
    slp->choice = SEQLOC_PACKED_PNT;
    slp->data.ptrvalue = (Pointer) pspp;
    loc = slp;
  }

  if (loc != NULL && loc->choice == SEQLOC_PACKED_PNT) {
    pspp = (PackSeqPntPtr) loc->data.ptrvalue;
    if (pspp != NULL) {
      PackSeqPntPut (pspp, pt);
    }
  }

  return loc;
}

static SeqAnnotPtr ReadRestrictionSiteTable (FILE *fp, CharPtr seqid, CharPtr annotname)

{
  DbtagPtr       dbt;
  AnnotDescrPtr  desc;
  Char           line [1023];
  Char           name [64];
  ObjectIdPtr    oip;
  Int4           pos;
  SeqFeatPtr     prev;
  Int4           pt;
  RsiteRefPtr    rrp;
  SeqAnnotPtr    sap = NULL;
  SeqFeatPtr     sfp = NULL;
  SeqIdPtr       sip;
  CharPtr        str;
  long int       val;

  if (fp == NULL || seqid == NULL) return NULL;
  sip = SeqIdFindBest (MakeSeqID (seqid), 0);
  if (sip == NULL) return NULL;
  name [0] = '\0';

  pos = ftell (fp);
  str = ReadALine (line, sizeof (line), fp);
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        fseek (fp, pos, SEEK_SET);
        SeqIdFree (sip);
        return sap;
      } else if (StringNCmp (line, "//", 2) == 0) {
        SeqIdFree (sip);
        return sap;
      }

      if (sscanf (line, "%ld", &val) == 1) {
        pt = (Int4) val;

        if (! HasNoText (name)) {

          if (sap == NULL) {
            sap = SeqAnnotNew ();
            if (sap != NULL) {
              sap->type = 1;
              if (! HasNoText (annotname)) {
                desc = ValNodeNew (NULL);
                if (desc != NULL) {
                  desc->choice = Annot_descr_name;
                  desc->data.ptrvalue = StringSave (annotname);
                  sap->desc = desc;
                }
              }
            }
          }

          if (sfp == NULL) {
            sfp = SeqFeatNew ();
            if (sfp != NULL) {
              sfp->data.choice = SEQFEAT_RSITE;
              dbt = DbtagNew ();
              if (dbt != NULL) {
                dbt->db = StringSave ("REBASE");
                oip = ObjectIdNew ();
                if (oip != NULL) {
                  oip->str = StringSave (name);
                }
                dbt->tag = oip;
              }
              rrp = ValNodeNew (NULL);
              if (rrp != NULL) {
                rrp->choice = 2;
                rrp->data.ptrvalue = dbt;
              }
              sfp->data.value.ptrvalue = (Pointer) rrp;

              if (sap != NULL) {
                if (sap->data != NULL) {
                  prev = sap->data;
                  while (prev->next != NULL) {
                    prev = prev->next;
                  }
                  prev->next = sfp;
                } else {
                  sap->data = (Pointer) sfp;
                }
              }
            }
          }

          sfp->location = AddPointToLocation (sfp->location, sip, pt);

        }

      } else {
        StringNCpy_0 (name, line, sizeof (name));
        sfp = NULL;
      }

    }

    pos = ftell (fp);
    str = ReadALine (line, sizeof (line), fp);
  }

  SeqIdFree (sip);
  return sap;
}

/* ReadMessageStrings allows retired services to announce replacement URLs. */

static void ReadMessageStrings (FILE *fp)

{
  Boolean     done = FALSE;
  ValNodePtr  head = NULL;
  size_t      len;
  Char        line [1023];
  Int4        pos;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;

  if (fp == NULL) return;

  pos = ftell (fp);
  str = ReadALine (line, sizeof (line), fp);
  while (str != NULL && (! done)) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        fseek (fp, pos, SEEK_SET);
        done = TRUE;
      } else if (StringNCmp (line, "//", 2) == 0) {
        done = TRUE;
      }

      if (! done) {
        ValNodeCopyStr (&head, 0, line);
      }
      /* Message (MSG_POST, "%s\n", line); */
    }

    if (! done) {
      pos = ftell (fp);
      str = ReadALine (line, sizeof (line), fp);
    }
  }

  for (vnp = head, len = 0; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (str != NULL) {
      len += StringLen (str) + 1;
    }
  }
  if (len > 0) {
    ptr = MemNew (sizeof (Char) * (len + 2));
    if (ptr != NULL) {
      for (vnp = head, tmp = NULL; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (str != NULL) {
          if (tmp == NULL) {
            tmp = ptr;
          } else {
            tmp = StringMove (tmp, "\n");
          }
          tmp = StringMove (tmp, str);
        }
      }
      Message (MSG_POST, "%s\n", ptr);
      MemFree (ptr);
    }
  }

  ValNodeFreeData (head);
}

/* ReadUidList reads lines of uids (or accessions) into a byte store. */

static ByteStorePtr ReadUidList (FILE *fp, Boolean nucdb, Boolean lastResortSeqIDs)

{
  Boolean       allDigits;
  Boolean       abort = FALSE;
  ByteStorePtr  bs;
  Char          ch;
  Char          line [1023];
  Int4          pos;
  CharPtr       ptr;
  CharPtr       str;
  TextSeqId     tsid;
  Int4          uid;
  long int      val;
  ValNode       vn;

  if (fp == NULL) return NULL;
  bs = BSNew (128);
  if (bs == NULL) return NULL;

  pos = ftell (fp);
  str = ReadALine (line, sizeof (line), fp);
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringNCmp (line, ">", 1) == 0 ||
          StringNCmp (line, "LOCUS ", 6) == 0 ||
          StringNCmp (line, "ID ", 3) == 0 ||
          StringStr (line, "::=") != NULL) {
        fseek (fp, pos, SEEK_SET);
        if (abort) {
          bs = BSFree (bs);
        }
        return bs;
      } else if (StringNCmp (line, "//", 2) == 0) {
        if (abort) {
          bs = BSFree (bs);
        }
        return bs;
      }

      allDigits = TRUE;
      ptr = line;
      ch = *ptr;
      while (ch != '\0' && allDigits) {
        if (! IS_DIGIT (ch)) {
          allDigits = FALSE;
        }
        ptr++;
        ch = *ptr;
      }
      if (allDigits && sscanf (line, "%ld", &val) == 1) {
        uid = (Int4) val;
        BSWrite (bs, &uid, sizeof (Int4));
      } else if (nucdb) {
        tsid.name = NULL;
        tsid.accession = line;
        tsid.release = NULL;
        tsid.version = INT2_MIN;
        vn.choice = (Uint1) SEQID_GENBANK;
        vn.data.ptrvalue = (Pointer) (&tsid);
        uid = GetGIForSeqId (&vn);
        if (uid > 0) {
          BSWrite (bs, &uid, sizeof (Int4));
        } else if (lastResortSeqIDs) {
          abort = TRUE;
        }
      }

    }

    pos = ftell (fp);
    str = ReadALine (line, sizeof (line), fp);
  }

  if (abort) {
    bs = BSFree (bs);
  }
  return bs;
}

/* ReadAsnFastaOrFlatFile reads lines, looking for starts of ASN.1, FASTA, GenBank, EMBL,
or GenPept files.  It then calls the appropriate read function, which is responsible for
reading the sequence (or object) and restoring the file pointer to the beginning of the
next record. */

NLM_EXTERN Pointer ReadAsnFastaOrFlatFile (FILE *fp, Uint2Ptr datatypeptr, Uint2Ptr entityIDptr,
                                           Boolean forceNuc, Boolean forceProt, Boolean parseFastaSeqId)

{
  AsnIoPtr       aip;
  CharPtr        annotname;
  Int4           begin;
  ByteStorePtr   bs = NULL;
  BioseqPtr      bsp = NULL;
  Char           ch;
  Uint1          choice = 0;
  Int2           db = -1;
  Boolean        inLetters;
  Boolean        isProt = FALSE;
  Int4           j;
  Char           line [1023];
  Boolean        mayBeAccessionList = TRUE;
  Boolean        mayBePlainFasta = TRUE;
  Int2           numDigits;
  Int2           numLetters;
  Int4           numLinks;
  ObjectIdPtr    oip;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp = NULL;
  Int4           pos;
  ValNodePtr     pip;
  ProjectPtr     proj = NULL;
  BoolPtr        protPtr;
  Pointer        ptr = NULL;
  SeqAnnotPtr    sap = NULL;
  SeqEntryPtr    sep;
  Char           seqid [128];
  CharPtr        str;
  CharPtr        title = NULL;
  CharPtr        tmp;
  UserFieldPtr   ufp;
  Int4           uid;
  UserObjectPtr  uop;
  ValNodePtr     vnp;

  if (fp == NULL) return NULL;

  if (datatypeptr != NULL) *datatypeptr = 0;
  if (entityIDptr != NULL) *entityIDptr = 0;

  if (forceNuc) {
    isProt = FALSE;
    protPtr = NULL;
  } else if (forceProt) {
    isProt = TRUE;
    protPtr = NULL;
  } else {
    protPtr = &isProt;
  }

  seqid [0] = '\0';
  pos = ftell (fp);
  begin = pos;
  str = ReadALine (line, sizeof (line), fp);
  while (str != NULL) {

    if (! HasNoText (line)) {

      if (StringStr (line, "::=") != NULL) {

        mayBePlainFasta = FALSE;
        mayBeAccessionList = FALSE;
        tmp = line;
        ch = *tmp;
        while (ch != '\0' && (! IS_WHITESP (ch))) {
          tmp++;
          ch = *tmp;
        }
        *tmp = '\0';

        omp = ObjMgrReadLock ();
        omtp = ObjMgrTypeFind (omp, 0, line, NULL);
        ObjMgrUnlock ();

        if (omtp != NULL) {
          fseek (fp, pos, SEEK_SET);
          aip = AsnIoNew (ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
          ptr = (*(omtp->asnread)) (aip, NULL);
          pos = AsnIoTell (aip);
          AsnIoFree (aip, FALSE);
          fseek (fp, pos, SEEK_SET);

          if (ptr == NULL) {
            ErrPostEx (SEV_ERROR, 0, 0, "Couldn't read type [%s]", omtp->asnname);
          } else {
            if (datatypeptr != NULL) {
              *datatypeptr = omtp->datatype;
            }
            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (omtp->datatype, ptr);
            }
          }
        } else {
          ErrPostEx (SEV_ERROR, 0, 0, "Couldn't read type [%s]", line);
        }
        return ptr;

      } else if (line [0] == '>') {

        mayBePlainFasta = FALSE;
        mayBeAccessionList = FALSE;
        db = -1;
        if (StringNCmp (line, ">PubMed", 7) == 0) {
          db = 0;
        } else if (StringNCmp (line, ">Protein", 8) == 0) {
          db = 1;
        } else if (StringNCmp (line, ">Nucleotide", 11) == 0) {
          db = 2;
        } else if (StringNCmp (line, ">Structure", 10) == 0) {
          db = 3;
        } else if (StringNCmp (line, ">Genome", 7) == 0) {
          db = 4;
        }
        if (db != -1) {
          bs = ReadUidList (fp, (Boolean) (db == 2), FALSE);
          if (bs != NULL) {
            proj = ProjectNew ();
            if (proj != NULL) {
              pip = ValNodeNew (NULL);
              if (pip != NULL) {
                switch (db) {
                  case 0 :
                    choice = ProjectItem_pmuid;
                    break;
                  case 1 :
                    choice = ProjectItem_protuid;
                    break;
                  case 2 :
                    choice = ProjectItem_nucuid;
                    break;
                  case 3 :
                    choice = ProjectItem_genomeuid;
                    break;
                  case 4 :
                    choice = ProjectItem_structuid;
                    break;
                  default :
                    choice = 0;
                    break;
                }
                pip->choice = choice;
                proj->data = pip;
                numLinks = BSLen (bs) / sizeof (Int4);
                BSSeek (bs, 0L, 0);
                for (j = 0; j < numLinks; j++) {
                  BSRead (bs, &uid, sizeof (Int4));
                  ValNodeAddInt ((ValNodePtr PNTR) &(pip->data.ptrvalue), choice, uid);
                  /*
                  switch (db) {
                    case 0 :
                      ValNodeAddInt (&(pip->data.ptrvalue), ProjectItem_pmid, uid);
                      break;
                    case 1 :
                    case 2 :
                    case 3 :
                      sip = ValNodeNew (NULL);
                      if (sip != NULL) {
                        sip->choice = SEQID_GI;
                        sip->data.intvalue = uid;
                      }
                      break;
                    case 4 :
                      break;
                    default :
                      break;
                  }
                  */
                }
              }
            }
            bs = BSFree (bs);

            if (datatypeptr != NULL) {
              *datatypeptr = OBJ_PROJECT;
            }
            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (OBJ_PROJECT, (Pointer) proj);
            }

            return (Pointer) proj;
          }

        } else if (StringNICmp (line, ">Feature", 8) == 0) {

          annotname = GetSeqId (seqid, line, sizeof (seqid), TRUE);
          if (! HasNoText (seqid)) {
            sap = ReadFeatureTable (fp, seqid, annotname);
            if (sap != NULL) {
              if (datatypeptr != NULL) {
                *datatypeptr = OBJ_SEQANNOT;
              }
              if (entityIDptr != NULL) {
                *entityIDptr = ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
              }
              return (Pointer) sap;
            }
          }

        } else if (StringNICmp (line, ">Restriction", 12) == 0) {

          annotname = GetSeqId (seqid, line, sizeof (seqid), TRUE);
          if (! HasNoText (seqid)) {
            sap = ReadRestrictionSiteTable (fp, seqid, annotname);
            if (sap != NULL) {
              if (datatypeptr != NULL) {
                *datatypeptr = OBJ_SEQANNOT;
              }
              if (entityIDptr != NULL) {
                *entityIDptr = ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
              }
              return (Pointer) sap;
            }
          }

        } else if (StringNICmp (line, ">Contig", 7) == 0) {

          sep = ReadContigList (fp, FALSE);
          if (sep != NULL && IS_Bioseq (sep)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {

              oip = ObjectIdNew ();
              oip->str = StringSave ("info");
              uop = UserObjectNew ();
              uop->type = oip;
              uop->_class = StringSave ("Genomes");

              oip = ObjectIdNew ();
              oip->id = 0;
              ufp = UserFieldNew ();
              ufp->choice = 2;
              ufp->data.intvalue = 0;
              ufp->label = oip;

              uop->data = ufp;

              vnp = ValNodeNew (NULL);
              vnp->choice = Seq_descr_user;
              vnp->data.ptrvalue = (Pointer) uop;
              vnp->next = bsp->descr;
              bsp->descr = vnp;

              if (datatypeptr != NULL) {
                *datatypeptr = OBJ_BIOSEQ;
              }
              if (entityIDptr != NULL) {
                *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
              }
            }
            return (Pointer) bsp;
          }

        } else if (StringNICmp (line, ">Message", 8) == 0) {

          ReadMessageStrings (fp);

        } else {

          title = NULL;
          tmp = StringChr (line + 1, '[');
          if (tmp != NULL) {
            if (StringStr (tmp, "[prot=") != NULL || StringStr (tmp, "[gene=") != NULL) {
              title = StringSave (tmp);
            }
          }
          if (parseFastaSeqId) {
            GetSeqId (seqid, line + 1, sizeof (seqid), FALSE);
            if (! HasNoText (seqid)) {
              bs = ReadFlatFileDNA (fp, protPtr, forceNuc, forceProt);
            }
          } else {
            bs = ReadFlatFileDNA (fp, protPtr, forceNuc, forceProt);
          }
          if (bs == NULL && title != NULL) {
            title = MemFree (title);
          }
        }

      } else if (StringNCmp (line, "LOCUS ", 6) == 0 || StringNCmp (line, "ID ", 3) == 0) {

        mayBePlainFasta = FALSE;
        mayBeAccessionList = FALSE;
        GetSeqId (seqid, line, sizeof (seqid), TRUE);

      } else if (StringNCmp (line, "ORIGIN ", 7) == 0 || StringNCmp (line, "SQ ", 3) == 0) {

        mayBePlainFasta = FALSE;
        mayBeAccessionList = FALSE;
        if (! HasNoText (seqid)) {
          bs = ReadFlatFileDNA (fp, protPtr, forceNuc, forceProt);
        }

      } else {
        if (mayBePlainFasta) {
          tmp = line;
          ch = *tmp;
          while (ch != '\0') {
            if (IS_WHITESP (ch)) {
            } else if (! (IS_ALPHA (ch) || ch == '*' || ch == '-')) {
              mayBePlainFasta = FALSE;
            } else if (protPtr != NULL) {
              ch = TO_UPPER (ch);
              if (StringChr ("EFILPQZ", ch) != NULL) {
                isProt = TRUE;
              }
            }
            tmp++;
            ch = *tmp;
          }
        }
        if (mayBeAccessionList) {
          inLetters = TRUE;
          numLetters = 0;
          numDigits = 0;
          tmp = line;
          ch = *tmp;
          while (ch != '\0') {
            if (IS_WHITESP (ch)) {
            } else if (IS_ALPHA (ch)) {
              if (! inLetters) {
                mayBeAccessionList = FALSE;
                numLetters++;
              }
            } else if (IS_DIGIT (ch)) {
              inLetters = FALSE;
              numDigits++;
            } else {
              mayBeAccessionList = FALSE;
            }
            tmp++;
            ch = *tmp;
          }
          if (numLetters == 1 && numDigits == 5) {
          } else if (numLetters == 2 && numDigits == 6) {
          } else {
            mayBeAccessionList = FALSE;
          }
        }
      }

      if (bs != NULL) {
        sep = SeqEntryNew ();
        if (sep != NULL) {
          bsp = BioseqNew ();
          if (bsp != NULL) {
            sep->choice = 1;
            sep->data.ptrvalue = (Pointer) bsp;
            SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

            if (isProt) {
              bsp->mol = Seq_mol_aa;
              bsp->seq_data_type = Seq_code_ncbieaa;
            } else {
              bsp->mol = Seq_mol_na;
              bsp->seq_data_type = Seq_code_iupacna;
            }

            bsp->repr = Seq_repr_raw;
            bsp->length = 0;
            if (parseFastaSeqId) {
              bsp->id = MakeSeqID (seqid);
            } else {
              bsp->id = MakeNewProteinSeqId (NULL, NULL);
            }
            SeqMgrAddToBioseqIndex (bsp);

            bsp->seq_data = bs;
            bsp->length = BSLen (bs);

            BioseqPack (bsp);

            if (title != NULL) {
              vnp = CreateNewDescriptor (sep, Seq_descr_title);
              if (vnp != NULL) {
                vnp->data.ptrvalue = (Pointer) title;
                title = NULL;
              }
              bsp->descr = vnp;
            }

            if (datatypeptr != NULL) {
              *datatypeptr = OBJ_BIOSEQ;
            }
            if (entityIDptr != NULL) {
              *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
            }
          }
        }
        return (Pointer) bsp;
      }

    }

    pos = ftell (fp);
    str = ReadALine (line, sizeof (line), fp);
  }

  if (mayBePlainFasta) {

    fseek (fp, begin, SEEK_SET);
    sep = FastaToSeqEntryEx (fp, (Boolean) (! isProt), NULL, FALSE);
    if (sep != NULL && IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != NULL) {
        if (datatypeptr != NULL) {
          *datatypeptr = OBJ_BIOSEQ;
        }
        if (entityIDptr != NULL) {
          *entityIDptr = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
        }
      }
      return (Pointer) bsp;
    }

  }

  if (mayBeAccessionList) {

    fseek (fp, begin, SEEK_SET);
    bs = ReadUidList (fp, TRUE, TRUE);
    if (bs != NULL) {
      numLinks = BSLen (bs) / sizeof (Int4);
      if (numLinks < 1) {
        bs = BSFree (bs);
        return NULL;
      }
      proj = ProjectNew ();
      if (proj != NULL) {
        pip = ValNodeNew (NULL);
        if (pip != NULL) {
          choice = ProjectItem_nucuid;
          pip->choice = choice;
          proj->data = pip;
          BSSeek (bs, 0L, 0);
          for (j = 0; j < numLinks; j++) {
            BSRead (bs, &uid, sizeof (Int4));
            ValNodeAddInt ((ValNodePtr PNTR) &(pip->data.ptrvalue), choice, uid);
          }
        }
      }
      bs = BSFree (bs);

      if (datatypeptr != NULL) {
        *datatypeptr = OBJ_PROJECT;
      }
      if (entityIDptr != NULL) {
        *entityIDptr = ObjMgrRegister (OBJ_PROJECT, (Pointer) proj);
      }

      return (Pointer) proj;
    }

  }

  return NULL;
}

