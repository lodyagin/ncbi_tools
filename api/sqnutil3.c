/*   sqnutil3.c
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
* File Name:  sqnutil3.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   2/7/00
*
* $Revision: 6.231 $
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
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <tofasta.h>
#include <parsegb.h>
#include <utilpars.h>
#include <validatr.h>
#include <explore.h>
#include <salsap.h>
#include <salutil.h>
#include <salpedit.h>
#include <alignmgr2.h>
#include <actutils.h>
#include <utilpub.h>
/* included for discrepancy report */
#include <asn2gnbk.h>
#include <asn2gnbp.h>
#include <valid.h>
#include <findrepl.h>

/* functions for associating CDS and parent mRNA using featureIDs */

NLM_EXTERN void ClearFeatIDs (
  SeqFeatPtr sfp
)

{
  if (sfp == NULL) return;
  SeqFeatIdFree (&sfp->id);
  sfp->id.choice = 0;
}

NLM_EXTERN void ClearFeatIDXrefs (
  SeqFeatPtr sfp
)

{
  SeqFeatXrefPtr  xref, next, PNTR prevlink;

  if (sfp == NULL) return;

  prevlink = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  xref = sfp->xref;
  while (xref != NULL) {
    next = xref->next;

    if (xref->id.choice != 0) {
      SeqFeatIdFree (&xref->id);
      xref->id.choice = 0;
    }
    if (xref->id.choice == 0 && xref->data.choice == 0) {
      *prevlink = xref->next;
      xref->next = NULL;
      MemFree (xref);
    } else {
      prevlink = (SeqFeatXrefPtr PNTR) &(xref->next);
    }

    xref = next;
  }
}

static void SfpClearFeatIDs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL) return;
  ClearFeatIDs (sfp);
  ClearFeatIDXrefs (sfp);
}

NLM_EXTERN void ClearFeatureIDs (
  SeqEntryPtr sep
)

{
  VisitFeaturesInSep (sep, NULL, SfpClearFeatIDs);
}

typedef struct idpair {
  Int4  before;
  Int4  after;
} IdPairData, PNTR IdPairPtr;

typedef struct fiddata {
  Int4       highestID;
  Int4       highestRef;
  Int4       offset;
  Int4       count;
  IdPairPtr  pairs;
} FidData, PNTR FidDataPtr;

static void FindHighestFeatID (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr      fip;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;

  if (sfp->id.choice == 3) {
    oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        if (oip->id >= fip->highestID) {
          fip->highestID = oip->id;
        }
      }
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 3) continue;
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        if (oip->id >= fip->highestRef) {
          fip->highestRef = oip->id;
        }
      }
    }
  }
}

NLM_EXTERN Int4 FindHighestFeatureID (
  SeqEntryPtr sep
)

{
  FidData  fd;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.highestID = 0;
  fd.highestRef = 0;
  VisitFeaturesInSep (sep, (Pointer) &fd, FindHighestFeatID);
  return fd.highestID;
}

static void SfpAssignFeatIDs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr   fip;
  ObjectIdPtr  oip;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;

  if (sfp->id.choice == 3) return;
  oip = ObjectIdNew ();
  if (oip == NULL) return;

  (fip->highestID)++;
  oip->id = fip->highestID;

  sfp->id.value.ptrvalue = (Pointer) oip;
  sfp->id.choice = 3;
}

NLM_EXTERN void AssignFeatureIDs (
  SeqEntryPtr sep
)

{
  FidData  fd;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.highestID = 0;
  fd.highestRef = 0;
  VisitFeaturesInSep (sep, (Pointer) &fd, FindHighestFeatID);
  VisitFeaturesInSep (sep, (Pointer) &fd, SfpAssignFeatIDs);
}

static void SfpOffsetFeatIDs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr   fip;
  ObjectIdPtr  oip;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;

  if (sfp->id.choice == 3) {
    oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        oip->id += fip->offset;
      }
    }
  }
}

NLM_EXTERN void OffsetFeatureIDs (
  SeqEntryPtr sep,
  Int4 offset
)

{
  FidData  fd;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.offset = offset;
  VisitFeaturesInSep (sep, (Pointer) &fd, SfpOffsetFeatIDs);
}

static void SfpOffsetFeatIDXrefs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr      fip;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 3) continue;
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        oip->id += fip->offset;
      }
    }
  }
}

NLM_EXTERN void OffsetFeatureIDXrefs (
  SeqEntryPtr sep,
  Int4 offset
)

{
  FidData  fd;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.offset = offset;
  VisitFeaturesInSep (sep, (Pointer) &fd, SfpOffsetFeatIDXrefs);
}

static void SfpMakePairList (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr   fip;
  Int4         idx;
  IdPairPtr    ipp;
  ObjectIdPtr  oip;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;
  if (fip->pairs == NULL) return;

  if (sfp->id.choice != 3) return;
  oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
  if (oip == NULL) return;

  idx = fip->highestID;
  ipp = &(fip->pairs [idx]);

  (fip->highestID)++;
  ipp->before = oip->id;
  ipp->after = fip->highestID;
}

static int LIBCALLBACK SortPairList (VoidPtr ptr1, VoidPtr ptr2)

{
  IdPairPtr  ipp1 = (IdPairPtr) ptr1;
  IdPairPtr  ipp2 = (IdPairPtr) ptr2;

  if (ipp1 == NULL || ipp2 == NULL) return 0;
  if (ipp1->before > ipp2->before) return 1;
  if (ipp1->before < ipp2->before) return -1;
  return 0;
}

static Int4 LookupNewFeatID (
  FidDataPtr fip,
  Int4 before
)

{
  IdPairPtr  ipp;
  Int4       L;
  Int4       mid;
  Int4       R;

  if (fip == NULL || fip->pairs == NULL || fip->count < 1) return 0;

  L = 0;
  R = fip->count - 1;
  while (L < R) {
    mid = (L + R) / 2;
    ipp = &(fip->pairs [mid]);
    if (ipp->before < before) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (R < fip->count) {
    ipp = &(fip->pairs [R]);
    if (ipp->before == before) return ipp->after;
  }

  return 0;
}

static void SfpReassignPairList (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  FidDataPtr      fip;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  fip = (FidDataPtr) userdata;
  if (fip == NULL) return;
  if (fip->pairs == NULL) return;

  if (sfp->id.choice == 3) {
    oip = (ObjectIdPtr) sfp->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        oip->id = LookupNewFeatID (fip, oip->id);
      }
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice != 3) continue;
    oip = (ObjectIdPtr) xref->id.value.ptrvalue;
    if (oip != NULL) {
      if (oip->str == NULL) {
        oip->id = LookupNewFeatID (fip, oip->id);
      }
    }
  }
}

NLM_EXTERN void ReassignFeatureIDs (
  SeqEntryPtr sep
)

{
  Int4     count;
  FidData  fd;

  count = VisitFeaturesInSep (sep, NULL, NULL);
  if (count < 1) return;

  MemSet ((Pointer) &fd, 0, sizeof (FidData));
  fd.highestID = 0;
  fd.highestRef = 0;
  fd.count = count;
  fd.pairs = (IdPairPtr) MemNew (sizeof (IdPairData) * (count + 1));
  if (fd.pairs == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) &fd, SfpMakePairList);

  HeapSort (fd.pairs, (size_t) count, sizeof (IdPairData), SortPairList);

  VisitFeaturesInSep (sep, (Pointer) &fd, SfpReassignPairList);

  MemFree (fd.pairs);
}

typedef struct vcmdata {
  Boolean     accounted_for;
  SeqFeatPtr  cds;
  SeqFeatPtr  mrna;
  SeqFeatPtr  partner;
} VcmData, PNTR VcmDataPtr;

typedef struct loopdata {
  Int2        count;
  SeqFeatPtr  cds;
  SeqFeatPtr  mrna;
} LoopData, PNTR LoopDataPtr;

static Boolean LIBCALLBACK GetSingleMrnaProc (
  SeqFeatPtr mrna,
  SeqMgrFeatContextPtr context
)

{
  LoopDataPtr  ldp;
  VcmDataPtr   vdp;

  ldp = (LoopDataPtr) context->userdata;

  vdp = (VcmDataPtr) mrna->idx.scratch;
  if (vdp != NULL && vdp->accounted_for) return TRUE;

  (ldp->count)++;
  ldp->mrna = mrna;

  return TRUE;
}

static void BspLinkCDSmRNAbyOverlap (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Int2               count;
  SeqMgrFeatContext  fcontext;
  Boolean            goOn;
  Int4               id;
  LoopData           ld;
  ObjectIdPtr        oip;
  SeqFeatPtr         partner, sfp;
  VcmDataPtr         vdp;
  SeqFeatXrefPtr     xref;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;

  /* add scratch structure to CDS and mRNA features */

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
  while (sfp != NULL) {
    sfp->idx.scratch = (Pointer) MemNew (sizeof (VcmData));
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &fcontext);
  while (sfp != NULL) {
    sfp->idx.scratch = (Pointer) MemNew (sizeof (VcmData));
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_mRNA, &fcontext);
  }

  /* loop through CDS features, finding single unused mRNA partner */

  goOn = TRUE;
  while (goOn) {
    goOn = FALSE;
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
    while (sfp != NULL) {
      vdp = (VcmDataPtr) sfp->idx.scratch;
      if (vdp != NULL && (! vdp->accounted_for)) {
        ld.count = 0;
        ld.cds = sfp;
        ld.mrna = NULL;
        if (sfp->excpt &&
            (StringISearch (sfp->except_text, "ribosomal slippage") != NULL ||
             StringISearch (sfp->except_text, "trans-splicing") != NULL)) {
          count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                   LOCATION_SUBSET, (Pointer) &ld,
                                                   GetSingleMrnaProc);
        } else {
          count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_mRNA, NULL, 0,
                                                   CHECK_INTERVALS, (Pointer) &ld,
                                                   GetSingleMrnaProc);
        }
        if (ld.count == 1 && ld.mrna != NULL) {
          vdp->accounted_for = TRUE;
          vdp->cds = ld.cds;
          vdp->mrna = ld.mrna;
          vdp->partner = ld.mrna;
          vdp = (VcmDataPtr) ld.mrna->idx.scratch;
          if (vdp != NULL) {
            vdp->accounted_for = TRUE;
            vdp->cds = ld.cds;
            vdp->mrna = ld.mrna;
            vdp->partner = ld.cds;
            goOn = TRUE;
          }
        }
      }
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext);
    }
  }

  /* assign xrefs between CDS and mRNA features */

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    vdp = (VcmDataPtr) sfp->idx.scratch;
    if (vdp != NULL && vdp->accounted_for) {
      partner = vdp->partner;
      if (partner != NULL && partner->id.choice == 3) {
        oip = (ObjectIdPtr) partner->id.value.ptrvalue;
        if (oip != NULL && oip->str == NULL) {
          id = oip->id;
          if (id > 0) {
            for (xref = sfp->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
            if (xref != NULL) {
              oip = (ObjectIdPtr) xref->id.value.ptrvalue;
              if (oip != NULL) {
                if (oip->str != NULL) {
                  oip->str = MemFree (oip->str);
                }
                oip->id = id;
              }
            } else {
              xref = SeqFeatXrefNew ();
              if (xref != NULL) {
                oip = ObjectIdNew ();
                if (oip != NULL) {
                  oip->id = id;
                  xref->id.choice = 3;
                  xref->id.value.ptrvalue = (Pointer) oip;
                  xref->next = sfp->xref;
                  sfp->xref = xref;
                }
              }
            }
          }
        }
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  /* free scratch structure in CDS and mRNA features */

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (sfp->idx.scratch != NULL) {
      sfp->idx.scratch = MemFree (sfp->idx.scratch);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
}

NLM_EXTERN void LinkCDSmRNAbyOverlap (
  SeqEntryPtr sep
)

{
  AssignFeatureIDs (sep);
  VisitBioseqsInSep (sep, NULL, BspLinkCDSmRNAbyOverlap);
}

static void BspLinkCDSmRNAbyLabel (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqFeatPtr         cds, mrna;
  SeqMgrFeatContext  ccontext;
  SeqMgrFeatContext  mcontext;
  Int4               id;
  ObjectIdPtr        oip;
  SeqFeatXrefPtr     xref;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;

  /* loop through CDS features, finding mRNA partner by label */

  cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &ccontext);
  while (cds != NULL) {
    if (StringDoesHaveText (ccontext.label)) {
      mrna = SeqMgrGetFeatureByLabel (bsp, ccontext.label, 0, FEATDEF_mRNA, &mcontext);
      if (mrna != NULL && StringCmp (ccontext.label, mcontext.label) == 0) {
        if (cds->id.choice == 3 && mrna->id.choice == 3) {

          /* assign xrefs between CDS and mRNA features */

          oip = (ObjectIdPtr) mrna->id.value.ptrvalue;
          if (oip != NULL && oip->str == NULL) {
            id = oip->id;
            if (id > 0) {
              for (xref = cds->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
              if (xref != NULL) {
                oip = (ObjectIdPtr) xref->id.value.ptrvalue;
                if (oip != NULL) {
                  if (oip->str != NULL) {
                    oip->str = MemFree (oip->str);
                  }
                  oip->id = id;
                }
              } else {
                xref = SeqFeatXrefNew ();
                if (xref != NULL) {
                  oip = ObjectIdNew ();
                  if (oip != NULL) {
                    oip->id = id;
                    xref->id.choice = 3;
                    xref->id.value.ptrvalue = (Pointer) oip;
                    xref->next = cds->xref;
                    cds->xref = xref;
                  }
                }
              }
            }
          }

          oip = (ObjectIdPtr) cds->id.value.ptrvalue;
          if (oip != NULL && oip->str == NULL) {
            id = oip->id;
            if (id > 0) {
              for (xref = mrna->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
              if (xref != NULL) {
                oip = (ObjectIdPtr) xref->id.value.ptrvalue;
                if (oip != NULL) {
                  if (oip->str != NULL) {
                    oip->str = MemFree (oip->str);
                  }
                  oip->id = id;
                }
              } else {
                xref = SeqFeatXrefNew ();
                if (xref != NULL) {
                  oip = ObjectIdNew ();
                  if (oip != NULL) {
                    oip->id = id;
                    xref->id.choice = 3;
                    xref->id.value.ptrvalue = (Pointer) oip;
                    xref->next = mrna->xref;
                    mrna->xref = xref;
                  }
                }
              }
            }
          }
        }
      }
    }
    cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, 0, &ccontext);
  }
}

NLM_EXTERN void LinkCDSmRNAbyLabel (
  SeqEntryPtr sep
)

{
  AssignFeatureIDs (sep);
  VisitBioseqsInSep (sep, NULL, BspLinkCDSmRNAbyLabel);
}

typedef struct ovpdata {
  SeqFeatPtr  sfp;
  Char        revstr [42];
} OvpData, PNTR OvpDataPtr;

static int LIBCALLBACK SortOvpByString (VoidPtr ptr1, VoidPtr ptr2)

{
  OvpDataPtr  odp1;
  OvpDataPtr  odp2;
  CharPtr     str1;
  CharPtr     str2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  odp1 = *((OvpDataPtr PNTR) ptr1);
  odp2 = *((OvpDataPtr PNTR) ptr2);
  if (odp1 == NULL || odp2 == NULL) return 0;
  str1 = odp1->revstr;
  str2 = odp2->revstr;
  if (str1 == NULL || str2 == NULL) return 0;
  return StringICmp (str1, str2);
}

static void FindProtBsp (BioseqPtr bsp, Pointer userdata)

{
  BioseqPtr PNTR  protP;

  if (bsp == NULL || ! (ISA_aa (bsp->mol))) return;
  protP = (BioseqPtr PNTR) userdata;
  *protP = bsp;
}

static void BspLinkCDSmRNAbyProduct (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioseqSetPtr       bssp;
  Char               buf [42];
  BioseqPtr          cdna, prot;
  SeqFeatPtr         cds, mrna, sfp;
  OvpDataPtr         PNTR cdsarray = NULL, PNTR mrnaarray = NULL;
  ValNodePtr         cdshead = NULL, mrnahead = NULL, vnp;
  int                compare;
  Uint2              entityID;
  SeqMgrFeatContext  fcontext;
  Int2               i, numcds, nummrna, L, R, mid;
  Int4               id;
  OvpDataPtr         odp;
  ObjectIdPtr        oip;
  SeqEntryPtr        sep;
  SeqIdPtr           sip;
  SeqFeatXrefPtr     xref;

  if (bsp == NULL || ISA_aa (bsp->mol)) return;

  numcds = 0;
  nummrna = 0;

  /* count CDS and mRNA features, make revstr from product SeqId */

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    switch (sfp->idx.subtype) {
      case FEATDEF_CDS :
        if (sfp->product != NULL) {
          numcds++;
          sip = SeqLocId (sfp->product);
          if (sip == NULL) break;
          MakeReversedSeqIdString (sip, buf, sizeof (buf) - 1);
          if (StringHasNoText (buf)) break;
          odp = (OvpDataPtr) MemNew (sizeof (OvpData));
          if (odp == NULL) break;
          odp->sfp = sfp;
          StringCpy (odp->revstr, buf);
          vnp = ValNodeAddPointer (NULL, 0, (Pointer) odp);
          if (vnp == NULL) break;
          vnp->next = cdshead;
          cdshead = vnp;
        }
        break;
      case FEATDEF_mRNA :
        if (sfp->product != NULL) {
          nummrna++;
          sip = SeqLocId (sfp->product);
          if (sip == NULL) break;
          MakeReversedSeqIdString (sip, buf, sizeof (buf) - 1);
          if (StringHasNoText (buf)) break;
          odp = (OvpDataPtr) MemNew (sizeof (OvpData));
          if (odp == NULL) break;
          odp->sfp = sfp;
          StringCpy (odp->revstr, buf);
          vnp = ValNodeAddPointer (NULL, 0, (Pointer) odp);
          if (vnp == NULL) break;
          vnp->next = mrnahead;
          mrnahead = vnp;
        }
        break;
      default :
        break;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }

  if (numcds > 0 && nummrna > 0) {
    cdsarray = (OvpDataPtr PNTR) MemNew (sizeof (OvpDataPtr) * (numcds + 1));
    mrnaarray = (OvpDataPtr PNTR) MemNew (sizeof (OvpDataPtr) * (nummrna + 1));

    /* populate and sort arrays to search for feature by product SeqId */

    if (cdsarray != NULL && mrnaarray != NULL) {
      for (vnp = cdshead, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        cdsarray [i] = (OvpDataPtr) vnp->data.ptrvalue;
      }
      for (vnp = mrnahead, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        mrnaarray [i] = (OvpDataPtr) vnp->data.ptrvalue;
      }

      HeapSort (cdsarray, (size_t) numcds, sizeof (OvpDataPtr), SortOvpByString);
      HeapSort (mrnaarray, (size_t) nummrna, sizeof (OvpDataPtr), SortOvpByString);

      for (i = 0; i < nummrna; i++) {
        odp = (OvpDataPtr) mrnaarray [i];
        if (odp == NULL) continue;
        mrna = odp->sfp;
        if (mrna == NULL || mrna->product == NULL) continue;
        sip = SeqLocId (mrna->product);
        if (sip == NULL) continue;

        cdna = BioseqLockById (sip);
        if (cdna == NULL) continue;
        entityID = ObjMgrGetEntityIDForPointer (cdna);
        if (entityID < 1) continue;
        if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
          sep = GetTopSeqEntryForEntityID (entityID);
          if (sep == NULL) continue;
          AssignIDsInEntity (entityID, 0, NULL);
        }
        if (cdna->idx.parenttype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) cdna->idx.parentptr;
          if (bssp == NULL) continue;
          if (bssp->_class == BioseqseqSet_class_nuc_prot) {
            prot = NULL;
            if (VisitBioseqsInSet (bssp, (Pointer) &prot, FindProtBsp) == 2) {
              for (sip = prot->id; sip != NULL; sip = sip->next) {
                MakeReversedSeqIdString (sip, buf, sizeof (buf) - 1);
    
                /* binary search */
    
                L = 0;
                R = numcds - 1;
                while (L < R) {
                  mid = (L + R) / 2;
                  odp = cdsarray [mid];
                  compare = StringCmp (odp->revstr, buf);
                  if (compare < 0) {
                    L = mid + 1;
                  } else {
                    R = mid;
                  }
                }
                odp = cdsarray [R];
                if (odp != NULL && StringCmp (odp->revstr, buf) == 0) {
                  cds = odp->sfp;
                  if (cds == NULL) continue;
    
                  /* make reciprocal feature ID xrefs */
    
                  if (cds->id.choice == 3) {
                    oip = (ObjectIdPtr) cds->id.value.ptrvalue;
                    if (oip != NULL && oip->str == NULL) {
                      id = oip->id;
                      if (id > 0) {
                        for (xref = mrna->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
                        if (xref != NULL) {
                          oip = (ObjectIdPtr) xref->id.value.ptrvalue;
                          if (oip != NULL) {
                            if (oip->str != NULL) {
                              oip->str = MemFree (oip->str);
                            }
                            oip->id = id;
                          }
                        } else {
                          xref = SeqFeatXrefNew ();
                          if (xref != NULL) {
                            oip = ObjectIdNew ();
                            if (oip != NULL) {
                              oip->id = id;
                              xref->id.choice = 3;
                              xref->id.value.ptrvalue = (Pointer) oip;
                              xref->next = mrna->xref;
                              mrna->xref = xref;
                            }
                          }
                        }
                      }
                    }
                  }
    
                  if (mrna->id.choice == 3) {
                    oip = (ObjectIdPtr) mrna->id.value.ptrvalue;
                    if (oip != NULL && oip->str == NULL) {
                      id = oip->id;
                      if (id > 0) {
                        for (xref = cds->xref; xref != NULL && xref->id.choice != 3; xref = xref->next) continue;
                        if (xref != NULL) {
                          oip = (ObjectIdPtr) xref->id.value.ptrvalue;
                          if (oip != NULL) {
                            if (oip->str != NULL) {
                              oip->str = MemFree (oip->str);
                            }
                            oip->id = id;
                          }
                        } else {
                          xref = SeqFeatXrefNew ();
                          if (xref != NULL) {
                            oip = ObjectIdNew ();
                            if (oip != NULL) {
                              oip->id = id;
                              xref->id.choice = 3;
                              xref->id.value.ptrvalue = (Pointer) oip;
                              xref->next = cds->xref;
                              cds->xref = xref;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        BioseqUnlock (cdna);
      }
    }

    /* clean up */

    MemFree (cdsarray);
    MemFree (mrnaarray);
  }

  /* more cleanup */

  ValNodeFreeData (cdshead);
  ValNodeFreeData (mrnahead);
}

NLM_EXTERN void LinkCDSmRNAbyProduct (
  SeqEntryPtr sep
)

{
  AssignFeatureIDs (sep);
  VisitBioseqsInSep (sep, NULL, BspLinkCDSmRNAbyProduct);
}

NLM_EXTERN void StripFeatIDXrefAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr    amp;
  AsnTypePtr      atp, atp_se, atp_sfx, atp_sfxe;
  DataVal         dv;
  Boolean         inxrefs;
  SeqFeatXrefPtr  xref;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_sfx = AsnFind ("Seq-feat.xref");
  atp_sfxe = AsnFind ("Seq-feat.xref.E");
  if (atp_se == NULL || atp_sfx == NULL || atp_sfxe == NULL) return;

  inxrefs = FALSE;
  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_sfxe) {
      xref = SeqFeatXrefAsnRead (aip, atp);
      if (xref->data.choice != 0) {
        if (! inxrefs) {
          inxrefs = TRUE;
          AsnOpenStruct (aop, atp_sfx, (Pointer) NULL);
        }
        SeqFeatXrefAsnWrite (xref, aop, atp);
      }
      SeqFeatXrefFree (xref);
    } else if (atp == atp_sfx) {
      AsnReadVal (aip, atp, &dv);
      /* only send struct as open and close item */
      AsnKillValue (atp, &dv);
    } else {
      if (inxrefs) {
        AsnCloseStruct (aop, atp_sfx, (Pointer) NULL);
        inxrefs = FALSE;
      }
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

NLM_EXTERN void StripSeqDataGapAsnFilter (
  AsnIoPtr aip,
  AsnIoPtr aop
)

{
  AsnModulePtr  amp;
  AsnTypePtr    atp, atp_se, atp_dsl;
  DataVal       dv;
  SeqLitPtr     slp;

  if (aip == NULL || aop == NULL) return;

  amp = AsnAllModPtr ();
  if (amp == NULL) return;
  atp_se = AsnFind ("Seq-entry");
  atp_dsl = AsnFind ("Delta-seq.literal");
  if (atp_se == NULL || atp_dsl == NULL) return;

  atp = atp_se;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_dsl) {
      slp = SeqLitAsnRead (aip, atp);
      if (slp != NULL && slp->seq_data != NULL && slp->seq_data_type == Seq_code_gap) {
        slp->seq_data = SeqDataFree (slp->seq_data, slp->seq_data_type);
      }
      SeqLitAsnWrite (slp, aop, atp);
      SeqLitFree (slp);
    } else {
      AsnReadVal (aip, atp, &dv);
      AsnWrite (aop, atp, &dv);
      AsnKillValue (atp, &dv);
    }
  }
}

/* CautiousSeqEntryCleanup section */

static Boolean EmptyOrNullString (CharPtr str)

{
  Char  ch;

  if (str == NULL) return TRUE;
  ch = *str;
  while (ch != '\0') {
    if (ch > ' ' && ch <= '~') return FALSE;
    str++;
    ch = *str;
  }
  return TRUE;
}

/* RemoveMultipleTitles currently removes FIRST title in chain */

static void RemoveMultipleTitles (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqDescrPtr    descr = NULL;
  SeqDescrPtr    lasttitle = NULL;
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    descr = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    descr = bssp->descr;
  } else return;
  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) continue;
    if (lasttitle != NULL) {
      if (lasttitle->extended != 0) {
        ovp = (ObjValNodePtr) lasttitle;
        ovp->idx.deleteme = TRUE;
      }
      lasttitle = sdp;
    } else {
      lasttitle = sdp;
    }
  }
}

static void MakeBioSourceCopy (SeqEntryPtr sep, Pointer userdata)

{
  BioSourcePtr  biop;
  OrgRefPtr     master;
  OrgRefPtr     orp;
  SeqDescrPtr   sdp;

  master = (OrgRefPtr) userdata;
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
  if (sdp != NULL) return;
  biop = BioSourceNew ();
  if (biop == NULL) return;
  orp = OrgRefNew ();
  if (orp == NULL) return;
  biop->org = orp;
  orp->taxname = StringSave (master->taxname);
  orp->common = StringSave (master->common);
  sdp = CreateNewDescriptor (sep, Seq_descr_source);
  if (sdp == NULL) return;
  sdp->data.ptrvalue = (Pointer) biop;
}

static void ReplicatePopPhyMutSetBioSource (SeqEntryPtr sep)

{
  BioSourcePtr   biop;
  BioseqSetPtr   bssp;
  OrgRefPtr      orp;
  ObjValNodePtr  ovp;
  SeqDescrPtr    sdp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (bssp->_class == 7 ||
      (bssp->_class >= 13 && bssp->_class <= 16)) {
    sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
    if (sdp == NULL) return;
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop == NULL) return;
    orp = biop->org;
    if (orp == NULL) return;
    VisitElementsInSep (sep, (Pointer) orp, MakeBioSourceCopy);
    if (sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
    }
  }
}

static SeqFeatPtr BestCDS (SeqLocPtr loc, ValNodePtr cdslist)

{
  SeqFeatPtr  best_cds = NULL;
  SeqFeatPtr  cds;
  Int4        diff;
  Int4        min = INT4_MAX;
  ValNodePtr  vnp;

  if (loc == NULL || cdslist == NULL) return NULL;
  for (vnp = cdslist; vnp != NULL; vnp = vnp->next) {
    cds = (SeqFeatPtr) vnp->data.ptrvalue;
    diff = SeqLocAinB (loc, cds->location);
    if (diff >= 0) {
      if (diff < min) {
        min = diff;
        best_cds = cds;
      }
    }
  }
  return best_cds;
}

#define num_bond 5
static CharPtr feat_bond [num_bond] = {
  NULL,
  "disulfide bond",
  "thiolester bond",
  "xlink bond",
  "thioether bond"
};

#define num_site 27
static CharPtr feat_site [num_site] = {
  NULL, 
  "active", 
  "binding",
  "cleavage",
  "inhibit",
  "modified",
  "glycosylation",
  "myristoylation",
  "mutagenized",
  "metal-binding",
  "phosphorylation",
  "acetylation",
  "amidation",
  "methylation",
  "hydroxylation",
  "sulfatation",
  "oxidative-deamination",
  "pyrrolidone-carboxylic-acid",
  "gamma-carboxyglutamic-acid",
  "blocked",
  "lipid-binding",
  "np-binding",
  "dna-binding",
  "signal-peptide",
  "transit-peptide",
  "transmembrane-region",
  "nitrosylation"
};

static Int2 FindStr (CharPtr PNTR array, Int2 array_num, CharPtr str)

{
  Int2 i;

  for (i = 0; i < array_num; i++) {
    if (array [i] == NULL) continue;
    if (StringNCmp (str, array [i], StringLen (array [i])) == 0) return i;
  }
  return -1;
}

static SeqLocPtr fake_bond_loc (SeqLocPtr slp)

{
  SeqLocPtr loc, l, lnext, ldata;

  if (slp == NULL) return NULL;
  loc = MemNew (sizeof (SeqLoc));
  MemCopy (loc, slp, sizeof (SeqLoc));
  ldata = (SeqLocPtr) loc->data.ptrvalue;
  if (slp->choice != SEQLOC_MIX) return loc;
  for (l = ldata; l != NULL; l = lnext) {
    lnext = l->next;
    if (l->choice == SEQLOC_NULL) {
      ldata = remove_node (ldata, l);
    }
  }
  return loc;
}

static void ConvertImpFeatToProt (SeqFeatPtr feat, Pointer userdata)

{
  SeqFeatPtr  best_cds = NULL;
  Int2        bond = 0;
  BioseqPtr   bsp;
  ValNodePtr  cdslist;
  Uint1       choice = 0;
  Int4        frame;
  ImpFeatPtr  ifp;
  SeqLocPtr   loc;
  Uint1       processed = 0;
  ProtRefPtr  prp;
  SeqFeatPtr  sfp;
  SeqIdPtr    sip;
  Int2        site = 0;
  SeqLocPtr   slp;
  Uint1       subtype = 0;

  if (feat == NULL || feat->data.choice != SEQFEAT_IMP) return;
  ifp = (ImpFeatPtr) feat->data.value.ptrvalue;
  if (ifp == NULL) return;
  cdslist = (ValNodePtr) userdata;
  if (StringCmp (ifp->key, "mat_peptide") == 0) {
    processed = 2;
    choice = SEQFEAT_PROT;
    subtype = FEATDEF_mat_peptide_aa;
  } else if (StringCmp (ifp->key, "sig_peptide") == 0) {
    processed = 3;
    choice = SEQFEAT_PROT;
    subtype = FEATDEF_sig_peptide_aa;
  } else if (StringCmp (ifp->key, "transit_peptide") == 0) {
    processed = 4;
    choice = SEQFEAT_PROT;
    subtype = FEATDEF_transit_peptide_aa;
  } else if (StringCmp (ifp->key, "misc_feature") == 0 && feat->comment != NULL) {
    site = FindStr (feat_site, num_site, feat->comment);
    if (site != -1) {
      choice = SEQFEAT_SITE;
      subtype = FEATDEF_SITE;
    } else {
      bond = FindStr (feat_bond, num_bond, feat->comment);
      if (bond != -1) {
        choice = SEQFEAT_BOND;
        subtype = FEATDEF_BOND;
      }
    }
  }
  if (choice == 0) return;

  if (processed != 0 || site != 0) {
    best_cds = BestCDS (feat->location, cdslist);
  } else if (bond != 0) {
    loc = fake_bond_loc (feat->location);
    best_cds = BestCDS (loc, cdslist);
    SeqLocFree (loc);
  }
  if (best_cds == NULL) return;
  slp = dnaLoc_to_aaLoc (best_cds, feat->location, TRUE, &frame, FALSE);
  if (slp == NULL) return;
  sip = SeqLocId (best_cds->product);
  if (sip == NULL) return;
  bsp = BioseqLockById (sip);
  if (bsp == NULL) return;
  sfp = CreateNewFeatureOnBioseq (bsp, choice, slp);
  BioseqUnlock (bsp);
  if (sfp == NULL) return;

  sfp->partial = feat->partial;
  sfp->excpt = feat->excpt;
  sfp->exp_ev = feat->exp_ev;
  sfp->pseudo = feat->pseudo;

  sfp->comment = feat->comment;
  feat->comment = NULL;
  sfp->qual = feat->qual;
  feat->qual = NULL;
  sfp->title = feat->title;
  feat->title = NULL;
  sfp->ext = feat->ext;
  feat->ext = NULL;
  sfp->cit = feat->cit;
  feat->cit = NULL;

  sfp->xref = feat->xref;
  feat->xref = NULL;
  sfp->dbxref = feat->dbxref;
  feat->dbxref = NULL;
  sfp->except_text = feat->except_text;
  feat->except_text = NULL;

  if (choice == SEQFEAT_PROT) {
    prp = ProtRefNew ();
    sfp->data.value.ptrvalue = (Pointer) prp;
    if (prp != NULL) {
      prp->processed = processed;
    }
    switch (processed) {
    }
  } else if (choice == SEQFEAT_SITE) {
    sfp->data.value.intvalue = site;
  } else if (choice == SEQFEAT_BOND) {
    sfp->data.value.intvalue = bond;
  }
  sfp->idx.subtype = subtype;

  feat->idx.deleteme = TRUE;
}

static void GetListOfCDSs (SeqFeatPtr sfp, Pointer userdata)

{
  ValNodePtr PNTR  head;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  head = (ValNodePtr PNTR) userdata;
  ValNodeAddPointer (head, 0, sfp->data.value.ptrvalue);
}

static void ChangeImpFeatToProt (SeqEntryPtr sep)

{
  ValNodePtr  cdslist = NULL;

  VisitFeaturesInSep (sep, (Pointer) &cdslist, GetListOfCDSs);
  VisitFeaturesInSep (sep, (Pointer) cdslist, ConvertImpFeatToProt);
  ValNodeFree (cdslist);
}

static void MergeAdjacentAnnotsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   nextsap;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1 && nextsap != NULL && nextsap->type == 1) {
      if (sap->id == NULL && nextsap->id == NULL &&
          sap->name == NULL && nextsap->name == NULL &&
          sap->db == 0 && nextsap->db == 0 &&
          sap->desc == NULL && nextsap->desc == NULL &&
          sap->data != NULL && nextsap->data != NULL) {
        sfp = (SeqFeatPtr) sap->data;
        while (sfp->next != NULL) {
          sfp = sfp->next;
        }
        sfp->next = (SeqFeatPtr) nextsap->data;
        nextsap->data = NULL;
        sap->next = nextsap->next;
        SeqAnnotFree (nextsap);
        nextsap = sap->next;
      }
    }
    sap = nextsap;
  }
}

NLM_EXTERN Boolean PubIsEffectivelyEmpty (PubdescPtr pdp)

{
  ValNodePtr  vnp;

  if (pdp == NULL) return FALSE;
  vnp = pdp->pub;
  if (vnp != NULL && vnp->next == NULL && vnp->choice == PUB_Gen) {
    if (empty_citgen ((CitGenPtr) vnp->data.ptrvalue)) {
      return TRUE;
    }
  }
  return FALSE;
}

static void MarkEmptyDescsForCleanup (SeqDescrPtr sdp, Pointer userdata)

{
  GBBlockPtr     gbp;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;
  CharPtr        str;

  if (sdp == NULL || sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  if (sdp->choice == Seq_descr_title) {
    str = (CharPtr) sdp->data.ptrvalue;
    if (StringHasNoText (str)) {
      ovp->idx.deleteme = TRUE;
    }
  } else if (sdp->choice == Seq_descr_pub) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp == NULL) return;
    if (PubIsEffectivelyEmpty (pdp)) {
      ovp->idx.deleteme = TRUE;
    }
  } else if (sdp->choice == Seq_descr_genbank) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp == NULL) return;
    /* gbp->source = MemFree (gbp->source); */
    /* gbp->origin = MemFree (gbp->origin); */
    gbp->taxonomy = MemFree (gbp->taxonomy);
    if (gbp->extra_accessions == NULL && gbp->source == NULL &&
        gbp->keywords == NULL && gbp->origin == NULL &&
        gbp->date == NULL && gbp->entry_date == NULL &&
        gbp->div == NULL && gbp->taxonomy == NULL) {
      ovp->idx.deleteme = TRUE;
    }
  }
}

static void MarkEmptyFeatsForCleanup (SeqFeatPtr sfp, Pointer userdata)

{
  GeneRefPtr  grp;
  PubdescPtr  pdp;
  ProtRefPtr  prp;
  ValNodePtr  vnp;

  if (sfp == NULL) return;
  if (sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (EmptyOrNullString (grp->locus)) {
      grp->locus = MemFree (grp->locus);
    }
    if (EmptyOrNullString (grp->allele)) {
      grp->allele = MemFree (grp->allele);
    }
    if (EmptyOrNullString (grp->desc)) {
      grp->desc = MemFree (grp->desc);
    }
    if (EmptyOrNullString (grp->maploc)) {
      grp->maploc = MemFree (grp->maploc);
    }
    if (EmptyOrNullString (grp->locus_tag)) {
      grp->locus_tag = MemFree (grp->locus_tag);
    }
    if (EmptyOrNullString (grp->locus) &&
        EmptyOrNullString (grp->allele) &&
        EmptyOrNullString (grp->desc) &&
        EmptyOrNullString (grp->maploc) &&
        EmptyOrNullString (grp->locus_tag) &&
        grp->db == NULL && grp->syn == NULL) {
      sfp->idx.deleteme = TRUE;
    }
  } else if (sfp->data.choice == SEQFEAT_PROT && sfp->data.value.ptrvalue != NULL) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp->processed != 3 && prp->processed != 4) {
      vnp = prp->name;
      if ((vnp == NULL || EmptyOrNullString ((CharPtr) vnp->data.ptrvalue)) &&
          EmptyOrNullString (prp->desc) &&
          prp->ec == NULL && prp->activity == NULL && prp->db == NULL) {
        sfp->idx.deleteme = TRUE;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_PUB && sfp->data.value.ptrvalue != NULL) {
    pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    if (PubIsEffectivelyEmpty (pdp)) {
      sfp->idx.deleteme = TRUE;
    }
  } else if (sfp->data.choice == SEQFEAT_COMMENT && EmptyOrNullString (sfp->comment)) {
    sfp->idx.deleteme = TRUE;
  }
}

static void ConvertPubFeatDescProc (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr      bsp;
  size_t         len;
  ObjValNodePtr  ovp;
  PubdescPtr     pdp;
  SeqDescPtr     sdp;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;
  CharPtr        str;
  ValNode        vn;

  /* look for publication features */
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB) return;
  /* get bioseq by feature location */
  sip = SeqLocId (sfp->location);
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  sip = SeqIdFindBest(bsp->id, 0);
  if (sip == NULL) return;
  vn.choice = SEQLOC_WHOLE;
  vn.extended = 0;
  vn.data.ptrvalue = (Pointer) sip;
  vn.next = NULL;
  /* is feature full length? */
  if (SeqLocCompare (sfp->location, &vn) != SLC_A_EQ_B) return;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;
  sdp = CreateNewDescriptor (sep, Seq_descr_pub);
  if (sdp == NULL) return;
  /* move publication from feature to descriptor */
  sdp->data.ptrvalue = sfp->data.value.ptrvalue;
  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.subtype = Seq_descr_pub;
  }
  sfp->data.value.ptrvalue = NULL;
  /* flag old feature for removal */
  sfp->idx.deleteme = TRUE;
  /* move comment to remark */
  if (sfp->comment == NULL) return;
  pdp = (PubdescPtr) sdp->data.ptrvalue;
  if (pdp == NULL) return;
  if (pdp->comment == NULL) {
    pdp->comment = sfp->comment;
    sfp->comment = NULL;
  } else {
    len = StringLen (pdp->comment) + StringLen (sfp->comment) + 5;
    str = MemNew (sizeof (Char) * len);
    StringCpy (str, pdp->comment);
    StringCat (str, "; ");
    StringCat (str, sfp->comment);
    pdp->comment = MemFree (pdp->comment);
    pdp->comment = str;
  }
}

extern void ConvertSourceFeatDescProc (SeqFeatPtr sfp, Pointer userdata)

{
  BioSourcePtr   biop;
  BioseqPtr      bsp;
  SubSourcePtr   lastssp;
  ObjValNodePtr  ovp;
  SeqDescPtr     sdp;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;
  SubSourcePtr   ssp;
  ValNode        vn;
  ValNodePtr     last_dbxref;

  /* look for biosource features */
  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC) return;
  /* get bioseq by feature location */
  sip = SeqLocId (sfp->location);
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  sip = SeqIdFindBest(bsp->id, 0);
  if (sip == NULL) return;
  vn.choice = SEQLOC_WHOLE;
  vn.extended = 0;
  vn.data.ptrvalue = (Pointer) sip;
  vn.next = NULL;
  /* is feature full length? */
  if (SeqLocCompare (sfp->location, &vn) != SLC_A_EQ_B) return;
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;
  sdp = CreateNewDescriptor (sep, Seq_descr_source);
  if (sdp == NULL) return;
  /* move biosource from feature to descriptor */
  sdp->data.ptrvalue = sfp->data.value.ptrvalue;
  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.subtype = Seq_descr_source;
  }
  sfp->data.value.ptrvalue = NULL;
  /* flag old feature for removal */
  sfp->idx.deleteme = TRUE;
  /* move comment to subsource note */
  if (sfp->comment == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  ssp = SubSourceNew ();
  if (ssp == NULL) return;
  ssp->subtype = SUBSRC_other;
  ssp->name = sfp->comment;
  sfp->comment = NULL;
  /* link in at end, since BasicSeqEntry will have sorted this list */
  if (biop->subtype == NULL) {
    biop->subtype = ssp;
  } else {
    lastssp = biop->subtype;
    while (lastssp->next != NULL) {
      lastssp = lastssp->next;
    }
    lastssp->next = ssp;
  }

  /* move dbxrefs on feature to source */
  if (sfp->dbxref != NULL) {
    if (biop->org == NULL) {
      biop->org = OrgRefNew();
    }
    last_dbxref = biop->org->db;
    while (last_dbxref != NULL && last_dbxref->next != NULL) {
      last_dbxref = last_dbxref->next;
    }
    if (last_dbxref == NULL) {    
      biop->org->db = sfp->dbxref;
    } else {
      last_dbxref->next = sfp->dbxref;
    }
    sfp->dbxref = NULL;
  }
}

static void PromoteOrgRefDescToBioSource (SeqDescrPtr sdp, Pointer userdata)

{
  BioSourcePtr   biop;
  OrgRefPtr      orp;
  ObjValNodePtr  ovp;

  if (sdp->choice != Seq_descr_org) return;
  orp = (OrgRefPtr) sdp->data.ptrvalue;
  if (orp == NULL) return;
  biop = BioSourceNew ();
  if (biop == NULL) return;
  biop->org = orp;
  sdp->choice = Seq_descr_source;
  sdp->data.ptrvalue = (Pointer) biop;
  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    ovp->idx.subtype = Seq_descr_source;
  }
}

static void PromoteOrgRefFeatToBioSource (SeqFeatPtr sfp, Pointer userdata)

{
  BioSourcePtr  biop;
  OrgRefPtr     orp;

  if (sfp->data.choice != SEQFEAT_ORG) return;
  orp = (OrgRefPtr) sfp->data.value.ptrvalue;
  if (orp == NULL) return;
  biop = BioSourceNew ();
  if (biop == NULL) return;
  biop->org = orp;
  sfp->data.choice = SEQFEAT_BIOSRC;
  sfp->data.value.ptrvalue = (Pointer) biop;
  sfp->idx.subtype = FEATDEF_BIOSRC;
}

static void DeleteBadMarkedGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatXrefPtr       nextxref;
  SeqFeatXrefPtr PNTR  prevxref;
  Boolean              unlink;
  SeqFeatXrefPtr       xref;

  if (sfp == NULL) return;
  xref = sfp->xref;
  prevxref = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  while (xref != NULL) {
    nextxref = xref->next;
    unlink = FALSE;
    if (xref->specialCleanupFlag && xref->data.choice == SEQFEAT_GENE) {
      if (SeqMgrGetOverlappingGene (sfp->location, NULL) != NULL) {
        unlink = TRUE;
      }
    }
    xref->specialCleanupFlag = FALSE;
    if (unlink) {
      *(prevxref) = xref->next;
      xref->next = NULL;
      SeqFeatXrefFree (xref);
    } else {
      prevxref = (SeqFeatXrefPtr PNTR) &(xref->next);
    }
    xref = nextxref;
  }
}

static void LookForMarkedGeneXrefs (SeqFeatPtr sfp, Pointer userdata)

{
  BoolPtr         hasMarkedGenes;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || sfp->xref == NULL) return;
  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->specialCleanupFlag) {
      hasMarkedGenes = (BoolPtr) userdata;
      *hasMarkedGenes = TRUE;
      return;
    }
  }
}

NLM_EXTERN void CautiousSeqEntryCleanup (SeqEntryPtr sep, SeqEntryFunc taxfun, SeqEntryFunc taxmerge)

{
  /*
  Boolean      correct = FALSE;
  */
  Uint2        entityID;
  Boolean      hasMarkedGenes;
  ErrSev       lsev;
  ErrSev       msev;
  SeqEntryPtr  oldscope;
  /*
  Boolean      strip = TRUE;
  */
  Boolean      taxserver;

  if (sep == NULL) return;
  msev = ErrSetMessageLevel (SEV_MAX);
  lsev = ErrSetLogLevel (SEV_MAX);
  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  BasicSeqEntryCleanup (sep);

  VisitFeaturesInSep (sep, NULL, PromoteOrgRefFeatToBioSource);
  VisitDescriptorsInSep (sep, NULL, PromoteOrgRefDescToBioSource);

  oldscope = SeqEntrySetScope (sep);
  VisitFeaturesInSep (sep, NULL, ConvertSourceFeatDescProc);
  VisitFeaturesInSep (sep, NULL, ConvertPubFeatDescProc);
  SeqEntrySetScope (oldscope);

  VisitFeaturesInSep (sep, NULL, MarkEmptyFeatsForCleanup);
  VisitDescriptorsInSep (sep, NULL, MarkEmptyDescsForCleanup);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);

  SeqEntryExplore (sep, NULL, MergeAdjacentAnnotsCallback);

  ChangeImpFeatToProt (sep);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);

  VisitBioseqsInSep (sep, NULL, ExtendSingleGeneOnMRNA);

  ReplicatePopPhyMutSetBioSource (sep);
  SeqEntryExplore (sep, NULL, RemoveMultipleTitles);

  /* LoopSeqEntryToAsn3 section here */
  taxserver = (Boolean) (taxfun != NULL || taxmerge != NULL);

  /*
  if (correct) {
    SeqEntryExplore(sep, (Pointer)(&porg), CorrectSourceFeat);
  }
  */







  /* a few more things to do here */

  hasMarkedGenes = FALSE;
  VisitFeaturesInSep (sep, (Pointer) &hasMarkedGenes, LookForMarkedGeneXrefs);
  if (hasMarkedGenes) {
    SeqMgrIndexFeatures (entityID, NULL);
    VisitFeaturesInSep (sep, NULL, DeleteBadMarkedGeneXrefs);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }

  BasicSeqEntryCleanup (sep);

  AssignIDsInEntity (entityID, 0, NULL);

  ErrSetMessageLevel (msev);
  ErrSetLogLevel (lsev);
}

/*
static Int4 LoopSeqEntryToAsn3 (SeqEntryPtr sep, Boolean strip, Boolean correct, SeqEntryFunc taxfun, SeqEntryFunc taxmerge)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   oldscope;
  Int4          rsult;
  Boolean       taxserver;

  rsult = 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        rsult += LoopSeqEntryToAsn3 (sep, strip, correct, taxfun, taxmerge);
      }
      return rsult;
    }
  }
  oldscope = SeqEntrySetScope (sep);
  taxserver = (Boolean) (taxfun != NULL || taxmerge != NULL);
  rsult = SeqEntryToAsn3Ex (sep, strip, correct, taxserver, taxfun, taxmerge);
  SeqEntrySetScope (oldscope);
  return rsult;
}
  LoopSeqEntryToAsn3 (sep, TRUE, FALSE, taxfun, taxmerge);

*/

typedef struct featdefNameStruct {
  Uint1    type;
  CharPtr  name;
} FeatdefNameData, PNTR FeatdefNamePtr;

static FeatdefNameData featdefWithName [] = {
  { FEATDEF_10_signal ,          "-10_signal"         },
  { FEATDEF_35_signal ,          "-35_signal"         },
  { FEATDEF_3clip ,              "3'clip"             },
  { FEATDEF_3UTR ,               "3'UTR"              },
  { FEATDEF_5clip ,              "5'clip"             },
  { FEATDEF_5UTR ,               "5'UTR"              },
  { FEATDEF_attenuator ,         "attenuator"         },
  { FEATDEF_BOND ,               "Bond"               },
  { FEATDEF_CAAT_signal ,        "CAAT_signal"        },
  { FEATDEF_CDS ,                "CDS"                },
  { FEATDEF_PUB ,                "Cit"                },
  { FEATDEF_COMMENT ,            "Comment"            },
  { FEATDEF_conflict ,           "conflict"           },
  { FEATDEF_C_region ,           "C_region"           },
  { FEATDEF_D_loop ,             "D-loop"             },
  { FEATDEF_D_segment ,          "D_segment"          },
  { FEATDEF_enhancer ,           "enhancer"           },
  { FEATDEF_exon ,               "exon"               },
  { FEATDEF_gap ,                "gap"                },
  { FEATDEF_GC_signal ,          "GC_signal"          },
  { FEATDEF_GENE ,               "Gene"               },
  { FEATDEF_HET ,                "Het"                },
  { FEATDEF_iDNA ,               "iDNA"               },
  { FEATDEF_IMP ,                "Import"             },
  { FEATDEF_Imp_CDS ,            "Imp_CDS"            },
  { FEATDEF_intron ,             "intron"             },
  { FEATDEF_J_segment ,          "J_segment"          },
  { FEATDEF_LTR ,                "LTR"                },
  { FEATDEF_mat_peptide_aa ,     "mat_peptide"        },
  { FEATDEF_mat_peptide ,        "mat_peptide_nt"     },
  { FEATDEF_misc_binding ,       "misc_binding"       },
  { FEATDEF_misc_difference ,    "misc_difference"    },
  { FEATDEF_misc_feature ,       "misc_feature"       },
  { FEATDEF_misc_recomb ,        "misc_recomb"        },
  { FEATDEF_otherRNA ,           "misc_RNA"           },
  { FEATDEF_misc_signal ,        "misc_signal"        },
  { FEATDEF_misc_structure ,     "misc_structure"     },
  { FEATDEF_modified_base ,      "modified_base"      },
  { FEATDEF_mRNA ,               "mRNA"               },
  { FEATDEF_NON_STD_RESIDUE ,    "NonStdRes"          },
  { FEATDEF_NUM ,                "Num"                },
  { FEATDEF_N_region ,           "N_region"           },
  { FEATDEF_ncRNA ,              "ncRNA"              },
  { FEATDEF_old_sequence ,       "old_sequence"       },
  { FEATDEF_operon ,             "operon"             },
  { FEATDEF_oriT ,               "oriT"               },
  { FEATDEF_polyA_signal ,       "polyA_signal"       },
  { FEATDEF_polyA_site ,         "polyA_site"         },
  { FEATDEF_preRNA ,             "precursor_RNA"      },
  { FEATDEF_preprotein ,         "preprotein"         },
  { FEATDEF_primer_bind ,        "primer_bind"        },
  { FEATDEF_prim_transcript ,    "prim_transcript"    },
  { FEATDEF_promoter ,           "promoter"           },
  { FEATDEF_PROT ,               "Protein"            },
  { FEATDEF_protein_bind ,       "protein_bind"       },
  { FEATDEF_RBS ,                "RBS"                },
  { FEATDEF_REGION ,             "Region"             },
  { FEATDEF_repeat_region ,      "repeat_region"      },
  { FEATDEF_repeat_unit ,        "repeat_unit"        },
  { FEATDEF_rep_origin ,         "rep_origin"         },
  { FEATDEF_rRNA ,               "rRNA"               },
  { FEATDEF_RSITE ,              "Rsite"              },
  { FEATDEF_satellite ,          "satellite"          },
  { FEATDEF_scRNA ,              "scRNA"              },
  { FEATDEF_PSEC_STR ,           "SecStr"             },
  { FEATDEF_sig_peptide_aa ,     "sig_peptide"        },
  { FEATDEF_sig_peptide ,        "sig_peptide_nt"     },
  { FEATDEF_SITE ,               "Site"               },
  { FEATDEF_site_ref ,           "Site-ref"           },
  { FEATDEF_snoRNA ,             "snoRNA"             },
  { FEATDEF_snRNA ,              "snRNA"              },
  { FEATDEF_source ,             "source"             },
  { FEATDEF_BIOSRC ,             "Src"                },
  { FEATDEF_stem_loop ,          "stem_loop"          },
  { FEATDEF_STS ,                "STS"                },
  { FEATDEF_S_region ,           "S_region"           },
  { FEATDEF_TATA_signal ,        "TATA_signal"        },
  { FEATDEF_terminator ,         "terminator"         },
  { FEATDEF_tmRNA ,              "tmRNA"              },
  { FEATDEF_transit_peptide_aa , "transit_peptide"    },
  { FEATDEF_transit_peptide ,    "transit_peptide_nt" },
  { FEATDEF_tRNA ,               "tRNA"               },
  { FEATDEF_TXINIT ,             "TxInit"             },
  { FEATDEF_unsure ,             "unsure"             },
  { FEATDEF_USER ,               "User"               },
  { FEATDEF_variation ,          "variation"          },
  { FEATDEF_virion ,             "virion"             },
  { FEATDEF_V_region ,           "V_region"           },
  { FEATDEF_V_segment ,          "V_segment"          },
  { FEATDEF_SEQ ,                "Xref"               }
};

NLM_EXTERN Uint1 FindFeatDefTypeFromKey (CharPtr key)

{
  Int2  L, R, mid;

  if (key == NULL || *key == '\0') return FEATDEF_BAD;

  L = 0;
  R = (sizeof (featdefWithName) / sizeof (FeatdefNameData)) - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (featdefWithName [mid].name, key) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (featdefWithName [R].name, key) == 0) {
    return featdefWithName [R].type;
  }

  return FEATDEF_BAD;
}

static CharPtr featurekeys [] = {
  "???" ,
  "Gene" ,
  "Org" ,
  "CDS" ,
  "Protein" ,
  "precursor_RNA" ,
  "mRNA" ,
  "tRNA" ,
  "rRNA" ,
  "snRNA" ,
  "scRNA" ,
  "misc_RNA" ,
  "Cit" ,
  "Xref" ,
  "Import" ,
  "allele" ,
  "attenuator" ,
  "C_region" ,
  "CAAT_signal" ,
  "CDS" ,
  "conflict" ,
  "D-loop" ,
  "D_segment" ,
  "enhancer" ,
  "exon" ,
  "GC_signal" ,
  "iDNA" ,
  "intron" ,
  "J_segment" ,
  "LTR" ,
  "mat_peptide" ,
  "misc_binding" ,
  "misc_difference" ,
  "misc_feature" ,
  "misc_recomb" ,
  "misc_RNA" ,
  "misc_signal" ,
  "misc_structure" ,
  "modified_base" ,
  "mutation" ,
  "N_region" ,
  "old_sequence" ,
  "polyA_signal" ,
  "polyA_site" ,
  "precursor_RNA" ,
  "prim_transcript" ,
  "primer_bind" ,
  "promoter" ,
  "protein_bind" ,
  "RBS" ,
  "repeat_region" ,
  "repeat_unit" ,
  "rep_origin" ,
  "S_region" ,
  "satellite" ,
  "sig_peptide" ,
  "source" ,
  "stem_loop" ,
  "STS" ,
  "TATA_signal" ,
  "terminator" ,
  "transit_peptide" ,
  "unsure" ,
  "V_region" ,
  "V_segment" ,
  "variation" ,
  "virion" ,
  "3'clip" ,
  "3'UTR" ,
  "5'clip" ,
  "5'UTR" ,
  "-10_signal" ,
  "-35_signal" ,
  "Site-ref" ,
  "Region" ,
  "Comment" ,
  "Bond" ,
  "Site" ,
  "Rsite" ,
  "User" ,
  "TxInit" ,
  "Num" ,
  "SecStr" ,
  "NonStdRes" ,
  "Het" ,
  "Src" ,
  "proprotein" ,
  "mat_peptide" ,
  "sig_peptide" ,
  "transit_peptide",
  "snoRNA",
  "gap",
  "operon",
  "oriT",
  "ncRNA",
  "tmRNA"
};

NLM_EXTERN CharPtr FindKeyFromFeatDefType (Uint1 type, Boolean forGBFF)

{
  CharPtr  key;

  if (type < FEATDEF_GENE || type >= FEATDEF_MAX) {
    type = FEATDEF_BAD;
  }
  key = featurekeys [type];

  if (forGBFF) {
    if (type == FEATDEF_GENE) {
      key = "gene";
    } else if (type == FEATDEF_REGION ||
               type == FEATDEF_COMMENT ||
               type == FEATDEF_BOND ||
               type == FEATDEF_SITE) {
      key = "misc_feature";
    }
  }

  return key;
}

/* tRNA codon index to codon string lookup table functions */

typedef struct gcCodonStruct {
  Uint1    index;
  CharPtr  codon;
} GcCodonData, PNTR GcCodonPtr;

static CharPtr    gcCodonStrings = NULL;
static GcCodonPtr codonGcIndex = NULL;

/* mapping from NCBI2na to codon codes */

static Uint1 codon_xref [4] = {
  2,  /* A */
  1,  /* C */
  3,  /* G */
  0   /* T */
};

static int LIBCALLBACK SortCodonByString (
  VoidPtr vp1,
  VoidPtr vp2
)

{
  int         compare;
  GcCodonPtr  gcp1 = vp1;
  GcCodonPtr  gcp2 = vp2;

  if (gcp1 == NULL || gcp2 == NULL) return 0;

  compare = StringICmp (gcp1->codon, gcp2->codon);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  return 0;
}

static void InitGcCodons (void)

{
  Uint1           codon [4], index;
  GcCodonPtr      codonGcIdx;
  CharPtr         gcCodonStr;
  Int2            i, j, k;
  int             idx, offset;
  CharPtr         ptr;
  Uint1           residue;
  SeqMapTablePtr  smtp;

  if (codonGcIndex != NULL && gcCodonStrings != NULL) return;

  gcCodonStr = (CharPtr) MemNew (sizeof (Char) * 256);
  if (gcCodonStr == NULL) return;
  codonGcIdx = (GcCodonPtr) MemNew (sizeof (GcCodonData) * 64);
  if (codonGcIdx == NULL) return;

  smtp = SeqMapTableFind (Seq_code_iupacna, Seq_code_ncbi2na);
  if (smtp == NULL) return;

  for (idx = 0; idx < 64; idx++) {
    index = (Uint1) idx;

    for (i = 0, j = 16; i < 3; i++, j /= 4) {
      residue = (Uint1) ((Int2) index / j);
      index -= (Uint1) (residue * j);
      for (k = 0; k < 4; k++) {
        if (codon_xref [k] == residue) {
          residue = (Uint1) k;
          break;
        }
      }
      residue = SeqMapTableConvert (smtp, residue);
      codon [i] = residue;
    }
    codon [3] = 0;

    offset = 4 * idx;
    ptr = gcCodonStr + offset;
    StringCpy (ptr, (CharPtr) codon);

    codonGcIdx [idx].index = (Uint1) idx;
    codonGcIdx [idx].codon = ptr;
  }

  HeapSort (codonGcIdx, (size_t) 64, sizeof (GcCodonData), SortCodonByString);

  gcCodonStrings = gcCodonStr;
  codonGcIndex = codonGcIdx;
}

NLM_EXTERN Uint1 CodonToGcIndex (CharPtr codon)

{
  Char  ch;
  Int2  i, L, R, mid;
  Char  tmp [4];

  if (codonGcIndex == NULL) {
    InitGcCodons ();
  }
  if (codonGcIndex == NULL) return 255;
  if (StringLen (codon) != 3) return 255;
  StringNCpy_0 (tmp, codon, sizeof (tmp));

  for (i = 0; i < 3; i++) {
    ch = tmp [i];
    ch = TO_UPPER (ch);
    if (ch == 'U') {
       ch = 'T';
    }
    tmp [i] = ch;
  }

  L = 0;
  R = 63;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (codonGcIndex [mid].codon, tmp) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (codonGcIndex [R].codon, tmp) == 0) {
    return codonGcIndex [R].index;
  }

  return 255;
}

NLM_EXTERN CharPtr GcIndextoCodon (Uint1 index)

{
  int      offset;
  CharPtr  ptr;

  if (gcCodonStrings == NULL) {
    InitGcCodons ();
  }
  if (gcCodonStrings == NULL) return NULL;
  if (index > 63) return NULL;

  offset = 4 * index;
  ptr = gcCodonStrings + offset;

  return ptr;
}

static FloatHi GetCddBitScore (SeqFeatPtr sfp)

{
  ObjectIdPtr    oip;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;

  if (sfp == NULL) return 0.0;
  uop = sfp->ext;
  if (uop == NULL) return 0.0;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "cddScoreData") != 0) return 0.0;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip != NULL && StringICmp (oip->str, "bit_score") == 0) {
      if (ufp->choice == 3) {
        return ufp->data.realvalue;
      }
    }
  }
  return 0.0;
}

static Boolean FeatIsCDD (
  SeqFeatPtr sfp,
  FloatHi PNTR scoreP
)

{
  DbtagPtr    dbt;
  ValNodePtr  vnp;

  if (scoreP != NULL) {
    *scoreP = 0.0;
  }
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (StringCmp (dbt->db, "CDD") == 0 || StringCmp (dbt->db, "cdd") == 0) {
        if (scoreP != NULL) {
          *scoreP = GetCddBitScore (sfp);
        }
        return TRUE;
      }
    }
  }

  return FALSE;
}
static void BestCDDperBioseq (BioseqPtr bsp, Pointer userdata)

{
  SeqFeatPtr         best;
  SeqMgrFeatContext  context;
  FloatHi            currscore;
  Int4               right;
  SeqFeatPtr         sfp;
  FloatHi            topscore;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
  while (sfp != NULL) {
    if (context.featdeftype == FEATDEF_REGION && FeatIsCDD (sfp, &currscore)) {
      best = sfp;
      right = context.right;
      topscore = currscore;
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
      while (sfp != NULL && context.featdeftype == FEATDEF_REGION &&
             FeatIsCDD (sfp, &currscore) && context.left < right) {
        right = MAX (context.right, right);
        if (currscore <= topscore) {
          sfp->idx.deleteme = TRUE;
        } else {
          best->idx.deleteme = TRUE;
          best = sfp;
          topscore = currscore;
        }
        sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
      }
    } else {
      sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
    }
  }
}

NLM_EXTERN void LeaveBestCDD (SeqEntryPtr sep)

{
  Uint2  entityID;

  if (sep == NULL) return;
  entityID = ObjMgrGetEntityIDForChoice (sep);
  if (entityID < 1) return;

  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  VisitBioseqsInSep (sep, NULL, BestCDDperBioseq);
  DeleteMarkedObjects (entityID, 0, NULL);

  SeqMgrClearFeatureIndexes (entityID, NULL);
}

static CharPtr CompressNonBases (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_ALPHA (ch)) {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

static void LIBCALLBACK SPStreamToRaw (
  CharPtr sequence,
  Pointer userdata
)

{
  ByteStorePtr  bs;
  Char          ch;
  size_t        len;
  CharPtr       tmp;

  bs = (ByteStorePtr) userdata;
  tmp = sequence;
  ch = *tmp;
  while (ch != '\0') {
    if (ch == '\n' || ch == '\r' || ch == '\t') {
      *tmp = ' ';
    } else {
      *tmp = TO_UPPER (ch);
    }
    tmp++;
    ch = *tmp;
  }
  TrimSpacesAroundString (sequence);
  CompressNonBases (sequence);

  len = StringLen (sequence);
  if (len < 1) return;
  BSWrite (bs, sequence, len * sizeof (Char));
}

NLM_EXTERN void SegOrDeltaBioseqToRaw (BioseqPtr bsp)

{
  ByteStorePtr  bs;

  if (bsp == NULL || (bsp->repr != Seq_repr_seg && bsp->repr != Seq_repr_delta)) return;
  if (! ISA_na (bsp->mol)) return;
  bs = BSNew (bsp->length);
  if (bs == NULL) return;

  SeqPortStream (bsp, STREAM_EXPAND_GAPS | STREAM_CORRECT_INVAL, (Pointer) bs, SPStreamToRaw);

  if (bsp->repr == Seq_repr_seg && bsp->seq_ext_type == 1) {
    bsp->seq_ext = SeqLocSetFree ((ValNodePtr) bsp->seq_ext);
    bsp->seq_ext_type = 0;
  } else if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
    bsp->seq_ext = NULL; /* for now just NULL out */
    bsp->seq_ext_type = 0;
  }
  bsp->seq_data = SeqDataFree (bsp->seq_data, bsp->seq_data_type);
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bsp->repr = Seq_repr_raw;
  bsp->seq_data_type = Seq_code_iupacna;
}

typedef struct segtodelta
{
  ValNodePtr seq_ext;
  Int4       len;
  SeqIdPtr   master_sip;
  BioseqPtr  master_bsp;  
  Int4       num_segs_converted;
} SegToDeltaData, PNTR SegToDeltaPtr;


static ValNodePtr CombineDescriptorLists (ValNodePtr target, ValNodePtr insert)
{
  ValNodePtr combined_list = NULL;
  ValNodePtr vnp, vnp_next;
  ValNodePtr title_descr = NULL, prev_descr = NULL;
  CharPtr    combined_title;
  Int4       combined_title_len;
  
  if (target == NULL)
  {
    combined_list = insert;
  }
  else if (insert == NULL)
  {
    combined_list = target;
  }
  else
  {
    combined_list = target;
      for (vnp = target; vnp->next != NULL; vnp = vnp->next)
      {
        if (vnp->choice == Seq_descr_title)
        {
          title_descr = vnp;
        }
      }
      prev_descr = vnp;
      if (title_descr == NULL)
      {
        prev_descr->next = insert;
      }
      else
      {
        for (vnp = insert; vnp != NULL; vnp = vnp_next)
        {
          vnp_next = vnp->next;
          vnp->next = NULL;
          if (vnp->choice == Seq_descr_title)
          {
            /* combine with previous title */
            combined_title_len = StringLen (title_descr->data.ptrvalue)
                                + StringLen (vnp->data.ptrvalue)
                                + 3;
            combined_title = (CharPtr) MemNew (sizeof (Char) * combined_title_len);
            if (combined_title != NULL)
            {
              StringCpy (combined_title, title_descr->data.ptrvalue);
              StringCat (combined_title, "; ");
              StringCat (combined_title, vnp->data.ptrvalue);
              title_descr->data.ptrvalue = MemFree (title_descr->data.ptrvalue);
              title_descr->data.ptrvalue = combined_title;
            }
            ValNodeFreeData (vnp);
          }
          else
          {
            /* add to master list */
            prev_descr->next = vnp;
            prev_descr = vnp;
          }
        } 
      }
  }
  return combined_list;
}

static void MoveSegmentLocToMaster (SeqLocPtr slp, SegToDeltaPtr sdp)
{
  SeqIntPtr     sintp;
  SeqLocPtr     slp2;
  SeqPntPtr     spp;
  PackSeqPntPtr pspp;
  Int4          i;
  
  if (slp == NULL || sdp == NULL) return;
  
  switch (slp->choice)
  {
    case SEQLOC_WHOLE:
    case SEQLOC_EMPTY:
      slp->data.ptrvalue = SeqIdFree (slp->data.ptrvalue);
      slp->data.ptrvalue = SeqIdDup (sdp->master_sip);
      break;
    case SEQLOC_INT:
      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL)
      {
        sintp->id = SeqIdFree (sintp->id);
        sintp->id = SeqIdDup (sdp->master_sip);
        sintp->from += sdp->len;
        sintp->to += sdp->len;
        /* strand stays the same */
      }
      break;
    case SEQLOC_PACKED_INT:
    case SEQLOC_MIX:    
    case SEQLOC_EQUIV:
            slp2 = (SeqLocPtr)slp->data.ptrvalue;
            while (slp2 != NULL)
            {
                MoveSegmentLocToMaster (slp2, sdp);
                slp2 = slp2->next;
            }
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL)
      {
        spp->id = SeqIdFree (spp->id);
        spp->id = SeqIdDup (sdp->master_sip);
        spp->point += sdp->len;
      }
      break;
    case SEQLOC_PACKED_PNT:
            pspp = (PackSeqPntPtr)slp->data.ptrvalue;
            while (pspp != NULL)
      {
        for (i = 0; i < pspp->used; i++)
        {
          pspp->pnts[i] += sdp->len;
        }
        pspp->id = SeqIdFree (pspp->id);
        pspp->id = SeqIdDup (sdp->master_sip);
        pspp = pspp->next;
      }
      break;
  }
}

static void MoveSegmentFeaturesToMaster (SeqFeatPtr sfp, Pointer userdata)

{
  SegToDeltaPtr   segdeltptr;

  if (sfp == NULL || userdata == NULL) return;
  
  segdeltptr = (SegToDeltaPtr) userdata;

  MoveSegmentLocToMaster (sfp->location, segdeltptr);
}

#if 0
static void AdjustAlignmentOffsetsForDeltaConversion (SeqAlignPtr salp, Int4Ptr offsets, BoolPtr is_gap, Int4 num_sets)
{
  DenseSegPtr dsp;
  Int4        aln_seg_num, j, index;
  
  if (salp == NULL || offsets == NULL) return;

  /* adjust alignment starts to match delta sequence coordinates */
  if (salp->segtype == 2)
  {
    dsp = (DenseSegPtr) (salp->segs);
    aln_seg_num = 0;
    for (j = 0; j < num_sets; j++)
    {
      if (!is_gap [j])
      {
        for (index = 0; index < dsp->numseg; index++)
        {
          if (dsp->starts [dsp->dim * index + aln_seg_num] != -1)
          {
            dsp->starts [dsp->dim * index + aln_seg_num] += offsets [j];  
          }
        }
        aln_seg_num++;
      }
    }
  }      
}
#endif

static SeqAnnotPtr CombineAnnots (SeqAnnotPtr target, SeqAnnotPtr insert, Int4 offset)
{
  SeqAnnotPtr combined_list = NULL;
  SeqAnnotPtr feature_sap = NULL;
  SeqAnnotPtr prev_sap = NULL;
  SeqAnnotPtr sap, next_sap;
  SeqFeatPtr  last_feat, first_feat;
  
  if (target == NULL)
  {
    combined_list = insert;
  }
  else if (insert == NULL)
  {
    combined_list = target;
  }
  else
  {
    combined_list = target;
    for (sap = target; sap != NULL; sap = sap->next)
    {
      if (sap->type == 1 && sap->name == NULL && sap->desc == NULL)
      {
        feature_sap = sap;
      }
      prev_sap = sap;
    }
    for (sap = insert; sap != NULL; sap = next_sap)
    {
      next_sap = sap->next;
      sap->next = NULL;
      if (sap->type == 1 && sap->name == NULL && sap->desc == NULL && feature_sap != NULL)
      {
        first_feat = (SeqFeatPtr) sap->data;
        if (first_feat != NULL)
        {
          for (last_feat = (SeqFeatPtr) feature_sap->data;
               last_feat != NULL && last_feat->next != NULL;
               last_feat = last_feat->next)
          {  
          }
          if (last_feat == NULL)
          {
            feature_sap->data = first_feat;    
          }
          else
          {
            last_feat->next = first_feat;
          }
        }
        sap->data = NULL;
        SeqAnnotFree (sap);
      }
      else
      {
        prev_sap->next = sap;
        prev_sap = sap;
      }
    }
  }
  return combined_list;
}

static Int4 AddGapSeqLit (ValNodePtr PNTR seq_ext)
{
  SeqLitPtr       slip;
  IntFuzzPtr      ifp;
  CharPtr         gap_chars = "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN"
                              "NNNNNNNNNN";
                              
  if (seq_ext == NULL) return 0;
                                
  slip = (SeqLitPtr) MemNew (sizeof (SeqLit));
  if (slip != NULL) {
    slip->length = 100;
    ValNodeAddPointer (seq_ext, (Int2) 2, (Pointer) slip);
    ifp = IntFuzzNew ();
    ifp->choice = 4;
      
    slip->fuzz = ifp;
    slip->seq_data = (SeqDataPtr) BSNew (slip->length);
    slip->seq_data_type = Seq_code_iupacna;
    AddBasesToByteStore ((ByteStorePtr) slip->seq_data, gap_chars);
    return 100;
  }
  return 0;
}

static Boolean LIBCALLBACK 
AddSegmentToDeltaSeq 
(SeqLocPtr slp,
 SeqMgrSegmentContextPtr context)

{
  SegToDeltaPtr   segdeltptr;
  SeqIdPtr        sip;
  BioseqPtr       bsp;
  CharPtr         bases;
  SeqLitPtr       slip;

  SeqLocPtr         loc;

  if (slp == NULL || context == NULL) return FALSE;
  segdeltptr = (SegToDeltaPtr) context->userdata;
  if (segdeltptr == NULL) return FALSE;

  sip = SeqLocId (slp);
  
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) 
  {
    return TRUE;
  }

  bsp = BioseqFind (sip);

  if (bsp == NULL)
  {
    return TRUE;
  }
  
  bases = GetSequenceByBsp (bsp);
  if (bases == NULL) 
  {
    bsp->idx.deleteme = TRUE;
    return TRUE;    
  }
  
  if (segdeltptr->seq_ext != NULL)
  {
    /* insert gap of unknown length between the previous segment
     * and this one.
     */
    segdeltptr->len += AddGapSeqLit (&(segdeltptr->seq_ext));
  }

  /* move descriptors to master_bsp */
  segdeltptr->master_bsp->descr = CombineDescriptorLists (segdeltptr->master_bsp->descr, bsp->descr);
  bsp->descr = NULL;
  
  /* move features to master_bsp */
  VisitFeaturesOnBsp (bsp, segdeltptr, MoveSegmentFeaturesToMaster);
  segdeltptr->master_bsp->annot = CombineAnnots (segdeltptr->master_bsp->annot, bsp->annot, segdeltptr->len);
  bsp->annot = NULL;
  
  slip = (SeqLitPtr) MemNew (sizeof (SeqLit));
  if (slip != NULL) 
  {
    slip->length = StringLen (bases);
    ValNodeAddPointer (&(segdeltptr->seq_ext), (Int2) 2, (Pointer) slip);
    slip->seq_data = (SeqDataPtr) BSNew (slip->length);
    slip->seq_data_type = Seq_code_iupacna;
    AddBasesToByteStore ((ByteStorePtr) slip->seq_data, bases);
    segdeltptr->len += slip->length;
  }

  segdeltptr->num_segs_converted ++;
  return TRUE;
}

static BioseqPtr GetDeltaSeqFromMasterSeg (BioseqPtr bsp)
{
  BioseqPtr      new_bsp;
  SegToDeltaData sdd;
  BioseqSetPtr   segset;
  
  if (bsp == NULL || bsp->repr != Seq_repr_seg 
      || bsp->seq_ext == NULL || bsp->seq_ext_type != 1) 
  {
    return NULL;
  }
  
  if (! ISA_na (bsp->mol)) return NULL;

  /* use SeqMgrExploreSegments to build a list of SeqLitPtr */
  sdd.seq_ext = NULL;
  sdd.len = 0;
  sdd.master_bsp = bsp;
  sdd.master_sip = bsp->id;
  sdd.num_segs_converted = 0;
  
  /* move descriptors and features from segset to master seg */
  if (bsp->idx.parenttype == OBJ_BIOSEQSET)
  {
    segset = (BioseqSetPtr) bsp->idx.parentptr;
    if (segset != NULL)
    {
      bsp->descr = CombineDescriptorLists (bsp->descr, segset->descr);
      segset->descr = NULL;
    }
  }  

  SeqMgrExploreSegments (bsp, (Pointer) &sdd, AddSegmentToDeltaSeq);
  
  new_bsp = BioseqNew ();
  new_bsp->descr = bsp->descr;
  bsp->descr = NULL;
  new_bsp->annot = bsp->annot;
  bsp->annot = NULL;
  new_bsp->seq_data = NULL;
  new_bsp->seq_data_type = 0;
  new_bsp->repr = Seq_repr_delta;
  new_bsp->seq_ext_type = 4;
  new_bsp->seq_ext = sdd.seq_ext;
  new_bsp->length = sdd.len;
  new_bsp->id = SeqIdDup (bsp->id); 
/*  new_bsp->id = MakeUniqueSeqID ("delta_"); */
  new_bsp->mol = bsp->mol;

  BioseqPack (new_bsp);  
  return new_bsp;
}

NLM_EXTERN void ConvertSegSetsToDeltaSequences (SeqEntryPtr sep)
{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sub_sep, prev_sep, next_sep;
  ObjMgrDataPtr omdptop;
  ObjMgrData    omdata;
  Uint2         parenttype;
  Pointer       parentptr;
  SeqEntryPtr   new_sep;
  BioseqPtr     bsp, new_bsp = NULL;
  BioseqSetPtr  parent_set;
  
  if (sep == NULL || !IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class == 2)
  {
    SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
    GetSeqEntryParent (sep, &parentptr, &parenttype);
  
    parent_set = (BioseqSetPtr)(bssp->idx.parentptr);
    prev_sep = NULL;
    for (sub_sep = bssp->seq_set; sub_sep != NULL && !IS_Bioseq (sub_sep); sub_sep = sub_sep->next)
    {
      prev_sep = sub_sep;
    }
    if (sub_sep != NULL)
    {
      bsp = sub_sep->data.ptrvalue;
      new_bsp = GetDeltaSeqFromMasterSeg (sub_sep->data.ptrvalue);
      new_sep = SeqEntryNew();
      new_sep->choice = 1;
      new_sep->data.ptrvalue = new_bsp;
            
      /* add new seq entry to parent set */
      AddSeqEntryToSeqEntry (parent_set->seqentry, new_sep, TRUE);

      /* remove segset */      
      bssp->idx.deleteme = TRUE;
    }
    SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
    DeleteMarkedObjects (0, OBJ_BIOSEQSET, parent_set);
    SeqMgrReplaceInBioseqIndex (new_bsp); 
  }
  else
  {
    for (sub_sep = bssp->seq_set; sub_sep != NULL; sub_sep = next_sep)
    {
      next_sep = sub_sep->next;
      ConvertSegSetsToDeltaSequences (sub_sep);
    }
  }
}

static PubMedFetchFunc pmf_pubfetch = NULL;

NLM_EXTERN void LIBCALL PubMedSetFetchFunc (PubMedFetchFunc func)

{
  pmf_pubfetch = func;
}

NLM_EXTERN PubmedEntryPtr LIBCALL GetPubMedForUid (Int4 uid)

{
  PubMedFetchFunc  func;

  if (uid < 1) return NULL;
  func = pmf_pubfetch;
  if (func == NULL) return NULL;
  return func (uid);
}

static Boolean IsTerminator (int c)
{
  if (c == '\n' || c == '\r') {
    return TRUE;
  } else {
    return FALSE;
  }
}

typedef struct bufferedread {
  CharPtr data;
  Int4    len;
  Int4    offset;
} BufferedReadData, PNTR BufferedReadPtr;

static BufferedReadPtr BufferedReadFree (BufferedReadPtr brp)
{
  if (brp == NULL) return NULL;
  if (brp->data != NULL) {
    MemFree (brp->data);
    brp->data = NULL;
  }
  brp->offset = 0;
  brp->len = 0;
  return NULL;
}

extern void FreeBufferedReadList (ValNodePtr vnp)
{
  if (vnp == NULL) return;
  FreeBufferedReadList (vnp->next);
  vnp->next = NULL;
  vnp->data.ptrvalue = BufferedReadFree ( (BufferedReadPtr)vnp->data.ptrvalue); 
  ValNodeFree (vnp);
}

/* three possible return codes:
 * 0 = no terminators seen at all
 * 1 = have terminator plus one character
 * 2 = last is terminator - need more characters
 */
static Int4 HasTerminator (ValNodePtr list, Int4 PNTR len)
{
  CharPtr      cp;
  ValNodePtr   vnp;
  BufferedReadPtr brp;

  if (len == NULL) return 0;
  *len = 0;
  if (list == NULL) return 0;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    if (vnp->data.ptrvalue == NULL) continue;
    brp = (BufferedReadPtr) vnp->data.ptrvalue;
    if (brp->data == NULL) continue;
    for (cp = brp->data + brp->offset; *cp != 0; cp++) {
      if (IsTerminator (*cp)) {
        if (* (cp + 1) != 0 || vnp->next != NULL) {
          return 1;
        } else {
          return 2;
        }
      } else { 
        (*len) ++;
      }
    }
  }
  return 0;
}

static CharPtr GetLineFromBuffer (ValNodePtr PNTR current_data, Int4 len)
{
  ValNodePtr      vnp, next_vnp;
  BufferedReadPtr brp;
  CharPtr         cp;
  CharPtr         new_line;
  Int4            ctr;
  Char            this_terminator;
  CharPtr         next_char;

  if (current_data == NULL || *current_data == NULL) return NULL;

  new_line = MemNew (len + 1);
  if (new_line == NULL) return NULL;

  ctr = 0;
  vnp = *current_data;
  while (vnp != NULL && ctr < len) {
    if ((brp = (BufferedReadPtr)vnp->data.ptrvalue) == NULL || brp->data == NULL) {
      next_vnp = vnp->next;
      vnp->next = NULL;
      vnp->data.ptrvalue = BufferedReadFree (brp);
      ValNodeFree (vnp);
      vnp = next_vnp;
    } else {
      if (ctr + brp->len <= len) {
        MemCpy (new_line + ctr, brp->data + brp->offset, brp->len);
        ctr += brp->len;
        next_vnp = vnp->next;
        vnp->next = NULL;
        vnp->data.ptrvalue = BufferedReadFree (brp);
        ValNodeFree (vnp);
        vnp = next_vnp;
      } else {
        MemCpy (new_line + ctr, brp->data + brp->offset, len - ctr);
        brp->offset += len - ctr;
        brp->len -= (len - ctr);
        ctr = len;
      }
    }
  }
  if (vnp != NULL) {
    brp = (BufferedReadPtr)vnp->data.ptrvalue;
    if (brp->len >= 0) {
      cp = brp->data + brp->offset;
      this_terminator = *cp;
      /* handle condition when last character in data is terminator */
      if (* (cp + 1) == 0) {
        next_vnp = vnp->next;
        vnp->next = NULL;
        vnp->data.ptrvalue = BufferedReadFree (brp);
        ValNodeFree (vnp);
        vnp = next_vnp;
        while (vnp != NULL && (brp = (BufferedReadPtr)vnp->data.ptrvalue) == NULL) {
          next_vnp = vnp->next;
          vnp->next = NULL;
          vnp->data.ptrvalue = BufferedReadFree (brp);
          ValNodeFree (vnp);
          vnp = next_vnp;
        }
        if (vnp == NULL) {
          *current_data = NULL;
          new_line [len] = 0;
          return new_line;
        } else {
          next_char = brp->data + brp->offset;
          if (IsTerminator (*next_char) && *next_char != this_terminator) {
            brp->offset ++;
            brp->len --;
            if (brp->len == 0) {
              next_vnp = vnp->next;
              vnp->next = NULL;
              vnp->data.ptrvalue = BufferedReadFree (brp);
              ValNodeFree (vnp);
              vnp = next_vnp;
            }
          }
        }
      } else {
        next_char = cp + 1;
        if (IsTerminator (*next_char) && *next_char != this_terminator) {
          brp->offset += 2;
          brp->len -= 2;
        } else {
          brp->offset ++;
          brp->len --;
        }
      }
      if (brp->len <= 0) {
        next_vnp = vnp->next;
        vnp->next = NULL;
        vnp->data.ptrvalue = BufferedReadFree (brp);
        ValNodeFree (vnp);
        vnp = next_vnp;
      }
    }
  }
  *current_data = vnp;
  new_line [len] = 0;
  return new_line;
}

#define READ_BUFFER_SIZE 5000

static ValNodePtr AddToBuffer (ValNodePtr current_data, FILE *fp)
{
  ValNodePtr vnp;
  BufferedReadPtr brp;

  vnp = ValNodeNew (current_data);
  if (vnp == NULL) return NULL;
 
  brp = (BufferedReadPtr) MemNew (sizeof (BufferedReadData));
  if (brp == NULL) return NULL;
  brp->data = MemNew (READ_BUFFER_SIZE);
  if (brp->data == NULL) return NULL;
  brp->offset = 0;
 
  brp->len = fread (brp->data, 1, READ_BUFFER_SIZE - 1, fp);
  *(char *)(brp->data + brp->len) = 0; 

  vnp->data.ptrvalue = brp;
  return vnp;
}

extern CharPtr MyFGetLine (FILE *fp, ValNodePtr PNTR current_data)
{
  Int4       terminator_status;
  Int4       data_len;
  ValNodePtr last_vnp;

  terminator_status = HasTerminator (*current_data, &data_len);
  while (!feof (fp) && terminator_status == 0) {
    last_vnp = AddToBuffer (*current_data, fp);
    if (*current_data == NULL) {
      *current_data = last_vnp;
    }
    terminator_status = HasTerminator (*current_data, &data_len);
  }

  if (!feof (fp) && terminator_status == 2) {
    AddToBuffer (*current_data, fp);
  }
  return GetLineFromBuffer (current_data, data_len);
} 

/* PCR_primer manipulation functions */

static ValNodePtr ParsePCRComponent (
  CharPtr strs
)

{
  ValNodePtr  head = NULL;
  size_t      len;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  if (tmp == NULL) return NULL;

  str = tmp;
  len = StringLen (str);
  if (len > 1 && *str == '(' && str [len - 1] == ')' && StringChr (str + 1, '(') == NULL) {
    str [len - 1] = '\0';
    str++;
  }

  while (StringDoesHaveText (str)) {
    ptr = StringChr (str, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }

    TrimSpacesAroundString (str);
    ValNodeCopyStr (&head, 0, str);

    str = ptr;
  }

  MemFree (tmp);
  return head;
}

NLM_EXTERN ValNodePtr ParsePCRStrings (
  CharPtr fwd_primer_seq,
  CharPtr rev_primer_seq,
  CharPtr fwd_primer_name,
  CharPtr rev_primer_name
)

{
  ValNodePtr  curr_fwd_name;
  ValNodePtr  curr_fwd_seq;
  ValNodePtr  curr_rev_name;
  ValNodePtr  curr_rev_seq;
  CharPtr     fwd_name;
  CharPtr     fwd_seq;
  CharPtr     rev_name;
  CharPtr     rev_seq;
  ValNodePtr  fwd_name_list = NULL;
  ValNodePtr  fwd_seq_list = NULL;
  ValNodePtr  rev_name_list = NULL;
  ValNodePtr  rev_seq_list = NULL;
  ValNodePtr  head = NULL;
  Boolean     okay;
  Int2        orig_order = 0;
  PcrSetPtr   psp;

  fwd_seq_list = ParsePCRComponent (fwd_primer_seq);
  rev_seq_list = ParsePCRComponent (rev_primer_seq);
  fwd_name_list = ParsePCRComponent (fwd_primer_name);
  rev_name_list = ParsePCRComponent (rev_primer_name);
    
  curr_fwd_seq = fwd_seq_list;
  curr_rev_seq = rev_seq_list;
  curr_fwd_name = fwd_name_list;
  curr_rev_name = rev_name_list;

  while (curr_fwd_seq != NULL || curr_rev_seq != NULL || curr_fwd_name != NULL || curr_rev_name != NULL) {
    fwd_seq = NULL;
    rev_seq = NULL;
    fwd_name = NULL;
    rev_name = NULL;
    okay = FALSE;

    if (curr_fwd_seq != NULL) {
      fwd_seq = (CharPtr) curr_fwd_seq->data.ptrvalue;
      curr_fwd_seq = curr_fwd_seq->next;
      okay = TRUE;
    }

    if (curr_rev_seq != NULL) {
      rev_seq = (CharPtr) curr_rev_seq->data.ptrvalue;
      curr_rev_seq = curr_rev_seq->next;
      okay = TRUE;
    }

    if (curr_fwd_name != NULL) {
      fwd_name = (CharPtr) curr_fwd_name->data.ptrvalue;
      curr_fwd_name = curr_fwd_name->next;
      okay = TRUE;
    }

    if (curr_rev_name != NULL) {
      rev_name = (CharPtr) curr_rev_name->data.ptrvalue;
      curr_rev_name = curr_rev_name->next;
      okay = TRUE;
    }

    if (okay) {
      psp = (PcrSetPtr) MemNew (sizeof (PcrSet));
      if (psp != NULL) {
        psp->fwd_seq = StringSaveNoNull (fwd_seq);
        psp->rev_seq = StringSaveNoNull (rev_seq);
        psp->fwd_name = StringSaveNoNull (fwd_name);
        psp->rev_name = StringSaveNoNull (rev_name);
        orig_order++;
        psp->orig_order = orig_order;
        ValNodeAddPointer (&head, 0, (Pointer) psp);
      }
    }
  }

  ValNodeFreeData (fwd_seq_list);
  ValNodeFreeData (rev_seq_list);
  ValNodeFreeData (fwd_name_list);
  ValNodeFreeData (rev_name_list);

  return head;
}

NLM_EXTERN ValNodePtr ParsePCRSet (
  BioSourcePtr biop
)

{
  CharPtr       fwd_primer_seq = NULL;
  CharPtr       rev_primer_seq = NULL;
  CharPtr       fwd_primer_name = NULL;
  CharPtr       rev_primer_name = NULL;
  SubSourcePtr  ssp;

  if (biop == NULL) return NULL;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      fwd_primer_seq = ssp->name;
    } else if (ssp->subtype == SUBSRC_rev_primer_seq) {
      rev_primer_seq = ssp->name;
    } else if (ssp->subtype == SUBSRC_fwd_primer_name) {
      fwd_primer_name = ssp->name;
    } else if (ssp->subtype == SUBSRC_rev_primer_name) {
      rev_primer_name = ssp->name;
    }
  }

  return ParsePCRStrings (fwd_primer_seq, rev_primer_seq, fwd_primer_name, rev_primer_name);
}

NLM_EXTERN int LIBCALLBACK SortVnpByPCRSetSeq (VoidPtr ptr1, VoidPtr ptr2)

{
  int         compare;
  PcrSetPtr   psp1, psp2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  psp1 = (PcrSetPtr) vnp1->data.ptrvalue;
  psp2 = (PcrSetPtr) vnp2->data.ptrvalue;
  if (psp1 == NULL || psp2 == NULL) return 0;

  compare = StringICmp (psp1->fwd_seq, psp2->fwd_seq);
  if (compare != 0) return compare;

  compare = StringICmp (psp1->rev_seq, psp2->rev_seq);
  if (compare != 0) return compare;

  compare = StringICmp (psp1->fwd_name, psp2->fwd_name);
  if (compare != 0) return compare;

  compare = StringICmp (psp1->rev_name, psp2->rev_name);
  if (compare != 0) return compare;

  if (psp1->orig_order > psp2->orig_order) {
    return 1;
  } else if (psp1->orig_order < psp2->orig_order) {
    return -1;
  }

  return 0;
}

NLM_EXTERN ValNodePtr UniqueVnpByPCRSetSeq (ValNodePtr pset)

{
  PcrSetPtr     last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  PcrSetPtr     psp;
  ValNodePtr    vnp;

  if (pset == NULL) return NULL;
  last = (PcrSetPtr) pset->data.ptrvalue;
  vnp = pset->next;
  prev = (Pointer PNTR) &(pset->next);
  while (vnp != NULL) {
    next = vnp->next;
    psp = (PcrSetPtr) vnp->data.ptrvalue;
    if (last != NULL && psp != NULL &&
        StringICmp (last->fwd_seq, psp->fwd_seq) == 0 &&
        StringICmp (last->rev_seq, psp->rev_seq) == 0 &&
        StringICmp (last->fwd_name, psp->fwd_name) == 0 &&
        StringICmp (last->rev_name, psp->rev_name) == 0) {
      vnp->next = NULL;
      *prev = next;
      MemFree (psp->fwd_seq);
      MemFree (psp->rev_seq);
      MemFree (psp->fwd_name);
      MemFree (psp->rev_name);
      ValNodeFreeData (vnp);
    } else {
      last = (PcrSetPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return pset;
}

NLM_EXTERN int LIBCALLBACK SortVnpByPCRSetOrder (VoidPtr ptr1, VoidPtr ptr2)

{
  PcrSetPtr   psp1, psp2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  psp1 = (PcrSetPtr) vnp1->data.ptrvalue;
  psp2 = (PcrSetPtr) vnp2->data.ptrvalue;
  if (psp1 == NULL || psp2 == NULL) return 0;

  if (psp1->orig_order > psp2->orig_order) {
    return 1;
  } else if (psp1->orig_order < psp2->orig_order) {
    return -1;
  }

  return 0;
}

static CharPtr CombinePCRItems (
  ValNodePtr list
)

{
  Int4        count;
  size_t      len;
  CharPtr     ptr;
  CharPtr     str;
  ValNodePtr  vnp;

  if (list == NULL) return NULL;
  count = ValNodeLen (list);
  if (count == 1) {
    ptr = (CharPtr) list->data.ptrvalue;
    return StringSaveNoNull (ptr);
  }

  len = 0;
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    ptr = (CharPtr) vnp->data.ptrvalue;
    if (ptr == NULL) continue;
    len += StringLen (ptr) + 1;
  }
  str = (CharPtr) MemNew (sizeof (Char) * (len + 4));
  if (str == NULL) return NULL;
  StringCpy (str, "(");

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    ptr = (CharPtr) vnp->data.ptrvalue;
    if (ptr == NULL) continue;
    StringCat (str, ptr);
    if (vnp->next != NULL) {
      StringCat (str, ",");
    }
  }

  StringCat (str, ")");
  return str;
}

NLM_EXTERN SubSourcePtr WritePCRSet (
  ValNodePtr pset
)

{
  ValNodePtr    fwd_name_list = NULL;
  ValNodePtr    fwd_seq_list = NULL;
  ValNodePtr    rev_name_list = NULL;
  ValNodePtr    rev_seq_list = NULL;
  SubSourcePtr  head = NULL;
  SubSourcePtr  last = NULL;
  PcrSetPtr     psp;
  SubSourcePtr  ssp;
  CharPtr       str;
  ValNodePtr    vnp;

  if (pset == NULL) return NULL;

  for (vnp = pset; vnp != NULL; vnp = vnp->next) {
    psp = (PcrSetPtr) vnp->data.ptrvalue;
    if (psp == NULL) continue;
    if (StringDoesHaveText (psp->fwd_seq)) {
      ValNodeCopyStr (&fwd_seq_list, 0, psp->fwd_seq);
    }
    if (StringDoesHaveText (psp->rev_seq)) {
      ValNodeCopyStr (&rev_seq_list, 0, psp->rev_seq);
    }
    if (StringDoesHaveText (psp->fwd_name)) {
      ValNodeCopyStr (&fwd_name_list, 0, psp->fwd_name);
    }
    if (StringDoesHaveText (psp->rev_name)) {
      ValNodeCopyStr (&rev_name_list, 0, psp->rev_name);
    }
  }

  str = CombinePCRItems (fwd_seq_list);
  if (str != NULL) {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_fwd_primer_seq;
    ssp->name = str;
    if (head == NULL) {
      head = ssp;
    }
    if (last != NULL) {
      last->next = ssp;
    }
    last = ssp;
  }

  str = CombinePCRItems (rev_seq_list);
  if (str != NULL) {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_rev_primer_seq;
    ssp->name = str;
    if (head == NULL) {
      head = ssp;
    }
    if (last != NULL) {
      last->next = ssp;
    }
    last = ssp;
  }

  str = CombinePCRItems (fwd_name_list);
  if (str != NULL) {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_fwd_primer_name;
    ssp->name = str;
    if (head == NULL) {
      head = ssp;
    }
    if (last != NULL) {
      last->next = ssp;
    }
    last = ssp;
  }

  str = CombinePCRItems (rev_name_list);
  if (str != NULL) {
    ssp = SubSourceNew ();
    ssp->subtype = SUBSRC_rev_primer_name;
    ssp->name = str;
    if (head == NULL) {
      head = ssp;
    }
    if (last != NULL) {
      last->next = ssp;
    }
    last = ssp;
  }

  return head;
}

NLM_EXTERN ValNodePtr FreePCRSet (
  ValNodePtr pset
)

{
  PcrSetPtr   psp;
  ValNodePtr  vnp;

  if (pset == NULL) return NULL;

  for (vnp = pset; vnp != NULL; vnp = vnp->next) {
    psp = (PcrSetPtr) vnp->data.ptrvalue;
    if (psp == NULL) continue;
    MemFree (psp->fwd_seq);
    MemFree (psp->rev_seq);
    MemFree (psp->fwd_name);
    MemFree (psp->rev_name);
  }

  return ValNodeFreeData (pset);
}




static void AddDefLinesToAlignmentSequences 
(TAlignmentFilePtr afp,
 SeqEntryPtr sep_head)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  sep;
  Int4         index;
  ValNodePtr   sdp;
  CharPtr      new_title;
  Uint4         new_title_len;
  Int4         curr_seg;
  Int4         num_sets = 1;
  Boolean      one_defline_per_sequence = TRUE;
  Boolean      all_extra_empty;

 
  if (afp == NULL || sep_head == NULL || ! IS_Bioseq_set (sep_head))
  {
    return;
  }
  if ((afp->num_deflines == 0 || afp->deflines == NULL)
    && (afp->num_organisms == 0 || afp->organisms == NULL))
  {
    return;
  }
  bssp = sep_head->data.ptrvalue;
  
  /* find out if all of our deflines are real */
  if (afp->num_segments > 1 && afp->num_deflines == afp->num_sequences)
  {
    one_defline_per_sequence = FALSE;
    num_sets = afp->num_sequences / afp->num_segments;
    all_extra_empty = TRUE;
    for (curr_seg = num_sets; curr_seg < afp->num_deflines && all_extra_empty; curr_seg ++)
    {
      if (afp->deflines [curr_seg] != NULL)
      {
        all_extra_empty = FALSE;
      }
    }
      if (all_extra_empty)
      {
          one_defline_per_sequence = TRUE;
      }
  }
  
  for (sep = bssp->seq_set, index = 0;
       sep != NULL && (index < afp->num_deflines || index < afp->num_organisms);
       sep = sep->next, index++)
  {
    new_title_len = 0;
    /* get lengths for organisms for this sequence */
    
    if (afp->num_segments > 1 && afp->num_organisms == afp->num_sequences)
    {
      /* have one organism per segment, in which case use only the first one */
      curr_seg = index * afp->num_segments;
    }
    else
    { /* otherwise one organism per sequence */
      curr_seg = index;
    }
    if (curr_seg < afp->num_organisms)
    {
      new_title_len += StringLen (afp->organisms [curr_seg]) + 1;
    }

    /* get lengths for deflines for this sequence */
    if (! one_defline_per_sequence)
    { /* have one defline per segment, in which use only the first one */
      curr_seg = index * afp->num_segments;
    }
    else
    { /* otherwise one defline per sequence */
      curr_seg = index;    
    }
    if (curr_seg < afp->num_deflines) 
    {
      new_title_len += StringLen (afp->deflines [curr_seg]) + 1;
    }
      
    if (new_title_len > 0) {
      new_title = (CharPtr) MemNew (new_title_len);
      if (new_title == NULL) return;
      new_title [0] = 0;
      
      /* list organisms at beginning of new defline */
      if (afp->num_segments > 1 && afp->num_organisms == afp->num_sequences)
      { /* have one organism per segment, in which case use only first one */
        curr_seg = index * afp->num_segments;
      }
      else
      { /* otherwise one organism per sequence */
          curr_seg = index;
      }

      if (curr_seg < afp->num_organisms) {
        StringCat (new_title, afp->organisms [curr_seg]);
        if (new_title_len > StringLen (new_title) + 1)
        {
          StringCat (new_title, " ");
        }
      }
      
      if (!one_defline_per_sequence)
      { /* have one defline per segment, in which case all go to same sequence */
        curr_seg = index * afp->num_segments;
      }
      else 
      {
          curr_seg = index;
      }
      if (curr_seg < afp->num_deflines) 
      {
        StringCat (new_title, afp->deflines [curr_seg]);
      }

      sdp = CreateNewDescriptor (sep, Seq_descr_title);
      if (sdp != NULL) {
        sdp->data.ptrvalue = new_title;
      } else {
        MemFree (new_title);
      }
    }
  }
}

#if 0
static SeqEntryPtr 
MakeDeltaSetFromAlignment 
(SeqEntryPtr sep_list,
 TAlignmentFilePtr afp,
 Uint1 moltype,
 Int4  gap_length
 )
{
  BioseqPtr    bsp, deltabsp;
  SeqEntryPtr  this_list, last_sep, next_list, sep, nextsep;
  SeqEntryPtr  topsep, last_delta_sep;
  SeqIdPtr     sip;
  Int4         curr_seg;
  CharPtr      seqbuf;
  ValNodePtr   vnp;
  SeqLitPtr    slp;
  IntFuzzPtr   ifp;
  SeqEntryPtr  delta_list = NULL;
  
  delta_list = NULL;
  last_delta_sep = NULL;
  this_list = sep_list;
  while (this_list != NULL)
  {
    last_sep = this_list;
    curr_seg = 0;
    while (last_sep != NULL && curr_seg < afp->num_segments - 1)
    {
        last_sep = last_sep->next;
        curr_seg++;
    }
    if (last_sep == NULL) return NULL;
    next_list = last_sep->next;
    last_sep->next = NULL;
  
    bsp = (BioseqPtr)this_list->data.ptrvalue;
    if (bsp == NULL) return NULL;

    sip = SeqIdDup (bsp->id);
    vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);

    deltabsp = BioseqNew ();
    if (deltabsp == NULL) return NULL;
    deltabsp->repr = Seq_repr_delta;
    deltabsp->seq_ext_type = 4;
    deltabsp->mol = moltype;
    deltabsp->length = 0;

    topsep = SeqEntryNew ();
    if (topsep == NULL) return NULL;
    topsep->choice = 1;
    topsep->data.ptrvalue = (Pointer) deltabsp;

    for (sep = this_list; sep != NULL; sep = nextsep) {
      nextsep = sep->next;
      sep->next = NULL;

      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp == NULL) continue;

      if (bsp->repr == Seq_repr_raw) {
        BioseqRawConvert (bsp, Seq_code_iupacna);
        seqbuf = BSMerge ((ByteStorePtr) bsp->seq_data, NULL);
        slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slp == NULL) continue;

        slp->length = bsp->length;
        ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);
        slp->seq_data = BSNew (slp->length);
        slp->seq_data_type = Seq_code_iupacna;
        AddBasesToByteStore (slp->seq_data, seqbuf);
        MemFree(seqbuf);

        deltabsp->length += slp->length;

      } else if (bsp->repr == Seq_repr_virtual) {
        slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slp == NULL) continue;
        slp->length = bsp->length;
        if (slp == NULL) continue;

        slp->length = bsp->length;
        ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);
        if (slp->length < 1) {
          slp->length = 0;
          ifp = IntFuzzNew ();
          ifp->choice = 4;
          slp->fuzz = ifp;
        }

        deltabsp->length += slp->length;
      }
      SeqEntryFree (sep);
      
      if (nextsep != NULL)
      {
        /* add gap */
        slp = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slp == NULL) continue;
        slp->length = gap_length;
        ValNodeAddPointer ((ValNodePtr PNTR) &(deltabsp->seq_ext), (Int2) 2, (Pointer) slp);
        deltabsp->length += slp->length;        
      }
    }

    ValNodeLink (&(deltabsp->descr), vnp);
    deltabsp->id = sip;
    
    if (last_delta_sep == NULL)
    {
        delta_list = topsep;
    }
    else 
    {
        last_delta_sep->next = topsep;    
    }
    last_delta_sep = topsep;
    
    this_list = next_list;    
  }
  return delta_list;
}
#endif

static void RenameSegSet (SeqEntryPtr sep)
{
  BioseqSetPtr bssp, seg_bssp;
  SeqEntryPtr  seg_sep;
  BioseqPtr    main_bsp = NULL;
  BioseqPtr    seg_bsp = NULL;
  Char         new_id_str [255];
  
  if (sep == NULL || !IS_Bioseq_set (sep) || (bssp = sep->data.ptrvalue) == NULL
      || bssp->_class != BioseqseqSet_class_segset)
  {
      return;
  }
  
  sep = bssp->seq_set;
  while (sep != NULL && (seg_bsp == NULL || main_bsp == NULL))
  {
      if (IS_Bioseq (sep))
      {
        main_bsp = (BioseqPtr) sep->data.ptrvalue;
      }
      else if (IS_Bioseq_set (sep))
      {
        seg_bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (seg_bssp != NULL && seg_bssp->_class == BioseqseqSet_class_parts)
        {
            seg_sep = seg_bssp->seq_set;
            while (seg_sep != NULL && seg_bsp == NULL)
            {
              if (IS_Bioseq (seg_sep))
              {
                seg_bsp = seg_sep->data.ptrvalue;
            }
            seg_sep = seg_sep->next;
          }
        }
      }
      sep = sep->next;
  }
  if (main_bsp == NULL || seg_bsp == NULL)
  {
      return;
  }
  SeqIdWrite (seg_bsp->id, new_id_str, PRINTID_FASTA_SHORT, sizeof (new_id_str) - 7);
  StringCat (new_id_str, "_master");
  SeqIdFree (main_bsp->id);
  main_bsp->id = MakeSeqID (new_id_str);
}

static SeqEntryPtr 
MakeSegmentedSetFromAlignment 
(SeqEntryPtr       sep_list,
 TAlignmentFilePtr afp,
 Uint1             moltype,
 Int4Ptr           segs_per_set)
{
  SeqEntryPtr  this_list, last_sep, next_list, nextsep, last_segset;
  Int4         curr_seg;
  Int4         set_index = 0;
  
  this_list = sep_list;
  sep_list = NULL;
  last_segset = NULL;
  while (this_list != NULL)
  {
    last_sep = this_list;
    curr_seg = 0;
    while (last_sep != NULL && curr_seg < segs_per_set [set_index] - 1)
    {
      if (!IS_Bioseq (last_sep)) return NULL;
      last_sep = last_sep->next;
      curr_seg++;
    }
    if (last_sep == NULL) return NULL;
    next_list = last_sep->next;
    last_sep->next = NULL;
    
    last_sep = this_list->next;
    this_list->next = NULL;
    while (last_sep != NULL)
    {
      nextsep = last_sep->next;
      last_sep->next = NULL;
      AddSeqEntryToSeqEntry (this_list, last_sep, FALSE);
      last_sep = nextsep;
    }

    /* fix IDs for seg sets */    
    RenameSegSet (this_list);
    
    if (sep_list == NULL) 
    {
      sep_list = this_list;
    }
    else
    {
      last_segset->next = this_list;
    }
    last_segset = this_list;
    
    this_list = next_list;
    set_index++;
  }
  return sep_list;     
}


extern CharPtr AlignmentStringToSequenceString (CharPtr aln_str, Uint1 moltype)
{
  CharPtr cp_aln, cp_seq;
  Char    ch;
  CharPtr seq_str;
  
  if (aln_str == NULL) return NULL;
  seq_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (aln_str) + 1));
  if (seq_str == NULL) return NULL;
  cp_seq = seq_str;
  for (cp_aln = aln_str; *cp_aln != 0; cp_aln++)
  {
    ch = *cp_aln;
    ch = TO_UPPER (ch); 
    if ( ISA_na (moltype) ) 
    {
      if (ch == 'U') ch = 'T';
      if (ch == 'X') ch = 'N';
      if ( StringChr ("EFIJLOPQXZ-.*", ch) == NULL )  
      { 
        *cp_seq = ch;
        cp_seq++;
      }
    }
    else 
    {
      if ( StringChr("JO-.", ch) == NULL ) 
      {
        *cp_seq = ch;
        cp_seq++;
      } 
    }
  }
  *cp_seq = 0;
  return seq_str;
}

static SeqEntryPtr SequenceStringToSeqEntry (CharPtr str, SeqIdPtr sip, Uint1 mol_type)
{
  SeqEntryPtr  sep;
  BioseqPtr    bsp;
  ByteStorePtr bs;

  if (str == NULL || sip == NULL) return NULL;
  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;
  bsp = BioseqNew ();
  if (bsp == NULL) 
  { 
    ValNodeFree (sep); 
    return NULL; 
  }
  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) bsp;
  bsp->id = SeqIdDup (sip);
  bsp->id->next = NULL;
  SeqMgrAddToBioseqIndex (bsp);
  bsp->repr = Seq_repr_raw;
  if ( ISA_na (mol_type) ) 
  {
    bsp->mol = Seq_mol_na;
    bsp->seq_data_type = Seq_code_iupacna;
  } 
  else
  {
    bsp->mol = Seq_mol_aa;
    bsp->seq_data_type = Seq_code_ncbieaa;
  }
  bsp->length = StringLen (str);
  if ( bsp->length == 0 ) 
  {
    BioseqFree (bsp);
    ValNodeFree (sep);
    return NULL;
  }
  bs = BSNew (bsp->length);
  bsp->seq_data = (SeqDataPtr) bs;
  BSWrite (bs, str, bsp->length);

  return sep;
}

#if 0
static SeqEntryPtr MakeDeltaSeqsFromAlignmentSequences (TAlignmentFilePtr afp, Uint1 moltype, CharPtr PNTR seq_str)
{
  Int4            num_sets, next_start, k, index;
  SeqIdPtr        sip;
  SeqLitPtr       slip;
  SeqEntryPtr     sep_list = NULL, sep, sep_last = NULL;
  BioseqPtr       new_bsp;
  ValNodePtr      seq_ext = NULL;
  
  if (afp == NULL || seq_str == NULL) return NULL;
  
  num_sets = afp->num_sequences / afp->num_segments;
  for (k = 0; k < num_sets; k++)
  {
    sep = SeqEntryNew ();
    if (sep == NULL) return NULL;
    new_bsp = BioseqNew ();
    if (new_bsp == NULL) return NULL;
    sip = MakeSeqID (afp->ids [k * afp->num_segments]);
    new_bsp->id = sip;
    sep->choice = 1;
    sep->data.ptrvalue = new_bsp;
    SeqMgrAddToBioseqIndex (new_bsp);
    
    if (sep_last == NULL)
    {
      sep_list = sep;
    }
    else
    {
      sep_last->next = sep;
    }
    sep_last = sep;
    
    new_bsp->seq_data = NULL;
    new_bsp->seq_data_type = 0;
    new_bsp->repr = Seq_repr_delta;
    new_bsp->seq_ext_type = 4;
    new_bsp->mol = moltype;
    new_bsp->seq_ext = NULL;
    new_bsp->length = 0;
    next_start = (k + 1) * afp->num_segments;
    seq_ext = NULL;
    for (index = k * afp->num_segments; index < next_start; index++)
    {
      if (seq_ext != NULL)
      {
        /* insert gap of unknown length between the previous segment
         * and this one.
         */
        new_bsp->length += AddGapSeqLit (&seq_ext);
      }

      if (StringHasNoText (seq_str [index]))
      {
        /* add gap to represent missing sequence */
        new_bsp->length += AddGapSeqLit (&seq_ext);        
      }
      else
      {
        slip = (SeqLitPtr) MemNew (sizeof (SeqLit));
        if (slip != NULL) 
        {
          slip->length = StringLen (seq_str [index]);
          ValNodeAddPointer (&seq_ext, (Int2) 2, (Pointer) slip);
          slip->seq_data = BSNew (slip->length);
          slip->seq_data_type = Seq_code_iupacna;
          AddBasesToByteStore (slip->seq_data, seq_str [index]);
          new_bsp->length += slip->length;
        } 
      }
    }
    new_bsp->seq_ext = seq_ext;
    BioseqPack (new_bsp);        
  }
    
  return sep_list;
}
#endif

static SeqIdPtr GetFarPointerID (CharPtr id_str)
{
  CharPtr  tmp_id_str;
  CharPtr  cp_start, cp_end;
  Int4     len;
  SeqIdPtr sip;
  
  if (id_str == NULL)
  {
    return NULL;
  }

  cp_start = StringChr (id_str, '|');
  if (cp_start == NULL)
  {
    cp_start = id_str;
    len = StringLen (id_str);
  }
  else
  {
    cp_start++;
    cp_end = StringChr (cp_start, '|');
    if (cp_end == NULL)
    {
      len = StringLen (cp_start);
    }
    else
    {
      len = cp_end - cp_start;
    }
  }
  if (len == 0)
  {
    return NULL;
  }
  tmp_id_str = (CharPtr) MemNew ((len + 4) * sizeof (Char));
  if (tmp_id_str == NULL)
  {
    return NULL;
  }
  StringCpy (tmp_id_str, "acc");
  StringNCat (tmp_id_str, cp_start, len);
  tmp_id_str [len + 3] = 0;
  sip = MakeSeqID (tmp_id_str);
  MemFree (tmp_id_str);
  return sip;
}

static void ReplacePipesWithUnderscores (CharPtr seqid_str)
{
  CharPtr cp;
  
  if (seqid_str == NULL)
  {
    return;
  }
  
  cp = seqid_str;
  while (*cp != 0)
  {
    if (*cp == '|')
    {
      *cp = '_';
    }
    cp++;
  }
}

extern SeqEntryPtr MakeSequinDataFromAlignmentEx (TAlignmentFilePtr afp, Uint1 moltype, Boolean check_ids) 
{
  SeqIdPtr    PNTR sip_list;
  SeqIdPtr    PNTR sip_prev;
  SeqAnnotPtr sap = NULL;
  SeqAlignPtr salp_list, salp_last;
  ValNodePtr  PNTR seqvnp;
  SeqEntryPtr sep_list;
  SeqEntryPtr sep, sep_prev;
  SeqIdPtr    sip;
  ValNodePtr  vnp;
  Int4        index, curr_seg, num_sets;
  BioseqPtr   bsp;
  CharPtr     tmp_id_str;
  MsgAnswer   ans;
  Int4Ptr      segs_per_set = NULL;
  Int4Ptr      segs_per_aln = NULL;
  Boolean      found_empty_seg = FALSE;
  CharPtr      seq_data = NULL;

  if (afp == NULL) return NULL;
  
  if (afp->num_sequences == 0) return NULL;
  if (afp->num_segments < 1) return NULL;
  
  sip_list = (SeqIdPtr PNTR) MemNew (afp->num_segments * sizeof (SeqIdPtr));  
  sip_prev = (SeqIdPtr PNTR) MemNew (afp->num_segments * sizeof (SeqIdPtr));  
  seqvnp = (ValNodePtr PNTR) MemNew (afp->num_segments * sizeof (ValNodePtr));  
  segs_per_set = (Int4Ptr) MemNew (sizeof (Int4Ptr) * afp->num_sequences);
  segs_per_aln = (Int4Ptr) MemNew (sizeof (Int4Ptr) * afp->num_segments);
  if (sip_list == NULL || sip_prev == NULL || seqvnp == NULL
      || segs_per_set == NULL || segs_per_aln == NULL)
  {
    MemFree (sip_list);
    MemFree (sip_prev);
      MemFree (seqvnp);
      MemFree (segs_per_set);
      MemFree (segs_per_aln);
      return NULL;
  }
 
  for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg ++)
  {
    sip_list [curr_seg] = NULL;
    sip_prev [curr_seg] = NULL;
      seqvnp [curr_seg] = NULL;
      segs_per_aln [curr_seg] = 0;
  }

  sep_list = NULL;
  sep_prev = NULL;
  curr_seg = 0;

  for (index = 0; index < afp->num_sequences; index++) {
    seq_data = AlignmentStringToSequenceString (afp->sequences [index], moltype);
    if (StringHasNoText (seq_data))
    {
      found_empty_seg = TRUE;
    }
    else
    {
      sip = MakeSeqID (afp->ids [index]);
      if (sip == NULL && StringChr (afp->ids [index], '|') != NULL)
      {
        ReplacePipesWithUnderscores (afp->ids [index]);
        sip = MakeSeqID (afp->ids [index]);
      }
      if (sip != NULL)
      {
        sip->next = SeqIdFree (sip->next);
      }
      if (check_ids && StringNCmp (afp->ids[index], "acc", 3) != 0)
      {
        bsp = BioseqFind (sip);
        if (bsp == NULL)
        {
          sip = SeqIdFree (sip);
          tmp_id_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (afp->ids [index]) + 4));
          sprintf (tmp_id_str, "gb|%s", afp->ids [index]);
          sip = MakeSeqID (tmp_id_str);
          MemFree (tmp_id_str);
          bsp = BioseqFind (sip);
        }
        if (bsp == NULL)
        {
          ans = Message (MSG_YN, "Can't find sequence %s in set - is this a far pointer?", afp->ids[index]);
          if (ans == ANS_YES)
          {
            sip = SeqIdFree (sip);
            sip = GetFarPointerID (afp->ids [index]);
          }
          else
          {
            sip = SeqIdFree (sip);
            sip = MakeSeqID (afp->ids [index]);
          }
          if (sip != NULL)
          {
            sip->next = SeqIdFree (sip->next);
          }
        }
      }

      sep = SequenceStringToSeqEntry (seq_data, sip, moltype);
      if (sep != NULL) {
        if (sep_list == NULL) {
          sep_list = sep;
        } else {
          sep_prev->next = sep;
        }
        sep_prev = sep;
        vnp = ValNodeNew (seqvnp[curr_seg]);
        if (seqvnp[curr_seg] == NULL) seqvnp[curr_seg] = vnp;
        vnp->data.ptrvalue = afp->sequences [index];
      
        /* only add SeqID to list if adding segment */
        if (sip_prev[curr_seg] == NULL) {
          sip_list[curr_seg] = sip;
        } else {
          sip_prev[curr_seg]->next = sip;
        }
        sip_prev[curr_seg] = sip;
        
        /* add to totals for this set and for this alignment */
        segs_per_set [index / afp->num_segments] ++;
        segs_per_aln [index % afp->num_segments] ++;
      }
    }
    seq_data = MemFree (seq_data);
    curr_seg ++;
    if (curr_seg >= afp->num_segments) 
    {
      curr_seg = 0;
    }
  }

  if (found_empty_seg)
  {
    Boolean   indexerVersion;
    MsgAnswer ans = ANS_YES;
    
    if (afp->num_segments > 1)
    {
      indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
      if (indexerVersion)
      {
        ans = Message (MSG_YN, "This alignment of segmented sets contains a segment that is all gaps - do you wish to continue?");
      }
    }
    else
    {
      Message (MSG_ERROR, "This alignment contains a sequence that is all gaps.");
      ans = ANS_NO;
    }
    if (ans == ANS_NO)
    {
      for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg ++)
      {
        ValNodeFree (seqvnp [curr_seg]);
      }
      MemFree (seqvnp);
      MemFree (sip_list);
      MemFree (sip_prev);
      MemFree (segs_per_set);
      MemFree (segs_per_aln);
      sep_list = SeqEntryFree (sep_list);
      return NULL;
    }
  }

  
  if (afp->num_segments == 1) 
  {
    sap = LocalAlignToSeqAnnotDimn (seqvnp[0], sip_list[0], NULL, afp->num_sequences,
                                    0, NULL, FALSE);
    sep_list = make_seqentry_for_seqentry (sep_list);
    SeqAlignAddInSeqEntry (sep_list, sap);      
  } 
  else 
  {
    sep_list = MakeSegmentedSetFromAlignment (sep_list, afp, moltype, segs_per_set);
    sep_list = make_seqentry_for_seqentry (sep_list);
    num_sets = afp->num_sequences / afp->num_segments;
    salp_list = NULL;
    salp_last = NULL;

    for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg++)
    {      
      sap = LocalAlignToSeqAnnotDimn (seqvnp[curr_seg], sip_list[curr_seg], NULL, segs_per_aln [curr_seg],
                                    0, NULL, FALSE);
      if (sap != NULL)
      {
        SeqAlignAddInSeqEntry (sep_list, sap);
      }
    }
  }

  for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg ++)
  {
    ValNodeFree (seqvnp [curr_seg]);
  }
  MemFree (seqvnp);
  MemFree (sip_list);
  MemFree (sip_prev);
  MemFree (segs_per_set);
  MemFree (segs_per_aln);

  AddDefLinesToAlignmentSequences (afp, sep_list);

  return sep_list;
}

extern SeqEntryPtr MakeSequinDataFromAlignment (TAlignmentFilePtr afp, Uint1 moltype) 
{
  return MakeSequinDataFromAlignmentEx (afp, moltype, FALSE);
}

/* Create sequences and alignment annotation */

/**********************************************************/
extern SeqEntryPtr make_seqentry_for_seqentry (SeqEntryPtr sep)
{
  SeqEntryPtr  sep1 = NULL,
               tmp;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
 
  if (sep == NULL) return NULL;

  if (! IS_Bioseq (sep) && ! IS_Bioseq_set (sep)) {
    return sep;
  } else if (sep->next == NULL) {
    return sep;
  } else if ((bssp = BioseqSetNew ()) == NULL) {
    return sep;
  } else {
    bssp->_class = 14;
    bssp->seq_set = sep;
    sep1 = SeqEntryNew ();
    sep1->choice = 2;
    sep1->data.ptrvalue = bssp;
    SeqMgrLinkSeqEntry (sep1, 0, NULL);
          
    for (tmp = bssp->seq_set; tmp!=NULL; tmp=tmp->next) {
      if (IS_Bioseq(tmp)) {
        bsp = (BioseqPtr) tmp->data.ptrvalue;
        ObjMgrConnect (OBJ_BIOSEQ, (Pointer) bsp, OBJ_BIOSEQSET, (Pointer) bssp);
      }
    }
  }
  return sep1;
}

/* These two functions are used for removing mRNAs that overlap pseudo misc_feats
 * and marking genes that overlap pseudo misc_feats as pseudo.
 */
static void PseudoMiscFeatProcessingCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SeqFeatPtr        gene, mRNA;
  SeqMgrFeatContext gcontext, mcontext;
  
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_misc_feature) return;
  /* we only want to process misc_feats if the pseudo flag is set or the 
   * comment contains the word "pseudogene".
   */
#if 0
  if (!sfp->pseudo && StringISearch (sfp->comment, "pseudogene") == NULL) return;
#endif

  gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
  mRNA = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL,
                                      RANGE_MATCH, &mcontext);
  if (gene != NULL)
  {
      gene->pseudo = TRUE;
  }
  if (mRNA != NULL && mRNA->product == NULL) /* only delete mRNAs without products */
  {
      mRNA->idx.deleteme = TRUE;
  }  
}

extern void ProcessPseudoMiscFeatsForEntityID (Uint2 entityID)
{
  SeqEntryPtr sep;
  
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  
  VisitFeaturesInSep (sep, (Pointer) NULL, PseudoMiscFeatProcessingCallback);
  DeleteMarkedObjects (entityID, 0, NULL);      
}

/* These three functions are used for converting pseudo CDSs to misc_features. */
extern Boolean ConvertOnePseudoCDSToMiscFeat (SeqFeatPtr sfp)
{
  BioseqPtr  bsp;
  SeqFeatPtr new_sfp;
  ImpFeatPtr ifp;
  
  if (sfp == NULL || (sfp->data.choice != SEQFEAT_CDREGION) || (! sfp->pseudo)) return FALSE;
  
  bsp = BioseqFindFromSeqLoc (sfp->location);  
  if (bsp == NULL) return FALSE;
  ifp = ImpFeatNew ();
  if (ifp == NULL) return FALSE;
  new_sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_IMP, sfp->location);
  if (new_sfp == NULL) 
  {
      ImpFeatFree (ifp);
      return FALSE;
  }
  new_sfp->data.value.ptrvalue = (Pointer) ifp;
  ifp->key = StringSave ("misc_feature");
  new_sfp->comment = sfp->comment;
  sfp->comment = NULL;
  new_sfp->qual = sfp->qual;
  sfp->qual = NULL;
  
  if (sfp->product != NULL)
  {
      bsp = BioseqFindFromSeqLoc (sfp->product);
      sfp->product = SeqLocFree (sfp->product);
      bsp->idx.deleteme = TRUE;
  }
  sfp->idx.deleteme = TRUE;
  return TRUE;
}

static void ConvertPseudoCDSToMiscFeatCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ConvertOnePseudoCDSToMiscFeat (sfp);
}

extern void ConvertPseudoCDSToMiscFeatsForEntityID (Uint2 entityID)
{
  SeqEntryPtr sep;
  
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  
  VisitFeaturesInSep (sep, (Pointer) NULL, ConvertPseudoCDSToMiscFeatCallback);
  DeleteMarkedObjects (entityID, 0, NULL);      
}

typedef struct alignmentforbsp
{
  BioseqPtr   bsp;
  SeqAlignPtr salp_list;
  SeqAlignPtr salp_last;
  ValNodePtr  seq_annot_list;
} AlignmentForBspData, PNTR AlignmentForBspPtr;

static void FindAlignmentsForBioseqCallback (SeqAnnotPtr sap, Pointer userdata)
{
  AlignmentForBspPtr   afbp;
  SeqAlignPtr          salp;
  SeqIdPtr             sip;
  Boolean              found = FALSE;

  if (sap == NULL || sap->type != 2 || userdata == NULL) 
  {
    return;
  }
  afbp = (AlignmentForBspPtr) userdata;
  if (afbp->bsp == NULL)
  {
    return;
  }
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  for (sip = afbp->bsp->id; sip != NULL && !found; sip = sip->next)
  {
    if (SeqAlignFindSeqId (salp, sip))
    {
      salp = SeqAlignListDup(salp);
      AlnMgr2IndexSeqAlign(salp);
      if (afbp->salp_last == NULL)
      {
        afbp->salp_list = salp; 
      }
      else
      {
        afbp->salp_last->next = salp;
      }
      afbp->salp_last = salp;
      found = TRUE;
    }
  }
}

extern SeqAlignPtr FindAlignmentsForBioseq (BioseqPtr bsp)
{
  SeqEntryPtr         topsep;
  AlignmentForBspData afbd;
  SeqLocPtr           slp;
  SeqIdPtr            sip;
  
  if (bsp == NULL) return NULL;
  topsep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  afbd.salp_list = NULL;
  afbd.salp_last = NULL;
  if (bsp->repr == Seq_repr_seg)
  {
    for (slp = bsp->seq_ext; slp != NULL; slp = slp->next)
    {
      sip = SeqLocId (slp);
      afbd.bsp = BioseqFind (sip);
      VisitAnnotsInSep (topsep, &afbd, FindAlignmentsForBioseqCallback);
    }
  }
  else
  {
    afbd.bsp = bsp;
    VisitAnnotsInSep (topsep, &afbd, FindAlignmentsForBioseqCallback);
  }
  
  return afbd.salp_list;
}

static void FindAlignSeqAnnotsForBioseqCallback (SeqAnnotPtr sap, Pointer userdata)
{
  AlignmentForBspPtr   afbp;
  SeqAlignPtr          salp;
  SeqIdPtr             sip;
  Boolean              found = FALSE;

  if (sap == NULL || sap->type != 2 || userdata == NULL) 
  {
    return;
  }
  afbp = (AlignmentForBspPtr) userdata;
  if (afbp->bsp == NULL)
  {
    return;
  }
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  for (sip = afbp->bsp->id; sip != NULL && !found; sip = sip->next)
  {
    if (SeqAlignFindSeqId (salp, sip))
    {
      ValNodeAddPointer (&(afbp->seq_annot_list), 0, sap);
      found = TRUE;
    }
  }
}

extern ValNodePtr FindAlignSeqAnnotsForBioseq (BioseqPtr bsp)
{
  SeqEntryPtr         topsep;
  AlignmentForBspData afbd;
  SeqLocPtr           slp;
  SeqIdPtr            sip;
  
  if (bsp == NULL) return NULL;
  topsep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  afbd.salp_list = NULL;
  afbd.salp_last = NULL;
  afbd.seq_annot_list = NULL;
  if (bsp->repr == Seq_repr_seg)
  {
    for (slp = bsp->seq_ext; slp != NULL; slp = slp->next)
    {
      sip = SeqLocId (slp);
      afbd.bsp = BioseqFind (sip);
      VisitAnnotsInSep (topsep, &afbd, FindAlignSeqAnnotsForBioseqCallback);
    }
  }
  else
  {
    afbd.bsp = bsp;
    VisitAnnotsInSep (topsep, &afbd, FindAlignSeqAnnotsForBioseqCallback);
  }
  
  return afbd.seq_annot_list;
}

NLM_EXTERN void ChangeSeqIdToWorstID (SeqIdPtr sip)
{
  BioseqPtr       bsp;
  SeqIdPtr        id;
  Pointer         pnt;

  if (sip == NULL)
    return;
  bsp = BioseqFindCore (sip);
  if (bsp == NULL)
    return;
  id = SeqIdDup (SeqIdFindWorst (bsp->id));
  if (id == NULL)
    return;
  /* now remove SeqId contents to reuse SeqId valnode */
  pnt = sip->data.ptrvalue;
  switch (sip->choice) {
  case SEQID_LOCAL:            /* local */
    ObjectIdFree ((ObjectIdPtr) pnt);
    break;
  case SEQID_GIBBSQ:           /* gibbseq */
  case SEQID_GIBBMT:           /* gibbmt */
    break;
  case SEQID_GIIM:             /* giimid */
    GiimFree ((GiimPtr) pnt);
    break;
  case SEQID_GENBANK:          /* genbank */
  case SEQID_EMBL:             /* embl */
  case SEQID_PIR:              /* pir   */
  case SEQID_SWISSPROT:        /* swissprot */
  case SEQID_OTHER:            /* other */
  case SEQID_DDBJ:
  case SEQID_PRF:
  case SEQID_TPG:
  case SEQID_TPE:
  case SEQID_TPD:
  case SEQID_GPIPE:
    TextSeqIdFree ((TextSeqIdPtr) pnt);
    break;
  case SEQID_PATENT:           /* patent seq id */
    PatentSeqIdFree ((PatentSeqIdPtr) pnt);
    break;
  case SEQID_GENERAL:          /* general */
    DbtagFree ((DbtagPtr) pnt);
    break;
  case SEQID_GI:               /* gi */
    break;
  case SEQID_PDB:
    PDBSeqIdFree ((PDBSeqIdPtr) pnt);
    break;
  }
  sip->choice = id->choice;
  sip->data.ptrvalue = id->data.ptrvalue;
  SeqIdStripLocus (sip);
}

NLM_EXTERN void ChangeSeqLocToWorstID (SeqLocPtr slp)
{
  SeqLocPtr       loc;
  PackSeqPntPtr   psp;
  SeqBondPtr      sbp;
  SeqIntPtr       sinp;
  SeqIdPtr        sip;
  SeqPntPtr       spp;

  while (slp != NULL) {
    switch (slp->choice) {
    case SEQLOC_NULL:
      break;
    case SEQLOC_EMPTY:
    case SEQLOC_WHOLE:
      sip = (SeqIdPtr) slp->data.ptrvalue;
      ChangeSeqIdToWorstID (sip);
      break;
    case SEQLOC_INT:
      sinp = (SeqIntPtr) slp->data.ptrvalue;
      if (sinp != NULL) {
        sip = sinp->id;
        ChangeSeqIdToWorstID (sip);
      }
      break;
    case SEQLOC_PNT:
      spp = (SeqPntPtr) slp->data.ptrvalue;
      if (spp != NULL) {
        sip = spp->id;
        ChangeSeqIdToWorstID (sip);
      }
      break;
    case SEQLOC_PACKED_PNT:
      psp = (PackSeqPntPtr) slp->data.ptrvalue;
      if (psp != NULL) {
        sip = psp->id;
        ChangeSeqIdToWorstID (sip);
      }
      break;
    case SEQLOC_PACKED_INT:
    case SEQLOC_MIX:
    case SEQLOC_EQUIV:
      loc = (SeqLocPtr) slp->data.ptrvalue;
      while (loc != NULL) {
        ChangeSeqLocToWorstID (loc);
        loc = loc->next;
      }
      break;
    case SEQLOC_BOND:
      sbp = (SeqBondPtr) slp->data.ptrvalue;
      if (sbp != NULL) {
        spp = (SeqPntPtr) sbp->a;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToWorstID (sip);
        }
        spp = (SeqPntPtr) sbp->b;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToWorstID (sip);
        }
      }
      break;
    case SEQLOC_FEAT:
      break;
    default:
      break;
    }
    slp = slp->next;
  }
}

/* This function will remove DenDiag and pairwise alignments if they contain
 * the sequence identified by sip, otherwise it will remove the sequence from
 * the alignment.
 */
static SeqAlignPtr RemoveOneSequenceFromAlignment (SeqIdPtr sip, SeqAlignPtr salphead)
{
  Uint4       seqid_order;
  SeqIdPtr    tmpsip;
  SeqAlignPtr salp, salp_next, prev_salp, remove_salp, last_remove;
  
  if (!FindSeqIdinSeqAlign (salphead, sip)) return salphead;
  
  salp = salphead;
  prev_salp = NULL;
  remove_salp = NULL;
  last_remove = NULL;
  while (salp != NULL)
  {
    salp_next = salp->next;
    tmpsip = SeqIdPtrFromSeqAlign (salp);
    seqid_order = SeqIdOrderInBioseqIdList(sip, tmpsip);
    if (seqid_order == 0)
    {
      /* do nothing for this subalignment */
      prev_salp = salp;
    }
    else if (salp->dim == 2 || salphead->segtype ==1)
    {
      /* This is for a pairwise alignment or a DENDIAG alignment */
      if (prev_salp == NULL)
      {
          salphead = salp->next;
      }
      else
      {
          prev_salp->next = salp->next;
      }
      /* save the alignments that we want to free in a list and get rid of them
       * at the end - freeing them beforehand causes problems with listing the
       * IDs in the alignment.
       */
      salp->next = NULL;
      if (remove_salp == NULL)
      {
          remove_salp = salp;
      }
      else
      {
          last_remove->next = salp;
      }
      last_remove = salp;
    }
    else 
    {
      SeqAlignBioseqDeleteById (salphead, sip);  
      prev_salp = salp;
    }
    salp = salp_next;
  }
  /* Now we can free the alignment */
  SeqAlignFree (remove_salp);
  return salphead;
}

static void RemoveSequenceFromAlignmentsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  SeqIdPtr    sip;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  sip = (SeqIdPtr) userdata;
  sap->data = RemoveOneSequenceFromAlignment (sip, salp);
  /* if we've deleted all of the alignments, get rid of the annotation as well */
  if (sap->data == NULL)
  {
      sap->idx.deleteme = TRUE;
  }
}

typedef struct checkforremovesequencefromalignments
{
  Boolean  found_problem;
  SeqIdPtr sip;
} CheckForRemoveSequenceFromAlignmentsData, PNTR CheckForRemoveSequenceFromAlignmentsPtr;

/* This is the callback function for looking for pairwise alignments.
 * If we delete the first sequence in a pairwise alignment, we end up deleting
 * the entire alignment because that sequence is paired with every other sequence.
 */
static void CheckForRemoveSequenceFromAlignmentsProblemsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  CheckForRemoveSequenceFromAlignmentsPtr p;
  SeqAlignPtr salphead, salp;
  Uint4       seqid_order;
  SeqIdPtr    tmpsip;
  
  if (sap == NULL || sap->type != 2
      || (p = (CheckForRemoveSequenceFromAlignmentsPtr)userdata) == NULL
      || p->found_problem)
  {
      return;
  }
  salphead = (SeqAlignPtr) sap->data;
  if (salphead == NULL) return;
  
  if (!FindSeqIdinSeqAlign (salphead, p->sip))
  {
      return;
  }
  for (salp = salphead; salp != NULL; salp = salp->next)
  {
    tmpsip = SeqIdPtrFromSeqAlign (salp);
    seqid_order = SeqIdOrderInBioseqIdList(p->sip, tmpsip);
    if (seqid_order == 0)
    {
      continue;
    }
    else if (seqid_order == 1 && salp->dim == 2)
    {
      p->found_problem = TRUE;      
    }
  }
}

extern Boolean IsSequenceFirstInPairwise (SeqEntryPtr sep, SeqIdPtr sip)
{
  CheckForRemoveSequenceFromAlignmentsData data;
  
  if (sep == NULL || sip == NULL)
  {
    return FALSE;
  }
  
    data.sip = sip;
    data.found_problem = FALSE;
  
  VisitAnnotsInSep (sep, (Pointer) &data, CheckForRemoveSequenceFromAlignmentsProblemsCallback);
  return data.found_problem;
}

extern Boolean RemoveSequenceFromAlignments (SeqEntryPtr sep, SeqIdPtr sip)
{
  if (sep == NULL || sip == NULL)
  {
    return FALSE;
  }
  if (IsSequenceFirstInPairwise (sep, sip))
  {
    return FALSE;
  }
  VisitAnnotsInSep (sep, (Pointer) sip, RemoveSequenceFromAlignmentsCallback);
  return TRUE;
}


static CharPtr inferencePrefix [] = {
  "",
  "similar to sequence",
  "similar to AA sequence",
  "similar to DNA sequence",
  "similar to RNA sequence",
  "similar to RNA sequence, mRNA",
  "similar to RNA sequence, EST",
  "similar to RNA sequence, other RNA",
  "profile",
  "nucleotide motif",
  "protein motif",
  "ab initio prediction",
  "alignment",
  NULL
};

static Int2 ValidateInferenceAccession (CharPtr str, Char chr, Boolean fetchAccn, Boolean has_fetch_function)

{
  Int2      accnv, rsult;
  ErrSev    sev;
  SeqIdPtr  sip;
  CharPtr   tmp;

  if (StringHasNoText (str)) return EMPTY_INFERENCE_STRING;

  rsult = VALID_INFERENCE;

  tmp = StringChr (str, chr);
  if (tmp != NULL) {
    *tmp = '\0';
    tmp++;
    TrimSpacesAroundString (str);
    TrimSpacesAroundString (tmp);
    if (StringDoesHaveText (tmp)) {
      if (StringICmp (str, "INSD") == 0 || StringICmp (str, "RefSeq") == 0) {
        accnv = ValidateAccnDotVer (tmp);
        if (accnv == -5 || accnv == -6) {
          rsult = BAD_INFERENCE_ACC_VERSION;
        } else if (accnv != 0) {
          rsult = BAD_INFERENCE_ACCESSION;
        } else if (fetchAccn) {
          sip = SeqIdFromAccessionDotVersion (tmp);
          sev = ErrGetMessageLevel ();
          ErrSetMessageLevel (SEV_ERROR);
          if (has_fetch_function && GetGIForSeqId (sip) == 0) {
            rsult = ACC_VERSION_NOT_PUBLIC;
          }
          ErrSetMessageLevel (sev);
          SeqIdFree (sip);
        }
      }
    }
    if (StringChr (tmp, ' ') != NULL) rsult = SPACES_IN_INFERENCE;
  } else {
    rsult = SINGLE_INFERENCE_FIELD;
  }

  return rsult;
}

NLM_EXTERN Int2 ValidateInferenceQualifier (CharPtr val, Boolean fetchAccn)

{
  Int2           best, j, rsult, tmprsult;
  Char           ch;
  Boolean        has_fetch_function, same_species;
  size_t         len;
  CharPtr        nxt, ptr, rest, str;
  ObjMgrProcPtr  ompp = NULL;

  if (StringHasNoText (val)) return EMPTY_INFERENCE_STRING;

  rest = NULL;
  best = -1;
  for (j = 0; inferencePrefix [j] != NULL; j++) {
    len = StringLen (inferencePrefix [j]);
    if (StringNICmp (val, inferencePrefix [j], len) != 0) continue;
    rest = val + len;
    best = j;
  }

  if (best < 0 || inferencePrefix [best] == NULL) return BAD_INFERENCE_PREFIX;

  if (rest == NULL) return BAD_INFERENCE_BODY;

  same_species = FALSE;
  ch = *rest;
  while (IS_WHITESP (ch)) {
    rest++;
    ch = *rest;
  }
  if (StringNICmp (rest, "(same species)", 14) == 0) {
    same_species = TRUE;
    rest += 14;
  }
  ch = *rest;
  while (IS_WHITESP (ch) || ch == ':') {
    rest++;
    ch = *rest;
  }

  if (StringHasNoText (rest)) return BAD_INFERENCE_BODY;

  rsult = VALID_INFERENCE;
  if (same_species && best > 7) {
    rsult = SAME_SPECIES_MISUSED;
  }

  has_fetch_function = FALSE;
  while ((ompp = ObjMgrProcFindNext(NULL, OMPROC_FETCH, OBJ_SEQID, OBJ_SEQID, ompp)) != NULL) {
    if ((ompp->subinputtype == 0) && (ompp->suboutputtype == SEQID_GI)) {
      has_fetch_function = TRUE;
    }
  }

  str = StringSave (rest);

  if (best >= 1 && best <= 7) {
    tmprsult = ValidateInferenceAccession (str, ':', fetchAccn, has_fetch_function);
    if (tmprsult != VALID_INFERENCE) {
      rsult = tmprsult;
    }
  } else if (best == 12) {
    ptr = StringRChr (str, ':');
    while (ptr != NULL) {
      *ptr = '\0';
      ptr++;
      nxt = StringChr (ptr, ',');
      if (nxt != NULL) {
        *nxt = '\0';
      }
      tmprsult = ValidateInferenceAccession (ptr, '|', fetchAccn, has_fetch_function);
      if (tmprsult != VALID_INFERENCE) {
        rsult = tmprsult;
      }
      ptr = nxt;
    }
  }

  MemFree (str);

  return rsult;
}

extern void MergeFeatureIntervalsToParts (SeqFeatPtr sfp, Boolean ordered)
{
  BioseqPtr     bsp;
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  Boolean       noLeft;
  Boolean       noRight;
  RnaRefPtr     rrp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  tRNAPtr       trna;
  
  if (sfp == NULL || sfp->location == NULL) return;
  sip = SeqLocId (sfp->location);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return;
  if (ISA_aa (bsp->mol)) return;
  CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
  slp = SegLocToPartsEx (bsp, sfp->location, ordered);
  if (slp == NULL) return;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = slp;
  FreeAllFuzz (sfp->location);
  SetSeqLocPartial (sfp->location, noLeft, noRight);
  sfp->partial = (sfp->partial || noLeft || noRight);
  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL && crp->code_break != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          slp = SegLocToPartsEx (bsp, cbp->loc, ordered);
          if (slp != NULL) {
            cbp->loc = SeqLocFree (cbp->loc);
            cbp->loc = slp;
            FreeAllFuzz (cbp->loc);
          }
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 3 && rrp->ext.choice == 2) {
        trna = rrp->ext.value.ptrvalue;
        if (trna != NULL && trna->anticodon != NULL) {
          slp = SegLocToPartsEx (bsp, trna->anticodon, ordered);
          if (slp != NULL) {
            trna->anticodon = SeqLocFree (trna->anticodon);
            trna->anticodon = slp;
            FreeAllFuzz (trna->anticodon);
          }
        }
      }
      break;
    default :
      break;
  }  
}

extern void ExtendSingleGeneOnMRNA (BioseqPtr bsp, Pointer userdata)
{
  MolInfoPtr        mip;
  SeqDescrPtr       sdp;
  Boolean           is_mrna = FALSE, is_master_seq = FALSE;
  SeqFeatPtr        gene = NULL;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  Int4              num_cds = 0;
  SeqLocPtr         slp;
  Boolean           partial5, partial3;
  BioSourcePtr      biop;
  OrgRefPtr         orp;
  BioseqSetPtr      bssp;

  if (bsp == NULL || bsp->length == 0
      || !ISA_na (bsp->mol)) {
    return;
  }

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->biomol == MOLECULE_TYPE_MRNA) {
      is_mrna = TRUE;
    }
  }
  if (!is_mrna) {
    return;
  }

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      if (biop->origin == ORG_ARTIFICIAL) {
        orp = biop->org;
        if (orp != NULL) {
          if (StringICmp (orp->taxname, "synthetic construct") == 0) return;
        }
      }
    }
  }

  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset) {
      is_master_seq = TRUE;
    }
  }
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context)) {
    if (sfp->data.choice == SEQFEAT_GENE) {
      /* skip this sequence if it has more than one gene */
      if (gene == NULL) {
        gene = sfp;
      } else {
        return; 
      }
    } else if (sfp->data.choice == SEQFEAT_CDREGION) {
      num_cds++;
      /* skip this seuqence if it has more than one coding region */
      if (num_cds > 1 && !is_master_seq) {
        return;
      }
    }
  }
  if (gene != NULL && BioseqFindFromSeqLoc (gene->location) == bsp) {
    CheckSeqLocForPartial (gene->location, &partial5, &partial3);
    /* gene should cover entire length of sequence */
    slp = SeqLocIntNew (0, bsp->length - 1, SeqLocStrand (gene->location), SeqIdFindBest (bsp->id, 0));
    SetSeqLocPartial (slp, partial5, partial3);
    gene->location = SeqLocFree (gene->location);
    gene->location = slp;
    if (is_master_seq) {
      MergeFeatureIntervalsToParts (gene, FALSE);
    }
  }
}

/* Functions for the Discrepancy Report */


static Boolean IsmRNASequenceInGenProdSet (BioseqPtr bsp)
{
  SeqMgrDescContext dcontext;
  BioseqSetPtr      bssp;
  SeqDescrPtr sdp;
  MolInfoPtr  mip;
  Boolean rval = FALSE;

  if (bsp != NULL && bsp->mol == Seq_mol_rna && bsp->idx.parentptr != NULL && bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
      if (sdp != NULL && sdp->data.ptrvalue != NULL && sdp->choice == Seq_descr_molinfo) {
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip->biomol == MOLECULE_TYPE_MRNA) {
          rval = TRUE;
        }
      }
    } else if (bssp->_class == BioseqseqSet_class_nuc_prot && bssp->idx.parentptr != NULL && bssp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bssp->idx.parentptr;
      if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
        sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
        if (sdp != NULL && sdp->data.ptrvalue != NULL && sdp->choice == Seq_descr_molinfo) {
          mip = (MolInfoPtr) sdp->data.ptrvalue;
          if (mip->biomol == MOLECULE_TYPE_MRNA) {
            rval = TRUE;
          }
        }
      }
    }
  }
  return rval;
}


typedef struct skipmrnafeaturesingenprodset {
  Pointer userdata;
  VisitFeaturesFunc callback;
} SkipmRNAFeaturesInGenProdSetData, PNTR SkipmRNAFeaturesInGenProdSetPtr;

static void VisitGenProdSetFeaturesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  SkipmRNAFeaturesInGenProdSetPtr p;

  if (sfp == NULL || userdata == NULL) {
    return;
  }

  p = (SkipmRNAFeaturesInGenProdSetPtr) userdata;
  if (p->callback == NULL) {
    return;
  }

  if (!IsmRNASequenceInGenProdSet(BioseqFindFromSeqLoc (sfp->location))) {
    (p->callback) (sfp, p->userdata);
  }
}


extern void VisitGenProdSetFeatures (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback)
{
  SkipmRNAFeaturesInGenProdSetData d;

  d.callback = callback;
  d.userdata = userdata;
  VisitFeaturesInSep (sep, &d, VisitGenProdSetFeaturesCallback);
}


extern ClickableItemPtr 
NewClickableItem 
(Uint4           clickable_item_type,
 CharPtr         description_fmt,
 ValNodePtr      item_list)
{
  ClickableItemPtr dip;
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = clickable_item_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (description_fmt) + 15));
    sprintf (dip->description, description_fmt, ValNodeLen (item_list));
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = item_list;
    dip->subcategories = NULL;
    dip->expanded = FALSE;
    dip->level = 0;
  }
  return dip;  
}


extern ClickableItemPtr ClickableItemFree (ClickableItemPtr cip)
{
  if (cip != NULL)
  {
    cip->description = MemFree (cip->description);
    if (cip->datafree_func != NULL)
    {
      (cip->datafree_func) (cip->callback_data);
    }
    cip->item_list = ValNodeFree (cip->item_list);
  
    cip->subcategories = FreeClickableList (cip->subcategories);
    cip = MemFree (cip);
  }
  return cip;
}


extern ValNodePtr FreeClickableList (ValNodePtr list)
{
  ValNodePtr       list_next;
  
  while (list != NULL) {
    list_next = list->next;
    list->next = NULL;
    list->data.ptrvalue = ClickableItemFree (list->data.ptrvalue);
    list = ValNodeFree (list);
    list = list_next;
  }
  return list;
}

/* utility functions for the discrepancy report tests */
static void ValNodeLinkCopy (ValNodePtr PNTR list1, ValNodePtr list2)
{
  if (list1 == NULL) return;
  while (list2 != NULL)
  {
    ValNodeAddPointer (list1, list2->choice, list2->data.ptrvalue);
    list2 = list2->next;
  }
}


static Boolean ValNodeStringListMatch (ValNodePtr vnp1, ValNodePtr vnp2)
{
  if (vnp1 == NULL && vnp2 == NULL)
  {
    return TRUE;
  }
  else if (vnp1 == NULL || vnp2 == NULL)
  {
    return FALSE;
  }
  else if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0)
  {
    return FALSE;
  }
  else
  {
    return ValNodeStringListMatch (vnp1->next, vnp2->next);
  }
}


static ValNodePtr ItemListFromSubcategories (ValNodePtr subcategories)
{
  ValNodePtr       vnp;
  ClickableItemPtr cip;
  ValNodePtr       item_list = NULL;

  for (vnp = subcategories; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      ValNodeLinkCopy (&item_list, cip->item_list);
    }
  }
  return item_list;
}


extern GlobalDiscrepancyPtr 
GlobalDiscrepancyNew (CharPtr str, Uint1 data_choice, Pointer data)
{
  GlobalDiscrepancyPtr g;

  g = (GlobalDiscrepancyPtr) MemNew (sizeof (GlobalDiscrepancyData));
  g->str = StringSave (str);
  g->data_choice = data_choice;
  if (g->data_choice == 0) {
    g->data = StringSave (data);
  } else {
    g->data = data;
  }
  return g;
}


extern GlobalDiscrepancyPtr GlobalDiscrepancyFree (GlobalDiscrepancyPtr g)
{
  if (g != NULL) {
    g->str = MemFree (g->str);
    if (g->data_choice == 0) {
      g->data = MemFree (g->data);
    }
    g = MemFree (g);
  }
  return g;
}


extern ValNodePtr FreeGlobalDiscrepancyList (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = GlobalDiscrepancyFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


extern void ConvertGlobalDiscrepancyToText (GlobalDiscrepancyPtr g, Boolean use_feature_fmt)
{
  ValNode vn;
  ValNodePtr list_copy;

  if (g == NULL || g->data_choice == 0) return;
  
  vn.choice = g->data_choice;
  vn.data.ptrvalue = g->data;
  vn.next = NULL;

  g->data_choice = 0;
  if (use_feature_fmt) {  
    list_copy = ReplaceDiscrepancyItemWithFeatureTableStrings (&vn);
    g->data = list_copy->data.ptrvalue;
    list_copy = ValNodeFree (list_copy);
  } else {
    g->data = GetDiscrepancyItemText (&vn);
  }
}


extern void ConvertGlobalDiscrepancyListToText (ValNodePtr vnp, Boolean use_feature_fmt)
{
  while (vnp != NULL) {
    ConvertGlobalDiscrepancyToText (vnp->data.ptrvalue, use_feature_fmt);
    vnp = vnp->next;
  }
}


extern ValNodePtr GetGlobalDiscrepancyItem (GlobalDiscrepancyPtr g)
{
  ValNodePtr rval = NULL;
  if (g != NULL) {
    rval = ValNodeNew (NULL);
    rval->choice = g->data_choice;
    if (rval->choice == 0) {
      rval->data.ptrvalue = StringSave (g->data);
    } else {
      rval->data.ptrvalue = g->data;
    }
  }
  return rval;
}


extern CharPtr GetGlobalDiscrepancyStr (GlobalDiscrepancyPtr g)
{
  CharPtr rval = NULL;
  if (g != NULL) {
    rval = g->str;
  }
  return rval;
}


NLM_EXTERN int LIBCALLBACK SortVnpByGlobalDiscrepancyString (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  GlobalDiscrepancyPtr g1, g2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      g1 = (GlobalDiscrepancyPtr) vnp1->data.ptrvalue;
      g2 = (GlobalDiscrepancyPtr) vnp2->data.ptrvalue;
      if (g1 != NULL && g2 != NULL && g1->str != NULL && g2->str != NULL) {
        return StringICmp (g1->str, g2->str);
      }
    }
  }
  return 0;
}


static Int4 CountDupGlobalDiscrepancy (ValNodePtr vnp)
{
  GlobalDiscrepancyPtr g1, g2;
  Int4                 num_dup = 1;

  if (vnp == NULL 
      || (g1 = (GlobalDiscrepancyPtr) vnp->data.ptrvalue) == NULL
      || StringHasNoText (g1->str)) {
    return 0;
  } else if (vnp->next == NULL) {
    return 1;
  }
  vnp = vnp->next;
  while (vnp != NULL
         && (g2 = (GlobalDiscrepancyPtr) vnp->data.ptrvalue) != NULL 
         && StringICmp (g1->str, g2->str) == 0) {
    num_dup++;
    vnp = vnp->next;
  }
  return num_dup;
}


extern ClickableItemPtr
ReportNonUniqueGlobalDiscrepancy 
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 CharPtr    ind_cat_fmt,
 Uint4      clickable_item_type,
 Boolean    keep_top_category)

{
  Boolean       print_heading = TRUE;
  Int4          num_dup, total_dup = 0, i;
  ValNodePtr       item_list;
  ClickableItemPtr cip = NULL;
  ValNodePtr       subcategories = NULL;
  CharPtr          str;

  while (vnp != NULL) {
    num_dup = CountDupGlobalDiscrepancy (vnp);
    if (num_dup > 1) {
      total_dup += num_dup;
      str = GetGlobalDiscrepancyStr (vnp->data.ptrvalue);
      if (str == NULL) str = "";
      item_list = NULL;
      i = num_dup;
      while (i > 0) {
        ValNodeLink (&item_list, GetGlobalDiscrepancyItem (vnp->data.ptrvalue));
        vnp = vnp->next;
        i--;
      }
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (ind_cat_fmt) + StringLen (str) + 15));
      sprintf (cip->description, ind_cat_fmt, num_dup, str);
      cip->clickable_item_type = clickable_item_type;
      cip->item_list = item_list;
      ValNodeAddPointer (&subcategories, 0, cip);
    } else {
      vnp = vnp->next;
    }
  }
  if (subcategories != NULL) {
    if (subcategories->next == NULL && !keep_top_category) {
      cip = subcategories->data.ptrvalue;
      subcategories = ValNodeFree (subcategories);
    } else {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label_fmt) + 15));
      sprintf (cip->description, label_fmt, total_dup);
      cip->clickable_item_type = clickable_item_type;
      cip->subcategories = subcategories;
    }
  }
  return cip;
}


static Boolean IsLocusTagFormatBad (CharPtr locus_tag)
{
  CharPtr cp;
  Boolean after_underscore = FALSE;
  
  if (StringHasNoText (locus_tag))
  {
    return FALSE;
  }
  
  cp = locus_tag;
  if (!isalpha (*cp))
  {
    return TRUE;
  }
  cp++;
  while (*cp != 0)
  {
    if (*cp == '_')
    {
      if (after_underscore)
      {
        return TRUE;
      }
      else
      {
        after_underscore = TRUE;
        if (*(cp + 1) == 0)
        {
          return TRUE;
        }
      }
    }
    else if (!isalpha (*cp) && !isdigit (*cp))
    {
      return TRUE;
    }
    cp++;
  }
  if (after_underscore)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


extern ClickableItemPtr ReportBadLocusTagFormat (ValNodePtr list)
{
  ValNodePtr vnp, item_list = NULL;
  ClickableItemPtr cip = NULL;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    if (IsLocusTagFormatBad (GetGlobalDiscrepancyStr (vnp->data.ptrvalue))) {
      ValNodeLink (&item_list, GetGlobalDiscrepancyItem (vnp->data.ptrvalue));
    }
  }
  if (item_list != NULL) {
    cip = NewClickableItem (DISC_GENE_LOCUS_TAG_BAD_FORMAT, "%d locus tags are incorrectly formatted.", item_list);
  }
  return cip;
}


static CharPtr GetGlobalDiscrepancyPrefix (GlobalDiscrepancyPtr g)
{
  CharPtr cp, prefix = NULL;
  Int4    len;

  if (g == NULL) return NULL;
  cp = StringChr (g->str, '_');
  if (cp != NULL) {
    len = cp - g->str;
    prefix = MemNew (sizeof (Char) * (len + 1));
    StringNCpy (prefix, g->str, len);
    prefix[len] = 0;
  }
  return prefix;
}


static Int4 CountDupGlobalDiscrepancyPrefix (ValNodePtr vnp)
{
  GlobalDiscrepancyPtr g1, g2;
  CharPtr              cp;
  Int4                 len;
  Int4                 num_dup = 1;

  if (vnp == NULL 
      || (g1 = (GlobalDiscrepancyPtr) vnp->data.ptrvalue) == NULL
      || StringHasNoText (g1->str)
      || (cp = StringChr (g1->str, '_')) == NULL) {
    return 0;
  } else if (vnp->next == NULL) {
    return 1;
  }
  len = cp - g1->str + 1;
  vnp = vnp->next;
  while (vnp != NULL
         && (g2 = (GlobalDiscrepancyPtr) vnp->data.ptrvalue) != NULL 
         && StringNCmp (g1->str, g2->str, len) == 0) {
    num_dup++;
    vnp = vnp->next;
  }
  return num_dup;
}


extern ValNodePtr ReportInconsistentGlobalDiscrepancyPrefixes
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 Uint4      clickable_item_type)

{
  Boolean       print_heading = TRUE;
  Int4          num_dup;
  CharPtr       prefix;
  ValNodePtr    disc_list = NULL;
  ClickableItemPtr cip;

  if (vnp == NULL) return NULL;

  num_dup = CountDupGlobalDiscrepancyPrefix (vnp);
  if (num_dup < ValNodeLen (vnp)) {
    while (vnp != NULL) {
      prefix = GetGlobalDiscrepancyPrefix (vnp->data.ptrvalue);
      num_dup = CountDupGlobalDiscrepancyPrefix (vnp);
      if (num_dup < 1) {
        vnp = vnp->next;
      } else if (prefix != NULL) {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        cip->clickable_item_type = clickable_item_type;
        cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label_fmt) + StringLen (prefix) + 15));
        sprintf (cip->description, label_fmt, num_dup, prefix);
        /* skip duplicates without printing */
        while (num_dup > 0) {
          ValNodeLink (&cip->item_list, GetGlobalDiscrepancyItem (vnp->data.ptrvalue));
          vnp = vnp->next;
          num_dup--;
        }
        prefix = MemFree (prefix);
        ValNodeAddPointer (&disc_list, 0, cip);
      } else {
        /* skip items without prefix */
        while (num_dup > 0) {
          vnp = vnp->next;
          num_dup--;
        }
      }
    }
  } 
  return disc_list;
}


extern ValNodePtr ReportInconsistentGlobalDiscrepancyStrings
(ValNodePtr vnp, 
 CharPtr    label_fmt,
 Uint4      clickable_item_type)

{
  Boolean       print_heading = TRUE;
  Int4          num_dup;
  CharPtr       prefix;
  ValNodePtr    disc_list = NULL;
  ClickableItemPtr cip;

  if (vnp == NULL) return NULL;

  num_dup = CountDupGlobalDiscrepancy (vnp);
  if (num_dup < ValNodeLen (vnp)) {
    while (vnp != NULL) {
      prefix = GetGlobalDiscrepancyStr (vnp->data.ptrvalue);
      num_dup = CountDupGlobalDiscrepancy (vnp);
      if (num_dup < 1) {
        vnp = vnp->next;
      } else if (prefix != NULL) {
        cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
        cip->clickable_item_type = clickable_item_type;
        cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label_fmt) + StringLen (prefix) + 15));
        sprintf (cip->description, label_fmt, num_dup, prefix);
        /* skip duplicates without printing */
        while (num_dup > 0) {
          ValNodeLink (&cip->item_list, GetGlobalDiscrepancyItem (vnp->data.ptrvalue));
          vnp = vnp->next;
          num_dup--;
        }
        ValNodeAddPointer (&disc_list, 0, cip);
      } else {
        /* skip items without prefix */
        while (num_dup > 0) {
          vnp = vnp->next;
          num_dup--;
        }
      }
    }
  } 
  return disc_list;
}


extern ClickableItemPtr ReportMissingFields (ValNodePtr list, CharPtr label_fmt, Uint4 clickable_item_type)
{
  ClickableItemPtr cip;

  if (list == NULL) return NULL;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (label_fmt) + 15));
  sprintf (cip->description, label_fmt, ValNodeLen (list));
  cip->clickable_item_type = clickable_item_type;
  while (list != NULL) {
    ValNodeLink (&(cip->item_list), GetGlobalDiscrepancyItem (list->data.ptrvalue));
    list = list->next;
  }
  return cip;  
}





/* declarations for discrepancy tests */
extern void AddMissingAndSuperfluousGeneDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddDiscrepanciesForNonGeneLocusTags (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindMissingProteinIDs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindCDSmRNAGeneLocationDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindCDSGeneProductConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindDuplicateGeneLocus (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddECNumberNoteDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindPseudoDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddJoinedFeatureDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddOverlappingGeneDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddOverlappingCodingRegionDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddContainedCodingRegionDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void AddRNACDSOverlapDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindShortContigs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindNonmatchingContigSources (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindSuspectProductNames (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindSuspectPhrases (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindInconsistentSourceAndDefline (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindParticalCDSsInCompleteSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindUnknownProteinsWithECNumbers (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindShortSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void tRNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void tRNAFindBadLength (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void rRNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindRNAsWithoutProducts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindTranslExceptNotes (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindCDSOverlappingtRNAs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern int LIBCALLBACK SortVnpByClickableItemDescription (VoidPtr ptr1, VoidPtr ptr2);
extern void CountProteins (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);
extern void FindFeaturesOverlappingSrcFeatures (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list);

/* functions for the missing and superfluous gene tests */
extern Boolean GeneRefMatch (GeneRefPtr grp1, GeneRefPtr grp2)
{
  if (grp1 == NULL && grp2 == NULL)
  {
    return TRUE;
  }
  else if (grp1 == NULL || grp2 == NULL)
  {
    return FALSE;
  }
  else if (StringCmp (grp1->locus, grp2->locus) != 0
           || StringCmp (grp1->allele, grp2->allele) != 0
           || StringCmp (grp1->desc, grp2->desc) != 0
           || StringCmp (grp1->maploc, grp2->maploc) != 0
           || StringCmp (grp1->locus_tag, grp2->locus_tag) != 0
           || (grp1->pseudo && !grp2->pseudo)
           || (!grp1->pseudo && grp2->pseudo)
           || !ValNodeStringListMatch (grp1->db, grp2->db)
           || !ValNodeStringListMatch (grp1->syn, grp2->syn))
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static void ExtractGeneFromListByGeneRef (ValNodePtr PNTR list, GeneRefPtr grp)
{
  ValNodePtr prev = NULL, this_vnp, next_vnp;
  SeqFeatPtr gene_feat;
  
  if (list == NULL || grp == NULL)
  {
    return;
  }
  
  this_vnp = *list;
  while (this_vnp != NULL)
  {
    next_vnp = this_vnp->next;
    gene_feat = (SeqFeatPtr) this_vnp->data.ptrvalue;
    if (gene_feat != NULL && GeneRefMatch (gene_feat->data.value.ptrvalue, grp))
    {
      if (prev == NULL)
      {
        *list = next_vnp;
      }
      else
      {
        prev->next = next_vnp;
      }
      this_vnp->next = NULL;
      ValNodeFree (this_vnp);
    }
    else
    {
      prev = this_vnp;
    }
    this_vnp = next_vnp;    
  }
}


static void ExtractGeneFromListByGene (ValNodePtr PNTR list, SeqFeatPtr gene)
{
  ValNodePtr prev = NULL, this_vnp, next_vnp;
  
  if (list == NULL || gene == NULL)
  {
    return;
  }
  
  this_vnp = *list;
  while (this_vnp != NULL)
  {
    next_vnp = this_vnp->next;
    if (this_vnp->data.ptrvalue == gene)
    {
      if (prev == NULL)
      {
        *list = next_vnp;
      }
      else
      {
        prev->next = next_vnp;
      }
      this_vnp->next = NULL;
      ValNodeFree (this_vnp);
    }
    else
    {
      prev = this_vnp;
    }
    this_vnp = next_vnp;    
  }
}


static void 
CheckGenesForFeatureType 
(ValNodePtr PNTR features_without_genes,
 ValNodePtr PNTR superfluous_genes,
 BioseqPtr  bsp,
 Uint1      feature_type,
 Uint1      feature_subtype,
 Boolean    makes_gene_not_superfluous)
{
  SeqFeatPtr         sfp, gene_sfp;
  GeneRefPtr         grp;
  SeqMgrFeatContext  context;
  
  if (features_without_genes == NULL
      || superfluous_genes == NULL
      || bsp == NULL)
  {
    return;
  }
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, feature_type, feature_subtype, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, feature_type, feature_subtype, &context))
  {
    /* check for gene xref */
    grp = SeqMgrGetGeneXref (sfp);
    if (grp != NULL)
    {
      if (SeqMgrGeneIsSuppressed (grp))
      {
        ValNodeAddPointer (features_without_genes, OBJ_SEQFEAT, sfp);
      }
      else
      {
        ExtractGeneFromListByGeneRef (superfluous_genes, grp);
      }
    }
    else
    {
      gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
      if (gene_sfp == NULL)
      {
        ValNodeAddPointer (features_without_genes, OBJ_SEQFEAT, sfp);
      }
      else if (makes_gene_not_superfluous)
      {
        ExtractGeneFromListByGene (superfluous_genes, gene_sfp);
      }
    }  
  }  
}

typedef struct misssupergenes
{
  ValNodePtr missing_list;
  ValNodePtr super_list;
} MissSuperGenesData, PNTR MissSuperGenesPtr;


static void FindMissingGenes (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;
  ValNodePtr         features_without_genes = NULL;
  ValNodePtr         superfluous_genes = NULL;
  MissSuperGenesPtr  msgp;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  msgp = (MissSuperGenesPtr) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &context))
  {
    ValNodeAddPointer (&superfluous_genes, OBJ_SEQFEAT, sfp);
  }
  
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_CDREGION, 0, TRUE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_RNA, 0, TRUE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_IMP, FEATDEF_RBS, FALSE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_IMP, FEATDEF_exon, FALSE);
  CheckGenesForFeatureType (&features_without_genes, &superfluous_genes, bsp, 
                            SEQFEAT_IMP, FEATDEF_intron, FALSE);

  ValNodeLink (&(msgp->missing_list), features_without_genes);
  ValNodeLink (&(msgp->super_list), superfluous_genes);
}


static void 
GetPseudoAndNonPseudoGeneList 
(ValNodePtr      super_list,
 ValNodePtr PNTR pseudo_list, 
 ValNodePtr PNTR non_pseudo_list)
{
  ValNodePtr vnp;
  SeqFeatPtr gene;
  GeneRefPtr grp;
  
  if (pseudo_list == NULL || non_pseudo_list == NULL)
  {
    return;
  }
  *pseudo_list = NULL;
  *non_pseudo_list = NULL;
  
  for (vnp = super_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_SEQFEAT)
    {
      gene = (SeqFeatPtr) vnp->data.ptrvalue;
      if (gene != NULL && gene->data.choice == SEQFEAT_GENE)
      {
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
        if (gene->pseudo || (grp != NULL && grp->pseudo))
        {
          ValNodeAddPointer (pseudo_list, OBJ_SEQFEAT, gene);
        }
        else
        {
          ValNodeAddPointer (non_pseudo_list, OBJ_SEQFEAT, gene);
        }
      }
    }
  }
}


static void 
GetFrameshiftAndNonFrameshiftGeneList 
(ValNodePtr      super_list,
 ValNodePtr PNTR frameshift_list, 
 ValNodePtr PNTR non_frameshift_list)
{
  ValNodePtr vnp;
  SeqFeatPtr gene;
  
  if (frameshift_list == NULL || non_frameshift_list == NULL)
  {
    return;
  }
  *frameshift_list = NULL;
  *non_frameshift_list = NULL;
  
  for (vnp = super_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice == OBJ_SEQFEAT)
    {
      gene = (SeqFeatPtr) vnp->data.ptrvalue;
      if (gene != NULL 
          && (StringISearch (gene->comment, "frameshift") != NULL
              || StringISearch (gene->comment, "frame shift") != NULL)) 
      {
        ValNodeAddPointer (frameshift_list, OBJ_SEQFEAT, gene);
      }
      else
      {
        ValNodeAddPointer (non_frameshift_list, OBJ_SEQFEAT, gene);
      }
    }
  }
}


extern void AddMissingAndSuperfluousGeneDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip, pseudo_dip, non_pseudo_dip;
  CharPtr            missing_genes_fmt = "%d features have no genes.";
  CharPtr            extra_genes_fmt = "%d gene features are not associated with a CDS or RNA feature.";
  CharPtr            pseudo_extra_genes_fmt = "%d pseudo gene features are not associated with a CDS or RNA feature.";
  CharPtr            non_pseudo_frameshift_extra_genes_fmt = "%d non-pseudo gene features are not associated with a CDS or RNA feature and have frameshift in the comment.";
  CharPtr            non_pseudo_non_frameshift_extra_genes_fmt = "%d non-pseudo gene features are not associated with a CDS or RNA feature and do not have frameshift in the comment.";
  MissSuperGenesData msgd;
  ValNodePtr         non_pseudo_list = NULL, pseudo_list = NULL, vnp;
  ValNodePtr         non_frameshift_list = NULL, frameshift_list = NULL;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  msgd.missing_list = NULL;
  msgd.super_list = NULL;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &msgd, FindMissingGenes);
  }
  
  if (msgd.missing_list != NULL)
  {
    dip = NewClickableItem (DISC_GENE_MISSING, missing_genes_fmt, msgd.missing_list);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
  if (msgd.super_list != NULL)
  {
    GetPseudoAndNonPseudoGeneList (msgd.super_list, &pseudo_list, &non_pseudo_list);
    GetFrameshiftAndNonFrameshiftGeneList (non_pseudo_list, &frameshift_list, &non_frameshift_list);
    non_pseudo_list = ValNodeFree (non_pseudo_list);
    dip = NewClickableItem (DISC_SUPERFLUOUS_GENE, extra_genes_fmt, msgd.super_list);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);     
 
      if (frameshift_list != NULL) 
      {
        non_pseudo_dip = NewClickableItem (DISC_SUPERFLUOUS_GENE, non_pseudo_frameshift_extra_genes_fmt, frameshift_list);
        non_pseudo_dip->level = 1;
        ValNodeAddPointer (&(dip->subcategories), 0, non_pseudo_dip);
      }
      if (non_frameshift_list != NULL) 
      {
        non_pseudo_dip = NewClickableItem (DISC_SUPERFLUOUS_GENE, non_pseudo_non_frameshift_extra_genes_fmt, non_frameshift_list);
        non_pseudo_dip->level = 1;
        ValNodeAddPointer (&(dip->subcategories), 0, non_pseudo_dip);
      }
      if (pseudo_list != NULL) 
      {
        pseudo_dip = NewClickableItem (DISC_SUPERFLUOUS_GENE, pseudo_extra_genes_fmt, pseudo_list);
        pseudo_dip->level = 1;
        ValNodeAddPointer (&(dip->subcategories), 0, pseudo_dip);
      }
    }
  }  
}


/* test for missing or inconsistent protein IDs */
typedef struct prefixcheck 
{
  CharPtr prefix;
  ValNodePtr feature_list;
} PrefixCheckData, PNTR PrefixCheckPtr;


static ValNodePtr FreePrefixCheckList (ValNodePtr prefix_list)
{
  PrefixCheckPtr pcp;
  
  if (prefix_list == NULL)
  {
    return NULL;
  }
  
  prefix_list->next = FreePrefixCheckList (prefix_list->next);
  
  pcp = (PrefixCheckPtr) prefix_list->data.ptrvalue;
  if (pcp != NULL)
  {
    pcp->prefix = MemFree (pcp->prefix);
    pcp->feature_list = ValNodeFree (pcp->feature_list);
    pcp = MemFree (pcp);
  }
  prefix_list = ValNodeFree (prefix_list);
  return NULL;
}


static ClickableItemPtr InconsistentPrefix (PrefixCheckPtr pcp, CharPtr bad_fmt, DiscrepancyType disc_type)
{
  ClickableItemPtr dip = NULL;

  if (pcp == NULL || StringHasNoText (pcp->prefix) || pcp->feature_list == NULL)
  {
    return NULL;
  }
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = disc_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (pcp->prefix)+ 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (pcp->feature_list), pcp->prefix);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = pcp->feature_list;
    pcp->feature_list = NULL;
  }      
  return dip;
}


extern CharPtr discReportInconsistentLocusTagPrefixFmt = "%d features have locus tag prefix %s.";
extern CharPtr discReportInconsistentProteinIDPrefixFmt = "%d sequences have protein ID prefix %s.";
extern CharPtr discReportBadProteinIdFmt = "%d proteins have invalid IDs.";

static ClickableItemPtr InconsistentLocusTagPrefix (PrefixCheckPtr pcp)
{
  return InconsistentPrefix (pcp, discReportInconsistentLocusTagPrefixFmt, DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX);
}


extern void FindProteinIDCallback (BioseqPtr bsp, Pointer userdata)
{
  ProtIdListsPtr pip;
  SeqIdPtr       sip;
  DbtagPtr       dbt = NULL;

  if (bsp == NULL || ! ISA_aa (bsp->mol) || userdata == NULL)
  {
    return;
  }

  pip = (ProtIdListsPtr) userdata;

  for (sip = bsp->id; sip != NULL && dbt == NULL; sip = sip->next)
  {
    if (sip->choice == SEQID_GENERAL)
    {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (IsSkippableDbtag(dbt))
      {
        dbt = NULL;
      }
    }
  }
  if (dbt == NULL)
  {
    ValNodeAddPointer (&(pip->missing_gnl_list), 0, GlobalDiscrepancyNew (NULL, OBJ_BIOSEQ, bsp));
  }
  else
  {  
    ValNodeAddPointer (&(pip->gnl_list), 0, GlobalDiscrepancyNew (dbt->db, OBJ_BIOSEQ, bsp));
  }
}



extern void FindMissingProteinIDs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  ProtIdListsData  pid;
  ValNodePtr       vnp;
  
  if (discrepancy_list == NULL) return;
  
  MemSet (&pid, 0, sizeof (ProtIdListsData));
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &pid, FindProteinIDCallback);
  }
  
  if (pid.missing_gnl_list != NULL)
  {
    dip = ReportMissingFields (pid.missing_gnl_list, discReportBadProteinIdFmt, DISC_MISSING_PROTEIN_ID);
    if (dip != NULL) {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
    pid.missing_gnl_list = FreeGlobalDiscrepancyList (pid.missing_gnl_list);
  }
  if (pid.gnl_list != NULL)
  {
    pid.gnl_list = ValNodeSort (pid.gnl_list, SortVnpByGlobalDiscrepancyString);
    ValNodeLink (discrepancy_list, 
                 ReportInconsistentGlobalDiscrepancyStrings (pid.gnl_list,
                                                             discReportInconsistentProteinIDPrefixFmt,
                                                             DISC_INCONSISTENT_PROTEIN_ID_PREFIX));
    pid.gnl_list = FreeGlobalDiscrepancyList (pid.gnl_list);
  } 
}


typedef struct locustagcheck
{
  ValNodePtr locus_tags_list;
  ValNodePtr missing_list;
} LocusTagCheckData, PNTR LocusTagCheckPtr;

static void GeneLocusTagDiscrepancyCallback (ValNodePtr item_list, Pointer userdata)
{
  Message (MSG_OK, "I could launch the editor for the individual gene...");
}

static void CheckGeneLocusTag (SeqFeatPtr sfp, Pointer userdata)
{
  GeneRefPtr         grp;
  LocusTagCheckPtr   ltcp;
  
  if (sfp == NULL || userdata == NULL || sfp->data.choice != SEQFEAT_GENE || sfp->data.value.ptrvalue == NULL)
  {
    return;
  }
  
  ltcp = (LocusTagCheckPtr) userdata;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp != NULL) {
    if (grp->pseudo) return;
    if (StringDoesHaveText (grp->locus_tag)) {
      ValNodeAddPointer (&(ltcp->locus_tags_list), 0, 
                          GlobalDiscrepancyNew (grp->locus_tag, OBJ_SEQFEAT, sfp));
    } else {
      ValNodeAddPointer (&(ltcp->missing_list), 0,
                          GlobalDiscrepancyNew (NULL, OBJ_SEQFEAT, sfp));
    }
  }

}

static Boolean AlreadyInList (ValNodePtr vnp, SeqFeatPtr sfp)
{
  while (vnp != NULL && vnp->data.ptrvalue != sfp)
  {
    vnp = vnp->next;
  }
  if (vnp == NULL)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static SeqFeatPtr GetNextGene (SeqFeatPtr sfp)
{
  BioseqPtr bsp;
  SeqFeatPtr sfp_next;
  SeqMgrFeatContext fcontext;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return NULL;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return NULL;
  /* initialize fcontext for search */
  sfp_next = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &fcontext);
  while (sfp_next != sfp && sfp_next != NULL)
  {
    sfp_next = SeqMgrGetNextFeature (bsp, sfp_next, SEQFEAT_GENE, FEATDEF_GENE, &fcontext);
  }
  if (sfp_next != sfp) return NULL;

  sfp_next = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &fcontext);
  return sfp_next;
}


static ValNodePtr FindValNodeForGlobalDiscrepancyFeature (ValNodePtr start_search, Int4 len_search, SeqFeatPtr sfp)
{
  GlobalDiscrepancyPtr g;

  while (start_search != NULL && len_search > 0) {
    g = (GlobalDiscrepancyPtr) start_search->data.ptrvalue;
    if (g != NULL && g->data_choice == OBJ_SEQFEAT && g->data == sfp) {
      return start_search;
    } else {
      start_search = start_search->next;
      len_search--;
    }
  }
  return NULL;
}


static ValNodePtr FindAdjacentGenesInSubList (ValNodePtr sub_list, Int4 list_len)
{
  GlobalDiscrepancyPtr g;
  SeqFeatPtr           sfp, sfp_next;
  ValNodePtr           vnp, found_match, adj_list = NULL;
  Int4                 len;
  
  vnp = sub_list;
  len = list_len;
  while (vnp != NULL && len > 0) {
    g = (GlobalDiscrepancyPtr) vnp->data.ptrvalue;
    if (g->data_choice == OBJ_SEQFEAT && g->data != NULL) {
      sfp = g->data;
      sfp_next = GetNextGene (sfp);
      if (sfp_next != NULL) {
        found_match = FindValNodeForGlobalDiscrepancyFeature (sub_list, list_len, sfp_next);
        if (found_match != NULL) {
          if (vnp->choice == 0) {
            ValNodeAddPointer (&adj_list, OBJ_SEQFEAT, sfp);
            vnp->choice = 1;
          }
          if (found_match->choice == 0) {
            ValNodeAddPointer (&adj_list, OBJ_SEQFEAT, sfp_next);
            found_match->choice = 1;
          }
        }
      }
    }
    vnp = vnp->next;
    len --;
  }
  /* set choices back to zero */
  for (vnp = sub_list, len = 0; vnp != NULL && len < list_len; vnp = vnp->next, len++) {
    vnp->choice = 0;
  }
  return adj_list;
}
 

extern ClickableItemPtr FindAdjacentDuplicateLocusTagGenes (ValNodePtr locus_tag_list)
{
  ValNodePtr       vnp, adjacent_list = NULL;
  ClickableItemPtr cip = NULL;
  Int4             num_dup;
  CharPtr          duplicate_adjacent_fmt = "%d genes are adjacent to another gene with the same locus tag.";

  vnp = locus_tag_list;
  while (vnp != NULL) {
    num_dup = CountDupGlobalDiscrepancy (vnp);
    if (num_dup > 1) {
      ValNodeLink (&adjacent_list, FindAdjacentGenesInSubList (vnp, num_dup));      
      while (num_dup > 0) {
        vnp = vnp->next;
        num_dup--;
      }
    } else {
      vnp = vnp->next;
    }
  }

  if (adjacent_list != NULL) {
    cip = NewClickableItem (DISC_GENE_DUPLICATE_LOCUS_TAG, duplicate_adjacent_fmt, adjacent_list);
  }
  return cip;
}


NLM_EXTERN ValNodePtr ValNodeDupStringList (ValNodePtr vnp)
{
  ValNodePtr cpy = NULL, last = NULL, tmp;

  while (vnp != NULL)
  {
    tmp = ValNodeNew (NULL);
    tmp->choice = vnp->choice;
    tmp->data.ptrvalue = StringSave (vnp->data.ptrvalue);
    if (last == NULL)
    {
      cpy = tmp;
    }
    else
    {
      last->next = tmp;
    }
    last = tmp;
    vnp = vnp->next;
  }
  return cpy;
}
      

NLM_EXTERN ValNodePtr FindBadLocusTagsInList (ValNodePtr list)
{
  ValNodePtr       bad_list = NULL, list_copy;
  Boolean          reported_last = FALSE;
  ValNodePtr       vnp;
  
  list_copy = ValNodeDupStringList (list);
  list_copy = ValNodeSort (list_copy, SortVnpByString);

  for (vnp = list_copy; vnp != NULL; vnp = vnp->next) {
    /* look for badly formatted locus tags */
    if (IsLocusTagFormatBad (vnp->data.ptrvalue)) {
      ValNodeAddPointer(&bad_list, eLocusTagErrorBadFormat, StringSave (vnp->data.ptrvalue));
    }
  }
  list_copy = ValNodeFreeData (list_copy);
  return bad_list;  
}


extern CharPtr discReportDuplicateLocusTagFmt = "%d genes have duplicate locus tags.";
extern CharPtr discReportOneDuplicateLocusTagFmt = "%d genes have locus tag %s.";
extern CharPtr discReportMissingLocusTags = "%d genes have no locus tags.";

extern void AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  LocusTagCheckData  ltcd;
  ClickableItemPtr dip = NULL, dip_sub;
  ValNodePtr         duplicate_list = NULL;
  CharPtr            bad_fmt = "%d locus tags are incorrectly formatted.";
  ValNodePtr         vnp;
  
  if (discrepancy_list == NULL) return;
  ltcd.locus_tags_list = NULL;
  ltcd.missing_list = NULL;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitGenProdSetFeatures (vnp->data.ptrvalue, &ltcd, CheckGeneLocusTag);
  }
  
  if (ltcd.locus_tags_list != NULL) {
    ltcd.locus_tags_list = ValNodeSort (ltcd.locus_tags_list, SortVnpByGlobalDiscrepancyString);
    ltcd.missing_list = ValNodeSort (ltcd.missing_list, SortVnpByGlobalDiscrepancyString);
    
    if (ltcd.missing_list != NULL) {
      dip = ReportMissingFields (ltcd.missing_list, discReportMissingLocusTags, DISC_GENE_MISSING_LOCUS_TAG);
      if (dip != NULL) {
        ValNodeAddPointer (discrepancy_list, 0, dip);
      }
    }
    dip = ReportNonUniqueGlobalDiscrepancy (ltcd.locus_tags_list, 
                                            discReportDuplicateLocusTagFmt, 
                                            discReportOneDuplicateLocusTagFmt,
                                            DISC_GENE_DUPLICATE_LOCUS_TAG,
                                            FALSE);
    if (dip != NULL) {
      dip_sub = FindAdjacentDuplicateLocusTagGenes (ltcd.locus_tags_list);
      if (dip_sub != NULL) {
        ValNodeAddPointer (&(dip->subcategories), 0, dip_sub);
      }
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }

    /* inconsistent locus tags */
    ValNodeLink (discrepancy_list,
                 ReportInconsistentGlobalDiscrepancyPrefixes (ltcd.locus_tags_list,
                                                              discReportInconsistentLocusTagPrefixFmt,
                                                              DISC_GENE_LOCUS_TAG_INCONSISTENT_PREFIX));
    /* bad formats */
    dip = ReportBadLocusTagFormat (ltcd.locus_tags_list);
    if (dip != NULL) {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  ltcd.locus_tags_list = FreeGlobalDiscrepancyList (ltcd.locus_tags_list);
  ltcd.missing_list = FreeGlobalDiscrepancyList (ltcd.missing_list);
}

static void AddDiscrepancyForNonGeneLocusTag (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR    locus_tag_list;
  GBQualPtr          qual;
  
  if (sfp == NULL || userdata == NULL || sfp->data.choice == SEQFEAT_GENE)
  {
    return;
  }
  
  locus_tag_list = (ValNodePtr PNTR) userdata;
  
  for (qual = sfp->qual; qual != NULL; qual = qual->next)
  {
    if (StringICmp(qual->qual, "locus_tag") == 0) 
    {
      ValNodeAddPointer (locus_tag_list, OBJ_SEQFEAT, sfp);
      return;
    }
  }
}

extern void AddDiscrepanciesForNonGeneLocusTags (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr locus_tag_list = NULL, vnp;
  CharPtr    bad_fmt = "%d non-gene features have locus tags.";
  ClickableItemPtr dip;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &locus_tag_list, AddDiscrepancyForNonGeneLocusTag);
  }
  
  if (locus_tag_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_NON_GENE_LOCUS_TAG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (locus_tag_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = locus_tag_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


static Boolean 
IsGeneLocationOk 
(SeqMgrFeatContextPtr feat_context, 
 SeqMgrFeatContextPtr gene_context,
 BioseqPtr            bsp)
{
  SeqFeatPtr        rbs_sfp;
  SeqMgrFeatContext rbs_context;
  
  if (feat_context == NULL || gene_context == NULL)
  {
    return FALSE;
  }  
  else if ((feat_context->strand == Seq_strand_minus && gene_context->strand != Seq_strand_minus)
           || (feat_context->strand != Seq_strand_minus && gene_context->strand == Seq_strand_minus))
  {
    return FALSE;
  }
  else if (gene_context->left == feat_context->left && gene_context->right == feat_context->right)
  {
    return TRUE;
  }
  else if ((gene_context->strand == Seq_strand_minus && gene_context->left == feat_context->left)
           || (gene_context->strand != Seq_strand_minus && gene_context->right == feat_context->right))
  {
    /* find RBS to extend gene on 5' end */
    for (rbs_sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_RBS, &rbs_context);
         rbs_sfp != NULL;
         rbs_sfp = SeqMgrGetNextFeature (bsp, rbs_sfp, 0, FEATDEF_RBS, &rbs_context))
    {
      if (rbs_context.strand != gene_context->strand)
      {
        continue;
      }
      if (rbs_context.strand == Seq_strand_minus)
      {
        if (rbs_context.right == gene_context->right 
            && rbs_context.left >= feat_context->right)
        {
          return TRUE;
        }
      }
      else
      {
        if (rbs_context.left == gene_context->left
            && rbs_context.right <= feat_context->left)
        {
          return  TRUE;
        }
      }
    }
  }
  return FALSE;
}

static ClickableItemPtr GeneLocationDiscrepancy (Uint1 feature_type, SeqFeatPtr gene, SeqFeatPtr sfp)
{
  ClickableItemPtr cip;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));
  cip->clickable_item_type = DISC_GENE_CDS_mRNA_LOCATION_CONFLICT;
  if (feature_type == SEQFEAT_CDREGION) {
    cip->description = StringSave ("Coding region location does not match gene location");
  } else if (feature_type == SEQFEAT_RNA) {
    cip->description = StringSave ("RNA feature location does not match gene location");
  } else {
    cip->description = StringSave ("Feature location does not match gene location");
  }
  ValNodeAddPointer (&cip->item_list, OBJ_SEQFEAT, sfp);
  ValNodeAddPointer (&cip->item_list, OBJ_SEQFEAT, gene);
  return cip;
}

static ClickableItemPtr MissingGeneXrefDiscrepancy (Uint1 feature_type, SeqFeatPtr sfp)
{
  ClickableItemPtr cip;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));
  cip->clickable_item_type = DISC_GENE_CDS_mRNA_LOCATION_CONFLICT;
  if (feature_type == SEQFEAT_CDREGION) {
    cip->description = StringSave ("Coding region xref gene does not exist");
  } else if (feature_type == SEQFEAT_RNA) {
    cip->description = StringSave ("RNA feature xref gene does not exist");
  } else {
    cip->description = StringSave ("Feature xref gene does not exist");
  }
  ValNodeAddPointer (&cip->item_list, OBJ_SEQFEAT, sfp);
  return cip;
}


static void
CheckFeatureTypeForLocationDiscrepancies 
(BioseqPtr       bsp, 
 Uint1           feature_type,
 ValNodePtr PNTR discrepancy_list)
{
  SeqMgrFeatContext context, gene_context;
  GeneRefPtr        grp;
  SeqFeatPtr        sfp, gene_sfp;
  Boolean           found_match;
  
  if (bsp == NULL || ISA_aa (bsp->mol) || discrepancy_list == NULL || IsmRNASequenceInGenProdSet(bsp))
  {
    return;
  }

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, feature_type, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, feature_type, 0, &context))
  {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL)
    {
      gene_sfp = SeqMgrGetOverlappingGene (sfp->location, &gene_context);
      if (gene_sfp != NULL && !IsGeneLocationOk (&context, &gene_context, bsp))
      {
        ValNodeAddPointer (discrepancy_list, 0, GeneLocationDiscrepancy(feature_type, gene_sfp, sfp));
      }
    }
    else if (!SeqMgrGeneIsSuppressed (grp))
    {
      found_match = FALSE;
      for (gene_sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &gene_context);
           gene_sfp != NULL && ! found_match;
           gene_sfp = SeqMgrGetNextFeature (bsp, gene_sfp, SEQFEAT_GENE, FEATDEF_GENE, &gene_context))
      {
        if (GeneRefMatch (gene_sfp->data.value.ptrvalue, grp) && gene_context.strand == context.strand)
        {
          if (IsGeneLocationOk (&context, &gene_context, bsp))
          {
            found_match = TRUE;
          }
          else
          {
            ValNodeAddPointer (discrepancy_list, 0, GeneLocationDiscrepancy(feature_type, gene_sfp, sfp));
          }
        }
      }
      if (!found_match) {
        ValNodeAddPointer (discrepancy_list, 0, MissingGeneXrefDiscrepancy(feature_type, sfp));
      }
    }
  }
  
}

static void CDSmRNAGeneLocationDiscrepanciesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR   discrepancy_list;

  if (bsp == NULL || ! ISA_na (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  discrepancy_list = (ValNodePtr PNTR) userdata;
  
  CheckFeatureTypeForLocationDiscrepancies (bsp, SEQFEAT_CDREGION, discrepancy_list);
  CheckFeatureTypeForLocationDiscrepancies (bsp, SEQFEAT_RNA, discrepancy_list);
}

static ValNodePtr ValNodePointerDup (ValNodePtr vnp);

extern void FindCDSmRNAGeneLocationDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr         subcategories = NULL, feature_list = NULL, vnp;
  CharPtr            bad_fmt = "%d coding regions or mRNAs have inconsistent gene locations.";
  ClickableItemPtr dip;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &subcategories, CDSmRNAGeneLocationDiscrepanciesCallback);
  }
  
  if (subcategories != NULL)
  {
    for (vnp = subcategories; vnp != NULL; vnp = vnp->next) {
      dip = vnp->data.ptrvalue;
      if (dip != NULL && dip->item_list != NULL) {
        ValNodeLink (&feature_list, ValNodePointerDup (dip->item_list));
      }
    }
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_GENE_CDS_mRNA_LOCATION_CONFLICT;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (subcategories));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->subcategories = subcategories;

      dip->item_list = feature_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


typedef struct cdsgeneproduct 
{
  ValNodePtr cds_list;
  CharPtr    gene_locus;
  CharPtr    product_name;
} CDSGeneProductData, PNTR CDSGeneProductPtr;


static ValNodePtr CDSGeneProductListFree (ValNodePtr cds_list)
{
  CDSGeneProductPtr cgpp;
  
  if (cds_list == NULL) {
    return cds_list;
  }
  
  cds_list->next = CDSGeneProductListFree(cds_list->next);
  
  cgpp = (CDSGeneProductPtr) cds_list->data.ptrvalue;
  if (cgpp != NULL) {
    cgpp->cds_list = ValNodeFree (cgpp->cds_list);
  }
  ValNodeFreeData (cds_list);
  return NULL;
}


static CharPtr GetGeneLabel (SeqFeatPtr sfp)
{
  GeneRefPtr grp;
  SeqFeatPtr gene_sfp;
  
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL)
  {
    gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
    if (gene_sfp != NULL)
    {
      grp = gene_sfp->data.value.ptrvalue;
    }
  }
  if (grp != NULL)
  {
    if (!StringHasNoText (grp->locus))
    {
      return grp->locus;
    }
  }
  return NULL;
}

static void FindCDSGeneProductConflictsCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CDSGeneProductPtr cgpp, cgpp_compare;
  SeqMgrFeatContext context;
  ValNodePtr PNTR   cds_list;
  ValNodePtr        vnp;
  Boolean           found_match = FALSE;
  CharPtr           gene_label;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || userdata == NULL)
  {
    return;
  }
  
  sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context);
  if (sfp == NULL)
  {
    return;
  }
  
  gene_label = GetGeneLabel (sfp);
  if (StringHasNoText (gene_label)) return;

  cds_list = (ValNodePtr PNTR) userdata;

  if (*cds_list == NULL) {
    cgpp = (CDSGeneProductPtr) MemNew (sizeof (CDSGeneProductData));
    if (cgpp != NULL)
    {
      ValNodeAddPointer (&(cgpp->cds_list), OBJ_SEQFEAT, sfp);
      cgpp->gene_locus = gene_label;
      cgpp->product_name = StringSave (context.label);
      ValNodeAddPointer (cds_list, 0, cgpp);
    }
  } else {
    vnp = *cds_list;
    while (vnp != NULL && !found_match)
    {
      cgpp_compare = (CDSGeneProductPtr) vnp->data.ptrvalue;
      if (cgpp_compare != NULL 
          && StringCmp (cgpp_compare->gene_locus, gene_label) == 0
          && StringCmp (cgpp_compare->product_name, context.label) != 0)
      {
        found_match = TRUE;
        vnp->choice = 1;
        ValNodeAddPointer (&(cgpp_compare->cds_list), OBJ_SEQFEAT, sfp);
      }
      vnp = vnp->next;
    }
    if (!found_match) {
      cgpp = (CDSGeneProductPtr) MemNew (sizeof (CDSGeneProductData));
      if (cgpp != NULL)
      {
        ValNodeAddPointer (&(cgpp->cds_list), OBJ_SEQFEAT, sfp);
        cgpp->gene_locus = gene_label;
        cgpp->product_name = StringSave (context.label);
        ValNodeAddPointer (cds_list, 0, cgpp);
      }
    }
  }  
}

extern void FindCDSGeneProductConflicts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr         cds_list = NULL, non_conflict = NULL, vnp;
  CDSGeneProductPtr  cgpp;
  CharPtr            bad_fmt = "%d coding regions have the same gene name as another coding region but a different product.";
  CharPtr            bad_cat_fmt = "%d coding regions have the same gene name(%s) as another coding region but a different product.";
  ClickableItemPtr dip;
  ValNodePtr         item_list = NULL, cds_vnp;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitGenProdSetFeatures (vnp->data.ptrvalue, &cds_list, FindCDSGeneProductConflictsCallback);
  }

  /* remove CDSs without conflicts */
  non_conflict = ValNodeExtractList (&cds_list, 0);
  non_conflict = CDSGeneProductListFree (non_conflict);
  
  /* for each item, replace structure used for search with just the feature */
  for (vnp = cds_list; vnp != NULL; vnp = vnp->next)
  {
    cgpp = (CDSGeneProductPtr) vnp->data.ptrvalue;
    if (cgpp != NULL)
    {
      for (cds_vnp = cgpp->cds_list; cds_vnp != NULL; cds_vnp = cds_vnp->next) {
          ValNodeAddPointer (&item_list, OBJ_SEQFEAT, cds_vnp->data.ptrvalue);
      }
      
      cgpp->product_name = MemFree (cgpp->product_name);
      
      dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      if (dip != NULL)
      {
        dip->clickable_item_type = DISC_GENE_PRODUCT_CONFLICT;
        dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (cgpp->gene_locus) + 15));
        sprintf (dip->description, bad_cat_fmt, ValNodeLen (cgpp->cds_list), cgpp->gene_locus == NULL ? "" : cgpp->gene_locus);
        dip->item_list = cgpp->cds_list;
        cgpp->cds_list = NULL;

        vnp->choice = 0;
        vnp->data.ptrvalue = dip;      
      } else {
        cgpp->cds_list = ValNodeFree (cgpp->cds_list);
        vnp->choice = 0;
        vnp->data.ptrvalue = NULL;
      }
      /* note - we are not freeing gene_locus because we didn't make a copy */
      cgpp = MemFree (cgpp);
    }
  }
    
  if (cds_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_GENE_PRODUCT_CONFLICT;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (item_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = item_list;
      dip->subcategories = cds_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
}


static void DuplicateGeneLocusCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR gene_list;
  SeqFeatPtr      sfp_compare;
  GeneRefPtr      grp, grp_compare;
  ValNodePtr      prev = NULL, vnp;
  Boolean         found_match = FALSE;
  Uint1           new_choice = 0;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE || sfp->data.value.ptrvalue == NULL
      || userdata == NULL)
  {
    return;
  }
  
  gene_list = (ValNodePtr PNTR) userdata;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (StringHasNoText (grp->locus))
  {
    return;
  }
  
  if (*gene_list == NULL)
  {
    ValNodeAddPointer (gene_list, 0, sfp);
  }
  else
  {
    vnp = *gene_list;
    while (vnp != NULL && !found_match)
    {
      sfp_compare = (SeqFeatPtr) vnp->data.ptrvalue;
      grp_compare = (GeneRefPtr) sfp_compare->data.value.ptrvalue;
      if (StringCmp (grp_compare->locus, grp->locus) == 0)
      {
        found_match = TRUE;
        vnp->choice = OBJ_SEQFEAT;
        new_choice = OBJ_SEQFEAT;
      }
      prev = vnp;
      vnp = vnp->next;
    }
    
    if (found_match)
    {
      vnp = prev;
      /* insert at end of matches */
      while (found_match && vnp != NULL)
      {
        sfp_compare = (SeqFeatPtr) vnp->data.ptrvalue;
        grp_compare = (GeneRefPtr) sfp_compare->data.value.ptrvalue;
        if (StringCmp (grp_compare->locus, grp->locus) != 0)
        {
          found_match = FALSE;
        }
        else
        {
          prev = vnp;
        }
        vnp = vnp->next;
      }
    }

    /* add to list */
    vnp = ValNodeNew (NULL);
    vnp->choice = new_choice;
    vnp->data.ptrvalue = sfp;
    vnp->next = prev->next;
    prev->next = vnp;
      
  }
  
}


extern void FindDuplicateGeneLocus (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr         gene_list = NULL, non_conflict = NULL, vnp;
  CharPtr            bad_fmt = "%d genes have the same locus as another gene.";
  ClickableItemPtr dip;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitGenProdSetFeatures (vnp->data.ptrvalue, &gene_list, DuplicateGeneLocusCallback);
  }

  /* remove Genes without conflicts */
  non_conflict = ValNodeExtractList (&gene_list, 0);
  non_conflict = ValNodeFree (non_conflict);
  
  if (gene_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_GENE_DUPLICATE_LOCUS;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (gene_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = gene_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }  
  
}


static void FindECNumberNotes (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR    ec_number_features;
  BioseqPtr          prot_bsp;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         prot_sfp;
  ProtRefPtr         prp;
  ValNodePtr         vnp;
  
  if (sfp == NULL || userdata == NULL || StringHasNoText (sfp->comment))
  {
    return;
  }
  
  ec_number_features = (ValNodePtr PNTR) userdata;
  
  if (LookForECnumberPattern (sfp->comment))
  {
    ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
  }
  else if (sfp->data.choice == SEQFEAT_CDREGION && sfp->product != NULL) 
  {
    prot_bsp = BioseqFindFromSeqLoc(sfp->product);
    prot_sfp = SeqMgrGetNextFeature(prot_bsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
    if (prot_sfp != NULL && prot_sfp->data.value.ptrvalue != NULL) {
      prp = (ProtRefPtr) prot_sfp->data.value.ptrvalue;
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        if (LookForECnumberPattern (vnp->data.ptrvalue)) {
          ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
          return;
        }
      }
      if (LookForECnumberPattern (prp->desc)) {
        ValNodeAddPointer (ec_number_features, OBJ_SEQFEAT, sfp);
        return;
      }
    }
  }  
}

extern void AddECNumberNoteDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr ec_number_features = NULL, vnp;
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d features have EC numbers in notes or products.";
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &ec_number_features, FindECNumberNotes);
  }
  
  if (ec_number_features != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_EC_NUMBER_NOTE;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (ec_number_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = ec_number_features;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
}


static void FindPseudoDiscrepanciesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR pseudo_features;
  GeneRefPtr      grp;
  SeqFeatPtr      gene_sfp = NULL;
  
  if (sfp == NULL || (sfp->data.choice != SEQFEAT_CDREGION && sfp->data.choice != SEQFEAT_RNA)
      || userdata == NULL)
  {
    return;
  }
  
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL)
  {
    return;
  }
  
  gene_sfp = SeqMgrGetOverlappingGene (sfp->location, NULL);
  if (gene_sfp == NULL)
  {
    return;
  }
  
  if (sfp->pseudo && ! gene_sfp->pseudo)
  {
    pseudo_features = (ValNodePtr PNTR) userdata;
    ValNodeAddPointer (pseudo_features, OBJ_SEQFEAT, sfp);
    if (gene_sfp != NULL)
    {
      ValNodeAddPointer (pseudo_features, OBJ_SEQFEAT, gene_sfp);
    }
  }
}


extern void FindPseudoDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr pseudo_features = NULL, vnp;
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d CDSs, RNAs, and genes have mismatching pseudos.";
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &pseudo_features, FindPseudoDiscrepanciesCallback);
  }
  
  if (pseudo_features != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_PSEUDO_MISMATCH;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (pseudo_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = pseudo_features;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
  
}


static void FindJoinedLocations (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR joined_features;
  
  if (sfp == NULL || userdata == NULL || sfp->location == NULL)
  {
    return;
  }
  
  joined_features = (ValNodePtr PNTR) userdata;
  if (sfp->location->choice == SEQLOC_MIX)
  {
    ValNodeAddPointer (joined_features, OBJ_SEQFEAT, sfp);
  }
}

static Boolean HasRibosomalSlippageException (SeqFeatPtr sfp)
{
  if (sfp == NULL || !sfp->excpt || StringISearch (sfp->except_text, "ribosomal slippage") == NULL)
  {
    return FALSE;
  }
  else 
  {
    return TRUE;
  }
}

extern void AddJoinedFeatureDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr joined_features = NULL, vnp;
  ValNodePtr slippage = NULL, no_slippage = NULL;
  SeqFeatPtr sfp;
  
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d features have joined locations.";
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &joined_features, FindJoinedLocations);
  }
  
  if (joined_features != NULL)
  {
    /* get list of features with and without ribosomal slippage exception */
    for (vnp = joined_features; vnp != NULL; vnp = vnp->next)
    {
      if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) 
      {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (HasRibosomalSlippageException (sfp))
        {
          ValNodeAddPointer (&slippage, OBJ_SEQFEAT, sfp);
        }
        else
        {
          ValNodeAddPointer (&no_slippage, OBJ_SEQFEAT, sfp);
        }
      }
    }

    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_JOINED_FEATURES;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (joined_features));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = joined_features;

      if (slippage != NULL) {
        ValNodeAddPointer (&(dip->subcategories), 0, NewClickableItem (DISC_JOINED_FEATURES, "%d features have joined location and ribosomal slippage exception", slippage));
      }
      if (no_slippage != NULL) {
        ValNodeAddPointer (&(dip->subcategories), 0, NewClickableItem (DISC_JOINED_FEATURES, "%d features have joined location but no ribosomal slippage exception", no_slippage));
      }
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


static void FindOverlappingGenes (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp, sfp_compare;
  SeqMgrFeatContext  context;
  ValNodePtr PNTR    overlapping_genes = NULL, non_overlap;
  ValNodePtr         gene_list = NULL, vnp, vnp_next;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  overlapping_genes = (ValNodePtr PNTR) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_GENE, FEATDEF_GENE, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_GENE, FEATDEF_GENE, &context))
  {
    ValNodeAddPointer (&gene_list, 0, sfp);
  }
  
  for (vnp = gene_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      sfp_compare = (SeqFeatPtr) vnp_next->data.ptrvalue;
      
      if (SeqLocStrand (sfp->location) != SeqLocStrand (sfp_compare->location))
      {
        continue;
      }
      
      if (SeqLocCompare (sfp->location, sfp_compare->location) != SLC_NO_MATCH)
      {
        vnp->choice = OBJ_SEQFEAT;
        vnp_next->choice = OBJ_SEQFEAT;
      }
    }
  }
  
  non_overlap = ValNodeExtractList (&gene_list, 0);
  non_overlap = ValNodeFree (non_overlap);
  ValNodeLink (overlapping_genes, gene_list);
  
}


extern void AddOverlappingGeneDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d genes overlap another gene on the same strand.";
  ValNodePtr         overlapping_genes = NULL, vnp;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &overlapping_genes, FindOverlappingGenes);
  }
  
  if (overlapping_genes != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_OVERLAPPING_GENES;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (overlapping_genes));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = overlapping_genes;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


typedef struct cdsoverlap 
{
  CharPtr    product_name;
  SeqFeatPtr sfp;  
  Int4       left;
  Int4       right;
} CDSOverlapData, PNTR CDSOverlapPtr;


static CDSOverlapPtr CDSOverlapNew (SeqFeatPtr sfp, CharPtr product_name, Int4 left, Int4 right)
{
  CDSOverlapPtr cop;
  
  cop = (CDSOverlapPtr) MemNew (sizeof (CDSOverlapData));
  if (cop != NULL)
  {
    cop->product_name = StringSave (product_name);
    cop->sfp = sfp;
    cop->left = left;
    cop->right = right;
  }
  return cop;
}


static ValNodePtr FreeCDSOverlapList (ValNodePtr vnp)
{
  CDSOverlapPtr cop;
  
  if (vnp != NULL)  
  {
    vnp->next = FreeCDSOverlapList (vnp->next);
    cop = (CDSOverlapPtr) vnp->data.ptrvalue;
    if (cop != NULL)
    {
      cop->product_name = MemFree (cop->product_name);
      cop = MemFree (cop);
      vnp->data.ptrvalue = NULL;
    }
    vnp = ValNodeFree (vnp);
  }
  return vnp;
}


static ValNodePtr FeatureListFromOverlapList (ValNodePtr vnp)
{
  ValNodePtr     feat_list = NULL;
  CDSOverlapPtr cop;
  
  while (vnp != NULL)
  {
    if (vnp->choice != 0 && vnp->data.ptrvalue != NULL)
    {
      cop = (CDSOverlapPtr) vnp->data.ptrvalue;
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, cop->sfp);
    }
    vnp = vnp->next;
  }
  return feat_list;
}


static CharPtr similar_product_words[] = 
{ "transposase",
  "integrase"
};

const int num_similar_product_words = sizeof (similar_product_words) / sizeof (CharPtr);

static CharPtr ignore_similar_product_words[] = 
{ "hypothetical protein",
  "phage"
};

const int num_ignore_similar_product_words = sizeof (ignore_similar_product_words) / sizeof (CharPtr);


static Boolean OverlappingProductNameSimilar (CharPtr str1, CharPtr str2)
{
  Int4 i;
  Boolean str1_has_similarity_word = FALSE, str2_has_similarity_word = FALSE;
  
  if (StringHasNoText (str1) && StringHasNoText (str2))
  {
    return TRUE;
  }
  else if (StringHasNoText (str1) || StringHasNoText (str2))
  {
    return FALSE;
  }
  
  /* if both product names contain one of the special case similarity words,
   * the product names are similar. */
  for (i = 0; i < num_similar_product_words; i++)
  {
    if (StringISearch (str1, similar_product_words [i]) != NULL)
    {
      str1_has_similarity_word = TRUE;
    }
    if (StringISearch (str2, similar_product_words [i]) != NULL)
    {
      str2_has_similarity_word = TRUE;
    }
  }
  if (str1_has_similarity_word && str2_has_similarity_word)
  {
    return TRUE;
  }
  
  /* otherwise, if one of the product names contains one of special ignore similarity
   * words, the product names are not similar.
   */
  for (i = 0; i < num_ignore_similar_product_words; i++)
  {
    if (StringISearch (str1, ignore_similar_product_words[i]) != NULL
        || StringISearch (str2, ignore_similar_product_words[i]) != NULL)
    {
      return FALSE;
    }
  }
  
  if (StringICmp (str1, str2) == 0)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static void FindOverlappingCDSs (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;
  ValNodePtr PNTR    overlapping_cds = NULL, cds_list;
  ValNodePtr         overlap_list = NULL, vnp, vnp_next;
  CDSOverlapPtr      cop, cop_compare;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  overlapping_cds = (ValNodePtr PNTR) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &context))
  {
    ValNodeAddPointer (&overlap_list, 0, CDSOverlapNew (sfp, context.label, context.left, context.right));
  }
  
  for (vnp = overlap_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {
    cop = (CDSOverlapPtr) vnp->data.ptrvalue;
    if (cop == NULL)
    {
      continue;
    }
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      cop_compare = (CDSOverlapPtr) vnp_next->data.ptrvalue;
      if (cop_compare == NULL)
      {
        continue;
      }
      else if (cop_compare->left > cop->right)
      {
        break;
      }
      if (!OverlappingProductNameSimilar (cop->product_name, cop_compare->product_name))
      {
        continue;
      }
      if (SeqLocStrand (cop->sfp->location) != SeqLocStrand (cop_compare->sfp->location))
      {
        continue;
      }
      
      if (SeqLocCompare (cop->sfp->location, cop_compare->sfp->location) != SLC_NO_MATCH)
      {
        vnp->choice = OBJ_SEQFEAT;
        vnp_next->choice = OBJ_SEQFEAT;
      }
    }
  }
  
  cds_list = FeatureListFromOverlapList(overlap_list);
  if (cds_list != NULL)
  {
    ValNodeLink (overlapping_cds, cds_list);
  }
  overlap_list = FreeCDSOverlapList (overlap_list);
}


extern void AddOverlappingCodingRegionDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d coding regions overlap another coding region with a similar or identical name.";
  ValNodePtr         overlapping_cds = NULL, vnp;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &overlapping_cds, FindOverlappingCDSs);
  }
  
  if (overlapping_cds != NULL)
  {
    dip = NewClickableItem (DISC_OVERLAPPING_CDS, bad_fmt, overlapping_cds);

    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}

static void FindContainedCDSs (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp, sfp_compare;
  SeqMgrFeatContext  context;
  ValNodePtr PNTR    contained_cds = NULL;
  ValNodePtr         cds_list = NULL;
  ValNodePtr         contained_list = NULL, vnp, vnp_next;
  Int2               loc_compare;
  Uint1              strand, strand_compare;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  contained_cds = (ValNodePtr PNTR) userdata;
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &context))
  {
    ValNodeAddPointer (&cds_list, OBJ_SEQFEAT, sfp);
  }
  
  for (vnp = cds_list; vnp != NULL && vnp->next != NULL; vnp = vnp->next)
  {    
    sfp = vnp->data.ptrvalue;
    strand = SeqLocStrand (sfp->location);
    for (vnp_next = vnp->next; vnp_next != NULL; vnp_next = vnp_next->next)
    {
      sfp_compare = vnp_next->data.ptrvalue;
      strand_compare = SeqLocStrand (sfp_compare->location);
      if ((strand == Seq_strand_minus && strand_compare != Seq_strand_minus)
          || (strand != Seq_strand_minus && strand_compare == Seq_strand_minus)) 
      {
        continue;
      }
      loc_compare = SeqLocCompare (sfp->location, sfp_compare->location);
      if (loc_compare == SLC_A_IN_B || loc_compare == SLC_B_IN_A || loc_compare == SLC_A_EQ_B)
      {
        if (!AlreadyInList (contained_list, sfp)) 
        {
          ValNodeAddPointer (&contained_list, OBJ_SEQFEAT, sfp);
        }
        if (!AlreadyInList (contained_list, sfp_compare)) 
        {
          ValNodeAddPointer (&contained_list, OBJ_SEQFEAT, sfp_compare);
        }
      }
    }
  }
  cds_list = ValNodeFree (cds_list);
  
  if (contained_list != NULL)
  {
    ValNodeLink (contained_cds, contained_list);
  }
}


extern void AddContainedCodingRegionDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d coding regions are completely contained in another coding region.";
  ValNodePtr         contained_cds = NULL, vnp;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &contained_cds, FindContainedCDSs);
  }
  
  if (contained_cds != NULL)
  {
    dip = NewClickableItem (DISC_CONTAINED_CDS, bad_fmt, contained_cds);

    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


typedef struct cdsrnaoverlap {
  ValNodePtr cds_in_rna;
  ValNodePtr rna_in_cds;
  ValNodePtr exact_match;
  ValNodePtr overlap_same_strand;
  ValNodePtr overlap_opp_strand;
  ValNodePtr overlap;
  ValNodePtr all;
} CDSRNAOverlapData, PNTR CDSRNAOverlapPtr;

static void FindCDSRNAOverlaps (BioseqPtr bsp, Pointer data)
{
  CDSRNAOverlapPtr  p;
  SeqFeatPtr        sfp, rna;
  ValNodePtr        rna_list = NULL, vnp;
  SeqMgrFeatContext fcontext;
  Int2              cmp;
  Uint1             strand1, strand2;

  if (bsp == NULL || data == NULL) return;

  p = (CDSRNAOverlapPtr) data;

  for (rna = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, 0, &fcontext);
       rna != NULL;
       rna = SeqMgrGetNextFeature (bsp, rna, SEQFEAT_RNA, 0, &fcontext))
  {
    if (rna->idx.subtype == FEATDEF_mRNA) continue;
    ValNodeAddPointer (&rna_list, OBJ_SEQFEAT, rna);
  }

  if (rna_list == NULL) return;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext))
  {
    for (vnp = rna_list; vnp != NULL; vnp = vnp->next)
    {
      rna = (SeqFeatPtr) vnp->data.ptrvalue;
      cmp = SeqLocCompare (sfp->location, rna->location);
      if (cmp == SLC_A_EQ_B)
      {
        ValNodeAddPointer (&(p->exact_match), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->exact_match), OBJ_SEQFEAT, rna);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, rna);
      }
      else if (cmp == SLC_A_IN_B)
      {
        ValNodeAddPointer (&(p->cds_in_rna), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->cds_in_rna), OBJ_SEQFEAT, rna);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, rna);
      }
      else if (cmp == SLC_B_IN_A)
      {
        ValNodeAddPointer (&(p->rna_in_cds), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->rna_in_cds), OBJ_SEQFEAT, rna);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, rna);
      }
      else if (cmp != SLC_NO_MATCH)
      {
        strand1 = SeqLocStrand (sfp->location);
        strand2 = SeqLocStrand (rna->location);
        if ((strand1 == Seq_strand_minus && strand2 != Seq_strand_minus)
          || (strand2 == Seq_strand_minus && strand1 != Seq_strand_minus))
        {
          ValNodeAddPointer (&(p->overlap_opp_strand), OBJ_SEQFEAT, sfp);
          ValNodeAddPointer (&(p->overlap_opp_strand), OBJ_SEQFEAT, rna);
        }
        else
        {
          ValNodeAddPointer (&(p->overlap_same_strand), OBJ_SEQFEAT, sfp);
          ValNodeAddPointer (&(p->overlap_same_strand), OBJ_SEQFEAT, rna);
        }
        ValNodeAddPointer (&(p->overlap), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->overlap), OBJ_SEQFEAT, rna);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, sfp);
        ValNodeAddPointer (&(p->all), OBJ_SEQFEAT, rna);
      }
    }
  }
  rna_list = ValNodeFree (rna_list);
}


static ClickableItemPtr DiscrepancyForPairs (Uint4 item_type, CharPtr bad_fmt, ValNodePtr item_list)
{
  ClickableItemPtr dip;
  Int4             num_feat;

  if (StringHasNoText (bad_fmt)) return NULL;
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  dip->clickable_item_type = item_type;
  dip->item_list = item_list;
  num_feat = ValNodeLen (dip->item_list) / 2;
  dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
  sprintf (dip->description, bad_fmt, num_feat);
  return dip;
}


extern void AddRNACDSOverlapDiscrepancies (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip, overlap_dip;
  ValNodePtr       vnp;
  CDSRNAOverlapData d;

  if (discrepancy_list == NULL)
  {
    return;
  }
  
  MemSet (&d, 0, sizeof (CDSRNAOverlapData));
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &d, FindCDSRNAOverlaps);
  }
  
  if (d.all != NULL)
  {
    dip = DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, "%d coding regions overlap RNA features", d.all);
    if (d.exact_match != NULL)
    {
      ValNodeAddPointer (&(dip->subcategories), 0, 
                         DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                              "%d coding region locations exactly match an RNA location",
                                              d.exact_match));
    }
    if (d.cds_in_rna != NULL)
    {
      ValNodeAddPointer (&(dip->subcategories), 0, 
                         DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP,
                                              "%d coding regions are completely contained in RNAs",
                                              d.cds_in_rna));
    }
    if (d.rna_in_cds != NULL)
    {
      ValNodeAddPointer (&(dip->subcategories), 0, 
                         DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                              "%d coding regions completely contain RNAs",
                                              d.rna_in_cds));
    }
    if (d.overlap != NULL)
    {
      overlap_dip = DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                          "%d coding regions overlap RNAs (no containment)",
                                          d.overlap);
      if (d.overlap_same_strand != NULL) {
        ValNodeAddPointer (&(overlap_dip->subcategories), 0, 
                            DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                              "%d coding regions overlap RNAs on the same strand (no containment)",
                                              d.overlap_same_strand));
      }
      if (d.overlap_same_strand != NULL) {
        ValNodeAddPointer (&(overlap_dip->subcategories), 0, 
                            DiscrepancyForPairs (DISC_RNA_CDS_OVERLAP, 
                                              "%d coding regions overlap RNAs on the opposite strand (no containment)",
                                              d.overlap_opp_strand));                                                   
      }
      ValNodeAddPointer (&(dip->subcategories), 0, overlap_dip);
    }
    
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}


static void FindShortContigsCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR bioseq_list;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL || bsp->length >= 200) {
    return;
  }

  if (IsmRNASequenceInGenProdSet (bsp)) {
    return;
  }
  
  bioseq_list = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (bioseq_list, OBJ_BIOSEQ, bsp);
}

extern void FindShortContigs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d contigs are shorter than 200 nt.";
  ValNodePtr         bioseq_list = NULL, vnp;
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &bioseq_list, FindShortContigsCallback);
  }
  
  if (bioseq_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_SHORT_CONTIG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (bioseq_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = bioseq_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}

static void FindShortSequencesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR bioseq_list;
  BioseqSetPtr    bssp;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL || bsp->length >= 50
      || IsmRNASequenceInGenProdSet (bsp))
  {
    return;
  }
  
  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) {
      return;
    }
  }
  
  bioseq_list = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (bioseq_list, OBJ_BIOSEQ, bsp);
}

extern void FindShortSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  CharPtr            bad_fmt = "%d sequences are shorter than 50 nt.";
  ValNodePtr         bioseq_list = NULL, vnp;
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &bioseq_list, FindShortSequencesCallback);
  }
  
  if (bioseq_list != NULL)
  {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    if (dip != NULL)
    {
      dip->clickable_item_type = DISC_SHORT_CONTIG;
      dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
      sprintf (dip->description, bad_fmt, ValNodeLen (bioseq_list));
      dip->callback_func = NULL;
      dip->datafree_func = NULL;
      dip->callback_data = NULL;
      dip->item_list = bioseq_list;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}

typedef struct sdpandbsp {
  SeqDescrPtr sdp;
  BioseqPtr bsp;
} SdpAndBspData, PNTR SdpAndBspPtr;

typedef struct biosrccheck 
{
  BioSourcePtr biop;
  ValNodePtr   sdp_list;
} BioSrcCheckData, PNTR BioSrcCheckPtr;

static ValNodePtr FreeBioSrcCheckList (ValNodePtr biosrc_list)
{
  BioSrcCheckPtr  bscp;
  
  if (biosrc_list == NULL)
  {
    return NULL;
  }
  
  biosrc_list->next = FreeBioSrcCheckList (biosrc_list->next);
  
  bscp = (BioSrcCheckPtr) biosrc_list->data.ptrvalue;
  if (bscp != NULL)
  {
    bscp->sdp_list = ValNodeFreeData (bscp->sdp_list);
    bscp = MemFree (bscp);
  }
  biosrc_list = ValNodeFree (biosrc_list);
  return NULL;
}


static void FindInconsistentSourcesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR biosrc_list, vnp;
  SeqDescrPtr     sdp;
  BioSrcCheckPtr  bscp;
  Boolean         found = FALSE;
  SeqMgrDescContext context;
  SdpAndBspPtr      sabp;
  
  if (bsp == NULL || !ISA_na (bsp->mol) || userdata == NULL)
  {
    return;
  }
  
  biosrc_list = (ValNodePtr PNTR) userdata;
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp != NULL)
  {
    sabp = (SdpAndBspPtr) MemNew (sizeof (SdpAndBspData));
    sabp->sdp = sdp;
    sabp->bsp = bsp;

    for (vnp = *biosrc_list; vnp != NULL && !found; vnp = vnp->next)
    {
      bscp = (BioSrcCheckPtr) vnp->data.ptrvalue;
      if (bscp != NULL && BioSourceMatch (sdp->data.ptrvalue, bscp->biop))
      {
        ValNodeAddPointer (&(bscp->sdp_list), 0, sabp);
        found = TRUE;
      }
    }
    if (!found)
    {
      bscp = (BioSrcCheckPtr) MemNew (sizeof (BioSrcCheckData));
      if (bscp != NULL)
      {
        bscp->biop = sdp->data.ptrvalue;
        ValNodeAddPointer (&(bscp->sdp_list), 0, sabp);
        ValNodeAddPointer (biosrc_list, 0, bscp);
      }
    }
  }
}


static ClickableItemPtr InconsistentBiosrc (BioSrcCheckPtr bscp)
{
  ClickableItemPtr dip = NULL;
  CharPtr          bad_fmt = "%d contigs have identical sources that do not match another contig source.";
  ValNodePtr       vnp;
  SdpAndBspPtr     sabp;

  if (bscp == NULL || bscp->sdp_list == NULL)
  {
    return NULL;
  }
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = DISC_INCONSISTENT_BIOSRC;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (bscp->sdp_list));
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    for (vnp = bscp->sdp_list; vnp != NULL; vnp = vnp->next) {
      sabp = (SdpAndBspPtr) vnp->data.ptrvalue;
      ValNodeAddPointer (&(dip->item_list), OBJ_SEQDESC, sabp->sdp);
      ValNodeAddPointer (&(dip->item_list), OBJ_BIOSEQ, sabp->bsp);
    }
  }      
  return dip;
}


extern void FindNonmatchingContigSources (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ClickableItemPtr dip;
  ValNodePtr       biosrc_list = NULL, vnp, vnp_s, sub_list = NULL, item_list = NULL;
  BioSrcCheckPtr   bscp;
  SdpAndBspPtr     sabp;
  CharPtr          disc_fmt = "%d inconsistent contig sources";

  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &biosrc_list, FindInconsistentSourcesCallback);
  }
  
  if (biosrc_list != NULL && biosrc_list->next != NULL)
  {
    for (vnp = biosrc_list; vnp != NULL; vnp = vnp->next)
    {
      bscp = (BioSrcCheckPtr) vnp->data.ptrvalue;
      dip = InconsistentBiosrc (bscp);
      ValNodeAddPointer (&sub_list, 0, dip);

      /* add sdp and bsp to item list */
      for (vnp_s = bscp->sdp_list; vnp_s != NULL; vnp_s = vnp_s->next) {
        sabp = (SdpAndBspPtr) vnp_s->data.ptrvalue;
        ValNodeAddPointer (&item_list, OBJ_SEQDESC, sabp->sdp);
        ValNodeAddPointer (&item_list, OBJ_BIOSEQ, sabp->bsp);
      }
    }
  }
  biosrc_list = FreeBioSrcCheckList (biosrc_list);
  if (item_list != NULL) {
    dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (dip, 0, sizeof (ClickableItemData));
    dip->clickable_item_type = DISC_INCONSISTENT_BIOSRC;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (disc_fmt) + 15));
    sprintf (dip->description, disc_fmt, ValNodeLen (item_list) / 2);
    dip->item_list = item_list;
    dip->subcategories = sub_list;
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}

static CharPtr suspect_product_names[] = 
{
"Similar to",
"Related to",
"interrupt",
"Homolog",
"Homologue",
"Fragment",
"Frameshift",
"Intein",
"N-term",
"N term",
"C-term",
"C term",
"Chloroplast",
"Mitochondrial",
"may contain a plural",
"Brackets or parenthesis [] ()",
"ending with period, comma, or hyphen",
"beginning with period, comma, or hyphen",
"ending with like",
"COG",     /* wholeword start */
"Subtilis",
"coli",
"pseudo",
"Yersinia",
"gene",
"genes",
"protein", /* singleword start */
"putative",
"probable", /* singleword end */
"unknown", /* unknown */
"ortholog",
"orthologue",
"paralog",
"paralogue",
"bifunctional protein",
"pseudogene",
"frame shift",
"protien",
"partial",
"B.subtilis",
"E.coli",
"Escherichia",
"Bacillus",
"Staphlococcus",
"aureus",
"Salmonella",
"Streptococcus",
"Staphlococcal",
"streptococcal",
"Helicobacter",
"pylori",
"Campylobacter",
"Jejuni",
"Pestis",
"Rhodobacter",
"sphaeroides",
"or related",
"authentic point mutation",
"novel protein",
"ttg start",
"domain protein domain protein",
"deletion",
"truncat"
};

const int num_suspect_product_names = sizeof (suspect_product_names) / sizeof (CharPtr);

const int term_start = 8;
const int term_end = 11;

const int plural_name = 14;
const int brackets_name = 15;
const int end_with_punct_name = 16;
const int begin_with_punct_name = 17;
const int end_with_like = 18;

const int wholeword_start = 19;
const int wholeword_end = 25;

const int singleword_start = 26;
const int singleword_end = 28;

const int unknown_name = 29;

static Boolean MayContainPlural (CharPtr str)
{
  CharPtr cp;
  Int4    word_len = 0;
  Boolean may_contain_plural = FALSE;

  if (str == NULL) return FALSE;
  cp = str;
  while (*cp != 0 && !may_contain_plural) {
    if (isalpha (*cp)) {
      if (word_len > 2 
        && *cp == 's'
        && *(cp - 1) != 's'
        && *(cp - 1) != 'i'
        && *(cp - 1) != 'u'
        && (*(cp + 1) == ',' || *(cp + 1) == 0)) {
        may_contain_plural = TRUE;
      } else {
        word_len++;
      }
    } else {
      word_len = 0;
    }
    cp++;
  }
  return may_contain_plural;
}

static void FindSuspectProductNamesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR feature_list;
  Int4            k, len;
  ProtRefPtr      prp;
  ValNodePtr      vnp;
  CharPtr         str;
  BioseqPtr       bsp;
  SeqFeatPtr      cds;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL
      || userdata == NULL)
  {
    return;
  }
  
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  feature_list = (ValNodePtr PNTR) userdata;

  /* add coding region rather than protein */
  if (sfp->idx.subtype == FEATDEF_PROT) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
      if (cds != NULL) {
        sfp = cds;
      }
    }
  }
  
  for (k = 0; k < num_suspect_product_names; k++)
  {
    if (k == plural_name) 
    {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      {
        if (vnp->data.ptrvalue == NULL) continue;
        if (MayContainPlural(vnp->data.ptrvalue))
        {
          ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
          break;
        }
      }
    }
    else if (k == brackets_name)
    {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      {
        if (StringChr (vnp->data.ptrvalue, '[') != NULL
            || StringChr (vnp->data.ptrvalue, ']') != NULL
            || StringChr (vnp->data.ptrvalue, '(') != NULL
            || StringChr (vnp->data.ptrvalue, ')') != NULL)
        {
          ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
          break;
        }
      }
    }
    else if (k == end_with_punct_name)
    {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      {
        if (vnp->data.ptrvalue == NULL) continue;
        len = StringLen (vnp->data.ptrvalue);
        str = (CharPtr) vnp->data.ptrvalue;
        if (str[len - 1] == '.' || str[len - 1] == ',' || str[len - 1] == '-')
        {
          ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
          break;
        }
      }
    }
    else if (k == begin_with_punct_name)
    {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      {
        if (vnp->data.ptrvalue == NULL) continue;
        str = (CharPtr) vnp->data.ptrvalue;
        if (str[0] == '.' || str[0] == ',' || str[0] == '-')
        {
          ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
          break;
        }
      }
    }
    else if (k == end_with_like) 
    {
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
      {
        if (vnp->data.ptrvalue == NULL) continue;
        str = (CharPtr) vnp->data.ptrvalue;
        len = StringLen (str);
        if (len > 4 && StringCmp (str + len - 4, "like") == 0) 
        {
          ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
          break;
        }
      }
    }
    else
    {
      if (k >= term_start && k <= term_end) {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
        {
          /* don't bother searching for c-term or n-term if product name contains "domain" */
          if (StringISearch (vnp->data.ptrvalue, "domain") != NULL) continue;
          str = StringISearch(vnp->data.ptrvalue, suspect_product_names[k]);
          /* c-term and n-term must be either first word or separated from other word by space, num, or punct */
          if (str != NULL && (str == vnp->data.ptrvalue || !isalpha (*(str - 1))))
          {
            ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
            break;
          }
        }
      }
      else if (k >= wholeword_start && k <= wholeword_end)
      {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
        {
          str = StringISearch(vnp->data.ptrvalue, suspect_product_names[k]);
          len = StringLen (suspect_product_names[k]);
          if (str != NULL 
              && (str == vnp->data.ptrvalue || isalpha (*(str - 1)))
              && (str [len - 1] == 0 || !isalpha (str [len - 1])))
          {
            ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
            break;
          }
        }
      }
      else if (k >= singleword_start && k <= singleword_end)
      {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
        {
          if (StringICmp (prp->name->data.ptrvalue, suspect_product_names[k]) == 0) {
            ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
          }
        }
      }
      else
      {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
        {
          if (StringISearch(vnp->data.ptrvalue, suspect_product_names[k]) != NULL)
          {
            if (k != unknown_name 
                || (StringISearch (vnp->data.ptrvalue, "protein of unknown function") == NULL
                    && StringISearch (vnp->data.ptrvalue, "domain of unknown function") == NULL)) {
              ValNodeAddPointer (&(feature_list[k]), OBJ_SEQFEAT, sfp);
              break;
            }
          }
        }
      }
    }
  }
}


static ClickableItemPtr SuspectPhrase (Uint4 clickable_item_type, CharPtr phrase, CharPtr feat_type, ValNodePtr feature_list)
{
  ClickableItemPtr dip = NULL;
  CharPtr            bad_fmt = "%d %ss contain %s";

  if (feature_list == NULL || StringHasNoText (phrase) || StringHasNoText (feat_type))
  {
    return NULL;
  }
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = clickable_item_type;
    dip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_fmt) + StringLen (phrase) + StringLen (feat_type) + 15));
    sprintf (dip->description, bad_fmt, ValNodeLen (feature_list), feat_type, phrase);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = feature_list;
  }      
  return dip;
}

extern void FindSuspectProductNames (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr PNTR   feature_list = NULL;
  ValNodePtr         master_list = NULL, vnp;
  Int4               k;
  ClickableItemPtr dip;
  ValNodePtr         subcategories = NULL;
  
  if (discrepancy_list == NULL) return;

  feature_list = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * num_suspect_product_names);
  if (feature_list == NULL) return;
  
  /* initialize array for suspicious product names */
  for (k = 0; k < num_suspect_product_names; k++)
  {
    feature_list[k] = NULL;
  }
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) 
  {
    VisitGenProdSetFeatures (vnp->data.ptrvalue, feature_list, FindSuspectProductNamesCallback);
  }
  
  for (k = 0; k < num_suspect_product_names; k++)
  {
    if (feature_list[k] != NULL)
    {
      dip = SuspectPhrase (DISC_SUSPECT_PRODUCT_NAME, suspect_product_names[k], "product name", feature_list[k]);
      if (dip != NULL)
      {
        ValNodeAddPointer (&subcategories, 0, dip);
      }
      ValNodeLinkCopy (&master_list, feature_list[k]);
    }
  }
  
  if (master_list != NULL)
  {
    dip = SuspectPhrase (DISC_SUSPECT_PRODUCT_NAME, "suspect phrase or characters", "product_name", master_list);
    if (dip != NULL)
    {
      dip->subcategories = subcategories;
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }

  MemFree (feature_list);
}


static CharPtr suspect_phrases[] = 
{
"fragment",
"frameshift"
};

const int num_suspect_phrases = sizeof (suspect_phrases) / sizeof (CharPtr);


static void FindSuspectPhrasesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR feature_list;
  Int4            k;
  ProtRefPtr      prp;
  CharPtr         check_str = NULL;
  
  if (sfp == NULL || (sfp->data.choice != SEQFEAT_PROT && sfp->data.choice != SEQFEAT_CDREGION) || sfp->data.value.ptrvalue == NULL
      || userdata == NULL)
  {
    return;
  }
  
  if (sfp->data.choice == SEQFEAT_PROT) {  
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    check_str = prp->desc;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    check_str = sfp->comment;
  }
  if (StringHasNoText (check_str)) return;
  
  feature_list = (ValNodePtr PNTR) userdata;
  
  for (k = 0; k < num_suspect_phrases; k++)
  {
    if (StringISearch(check_str, suspect_phrases[k]) != NULL)
    {
      ValNodeAddPointer (feature_list, OBJ_SEQFEAT, sfp);
      break;
    }
  }  
}

extern void FindSuspectPhrases (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr     feature_list = NULL, vnp;
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &feature_list, FindSuspectPhrasesCallback);
  }
  
  dip = SuspectPhrase (DISC_SUSPECT_PHRASES, "fragment or frameshift", "cds comments or protein description", feature_list);
  if (dip != NULL)
  {
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }

}


static void FindUnknownProteinsWithECNumbersCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ProtRefPtr prp;
  ValNodePtr PNTR feature_list;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT || sfp->data.value.ptrvalue == NULL || userdata == NULL) 
  {
    return;
  }
  
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp->name == NULL || prp->ec == NULL) return;

  if (StringISearch (prp->name->data.ptrvalue, "hypothetical protein") != NULL
      || StringISearch (prp->name->data.ptrvalue, "unknown protein") != NULL) 
  {
    feature_list = (ValNodePtr PNTR) userdata;
    ValNodeAddPointer (feature_list, OBJ_SEQFEAT, sfp);
  }
}


extern void FindUnknownProteinsWithECNumbers (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr         feature_list = NULL, vnp;
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;
  
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitFeaturesInSep (vnp->data.ptrvalue, &feature_list, FindUnknownProteinsWithECNumbersCallback);
  }
  
  if (feature_list != NULL) {
    dip = NewClickableItem (DISC_EC_NUMBER_ON_HYPOTHETICAL_PROTEIN, "%d protein features have an EC number and a protein name of 'unknown protein' or 'hypothetical protein'", feature_list);
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }  
}

static ClickableItemPtr InconsistentSourceDefline (SeqDescrPtr biop_sdp, SeqDescrPtr title_sdp)
{
  ClickableItemPtr dip = NULL;
  CharPtr            bad_fmt = "Organism description not found in definition line: %s.";
  BioSourcePtr       biop;
  CharPtr            desc = NULL;

  if (biop_sdp == NULL || title_sdp == NULL)
  {
    return NULL;
  }
  
  biop = (BioSourcePtr) biop_sdp->data.ptrvalue;
  if (biop != NULL && biop->org != NULL && !StringHasNoText (biop->org->taxname))
  {
    desc = biop->org->taxname;
  }
  else
  {
    desc = title_sdp->data.ptrvalue;
  }
  if (StringHasNoText (desc)) {
    return NULL;
  }
  
  dip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  if (dip != NULL)
  {
    dip->clickable_item_type = DISC_INCONSISTENT_BIOSRC_DEFLINE;
    dip->description = (CharPtr)MemNew (StringLen (bad_fmt) + StringLen (desc));
    sprintf (dip->description, bad_fmt, desc);
    dip->callback_func = NULL;
    dip->datafree_func = NULL;
    dip->callback_data = NULL;
    dip->item_list = NULL;
    ValNodeAddPointer (&(dip->item_list), OBJ_SEQDESC, biop_sdp);
    ValNodeAddPointer (&(dip->item_list), OBJ_SEQDESC, title_sdp);
  }      
  return dip;
}


static void FindInconsistentSourceAndDeflineCallback (BioseqPtr bsp, Pointer userdata)
{
  ClickableItemPtr dip;
  ValNodePtr PNTR discrepancy_list;
  SeqDescrPtr        biop_sdp, title_sdp;
  SeqMgrDescContext  context;
  BioSourcePtr       biop;
  
  discrepancy_list = (ValNodePtr PNTR) userdata;
  if (bsp == NULL || discrepancy_list == NULL) return;
  
  biop_sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_source, &context);
  if (biop_sdp == NULL || biop_sdp->data.ptrvalue == NULL)
  {
    return;
  }
  biop = (BioSourcePtr) biop_sdp->data.ptrvalue;
  if (biop->org == NULL)
  {
    return;
  }
  if (StringHasNoText (biop->org->taxname)) 
  {
    return;
  }
  
  title_sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_title, &context);
  if (title_sdp == NULL) return;
  
  if (StringStr (title_sdp->data.ptrvalue, biop->org->taxname) == NULL)
  {
    dip = InconsistentSourceDefline (biop_sdp, title_sdp);
    if (dip != NULL)
    {
      ValNodeAddPointer (discrepancy_list, 0, dip);
    }
  }
}


extern void FindInconsistentSourceAndDefline (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{  
  ValNodePtr disc_pairs = NULL, vnp;
  CharPtr    bad_fmt = "%d sources do not match definition lines.";
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &disc_pairs, FindInconsistentSourceAndDeflineCallback);
  }

  if (disc_pairs == NULL) 
  {
    return;
  }
  else if (disc_pairs->next == NULL)
  {
    ValNodeLink (discrepancy_list, disc_pairs);
  }
  else
  {
    dip = NewClickableItem (DISC_INCONSISTENT_BIOSRC_DEFLINE, bad_fmt, disc_pairs);
    dip->item_list = NULL;
    dip->subcategories = disc_pairs;
    
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}


static void FindParticalCDSsInCompleteSequencesCallback (BioseqPtr bsp, Pointer userdata)
{
  ValNodePtr PNTR    cds_list;
  SeqDescrPtr        molinfo_sdp;
  SeqMgrDescContext  context;
  SeqFeatPtr         cds;
  SeqMgrFeatContext  fcontext;
  MolInfoPtr         mip;
  Boolean            partial5, partial3;
  
  cds_list = (ValNodePtr PNTR) userdata;
  if (bsp == NULL || cds_list == NULL) return;
  
  molinfo_sdp = SeqMgrGetNextDescriptor(bsp, NULL, Seq_descr_molinfo, &context);
  if (molinfo_sdp == NULL || molinfo_sdp->data.ptrvalue == NULL)
  {
    return;
  }
  mip = (MolInfoPtr) molinfo_sdp->data.ptrvalue;
  if (mip->completeness != 1)
  {
    return;
  }
  
  cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
  while (cds != NULL) {
      CheckSeqLocForPartial (cds->location, &partial5, &partial3);
      if (cds->partial || partial5 || partial3) {
          ValNodeAddPointer (cds_list, OBJ_SEQFEAT, cds);
      }
      cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, FEATDEF_CDS, &fcontext);
  }
}


extern void FindParticalCDSsInCompleteSequences (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{  
  ValNodePtr cds_list = NULL, vnp;
  CharPtr    bad_fmt = "%d partial CDSs in complete sequences.";
  ClickableItemPtr dip;
  
  if (discrepancy_list == NULL) return;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    VisitBioseqsInSep (vnp->data.ptrvalue, &cds_list, FindParticalCDSsInCompleteSequencesCallback);
  }

  if (cds_list == NULL) 
  {
    return;
  }
  else
  {
    dip = NewClickableItem (DISC_PARTIAL_CDS_IN_COMPLETE_SEQUENCE, bad_fmt, cds_list);
    dip->subcategories = NULL;
    
    ValNodeAddPointer (discrepancy_list, 0, dip);
  }
}

static Boolean RnaRefMatch (RnaRefPtr rrp1, RnaRefPtr rrp2)
{
  tRNAPtr tp1, tp2;
  Boolean rval = FALSE;

  if (rrp1 == NULL && rrp2 == NULL) {
    rval = TRUE;
  } else if (rrp1 == NULL || rrp2 == NULL) {
    rval = FALSE;
  } else if (rrp1->type != rrp2->type) {
    rval = FALSE;
  } else if (rrp1->ext.choice != rrp2->ext.choice) {
    return FALSE;
  } else {
    switch (rrp1->ext.choice) {
      case 0:
        rval = TRUE;
        break;
      case 1:
        if (StringCmp (rrp1->ext.value.ptrvalue, rrp2->ext.value.ptrvalue) == 0) {
          rval = TRUE;
        } else {
          rval = FALSE;
        }
        break;
      case 2:
        tp1 = rrp1->ext.value.ptrvalue;
        tp2 = rrp2->ext.value.ptrvalue;
        if (tp1 == NULL && tp2 == NULL) {
          rval = TRUE;
        } else if (tp1 == NULL || tp2 == NULL) {
          rval = FALSE;
        } else if (tp1->aa == tp2->aa) {
          rval = TRUE;
        } else {
          rval = FALSE;
        }
        break;
      default:
        rval = FALSE;
        break;
    }
  }
  return rval;
}

static void AddRNAMatch (SeqFeatPtr sfp, ValNodePtr PNTR puniq_list)
{
  ValNodePtr vnp, uniq_vnp;
  SeqFeatPtr sfp_match;
  RnaRefPtr rrp_match, rrp_find;
  Boolean   found_match = FALSE;

  if (sfp == NULL || sfp->data.value.ptrvalue == NULL || puniq_list == NULL) return;
  rrp_find = (RnaRefPtr) sfp->data.value.ptrvalue;

  uniq_vnp = *puniq_list;

  if (uniq_vnp == NULL) {
    vnp = ValNodeNew(NULL);
    vnp->choice = OBJ_SEQFEAT;
    vnp->data.ptrvalue = sfp;
    vnp->next = NULL;
    ValNodeAddPointer (puniq_list, 0, vnp);
    found_match = TRUE;
  }
  while (uniq_vnp != NULL && !found_match) {
    vnp = uniq_vnp->data.ptrvalue;
    if (vnp == NULL) {
      /* fill in empty list */
      ValNodeAddPointer (&vnp, OBJ_SEQFEAT, sfp);
      uniq_vnp->data.ptrvalue = vnp;
      found_match = TRUE;
    } else {
      sfp_match = vnp->data.ptrvalue;
      if (sfp_match != NULL && sfp_match->data.choice == SEQFEAT_RNA && sfp_match->data.value.ptrvalue != NULL) {
        rrp_match = sfp_match->data.value.ptrvalue;
        if (RnaRefMatch(rrp_match, rrp_find)) {
          ValNodeAddPointer (&vnp, OBJ_SEQFEAT, sfp);
          found_match = TRUE;
          /* set flag so we know this list has duplicates */
          uniq_vnp->choice = 1;
        } 
      }
      if (!found_match) {
        if (uniq_vnp->next == NULL) {
          /* add to end of list */
          uniq_vnp->next = ValNodeNew(NULL);
          uniq_vnp->next->next = NULL;
          uniq_vnp->next->choice = 0;
          vnp = ValNodeNew(NULL);
          vnp->choice = OBJ_SEQFEAT;
          vnp->data.ptrvalue = sfp;
          vnp->next = NULL;
          uniq_vnp->next->data.ptrvalue = vnp;
          found_match = TRUE;
        } else {
          uniq_vnp = uniq_vnp->next;
        }
      }
    }
  }
}
  

static void FindDupRNAsInList (ValNodePtr rna_list, ValNodePtr PNTR discrepancy_list, CharPtr label, CharPtr id_str)
{
  ValNodePtr vnp, uniq_list = NULL;
  ValNodePtr dup_list = NULL;
  CharPtr          dup_fmt = "%d %s features on %s have the same name (%s)";
  ClickableItemPtr cip;
  SeqFeatPtr       sfp;
  SeqMgrFeatContext fcontext;

  for (vnp = rna_list; vnp != NULL; vnp = vnp->next) {
    AddRNAMatch (vnp->data.ptrvalue, &uniq_list);
  }

  dup_list = ValNodeExtractList (&uniq_list, 1);

  for (vnp = uniq_list; vnp != NULL; vnp = vnp->next) {
    uniq_list->data.ptrvalue = ValNodeFree (uniq_list->data.ptrvalue);
  }
  uniq_list = ValNodeFree (uniq_list);

  for (vnp = dup_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->item_list = vnp->data.ptrvalue;
    sfp = cip->item_list->data.ptrvalue;
    sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (dup_fmt) + StringLen (label) + StringLen (id_str) + StringLen (fcontext.label) + 15));
    sprintf (cip->description, dup_fmt, ValNodeLen (cip->item_list), label, id_str, fcontext.label);
    if (sfp->idx.subtype == FEATDEF_tRNA) {
      cip->clickable_item_type = DISC_DUP_TRNA;
    } else {
      cip->clickable_item_type = DISC_DUP_RRNA;
    }
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
  dup_list = ValNodeFree (dup_list);

}

typedef struct desiredaa {
  Char    short_symbol;
  CharPtr long_symbol;
  Int4    num_expected;
} DesiredAAData, PNTR DesiredAAPtr;

static DesiredAAData desired_aaList [] = {
{'A', "Ala", 1 },
{'B', "Asx", 0 },
{'C', "Cys", 1 },
{'D', "Asp", 1 },
{'E', "Glu", 1 },
{'F', "Phe", 1 },
{'G', "Gly", 1 },
{'H', "His", 1 },
{'I', "Ile", 1 },
{'J', "Xle", 0 },
{'K', "Lys", 1 },
{'L', "Leu", 2 },
{'M', "Met", 1 },
{'N', "Asn", 1 },
{'P', "Pro", 1 },
{'Q', "Gln", 1 },
{'R', "Arg", 1 },
{'S', "Ser", 2 },
{'T', "Thr", 1 },
{'V', "Val", 1 },
{'W', "Trp", 1 },
{'X', "Xxx", 0 },
{'Y', "Tyr", 1 },
{'Z', "Glx", 0 },
{'U', "Sec", 0 },
{'O', "Pyl", 0 },
{'*', "Ter", 0 }
};

static void AddMissingtRNADiscrepancy (CharPtr str, ValNodePtr PNTR discrepancy_list, CharPtr id_str, BioseqPtr bsp)
{
  ClickableItemPtr cip;
  CharPtr          desc_fmt = "Sequence %s is missing trna-%s";

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));

  cip->clickable_item_type = DISC_COUNT_TRNA;
  ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (desc_fmt) + StringLen (id_str) + StringLen (str)));
  sprintf (cip->description, desc_fmt, id_str, str);
  ValNodeAddPointer (discrepancy_list, 0, cip);
}

static void 
AddExtratRNADiscrepancy 
(CharPtr         str, 
 Int4            num, 
 ValNodePtr PNTR discrepancy_list, 
 CharPtr         id_str, 
 BioseqPtr       bsp,
 ValNodePtr      rna_list)
{
  ClickableItemPtr cip;
  SeqMgrFeatContext fcontext;
  CharPtr          desc_fmt = "Sequence %s has %d trna-%s features";
  ValNodePtr       vnp;
  SeqFeatPtr       sfp;

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));

  cip->clickable_item_type = DISC_COUNT_TRNA;
  ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
  for (vnp = rna_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    if (StringSearch (fcontext.label, str) != NULL) {
      ValNodeAddPointer (&(cip->item_list), OBJ_SEQFEAT, sfp);
    }
  }
  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (desc_fmt) + StringLen (id_str) + StringLen (str) + 15));
  sprintf (cip->description, desc_fmt, id_str, num, str);
  ValNodeAddPointer (discrepancy_list, 0, cip);
}

static void FindMissingRNAsInList (ValNodePtr rna_list, ValNodePtr PNTR discrepancy_list, CharPtr id_str, BioseqPtr bsp)
{
  ValNodePtr       vnp;
  SeqFeatPtr       sfp;
  SeqMgrFeatContext fcontext;
  Uint1            num;
  Int4Ptr          num_present;
  Uint1            i;

  num = sizeof (desired_aaList) / sizeof (DesiredAAData);

  num_present = (Int4Ptr) MemNew (sizeof (Int4) * num);
  MemSet (num_present, 0, sizeof (Int4) * num);

  for (vnp = rna_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
    for (i = 0; i < num; i++) {
      if (StringSearch (fcontext.label, desired_aaList[i].long_symbol) != NULL) {
        num_present[i] ++;
        break;
      }
    }
  }
  for (i = 0; i < num; i++) {
    if (num_present[i] < desired_aaList[i].num_expected) {
      AddMissingtRNADiscrepancy (desired_aaList[i].long_symbol, discrepancy_list, id_str, bsp);
    } else if (num_present[i] > desired_aaList[i].num_expected) {
      AddExtratRNADiscrepancy (desired_aaList[i].long_symbol, num_present[i], discrepancy_list, id_str, bsp, rna_list);
    }
  }
}

typedef struct featcount {
  Uint1      featdeftype;
  ValNodePtr discrepancy_list;
} FeatCountData, PNTR FeatCountPtr;

static void RNACountFeaturesBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  ValNodePtr         feat_list = NULL;
  BioSourcePtr       biop;
  Boolean            run_test = FALSE;
  FeatCountPtr       fcp;
  CharPtr            count_fmt = "%d %s features found on %s";
  CharPtr            label;
  ClickableItemPtr   cip;
  Char        id_str[45];

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }

  fcp = (FeatCountPtr) userdata;

  /* look for Bioseq with organelle */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  while (sdp != NULL && !run_test) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL 
        && (biop->genome == GENOME_plastid
            || biop->genome == GENOME_mitochondrion
            || biop->genome == GENOME_chloroplast)) {
      run_test = TRUE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
  }
  if (!run_test) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, fcp->featdeftype, &fcontext);
  while (sfp != NULL) {
    ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, fcp->featdeftype, &fcontext);
  }

  if (feat_list != NULL) {
    label = (CharPtr) FeatDefTypeLabel(feat_list->data.ptrvalue);
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, sizeof (id_str) - 1);
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (count_fmt) + StringLen (label) + StringLen (id_str) + 15));
    sprintf (cip->description, count_fmt, ValNodeLen (feat_list), label, id_str);
    cip->item_list = feat_list;
    if (fcp->featdeftype == FEATDEF_tRNA) {
      cip->clickable_item_type = DISC_COUNT_TRNA;
    } else {
      cip->clickable_item_type = DISC_COUNT_RRNA;
    }
    ValNodeAddPointer (&(fcp->discrepancy_list), 0, cip);
    if (fcp->featdeftype == FEATDEF_tRNA) {
      FindMissingRNAsInList (feat_list, &(fcp->discrepancy_list), id_str, bsp);
    } else {
      FindDupRNAsInList (feat_list, &(fcp->discrepancy_list), label, id_str);
    }
  }
}


static CharPtr GetRNACompareString (ValNodePtr vnp)
{
  ClickableItemPtr cip;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;
  cip = (ClickableItemPtr) vnp->data.ptrvalue;
  return cip->description;
}

static BioseqPtr GetRNATestBioseq (ValNodePtr vp)
{
  ClickableItemPtr cip;
  ValNodePtr       vnp;
  BioseqPtr        bsp = NULL;
  SeqFeatPtr       sfp;

  if (vp == NULL || vp->data.ptrvalue == NULL) return NULL;
  cip = (ClickableItemPtr) vp->data.ptrvalue;
  for (vnp = cip->item_list; vnp != NULL && bsp == NULL; vnp = vnp->next) {
    if (vnp->data.ptrvalue == NULL) continue;
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
    }
  }
  return bsp;
}

static void AddRNANumList (ValNodePtr PNTR discrepancy_list, ValNodePtr list_start)
{
  ClickableItemPtr cip;
  CharPtr          cp;
  CharPtr          desc_fmt = "%d sequences have ";
  CharPtr          desc_str;
  Int4             copy_len, orig_len;
  ValNodePtr       vnp;
  BioseqPtr        bsp;

  if (discrepancy_list == NULL || list_start == NULL) return;
  desc_str = GetRNACompareString (list_start);
  cp = StringSearch (desc_str, " found on");
  if (cp != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->subcategories = list_start;
    copy_len = cp - desc_str;
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (desc_fmt) + 15 + copy_len));
    sprintf (cip->description, desc_fmt, ValNodeLen (list_start));
    orig_len = StringLen (cip->description);
    StringNCat (cip->description, desc_str, copy_len);
    cip->description [orig_len + copy_len] = 0;

    for (vnp = list_start; vnp != NULL; vnp = vnp->next) {
      bsp = GetRNATestBioseq (vnp);
      if (bsp != NULL) {
        ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, bsp);
      }
    }
    ValNodeAddPointer (discrepancy_list, 0, cip);
  } else {
    ValNodeLink (discrepancy_list, list_start);
  }
}

extern int LIBCALLBACK SortVnpByClickableItemDescription (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = GetRNACompareString (vnp1);
      str2 = GetRNACompareString (vnp2);
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
      }
    }
  }
  return 0;
}

static void RNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list, Uint1 featdeftype)
{
  ValNodePtr       vnp, list_start, list_last;
  FeatCountData    fcd;
  SeqEntryPtr      sep;
  CharPtr          cp, compare1;
  Int4             compare_len;

  fcd.featdeftype = featdeftype;
  fcd.discrepancy_list = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &fcd, RNACountFeaturesBioseqCallback);
  }

  /* count how many Bioseqs have different numbers of features */
  fcd.discrepancy_list = ValNodeSort (fcd.discrepancy_list, SortVnpByClickableItemDescription);

  while (fcd.discrepancy_list != NULL) {
    list_start = fcd.discrepancy_list;
    compare1 = GetRNACompareString (list_start);
    cp = StringSearch (compare1, "found on");
    if (cp == NULL) {
      fcd.discrepancy_list = fcd.discrepancy_list->next;
      list_start->next = NULL;
      AddRNANumList (discrepancy_list, list_start);
    } else {
      compare_len = cp - compare1;
      list_last = list_start;
      vnp = list_start->next;
      while (vnp != NULL && StringNCmp (compare1, GetRNACompareString (vnp), compare_len) == 0) {
        list_last = vnp;
        vnp = vnp->next;
      }
      
      list_last->next = NULL;
      fcd.discrepancy_list = vnp;
      AddRNANumList (discrepancy_list, list_start);
    }
  }
}

extern void tRNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  RNACountFeaturesAndFindDups (discrepancy_list, sep_list, FEATDEF_tRNA);
}

extern void rRNACountFeaturesAndFindDups (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  RNACountFeaturesAndFindDups (discrepancy_list, sep_list, FEATDEF_rRNA);
}

static void CountShorttRNA (SeqFeatPtr sfp, Pointer data)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_tRNA || data == NULL) return;

  if (SeqLocLen (sfp->location) < 50) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}

static void CountLongtRNA (SeqFeatPtr sfp, Pointer data)
{
  if (sfp == NULL || sfp->idx.subtype != FEATDEF_tRNA || data == NULL) return;

  if (SeqLocLen (sfp->location) > 90) {
    ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
  }
}

extern void tRNAFindBadLength (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr  vnp;
  SeqEntryPtr sep;
  ValNodePtr  too_short = NULL, too_long = NULL;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &too_short, CountShorttRNA);
  }
  if (too_short != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BADLEN_TRNA, "%d tRNAs are too short", too_short));
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &too_long, CountLongtRNA);
  }
  if (too_short != NULL) {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (DISC_BADLEN_TRNA, "%d tRNAs are too long", too_long));
  }

}

static void FindRNAsWithoutProductsCallback (SeqFeatPtr sfp, Pointer data)
{
  RnaRefPtr rp;

  if (sfp != NULL && sfp->data.choice == SEQFEAT_RNA && sfp->data.value.ptrvalue != NULL) {
    rp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rp->ext.choice == 0 || (rp->ext.choice == 1 && StringHasNoText (rp->ext.value.ptrvalue))) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, OBJ_SEQFEAT, sfp);
    }
  }
}

extern void FindRNAsWithoutProducts (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;
  ValNodePtr       rna_list = NULL;
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitFeaturesInSep (sep, &rna_list, FindRNAsWithoutProductsCallback);
  }
  if (rna_list != NULL) {
   cip = NewClickableItem (DISC_RNA_NO_PRODUCT, "%d RNA features have no product", rna_list);
   ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void tRNASameStrandBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqMgrFeatContext fcontext;
  SeqMgrDescContext dcontext;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  ValNodePtr         feat_list = NULL;
  BioSourcePtr       biop;
  Boolean            run_test = FALSE;
  ValNodePtr PNTR    discrepancy_list;
  ClickableItemPtr   cip;

  Uint1              strand, this_strand;
  Boolean            mixed_strand = FALSE;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }

  discrepancy_list = (ValNodePtr PNTR) userdata;

  /* look for Bioseq with organelle */
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  while (sdp != NULL && !run_test) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL 
        && (biop->genome == GENOME_plastid
            || biop->genome == GENOME_mitochondrion
            || biop->genome == GENOME_chloroplast)) {
      run_test = TRUE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
  }
  if (!run_test) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, FEATDEF_tRNA, &fcontext);
  while (sfp != NULL && !mixed_strand) {
    if (feat_list == NULL) {
      strand = SeqLocStrand (sfp->location);
    } else {
      this_strand = SeqLocStrand (sfp->location);
      if ((strand == Seq_strand_minus && this_strand != Seq_strand_minus)
          || (strand != Seq_strand_minus && this_strand == Seq_strand_minus)) {
        mixed_strand = TRUE;
      }
    }
    ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, FEATDEF_tRNA, &fcontext);
  }

  if (mixed_strand) {
    feat_list = ValNodeFree (feat_list);
  } else if (feat_list != NULL) {
    if (strand == Seq_strand_minus) {
      cip = NewClickableItem (DISC_STRAND_TRNA, "%d tRNAs on minus strand", feat_list);
    } else {
      cip = NewClickableItem (DISC_STRAND_TRNA, "%d tRNAs on plus strand", feat_list);
    }
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


extern void FindtRNAsOnSameStrand (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, discrepancy_list, tRNASameStrandBioseqCallback);
  }
}


typedef struct translnote {
  ValNodePtr transl_no_note;
  ValNodePtr note_no_transl;
  ValNodePtr transl_too_long;
} TranslNoteData, PNTR TranslNotePtr;

static Boolean CodingRegionHasTranslExcept (SeqFeatPtr sfp)
{
  CodeBreakPtr cbp;
  Int4         len;
  CdRegionPtr  crp;
  SeqLocPtr    slp;
  Int4         codon_start, codon_stop, pos, codon_length;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue) == NULL
      || crp->code_break == NULL)
  {
      return FALSE;
  }

  len = SeqLocLen (sfp->location);
  if (crp->frame > 1) 
  {
      len -= crp->frame - 1;
  }
  if (len % 3 == 0) 
  {
      return FALSE;
  }
  for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next)
  {
    if (cbp->aa.choice != 1 || cbp->aa.value.intvalue != 42) {
      continue;
    }
    codon_start = INT4_MAX;
    codon_stop = -10;
    slp = NULL;
    while ((slp = SeqLocFindNext (cbp->loc, slp)) != NULL) {
      pos = GetOffsetInLoc (slp, sfp->location, SEQLOC_START);
      if (pos < codon_start)
      {
        codon_start = pos;
        pos = GetOffsetInLoc (slp, sfp->location, SEQLOC_STOP);
        if (pos > codon_stop)
        {
          codon_stop = pos;
        }
        codon_length = codon_stop - codon_start;      /* codon length */
        if (codon_length >= 0 && codon_length <= 1 && codon_stop == len - 1)
        {                       /*  a codon */
          /* allowing a partial codon at the end */
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}

static Boolean TranslTooLong (SeqFeatPtr sfp)
{
  CodeBreakPtr cbp;
  CdRegionPtr  crp;
  
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION
      || (crp = (CdRegionPtr)sfp->data.value.ptrvalue) == NULL
      || crp->code_break == NULL)
  {
      return FALSE;
  }

  for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next)
  {
    if (cbp->aa.choice == 1 
        && cbp->aa.value.intvalue == 42
        && SeqLocLen (cbp->loc) > 3) {
      return TRUE;
    }
  }
  return FALSE;
}

static void FindTranslNoNote (SeqFeatPtr sfp, Pointer userdata)
{
  TranslNotePtr tnp;
  CharPtr       note_txt = "TAA stop codon is completed by the addition of 3' A residues to the mRNA";

  if (sfp != NULL && userdata != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
    tnp = (TranslNotePtr) userdata;
    if (CodingRegionHasTranslExcept (sfp)) {
      if (StringStr (sfp->comment, note_txt) == NULL) {
        ValNodeAddPointer (&(tnp->transl_no_note), OBJ_SEQFEAT, sfp);
      }
    } else if (StringStr (sfp->comment, note_txt) != NULL) {
      ValNodeAddPointer (&(tnp->note_no_transl), OBJ_SEQFEAT, sfp);
    }
    if (TranslTooLong(sfp)) {
      ValNodeAddPointer (&(tnp->transl_too_long), OBJ_SEQFEAT, sfp);
    }
  }
}

extern void FindTranslExceptNotes (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  TranslNoteData   tnd;
  SeqEntryPtr      sep;
  ClickableItemPtr cip;
  CharPtr          transl_no_note_fmt = "%d features have a translation exception but no note";
  CharPtr          note_no_transl_fmt = "%d features have a note but not translation exception";
  CharPtr          transl_too_long_fmt = "%d features have translation exceptions longer than 3 bp";

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    tnd.transl_no_note = NULL;
    tnd.note_no_transl = NULL;
    tnd.transl_too_long = NULL;
    VisitFeaturesInSep (sep, &tnd, FindTranslNoNote);
    if (tnd.transl_no_note != NULL) {
      cip = NewClickableItem (DISC_TRANSL_NO_NOTE, transl_no_note_fmt, tnd.transl_no_note);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
    if (tnd.note_no_transl != NULL) {
      cip = NewClickableItem (DISC_NOTE_NO_TRANSL, note_no_transl_fmt, tnd.note_no_transl);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
    if (tnd.transl_too_long != NULL) {
      cip = NewClickableItem (DISC_TRANSL_TOO_LONG, transl_too_long_fmt, tnd.note_no_transl);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
}

static Boolean GetOverlappingTRNAs (BioseqPtr bsp, SeqLocPtr slp, Int4 loc_right, ValNodePtr PNTR list)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext  context;
  Boolean            found_any = FALSE;

  if (bsp == NULL || slp == NULL || list == NULL) return FALSE;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, FEATDEF_tRNA, &context);
       sfp != NULL && context.left <= loc_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_RNA, FEATDEF_tRNA, &context))
  {
    if (SeqLocStrand (sfp->location) == SeqLocStrand (slp) && SeqLocCompare (sfp->location, slp) != SLC_NO_MATCH) {
      ValNodeAddPointer (list, OBJ_SEQFEAT, sfp);
      found_any = TRUE;
    }
  }
  return found_any;
}

static void FindCDSOverlappingtRNAsBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  context;
  ValNodePtr         subcategories = NULL;
  ValNodePtr PNTR    discrepancy_list;
  ValNodePtr         item_list, all_item_list = NULL;
  ValNodePtr         trna_list = NULL;
  ClickableItemPtr   cip;
  CharPtr            list_fmt = "%d coding regions have overlapping tRNAs";
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  discrepancy_list = (ValNodePtr PNTR) userdata;
  
  
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, FEATDEF_CDS, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, FEATDEF_CDS, &context))
  {
    
    item_list = NULL;
    trna_list = NULL;
    if (GetOverlappingTRNAs (bsp, sfp->location, context.right, &trna_list)) {
      ValNodeAddPointer (&item_list, OBJ_SEQFEAT, sfp);
      ValNodeLink (&item_list, trna_list);
      ValNodeLink (&all_item_list, ValNodePointerDup(item_list));
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      MemSet (cip, 0, sizeof (ClickableItemData));
      cip->item_list = item_list;
      cip->description = StringSave ("Coding region overlaps tRNAs");
      cip->clickable_item_type = DISC_CDS_OVERLAP_TRNA;
      ValNodeAddPointer (&subcategories, 0, cip);
    }
  }  
  if (subcategories != NULL) {
    cip = (ClickableItemPtr) MemNew (sizeof(ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->clickable_item_type = DISC_CDS_OVERLAP_TRNA;
    cip->item_list = all_item_list;
    cip->subcategories = subcategories;
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (list_fmt) + 15));
    sprintf (cip->description, list_fmt, ValNodeLen (subcategories));
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


extern void FindCDSOverlappingtRNAs (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp, this_list = NULL;
  SeqEntryPtr      sep;
  ClickableItemPtr cip;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, &this_list, FindCDSOverlappingtRNAsBioseqCallback);
  }
  if (this_list != NULL) {
    cip = NewClickableItem (DISC_CDS_OVERLAP_TRNA, "%d Bioseqs have coding regions that overlap tRNAs", this_list);
    cip->subcategories = this_list;
    cip->item_list = ItemListFromSubcategories (this_list);
    ValNodeAddPointer (discrepancy_list, 0, cip);
  }
}


static void FindFeaturesOverlappingSrcFeaturesBioseqCallback (BioseqPtr bsp, Pointer data)
{
  ClickableItemPtr  cip;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  ValNodePtr        this_list, src_vnp;

  if (bsp == NULL || data == NULL) return;


  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_BIOSRC, &fcontext);
  while (sfp != NULL)
  {
    this_list = ListFeaturesOverlappingLocation (bsp, sfp->location, 0, 0);
    if (this_list != NULL)
    {
      cip = NewClickableItem (DISC_FEAT_OVERLAP_SRCFEAT, "%d features overlap a source feature", this_list);
      /* insert source feature at beginning of item list */
      src_vnp = ValNodeNew (NULL);     
      src_vnp->choice = OBJ_SEQFEAT;
      src_vnp->data.ptrvalue = sfp;
      src_vnp->next = cip->item_list;
      cip->item_list = src_vnp;
      ValNodeAddPointer ((ValNodePtr PNTR) data, 0, cip);
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_BIOSRC, &fcontext);
  }
}


extern void FindFeaturesOverlappingSrcFeatures (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    VisitBioseqsInSep (sep, discrepancy_list, FindFeaturesOverlappingSrcFeaturesBioseqCallback);
  }
}


extern CharPtr discReportDuplicateProteinIDFmt = "%d coding regions have non-unique protein IDs";
extern CharPtr discReportOneDuplicateProteinIDFmt = "%d coding regions have protein ID %s";
extern CharPtr discReportMissingProteinIDFmt = "%d coding regions have missing protein IDs";
extern CharPtr discReportDuplicateTranscriptIdFmt = "%d mRNAs have non-unique transcript IDs";
extern CharPtr discReportOneDuplicateTranscriptIdFmt = "%d mRNAs have non-unique transcript ID %s";
extern CharPtr discReportMissingTranscriptIDFmt = "%d mRNAs have missing transcript IDs";


/* look for duplicate protein IDs and duplicate transcript IDs */
/* every coding region should have a protein ID and a transcript ID */
/* RNA should have a transcript ID to match. */
static void CheckGenProdSetBioseq (BioseqPtr bsp, GenProdSetDiscrepancyListsPtr lists)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  SeqIdPtr          sip;
  Char              buf [96];

  if (bsp == NULL || !ISA_na (bsp->mol) || lists == NULL) {
    return;
  }

  /* look for missing protein IDs and duplicate protein IDs on coding regions */
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_CDREGION, 0, &fcontext)) {
    if (sfp->product == NULL) {
      if (!sfp->pseudo) {
        ValNodeAddPointer (&(lists->missing_protein_id), 0,
                           GlobalDiscrepancyNew (NULL, OBJ_SEQFEAT, sfp));
      }
    } else {
      sip = SeqLocId (sfp->product);
      SeqIdWrite (sip, buf, PRINTID_REPORT, sizeof (buf) - 1);
      ValNodeAddPointer (&(lists->cds_product_list), 0,
                         GlobalDiscrepancyNew (buf, OBJ_SEQFEAT, sfp));
    }
  }

  /* look for missing transcript IDs and duplicate transcript IDs on mRNAs */
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_mRNA, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_mRNA, &fcontext)) {
    if (sfp->product == NULL) {
      ValNodeAddPointer (&(lists->missing_mrna_product), 0, 
                         GlobalDiscrepancyNew (NULL, OBJ_SEQFEAT, sfp));
    } else {
      sip = SeqLocId (sfp->product);
      SeqIdWrite (sip, buf, PRINTID_REPORT, sizeof (buf) - 1);
      ValNodeAddPointer (&(lists->mrna_product_list), 0,
                         GlobalDiscrepancyNew (buf, OBJ_SEQFEAT, sfp));
    }
  }
}


extern void CheckGenProdSetsInSeqEntry (SeqEntryPtr sep, GenProdSetDiscrepancyListsPtr lists)
{
  BioseqSetPtr bssp;

  if (sep == NULL || !IS_Bioseq_set (sep) || sep->data.ptrvalue == NULL || lists == NULL) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
    if (IS_Bioseq (bssp->seq_set)) {
      CheckGenProdSetBioseq(bssp->seq_set->data.ptrvalue, lists);
    }
  } else {
    sep = bssp->seq_set;
    while (sep != NULL) {
      CheckGenProdSetsInSeqEntry (sep, lists);
      sep = sep->next;
    }
  }
}


static void CheckListForGenProdSets (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp, disc_list = NULL;
  ClickableItemPtr cip;
  GenProdSetDiscrepancyListsData lists;

  MemSet (&lists, 0, sizeof (GenProdSetDiscrepancyListsData));
  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    CheckGenProdSetsInSeqEntry (vnp->data.ptrvalue, &lists);
  }

  if (lists.missing_protein_id != NULL) {
    cip = ReportMissingFields (lists.missing_protein_id, discReportMissingProteinIDFmt, DISC_MISSING_GENPRODSET_PROTEIN);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    lists.missing_protein_id = FreeGlobalDiscrepancyList (lists.missing_protein_id);
  }

  if (lists.cds_product_list != NULL) {
    lists.cds_product_list = ValNodeSort (lists.cds_product_list, SortVnpByGlobalDiscrepancyString);
    cip = ReportNonUniqueGlobalDiscrepancy (lists.cds_product_list, 
                                            discReportDuplicateProteinIDFmt,
                                            discReportOneDuplicateProteinIDFmt,
                                            DISC_DUP_GENPRODSET_PROTEIN,
                                            FALSE);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    lists.cds_product_list = FreeGlobalDiscrepancyList (lists.cds_product_list);
  }
  
  
  if (lists.missing_mrna_product != NULL) {
    cip = ReportMissingFields (lists.missing_mrna_product, discReportMissingTranscriptIDFmt, DISC_MISSING_GENPRODSET_TRANSCRIPT_ID);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    lists.missing_mrna_product = FreeGlobalDiscrepancyList (lists.missing_mrna_product);
  }


  if (lists.mrna_product_list != NULL) {
    lists.mrna_product_list = ValNodeSort (lists.mrna_product_list, SortVnpByGlobalDiscrepancyString);
    cip = ReportNonUniqueGlobalDiscrepancy (lists.mrna_product_list, 
                                            discReportDuplicateTranscriptIdFmt,
                                            discReportOneDuplicateTranscriptIdFmt,
                                            DISC_DUP_GENPRODSET_TRANSCRIPT_ID,
                                            FALSE);
    if (cip != NULL) {
      ValNodeAddPointer (&disc_list, 0, cip);
    }
    lists.mrna_product_list = FreeGlobalDiscrepancyList (lists.mrna_product_list);
  }
  
    


  if (disc_list != NULL) {
    if (disc_list->next == NULL) {
      ValNodeLink (discrepancy_list, disc_list);
    } else {
      cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
      cip->description = StringSave ("GenProdSet Errors");
      cip->subcategories = disc_list;
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
}


static void CountProteinsBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  if (bsp != NULL && ISA_aa (bsp->mol) && userdata != NULL) {
    ValNodeAddPointer ((ValNodePtr PNTR) userdata, OBJ_BIOSEQ, bsp);
  }
}

extern void CountProteins (ValNodePtr PNTR discrepancy_list, ValNodePtr sep_list)
{
  ValNodePtr       vnp;
  SeqEntryPtr      sep;
  ValNodePtr       proteins;
  ClickableItemPtr cip;
  CharPtr          prot_count_fmt = "%d protein sequences in record";

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next) {
    sep = vnp->data.ptrvalue;
    proteins = NULL;
    VisitBioseqsInSep (sep, &proteins, CountProteinsBioseqCallback);
    if (proteins != NULL) {
      cip = NewClickableItem (DISC_COUNT_PROTEINS, prot_count_fmt, proteins);
      ValNodeAddPointer (discrepancy_list, 0, cip);
    }
  }
}


static void 
RemoveUnwantedDiscrepancyItems 
(ValNodePtr PNTR      discrepancy_list,
 DiscrepancyConfigPtr dcp)
{
  ValNodePtr         vnp, prev = NULL, vnp_next;
  ClickableItemPtr dip;
  
  if (dcp == NULL || discrepancy_list == NULL || *discrepancy_list == NULL)
  {
    return;
  }
  
  for (vnp = *discrepancy_list; vnp != NULL; vnp = vnp_next)
  {
    vnp_next = vnp->next;
    dip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (dip == NULL || ! dcp->conf_list[dip->clickable_item_type])
    {
      if (prev == NULL)
      {
        *discrepancy_list = vnp_next;
      }
      else
      {
        prev->next = vnp_next;
      }
      vnp->next = NULL;
      vnp = FreeClickableList (vnp);
    }
    else
    {
      prev = vnp;
    }
  }
  
}


static void SetDiscrepancyLevels (ValNodePtr discrepancy_list, Int4 level)
{
  ClickableItemPtr dip;
  
  while (discrepancy_list != NULL)
  {
    dip = (ClickableItemPtr) discrepancy_list->data.ptrvalue;
    if (dip != NULL)
    {
      dip->level = level;
      SetDiscrepancyLevels (dip->subcategories, level + 1);
    }
    discrepancy_list = discrepancy_list->next;
  }
}



typedef struct discrepancyinfo 
{
  CharPtr                conf_name;
  CharPtr                setting_name;
  PerformDiscrepancyTest test_func;
} DiscrepancyInfoData, PNTR DiscrepancyInfoPtr;

static DiscrepancyInfoData discrepancy_info_list[] = 
{
  { "Missing Genes", "MISSING_GENES", AddMissingAndSuperfluousGeneDiscrepancies },
  { "Extra Genes", "EXTRA_GENES", AddMissingAndSuperfluousGeneDiscrepancies },
  { "Missing Locus Tags", "MISSING_LOCUS_TAGS", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Duplicate Locus Tags", "DUPLICATE_LOCUS_TAGS", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Bad Locus Tag Format", "BAD_LOCUS_TAG_FORMAT", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Inconsistent Locus Tag Prefix", "INCONSISTENT_LOCUS_TAG_PREFIX", AddDiscrepanciesForMissingOrNonUniqueGeneLocusTags },
  { "Nongene Locus Tag", "NON_GENE_LOCUS_TAG", AddDiscrepanciesForNonGeneLocusTags },
  { "Missing Protein ID", "MISSING_PROTEIN_ID", FindMissingProteinIDs },
  { "Inconsistent Protein ID", "INCONSISTENT_PROTEIN_ID", FindMissingProteinIDs },
  { "Feature Location Conflict", "FEATURE_LOCATION_CONFLICT", FindCDSmRNAGeneLocationDiscrepancies },
  { "Gene Product Conflict", "GENE_PRODUCT_CONFLICT", FindCDSGeneProductConflicts },
  { "Duplicate Gene Locus", "DUPLICATE_GENE_LOCUS", FindDuplicateGeneLocus },
  { "EC Number Note", "EC_NUMBER_NOTE", AddECNumberNoteDiscrepancies },
  { "Pseudo Mismatch", "PSEUDO_MISMATCH", FindPseudoDiscrepancies },
  { "Joined Features", "JOINED_FEATURES", AddJoinedFeatureDiscrepancies },
  { "Overlapping Genes", "OVERLAPPING_GENES", AddOverlappingGeneDiscrepancies },
  { "Overlapping CDS", "OVERLAPPING_CDS", AddOverlappingCodingRegionDiscrepancies },
  { "Contained CDS", "CONTAINED_CDS", AddContainedCodingRegionDiscrepancies },
  { "CDS RNA Overlap", "RNA_CDS_OVERLAP", AddRNACDSOverlapDiscrepancies },
  { "Short Contig", "SHORT_CONTIG", FindShortContigs },
  { "Inconsistent BioSource", "INCONSISTENT_BIOSOURCE", FindNonmatchingContigSources },
  { "Suspect Product Name", "SUSPECT_PRODUCT_NAMES", FindSuspectProductNames },
  { "Inconsistent Source And Definition Line", "INCONSISTENT_SOURCE_DEFLINE", FindInconsistentSourceAndDefline },
  { "Partial CDSs in Complete Sequences", "PARTIAL_CDS_COMPLETE_SEQUENCE", FindParticalCDSsInCompleteSequences },
  { "Hypothetical or Unknown Protein with EC Number", "EC_NUMBER_ON_UNKNOWN_PROTEIN", FindUnknownProteinsWithECNumbers },
  { "Find Missing Tax Lookups", "TAX_LOOKUP_MISSING", NULL } ,
  { "Find Tax Lookup Mismatches", "TAX_LOOKUP_MISMATCH", NULL },
  { "Find Short Sequences", "SHORT_SEQUENCES", FindShortSequences },
  { "Suspect Phrases", "SUSPECT_PHRASES", FindSuspectPhrases },
  { "Count tRNAs", "COUNT_TRNAS", tRNACountFeaturesAndFindDups },
  { "Find Duplicate tRNAs", "FIND_DUP_TRNAS", tRNACountFeaturesAndFindDups },
  { "Find short and long tRNAs", "FIND_BADLEN_TRNAS", tRNAFindBadLength},
  { "Find tRNAs on the same strand", "FIND_STRAND_TRNAS", FindtRNAsOnSameStrand},
  { "Count rRNAs", "COUNT_RRNAS", rRNACountFeaturesAndFindDups },
  { "Find Duplicate rRNAs", "FIND_DUP_RRNAS", rRNACountFeaturesAndFindDups },
  { "Find RNAs without Products", "RNA_NO_PRODUCT", FindRNAsWithoutProducts },
  { "Transl_except without Note", "TRANSL_NO_NOTE", FindTranslExceptNotes },
  { "Note without Transl_except", "NOTE_NO_TRANSL", FindTranslExceptNotes },
  { "Transl_except longer than 3", "TRANSL_TOO_LONG", FindTranslExceptNotes },
  { "CDS tRNA overlaps", "CDS_TRNA_OVERLAP", FindCDSOverlappingtRNAs },
  { "Count Proteins", "COUNT_PROTEINS", CountProteins },
  { "Features Intersecting Source Features", "DISC_FEAT_OVERLAP_SRCFEAT", FindFeaturesOverlappingSrcFeatures },
  { "CDS on GenProdSet without protein", "MISSING_GENPRODSET_PROTEIN", CheckListForGenProdSets},
  { "Multiple CDS on GenProdSet, same protein", "DUP_GENPRODSET_PROTEIN", CheckListForGenProdSets},
  { "mRNA on GenProdSet without transcript ID", "MISSING_GENPRODSET_TRANSCRIPT_ID", CheckListForGenProdSets},
  { "mRNA on GenProdSet with duplicate ID", "DISC_DUP_GENPRODSET_TRANSCRIPT_ID", CheckListForGenProdSets}

};

extern CharPtr GetDiscrepancyTestConfName (DiscrepancyType dtype) 
{
  return discrepancy_info_list[dtype].conf_name;
}

extern CharPtr GetDiscrepancyTestSettingName (DiscrepancyType dtype) 
{
  return discrepancy_info_list[dtype].setting_name;
}

extern DiscrepancyType GetDiscrepancyTypeFromSettingName (CharPtr setting_name)
{
  Int4 i;

  if (StringHasNoText (setting_name)) {
    return MAX_DISC_TYPE;
  }
  for (i = 0; i < MAX_DISC_TYPE; i++) {
    if (StringICmp (setting_name, discrepancy_info_list[i].setting_name) == 0) {
      return (DiscrepancyType) i;
    }
  }
  return MAX_DISC_TYPE;
}


/* Note that this function contains a hack - it assumes that all of the
 * test types that use the same collection function are listed together.
 */
extern ValNodePtr CollectDiscrepancies (DiscrepancyConfigPtr dcp, ValNodePtr sep_list, PerformDiscrepancyTest taxlookup)
{
  ValNodePtr             discrepancy_list = NULL;
  Int4                   i;
  PerformDiscrepancyTest last_test_func = NULL;

  discrepancy_info_list[DISC_NO_TAXLOOKUP].test_func = taxlookup;
  discrepancy_info_list[DISC_BAD_TAXLOOKUP].test_func = taxlookup;

  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    if ((dcp == NULL || dcp->conf_list[i])
        && discrepancy_info_list[i].test_func != NULL
        && discrepancy_info_list[i].test_func != last_test_func)
    {
      discrepancy_info_list[i].test_func (&discrepancy_list, sep_list);
      last_test_func = discrepancy_info_list[i].test_func;
    }
  }
  
  /* because some tests are run together, need to remove unwanted results */
  RemoveUnwantedDiscrepancyItems (&discrepancy_list, dcp);

  /* normalize the discrepancy levels so that they will be correctly displayed */
  SetDiscrepancyLevels (discrepancy_list, 0);
  return discrepancy_list;  
}


static CharPtr GetLocusTagForFeature (SeqFeatPtr sfp)
{
  GeneRefPtr grp = NULL;
  SeqFeatPtr gene;
  
  if (sfp == NULL) {
    return NULL;
  }
  if (sfp->idx.subtype == FEATDEF_GENE) {
    grp = sfp->data.value.ptrvalue;
  } else {
    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL) {
      gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
      if (gene != NULL) {
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
      }
    }
  }

  if (grp == NULL) {
    return NULL;
  } else {
    return grp->locus_tag;
  }
}


static CharPtr GetBioseqLabel (BioseqPtr bsp)
{
  Char        id_str[45];

  if (bsp == NULL) {
    return NULL;
  }
  
  SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str, PRINTID_REPORT, 39);
  return StringSave (id_str);
}


extern CharPtr GetBioseqSetLabel (BioseqSetPtr bssp)
{
  SeqEntryPtr sep;
  Char        id_str[45];
  BioseqPtr   bsp;

  if (bssp == NULL) {
    return NULL;
  }
  if (bssp->_class == BioseqseqSet_class_segset) {
    sprintf (id_str, "ss|");
  } else if (bssp->_class == BioseqseqSet_class_nuc_prot) {
    sprintf (id_str, "np|");
  } else {
    return NULL;
  }
  sep = bssp->seq_set;
  if (IS_Bioseq(sep) && sep->data.ptrvalue != NULL) {
    bsp = sep->data.ptrvalue;
    SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_str + 3, PRINTID_REPORT, 39);
    return StringSave (id_str);
  }
  return NULL;
}

extern CharPtr GetDiscrepancyItemText (ValNodePtr vnp)
{
  CharPtr           row_text = NULL, tmp;
  SeqFeatPtr        sfp, cds;
  BioseqPtr         bsp;
  SeqMgrFeatContext context;
  CharPtr           location;
  CharPtr           label;
  SeqDescrPtr       sdp;
  CharPtr           locus_tag = "";
  ObjValNodePtr     ovn;
  SeqEntryPtr       sep;
  
  if (vnp == NULL)
  {
    return NULL;
  }
  if (vnp->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
      if (SeqMgrFeaturesAreIndexed(sfp->idx.entityID) == 0) {
        SeqMgrIndexFeatures (sfp->idx.entityID, NULL);
      }

      sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &context);
      if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        if (bsp != NULL) {
          cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
          if (cds != NULL) {
            sfp = cds;
          }
        }
      }
      if (sfp != NULL)
      {
        location = SeqLocPrintUseBestID (sfp->location);
        label = (CharPtr) FeatDefTypeLabel(sfp);
        locus_tag = GetLocusTagForFeature (sfp);

        row_text = (CharPtr) MemNew (sizeof (Char) * 
                                     (StringLen (label) 
                                      + StringLen (context.label) 
                                      + StringLen (location) 
                                      + StringLen (locus_tag)
                                      + 6));
        sprintf (row_text, "%s\t%s\t%s\t%s\n", label, context.label, location, locus_tag == NULL ? "" : locus_tag);
        location = MemFree (location);
      }
    }
  }
  else if (vnp->choice == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL)
    {
      tmp = GetBioseqLabel (vnp->data.ptrvalue);
      row_text = (CharPtr) MemNew (sizeof(Char) * (StringLen (tmp) + 2));
      sprintf (row_text, "%s\n", tmp);
      tmp = MemFree (tmp);
    }
  }
  else if (vnp->choice == OBJ_BIOSEQSET) 
  {
    tmp = GetBioseqSetLabel (vnp->data.ptrvalue);
    row_text = (CharPtr) MemNew (sizeof(Char) * (StringLen (tmp) + 2));
    sprintf (row_text, "%s\n", tmp);
    tmp = MemFree (tmp);
  }
  else if (vnp->choice == OBJ_SEQENTRY)
  {
    sep = (SeqEntryPtr) vnp->data.ptrvalue;
    if (sep != NULL && sep->data.ptrvalue != NULL) {
      tmp = NULL;
      if (IS_Bioseq(sep)) {
        tmp = GetBioseqLabel (sep->data.ptrvalue);
      } else if (IS_Bioseq_set (sep)) {
        tmp = GetBioseqSetLabel (sep->data.ptrvalue);
      }
      if (tmp != NULL) {
        row_text = (CharPtr) MemNew (sizeof(Char) * (StringLen (tmp) + 2));
        sprintf (row_text, "%s\n", tmp);
        tmp = MemFree (tmp);
      }
    }
  }
  else if (vnp->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) vnp->data.ptrvalue;
    if (sdp != NULL)
    {
      bsp = NULL;
      if (sdp->extended != 0) {
        ovn = (ObjValNodePtr) sdp;
        if (ovn->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovn->idx.parentptr;
        }
      }
      if (bsp == NULL) {
        row_text = (CharPtr) MemNew (sizeof (Char) * 61);
        SeqDescLabel (sdp, row_text, 59, TRUE);
      } else {
        row_text = (CharPtr) MemNew (sizeof (Char) * (61 + 41));
        SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), row_text, PRINTID_REPORT, 39);
        row_text[39] = 0;
        StringCat (row_text, ":");
        SeqDescLabel (sdp, row_text + StringLen (row_text), 59, TRUE);
      }
      StringCat (row_text, "\n");
    }
  }
  
  return row_text;
}

extern CharPtr GetParentLabelForDiscrepancyItem (ValNodePtr vnp)
{
  CharPtr label = NULL;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovn;
  BioseqPtr  bsp;

  if (vnp == NULL || vnp->data.ptrvalue == NULL) return NULL;

  switch (vnp->choice)
  {
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) vnp->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      label = GetBioseqLabel (bsp);
      break;
    case OBJ_SEQDESC:
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;
      if (sdp != NULL)
      {
        if (sdp->extended != 0) {
          ovn = (ObjValNodePtr) sdp;
          if (ovn->idx.parenttype == OBJ_BIOSEQ) {
            label = GetBioseqLabel ((BioseqPtr) ovn->idx.parentptr);
          } else if (ovn->idx.parenttype == OBJ_BIOSEQSET) {
            label = GetBioseqSetLabel ((BioseqSetPtr) ovn->idx.parentptr);
          }
        }
      }
      break;
    case OBJ_BIOSEQ:
      label = GetBioseqLabel (vnp->data.ptrvalue);
      break;
    case OBJ_BIOSEQSET:
      label = GetBioseqSetLabel (vnp->data.ptrvalue);
      break;
  }
  return label;
}


static ValNodePtr ValNodePointerDup (ValNodePtr vnp)
{
  ValNodePtr vnp_new = NULL;
  
  if (vnp != NULL)
  {
    vnp_new = ValNodeNew (NULL);
    vnp_new->choice = vnp->choice;
    vnp_new->data.ptrvalue = vnp->data.ptrvalue;
    vnp_new->next = ValNodePointerDup (vnp->next);
  }
  return vnp_new;
}


typedef struct ftstrings {
  CharPtr header;
  CharPtr desc;
} FTStringsData, PNTR FTStringsPtr;


static FTStringsPtr FTStringsNew (CharPtr header, CharPtr desc)
{
  FTStringsPtr f;

  f = (FTStringsPtr) MemNew (sizeof (FTStringsData));
  f->header = header;
  f->desc = desc;
  return f;
}


static FTStringsPtr FTStringsFree (FTStringsPtr f)
{
  if (f != NULL) {
    f->header = MemFree (f->header);
    f->desc = MemFree (f->desc);
    f = MemFree (f);
  }
  return f;
}


extern ValNodePtr ReplaceDiscrepancyItemWithFeatureTableStrings (ValNodePtr feat_list)
{
  BioseqPtr       bsp, prot_bsp;
  CstType         custom_flags = 0;
  Asn2gbJobPtr    ajp;
  BaseBlockPtr    bbp;
  Int4            index;
  SeqFeatPtr      sfp, cds;
  ValNodePtr      vnp, list_copy = NULL, list_vnp;
  CharPtr         feature_table_header = NULL, feat_desc;
  FTStringsPtr    fts;
  
  if (feat_list == NULL) return NULL;
  
  list_copy = ValNodePointerDup (feat_list);
  for (vnp = list_copy; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->choice != OBJ_SEQFEAT || vnp->data.ptrvalue == NULL) continue;

    sfp = (SeqFeatPtr) vnp->data.ptrvalue;

    if (sfp->idx.subtype == FEATDEF_PROT) {
      prot_bsp = BioseqFindFromSeqLoc (sfp->location);
      if (prot_bsp != NULL) {
        cds = SeqMgrGetCDSgivenProduct (prot_bsp, NULL);
        if (cds != NULL) {
          sfp = cds;
        }
      }      
    }      
    bsp = BioseqFindFromSeqLoc (sfp->location);
    feature_table_header = NULL;
    ajp = asn2gnbk_setup (bsp, NULL, NULL, FTABLE_FMT, DUMP_MODE, NORMAL_STYLE,
                          0, 0, custom_flags, NULL);
    if (ajp == NULL) {
      continue;
    }

    for (index = 0; index < ajp->numParagraphs; index++) 
    {
      bbp = ajp->paragraphArray [index];
      if (bbp->blocktype == FEATHEADER_BLOCK) {
        feature_table_header = asn2gnbk_format (ajp, (Int4) index);
      } else if (bbp->blocktype == FEATURE_BLOCK) {
        for (list_vnp = vnp; list_vnp != NULL; list_vnp = list_vnp->next)
        {
          if (list_vnp->choice != OBJ_SEQFEAT || list_vnp->data.ptrvalue == NULL) continue;
          sfp = (SeqFeatPtr) list_vnp->data.ptrvalue;
          if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT) {
            prot_bsp = BioseqFindFromSeqLoc (sfp->location);
            if (prot_bsp != NULL) {
              cds = SeqMgrGetCDSgivenProduct (prot_bsp, NULL);
              if (cds != NULL) {
                sfp = cds;
              }
            }
          }      
           
          if (sfp != NULL 
              && bbp->entityID == sfp->idx.entityID
              && bbp->itemtype == sfp->idx.itemtype
              && bbp->itemID == sfp->idx.itemID)
          {
            /* replace list feature with description, change choice */
            list_vnp->choice = 0;
            feat_desc = asn2gnbk_format (ajp, (Int4) index);
            list_vnp->data.ptrvalue = FTStringsNew (StringSave (feature_table_header), feat_desc);
          }
        }
      }
    }
    asn2gnbk_cleanup (ajp);
    feature_table_header = MemFree (feature_table_header);
  }

  /* now remove redundant headers */
  for (list_vnp = list_copy; list_vnp != NULL; list_vnp = list_vnp->next) {
    if (list_vnp->choice != 0) continue;
    fts = (FTStringsPtr) list_vnp->data.ptrvalue;
    if (feature_table_header == NULL
        || StringCmp (feature_table_header, fts->header) != 0) {
      feature_table_header = MemFree (feature_table_header);
      feature_table_header = fts->header;
      fts->header = NULL;
      list_vnp->data.ptrvalue = (CharPtr) MemNew (sizeof (Char) * (StringLen (feature_table_header) + StringLen (fts->desc) + 2));
      StringCpy (list_vnp->data.ptrvalue, feature_table_header);
      StringCat (list_vnp->data.ptrvalue, fts->desc);
    } else {
      list_vnp->data.ptrvalue = fts->desc;
      fts->desc = NULL;
    }
    fts = FTStringsFree (fts);
  }
  feature_table_header = MemFree (feature_table_header);
  return list_copy;
}

static CharPtr GetItemFilename (ValNodePtr item, ValNodePtr filename_list)
{
  ValNodePtr       f_vnp;
  SeqFeatPtr       sfp;
  BioseqPtr        bsp;
  Uint2            entityID = 0;
  SeqDescrPtr      sdp;
  ObjValNodePtr    ovn;

  if (item == NULL || filename_list == NULL) return NULL;

  if (item->choice == OBJ_SEQFEAT)
  {
    sfp = (SeqFeatPtr) item->data.ptrvalue;
    if (sfp != NULL)
    {
      entityID = sfp->idx.entityID;
    }
  }
  else if (item->choice == OBJ_BIOSEQ)
  {
    bsp = (BioseqPtr) item->data.ptrvalue;
    if (bsp != NULL)
    {
      entityID = bsp->idx.entityID;
    }
  }
  else if (item->choice == OBJ_SEQDESC)
  {
    sdp = (SeqDescrPtr) item->data.ptrvalue;
    if (sdp != NULL)
    {
      bsp = NULL;
      if (sdp->extended != 0) {
        ovn = (ObjValNodePtr) sdp;
        entityID = ovn->idx.entityID;
      }
    }
  }
  
  for (f_vnp = filename_list; f_vnp != NULL; f_vnp = f_vnp->next) {
    if (f_vnp->choice == FILENAME_LIST_ENTITY_ID_ITEM && f_vnp->data.intvalue == entityID) {
      f_vnp = f_vnp->next;
      if (f_vnp != NULL && f_vnp->choice == FILENAME_LIST_FILENAME_ITEM && f_vnp->data.ptrvalue != NULL) {
        return f_vnp->data.ptrvalue;
      }
      break;
    }
  }
  return NULL;
}


extern void WriteDiscrepancyEx (FILE *fp, ClickableItemPtr dip, Boolean use_feature_table_fmt, ValNodePtr filename_list, CharPtr descr_prefix)
{
  ValNodePtr vnp, list_copy = NULL;
  CharPtr    row_text, filename;
  
  if (fp == NULL || dip == NULL)
  {
    return;
  }
  
  if (!StringHasNoText (descr_prefix)) {
    fprintf (fp, "%s:", descr_prefix);
  }
  fprintf (fp, "%s\n", dip->description);
  vnp = dip->item_list;
  
  if (use_feature_table_fmt)
  {
    list_copy = ReplaceDiscrepancyItemWithFeatureTableStrings (vnp);
    vnp = list_copy;
  }

  while (vnp != NULL)
  {
    filename = NULL;
    if (vnp->choice == 0)
    {
      row_text = vnp->data.ptrvalue;
    }
    else
    {
      row_text = GetDiscrepancyItemText (vnp);
      filename = GetItemFilename (vnp, filename_list);
    }
    if (row_text != NULL)
    {
      if (filename != NULL) {
        fprintf (fp, "%s:", filename);
      }
      fprintf (fp, row_text);
      row_text = MemFree (row_text);
    }
    vnp = vnp->next;
  }
  
  fprintf (fp, "\n");
}

extern void WriteDiscrepancy (FILE *fp, ClickableItemPtr dip, Boolean use_feature_table_fmt)
{
  WriteDiscrepancyEx (fp, dip, use_feature_table_fmt, NULL, NULL);
}









/* DiscrepancyConfig functions */
extern DiscrepancyConfigPtr DiscrepancyConfigFree (DiscrepancyConfigPtr dcp)
{
  return MemFree (dcp);  
}

extern void DisableTRNATests (DiscrepancyConfigPtr dcp)
{
  if (dcp != NULL) {
    dcp->conf_list[DISC_COUNT_TRNA] = FALSE;
    dcp->conf_list[DISC_DUP_TRNA] = FALSE;
    dcp->conf_list[DISC_BADLEN_TRNA] = FALSE;
    dcp->conf_list[DISC_COUNT_RRNA] = FALSE;
    dcp->conf_list[DISC_DUP_RRNA] = FALSE;
    dcp->conf_list[DISC_TRANSL_NO_NOTE] = FALSE;
    dcp->conf_list[DISC_NOTE_NO_TRANSL] = FALSE;
    dcp->conf_list[DISC_TRANSL_TOO_LONG] = FALSE;
    dcp->conf_list[DISC_CDS_OVERLAP_TRNA] = FALSE;
    dcp->conf_list[DISC_COUNT_PROTEINS] = FALSE;
  }
}

extern DiscrepancyConfigPtr DiscrepancyConfigNew (void)
{
  DiscrepancyConfigPtr dcp;
  Int4                 i;
  
  dcp = (DiscrepancyConfigPtr) MemNew (sizeof (DiscrepancyConfigData));
  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    dcp->conf_list[i] = TRUE;
  }

  /* by default, tRNA and rRNA tests are off */
  DisableTRNATests (dcp);

  dcp->use_feature_table_format = FALSE;
  return dcp;
}

extern DiscrepancyConfigPtr ReadDiscrepancyConfig (void)
{
  DiscrepancyConfigPtr dcp;
  Int4                 i;
  Char                 str[20];
  
  dcp = DiscrepancyConfigNew();
  if (dcp != NULL)
  {
    for (i = 0; i < MAX_DISC_TYPE; i++)
    {
      if (GetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", discrepancy_info_list[i].setting_name, NULL, str, sizeof (str))) {
        if (StringICmp (str, "FALSE") == 0) {
          dcp->conf_list[i] = FALSE;
        } else if (StringICmp (str, "TRUE") == 0) {
          dcp->conf_list[i] = TRUE;
        }
      }
    }
    if (GetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", "USE_FEATURE_TABLE_FORMAT", NULL, str, sizeof (str))) {
      if (StringICmp (str, "TRUE") == 0) {
        dcp->use_feature_table_format = TRUE;
      }
    }
  }
  return dcp;
}


extern void SaveDiscrepancyConfig (DiscrepancyConfigPtr dcp)
{
  Int4 i;
  
  if (dcp == NULL)
  {
    return;
  }
  
  for (i = 0; i < MAX_DISC_TYPE; i++)
  {
    if (dcp->conf_list[i])
    {
      SetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", discrepancy_info_list[i].setting_name, "TRUE");
    }
    else
    {
      SetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", discrepancy_info_list[i].setting_name, "FALSE");
    }
  }
  if (dcp->use_feature_table_format)
  {
    SetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", "USE_FEATURE_TABLE_FORMAT", "TRUE");
  }
  else
  {
    SetAppParam ("SEQUINCUSTOM", "DISCREPANCY_REPORT", "USE_FEATURE_TABLE_FORMAT", "FALSE");
  }
}


extern CharPtr SetDiscrepancyReportTestsFromString (CharPtr list, Boolean enable, DiscrepancyConfigPtr dcp)
{
  CharPtr         ptr, tmp, name_start, err_msg;
  DiscrepancyType test_type;
  CharPtr         err_fmt = "%s is an unrecognized test name";
  
  if (dcp == NULL) return StringSave ("Unable to configure");

  if (!StringDoesHaveText (list)) {
      return StringSave ("No tests specified!");
  }

  tmp = StringSave (list);
  name_start = tmp;
  while (name_start != NULL && StringDoesHaveText (name_start)) {
    ptr = StringChr (name_start, ',');
    if (ptr != NULL) {
      *ptr = 0;
    }
    TrimSpacesAroundString (name_start);
    test_type = GetDiscrepancyTypeFromSettingName (name_start);
    if (test_type == MAX_DISC_TYPE) {
      err_msg = (CharPtr) MemNew (StringLen (err_fmt) + StringLen (name_start));
      sprintf (err_msg, err_fmt, name_start);
      tmp = MemFree (tmp);
      return err_msg;
    }
    dcp->conf_list[test_type] = enable;
    if (ptr == NULL) {
      name_start = NULL;
    } else {
      name_start = ptr + 1;
    }
  }
  tmp = MemFree (tmp);
  return NULL;  
}


/* functions for writing discrepancy report to file */
extern void WriteAsnDiscReport (ValNodePtr discrepancy_list, FILE *ofp, DiscReportOutputConfigPtr oc, Boolean use_flag)
{
  ValNodePtr       vnp;
  ClickableItemPtr cip;
  CharPtr          setting_name, prefix;
  CharPtr          prefix_fmt = "DiscRep:%s:";

  if (ofp == NULL || oc == NULL) return;

  for (vnp = discrepancy_list; vnp != NULL; vnp = vnp->next) {
    cip = (ClickableItemPtr) vnp->data.ptrvalue;
    if (cip != NULL) {
      prefix = NULL;
      if (use_flag) {
        setting_name = GetDiscrepancyTestSettingName (cip->clickable_item_type);
        if (StringHasNoText (setting_name)) {
          prefix = StringSave ("DiscRep:");
        } else {
          prefix = (CharPtr) MemNew (sizeof (Char) * (StringLen (prefix_fmt) + StringLen (setting_name)));
          sprintf (prefix, prefix_fmt, setting_name);
        }
      }
      if (oc->summary_report) {
        fprintf (ofp, "%s%s\n", prefix == NULL ? "" : prefix, cip->description);           
      } else {
        WriteDiscrepancyEx (ofp, cip, oc->use_feature_table_format, oc->filename_list, prefix);
      }
      prefix = MemFree (prefix);
      if ((cip->item_list == NULL || oc->expand_report_categories[cip->clickable_item_type]) && cip->subcategories != NULL) {
        if (use_flag && cip->clickable_item_type == DISC_INCONSISTENT_BIOSRC_DEFLINE) {
          WriteAsnDiscReport (cip->subcategories, ofp, oc, FALSE);
        } else {
          WriteAsnDiscReport (cip->subcategories, ofp, oc, use_flag);
        }
      }
    }
  }

}


static int LIBCALLBACK SortVnpByDiscrepancyType (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  ClickableItemPtr c1, c2;
  CharPtr          cp1, cp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      c1 = (ClickableItemPtr) vnp1->data.ptrvalue;
      c2 = (ClickableItemPtr) vnp2->data.ptrvalue;
      if (c1 != NULL && c2 != NULL) {
        if (c1->clickable_item_type < c2->clickable_item_type) {
          return -1;
        } else if (c1->clickable_item_type > c2->clickable_item_type) {
          return 1;
        } else {
          if (c1->description == NULL && c2->description == NULL) {
            return 0;
          } else if (c1->description == NULL) {
            return -1;
          } else if (c2->description == NULL) {
            return 1;
          } else {
            cp1 = c1->description;
            while (isdigit (*cp1)) {
              cp1++;
            }
            cp2 = c2->description;
            while (isdigit (*cp2)) {
              cp2++;
            }
            return StringCmp (cp1, cp2);
          }
        }
      }
    }
  }
  return 0;
}


static ClickableItemPtr CombineDiscrepancyReports (ClickableItemPtr cip1, ClickableItemPtr cip2)
{
  CharPtr cp1, cp2, num_start1, num_start2, num_buf;
  Char    fixed_buf[15];
  Int4    common_start_len = 0, common_len_end = 0;
  Int4    num_len1, num_len2, num_items1, num_items2;
  ClickableItemPtr combined = NULL;
  

  if (cip1 == NULL || cip2 == NULL || cip1->clickable_item_type != cip2->clickable_item_type
      || StringHasNoText (cip1->description) || StringHasNoText (cip2->description)) {
    return NULL;
  }

  cp1 = cip1->description;
  cp2 = cip2->description;
    
  while (*cp1 == *cp2 && *cp1 != 0 && *cp2 != 0 && !isdigit (*cp1)) {
    cp1++;
    cp2++;
    common_start_len++;
  }
  if (*cp1 == 0 && *cp2 == 0) {
    /* entire description matches */
    combined = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    combined->clickable_item_type = cip1->clickable_item_type;
    combined->description = StringSave (cip1->description);
    combined->item_list = cip1->item_list;
    cip1->item_list = NULL;
    combined->subcategories = cip1->subcategories;
    cip1->subcategories = NULL;
    ValNodeLink (&(combined->item_list), cip2->item_list);
    cip2->item_list = NULL;
    ValNodeLink (&(combined->subcategories), cip2->subcategories);
    cip2->subcategories = NULL;
  } else if (isdigit (*cp1) && isdigit (*cp2) && (cp1 == cip1->description || isspace (*(cp1 - 1)))) {
    num_start1 = cp1;
    num_len1 = 0;
    while (isdigit (*cp1)) {
      cp1++;
      num_len1++;
    }
    num_start2 = cp2;
    num_len2 = 0;
    while (isdigit (*cp2)) {
      cp2++;
      num_len2++;
    }
    if ((*cp1 == 0 || isspace (*cp1)) && StringCmp (cp1, cp2) == 0) {
      /* matches on the other side of the number */
      /* build combined description */
      if (num_len1 < sizeof (fixed_buf)) {
        StringNCpy (fixed_buf, num_start1, num_len1);
        fixed_buf[num_len1] = 0;
        num_items1 = atoi(fixed_buf);
      } else {
        num_buf = (CharPtr) MemNew (sizeof (Char) * (num_len1 + 1));
        StringNCpy (num_buf, num_start1, num_len1);
        num_buf[num_len1] = 0;
        num_items1 = atoi (num_buf);
        num_buf = MemFree (num_buf);
      }
      if (num_len2 < sizeof (fixed_buf) - 1) {
        StringNCpy (fixed_buf, num_start2, num_len2);
        fixed_buf[num_len2] = 0;
        num_items2 = atoi(fixed_buf);
      } else {
        num_buf = (CharPtr) MemNew (sizeof (Char) * (num_len2 + 1));
        StringNCpy (num_buf, num_start2, num_len2);
        num_buf[num_len2] = 0;
        num_items2 = atoi (num_buf);
        num_buf = MemFree (num_buf);
      }
      
      combined = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));

      combined->description = (CharPtr) MemNew (sizeof (Char) * (common_start_len + sizeof (fixed_buf) + StringLen (cp1) + 1));
      StringNCpy (combined->description, cip1->description, common_start_len);
      sprintf (fixed_buf, "%d", num_items1 + num_items2);
      StringCat (combined->description, fixed_buf);
      StringCat (combined->description, cp1);

      combined->clickable_item_type = cip1->clickable_item_type;  
      combined->item_list = cip1->item_list;
      cip1->item_list = NULL;
      combined->subcategories = cip1->subcategories;
      cip1->subcategories = NULL;
      ValNodeLink (&(combined->item_list), cip2->item_list);
      cip2->item_list = NULL;
      ValNodeLink (&(combined->subcategories), cip2->subcategories);
      cip2->subcategories = NULL;
    } else {
      combined = NULL;
    }
  }
  if (combined != NULL && combined->subcategories != NULL) {
    CollateDiscrepancyReports (&(combined->subcategories));
  }
  return combined;
}


extern void CollateDiscrepancyReports (ValNodePtr PNTR discrepancy_reports)
{
  ValNodePtr vnp, tmp;
  ClickableItemPtr combined;

  *discrepancy_reports = ValNodeSort (*discrepancy_reports, SortVnpByDiscrepancyType);

  vnp = *discrepancy_reports;
  while (vnp != NULL && vnp->next != NULL) {
    combined = CombineDiscrepancyReports (vnp->data.ptrvalue, vnp->next->data.ptrvalue);
    if (combined != NULL) {
      vnp->data.ptrvalue = ClickableItemFree (vnp->data.ptrvalue);
      vnp->next->data.ptrvalue = ClickableItemFree (vnp->next->data.ptrvalue);
      tmp = vnp->next;
      vnp->next = vnp->next->next;
      tmp->next = NULL;
      tmp = ValNodeFree (tmp);
      vnp->data.ptrvalue = combined;
    } else {
      vnp = vnp->next;
    }
  }
}


extern CharPtr ExpandDiscrepancyReportTestsFromString (CharPtr list, Boolean expand, DiscReportOutputConfigPtr dcp)
{
  CharPtr         ptr, tmp, name_start, err_msg;
  Int4            i;
  DiscrepancyType test_type;
  CharPtr         err_fmt = "%s is an unrecognized test name";
  
  if (dcp == NULL) return StringSave ("Unable to configure");

  if (!StringDoesHaveText (list)) {
    return NULL;
  } else if (StringICmp (list, "all") == 0) {
    for (i = 0; i < MAX_DISC_TYPE; i++) {
      dcp->expand_report_categories[i] = expand;
    }
  } else {
    tmp = StringSave (list);
    name_start = tmp;
    while (name_start != NULL && StringDoesHaveText (name_start)) {
      ptr = StringChr (name_start, ',');
      if (ptr != NULL) {
        *ptr = 0;
      }
      TrimSpacesAroundString (name_start);
      test_type = GetDiscrepancyTypeFromSettingName (name_start);
      if (test_type == MAX_DISC_TYPE) {
        err_msg = (CharPtr) MemNew (StringLen (err_fmt) + StringLen (name_start));
        sprintf (err_msg, err_fmt, name_start);
        tmp = MemFree (tmp);
        return err_msg;
      }
      dcp->expand_report_categories[test_type] = expand;
      if (ptr == NULL) {
        name_start = NULL;
      } else {
        name_start = ptr + 1;
      }
    }
    tmp = MemFree (tmp);
  }
  return NULL;  
}


extern ValNodePtr FreeFilenameList (ValNodePtr filename_list)
{
  ValNodePtr vnp_next;
  
  while (filename_list != NULL) {
    vnp_next = filename_list->next;
    filename_list->next = NULL;
    if (filename_list->choice == FILENAME_LIST_FILENAME_ITEM) {
      filename_list = ValNodeFreeData (filename_list);
    } else {
      filename_list = ValNodeFree (filename_list);
    }
    filename_list = vnp_next;
  }
  return filename_list;
}


/* Barcode Discrepancy Function */

/* 
 * list of names for the individual tests.
 * Note - this array should have eBarcodeTest_LAST elements (see sqnutils.h for value of eBarcodeTest_LAST).
 */
static CharPtr BarcodeTestNames[] = 
{ "Too Short",
  "Missing Primers",
  "Missing Country",
  "Missing Specimen Voucher",
  "Too Many Ns"
};


extern CharPtr GetBarcodeTestName (Int4 i)
{
  if (i < 0 || i >= sizeof (BarcodeTestNames) / sizeof (CharPtr)) 
  {
    return NULL;
  } 
  else 
  {
    return BarcodeTestNames[i];
  }
}


extern Int4 GetBarcodeTestNumFromBarcodeTestName (CharPtr test_name)
{
  Int4 i;

  if (StringHasNoText (test_name)) {
    return eBarcodeTest_LAST;
  }
  for (i = 0; i < eBarcodeTest_LAST; i++) {
    if (StringICmp (test_name, BarcodeTestNames[i]) == 0) {
      return i;
    }
  }
  return eBarcodeTest_LAST;
}


/* Functions for creating and freeing configurations for the Barcode Tests. */

extern BarcodeTestConfigPtr BarcodeTestConfigNew()
{
  BarcodeTestConfigPtr cfg;
  Int4                 i;

  cfg = (BarcodeTestConfigPtr) MemNew (sizeof (BarcodeTestConfigData));
  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    cfg->conf_list[i] = TRUE;
  }
  cfg->min_length = 500;
  cfg->min_n_percent = 1.0;
  return cfg;
}


extern BarcodeTestConfigPtr BarcodeTestConfigFree (BarcodeTestConfigPtr cfg)
{
  if (cfg != NULL)
  {
    cfg = MemFree (cfg);
  }
  return cfg;
}


/* A BarcodeTestResults lists the Bioseq that the test was performed on,
 * indicates whether each test passed, and give the percentage of Ns
 * (if the value is above the minimum in the configuration).
 */
extern BarcodeTestResultsPtr BarcodeTestResultsNew ()
{
  BarcodeTestResultsPtr res;

  res = (BarcodeTestResultsPtr) MemNew (sizeof (BarcodeTestResultsData));
  MemSet (res, 0, sizeof (BarcodeTestResultsData));
  return res;
}


extern BarcodeTestResultsPtr BarcodeTestResultsFree (BarcodeTestResultsPtr res)
{
  if (res != NULL)
  {
    res = MemFree (res);
  }
  return res;
}


extern BarcodeTestResultsPtr BarcodeTestResultsCopy (BarcodeTestResultsPtr res)
{
  BarcodeTestResultsPtr res_new = NULL;

  if (res != NULL)
  {
    res_new = BarcodeTestResultsNew();
    MemCopy (res_new, res, sizeof (BarcodeTestResultsData));
  }
  return res_new;
}


extern ValNodePtr BarcodeTestResultsListFree (ValNodePtr res_list)
{
  ValNodePtr vnp;

  if (res_list != NULL) 
  {
    vnp = res_list->next;
    res_list->next = NULL;
    res_list->data.ptrvalue = BarcodeTestResultsFree (res_list->data.ptrvalue);
    ValNodeFree (res_list);
    BarcodeTestResultsListFree (vnp);
  }
  return res_list;
}


/* determines whether barcode tests should be performed on a sequence -
 * no barcode keyword, no barcode tests needed.
 */
static Boolean HasBARCODEKeyword (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  MolInfoPtr        mip;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_molinfo, &dcontext))
  {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->tech == MI_TECH_barcode)
    {
      found = TRUE;
    }
  }
  return found;
}

/*
 * Finds the MolInfo descriptor for the Bioseq and removes the BARCODE technique.
 * Returns true if the BARCODE technique was present before it was removed.
 */
static Boolean RemoveBarcodeTechFromBioseq (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  MolInfoPtr        mip;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext))
  {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->tech == MI_TECH_barcode)
    {
      mip->tech = MI_TECH_unknown;
      found = TRUE;
    }
  }
  return found;
}


/* 
 * Adds the BARCODE technique to the MolInfo descriptor for the Bioseq.
 * Will create a new MolInfo descriptor for the Bioseq if it doesn't 
 * find one already there.
 */
static void ApplyBarcodeTechToBioseq (BioseqPtr bsp)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext dcontext;
  Boolean           found = FALSE;
  MolInfoPtr        mip;
  SeqEntryPtr       sep;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext))
  {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip == NULL)
    {
      mip = MolInfoNew ();
      sdp->data.ptrvalue = mip;
    }
    mip->tech = MI_TECH_barcode;
    found = TRUE;
  }

  if (!found) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
    mip = MolInfoNew();
    mip->tech = MI_TECH_barcode;
    sdp->data.ptrvalue = mip;
  }
}


/* Used for generating Discrepancy Report style data,
 * where Bioseqs are listed separately for each test they fail.
 */
typedef struct barcodesearch {
  ValNodePtr           bioseq_list;
  BarcodeTestConfigPtr cfg;
} BarcodeSearchData, PNTR BarcodeSearchPtr;

static void FindShortBarcodeSequencesCallback (BioseqPtr bsp, Pointer userdata)
{
  BarcodeSearchPtr bsd;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL || !HasBARCODEKeyword(bsp)) return;

  bsd = (BarcodeSearchPtr) userdata;
  if (bsp->length < bsd->cfg->min_length) 
  { 
    ValNodeAddPointer (&(bsd->bioseq_list), OBJ_BIOSEQ, bsp);
  }
}  

typedef Boolean (*BarcodeBioSourceTestFunc) PROTO ((BioSourcePtr));

static Boolean HasForwardAndReversePrimers (BioSourcePtr biop)
{
  Boolean           found_fwd_seq = FALSE;
  Boolean           found_rev_seq = FALSE;
  SubSourcePtr      ssp;

  if (biop == NULL || biop->subtype == NULL) return FALSE;

  for (ssp = biop->subtype; ssp != NULL && (!found_fwd_seq || !found_rev_seq); ssp = ssp->next)
  {
    if (ssp->subtype == SUBSRC_fwd_primer_seq)
    {
      found_fwd_seq = TRUE;
    }
    else if (ssp->subtype == SUBSRC_rev_primer_seq)
    {
      found_rev_seq = TRUE;
    }
  }

  return found_fwd_seq && found_rev_seq;
}


static Boolean HasCountry (BioSourcePtr biop)
{
  SubSourcePtr      ssp;
  Boolean           found = FALSE;

  if (biop == NULL || biop->subtype == NULL) return FALSE;

  for (ssp = biop->subtype; ssp != NULL && !found; ssp = ssp->next)
  {
    if (ssp->subtype == SUBSRC_country)
    {
      found = TRUE;
    }
  }

  return found;
}

static Boolean HasSpecimenVoucher (BioSourcePtr biop)
{
  OrgModPtr mod;

  if (biop == NULL || biop->org == NULL || biop->org->orgname == NULL) return FALSE;

  for (mod = biop->org->orgname->mod; 
       mod != NULL && mod->subtype != ORGMOD_specimen_voucher; 
       mod = mod->next)
  {}

  if (mod == NULL) 
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static Boolean BarcodeBioSourceTest (BioseqPtr bsp, BarcodeBioSourceTestFunc test_func)
{
  SeqDescrPtr       sdp;
  BioSourcePtr      biop;
  SeqMgrDescContext context;
  Boolean           found = FALSE;
  
  if (bsp == NULL || ISA_aa (bsp->mol) || !HasBARCODEKeyword (bsp) || test_func == NULL)
  {
    return FALSE;
  }
    
  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
       sdp != NULL && !found;
       sdp = SeqMgrGetNextDescriptor (bsp,sdp, Seq_descr_source, &context))
  {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    found = test_func(biop);
  }

  return !found;
}


static void BarcodeBioSourceTestCallback (BioseqPtr bsp, Pointer userdata, BarcodeBioSourceTestFunc test_func)
{  
  BarcodeSearchPtr bsd;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL || test_func == NULL)
  {
    return;
  }

  bsd = (BarcodeSearchPtr) userdata;
  
  if (BarcodeBioSourceTest (bsp, test_func))
  {
    ValNodeAddPointer (&(bsd->bioseq_list), OBJ_BIOSEQ, bsp);
  }
}


static void FindMissingForwardAndReversePrimers (BioseqPtr bsp, Pointer userdata)
{  
  BarcodeBioSourceTestCallback (bsp, userdata, HasForwardAndReversePrimers);
}


static void FindMissingCountryAndLatLon (BioseqPtr bsp, Pointer userdata)
{  
  BarcodeBioSourceTestCallback (bsp, userdata, HasCountry);
}


static void FindMissingSpecimenVoucher (BioseqPtr bsp, Pointer userdata)
{  
  BarcodeBioSourceTestCallback (bsp, userdata, HasSpecimenVoucher);
}


static void 
BarcodeTestForSeqEntry 
(SeqEntryPtr          sep,
 ValNodePtr PNTR      discrepancy_list, 
 VisitBioseqsFunc     callback, 
 CharPtr              fmt,
 BarcodeTestConfigPtr cfg)
{
  BarcodeSearchData bsd;

  bsd.bioseq_list = NULL;
  bsd.cfg = cfg;
  if (bsd.cfg == NULL)
  {
    bsd.cfg = BarcodeTestConfigNew ();
  }
  VisitBioseqsInSep (sep, &bsd, callback);

  if (bsd.bioseq_list != NULL)
  {
    ValNodeAddPointer (discrepancy_list, 0, NewClickableItem (0, fmt, bsd.bioseq_list));
  }
  if (cfg != bsd.cfg)
  {
    bsd.cfg = BarcodeTestConfigFree (bsd.cfg);
  }
}

 
static void LIBCALLBACK CountNProc (CharPtr sequence, Pointer userdata)
{
  Int4Ptr p_i;
  CharPtr cp;

  if (sequence == NULL || userdata == NULL) return;
  p_i = (Int4Ptr) userdata;

  for (cp = sequence; *cp != 0; cp++)
  {
    if (*cp == 'N')
    {
      (*p_i) ++;
    }
  }
}


static FloatLo PercentNInBioseq (BioseqPtr bsp)
{
  Int4 num_n = 0;
  
  if (bsp->length == 0) return 0;

  SeqPortStream (bsp, STREAM_EXPAND_GAPS, (Pointer) &num_n, CountNProc);

  return ((FloatLo)num_n * 100) / (FloatLo) bsp->length;
}


static void PercentNDiscrepancy (BioseqPtr bsp, Pointer userdata)
{
  FloatLo pct;
  BarcodeSearchPtr bs;

  if (bsp == NULL || ISA_aa (bsp->mol) || HasBARCODEKeyword (bsp) || userdata == NULL)
  {
    return;
  }

  bs = (BarcodeSearchPtr) userdata;

  pct = PercentNInBioseq (bsp);
  if (pct > bs->cfg->min_n_percent) 
  {
    ValNodeAddPointer (&(bs->bioseq_list), OBJ_BIOSEQ, bsp);
  }
}

static void PercentNDiscrepanciesForSeqEntry (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list, BarcodeTestConfigPtr cfg)
{
  BarcodeSearchData bsd;
  ValNodePtr subcategories = NULL, vnp;
  ClickableItemPtr cip;
  CharPtr fmt = "Sequence has %.1f percent Ns";
  CharPtr top_fmt = "%d sequences have > %.1f%% Ns";
  FloatLo pct;

  bsd.bioseq_list = NULL;
  bsd.cfg = cfg;
  if (bsd.cfg == NULL)
  {
    bsd.cfg = BarcodeTestConfigNew ();
  }
  VisitBioseqsInSep (sep, &bsd, PercentNDiscrepancy);

  if (bsd.bioseq_list == NULL) return;

  for (vnp = bsd.bioseq_list; vnp != NULL; vnp = vnp->next)
  {
    cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
    MemSet (cip, 0, sizeof (ClickableItemData));
    cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + 5));
    pct = PercentNInBioseq (vnp->data.ptrvalue);
    sprintf (cip->description, fmt, pct);
    ValNodeAddPointer (&(cip->item_list), OBJ_BIOSEQ, vnp->data.ptrvalue);
    ValNodeAddPointer (&subcategories, 0, cip);
  }

  cip = (ClickableItemPtr) MemNew (sizeof (ClickableItemData));
  MemSet (cip, 0, sizeof (ClickableItemData));
  cip->description = (CharPtr) MemNew (sizeof (Char) * (StringLen (top_fmt) + 10));
  sprintf (cip->description, fmt, ValNodeLen (bsd.bioseq_list), bsd.cfg->min_n_percent);
  cip->item_list = bsd.bioseq_list;
  cip->subcategories = subcategories;
  ValNodeAddPointer (discrepancy_list, 0, cip);

  if (bsd.cfg != cfg)
  {
    bsd.cfg = BarcodeTestConfigFree (bsd.cfg);
  }
}


static void GetBarcodeDiscrepanciesForSeqEntry (SeqEntryPtr sep, ValNodePtr PNTR discrepancy_list, BarcodeTestConfigPtr cfg)
{
  if (cfg == NULL) return;

  if (cfg->conf_list[eBarcodeTest_Length]) 
  {
    BarcodeTestForSeqEntry (sep, discrepancy_list, FindShortBarcodeSequencesCallback, "%d sequences are shorter than 500 nucleotides", cfg);
  }
  if (cfg->conf_list[eBarcodeTest_Primers])
  {
    BarcodeTestForSeqEntry (sep, discrepancy_list, FindMissingForwardAndReversePrimers, "%d sequences are missing forward and/or reverse primers", cfg);
  }
  if (cfg->conf_list[eBarcodeTest_Country])
  {
    BarcodeTestForSeqEntry (sep, discrepancy_list, FindMissingCountryAndLatLon, "%d sequences are missing country", cfg);
  }
  if (cfg->conf_list[eBarcodeTest_SpecimenVoucher])
  {
    BarcodeTestForSeqEntry (sep, discrepancy_list, FindMissingSpecimenVoucher, "%d sequences are missing specimen voucher", cfg);
  }
  if (cfg->conf_list[eBarcodeTest_PercentN])
  {
    PercentNDiscrepanciesForSeqEntry (sep, discrepancy_list, cfg);
  }
}


extern ValNodePtr GetBarcodeDiscrepancies (ValNodePtr sep_list, BarcodeTestConfigPtr cfg)
{
  ValNodePtr    vnp, discrepancy_list = NULL;
  SeqEntryPtr   sep;
  BarcodeTestConfigPtr local_cfg;

  if (cfg == NULL) 
  {
    local_cfg = BarcodeTestConfigNew();
  }
  else
  {
    local_cfg = cfg;
  }

  for (vnp = sep_list; vnp != NULL; vnp = vnp->next)
  {
    sep = vnp->data.ptrvalue;
    GetBarcodeDiscrepanciesForSeqEntry (sep, &discrepancy_list, local_cfg);
  }

  if (local_cfg != cfg)
  {
    local_cfg = BarcodeTestConfigFree (local_cfg);
  }

  /* normalize the discrepancy levels so that they will be correctly displayed */
  SetDiscrepancyLevels (discrepancy_list, 0);
  return discrepancy_list;  
}


/* This section is used for generating the "Failure Report" and "Compliance Report". */

typedef struct barcodebioseqsearch {
  ValNodePtr           results_list;
  BarcodeTestConfigPtr cfg;
  Boolean              collect_positives;
} BarcodeBioseqSearchData, PNTR BarcodeBioseqSearchPtr;


extern Boolean IsBarcodeID (SeqIdPtr sip)
{
  DbtagPtr dbt;

  if (sip == NULL) return FALSE;
  if (sip->choice != SEQID_GENERAL) return FALSE;
  dbt = (DbtagPtr) sip->data.ptrvalue;
  if (dbt == NULL) return FALSE;
  if (StringICmp (dbt->db, "uoguelph") == 0) return TRUE;
  return FALSE;
}

#define cMaxBarcodeIDStringLen 200
#define cMaxGenbankIDStringLen 200

extern CharPtr BarcodeTestBarcodeIdString (BioseqPtr bsp)
{
  SeqIdPtr barcode_id;
  Char     barcode_id_str[cMaxBarcodeIDStringLen];

  if (bsp == NULL) return NULL;

  barcode_id = bsp->id;
  while (barcode_id != NULL && !IsBarcodeID (barcode_id))
  {
    barcode_id = barcode_id->next;
  }

  if (barcode_id == NULL) 
  {
    barcode_id = bsp->id;
    while (barcode_id != NULL && barcode_id->choice != SEQID_LOCAL)
    {
      barcode_id = barcode_id->next;
    }
  }

  if (barcode_id == NULL) 
  {
    sprintf (barcode_id_str, "NO");
  }
  else
  {
    SeqIdWrite (barcode_id, barcode_id_str, PRINTID_FASTA_SHORT, sizeof (barcode_id_str) - 1);
  }
  return StringSave (barcode_id_str);
}

extern CharPtr BarcodeTestGenbankIdString (BioseqPtr bsp)
{
  SeqIdPtr genbank_id;
  Char     genbank_id_str[cMaxGenbankIDStringLen];
  CharPtr  src, dst;

  genbank_id = bsp->id;
  while (genbank_id != NULL && genbank_id->choice != SEQID_GENBANK)
  {
    genbank_id = genbank_id->next;
  }
  if (genbank_id == NULL) 
  {
    sprintf (genbank_id_str, "NO");
  }
  else
  {
    SeqIdWrite (genbank_id, genbank_id_str, PRINTID_FASTA_SHORT, sizeof (genbank_id_str) - 1);
    if (StringNICmp (genbank_id_str, "gb|", 3) == 0) {
      src = genbank_id_str + 3;
      dst = genbank_id_str;
      while (*src != 0) {
        *dst = *src;
        dst++;
        src++;
      }
      dst[0] = 0;
    }
    if (genbank_id_str[StringLen (genbank_id_str) - 1] == '|') {
      genbank_id_str[StringLen (genbank_id_str) - 1] = 0;
    }
  }
  return StringSave (genbank_id_str);
}



static CharPtr GetBarcodeTestFailureReasons (BarcodeTestResultsPtr res)
{
  Int4             i, msg_len = 0;
  Boolean          any_failed = FALSE;
  CharPtr          msg;
  Char             pct[5];

  if (res == NULL || res->bsp == NULL) return NULL;

  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i]) 
    {
      msg_len += StringLen (GetBarcodeTestName (i)) + 2;
      if (i == eBarcodeTest_PercentN)
      {
        msg_len += 5;
      }
      any_failed = TRUE;
    }
  }
  if (!any_failed) return NULL;
   
  msg = (CharPtr) MemNew (sizeof (Char) * msg_len);
  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i])
    {
      StringCat (msg, GetBarcodeTestName(i));
      if (i == eBarcodeTest_PercentN) 
      {
        sprintf (pct, ":%.1f", res->n_percent);
        StringCat (msg, pct);
      }
      StringCat (msg, ",");
    }
  }
  /* remove trailing comma */
  msg[StringLen(msg) - 1] = 0;
          
  return msg;

}


static CharPtr SummaryTextFromBarcodeTestResults (BarcodeTestResultsPtr res)
{
  Int4             i, msg_len = 0;
  Boolean          any_failed = FALSE;
  CharPtr          msg, genbank_id, barcode_id;
  Char             pct[5];

  if (res == NULL || res->bsp == NULL) return NULL;

  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i]) 
    {
      msg_len += StringLen (GetBarcodeTestName (i)) + 2;
      if (i == eBarcodeTest_PercentN)
      {
        msg_len += 5;
      }
      any_failed = TRUE;
    }
  }
  if (!any_failed) return NULL;

  genbank_id = BarcodeTestGenbankIdString (res->bsp);
  barcode_id = BarcodeTestBarcodeIdString (res->bsp);
  
  msg_len += StringLen (genbank_id) + StringLen (barcode_id) + 2;
 
  msg = (CharPtr) MemNew (sizeof (Char) * msg_len);
  sprintf (msg, "%s\t%s\t", barcode_id, genbank_id);
  for (i = 0; i < eBarcodeTest_LAST; i++)
  {
    if (res->failed_tests[i])
    {
      StringCat (msg, GetBarcodeTestName(i));
      if (i == eBarcodeTest_PercentN) 
      {
        sprintf (pct, ":%.1f", res->n_percent);
        StringCat (msg, pct);
      }
      StringCat (msg, ",");
    }
  }
  /* remove trailing comma */
  msg[StringLen(msg) - 1] = 0;
          
  return msg;
}


extern Boolean PassBarcodeTests (BarcodeTestResultsPtr res)
{
  if (res == NULL
      || res->failed_tests[eBarcodeTest_Length]
      || res->failed_tests[eBarcodeTest_Primers]
      || res->failed_tests[eBarcodeTest_Country]
      || res->failed_tests[eBarcodeTest_SpecimenVoucher]
      || res->failed_tests[eBarcodeTest_PercentN])
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static BarcodeTestResultsPtr BarcodeTestResultsForBioseq (BioseqPtr bsp, BarcodeTestConfigPtr cfg)
{
  BarcodeTestResultsPtr res = NULL;

  if (bsp == NULL || ISA_aa (bsp->mol) || !HasBARCODEKeyword (bsp) || cfg == NULL)
  {
    return NULL;
  }

  res = BarcodeTestResultsNew ();

  res->bsp = bsp;

  if (bsp->length < cfg->min_length && cfg->conf_list[eBarcodeTest_Length])
  {
    res->failed_tests[eBarcodeTest_Length] = TRUE;
  } 

  if (cfg->conf_list[eBarcodeTest_Primers])
  {
    res->failed_tests[eBarcodeTest_Primers] = BarcodeBioSourceTest(bsp, HasForwardAndReversePrimers);
  }
  if (cfg->conf_list[eBarcodeTest_Country])
  {
    res->failed_tests[eBarcodeTest_Country] = BarcodeBioSourceTest(bsp, HasCountry);
  }
  if (cfg->conf_list[eBarcodeTest_SpecimenVoucher])
  {
    res->failed_tests[eBarcodeTest_SpecimenVoucher] = BarcodeBioSourceTest(bsp, HasSpecimenVoucher);
  }

  if (cfg->conf_list[eBarcodeTest_PercentN])
  {
    res->n_percent = PercentNInBioseq (bsp);
    res->failed_tests[eBarcodeTest_PercentN] = (Boolean)(res->n_percent > cfg->min_n_percent);
  }
 
  return res;
}


static void FailedBarcodeTests (BioseqPtr bsp, Pointer userdata)
{
  BarcodeBioseqSearchPtr sp;
  BarcodeTestResultsPtr  res = NULL;

  if (bsp == NULL || ISA_aa (bsp->mol) || !HasBARCODEKeyword (bsp) || userdata == NULL)
  {
    return;
  }

  sp = (BarcodeBioseqSearchPtr) userdata;
  if (sp->cfg == NULL) return;

  res = BarcodeTestResultsForBioseq (bsp, sp->cfg);
  if (res == NULL) return;

  if (sp->collect_positives
      || !PassBarcodeTests(res))
  {
    ValNodeAddPointer (&(sp->results_list), 0, res);
  }
  else
  {
    res = BarcodeTestResultsFree (res);
  }
}


extern ValNodePtr GetBarcodeFailedAccessionList (SeqEntryPtr sep, BarcodeTestConfigPtr cfg)
{
  BarcodeBioseqSearchData sd;

  if (cfg == NULL)
  {
    sd.cfg = BarcodeTestConfigNew();
  }
  else
  {
    sd.cfg = cfg;
  }

  sd.results_list = NULL;
  sd.collect_positives = FALSE;

  VisitBioseqsInSep (sep, &sd, FailedBarcodeTests);

  if (sd.cfg != cfg)
  {
    sd.cfg = BarcodeTestConfigFree (sd.cfg);
  }
  return sd.results_list;
}


extern ValNodePtr GetBarcodePassFail (SeqEntryPtr sep, BarcodeTestConfigPtr cfg)
{
  BarcodeBioseqSearchData sd;

  if (cfg == NULL)
  {
    sd.cfg = BarcodeTestConfigNew();
  }
  else
  {
    sd.cfg = cfg;
  }

  sd.results_list = NULL;
  sd.collect_positives = TRUE;

  VisitBioseqsInSep (sep, &sd, FailedBarcodeTests);

  if (sd.cfg != cfg)
  {
    sd.cfg = BarcodeTestConfigFree (sd.cfg);
  }
  return sd.results_list;
}


/* Report lists each Bioseq and whether the Bioseq passed all tests
 * or failed at least one.
 */
extern void WriteBarcodeTestCompliance (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               barcode_id, genbank_id;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    barcode_id = BarcodeTestBarcodeIdString (res->bsp);
    genbank_id = BarcodeTestGenbankIdString (res->bsp);
    fprintf (fp, "%s\t%s\t%s\n", barcode_id, genbank_id, 
                                 PassBarcodeTests (res) ? "PASS" : "FAIL");
    barcode_id = MemFree (barcode_id);
    genbank_id = MemFree (genbank_id);
  }
}


/* Report lists each Bioseq and whether the Bioseq passed all tests
 * or failed at least one.
 */
extern void WriteBarcodeTestComprehensive (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               barcode_id, genbank_id, reason;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    barcode_id = BarcodeTestBarcodeIdString (res->bsp);
    genbank_id = BarcodeTestGenbankIdString (res->bsp);
    reason = GetBarcodeTestFailureReasons (res);
    fprintf (fp, "%s\t%s\t%s\t%s\n", barcode_id, genbank_id, 
                                 PassBarcodeTests (res) ? "PASS" : "FAIL",
                                 reason == NULL ? "" : reason);
    barcode_id = MemFree (barcode_id);
    genbank_id = MemFree (genbank_id);
    reason = MemFree (reason);
  }
}


/* Create a tag table for updates */
extern void WriteBarcodeTagTable (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               barcode_id, genbank_id;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    barcode_id = BarcodeTestBarcodeIdString (res->bsp);
    genbank_id = BarcodeTestGenbankIdString (res->bsp);
    fprintf (fp, "%s\t%s\t\n", genbank_id, barcode_id);
    barcode_id = MemFree (barcode_id);
    genbank_id = MemFree (genbank_id);
  }
}


/* Report lists the individual tests that each Bioseq failed. */
extern void WriteBarcodeDiscrepancies (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               msg;

  if (fp == NULL) return;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    
    msg = SummaryTextFromBarcodeTestResults (res);
    fprintf (fp, "%s\n", msg);
    msg = MemFree (msg);
  }
}


/* Removes Barcode tech from all Bioseqs in supplied list */
extern void RemoveBarcodeKeywords (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               msg;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (res != NULL && res->bsp != NULL)
    {
      if (RemoveBarcodeTechFromBioseq (res->bsp))
      {
        if (fp != NULL)
        {
          msg = SummaryTextFromBarcodeTestResults (res);
          fprintf (fp, "%s\n", msg);
          msg = MemFree (msg);
        }
      }
    }
  }
}  


/* Applies Barcode technique to all Bioseqs in supplied list */
/* Used by Barcode Discrepancy Tool for the UNDO button.     */
extern void ApplyBarcodeKeywords (FILE *fp, ValNodePtr results_list)
{
  BarcodeTestResultsPtr res;
  ValNodePtr            vnp;
  CharPtr               msg;

  for (vnp = results_list; vnp != NULL; vnp = vnp->next)
  {
    res = (BarcodeTestResultsPtr) vnp->data.ptrvalue;
    if (res != NULL && res->bsp != NULL)
    {
      ApplyBarcodeTechToBioseq (res->bsp);
      if (fp != NULL)
      {
        msg = SummaryTextFromBarcodeTestResults (res);
        fprintf (fp, "%s\n", msg);
        msg = MemFree (msg);
      }
    }
  }
}  
  

#if defined (WIN32)
extern char * __stdcall AbstractReadFunction (Pointer userdata)
#else
extern char * AbstractReadFunction (Pointer userdata)
#endif
{
  ReadBufferPtr rbp;

  if (userdata == NULL) return NULL;

  rbp = (ReadBufferPtr) userdata;

  return MyFGetLine (rbp->fp, &(rbp->current_data));
}

#if defined (WIN32)
extern void __stdcall AbstractReportError (
#else
extern void AbstractReportError (
#endif
  TErrorInfoPtr err_ptr,
  Pointer      userdata
)
{
  TErrorInfoPtr PNTR list;
  TErrorInfoPtr last;

  if (err_ptr == NULL || userdata == NULL) return;

  list = (TErrorInfoPtr PNTR) userdata;

  if (*list == NULL)
  {
    *list = err_ptr;
  }
  else
  {
    for (last = *list; last != NULL && last->next != NULL; last = last->next)
    {}
    last->next = err_ptr;
  }

}



extern Boolean ParseLatLon (
  CharPtr lat_lon,
  FloatHi PNTR latP,
  FloatHi PNTR lonP
)

{
  char    ew;
  double  lat;
  double  lon;
  char    ns;

  if (latP != NULL) {
    *latP = 0.0;
  }
  if (lonP != NULL) {
    *lonP = 0.0;
  }

  if (StringHasNoText (lat_lon)) return FALSE;

  if (sscanf (lat_lon, "%lf %c %lf %c", &lat, &ns, &lon, &ew) == 4) {
    if (lon < -180.0) {
      lon = -180.0;
    }
    if (lat < -90.0) {
      lat = -90.0;
    }
    if (lon > 180.0) {
      lon = 180.0;
    }
    if (lat > 90.0) {
      lat = 90.0;
    }
    if (ew == 'W') {
      lon = -lon;
    }
    if (ns == 'S') {
      lat = -lat;
    }

    if (latP != NULL) {
      *latP = (FloatHi) lat;
    }
    if (lonP != NULL) {
      *lonP = (FloatHi) lon;
    }

    return TRUE;
  }

  return FALSE;
}

static CharPtr print_lat_lon_fmt = "%.*f %c %.*f %c";

extern void IsCorrectLatLonFormat (CharPtr lat_lon, BoolPtr format_correct, BoolPtr lat_in_range, BoolPtr lon_in_range)
{
  FloatLo  ns, ew;
  Char     lon, lat;
  Boolean  format_ok = FALSE, lat_ok = FALSE, lon_ok = FALSE;
  Int4     processed;

  if (StringHasNoText (lat_lon))
  {
    format_ok = FALSE;
  }
  else if (sscanf (lat_lon, "%f %c %f %c%n", &ns, &lat, &ew, &lon, &processed) != 4
           || processed != StringLen (lat_lon))
  {
    format_ok = FALSE;
  }
  else if ((lat != 'N' && lat != 'S') || (lon != 'E' && lon != 'W'))
  {
    format_ok = FALSE;
  }
  else 
  {
    format_ok = TRUE;
    if (ns <= 90 && ns >= 0)
    {
      lat_ok = TRUE;
    }
    if (ew <= 180 && ew >= 0)
    {
      lon_ok = TRUE;
    }
  }

  if (format_correct != NULL)
  {
    *format_correct = format_ok;
  }
  if (lat_in_range != NULL)
  {
    *lat_in_range = lat_ok;
  }
  if (lon_in_range != NULL)
  {
    *lon_in_range = lon_ok;
  }
}


static Boolean IsDirectionChar (Char dir)
{
  if (dir == 'E' || dir == 'W' || dir == 'N' || dir == 'S')
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static CharPtr MakeLatLonFromParts (FloatHi lat, Char ns, Int4 prec1, FloatHi lon, Char ew, Int4 prec2)
{
  CharPtr buf;

  buf = (CharPtr) MemNew (sizeof(Char) * 50);

  /* choose default directions when none supplied */
  if (ns == 0 && ew == 0)
  {
    ns = 'N';
    ew = 'E';
  }
  else if (ns == 0)
  {
    if (ew == 'E' || ew == 'W') 
    {
      ns = 'N';
    }
    else
    {
      ns = 'E';
    }
  }
  else if (ew == 0)
  {
    if (ns == 'N' || ns == 'S') 
    {
      ew = 'E';
    }
    else 
    {
      ew = 'N';
    }
  }

  /* correct -E to +W, -W to +W, -N to +S, -S to +S */
  if (lat < 0.0) 
  {
    if (ns == 'E') 
    {
      ns = 'W';
    }
    else if (ns == 'N')
    {
      ns = 'S';
    }
    lat = 0.0 - lat;
  }

  if (lon < 0.0) 
  {
    if (ew == 'E') 
    {
      ew = 'W';
    }
    else if (ew == 'N')
    {
      ew = 'S';
    }
    lon = 0.0 - lon;
  }

  if (ns == 'E' || ns == 'W')
  {
    sprintf (buf, print_lat_lon_fmt, prec2, lon, ew, prec1, lat, ns);
  }
  else
  {
    sprintf (buf, print_lat_lon_fmt, prec1, lat, ns, prec2, lon, ew);
  }
  return buf;
}


static Boolean ParseNumericFromDToken (CharPtr dtoken, FloatHiPtr val, Int4Ptr prec)
{
  FloatLo  a, b, c;
  FloatHi  f = 0.0;
  Int4     i, j, k;
  Boolean  rval = FALSE;
  Int4     processed, len, dec_size;
  CharPtr  cp;

  if (StringHasNoText (dtoken) || val == NULL || prec == NULL)
  {
    return FALSE;
  }

  len = StringLen (dtoken);
  if ((sscanf (dtoken, "%d.%d.%d%n", &i, &j, &k, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%d.%d.%d'%n", &i, &j, &k, &processed) == 3 && processed == len))
  {
    f = (FloatHi) i + (FloatHi)j / (FloatHi)60.0 + (FloatHi)k / (FloatHi)3600.0;
    *prec = 4;
    rval = TRUE;
  }
  else if ((sscanf (dtoken, "%f:%f:%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f:%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f %f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f'%f''%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f'%f\"%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f'%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f'%f'%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f'%f'%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f'%f'%f'%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f-%f-%f%n", &a, &b, &c, &processed) == 3 && processed == len)
      || (sscanf (dtoken, "%f %f-%f%n", &a, &b, &c, &processed) == 3 && processed == len))
  {
    f = a +  b / (FloatHi)60.0 + c / (FloatHi)3600.0;
    *prec = 4;
    rval = TRUE;
  }
  else if ((sscanf (dtoken, "%f %f%n", &a, &b, &processed) == 2 && processed == len)
           || (sscanf (dtoken, "%f:%f%n", &a, &b, &processed) == 2 && processed == len)
           || (sscanf (dtoken, "%f %f'%n", &a, &b, &processed) == 2 && processed == len)
           || (sscanf (dtoken, "%f'%f'%n", &a, &b, &processed) == 2 && processed == len)
           || (sscanf (dtoken, "%f'%f%n", &a, &b, &processed) == 2 && processed == len))
  {
    if (a < 0)
    {
      f = (FloatHi) a - b / (FloatHi) 60.0;
    }
    else
    {
      f = (FloatHi) a + b / (FloatHi) 60.0;
    }
    cp = StringRChr (dtoken, '.');
    if (cp == NULL)
    {
      dec_size = 0;
    }
    else
    {
      dec_size = StringLen (StringChr (dtoken, '.') + 1);
      if (dtoken[StringLen(dtoken) - 1] == '\'')
      {
        dec_size--;
      }
    }

    *prec = 2 + dec_size;
    rval = TRUE;
  }
  else if ((sscanf (dtoken, "%f%n", &a, &processed) == 1 && processed == len)
           || (sscanf (dtoken, "%f'%n", &a, &processed) == 1 && processed == len))
  {
    rval = TRUE;
    cp = StringChr (dtoken, '.');
    if (cp == NULL)
    {
      *prec = 2;
    }
    else
    {
      *prec = MAX (2, StringLen (cp + 1));
    }
    f = (FloatHi) a;
  }

  if (rval)
  {
    *val = f;
  } 
  return rval;
}

static Boolean ParseFromDToken (CharPtr dtoken, FloatHiPtr val, CharPtr d, Int4Ptr prec)
{
  FloatHi f;
  Char    dir = 0;
  Boolean rval = FALSE;
  Int4    token_len;

  if (StringHasNoText (dtoken) || val == NULL || d == NULL) 
  {
    return FALSE;
  }

  token_len = StringLen (dtoken);

  if (IsDirectionChar (dtoken[0])) 
  {
    dir = dtoken[0];
    rval = ParseNumericFromDToken (dtoken + 1, &f, prec);
  }
  else if (IsDirectionChar (dtoken[token_len - 1]))
  {
    dir = dtoken[token_len - 1];
    dtoken[token_len - 1] = 0;
    token_len--;
    while (token_len > 0 && isspace (dtoken[token_len - 1]))
    {
      dtoken[token_len - 1] = 0;
      token_len --; 
    }
    rval = ParseNumericFromDToken (dtoken, &f, prec);
  }
  else
  {
    rval = ParseNumericFromDToken (dtoken, &f, prec);
  }
  if (rval)
  {
    *val = f;
    *d = dir;
  }
  return rval;
}


static Boolean ParseFromLToken (CharPtr ltoken, Boolean first, FloatHiPtr val, CharPtr d, Int4Ptr prec)
{
  CharPtr  dtoken;
  Boolean  rval = FALSE;
  FloatHi  f;
  Char     dir;
  Char     plus_dir, minus_dir;
  Int4     len;

  if (StringHasNoText (ltoken) || val == NULL || d == NULL)
  {
    return rval;
  }
  len = StringLen (ltoken);
  if (StringNCmp (ltoken, "LAT", 3) == 0)
  {
    dtoken = ltoken + 3;
    plus_dir = 'N';
    minus_dir = 'S';
  } 
  else if (StringNCmp (ltoken, "LONG", 4) == 0)
  {
    dtoken = ltoken + 4;
    plus_dir = 'E';
    minus_dir = 'W';
  }
  else if (len > 3 && StringCmp (ltoken + len - 3, "LAT") == 0)
  {
    ltoken[len - 3] = 0;
    dtoken = ltoken;
    plus_dir = 'N';
    minus_dir = 'S';
  }
  else if (len > 4 && StringCmp (ltoken + len - 4, "LONG") == 0)
  {
    ltoken[len - 4] = 0;
    dtoken = ltoken;
    plus_dir = 'E';
    minus_dir = 'W';
  }
  else if (first)
  {
    dtoken = ltoken;
    plus_dir = 'N';
    minus_dir = 'S';
  }
  else
  {
    dtoken = ltoken;
    plus_dir = 'E';
    minus_dir = 'W';
  }
  /* trim space and punctuation from beginning */
  while (isspace (*dtoken) || (*dtoken != '-' && ispunct(*dtoken))) 
  {
    dtoken++;
  }
  /* trim space from end */
  len = StringLen (dtoken);
  while (len > 0 && isspace (dtoken[len - 1])) {
    dtoken[len - 1] = 0;
    len--;
  }
  if (ParseFromDToken (dtoken, &f, &dir, prec))
  {
    if (dir == 0)
    {
       if (f < 0)
       {
         dir = minus_dir;
         f = 0 - f;
       }
       else
       {
         dir = plus_dir;
       }
       rval = TRUE;
    }
    else if (dir == plus_dir || dir == minus_dir)
    {
      rval = TRUE;
    }
  }

  if (rval)
  {
    *val = f;
    *d = dir;
  }
  return rval;

}


static CharPtr MakeToken(CharPtr token1, CharPtr token2)
{
  Int4    token_len;
  CharPtr token;

  if (StringHasNoText (token1)) return NULL;
  while (isspace (*token1) || (ispunct (*token1) && *token1 != '-'))
  {
    token1++;
  }
  if (*token1 == 0)
  {
    return NULL;
  }
  if (token2 == NULL)
  {
    token_len = StringLen (token1) + 1;
  }
  else
  {
    token_len = token2 - token1 + 1;
  }
  token = (CharPtr) MemNew (sizeof (Char) * token_len);
  strncpy (token, token1, token_len - 1);
  token[token_len - 1] = 0;
  while ((isspace (token[token_len - 2]) || ispunct (token[token_len - 2])) && token_len > 2)
  {
    token[token_len - 2] = 0;
    token_len--;
  }
  return token;
}


typedef struct replacepair {
  CharPtr find;
  CharPtr replace;
} ReplacePairData, PNTR ReplacePairPtr;

static ReplacePairData latlon_replace_list[] = {
 { "LONGITUDE", "LONG" },
 { "LONG.",     "LONG" },
 { "LON.",      "LONG" },
 { "LATITUDE",  "LAT"  },
 { "LAT.",      "LAT"  },
 { "DEGREES",   " " },
 { "DEGREE",    " " },
 { "DEG.",      " " },
 { "DEG",       " " },
 { "MIN.",      "'" },
 { "MIN",       "'" },
 { "SEC.",      "''" },
 { "SEC",       "''" },
 { "NORTH",     "N"    },
 { "SOUTH",     "S"    },
 { "EAST",      "E"    },
 { "WEST",      "W"    },
};

static Int4 num_latlon_replace = sizeof (latlon_replace_list) / sizeof (ReplacePairData);


static Boolean CommaShouldBePeriod (CharPtr pStr, CharPtr pComma)
{
  CharPtr cp;
  Boolean rval = FALSE;

  if (StringHasNoText (pStr) || pComma == NULL || pComma <= pStr) return FALSE;

  cp = pComma - 1;
  if (!isdigit (*cp)) return FALSE;

  while (cp > pStr && isdigit (*cp)) 
  {
    cp--;
  }
  if (*cp != '.' && isdigit (*(pComma + 1)))
  {
    rval = TRUE;
  }
  return rval;
}

extern CharPtr FixLatLonFormat (CharPtr orig_lat_lon)
{
  FloatHi lon, lat;
  Char    ns, ew;
  CharPtr cpy, cp, dst;
  CharPtr first_dash = NULL, second_dash = NULL, first_space = NULL, second_space = NULL;
  CharPtr rval = NULL;
  CharPtr word1 = NULL, word2 = NULL;
  Boolean bad_letter_found = FALSE, replace_found, comma_sep = FALSE;
  Int4    i;
  CharPtr ltoken1 = NULL, ltoken2 = NULL;
  CharPtr dtoken1 = NULL, dtoken2 = NULL;
  Int4    find_len, replace_len, jump_len;
  Int4    prec1, prec2;
  CharPtr extra_text = NULL;

  if (StringHasNoText (orig_lat_lon))
  {
    return NULL;
  }

  cpy = StringSave (orig_lat_lon);

  cp = cpy;
  while (*cp != 0)
  {
    if (*cp == 'O' 
        && (cp == cpy || !isalpha (*(cp - 1))))
    {
      *cp = '0';
    }
    else if (*cp == 'o' && cp != cpy && !isalpha (*(cp - 1)) && !isalpha (*(cp + 1)))
    {
      *cp = ' ';
    }
    else if (*cp == '#') /* # is sometimes used for degree, sometimes separator */
    {
      *cp = ' ';
    }
    else if (isalpha (*cp)) 
    {
      *cp = toupper (*cp);
    }
    else if (*cp == ',')
    {
      if (CommaShouldBePeriod (cpy, cp))
      {
        *cp = '.';
      }
    }
    cp++;
  }

  /* fix words */
  cp = cpy;
  dst = cpy;
  while (*cp != 0 && !bad_letter_found && extra_text == NULL)
  {
    if (isalpha (*cp)) 
    {
      replace_found = FALSE;
      for (i = 0; i < num_latlon_replace && !replace_found; i++) 
      {
        find_len = StringLen (latlon_replace_list[i].find);
        replace_len = StringLen (latlon_replace_list[i].replace);
        jump_len = 0;
        if (StringNICmp (cp, latlon_replace_list[i].find, find_len) == 0)
        {
          jump_len = find_len;
        }
        else if (StringNCmp (cp, latlon_replace_list[i].replace, replace_len) == 0)
        {
          jump_len = replace_len;
        }
        else
        {
          continue;
        }

        if (i < 5)
        {
          if (ltoken1 == NULL)
          {
            ltoken1 = dst;
          }
          else if (ltoken2 == NULL)
          {
            ltoken2 = dst;
          }
          else
          {
            bad_letter_found = TRUE;
          }
        }
        else if (i >= 13)
        {
          if (dtoken1 == NULL)
          {
            dtoken1 = dst;
          }
          else if (dtoken2 == NULL)
          {
            dtoken2 = dst;
          }
          else if (!comma_sep)
          {
            bad_letter_found = TRUE;
          }
        }

        if ((latlon_replace_list[i].replace[0] == '\'' || latlon_replace_list[i].replace[0] == ' ')
            && dst > cpy && isspace (*(dst - 1)))
        {
          /* no double spaces, no spaces before tick marks */ 
          dst--;
        }
        if (replace_len == 1 && latlon_replace_list[i].replace[0] == ' ')
        {
          if (isspace (*(cp + jump_len)) || *(cp + jump_len) == 0 || *(cp + jump_len) == '\'')
          {
            /* no double spaces */
          }
          else 
          {
            *dst = ' ';
            dst++;
          }
        }          
        else 
        {
          StringNCpy (dst, latlon_replace_list[i].replace, replace_len);
          dst += replace_len;
        }
        cp += jump_len;
        replace_found = TRUE;
      } 
      if (!replace_found)
      {
        bad_letter_found = 1;
      }
    }
    else if (isspace (*cp))
    {
      if (isspace (*(cp + 1)) || *(cp + 1) == 0 || *(cp + 1) == '\''
          || (dst > cpy && isspace (*(dst - 1))))
      {
        cp++;
      }
      else
      {
        *dst = ' ';
        if (first_space == NULL)
        {
          first_space = dst;
        }
        else if (second_space == NULL)
        {
          second_space = dst;
        }        
        dst++;
        cp++;
      }
    }
    else if (*cp == '-')
    {
      *dst = '-';
      if (first_dash == NULL)
      {
        first_dash = dst;
      }
      else if (second_dash == NULL)
      {
        second_dash = dst;
      }
      dst++;
      cp++;
    }
    else if (*cp == ',')
    {
      if (!comma_sep && ((ltoken1 != NULL && ltoken2 != NULL) || (dtoken1 != NULL && dtoken2 != NULL)))
      {
        extra_text = orig_lat_lon + (cp - cpy);
      } else {  
        *dst = ' ';
        dst++;
        cp++;
        if (dtoken1 != NULL && dtoken2 == NULL)
        {
          dtoken2 = dst;
          comma_sep = TRUE;
        }
        else if (ltoken1 == NULL && ltoken2 == NULL)
        {
          ltoken1 = cpy;
          ltoken2 = dst;
          comma_sep = TRUE;
        }
      }
    }
    else if (*cp == ';' && cp > cpy && isdigit(*(cp - 1)) && isdigit(*(cp + 1)))
    {
      /* replace typo semicolon with colon */
      *dst = ':';
      dst++;
      cp++;
    }
    else
    {
      *dst = *cp;
      dst++;
      cp++;
    }
  }

  *dst = 0;

  /* have to have both ltokens or none */
  if (ltoken1 != NULL && ltoken2 == NULL)
  {
    bad_letter_found = 1;
  }
  /* if no ltokens, must have both dtokens */
  else if (ltoken1 == NULL && (dtoken1 == NULL || dtoken2 == NULL))
  { 
    /* use space to separate the two tokens */
    if (first_space != NULL && second_space == NULL)
    {
      ltoken1 = cpy;
      ltoken2 = first_space + 1;
    }
    /* allow a dash to separate the two tokens if no spaces and only one dash */
    else if (first_space == NULL && second_space == NULL && first_dash != NULL && second_dash == NULL)
    {
      ltoken1 = cpy;
      *first_dash = ' ';
      ltoken2 = first_dash + 1;
    }
    else if (dtoken1 != NULL && dtoken2 == NULL && dtoken1 > cpy && dtoken1 < cpy + StringLen (cpy) - 1)
    {
      word1 = MakeToken (cpy, dtoken1 + 1);
      if (ParseFromDToken (word1, &lat, &ns, &prec1))
      {
        /* first portion parses ok, assume user just left off direction for second token */
        /* letters end tokens */
        dtoken2 = dtoken1 + 1;
        dtoken1 = cpy;
      }
      else
      {
        bad_letter_found = 1;
      }
      word1 = MemFree (word1);
    }
    else
    {
      bad_letter_found = 1;
    }
  }
  if (first_space == NULL && first_dash != NULL && second_dash == NULL && !comma_sep)
  {
    /* don't let the dash dividing the tokens be used as minus sign */
    *first_dash = ' ';
  }

  if (bad_letter_found)
  {
  }
  else if (ltoken1 != NULL)
  {
    /* if latitude and longitude are at end of token, change start */
    if (ltoken1 != cpy)
    {
      ltoken2 = ltoken1 + 3;
      if (*ltoken2 == 'G') 
      {
        ltoken2++;
      }
      ltoken1 = cpy;
    }
    word1 = MakeToken(ltoken1, ltoken2);
    word2 = MakeToken(ltoken2, NULL);
    if (ParseFromLToken (word1, TRUE, &lat, &ns, &prec1)
        && ParseFromLToken (word2, FALSE, &lon, &ew, &prec2))
    {
      rval = MakeLatLonFromParts (lat, ns, prec1, lon, ew, prec2);
    }
  }
  else
  {
    if (dtoken1 != cpy) 
    {
      /* letters end tokens */
      dtoken2 = dtoken1 + 1;
      dtoken1 = cpy;
    }
    word1 = MakeToken (dtoken1, dtoken2);
    word2 = MakeToken (dtoken2, NULL);
    if (ParseFromDToken (word1, &lat, &ns, &prec1)
        && ParseFromDToken (word2, &lon, &ew, &prec2))
    {
      rval = MakeLatLonFromParts (lat, ns, prec1, lon, ew, prec2);
    }
  }
      
  word1 = MemFree (word1);
  word2 = MemFree (word2);
  cpy = MemFree (cpy);
  
  if (rval != NULL && extra_text != NULL)
  {
    cpy = (CharPtr) MemNew (sizeof (Char) * (StringLen (rval) + StringLen (extra_text) + 1));
    sprintf (cpy, "%s%s", rval, extra_text);
    rval = MemFree (rval);
    rval = cpy;
  }
  return rval;
}


static void TestLatLonFormatting (FILE *fp)
{
  CharPtr tests[]  = 
  { "100.12 N 200.12 E",     /* already correct */
    "100 N 200 E",           /* correctable */
    "100.1 N 200.2 E",       /* correctable */
    "1OO.1 N 200.2 E",       /* correctable (replace capital o with zero) */
    "100.1 N, 200.2 E",      /* correctable (remove comma) */
    "E 100, S 120",          /* correctable (remove comma, reverse order, letters before numbers */
    "latitude: 200 N longitude: 100 E",
    "latitude: 200 E longitude: 100 N", /* NOT correctable */
    "N 37 45.403', 119 1.456' W",
    "38 52 56 N 84 44 53 W",
    "49 29 50 N 80 25 52 W",
    "39N 93W",
    "42:43:13N 01:0015W",
    "02deg 33min 00.7sec S 45deg 01min 38.8sec W",
    "42:24:37.9 N 85:22:11.7 W",
    "10 N 124 E",
    "41deg30'' S 145deg37' E",
    "59.30deg N 22.40deg E",
    "35 N 134 E",
    "2 S 114 E",
    "24deg 24.377' N 101deg 23.073' W'",
    "26deg 57.9' N 102deg 08.3 W'",
    "38 11 44.66 North 0 35 01.93 West",
    "62.08 N 129.682",
    "64.444 N -164.973",
    "62.033 N -146.533",
    "67 N -51",
    "69.107 N 124.195",
    "2:46:00-59:41:00",
    "64 degree 55 N 25 degree 05 E",
    "64.907 N -166.18",
    "2:46:00-59:41:00",
    "66 degree 21 N 29 degree 21 E",
    "37deg27N 121deg52'W",
    "01deg31'25''N 66''33'31''W",
    "07deg33'30''N 69deg20'W",
    "10.8439,-85.6138",
    "11.03,-85.527",
    "8 deg 45 min S, 63 deg 26 min W",
    "29deg 49' 23.7' N; 106deg 23' 15.8'W",
    "7:46S, 12:30E",
    "35deg48'50'' N; 82deg5658'' W",
    "45deg34.18''N, 122deg12.00 'W",
    "37deg27N, 121deg52'W",
    "41:00;00N 20:45:00E",
    "02 deg 28' 29# S, 56 deg 6' 31# W"
};
  Int4 test_num, num_tests = sizeof (tests) / sizeof (char *);
  CharPtr fix;
  Int4 num_pass = 0, num_formatted = 0;
  Boolean format_ok, lat_in_range, lon_in_range;

  if (fp == NULL) return;

  for (test_num = 0; test_num < num_tests; test_num++)
  {
    fprintf (fp, "Test %d: %s\n", test_num, tests[test_num]);
    fix = FixLatLonFormat (tests[test_num]);
    if (fix == NULL) 
    {
      fprintf (fp, "Unable to correct format\n");
    }
    else
    {
      IsCorrectLatLonFormat (fix, &format_ok, &lat_in_range, &lon_in_range);
      if (format_ok)
      {
        num_formatted ++;
        fprintf (fp, "Correction succeeded:%s\n", fix);
        num_pass++;
      }
      else
      {
        num_formatted ++;
        fprintf (fp, "Correction failed:%s\n", fix);
      }
    }
  }
  fprintf (fp, "Formats %d out of %d, %d succeed\n", num_formatted, num_tests, num_pass);
}


static CharPtr StringFromObjectID (ObjectIdPtr oip)
{
  CharPtr    str;
  if (oip == NULL) return NULL;

  if (oip->id > 0)
  {
    str = (CharPtr) MemNew (sizeof (Char) * 20);
    sprintf (str, "%d", oip->id);
  }
  else
  {
    str = StringSave (oip->str);
  }
  return str;
}

extern void ApplyBarcodeDbxrefToBioSource (BioSourcePtr biop, ObjectIdPtr oip)
{
  ValNodePtr vnp;
  DbtagPtr   dbt;
  CharPtr    str, cmp;
  Boolean    found = FALSE;

  if (biop == NULL || oip == NULL) return;

  if (biop->org == NULL)
  {
    biop->org = OrgRefNew();
  }

  str = StringFromObjectID (oip);

  for (vnp = biop->org->db; vnp != NULL && !found; vnp = vnp->next)
  {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL || dbt->tag == NULL) continue;
    if (StringCmp (dbt->db, "BOLD") != 0) continue;
    cmp = StringFromObjectID (dbt->tag);
    if (StringCmp (str, cmp) == 0) found = TRUE;
    cmp = MemFree (cmp);
  }
  if (found) 
  {
    str = MemFree (str);
  }
  else
  {
    dbt = DbtagNew ();
    dbt->db = StringSave ("BOLD");
    dbt->tag = ObjectIdNew();
    dbt->tag->str = str;
    ValNodeAddPointer (&(biop->org->db), 0, dbt);
  }
}


extern void ApplyBarcodeDbxrefsToBioseq (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr       sdp;
  SeqMgrDescContext context;
  SeqIdPtr          sip;
  DbtagPtr          dbt;

  if (bsp == NULL) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next)
  {
    if (IsBarcodeID (sip) && sip->choice == SEQID_GENERAL && sip->data.ptrvalue != NULL) 
    {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
      if (sdp != NULL)
      {
        ApplyBarcodeDbxrefToBioSource ((BioSourcePtr) sdp->data.ptrvalue, dbt->tag);
      }
    }
  }
}


/* Code for Country Fixup */
typedef struct  countrystatelist {
  CharPtr PNTR state_list;
  CharPtr country_name;
} CountryStateListData, PNTR CountryStateListPtr;

static Boolean IsMatchInSecondChoiceLists (CharPtr find_str, Int4 match_len, CountryStateListPtr second_choice_lists, CharPtr whole_string)
{
  Int4    i, j;
  Boolean in_lists = FALSE;
  CharPtr cp;
  Int4 len_second_choice;

  if (StringHasNoText (find_str) || match_len < 1 || second_choice_lists == NULL) return FALSE;

  for (i = 0; second_choice_lists[i].state_list != NULL && !in_lists; i++)
  {
    for (j = 0; second_choice_lists[i].state_list[j] != NULL && !in_lists; j++) 
    {
      len_second_choice = StringLen (second_choice_lists[i].state_list[j]);      
      if (len_second_choice == match_len 
          &&StringNCmp (find_str, second_choice_lists[i].state_list[j], match_len) == 0) 
      {
        in_lists = TRUE;
      }
      else if ((cp = StringSearch (second_choice_lists[i].state_list[j], find_str)) != NULL
               && StringSearch (whole_string, second_choice_lists[i].state_list[j]) != NULL)
      {
        in_lists = TRUE;
      }
    }
  }
  return in_lists;
}


static Boolean IsBodyOfWater (CharPtr str)
{
  if (StringHasNoText (str)) return FALSE;
  if (StringSearch (str, "Ocean") != NULL) return TRUE;
  if (StringSearch (str, "Gulf") != NULL) return TRUE;
  if (StringSearch (str, "Sea") != NULL) return TRUE;
  return FALSE;
}


static Boolean IsSubstringOfStringInList (CharPtr whole_str, CharPtr match_p, CharPtr match_str, CharPtr PNTR list)
{
  CharPtr cp;
  Int4    context_len, find_len;
  Boolean rval = FALSE;

  if (list == NULL || StringHasNoText (whole_str) || match_p == NULL || match_p < whole_str) {
    return FALSE;
  }
  find_len = StringLen (match_str);
  while (*list != NULL && !rval) {
    context_len = StringLen (*list);
    if (find_len < context_len) {
      cp = StringSearch (whole_str, *list);
      while (cp != NULL && !rval) {
        if (match_p < cp) {
          cp = NULL;
        } else if (cp + context_len > match_p) {
          rval = TRUE;
        } else {
          cp = StringSearch (cp + 1, *list);
        }
      }
    }
    list++;
  }
  return rval;
}


static CharPtr bad_context_names[] = {
 "Gibraltar Range National Park",
 "Western Australia",
 "WSW Chihuahua",
 "Mississippi River",
 NULL };

static Boolean IsBadContextName (CharPtr whole_str, CharPtr match_p, CharPtr match_str)
{
  Boolean rval;
  Int4    len;

  rval = IsSubstringOfStringInList (whole_str, match_p, match_str, bad_context_names);
  if (!rval) {
    len = StringLen (match_str);
    if (StringCmp (match_p + len, " River") == 0) {
      rval = TRUE;
    } else if (StringCmp (match_p + len, " State University") == 0) {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean IsPartOfStateName (CharPtr whole_str, CharPtr match_p, CharPtr match_str, CountryStateListPtr second_choice_lists)
{
  Boolean rval = FALSE;
  Int4    i;

  if (second_choice_lists == NULL) return FALSE;

  for (i = 0; second_choice_lists[i].state_list != NULL && !rval; i++) {    
    rval = IsSubstringOfStringInList (whole_str, match_p, match_str, second_choice_lists[i].state_list);
  }
  return rval;
}


static CharPtr
FindStringInStringWithContext
(CharPtr search_str, CharPtr look_for, CharPtr PNTR list, CountryStateListPtr second_choice_lists)
{
  Int4          len_match;
  CharPtr       cp;

  if (StringHasNoText (search_str) || StringHasNoText (look_for)) {
    return NULL;
  }

  cp = StringISearch (search_str, look_for);
  len_match = StringLen (look_for);
  while (cp != NULL) {
    /* if character after match is alpha, continue */
    if (isalpha ((Int4)(cp [len_match]))
        /* if character before match is alpha, continue */
        || (cp > search_str && isalpha ((Int4)(*(cp - 1))))
        /* if match is part of a known "bad context", continue */
        || IsBadContextName (search_str, cp, look_for)
        /* if is shorter match for other item, continue */
        || IsSubstringOfStringInList (search_str, cp, look_for, list)
        || IsPartOfStateName (search_str, cp, look_for, second_choice_lists)) {
      cp = StringSearch (cp + len_match, look_for);
    } else {
      return cp;
    }
  }
  return cp;
}
 

static ValNodePtr FindBestStringMatch (CharPtr PNTR list, CharPtr find_str, CountryStateListPtr second_choice_lists)
{
  CharPtr PNTR  ptr;
  Int4          len_find;
  CharPtr       cp;
  Boolean       ocean_best, ocean_this;
  Boolean       best_in_second, this_in_second;
  ValNodePtr    match_list = NULL, vnp, best_vnp;
  
  if (list == NULL || find_str == NULL) return NULL;
  
  len_find = StringLen (find_str);
  
  /* first, find all matches */
  for (ptr = list; ptr != NULL && *ptr != NULL; ptr++)
  {
    cp = FindStringInStringWithContext (find_str, *ptr, list, second_choice_lists);
    if (cp != NULL) {
      ValNodeAddPointer (&match_list, 1, *ptr);
    }
  } 

  if (match_list == NULL) return NULL;

  /* now eliminate matches where we know we have a preference */
  best_vnp = match_list;
  for (vnp = match_list->next; vnp != NULL; vnp = vnp->next) 
  {
    if (StringSearch (vnp->data.ptrvalue, best_vnp->data.ptrvalue) != NULL)
    {
      /* best is inside this one */
      best_vnp->choice = 0;
      best_vnp = vnp;
    } else if (StringSearch (best_vnp->data.ptrvalue, vnp->data.ptrvalue) != NULL) {
      /* this is inside best */
      vnp->choice = 0;
    } else {
      /* prefer non-ocean to ocean */
      ocean_best = IsBodyOfWater (best_vnp->data.ptrvalue);
      ocean_this = IsBodyOfWater (vnp->data.ptrvalue);
      if (ocean_this && !ocean_best)
      {
        /* disregard this one */
        vnp->choice = 0;
        continue;
      }
      else if (!ocean_this && ocean_best) 
      {
        /* definitely take this to replace best */
        best_vnp->choice = 0;
        best_vnp = vnp;     
      }
      else if (second_choice_lists != NULL)
      {
        best_in_second = IsMatchInSecondChoiceLists (best_vnp->data.ptrvalue, StringLen (best_vnp->data.ptrvalue), second_choice_lists, find_str);
        this_in_second = IsMatchInSecondChoiceLists (vnp->data.ptrvalue, StringLen (vnp->data.ptrvalue), second_choice_lists, find_str);
        /* if this choice is a second choice, but the previous best wasn't, don't bother with this */
        if (this_in_second && !best_in_second) {
          vnp->choice = 0;
        } else if (!this_in_second && best_in_second) {
          /* if previous choice was in the secondary lists, prefer this and ignore previous */
          best_vnp->choice = 0;
          best_vnp = vnp;
        }
      }
    }
  }
  vnp = ValNodeExtract (&match_list, 0);
  vnp = ValNodeFree (vnp);
  return match_list;
}

static CharPtr usa_state_list[] = 
{
  "Alabama",
  "Alaska",
  "Arizona",
  "Arkansas",
  "California",
  "Colorado",
  "Connecticut",
  "Delaware",
  "Florida",
  "Georgia",
  "Hawaii",
  "Idaho",
  "Illinois",
  "Indiana",
  "Iowa",
  "Kansas",
  "Kentucky",
  "Louisiana",
  "Maine",
  "Maryland",
  "Massachusetts",
  "Michigan",
  "Minnesota",
  "Mississippi",
  "Missouri",
  "Montana",
  "Nebraska",
  "Nevada",
  "New Hampshire",
  "New Jersey",
  "New Mexico",
  "New York",
  "North Carolina",
  "North Dakota",
  "Ohio",
  "Oklahoma",
  "Oregon",
  "Pennsylvania",
  "Rhode Island",
  "South Carolina",
  "South Dakota",
  "Tennessee",
  "Texas",
  "Utah",
  "Vermont",
  "Virginia",
  "Washington", 
  "Washington, DC",
  "West Virginia",
  "Wisconsin",
  "Wyoming",
  NULL
};

static CharPtr uk_state_list[] = {
  "England",
  "Scotland",
  NULL
};

static CharPtr canada_province_list[] = {
  "Alberta",
  "British Columbia",
  "Manitoba",
  "New Brunswick",
  "Newfoundland and Labrador",
  "Northwest Territories",
  "Nova Scotia",
  "Nunavut",
  "Ontario",
  "Prince Edward Island",
  "Quebec",
  "Saskatchewan",
  "Yukon",
  NULL
};

static CharPtr australia_state_list[] = {
  "Australian Capital Territory",
  "Jervis Bay Territory",
  "New South Wales",
  "Northern Territory",
  "Queensland",
  "South Australia",
  "Tasmania",
  "Victoria",
  "Western Australia",
  NULL
};

static CharPtr mx_state_list[] = 
{
  "Aguascalientes",
  "Baja California",
  "Baja California Sur",
  "Campeche",
  "Chiapas",
  "Chihuahua",
  "Coahuila",
  "Colima",
  "Distrito Federal",
  "Durango",
  "Estado de Mexico",
  "Guanajuato",
  "Guerrero",
  "Hidalgo",
  "Jalisco",
  "Michoacan",
  "Morelos",
  "Nayarit",
  "Nuevo Leon",
  "Oaxaca",
  "Puebla",
  "Queretaro",
  "Quintana Roo",
  "San Luis Potosi",
  "Sinaloa",
  "Sonora",
  "Tabasco",
  "Tamaulipas",
  "Tlaxcala",
  "Veracruz",
  "Yucatan",
  "Zacatecas",
  NULL
};

static CharPtr portugal_state_list[] = {
  "Azores",
  NULL
};

static CharPtr ecuador_state_list[] = {
  "Galapagos",
  NULL
};


static CountryStateListData country_state_list[] = {
{ usa_state_list, "USA" },
{ uk_state_list, "United Kingdom"},
{ canada_province_list, "Canada"},
{ australia_state_list, "Australia"},
{ mx_state_list, "Mexico"},
{ portugal_state_list, "Portugal"},
{ ecuador_state_list, "Ecuador"},
{ NULL, NULL}
};


static CharPtr FindStateMatch (CharPtr search, CharPtr PNTR country, Int4Ptr state_len, BoolPtr pMulti)
{
  CharPtr best_match = NULL;
  Int4    i;
  ValNodePtr state_matches;

  if (StringHasNoText (search)) return NULL;
  if (country != NULL) {
    *country = NULL;
  }

  for (i = 0; country_state_list[i].state_list != NULL; i++) {    
    state_matches = FindBestStringMatch (country_state_list[i].state_list, search, NULL);
    if (state_matches != NULL) {
      if (state_matches->next == NULL && best_match == NULL) {
        best_match = state_matches->data.ptrvalue;
        *country = country_state_list[i].country_name;
      } else {
        *pMulti = TRUE;
        return NULL;
      }
      state_matches = ValNodeFree (state_matches);
    }
  }
  if (best_match != NULL && state_len != NULL) {
    *state_len = StringLen (best_match);
  }
  return best_match;
}


static CharPtr FindStateMatchForCountry (CharPtr search, CharPtr country, BoolPtr pMulti)
{
  ValNodePtr state_matches;
  CharPtr state_match = NULL;
  Int4    i;

  if (StringHasNoText (search) || StringHasNoText (country)) return NULL;

  for (i = 0; country_state_list[i].state_list != NULL; i++) {    
    if (StringCmp (country, country_state_list[i].country_name) == 0) {
      state_matches = FindBestStringMatch (country_state_list[i].state_list, search, NULL);
      if (state_matches != NULL) {
        if (state_matches->next == NULL) {
          state_match = state_matches->data.ptrvalue;
        } else {
          *pMulti = TRUE;
        }
        state_matches = ValNodeFree (state_matches);
        return state_match;
      }
    }
  }
  return NULL;
}


static void FixCountryStringForStateName (CharPtr PNTR pCountry, CharPtr state_name, CharPtr country_name)
{
  CharPtr cp;
  Int4    len_state, len_country, len_qual, len_name, len_after;
  CharPtr before, newname;
  
  if (pCountry == NULL
      || StringHasNoText (*pCountry)
      || StringHasNoText (state_name)
      || StringHasNoText (country_name))
  {
    return;
  }
  
  cp = StringStr (*pCountry, state_name);
  if (cp == NULL)
  {
    return;
  }
  len_state = StringLen (state_name);
  if (isalpha ((Int4)(cp [len_state])))
  {
    return;
  }
  
  len_country = StringLen (country_name);
  
  len_qual = StringLen (*pCountry);
 	if (cp == *pCountry)
  {
    len_after = len_qual - len_state;
    newname = (CharPtr) MemNew ((5 + len_country + len_state + len_after) * sizeof (Char));
    sprintf (newname, "%s: %s", country_name, state_name);
    if (len_after > 0)
    {
      StringCat (newname, ", ");
      StringCat (newname, *pCountry + len_state);
    }
    *pCountry = MemFree (*pCountry);
    *pCountry = newname;
  }
  else
  {
    newname = (CharPtr) MemNew (len_qual + 5 + len_country);
    *(cp - 1) = 0;
    before = StringSave (*pCountry);
    sprintf (newname, "%s: %s, ", country_name, state_name);
    StringNCpy (newname + 4 + len_country + len_state, before, StringLen (before));
    StringCpy (newname + 4 + len_country + len_state + StringLen (before), cp + len_state);
    len_name = StringLen (newname);
    while (isspace ((Int4)(newname[len_name - 1])) || ispunct ((Int4)(newname [len_name - 1])))
    {
    	newname [len_name - 1] = 0;
   	  len_name --;
    }
    /* get rid of trailing comma if necessary */
    if (len_name == StringLen (country_name) + 3 + StringLen (state_name)
        && newname [len_name - 1] == ',')
    {
    	newname [len_name - 1] = 0;
   	  len_name --;
    }
    before = MemFree (before);
    MemFree (*pCountry);
    *pCountry = newname;
  }
}


static CharPtr FindCountryMatch (CharPtr search_str, CharPtr PNTR country_list, BoolPtr isMulti)
{
  ValNodePtr match_list;
  CharPtr best_match = NULL;

  if (StringSearch (search_str, "Yugoslavia")) {
    *isMulti = TRUE;
    return NULL;
  }

  match_list = FindBestStringMatch (country_list, search_str, country_state_list);
  if (match_list != NULL) {
    if (match_list->next == NULL) {
      best_match = match_list->data.ptrvalue;
    } else {
      *isMulti = TRUE;
    }
    match_list = ValNodeFree (match_list);
  }
  return best_match;
}


static ReplacePairData country_name_fixes[] = {
 {"Vietnam", "Viet Nam"},
 {"Ivory Coast", "Cote d'Ivoire"},
 {"United States of America", "USA"},
 {"U.S.A.", "USA"},
 {"The Netherlands", "Netherlands"},
 {NULL, NULL}
};

static void FixCountryNames (CharPtr PNTR pCountry)
{
  ReplacePairPtr fix;

  if (pCountry == NULL || StringHasNoText (*pCountry))
  {
    return;
  }

  fix = country_name_fixes;
  while (fix->find != NULL) 
  {
    if (StringStr (*pCountry, fix->replace) == NULL || StringSearch (fix->find, fix->replace) != NULL) {
      FindReplaceString (pCountry, fix->find, fix->replace, FALSE, TRUE);
    }
    fix++;
  }
}


static ReplacePairData us_state_abbrev_fixes[] = {
 {"AL", "Alabama"},
 {"AK", "Alaska"},
 {"AZ", "Arizona"},
 {"AR", "Arkansas"},
 {"CA", "California"},
 {"CO", "Colorado"},
 {"CN", "Connecticut"},
 {"DE", "Delaware"},
 {"FL", "Florida"},
 {"GA", "Georgia"},
 {"HI", "Hawaii"},
 {"ID", "Idaho"},
 {"IL", "Illinois"},
 {"IN", "Indiana"},
 {"IA", "Iowa"},
 {"KS", "Kansas"},
 {"KY", "Kentucky"},
 {"LA", "Louisiana"},
 {"ME", "Maine"},
 {"MD", "Maryland"},
 {"MA", "Massachusetts"},
 {"MI", "Michigan"},
 {"MN", "Minnesota"},
 {"MS", "Mississippi"},
 {"MO", "Missouri"},
 {"MT", "Montana"},
 {"NE", "Nebraska"},
 {"NV", "Nevada"},
 {"NH", "New Hampshire"},
 {"NJ", "New Jersey"},
 {"NM", "New Mexico"},
 {"NY", "New York"},
 {"NC", "North Carolina"},
 {"ND", "North Dakota"},
 {"OH", "Ohio"},
 {"OK", "Oklahoma"},
 {"OR", "Oregon"},
 {"PA", "Pennsylvania"},
 {"RI", "Rhode Island"},
 {"SC", "South Carolina"},
 {"SD", "South Dakota"},
 {"TN", "Tennessee"},
 {"TX", "Texas"},
 {"UT", "Utah"},
 {"VT", "Vermont"},
 {"VA", "Virginia"},
 {"WA", "Washington"}, 
 {"WV", "West Virginia"},
 {"WI", "Wisconsin"},
 {"WY", "Wyoming"},
 {NULL, NULL}
};


NLM_EXTERN CharPtr GetStateAbbreviation (CharPtr state)
{
  ReplacePairPtr fix;
  CharPtr        abbrev = NULL;

  fix = us_state_abbrev_fixes;
  while (fix->find != NULL && abbrev == NULL) {
    if (StringICmp (fix->replace, state) == 0) {
      abbrev = fix->find;
    } 
    fix++;
  }
  return abbrev;
}


static void FixUSStateAbbreviations (CharPtr PNTR pCountry)
{
  ReplacePairPtr fix;

  if (pCountry == NULL || StringHasNoText (*pCountry))
  {
    return;
  }

  fix = us_state_abbrev_fixes;
  while (fix->find != NULL) 
  {
    FindReplaceString (pCountry, fix->find, fix->replace, TRUE, TRUE);
    fix++;
  }
}


static CharPtr MoveStateAndAddComma (CharPtr cntry_str, CharPtr state_match, Int4 len_cntry)
{
  CharPtr newname = NULL, cp;
  Int4    len_state, len_qual, len_after, len_before;
  
  if (StringHasNoText (cntry_str) || StringHasNoText (state_match) || len_cntry < 1)
  {
    return cntry_str;
  }
  
  cp = StringISearch (cntry_str + len_cntry + 2, state_match);
  if (cp != NULL)
  {
    len_state = StringLen (state_match);
    len_qual = StringLen (cntry_str);
    
    if (cp == cntry_str + len_cntry + 2)
    {
      /* state is at beginning of string */
      len_after = len_qual - len_cntry - 2 - len_state;
      if (len_after == 0 || cntry_str [len_cntry + 2 + len_state] == ',')
      {
        /* already in correct format, nothing after state name */
        /* just copy in state name, in case we are correcting case */
        StringNCpy (cp, state_match, len_state);
        return cntry_str;
      }
      else
      {
        /* insert comma */
        newname = (CharPtr) MemNew (StringLen (cntry_str) + 3);
        StringNCpy (newname, cntry_str, len_cntry + 2 + len_state);
        newname [len_cntry + 2 + len_state] = 0;
        StringCat (newname, ",");
        StringCat (newname, cntry_str + len_cntry + 2 + len_state);
        cntry_str = MemFree (cntry_str);
        cntry_str = newname;
      }
    }
    else
    {
      newname = (CharPtr) MemNew (StringLen (cntry_str) + 3);
      StringNCpy (newname, cntry_str, len_cntry + 2);
      newname [len_cntry + 2] = 0;
      StringCat (newname, state_match);
      StringCat (newname, ", ");
      len_before = cp - cntry_str - 3 - len_cntry;
      StringNCpy (newname + len_cntry + 2 + len_state + 2,
                  cntry_str + len_cntry + 2,
                  len_before);
      newname [len_cntry + 2 + len_state + 2 + len_before] = 0;
      StringCat (newname, cp + len_state);
      cntry_str = MemFree (cntry_str);
      cntry_str = newname;
    }
  }
  return cntry_str;
}

typedef struct namedregion {
  CharPtr country;
  CharPtr state;
  CharPtr region;
} NamedRegionData, PNTR NamedRegionPtr;


static NamedRegionData named_regions[] = {
{ "USA", "Alaska", "Aleutian Islands" }
};

static Int4 num_named_regions = sizeof (named_regions) / sizeof (NamedRegionData);

static void TrimInternalSpacesAndLeadingPunct (CharPtr str)
{
  CharPtr src, dst;

  src = str;
  dst = str;

  while (*src != 0) {
    if (isspace (*src)) {
      if (dst > str && !isspace (*(dst - 1))) {
        *dst = ' ';
        dst++;
      }
    } else if (ispunct (*src)) {
      if (dst > str) {
        *dst = *src;
        dst++;
      }
    } else {
      *dst = *src;
      dst++;
    }
    src++;
  }
  if (dst > src && (isspace (*(dst - 1)))) {
    *(dst - 1) = 0;
  } else {
    *dst = 0;
  }
}

static void FixForNamedRegions (CharPtr PNTR pCountry)
{
  Int4 i, country_len, state_len, region_len;
  CharPtr region = NULL, country, state, new_str;

  if (pCountry == NULL || StringHasNoText (*pCountry)) return;

  for (i = 0; i < num_named_regions && region == NULL; i++) {
    region = StringSearch (*pCountry, named_regions[i].region);
    if (region != NULL) {
      country_len = StringLen (named_regions[i].country);
      state_len = StringLen (named_regions[i].state);
      country = StringSearch (*pCountry, named_regions[i].country);
      region_len = StringLen (named_regions[i].region);
      if (country != NULL) {
        MemSet (country, ' ', country_len);
        if (*(country + country_len) == ':') {
          *(country + country_len) = ' ';
        }
      }
      state = StringSearch (*pCountry, named_regions[i].state);
      if (state != NULL) {
        MemSet (state, ' ', state_len);
        if (*(state + state_len) == ',') {
          *(state + state_len) = ' ';
        }
      }
      MemSet (region, ' ', region_len);
      if (ispunct (*(region + region_len))) {
        *(region + region_len) = ' ';
      }
      TrimInternalSpacesAndLeadingPunct (*pCountry);
      new_str = (CharPtr) MemNew (sizeof (Char) * (country_len + state_len + region_len + 7 + StringLen (*pCountry)));
      sprintf (new_str, "%s: %s, %s", named_regions[i].country, named_regions[i].state, named_regions[i].region);
      if (!StringHasNoText (*pCountry)) {
        StringCat (new_str, ", ");
        StringCat (new_str, *pCountry);
      }
      *pCountry = MemFree (*pCountry);
      *pCountry = new_str;
    }
  }
}


static void FindCountryName (CharPtr PNTR pCountry, CharPtr PNTR country_list)
{
  CharPtr       best_match = NULL, state_match, state_country = NULL;
  CharPtr       cp, before, newname, after;
  Int4          len_cntry = 0, len_state = 0, len_qual, len_name;
  Boolean       state_multi = FALSE, country_multi = FALSE;

  if (pCountry == NULL || StringHasNoText (*pCountry))
  {
    return;
  }
  
  best_match = FindCountryMatch (*pCountry, country_list, &country_multi);
  if (country_multi) {
    *pCountry = MemFree (*pCountry);
    return;
  }
  state_match = FindStateMatch (*pCountry, &state_country, &len_state, &state_multi);

  if ((best_match == NULL && state_match == NULL) || (best_match == NULL && state_multi)) {
    *pCountry = MemFree (*pCountry);
    return;
  } else if (best_match != NULL && state_match != NULL && StringCmp (best_match, state_country) != 0) {
    state_match = NULL;
  }

  /* if match could be a country or a state, treat it as a country */
  if (StringCmp (best_match, state_match) == 0) {
    state_match = NULL;
  }

  if (IsBodyOfWater (best_match) && state_match != NULL)
  {
    /* prefer state to body of water */
    best_match = NULL;
  }

  /* if we have a country and a state, but the state is for a different country, drop the state */
  if (state_match != NULL && best_match != NULL && StringNCmp (state_country, best_match, len_cntry) != 0)
  {
    state_multi = FALSE;
    state_match = FindStateMatchForCountry (*pCountry, best_match, &state_multi);
    if (state_multi) {
      *pCountry = MemFree (*pCountry);
      return;
    }
  }

  if (best_match != NULL && StringCmp (best_match, "USA") == 0 && StringLen (*pCountry) > 3 && state_match == NULL) 
  {
    FixUSStateAbbreviations (pCountry);
    state_multi = FALSE;
    state_match = FindStateMatchForCountry (*pCountry, best_match, &state_multi);
    if (state_multi)
    {
      *pCountry = MemFree (*pCountry);
      return;
    }
    if (state_match != NULL) {
      FindReplaceString (pCountry, "USA:", "", TRUE, TRUE);
      FindReplaceString (pCountry, "USA", "", TRUE, TRUE);
      best_match = NULL;
      state_country = "USA";
    }
  }
  
  if (best_match == NULL && state_match == NULL) {
    *pCountry = MemFree (*pCountry);
    return;
  }
  else if (best_match == NULL && state_match != NULL)
  {
    FixCountryStringForStateName (pCountry, state_match, state_country);
  }
  else
  { 	
  	cp = StringISearch (*pCountry, best_match);
  	len_cntry = StringLen (best_match);
  	after = cp + len_cntry;
  	while (isspace (*after) || ispunct(*after)) 
  	{
  	  after++;
  	}
  	
    if (cp != NULL && !isalpha ((Int4)(cp [len_cntry])))
    {
      len_qual = StringLen (*pCountry);
    	if (cp == *pCountry)
    	{
     	  newname = (CharPtr) MemNew (len_cntry + StringLen (after) + 3);
     	  sprintf (newname, "%s: %s", best_match, after);
      }
      else
      {
        /* strip spaces and punctuation from before */
       	*(cp - 1) = 0;
       	before = cp - 2;
       	while (before >= *pCountry  
       	       && (isspace (*before) || ispunct (*before)))
       	{
       	  *before = 0;
       	  before--;
       	}
       	before = *pCountry;
       	while (isspace (*before) || ispunct(*before))
       	{
       	  before++;
       	}
       	
       	newname = (CharPtr) MemNew (len_cntry + StringLen (before) + StringLen (after) + 4);
       	sprintf (newname, "%s: %s%s%s", best_match, before,
       	         StringHasNoText (before) || StringHasNoText(after) ? "" : " ",
       	         after);
      }
      if (state_match != NULL)
      {
        newname = MoveStateAndAddComma (newname, state_match, len_cntry);          
      }
        
      /* remove trailing spaces and punctuation */
      len_name = StringLen (newname);        
      while (isspace ((Int4)(newname[len_name - 1])) 
             || newname [len_name - 1] == ','
             || newname [len_name - 1] == ':'
             || newname [len_name - 1] == ';')
      {
        newname [len_name - 1] = 0;
        len_name --;
      }                        	      	
      MemFree (*pCountry);
      *pCountry = newname;
    }
  }    
}


static void CountryColonToComma (CharPtr PNTR country_str)
{
  CharPtr cp, cp1, cp2, new_name;
  Int4    pre_len;
  
  if (country_str == NULL || *country_str == NULL) {
    return;
  }
  
  cp = StringChr (*country_str, ':');
  if (cp == NULL) return;
  cp = StringChr (cp + 1, ':');
  while (cp != NULL) {
    cp1 = cp;
    while (cp1 > *country_str && (isspace (*(cp1 - 1)) || *(cp1 - 1) == ',')) {
      cp1--;
    }
    pre_len = cp1 - *country_str;
    cp2 = cp + 1;
    while (isspace (*cp2) || *cp2 == ',') {
      cp2++;
    }
    new_name = (CharPtr) MemNew ((pre_len + StringLen (cp2) + 3) * sizeof (Char));
    StringNCpy (new_name, *country_str, pre_len);
    StringCat (new_name, ", ");
    StringCat (new_name, cp2);
    *country_str = MemFree (*country_str);
    *country_str = new_name;
    cp = StringChr ((*country_str) + pre_len, ':');
  } 
}

static void RemoveDoubleCommas (CharPtr PNTR country_str)
{
  CharPtr cp, cp1, cp2, new_name;
  Int4    pre_len;
  Boolean found_second_comma;
  
  if (country_str == NULL || *country_str == NULL) {
    return;
  }
  
  cp = StringChr (*country_str, ',');
  while (cp != NULL) {
    cp1 = cp;
    while (cp1 > *country_str && (isspace (*(cp1 - 1)) || *(cp1 - 1) == ',')) {
      cp1--;
    }
    pre_len = cp1 - *country_str;
    cp2 = cp + 1;
    found_second_comma = FALSE;
    while (isspace (*cp2) || *cp2 == ',') {
      if (*cp2 == ',') {
        found_second_comma = TRUE;
      }
      cp2++;
    }
    
    if (cp1 < cp || found_second_comma || cp2 > cp + 2) {
      new_name = (CharPtr) MemNew ((pre_len + StringLen (cp2) + 3) * sizeof (Char));
      StringNCpy (new_name, *country_str, pre_len);
      StringCat (new_name, ", ");
      StringCat (new_name, cp2);
      *country_str = MemFree (*country_str);
      *country_str = new_name;
      cp = StringChr ((*country_str) + pre_len + 1, ',');
    } else {
      cp = StringChr (cp2, ',');
    }
  }   
}


static Boolean ContainsMultipleCountryNames (CharPtr PNTR list, CharPtr search_str)
{
  CharPtr PNTR  ptr;
  Int4          len_match;
  CharPtr       cp;
  Boolean       found_one = FALSE;
  
  if (list == NULL || search_str == NULL) return FALSE;
  
  for (ptr = list; ptr != NULL && *ptr != NULL; ptr++)
  {
    cp = StringISearch (search_str, *ptr);
    len_match = StringLen (*ptr);
    while (cp != NULL) {
      /* if character after match is alpha, continue */
      if (isalpha ((Int4)(cp [len_match]))
          /* if character before match is alpha, continue */
          || (cp > search_str && isalpha ((Int4)(*(cp - 1))))
        /* if is shorter match for other item, continue */
        || IsSubstringOfStringInList (search_str, cp, *ptr, list)) {
        cp = StringSearch (cp + len_match, *ptr);
      } else if (found_one) {
        return TRUE;
      } else {
        found_one = TRUE;
        cp = StringSearch (cp + len_match, *ptr);
      }
    }
  }
  return FALSE;
}


static CharPtr NewFixCountry (CharPtr country, CharPtr PNTR country_list)
{
  CharPtr cp, next_sep, start_after;
  CharPtr valid_country = NULL, new_country = NULL, tmp;
  Char    ch;
  CharPtr separator_list = ",:";
  Boolean too_many_countries = FALSE;
  Int4    len_country, len_before, len_after, len_diff;
  ReplacePairPtr fix;
  Boolean fix_found;

  country = StringSave (country);
  cp = country;
  while (*cp != 0 && !too_many_countries) {
    next_sep = cp + StringCSpn (cp, separator_list);
    ch = *next_sep;
    *next_sep = 0;
    
    if (CountryIsValid (cp, NULL)) {
      if (valid_country == NULL) {
        valid_country = cp;
      } else {
        too_many_countries = TRUE;
      }
    } else {
      /* see if this is a fixable country */
      fix = country_name_fixes;
      fix_found = FALSE;
      while (fix->find != NULL && !fix_found) {
        if (StringCmp (fix->find, cp) == 0) {
          fix_found = TRUE;
          if (valid_country == NULL) {
            len_before = cp - country;
            if (ch == 0) {
              len_after = 0;
            } else {
              len_after = StringLen (next_sep + 1) + 1;
            }
            len_diff = StringLen (fix->replace) - StringLen (fix->find);
            len_country = StringLen (country) + len_diff + len_after + 1;       
            tmp = (CharPtr) MemNew (sizeof (Char) * len_country);
            if (len_before > 0) {
              StringNCpy (tmp, country, len_before);
            }
            StringCpy (tmp + len_before, fix->replace);
            if (len_after > 0) {
              StringCpy (tmp + len_before + StringLen (fix->replace) + 1, next_sep + 1);
            }
            cp = tmp + len_before;
            valid_country = cp;
            next_sep = tmp + (next_sep - country) + len_diff;
            country = MemFree (country);
            country = tmp;
          } else {
            too_many_countries = TRUE;
          }
        }
        fix++;
      }
    }

    *next_sep = ch;
    if (*next_sep == 0) {
      cp = next_sep;
    } else {
      cp = next_sep + 1;
      while (isspace (*cp)) {
        cp++;
      }
    }
  }
  if (valid_country != NULL && !too_many_countries) {
    too_many_countries = ContainsMultipleCountryNames (country_list, country);
  }

  if (valid_country != NULL && !too_many_countries) {
    len_country = StringCSpn (valid_country, separator_list);
    len_before = valid_country - country;

    while (len_before > 0 
           && (isspace (country [len_before - 1]) 
               || StringChr (separator_list, country [len_before - 1]) != NULL)) {
      len_before--;
    }
    start_after = valid_country + len_country;
    while (*start_after != 0 
           && (isspace (*start_after)
               || StringChr (separator_list, *start_after) != NULL)) {
      start_after++;
    }

    len_after = StringLen (start_after);

    new_country = MemNew (sizeof (Char) * (len_country + len_before + len_after + 5));
    
    StringNCpy (new_country, valid_country, len_country);
    if (len_before > 0 || len_after > 0) {
      StringCat (new_country, ": ");
      if (len_before > 0) {
        StringNCat (new_country, country, len_before);
        if (len_after > 0) {
          StringCat (new_country, ", ");
        }
      }
      if (len_after > 0) {
        StringCat (new_country, start_after);
      }
    }
  }
  country = MemFree (country);
  return new_country;
}


extern CharPtr GetCountryFix (CharPtr country, CharPtr PNTR country_list)
{
  CharPtr new_country;

  if (StringHasNoText (country)) return NULL;
#if 1
  new_country = NewFixCountry (country, country_list);
#else
  new_country = StringSave (country);
  FixCountryNames (&new_country);  	
  FindCountryName (&new_country, country_list);
  CountryColonToComma (&new_country);
  RemoveDoubleCommas (&new_country);
  FixForNamedRegions (&new_country);
#endif
  return new_country;
}


extern ValNodePtr ListFeaturesInLocation (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice)
{
  ValNodePtr        feat_list = NULL;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  Int4              loc_left, loc_right, tmp;

  if (bsp == NULL || slp == NULL) return NULL;

  loc_left = SeqLocStart (slp);
  loc_right = SeqLocStop (slp);
  if (loc_left > loc_right) {
    tmp = loc_left;
    loc_left = loc_right;
    loc_right = tmp;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, seqfeatChoice, featdefChoice, &fcontext);
       sfp != NULL && fcontext.left <= loc_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, seqfeatChoice, featdefChoice, &fcontext))
  {
    if (fcontext.right < loc_left) continue;
    if (SeqLocCompare (sfp->location, slp) == SLC_A_IN_B)
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


extern ValNodePtr ListFeaturesOverlappingLocation (BioseqPtr bsp, SeqLocPtr slp, Uint1 seqfeatChoice, Uint1 featdefChoice)
{
  ValNodePtr        feat_list = NULL;
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;
  Int4              loc_left, loc_right, tmp;
  Int4              cmp;

  if (bsp == NULL || slp == NULL) return NULL;

  loc_left = SeqLocStart (slp);
  loc_right = SeqLocStop (slp);
  if (loc_left > loc_right) {
    tmp = loc_left;
    loc_left = loc_right;
    loc_right = tmp;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, seqfeatChoice, featdefChoice, &fcontext);
       sfp != NULL && fcontext.left <= loc_right;
       sfp = SeqMgrGetNextFeature (bsp, sfp, seqfeatChoice, featdefChoice, &fcontext))
  {
    cmp = SeqLocCompare (sfp->location, slp);
    if (cmp != SLC_NO_MATCH)
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


static void CDSInSrcFeatCallback (BioseqPtr bsp, Pointer data)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        sfp;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_BIOSRC, &fcontext);
  while (sfp != NULL)
  {
    ValNodeLink ((ValNodePtr PNTR) data, ListFeaturesInLocation (bsp, sfp->location, 0, FEATDEF_CDS));
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_BIOSRC, &fcontext);
  }
}

extern ValNodePtr ListCodingRegionsContainedInSourceFeatures (SeqEntryPtr sep)
{
  ValNodePtr feat_list = NULL;

  VisitBioseqsInSep (sep, &feat_list, CDSInSrcFeatCallback);
  return feat_list;
}


