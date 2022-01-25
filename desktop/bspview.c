/*   bspview.c
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
* File Name:  bspview.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/30/95
*
* $Revision: 6.54 $
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

#include <bspview.h>
#include <document.h>
#include <picture.h>
#include <viewer.h>
#include <drawseq.h>
#include <objfdef.h>
#include <gather.h>
#include <subutil.h>
#include <asn2ff.h>
#include <tofasta.h>
#include <txalign.h>
#include <fstyle.h>
#include <picturep.h>
#include <drawingp.h>
#include <viewerp.h>
#include <objentr.h>
#include <accentr.h>
#include <mapgene.h>
#include <saledit.h>
#include <fea2seg.h>
#include <blast.h>
#include <blastpri.h>
#include <explore.h>

extern ForM smartBioseqViewForm;
ForM smartBioseqViewForm = NULL;


typedef struct bioseqviewform {
  FORM_MESSAGE_BLOCK

  BioseqPagePtr   bioseqNucPageList;
  BioseqPagePtr   bioseqProtPageList;
  BioseqPagePtr   currentBioseqPage;

  Int2            currentNucPage;
  Int2            currentProtPage;

  Handle          nucViewControl;
  Handle          protViewControl;
  Handle          targetControl;
  EnumFieldAssoc  PNTR targetAlist;
  Boolean         usePopupForTarget;
  Int4            numTargets;
  Int2            targetScratchSpace;
  GrouP           controls;
  GrpActnProc     updateControls;
  GrouP           retrieveAlignments;
  SeqViewUpdateFetchCounts  updateCounts;
  Boolean         hasaligns;

  EnumFieldAssoc  PNTR workingAlist;
  Int2            workingCount;
  Int2            workingTargets;

  BioseqViewData  bvd;

  Boolean         cleanupObjectPtr;
  WndActnProc     activateForm;

  ForM            toolForm;
} BioseqViewForm, PNTR BioseqViewFormPtr;

static void LookForGenomeTag (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     descr;
  BoolPtr        rsltptr;
  UserObjectPtr  uop;

  rsltptr = (BoolPtr) mydata;
  if (rsltptr == NULL) return;
  descr = NULL;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    descr = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    descr = bssp->descr;
  } else return;
  while (descr != NULL) {
    if (descr->choice == Seq_descr_user) {
      uop = (UserObjectPtr) descr->data.ptrvalue;
      if (uop != NULL && StringICmp (uop->_class, "Genomes") == 0) {
        *rsltptr = TRUE;
      }
    }
    descr = descr->next;
  }
}

extern Boolean LIBCALL IsAGenomeRecord (SeqEntryPtr sep)

{
  Boolean  rsult;

  rsult = FALSE;
  SeqEntryExplore (sep, &rsult, LookForGenomeTag);
  return rsult;
}

typedef struct updatesegstruc {
  BioseqSetPtr      parts;
  BioseqPtr         segseq;
  BioseqSetPtr      segset;
} UpdateSegStruc, PNTR UpdateSegStrucPtr;

static void FindSegSetComponentsCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)

{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  UpdateSegStrucPtr  ussp;

  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
    ussp = (UpdateSegStrucPtr) mydata;
    if (sep->choice == 1) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (ISA_na (bsp->mol) && bsp->repr == Seq_repr_seg) {
        ussp->segseq = bsp;
      }
    } else if (sep->choice == 2) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp->_class == 2) {
        ussp->segset = bssp;
      } else if (bssp->_class == 4) {
        ussp->parts = bssp;
      }
    }
  }
}

static Int4 UpdateSegList (SeqEntryPtr sep, Pointer mydata,
                           SeqEntryFunc mycallback,
                           Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return index;
  if (mycallback != NULL)
    (*mycallback) (sep, mydata, index, indent);
  index++;
  if (IS_Bioseq (sep)) return index;
  if (Bioseq_set_class (sep) == 4) return index;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  sep = bssp->seq_set;
  indent++;
  while (sep != NULL) {
    index = UpdateSegList (sep, mydata, mycallback, index, indent);
    sep = sep->next;
  }
  return index;
}

#define UpdateSegExplore(a,b,c) UpdateSegList(a, b, c, 0L, 0);

extern Boolean IsSegmentedBioseqWithoutParts (SeqEntryPtr sep)

{
  UpdateSegStruc  uss;

  if (sep == NULL) return FALSE;
  uss.segseq = NULL;
  uss.parts = NULL;
  uss.segset = NULL;
  UpdateSegExplore (sep, (Pointer) &uss, FindSegSetComponentsCallback);
  if (uss.segseq != NULL && uss.parts == NULL) return TRUE;
  return FALSE;
}

static void LookForDeltaBioseq (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  BoolPtr    rsltptr;

  rsltptr = (BoolPtr) mydata;
  if (rsltptr == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL && bsp->repr == Seq_repr_delta) {
      *rsltptr = TRUE;
    }
  }
}

extern Boolean IsADeltaBioseq (SeqEntryPtr sep)

{
  Boolean  rsult;

  rsult = FALSE;
  SeqEntryExplore (sep, &rsult, LookForDeltaBioseq);
  return rsult;
}

static SeqIdPtr get_db_sip (SeqFeatPtr sfp)
{
	GeneRefPtr grp;
	ValNodePtr db;
	DbtagPtr db_tag;
	ObjectIdPtr obj_id;
	CharPtr acc;

	if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return NULL;
	grp = sfp->data.value.ptrvalue;
	if (grp == NULL) return NULL;
	for (db=grp->db; db!=NULL; db=db->next)
	{
		db_tag = db->data.ptrvalue;
		if(StringICmp(db_tag->db, "GenBank") == 0)
		{
			obj_id = db_tag->tag;
			acc = obj_id->str;
			return gb_id_make(NULL, acc);
		}
	}

	return NULL;
}
	
static Boolean NamedAlignmentProc (GatherContextPtr gcp)

{
  AnnotDescrPtr  adp;
  Uint1          extra_type;
  ObjectIdPtr    oip;
  BoolPtr        rsult;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  SeqIdPtr       sip;
  UserObjectPtr  uop;

  if (gcp == NULL || gcp->thisitem == NULL) {
    return FALSE;
  }
  rsult = (BoolPtr) gcp->userdata;
  if (rsult == NULL) return FALSE;
  switch (gcp->thistype) {
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      if (gcp->parenttype == OBJ_SEQANNOT) {
        sap = (SeqAnnotPtr) gcp->parentitem;
        if (sap != NULL) {
          adp = sap->desc;
          while (adp != NULL) {
            if (adp->choice == Annot_descr_user) {
              uop = (UserObjectPtr) adp->data.ptrvalue;
              if (uop != NULL) {
                oip = uop->type;
                if (oip != NULL) {
                  if (StringICmp (oip->str, "Hist Seqalign") == 0) {
                    *rsult = TRUE;
                  }
                }
              }
            }
            adp = adp->next;
          }
        }
      }
      return TRUE;
    case OBJ_BIOSEQ_MAPFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      extra_type = ck_seqfeat_extra (sfp);
      if (extra_type & EXTRA_GENBANK)
      {
        sip = get_db_sip (sfp);
        if (sip != NULL) {
          *rsult = TRUE;
        }
        SeqIdFree (sip);
        return TRUE;
      }
      break;
    default :
      break;
  }
  return FALSE;
}

extern Boolean LIBCALL IsANamedAlignment (Uint2 entityID, Uint2 itemID, Uint2 itemtype)

{
  Boolean  rsult;

  rsult = FALSE;
  GatherItem (entityID, itemID, itemtype, (Pointer) (&rsult), NamedAlignmentProc);
  return rsult;
}

extern Boolean LIBCALL LaunchViewerNotEditor (BioseqViewPtr bvp, SeqEntryPtr sep,
                                              Uint2 entityID, Uint2 itemID, Uint2 itemtype)

{
  if (bvp == NULL) return FALSE;
  if (! bvp->launchEditors) return FALSE;
  if (IsAGenomeRecord (sep) ||
      IsSegmentedBioseqWithoutParts (sep) ||
      IsADeltaBioseq (sep) ||
      IsANamedAlignment (entityID, itemID, itemtype)) {
    if (itemtype == OBJ_BIOSEQ_SEG || itemtype == OBJ_BIOSEQ_DELTA ||
        itemtype == OBJ_SEQALIGN || itemtype == OBJ_SEQHIST_ALIGN ||
        itemtype == OBJ_BIOSEQ_MAPFEAT) {
      return TRUE;
    }
  }
  return FALSE;
}

static void AddOneBlastAlignment (BioseqPtr subject, BioseqPtr query)

{
  Uint1                align_type = 0;
  BioseqPtr            bsp;
  BioseqSetPtr         bssp;
  SeqAnnotPtr          curr;
  BLAST_OptionsBlkPtr  options = NULL;
  SeqAlignPtr          prev;
  SeqAlignPtr          salp;
  SeqAnnotPtr          sap;
  SeqAnnotPtr PNTR     sapp;
  BlastSearchBlkPtr    search;
  SeqEntryPtr          sep;

  if (subject == NULL || query == NULL) return;
  sap = NULL;
  salp = NULL;
  if (ISA_na (subject->mol)) {
    if (! ISA_na (query->mol)) return;
    align_type = 1;
    options = BLASTOptionNew ("blastn", TRUE);
    if (options != NULL) {
      options->gapped_calculation = TRUE;
      options->db_length = 100000000;
#ifdef WIN16
      options->wordsize = 10;
#else
      options->wordsize = 12;
#endif
    }
  } else if (ISA_aa (subject->mol)) {
    if (! ISA_aa (query->mol)) return;
    align_type = 2;
    options = BLASTOptionNew ("blastp", TRUE);
    if (options != NULL) {
      options->gapped_calculation = TRUE;
      options->db_length = 20000000;
      options->threshold_second = 12;
    }
  } else return;
  search = BLASTSetUpSearch (subject, options->program_name, 0, 0, NULL, options, NULL);

  salp = BlastSequencesOnTheFly (search, query);
  if (salp != NULL) {
    if (sap == NULL) {
      sap = SeqAnnotNew ();
      if (sap != NULL) {
        sap->type = 2;
      }
    }
    if (sap != NULL) {
      if (sap->data != NULL) {
        prev = sap->data;
        while (prev->next != NULL) {
          prev = prev->next;
        }
        prev->next = salp;
      } else {
        sap->data = (Pointer) salp;
      }
    }
  }

  BLASTOptionDelete (options);
  BlastSearchBlkDestruct (search);

  if (sap == NULL) return;

  AddAlignInfoToSeqAnnot (sap, align_type);
  /*
  ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
  */
  sapp = NULL;
  sep = SeqMgrGetSeqEntryForData (subject);
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sapp = &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sapp = &(bssp->annot);
  }
  if (sapp != NULL) {
    if (*sapp != NULL) {
      curr = *sapp;
      while (curr->next != NULL) {
        curr = curr->next;
      }
      curr->next = sap;
    } else {
      *sapp = sap;
    }
  }
}

static Boolean LaunchSequenceViewer (SeqIdPtr sip, BioseqPtr query)

{
  BioseqPtr        bsp;
  Uint2            entityID;
  Int2             handled;
  Uint2            itemID;
  SeqViewProcsPtr  svpp;

  if (sip == NULL) return FALSE;
  bsp = BioseqLockById (sip);
  if (bsp == NULL) return FALSE;
  entityID = BioseqFindEntity (sip, &itemID);
  if (entityID == 0) return FALSE;
  WatchCursor ();
  svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
  if (svpp != NULL) {
    svpp->forceSeparateViewer = TRUE;
    if (query != NULL && query->repr == Seq_repr_raw && query->length < 100000 &&
        bsp->repr == Seq_repr_raw && bsp->length < 100000) {
      if (svpp->alignWithChecked != NULL) {
        if (GetStatus (svpp->alignWithChecked)) {
          AddOneBlastAlignment (bsp, query);
        }
      } else if (svpp->alignDefault) {
        AddOneBlastAlignment (bsp, query);
      }
    }
  }
  handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, itemID,
                              OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
  ArrowCursor ();
  if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
    Message (MSG_FATAL, "Unable to launch viewer.");
  } else {
    ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
  }
  BioseqUnlockById (sip);

  return TRUE;
}

static Boolean LaunchPrim (GatherContextPtr gcp)

{
  SeqAlignPtr   align;
  BioseqPtr     bsp;
  DenseDiagPtr  ddp;
  DeltaSeqPtr   delsp;
  DenseSegPtr   dsp;
  Uint1         extra_type;
  SeqLocPtr     seg_loc;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  StdSegPtr     ssp;
/* SeqLocPtr     tloc; */
  ValNode       vn;

  if (gcp == NULL || gcp->thisitem == NULL) {
    Beep ();
    return FALSE;
  }
  bsp = (BioseqPtr) gcp->userdata;
  if (bsp == NULL) return FALSE;
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = (Pointer) bsp->id;
  slp = &vn;
  switch (gcp->thistype) {
    case OBJ_BIOSEQ_SEG :
      seg_loc = (SeqLocPtr) gcp->thisitem;
      sip = SeqLocId (seg_loc);
      if (! LaunchSequenceViewer (sip, bsp)) {
        Beep ();
        return FALSE;
      }
      return TRUE;
    case OBJ_BIOSEQ_DELTA :
      delsp = (DeltaSeqPtr) gcp->thisitem;
      if (delsp != NULL && delsp->choice == 1) {
        seg_loc = (SeqLocPtr) delsp->data.ptrvalue;
        sip = SeqLocId (seg_loc);
        if (! LaunchSequenceViewer (sip, bsp)) {
          Beep ();
          return FALSE;
        }
      }
      return TRUE;
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      align = (SeqAlignPtr) gcp->thisitem;
      sip = NULL;
      if (align->segtype == 1) {
        ddp = (DenseDiagPtr) align->segs;
        if (ddp != NULL) {
          for (sip = ddp->id; sip != NULL; sip = sip->next) {
            if (! SeqIdForSameBioseq (sip, SeqLocId (slp)))
              break;
          }
        }
      } else if (align->segtype == 2) {
        dsp = (DenseSegPtr) align->segs;
        if (dsp != NULL) {
          if (dsp->ids != NULL) {
            sip = dsp->ids->next;
          }
          /*
          for (sip = dsp->ids; sip != NULL; sip = sip->next) {
            if (! SeqIdForSameBioseq (sip, SeqLocId (slp)))
              break;
          }
          */
        }
      } else if (align->segtype == 3) {
        ssp = (StdSegPtr) align->segs;
        if (ssp != NULL && ssp->loc != NULL) {
          if (ssp->loc->next != NULL) {
            sip = SeqLocId (ssp->loc->next);
          }
          /*
          for (tloc = ssp->loc; tloc != NULL; tloc = tloc->next) {
            if (! SeqIdForSameBioseq (SeqLocId (tloc), SeqLocId (slp))) {
              sip = SeqLocId (tloc);
              break;
            }
          }
          */
        }
      }
      if (sip != NULL) {
        if (! LaunchSequenceViewer (sip, bsp)) {
          Beep ();
          return FALSE;
        }
      }
      return TRUE;
    case OBJ_BIOSEQ_MAPFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      extra_type = ck_seqfeat_extra (sfp);
      if (extra_type & EXTRA_GENBANK)
      {
        sip = get_db_sip (sfp);
        if (! LaunchSequenceViewer (sip, bsp)) {
          SeqIdFree (sip);
          Beep ();
          return FALSE;
        }
        SeqIdFree (sip);
        return TRUE;
      }
      if (sfp->product != NULL) {
        sip = SeqLocId (sfp->product);
        if (sip != NULL) {
          if (! LaunchSequenceViewer (sip, bsp)) {
            Beep ();
            return FALSE;
          }
        }
      }
      break;
    default :
      break;
  }
  return FALSE;
}

void LIBCALL LaunchNewBioseqViewer (BioseqPtr bsp, Uint2 entityID, Uint2 itemID, Uint2 itemtype)

{
  GatherItem (entityID, itemID, itemtype, (Pointer) bsp, LaunchPrim);
}

extern void LIBCALL AddBioseqPageToList (BioseqPagePtr PNTR head, BioseqPagePtr bpp)

{
  BioseqPagePtr  newbpp;
  BioseqPagePtr  tmpbpp;

  if (head == NULL || bpp == NULL) return;
  newbpp = MemNew (sizeof (BioseqPageData));
  if (newbpp == NULL) return;
  if (*head != NULL) {
    tmpbpp = *head;
    while (tmpbpp->next != NULL) {
      tmpbpp = tmpbpp->next;
    }
    tmpbpp->next = newbpp;
  } else {
    *head = newbpp;
  }
  newbpp->label = StringSaveNoNull (bpp->label);
  newbpp->nucOK = bpp->nucOK;
  newbpp->protOK = bpp->protOK;
  newbpp->genomeOK = bpp->genomeOK;
  newbpp->needAlignment = bpp->needAlignment;
  newbpp->maxLength = bpp->maxLength;
  newbpp->populate = bpp->populate;
  newbpp->show = bpp->show;
  newbpp->highlight = bpp->highlight;
  newbpp->toClipboard = bpp->toClipboard;
  newbpp->print = bpp->print;
  newbpp->export = bpp->export;
  newbpp->togif = bpp->togif;
  newbpp->resize = bpp->resize;
  newbpp->next = NULL;
}

extern BioseqPagePtr LIBCALL BioseqPageListFree (BioseqPagePtr bpp)

{
  BioseqPagePtr  next;

  while (bpp != NULL) {
    next = bpp->next;
    bpp->label = MemFree (bpp->label);
    MemFree (bpp);
    bpp = next;
  }
  return NULL;
}

static void LookInSeqIdList (SeqIdPtr sip, ValNodePtr PNTR vnpp, Uint1 align_type, Boolean useUids)

{
  Char  str [64];
  Int4  uid;

  while (sip != NULL) {
    switch (sip->choice) {
      case SEQID_NOT_SET :
      case SEQID_LOCAL :
        break;
      case SEQID_GI :
        if (useUids) {
          uid = (Int4) sip->data.intvalue;
          ValNodeAddInt (vnpp, align_type, uid);
        } else {
          SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
          ValNodeCopyStr (vnpp, align_type, str);
        }
        break;
      default :
        if (useUids) {
          uid = GetGIForSeqId (sip);
          if (uid > 0) {
            ValNodeAddInt (vnpp, align_type, uid);
          }
        } else {
          SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
          ValNodeCopyStr (vnpp, align_type, str);
        }
        break;
    }
    sip = sip->next;
  }
}

static void LookInSeqLocList (SeqLocPtr slp, ValNodePtr PNTR vnpp, Uint1 align_type, Boolean useUids)

{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        LookInSeqIdList (sip, vnpp, align_type, useUids);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          LookInSeqIdList (sip, vnpp, align_type, useUids);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          LookInSeqIdList (sip, vnpp, align_type, useUids);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          LookInSeqIdList (sip, vnpp, align_type, useUids);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          LookInSeqIdList (loc, vnpp, align_type, useUids);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            LookInSeqIdList (sip, vnpp, align_type, useUids);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            LookInSeqIdList (sip, vnpp, align_type, useUids);
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }
}

static void LIBCALL GetSeqIdsForOneSeqAnnot (SeqAnnotPtr sap, ValNodePtr PNTR vnpp, Uint1 align_type, Boolean useUids)

{
  DenseDiagPtr   ddp;
  DenseSegPtr    dsp;
  SeqAlignPtr    sal;
  StdSegPtr      ssp;

  if (sap == NULL || vnpp == NULL) return;
  if (sap->type == 2) {
    sal = (SeqAlignPtr) sap->data;
    while (sal != NULL) {
      if (sal->segtype == 1) {
        ddp = (DenseDiagPtr) sal->segs;
        if (ddp != NULL) {
          LookInSeqIdList (ddp->id, vnpp, align_type, useUids);
        }
      } else if (sal->segtype == 2) {
        dsp = (DenseSegPtr) sal->segs;
        if (dsp != NULL) {
          LookInSeqIdList (dsp->ids, vnpp, align_type, useUids);
        }
      } else if (sal->segtype == 3) {
        ssp = (StdSegPtr) sal->segs;
        if (ssp != NULL) {
           LookInSeqLocList (ssp->loc, vnpp, align_type, useUids);
        }
      }
      sal = sal->next;
    }
  }
}

extern void LIBCALL GetUidsForOneSeqAnnot (SeqAnnotPtr sap, ValNodePtr PNTR vnpp, Uint1 align_type)

{
  GetSeqIdsForOneSeqAnnot (sap, vnpp, align_type, TRUE);
}

static void GetAlignmentsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent, Boolean useUids)

{
  Uint1          align_type;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ObjectIdPtr    oip;
  SeqAnnotPtr    sap;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;
  ValNodePtr     vnp;

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
    if (sap->type == 2) {
      align_type = 0;
      for (vnp = sap->desc; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == Annot_descr_user) {
          uop = (UserObjectPtr) vnp->data.ptrvalue;
          if (uop != NULL) {
            oip = uop->type;
            if (oip != NULL && StringICmp (oip->str, "Blast Type") == 0) {
              ufp = uop->data;
              if (ufp != NULL && ufp->choice == 2) {
                align_type = ufp->data.intvalue;
              }
            }
          }
        }
      }
      GetSeqIdsForOneSeqAnnot (sap, (ValNodePtr PNTR) mydata, align_type, useUids);
    }
    sap = sap->next;
  }
}

static void GetUidsAlignmentsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  GetAlignmentsCallback (sep, mydata, index, indent, TRUE);
}

static void GetStrIdsAlignmentsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  GetAlignmentsCallback (sep, mydata, index, indent, FALSE);
}

extern int LIBCALLBACK SortByVnpDataIntvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->data.intvalue > vnp2->data.intvalue) {
        return 1;
      } else if (vnp1->data.intvalue < vnp2->data.intvalue) {
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

static int LIBCALLBACK SortByName (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (CharPtr) vnp1->data.ptrvalue;
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
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

extern ValNodePtr LIBCALL GetUidsForSeqEntryAligns (SeqEntryPtr sep)

{
  Uint1            choice;
  ValNodePtr       head;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  Int4             uid;
  ValNodePtr       vnp;

  head = NULL;
  if (sep == NULL) return NULL;
  SeqEntryExplore (sep, (Pointer) (&head), GetUidsAlignmentsCallback);
  if (head == NULL) return NULL;
  head = SortValNode (head, SortByVnpDataIntvalue);
  uid = 0;
  choice = 0;
  prev = (ValNodePtr PNTR) &head;
  vnp = head;
  while (vnp != NULL) {
    next = vnp->next;
    if (vnp->data.intvalue == uid && vnp->choice == choice) {
      *(prev) = vnp->next;
      vnp->next = NULL;
      ValNodeFree (vnp);
    } else {
      uid = vnp->data.intvalue;
      choice = vnp->choice;
      prev = (ValNodePtr PNTR) &(vnp->next);
    }
    vnp = next;
  }
  return head;
}

extern ValNodePtr LIBCALL GetIdStringsForSeqEntryAligns (SeqEntryPtr sep)

{
  Uint1            choice;
  ValNodePtr       head;
  CharPtr          last;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  head = NULL;
  if (sep == NULL) return NULL;
  SeqEntryExplore (sep, (Pointer) (&head), GetStrIdsAlignmentsCallback);
  if (head == NULL) return NULL;
  head = SortValNode (head, SortByName);
  last = NULL;
  choice = 0;
  prev = (ValNodePtr PNTR) &head;
  vnp = head;
  while (vnp != NULL) {
    next = vnp->next;
    if (StringICmp ((CharPtr) vnp->data.ptrvalue, last) == 0 && vnp->choice == choice) {
      *(prev) = vnp->next;
      vnp->next = NULL;
      ValNodeFree (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      choice = vnp->choice;
      prev = (ValNodePtr PNTR) &(vnp->next);
    }
    vnp = next;
  }
  return head;
}

static void ResizeViewForm (WindoW w)

{
  BioseqViewFormPtr  bfp;
  Int2               bottom;
  BioseqPagePtr      bpp;
  RecT               c;
  Int2               diff = 0;
  Int2               gap = 0;
  Int2               height;
  RecT               r;
  RecT               s;
  RecT               t;
  Int2               width;

  bfp = (BioseqViewFormPtr) GetObjectExtra (w);
  if (bfp == NULL) return;
  WatchCursor ();
  ObjectRect (w, &r);
  width = r.right - r.left;
  height = r.bottom - r.top;
  bottom = height - 10;
  SafeHide (bfp->controls);
  SafeHide (bfp->retrieveAlignments);

  if (bfp->controls != NULL) {
    GetPosition (bfp->controls, &c);
    LoadRect (&t, c.left, c.top, c.right, c.bottom);
    diff = t.bottom - t.top;
    gap = 10;
    t.bottom = height - 10;
    t.top = t.bottom - diff;
    t.right = width - 10;
    SetPosition (bfp->controls, &t);
    AdjustPrnt (bfp->controls, &t, FALSE);
    bottom = t.top - gap;
  } else if (bfp->retrieveAlignments != NULL) {
    GetPosition (bfp->retrieveAlignments, &c);
    LoadRect (&t, c.left, c.top, c.right, c.bottom);
    diff = t.bottom - t.top;
    gap = 10;
    t.bottom = height - 10;
    t.top = t.bottom - diff;
    SetPosition (bfp->retrieveAlignments, &t);
    AdjustPrnt (bfp->retrieveAlignments, &t, FALSE);
    if (bfp->hasaligns) {
      bottom = t.top - gap;
    }
  }

  if (bfp->bvd.vwr != NULL) {
    GetPosition (bfp->bvd.vwr, &s);
    s.right = width - 10;
    s.bottom = bottom;
    SetPosition (bfp->bvd.vwr, &s);
    AdjustPrnt (bfp->bvd.vwr, &s, FALSE);
  }
  if (bfp->bvd.doc != NULL) {
    GetPosition (bfp->bvd.doc, &s);
    s.right = width - 10;
    s.bottom = bottom;
    SetPosition (bfp->bvd.doc, &s);
    AdjustPrnt (bfp->bvd.doc, &s, FALSE);
  }
  if (bfp->bvd.text != NULL) {
    GetPosition (bfp->bvd.text, &s);
    s.right = width - 10;
    s.bottom = bottom;
    SetPosition (bfp->bvd.text, &s);
    AdjustPrnt (bfp->bvd.text, &s, FALSE);
  }
  if (bfp->bvd.pnl != NULL) {
    GetPosition (bfp->bvd.pnl, &s);
    s.right = width - 10;
    s.bottom = bottom;
    SetPosition (bfp->bvd.pnl, &s);
    AdjustPrnt (bfp->bvd.pnl, &s, FALSE);
  }
  /*
  if (bfp->bvd.vwr != NULL) {
    if (Visible (bfp->bvd.vwr) && AllParentsVisible (bfp->bvd.vwr)) {
      ViewerWasResized (bfp->bvd.vwr);
    }
  }
  if (bfp->bvd.doc != NULL) {
    if (Visible (bfp->bvd.doc) && AllParentsVisible (bfp->bvd.doc)) {
      UpdateDocument (bfp->bvd.doc, 0, 0);
    }
  }
  */
  bpp = bfp->currentBioseqPage;
  if (bpp != NULL) {
    if (bpp->resize != NULL) {
      bpp->resize (&(bfp->bvd));
    }
  }
  if (bfp->controls != NULL) {
    SafeShow (bfp->controls);
  } else if (bfp->retrieveAlignments != NULL && bfp->hasaligns) {
    SafeShow (bfp->retrieveAlignments);
  }
  ArrowCursor ();
  Update ();
}

static void PopTargetAlistProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqViewFormPtr  bfp;
  BioseqPtr          bsp;
  CharPtr            ptr;
  SeqIdPtr           sip;
  Char               str [128];

  bfp = (BioseqViewFormPtr) mydata;
  if (bfp != NULL && sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sip = SeqIdFindWorst (bsp->id);
    SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
    ptr = StringChr (str, '|');
    if (ptr == NULL) {
      ptr = str;
    } else {
      ptr++;
    }
    bfp->workingAlist [bfp->workingCount].name = StringSave (ptr);
    (bfp->workingCount)++;
    if (bsp == bfp->bvd.bsp) {
      bfp->targetScratchSpace = index + 1;
    }
  }
}

static EnumFieldAssocPtr MakeTargetAlist (BioseqViewFormPtr bfp, SeqEntryPtr sep)

{
  EnumFieldAssocPtr  alist;
  Int2               num;

  if (bfp == NULL || sep == NULL) return NULL;
  bfp->workingAlist = NULL;
  bfp->workingCount = 0;
  num = SequinEntryCount (sep);
  bfp->workingAlist = MemNew (sizeof (EnumFieldAssoc) * (num + 3));
  if (bfp->workingAlist == NULL) return NULL;
  bfp->workingAlist [bfp->workingCount].name = StringSave ("ALL SEQUENCES");
  (bfp->workingCount)++;
  bfp->workingTargets = SequinEntryExplore (sep, (Pointer) bfp, PopTargetAlistProc);
  alist = bfp->workingAlist;
  bfp->workingAlist = NULL;
  bfp->workingCount = 0;
  return alist;
}

static Int2 PopulateTarget (BioseqViewFormPtr bfp)

{
  EnumFieldAssocPtr  ap;
  Int2               count;
  Uint2              entityID;
  SeqEntryPtr        sep;
  Int2               val;

  val = 0;
  if (bfp != NULL && bfp->bvd.bsp != NULL) {
    bfp->targetAlist = FreeEnumFieldAlist (bfp->targetAlist);
    entityID = ObjMgrGetEntityIDForPointer (bfp->bvd.bsp);
    sep = GetTopSeqEntryForEntityID (entityID);
    if (sep != NULL) {
      bfp->targetScratchSpace = 0;
      bfp->workingTargets = 0;
      bfp->targetAlist = MakeTargetAlist (bfp, sep);
      for (ap = bfp->targetAlist, count = 0; ap != NULL && ap->name != NULL; ap++, count++) {
        if (bfp->usePopupForTarget) {
          if (count < 32) {
            PopupItem (bfp->targetControl, ap->name);
          }
        } else {
          ListItem (bfp->targetControl, ap->name);
        }
      }
      bfp->numTargets = bfp->workingTargets;
      val = bfp->targetScratchSpace;
    }
  }
  return val;
}

static Int4 GetUidFromBsp (BioseqPtr bsp)

{
  SeqIdPtr  sip;

  if (bsp == NULL) return 0;
  sip = bsp->id;
  while (sip != NULL) {
    if (sip->choice == SEQID_GI) {
      return (Int4) sip->data.intvalue;
    }
    sip = sip->next;
  }
  return 0;
}

static void CheckForCookedBioseqs (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BoolPtr    bp;
  BioseqPtr  bsp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bp = (BoolPtr) mydata;
  if (bp == NULL) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (bsp->repr != Seq_repr_raw && bsp->repr != Seq_repr_seg) {
    *bp = FALSE;
  }
}

static void BioseqPtrToBioseqForm (ForM f, Pointer data)

{
  EnumFieldAssocPtr  alist;
  Boolean            allRawOrSeg = TRUE;
  EnumFieldAssocPtr  ap1, ap2;
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;
  BioseqPtr          bsp;
  Int4               count;
  Uint2              entityID;
  SeqEntryPtr        sep;
  Int2               val;

  bfp = (BioseqViewFormPtr) GetObjectExtra (f);
  bsp = (BioseqPtr) data;
  if (bfp != NULL && bsp != NULL) {
    /*
    if (bfp->bvd.bsp != bsp) {
      BioseqLock (bsp);
      BioseqUnlock (bfp->bvd.bsp);
    }
    */
    bfp->bvd.bsp = bsp;
    bfp->docuid = GetUidFromBsp (bsp);
    if (ISA_na (bsp->mol)) {
      bfp->doctype = TYP_NT;
    } else if (ISA_aa (bsp->mol)) {
      bfp->doctype = TYP_AA;
    }
    entityID = ObjMgrGetEntityIDForPointer (bsp);
    sep = GetTopSeqEntryForEntityID (entityID);
    count = SequinEntryCount (sep);
    if (bfp->numTargets == count && sep != NULL) {
      alist = MakeTargetAlist (bfp, sep);
      if (alist != NULL && bfp->targetAlist) {
        for (ap1 = alist, ap2 = bfp->targetAlist;
             ap1->name != NULL && ap2->name != NULL;
             ap1++, ap2++) {
          if (StringICmp (ap1->name, ap2->name) != 0) {
            count = bfp->numTargets + 1; /* seqIDs have changed, so force repopulation */
          }
        }
      }
      alist = FreeEnumFieldAlist (alist);
    }
    if (bfp->numTargets != count) {
      val = GetValue (bfp->targetControl);
      Hide (bfp->targetControl);
      Update ();
      Reset (bfp->targetControl);
      PopulateTarget (bfp);
      SetValue (bfp->targetControl, val);
      Show (bfp->targetControl);
    }
    Update ();
    bpp = bfp->currentBioseqPage;
    if (bpp != NULL) {
      if (bpp->populate != NULL) {
        /* oldErrSev = ErrSetMessageLevel (SEV_FATAL); */
        if (bfp->bvd.hasTargetControl) { /* just sequin, for now */
          sep = GetTopSeqEntryForEntityID (entityID);
          SeqEntryExplore (sep, (Pointer) (&allRawOrSeg), CheckForCookedBioseqs);
          if (allRawOrSeg) {
            SeqMgrIndexFeatures (entityID, NULL, OM_LABEL_CONTENT);
          }
        }
        BioseqLock (bsp);
        bpp->populate (&(bfp->bvd));
        BioseqUnlock (bsp);
        Update ();
        /* ErrSetMessageLevel (oldErrSev); */
      }
      if (bfp->retrieveAlignments != NULL && bfp->updateCounts != NULL) {
        entityID = ObjMgrGetEntityIDForPointer (bsp);
        sep = GetTopSeqEntryForEntityID (entityID);
        bfp->hasaligns = bfp->updateCounts (bfp->retrieveAlignments, sep);
        if (Visible (bfp->retrieveAlignments)) {
          if (! bfp->hasaligns) {
            ResizeViewForm ((WindoW) f);
          }
        } else {
          if (bfp->hasaligns) {
            ResizeViewForm ((WindoW) f);
          }
        }
      }
      if (bpp->show != NULL) {
        bpp->show (&(bfp->bvd), TRUE);
      }
    }
  }
}

static void LIBCALL CopyBioseqViewFormToClipboard (Pointer formDataPtr)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;

  bfp = (BioseqViewFormPtr) formDataPtr;
  if (bfp == NULL || bfp->bvd.bsp == NULL) return;
  bpp = bfp->currentBioseqPage;
  if (bpp == NULL) return;
  if (bpp->toClipboard != NULL) {
    bpp->toClipboard (&(bfp->bvd));
  }
}

static void LIBCALL ExportBioseqViewFormToFile (Pointer formDataPtr, CharPtr filename)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;
  Char               dfault [32];

  bfp = (BioseqViewFormPtr) formDataPtr;
  if (bfp == NULL || bfp->bvd.bsp == NULL) return;
  bpp = bfp->currentBioseqPage;
  if (bpp == NULL) return;
  if (bpp->export != NULL) {
    GetTitle (bfp->form, dfault, sizeof (dfault));
    bpp->export (&(bfp->bvd), NULL, dfault);
  }
}

static Boolean ShortDefFastaFileFunc (BioseqPtr bsp, Int2 key, CharPtr buf, Uint4 buflen, Pointer data)

{
	Int2  code;
	CharPtr  buffer;
	Uint2  entityID;
	FILE * fp;
	size_t  len;
	Char  org [200];
	CharPtr  ptr;
	SeqEntryPtr  sep;
	Char  tmp [16];

	fp = (FILE *)data;

	switch (key)
	{
		case FASTA_ID:
			fprintf(fp, ">%s ", buf);
			break;
		case FASTA_DEFLINE:
			entityID = ObjMgrGetEntityIDForPointer (bsp);
			sep = GetBestTopParentForData (entityID, bsp);
			code = SeqEntryToGeneticCode (sep, NULL, org, sizeof (org) - 21);
			if (! StringHasNoText (org)) {
				StringCat (org, "]");
				if (code > 0) {
			 		sprintf (tmp, " [code=%d]", (int) code);
			 		StringCat (org, tmp);
				}
			}
			len = StringLen (buf) + StringLen (org) + 20;
			buffer = MemNew (len);
			ptr = StringChr (buf, ' ');
			if (ptr != NULL) {
			  *ptr = '\0';
			  ptr++;
			}
			StringCpy (buffer, buf);
			if (org [0] != '\0') {
			  StringCat (buffer, " [org=");
			  StringCat (buffer, org);
			}
			if (! StringHasNoText (ptr)) {
			  StringCat (buffer, " ");
			  StringCat (buffer, ptr);
			}
			if (StringLen (buffer) > 253) {
			  buffer [251] = '.';
			  buffer [252] = '.';
			  buffer [253] = '.';
			  buffer [254] = '\0';
			}
			fprintf(fp, ">%s\n", buffer);
			MemFree (buffer);
			break;
		case FASTA_SEQLINE:
			fprintf(fp, "%s\n", buf);
			break;
		case FASTA_EOS:   /* end of sequence */
			break;
		default:
			break;
	}
	return TRUE;
}

extern Boolean SeqnSeqEntrysToFasta (SeqEntryPtr sep, FILE *fp, Boolean is_na, Uint1 group_segs)

{
    FastaDat tfa;
    MyFsa mfa;
    Char buf[255];

    if ((sep == NULL) || (fp == NULL))
        return FALSE;

    MemSet ((Pointer) (&mfa), 0, sizeof (MyFsa));
    mfa.buf = buf;
    mfa.buflen = 254;
    mfa.seqlen = 70;
    mfa.mydata = (Pointer)fp;
    mfa.myfunc = ShortDefFastaFileFunc;
    mfa.bad_asn1 = FALSE;
    mfa.order = 0;
    mfa.accession = NULL;
    mfa.organism = NULL;
    mfa.do_virtual = FALSE;
    mfa.tech = 0;
    mfa.no_sequence = FALSE;
    mfa.formatdb	= FALSE;

    tfa.mfp = &mfa;
    tfa.is_na = is_na;

    if (is_na)
        mfa.code = Seq_code_iupacna;
    else
        mfa.code = Seq_code_ncbieaa;

    if (group_segs == 3)  /* do 2 things */
    {
        mfa.do_virtual = TRUE;
        group_segs = 1;
    }
    tfa.group_segs = group_segs;
    tfa.last_indent = -1;
    tfa.parts = -1;
    tfa.seg = -1;
    tfa.got_one = FALSE;
    SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryFasta);
    return tfa.got_one;
}

extern Boolean BioseqViewCanSaveFasta (ForM f, Boolean nucs, Boolean prots, Boolean onlyTarget)

{
  BioseqViewFormPtr  bfp;
  BioseqPtr          bsp;
  BioseqViewPtr      bvp;
  Uint2              entityID;
  SeqEntryPtr        sep;

  bfp = (BioseqViewFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    bvp = (&(bfp->bvd));
    if (bvp == NULL) return FALSE;
    bsp = bvp->bsp;
    if (bsp == NULL) return FALSE;
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (! onlyTarget) {
      entityID = ObjMgrGetEntityIDForChoice (sep);
      sep = GetTopSeqEntryForEntityID (entityID);
    }
    if (sep == NULL) return FALSE;
    if (nucs && SeqEntryHasNucs (sep)) return TRUE;
    if (prots && SeqEntryHasProts (sep)) return TRUE;
  }
  return FALSE;
}

extern Boolean ExportBioseqViewFasta (ForM f, CharPtr filename, Boolean nucs, Boolean prots, Boolean onlyTarget)

{
  BioseqViewFormPtr  bfp;
  BioseqPtr          bsp;
  BioseqViewPtr      bvp;
  Uint2              entityID;
  FILE               *fp;
  Uint1              group_segs;
  Char               path [PATH_MAX];
  SeqEntryPtr        sep;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  bfp = (BioseqViewFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    bvp = (&(bfp->bvd));
    if (bvp == NULL) return FALSE;
    bsp = bvp->bsp;
    if (bsp == NULL) return FALSE;
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (! onlyTarget) {
      entityID = ObjMgrGetEntityIDForChoice (sep);
      sep = GetTopSeqEntryForEntityID (entityID);
    }
    if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
#ifdef WIN_MAC
      fp = FileOpen (path, "r");
      if (fp != NULL) {
        FileClose (fp);
      } else {
        FileCreate (path, "TEXT", "ttxt");
      }
#endif
      fp = FileOpen (path, "w");
      if (fp != NULL) {
        WatchCursor ();
        Update ();
        group_segs = 0;
        if (bsp->repr == Seq_repr_seg) {
          group_segs = 2;
        }
        if (nucs) {
          SeqnSeqEntrysToFasta (sep, fp, TRUE, group_segs);
        }
        if (prots) {
          SeqnSeqEntrysToFasta (sep, fp, FALSE, 0);
        }
        FileClose (fp);
        ArrowCursor ();
        Update ();
        return TRUE;
      }
    }
  }
  return FALSE;
}

static void LIBCALL PrintBioseqViewForm (Pointer formDataPtr)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;

  bfp = (BioseqViewFormPtr) formDataPtr;
  if (bfp == NULL || bfp->bvd.bsp == NULL) return;
  bpp = bfp->currentBioseqPage;
  if (bpp == NULL) return;
  if (bpp->print != NULL) {
    bpp->print (&(bfp->bvd));
  }
}

extern void LIBCALL NewSaveBioseqViewFormGifItemTable (Pointer formDataPtr, CharPtr filename)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;

  bfp = (BioseqViewFormPtr) formDataPtr;
  if (bfp == NULL || bfp->bvd.bsp == NULL) return;
  bpp = bfp->currentBioseqPage;
  if (bpp == NULL) return;
  if (bpp->togif != NULL) {
    bpp->togif (&(bfp->bvd), filename, NULL);
  }
}

static void SetBioseqImportExportItems (BioseqViewFormPtr bfp)

{
  BioseqPagePtr  bpp;
  IteM           exportItm;
  Char           str [64];
  CharPtr        tmp;

  if (bfp == NULL || bfp->bvd.bsp == NULL) return;
  bpp = bfp->currentBioseqPage;
  if (bpp == NULL) return;
  exportItm = FindFormMenuItem ((BaseFormPtr) bfp, VIB_MSG_EXPORT);
  if (bpp->export != NULL) {
    tmp = StringMove (str, "Export ");
    StringNCpy_0 (tmp, bpp->label, sizeof (str) - 12);
    StringCat (tmp, "...");
    SafeSetTitle (exportItm, str);
    SafeEnable (exportItm);
  } else {
    SafeSetTitle (exportItm, "Export...");
    SafeDisable (exportItm);
  }
}

static void SetCurrentPagePointers (BioseqViewFormPtr bfp)

{
  BioseqPagePtr  bpp;
  Int2           page;

  if (bfp == NULL || bfp->bvd.bsp == NULL) return;
  bpp = NULL;
  page = 0;
  if (ISA_na (bfp->bvd.bsp->mol)) {
    bpp = bfp->bioseqNucPageList;
    page = bfp->currentNucPage;
  } else if (ISA_aa (bfp->bvd.bsp->mol)) {
    bpp = bfp->bioseqProtPageList;
    page = bfp->currentProtPage;
  }
  while (bpp != NULL && page > 0) {
    bpp = bpp->next;
    page--;
  }
  bfp->currentBioseqPage = bpp;
}

static void AdjustDynamicGraphicViewer (BioseqViewPtr bvp)

{
  if (bvp == NULL) return;
  if (Visible (bvp->vwr)) {
    if (PictureGrew (bvp->vwr)) {
      PictureHasEnlarged (bvp->vwr);
    }
  }
}

static void ChangeBioseqViewTabs (VoidPtr data, Int2 newval, Int2 oldval)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;
  BioseqPtr          bsp;

  bfp = (BioseqViewFormPtr) data;
  if (bfp != NULL && bfp->bvd.bsp != NULL) {
    WatchCursor ();
    bsp = bfp->bvd.bsp;
    bfp->bvd.scaleNotCalculated = TRUE;
    bfp->bvd.moveToOldPos = FALSE;
    bpp = bfp->currentBioseqPage;
    if (bpp != NULL && bpp->show != NULL) {
      bpp->show (&(bfp->bvd), FALSE);
    }
    Update ();
    bfp->bvd.old_rect_shown = FALSE;
    if (ISA_na (bsp->mol)) {
      bfp->currentNucPage = newval;
    } else if (ISA_aa (bsp->mol)) {
      bfp->currentProtPage = newval;
    }
    SetCurrentPagePointers (bfp);
    PointerToForm (bfp->form, (Pointer) bfp->bvd.bsp);
    SetBioseqImportExportItems (bfp);
    ArrowCursor ();
    Update ();
    AdjustDynamicGraphicViewer (&(bfp->bvd));
  }
}

static void ChangeBioseqDocText (PopuP p)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (p);
  if (bfp != NULL && bfp->bvd.bsp != NULL) {
    WatchCursor ();
    bfp->bvd.moveToOldPos = FALSE;
    bpp = bfp->currentBioseqPage;
    if (bpp != NULL && bpp->show != NULL) {
      bpp->show (&(bfp->bvd), FALSE);
    }
    Update ();
    bfp->bvd.useScrollText = (Boolean) (! bfp->bvd.useScrollText);
    PointerToForm (bfp->form, (Pointer) bfp->bvd.bsp);
    SetBioseqImportExportItems (bfp);
    ArrowCursor ();
    Update ();
    AdjustDynamicGraphicViewer (&(bfp->bvd));
  }
}

static void ChangeBioseqViewGroup (GrouP g)

{
  BioseqViewFormPtr  bfp;
  Int2               val;

  bfp = (BioseqViewFormPtr) GetObjectExtra (g);
  if (bfp != NULL) {
    val = GetValue (g);
    ChangeBioseqViewTabs ((VoidPtr) bfp, val - 1, 0);
  }
}

static void ChangeBioseqViewPopup (PopuP p)

{
  BioseqViewFormPtr  bfp;
  Int2               val;

  bfp = (BioseqViewFormPtr) GetObjectExtra (p);
  if (bfp != NULL) {
    val = GetValue (p);
    ChangeBioseqViewTabs ((VoidPtr) bfp, val - 1, 0);
  }
}

static Boolean ChangeBioseqItemID (GatherContextPtr gcp)

{
  BioseqViewFormPtr  bfp;

  if (gcp == NULL) return TRUE;
  bfp = (BioseqViewFormPtr) gcp->userdata;
  if (bfp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ) {
    if (bfp->bvd.bsp == (BioseqPtr) gcp->thisitem) {
      bfp->input_itemID = gcp->itemID;
      return FALSE;
    }
  }
  return TRUE;
}

static void ChangeTarget (Handle targ)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;
  BioseqPtr          bsp;
  GatherScope        gs;
  SeqEntryPtr        sep;
  Int2               val;

  bfp = (BioseqViewFormPtr) GetObjectExtra (targ);
  if (bfp != NULL) {
    bfp->bvd.viewWholeEntity = FALSE;
    val = GetValue (targ);
    if (val == 1) {
      bfp->bvd.viewWholeEntity = TRUE;
    } else {
      val--;
    }
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep != NULL) {
      sep = FindNthSequinEntry (sep, val);
      if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        /*
        if (bfp->bvd.bsp != bsp) {
          BioseqLock (bsp);
          BioseqUnlock (bfp->bvd.bsp);
        }
        */
        bfp->bvd.bsp = bsp;
        bfp->bvd.scaleNotCalculated = TRUE;
        bfp->bvd.moveToOldPos = FALSE;
        bpp = bfp->currentBioseqPage;
        if (bpp != NULL && bpp->show != NULL) {
          bpp->show (&(bfp->bvd), FALSE);
        }
        Update ();
        if (ISA_na (bsp->mol)) {
          SafeHide (bfp->protViewControl);
          SafeShow (bfp->nucViewControl);
        } else if (ISA_aa (bsp->mol)) {
          SafeHide (bfp->nucViewControl);
          SafeShow (bfp->protViewControl);
        }
        SetCurrentPagePointers (bfp);
        if (bfp->bvd.slp_list != NULL) {
          bfp->bvd.slp_list = free_slp_list (bfp->bvd.slp_list);
        }
        if (bfp->bvd.anp_node != NULL) {
          bfp->bvd.anp_node = FreeAlignNode (bfp->bvd.anp_node);
        }
        if (bfp->bvd.g_list != NULL) {
          bfp->bvd.g_list = ValNodeFreeData (bfp->bvd.g_list);
        }
        if (bfp->bvd.ftype_list != NULL) {
          bfp->bvd.ftype_list = ValNodeFree (bfp->bvd.ftype_list);
        }
        if (bfp->bvd.sentinelList != NULL) {
          bfp->bvd.sentinelList = ValNodeFreeData (bfp->bvd.sentinelList);
        }
        if (bfp->bvd.entityList != NULL) {
          bfp->bvd.entityList = ValNodeFree (bfp->bvd.entityList);
        }
        if (bfp->bvd.viewWholeEntity) {
          bfp->input_itemID = 0;
          PointerToForm (bfp->form, (Pointer) bfp->bvd.bsp);
          SetBioseqImportExportItems (bfp);
          SendMessageToForm (bfp->form, VIB_MSG_CHANGE);
          Update ();
          AdjustDynamicGraphicViewer (&(bfp->bvd));
          if (bfp->controls != NULL && bfp->updateControls != NULL) {
            bfp->updateControls (bfp->controls);
          }
          return;
        }
        MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
        gs.seglevels = 1;
        MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
        gs.ignore[OBJ_BIOSEQ] = FALSE;
        gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
        GatherEntity (bfp->input_entityID, (Pointer) bfp, ChangeBioseqItemID, &gs);
        PointerToForm (bfp->form, (Pointer) bfp->bvd.bsp);
        SetBioseqImportExportItems (bfp);
        SendMessageToForm (bfp->form, VIB_MSG_CHANGE);
        Update ();
        AdjustDynamicGraphicViewer (&(bfp->bvd));
        if (bfp->controls != NULL && bfp->updateControls != NULL) {
          bfp->updateControls (bfp->controls);
        }
      }
    }
  }
}

extern void SetBioseqViewTarget (BaseFormPtr fp, CharPtr seqId)

{
  EnumFieldAssocPtr  ap;
  BioseqViewFormPtr  bfp;
  Int2               val;

  bfp = (BioseqViewFormPtr) fp;
  if (bfp == NULL || StringHasNoText (seqId) || bfp->targetAlist == NULL) return;
  for (ap = bfp->targetAlist, val = 1; ap != NULL && ap->name != NULL; ap++, val++) {
    if (StringICmp (ap->name, seqId) == 0) {
      SetValue (bfp->targetControl, val);
      ChangeTarget ((Handle) bfp->targetControl);
      return;
    }
  }
}

extern BioseqViewPtr GetBioseqViewPtrFromBaseFormPtr (BaseFormPtr fp)

{
  BioseqViewFormPtr  bfp;
  bfp = (BioseqViewFormPtr) fp;
  if (bfp == NULL) return NULL;
  return (&(bfp->bvd));
}

static void ChangeStyle (PopuP p)

{
  BioseqViewFormPtr  bfp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (p);
  if (bfp != NULL) {
    bfp->bvd.moveToOldPos = TRUE;
    PointerToForm (bfp->form, bfp->bvd.bsp);
    Update ();
    AdjustDynamicGraphicViewer (&(bfp->bvd));
  }
}

static void ChangeScale (PopuP p)

{
  BioseqViewFormPtr  bfp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (p);
  if (bfp != NULL) {
    bfp->bvd.moveToOldPos = TRUE;
    PointerToForm (bfp->form, bfp->bvd.bsp);
    Update ();
    AdjustDynamicGraphicViewer (&(bfp->bvd));
  }
}

static void ChangeSalsaControls (PopuP p)

{
  BioseqViewFormPtr  bfp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (p);
  if (bfp != NULL) {
    bfp->bvd.moveToOldPos = TRUE;
    PointerToForm (bfp->form, bfp->bvd.bsp);
    Update ();
    AdjustDynamicGraphicViewer (&(bfp->bvd));
  }
}

static void ChangeFlatFileMode (PopuP p)

{
  BioseqViewFormPtr  bfp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (p);
  if (bfp != NULL) {
    bfp->bvd.moveToOldPos = TRUE;
    PointerToForm (bfp->form, bfp->bvd.bsp);
    Update ();
    AdjustDynamicGraphicViewer (&(bfp->bvd));
  }
}

static void DuplicateViewProc (ButtoN b)

{
  BioseqViewFormPtr  bfp;
  Int2               handled;
  Uint2              itemID;
  SeqViewProcsPtr    svpp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  if (bfp->input_itemtype == OBJ_BIOSEQ) {
    WatchCursor ();
    itemID = bfp->input_itemID;
    if (itemID == 0) {
      itemID = 1;
    }
    svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
    if (svpp != NULL) {
      svpp->forceSeparateViewer = TRUE;
    }
    handled = GatherProcLaunch (OMPROC_VIEW, FALSE, bfp->input_entityID, itemID,
                                OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
    ArrowCursor ();
    if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
      Message (MSG_ERROR, "Unable to launch additional viewer.");
    }
  } else if (bfp->input_itemtype == OBJ_MEDLINE_ENTRY) {
    WatchCursor ();
    itemID = bfp->input_itemID;
    if (itemID == 0) {
      itemID = 1;
    }
    handled = GatherProcLaunch (OMPROC_VIEW, FALSE, bfp->input_entityID, itemID,
                                OBJ_MEDLINE_ENTRY, 0, 0, OBJ_MEDLINE_ENTRY, 0);
    ArrowCursor ();
    if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
      Message (MSG_ERROR, "Unable to launch additional viewer.");
    }
  }
}

static void CleanupBioseqForm (GraphiC g, VoidPtr data)

{
  BioseqViewFormPtr  bfp;
  Uint2              userkey;
  ValNodePtr         vnp;

  bfp = (BioseqViewFormPtr) data;
  if (bfp != NULL) {
    WatchCursor ();
    /*
    BioseqUnlock (bfp->bvd.bsp);
    */
    if (bfp->input_entityID > 0 && bfp->userkey > 0) {
      userkey = bfp->userkey;
      bfp->userkey = 0;
      ObjMgrFreeUserData (bfp->input_entityID, bfp->procid, bfp->proctype, userkey);
      /* this may trigger another remove, hence bfp->userkey first set to 0 */
      for (vnp = bfp->bvd.entityList; vnp != NULL; vnp = vnp->next) {
        if (bfp->input_entityID != (Uint2) vnp->data.intvalue) {
          ObjMgrFreeUserData ((Uint2) vnp->data.intvalue, bfp->procid, bfp->proctype, userkey);
        }
      }
    }
    bfp->bvd.pict = DeletePicture (bfp->bvd.pict);
    if (bfp->bvd.slp_list != NULL) {
      bfp->bvd.slp_list = free_slp_list (bfp->bvd.slp_list);
    }
    if (bfp->bvd.anp_node != NULL) {
      bfp->bvd.anp_node = FreeAlignNode (bfp->bvd.anp_node);
    }
    if (bfp->bvd.g_list != NULL) {
      bfp->bvd.g_list = ValNodeFreeData (bfp->bvd.g_list);
    }
    if (bfp->bvd.ftype_list != NULL) {
      bfp->bvd.ftype_list = ValNodeFree (bfp->bvd.ftype_list);
    }
    if (bfp->bvd.sentinelList != NULL) {
      bfp->bvd.sentinelList = ValNodeFreeData (bfp->bvd.sentinelList);
    }
    if (bfp->bvd.entityList != NULL) {
      bfp->bvd.entityList = ValNodeFree (bfp->bvd.entityList);
    }
    if (bfp->cleanupObjectPtr && bfp->objectDataPtr != NULL) {
      SeqEntryFree ((SeqEntryPtr) bfp->objectDataPtr);
    }
    bfp->bioseqNucPageList = BioseqPageListFree (bfp->bioseqNucPageList);
    bfp->bioseqProtPageList = BioseqPageListFree (bfp->bioseqProtPageList);
    bfp->targetAlist = FreeEnumFieldAlist (bfp->targetAlist);
    ArrowCursor ();
  }
  StdCleanupFormProc (g, data);
}

static void WriteEntityProc (Uint2 entityID, CharPtr path)

{
  if (entityID > 0 && path != NULL && *path != '\0') {
    ObjMgrSetDirtyFlag (entityID, FALSE);
  }
}

static void CloseViewFormProc (WindoW w)

{
  BioseqViewFormPtr  bfp;
  Char               path [PATH_MAX];
  MsgAnswer          response;

  bfp = (BioseqViewFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    if (ObjMgrGetDirtyFlag (bfp->input_entityID)) {
      response = Message (MSG_YNC, "Do you wish to save the changes?");
      if (response == ANS_YES) {
        if (GetOutputFileName (path, sizeof (path), NULL)) {
          WriteEntityProc (bfp->input_entityID, path);
        } else {
          return;
        }
      } else if (response == ANS_CANCEL) {
        return;
      }
    }
  }
  RemoveSeqEntryViewer ((ForM) w);
}

static Boolean ExportBioseqViewForm (ForM f, CharPtr filename)

{
  ExportBioseqViewFormToFile (GetObjectExtra (f), filename);
  return TRUE;
}

static void BioseqFormMessage (ForM f, Int2 mssg)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;
  SeqViewProcsPtr    svpp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_REDRAW :
        svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
        if (svpp != NULL) {
          bfp->bvd.displayFont = svpp->displayFont;
          bpp = bfp->currentBioseqPage;
          if (bpp != NULL && bpp->show != NULL) {
            bpp->show (&(bfp->bvd), FALSE);
          }
          PointerToForm (bfp->form, bfp->bvd.bsp);
          Update ();
          AdjustDynamicGraphicViewer (&(bfp->bvd));
        }
        break;
      case VIB_MSG_COPY :
        CopyBioseqViewFormToClipboard (bfp);
        break;
      case VIB_MSG_EXPORT :
        ExportBioseqViewForm (f, NULL);
        break;
      case VIB_MSG_PRINT :
        PrintBioseqViewForm (bfp);
        break;
      default :
        if (bfp->appmessage != NULL) {
          bfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static Handle CreateViewControl (GrouP g, BioseqViewFormPtr bfp,
                                 SeqViewProcsPtr svpp, BioseqPagePtr bpp,
                                 Int2 page, Int2 pixwidth)


{
  Int2          count;
  Int2          j;
  GrouP         k;
  PopuP         pops;
  PrompT        ppt;
  Int2          radiowidth;
  GrouP         rads;
  DialoG        tbs;
  CharPtr PNTR  titles;
  BioseqPagePtr tmp;
  Int2          wid;

  if (bfp != NULL & svpp != NULL && bpp != NULL) {
    count = 0;
    for (tmp = bpp; tmp != NULL; tmp = tmp->next) {
      count++;
    }
    titles = MemNew (sizeof (CharPtr PNTR) * (count + 1));
    count = 0;
    for (tmp = bpp; tmp != NULL; tmp = tmp->next) {
      titles [count] = tmp->label;
      count++;
    }
    switch (svpp->useFolderTabs) {
      case CHANGE_VIEW_NOTABS :
        return NULL;
      case CHANGE_VIEW_FOLDERTABS :
        tbs = CreateFolderTabs (g, titles, page,
                                0, 0, SYSTEM_FOLDER_TAB,
                                ChangeBioseqViewTabs, (Pointer) bfp);
        return (Handle) tbs;
      case CHANGE_VIEW_TEXTTABS :
        tbs = CreateTextTabs (g, titles, page,
                              0, 0, SYSTEM_TEXT_TAB,
                              ChangeBioseqViewTabs, (Pointer) bfp);
        return (Handle) tbs;
      case CHANGE_VIEW_RADIOBUTTONS :
        k = HiddenGroup (g, -2, 0, NULL);
        ppt = StaticPrompt (k, "Format:", 0, 0, programFont, 'l');
        SelectFont (programFont);
        pixwidth -= StringWidth ("Format:");
        SelectFont (systemFont);
        wid = 0;
        radiowidth = 0;
        for (j = 0; titles [j] != NULL; j++) {
          wid += StringWidth (titles [j]);
#ifdef WIN_MOTIF
          pixwidth -= 25;
#else
          pixwidth -= 20;
#endif
          radiowidth++;
        }
        if (wid > pixwidth) {
          radiowidth = (radiowidth + 1) / 2;
        }
        rads = HiddenGroup (k, radiowidth, 0, ChangeBioseqViewGroup);
        SetObjectExtra (rads, (Pointer) bfp, NULL);
        for (j = 0; titles [j] != NULL; j++) {
          RadioButton (rads, titles [j]);
        }
        SetValue (rads, page + 1);
        AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) rads, NULL);
        return (Handle) rads;
      case CHANGE_VIEW_POPUP :
        k = HiddenGroup (g, -2, 0, NULL);
        StaticPrompt (k, "Display Format", 0, popupMenuHeight, programFont, 'l');
        pops = PopupList (k, TRUE, ChangeBioseqViewPopup);
        SetObjectExtra (pops, (Pointer) bfp, NULL);
        for (j = 0; titles [j] != NULL; j++) {
          PopupItem (pops, titles [j]);
        }
        SetValue (pops, page + 1);
        return (Handle) pops;
      default :
        break;
    }
    MemFree (titles);
  }
  return NULL;
}

static void BioseqViewFormActivate (WindoW w)

{
  BioseqViewFormPtr  bfp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    if (bfp->activate != NULL) {
      bfp->activate (w);
    }
    SetBioseqImportExportItems (bfp);
    if (bfp->bvd.legendOK) {
      SafeEnable (bfp->bvd.legendItem);
    } else {
      SafeDisable (bfp->bvd.legendItem);
    }
  }
}

static void CopyBioseqReportSpecs (BioseqPagePtr bpp, BioseqPagePtr PNTR head,
                                   Int4 length, Boolean hasAlignments,
                                   Boolean isGenome, Boolean nucOK, Boolean protOK)

{
  Boolean  okay;

  if (bpp == NULL || head == NULL) return;
  while (bpp != NULL) {
    okay = FALSE;
    if (isGenome) {
      if (bpp->genomeOK) {
        okay = TRUE;
      }
    } else if (nucOK && bpp->nucOK) {
      okay = TRUE;
    } else if (protOK && bpp->protOK) {
      okay = TRUE;
    }
    if (okay) {
      if (bpp->maxLength < 0 || length <= bpp->maxLength) {
        if (hasAlignments || (! bpp->needAlignment)) {
          AddBioseqPageToList (head, bpp);
        }
      }
    }
    bpp = bpp->next;
  }
}

static Boolean HasAlignments (BioseqPtr bsp)

{
  return TRUE;
}

/*
extern void Nlm_DisplayEnvironmentVariables (void);
static void EnviroProc (ButtoN b)

{
  Nlm_DisplayEnvironmentVariables ();
}
*/

static ForM LIBCALL CreateNewSeqEntryViewFormEx (Int2 left, Int2 top, CharPtr title,
                                                 BioseqPtr bsp, SeqViewProcsPtr svpp,
                                                 Uint2 entD, Uint2 itemID, Uint2 itemtype,
                                                 Boolean smart)

{
  ButtoN               b;
  BioseqViewFormPtr    bfp;
  BioseqPagePtr        bpp;
  WndActnProc          close;
  Int4                 count;
  GrouP                d;
  ButtoN               dp;
  Int2                 delta;
  Uint2                entityID;
  FonT                 fnt;
  GrouP                g;
  GrouP                h;
  Boolean              hasAlignments;
  Int2                 i;
  Int2                 j;
  GrouP                k;
  Int4                 length;
  SeqViewFetchAlignsProc  makeAlignBtn = NULL;
  SeqViewControlsProc  makeControls = NULL;
  Int2                 mssg;
  Int2                 numStyles;
  Int2                 pixheight;
  Int2                 pixwidth;
  GrouP                pnlGrp;
  PrompT               ppt;
  PoinT                pt;
  Int2                 pty;
  GrouP                q;
  RecT                 r;
  RecT                 r1;
  RecT                 r2;
  RecT                 r3;
  GrouP                s;
  SeqEntryPtr          sep;
  CharPtr              str;
  CharPtr              styleName;
  Int2                 val;
  WindoW               w;
  PopuP                x;
  GrouP                y;
  ButtoN               z;

  w = NULL;
  bfp = (BioseqViewFormPtr) MemNew (sizeof (BioseqViewForm));
  if (bfp != NULL /* && bsp != NULL */ ) {
    close = CloseViewFormProc;
    if (svpp == NULL) {
      svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
    }
    if (svpp != NULL && svpp->closeForm != NULL) {
      close = svpp->closeForm;
    }
    w = DocumentWindow (left, top, -10, -10, title,
                        close, ResizeViewForm);
    SetObjectExtra (w, bfp, CleanupBioseqForm);
    bfp->form = (ForM) w;
    bfp->actproc = NULL;
    bfp->toform = BioseqPtrToBioseqForm;
    bfp->exportform = ExportBioseqViewForm;
    bfp->formmessage = BioseqFormMessage;

    bfp->input_entityID = entD;
    bfp->input_itemID = itemID;
    bfp->input_itemtype = itemtype;
    bfp->this_itemtype = OBJ_BIOSEQ;
    if (bsp != NULL) {
      bfp->this_subtype = bsp->repr;
    } else {
      bfp->this_subtype = Seq_repr_raw;
    }

    bfp->bvd.form = bfp->form;
    bfp->bvd.bsp = bsp;
    /*
    BioseqLock (bsp);
    */
  
    if (svpp != NULL) {
      bfp->cleanupObjectPtr = svpp->cleanupObjectPtr;
      bfp->appmessage = svpp->handleMessages;
      bfp->activateForm = svpp->activateForm;

      if (svpp->allowScrollText) {
        bfp->bvd.useScrollText = svpp->startInScrollText;
      }
      bfp->bvd.launchEditors = svpp->launchEditors;
      bfp->bvd.launchSubviewers = svpp->launchSubviewers;
      bfp->bvd.sendSelectMessages = svpp->sendSelectMessages;
      bfp->bvd.highlightSelections = svpp->highlightSelections;
      bfp->bvd.hasTargetControl = svpp->hasTargetControl;
      bfp->bvd.displayFont = svpp->displayFont;

      if (! StringHasNoText (svpp->filepath)) {
        bfp->filepath = StringSave (svpp->filepath);
      }

      if (bsp != NULL) {
        hasAlignments = HasAlignments (bsp);
        entityID = ObjMgrGetEntityIDForPointer (bsp);
        sep = GetTopSeqEntryForEntityID (entityID);
        bfp->bvd.isGenome = IsAGenomeRecord (sep);
        length = bsp->length;
      } else {
        hasAlignments = TRUE;
        entityID = 0;
        sep = NULL;
        bfp->bvd.isGenome = FALSE;
        length = 0;
      }
      CopyBioseqReportSpecs (svpp->pageSpecs, &(bfp->bioseqNucPageList),
                             length, hasAlignments, bfp->bvd.isGenome,
                             TRUE, FALSE);
      CopyBioseqReportSpecs (svpp->pageSpecs, &(bfp->bioseqProtPageList),
                             length, hasAlignments, bfp->bvd.isGenome,
                             FALSE, TRUE);
      if (bfp->bioseqNucPageList == NULL && bfp->bioseqProtPageList == NULL) {
        Message (MSG_ERROR, "No acceptable report forms are currently registered");
      }
      makeControls = svpp->makeControls;
      bfp->updateControls = svpp->updateControls;
      makeAlignBtn = svpp->makeAlignBtn;
      bfp->updateCounts = svpp->updateCounts;
    }

    if (svpp != NULL && svpp->createMenus != NULL) {
      svpp->createMenus (w);
    }

    g = HiddenGroup (w, -1, 0, NULL);
    SetObjectExtra (g, bfp, NULL);
    SetGroupSpacing (g, 3, 10);

    /*
    if (svpp != NULL && svpp->createToolBar != NULL) {
      svpp->createToolBar (g);
    }
    */

    fnt = programFont;
    if (bfp->bvd.displayFont != NULL) {
      fnt = bfp->bvd.displayFont;
    }
    SelectFont (fnt);
    pixwidth = MIN ((Int2) (80 * CharWidth ('0') + 18),
                    (Int2) (screenRect.right - screenRect.left - 20));
    SelectFont (systemFont);
    pixheight = 14 * stdLineHeight;
    bfp->currentBioseqPage = NULL;
    bfp->currentNucPage = 0;
    bfp->currentProtPage = 0;
    bfp->targetAlist = NULL;

    /*
    bfp->ffstyle = GENBANK_STYLE;
    bfp->anp_node = NULL;
    if (Nlm_GetAppProperty ("SequinUseEMBLStyle") != NULL) {
      bfp->ffstyle = EMBL_STYLE;
    } else if (Nlm_GetAppProperty ("SequinUseDDBJStyle") != NULL) {
      bfp->ffstyle = DDBJ_STYLE;
    }
    */
    b = NULL;
    k = NULL;
    if (svpp != NULL) {
      if (bfp->bvd.hasTargetControl) {
        k = HiddenGroup (g, -4, 0, NULL);
        ppt = StaticPrompt (k, "Target Sequence", 0, 0, programFont, 'l');
        count = 0;
        if (bsp != NULL) {
          entityID = ObjMgrGetEntityIDForPointer (bfp->bvd.bsp);
          sep = GetTopSeqEntryForEntityID (entityID);
          count = SequinEntryCount (sep);
        }
        if (bsp != NULL && count < 32 && (! smart)) {
          bfp->usePopupForTarget = TRUE;
          bfp->targetControl = (Handle) PopupList (k, TRUE, (PupActnProc) ChangeTarget);
        } else {
          bfp->usePopupForTarget = FALSE;
          bfp->targetControl = (Handle) SingleList (k, 14, 3, (LstActnProc) ChangeTarget);
        }
        SetObjectExtra (bfp->targetControl, (Pointer) bfp, NULL);
        val = PopulateTarget (bfp);
        SetValue (bfp->targetControl, val + 1);
        if (svpp->hasDoneButton) {
          b = PushButton (k, "Done", StdSendAcceptButtonMessageProc);
          SetObjectExtra (b, bfp, NULL);
        }
        AlignObjects (ALIGN_VERTICAL, (HANDLE) ppt, (HANDLE) bfp->targetControl, NULL);
        GetPosition (bfp->targetControl, &r1);
        GetPosition (b, &r2);
        delta = (r1.bottom - r1.top) - (r2.bottom - r2.top);
        if (delta > 0) {
          OffsetRect (&r2, 0, delta / 2);
          SetPosition (b, &r2);
          AdjustPrnt (b, &r2, FALSE);
        }
        /* PushButton (k, "Enviro", EnviroProc); */
      }
      bfp->currentNucPage = svpp->initNucPage;
      bfp->currentProtPage = svpp->initProtPage;

      str = NULL;
      bpp = NULL;
      if (bfp->bvd.isGenome) {
        if (svpp->initGenomeLabel != NULL) {
          bpp = bfp->bioseqNucPageList;
          str = svpp->initGenomeLabel;
        }
      } else {
        if (svpp->initNucLabel != NULL) {
          bpp = bfp->bioseqNucPageList;
          str = svpp->initNucLabel;
        }
      }
      if (str != NULL && bpp != NULL) {
        j = 0;
        while (bpp != NULL) {
          if (StringICmp (bpp->label, str) == 0) {
            bfp->currentNucPage = j;
          }
          bpp = bpp->next;
          j++;
        }
      }
      if (svpp->initProtLabel != NULL) {
        bpp = bfp->bioseqProtPageList;
        str = svpp->initProtLabel;
        j = 0;
        while (bpp != NULL) {
          if (StringICmp (bpp->label, str) == 0) {
            bfp->currentProtPage = j;
          }
          bpp = bpp->next;
          j++;
        }
      }

      k = HiddenGroup (g, -2, 0, NULL);
      if (bfp->bvd.hasTargetControl || bsp == NULL) {
        q = HiddenGroup (k, 0, 0, NULL);
        bfp->nucViewControl = CreateViewControl (q, bfp, svpp, bfp->bioseqNucPageList,
                                                 bfp->currentNucPage, pixwidth);
        bfp->protViewControl = CreateViewControl (q, bfp, svpp, bfp->bioseqProtPageList,
                                                  bfp->currentProtPage, pixwidth);
      } else if (ISA_na (bsp->mol)) {
        bfp->nucViewControl = CreateViewControl (k, bfp, svpp, bfp->bioseqNucPageList,
                                                 bfp->currentNucPage, pixwidth);
      } else if (ISA_aa (bsp->mol)) {
        bfp->protViewControl = CreateViewControl (k, bfp, svpp, bfp->bioseqProtPageList,
                                                  bfp->currentProtPage, pixwidth);
      }
      if (bsp == NULL || ISA_na (bsp->mol)) {
        Hide (bfp->protViewControl);
      } else {
        Hide (bfp->nucViewControl);
      }
      if (bfp->nucViewControl != NULL) {
        ObjectRect (bfp->nucViewControl, &r);
        pixwidth = MAX (pixwidth, r.right - 2 * r.left - Nlm_vScrollBarWidth);
      } else if (bfp->protViewControl != NULL) {
        ObjectRect (bfp->protViewControl, &r);
        pixwidth = MAX (pixwidth, r.right - 2 * r.left - Nlm_vScrollBarWidth);
      }
      pixwidth = MAX (pixwidth, svpp->minPixelWidth);
      pixheight = MAX (pixheight, svpp->minPixelHeight);
      if (bfp->nucViewControl != NULL && Visible (bfp->nucViewControl)) {
        bfp->currentBioseqPage = bfp->bioseqNucPageList;
      } else if (bfp->protViewControl != NULL && Visible (bfp->protViewControl)) {
        bfp->currentBioseqPage = bfp->bioseqProtPageList;
      } else if (svpp != NULL) {
        bfp->currentBioseqPage = svpp->pageSpecs;
      }
    }

    if (svpp != NULL && svpp->useFolderTabs == CHANGE_VIEW_POPUP && k != NULL) {
      d = HiddenGroup (k, -2, 0, NULL);
      if (svpp->hasDuplicateButton) {
        dp = PushButton (d, "Duplicate", DuplicateViewProc);
        SetObjectExtra (dp, bfp, NULL);
      }
      h = HiddenGroup (d, 0, 0, NULL);
    } else {
      d = HiddenGroup (g, -2, 0, NULL);
      if (svpp->hasDuplicateButton) {
        dp = PushButton (d, "Duplicate", DuplicateViewProc);
        SetObjectExtra (dp, bfp, NULL);
      }
      h = HiddenGroup (d, 0, 0, NULL);
    }

    x = NULL;
    bfp->bvd.docTxtControlGrp = HiddenGroup (h, -4, 0, NULL);
    if (svpp != NULL && svpp->allowScrollText) {
      StaticPrompt (bfp->bvd.docTxtControlGrp, "Display Type", 0, popupMenuHeight, programFont, 'l');
      x = PopupList (bfp->bvd.docTxtControlGrp, TRUE, ChangeBioseqDocText);
      SetObjectExtra (x, bfp, NULL);
      PopupItem (x, "Document");
      PopupItem (x, "Text");
    }
    if (bfp->bvd.hasTargetControl && GetAppProperty ("InternalNcbiSequin") != NULL) {
      bfp->bvd.modeControlGrp = HiddenGroup (bfp->bvd.docTxtControlGrp, -4, 0, NULL);
      StaticPrompt (bfp->bvd.modeControlGrp, "Mode", 0, popupMenuHeight, programFont, 'l');
      bfp->bvd.ffModeCtrl = PopupList (bfp->bvd.modeControlGrp, TRUE, ChangeFlatFileMode);
      SetObjectExtra (bfp->bvd.ffModeCtrl, bfp, NULL);
      PopupItem (bfp->bvd.ffModeCtrl, "Sequin");
      PopupItem (bfp->bvd.ffModeCtrl, "Release");
      SetValue (bfp->bvd.ffModeCtrl, 1);
    }
    Hide (bfp->bvd.docTxtControlGrp);

    s = HiddenGroup (h, -4, 0, NULL);
    bfp->bvd.styleControlGrp = HiddenGroup (s, -6, 0, NULL);
    StaticPrompt (bfp->bvd.styleControlGrp, "Style", 0, popupMenuHeight, programFont, 'l');
    bfp->bvd.style = PopupList (bfp->bvd.styleControlGrp, TRUE, ChangeStyle);
    SetObjectExtra (bfp->bvd.style, bfp, NULL);
    numStyles = GetMuskTotalSt ();
    for (i = 0; i < numStyles; i++) {
      styleName = GetMuskStyleName (i);
      if (StringCmp (styleName, "StyleX") != 0) {
        PopupItem (bfp->bvd.style, styleName);
      }
    }
    SetValue (bfp->bvd.style, GetMuskCurrentSt () + 1);

    y = HiddenGroup (s, 0, 0, NULL);
    bfp->bvd.scaleControlGrp = HiddenGroup (y, -6, 0, NULL);
    StaticPrompt (bfp->bvd.scaleControlGrp, "Scale", 0, popupMenuHeight, programFont, 'l');
    bfp->bvd.scale = PopupList (bfp->bvd.scaleControlGrp, TRUE, ChangeScale);
    SetObjectExtra (bfp->bvd.scale, bfp, NULL);

    bfp->bvd.findGeneGrp = HiddenGroup (y, -4, 0, NULL);
    z = PushButton (bfp->bvd.findGeneGrp, "Find by Gene or Product", ShowGeneList);
    SetObjectExtra (z, (Pointer) &(bfp->bvd), NULL);

    Hide (bfp->bvd.styleControlGrp);
    Hide (bfp->bvd.scaleControlGrp);
    Hide (bfp->bvd.findGeneGrp);

    if (b != NULL) {
      GetPosition (b, &r1);
      GetPosition (bfp->bvd.scale, &r2);
      delta = r2.left - r1.left;
      if (delta > 0) {
        OffsetRect (&r1, delta, 0);
        SetPosition (b, &r1);
        AdjustPrnt (b, &r1, FALSE);
      }
    }

    if (bfp->bvd.launchEditors) {
      bfp->bvd.clickMe = StaticPrompt (g, "Double click on an item to launch the appropriate editor.",
                                       0, 0, programFont, 'l');
      Hide (bfp->bvd.clickMe);
    }

    h = HiddenGroup (g, 0, 0, NULL);
    SetGroupMargins (h, 1, 1);

    bfp->bvd.vwr = CreateViewer (h, pixwidth, pixheight, TRUE, TRUE);
    SetObjectExtra (bfp->bvd.vwr, (Pointer) &(bfp->bvd), NULL);
    bfp->bvd.pict = NULL;
    bfp->bvd.scaleNotCalculated = TRUE;
    bfp->bvd.moveToOldPos = FALSE;
    bfp->bvd.minIndex = 1;
    bfp->bvd.expansion = 1;
    Hide (bfp->bvd.vwr);

    fnt = programFont;
    if (bfp->bvd.displayFont != NULL) {
      fnt = bfp->bvd.displayFont;
    }
    bfp->bvd.text = ScrollText (h, pixwidth / stdCharWidth,
                                pixheight / stdLineHeight, fnt, FALSE, NULL);
    SetObjectExtra (bfp->bvd.text, (Pointer) &(bfp->bvd), NULL);
    Hide (bfp->bvd.text);

    bfp->bvd.doc = DocumentPanel (h, pixwidth, pixheight);
    SetObjectExtra (bfp->bvd.doc, (Pointer) &(bfp->bvd), NULL);
    Hide (bfp->bvd.doc);

    bfp->bvd.pnlParentGrp = HiddenGroup (h, -1, 0, NULL);
    GetNextPosition (bfp->bvd.pnlParentGrp, &pt);
    pty = pt.y;
    pnlGrp = HiddenGroup (bfp->bvd.pnlParentGrp, 6, 0, NULL);
    StaticPrompt (pnlGrp, "Sequences", 0, popupMenuHeight, programFont, 'l');
    bfp->bvd.seqControl = PopupList (pnlGrp, TRUE, ChangeSalsaControls);
    SetObjectExtra (bfp->bvd.seqControl, bfp, NULL);
    PopupItem (bfp->bvd.seqControl, "Target");
    PopupItem (bfp->bvd.seqControl, "Aligned");
    SetValue (bfp->bvd.seqControl, 1);
    StaticPrompt (pnlGrp, "Features", 0, popupMenuHeight, programFont, 'l');
    bfp->bvd.featControl = PopupList (pnlGrp, TRUE, ChangeSalsaControls);
    SetObjectExtra (bfp->bvd.featControl, bfp, NULL);
    PopupItem (bfp->bvd.featControl, "None");
    PopupItem (bfp->bvd.featControl, "Target");
    PopupItem (bfp->bvd.featControl, "Aligned");
    SetValue (bfp->bvd.featControl, 2);
    StaticPrompt (pnlGrp, "Numbering", 0, popupMenuHeight, programFont, 'l');
    bfp->bvd.numControl = PopupList (pnlGrp, TRUE, ChangeSalsaControls);
    SetObjectExtra (bfp->bvd.numControl, bfp, NULL);
    PopupItem (bfp->bvd.numControl, "None");
    PopupItem (bfp->bvd.numControl, "Side");
    PopupItem (bfp->bvd.numControl, "Top");
    SetValue (bfp->bvd.numControl, 3);
    GetNextPosition (bfp->bvd.pnlParentGrp, &pt);
    bfp->bvd.pnl = SalsaTextPanel (bfp->bvd.pnlParentGrp, pixwidth, pixheight - (pt.y - pty));
    SetObjectExtra (bfp->bvd.pnl, (Pointer) &(bfp->bvd), NULL);
    Hide (bfp->bvd.pnl);
    Hide (bfp->bvd.pnlParentGrp);

    if (makeControls != NULL) {
      if (bsp == NULL) {
        bfp->docuid = 0;
        bfp->doctype = TYP_NT;
        bfp->controls = makeControls (g, (BaseFormPtr) bfp, bfp->doctype, bfp->docuid);
      } else {
        bfp->docuid = GetUidFromBsp (bsp);
        if (bfp->docuid > 0) {
          if (ISA_na (bsp->mol)) {
            bfp->doctype = TYP_NT;
          } else if (ISA_aa (bsp->mol)) {
            bfp->doctype = TYP_AA;
          }
          bfp->controls = makeControls (g, (BaseFormPtr) bfp, bfp->doctype, bfp->docuid);
        } else if (makeAlignBtn != NULL) {
          bfp->retrieveAlignments = makeAlignBtn (g, (BaseFormPtr) bfp);
          if (bfp->retrieveAlignments != NULL && bfp->updateCounts != NULL) {
            entityID = ObjMgrGetEntityIDForPointer (bsp);
            sep = GetTopSeqEntryForEntityID (entityID);
            bfp->hasaligns = bfp->updateCounts (bfp->retrieveAlignments, sep);
          }
        }
      }
    }

    AlignObjects (ALIGN_CENTER, (HANDLE) bfp->bvd.text, (HANDLE) bfp->bvd.doc, NULL);
    AlignObjects (ALIGN_RIGHT, (HANDLE) bfp->bvd.vwr, (HANDLE)  bfp->bvd.text,
                  (HANDLE) bfp->bvd.doc, (HANDLE) bfp->bvd.pnl, NULL);
    AlignObjects (ALIGN_LOWER, (HANDLE) bfp->bvd.vwr, (HANDLE)  bfp->bvd.text,
                  (HANDLE) bfp->bvd.doc, (HANDLE) bfp->bvd.pnl, NULL);

    GetPosition (bfp->bvd.vwr, &r3);
    AdjustPrnt (bfp->bvd.vwr, &r3, FALSE);
    GetPosition (bfp->bvd.text, &r3);
    AdjustPrnt (bfp->bvd.text, &r3, FALSE);
    GetPosition (bfp->bvd.doc, &r3);
    AdjustPrnt (bfp->bvd.doc, &r3, FALSE);
    GetPosition (bfp->bvd.pnl, &r3);
    AdjustPrnt (bfp->bvd.pnl, &r3, FALSE);

    if (bfp->bvd.useScrollText) {
      SetValue (x, 2);
    } else {
      SetValue (x, 1);
    }

    RealizeWindow (w);

    mssg = RegisterFormMenuItemName ("SequinLegendItem");
    bfp->bvd.legendItem = FindFormMenuItem ((BaseFormPtr) bfp, mssg);
    bfp->activate = NULL;
    if (svpp != NULL) {
      bfp->activate = svpp->activateForm;
    }
    SetActivate (w, BioseqViewFormActivate);
    Update ();
    BioseqViewFormActivate ((WindoW) bfp->form);

    SendMessageToForm (bfp->form, VIB_MSG_INIT);
    SetCurrentPagePointers (bfp);
    if (bfp->input_entityID > 0 && bfp->bvd.hasTargetControl) {
      ChangeTarget ((Handle) bfp->targetControl); /* shows correct page and populates */
    } else {
      PointerToForm (bfp->form, bfp->bvd.bsp); /* shows correct page and populates */
    }
    Update ();
    AdjustDynamicGraphicViewer (&(bfp->bvd));
  }
  return (ForM) w;
}

extern ForM LIBCALL CreateNewSeqEntryViewForm (Int2 left, Int2 top, CharPtr title,
                                               BioseqPtr bsp, SeqViewProcsPtr svpp)

{
  return CreateNewSeqEntryViewFormEx (left, top, title, bsp, svpp, 0, 0, 0, FALSE);
}

static void ShowFeatLegend (IteM i)

{
  BioseqViewFormPtr  bfp;
  WindoW nw;
  SegmenT pic;
  Int2 style;
  CharPtr style_name;
  VieweR viewer;
  BoxInfo box_i;
  Int2 width, height;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = (BioseqViewFormPtr) GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  if (bfp->bvd.ftype_list == NULL) return;
  style = GetValue (bfp->bvd.style) - 1;
  style_name = GetMuskStyleName (style);

  pic = pic_for_f_legend (bfp->bvd.ftype_list, style_name, 120);
  if (pic == NULL) return;
  SegmentBox (pic, &box_i);
  width = (box_i.right - box_i.left) + 20;
  height = (box_i.top - box_i.bottom) + 10;
  nw = FixedWindow (-50, -33, -10, -10, "Feature Legend", StdCloseWindowProc);
  viewer = CreateViewer (nw, width, height, FALSE, FALSE);
  AttachPicture (viewer, pic, 0, 0, UPPER_LEFT, 1, 1, NULL);
  Show (nw);
  Select (nw);
}

extern IteM CreateLegendItem (MenU m, BaseFormPtr bfp)

{
  IteM  i;
  Int2  mssg;

  i = CommandItem (m, "Legend...", ShowFeatLegend);
  SetObjectExtra (i, bfp, NULL);
  mssg = RegisterFormMenuItemName ("SequinLegendItem");
  SetFormMenuItem (bfp, mssg, i);
  return i;
}

extern void EnableDisableLegendItem (BioseqViewPtr bvp, Boolean enable)

{
  if (bvp == NULL) return;
  bvp->legendOK = enable;
  if (enable) {
    SafeEnable (bvp->legendItem);
  } else {
    SafeDisable (bvp->legendItem);
  }
}

static Boolean SetClickmeTitle (GatherContextPtr gcp)

{
  BioseqViewFormPtr  bfp;
  Char               buf [80];
  ObjMgrPtr          omp;
  ObjMgrTypePtr      omtp;

  if (gcp == NULL || gcp->thisitem == NULL) {
    return FALSE;
  }
  bfp = (BioseqViewFormPtr) gcp->userdata;
  if (bfp == NULL) return FALSE;
  buf [0] = '\0';
  omp = ObjMgrGet ();
  if (omp != NULL) {
    omtp = ObjMgrTypeFind (omp, gcp->thistype, NULL, NULL);
    if (omtp != NULL && omtp->labelfunc != NULL) {
      (*(omtp->labelfunc)) (gcp->thisitem, buf, sizeof (buf) - 1, OM_LABEL_BOTH);
    } else if (gcp->thistype == OBJ_BIOSEQ_SEG) {
      SeqLocLabel (gcp->thisitem, buf, sizeof (buf) - 1, OM_LABEL_BOTH);
    }
  }
  if (! StringHasNoText (buf)) {
     SafeSetTitle (bfp->bvd.clickMe, buf);
    return TRUE;
  }
  /* SafeSetTitle (bfp->bvd.clickMe, "Double click on an item to launch the appropriate editor."); */
  SafeSetTitle (bfp->bvd.clickMe, "");
  return FALSE;
}

extern Boolean InBioseqViewEntityList (Uint2 entityID, BioseqViewPtr bvp)

{
  ValNodePtr  vnp;

  if (entityID == 0 || bvp == NULL) return FALSE;
  for (vnp = bvp->entityList; vnp != NULL; vnp = vnp->next) {
    if (entityID == (Uint2) vnp->data.intvalue) {
      return TRUE;
    }
  }
  return FALSE;
}

extern Int2 LIBCALLBACK BioseqViewMsgFunc (OMMsgStructPtr ommsp)

{
  BioseqViewFormPtr  bfp;
  BioseqPagePtr      bpp;
  BioseqPtr          bsp;
  Boolean            do_refresh;
  ObjMgrDataPtr      omdp;
  OMUserDataPtr      omudp;
  SelStructPtr       sel;
  SeqEntryPtr        sep;
  SeqLocPtr          slp;
  Int2               val;

  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  bfp = (BioseqViewFormPtr) omudp->userdata.ptrvalue;
  if (bfp == NULL) return OM_MSG_RET_ERROR;
  bsp = bfp->bvd.bsp;
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  bpp = bfp->currentBioseqPage;
  if (bpp == NULL) return OM_MSG_RET_ERROR;

  if (ommsp->entityID != bfp->input_entityID) {
    if (! InBioseqViewEntityList (ommsp->entityID, &(bfp->bvd))) return OM_MSG_RET_OK;
  }

  do_refresh = FALSE;
  switch (ommsp->message) {
    case OM_MSG_DEL:
      if (ommsp->entityID == bfp->input_entityID) {
        do_refresh = TRUE;
        omdp = ObjMgrGetData (ommsp->entityID);
        if (omdp != NULL) {
          if (ObjMgrWholeEntity (omdp, ommsp->itemID, ommsp->itemtype)) {
            if (bfp != NULL) {
              RemoveSeqEntryViewer (bfp->form);
            }
            return OM_MSG_RET_OK;
          }
        }
      }
      break;
    case OM_MSG_CREATE:
      if (ommsp->entityID == bfp->input_entityID) {
        do_refresh = TRUE;
      }
      break;
    case OM_MSG_UPDATE:
      do_refresh = TRUE;
      break;
    case OM_MSG_SELECT:
      if (! bfp->bvd.highlightSelections) return OM_MSG_RET_OK;
      if (bpp->highlight == NULL) return OM_MSG_RET_OK;
      slp = NULL;
      if (ommsp->regiontype == OM_REGION_SEQLOC) {
        slp = (SeqLocPtr) ommsp->region;
      }
      bpp->highlight (&(bfp->bvd), ommsp->entityID, ommsp->itemID, ommsp->itemtype, slp, TRUE, TRUE);
      SendMessageToForm (bfp->form, VIB_MSG_SELECT);
      break;
    case OM_MSG_DESELECT:
      if (! bfp->bvd.highlightSelections) return OM_MSG_RET_OK;
      if (bpp->highlight == NULL) return OM_MSG_RET_OK;
      slp = NULL;
      if (ommsp->regiontype == OM_REGION_SEQLOC) {
        slp = (SeqLocPtr) ommsp->region;
      }
      bpp->highlight (&(bfp->bvd), ommsp->entityID, ommsp->itemID, ommsp->itemtype, slp, FALSE, FALSE);
      SendMessageToForm (bfp->form, VIB_MSG_SELECT);
      break;
    case OM_MSG_CACHED:
      break;
    case OM_MSG_UNCACHED:
      break;
    case OM_MSG_TO_CLIPBOARD:
      break;
    default :
      break;
  }
  if (do_refresh) {
    if (bfp->bvd.anp_node != NULL) {
      bfp->bvd.anp_node = FreeAlignNode (bfp->bvd.anp_node);
    }
    if (bfp->targetControl) {
      bfp->bvd.viewWholeEntity = FALSE;
      val = GetValue (bfp->targetControl);
      if (val == 1) {
        bfp->bvd.viewWholeEntity = TRUE;
      } else {
        val--;
      }
      sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
      if (sep != NULL) {
        sep = FindNthSequinEntry (sep, val);
        if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          /*
          bfp->scaleNotCalculated = TRUE;
          */
          bfp->bvd.moveToOldPos = TRUE;
          /*
          if (bfp->bvd.bsp != bsp) {
            BioseqLock (bsp);
            BioseqUnlock (bfp->bvd.bsp);
          }
          */
          bfp->bvd.bsp = bsp;
        }
      }
    }
    PointerToForm (bfp->form, (Pointer) bfp->bvd.bsp);
    AdjustDynamicGraphicViewer (&(bfp->bvd));
  }
  if (ommsp->message == OM_MSG_SELECT || ommsp->message == OM_MSG_DESELECT) {
    ResetClip ();
    sel = ObjMgrGetSelected ();
    if (sel != NULL && sel->next == NULL &&
        (sel->entityID == bfp->input_entityID || InBioseqViewEntityList (sel->entityID, &(bfp->bvd)))) {
      GatherItem (sel->entityID, sel->itemID, sel->itemtype, (Pointer) bfp, SetClickmeTitle);
    } else {
      /* SafeSetTitle (bfp->bvd.clickMe, "Double click on an item to launch the appropriate editor."); */
      SafeSetTitle (bfp->bvd.clickMe, "");
    }
    if (ommsp->message == OM_MSG_SELECT) {
      Update ();
    }
  }
  return OM_MSG_RET_OK;
}

static void CleanSmartViewer (BioseqViewFormPtr bfp)

{
  BioseqPagePtr  bpp;
  Uint2          userkey;
  ValNodePtr     vnp;

  if (bfp == NULL) return;
  bpp = bfp->currentBioseqPage;
  if (bpp != NULL && bpp->show != NULL) {
    bpp->show (&(bfp->bvd), FALSE);
    Update ();
  }
  if (bfp->input_entityID > 0) {
    if (bfp->userkey > 0) {
      userkey = bfp->userkey;
      bfp->userkey = 0;
      ObjMgrFreeUserData (bfp->input_entityID, bfp->procid, bfp->proctype, userkey);
      /* this may trigger another remove, hence bfp->userkey first set to 0 */
      for (vnp = bfp->bvd.entityList; vnp != NULL; vnp = vnp->next) {
        if (bfp->input_entityID != (Uint2) vnp->data.intvalue) {
          ObjMgrFreeUserData ((Uint2) vnp->data.intvalue, bfp->procid, bfp->proctype, userkey);
        }
      }
    }
    bfp->bvd.pict = DeletePicture (bfp->bvd.pict);
    if (bfp->bvd.slp_list != NULL) {
      bfp->bvd.slp_list = free_slp_list (bfp->bvd.slp_list);
    }
    if (bfp->bvd.anp_node != NULL) {
      bfp->bvd.anp_node = FreeAlignNode (bfp->bvd.anp_node);
    }
    if (bfp->bvd.g_list != NULL) {
      bfp->bvd.g_list = ValNodeFreeData (bfp->bvd.g_list);
    }
    if (bfp->bvd.ftype_list != NULL) {
      bfp->bvd.ftype_list = ValNodeFree (bfp->bvd.ftype_list);
    }
    if (bfp->bvd.sentinelList != NULL) {
      bfp->bvd.sentinelList = ValNodeFreeData (bfp->bvd.sentinelList);
    }
    if (bfp->bvd.entityList != NULL) {
      bfp->bvd.entityList = ValNodeFree (bfp->bvd.entityList);
    }
  }
  if (bfp->toolForm != NULL) {
    /* Hide (bfp->toolForm); */
    bfp->toolForm = Remove (bfp->toolForm);
  }
}

static void ToolFormHideWindowProc (WindoW w)

{
  Hide (w);
}

extern ForM MakeToolFormForBioseqView (BaseFormPtr bafp, GrpActnProc createToolBar)

{
  BioseqViewFormPtr  bfp;
  GrouP              g;
  Char               str [256];
  WindoW             w;

  bfp = (BioseqViewFormPtr) bafp;
  if (bfp == NULL || createToolBar == NULL) return NULL;
  if (bfp->toolForm != NULL) return bfp->toolForm;
  GetTitle (bfp->form, str, sizeof (str));
  TrimSpacesAroundString (str);
  if (StringHasNoText (str)) {
    StringCpy (str, "ToolBar");
  }
  w = FixedWindow (-5, -50, -10, -10, str, ToolFormHideWindowProc);
  if (w == NULL) return NULL;
  g = HiddenGroup (w, -1, 0, NULL);
  SetObjectExtra (g, bfp, NULL);
  createToolBar (g);
  RealizeWindow (w);
  bfp->toolForm = (ForM) w;
  return (ForM) w;
}

extern ForM RemoveSeqEntryViewer (ForM f)

{
  BioseqViewFormPtr  bfp;
  SeqViewProcsPtr    svpp;

  bfp = (BioseqViewFormPtr) GetObjectExtra (f);
  if (f == smartBioseqViewForm) {
    svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
    if (svpp != NULL && (! svpp->keepSmartViewerVisible)) {
      Hide (f);
      if (bfp != NULL) {
        bfp->toolForm = Remove (bfp->toolForm);
      }
    }
    CleanSmartViewer (bfp);
  } else {
    if (bfp != NULL) {
      bfp->toolForm = Remove (bfp->toolForm);
    }
    Remove (f);
  }
  return NULL;
}

extern Int2 LIBCALLBACK NewSeqEntryViewGenFunc (Pointer data)

{
  BioseqContextPtr   bcp;
  BioseqViewFormPtr  bfp;
  BioseqPtr          bsp;
  ForM               f;
  OMProcControlPtr   ompcp;
  OMUserDataPtr      omudp;
  ValNodePtr         sdp;
  SeqEntryPtr        sep;
  SeqIdPtr           sip;
  Char               str [64];
  SeqViewProcsPtr    svpp;
  ValNodePtr         ttl;
  WindoW             w;

  ompcp = (OMProcControlPtr) data;
  bsp = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    case OBJ_BIOSEQSET :
      return OM_MSG_RET_ERROR;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  sip = SeqIdFindWorst (bsp->id);
  if (sip == NULL) return OM_MSG_RET_ERROR;
  SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (IsAGenomeRecord (sep)) {
    bcp = BioseqContextNew (bsp);
    ttl = NULL;
    sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_title, NULL, NULL);
    while (sdp != NULL) {
      ttl = sdp;
      sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_title, sdp, NULL);
    }
    BioseqContextFree (bcp);
    if (ttl != NULL && (! StringHasNoText ((CharPtr) ttl->data.ptrvalue))) {
      StringNCpy_0 (str, (CharPtr) ttl->data.ptrvalue, sizeof (str));
    }
  }
  w = (WindoW) CreateNewSeqEntryViewFormEx (-50, -33, str, bsp, NULL, ompcp->input_entityID,
                                            ompcp->input_itemID, ompcp->input_itemtype, FALSE);
  bfp = (BioseqViewFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    bfp->input_entityID = ompcp->input_entityID;
    bfp->input_itemID = ompcp->input_itemID;
    bfp->input_itemtype = ompcp->input_itemtype;
    bfp->this_itemtype = OBJ_BIOSEQ;
    bfp->this_subtype = bsp->repr;
    bfp->procid = ompcp->proc->procid;
    bfp->proctype = ompcp->proc->proctype;
    bfp->docuid = GetUidFromBsp (bsp);
    if (ISA_na (bsp->mol)) {
      bfp->doctype = TYP_NT;
    } else if (ISA_aa (bsp->mol)) {
      bfp->doctype = TYP_AA;
    }
    bfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_VIEW, bfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) bfp;
      omudp->messagefunc = BioseqViewMsgFunc;
    }
    SendMessageToForm (bfp->form, VIB_MSG_CHANGE);
  }
  Show (w);
  Select (w);
  if (bfp != NULL) {
    svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
    if (svpp != NULL && svpp->createToolBar != NULL) {
      f = MakeToolFormForBioseqView ((BaseFormPtr) bfp, svpp->createToolBar);
      Show (f);
    }
  }
  return OM_MSG_RET_DONE;
}

extern Int2 LIBCALLBACK SmartSeqEntryViewGenFunc (Pointer data)

{
  BioseqContextPtr   bcp;
  BioseqViewFormPtr  bfp;
  BioseqPtr          bsp;
  ForM               f;
/* ObjMgrDataPtr      omdp; */
  OMProcControlPtr   ompcp;
  OMUserDataPtr      omudp;
  Boolean            reusing;
  ValNodePtr         sdp;
  SeqEntryPtr        sep;
  SeqIdPtr           sip;
  Char               str [64];
  SeqViewProcsPtr    svpp;
  ValNodePtr         ttl;
  Int2               val;
  WindoW             w;

  ompcp = (OMProcControlPtr) data;
  bsp = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    case OBJ_BIOSEQSET :
      return OM_MSG_RET_ERROR;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  /* if (bsp == NULL) return OM_MSG_RET_ERROR; */
  StringCpy (str, "no record");
  if (bsp != NULL) {
    sip = SeqIdFindWorst (bsp->id);
    if (sip == NULL) return OM_MSG_RET_ERROR;
    SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
    sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
    if (IsAGenomeRecord (sep)) {
      bcp = BioseqContextNew (bsp);
      ttl = NULL;
      sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_title, NULL, NULL);
      while (sdp != NULL) {
        ttl = sdp;
        sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_title, sdp, NULL);
      }
      BioseqContextFree (bcp);
      if (ttl != NULL && (! StringHasNoText ((CharPtr) ttl->data.ptrvalue))) {
        StringNCpy_0 (str, (CharPtr) ttl->data.ptrvalue, sizeof (str));
      }
    }
  }
  w = NULL;
  reusing = FALSE;
  if (smartBioseqViewForm != NULL) {
    svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
    if (svpp != NULL && svpp->forceSeparateViewer) {
      w = (WindoW) CreateNewSeqEntryViewFormEx (-50, -33, str, bsp, NULL, ompcp->input_entityID,
                                                ompcp->input_itemID, ompcp->input_itemtype, TRUE);
      if (smartBioseqViewForm == NULL) {
        smartBioseqViewForm = (ForM) w;
      }
    } else {
      /*
      bfp = (BioseqViewFormPtr) GetObjectExtra (smartBioseqViewForm);

      if (bfp != NULL) {
          if (bfp->input_entityID > 0) {
            omdp = ObjMgrGetData (bfp->input_entityID);
            if (omdp != NULL && omdp->dirty) {
                SendMessageToForm (smartBioseqViewForm, VIB_MSG_ACCEPT);
                ObjMgrSetDirtyFlag (bfp->input_entityID, FALSE);
            } else if(Visible((WindoW)smartBioseqViewForm) ||
                      Nlm_IconicWindow((WindoW)smartBioseqViewForm)) {
                SendMessageToForm (smartBioseqViewForm, VIB_MSG_RESET);
            }
          }
      }

      CleanSmartViewer (bfp);
      Update ();
      */
      w = (WindoW) smartBioseqViewForm;
      SetTitle (w, str);
      reusing = TRUE;
    }
  } else {
    w = (WindoW) CreateNewSeqEntryViewFormEx (-50, -33, str, bsp, NULL, ompcp->input_entityID,
                                              ompcp->input_itemID, ompcp->input_itemtype, TRUE);
    if (smartBioseqViewForm == NULL) {
      smartBioseqViewForm = (ForM) w;
    }
  }
  bfp = (BioseqViewFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    bfp->bvd.bsp = bsp;
    SafeHide (bfp->targetControl);
    Update ();
    Reset (bfp->targetControl);
    val = PopulateTarget (bfp);
    SetValue (bfp->targetControl, val + 1);
    SafeShow (bfp->targetControl);
    bfp->input_entityID = ompcp->input_entityID;
    bfp->input_itemID = ompcp->input_itemID;
    bfp->input_itemtype = ompcp->input_itemtype;
    bfp->this_itemtype = OBJ_BIOSEQ;
    bfp->procid = ompcp->proc->procid;
    bfp->proctype = ompcp->proc->proctype;
    if (bsp != NULL) {
      bfp->this_subtype = bsp->repr;
      bfp->docuid = GetUidFromBsp (bsp);
    } else {
      bfp->this_subtype = Seq_repr_raw;
      bfp->docuid = 0;
    }
    if (bsp == NULL) {
    } else if (ISA_na (bsp->mol)) {
      bfp->doctype = TYP_NT;
    } else if (ISA_aa (bsp->mol)) {
      bfp->doctype = TYP_AA;
    }
    bfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_VIEW, bfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) bfp;
      omudp->messagefunc = BioseqViewMsgFunc;
    }
    if (reusing) {
      ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    }
    SendMessageToForm (bfp->form, VIB_MSG_CHANGE);
  }
  if (bsp != NULL) {
    Show (w);
    Select (w);
  }
  if (bfp != NULL) {
    svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
    if (svpp != NULL && svpp->createToolBar != NULL) {
      f = MakeToolFormForBioseqView ((BaseFormPtr) bfp, svpp->createToolBar);
      Show (f);
    }
  }
  return OM_MSG_RET_DONE;
}

