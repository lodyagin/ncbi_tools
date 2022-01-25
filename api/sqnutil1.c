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
* $Revision: 6.913 $
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
#include <subutil.h>
#include <asn2gnbi.h>
#include <salpacc.h>
#include <alignmgr2.h>
#include <valid.h>
#include <objvalid.h>
#include <valapi.h>
#include <findrepl.h>


#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

static int descr_insert_order [] = {
  Seq_descr_title,
  Seq_descr_source,
  Seq_descr_molinfo,
  Seq_descr_het,
  Seq_descr_pub,
  Seq_descr_comment,
  Seq_descr_name,
  Seq_descr_user,
  Seq_descr_maploc,
  Seq_descr_region,
  Seq_descr_num,
  Seq_descr_dbxref,
  Seq_descr_mol_type,
  Seq_descr_modif,
  Seq_descr_method,
  Seq_descr_org,
  Seq_descr_sp,
  Seq_descr_pir,
  Seq_descr_prf,
  Seq_descr_pdb,
  Seq_descr_embl,
  Seq_descr_genbank,
  Seq_descr_modelev,
  Seq_descr_create_date,
  Seq_descr_update_date,
  0
};

static void NormalizeDescriptorProc (
  SeqEntryPtr sep,
  Pointer data,
  Int4 index,
  Int2 indent
)

{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  /* arrays are SEQDESCR_MAX + 1, last slot stores unexpected descriptor numbers */
  SeqDescrPtr       first [SEQDESCR_MAX + 1];
  SeqDescrPtr       last [SEQDESCR_MAX + 1];
  int               i;
  int               idx;
  SeqDescrPtr PNTR  head = NULL;
  SeqDescrPtr PNTR  prev = NULL;
  SeqDescrPtr       next;
  SeqDescrPtr       sdp;

  if (sep == NULL) return;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    head = &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    head = &(bssp->descr);
  }
  if (head == NULL) return;

  MemSet ((Pointer) &first, 0, sizeof (first));
  MemSet ((Pointer) &last, 0, sizeof (last));

  prev = head;
  sdp = *prev;
  while (sdp != NULL) {
    next = sdp->next;

    *prev = sdp->next;
    sdp->next = NULL;

    idx = (int) sdp->choice;
    /* unexpected descriptor numbers go into last slot */
    if (idx <= 0 || idx >= SEQDESCR_MAX) {
      idx = SEQDESCR_MAX;
    }
    if (idx > 0 && idx <= SEQDESCR_MAX) {
      if (first [idx] == NULL) {
        first [idx] = sdp;
      }
      if (last [idx] != NULL) {
        (last [idx])->next = sdp;
      }
      last [idx] = sdp;
    }

    sdp = next;
  }

  for (i = 0; descr_insert_order [i] != 0; i++) {
    idx = descr_insert_order [i];
    sdp = first [idx];
    if (sdp == NULL) continue;
    ValNodeLink (head, sdp);
  }
}

NLM_EXTERN void NormalizeDescriptorOrder (
  SeqEntryPtr sep
)

{
  SeqEntryExplore (sep, NULL, NormalizeDescriptorProc);
}

typedef struct orgscan {
  ObjMgrPtr     omp;
  Int2          nuclCode;
  Int2          mitoCode;
  Int2          pstdCode;
  Boolean       mito;
  Boolean       plastid;
  Char          taxname [196];
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
  Int2           pstdCode = 0;
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
          //LCOV_EXCL_START
          //org features are converted to biosrc features in BasicCleanup
          orp = (OrgRefPtr) sfp->data.value.ptrvalue;
          break;
          //LCOV_EXCL_STOP
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
          //LCOV_EXCL_START
          // org descriptors are converted to biosrc descriptors in basiccleanup
          orp = (OrgRefPtr) sdp->data.ptrvalue;
          break;
          //LCOV_EXCL_STOP
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
    mito = (Boolean) (biop->genome == GENOME_kinetoplast ||
                      biop->genome == GENOME_mitochondrion ||
                      biop->genome == GENOME_hydrogenosome);
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
      pstdCode = onp->pgcode;
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
      osp->plastid = (Boolean) (biop->genome == GENOME_chloroplast ||
                                biop->genome == GENOME_chromoplast ||
                                biop->genome == GENOME_plastid ||
                                biop->genome == GENOME_cyanelle ||
                                biop->genome == GENOME_apicoplast ||
                                biop->genome == GENOME_leucoplast ||
                                biop->genome == GENOME_proplastid ||
                                biop->genome == GENOME_chromatophore);
      if (doCodes) {
        osp->nuclCode = nuclCode;
        osp->mitoCode = mitoCode;
        osp->pstdCode = pstdCode;
      }
      if (doTaxname) {
        StringNCpy_0 (osp->taxname, taxname, sizeof (osp->taxname));
      }
    }
  }

  return TRUE;
}

//LCOV_EXCL_START
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
  osp.plastid = FALSE;
  osp.nuclCode = 0;
  osp.mitoCode = 0;
  osp.pstdCode = 0;
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
  if (osp.plastid) {
    if (osp.pstdCode > 0) {
      return osp.pstdCode;
    } else {
      return 11;
    }
  } else if (osp.mito) {
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


NLM_EXTERN Boolean BioseqToGeneticCode (
  BioseqPtr bsp,
  Int2Ptr gencodep,
  BoolPtr mitop,
  BoolPtr plastidp,
  CharPtr taxnamep,
  size_t maxsize,
  BioSourcePtr PNTR biopp
)

{
  BioSourcePtr       biop = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Int2               gencode = 0;
  Boolean            mito = FALSE;
  Int2               mitoCode = 0;
  Int2               nuclCode = 0;
  Int2               pstdCode = 0;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  Boolean            plastid = FALSE;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  CharPtr            taxname = NULL;

  if (bsp == NULL) return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }

  if (biop == NULL) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }

  if (biop == NULL) return FALSE;
  orp = biop->org;
  if (orp == NULL) return FALSE;

  taxname = orp->taxname;
  if (StringHasNoText (taxname)) return FALSE;

  onp = orp->orgname;
  if (onp != NULL) {
    nuclCode = onp->gcode;
    mitoCode = onp->mgcode;
    pstdCode = onp->pgcode;
  }

  mito = (Boolean) (biop->genome == GENOME_kinetoplast ||
                    biop->genome == GENOME_mitochondrion ||
                    biop->genome == GENOME_hydrogenosome);

  plastid = (Boolean) (biop->genome == GENOME_chloroplast ||
                       biop->genome == GENOME_chromoplast ||
                       biop->genome == GENOME_plastid ||
                       biop->genome == GENOME_cyanelle ||
                       biop->genome == GENOME_apicoplast ||
                       biop->genome == GENOME_leucoplast ||
                       biop->genome == GENOME_proplastid ||
                       biop->genome == GENOME_chromatophore);

  if (plastid) {
    if (pstdCode > 0) {
      gencode = pstdCode;
    } else {
      gencode = 11;
    }
  } else if (mito) {
    gencode = mitoCode;
  } else {
    gencode = nuclCode;
  }

  if (gencodep != NULL) {
    *gencodep = gencode;
  }
  if (mitop != NULL) {
    *mitop = mito;
  }
  if (plastidp != NULL) {
    *plastidp = plastid;
  }
  if (taxnamep != NULL && maxsize > 0) {
    StringNCpy_0 (taxnamep, taxname, maxsize);
  }
  if (biopp != NULL) {
    *biopp = biop;
  }

  return TRUE;
}


typedef struct commontitle {
  BioseqPtr bsp;
  SeqDescPtr sdp;
} CommonTitleData, PNTR CommonTitlePtr;


static CommonTitlePtr CommonTitleNew (BioseqPtr bsp, SeqDescPtr sdp)
{
  CommonTitlePtr c = (CommonTitlePtr) MemNew (sizeof (CommonTitleData));
  c->bsp = bsp;
  c->sdp = sdp;
  return c;
}


static CommonTitlePtr CommonTitleFree (CommonTitlePtr c)
{
  if (c != NULL) {
    c = MemFree (c);
  }
  return c;
}


static ValNodePtr CommonTitleListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = CommonTitleFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static void RemoveCommonTitles (ValNodePtr vnp, CharPtr common_title)
{
  CommonTitlePtr c;
  ObjValNodePtr  ovp;

  while (vnp != NULL) {
    c = vnp->data.ptrvalue;
    if (StringCmp (c->sdp->data.ptrvalue, common_title) == 0 && c->sdp->extended > 0) {
      ovp = (ObjValNodePtr) c->sdp;
      ovp->idx.deleteme = TRUE;
    }
    vnp = vnp->next;
  }
}


static int LIBCALLBACK SortCommonTitle (VoidPtr ptr1, VoidPtr ptr2)

{
  CommonTitlePtr c1;
  CommonTitlePtr c2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      c1 = (CommonTitlePtr) vnp1->data.ptrvalue;
      c2 = (CommonTitlePtr) vnp2->data.ptrvalue;
      if (c1 != NULL && c2 != NULL && c1->sdp != NULL && c2->sdp != NULL
          && c1->sdp->data.ptrvalue != NULL && c2->sdp->data.ptrvalue != NULL) {
        return StringCmp (c1->sdp->data.ptrvalue, c2->sdp->data.ptrvalue);
      }
    }
  }
  return 0;
}


static void CollectCommonTitle (BioseqPtr bsp, Pointer data)
{
  SeqDescPtr sdp;

  if (bsp == NULL || ISA_aa (bsp->mol) || data == NULL) {
    return;
  }

  sdp = bsp->descr;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_title) {
      ValNodeAddPointer ((ValNodePtr PNTR) data, 0, CommonTitleNew (bsp, sdp));
    }
    sdp = sdp->next;
  }
}


static CharPtr FindCommonTitleFromList (ValNodePtr list)
{
  ValNodePtr vnp;
  CommonTitlePtr c;
  Int4 num_common = 0, num_total, num_expected;
  CharPtr common_title;

  if (list == NULL) {
    return NULL;
  }
  num_total = ValNodeLen (list);
  if (num_total % 2 != 0 || num_total < 4) {
    return NULL;
  }
  num_expected = num_total / 2;

  c = list->data.ptrvalue;
  common_title = c->sdp->data.ptrvalue;
  num_common = 1;

  for (vnp = list->next; vnp != NULL; vnp = vnp->next) {
    c = (CommonTitlePtr) vnp->data.ptrvalue;
    if (StringCmp (common_title, c->sdp->data.ptrvalue) == 0) {
      num_common++;
    } else {
      num_common = 1;
      common_title = c->sdp->data.ptrvalue;
    }
  }
  if (num_common == num_expected) {
    return StringSave (common_title);
  }

  return NULL;
}


static void PromoteCommonTitlesSetCallback (BioseqSetPtr bssp, Pointer data)
{
  ValNodePtr list = NULL;
  CharPtr common_title = NULL;
  SeqDescrPtr sdp;
  Int4        num_member = 0;
  SeqEntryPtr s;
  CharPtr     set_title = NULL;

  if (bssp == NULL || !GetsDocsumTitle (bssp->_class)) {
    return;
  }

  VisitBioseqsInSet (bssp, &list, CollectCommonTitle);
  list = ValNodeSort (list, SortCommonTitle);

  common_title = FindCommonTitleFromList (list);
  if (common_title != NULL) {
    s = bssp->seq_set;
    while (s != NULL) {
      num_member++;
      s = s->next;
    }
    if (ValNodeLen (list) == num_member) {
      for (sdp = bssp->descr; sdp != NULL && set_title == NULL; sdp = sdp->next) {
        if (sdp->choice == Seq_descr_title) {
          set_title = sdp->data.ptrvalue;
        }
      }
      if (set_title != NULL
          && StringCmp (set_title, common_title) != 0) {
        /* don't remove, the seq titles just happen to be identical */
        common_title = MemFree (common_title);
      }
    }
  }
  if (common_title != NULL) {
    sdp = SeqDescrNew (NULL);
    sdp->choice = Seq_descr_title;
    sdp->data.ptrvalue = common_title;
    sdp->next = bssp->descr;
    bssp->descr = sdp;
    RemoveCommonTitles (list, common_title);
  }
  list = CommonTitleListFree(list);
}


NLM_EXTERN void PromoteCommonTitlesToSet (SeqEntryPtr sep)
{
  VisitSetsInSep (sep, NULL, PromoteCommonTitlesSetCallback);
}
//LCOV_EXCL_STOP

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
      //LCOV_EXCL_START
      //cleanup functions only call this during RenormalizeNucProtSets,
      //and only for Bioseqs
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
    //LCOV_EXCL_STOP
  } else return;
  hastitle = FALSE;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_title) {
      if (hastitle) {
        //LCOV_EXCL_START
        //when called from RenormalizeNucProtSets, 
        //extra titles are already gone
        *(prevsdp) = sdp->next;
        sdp->next = NULL;
        SeqDescFree (sdp);
        //LCOV_EXCL_STOP
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

NLM_EXTERN Int4 RenormalizeNucProtSets (SeqEntryPtr sep, Boolean relink)

{
  SeqAnnotPtr    annot;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     descr;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqAnnotPtr    sap, tmp_sap;
  SeqEntryPtr    seqentry;
  Int4           num_renormalized = 0;

  if (sep == NULL) return 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (bssp->_class >= 13 && bssp->_class <= 16) ||
                         bssp->_class == BioseqseqSet_class_wgs_set ||
                         bssp->_class == BioseqseqSet_class_gen_prod_set ||
                         bssp->_class == BioseqseqSet_class_small_genome_set)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        num_renormalized += RenormalizeNucProtSets (sep, relink);
      }
      return num_renormalized;
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
        bssp->seqentry = NULL;
        MemFree (seqentry);
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
            //LCOV_EXCL_START
            //should not have set inside nuc-prot set
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          ValNodeLink (&(bssp->descr), descr);
          if (bssp->annot == NULL) {
            bssp->annot = annot;
            annot = NULL;
          } else {
            sap = bssp->annot;
          }
          //LCOV_EXCL_STOP
        }
        if (sap != NULL) {
          tmp_sap = sap;
          while (tmp_sap->next != NULL) {
            tmp_sap = tmp_sap->next;
          }
          tmp_sap->next = annot;
          MergeAdjacentAnnotsInList (sap);
        }

        DeleteMultipleTitles (sep, NULL, 0, 0);

        if (relink) {
          SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
          RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
        }
        num_renormalized++;
      }
    }
  }
  return num_renormalized;
}


//LCOV_EXCL_START
//only used by RemoveSingleItemSet, which is not used by cleanup
static Boolean SetHasAlignments (BioseqSetPtr bssp)
{
  SeqAnnotPtr sap;
  Boolean     rval = FALSE;

  if (bssp == NULL) {
    return FALSE;
  }
  for (sap = bssp->annot; sap != NULL && !rval; sap = sap->next) {
    if (sap->type == 2) {
      rval = TRUE;
    }
  }
  return rval;
}


//not used by cleanup
NLM_EXTERN Int4 RemoveSingleItemSet(SeqEntryPtr sep, Boolean relink)
{
  SeqAnnotPtr    annot;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     descr;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqAnnotPtr    sap, tmp_sap;
  SeqEntryPtr    seqentry, sep_next;
  Int4           num_renormalized = 0;

  if (sep == NULL 
      || !IS_Bioseq_set (sep) 
      || (bssp = (BioseqSetPtr) sep->data.ptrvalue) == NULL) {
    return 0;
  }

  if ((bssp->_class == BioseqseqSet_class_pop_set
       || bssp->_class == BioseqseqSet_class_phy_set
       || bssp->_class == BioseqseqSet_class_mut_set
       || bssp->_class == BioseqseqSet_class_eco_set)
      && bssp->seq_set != NULL
      && bssp->seq_set->next == NULL
      && !SetHasAlignments(bssp)) {

    seqentry = bssp->seq_set;

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
    bssp->seqentry = NULL;
    MemFree (seqentry);
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
      tmp_sap = sap;
      while (tmp_sap->next != NULL) {
        tmp_sap = tmp_sap->next;
      }
      tmp_sap->next = annot;
      MergeAdjacentAnnotsInList (sap);
    }

    DeleteMultipleTitles (sep, NULL, 0, 0);

    if (relink) {
      SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
      RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
    }
    num_renormalized++;
  } else {
    for (sep = bssp->seq_set; sep != NULL; sep = sep_next) {
      sep_next = sep->next;
      num_renormalized += RemoveSingleItemSet (sep, relink);
    }
  }

  return num_renormalized;
}


static Boolean IsExtractableDescriptor (SeqDescPtr sdp)
{
  if (sdp == NULL) {
    return FALSE;
  }
  if (sdp->choice == Seq_descr_pub || sdp->choice == Seq_descr_source) {
    return TRUE;
  } else if (sdp->choice == Seq_descr_user && IsDBLinkObject(sdp->data.ptrvalue)) {
    return TRUE;
  } else {
    return FALSE;
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
    if (IsExtractableDescriptor(sdp)) {
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
          if (bssp->_class != BioseqseqSet_class_genbank &&
              bssp->_class != BioseqseqSet_class_mut_set &&
              bssp->_class != BioseqseqSet_class_pop_set &&
              bssp->_class != BioseqseqSet_class_phy_set &&
              bssp->_class != BioseqseqSet_class_eco_set &&
              bssp->_class != BioseqseqSet_class_wgs_set &&
              bssp->_class != BioseqseqSet_class_small_genome_set &&
              (bssp->_class != BioseqseqSet_class_gen_prod_set ||
           (! tdp->skipGenProdSet))) {
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
  BioseqSetPtr  bssp;
  BioseqSetPtr  parent;
  GatherScope   gs;
  TargetData    td;

  td.bsp = bsp;
  td.nps = NULL;
  td.skipGenProdSet = skipGenProdSet;
  if (entityID > 0 && bsp != NULL) {
    if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bsp->idx.parentptr;
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts && bssp->idx.parenttype == OBJ_BIOSEQSET) {
        parent = (BioseqSetPtr) bssp->idx.parentptr;
        if (parent != NULL && parent->_class == BioseqseqSet_class_segset) {
          bssp = parent;
        }
      }
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset && bssp->idx.parenttype == OBJ_BIOSEQSET) {
        parent = (BioseqSetPtr) bssp->idx.parentptr;
        if (parent != NULL && parent->_class == BioseqseqSet_class_nuc_prot) {
          bssp = parent;
        }
      }
      if (bssp != NULL && bssp->seqentry != NULL) {
        if (bssp->_class == BioseqseqSet_class_nuc_prot ||
            bssp->_class == BioseqseqSet_class_segset ||
            bssp->_class == BioseqseqSet_class_parts) {
          return bssp->seqentry;
        }
      }
      if (bsp->seqentry != NULL) {
        return bsp->seqentry;
      }
    }
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

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForItemIDEx (Uint2 entityID, Uint4 itemID, Uint2 itemtype, Boolean skipGenProdSet)

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

NLM_EXTERN SeqEntryPtr LIBCALL GetBestTopParentForItemID (Uint2 entityID, Uint4 itemID, Uint2 itemtype)

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
//LCOV_EXCL_STOP

NLM_EXTERN Boolean CheckSeqLocForPartialEx (SeqLocPtr location, BoolPtr p5ptr, BoolPtr p3ptr, Int4Ptr limptr)

{
  SeqLocPtr   firstSlp;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  Int4        lim;
  Boolean     partial5;
  Boolean     partial3;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  SeqPntPtr   spp;

  partial5 = FALSE;
  partial3 = FALSE;
  lim = -1;
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
        ifp = spp->fuzz;
        if (ifp != NULL && ifp->choice == 4) {
          lim = ifp->a;
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
        ifp = spp->fuzz;
        if (ifp != NULL && ifp->choice == 4) {
          lim = ifp->a;
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
  if (limptr != NULL) {
    *limptr = lim;
  }
  return (Boolean) (partial5 || partial3 || lim == 3 || lim == 4);
}

NLM_EXTERN Boolean CheckSeqLocForPartial (SeqLocPtr location, BoolPtr p5ptr, BoolPtr p3ptr)

{
  return CheckSeqLocForPartialEx (location, p5ptr, p3ptr, NULL);
}

static void ConvertWholeToIntLoc (SeqLocPtr slp)
{
  BioseqPtr bsp;
  SeqIntPtr sip;
  
  if (slp == NULL || slp->choice != SEQLOC_WHOLE || slp->data.ptrvalue == NULL)
  {
    return;
  }
  bsp = BioseqFind (slp->data.ptrvalue);
  if (bsp == NULL)
  {
    return;
  }
  
  sip = SeqIntNew ();
  if (sip != NULL)
  {
    sip->from = 0;
    sip->to = bsp->length - 1;
    sip->id = SeqIdDup (SeqIdFindBest (bsp->id, 0));
    sip->strand = bsp->strand;
    slp->data.ptrvalue = SeqIdFree (slp->data.ptrvalue);
    slp->data.ptrvalue = sip;
    slp->choice = SEQLOC_INT;
  } 
}

NLM_EXTERN void SetSeqLocPartialEx (SeqLocPtr location, Boolean partial5, Boolean partial3, Int4 lim)

{
  SeqLocPtr   firstSlp;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  SeqPntPtr   spp;

  if (location != NULL) {
    /* if whole, need to convert to int */
    if (partial5 || partial3)
    {
      ConvertWholeToIntLoc (location);
    }
  
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
              sip->if_to = IntFuzzFree (sip->if_to);
              sip->if_to = ifp;
              ifp->a = 1;
            } else {
              sip->if_from = IntFuzzFree (sip->if_from);
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
              spp->fuzz = IntFuzzFree (spp->fuzz);
              spp->fuzz = ifp;
              ifp->a = 1;
            } else {
              spp->fuzz = IntFuzzFree (spp->fuzz);
              spp->fuzz = ifp;
              ifp->a = 2;
            }
          }
        } else if (lim == 3 || lim == 4) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            spp->fuzz = IntFuzzFree (spp->fuzz);
            spp->fuzz = ifp;
            ifp->a = lim;
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
              sip->if_from = IntFuzzFree (sip->if_from);
              sip->if_from = ifp;
              ifp->a = 2;
            } else {
              sip->if_to = IntFuzzFree (sip->if_to);
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
              spp->fuzz = IntFuzzFree (spp->fuzz);
              spp->fuzz = ifp;
              ifp->a = 2;
            } else {
              spp->fuzz = IntFuzzFree (spp->fuzz);
              spp->fuzz = ifp;
              ifp->a = 1;
            }
          }
        } else if (lim == 3 || lim == 4) {
          ifp = IntFuzzNew ();
          if (ifp != NULL) {
            ifp->choice = 4;
            spp->fuzz = IntFuzzFree (spp->fuzz);
            spp->fuzz = ifp;
            ifp->a = lim;
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

NLM_EXTERN void SetSeqLocPartial (SeqLocPtr location, Boolean partial5, Boolean partial3)

{
  SetSeqLocPartialEx (location, partial5, partial3, -1);
}

//LCOV_EXCL_START
NLM_EXTERN ValNodePtr GetSeqLocPartialSet (SeqLocPtr location)

{
  ValNodePtr  head = NULL, last = NULL, vnp;
  Int4        lim;
  Boolean     noLeft;
  Boolean     noRight;
  SeqLocPtr   slp;
  Int4        val;

  if (location == NULL) return NULL;

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    CheckSeqLocForPartialEx (slp, &noLeft, &noRight, &lim);
    val = 0;
    if (noLeft) {
      val |= 2;
    }
    if (noRight) {
      val |= 1;
    }
    if (lim == 3) {
      val |= 4;
    } else if (lim == 4) {
      val |= 8;
    }
    vnp = ValNodeAddInt (&last, 0, val);
    if (head == NULL) {
      head = vnp;
    }
    last = vnp;
    slp = SeqLocFindNext (location, slp);
  }

  return head;
}

NLM_EXTERN void SetSeqLocPartialSet (SeqLocPtr location, ValNodePtr vnp)

{
  Int4        lim;
  Boolean     noLeft;
  Boolean     noRight;
  SeqLocPtr   slp;
  Int4        val;

  if (location == NULL || vnp == NULL) return;

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL && vnp != NULL) {
    val = (Int4) vnp->data.intvalue;
    noLeft = (Boolean) ((val & 2) != 0);
    noRight = (Boolean) ((val & 1) != 0);
    lim = -1;
    if ((val & 4) != 0) {
      lim = 3;
    } else if ((val & 8) != 0) {
      lim = 4;
    }
    SetSeqLocPartialEx (slp, noLeft, noRight, lim);
    slp = SeqLocFindNext (location, slp);
    vnp = vnp->next;
  }
}

/* KeyTag section */

NLM_EXTERN int LIBCALLBACK SortVnpByString (VoidPtr ptr1, VoidPtr ptr2)

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
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByStringCS (VoidPtr ptr1, VoidPtr ptr2)

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
        return StringCmp (str1, str2);
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByStringCI (VoidPtr ptr1, VoidPtr ptr2)

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
        return StringCmp (str1, str2);
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByStringCIUCFirst (VoidPtr ptr1, VoidPtr ptr2)

{
  int         comp;
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
        comp = StringICmp (str1, str2);
        if (comp != 0) return comp;
        return StringCmp (str1, str2);
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByStringCILCFirst (VoidPtr ptr1, VoidPtr ptr2)

{
  int         comp;
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
        comp = StringICmp (str1, str2);
        if (comp != 0) return comp;
        return StringCmp (str2, str1);
      }
    }
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortVnpByNaturalCS (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  str1 = (CharPtr) vnp1->data.ptrvalue;
  str2 = (CharPtr) vnp2->data.ptrvalue;
  if (str1 == NULL || str2 == NULL) return 0;

  return NaturalStringCmp (str1, str2);
}

NLM_EXTERN int LIBCALLBACK SortVnpByNaturalCI (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  str1 = (CharPtr) vnp1->data.ptrvalue;
  str2 = (CharPtr) vnp2->data.ptrvalue;
  if (str1 == NULL || str2 == NULL) return 0;

  return NaturalStringICmp (str1, str2);
}
//LCOV_EXCL_STOP

NLM_EXTERN ValNodePtr UniqueValNode (ValNodePtr list)

{
  CharPtr       last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       str;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringICmp (last, str) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFreeData (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

//LCOV_EXCL_START
NLM_EXTERN ValNodePtr UniqueStringValNodeCS (ValNodePtr list)

{
  CharPtr       last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       str;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringCmp (last, str) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFreeData (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN ValNodePtr UniqueStringValNodeCI (ValNodePtr list)

{
  CharPtr       last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  CharPtr       str;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (CharPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringICmp (last, str) == 0) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFreeData (vnp);
    } else {
      last = (CharPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN int LIBCALLBACK SortByChoice (VoidPtr ptr1, VoidPtr ptr2)

{
  Uint1       chs1;
  Uint1       chs2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  chs1 = (Uint1) vnp1->choice;
  chs2 = (Uint1) vnp2->choice;
  if (chs1 > chs2) {
    return 1;
  } else if (chs1 < chs2) {
    return -1;
  }
  return 0;
}

NLM_EXTERN int LIBCALLBACK SortByIntvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  Int4        val1;
  Int4        val2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  val1 = (Int4) vnp1->data.intvalue;
  val2 = (Int4) vnp2->data.intvalue;
  if (val1 > val2) {
    return 1;
  } else if (val1 < val2) {
    return -1;
  }
  return 0;
}

NLM_EXTERN ValNodePtr UniqueIntValNode (ValNodePtr list)

{
  Int4          curr, last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (Int4) list->data.intvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    curr = (Int4) vnp->data.intvalue;
    if (last == curr) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFree (vnp);
    } else {
      last = (Int4) vnp->data.intvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN int LIBCALLBACK SortByPtrvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  VoidPtr     val1;
  VoidPtr     val2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  val1 = (VoidPtr) vnp1->data.ptrvalue;
  val2 = (VoidPtr) vnp2->data.ptrvalue;
  if (val1 > val2) {
    return 1;
  } else if (val1 < val2) {
    return -1;
  }
  return 0;
}

NLM_EXTERN ValNodePtr UniquePtrValNode (ValNodePtr list)

{
  VoidPtr       curr, last;
  ValNodePtr    next;
  Pointer PNTR  prev;
  ValNodePtr    vnp;

  if (list == NULL) return NULL;
  last = (VoidPtr) list->data.ptrvalue;
  vnp = list->next;
  prev = (Pointer PNTR) &(list->next);
  while (vnp != NULL) {
    next = vnp->next;
    curr = (VoidPtr) vnp->data.ptrvalue;
    if (last == curr) {
      vnp->next = NULL;
      *prev = next;
      ValNodeFree (vnp);
    } else {
      last = (VoidPtr) vnp->data.ptrvalue;
      prev = (Pointer PNTR) &(vnp->next);
    }
    vnp = next;
  }

  return list;
}

NLM_EXTERN void KeyTagInit (KeyTag PNTR ktp, ValNodePtr list)

{
  Int2          i;
  CharPtr PNTR  index;
  Int2          num;
  ValNodePtr    vnp;

  if (ktp == NULL || list == NULL) return;
  list = ValNodeSort (list, SortVnpByString);
  list = UniqueValNode (list);
  num = ValNodeLen (list);
  index = MemNew (sizeof (CharPtr) * (num + 1));

  for (vnp = list, i = 0; vnp != NULL && i < num; vnp = vnp->next, i++) {
    index [i] = (CharPtr) vnp->data.ptrvalue;
  }

  ktp->num = num;
  ktp->list = list;
  ktp->index = index;
}

NLM_EXTERN void KeyTagClear (KeyTag PNTR ktp)

{
  if (ktp == NULL) return;
  ktp->num = 0;
  ktp->list = ValNodeFreeData (ktp->list);
  ktp->index = MemFree (ktp->index);
}

NLM_EXTERN Int2 KeyFromTag (KeyTag PNTR ktp, CharPtr tag)

{
  Int2  L, R, mid, compare;

  if (ktp == NULL || ktp->list == NULL || ktp->index == NULL) return 0;
  if (tag == NULL) return 0;

  L = 0;
  R = ktp->num - 1;
  while (L < R) {
    mid = (L + R) / 2;
    compare = StringICmp (ktp->index [mid], tag);
    if (compare < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }
  if (StringICmp (ktp->index [R], tag) == 0) {
    return (R + 1);
  }

  return 0;
}

NLM_EXTERN CharPtr TagFromKey (KeyTag PNTR ktp, Int2 key)

{
  if (ktp == NULL || ktp->list == NULL || ktp->index == NULL) return 0;
  if (key < 1 || key > ktp->num) return 0;
  key--;
  return ktp->index [key];
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

/*
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
*/

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

NLM_EXTERN void PromoteXrefsExEx (
  SeqFeatPtr sfp,
  BioseqPtr bsp,
  Uint2 entityID,
  Boolean include_stop,
  Boolean remove_trailingX,
  Boolean gen_prod_set,
  Boolean force_local_id,
  BoolPtr seq_fetch_failP
)

{
  Int2                 adv;
  ByteStorePtr         bs;
  BioseqSetPtr         bssp;
  Char                 ch;
  CharPtr              comment;
  CdRegionPtr          crp;
  Int2                 ctr = 1;
  ValNodePtr           descr;
  SeqFeatPtr           first;
  GBQualPtr            gbq;
  Int4                 i;
  Char                 id [128];
  SeqEntryPtr          last;
  Char                 lcl [128];
  BioseqPtr            mbsp;
  MolInfoPtr           mip;
  SeqEntryPtr          msep;
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
  ProtRefPtr           prp, prp2;
  SeqEntryPtr          psep;
  CharPtr              ptr;
  CharPtr              rnaseq;
  SeqEntryPtr          sep;
  SeqHistPtr           shp;
  SeqIdPtr             sip;
  SeqEntryPtr          target = NULL;
  Uint4                version = 0;
  long int             val;
  ValNodePtr           vnp;
  SeqFeatXrefPtr       xref;
  Boolean              ok_to_remove;
  /*
  DbtagPtr             dbt;
  SeqFeatPtr           gene;
  GeneRefPtr           grp;
  */

  if (seq_fetch_failP != NULL) {
    *seq_fetch_failP = FALSE;
  }

  if (sfp == NULL || bsp == NULL) return;

  /* set subtypes, used to find mRNA features for genomic product sets */

  first = sfp;
  while (sfp != NULL) {
    if (sfp->idx.subtype == 0) {
      sfp->idx.subtype = FindFeatDefType (sfp);
    }
    sfp = sfp->next;
  }

  /* no longer expand genes specified by qualifiers on other features (except repeat_region) */

  /*
  sfp = first;
  while (sfp != NULL) {
    prev = &(sfp->xref);
    xref = sfp->xref;
    while (xref != NULL) {
      next = xref->next;
      if (xref->data.choice == SEQFEAT_GENE &&
          sfp->data.choice != SEQFEAT_GENE &&
          sfp->idx.subtype != FEATDEF_repeat_region) {
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
                for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
                  dbt = (DbtagPtr) vnp->data.ptrvalue;
                  if (dbt == NULL) continue;
                  ValNodeAddPointer (&(gene->dbxref), 0, (Pointer) DbtagDup (dbt));
                }
              }
            }
          }
          *(prev) = next;
          xref->next = NULL;
          xref->data.choice = 0;
          SeqFeatXrefFree (xref);
        }
      } else {
        prev = &(xref->next);
      }
      xref = next;
    }
    sfp = sfp->next;
  }
  */

  /* expand mRNA features into cDNA product sequences */

  bssp = NULL;
  sep = NULL;
  last = NULL;
  if (gen_prod_set) {
    sep = GetTopSeqEntryForEntityID (entityID);
    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && bssp->seq_set != NULL) {
        last = bssp->seq_set;
        while (last->next != NULL) {
          last = last->next;
        }
      }
    }
  }

  if (gen_prod_set && sep != NULL && bssp != NULL && last != NULL) {
    target = sep;
    SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
    GetSeqEntryParent (target, &parentptr, &parenttype);
    sfp = first;
    while (sfp != NULL) {
      if (sfp->data.choice == SEQFEAT_RNA &&
          /* sfp->idx.subtype != FEATDEF_tRNA && */
          sfp->product == NULL && (! sfp->pseudo)) {
        gbq = sfp->qual;
        prevqual = (GBQualPtr PNTR) &(sfp->qual);
        id [0] = '\0';
        sip = NULL;
        comment = NULL;
        while (gbq != NULL) {
          nextqual = gbq->next;
          if (StringICmp (gbq->qual, "transcript_id") == 0) {
            if (StringDoesHaveText (id) && StringDoesHaveText (gbq->val)) {
              ErrPostEx (SEV_WARNING, ERR_FEATURE_QualWrongThisFeat,
                         "RNA transcript_id %s replacing %s", gbq->val, id);
            }
            *(prevqual) = gbq->next;
            gbq->next = NULL;
            StringNCpy_0 (id, gbq->val, sizeof (id));
            GBQualFree (gbq);
          } else if (StringICmp (gbq->qual, "comment") == 0 &&
                     StringDoesHaveText (gbq->val)) {
            *(prevqual) = gbq->next;
            gbq->next = NULL;
            comment = StringSave (gbq->val);
            GBQualFree (gbq);
          } else {
            prevqual = (GBQualPtr PNTR) &(gbq->next);
          }
          gbq = nextqual;
        }
        if (! StringHasNoText (id)) {
          if (StringChr (id, '|') != NULL) {
            sip = SeqIdParse (id);
          } else if (force_local_id) {
            sprintf (lcl, "lcl|%s", id);
            sip = SeqIdParse (lcl);
          } else {
            adv = ValidateAccnDotVer (id);
            if (adv == 0 || adv == -5) {
              ptr = StringChr (id, '.');
              if (ptr != NULL) {
                *ptr = '\0';
                ptr++;
                if (sscanf (ptr, "%ld", &val) == 1) {
                  version = (Uint4) val;
                }
              }
              sip = SeqIdFromAccession (id, version, NULL);
            } else {
              sprintf (lcl, "lcl|%s", id);
              sip = SeqIdParse (lcl);
            }
          }
        }
        if (sip != NULL || sfp->idx.subtype == FEATDEF_mRNA) {
          rnaseq = GetSequenceByFeature (sfp);
          if (rnaseq == NULL && seq_fetch_failP != NULL) {
            *seq_fetch_failP = TRUE;
          }
          if (rnaseq != NULL) {
            i = (Int4) StringLen (rnaseq);
            bs = BSNew (i + 2);
            if (bs != NULL) {
              BSWrite (bs, (VoidPtr) rnaseq, (Int4) StringLen (rnaseq));
              mbsp = BioseqNew ();
              if (mbsp != NULL) {
                mbsp->repr = Seq_repr_raw;
                mbsp->mol = Seq_mol_rna;
                mbsp->seq_data_type = Seq_code_iupacna;
                mbsp->seq_data = (SeqDataPtr) bs;
                mbsp->length = BSLen (bs);
                BioseqPack (mbsp);
                bs = NULL;
                /*
                sep = GetTopSeqEntryForEntityID (entityID);
                */
                old = SeqEntrySetScope (sep);
                if (sip != NULL) {
                  mbsp->id = sip;
                } else if (sfp->idx.subtype == FEATDEF_mRNA) {
                  /* actually just making rapid unique ID for mRNA */
                  mbsp->id = MakeNewProteinSeqIdEx (sfp->location, NULL, NULL, &ctr);
                }
                CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
                SeqMgrAddToBioseqIndex (mbsp);
                SeqEntrySetScope (old);
                msep = SeqEntryNew ();
                if (msep != NULL) {
                  msep->choice = 1;
                  msep->data.ptrvalue = (Pointer) mbsp;
                  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) mbsp, msep);
                  mip = MolInfoNew ();
                  if (mip != NULL) {
                    switch (sfp->idx.subtype) {
                      case FEATDEF_preRNA :
                        mip->biomol = MOLECULE_TYPE_PRE_MRNA;
                        break;
                      case FEATDEF_mRNA :
                        mip->biomol = MOLECULE_TYPE_MRNA;
                        break;
                      case FEATDEF_tRNA :
                        mip->biomol = MOLECULE_TYPE_TRNA;
                        break;
                      case FEATDEF_rRNA :
                        mip->biomol = MOLECULE_TYPE_RRNA;
                        break;
                      case FEATDEF_snRNA :
                        mip->biomol = MOLECULE_TYPE_SNRNA;
                        break;
                      case FEATDEF_scRNA :
                        mip->biomol = MOLECULE_TYPE_SCRNA;
                        break;
                      case FEATDEF_otherRNA :
                        mip->biomol = MOLECULE_TYPE_TRANSCRIBED_RNA;
                        break;
                      case FEATDEF_snoRNA :
                        mip->biomol = MOLECULE_TYPE_SNORNA;
                        break;
                      case FEATDEF_ncRNA :
                        mip->biomol = MOLECULE_TYPE_NCRNA;
                        break;
                      case FEATDEF_tmRNA :
                        mip->biomol = MOLECULE_TYPE_TMRNA;
                        break;
                      default :
                        mip->biomol = 0;
                        break;
                    }
                    if (partial5 && partial3) {
                      mip->completeness = 5;
                    } else if (partial5) {
                      mip->completeness = 3;
                    } else if (partial3) {
                      mip->completeness = 4;
                    }
                    vnp = CreateNewDescriptor (msep, Seq_descr_molinfo);
                    if (vnp != NULL) {
                      vnp->data.ptrvalue = (Pointer) mip;
                    }
                  }
                  if (comment != NULL) {
                    vnp = CreateNewDescriptor (msep, Seq_descr_comment);
                    if (vnp != NULL) {
                      vnp->data.ptrvalue = (Pointer) comment;
                    }
                  }
                  /* add mRNA sequence to genomic product set */
                  last->next = msep;
                  last = msep;
                  SetSeqFeatProduct (sfp, mbsp);
                }
              }
            }
            rnaseq = MemFree (rnaseq);
          }
        }
      }
      sfp = sfp->next;
    }
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }

  /* expand coding region features into protein product sequences */

  last = NULL;
  sfp = first;
  while (sfp != NULL) {
    prev = &(sfp->xref);
    xref = sfp->xref;
    while (xref != NULL) {
      next = xref->next;
      if (xref->data.choice == SEQFEAT_PROT &&
          sfp->data.choice == SEQFEAT_CDREGION &&
          sfp->product == NULL && (! sfp->pseudo)) {
        prp = (ProtRefPtr) xref->data.value.ptrvalue;
        ok_to_remove = TRUE;
        if (prp != NULL) {
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
/**
            crp->frame = 0;
**/
            bs = ProteinFromCdRegionEx (sfp, include_stop, remove_trailingX);
            if (bs == NULL && seq_fetch_failP != NULL) {
              *seq_fetch_failP = TRUE;
            }
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
                i = (Int4) StringLen (protseq);
                if (i > 0 && protseq [i - 1] == '*') {
                  protseq [i - 1] = '\0';
                }
                bs = BSNew (i + 2);
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
                pbsp->seq_data = (SeqDataPtr) bs;
                pbsp->length = BSLen (bs);
                bs = NULL;
                sep = NULL;
                mbsp = NULL;
                if (gen_prod_set) {
                  gbq = sfp->qual;
                  prevqual = (GBQualPtr PNTR) &(sfp->qual);
                  id [0] = '\0';
                  sip = NULL;
                  while (gbq != NULL) {
                    nextqual = gbq->next;
                    if (StringICmp (gbq->qual, "transcript_id") == 0) {
                      if (StringDoesHaveText (id) && StringDoesHaveText (gbq->val)) {
                        ErrPostEx (SEV_WARNING, ERR_FEATURE_QualWrongThisFeat,
                                   "CDS transcript_id %s replacing %s", gbq->val, id);
                      }
                      *(prevqual) = gbq->next;
                      gbq->next = NULL;
                      StringNCpy_0 (id, gbq->val, sizeof (id));
                      GBQualFree (gbq);
                    } else if (StringICmp (gbq->qual, "secondary_accession") == 0) {
                      *(prevqual) = gbq->next;
                      gbq->next = NULL;
                      shp = ParseStringIntoSeqHist (NULL, gbq->val);
                      if (shp != NULL) {
                        pbsp->hist = shp;
                      }
                      GBQualFree (gbq);
                    } else {
                      prevqual = (GBQualPtr PNTR) &(gbq->next);
                    }
                    gbq = nextqual;
                  }
                  if (StringHasNoText (id)) {
                    Message (MSG_POSTERR, "No transcript_id on CDS - unable to create nuc-prot set");
                  } else {
                    if (StringChr (id, '|') != NULL) {
                      sip = SeqIdParse (id);
                    } else if (force_local_id) {
                      sprintf (lcl, "lcl|%s", id);
                      sip = SeqIdParse (lcl);
                    } else {
                      adv = ValidateAccnDotVer (id);
                      if (adv == 0 || adv == -5) {
                        ptr = StringChr (id, '.');
                        if (ptr != NULL) {
                          *ptr = '\0';
                          ptr++;
                          if (sscanf (ptr, "%ld", &val) == 1) {
                            version = (Uint4) val;
                          }
                        }
                        sip = SeqIdFromAccession (id, version, NULL);
                      } else {
                        sprintf (lcl, "lcl|%s", id);
                        sip = SeqIdParse (lcl);
                      }
                    }
                  }
                  mbsp = BioseqFind (sip);
                  SeqIdFree (sip);
                  if (mbsp != NULL) {
                    sep = SeqMgrGetSeqEntryForData (mbsp);
                  /*
                  } else {
                    sep = GetBestTopParentForDataEx (entityID, bsp, TRUE);
                  */
                  }
                } else {
                  sep = GetBestTopParentForData (entityID, bsp);
                }
                if (sep == NULL) {
                  Message (MSG_POSTERR, "No location for nuc-prot set for CDS - unable to create nuc-prot set");
                  pbsp = BioseqFree (pbsp);
                  ok_to_remove = FALSE;
                } else {
                  old = SeqEntrySetScope (sep);
                  gbq = sfp->qual;
                  prevqual = (GBQualPtr PNTR) &(sfp->qual);
                  id [0] = '\0';
                  sip = NULL;
                  while (gbq != NULL) {
                    nextqual = gbq->next;
                    if (StringICmp (gbq->qual, "protein_id") == 0) {
                      if (StringDoesHaveText (id) && StringDoesHaveText (gbq->val)) {
                                ErrPostEx (SEV_WARNING, ERR_FEATURE_QualWrongThisFeat,
                          "CDS protein_id %s replacing %s", gbq->val, id);
                      }
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
                    } else if (force_local_id) {
                      sprintf (lcl, "lcl|%s", id);
                      sip = SeqIdParse (lcl);
                    } else {
                      adv = ValidateAccnDotVer (id);
                      if (adv == 0 || adv == -5) {
                        ptr = StringChr (id, '.');
                        if (ptr != NULL) {
                          *ptr = '\0';
                          ptr++;
                          if (sscanf (ptr, "%ld", &val) == 1) {
                            version = (Uint4) val;
                          }
                        }
                        sip = SeqIdFromAccession (id, version, NULL);
                      } else {
                        sprintf (lcl, "lcl|%s", id);
                        sip = SeqIdParse (lcl);
                      }
                    }
                  }
                  if (sip != NULL) {
                    pbsp->id = sip;
                  } else {
                    pbsp->id = MakeNewProteinSeqIdEx (sfp->location, NULL, NULL, &ctr);
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
                    /* the first protein may change the set/seq structure,
                    so goes through AddSeqEntryToSeqEntry */

                    if (gen_prod_set || last == NULL) {
                      descr = ExtractBioSourceAndPubs (sep);
                      AddSeqEntryToSeqEntry (sep, psep, TRUE);
                      ReplaceBioSourceAndPubs (sep, descr);
                      last = psep;
                    } else {
                      last->next = psep;
                      last = psep;
                    }
                    if (target == NULL) {
                      target = sep;
                      SaveSeqEntryObjMgrData (target, &omdptop, &omdata);
                      GetSeqEntryParent (target, &parentptr, &parenttype);
                    }
                    SetSeqFeatProduct (sfp, pbsp);
                    psep = SeqMgrGetSeqEntryForData (pbsp);
                    if (psep != NULL) {
                      last = psep;
                      prot = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
                      if (prot != NULL) {
                        prot->data.value.ptrvalue = (Pointer) prp;
                        SetSeqLocPartial (prot->location, partial5, partial3);
                        prot->partial = (Boolean) (partial5 || partial3);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        if (ok_to_remove) {
          xref->data.value.ptrvalue = NULL;
          *(prev) = next;
          xref->next = NULL;
          xref->data.choice = 0;
          SeqFeatXrefFree (xref);
        } else {
          prev = &(xref->next);
        }
      } else {
        prev = &(xref->next);
      }
      xref = next;
    }
    sfp = sfp->next;
  }

  /* expand mat_peptide features with protein_id qualifiers into protein product sequences */

  last = NULL;
  sfp = first;
  while (sfp != NULL) {
    if (sfp->data.choice == SEQFEAT_PROT && sfp->product == NULL) {
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      gbq = sfp->qual;
      prevqual = (GBQualPtr PNTR) &(sfp->qual);
      id [0] = '\0';
      sip = NULL;
      while (gbq != NULL) {
        nextqual = gbq->next;
        if (StringICmp (gbq->qual, "protein_id") == 0) {
          if (StringDoesHaveText (id) && StringDoesHaveText (gbq->val)) {
            ErrPostEx (SEV_WARNING, ERR_FEATURE_QualWrongThisFeat,
                       "Protein protein_id %s replacing %s",
                       gbq->val, id);
          }
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
        } else if (force_local_id) {
          sprintf (lcl, "lcl|%s", id);
          sip = SeqIdParse (lcl);
        } else {
          adv = ValidateAccnDotVer (id);
          if (adv == 0 || adv == -5) {
            ptr = StringChr (id, '.');
            if (ptr != NULL) {
              *ptr = '\0';
              ptr++;
              if (sscanf (ptr, "%ld", &val) == 1) {
                version = (Uint4) val;
              }
            }
            sip = SeqIdFromAccession (id, version, NULL);
          } else {
            sprintf (lcl, "lcl|%s", id);
            sip = SeqIdParse (lcl);
          }
        }
      }
      if (sip != NULL) {
        protseq = GetSequenceByFeature (sfp);
        if (protseq == NULL && seq_fetch_failP != NULL) {
          *seq_fetch_failP = TRUE;
        }
        if (protseq != NULL) {
          i = (Int4) StringLen (protseq);
          bs = BSNew (i + 2);
          if (bs != NULL) {
            BSWrite (bs, (VoidPtr) protseq, (Int4) StringLen (protseq));
            pbsp = BioseqNew ();
            if (pbsp != NULL) {
              pbsp->repr = Seq_repr_raw;
              pbsp->mol = Seq_mol_aa;
              pbsp->seq_data_type = Seq_code_ncbieaa;
              pbsp->seq_data = (SeqDataPtr) bs;
              pbsp->length = BSLen (bs);
              bs = NULL;
              /*
              sep = GetTopSeqEntryForEntityID (entityID);
              */
              sep = GetBestTopParentForData (entityID, bsp);
              old = SeqEntrySetScope (sep);
              if (sip != NULL) {
                pbsp->id = sip;
              } else {
                pbsp->id = MakeNewProteinSeqIdEx (sfp->location, NULL, NULL, &ctr);
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
                  mip->biomol = MOLECULE_TYPE_PEPTIDE;
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
                  AddSeqEntryToSeqEntry (sep, psep, TRUE);
                  last = psep;
                } else {
                  last->next = psep;
                  last = psep;
                }
                SetSeqFeatProduct (sfp, pbsp);
                if (prp != NULL) {
                  prp2 = AsnIoMemCopy ((Pointer) prp,
                                       (AsnReadFunc) ProtRefAsnRead,
                                       (AsnWriteFunc) ProtRefAsnWrite);
                  if (prp2 != NULL) {
                    psep = SeqMgrGetSeqEntryForData (pbsp);
                    if (psep != NULL) {
                      prot = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
                      if (prot != NULL) {
                        prot->data.value.ptrvalue = prp2;
                        SetSeqLocPartial (prot->location, partial5, partial3);
                        prot->partial = (Boolean) (partial5 || partial3);
                      }
                    }
                  }
                }
              }
            }
          }
          protseq = MemFree (protseq);
        }
      }
    }
    sfp = sfp->next;
  }

  if (target != NULL) {
    SeqMgrLinkSeqEntry (target, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (target, omdptop, &omdata);
  }
}

NLM_EXTERN void PromoteXrefsEx (SeqFeatPtr sfp, BioseqPtr bsp, Uint2 entityID, Boolean include_stop,
                                Boolean remove_trailingX, Boolean gen_prod_set)

{
  PromoteXrefsExEx (sfp, bsp, entityID, include_stop, remove_trailingX, gen_prod_set, FALSE, NULL);
}

NLM_EXTERN void PromoteXrefs (SeqFeatPtr sfp, BioseqPtr bsp, Uint2 entityID)

{
  PromoteXrefsExEx (sfp, bsp, entityID, TRUE, FALSE, FALSE, FALSE, NULL);
}
//LCOV_EXCL_STOP

/* begin BasicSeqEntryCleanup section */

static Boolean HasNoText (CharPtr str)

{
  Uchar  ch;    /* to use 8bit characters in multibyte languages */

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

//LCOV_EXCL_START
NLM_EXTERN CharPtr TrimSpacesAndSemicolons (CharPtr str)

{
  CharPtr  amp;
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && (ch <= ' ' || ch == ';')) {
      while (ch != '\0' && (ch <= ' ' || ch == ';')) {
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
      } else if (ch <= ' ') {
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
//LCOV_EXCL_STOP

NLM_EXTERN CharPtr TrimSpacesAndJunkFromEnds (
  CharPtr str,
  Boolean allowEllipsis
)

{
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  Boolean  isPeriod;
  Boolean  isTilde;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && (ch <= ' ' || ch == ',' || ch == ';')) {
      while (ch != '\0' && (ch <= ' ' || ch == ',' || ch == ';')) {
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
    dst = NULL;
    ptr = str;
    ch = *ptr;
    isPeriod = FALSE;
    isTilde = FALSE;
    while (ch != '\0') {
      if (ch <= ' ' || ch == '.' || ch == ',' || ch == '~' || ch == ';') {
        if (dst == NULL) {
          dst = ptr;
        }
        isPeriod = (Boolean) (isPeriod || ch == '.');
        isTilde = (Boolean) (isTilde || ch == '~');
      } else {
        dst = NULL;
        isPeriod = FALSE;
        isTilde = FALSE;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      /* allow one period at end */
      if (isPeriod) {
        *dst = '.';
        dst++;
        /* ellipsis are now okay */
        if (allowEllipsis && *dst == '.' && dst [1] == '.') {
          dst += 2;
        }
      } else if (isTilde) {
        /* allow double tilde at end */
        if (*dst == '~' && dst [1] == '~') {
          dst += 2;
        }
      }
      *dst = '\0';
    }
  }
  return str;
}

static CharPtr TrimSpacesSemicolonsAndCommas (CharPtr str)

{
  CharPtr  amp;
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && (ch <= ' ' || ch == ';' || ch == ',')) {
      while (ch != '\0' && (ch <= ' ' || ch == ';' || ch == ',')) {
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
      } else if (ch <= ' ') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
      } else if (ch == ';') {
        if (dst == NULL && amp == NULL) {
          dst = ptr;
        }
      } else if (ch == ',') {
        if (dst == NULL) {
          dst = ptr;
        }
        amp = NULL;
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

static CharPtr TrimFlankingQuotes (CharPtr str)

{
  size_t  len;

  if (str != NULL && str [0] != '\0') {
    len = StringLen (str);
    while (len > 0) {
      if (str [0] == '"' && str [len - 1] == '"') {
        str [0] = ' ';
        str [len - 1] = ' ';
      } else if (str [0] == '\'' && str [len - 1] == '\'') {
        str [0] = ' ';
        str [len - 1] = ' ';
      } else {
        return str;
      }
      TrimSpacesAroundString (str);
      len = StringLen (str);
    }
  }
  return str;
}

static void RemoveFlankingQuotes (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimFlankingQuotes (*strp);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void RemoveFlankingQuotesList (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimFlankingQuotes (vnp->data.ptrvalue);
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

static void CleanVisString (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesSemicolonsAndCommas (*strp);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void CleanVisStringAndCompress (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesSemicolonsAndCommas (*strp);
  Asn2gnbkCompressSpaces (*strp);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void CleanVisStringJunk (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesAndJunkFromEnds (*strp, TRUE);
  if (HasNoText (*strp)) {
    *strp = MemFree (*strp);
  }
}

static void CleanVisStringJunkAndCompress (CharPtr PNTR strp)

{
  if (strp == NULL) return;
  if (*strp == NULL) return;
  TrimSpacesAndJunkFromEnds (*strp, TRUE);
  Asn2gnbkCompressSpaces (*strp);
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

static CharPtr RemoveSpacesBetweenTildes (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
 CharPtr  ptr;
  CharPtr  tmp;

  if (str == NULL || str [0] == '\0') return str;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    *dst = ch;
    dst++;
    ptr++;
    if (ch == '~') {
      tmp = ptr;
      ch = *tmp;
      while (ch != 0 && ch <= ' ') {
        tmp++;
        ch = *tmp;
      }
      if (ch == '~') {
        ptr = tmp;
      }
    }
    ch = *ptr;
  }
  *dst = '\0';

  return str;
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
    TrimSpacesSemicolonsAndCommas (vnp->data.ptrvalue);
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

static void CleanVisStringJunkListAndCompress (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimSpacesSemicolonsAndCommas (vnp->data.ptrvalue);
    TrimSpacesAndJunkFromEnds (vnp->data.ptrvalue, TRUE);
    Asn2gnbkCompressSpaces (vnp->data.ptrvalue);
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

static void CleanVisStringListAndCompress (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimSpacesSemicolonsAndCommas (vnp->data.ptrvalue);
    Asn2gnbkCompressSpaces (vnp->data.ptrvalue);
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

static Boolean AlreadyInVnpListCaseSensitive (ValNodePtr head, ValNodePtr curr)

{
  if (head == NULL || curr == NULL) return FALSE;
  /* since we cannot sort these lists, must check against all previous entries */
  while (head != curr && head != NULL) {
    if (StringCmp (head->data.ptrvalue, curr->data.ptrvalue) == 0) return TRUE;
    head = head->next;
  }
  return FALSE;
}

static void CleanVisStringListCaseSensitive (ValNodePtr PNTR vnpp)

{
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (vnpp == NULL) return;
  prev = vnpp;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    TrimSpacesSemicolonsAndCommas (vnp->data.ptrvalue);
    if (HasNoText (vnp->data.ptrvalue) || AlreadyInVnpListCaseSensitive (*vnpp, vnp)) {
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

  if (StringICmp (gbq->qual, "map") == 0) {
    choice = 2;
  } else if (StringICmp (gbq->qual, "allele") == 0) {
    choice = 3;
  } else if (StringICmp (gbq->qual, "locus_tag") == 0) {
    choice = 4;
  } else if (StringICmp (gbq->qual, "old_locus_tag") == 0) {
    choice = 5;
  } else if (StringICmp (gbq->qual, "gene_synonym") == 0) {
    choice = 6;
  }
  if (choice > 0) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp == NULL) return FALSE;
    switch (choice) {
      case 2 :
        if (grp->maploc != NULL) return FALSE;
        if (StringHasNoText (gbq->val)) return FALSE;
        grp->maploc = StringSave (gbq->val);
        break;
      case 3 :
        if (StringHasNoText (gbq->val)) return FALSE;
        if (grp->allele != NULL) {
          if (StringICmp (gbq->val, grp->allele) == 0) return TRUE;
          return FALSE;
        }
        grp->allele = StringSave (gbq->val);
        break;
      case 4 :
        if (grp->locus_tag != NULL) return FALSE;
        if (StringHasNoText (gbq->val)) return FALSE;
        grp->locus_tag = StringSave (gbq->val);
        break;
      case 5 :
/* removed by indexer request */
/*        if (StringHasNoText (gbq->val)) return FALSE;
 *       if (grp->locus_tag != NULL) {
 *         if (StringICmp (gbq->val, grp->locus_tag) == 0) return TRUE;
 *         return FALSE;
 *       }
 */
        return FALSE;
        break;
      case 6 :
        if (StringHasNoText (gbq->val)) return FALSE;
        ValNodeCopyStr (&(grp->syn), 0, gbq->val);
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
    } else {
        ErrPostEx (SEV_WARNING, ERR_QUALIFIER_InvalidDataFormat,
                   "bad transl_except %s", qval);
        str = StringStr(qval, ",");
        if (str != NULL) {
            str = StringStr(str, ":");
            if (str != NULL) {
              str++;
            }
        }
    }

    if (str == NULL) return (Uint1) 'X';

       while (*str == ' ')
           ++str;
       for (eptr = str; *eptr != ')' && *eptr != ' ' && *eptr != '\0';  eptr++) continue;

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
   eptr = StringStr (bptr, ",aa:");
   if (eptr == NULL) {
     for (eptr = bptr; *eptr != ',' && *eptr != '\0'; eptr++) continue;
   }
   if (eptr == NULL) return NULL;

   return (TextSave(bptr, eptr-bptr));
}

//LCOV_EXCL_START
extern Boolean ParseAnticodon (SeqFeatPtr sfp, CharPtr val, Int4 offset);
extern Boolean ParseAnticodon (SeqFeatPtr sfp, CharPtr val, Int4 offset)

{
  Int4       diff;
  Int2       j;
  Boolean    locmap;
  int        num_errs;
  CharPtr    pos;
  Boolean    pos_range = FALSE;
  RnaRefPtr  rrp;
  SeqIntPtr  sintp;
  SeqIdPtr   sip;
  Boolean    sitesmap;
  SeqLocPtr  slp;
  SeqPntPtr  spp;
  Uint1      strand;
  Int4       temp;
  tRNAPtr    trp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return FALSE;
  if (StringHasNoText (val)) return FALSE;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return FALSE;

  if (rrp->ext.choice == 0 && rrp->ext.value.ptrvalue == NULL) {
    rrp->ext.choice = 2;
    trp = (tRNAPtr) MemNew (sizeof (tRNA));
    rrp->ext.value.ptrvalue = (Pointer) trp;
    if (trp != NULL) {
      trp->aatype = 2;
      for (j = 0; j < 6; j++) {
        trp->codon [j] = 255;
      }
    }
  }
  if (rrp->ext.choice != 2) return FALSE;

  trp = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trp == NULL) return FALSE;
      
  /* find SeqId to use */
  sip = SeqLocId (sfp->location);
  if (sip == NULL) {
    slp = SeqLocFindNext (sfp->location, NULL);
    if (slp != NULL) {
      sip = SeqLocId (slp);
    }
  }
  if (sip == NULL) return FALSE;

  /* parse location */
  pos = SimpleValuePos (val);
  if (pos == NULL) {
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "anticodon parsing failed, %s, drop the anticodon", val);
    return FALSE;
  }

  trp->anticodon = Nlm_gbparseint (pos, &locmap, &sitesmap, &num_errs, sip);
  if (trp->anticodon == NULL) {
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "anticodon parsing failed, %s, drop the anticodon", pos);
    MemFree (pos);
    return FALSE;
  }

  if (trp->anticodon->choice == SEQLOC_PNT) {
    /* allow a single point */
    spp = trp->anticodon->data.ptrvalue;
    if (spp != NULL) {
      spp->point += offset;
    }
  }
  if (trp->anticodon->choice == SEQLOC_INT) {
    sintp = trp->anticodon->data.ptrvalue;
    if (sintp == NULL) {
      MemFree (pos);
      return FALSE;
    }
    sintp->from += offset;
    sintp->to += offset;
    if (sintp->from > sintp->to) {
      temp = sintp->from;
      sintp->from = sintp->to;
      sintp->to = temp;
    }
    sintp->strand = SeqLocStrand (sfp->location);
    strand = sintp->strand;
    diff = SeqLocStop(trp->anticodon) - SeqLocStart(trp->anticodon); /* SeqLocStop/Start does not do what you think */
    /*
    if ((diff != 2 && (strand != Seq_strand_minus)) ||
        (diff != -2 && (strand == Seq_strand_minus))) {
      pos_range = TRUE;
    }
    */
    if (diff != 2) {
      pos_range = TRUE;
    }
    if (num_errs > 0 || pos_range) {
      ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
                 "anticodon range is wrong, %s, drop the anticodon", pos);
      MemFree (pos);
      return FALSE;
    }
    if (SeqLocCompare (sfp->location, trp->anticodon) != SLC_B_IN_A) {
      ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
                 "/anticodon not in tRNA: %s", val);
      MemFree (pos);
      return FALSE;
    }
  }

  MemFree (pos);

  return TRUE;
}
//LCOV_EXCL_STOP

extern Boolean ParseCodeBreak (SeqFeatPtr sfp, CharPtr val, Int4 offset);
extern Boolean ParseCodeBreak (SeqFeatPtr sfp, CharPtr val, Int4 offset)

{
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  Int4          diff;
  CodeBreakPtr  lastcbp;
  Boolean       locmap;
  int           num_errs;
  Boolean       packed_int = TRUE;
  CharPtr       pos;
  Boolean       pos_range = FALSE;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  Boolean       sitesmap;
  SeqLocPtr     slp;
  SeqLocPtr     slp1, slp2;
  SeqPntPtr     spp;
  Uint1         strand;
  Int4          temp;
  CharPtr       tmp;

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

  cbp = CodeBreakNew ();
  if (cbp == NULL) return FALSE;
  cbp->aa.choice = 1; /* ncbieaa */
  cbp->aa.value.intvalue = (Int4) GetQualValueAa (val);

  /* parse location */
  pos = SimpleValuePos (val);
  if (pos == NULL) {
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "transl_except parsing failed, %s, drop the transl_except", val);
    return FALSE;
  }
  if (StringChr (pos, ',') != NULL) {
    tmp = (CharPtr) MemNew ((StringLen (pos) + 10) * sizeof (Char));
    if (tmp != NULL) {
      sprintf (tmp, "join(%s)", pos);
      MemFree (pos);
      pos = tmp;
    }
  }
  cbp->loc = Nlm_gbparseint (pos, &locmap, &sitesmap, &num_errs, sip);
  if (cbp->loc == NULL) {
    CodeBreakFree (cbp);
    ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
               "transl_except parsing failed, %s, drop the transl_except", pos);
    MemFree (pos);
    return FALSE;
  }
  if (cbp->loc->choice == SEQLOC_PNT) {
    /* allow a single point */
    spp = cbp->loc->data.ptrvalue;
    if (spp != NULL) {
      spp->point += offset;
    }
  } else if (cbp->loc->choice == SEQLOC_INT) {
    sintp = cbp->loc->data.ptrvalue;
    if (sintp == NULL) {
      MemFree (pos);
      return FALSE;
    }
    sintp->from += offset;
    sintp->to += offset;
    if (sintp->from > sintp->to) {
      temp = sintp->from;
      sintp->from = sintp->to;
      sintp->to = temp;
    }
    sintp->strand = SeqLocStrand (sfp->location);
    strand = sintp->strand;
    diff = SeqLocStop(cbp->loc) - SeqLocStart(cbp->loc); /* SeqLocStop/Start does not do what you think */
    /*
    if ((diff != 2 && (strand != Seq_strand_minus)) ||
        (diff != -2 && (strand == Seq_strand_minus))) {
      pos_range = TRUE;
    }
    */
    if (diff != 2) {
      pos_range = TRUE;
    }
    if (num_errs > 0 || pos_range) {
      CodeBreakFree (cbp);
      ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
                 "transl_except range is wrong, %s, drop the transl_except", pos);
      MemFree (pos);
      return FALSE;
    }
    if (SeqLocCompare (sfp->location, cbp->loc) != SLC_B_IN_A) {
      CodeBreakFree (cbp);
      ErrPostEx (SEV_WARNING, ERR_FEATURE_LocationParsing,
                 "/transl_except not in CDS: %s", val);
      MemFree (pos);
      return FALSE;
    }
  } else {
    slp1 = dnaLoc_to_aaLoc (sfp, cbp->loc, TRUE, NULL, TRUE);
    if (slp1 != NULL) {
      slp2 = aaLoc_to_dnaLoc (sfp, slp1);
      if (slp2 != NULL) {
        SeqLocFree (cbp->loc);
        cbp->loc = slp2;
      }
      SeqLocFree (slp1);
    }
    slp = SeqLocFindNext (cbp->loc, NULL);
    while (slp != NULL) {
      if (slp->choice == SEQLOC_PNT) {
        spp = slp->data.ptrvalue;
        if (spp != NULL) {
          sintp = SeqIntNew();
          if (sintp != NULL) {
            sintp->id = SeqIdDup (spp->id);
            sintp->from = spp->point;
            sintp->to = spp->point;
            sintp->strand = SeqLocStrand (sfp->location);
            slp->choice = SEQLOC_INT;
            slp->data.ptrvalue = sintp;
            SeqPntFree (spp);
          }
        }
      }
      if (slp->choice == SEQLOC_INT) {
        sintp = slp->data.ptrvalue;
        if (sintp == NULL) {
          MemFree (pos);
          return FALSE;
        }
        sintp->from += offset;
        sintp->to += offset;
        if (sintp->from > sintp->to) {
          temp = sintp->from;
          sintp->from = sintp->to;
          sintp->to = temp;
        }
        sintp->strand = SeqLocStrand (sfp->location);
      } else {
        packed_int = FALSE;
      }
      slp = SeqLocFindNext (cbp->loc, slp);
    }
    slp = cbp->loc;
    if (packed_int && slp->choice == SEQLOC_MIX) {
      slp->choice = SEQLOC_PACKED_INT;
    }
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
  MemFree (pos);
  return TRUE;
}

static Boolean CodonsAlreadyInOrder (tRNAPtr trp)

{
  Int2  i, j;

  if (trp == NULL) return TRUE;
  for (i = 0, j = 1; i < 5; i++, j++) {
    if (trp->codon [i] > trp->codon [j]) return FALSE;
  }
  return TRUE;
}

static int LIBCALLBACK SortCodons (VoidPtr ptr1, VoidPtr ptr2)

{
  Uint1  codon1, codon2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  codon1 = *((Uint1Ptr) ptr1);
  codon2 = *((Uint1Ptr) ptr2);
  if (codon1 > codon2) {
    return 1;
  } else if (codon1 < codon2) {
    return -1;
  }
  return 0;
}

static void UniqueCodons (tRNAPtr trp)

{
  Int2   i, j;
  Uint1  last = 255, next;

  if (trp == NULL) return;

  for (i = 0, j = 0; i < 6; i++) {
    next = trp->codon [i];
    if (next != last) {
      trp->codon [j] = next;
      last = next;
      j++;
    }
  }
  while (j < 6) {
    trp->codon [j] = 255;
    j++;
  }
}

static CharPtr  codonLetterExpand [] =
{
  "?", "A", "C", "AC",
  "G", "AG", "CG", "ACG",
  "T", "AT", "CT", "ACT",
  "GT", "AGT", "CGT", "ACGT",
  NULL
};

NLM_EXTERN Boolean ParseDegenerateCodon (tRNAPtr trp, Uint1Ptr codon)

{
  Uint1    ch;
  Uint1    chrToInt [256];
  Int2     k;
  Uint1    i, j;
  Uint1    idx;
  CharPtr  intToChr = "?ACMGRSVTWYHKDBN";
  CharPtr  ptr, str;

  if (trp == NULL || codon == NULL) return FALSE;

  for (i = 0; i < 2; i++) {
    ch = codon [i];
    if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T') return FALSE;
  }

  for (k = 0; k < 256; k++) {
    chrToInt [k] = 0;
  }
  for (i = 1; i < 16; i++) {
    ch = intToChr [i];
    chrToInt [(int) ch] = i;
  }

  idx = chrToInt [(int) codon [2]];
  if (idx > 15) return FALSE;

  str = codonLetterExpand [idx];
  ptr = str;
  ch = *ptr;
  j = 0;
  codon [3] = '\0';
  while (ch != '\0' && j < 6) {
    codon [2] = ch;
    trp->codon [j] = IndexForCodon (codon, Seq_code_iupacna);
    ptr++;
    ch = *ptr;
    j++;
  }

  return TRUE;
}

static void CleanupTrna (SeqFeatPtr sfp, tRNAPtr trp)

{
  Uint1           aa = 0;
  Uint1           curraa;
  Uint1           from = 0;
  Int2            j;
  Boolean         justTrnaText;
  SeqMapTablePtr  smtp;
  Uint1           trpcodon [6];
  /*
  Char            codon [16];
  Int2            i;
  Boolean         okayToFree = TRUE;
  CharPtr         str;
  */

  /* look for tRNA-OTHER with actual amino acid in comment */

  if (trp == NULL) return;

  /*
  if (sfp != NULL && sfp->comment != NULL && trp->codon [0] == 255) {
    codon [0] = '\0';
    if (StringNICmp (sfp->comment, "codon recognized: ", 18) == 0) {
      StringNCpy_0 (codon, sfp->comment + 18, sizeof (codon));
    } else if (StringNICmp (sfp->comment, "codons recognized: ", 19) == 0) {
      StringNCpy_0 (codon, sfp->comment + 19, sizeof (codon));
    }
    if (StringDoesHaveText (codon)) {
      if (StringLen (codon) > 3 && codon [3] == ';') {
        codon [3] = '\0';
        okayToFree = FALSE;
      }
      if (StringLen (codon) == 3) {
        for (i = 0; i < 3; i++) {
          if (codon [i] == 'U') {
            codon [i] = 'T';
          }
        }
        if (ParseDegenerateCodon (trp, (Uint1Ptr) codon)) {
          if (okayToFree) {
            sfp->comment = MemFree (sfp->comment);
          } else {
            str = StringSave (sfp->comment + 22);
            TrimSpacesAroundString (str);
            sfp->comment = MemFree (sfp->comment);
            if (StringHasNoText (str)) {
              str = MemFree (str);
            }
            sfp->comment = str;
          }
        }
      }
    }
  }
  */

  if (! CodonsAlreadyInOrder (trp)) {
    StableMergeSort ((VoidPtr) &(trp->codon), 6, sizeof (Uint1), SortCodons);
  }
  UniqueCodons (trp);

  /* now always switch iupacaa to ncbieaa (was just for selenocysteine) */

  if (trp->aatype == 1 /* && trp->aa == 'U' */) {
    trp->aatype = 2;
  }

  if (sfp == NULL || sfp->comment == NULL) return;

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
  if (aa != 'X') {
    curraa = ParseTRnaString (sfp->comment, &justTrnaText, trpcodon, TRUE);
    if (aa == 0 && curraa != 0) {
      aa = curraa;
      trp->aa = curraa;
      trp->aatype = 2;
    }
    if (aa != 0 && aa == curraa) {
      if (justTrnaText) {
        for (j = 0; j < 6; j++) {
          if (trp->codon [j] == 255) {
            trp->codon [j] = trpcodon [j];
          }
        }
        if (StringCmp (sfp->comment, "fMet") != 0 && StringCmp (sfp->comment, "iMet") != 0) {
          sfp->comment = MemFree (sfp->comment);
        }
      }
    }
    return;
  }
  aa = ParseTRnaString (sfp->comment, &justTrnaText, trpcodon, TRUE);
  if (aa == 0) return;
  trp->aa = aa;
  trp->aatype = 2;
  if (justTrnaText) {
    for (j = 0; j < 6; j++) {
      if (trp->codon [j] == 255) {
        trp->codon [j] = trpcodon [j];
      }
    }
    if (StringCmp (sfp->comment, "fMet") != 0 && StringCmp (sfp->comment, "iMet") != 0) {
      sfp->comment = MemFree (sfp->comment);
    }
  }
}

NLM_EXTERN SeqFeatPtr LIBCALL GetBestProteinFeatureUnindexed (SeqLocPtr product)

{
  BioseqPtr    bsp;
  SeqFeatPtr   prot = NULL;
  SeqAnnotPtr  sap;
  SeqFeatPtr   tmp;
  ValNode      vn;

  if (product == NULL) return NULL;
  bsp = BioseqFindFromSeqLoc (product);
  if (bsp == NULL || bsp->repr != Seq_repr_raw) return NULL;
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
  return prot;
}

static void CleanupECNumber (CharPtr str)

{
  size_t len;

  len = StringLen (str);
  if (len < 1) return;
  if (str [len - 1] == '.') {
    str [len - 1] = ' ';
  }
  if (StringNICmp (str, "EC ", 3) == 0) {
    str [0] = ' ';
    str [1] = ' ';
  } else if (StringNICmp (str, "EC:", 3) == 0) {
    str [0] = ' ';
    str [1] = ' ';
    str [2] = ' ';
  }
  TrimSpacesAroundString (str);
}

static Boolean ECNumberCanBeSplit (CharPtr str)

{
  Char     ch;
  CharPtr  ptr;

  if (StringHasNoText (str)) return FALSE;

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if ((! IS_DIGIT (ch)) && ch != '.' && ch !='-' && ch !='n' && ch != ' ' && ch !=';') return FALSE;
    ptr++;
    ch = *ptr;
  }

  return TRUE;
}

static Boolean HandledGBQualOnCDS (SeqFeatPtr sfp, GBQualPtr gbq, ValNodePtr PNTR afterMe)

{
  Int2            choice = 0;
  CdRegionPtr     crp;
  Uint1           frame;
  ValNodePtr      gcp;
  ValNodePtr      prev;
  SeqFeatPtr      prot;
  ProtRefPtr      prp = NULL;
  Char            str [16];
  Int4            transl_table;
  int             val;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (StringICmp (gbq->qual, "product") == 0) {
    choice = 1;
  } else if (StringICmp (gbq->qual, "function") == 0) {
    choice = 2;
  } else if (StringICmp (gbq->qual, "EC_number") == 0) {
    choice = 3;
  } else if (StringICmp (gbq->qual, "prot_note") == 0) {
    choice = 4;
  }
  if (choice > 0) {
    prot = GetBestProteinFeatureUnindexed (sfp->product);
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
        if (prot != NULL && prot->data.value.ptrvalue != NULL) {
          if (*afterMe == NULL) {
            /* if protein product exists, product gbqual becomes first name */
            vnp = ValNodeCopyStr (NULL, 0, gbq->val);
            if (vnp != NULL) {
              vnp->next = prp->name;
              prp->name = vnp;
            }
            *afterMe = vnp;
          } else {
            vnp = ValNodeCopyStr (NULL, 0, gbq->val);
            prev = *afterMe;
            if (vnp != NULL) {
              vnp->next = prev->next;
              prev->next = vnp;
            }
            *afterMe = vnp;
          }
        } else {
          /* if local xref, append to name */
          ValNodeCopyStr (&(prp->name), 0, gbq->val);
        }
        break;
      case 2 :
        ValNodeCopyStr (&(prp->activity), 0, gbq->val);
        break;
      case 3 :
        ValNodeCopyStr (&(prp->ec), 0, gbq->val);
        break;
      case 4 :
        if (prot == NULL) {
          return FALSE;
        } else {
          prot->comment = StringSave (gbq->val);
        }
        break;
      default :
        break;
    }
    return TRUE;
  }

  if (StringICmp (gbq->qual, "transl_except") == 0) {
    return ParseCodeBreak (sfp, gbq->val, 0);
  }

  if (StringICmp (gbq->qual, "codon_start") == 0) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      frame = crp->frame;
      if (frame == 0) {
        StringNCpy_0 (str, gbq->val, sizeof (str));
        if (sscanf (str, "%d", &val) == 1) {
          if (val > 0 && val < 4) {
            crp->frame = (Uint1) val;
            return TRUE;
          }
        }
        frame = 1;
      }
      sprintf (str, "%d", (int) frame);
      if (StringICmp (str, gbq->val) == 0) {
        return TRUE;
      } else if (sfp->pseudo && sfp->product == NULL) {
        StringNCpy_0 (str, gbq->val, sizeof (str));
        if (sscanf (str, "%d", &val) == 1) {
          if (val > 0 && val < 4) {
            crp->frame = (Uint1) val;
            return TRUE;
          }
        }
      }
    }
  }

  if (StringICmp (gbq->qual, "transl_table") == 0) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      transl_table = 0;
      gcp = crp->genetic_code;
      if (gcp != NULL) {
        for (vnp = gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
          if (vnp->choice == 2 && vnp->data.intvalue != 0) {
            transl_table = vnp->data.intvalue;
          }
        }
        if (transl_table == 0) {
          transl_table = 1;
        }
        sprintf (str, "%ld", (long) transl_table);
        if (StringICmp (str, gbq->val) == 0) {
          return TRUE;
        }
      } else {
        StringNCpy_0 (str, gbq->val, sizeof (str));
        if (sscanf (str, "%d", &val) == 1) {
          vnp = ValNodeNew (NULL);
          if (vnp != NULL) {
            vnp->choice = 2;
            vnp->data.intvalue = (Int4) val;
            gcp = GeneticCodeNew ();
            if (gcp != NULL) {
              gcp->data.ptrvalue = vnp;
              crp->genetic_code = gcp;
              return TRUE;
            }
          }
        }
      }
    }
  }

  if (StringICmp (gbq->qual, "translation") == 0) {
    return TRUE;
  }

  return FALSE;
}


static Boolean HandledGBQualOnRNA (SeqFeatPtr sfp, GBQualPtr gbq, Boolean isEmblOrDdbj)

{
  Uint1      aa;
  BioseqPtr  bsp;
  Uint1      codon [6];
  Boolean    emptyRNA;
  Int4       from;
  Boolean    is_fMet = FALSE;
  Boolean    is_iMet = FALSE;
  Boolean    is_std_name = FALSE;
  Int2       j;
  Boolean    justTrnaText;
  size_t     len;
  CharPtr    name;
  CharPtr    ptr;
  RNAGenPtr  rgp;
  RnaRefPtr  rrp;
  SeqIntPtr  sintp;
  SeqIdPtr   sip;
  CharPtr    str;
  Char       tmp [64];
  Int4       to;
  tRNAPtr    trp;
  long int   val;

  is_std_name = (Boolean) (StringICmp (gbq->qual, "standard_name") == 0);
  if (StringICmp (gbq->qual, "product") == 0 ||
      (is_std_name && (! isEmblOrDdbj) )) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp == NULL) return FALSE;
    if (rrp->type == 0) {
      rrp->type = 255;
    }
    if (rrp->type == 255 && is_std_name) return FALSE;
    if (rrp->ext.choice == 1) {
      name = (CharPtr) rrp->ext.value.ptrvalue;
      if (StringHasNoText (name)) {
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        rrp->ext.choice = 0;
      }
    }
    if (rrp->ext.choice == 2) {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL) {
        if (trp->aatype == 0 && trp->aa == 0 && trp->anticodon == NULL) {
          emptyRNA = TRUE;
          for (j = 0; j < 6; j++) {
            if (trp->codon [j] != 255) {
              emptyRNA = FALSE;
            }
          }
          if (emptyRNA) {
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
            rrp->ext.choice = 0;
          }
        }
      }
    }
    if (rrp->type == 3 && rrp->ext.choice == 1) {
      name = (CharPtr) rrp->ext.value.ptrvalue;
      aa = ParseTRnaString (name, &justTrnaText, codon, FALSE);
      if (aa != 0) {
        is_fMet = (Boolean) (StringStr (name, "fMet") != NULL);
        is_iMet = (Boolean) (StringStr (name, "iMet") != NULL);
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        trp = (tRNAPtr) MemNew (sizeof (tRNA));
        if (trp != NULL) {
          trp->aatype = 2;
          for (j = 0; j < 6; j++) {
            trp->codon [j] = 255;
          }
          if (justTrnaText) {
            for (j = 0; j < 6; j++) {
              trp->codon [j] = codon [j];
            }
          }
          trp->aa = aa;
          rrp->ext.choice = 2;
          rrp->ext.value.ptrvalue = (Pointer) trp;
          if (aa == 'M') {
            if (is_fMet) {
              if (sfp->comment == NULL) {
                sfp->comment = StringSave ("fMet");
              } else {
                len = StringLen (sfp->comment) + StringLen ("fMet") + 5;
                str = MemNew (sizeof (Char) * len);
                StringCpy (str, sfp->comment);
                StringCat (str, "; ");
                StringCat (str, "fMet");
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = str;
              }
            }
            if (is_iMet) {
              if (sfp->comment == NULL) {
                sfp->comment = StringSave ("iMet");
              } else {
                len = StringLen (sfp->comment) + StringLen ("iMet") + 5;
                str = MemNew (sizeof (Char) * len);
                StringCpy (str, sfp->comment);
                StringCat (str, "; ");
                StringCat (str, "iMet");
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = str;
              }
            }
          }
          CleanupTrna (sfp, trp);
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
        if (trp->aa == 77) {
          if (StringICmp (gbq->val, "tRNA-fMet") == 0 || StringICmp (gbq->val, "tRNA-iMet") == 0) return FALSE;
        }
        if (trp->aa == ParseTRnaString (gbq->val, NULL, NULL, FALSE)) {
          return TRUE;
        }
      }
    }
    if (rrp->ext.choice == 3) {
      rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
      if (rgp == NULL) return FALSE;
      if (StringHasNoText (rgp->product)) {
        rgp->product = StringSave (gbq->val);
        return TRUE;
      }
      return FALSE;
    }
    if (rrp->ext.choice != 0 && rrp->ext.choice != 1) return FALSE;
    name = (CharPtr) rrp->ext.value.ptrvalue;
    if (! HasNoText (name)) {
      if (StringICmp (name, gbq->val) == 0) {
        return TRUE;
      }
      str = StringStr (gbq->val, "rDNA");
      if (str != NULL) {
        str [1] = 'R';
        if (StringICmp (name, gbq->val) == 0) {
          return TRUE;
        }
      }
      if (rrp->type == 255 || rrp->type == 8 || rrp->type == 9 || rrp->type == 10) {
        /* new convention follows ASN.1 spec comments, allows new RNA types */
        return FALSE;
      }
      /* subsequent /product now added to comment */
      if (sfp->comment == NULL) {
        sfp->comment = gbq->val;
        gbq->val = NULL;
      } else if (StringStr (gbq->val, sfp->comment) == NULL) {
        len = StringLen (sfp->comment) + StringLen (gbq->val) + 5;
        str = MemNew (sizeof (Char) * len);
        StringCpy (str, sfp->comment);
        StringCat (str, "; ");
        StringCat (str, gbq->val);
        sfp->comment = MemFree (sfp->comment);
        sfp->comment = str;
      }
      /* return FALSE; */
      return TRUE;
    }
    if (rrp->type == 8 || rrp->type == 9 || rrp->type == 10) {
      /* new convention follows ASN.1 spec comments, allows new RNA types */
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
  } else if (StringICmp (gbq->qual, "anticodon") == 0) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp == NULL) return FALSE;
    if (rrp->type == 0) {
      rrp->type = 255;
    }
    if (rrp->type == 3 && rrp->ext.choice == 0) {
      trp = (tRNAPtr) MemNew (sizeof (tRNA));
      if (trp != NULL) {
        rrp->ext.choice = 2;
        rrp->ext.value.ptrvalue = trp;
        for (j = 0; j < 6; j++) {
          trp->codon [j] = 255;
        }
      }
    }
    if (rrp->type == 3 && rrp->ext.choice == 2) {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL) {
        StringNCpy_0 (tmp, gbq->val, sizeof (tmp));
        ptr = StringStr (tmp, "(");
        if (ptr != NULL) {
          ptr = StringStr (ptr + 1, "pos");
          if (ptr != NULL) {
            ptr = StringStr (ptr + 3, ":");
          }
        }
        if (ptr != NULL) {
          str = ptr + 1;
          ptr = StringStr (str, "..");
          if (ptr != NULL) {
            *ptr = '\0';
            if (sscanf (str, "%ld", &val) == 1) {
              from = val - 1;
              str = ptr + 2;
              ptr = StringStr (str, ",");
              if (ptr != NULL) {
                *ptr = '\0';
                if (sscanf (str, "%ld", &val) == 1) {
                  to = val - 1;
                  sip = SeqLocId (sfp->location);
                  if (sip != NULL) {
                    bsp = BioseqFind (sip);
                    if (bsp != NULL) {
                      if (from >= 0 && from < bsp->length - 1) {
                        if (to >= 0 && to < bsp->length - 1) {
                          sintp = SeqIntNew ();
                          if (sintp != NULL) {
                            if (from > to) {
                              sintp->from = to;
                              sintp->to = from;
                              sintp->strand = Seq_strand_minus;
                            } else {
                              sintp->from = from;
                              sintp->to = to;
                              sintp->strand = Seq_strand_plus;
                            }
                            sintp->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
                            trp->anticodon = ValNodeAddPointer (NULL, SEQLOC_INT, (Pointer) sintp);
                            if (trp->aatype == 0 && trp->aa == 0) {
                              ptr = StringStr (ptr + 1, "aa:");
                              if (ptr != NULL) {
                                str = ptr + 3;
                                ptr = StringStr (str, ")");
                                if (ptr != NULL) {
                                  *ptr = '\0';
                                  trp->aa = ParseTRnaString (str, NULL, NULL, FALSE);
                                  if (trp->aa != 0) {
                                    trp->aatype = 2;
                                  }
                                }
                              }
                            }
                            return TRUE;
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
  } else if (StringICmp (gbq->qual, "standard_name") == 0) {
    choice = 4;
  } else if (StringICmp (gbq->qual, "label") == 0) {
    choice = 5;
  } else if (StringICmp (gbq->qual, "allele") == 0) {
      choice = 6;
  }
  if (choice == 1 || choice == 4) {
    vnp = prp->name;
    if (vnp != NULL && (! HasNoText (vnp->data.ptrvalue))) return FALSE;
    ValNodeCopyStr (&(prp->name), 0, gbq->val);
    /*
    vnp = prp->name;
    if (vnp != NULL && prp->desc != NULL) {
      if (StringICmp (vnp->data.ptrvalue, prp->desc) == 0) {
        prp->desc = MemFree (prp->desc);
      }
    }
    */
    return TRUE;
  } else if (choice == 2) {
    ValNodeCopyStr (&(prp->activity), 0, gbq->val);
    return TRUE;
  } else if (choice == 3) {
    ValNodeCopyStr (&(prp->ec), 0, gbq->val);
    return TRUE;
  } else if (choice == 5) {
    return FALSE; /* keep label gbqual only */
  } else if (choice == 6) {
      return FALSE;
  }

  if (StringICmp (gbq->qual, "experiment") == 0 ||
      StringICmp (gbq->qual, "inference") == 0) {
    return FALSE;
  }

  if (StringICmp (gbq->qual, "UniProtKB_evidence") == 0) {
    return FALSE;
  }

  return TRUE; /* all other gbquals not appropriate on protein features */
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
    /* return TRUE; */
  }
  return FALSE;
}

static void CleanupRptUnit (GBQualPtr gbq)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;
  len = StringLen (gbq->val) * 2 + 1;
  str = MemNew (sizeof (Char) * len);
  if (str == NULL) return;
  ptr = str;
  tmp = gbq->val;
  ch = *tmp;
  while (ch != '\0') {
    while (ch == '(' || ch == ')' || ch == ',') {
      *ptr = ch;
      ptr++;
      tmp++;
      ch = *tmp;
    }
    while (IS_WHITESP (ch)) {
      tmp++;
      ch = *tmp;
    }
    while (IS_DIGIT (ch)) {
      *ptr = ch;
      ptr++;
      tmp++;
      ch = *tmp;
    }
    if (ch == '.' || ch == '-') {
      while (ch == '.' || ch == '-') {
        tmp++;
        ch = *tmp;
      }
      *ptr = '.';
      ptr++;
      *ptr = '.';
      ptr++;
    }
    while (IS_WHITESP (ch)) {
      tmp++;
      ch = *tmp;
    }
    while (IS_DIGIT (ch)) {
      *ptr = ch;
      ptr++;
      tmp++;
      ch = *tmp;
    }
    while (IS_WHITESP (ch)) {
      tmp++;
      ch = *tmp;
    }
    if (ch == '\0' || ch == '(' || ch == ')' || ch == ',' || ch == '.' || IS_WHITESP (ch) || IS_DIGIT (ch)) {
    } else {
      MemFree (str);
      /* lower case the contents */
      ptr = gbq->val;
      ch = *ptr;
      while (ch != '\0') {
        if (IS_UPPER (ch)) {
          *ptr = TO_LOWER (ch);
        }
        ptr++;
        ch = *ptr;
      }
      return;
    }
  }
  *ptr = '\0';
  gbq->val = MemFree (gbq->val);
  gbq->val = str;
  /* and lower case the contents */
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_UPPER (ch)) {
      *ptr = TO_LOWER (ch);
    }
    ptr++;
    ch = *ptr;
  }
}

static void CleanupRptUnitSeq (GBQualPtr gbq)

{
  Char     ch;
  CharPtr  ptr;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;

  /* do not clean if val contains non-sequence characters */
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (StringChr ("ACGTUacgtu", ch) == NULL) return;
    ptr++;
    ch = *ptr;
  }

  /* lower case, and convert U to T */
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_UPPER (ch)) {
      ch = TO_LOWER (ch);
      *ptr = ch;
    }
    if (ch == 'u') {
      ch = 't';
      *ptr = ch;
    }
    ptr++;
    ch = *ptr;
  }
}

static void CleanupRptUnitRange (GBQualPtr gbq)

{
  Char     ch;
  Int2     dashes = 0;
  Int2     dots = 0;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '-') {
      dashes++;
    } else if (ch == '.') {
      dots++;
    } else if (IS_DIGIT (ch)) {
      /* okay */
    } else return;
    ptr++;
    ch = *ptr;
  }

  if (dashes > 0 && dots == 0) {
    len = StringLen (gbq->val + dashes);
    str = (CharPtr) MemNew (sizeof (Char) * (len + 5));
    tmp = str;
    ptr = gbq->val;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '-') {
        *tmp = '.';
        tmp++;
        *tmp = '.';
        tmp++;
      } else {
        *tmp = ch;
        tmp++;
      }
      ptr++;
      ch = *ptr;
    }
    gbq->val = MemFree (gbq->val);
    gbq->val = str;
  }
}

static void CleanupReplace (GBQualPtr gbq)

{
  Char     ch;
  CharPtr  ptr;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (StringChr ("ACGTUacgtu", ch) == NULL) return;
    ptr++;
    ch = *ptr;
  }
  /* lower case, and convert U to T */
  ptr = gbq->val;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_UPPER (ch)) {
      ch = TO_LOWER (ch);
      *ptr = ch;
    }
    if (ch == 'u') {
      ch = 't';
      *ptr = ch;
    }
    ptr++;
    ch = *ptr;
  }
}

static CharPtr evCategoryPfx [] = {
  "",
  "COORDINATES: ",
  "DESCRIPTION: ",
  "EXISTENCE: ",
  NULL
};

static void CleanupInference (GBQualPtr gbq)

{
  Char     ch;
  CharPtr  colon;
  CharPtr  dst;
  Int2     j;
  size_t   len;
  CharPtr  ptr;
  CharPtr  skip;
  CharPtr  space;
  CharPtr  str;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;

  str = gbq->val;
  space = NULL;
  colon = NULL;

  skip = NULL;
  for (j = 0; evCategoryPfx [j] != NULL; j++) {
    len = StringLen (evCategoryPfx [j]);
    if (StringNICmp (str, evCategoryPfx [j], len) != 0) continue;
    skip = str + len;
  }
  if (skip != NULL) {
    str = skip;
  }

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    *dst = ch;
    if (ch == ' ') {
      if (space == NULL) {
        space = dst;
      }
    } else if (ch == ':') {
      if (space != NULL) {
        dst = space;
        *dst = ch;
      }
      space = NULL;
      colon = dst;
    } else {
      if (space != NULL && colon != NULL) {
        colon++;
        dst = colon;
        *dst = ch;
      }
      space = NULL;
      colon = NULL;
    }
    dst++;
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    *dst = ch;
    if ((ch == ':' || ch == ',') && *(ptr + 1) == '?' && *(ptr + 2) == '|') {
      ptr += 2;
    }
    dst++;
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';
}

static CharPtr evCategoryNoSpace [] = {
  "",
  "COORDINATES:",
  "DESCRIPTION:",
  "EXISTENCE:",
  NULL
};

static void RepairInference (GBQualPtr gbq)

{
  Int2     j;
  size_t   len;
  CharPtr  ptr;
  CharPtr  skip;
  CharPtr  str;

  if (gbq == NULL) return;
  if (StringHasNoText (gbq->val)) return;

  str = gbq->val;
  for (j = 0; evCategoryNoSpace [j] != NULL; j++) {
    len = StringLen (evCategoryNoSpace [j]);
    if (StringNICmp (str, evCategoryNoSpace [j], len) != 0) continue;
    if (StringNICmp (str, evCategoryPfx [j], len + 1) == 0) continue;
    /* need to repair */
    skip = str + len;
    ptr = MemNew (StringLen (skip) + 20);
    if (ptr == NULL) return;
    StringCpy (ptr, evCategoryPfx [j]);
    StringCat (ptr, skip);
    gbq->val = MemFree (gbq->val);
    gbq->val = ptr;
    return;
  }
}

static void CleanupConsSplice (GBQualPtr gbq)

{
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;

  if (StringNICmp (gbq->val, "(5'site:", 8) != 0) return;
  ptr = StringStr (gbq->val, ",3'site:");
  if (ptr == NULL) return;
  len = StringLen (gbq->val) + 5;
  str = (CharPtr) MemNew (len);
  if (str == NULL) return;
  *ptr = '\0';
  ptr++;
  StringCpy (str, gbq->val);
  StringCat (str, ", ");
  StringCat (str, ptr);
  gbq->val = MemFree (gbq->val);
  gbq->val = str;
}

static Boolean ExpandParenGroup (GBQualPtr headgbq)

{
  Char       ch;
  GBQualPtr  lastgbq;
  size_t     len;
  Int2       nesting;
  GBQualPtr  newgbq;
  GBQualPtr  nextqual;
  CharPtr    ptr;
  CharPtr    str;
  CharPtr    tmp;

  nextqual = headgbq->next;
  lastgbq = headgbq;
  ptr = headgbq->val;
  tmp = StringSave (ptr + 1);
  len = StringLen (tmp);
  if (len > 0 && tmp [len - 1] == ')') {
    tmp [len - 1] = '\0';
  }
  str = tmp;
  nesting = 0;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '(') {
      nesting++;
    } else if (ch == ')') {
      nesting--;
      if (nesting < 0) {
        MemFree (tmp);
        return FALSE;
      }
    } else if (ch == ',') {
      if (nesting < 0) {
        MemFree (tmp);
        return FALSE;
      }
    }
    ptr++;
    ch = *ptr;
  }
  while (! StringHasNoText (str)) {
    ptr = StringChr (str, ',');
    if (ptr == NULL) {
      ptr = StringRChr (str, ')');
    }
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    newgbq = GBQualNew ();
    if (newgbq != NULL) {
      newgbq->qual = StringSave (headgbq->qual);
      newgbq->val = StringSave (str);
      newgbq->next = nextqual;
      lastgbq->next = newgbq;
      lastgbq = newgbq;
    }
    str = ptr;
  }
  MemFree (tmp);
  return TRUE;
}

static Boolean IsBaseRange (CharPtr str)

{
  CharPtr   ptr;
  Char      tmp [32];
  long int  val;

  if (StringLen (str) > 25) return FALSE;
  StringNCpy_0 (tmp, str, sizeof (tmp));
  ptr = StringStr (tmp, "..");
  if (ptr == NULL) return FALSE;
  *ptr = '\0';
  if (StringHasNoText (tmp)) return FALSE;
  if (sscanf (tmp, "%ld", &val) != 1 || val < 1) return FALSE;
  ptr += 2;
  if (StringHasNoText (ptr)) return FALSE;
  if (sscanf (ptr, "%ld", &val) != 1 || val < 1) return FALSE;
  return TRUE;
}

static void ModernizeFeatureGBQuals (SeqFeatPtr sfp)

{
  GBQualPtr       gbq;
  size_t          len;
  GBQualPtr       nextqual;
  GBQualPtr PNTR  prevqual;
  CharPtr         str;
  Boolean         unlink;

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
    if (StringICmp (gbq->qual, "rpt_unit_seq") == 0) {
      str = gbq->val;
      len = StringLen (str);
      if (len > 1 && *str == '{' && str [len - 1] == '}') {
        *str = '(';
        str [len - 1] = ')';
      }
      if (len > 1 && *str == '(' && str [len - 1] == ')' /* && StringChr (str + 1, '(') == NULL */) {
        if (ExpandParenGroup (gbq)) {
          nextqual = gbq->next;
          /* individual parsed out (xxx,xxx) qualifiers will be processed next, now get rid of original */
          unlink = TRUE;
        } else {
          unlink = FALSE;
        }
      } else {
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "rpt_type") == 0 ||
        StringICmp (gbq->qual, "rpt_unit") == 0 ||
        StringICmp (gbq->qual, "rpt_unit_range") == 0 ||
        StringICmp (gbq->qual, "rpt_unit_seq") == 0 ||
        StringICmp (gbq->qual, "replace") == 0 ||
        StringICmp (gbq->qual, "compare") == 0 ||
        StringICmp (gbq->qual, "old_locus_tag") == 0 ||
        StringICmp (gbq->qual, "usedin") == 0) {
      str = gbq->val;
      len = StringLen (str);
      if (len > 1 && *str == '{' && str [len - 1] == '}') {
        *str = '(';
        str [len - 1] = ')';
      }
      if (len > 1 && *str == '(' && str [len - 1] == ')' && StringChr (str + 1, '(') == NULL) {
        if (ExpandParenGroup (gbq)) {
          nextqual = gbq->next;
          /* individual parsed out (xxx,xxx) qualifiers will be processed next, now get rid of original */
          unlink = TRUE;
        } else {
          unlink = FALSE;
        }
      } else {
        unlink = FALSE;
      }
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


static void MendSatelliteQualifier (CharPtr PNTR satellite)
{
  Int4 microsatellite_len = StringLen ("microsatellite");
  Int4 minisatellite_len = StringLen ("minisatellite");
  Int4 satellite_len = StringLen ("satellite");
  Int4 type_len = 0;
  CharPtr new_qual, colon, src, dst;

  if (satellite == NULL || StringHasNoText (*satellite)) {
    return;
  }

  if (StringNCmp (*satellite, "microsatellite", microsatellite_len) == 0) {
    type_len = microsatellite_len;
  } else if (StringNCmp (*satellite, "minisatellite", minisatellite_len) == 0) {
    type_len = minisatellite_len;
  } else if (StringNCmp (*satellite, "satellite", satellite_len) == 0) {
    type_len = satellite_len;
  }
  
  if (type_len == 0) {
    new_qual = (CharPtr) MemNew (sizeof (Char) * (StringLen (*satellite) + satellite_len + 3));
    sprintf (new_qual, "satellite:%s", *satellite);
    *satellite = MemFree (*satellite);
    *satellite = new_qual;
  } else if (*(*satellite + type_len) == ' ') {
    *(*satellite + type_len) = ':';
  }

  /* remove spaces after colon */
  colon = StringChr (*satellite, ':');
  if (colon != NULL) {
    src = colon + 1;
    dst = colon + 1;
    while (*src == ' ') {
      src++;
    }
    while (*src != 0) {
      *dst = *src;
      dst++;
      src++;
    }
    *dst = 0;
  }
}


static void CleanupFeatureGBQuals (SeqFeatPtr sfp, Boolean isEmblOrDdbj)

{
  ValNodePtr      afterMe = NULL;
  Boolean         all_digits;
  Char            ch;
  DbtagPtr        db;
  GBQualPtr       gbq;
  GeneRefPtr      grp;
  ImpFeatPtr      ifp;
  size_t          len;
  GBQualPtr       nextqual;
  ObjectIdPtr     oip;
  GBQualPtr PNTR  prevqual;
  CharPtr         ptr;
  GBQualPtr       rpt_unit_range = NULL;
  GBQualPtr       rpt_unit_seq = NULL;
  CharPtr         str;
  CharPtr         tag;
  Boolean         unlink;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    CleanVisString (&(gbq->qual));
    CleanVisStringAndCompress (&(gbq->val));
    if (gbq->qual == NULL) {
      gbq->qual = StringSave ("");
    }
    if (StringIsJustQuotes (gbq->val)) {
      gbq->val = MemFree (gbq->val);
    }
    if (gbq->val == NULL) {
      gbq->val = StringSave ("");
    }
    if (StringICmp (gbq->qual, "replace") == 0) {
      if (sfp->data.choice == SEQFEAT_IMP) {
        ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
        if (ifp != NULL) {
          if (StringICmp (ifp->key, "variation") == 0 && gbq->val != NULL) {
            ptr = gbq->val;
            ch = *ptr;
            while (ch != '\0') {
              *ptr = TO_LOWER (ch);
              ptr++;
              ch = *ptr;
            }
          }
        }
      }
    }
    nextqual = gbq->next;
    unlink = TRUE;
    if (StringICmp (gbq->qual, "partial") == 0) {
      sfp->partial = TRUE;
    } else if (StringICmp (gbq->qual, "evidence") == 0) {
      /*
      if (StringICmp (gbq->val, "experimental") == 0) {
        if (sfp->exp_ev != 2) {
          sfp->exp_ev = 1;
        }
      } else if (StringICmp (gbq->val, "not_experimental") == 0) {
        sfp->exp_ev = 2;
      }
      */
    } else if (StringICmp (gbq->qual, "exception") == 0) {
      sfp->excpt = TRUE;
      if (! HasNoText (gbq->val)) {
        if (StringICmp (gbq->val, "TRUE") != 0) {
          if (sfp->except_text == NULL) {
            sfp->except_text = StringSaveNoNull (gbq->val);
          }
        }
      }
    } else if (StringICmp (gbq->qual, "note") == 0 ||
               StringICmp (gbq->qual, "notes") == 0 ||
               StringICmp (gbq->qual, "comment") == 0) {
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
    } else if (StringICmp (gbq->qual, "label") == 0) {
      if (StringICmp (gbq->val, FindKeyFromFeatDefType (sfp->idx.subtype, FALSE)) == 0) {
        /* skip label that is simply the feature key */
      } else if (sfp->comment == NULL || StringISearch (sfp->comment, gbq->qual) == NULL) {
        /* if label is not already in comment, append */
        len = StringLen (sfp->comment) + StringLen (gbq->val) + StringLen ("label: ") + 5;
        str = MemNew (sizeof (Char) * len);
        if (sfp->comment == NULL) {
          StringCpy (str, "label: ");
          StringCat (str, gbq->val);
          sfp->comment = str;
        } else {
          StringCpy (str, sfp->comment);
          StringCat (str, "; ");
          StringCat (str, "label: ");
          StringCat (str, gbq->val);
          sfp->comment = MemFree (sfp->comment);
          sfp->comment = str;
        }
      }
    } else if (StringICmp (gbq->qual, "db_xref") == 0) {
      tag = gbq->val;
      ptr = StringChr (tag, ':');
      if (ptr != NULL) {
        vnp = ValNodeNew (NULL);
        db = DbtagNew ();
        vnp->data.ptrvalue = db;
        *ptr = '\0';
        ptr++;
        db->db = StringSave (tag);
        oip = ObjectIdNew ();
        oip->str = StringSave (ptr);
        db->tag = oip;
        vnp->next = sfp->dbxref;
        sfp->dbxref = vnp;
      } else {
        /*
        db->db = StringSave ("?");
        oip = ObjectIdNew ();
        oip->str = StringSave (tag);
        db->tag = oip;
        vnp->next = sfp->dbxref;
        sfp->dbxref = vnp;
        */
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "gdb_xref") == 0) {
      vnp = ValNodeNew (NULL);
      db = DbtagNew ();
      vnp->data.ptrvalue = db;
      db->db = StringSave ("GDB");
      oip = ObjectIdNew ();
      oip->str = StringSave (gbq->val);
      db->tag = oip;
      vnp->next = sfp->dbxref;
      sfp->dbxref = vnp;
    } else if (StringICmp (gbq->qual, "cons_splice") == 0) {
      /*
      CleanupConsSplice (gbq);
      unlink = FALSE;
      */
    } else if (StringICmp (gbq->qual, "replace") == 0) {
      CleanupReplace (gbq);
      unlink = FALSE;
    } else if (StringICmp (gbq->qual, "rpt_unit_seq") == 0) {
      if (IsBaseRange (gbq->val)) {
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("rpt_unit_range");
        CleanupRptUnitRange (gbq);
      } else {
        CleanupRptUnitSeq (gbq);
      }
      unlink = FALSE;
    } else if (StringICmp (gbq->qual, "rpt_unit_range") == 0) {
      if (! IsBaseRange (gbq->val)) {
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("rpt_unit_seq");
        CleanupRptUnitSeq (gbq);
      } else {
        CleanupRptUnitRange (gbq);
      }
      unlink = FALSE;
    } else if (sfp->data.choice == SEQFEAT_GENE && HandledGBQualOnGene (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_CDREGION && HandledGBQualOnCDS (sfp, gbq, &afterMe)) {
    } else if (sfp->data.choice == SEQFEAT_RNA && HandledGBQualOnRNA (sfp, gbq, isEmblOrDdbj)) {
    } else if (sfp->data.choice == SEQFEAT_PROT && HandledGBQualOnProt (sfp, gbq)) {
    } else if (sfp->data.choice == SEQFEAT_IMP && HandledGBQualOnImp (sfp, gbq)) {
    } else if (StringICmp (gbq->qual, "rpt_unit") == 0) {
      if (IsBaseRange (gbq->val)) {
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("rpt_unit_range");
        unlink = FALSE;
      } else {
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("rpt_unit_seq");
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "EC_number") == 0) {
      CleanupECNumber (gbq->val);
      unlink = FALSE;
    } else if (StringICmp (gbq->qual, "pseudo") == 0) {
      sfp->pseudo = TRUE;
    } else if (StringICmp (gbq->qual, "pseudogene") == 0) {
      str = gbq->val;
      if (StringICmp (str, "processed") == 0 ||
          StringICmp (str, "unprocessed") == 0 ||
          StringICmp (str, "unitary") == 0 ||
          StringICmp (str, "allelic") == 0 ||
          StringICmp (str, "unknown") == 0) {
        sfp->pseudo = TRUE;
        ptr = str;
        ch = *ptr;
        while (ch != '\0') {
          if (IS_UPPER (ch)) {
            *ptr = TO_LOWER (ch);
          }
          ptr++;
          ch = *ptr;
        }
      }
      unlink = FALSE;
    } else if (StringICmp (gbq->qual, "ribosomal_slippage") == 0 ||
               StringICmp (gbq->qual, "ribosomal-slippage") == 0 ||
               StringICmp (gbq->qual, "ribosomal slippage") == 0) {
      sfp->excpt = TRUE;
      if (HasNoText (gbq->val)) {
        if (sfp->except_text == NULL) {
          sfp->except_text = StringSaveNoNull ("ribosomal slippage");
        }
      }
    } else if (StringICmp (gbq->qual, "trans_splicing") == 0 ||
               StringICmp (gbq->qual, "trans-splicing") == 0 ||
               StringICmp (gbq->qual, "trans splicing") == 0) {
      sfp->excpt = TRUE;
      if (HasNoText (gbq->val)) {
        if (sfp->except_text == NULL) {
          sfp->except_text = StringSaveNoNull ("trans-splicing");
        }
      }
    } else if (StringICmp (gbq->qual, "artificial_location") == 0 ||
               StringICmp (gbq->qual, "artificial-location") == 0 ||
               StringICmp (gbq->qual, "artificial location") == 0) {
      sfp->excpt = TRUE;
      if (HasNoText (gbq->val)) {
        if (sfp->except_text == NULL) {
          sfp->except_text = StringSaveNoNull ("artificial location");
        }
      }
    } else if (StringICmp (gbq->qual, "gene") == 0 && (! StringHasNoText (gbq->val))) {
      grp = GeneRefNew ();
      grp->locus = StringSave (gbq->val);
      xref = SeqFeatXrefNew ();
      xref->data.choice = SEQFEAT_GENE;
      xref->data.value.ptrvalue = (Pointer) grp;
      xref->specialCleanupFlag = TRUE; /* flag to test for overlapping gene later */
      xref->next = sfp->xref;
      sfp->xref = xref;
    } else if (sfp->data.choice != SEQFEAT_CDREGION && StringICmp (gbq->qual, "codon_start") == 0) {
      /* not legal on anything but CDS, so remove it */
    } else if (StringICmp (gbq->qual, "experiment") == 0 &&
               StringICmp (gbq->val, "experimental evidence, no additional details recorded") == 0) {
      /* remove default experiment string if instantiated */
    } else if (StringICmp (gbq->qual, "inference") == 0) {
      if (StringICmp (gbq->val, "non-experimental evidence, no additional details recorded") == 0) {
        /* remove default inference string if instantiated */
      } else {
        CleanupInference (gbq);
        RepairInference (gbq);
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "transposon") == 0) {
      if (StringICmp (gbq->val, "class I integron") == 0 ||
          StringICmp (gbq->val, "class II integron") == 0 ||
          StringICmp (gbq->val, "class III integron") == 0 ||
          StringICmp (gbq->val, "class 1 integron") == 0 ||
          StringICmp (gbq->val, "class 2 integron") == 0 ||
          StringICmp (gbq->val, "class 3 integron") == 0) {
        len = StringLen ("integron") + StringLen (gbq->val) + 5;
        str = MemNew (sizeof (Char) * len);
        StringCpy (str, "integron");
        StringCat (str, ":");
        ptr = StringStr (gbq->val, " integron");
        if (ptr != NULL) {
          *ptr = '\0';
        }
        StringCat (str, gbq->val);
        gbq->val = MemFree (gbq->val);
        gbq->val = str;
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("mobile_element");
        unlink = FALSE;
      } else {
        len = StringLen ("transposon") + StringLen (gbq->val) + 5;
        str = MemNew (sizeof (Char) * len);
        StringCpy (str, "transposon");
        StringCat (str, ":");
        StringCat (str, gbq->val);
        gbq->val = MemFree (gbq->val);
        gbq->val = str;
        gbq->qual = MemFree (gbq->qual);
        gbq->qual = StringSave ("mobile_element");
        unlink = FALSE;
      }
    } else if (StringICmp (gbq->qual, "insertion_seq") == 0) {
      len = StringLen ("insertion sequence") + StringLen (gbq->val) + 5;
      str = MemNew (sizeof (Char) * len);
      StringCpy (str, "insertion sequence");
      StringCat (str, ":");
      StringCat (str, gbq->val);
      gbq->val = MemFree (gbq->val);
      gbq->val = str;
      gbq->qual = MemFree (gbq->qual);
      gbq->qual = StringSave ("mobile_element");
      unlink = FALSE;
    } else if (StringCmp (gbq->qual, "satellite") == 0) {
      MendSatelliteQualifier(&(gbq->val));
      unlink = FALSE;
    } else {
      unlink = FALSE;
    }
    
    if (StringICmp (gbq->qual, "mobile_element") == 0) {
      if (sfp->data.choice == SEQFEAT_IMP) {
        ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
        if (ifp != NULL) {
          if (StringICmp (ifp->key, "repeat_region") == 0 && gbq->val != NULL) {
            gbq->qual = MemFree (gbq->qual);
            gbq->qual = StringSave ("mobile_element_type");
            ifp->key = MemFree (ifp->key);
            ifp->key = StringSave ("mobile_element");
            sfp->idx.subtype = FEATDEF_mobile_element;
          }
        }
      }
    }
    if (StringICmp (gbq->qual, "mobile_element") == 0) {
      gbq->qual = MemFree (gbq->qual);
      gbq->qual = StringSave ("mobile_element_type");
    }
    if (StringICmp (gbq->qual, "mobile_element_type") == 0) {
      if (StringStr (gbq->val, " :") != NULL || StringStr (gbq->val, ": ") != NULL) {
        len = StringLen (gbq->val) + 5;
        ptr = StringChr (gbq->val, ':');
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
          TrimSpacesAroundString (gbq->val);
          TrimSpacesAroundString (ptr);
          str = MemNew (sizeof (Char) * len);
          StringCpy (str, gbq->val);
          StringCat (str, ":");
          StringCat (str, ptr);
          gbq->val = MemFree (gbq->val);
          gbq->val = str;
        }
      }
    }

    if (StringICmp (gbq->qual, "estimated_length") == 0) {
      all_digits = TRUE;
      ptr = gbq->val;
      if (ptr != NULL) {
        ch = *ptr;
        while (ch != '\0') {
          if (! IS_DIGIT (ch)) {
            all_digits = FALSE;
          }
          ptr++;
          ch = *ptr;
        }
      }
      if (! all_digits) {
        if (StringICmp (gbq->val, "unknown") != 0) {
          MemFree (gbq->val);
          gbq->val = StringSave ("unknown");
        }
      }
    }

    if (sfp->data.choice == SEQFEAT_IMP) {
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (ifp != NULL) {
        if (StringICmp (ifp->key, "conflict") == 0 ) {
          ifp->key = MemFree (ifp->key);
          ifp->key = StringSave ("misc_difference");
          sfp->idx.subtype = FEATDEF_misc_difference;
          len = StringLen (sfp->comment) + StringLen ("conflict") + 5;
          str = MemNew (sizeof (Char) * len);
          if (sfp->comment == NULL) {
            StringCpy (str, "conflict");
            sfp->comment = str;
          } else {
            StringCpy (str, "conflict; ");
            StringCat (str, sfp->comment);
            sfp->comment = MemFree (sfp->comment);
            sfp->comment = str;
          }
        }
      }
    }

    if (rpt_unit_seq != NULL) {
      CleanupRptUnit (rpt_unit_seq);
    }
    if (rpt_unit_range != NULL) {
      CleanupRptUnit (rpt_unit_range);
    }

    if (StringHasNoText (gbq->qual) && StringHasNoText (gbq->val)) {
      unlink = TRUE;
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

static int LIBCALLBACK SortByGBQualKeyAndVal (VoidPtr ptr1, VoidPtr ptr2)

{
  int        compare;
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
  compare = StringICmp (str1, str2);
  if (compare != 0) return compare;
  str1 = (CharPtr) gbq1->val;
  str2 = (CharPtr) gbq2->val;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  return compare;
}

static Boolean GBQualsAlreadyInOrder (GBQualPtr list)

{
  int        compare;
  GBQualPtr  curr;
  GBQualPtr  next;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    compare = StringICmp (curr->qual, next->qual);
    if (compare > 0) return FALSE;
    if (compare == 0) {
      compare = StringICmp (curr->val, next->val);
      if (compare > 0) return FALSE;
    }
    curr = next;
    next = curr->next;
  }
  return TRUE;
}

NLM_EXTERN GBQualPtr SortFeatureGBQuals (GBQualPtr list)

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

  StableMergeSort (head, count, sizeof (GBQualPtr), SortByGBQualKeyAndVal);

  for (i = 0; i < count; i++) {
    gbq = head [i];
    gbq->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}

NLM_EXTERN void CleanupDuplicateGBQuals (GBQualPtr PNTR prevgbq)

{
  GBQualPtr  gbq;
  GBQualPtr  last = NULL;
  GBQualPtr  next;
  Boolean    unlink;

  if (prevgbq == NULL) return;
  gbq = *prevgbq;
  while (gbq != NULL) {
    next = gbq->next;
    unlink = FALSE;
    if (last != NULL) {
      if (StringICmp (last->qual, gbq->qual) == 0 &&
          StringICmp (last->val, gbq->val) == 0) {
        unlink = TRUE;
      }
    } else {
      last = gbq;
    }
    if (unlink) {
      *prevgbq = gbq->next;
      gbq->next = NULL;
      GBQualFree (gbq);
    } else {
      last = gbq;
      prevgbq = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = next;
  }
}

/* this identifies gbquals that should have been placed into special fields */

#define NUM_ILLEGAL_QUALS 14

/* StringICmp use of TO_UPPER means translation should go before transl_XXX */

static CharPtr illegalGbqualList [NUM_ILLEGAL_QUALS] = {
  "anticodon",
  "citation",
  "codon_start",
  "db_xref",
  "evidence",
  "exception",
  "gene",
  "note",
  "protein_id",
  "pseudo",
  "transcript_id",
  "translation",
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

static Boolean IsSubString (CharPtr str1, CharPtr str2)

{
  Char    ch;
  size_t  len1, len2;

  len1 = StringLen (str1);
  len2 = StringLen (str2);
  if (len1 >= len2) return FALSE;
  if (StringNICmp (str1, str2, len1) != 0) return FALSE;
  ch = str2 [len1];
  if (IS_ALPHANUM (ch)) return FALSE;
  return TRUE;
}

static int LIBCALLBACK SortByOrgModSubtype (VoidPtr ptr1, VoidPtr ptr2)

{
  int        compare;
  OrgModPtr  omp1;
  OrgModPtr  omp2;
  CharPtr    str1;
  CharPtr    str2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  omp1 = *((OrgModPtr PNTR) ptr1);
  omp2 = *((OrgModPtr PNTR) ptr2);
  if (omp1 == NULL || omp2 == NULL) return 0;
  if (omp1->subtype > omp2->subtype) {
    return 1;
  } else if (omp1->subtype < omp2->subtype) {
    return -1;
  }
  str1 = (CharPtr) omp1->subname;
  str2 = (CharPtr) omp2->subname;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  return compare;
}

static Boolean OrgModsAlreadyInOrder (OrgModPtr list)

{
  int        compare;
  OrgModPtr  curr;
  OrgModPtr  next;
  CharPtr    str1;
  CharPtr    str2;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    if (curr->subtype > next->subtype) return FALSE;
    str1 = (CharPtr) curr->subname;
    str2 = (CharPtr) next->subname;
    compare = StringICmp (str1, str2);
    if (compare > 0) return FALSE;
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

  StableMergeSort (head, count, sizeof (OrgModPtr), SortByOrgModSubtype);

  for (i = 0; i < count; i++) {
    omp = head [i];
    omp->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}


static void RemoveSpaceBeforeAndAfterColon (CharPtr str)
{
  CharPtr pColon, cp, src, dst;

  if (StringHasNoText (str)) {
    return;
  }

  pColon = StringChr (str, ':');
  while (pColon != NULL) {
    cp = pColon - 1;
    while (cp > str && isspace (*cp)) {
      cp--;
    }
    if (cp < str || !isspace (*cp)) {
      cp++;
    }
    *cp = ':';
    dst = cp + 1;
    cp = pColon + 1;
    while (isspace (*cp)) {
      cp++;
    }
    src = cp;
    pColon = dst - 1;
    if (src != dst) {
      while (*src != 0) {
        *dst = *src;
        dst++; src++;
      }
      *dst = 0;
    }
    pColon = StringChr (pColon + 1, ':');
  }
}

static void CorrectTildes (
  CharPtr PNTR str
)

{
#ifndef OS_MSWIN
  FindReplaceString (str, "were ~25 cm in height (~3 weeks)", "were ~~25 cm in height (~~3 weeks)", FALSE, FALSE);
  FindReplaceString (str, "generally ~3 weeks", "generally ~~3 weeks", FALSE, FALSE);
  FindReplaceString (str, "sequencing (~4 96-well plates)", "sequencing (~~4 96-well plates)", FALSE, FALSE);
  FindReplaceString (str, "size distribution (~2 kb)", "size distribution (~~2 kb)", FALSE, FALSE);
  FindReplaceString (str, "sequencing (~3 96-well plates)", "sequencing (~~3 96-well plates)", FALSE, FALSE);
  FindReplaceString (str, "vector. 1~2 ul of ligated", "vector. 1~~2 ul of ligated", FALSE, FALSE);
  /*
  FindReplaceString (str, "Lambda FLC I.~Islet cells were provided", "Lambda FLC I.~~Islet cells were provided", FALSE, FALSE);
  */
  FindReplaceString (str, "different strains~of mice", "different strains of mice", FALSE, FALSE);
  FindReplaceString (str, "oligo-dT-NotI primer~(5'-biotin", "oligo-dT-NotI primer (5'-biotin", FALSE, FALSE);
  FindReplaceString (str, "sizes of 200~800 bp were purified", "sizes of 200~~800 bp were purified", FALSE, FALSE);
  FindReplaceString (str, "Tween 20 (~50 ml per tree)", "Tween 20 (~~50 ml per tree)", FALSE, FALSE);
  FindReplaceString (str, "the SMART approach (~http://www.evrogen.com", "the SMART approach (http://www.evrogen.com", FALSE, FALSE);
  FindReplaceString (str, "the morning (~10 am) with", "the morning (~~10 am) with", FALSE, FALSE);
  FindReplaceString (str, "(host) sequences (~10%)", "(host) sequences (~~10%)", FALSE, FALSE);
  /*
  FindReplaceString (str, "unidirectionally.~ High quality", "unidirectionally. High quality", FALSE, FALSE);
  FindReplaceString (str, "onlysubmitted.~ Average", "onlysubmitted. Average", FALSE, FALSE);
  */
  FindReplaceString (str, "Plasmid; ~The F03-1270", "Plasmid; The F03-1270", FALSE, FALSE);
  FindReplaceString (str, "using STS-PCR~from Eb", "using STS-PCR from Eb", FALSE, FALSE);
  FindReplaceString (str, "specific to~the Eb", "specific to the Eb", FALSE, FALSE);
  FindReplaceString (str, "side of insert); , M.F., Lennon", "side of insert); Bonaldo, M.F., Lennon", FALSE, FALSE);
  FindReplaceString (str, "Uni-ZAP XR vector. 1~2 ul of", "Uni-ZAP XR vector. 1~~2 ul of", FALSE, FALSE);
  FindReplaceString (str, "from diploid~Secale montanum", "from diploid Secale montanum", FALSE, FALSE);
  FindReplaceString (str, "homology with~U43516,", "homology with U43516,", FALSE, FALSE);
  /*
  FindReplaceString (str, "from http//www.biobase.dk/~ddbase", "from http//www.biobase.dk/~~ddbase", FALSE, FALSE);
  */
  FindReplaceString (str, "plasmid; ~Assembled EST", "plasmid; Assembled EST", FALSE, FALSE);
  FindReplaceString (str, "databases.~Different cDNA", "databases. Different cDNA", FALSE, FALSE);
  FindReplaceString (str, "enzyme PstI.~DH5-alpha", "enzyme PstI. DH5-alpha", FALSE, FALSE);
  FindReplaceString (str, "as they~were prepared", "as they were prepared", FALSE, FALSE);
  FindReplaceString (str, "loci in~the genome", "loci in the genome", FALSE, FALSE);
  FindReplaceString (str, "P{CaSpeR}Cp1~50C (FBti0004219)", "P{CaSpeR}Cp1~~50C (FBti0004219)", FALSE, FALSE);
  FindReplaceString (str, "seedlings with 2~4 leaves", "seedlings with 2~~4 leaves", FALSE, FALSE);
  FindReplaceString (str, "tween 20 (~50mLs per tree)", "tween 20 (~~50mLs per tree)", FALSE, FALSE);
#endif
}

static void FixStrainForPrefix (OrgModPtr omp)

{
  Char        ch;
  CharPtr     cpy;
  ValNodePtr  head = NULL;
  size_t      len;
  CharPtr     pfx;
  CharPtr     sfx;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;

  if (omp == NULL || omp->subtype != ORGMOD_strain) return;
  str = omp->subname;
  if (StringHasNoText (str)) return;

  head = SplitStringAtSemicolon (str);
  if (head == NULL) return;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    TrimSpacesAroundString (str);

    pfx = NULL;
    sfx = NULL;
    if (StringNICmp (str, "ATCC", 4) == 0) {
      pfx = "ATCC";
      sfx = str + 4;
    } else if (StringNICmp (str, "DSM", 3) == 0) {
      pfx = "DSM";
      sfx = str + 3;
    }
    if (pfx == NULL || sfx == NULL) continue;

    ch = *sfx;
    if (ch == ':' || ch == '/') {
      sfx++;
    }
    cpy = StringSave (sfx);
    TrimSpacesAroundString(cpy);
    if (! StringIsAllDigits (cpy)) {
      cpy = MemFree (cpy);
      continue;
    }

    len = StringLen (pfx) + StringLen (cpy) + 3;
    tmp = (CharPtr) MemNew (len);
    if (tmp == NULL) continue;
    StringCpy (tmp, pfx);
    StringCat (tmp, " ");
    StringCat (tmp, cpy);
    vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
    vnp->data.ptrvalue = tmp;
    cpy = MemFree (cpy);
  }

  tmp = ValNodeMergeStrsEx (head, "; ");
  if (tmp == NULL) return;

  omp->subname = MemFree (omp->subname);
  omp->subname = tmp;
}

static void CleanOrgModListEx (OrgModPtr PNTR ompp, CharPtr orpcommon)

{
  Char            ch;
  OrgModPtr       last = NULL;
  OrgModPtr       next;
  OrgModPtr       omp;
  OrgModPtr       omp_anamorph, omp_gb_anamorph, omp_other;
  OrgModPtr PNTR  prev;
  CharPtr         ptr;
  Boolean         redund;
  CharPtr         str;
  CharPtr         tmp;
  Boolean         unlink;

  if (ompp == NULL) return;
  prev = ompp;
  omp = *ompp;
  while (omp != NULL) {
    next = omp->next;
    unlink= FALSE;
    CleanVisStringAndCompress (&(omp->subname));
    TrimSpacesAndJunkFromEnds (omp->subname, FALSE);
    RemoveFlankingQuotes (&(omp->subname));
    CleanVisStringAndCompress (&(omp->attrib));
    if (omp->subtype == ORGMOD_other && StringDoesHaveText (omp->subname)) {
      CorrectTildes (&(omp->subname));
    }
    if (omp->subtype == ORGMOD_common && StringICmp (omp->subname, orpcommon) == 0) {
      /*
      unlink = TRUE;
      */
    } else if (last != NULL) {
      if (HasNoText (omp->subname)) {
        unlink = TRUE;
      } else if ((last->subtype == omp->subtype &&
                 StringICmp (last->subname, omp->subname) == 0) ||
                 (last->subtype == omp->subtype &&
                 last->subtype == ORGMOD_other &&
                  StringStr (last->subname, omp->subname) != NULL)) {
        unlink = TRUE;
      } else if (last->subtype == omp->subtype &&
                 last->subtype == ORGMOD_other &&
                 IsSubString (last->subname, omp->subname)) {
        last->subname = MemFree (last->subname);
        last->subname = omp->subname;
        omp->subname = NULL;
        unlink = TRUE;
      }
    } else if (HasNoText (omp->subname) ||
               StringCmp (omp->subname, ")") == 0 ||
               StringCmp (omp->subname, "(") == 0) {
      unlink = TRUE;
    } else {
      last = omp;
    }
    if (unlink) {
      *prev = omp->next;
      omp->next = NULL;
      OrgModFree (omp);
    } else {
      last = omp;
      prev = &(omp->next);
    }
    omp = next;
  }


  for (omp = *ompp; omp != NULL; omp = omp->next) {
    if (omp->subtype != ORGMOD_specimen_voucher &&
        omp->subtype != ORGMOD_culture_collection &&
        omp->subtype != ORGMOD_bio_material) continue;
    if (StringHasNoText (omp->subname)) continue;
    RemoveSpaceBeforeAndAfterColon (omp->subname);
    ptr = StringStr (omp->subname, "::");
    if (ptr == NULL) continue;
    ptr++;
    tmp = ptr;
    tmp++;
    ch = *tmp;
    while (ch != '\0') {
      *ptr = ch;
      ptr++;
      tmp++;
      ch = *tmp;
    }
    *ptr = '\0';
  }

  omp_anamorph = NULL;
  omp_gb_anamorph = NULL;
  omp_other = NULL;
  redund = FALSE;

  for (omp = *ompp; omp != NULL; omp = omp->next) {
    if (omp->subtype == ORGMOD_anamorph) {
      omp_anamorph = omp;
    } else if (omp->subtype == ORGMOD_gb_anamorph) {
      omp_gb_anamorph = omp;
    } else if (omp->subtype == ORGMOD_other) {
      omp_other = omp;
    } else if (omp->subtype == ORGMOD_nat_host) {
      if (StringICmp (omp->subname, "human") == 0) {
        omp->subname = MemFree (omp->subname);
        omp->subname = StringSave ("Homo sapiens");
      }
    } else if (omp->subtype == ORGMOD_strain) {
      FixStrainForPrefix (omp);
    }
  }
  if (omp_other != NULL && StringNICmp (omp_other->subname, "anamorph:", 9) == 0) {
    ptr = omp_other->subname + 9;
    ch = *ptr;
    while (ch == ' ') {
      ptr++;
      ch = *ptr;
    }
    if (omp_anamorph != NULL) {
      str = omp_anamorph->subname;
      if (StringCmp (ptr, str) == 0) {
        redund = TRUE;
      }
    } else if (omp_gb_anamorph != NULL) {
      str = omp_gb_anamorph->subname;
      if (StringCmp (ptr, str) == 0) {
        redund = TRUE;
      }
    }
  }
  if (redund) {
    prev = ompp;
    omp = *ompp;
    while (omp != NULL) {
      next = omp->next;
      unlink= FALSE;
      if (omp == omp_other) {
        unlink= TRUE;
      }
      if (unlink) {
        *prev = omp->next;
        omp->next = NULL;
        OrgModFree (omp);
      } else {
        prev = &(omp->next);
      }
      omp = next;
    }
  }
}

NLM_EXTERN void CleanOrgModList (OrgModPtr PNTR ompp)

{
  CleanOrgModListEx (ompp, NULL);
}

static Boolean IsNoNameSubSource (SubSourcePtr ssp)

{
  if (ssp == NULL) return FALSE;

  return (Boolean) (ssp->subtype == SUBSRC_germline ||
                    ssp->subtype == SUBSRC_rearranged ||
                    ssp->subtype == SUBSRC_transgenic ||
                    ssp->subtype == SUBSRC_environmental_sample ||
                    ssp->subtype == SUBSRC_metagenomic);
}

static int LIBCALLBACK SortBySubSourceSubtype (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  SubSourcePtr  ssp1;
  SubSourcePtr  ssp2;
  CharPtr       str1;
  CharPtr       str2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  ssp1 = *((SubSourcePtr PNTR) ptr1);
  ssp2 = *((SubSourcePtr PNTR) ptr2);
  if (ssp1 == NULL || ssp2 == NULL) return 0;
  if (ssp1->subtype > ssp2->subtype) {
    return 1;
  } else if (ssp1->subtype < ssp2->subtype) {
    return -1;
  }
  if (IsNoNameSubSource (ssp1)) return 0;
  str1 = (CharPtr) ssp1->name;
  str2 = (CharPtr) ssp2->name;
  if (str1 == NULL || str2 == NULL) return 0;
  compare = StringICmp (str1, str2);
  return compare;
}

static Boolean SubSourceAlreadyInOrder (SubSourcePtr list)

{
  int           compare;
  SubSourcePtr  curr;
  SubSourcePtr  next;
  CharPtr       str1;
  CharPtr       str2;

  if (list == NULL || list->next == NULL) return TRUE;
  curr = list;
  next = curr->next;
  while (next != NULL) {
    if (curr->subtype > next->subtype) return FALSE;
    if (curr->subtype == next->subtype) {
      if (! IsNoNameSubSource (curr)) {
        str1 = (CharPtr) curr->name;
        str2 = (CharPtr) next->name;
        compare = StringICmp (str1, str2);
        if (compare > 0) return FALSE;
      }
    }
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

  StableMergeSort (head, count, sizeof (SubSourcePtr), SortBySubSourceSubtype);

  for (i = 0; i < count; i++) {
    ssp = head [i];
    ssp->next = head [i + 1];
  }

  list = head [0];
  MemFree (head);

  return list;
}

//LCOV_EXCL_START
static CharPtr TrimParenthesesAndCommasAroundString (CharPtr str)

{
  Uchar    ch;    /* to use 8bit characters in multibyte languages */
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

  if (StringStr (origval, newval) != NULL) return origval;
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
//LCOV_EXCL_STOP

static Uint1 LocationForPlastidText (CharPtr plastid_name)
{
  if (StringICmp (plastid_name, "chloroplast") == 0) {
    return GENOME_chloroplast;
  } else if (StringICmp (plastid_name, "chromoplast") == 0) {
    return GENOME_chromoplast;
  } else if (StringICmp (plastid_name, "kinetoplast") == 0) {
    return GENOME_kinetoplast;
  } else if (StringICmp (plastid_name, "plastid") == 0) {
    return GENOME_plastid;
  } else if (StringICmp (plastid_name, "apicoplast") == 0) {
    return GENOME_apicoplast;
  } else if (StringICmp (plastid_name, "leucoplast") == 0) {
    return GENOME_leucoplast;
  } else if (StringICmp (plastid_name, "proplastid") == 0) {
    return GENOME_proplastid;
  } else if (StringICmp (plastid_name, "chromatophore") == 0) {
    return GENOME_chromatophore;
  } else {
    return 0;
  }
}

//LCOV_EXCL_START
NLM_EXTERN void StringToLower (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    *str = TO_LOWER (ch);
    str++;
    ch = *str;
  }
}
//LCOV_EXCL_STOP


static void CleanPCRPrimerSeq (CharPtr seq)
{
  CharPtr ptr, src, dst, tmp;
  Char    ch;
  Boolean in_brackets = FALSE;
  Int4    i;

  if (StringHasNoText (seq)) {
    return;
  }

  /* upper case sequence */
  ptr = seq;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_UPPER (ch)) {
      *ptr = TO_LOWER (ch);
    }
    ptr++;
    ch = *ptr;
  }
  /* remove any spaces in sequence outisde of <modified base> */
  src = seq;
  dst = seq;
  ch = *src;
  while (ch != '\0') {
    if (ch == '<') {
      in_brackets = TRUE;
      *dst = ch;
      dst++;
    } else if (ch == '>') {
      in_brackets = FALSE;
      *dst = ch;
      dst++;
    } else if (ch != ' ') {
      *dst = ch;
      dst++;
    } else if (in_brackets) {
      *dst = ch;
      dst++;
    }
    src++;
    ch = *src;
  }
  *dst = '\0';
  /* upper case modified base <OTHER> */
  ptr = seq;
  tmp = StringStr (ptr, "<other>");
  while (tmp != NULL) {
    ptr = tmp + 7;
    for (i = 1; i < 6; i++) {
      ch = tmp [i];
      tmp [i] = TO_UPPER (ch);
    }
    tmp = StringStr (ptr, "<other>");
  }
}


static void CleanupPCRPrimers (PCRPrimerPtr PNTR pppp)

{
  PCRPrimerPtr       next;
  PCRPrimerPtr PNTR  prev;
  PCRPrimerPtr       ppp;
  PCRPrimerPtr       pr1, pr2;

  if (pppp == NULL) return;

  ppp = *pppp;
  while (ppp != NULL) {
    CleanVisString (&(ppp->seq));
    CleanPCRPrimerSeq (ppp->seq);
    CleanVisString (&(ppp->name));
    Asn2gnbkCompressSpaces (ppp->name);
    StringToLower (ppp->seq);

    ppp = ppp->next;
  }

  ppp = *pppp;
  for (pr1 = ppp; pr1 != NULL; pr1 = pr1->next)  {
    for (pr2 = pr1->next; pr2 != NULL; pr2 = pr2->next) {
      if (StringCmp (pr1->seq, pr2->seq) == 0 && StringCmp (pr1->name, pr2->name) == 0) {
        pr2->seq = MemFree (pr2->seq);
        pr2->name = MemFree (pr2->name);
      } else if (StringCmp (pr1->name, pr2->name) == 0) {
        if (StringHasNoText (pr1->seq)) {
          pr1->seq = MemFree (pr1->seq);
          pr1->seq = pr2->seq;
          pr2->seq = NULL;
        } else if (StringHasNoText (pr2->seq)) {
          pr2->seq = MemFree (pr2->seq);
          pr2->name = MemFree (pr2->name);
        }
      }
    }
  }

  prev = pppp;
  ppp = *pppp;
  while (ppp != NULL) {
    next = ppp->next;

    CleanVisString (&(ppp->seq));
    CleanPCRPrimerSeq (ppp->seq);
    CleanVisString (&(ppp->name));

    if (ppp->seq == NULL && ppp->name == NULL) {
      *prev = next;
      ppp->next = NULL;
      PCRPrimerFree (ppp);
    } else {
      StringToLower (ppp->seq);
      prev = &(ppp->next);
    }

    ppp = next;
  }

  /* fix artifact caused by fwd/rev-primer-seq starting with colon, separating name and seq */

  ppp = *pppp;
  if (ppp == NULL) return;
  next = ppp->next;
  if (next == NULL) return;
  if (next->next != NULL) return;

  if (ppp->name != NULL && ppp->seq == NULL && next->name == NULL && next->seq != NULL) {
    ppp->seq = next->seq;
    next->seq = NULL;
    ppp->next = NULL;
    PCRPrimerFree (next);
  } else if (ppp->seq != NULL && ppp->name == NULL && next->seq == NULL && next->name != NULL) {
    ppp->name = next->name;
    next->name = NULL;
    ppp->next = NULL;
    PCRPrimerFree (next);
  }
}

static Boolean PCRPrimersMatch (PCRPrimerPtr ppp1, PCRPrimerPtr ppp2)

{
  Int2          len1 = 0, len2 = 0, matches = 0;
  PCRPrimerPtr  pr1, pr2;

  if (ppp1 == NULL || ppp2 == NULL) return FALSE;

  for (pr1 = ppp1; pr1 != NULL; pr1 = pr1->next) {
    len1++;
  }
  for (pr2 = ppp2; pr2 != NULL; pr2 = pr2->next) {
    len2++;
  }
  if (len1 != len2) return FALSE;

  for (pr1 = ppp1; pr1 != NULL; pr1 = pr1->next) {
    for (pr2 = ppp2; pr2 != NULL; pr2 = pr2->next) {
      if (StringCmp (pr1->seq, pr2->seq) == 0 && StringCmp (pr1->name, pr2->name) == 0) {
        matches++;
      }
    }
  }

  if (matches == len1) return TRUE;

  return FALSE;
}

static Boolean PCRReactionSetsMatch (PCRReactionSetPtr prp1, PCRReactionSetPtr prp2)

{
  if (prp1 == NULL || prp2 == NULL) return FALSE;

  if (! PCRPrimersMatch (prp1->forward, prp2->forward)) return FALSE;
  if (! PCRPrimersMatch (prp1->reverse, prp2->reverse)) return FALSE;

  return TRUE;
}

static void CleanupPCRReactionSet (PCRReactionSetPtr PNTR prpp)

{
  PCRReactionSetPtr       curr;
  PCRReactionSetPtr       next;
  PCRReactionSetPtr PNTR  prev;
  PCRReactionSetPtr       prp;

  if (prpp == NULL) return;

  prp = *prpp;
  while (prp != NULL) {
    CleanupPCRPrimers (&(prp->forward));
    CleanupPCRPrimers (&(prp->reverse));
    prp = prp->next;
  }

  prev = prpp;
  prp = *prpp;
  while (prp != NULL) {
    next = prp->next;

    curr = next;
    while (curr != NULL) {
      if (PCRReactionSetsMatch (prp, curr)) {
        curr->forward = PCRPrimerFree (curr->forward);
        curr->reverse = PCRPrimerFree (curr->reverse);
      }
      curr = curr->next;
    }

    if (prp->forward == NULL && prp->reverse == NULL) {
      *prev = next;
      prp->next = NULL;
      PCRReactionFree (prp);
    } else {
      prev = &(prp->next);
    }

    prp = next;
  }

}

static void CleanupAltitude (SubSourcePtr ssp)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;

  if (ssp == NULL || StringHasNoText (ssp->name)) return;
  len = StringLen (ssp->name);
  if (len < 1) return;

  ptr = ssp->name;
  ch = *ptr;

  if (len > 2 && ptr [len-1] == '.') {
    ptr [len-1] = '\0';
  }

  if (ch == '+' || ch == '-') {
    ptr++;
    ch = *ptr;
  }

  if (! IS_DIGIT (ch)) return;

  ptr++;
  ch = *ptr;
  while (IS_DIGIT (ch)) {
    ptr++;
    ch = *ptr;
  }

  if (ch == '.') {
    ptr++;
    ch = *ptr;
    if (! IS_DIGIT (ch)) return;
    ptr++;
    ch = *ptr;
    while (IS_DIGIT (ch)) {
      ptr++;
      ch = *ptr;
    }
  }

  if (StringCmp (ptr, "m") == 0 ||
      StringCmp (ptr, "m.") == 0 ||
      StringCmp (ptr, " m") == 0||
      StringCmp (ptr, " meters") == 0||
      StringCmp (ptr, " metres") == 0) {
    *ptr = '\0';
    ptr = (CharPtr) MemNew (len + 5);
    if (ptr == NULL) return;
    StringCpy (ptr, ssp->name);
    StringCat (ptr, " m");
    ssp->name = MemFree (ssp->name);
    ssp->name = ptr;
  }
}

static CharPtr coll_date_month_abbrevs [12] =
{
  "-Jan-", "-Feb-", "-Mar-", "-Apr-", "-May-", "-Jun-",
  "-Jul-", "-Aug-", "-Sep-", "-Oct-", "-Nov-", "-Dec-"
};

static void CorrectMonthCapitalization (CharPtr str)

{
  Int2     i;
  Int2     j;
  CharPtr  month;
  CharPtr  ptr;

  for (i = 0; i < 12; i++) {
    month = coll_date_month_abbrevs [i];
    ptr = StringISearch (str, month);
    if (ptr == NULL) continue;
    for (j = 0; j < 5; j++) {
      ptr [j] = month [j];
    }
    return;
  }
}

typedef struct stringpair {
  CharPtr  from;
  CharPtr  to;
} StringPair, PNTR StringPairPtr;

static StringPair sex_conv[] = {
  { "asexual female",        "asexual and female"       },
  { "asexual male",          "asexual and male"         },
  { "dioecious female",      "dioecious and female"     },
  { "dioecious male",        "dioecious and male"       },
  { "f and m mixed",         "female, male, and mixed"  },
  { "f",                     "female"                   },
  { "f/m",                   "female and male"          },
  { "female,male",           "female and male"          },
  { "female/hermaphrodite",  "female and hermaphrodite" },
  { "female/male mixed",     "female, male, and mixed"  },
  { "female/male",           "female and male"          },
  { "m and f mixed",         "male, female, and mixed"  },
  { "m",                     "male"                     },
  { "m/f",                   "male and female"          },
  { "male,female",           "male and female"          },
  { "male/female mixed",     "male, female, and mixed"  },
  { "male/female",           "male and female"          },
  { "male/hermaphrodite",    "male and hermaphrodite"   },
  { "mixed female and male", "mixed, female, and male"  },
  { "mixed female/male",     "mixed, female, and male"  },
  { "mixed male and female", "mixed, male, and female"  },
  { "mixed male/female",     "mixed, male, and female"  },
  { NULL,                    NULL                       }
};

extern void CleanSubSourceList (SubSourcePtr PNTR sspp, Uint1 location)

{
  Char               ch;
  CharPtr            dst;
  Int2               i;
  Boolean            in_brackets = FALSE;
  SubSourcePtr       last = NULL;
  size_t             len;
  SubSourcePtr       next;
  SubSourcePtr PNTR  prev;
  CharPtr            ptr;
  CharPtr            src;
  SubSourcePtr       ssp;
  CharPtr            str;
  CharPtr            tmp;
  Boolean            unlink;
  /*
  FloatHi            ns, ew;
  Char               lon, lat;
  Int4               processed;
  */
  /*
  SubSourcePtr       fwd_seq = NULL, rev_seq = NULL, fwd_name = NULL, rev_name = NULL;
  size_t             len;
  */

  if (sspp == NULL) return;
  prev = sspp;
  ssp = *sspp;
  while (ssp != NULL) {
    next = ssp->next;
    unlink= FALSE;
    if (! IsNoNameSubSource (ssp)) {
      CleanVisStringAndCompress (&(ssp->name));
      TrimSpacesAndJunkFromEnds (ssp->name, FALSE);
      RemoveFlankingQuotes (&(ssp->name));
    } else /* if (StringICmp (ssp->name, "TRUE") == 0) */ {
      ssp->name = MemFree (ssp->name);
      ssp->name = StringSave ("");
    }
    if (ssp->subtype == SUBSRC_country) {
      CleanVisStringJunk (&(ssp->name));
      len = StringLen (ssp->name);
      if (len > 2) {
        str = ssp->name;
        if (str [len - 1] == ':') {
          str [len - 1] = '\0';
        }
      }
      if (StringICmp (ssp->name, "United States") == 0 ||
          StringICmp (ssp->name, "United States of America") == 0 ||
          StringICmp (ssp->name, "U.S.A.") == 0) {
        ssp->name = MemFree (ssp->name);
        ssp->name = StringSave ("USA");
      }
      if (StringNICmp (ssp->name, "United States:", 14) == 0) {
        str = ssp->name;
        str [0] = ' ';
        str [1] = ' ';
        str [2] = ' ';
        str [3] = ' ';
        str [4] = ' ';
        str [5] = ' ';
        str [6] = ' ';
        str [7] = ' ';
        str [8] = ' ';
        str [9] = ' ';
        str [10] = 'U';
        str [11] = 'S';
        str [12] = 'A';
        TrimSpacesAroundString (ssp->name);
      }
    } else if (ssp->subtype == SUBSRC_clone) {
      CleanVisStringJunk (&(ssp->name));
    } else if (ssp->subtype == SUBSRC_altitude) {
      if (ssp->name != NULL && (! AltitudeIsValid (ssp->name))) {
        CleanupAltitude (ssp);
      }
    } else if (ssp->subtype == SUBSRC_lat_lon) {
      /*
      str = ssp->name;
      if (str != NULL) {
        ptr = StringStr (str, " N, ");
        if (ptr == NULL) {
          ptr = StringStr (str, " S, ");
        }
        if (ptr != NULL) {
          ptr += 2;
          *ptr = ' ';
          Asn2gnbkCompressSpaces (str);
        }
      }
      */
      /*
      if (str != NULL && sscanf (str, "%lf %c, %lf %c%n", &ns, &lat, &ew, &lon, &processed) == 4 && processed == StringLen (str)) {
        ptr = StringChr (str, ',');
        if (ptr != NULL) {
          *ptr = ' ';
          Asn2gnbkCompressSpaces (str);
        }
      }
      */
    } else if (ssp->subtype == SUBSRC_other && StringDoesHaveText (ssp->name)) {
      CorrectTildes (&(ssp->name));
    } else if (ssp->subtype == SUBSRC_sex) {
      ptr = ssp->name;
      if (StringDoesHaveText (ptr)) {
        ch = *ptr;
        while (ch != '\0') {
          ch = TO_LOWER(ch);
          *ptr = ch;
          ptr++;
          ch = *ptr;
        }
        ptr = ssp->name;
        for (i = 0; sex_conv[i].from != NULL; i++) {
          if (StringCmp (ptr, sex_conv[i].from) == 0) {
            ssp->name = MemFree (ssp->name);
            ssp->name = StringSave (sex_conv[i].to);
            break;
          }
        }
      }
    } else if (ssp->subtype == SUBSRC_collection_date) {
      ptr = ssp->name;
      if (StringDoesHaveText (ptr)) {
        CorrectMonthCapitalization (ptr);
      }
    }
    if (ssp->subtype == SUBSRC_fwd_primer_seq ||
        ssp->subtype == SUBSRC_rev_primer_seq) {
      if (ssp->name != NULL) {
        /* upper case sequence */
        ptr = ssp->name;
        ch = *ptr;
        while (ch != '\0') {
          if (IS_UPPER (ch)) {
            *ptr = TO_LOWER (ch);
          }
          ptr++;
          ch = *ptr;
        }
        /* remove any spaces in sequence outisde of <modified base> */
        src = ssp->name;
        dst = ssp->name;
        ch = *src;
        while (ch != '\0') {
          if (ch == '<') {
            in_brackets = TRUE;
            *dst = ch;
            dst++;
          } else if (ch == '>') {
            in_brackets = FALSE;
            *dst = ch;
            dst++;
          } else if (ch != ' ') {
            *dst = ch;
            dst++;
          } else if (in_brackets) {
            *dst = ch;
            dst++;
          }
          src++;
          ch = *src;
        }
        *dst = '\0';
        /* upper case modified base <OTHER> */
        ptr = ssp->name;
        tmp = StringStr (ptr, "<other>");
        while (tmp != NULL) {
          ptr = tmp + 7;
          for (i = 1; i < 6; i++) {
            ch = tmp [i];
            tmp [i] = TO_UPPER (ch);
          }
          tmp = StringStr (ptr, "<other>");
        }
      }
    }
    /*
    if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      if (fwd_seq == NULL) {
        fwd_seq = ssp;
      } else {
        fwd_seq->name = CombineSplitQual (fwd_seq->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_rev_primer_seq) {
      if (rev_seq == NULL) {
        rev_seq = ssp;
      } else {
        rev_seq->name = CombineSplitQual (rev_seq->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_fwd_primer_name) {
      if (fwd_name == NULL) {
        fwd_name = ssp;
      } else {
        fwd_name->name = CombineSplitQual (fwd_name->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_rev_primer_name) {
      if (rev_name == NULL) {
        rev_name = ssp;
      } else {
        rev_name->name = CombineSplitQual (rev_name->name, ssp->name);
        unlink = TRUE;
      }
    }
    */
    CleanVisString (&(ssp->attrib));
    if (last != NULL) {
      if (HasNoText (ssp->name) && (! IsNoNameSubSource (ssp))) {
        unlink = TRUE;
      } else if (last->subtype == ssp->subtype &&
                 (IsNoNameSubSource (ssp) ||
                  StringICmp (last->name, ssp->name) == 0 ||
                  (last->subtype == SUBSRC_other &&
                   StringStr (last->name, ssp->name) != NULL))) {
        unlink = TRUE;
      } else if (last->subtype == ssp->subtype &&
                 last->subtype == SUBSRC_other &&
                 IsSubString (last->name, ssp->name)) {
        last->name = MemFree (last->name);
        last->name = ssp->name;
        ssp->name = NULL;
        unlink = TRUE;
      } else if (ssp->subtype == SUBSRC_plastid_name &&
                 location != 0
                 && location == LocationForPlastidText (ssp->name)) {
        unlink = TRUE;
      }
    } else if (HasNoText (ssp->name) && (! IsNoNameSubSource (ssp))) {
      unlink = TRUE;
    } else if (ssp->subtype == SUBSRC_plastid_name &&
               location != 0
               && location == LocationForPlastidText (ssp->name)) {
      unlink = TRUE;
    } else {
      last = ssp;
    }
    if (unlink) {
      *prev = ssp->next;
      ssp->next = NULL;
      SubSourceFree (ssp);
    } else {
      last = ssp;
      prev = &(ssp->next);
    }
    ssp = next;
  }
  /*
  if (fwd_seq != NULL) {
    if (StringChr (fwd_seq->name, ',') != NULL) {
      ptr = fwd_seq->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (fwd_seq->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, fwd_seq->name);
          StringCat (str, ")");
          fwd_seq->name = MemFree (fwd_seq->name);
          fwd_seq->name = str;
        }
      }
    }
  }
  if (rev_seq != NULL) {
    if (StringChr (rev_seq->name, ',') != NULL) {
      ptr = rev_seq->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (rev_seq->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, rev_seq->name);
          StringCat (str, ")");
          rev_seq->name = MemFree (rev_seq->name);
          rev_seq->name = str;
        }
      }
    }
  }
  if (fwd_name != NULL) {
    if (StringChr (fwd_name->name, ',') != NULL) {
      ptr = fwd_name->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (fwd_name->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, fwd_name->name);
          StringCat (str, ")");
          fwd_name->name = MemFree (fwd_name->name);
          fwd_name->name = str;
        }
      }
    }
  }
  if (rev_name != NULL) {
    if (StringChr (rev_name->name, ',') != NULL) {
      ptr = rev_name->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (rev_name->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, rev_name->name);
          StringCat (str, ")");
          rev_name->name = MemFree (rev_name->name);
          rev_name->name = str;
        }
      }
    }
  }
  */
}

//LCOV_EXCL_START
extern void CleanSubSourcePrimers (SubSourcePtr PNTR sspp)

{
  SubSourcePtr       fwd_seq = NULL, rev_seq = NULL, fwd_name = NULL, rev_name = NULL;
  size_t             len;
  SubSourcePtr       next;
  SubSourcePtr PNTR  prev;
  CharPtr            ptr;
  SubSourcePtr       ssp;
  CharPtr            str;
  Boolean            unlink;

  if (sspp == NULL) return;
  prev = sspp;
  ssp = *sspp;
  while (ssp != NULL) {
    next = ssp->next;
    unlink= FALSE;
    if (ssp->subtype == SUBSRC_fwd_primer_seq) {
      if (fwd_seq == NULL) {
        fwd_seq = ssp;
      } else {
        fwd_seq->name = CombineSplitQual (fwd_seq->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_rev_primer_seq) {
      if (rev_seq == NULL) {
        rev_seq = ssp;
      } else {
        rev_seq->name = CombineSplitQual (rev_seq->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_fwd_primer_name) {
      if (fwd_name == NULL) {
        fwd_name = ssp;
      } else {
        fwd_name->name = CombineSplitQual (fwd_name->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (ssp->subtype == SUBSRC_rev_primer_name) {
      if (rev_name == NULL) {
        rev_name = ssp;
      } else {
        rev_name->name = CombineSplitQual (rev_name->name, ssp->name);
        unlink = TRUE;
      }
    }
    if (unlink) {
      *prev = ssp->next;
      ssp->next = NULL;
      SubSourceFree (ssp);
    } else {
      prev = &(ssp->next);
    }
    ssp = next;
  }
  if (fwd_seq != NULL) {
    if (StringChr (fwd_seq->name, ',') != NULL) {
      ptr = fwd_seq->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (fwd_seq->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, fwd_seq->name);
          StringCat (str, ")");
          fwd_seq->name = MemFree (fwd_seq->name);
          fwd_seq->name = str;
        }
      }
    }
  }
  if (rev_seq != NULL) {
    if (StringChr (rev_seq->name, ',') != NULL) {
      ptr = rev_seq->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (rev_seq->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, rev_seq->name);
          StringCat (str, ")");
          rev_seq->name = MemFree (rev_seq->name);
          rev_seq->name = str;
        }
      }
    }
  }
  if (fwd_name != NULL) {
    if (StringChr (fwd_name->name, ',') != NULL) {
      ptr = fwd_name->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (fwd_name->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, fwd_name->name);
          StringCat (str, ")");
          fwd_name->name = MemFree (fwd_name->name);
          fwd_name->name = str;
        }
      }
    }
  }
  if (rev_name != NULL) {
    if (StringChr (rev_name->name, ',') != NULL) {
      ptr = rev_name->name;
      len = StringLen (ptr);
      if (ptr [0] != '(' || ptr [len - 1] != ')') {
        TrimParenthesesAndCommasAroundString (rev_name->name);
        str = MemNew (sizeof (Char) * (len + 4));
        if (str != NULL) {
          StringCpy (str, "(");
          StringCat (str, rev_name->name);
          StringCat (str, ")");
          rev_name->name = MemFree (rev_name->name);
          rev_name->name = str;
        }
      }
    }
  }
}
//LCOV_EXCL_STOP

static void OrpModToOrgMod (ValNodePtr PNTR vnpp, OrgModPtr PNTR ompp)

{
  Char        ch;
  ValNodePtr  next;
  Int2        numcommas;
  Int2        numspaces;
  OrgModPtr   omp;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     val;
  ValNodePtr  vnp;
  Uint1       subtype;

  if (vnpp == NULL || ompp == NULL) return;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    val = NULL;
    subtype = 0;
    StringHasOrgModPrefix (str, &val, &subtype, TRUE);
    if (val != NULL) {
      numspaces = 0;
      numcommas = 0;
      ptr = str;
      ch = *ptr;
      while (ch != '\0') {
        if (ch == ' ') {
          numspaces++;
        } else if (ch == ',') {
          numcommas++;
        }
        ptr++;
        ch = *ptr;
      }
      if (numspaces > 4 || numcommas > 0) {
        val = NULL;
      }
    }
    if (val != NULL) {
      omp = OrgModNew ();
      if (omp != NULL) {
        omp->subtype = (Uint1) subtype;
        omp->subname = StringSave (val);
        omp->next = *ompp;
        *ompp = omp;
      }
      *vnpp = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      vnpp = &(vnp->next);
    }
    vnp = next;
  }
}

static void StringHasSubSourcePrefix (CharPtr str, CharPtr PNTR pval, Uint1Ptr p_subtypeval, Boolean skippref)
{
  Int2          i;
  CharPtr       val = NULL;
  Uint1         subtype_val = 0;
  
  for (i = 0; current_subsource_subtype_alist[i].name != NULL && subtype_val == 0; i++) {
    val = StringHasPrefix (str, current_subsource_subtype_alist [i].name,
                           (Boolean) (current_subsource_subtype_alist[i].value == SUBSRC_germline ||
                                      current_subsource_subtype_alist[i].value == SUBSRC_rearranged ||
                                      current_subsource_subtype_alist[i].value == SUBSRC_transgenic ||
                                      current_subsource_subtype_alist[i].value == SUBSRC_environmental_sample ||
                                      current_subsource_subtype_alist[i].value == SUBSRC_metagenomic),
                           skippref);
    if (val != NULL) {
      subtype_val = current_subsource_subtype_alist[i].value;
    }
  }  
  if (subtype_val == 0) {
    for (i = 0; subsource_aliases[i].name != NULL && subtype_val == 0; i++) {
      val = StringHasPrefix (str, subsource_aliases [i].alias,
                             (Boolean) (subsource_aliases[i].value == SUBSRC_germline ||
                                        subsource_aliases[i].value == SUBSRC_rearranged ||
                                        subsource_aliases[i].value == SUBSRC_transgenic ||
                                        subsource_aliases[i].value == SUBSRC_environmental_sample ||
                                        subsource_aliases[i].value == SUBSRC_metagenomic),
                             skippref);
      if (val != NULL) {
        subtype_val = subsource_aliases[i].value;
      }
    }
  }
  if (pval != NULL) {
    *pval = val;
  }
  if (p_subtypeval != NULL) {
    *p_subtypeval = subtype_val;
  }
}

static void OrpModToSubSource (ValNodePtr PNTR vnpp, SubSourcePtr PNTR sspp)

{
  Char          ch;
  ValNodePtr    next;
  Int2          numcommas;
  Int2          numspaces;
  CharPtr       ptr;
  SubSourcePtr  ssp;
  CharPtr       str;
  CharPtr       val;
  ValNodePtr    vnp;
  Uint1         subtype_val = 0;

  if (vnpp == NULL || sspp == NULL) return;
  vnp = *vnpp;
  while (vnp != NULL) {
    next = vnp->next;
    str = (CharPtr) vnp->data.ptrvalue;
    val = NULL;
    subtype_val = 0;
    StringHasSubSourcePrefix (str, &val, &subtype_val, TRUE);

    if (val != NULL) {
      numspaces = 0;
      numcommas = 0;
      ptr = str;
      ch = *ptr;
      while (ch != '\0') {
        if (ch == ' ') {
          numspaces++;
        } else if (ch == ',') {
          numcommas++;
        }
        ptr++;
        ch = *ptr;
      }
      if (numspaces > 4 || numcommas > 0) {
        val = NULL;
      }
    }
    if (val != NULL) {
      ssp = SubSourceNew ();
      if (ssp != NULL) {
        ssp->subtype = subtype_val;
        ssp->name = StringSave (val);
        ssp->next = *sspp;
        *sspp = ssp;
      }
      *vnpp = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      vnpp = &(vnp->next);
    }
    vnp = next;
  }
}

static void GbqualToOrpMod (GBQualPtr PNTR prevgbq, ValNodePtr PNTR vnpp)

{
  GBQualPtr  gbq;
  size_t     len;
  GBQualPtr  next;
  CharPtr    str;
  Boolean    unlink;
  CharPtr    val;
  Uint1      subtype_val;

  if (prevgbq == NULL) return;
  gbq = *prevgbq;
  while (gbq != NULL) {
    next = gbq->next;
    unlink = FALSE;
    str = gbq->qual;
    if (str != NULL) {
      val = NULL;
      subtype_val = 0;
      StringHasOrgModPrefix (str, &val, &subtype_val, FALSE);
      if (val == NULL) {
        subtype_val = 0;
        StringHasSubSourcePrefix (str, &val, &subtype_val, FALSE);

      }
      if (val != NULL) {
        len = StringLen (gbq->val);
        str = MemNew (sizeof (Char) * (len + 64));
        if (str != NULL) {
          StringCpy (str, val);
          StringCat (str, "=");
          StringCat (str, gbq->val);
          ValNodeAddStr (vnpp, 0, str);
          unlink = TRUE; 
        }
      }
    }
    if (unlink) {
      *prevgbq = gbq->next;
      gbq->next = NULL;
      GBQualFree (gbq);
    } else {
      prevgbq = (GBQualPtr PNTR) &(gbq->next);
    }
    gbq = next;
  }
}

#define IS_WHITESP(c) (((c) == ' ') || ((c) == '\n') || ((c) == '\r') || ((c) == '\t'))

static Boolean IsStringSingleToken (CharPtr str)

{
  Char  ch;

  if (StringHasNoText (str)) return FALSE;

  ch = *str;
  while (ch != '\0') {
    if (IS_WHITESP (ch)) return FALSE;
    str++;
    ch = *str;
  }

  return TRUE;
}

static CharPtr FindAnOrgMod (OrgNamePtr onp, Uint1 subtype)

{
  OrgModPtr  omp;

  if (onp == NULL || subtype == 0) return NULL;

  for (omp = onp->mod; omp != NULL; omp = omp->next) {
    if (omp->subtype != subtype) continue;
    if (StringHasNoText (omp->subname)) continue;
    return omp->subname;
  }

  return NULL;
}

static CharPtr FindASubSource (BioSourcePtr biop, Uint1 subtype)

{
  SubSourcePtr  ssp;

  if (biop == NULL || subtype == 0) return NULL;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype != subtype) continue;
    if (StringHasNoText (ssp->name)) continue;
    return ssp->name;
  }

  return NULL;
}

static CharPtr FindNextSingleTilde (CharPtr str)

{
  Char  ch;

  if (StringHasNoText (str)) return NULL;

  ch = *str;
  while (ch != '\0') {
    if (ch == ' ') {
      if (str [1] == '~') {
        str++;
        ch = *str;
        while (ch == '~') {
          str++;
          ch = *str;
        }
      } else {
        str++;
        ch = *str;
      }
    } else if (ch == '~') {
      if (str [1] != '~') return str;
      str++;
      ch = *str;
      while (ch == '~') {
        str++;
        ch = *str;
      }
    } else {
      str++;
      ch = *str;
    }
  }

  return NULL;
}

static ValNodePtr SplitAtSingleTilde (CharPtr strs)

{
  ValNodePtr  head = NULL;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  str = tmp;

  while (StringDoesHaveText (str)) {
    ptr = FindNextSingleTilde (str);
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

static CharPtr MergeTildeStrings (ValNodePtr head)

{
  size_t      len = 0;
  CharPtr     prefix = "", ptr, str;
  ValNodePtr  vnp;

  if (head == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    len += StringLen (str) + 1;
  }
  if (len < 1) return NULL;

  ptr = MemNew (sizeof (Char) * (len + 2));
  if (ptr == NULL) return NULL;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    StringCat (ptr, prefix);
    StringCat (ptr, str);
    prefix = "~";
  }

  return ptr;
}


static void CleanupOrgModOther (BioSourcePtr biop, OrgNamePtr onp)

{
  ValNodePtr      head, vnp;
  OrgModPtr       next;
  OrgModPtr       omp;
  OrgModPtr PNTR  prev;
  CharPtr         str;
  Uint1           subtype_val;
  CharPtr         tmp;
  Boolean         unlink;
  CharPtr         val;

  if (biop == NULL || onp == NULL) return;

  prev = &(onp->mod);
  omp = onp->mod;
  while (omp != NULL) {
    next = omp->next;
    unlink= FALSE;
    if (omp->subtype == ORGMOD_other) {
      str = omp->subname;
      head = SplitAtSingleTilde (str);
      for (vnp = head; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        val = NULL;
        subtype_val = 0;
        StringHasOrgModPrefix (str, &val, &subtype_val, TRUE);
        if (val != NULL) {
          tmp = FindAnOrgMod (onp, subtype_val);
          if (tmp != NULL && StringICmp (tmp, val) == 0) {
            vnp->data.ptrvalue = NULL;
          }
        } else {
          subtype_val = 0;
          StringHasSubSourcePrefix (str, &val, &subtype_val, TRUE);
          if (val != NULL) {
            tmp = FindASubSource (biop, subtype_val);
            if (tmp != NULL && StringICmp (tmp, val) == 0) {
              vnp->data.ptrvalue = NULL;
            }
          }
        }
      }
      str = MergeTildeStrings (head);
      ValNodeFreeData (head);
      omp->subname = MemFree (omp->subname);
      omp->subname = str;
      if (StringHasNoText (str)) {
        unlink = TRUE;
      }
    } else if (omp->subtype == ORGMOD_bio_material 
               || omp->subtype == ORGMOD_culture_collection
               || omp->subtype == ORGMOD_specimen_voucher) {
      /*
      FixOrgModVoucher (omp);
      */
    }
    if (unlink) {
      *prev = omp->next;
      omp->next = NULL;
      OrgModFree (omp);
    } else {
      prev = &(omp->next);
    }
    omp = next;
  }
}

static void CleanupSubSourceOther (BioSourcePtr biop, OrgNamePtr onp)

{
  ValNodePtr         head, vnp;
  SubSourcePtr       next;
  SubSourcePtr PNTR  prev;
  SubSourcePtr       ssp;
  CharPtr            str;
  Uint1              subtype_val;
  CharPtr            tmp;
  Boolean            unlink;
  CharPtr            val;

  if (biop == NULL /* || onp == NULL */ ) return;

  prev = &(biop->subtype);
  ssp = biop->subtype;
  while (ssp != NULL) {
    next = ssp->next;
    unlink = FALSE;
    if (ssp->subtype == SUBSRC_other) {
      str = ssp->name;
      head = SplitAtSingleTilde (str);
      for (vnp = head; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        val = NULL;
        subtype_val = 0;
        StringHasOrgModPrefix (str, &val, &subtype_val, TRUE);
        if (val != NULL) {
          tmp = FindAnOrgMod (onp, subtype_val);
          if (tmp != NULL && StringICmp (tmp, val) == 0) {
            vnp->data.ptrvalue = NULL;
          }
        } else {
          subtype_val = 0;
          StringHasSubSourcePrefix (str, &val, &subtype_val, TRUE);
          if (val != NULL) {
            tmp = FindASubSource (biop, subtype_val);
            if (tmp != NULL && StringICmp (tmp, val) == 0) {
              vnp->data.ptrvalue = NULL;
            }
          }
        }
      }
      str = MergeTildeStrings (head);
      ValNodeFreeData (head);
      ssp->name = MemFree (ssp->name);
      ssp->name = str;
      if (StringHasNoText (str)) {
        unlink = TRUE;
      }
    }
    if (unlink) {
      *prev = ssp->next;
      ssp->next = NULL;
      SubSourceFree (ssp);
    } else {
      prev = &(ssp->next);
    }
    ssp = next;
  }
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
  } else if (str1 != NULL) {
    return 1;
  } else if (str2 != NULL) {
    return -1;
  }
  return 0;
}

static void FixNumericDbxref (DbtagPtr dbt)

{
  size_t       len;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  long         val;

  if (dbt != NULL) {
    oip = dbt->tag;
    if (oip != NULL) {
      ptr = oip->str;
      if (ptr != NULL && *ptr != '0' && StringIsAllDigits(ptr)) {
        len = StringLen (ptr);
        if (len < 10 || (len == 10 && StringCmp (ptr, "2147483647") <= 0)) {
          if (sscanf (oip->str, "%ld", &val) == 1) {
            oip->id = (Int4) val;
            oip->str = MemFree (oip->str);
          }
        }
      }
    }
  }
}

static void FixNumericDbxrefs (ValNodePtr vnp)

{
  DbtagPtr  dbt;

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      FixNumericDbxref (dbt);
    }
    vnp = vnp->next;
  }
}

static void FixOldDbxref (DbtagPtr dbt)

{
  Boolean      all_digits;
  Char         buf [32];
  Char         ch;
  CharPtr      ident;
  size_t       len;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  CharPtr      str;

  if (dbt != NULL) {

    TrimSpacesAroundString (dbt->db);
    oip = dbt->tag;
    if (oip != NULL && oip->str != NULL) {
      /*
      TrimSpacesAroundString (oip->str);
      */
      TrimSpacesSemicolonsAndCommas (oip->str);
    }

    if (StringICmp (dbt->db, "SWISS-PROT") == 0 &&
        StringCmp (dbt->db, "Swiss-Prot") != 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("Swiss-Prot");
    } else if (StringICmp (dbt->db, "SPTREMBL") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("TrEMBL");
    } else if (StringICmp (dbt->db, "SUBTILIS") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("SubtiList");
    } else if (StringICmp (dbt->db, "MGD") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("MGI");
    } else if (StringCmp (dbt->db, "cdd") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("CDD");
    } else if (StringCmp (dbt->db, "FlyBase") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("FLYBASE");
    } else if (StringCmp (dbt->db, "GENEDB") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("GeneDB");
    } else if (StringCmp (dbt->db, "GreengenesID") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("Greengenes");
    } else if (StringCmp (dbt->db, "HMPID") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("HMP");
    }
    if (StringICmp (dbt->db, "HPRD") == 0) {
      oip = dbt->tag;
      if (oip != NULL && StringDoesHaveText (oip->str)) {
        str = oip->str;
        if (str != NULL && StringNICmp (str, "HPRD_", 5) == 0) {
          str [0] = ' ';
          str [1] = ' ';
          str [2] = ' ';
          str [3] = ' ';
          str [4] = ' ';
          TrimSpacesAroundString (str);
        }
      }
    } else if (StringICmp (dbt->db, "MGI") == 0) {
      oip = dbt->tag;
      if (oip != NULL && oip->str != NULL && StringDoesHaveText (oip->str)) {
        str = oip->str;
        if (StringNICmp (str, "MGI:", 4) == 0 || StringNICmp (str, "MGD:", 4) == 0) {
          str [0] = ' ';
          str [1] = ' ';
          str [2] = ' ';
          str [3] = ' ';
          TrimSpacesAroundString (str);
        } else if (StringNICmp (str, "J:", 2) == 0) {
          ptr = str + 2;
          ch = *ptr;
          all_digits = TRUE;
          while (ch != '\0') {
            if (! IS_DIGIT (ch)) {
              all_digits = FALSE;
            }
            ptr++;
            ch = *ptr;
          }
          if (all_digits) {
            oip->str = MemFree (oip->str);
            oip->str = StringSave ("");
          }
        }
      }
    }
    if (StringICmp (dbt->db, "Swiss-Prot") == 0 ||
        StringICmp (dbt->db, "SWISSPROT") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("UniProt/Swiss-Prot");
    } else if (StringICmp (dbt->db, "TrEMBL") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("UniProt/TrEMBL");
    } else if (StringICmp (dbt->db, "LocusID") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("GeneID");
    } else if (StringICmp (dbt->db, "MaizeDB") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("MaizeGDB");
    }
    if (StringICmp (dbt->db, "UniProt/Swiss-Prot") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("UniProtKB/Swiss-Prot");
    } else if (StringICmp (dbt->db, "UniProt/TrEMBL") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("UniProtKB/TrEMBL");
    } else if (StringICmp (dbt->db, "Genew") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("HGNC");
    } else if (StringICmp (dbt->db, "IFO") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("NBRC");
    } else if (StringICmp (dbt->db, "BHB") == 0 ||
        StringICmp (dbt->db, "BioHealthBase") == 0) {
      dbt->db = MemFree (dbt->db);
      dbt->db = StringSave ("IRD");
    }

    oip = dbt->tag;
    if (oip != NULL && oip->str != NULL) {
      ident = oip->str;
      if (StringCmp (dbt->db, "HGNC") == 0 && StringNCmp (ident, "HGNC:", 5) == 0 ) {
        ident += 5;
        ptr = StringSave (ident);
        oip->str = MemFree (oip->str);
        oip->str = ptr;
      } else if (StringCmp (dbt->db, "VGNC") == 0 && StringNCmp (ident, "VGNC:", 5) == 0 ) {
        ident += 5;
        ptr = StringSave (ident);
        oip->str = MemFree (oip->str);
        oip->str = ptr;
      } else if (StringCmp (dbt->db, "MGI") == 0 && StringNCmp (ident, "MGI:", 4) == 0 ) {
        ident += 4;
        ptr = StringSave (ident);
        oip->str = MemFree (oip->str);
        oip->str = ptr;
      } else if (StringCmp (dbt->db, "RGD") == 0 && StringNCmp (ident, "RGD:", 4) == 0 ) {
        ident += 4;
        ptr = StringSave (ident);
        oip->str = MemFree (oip->str);
        oip->str = ptr;
      }
    }
    if (oip != NULL) {
      if (StringCmp (dbt->db, "HGNC") == 0 || StringCmp (dbt->db, "VGNC") == 0 || StringCmp (dbt->db, "MGI") == 0) {
        if (oip->str == NULL && oip->id > 0) {
          sprintf (buf, "%ld", (long) oip->id);
          ptr = StringSave (buf);
          oip->id = 0;
          oip->str = ptr;
        }
        ident = oip->str;
        if (ident != NULL) {
          if (StringChr (ident, ':') == NULL) {
            len = StringLen (dbt->db) + StringLen (ident) + 5;
            ptr = (CharPtr) MemNew (sizeof (Char) * len);
            if (ptr != NULL) {
              sprintf (ptr, "%s:%s", dbt->db, ident);
              oip->str = MemFree (oip->str);
              oip->str = ptr;
            }
          }
        }
      }
    }
  }
}

static void FixOldDbxrefs (ValNodePtr vnp, Boolean isEmblOrDdbj)

{
  DbtagPtr     dbt;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  CharPtr      tmp;
  ValNodePtr   vp2;

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      FixOldDbxref (dbt);

      if (! isEmblOrDdbj) {
        if (StringCmp (dbt->db, "HGNC") != 0 && StringCmp (dbt->db, "VGNC") != 0 && StringCmp (dbt->db, "MGI") != 0) {
          /* expand db_xrefs with colons inside tags */
          oip = dbt->tag;
          if (oip != NULL && oip->str != NULL) {
            ptr = StringChr (oip->str, ':');
            if (ptr != NULL) {
              if (StringHasNoText (ptr + 1)) {
                *ptr = '\0';
              } else {
                tmp = dbt->db;
                dbt = DbtagNew ();
                if (dbt != NULL) {
                  oip = ObjectIdNew ();
                  if (oip != NULL) {
                    vp2 = ValNodeNew (NULL);
                    if (vp2 != NULL) {
                      *ptr = '\0';
                      ptr++;
                      TrimSpacesAroundString (ptr);
                      dbt->db = StringSave (tmp);
                      oip->str = StringSave (ptr);
                      dbt->tag = oip;
                      vp2->data.ptrvalue = (Pointer) dbt;
                      vp2->next = vnp->next;
                      vnp->next = vp2;
                    }
                  }
                }
              }
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

static void CleanupObsoleteDbxrefs (ValNodePtr PNTR prevvnp)

{
  DbtagPtr     dbt;
  ValNodePtr   nextvnp;
  ObjectIdPtr  oip;
  CharPtr      str;
  Boolean      unlink;
  ValNodePtr   vnp;

  if (prevvnp == NULL) return;
  vnp = *prevvnp;
  while (vnp != NULL) {
    nextvnp = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      unlink = FALSE;
      str = (CharPtr) dbt->db;
      if (StringHasNoText (str) ||
          StringICmp (str, "PID") == 0 ||
          StringICmp (str, "PIDg") == 0 ||
          /*
          StringICmp (str, "PIDe") == 0 ||
          StringICmp (str, "PIDd") == 0 ||
          */
          /*
          StringICmp (str, "GI") == 0 ||
          */
          StringICmp (str, "NID") == 0) {
        unlink = TRUE;
      }
      oip = dbt->tag;
      if (oip == NULL) {
        unlink = TRUE;
      } else if (oip->str != NULL) {
        if (StringHasNoText (oip->str)) {
          unlink = TRUE;
        }
      } else if (oip->id == 0) {
        unlink = TRUE;
      }
      if (unlink) {
        *prevvnp = vnp->next;
        vnp->next = NULL;
        DbtagFree (dbt);
        ValNodeFree (vnp);
      } else {
        prevvnp = (ValNodePtr PNTR) &(vnp->next);
      }
    }
    vnp = nextvnp;
  }
}

static void CleanupGoDbxrefs (ValNodePtr vnp)

{
  DbtagPtr     dbt;
  size_t       idx;
  size_t       len;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  Char         tmp [32];

  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (StringICmp (dbt->db, "GO") == 0) {
        oip = dbt->tag;
        if (oip != NULL) {
          if (oip->str == NULL && oip->id > 0) {
            sprintf (tmp, "%ld", (long) oip->id);
            oip->str = StringSave (tmp);
            oip->id = 0;
          }
          ptr = oip->str;
          if (ptr != NULL && StringIsAllDigits(ptr)) {
            len = StringLen (ptr);
            if (len < 7) {
              idx = 7 - len;
              StringCpy (tmp, "0000000");
              tmp [idx] = '\0';
              StringCat (tmp, ptr);
              oip->str = MemFree (oip->str);
              oip->str = StringSave (tmp);
            }
          }
        }
      }
    }
    vnp = vnp->next;
  }
}

static int LIBCALLBACK SortCits (VoidPtr ptr1, VoidPtr ptr2)

{
  int         compare;
  Char        label1 [128], label2 [128];
  ValNodePtr  ppr1, ppr2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  ppr1 = *((ValNodePtr PNTR) ptr1);
  ppr2 = *((ValNodePtr PNTR) ptr2);
  if (ppr1 == NULL || ppr2 == NULL) return 0;
  PubLabel (ppr1, label1, 127, OM_LABEL_CONTENT);
  PubLabel (ppr2, label2, 127, OM_LABEL_CONTENT);
  compare = StringICmp (label1, label2);
  return compare;
}

static Boolean CitGenTitlesMatch (ValNodePtr pub1, ValNodePtr pub2)

{
  CitGenPtr  cgp1, cgp2;

  if (pub1->choice == PUB_Gen) {
    cgp1 = (CitGenPtr) pub1->data.ptrvalue;
    if (cgp1->serial_number != -1 && pub1->next != NULL) {
      pub1 = pub1->next;
    }
  }
  if (pub2->choice == PUB_Gen) {
    cgp2 = (CitGenPtr) pub2->data.ptrvalue;
    if (cgp2->serial_number != -1 && pub2->next != NULL) {
      pub2 = pub2->next;
    }
  }

  if (pub1->choice != PUB_Gen || pub2->choice != PUB_Gen) return TRUE;
  cgp1 = (CitGenPtr) pub1->data.ptrvalue;
  cgp2 = (CitGenPtr) pub2->data.ptrvalue;
  if (cgp1->title == NULL || cgp2->title == NULL) return TRUE;
  if (StringCmp (cgp1->title, cgp2->title) != 0) return FALSE;
  return TRUE;
}

static void CleanupDuplicateCits (ValNodePtr PNTR prevvnp)

{
  Char        label1 [128], label2 [128];
  ValNodePtr  last = NULL;
  ValNodePtr  nextvnp;
  Boolean     unlink;
  ValNodePtr  vnp;

  if (prevvnp == NULL) return;
  vnp = *prevvnp;
  while (vnp != NULL) {
    nextvnp = vnp->next;
    unlink = FALSE;
    if (last != NULL) {
      PubLabelUnique (last, label1, 127, OM_LABEL_CONTENT, TRUE);
      PubLabelUnique (vnp, label2, 127, OM_LABEL_CONTENT, TRUE);
      if (StringCmp (label1, label2) == 0 && CitGenTitlesMatch (last, vnp)) {
        unlink = TRUE;
      }
    } else {
      last = vnp;
    }
    if (unlink) {
      *prevvnp = vnp->next;
      vnp->next = NULL;
      PubFree (vnp);
    } else {
      last = vnp;
      prevvnp = (ValNodePtr PNTR) &(vnp->next);
    }
    vnp = nextvnp;
  }
}

/* name processing code from Sequin editors */

NLM_EXTERN void FirstNameToInitials (CharPtr first, CharPtr inits, size_t maxsize)

{
  Char  ch;
  Uint2  i;

  if (inits != NULL && maxsize > 0) {
    inits [0] = '\0';
    if (first != NULL) {
      i = 0;
      ch = *first;
      while (ch != '\0' && i < maxsize) {
        while (ch != '\0' && (ch <= ' ' || ch == '-')) {
          first++;
          ch = *first;
        }
        if (IS_ALPHA (ch)) {
          inits [i] = ch;
          i++;
          first++;
          ch = *first;
        }
        while (ch != '\0' && ch > ' ' && ch != '-') {
          first++;
          ch = *first;
        }
        if (ch == '-') {
          inits [i] = ch;
          i++;
          first++;
          ch = *first;
        }
      }
      inits [i] = '\0';
    }
  }
}

static void StripPeriods (CharPtr str)

{
  Char     ch;
  CharPtr  dst;

  if (str != NULL) {
    dst = str;
    ch = *str;
    while (ch != '\0') {
      if (ch != '.') {
        *dst = ch;
        dst++;
      }
      str++;
      ch = *str;
    }
    *dst = '\0';
  }
}

static void TrimLeadingSpaces (CharPtr str)

{
  Char     ch;
  CharPtr  dst;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ch = *str;
    while (ch != '\0' && ch <= ' ') {
      str++;
      ch = *str;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      str++;
      ch = *str;
    }
    *dst = '\0';
  }
}

static void ExtractSuffixFromInitials (NameStdPtr nsp)

{
  Char     ch;
  Boolean  has_period = FALSE;
  size_t   len;
  CharPtr  str;

  str = nsp->names [4];
  ch = *str;
  while (ch != '\0') {
    if (ch == '.') {
      has_period = TRUE;
    }
    str++;
    ch = *str;
  }
  if (! has_period) return;
  str = nsp->names [4];
  len = StringLen (str);
  if (len >= 4 && StringCmp (str +  len - 3, "III") == 0) {
    str [len - 3] = '\0';
    nsp->names [5] = StringSave ("III");
  } else if (len >= 5 && StringCmp (str +  len - 4, "III.") == 0) {
    str [len - 4] = '\0';
    nsp->names [5] = StringSave ("III");
  } else if (len >= 3 && StringCmp (str +  len - 2, "Jr") == 0) {
    str [len - 2] = '\0';
    nsp->names [5] = StringSave ("Jr");
  } else if (len >= 4 && StringCmp (str +  len - 3, "2nd") == 0) {
    str [len - 3] = '\0';
    nsp->names [5] = StringSave ("II");
  } else if (len >= 3 && StringCmp (str +  len - 2, "IV") == 0) {
    str [len - 2] = '\0';
    nsp->names [5] = StringSave ("IV");
  } else if (len >= 4 && StringCmp (str +  len - 3, "IV.") == 0) {
    str [len - 3] = '\0';
    nsp->names [5] = StringSave ("IV");
  }
}

static CharPtr NameStdPtrToTabbedString (NameStdPtr nsp, Boolean fixInitials)

{
  Char   first [256];
  Char   frstinits [64];
  Char   initials [64];
  Int2   j;
  Char   last [256];
  Char   middle [128];
  Char   str [512];
  Char   suffix [64];

  if (nsp == NULL) return NULL;
  if (nsp->names [5] == NULL && nsp->names [4] != NULL) {
    ExtractSuffixFromInitials (nsp);
  }
  str [0] = '\0';
  StringNCpy_0 (first, nsp->names [1], sizeof (first));
  TrimSpacesAroundString (first);
  StringNCpy_0 (initials, nsp->names [4], sizeof (initials));
  StripPeriods (initials);
  TrimLeadingSpaces (initials);
  StringNCpy_0 (last, nsp->names [0], sizeof (last));
  TrimLeadingSpaces (last);
  StringNCpy_0 (middle, nsp->names [2], sizeof (middle));
  TrimLeadingSpaces (middle);
  if (StringCmp (initials, "al") == 0 &&
      StringCmp (last, "et") == 0 &&
      first [0] == '\0') {
    initials [0] = '\0';
    StringCpy (last, "et al.");
  }
  /*
  if (first [0] == '\0') {
    StringNCpy_0 (first, initials, sizeof (first));
    if (IS_ALPHA (first [0])) {
      if (first [1] == '-') {
        first [3] = '\0';
      } else {
        first [1] = '\0';
      }
    } else {
      first [0] = '\0';
    }
  }
  */
  frstinits [0] = '\0';
  FirstNameToInitials (first, frstinits, sizeof (frstinits) - 1);
  StripPeriods (first);
  TrimLeadingSpaces (first);
  if (first [0] != '\0') {
    StringCat (str, first);
  } else {
    /*
    StringCat (str, " ");
    */
  }
  StringCat (str, "\t");
  if (fixInitials) {
    j = 0;
    while (initials [j] != '\0' && TO_UPPER (initials [j]) == TO_UPPER (frstinits [j])) {
      j++;
    }
    if (initials [j] != '\0') {
      StringCat (str, initials + j);
    } else {
      /*
      StringCat (str, " ");
      */
    }
  } else if (initials [0] != '\0') {
    StringCat (str, initials);
  } else if (frstinits [0] != '\0') {
    StringCat (str, frstinits);
  }
  StringCat (str, "\t");
  StringCat (str, last);
  StringNCpy_0 (suffix, nsp->names [5], sizeof (suffix));
  StringCat (str, "\t");
  StripPeriods (suffix);
  TrimLeadingSpaces (suffix);
  if (suffix [0] != '\0') {
    StringCat (str, suffix);
  } else {
    /*
    StringCat (str, " ");
    */
  }
  StringCat (str, "\t");
  StringCat (str, middle);
  StringCat (str, "\n");
  return StringSave (str);
}

static CharPtr XtractTagListColumn (CharPtr source, Int2 col)

{
  Char     ch;
  size_t   count;
  CharPtr  ptr;
  CharPtr  str;

  if (source == NULL || source [0] == '\0' || col < 0) return NULL;

  ptr = source;
  ch = *ptr;
  while (col > 0 && ch != '\n' && ch != '\0') {
    while (ch != '\t' && ch != '\n' && ch != '\0') {
      ptr++;
      ch = *ptr;
    }
    if (ch == '\t') {
      ptr++;
      ch = *ptr;
    }
    col--;
  }

  count = 0;
  ch = ptr [count];
  while (ch != '\t' && ch != '\n' && ch != '\0') {
    count++;
    ch = ptr [count];
  }
  str = (CharPtr) MemNew(count + 1);
  if (str != NULL) {
    MemCpy (str, ptr, count);
  }
  return str;
}

static NameStdPtr TabbedStringToNameStdPtr (CharPtr txt, Boolean fixInitials)

{
  Char        ch;
  CharPtr     first;
  Char        initials [64];
  Int2        j;
  Int2        k;
  Char        last;
  Int2        len;
  NameStdPtr  nsp;
  Char        periods [128];
  CharPtr     str;
  Char        str1 [64];
  Char        suffix [80];

  if (txt == NULL) return NULL;
  nsp = NameStdNew ();
  if (nsp == NULL) return NULL;
  nsp->names [0] = XtractTagListColumn (txt, 2);
  TrimLeadingSpaces (nsp->names [0]);
  first = XtractTagListColumn (txt, 0);
  StripPeriods (first);
  nsp->names [1] = StringSave (first);
  TrimLeadingSpaces (nsp->names [1]);
  str1 [0] = '\0';
  if (fixInitials) {
    FirstNameToInitials (first, str1, sizeof (str1) - 1);
  }
  str = XtractTagListColumn (txt, 1);
  StringNCat (str1, str, sizeof (str1) - 1);
  MemFree (str);
  j = 0;
  k = 0;
  ch = str1 [j];
  while (ch != '\0') {
    if (ch != ' ') {
      initials [k] = ch;
      k++;
    }
    j++;
    ch = str1 [j];
  }
  initials [k] = '\0';
  periods [0] = '\0';
          j = 0;
          ch = initials [j];
          while (ch != '\0') {
            if (ch == ',') {
              initials [j] = '.';
            }
            j++;
            ch = initials [j];
          }
          str = StringStr (initials, ".ST.");
          if (str != NULL) {
            *(str + 2) = 't';
          }
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
              last = ch;
      k++;
      j++;
      ch = initials [j];
              if (ch == '\0') {
                if (! (IS_LOWER (last))) {
                  periods [k] = '.';
                  k++;
                }
              /* } else if (ch == '.' && initials [j + 1] == '\0') { */
              } else if (! (IS_LOWER (ch))) {
                periods [k] = '.';
                k++;
              }
    }
  }
  if (k > 0 && periods [k - 1] != '.') {
    periods [k] = '.';
    k++;
  }
  periods [k] = '\0';
  nsp->names [4] = StringSave (periods);
  TrimLeadingSpaces (nsp->names [4]);
  str = XtractTagListColumn (txt, 3);
  StringNCpy_0 (str1, str, sizeof (str1));
  MemFree (str);
  j = 0;
  k = 0;
  ch = str1 [j];
  while (ch != '\0') {
    if (ch != ' ') {
      suffix [k] = ch;
      k++;
    }
    j++;
    ch = str1 [j];
  }
  suffix [k] = '\0';
  if (suffix [0] != '\0') {
    len = StringLen (suffix);
    if (len > 0 && suffix [len - 1] == '.') {
      suffix [len - 1] = '\0';
    }
    if (StringICmp (suffix, "1d") == 0) {
      StringCpy (suffix, "I");
    } else if (StringICmp (suffix, "1st") == 0) {
      StringCpy (suffix, "I");
    } else if (StringICmp (suffix, "2d") == 0) {
      StringCpy (suffix, "2nd");
    } else if (StringICmp (suffix, "3d") == 0) {
      StringCpy (suffix, "3rd");
    } else if (StringICmp (suffix, "Sr") == 0) {
      StringCpy (suffix, "Sr.");
    } else if (StringICmp (suffix, "Jr") == 0) {
      StringCpy (suffix, "Jr.");
    }
    /*
    len = StringLen (suffix);
    if (len > 0 && suffix [len - 1] != '.') {
      StringCat (suffix, ".");
    }
    */
    nsp->names [5] = StringSave (suffix);
    TrimLeadingSpaces (nsp->names [5]);
  }
  if (StringCmp (nsp->names [0], "et al") == 0) {
    nsp->names [0] = MemFree (nsp->names [0]);
    nsp->names [0] = StringSave ("et al.");
  }
  nsp->names [2] = XtractTagListColumn (txt, 4);
  TrimLeadingSpaces (nsp->names [2]);
  if (StringHasNoText (nsp->names [2])) {
    nsp->names [2] = MemFree (nsp->names [2]);
  }
  MemFree (first);
  return nsp;
}

static AffilPtr CleanAffil (AffilPtr afp)

{
  if (afp == NULL) return NULL;
  CleanVisStringJunkAndCompress (&(afp->affil));
  if (afp->choice == 2) {
    CleanVisStringJunkAndCompress (&(afp->div));
    CleanVisStringJunkAndCompress (&(afp->city));
    CleanVisStringJunkAndCompress (&(afp->sub));
    CleanVisStringJunkAndCompress (&(afp->country));
    CleanVisStringJunkAndCompress (&(afp->street));
    CleanVisStringJunkAndCompress (&(afp->email));
    CleanVisStringJunkAndCompress (&(afp->fax));
    CleanVisStringJunkAndCompress (&(afp->phone));
    CleanVisStringJunkAndCompress (&(afp->postal_code));
    TrimSpacesSemicolonsAndCommas (afp->postal_code);
    if (StringICmp (afp->country, "U.S.A.") == 0) {
      afp->country = MemFree (afp->country);
      afp->country = StringSave ("USA");
    }
    if (StringICmp (afp->country, "USA") == 0 && StringCmp (afp->country, "USA") != 0) {
      afp->country = MemFree (afp->country);
      afp->country = StringSave ("USA");
    }
    if (StringCmp (afp->country, "USA") == 0 && afp->sub != NULL) {
      StripPeriods (afp->sub);
      TrimSpacesAroundString (afp->sub);
    }
  }
  if (afp->affil == NULL &&
      afp->div == NULL &&
      afp->city == NULL &&
      afp->sub == NULL &&
      afp->country == NULL &&
      afp->street == NULL &&
      afp->email == NULL &&
      afp->fax == NULL &&
      afp->phone == NULL &&
      afp->postal_code == NULL) {
    afp = MemFree (afp);
  }
  return afp;
}

static void NormalizeAuthors (AuthListPtr alp, Boolean fixInitials)

{
  AuthorPtr        ap;
  CharPtr          initials;
  size_t           len;
  ValNodePtr       names;
  ValNodePtr       next;
  NameStdPtr       nsp;
  PersonIdPtr      pid;
  ValNodePtr PNTR  prev;
  CharPtr          str;
  Boolean          upcaseinits;
  ValNodePtr       vnp;
  Boolean          zap;

  if (alp == NULL) return;
  alp->affil = CleanAffil (alp->affil);

  if (alp->choice == 2 || alp->choice == 3) {
    for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      TrimSpacesAroundString (str);
      TrimSpacesAndJunkFromEnds (str, FALSE);
      Asn2gnbkCompressSpaces (str);
    }
  }
  if (alp->choice != 1) return;

  prev = &(alp->names);
  names = alp->names;
  while (names != NULL) {
    next = names->next;
    zap = FALSE;
    ap = names->data.ptrvalue;
    if (ap != NULL) {
      pid = ap->name;
      if (pid == NULL) {
        /* continue */
      } else if (pid->choice == 2) {
        nsp = pid->data;
        if (nsp != NULL /* && nsp->names [4] != NULL */) {
          upcaseinits = FALSE;
          initials = nsp->names [4];
          if (StringLen (initials) > 0) {
            if (IS_UPPER (initials [0])) {
              upcaseinits = TRUE;
            }
          }
          str = NameStdPtrToTabbedString (nsp, fixInitials);
          pid->data = NameStdFree (nsp);
          nsp = TabbedStringToNameStdPtr (str, fixInitials);
          if (upcaseinits) {
            initials = nsp->names [4];
            if (StringLen (initials) > 0) {
              if (IS_LOWER (initials [0])) {
                initials [0] = TO_UPPER (initials [0]);
              }
            }
          }
          pid->data = nsp;
          MemFree (str);
          CleanVisString (&(nsp->names [0]));
          CleanVisString (&(nsp->names [1]));
          CleanVisString (&(nsp->names [2]));
          CleanVisString (&(nsp->names [3]));
          CleanVisString (&(nsp->names [4]));
          CleanVisString (&(nsp->names [5]));
          CleanVisString (&(nsp->names [6]));
          if (StringCmp (nsp->names [0], "et") == 0 &&
              (StringCmp (nsp->names [4], "al") == 0 ||
               StringCmp (nsp->names [4], "al.") == 0 ||
               StringCmp (nsp->names [4], "Al.") == 0) &&
              (StringHasNoText (nsp->names [1]) ||
               StringCmp (nsp->names [1], "a") == 0)) {
            nsp->names [4] = MemFree (nsp->names [4]);
            nsp->names [1] = MemFree (nsp->names [1]);
            nsp->names [0] = MemFree (nsp->names [0]);
            nsp->names [0] = StringSave ("et al.");
          }
          str = nsp->names [0];
          len = StringLen (str);
          if (len > 4 && StringHasNoText (nsp->names [5])) {
            if (StringCmp (str + len - 4, " Jr.") == 0 ||
                StringCmp (str + len - 4, " Sr.") == 0) {
              nsp->names [5] = StringSave (str + len - 3);
              str [len - 4] = '\0';
              TrimSpacesAroundString (str);
            }
          }
          str = nsp->names [4];
          len = StringLen (str);
          if (len > 4 && StringHasNoText (nsp->names [5])) {
            if (StringCmp (str + len - 4, ".Jr.") == 0 ||
                StringCmp (str + len - 4, ".Sr.") == 0) {
              nsp->names [5] = StringSave (str + len - 3);
              str [len - 3] = '\0';
              TrimSpacesAroundString (str);
            }
          }
          if (StringHasNoText (nsp->names [0]) &&
              StringHasNoText (nsp->names [1]) &&
              StringHasNoText (nsp->names [2]) &&
              StringHasNoText (nsp->names [3]) &&
              StringHasNoText (nsp->names [4]) &&
              StringHasNoText (nsp->names [5]) &&
              StringHasNoText (nsp->names [6])) {
            zap = TRUE;
          }
          /* last name is required, so zap if not present */
          if (StringHasNoText (nsp->names [0])) {
            zap = TRUE;
          }
        }
      } else if (pid->choice == 3 || pid->choice == 4 || pid->choice == 5) {
        TrimSpacesAroundString ((CharPtr) pid->data);
        if (StringHasNoText ((CharPtr) pid->data)) {
          zap = TRUE;
        }
      }
    }
    if (zap) {
      /* remove empty authors */
      *prev = names->next;
      names->next = NULL;
      AuthorFree (ap);
      ValNodeFree (names);
    } else {
      prev = &(names->next);
    }
    names = next;
  }
  /* if no remaining authors, put in default author for legal ASN.1 */
  if (alp->names == NULL) {
    names = ValNodeNew (NULL);
    if (names != NULL) {
      /*
      ap = AuthorNew ();
      if (ap != NULL) {
        pid = PersonIdNew ();
        if (pid != NULL) {
          pid->choice = 4;
          pid->data = (Pointer) StringSave ("?");
          ap->name = pid;
          names->choice = 1;
          names->data.ptrvalue = ap;
          alp->names = names;
        }
      }
      */
      names->choice = 3;
      names->data.ptrvalue = StringSave ("?");
      alp->names = names;
      alp->choice = 3;
    }
  }
}

static void StrStripSpaces (
  CharPtr str
)

{
  CharPtr  new_str;

  if (str == NULL) return;

  new_str = str;
  while (*str != '\0') {
    *new_str++ = *str;
    if (*str == ' ' || *str == '\t' || *str == '(') {
      for (str++; *str == ' ' || *str == '\t'; str++) continue;
      if (*str == ')' || *str == ',') {
        if( *(new_str - 1) != '(' ) { // this if handles the case "\([ \t]*\)"
          --new_str;
        }
      }
    } else {
      str++;
    }
  }
  *new_str = '\0';
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

static void NormalizePubAuthors (ValNodePtr vnp, Boolean stripSerial, Boolean fixInitials)

{
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitPatPtr    cpp;
  CitSubPtr    csp;

  if (vnp == NULL) return;
  if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) return;
  if (vnp->data.ptrvalue == NULL) return;
  switch (vnp->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cgp->authors, fixInitials);
      break;
    case PUB_Sub :
      csp = (CitSubPtr) vnp->data.ptrvalue;
      NormalizeAuthors (csp->authors, fixInitials);
      break;
    case PUB_Article :
      cap = (CitArtPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cap->authors, fixInitials);
      if (cap->from == 2 || cap->from == 3) {
        cbp = (CitBookPtr) cap->fromptr;
        if (cbp != NULL) {
          NormalizeAuthors (cbp->authors, fixInitials);
        }
      }
      break;
    case PUB_Book :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cbp->authors, fixInitials);
      break;
    case PUB_Man :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      if (cbp->othertype == 2 && cbp->let_type == 3) {
        NormalizeAuthors (cbp->authors, fixInitials);
      }
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) vnp->data.ptrvalue;
      NormalizeAuthors (cpp->authors, fixInitials);
      NormalizeAuthors (cpp->applicants, fixInitials);
      NormalizeAuthors (cpp->assignees, fixInitials);
      break;
    default :
      break;
  }
}

static void NormalizeAPub (ValNodePtr vnp, Boolean stripSerial, Boolean fixInitials)

{
  AffilPtr     affil;
  AuthListPtr  alp;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitJourPtr   cjp;
  CitPatPtr    cpp;
  CitSubPtr    csp;
  ImprintPtr   imp;
  CharPtr      str;
  CharPtr      tmp;
  ValNodePtr   ttl;

  if (vnp == NULL) return;
  if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) return;
  if (vnp->data.ptrvalue == NULL) return;
  imp = NULL;
  switch (vnp->choice) {
    case PUB_Gen :
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (stripSerial) {
        cgp->serial_number = -1; /* but does not remove if empty */
      }
      if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
        cgp->cit [0] = 'U';
        /* cgp->date = DateFree (cgp->date); */ /* remove date if unpublished */
        if (cgp->journal == NULL) {
          cgp->volume = MemFree (cgp->volume);
          cgp->issue = MemFree (cgp->issue);
          cgp->pages = MemFree (cgp->pages);
        }
      }
      TrimSpacesAroundString (cgp->cit);
      if (StringDoesHaveText (cgp->title)) {
        StrStripSpaces (cgp->title);
      }
      break;
    case PUB_Sub :
      csp = (CitSubPtr) vnp->data.ptrvalue;
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
      if (imp != NULL && imp->date == NULL) {
        csp->imp = ImprintFree (csp->imp);
      }
      if (alp != NULL && alp->affil != NULL) {
        affil = alp->affil;
        if (affil->choice == 1) {
          str = affil->affil;
          if (StringNICmp (str, "to the ", 7) == 0) {
            if (StringNICmp (str + 24, " databases", 10) == 0) {
              str += 34;
              if (*str == '.') {
                str++;
              }
              tmp = StringSaveNoNull (TrimSpacesAroundString (str));
              affil->affil = MemFree (affil->affil);
              affil->affil = tmp;
            }
          }
        }
        alp->affil = CleanAffil (alp->affil);
      }
      imp = csp->imp;
      break;
    case PUB_Article :
      cap = (CitArtPtr) vnp->data.ptrvalue;
      if (cap != NULL) {
        if (cap->from == 1) {
          cjp = (CitJourPtr) cap->fromptr;
          if (cjp != NULL) {
            imp = cjp->imp;
          }
        } else if (cap->from == 2 || cap->from == 3) {
          cbp = (CitBookPtr) cap->fromptr;
          if (cbp != NULL) {
            imp = cbp->imp;
          }
        }
        for (ttl = cap->title; ttl != NULL; ttl = ttl->next) {
          if (ttl->choice == Cit_title_name) {
            str = (CharPtr) ttl->data.ptrvalue;
            if (StringHasNoText (str)) continue;
            StrStripSpaces (str);
          }
        }
      }
      break;
    case PUB_Book :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      if (cbp != NULL) {
        imp = cbp->imp;
      }
      break;
    case PUB_Man :
      cbp = (CitBookPtr) vnp->data.ptrvalue;
      if (cbp != NULL) {
        imp = cbp->imp;
        if (imp != NULL) {
          affil = imp->pub;
          if (affil != NULL && affil->choice == 1) {
            CleanVisStringJunkAndCompress (&(affil->affil));
          }
        }
      }
      break;
    case PUB_Patent :
      cpp = (CitPatPtr) vnp->data.ptrvalue;
      if (cpp != NULL) {
        if (StringCmp (cpp->country, "USA") == 0) {
          cpp->country = MemFree (cpp->country);
          cpp->country = StringSave ("US");
        }
      }
      break;
    default :
      break;
  }
  if (imp != NULL) {
    CleanVisStringAndCompress (&(imp->volume));
    CleanVisStringAndCompress (&(imp->issue));
    CleanVisStringAndCompress (&(imp->pages));
    CleanVisStringAndCompress (&(imp->section));
    CleanVisStringAndCompress (&(imp->part_sup));
    CleanVisStringAndCompress (&(imp->language));
    CleanVisStringAndCompress (&(imp->part_supi));
  }
}

//LCOV_EXCL_START
NLM_EXTERN void CleanUpPubdescAuthors (PubdescPtr pdp)

{
  Char             buf1 [121];
  Boolean          fixInitials = TRUE;
  Boolean          hasArt = FALSE;
  Boolean          hasUid = FALSE;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       vnp;

  if (pdp == NULL) return;
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      if (vnp->data.intvalue > 0) {
        hasUid = TRUE;
      }
    } else if (vnp->choice == PUB_Article) {
      hasArt = TRUE;
    }
  }
  if (hasArt && hasUid) {
    fixInitials = FALSE;
  }
  prev = &(pdp->pub);
  vnp = pdp->pub;
  while (vnp != NULL) {
    next = vnp->next;
    PubLabelUnique (vnp, buf1, sizeof (buf1) - 1, OM_LABEL_CONTENT, TRUE);
    NormalizePubAuthors (vnp, TRUE, fixInitials);
    vnp = next;
  }
}
//LCOV_EXCL_STOP

static int pub_order [] = {
  0,
  3,
  4,
  13,
  2,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  1
};

static int LIBCALLBACK SortByPubType (VoidPtr ptr1, VoidPtr ptr2)

{
  Uint1       chs1;
  Uint1       chs2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  chs1 = (Uint1) vnp1->choice;
  chs2 = (Uint1) vnp2->choice;
  if (chs1 < 14 && chs2 < 14) {
    chs1 = pub_order [chs1];
    chs2 = pub_order [chs2];
  }
  if (chs1 > chs2) {
    return 1;
  } else if (chs1 < chs2) {
    return -1;
  }
  return 0;
}

static void NormalizePubdesc (PubdescPtr pdp, Boolean stripSerial, Boolean doAuthors, ValNodePtr PNTR publist)

{
  ArticleIdPtr      aip;
  Int4              artpmid = 0;
  Char              buf1 [121];
  Char              buf2 [121];
  CitArtPtr         cap = NULL;
  CitGenPtr         cgp;
  CitJourPtr        cjp;
  Boolean           fixInitials = TRUE;
  Boolean           hasArt = FALSE;
  Boolean           hasUid = FALSE;
  ImprintPtr        imp;
  Int4              lastartpmid = 0;
  Int4              muid = 0;
  ValNodePtr        next;
  ArticleIdPtr      nextaip;
  Int4              pmid = 0;
  ValNodePtr PNTR   prev;
  ArticleIdPtr PNTR prevaip;
  ValNodePtr        vnp;

  if (pdp == NULL) return;
  CleanVisString (&(pdp->comment));
  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Muid) {
      if (vnp->data.intvalue > 0) {
        muid = vnp->data.intvalue;
      }
    }
    if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      if (vnp->data.intvalue > 0) {
        hasUid = TRUE;
      }
    } else if (vnp->choice == PUB_Article) {
      hasArt = TRUE;
    }
  }
  if (hasArt && hasUid) {
    fixInitials = FALSE;
  }
  if (pdp->pub != NULL) {
    pdp->pub = ValNodeSort (pdp->pub, SortByPubType);
  }

  /* remove zero muid where there is also a non-zero muid */
  prev = &(pdp->pub);
  vnp = pdp->pub;
  while (vnp != NULL) {
    next = vnp->next;
    if (vnp->choice == PUB_Muid && vnp->data.intvalue == 0 && muid != 0) {
      *prev = vnp->next;
      vnp->next = NULL;
      PubFree (vnp);
    } else {
      prev = &(vnp->next);
    }
    vnp = next;
  }

  prev = &(pdp->pub);
  vnp = pdp->pub;
  if (vnp != NULL && vnp->next == NULL && vnp->choice == PUB_Gen) {
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    buf1 [0] = '\0';
    PubLabelUnique (vnp, buf1, sizeof (buf1) - 1, OM_LABEL_CONTENT, TRUE);
    if (doAuthors) {
      NormalizeAuthors (cgp->authors, fixInitials);
    }
    if (stripSerial) {
      cgp->serial_number = -1;
    }
    if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
      cgp->cit [0] = 'U';
      /* cgp->date = DateFree (cgp->date); */ /* remove date if unpublished */
      if (cgp->journal == NULL) {
        cgp->volume = MemFree (cgp->volume);
        cgp->issue = MemFree (cgp->issue);
        cgp->pages = MemFree (cgp->pages);
      }
    }
    TrimSpacesAroundString (cgp->cit);
    if (StringDoesHaveText (cgp->title)) {
      StrStripSpaces (cgp->title);
    }
    buf2 [0] = '\0';
    PubLabelUnique (vnp, buf2, sizeof (buf2) - 1, OM_LABEL_CONTENT, TRUE);
    if (StringCmp (buf1, buf2) != 0) {
      ValNodeCopyStr (publist, 1, buf1);
      ValNodeCopyStr (publist, 2, buf2);
    }
    return; /* but does not remove if empty and only element of Pub */
  }
  while (vnp != NULL) {
    next = vnp->next;
    buf1 [0] = '\0';
    PubLabelUnique (vnp, buf1, sizeof (buf1) - 1, OM_LABEL_CONTENT, TRUE);
    if (doAuthors) {
      NormalizePubAuthors (vnp, stripSerial, fixInitials);
    }
    NormalizeAPub (vnp, stripSerial, fixInitials);
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
      if (cap != NULL && cap->from == 1) {
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp != NULL) {
          imp = cjp->imp;
          if (imp != NULL) {
            if (imp->pubstatus == PUBSTATUS_aheadofprint && imp->prepub != 2) {
              if (StringHasNoText (imp->volume) || StringHasNoText (imp->pages)) {
                imp->prepub = 2;
              }
            }
            if (imp->pubstatus == PUBSTATUS_aheadofprint && imp->prepub == 2) {
              if (StringDoesHaveText (imp->volume) && StringDoesHaveText (imp->pages)) {
                imp->prepub = 0;
              }
            }
            if (imp->pubstatus == PUBSTATUS_epublish && imp->prepub == 2) {
              imp->prepub = 0;
            }
          }
        }
      }
      if (cap != NULL) {
        aip = cap->ids;
        prevaip = (ArticleIdPtr PNTR) &(cap->ids);
        lastartpmid = 0;
        while (aip != NULL) {
          nextaip = aip->next;
          if (aip->choice == ARTICLEID_PUBMED) {
            artpmid = aip->data.intvalue;
            if (lastartpmid != 0 && lastartpmid == artpmid) {
              aip->next = NULL;
              *prevaip = nextaip;
              ArticleIdFree (aip);
            } else {
              prevaip = (ArticleIdPtr PNTR) &(aip->next);
            }
            lastartpmid = artpmid;
          } else {
            prevaip = (ArticleIdPtr PNTR) &(aip->next);
          }
          aip = nextaip;
        }
      }
    } else if (vnp->choice == PUB_PMid) {
      pmid = vnp->data.intvalue;
    }
    if (vnp->choice == PUB_Gen && empty_citgen ((CitGenPtr) vnp->data.ptrvalue)) {
      *prev = vnp->next;
      vnp->next = NULL;
      PubFree (vnp);
    } else {
      prev = &(vnp->next);
      buf2 [0] = '\0';
      PubLabelUnique (vnp, buf2, sizeof (buf2) - 1, OM_LABEL_CONTENT, TRUE);
      if (StringCmp (buf1, buf2) != 0) {
        ValNodeCopyStr (publist, 1, buf1);
        ValNodeCopyStr (publist, 2, buf2);
      }
    }
    vnp = next;
  }
  if (pmid == 0 && artpmid > 0) {
    ValNodeAddInt (&(pdp->pub), PUB_PMid, artpmid);
  } else if (pmid > 0 && artpmid == 0 && cap != NULL) {
    ValNodeAddInt (&(cap->ids), ARTICLEID_PUBMED, pmid);
  }
}

//LCOV_EXCL_START
NLM_EXTERN void CleanUpPubdescBody (PubdescPtr pdp, Boolean stripSerial)

{
  if (pdp == NULL) return;
  NormalizePubdesc (pdp, stripSerial, FALSE, NULL);
}
//LCOV_EXCL_STOP

static Boolean KeywordAlreadyInList (ValNodePtr head, CharPtr kwd)

{
  ValNodePtr  vnp;

  if (head == NULL || kwd == NULL) return FALSE;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (StringICmp ((CharPtr) vnp->data.ptrvalue, kwd) == 0) return TRUE;
  }

  return FALSE;
}

static Boolean CopyGeneXrefToGeneFeat (GeneRefPtr grp, GeneRefPtr grx)

{
  if (grp == NULL || grx == NULL) return FALSE;
  if (grx->db != NULL) {
    ValNodeLink (&(grp->db), grx->db);
    grx->db = NULL;
  }
  if (grx->locus == NULL && grx->allele == NULL &&
      grx->desc == NULL && grx->maploc == NULL &&
      grx->locus_tag == NULL && grx->db == NULL &&
      grx->syn == NULL) return TRUE;
  return FALSE;
}

static void HandleXrefOnGene (SeqFeatPtr sfp)

{
  GeneRefPtr           grp;
  GeneRefPtr           grx;
  SeqFeatXrefPtr       next;
  SeqFeatXrefPtr PNTR  prev;
  SeqFeatXrefPtr       xref;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_GENE) return;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;
   prev = &(sfp->xref);
  xref = sfp->xref;
  while (xref != NULL) {
    next = xref->next;
    if (xref->data.choice == SEQFEAT_GENE) {
      grx = (GeneRefPtr) xref->data.value.ptrvalue;
      if (CopyGeneXrefToGeneFeat (grp, grx)) {
        *(prev) = next;
        xref->next = NULL;
        SeqFeatXrefFree (xref);
      } else {
        prev = &(xref->next);
      }
    } else {
      prev = &(xref->next);
    }
    xref = next;
  }
}

static void CopyProtXrefToProtFeat (ProtRefPtr prp, ProtRefPtr prx)

{
  ValNodePtr       curr;
  size_t           len;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  CharPtr          str;

  if (prp == NULL || prx == NULL) return;

  if (prx->db != NULL) {
    ValNodeLink (&(prp->db), prx->db);
    prx->db = NULL;
  }

  prev = &(prx->name);
  curr = prx->name;
  while (curr != NULL) {
    next = curr->next;
    str = (CharPtr) curr->data.ptrvalue;
    if (! KeywordAlreadyInList (prp->name, str)) {
      ValNodeCopyStr (&(prp->name), 0, str);
      *(prev) = next;
      curr->next = NULL;
      curr->data.ptrvalue = NULL;
      ValNodeFree (curr);
    } else {
      prev = &(curr->next);
    }
    curr = next;
  }

  if (prp->desc == NULL) {
    prp->desc = prx->desc;
    prx->desc = NULL;
  } else if (prx->desc != NULL) {
    if (StringCmp (prx->desc, prp->desc) != 0) {
      len = StringLen (prp->desc) + StringLen (prx->desc) + 6;
      str = MemNew (len);
      if (str != NULL) {
        StringCpy (str, prp->desc);
        StringCat (str, "; ");
        StringCat (str, prx->desc);
        prp->desc = MemFree (prp->desc);
        prp->desc = str;
      }
    }
  }

  prev = &(prx->ec);
  curr = prx->ec;
  while (curr != NULL) {
    next = curr->next;
    str = (CharPtr) curr->data.ptrvalue;
    if (! KeywordAlreadyInList (prp->ec, str)) {
      ValNodeCopyStr (&(prp->ec), 0, str);
      *(prev) = next;
      curr->next = NULL;
      curr->data.ptrvalue = NULL;
      ValNodeFree (curr);
    } else {
      prev = &(curr->next);
    }
    curr = next;
  }

  prev = &(prx->activity);
  curr = prx->activity;
  while (curr != NULL) {
    next = curr->next;
    str = (CharPtr) curr->data.ptrvalue;
    if (! KeywordAlreadyInList (prp->activity, str)) {
      ValNodeCopyStr (&(prp->activity), 0, str);
      curr->data.ptrvalue = NULL;
    }
    *(prev) = next;
    curr->next = NULL;
    curr->data.ptrvalue = NULL;
    ValNodeFree (curr);
    curr = next;
  }
}

static Boolean InGpsGenomic (SeqFeatPtr sfp)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;

  if (sfp == NULL) return FALSE;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return FALSE;
  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    while (bssp != NULL) {
      if (bssp->_class == BioseqseqSet_class_nuc_prot) return FALSE;
      if (bssp->_class == BioseqseqSet_class_gen_prod_set) return TRUE;
      if (bssp->idx.parenttype != OBJ_BIOSEQSET) return FALSE;
      bssp = (BioseqSetPtr) bssp->idx.parentptr;
    }
  }
  return FALSE;
}

static void HandleXrefOnCDS (SeqFeatPtr sfp)

{
  SeqFeatXrefPtr       next;
  SeqFeatXrefPtr PNTR  prev;
  SeqFeatPtr           prot;
  ProtRefPtr           prp;
  ProtRefPtr           prx;
  SeqFeatXrefPtr       xref;

  if (sfp != NULL && sfp->product != NULL) {
    if (InGpsGenomic (sfp)) return;
    prot = GetBestProteinFeatureUnindexed (sfp->product);
    if (prot != NULL) {
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
      if (prp != NULL) {
        prev = &(sfp->xref);
        xref = sfp->xref;
        while (xref != NULL) {
          next = xref->next;
          if (xref->data.choice == SEQFEAT_PROT) {
            prx = (ProtRefPtr) xref->data.value.ptrvalue;
            CopyProtXrefToProtFeat (prp, prx);
            *(prev) = next;
            xref->next = NULL;
            SeqFeatXrefFree (xref);
          } else {
            prev = &(xref->next);
          }
          xref = next;
        }
      }
    }
  }
}

static void CleanUserStrings (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  CharPtr PNTR  cpp;
  Int4          i;
  ObjectIdPtr   oip;

  oip = ufp->label;
  if (oip != NULL && oip->str != NULL) {
    if (! StringHasNoText (oip->str)) {
      CleanVisString (&(oip->str));
    }
  }
  if (ufp->choice == 1) {
    if (! StringHasNoText ((CharPtr) ufp->data.ptrvalue)) {
      CleanVisStringAndCompress ((CharPtr PNTR) &(ufp->data.ptrvalue));
    }
  } else if (ufp->choice == 7) {
    cpp = (CharPtr PNTR) ufp->data.ptrvalue;
    if (cpp != NULL) {
      for (i = 0; i < ufp->num; i++) {
        TrimSpacesSemicolonsAndCommas (cpp [i]);
        Asn2gnbkCompressSpaces (cpp [i]);
      }
    }
  }
}

static void CleanUserFields (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  oip = ufp->label;
  if (oip != NULL && oip->str != NULL) {
    if (! StringHasNoText (oip->str)) {
      CleanVisString (&(oip->str));
    }
  }
  VisitUserFieldsInUfp (ufp, userdata, CleanUserStrings);
}

//LCOV_EXCL_START
NLM_EXTERN UserFieldPtr LIBCALL UserFieldSort (UserFieldPtr list, int (LIBCALLBACK *compar ) PROTO((VoidPtr, VoidPtr)))

{
  Int4          count, i;
  UserFieldPtr  PNTR head;
  UserFieldPtr  tmp;

  if (list == NULL) return NULL;

  count = 0;
  for (tmp = list; tmp != NULL; tmp = tmp->next) {
    count++;
  }

  head = (UserFieldPtr *) MemNew (((size_t) count + 1) * sizeof (UserFieldPtr));

  for (tmp = list, i = 0; tmp != NULL && i < count; i++) {
    head [i] = tmp;
    tmp = tmp->next;
  }

  HeapSort (head, (size_t) count, sizeof (UserFieldPtr), compar);

  for (i = 0; i < count; i++) {
    tmp = head [i];
    tmp->next = head [i + 1];
  }
  list = head [0];

  MemFree (head);

  return list;
}
//LCOV_EXCL_STOP

/*
static CharPtr barcodeOrder [] = {
  "",
  "StructuredCommentPrefix",
  "Barcode Index Number",
  "Order Assignment",
  "iBOL Working Group",
  "iBOL Release Status",
  "Tentative Name",
  "StructuredCommentSuffix",
  NULL
};

static Int2 GetBarcodeOrder (CharPtr str)

{
  Int2  i;

  if (StringHasNoText (str)) return 0;

  for (i = 1; barcodeOrder [i] != NULL; i++) {
    if (StringCmp (str, barcodeOrder [i]) == 0) return i;
  }

  return 0;
}

static int LIBCALLBACK ReorderBarcodeFields (VoidPtr ptr1, VoidPtr ptr2)

{
  Int2          idx1, idx2;
  ObjectIdPtr   lbl1, lbl2;
  CharPtr       str1, str2;
  UserFieldPtr  ufp1, ufp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  ufp1 = *((UserFieldPtr PNTR) ptr1);
  ufp2 = *((UserFieldPtr PNTR) ptr2);
  if (ufp1 == NULL || ufp2 == NULL) return 0;

  lbl1 = (ObjectIdPtr) ufp1->label;
  lbl2 = (ObjectIdPtr) ufp2->label;
  if (lbl1 == NULL || lbl2 == NULL) return 0;

  str1 = (CharPtr) lbl1->str;
  str2 = (CharPtr) lbl2->str;
  if (str1 == NULL || str2 == NULL) return 0;

  idx1 = GetBarcodeOrder (str1);
  idx2 = GetBarcodeOrder (str2);

  if (idx1 > idx2) return 1;
  if (idx1 < idx2) return -1;

  return 0;
}
*/

NLM_EXTERN void CleanStructuredComment (
  UserObjectPtr uop
)

{
  Boolean      genome_assembly_data = FALSE, ibol_data = FALSE;
  UserFieldPtr ufp;
  CharPtr      str, core, new_str;

  if (uop == NULL || uop->type == NULL 
      || StringCmp (uop->type->str, "StructuredComment") != 0) {
    return;
  }

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->label != NULL 
        && ufp->choice == 1 
        && (str = (CharPtr) ufp->data.ptrvalue) != NULL) {
      if (StringCmp (ufp->label->str, "StructuredCommentPrefix") == 0) {
        core = StructuredCommentDbnameFromString(str);
        new_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (core) + 15));
        sprintf (new_str, "##%s-START##", core);
        str = MemFree (str);
        ufp->data.ptrvalue = new_str;
        if (StringCmp (core, "Genome-Assembly-Data") == 0) {
          genome_assembly_data = TRUE;
        } else if (StringCmp (core, "International Barcode of Life (iBOL)Data") == 0) {
          ibol_data = TRUE;
        }
        core = MemFree (core);
      } else if (StringCmp (ufp->label->str, "StructuredCommentSuffix") == 0) {
        core = StructuredCommentDbnameFromString(str);
        new_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (core) + 15));
        sprintf (new_str, "##%s-END##", core);
        str = MemFree (str);
        ufp->data.ptrvalue = new_str;
        if (StringCmp (core, "Genome-Assembly-Data") == 0) {
          genome_assembly_data = TRUE;
        } else if (StringCmp (core, "International Barcode of Life (iBOL)Data") == 0) {
          ibol_data = TRUE;
        }
        core = MemFree (core);
      }
    }
  }

  if (genome_assembly_data) {
    for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
      if (ufp->label != NULL 
          && ufp->choice == 1 
          && (str = (CharPtr) ufp->data.ptrvalue) != NULL) {
        if (StringCmp (ufp->label->str, "Finishing Goal") == 0 ||
            StringCmp (ufp->label->str, "Current Finishing Status") == 0) {
          if (StringCmp (str, "High Quality Draft") == 0) {
            ufp->data.ptrvalue = StringSave ("High-Quality Draft");
            str = MemFree (str);
          } else if (StringCmp (str, "Improved High Quality Draft") == 0) {
            ufp->data.ptrvalue = StringSave ("Improved High-Quality Draft");
            str = MemFree (str);
          } else if (StringCmp (str, "Annotation Directed") == 0) {
            ufp->data.ptrvalue = StringSave ("Annotation-Directed Improvement");
            str = MemFree (str);
          } else if (StringCmp (str, "Non-contiguous Finished") == 0) {
            ufp->data.ptrvalue = StringSave ("Noncontiguous Finished");
            str = MemFree (str);
          }
        } else if (StringCmp(ufp->label->str, "Assembly Date") == 0) {
          str = (CharPtr) ufp->data.ptrvalue;
          ReformatAssemblyDate(&str);
          ufp->data.ptrvalue = str;
        }
      }
    }
  }

  if (ibol_data) {
    /*
    uop->data = UserFieldSort (uop->data, ReorderBarcodeFields);
    */
    ReorderStructuredCommentFields (uop);
  }
}


//LCOV_EXCL_START
// change made as a result of SQD-2399, which will not be implemented for the C++ Toolkit
// going forward. bad data was generated internally, production process has been fixed.
static void CleanRefGeneTrackingUserObject (
  UserObjectPtr uop
)

{
  UserFieldPtr  asmbly = NULL, entry, tmp, ufp;
  ObjectIdPtr   oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "RefGeneTracking") != 0) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL) continue;
    if (StringCmp (oip->str, "Assembly") != 0) continue;
    asmbly = ufp;
    break;
  }

  if (asmbly == NULL || asmbly->choice != 11) return;
  tmp = asmbly->data.ptrvalue;
  if (tmp == NULL || tmp->choice == 11) return;

  entry = UserFieldNew ();
  if (entry == NULL) return;
  oip = ObjectIdNew ();
  if (oip == NULL) return;

  entry->data.ptrvalue = (Pointer) tmp;
  entry->choice = 11;
  entry->label = oip;
  oip->id = 0;

  asmbly->data.ptrvalue = (Pointer) entry;
  asmbly->choice = 11;
}
//LCOV_EXCL_STOP

static void CleanUserObject (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  oip = uop->type;
  if (oip != NULL && oip->str != NULL) {
    if (! StringHasNoText (oip->str)) {
      CleanVisString (&(oip->str));
    }
  }
  VisitUserFieldsInUop (uop, userdata, CleanUserFields);
  CleanStructuredComment (uop);
  CleanRefGeneTrackingUserObject (uop);
}

static CharPtr bsecSiteList [] = {
  "", "active", "binding", "cleavage", "inhibit", "modifi",
  "glycosylation", "myristoylation", "mutagenized", "metal-binding",
  "phosphorylation", "acetylation", "amidation", "methylation",
  "hydroxylation", "sulfatation", "oxidative-deamination",
  "pyrrolidone-carboxylic-acid", "gamma-carboxyglutamic-acid",
  "blocked", "lipid-binding", "np-binding", "DNA-binding",
  "signal-peptide", "transit-peptide", "transmembrane-region",
  "nitrosylation", NULL
};

static CharPtr uninfStrings [] = {
  "signal",
  "transit",
  "peptide",
  "signal peptide",
  "signal-peptide",
  "signal_peptide",
  "transit peptide",
  "transit-peptide",
  "transit_peptide",
  "unnamed",
  "unknown",
  "putative",
  NULL
};

static Boolean InformativeString (CharPtr str)

{
  Int2  i;

  if (StringHasNoText (str)) return FALSE;

  for (i = 0; uninfStrings [i] != NULL; i++) {
    if (StringICmp (str, uninfStrings [i]) == 0) return FALSE;
  }

  return TRUE;
}

static void CleanUpExceptText (SeqFeatPtr sfp)

{
  ValNodePtr  head, vnp;
  size_t      len;
  CharPtr     prefix, ptr, str, tmp;

  if (sfp == NULL || sfp->except_text == NULL) return;
  if (StringStr (sfp->except_text, "ribosome slippage") == NULL &&
      StringStr (sfp->except_text, "trans splicing") == NULL &&
      StringStr (sfp->except_text, "alternate processing") == NULL &&
      StringStr (sfp->except_text, "non-consensus splice site") == NULL &&
      StringStr (sfp->except_text, "adjusted for low quality genome") == NULL) return;

  head = NULL;
  str = sfp->except_text;
  tmp = str;
  while (! StringHasNoText (tmp)) {
    ptr = StringChr (tmp, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (tmp);
    ValNodeCopyStr (&head, 0, tmp);
    tmp = ptr;
  }
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    tmp = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (tmp)) continue;
    if (StringCmp (tmp, "ribosome slippage") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("ribosomal slippage");
    } else if (StringCmp (tmp, "trans splicing") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("trans-splicing");
    } else if (StringCmp (tmp, "alternate processing") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("alternative processing");
    } else if (StringCmp (tmp, "non-consensus splice site") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("nonconsensus splice site");
    } else if (StringCmp (tmp, "adjusted for low quality genome") == 0) {
      vnp->data.ptrvalue = MemFree (tmp);
      vnp->data.ptrvalue = StringSave ("adjusted for low-quality genome");
    }
  }

  len = 0;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    tmp = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (tmp)) continue;
    len += StringLen (tmp) + 2;
  }

  str = (CharPtr) MemNew (len + 2);
  if (str == NULL) return;

  prefix = "";
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    tmp = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (tmp)) continue;
    StringCat (str, prefix);
    StringCat (str, tmp);
    prefix = ", ";
  }

  sfp->except_text = MemFree (sfp->except_text);
  sfp->except_text = str;

  ValNodeFreeData (head);
}

static Boolean ExpandGeneSynCom (ValNodePtr headsyn)

{
  ValNodePtr  lastsyn;
  ValNodePtr  newsyn;
  ValNodePtr  nextsyn;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;

  str = (CharPtr) headsyn->data.ptrvalue;
  if (StringHasNoText (str)) return TRUE;
  if (StringChr (str, ',') == NULL) return FALSE;

  nextsyn = headsyn->next;
  lastsyn = headsyn;
  tmp = StringSave ((CharPtr) headsyn->data.ptrvalue);
  str = tmp;

  while (! StringHasNoText (str)) {
    ptr = StringChr (str, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    newsyn = ValNodeNew (NULL);
    if (newsyn != NULL) {
      newsyn->data.ptrvalue = StringSave (str);
      newsyn->next = nextsyn;
      lastsyn->next = newsyn;
      lastsyn = newsyn;
    }
    str = ptr;
  }

  MemFree (tmp);
  return TRUE;
}

static Boolean ExpandGeneSynSem (ValNodePtr headsyn)

{
  ValNodePtr  lastsyn;
  ValNodePtr  newsyn;
  ValNodePtr  nextsyn;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;

  str = (CharPtr) headsyn->data.ptrvalue;
  if (StringHasNoText (str)) return TRUE;
  if (StringStr (str, "; ") == NULL) return FALSE;

  nextsyn = headsyn->next;
  lastsyn = headsyn;
  tmp = StringSave ((CharPtr) headsyn->data.ptrvalue);
  str = tmp;

  while (! StringHasNoText (str)) {
    ptr = StringStr (str, "; ");
    if (ptr != NULL) {
      ptr++;
      *ptr = '\0';
      ptr++;
    }
    TrimSpacesAroundString (str);
    newsyn = ValNodeNew (NULL);
    if (newsyn != NULL) {
      newsyn->data.ptrvalue = StringSave (str);
      newsyn->next = nextsyn;
      lastsyn->next = newsyn;
      lastsyn = newsyn;
    }
    str = ptr;
  }

  MemFree (tmp);
  return TRUE;
}

static void ExpandGeneSynList (GeneRefPtr grp)

{
  ValNodePtr       currsyn;
  ValNodePtr       nextsyn;
  ValNodePtr PNTR  prevsyn;

  if (grp == NULL || grp->syn == NULL) return;

  currsyn = grp->syn;
  prevsyn = &(grp->syn);
  while (currsyn != NULL) {
    if (ExpandGeneSynCom (currsyn)) {
      nextsyn = currsyn->next;
      *(prevsyn) = currsyn->next;
      currsyn->next = NULL;
      ValNodeFreeData (currsyn);
    } else {
      nextsyn = currsyn->next;
      prevsyn = (ValNodePtr PNTR) &(currsyn->next);
    }
    currsyn = nextsyn;
  }

  currsyn = grp->syn;
  prevsyn = &(grp->syn);
  while (currsyn != NULL) {
    if (ExpandGeneSynSem (currsyn)) {
      nextsyn = currsyn->next;
      *(prevsyn) = currsyn->next;
      currsyn->next = NULL;
      ValNodeFreeData (currsyn);
    } else {
      nextsyn = currsyn->next;
      prevsyn = (ValNodePtr PNTR) &(currsyn->next);
    }
    currsyn = nextsyn;
  }
}

typedef struct gosstruc {
  CharPtr       term;
  Char          goid [32];
  CharPtr       evidence;
  Int4          pmid;
  CharPtr       goref;
  UserFieldPtr  ufp;
} GosStruc, PNTR GosStrucPtr;

static int LIBCALLBACK SortVnpByGssp (VoidPtr ptr1, VoidPtr ptr2)

{
  int           compare;
  GosStrucPtr   gsp1, gsp2;
  ValNodePtr    vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  gsp1 = (GosStrucPtr) vnp1->data.ptrvalue;
  gsp2 = (GosStrucPtr) vnp2->data.ptrvalue;
  if (gsp1 == NULL || gsp2 == NULL) return 0;

  compare = StringICmp (gsp1->goid, gsp2->goid);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->term, gsp2->term);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  compare = StringICmp (gsp1->evidence, gsp2->evidence);
  if (compare > 0) {
    return 1;
  } else if (compare < 0) {
    return -1;
  }

  if (gsp1->pmid == 0) return 1;
  if (gsp2->pmid == 0) return -1;
  if (gsp1->pmid > gsp2->pmid) {
    return 1;
  } else if (gsp1->pmid < gsp2->pmid) {
    return -1;
  }

  return 0;
}

static CharPtr bsecGoQualType [] = {
  "", "Process", "Component", "Function", NULL
};

static CharPtr bsecGoFieldType [] = {
  "", "text string", "go id", "pubmed id", "go ref", "evidence", NULL
};

static UserFieldPtr SortGoTerms (
  UserFieldPtr entryhead
)

{
  UserFieldPtr  entry, topufp, ufp, lastufp;
  CharPtr       evidence, goid, goref, textstr;
  Char          gid [32];
  GosStrucPtr   gsp, lastgsp;
  ValNodePtr    head = NULL, vnp;
  Int2          j;
  ObjectIdPtr   oip;
  Int4          pmid;

  if (entryhead == NULL) return entryhead;

  for (entry = entryhead; entry != NULL; entry = entry->next) {
    if (entry == NULL || entry->choice != 11) break;
    topufp = (UserFieldPtr)  entry->data.ptrvalue;
    if (topufp == NULL) continue;

    textstr = NULL;
    evidence = NULL;
    goid = NULL;
    goref = NULL;
    pmid = 0;
    for (ufp = topufp; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      for (j = 0; bsecGoFieldType [j] != NULL; j++) {
        if (StringICmp (oip->str, bsecGoFieldType [j]) == 0) break;
      }
      if (bsecGoFieldType [j] == NULL) continue;
      switch (j) {
        case 1 :
          if (ufp->choice == 1) {
            textstr = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 2 :
          if (ufp->choice == 1) {
            goid = (CharPtr) ufp->data.ptrvalue;
          } else if (ufp->choice == 2) {
            sprintf (gid, "%ld", (long) (Int4) ufp->data.intvalue);
            goid = (CharPtr) gid;
          }
          break;
        case 3 :
          if (ufp->choice == 2) {
            pmid = (Int4) ufp->data.intvalue;
          }
          break;
        case 4 :
          if (ufp->choice == 1) {
            goref = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        case 5 :
          if (ufp->choice == 1) {
            evidence = (CharPtr) ufp->data.ptrvalue;
          }
          break;
        default :
          break;
      }
    }

    if (StringDoesHaveText (textstr)) {
      gsp = (GosStrucPtr) MemNew (sizeof (GosStruc));
      if (gsp != NULL) {
        gsp->term = textstr;
        StringNCpy_0 (gsp->goid, goid, sizeof (gsp->goid));
        gsp->evidence = evidence;
        gsp->pmid = pmid;
        gsp->goref = goref;
        gsp->ufp = entry;
        ValNodeAddPointer (&head, 0, (Pointer) gsp);
      }
    }
  }

  if (head == NULL) return entryhead;
  head = ValNodeSort (head, SortVnpByGssp);

  entryhead = NULL;
  lastgsp = NULL;
  lastufp = NULL;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    gsp = (GosStrucPtr) vnp->data.ptrvalue;
    if (gsp == NULL || gsp->ufp == NULL) continue;
    if (lastgsp != NULL &&
        (StringICmp (gsp->term, lastgsp->term) == 0 || StringICmp (gsp->goid, lastgsp->goid) == 0) &&
         (gsp->pmid == lastgsp->pmid &&
          StringICmp (gsp->goref, lastgsp->goref) == 0 &&
          StringICmp (gsp->evidence, lastgsp->evidence) == 0)) {
      gsp->ufp->next = NULL;
      UserFieldFree (gsp->ufp);
    } else {
      if (lastufp != NULL) {
        lastufp->next = gsp->ufp;
      } else {
        entryhead = gsp->ufp;
      }
      lastufp = gsp->ufp;
      lastufp->next = NULL;
    }
    lastgsp = gsp;
  }

  ValNodeFreeData (head);

  return entryhead;
}

static void SortGoTermsUfp (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr  entry;
  Int2          i;
  ObjectIdPtr   oip;
 
  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; bsecGoQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, bsecGoQualType [i]) == 0) break;
  }
  if (bsecGoQualType [i] == NULL) return;

  entry = ufp->data.ptrvalue;
  if (entry == NULL || entry->choice != 11) return;

  ufp->data.ptrvalue = SortGoTerms (entry);
}

static void SortGoTermsSfp (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    VisitUserFieldsInUop (uop, userdata, SortGoTermsUfp);
  }
}

static void CleanupGoTerms (
  UserFieldPtr entryhead
)

{
  UserFieldPtr  entry, topufp, ufp;
  CharPtr       goid, goref, str;
  Int2          j;
  ObjectIdPtr   oip;

  if (entryhead == NULL) return;

  for (entry = entryhead; entry != NULL; entry = entry->next) {
    if (entry == NULL || entry->choice != 11) break;
    topufp = (UserFieldPtr)  entry->data.ptrvalue;
    if (topufp == NULL) continue;

    goid = NULL;
    goref = NULL;
    for (ufp = topufp; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      for (j = 0; bsecGoFieldType [j] != NULL; j++) {
        if (StringICmp (oip->str, bsecGoFieldType [j]) == 0) break;
      }
      if (bsecGoFieldType [j] == NULL) continue;
      switch (j) {
        case 2 :
          if (ufp->choice == 1) {
            goid = (CharPtr) ufp->data.ptrvalue;
            if (goid != NULL && *goid != '\0') {
              if (StringNICmp (goid, "GO:", 3) == 0) {
                str = StringSave (goid + 3);
                ufp->data.ptrvalue = (Pointer) str;
                MemFree (goid);
              }
            }
          }
          break;
        case 4 :
          if (ufp->choice == 1) {
            goref = (CharPtr) ufp->data.ptrvalue;
            if (goref != NULL && *goref != '\0') {
              if (StringNICmp (goref, "GO_REF:", 7) == 0) {
                str = StringSave (goref + 7);
                ufp->data.ptrvalue = (Pointer) str;
                MemFree (goref);
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

static void CleanupGoTermsUfp (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr  entry;
  Int2          i;
  ObjectIdPtr   oip;
 
  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; bsecGoQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, bsecGoQualType [i]) == 0) break;
  }
  if (bsecGoQualType [i] == NULL) return;

  entry = ufp->data.ptrvalue;
  if (entry == NULL || entry->choice != 11) return;

  CleanupGoTerms (entry);
}

static void CleanupGoTermsSfp (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "GeneOntology") == 0) {
    VisitUserFieldsInUop (uop, userdata, CleanupGoTermsUfp);
  }
}

static CharPtr CleanUpSgml (
  CharPtr str
)

{
  Int2     ascii_len;
  Char     buf [256];
  CharPtr  ptr;

  if (StringHasNoText (str)) return NULL;
  if (StringChr (str, '&') == NULL) return NULL;

  ascii_len = Sgml2AsciiLen (str);
  if (ascii_len + 2 >= sizeof (buf)) return NULL;

  buf [0] = '\0';
  Sgml2Ascii (str, buf, ascii_len + 1);
  if (StringHasNoText (buf)) return NULL;
  if (StringCmp (str, buf) == 0) return NULL;

  ptr = StringChr (buf, '<');
  if (ptr != NULL) {
    *ptr = ' ';
  }
  ptr = StringChr (buf, '>');
  if (ptr != NULL) {
    *ptr = ' ';
  }
  TrimSpacesAroundString (buf);
  Asn2gnbkCompressSpaces (buf);

  return StringSave (buf);
}

/* special exception for genome pipeline rRNA names */

static Boolean NotExceptedRibosomalName (
  CharPtr name
)

{
  Char     ch;
  CharPtr  str;

  str = StringStr (name, " ribosomal");
  if (str == NULL) return FALSE;

  str += 10;
  ch = *str;
  while (ch != '\0') {
    if (ch == ' ' || IS_DIGIT (ch)) {
      /* okay */
    } else {
      return TRUE;
    }
    str++;
    ch = *str;
  }

  return FALSE;
}

//LCOV_EXCL_START
NLM_EXTERN void CleanupSubSourceOrgModOtherFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioSourcePtr  biop;
  OrgNamePtr    onp = NULL;
  OrgRefPtr     orp;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_BIOSRC) return;
  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  if (biop == NULL) return;
  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (orp != NULL) {
      CleanupOrgModOther (biop, onp);
    }
  }
  CleanupSubSourceOther (biop, onp);
}

NLM_EXTERN void CleanupSubSourceOrgModOtherDesc (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BioSourcePtr  biop;
  OrgNamePtr    onp = NULL;
  OrgRefPtr     orp;

  if (sdp == NULL) return;
  if (sdp->choice != Seq_descr_source) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (orp != NULL) {
      CleanupOrgModOther (biop, onp);
    }
  }
  CleanupSubSourceOther (biop, onp);
}
//LCOV_EXCL_STOP


typedef struct xmltable {
  CharPtr  code;
  size_t   len;
  CharPtr  letter;
} XmlTable, PNTR XmlTablePtr;

static XmlTable xmlunicodes [] = {
  { "&amp",     4, "&"},
  { "&apos",    5, "\'"},
  { "&gt",      3, ">"},
  { "&lt",      3, "<"},
  { "&quot",    5, "\""},
  { "&#13&#10", 8, ""},
  { "&#916",    5, "Delta"},
  { "&#945",    5, "alpha"},
  { "&#946",    5, "beta"},
  { "&#947",    5, "gamma"},
  { "&#952",    5, "theta"},
  { "&#955",    5, "lambda"},
  { "&#956",    5, "mu"},
  { "&#957",    5, "nu"},
  { "&#8201",   6, " "},
  { "&#8206",   6, ""},
  { "&#8242",   6, "'"},
  { "&#8594",   6, "->"},
  { "&#8722",   6, "-"},
  { "&#8710",   6, "delta"},
  { "&#64257",  7, "fi"},
  { "&#64258",  7, "fl"},
  { "&#65292",  7, ","},
  { NULL,       0, ""}
};

static CharPtr BSECDecodeXml (
  CharPtr str
)

{
  Char         ch, nxt;
  CharPtr      dst, ptr, src;
  Int2         i;
  size_t       len;
  XmlTablePtr  xtp;

  if (StringHasNoText (str)) return str;

  src = str;
  dst = str;
  ch = *src;
  while (ch != '\0') {
    if (ch == '&') {
      xtp = NULL;
      len = 1;
      for (i = 0; xmlunicodes [i].code != NULL; i++) {
        if (StringNICmp (src, xmlunicodes [i].code, xmlunicodes [i].len) == 0) {
          nxt = *(src +xmlunicodes [i].len);
          if (nxt == ';') {
            xtp = &(xmlunicodes [i]);
            len = xtp->len + 1;
            break;
          } else if (nxt == ' ' || nxt == '\0') {
            xtp = &(xmlunicodes [i]);
            len = xtp->len;
            break;
          }
        }
      }
      if (xtp != NULL) {
        if (StringLen (xtp->letter) > 0) {
          ptr = xtp->letter;
          ch = *ptr;
          while (ch != '\0') {
            *dst = ch;
            dst++;
            ptr++;
            ch = *ptr;
          }
        }
        src += len;
      } else {
        *dst = ch;
        dst++;
        src++;
      }
    } else {
      *dst = ch;
      dst++;
      src++;
    }
    ch = *src;
  }
  *dst = '\0';

  return str;
}

static void CleanupFeatureStrings (
  SeqFeatPtr sfp,
  Boolean isJscan,
  Boolean isEmblOrDdbj,
  Boolean stripSerial,
  Boolean modernizeFeats,
  ValNodePtr PNTR publist
)

{
  Uint1                aa;
  BioSourcePtr         biop;
  Char                 ch;
  Uint1                codon [6];
  GeneNomenclaturePtr  gnp;
  GeneRefPtr           grp;
  ImpFeatPtr           ifp;
  Boolean              is_fMet = FALSE;
  Boolean              is_iMet = FALSE;
  Int2                 j;
  Boolean              justTrnaText;
  size_t               len;
  CharPtr              name;
  ObjectIdPtr          oip;
  OrgNamePtr           onp = NULL;
  OrgRefPtr            orp;
  PubdescPtr           pdp;
  ProtRefPtr           prp;
  CharPtr              ptr;
  RNAGenPtr            rgp;
  RNAQualPtr           rqp;
  RnaRefPtr            rrp;
  SubSourcePtr         ssp;
  CharPtr              str;
  CharPtr              suff;
  CharPtr              temp;
  Char                 tmp [64];
  Boolean              trimming_junk;
  tRNAPtr              trp;
  UserFieldPtr         ufp;
  UserObjectPtr        uop;
  CharPtr              val;
  ValNodePtr           vnp, vnp2;
  SeqFeatXrefPtr       xref;

  if (sfp == NULL) return;
  BSECDecodeXml (sfp->comment);
  CleanVisStringAndCompress (&(sfp->comment));
  len = StringLen (sfp->comment);
  if (len > 4) {
    if (StringCmp (sfp->comment + len - 3, ",..") == 0 ||
        StringCmp (sfp->comment + len - 3, ".,.") == 0 ||
        StringCmp (sfp->comment + len - 3, "..,") == 0 ||
        StringCmp (sfp->comment + len - 3, ",.,") == 0) {
      sfp->comment [len - 3] = '.';
      sfp->comment [len - 2] = '.';
      sfp->comment [len - 1] = '.';
    }
  }
  BSECDecodeXml (sfp->title);
  CleanVisString (&(sfp->title));
  CleanVisString (&(sfp->except_text));
  if (StringDoesHaveText (sfp->except_text)) {
    CleanUpExceptText (sfp);
  }
  CleanDoubleQuote (sfp->comment);
  if (StringCmp (sfp->comment, ".") == 0) {
    sfp->comment = MemFree (sfp->comment);
  }
  /*
  if (sfp->ext != NULL) {
    VisitUserObjectsInUop (sfp->ext, NULL, SortGoTermsSfp);
  }
  */
  if (sfp->ext != NULL) {
    VisitUserObjectsInUop (sfp->ext, NULL, CleanupGoTermsSfp);
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->data.choice != SEQFEAT_PROT) continue;
    prp = (ProtRefPtr) xref->data.value.ptrvalue;
    if (prp == NULL) continue;
    RemoveFlankingQuotes (&(prp->desc));
    RemoveFlankingQuotesList (&(prp->name));
    CleanVisStringAndCompress (&(prp->desc));
    CleanVisStringListAndCompress (&(prp->name));
  }

  switch (sfp->data.choice) {
    case SEQFEAT_BOND :
    case SEQFEAT_PSEC_STR :
    case SEQFEAT_COMMENT:
      return;
    case SEQFEAT_SITE :
      for (j = 0; bsecSiteList [j] != NULL; j++) {
        StringNCpy_0 (tmp, bsecSiteList [j], sizeof (tmp));
        len = StringLen (tmp);
        if (StringNICmp (sfp->comment, tmp, len) == 0) {
          if (sfp->data.value.intvalue == 0 || sfp->data.value.intvalue == 255) {
            sfp->data.value.intvalue = j;
            if (StringHasNoText (sfp->comment + len) || StringICmp (sfp->comment + len, " site") == 0) {
              sfp->comment = MemFree (sfp->comment);
            }
          }
        } else {
          val = tmp;
          ch = *val;
          while (ch != '\0') {
            if (ch == '-') {
              *val = ' ';
            }
            val++;
            ch = *val;
          }
          if (StringNICmp (sfp->comment, tmp, len) == 0) {
            if (sfp->data.value.intvalue == 0 || sfp->data.value.intvalue == 255) {
              sfp->data.value.intvalue = j;
              if (StringHasNoText (sfp->comment + len) || StringICmp (sfp->comment + len, " site") == 0) {
                sfp->comment = MemFree (sfp->comment);
              }
            }
          }
        }
      }
      break;
    default :
      break;
  }
  if (sfp->data.value.ptrvalue == NULL) return;

  biop = NULL;
  orp = NULL;
  switch (sfp->data.choice) {
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
      }
    default :
      break;
  }
  if (orp != NULL && sfp->qual != NULL) {
    GbqualToOrpMod (&(sfp->qual), &(orp->mod));
  }

  biop = NULL;
  orp = NULL;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (sfp->xref != NULL) {
        HandleXrefOnGene (sfp);
      }
      BSECDecodeXml (grp->locus);
      CleanVisStringAndCompress (&(grp->locus));
      /*
      if (isJscan && StringDoesHaveText (grp->locus)) {
        ptr = CleanUpSgml (grp->locus);
        if (ptr != NULL) {
          grp->locus = MemFree (grp->locus);
          grp->locus = StringSave (ptr);
        }
      }
      */
      CleanVisString (&(grp->allele));
      CleanVisStringAndCompress (&(grp->desc));
      CleanVisString (&(grp->maploc));
      CleanVisString (&(grp->locus_tag));
      ExpandGeneSynList (grp);
      /*
      if (isJscan && grp->syn != NULL) {
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          ptr = CleanUpSgml (str);
          if (ptr != NULL) {
            vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
            vnp->data.ptrvalue = StringSave (ptr);
          }
        }
      }
      */
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        BSECDecodeXml (str);
      }
      CleanVisStringListCaseSensitive (&(grp->syn));
      grp->syn = ValNodeSort (grp->syn, SortVnpByStringCS);
      grp->syn = UniqueStringValNodeCS (grp->syn);
      grp->syn = ValNodeSort (grp->syn, SortVnpByStringCILCFirst);
      CleanDoubleQuote (grp->locus);
      CleanDoubleQuote (grp->allele);
      CleanDoubleQuote (grp->desc);
      /*
      if (isJscan && StringDoesHaveText (grp->desc)) {
        ptr = CleanUpSgml (grp->desc);
        if (ptr != NULL) {
          grp->desc = MemFree (grp->desc);
          grp->desc = StringSave (ptr);
        }
      }
      */
      CleanDoubleQuote (grp->maploc);
      CleanDoubleQuote (grp->locus_tag);
      CleanDoubleQuoteList (grp->syn);
      FixOldDbxrefs (grp->db, isEmblOrDdbj);
      FixNumericDbxrefs (grp->db);
      grp->db = ValNodeSort (grp->db, SortDbxref);
      CleanupDuplicateDbxrefs (&(grp->db));
      CleanupObsoleteDbxrefs (&(grp->db));
      CleanupGoDbxrefs (grp->db);
      /* now move grp->dbxref to sfp->dbxref */
      vnp = grp->db;
      grp->db = NULL;
      ValNodeLink ((&sfp->dbxref), vnp);
      if (grp->locus != NULL && grp->syn != NULL) {
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringCmp (grp->locus, str) == 0) {
            vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
          }
        }
        CleanVisStringListCaseSensitive (&(grp->syn));
      }
      gnp = grp->formal_name;
      if (gnp != NULL) {
        FixOldDbxref (gnp->source);
        FixNumericDbxref (gnp->source);
      }
      /*
      if (grp->locus != NULL && sfp->comment != NULL) {
        if (StringCmp (grp->locus, sfp->comment) == 0) {
          sfp->comment = MemFree (sfp->comment);
        }
      }
      */
      break;
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      break;
    case SEQFEAT_CDREGION :
      if (sfp->xref != NULL && sfp->product != NULL) {
        HandleXrefOnCDS (sfp);
      }
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        CleanupECNumber (str);
        if (ECNumberCanBeSplit (str)) {
          ptr = str;
          ch = *ptr;
          while (ch != '\0' && ch != ' ' && ch != ';') {
            ptr++;
            ch = *ptr;
          }
          if (ch != '\0') {
            *ptr = '\0';
            ptr++;
            vnp2 = ValNodeCopyStr (NULL, 0, ptr);
            if (vnp2 != NULL) {
              vnp2->next = vnp->next;
              vnp->next = vnp2;
            }
          }
        }
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        BSECDecodeXml (str);
      }
      BSECDecodeXml (prp->desc);
      CleanVisStringAndCompress (&(prp->desc));
      CleanVisStringJunkListAndCompress (&(prp->name));
      CleanVisStringList (&(prp->ec));
      CleanVisStringJunkListAndCompress (&(prp->activity));
      CleanDoubleQuote (prp->desc);
      CleanDoubleQuoteList (prp->name);
      CleanDoubleQuoteList (prp->ec);
      CleanDoubleQuoteList (prp->activity);
      RemoveFlankingQuotes (&(prp->desc));
      RemoveFlankingQuotesList (&(prp->name));
      FixOldDbxrefs (prp->db, isEmblOrDdbj);
      FixNumericDbxrefs (prp->db);
      prp->db = ValNodeSort (prp->db, SortDbxref);
      CleanupDuplicateDbxrefs (&(prp->db));
      CleanupObsoleteDbxrefs (&(prp->db));
      CleanupGoDbxrefs (prp->db);
      /* now move prp->dbxref to sfp->dbxref */
      vnp = prp->db;
      prp->db = NULL;
      ValNodeLink ((&sfp->dbxref), vnp);
      if (prp->processed != 3 && prp->processed != 4 && prp->processed != 5 &&
          prp->name == NULL && sfp->comment != NULL) {
        if (StringICmp (sfp->comment, "putative") != 0) {
          ValNodeAddStr (&(prp->name), 0, sfp->comment);
          sfp->comment = NULL;
        }
      }
      if (prp->processed == 3 || prp->processed == 4 || prp->processed == 5) {
        if (prp->name != NULL) {
          str = (CharPtr) prp->name->data.ptrvalue;
          if ((StringStr (str, "putative") != NULL ||
               StringStr (str, "put. ") != NULL) &&
              sfp->comment == NULL) {
            sfp->comment = StringSave ("putative");
          }
          if (! InformativeString (str)) {
            prp->name = ValNodeFreeData (prp->name);
          }
        }
      }
      if ((prp->processed == 1 || prp->processed == 2) && prp->name == NULL) {
        ValNodeCopyStr (&(prp->name), 0, "unnamed");
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringICmp (str, "RbcL") == 0 || StringICmp (str, "rubisco large subunit") == 0) {
          vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit");
          MemFree (str);
        } else if (StringICmp (str, "RbcS") == 0 || StringICmp (str, "rubisco small subunit") == 0) {
          vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit");
          MemFree (str);
        }
      }
      /*
      if (StringDoesHaveText (prp->desc)) {
        for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          if (StringCmp (prp->desc, str) == 0) {
            prp->desc = MemFree (prp->desc);
          }
        }
      }
      */
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->ext.choice == 1) {
        BSECDecodeXml ((CharPtr) rrp->ext.value.ptrvalue);
        str = (CharPtr) rrp->ext.value.ptrvalue;
        CleanVisStringAndCompress ((CharPtr PNTR) &(rrp->ext.value.ptrvalue));
        CleanDoubleQuote ((CharPtr) rrp->ext.value.ptrvalue);
        RemoveFlankingQuotes ((CharPtr PNTR) &(rrp->ext.value.ptrvalue));
        if (rrp->ext.value.ptrvalue == NULL) {
          rrp->ext.choice = 0;
        } else if (rrp->type == 4) {
          name = (CharPtr) rrp->ext.value.ptrvalue;
          len = StringLen (name);
          if (len > 5) {
            if (len > 16 && StringNICmp (name + len - 16, " ribosomal RNA .", 14) == 0) {
              name [len-2] = '\0';
              len = StringLen (name);
            }
            if (len > 14 && StringNICmp (name + len - 14, " ribosomal rRNA", 14) == 0) {
            } else if (StringNICmp (name + len - 5, " rRNA", 5) == 0) {
              str = MemNew (len + 10);
              if (str != NULL) {
                StringNCpy (str, name, len - 5);
                StringCat (str, " ribosomal RNA");
                rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                rrp->ext.value.ptrvalue = (Pointer) str;
              }
            } else if (StringNICmp (name + len - 5, "_rRNA", 5) == 0) {
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
          aa = ParseTRnaString (name, &justTrnaText, codon, FALSE);
          if (aa != 0) {
            is_fMet = (Boolean) (StringStr (name, "fMet") != NULL);
            is_iMet = (Boolean) (StringStr (name, "iMet") != NULL);
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
            trp = (tRNAPtr) MemNew (sizeof (tRNA));
            if (trp != NULL) {
              trp->aatype = 2;
              for (j = 0; j < 6; j++) {
                trp->codon [j] = 255;
              }
              if (justTrnaText) {
                for (j = 0; j < 6; j++) {
                  trp->codon [j] = codon [j];
                }
              }
              trp->aa = aa;
              rrp->ext.choice = 2;
              rrp->ext.value.ptrvalue = (Pointer) trp;
              CleanupTrna (sfp, trp);
            }
            if (is_fMet) {
              if (sfp->comment == NULL) {
                sfp->comment = StringSave ("fMet");
              } else {
                len = StringLen (sfp->comment) + StringLen ("fMet") + 5;
                str = MemNew (sizeof (Char) * len);
                StringCpy (str, sfp->comment);
                StringCat (str, "; ");
                StringCat (str, "fMet");
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = str;
              }
            }
            if (is_iMet) {
              if (sfp->comment == NULL) {
                sfp->comment = StringSave ("iMet");
              } else {
                len = StringLen (sfp->comment) + StringLen ("iMet") + 5;
                str = MemNew (sizeof (Char) * len);
                StringCpy (str, sfp->comment);
                StringCat (str, "; ");
                StringCat (str, "iMet");
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = str;
              }
            }
          }
        }
      } else if (rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        CleanupTrna (sfp, trp);
      } else if (rrp->type == 3 && (! StringHasNoText (sfp->comment))) {
        aa = ParseTRnaString (sfp->comment, &justTrnaText, codon, TRUE);
        if (aa != 0) {
          trp = (tRNAPtr) MemNew (sizeof (tRNA));
          if (trp != NULL) {
            trp->aatype = 2;
            for (j = 0; j < 6; j++) {
              trp->codon [j] = 255;
            }
            if (justTrnaText) {
              for (j = 0; j < 6; j++) {
                trp->codon [j] = codon [j];
              }
            }
            trp->aa = aa;
            rrp->ext.choice = 2;
            rrp->ext.value.ptrvalue = (Pointer) trp;
            if (justTrnaText) {
              if (StringCmp (sfp->comment, "tRNA-fMet") != 0 &&
                  StringCmp (sfp->comment, "fMet") != 0 &&
                  StringCmp (sfp->comment, "fMet tRNA") != 0 &&
                  StringCmp (sfp->comment, "fMet-tRNA") != 0) {
                sfp->comment = MemFree (sfp->comment);
              } else {
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = StringSave ("fMet");
              }
              if (StringCmp (sfp->comment, "tRNA-iMet") != 0 &&
                  StringCmp (sfp->comment, "iMet") != 0 &&
                  StringCmp (sfp->comment, "iMet tRNA") != 0 &&
                  StringCmp (sfp->comment, "iMet-tRNA") != 0) {
                sfp->comment = MemFree (sfp->comment);
              } else {
                sfp->comment = MemFree (sfp->comment);
                sfp->comment = StringSave ("iMet");
              }
            }
          }
        }
      }
      if (rrp->ext.choice == 3) {
        rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
        if (rgp != NULL) {
          str = (CharPtr) rgp->product;
          CleanVisStringAndCompress (&(rgp->product));
          CleanDoubleQuote (rgp->product);
          RemoveFlankingQuotes (&(rgp->product));
          if (StringICmp (rgp->product, "internal transcribed spacer 1 (ITS1)") == 0) {
            rgp->product = MemFree (rgp->product);
            rgp->product = StringSave ("internal transcribed spacer 1");
          } else if (StringICmp (rgp->product, "internal transcribed spacer 2 (ITS2)") == 0) {
            rgp->product = MemFree (rgp->product);
            rgp->product = StringSave ("internal transcribed spacer 2");
          } else if (StringICmp (rgp->product, "internal transcribed spacer 3 (ITS3)") == 0) {
            rgp->product = MemFree (rgp->product);
            rgp->product = StringSave ("internal transcribed spacer 3");
          }
          CleanVisStringAndCompress (&(rgp->_class));
          CleanDoubleQuote (rgp->_class);
          for (rqp = rgp->quals; rqp != NULL; rqp = rqp->next) {
            CleanVisStringAndCompress (&(rqp->qual));
            CleanDoubleQuote (rqp->qual);
            CleanVisStringAndCompress (&(rqp->val));
            CleanDoubleQuote (rqp->val);
          }
        }
      }
      if (rrp->ext.choice == 0 && sfp->comment != NULL && rrp->type == 4) {
        len = StringLen (sfp->comment);
        if (len > 15 && len < 20) {
          if (StringNICmp (sfp->comment + len - 15, "S ribosomal RNA", 15) == 0) {
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = sfp->comment;
            sfp->comment = NULL;
          }
        } else if (len > 6 && len < 20) {
          if (StringNICmp (sfp->comment + len - 6, "S rRNA", 6) == 0) {
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = sfp->comment;
            sfp->comment = NULL;
          }
        }
      }
/*
 * This section has been commented out based on a request by DeAnne Cravaritis.
 * If left in, this causes unexpected results when RNA comments are copied to 
 * the product name or vice versa.      
      if (rrp->ext.choice == 1 && rrp->ext.value.ptrvalue != NULL) {
        if (StringICmp ((CharPtr) rrp->ext.value.ptrvalue, sfp->comment) == 0) {
          sfp->comment = MemFree (sfp->comment);
        }
      }
*/      
      if (rrp->type == 4 && rrp->ext.choice == 1 ) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        len = StringLen (name);
        if (len > 5 && NotExceptedRibosomalName (name)) {
          suff = NULL;
          str = StringStr (name, " ribosomal");
          if (str != NULL) {
            suff = str + 10;
            ch = *suff;
            if (ch != '\0' && ch != ' ') {
              suff = NULL;
              str = NULL;
            }
          }
          if (str == NULL) {
            str = StringStr (name, " rRNA");
            if (str != NULL) {
              suff = str + 5;
              ch = *suff;
              if (ch != '\0' && ch != ' ') {
                suff = NULL;
                str = NULL;
              }
            }
          }
          if (suff != NULL && StringNICmp (suff, " RNA", 4) == 0) {
            suff += 4;
          }
          if (suff != NULL && StringNICmp (suff, " DNA", 4) == 0) {
            suff += 4;
          }
          if (suff != NULL && StringNICmp (suff, " ribosomal", 10) == 0) {
            suff += 10;
          }
          TrimSpacesAroundString (suff);
          if (str != NULL) {
            *str = '\0';
            len = StringLen (name);
            if (StringHasNoText (suff)) {
              suff = NULL;
            }
            if (suff != NULL) {
              len += StringLen (suff) + 2;
            }
            str = MemNew (len + 15);
            if (str != NULL) {
              StringCpy (str, name);
              StringCat (str, " ribosomal RNA");
              if (suff != NULL) {
                ch = *suff;
                if (ch != ',' && ch != ';') {
                  StringCat (str, " ");
                }
                StringCat (str, suff);
              }
              rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
              rrp->ext.value.ptrvalue = (Pointer) str;
            }
          }
        }
        name = (CharPtr) rrp->ext.value.ptrvalue;
        len = StringLen (name);
        if (len > 5) {
          ch = *name;
          while (ch != '\0' && (ch == '.' || (IS_DIGIT (ch)))) {
            name++;
            ch = *name;
          }
          /*
          if (ch == 's' && StringCmp (name, "s ribosomal RNA") == 0) {
            *name = 'S';
          }
          */
          if (ch == 's' && name [1] == ' ') {
            *name = 'S';
          }
        }
        StrStripSpaces ((CharPtr) rrp->ext.value.ptrvalue);
        name = (CharPtr) rrp->ext.value.ptrvalue;
        len = StringLen (name);
        if (len > 17) {
          if (StringNICmp (name + len - 17, "ribosomal RNA RNA", 17) == 0) {
            *(name + len - 4) = '\0';
          }
        }
        trimming_junk = TRUE;
        while (trimming_junk) {
          StrStripSpaces ((CharPtr) rrp->ext.value.ptrvalue);
          name = (CharPtr) rrp->ext.value.ptrvalue;
          ptr = StringStr (name, "ribosomal ribosomal");
          if (ptr != NULL) {
            suff = ptr + 19;
            *(ptr + 10) = '\0';
            temp = MemNew (StringLen (name) + StringLen (suff) + 2);
            TrimSpacesAroundString (suff);
            StringCpy (temp, name);
            if (suff [0] != ' ' && suff [0] != '\0') {
              StringCat (temp, " ");
            }
            StringCat (temp, suff);
            rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
            rrp->ext.value.ptrvalue = (Pointer) temp;
          } else {
            ptr = StringStr (name, "RNA RNA");
            if (ptr != NULL) {
              suff = ptr + 7;
              *(ptr + 4) = '\0';
              temp = MemNew (StringLen (name) + StringLen (suff) + 2);
              TrimSpacesAroundString (suff);
              StringCpy (temp, name);
              if (suff [0] != ' ' && suff [0] != '\0') {
                StringCat (temp, " ");
              }
              StringCat (temp, suff);
              rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
              rrp->ext.value.ptrvalue = (Pointer) temp;
            } else {
              ptr = StringStr (name, "ribosomal RNA ribosomal");
              if (ptr != NULL) {
                suff = ptr + 23;
                *(ptr + 14) = '\0';
                temp = MemNew (StringLen (name) + StringLen (suff) + 2);
                TrimSpacesAroundString (suff);
                StringCpy (temp, name);
                if (suff [0] != ' ' && suff [0] != '\0') {
                  StringCat (temp, " ");
                }
                StringCat (temp, suff);
                rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                rrp->ext.value.ptrvalue = (Pointer) temp;
              } else {
                ptr = StringStr (name, "ribosomal rRNA");
                if (ptr != NULL) {
                  suff = ptr + 14;
                  *(ptr + 10) = '\0';
                  temp = MemNew (StringLen (name) + StringLen (" RNA") + StringLen (suff) + 2);
                  TrimSpacesAroundString (suff);
                  StringCpy (temp, name);
                  StringCat (temp, " RNA");
                  if (suff [0] != ' ' && suff [0] != '\0') {
                    StringCat (temp, " ");
                  }
                  StringCat (temp, suff);
                  rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                  rrp->ext.value.ptrvalue = (Pointer) temp;
                } else {
                  ptr = StringStr (name, "RNA rRNA");
                  if (ptr != NULL) {
                    suff = ptr + 8;
                    *(ptr + 3) = '\0';
                    temp = MemNew (StringLen (name) + StringLen (suff) + 2);
                    TrimSpacesAroundString (suff);
                    StringCpy (temp, name);
                    if (suff [0] != ' ' && suff [0] != '\0') {
                      StringCat (temp, " ");
                    }
                    StringCat (temp, suff);
                    rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
                    rrp->ext.value.ptrvalue = (Pointer) temp;
                  } else {
                    trimming_junk = FALSE;
                  }
                }
              }
            }
          }
        }
        TrimSpacesAroundString ((CharPtr) rrp->ext.value.ptrvalue);
        /*
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringICmp (name, "16S rRNA. Bacterial SSU") == 0) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.value.ptrvalue = StringSave ("16S ribosomal RNA");
        } else if (StringICmp (name, "23S rRNA. Bacterial LSU") == 0) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.value.ptrvalue = StringSave ("23S ribosomal RNA");
        } else if (StringICmp (name, "5S rRNA. Bacterial TSU") == 0) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.value.ptrvalue = StringSave ("5S ribosomal RNA");
        } else if (StringICmp (name, "Large Subunit Ribosomal RNA; lsuRNA; 23S ribosomal RNA") == 0) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.value.ptrvalue = StringSave ("23S ribosomal RNA");
        } else if (StringICmp (name, "Small Subunit Ribosomal RNA; ssuRNA; 16S ribosomal RNA") == 0) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.value.ptrvalue = StringSave ("16S ribosomal RNA");
        } else if (StringICmp (name, "Small Subunit Ribosomal RNA; ssuRNA; SSU ribosomal RNA") == 0) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.value.ptrvalue = StringSave ("small subunit ribosomal RNA");
        } else if (StringICmp (name, "Large Subunit Ribosomal RNA; lsuRNA; LSU ribosomal RNA") == 0) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.value.ptrvalue = StringSave ("large subunit ribosomal RNA");
        }
        */
      }
      /*
      if (rrp->type == 2 && rrp->ext.choice == 0 && sfp->comment != NULL) {
        rrp->ext.choice = 1;
        rrp->ext.value.ptrvalue = sfp->comment;
        sfp->comment = NULL;
      }
      */
      if (rrp->type == 2 && rrp->ext.choice == 0 && sfp->comment != NULL) {
        len = StringLen (sfp->comment);
        if (len > 5) {
          if (StringNICmp (sfp->comment + len - 4, " RNA", 4) == 0 ||
              StringNICmp (sfp->comment + len - 5, " mRNA", 5) == 0) {
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = sfp->comment;
            sfp->comment = NULL;
          }
        }
      }
      if (rrp->type == 255 || rrp->type == 10) {
        name = GetRNARefProductString (rrp, NULL);
        if (StringICmp (name, "its1") == 0 || StringICmp (name, "its 1") == 0) {
          SetRNARefProductString (rrp, NULL, "internal transcribed spacer 1", ExistingTextOption_replace_old);
        } else if (StringICmp (name, "its2") == 0 || StringICmp (name, "its 2") == 0) {
          SetRNARefProductString (rrp, NULL, "internal transcribed spacer 2", ExistingTextOption_replace_old);
        } else if (StringICmp (name, "its3") == 0 || StringICmp (name, "its 3") == 0) {
          SetRNARefProductString (rrp, NULL, "internal transcribed spacer 3", ExistingTextOption_replace_old);
        }
        name = MemFree (name);
      }
      if ((rrp->type == 255 || rrp->type == 10) && rrp->ext.choice == 0 && sfp->comment != NULL) {
        if (StringICmp (sfp->comment, "internal transcribed spacer 1") == 0 ||
            StringICmp (sfp->comment, "internal transcribed spacer 2") == 0 ||
            StringICmp (sfp->comment, "internal transcribed spacer 3") == 0) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = sfp->comment;
          sfp->comment = NULL;
        } else if (StringICmp (sfp->comment, "internal transcribed spacer 1 (ITS1)") == 0 ||
            StringICmp (sfp->comment, "internal transcribed spacer 2 (ITS2)") == 0 ||
            StringICmp (sfp->comment, "internal transcribed spacer 3 (ITS3)") == 0) {
          ptr = StringStr (sfp->comment, " (");
          if (ptr != NULL) {
            *ptr = '\0';
          }
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = sfp->comment;
          sfp->comment = NULL;
        } else if (StringICmp (sfp->comment, "ITS1") == 0 || StringICmp (sfp->comment, "ITS 1") == 0) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("internal transcribed spacer 1");
          sfp->comment = MemFree (sfp->comment);
        } else if (StringICmp (sfp->comment, "ITS2") == 0 || StringICmp (sfp->comment, "ITS 2") == 0) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("internal transcribed spacer 2");
          sfp->comment = MemFree (sfp->comment);
        } else if (StringICmp (sfp->comment, "ITS3") == 0 || StringICmp (sfp->comment, "ITS 3") == 0) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("internal transcribed spacer 3");
          sfp->comment = MemFree (sfp->comment);
        }
      }
      break;
    case SEQFEAT_PUB :
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      CleanDoubleQuote (pdp->comment);
      NormalizePubdesc (pdp, stripSerial, TRUE, publist);
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
      CleanVisStringAndCompress ((CharPtr PNTR) &(sfp->data.value.ptrvalue));
      CleanDoubleQuote ((CharPtr) sfp->data.value.ptrvalue);
      if (sfp->data.value.ptrvalue == NULL) {
        sfp->data.choice = SEQFEAT_COMMENT;
      } else {
        if (sfp->ext != NULL) {
          uop = FindUopByTag (sfp->ext, "cddScoreData");
          if (uop != NULL) {
            for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
              if (ufp->choice != 1) continue;
              oip = ufp->label;
              if (oip == NULL) continue;
              if (StringICmp (oip->str, "definition") == 0) {
                CleanVisStringAndCompress ((CharPtr PNTR) &(ufp->data.ptrvalue));
                CleanDoubleQuote ((CharPtr) ufp->data.ptrvalue);
              }
            }
          }
        }
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
      VisitAllUserObjectsInUop ((UserObjectPtr) sfp->data.value.ptrvalue, NULL, CleanUserObject);
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
      if (biop != NULL) {
        if (biop->genome == GENOME_virion) {
          biop->genome = GENOME_unknown;
        }
        orp = biop->org;
        if (orp != NULL) {
          CleanVisStringListAndCompress (&(orp->mod));
          OrpModToSubSource (&(orp->mod), &(biop->subtype));
          onp = orp->orgname;
          if (onp != NULL) {
            CleanupOrgModOther (biop, onp);
          }
        }
        biop->subtype = SortSubSourceList (biop->subtype);
        CleanSubSourceList (&(biop->subtype), biop->genome);
        CleanupSubSourceOther (biop, onp);
        biop->subtype = SortSubSourceList (biop->subtype);
        if (modernizeFeats) {
          ModernizePCRPrimers (biop);
        }
        CleanupPCRReactionSet (&(biop->pcr_primers));
        if (biop->genome == GENOME_unknown || biop->genome == GENOME_genomic) {
          for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
            if (ssp->subtype == SUBSRC_plasmid_name) {
              biop->genome = GENOME_plasmid;
            }
          }
        }
      }
      break;
    default :
      break;
  }
  if (orp != NULL) {
    CleanVisStringAndCompress (&(orp->taxname));
    CleanVisStringAndCompress (&(orp->common));
    CleanVisStringList (&(orp->mod));
    CleanVisStringList (&(orp->syn));
    FixOldDbxrefs (orp->db, isEmblOrDdbj);
    FixNumericDbxrefs (orp->db);
    orp->db = ValNodeSort (orp->db, SortDbxref);
    orp->syn = ValNodeSort (orp->syn, SortVnpByString);
    orp->syn = UniqueValNode (orp->syn);
    CleanupDuplicateDbxrefs (&(orp->db));
    CleanupObsoleteDbxrefs (&(orp->db));
    CleanupGoDbxrefs (orp->db);
    onp = orp->orgname;
    while (onp != NULL) {
      CleanVisString (&(onp->attrib));
      CleanVisString (&(onp->lineage));
      CleanVisString (&(onp->div));
      OrpModToOrgMod (&(orp->mod), &(onp->mod));
      onp->mod = SortOrgModList (onp->mod);
      CleanOrgModListEx (&(onp->mod), orp->common);
      onp->mod = SortOrgModList (onp->mod);
      onp = onp->next;
    }
  }
}

static ValNodePtr SplitStringsAtSemicolon (ValNodePtr PNTR head)

{
  ValNodePtr  curr, vnp;
  CharPtr     ptr, str;

  if (head == NULL || *head == NULL) return NULL;

  curr = *head;
  while (curr != NULL) {
    str = (CharPtr) curr->data.ptrvalue;
    ptr = StringChr (str, ';');
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
      vnp = ValNodeCopyStr (NULL, 0, ptr);
      if (vnp != NULL) {
        vnp->next = curr->next;
        curr->next = vnp;
      }
    }
    curr = curr->next;
  }

  return *head;
}


static void CleanupDescriptorStrings (
  ValNodePtr sdp,
  Boolean stripSerial,
  Boolean modernizeFeats,
  ValNodePtr PNTR publist,
  Boolean isEmblOrDdbj
)

{
  BioSourcePtr  biop;
  EMBLBlockPtr  ebp;
  GBBlockPtr    gbp;
  OrgNamePtr    onp = NULL;
  OrgRefPtr     orp;
  PubdescPtr    pdp;
  PirBlockPtr   pir;
  PrfBlockPtr   prf;
  SPBlockPtr    sp;
  SubSourcePtr  ssp;
  CharPtr       str;
  ValNodePtr    vnp;

  if (sdp == NULL) return;
  switch (sdp->choice) {
    case Seq_descr_mol_type :
    case Seq_descr_method :
      return;
    default :
      break;
  }
  if (sdp->data.ptrvalue == NULL) return;

  biop = NULL;
  orp = NULL;
  switch (sdp->choice) {
    case Seq_descr_mol_type :
      break;
    case Seq_descr_modif :
      break;
    case Seq_descr_method :
      break;
    case Seq_descr_name :
      CleanVisString ((CharPtr PNTR) &sdp->data.ptrvalue);
      if (sdp->data.ptrvalue == NULL) {
        sdp->data.ptrvalue = StringSave ("");
      }
      break;
    case Seq_descr_title :
      BSECDecodeXml ((CharPtr) sdp->data.ptrvalue);
      str = (CharPtr) sdp->data.ptrvalue;
      CleanVisStringAndCompress ((CharPtr PNTR) &sdp->data.ptrvalue);
      if (sdp->data.ptrvalue == NULL) {
        sdp->data.ptrvalue = StringSave ("");
      }
      break;
    case Seq_descr_org :
      orp = (OrgRefPtr) sdp->data.ptrvalue;
      break;
    case Seq_descr_comment :
      BSECDecodeXml ((CharPtr) sdp->data.ptrvalue);
      CleanVisStringJunk ((CharPtr PNTR) &sdp->data.ptrvalue);
      RemoveSpacesBetweenTildes ((CharPtr) sdp->data.ptrvalue);
      if (sdp->data.ptrvalue == NULL) {
        sdp->data.ptrvalue = StringSave ("");
      }
      break;
    case Seq_descr_num :
      break;
    case Seq_descr_maploc :
      break;
    case Seq_descr_pir :
      pir = (PirBlockPtr) sdp->data.ptrvalue;
      SplitStringsAtSemicolon (&(pir->keywords));
      break;
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      SplitStringsAtSemicolon (&(gbp->keywords));
      for (vnp = gbp->keywords; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringICmp (str, "TPA:reassembly") == 0) {
          vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
          vnp->data.ptrvalue = StringSave ("TPA:assembly");
        } else if (StringICmp (str, "TPA_reassembly") == 0) {
          vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
          vnp->data.ptrvalue = StringSave ("TPA:assembly");
        } else if (StringICmp (str, "TPA_assembly") == 0) {
          vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
          vnp->data.ptrvalue = StringSave ("TPA:assembly");
        }
      }
      CleanVisStringList (&(gbp->extra_accessions));
      gbp->extra_accessions = ValNodeSort (gbp->extra_accessions, SortVnpByString);
      gbp->extra_accessions = UniqueValNode (gbp->extra_accessions);
      if (isEmblOrDdbj) {
        CleanVisStringListCaseSensitive (&(gbp->keywords));
      } else {
        CleanVisStringList (&(gbp->keywords));
      }
      CleanVisStringJunk (&(gbp->source));
      if (StringCmp (gbp->source, ".") == 0) {
        gbp->source = MemFree (gbp->source);
      }
      CleanVisStringJunk (&(gbp->origin));
      if (StringCmp (gbp->origin, ".") == 0) {
        gbp->origin = MemFree (gbp->origin);
      }
      CleanVisString (&(gbp->date));
      CleanVisString (&(gbp->div));
      CleanVisString (&(gbp->taxonomy));
      break;
    case Seq_descr_pub :
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      CleanDoubleQuote (pdp->comment);
      NormalizePubdesc (pdp, stripSerial, TRUE, publist);
      break;
    case Seq_descr_region :
      CleanVisString ((CharPtr PNTR) &sdp->data.ptrvalue);
      if (sdp->data.ptrvalue == NULL) {
        sdp->data.ptrvalue = StringSave ("");
      }
      break;
    case Seq_descr_user :
      VisitAllUserObjectsInUop ((UserObjectPtr) sdp->data.ptrvalue, NULL, CleanUserObject);
      break;
    case Seq_descr_sp :
      sp = (SPBlockPtr) sdp->data.ptrvalue;
      SplitStringsAtSemicolon (&(sp->keywords));
      break;
    case Seq_descr_dbxref :
      break;
    case Seq_descr_embl :
      ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
      CleanVisStringList (&(ebp->extra_acc));
      ebp->extra_acc = ValNodeSort (ebp->extra_acc, SortVnpByString);
      SplitStringsAtSemicolon (&(ebp->keywords));
      CleanVisStringListCaseSensitive (&(ebp->keywords));
      break;
    case Seq_descr_create_date :
      break;
    case Seq_descr_update_date :
      break;
    case Seq_descr_prf :
      prf = (PrfBlockPtr) sdp->data.ptrvalue;
      SplitStringsAtSemicolon (&(prf->keywords));
      break;
    case Seq_descr_pdb :
      break;
    case Seq_descr_het :
      break;
    case Seq_descr_source :
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        if (biop->genome == GENOME_virion) {
          biop->genome = GENOME_unknown;
        }
        orp = biop->org;
        if (orp != NULL) {
          CleanVisStringList (&(orp->mod));
          OrpModToSubSource (&(orp->mod), &(biop->subtype));
          onp = orp->orgname;
          if (onp != NULL) {
            CleanupOrgModOther (biop, onp);
          }
        }
        biop->subtype = SortSubSourceList (biop->subtype);
        CleanSubSourceList (&(biop->subtype), biop->genome);
        CleanupSubSourceOther (biop, onp);
        biop->subtype = SortSubSourceList (biop->subtype);
        if (modernizeFeats) {
          ModernizePCRPrimers (biop);
        }
        CleanupPCRReactionSet (&(biop->pcr_primers));
        if (biop->genome == GENOME_unknown || biop->genome == GENOME_genomic) {
          for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
            if (ssp->subtype == SUBSRC_plasmid_name) {
              biop->genome = GENOME_plasmid;
            }
          }
        }
      }
      break;
    case Seq_descr_molinfo :
      break;
    default :
      break;
  }
  if (orp != NULL) {
    CleanVisStringAndCompress (&(orp->taxname));
    CleanVisStringAndCompress (&(orp->common));
    CleanVisStringList (&(orp->mod));
    CleanVisStringList (&(orp->syn));
    FixOldDbxrefs (orp->db, isEmblOrDdbj);
    FixNumericDbxrefs (orp->db);
    orp->db = ValNodeSort (orp->db, SortDbxref);
    orp->syn = ValNodeSort (orp->syn, SortVnpByString);
    orp->syn = UniqueValNode (orp->syn);
    CleanupDuplicateDbxrefs (&(orp->db));
    CleanupObsoleteDbxrefs (&(orp->db));
    CleanupGoDbxrefs (orp->db);
    onp = orp->orgname;
    while (onp != NULL) {
      CleanVisString (&(onp->attrib));
      CleanVisString (&(onp->lineage));
      CleanVisString (&(onp->div));
      OrpModToOrgMod (&(orp->mod), &(onp->mod));
      onp->mod = SortOrgModList (onp->mod);
      CleanOrgModListEx (&(onp->mod), orp->common);
      onp->mod = SortOrgModList (onp->mod);
      onp = onp->next;
    }
  }
}

static Int2 CheckForQual (GBQualPtr gbqual, CharPtr string_q, CharPtr string_v)

{
  GBQualPtr curq;

  for (curq = gbqual; curq; curq = curq->next) {
    if (StringCmp (string_q, curq->qual) == 0) {
      if (curq->val == NULL) {
        curq->val = StringSave (string_v);
        return 1;
      } 
      if (StringCmp (string_v, curq->val) == 0) return 1;
    }
  }
  return 0;
}
static GBQualPtr AddGBQual (GBQualPtr gbqual, CharPtr qual, CharPtr val)

{
  GBQualPtr curq;

  if (StringCmp (qual, "translation") == 0) {
    if (val == NULL)  return gbqual;
    if (*val == '\0') return gbqual;
  }
  if (gbqual) {
    if (CheckForQual (gbqual, qual, val) == 1) return gbqual;
    for (curq = gbqual; curq->next != NULL; curq = curq->next) continue;
    curq->next = GBQualNew ();
    curq = curq->next;
    if (val)
      curq->val = StringSave (val);
    curq->qual = StringSave (qual);
  } else {
    gbqual = GBQualNew ();
    gbqual->next = NULL;
    if (val)
      gbqual->val = StringSave (val);
    gbqual->qual = StringSave (qual);
  }
  return gbqual;
}

static void AddReplaceQual (SeqFeatPtr sfp, CharPtr p)

{
  CharPtr s, val;

  val = StringChr (p, '\"');
  if (val == NULL) return;
  val++;
  s = p + StringLen (p) - 1;
  if (*s != ')') return;
  for (s--; s > val && *s != '\"'; s--) continue;
  if (*s != '\"') return;
  *s = '\0';
  sfp->qual = (GBQualPtr) AddGBQual (sfp->qual, "replace", val);
  *s = '\"';
}

//LCOV_EXCL_START
NLM_EXTERN Boolean SerialNumberInString (CharPtr str)

{
  Char     ch;
  Boolean  hasdigits;
  CharPtr  ptr;
  Boolean  suspicious = FALSE;

  if (str == NULL || StringHasNoText (str)) return FALSE;
  ptr = StringChr (str, '[');

  /* bail if first digit after bracket is 0 */
  if (ptr != NULL && ptr [1] == '0') return FALSE;

  while ((! suspicious) && ptr != NULL) {
    hasdigits = FALSE;
    ptr++;
    ch = *ptr;
    while (IS_DIGIT (ch)) {
      hasdigits = TRUE;
      ptr++;
      ch = *ptr;
    }
    if (ch == ']' && hasdigits) {
      suspicious = TRUE;
    }
    if (! suspicious) {
      ptr = StringChr (ptr, '[');
    }
  }
  return suspicious;
}
//LCOV_EXCL_STOP

/* now only strips serials for local, general, refseq, and 2+6 genbank ids */
static void CheckForSwissProtID (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  SeqIdPtr      sip;
  BoolPtr       stripSerial;
  TextSeqIdPtr  tsip;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    stripSerial = (BoolPtr) mydata;
    if (stripSerial == NULL) return;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_GIBBSQ :
        case SEQID_GIBBMT :
          *stripSerial = FALSE;
          break;
        case SEQID_EMBL :
        case SEQID_PIR :
        case SEQID_SWISSPROT :
        case SEQID_PATENT :
        case SEQID_DDBJ :
        case SEQID_PRF :
        case SEQID_PDB :
        case SEQID_TPE:
        case SEQID_TPD:
        case SEQID_GPIPE:
          *stripSerial = FALSE;
          break;
        case SEQID_GENBANK :
        case SEQID_TPG:
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL) {
            if (StringLen (tsip->accession) == 6) {
              *stripSerial = FALSE;
            }
          }
          break;
        case SEQID_NOT_SET :
        case SEQID_LOCAL :
        case SEQID_OTHER :
        case SEQID_GENERAL :
          break;
        default :
          break;
      }
    }
  }
}

static void CheckForEmblDdbjID (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  BoolPtr    isEmblOrDdbj;
  SeqIdPtr   sip;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    isEmblOrDdbj = (BoolPtr) mydata;
    if (isEmblOrDdbj == NULL) return;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_EMBL :
        case SEQID_DDBJ :
          *isEmblOrDdbj = TRUE;
          break;
        default :
          break;
      }
    }
  }
}

static void CheckForJournalScanID (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  BoolPtr    isJScan;
  SeqIdPtr   sip;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    isJScan = (BoolPtr) mydata;
    if (isJScan == NULL) return;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_GIBBSQ :
        case SEQID_GIBBMT :
        case SEQID_GIIM :
          *isJScan = TRUE;
          break;
        default :
          break;
      }
    }
  }
}

NLM_EXTERN Boolean FixWrongFuzzOnPlusStrand (SeqLocPtr location)

{
  SeqLocPtr   firstSlp;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  Boolean     res = FALSE;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  SeqPntPtr   spp;

  if (location == NULL) return FALSE;

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

  if (firstSlp != NULL && firstSlp->choice == SEQLOC_INT && firstSlp->data.ptrvalue != NULL) {
    sip = (SeqIntPtr) firstSlp->data.ptrvalue;
    if (sip != NULL && (sip->strand == Seq_strand_plus || sip->strand == Seq_strand_unknown)) {
      if (sip->if_to != NULL && sip->if_from == NULL) {
        sip->if_from = IntFuzzFree (sip->if_from);
        ifp = IntFuzzNew ();
        if (ifp != NULL) {
          ifp->choice = 4;
          sip->if_from = ifp;
          ifp->a = 2;
          res = TRUE;
        }
      }
    }
  }

  if (lastSlp != NULL && lastSlp->choice == SEQLOC_INT && lastSlp->data.ptrvalue != NULL) {
    sip = (SeqIntPtr) lastSlp->data.ptrvalue;
    if (sip != NULL && (sip->strand == Seq_strand_plus || sip->strand == Seq_strand_unknown)) {
      if (sip->if_to == NULL && sip->if_from != NULL) {
        sip->if_to = IntFuzzFree (sip->if_to);
        ifp = IntFuzzNew ();
        if (ifp != NULL) {
          ifp->choice = 4;
          sip->if_to = ifp;
          ifp->a = 1;
          res = TRUE;
        }
      }
    }
  }

  return res;
}

NLM_EXTERN Boolean FixWrongFuzzOnMinusStrand (SeqLocPtr location)

{
  SeqLocPtr   firstSlp;
  IntFuzzPtr  ifp;
  SeqLocPtr   lastSlp;
  Boolean     res = FALSE;
  SeqIntPtr   sip;
  SeqLocPtr   slp;
  SeqPntPtr   spp;

  if (location == NULL) return FALSE;

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

  if (firstSlp != NULL && firstSlp->choice == SEQLOC_INT && firstSlp->data.ptrvalue != NULL) {
    sip = (SeqIntPtr) firstSlp->data.ptrvalue;
    if (sip != NULL && (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev)) {
      if (sip->if_to == NULL && sip->if_from != NULL) {
        sip->if_from = IntFuzzFree (sip->if_from);
        ifp = IntFuzzNew ();
        if (ifp != NULL) {
          ifp->choice = 4;
          sip->if_to = ifp;
          ifp->a = 1;
          res = TRUE;
        }
      }
    }
  }

  if (lastSlp != NULL && lastSlp->choice == SEQLOC_INT && lastSlp->data.ptrvalue != NULL) {
    sip = (SeqIntPtr) lastSlp->data.ptrvalue;
    if (sip != NULL && (sip->strand == Seq_strand_minus || sip->strand == Seq_strand_both_rev)) {
      if (sip->if_to != NULL && sip->if_from == NULL) {
        sip->if_to = IntFuzzFree (sip->if_to);
        ifp = IntFuzzNew ();
        if (ifp != NULL) {
          ifp->choice = 4;
          sip->if_from = ifp;
          ifp->a = 2;
          res = TRUE;
        }
      }
    }
  }

  return res;
}

NLM_EXTERN void CleanUpSeqLoc (SeqLocPtr slp)

{
  BioseqPtr  bsp;
  SeqLocPtr  curr;
  SeqLocPtr  head;
  SeqLocPtr  last;
  SeqLocPtr  loc;
  SeqLocPtr  next;
  SeqIdPtr   sip;
  SeqIntPtr  sintp;
  SeqPntPtr  spp;
  Int4       swp;
  SeqLocPtr  tail;

  if (slp == NULL) return;

  if (slp->choice == SEQLOC_WHOLE) {
    sip = (SeqIdPtr) slp->data.ptrvalue;
    if (sip != NULL) {
      bsp = BioseqFind (sip);
      if (bsp != NULL) {
        sintp = SeqIntNew ();
        if (sintp != NULL) {
          sintp->from = 0;
          sintp->to = bsp->length - 1;
          sintp->id = sip; /* reuse existing slp->data.ptrvalue, no need to free */
          slp->choice = SEQLOC_INT;
          slp->data.ptrvalue = (Pointer) sintp;
        }
      }
    }
  }

  /* from < to for all intervals */
  loc = SeqLocFindNext (slp, NULL);
  while (loc != NULL) {
    if (loc->choice == SEQLOC_INT) {
      sintp = (SeqIntPtr) loc->data.ptrvalue;
      if (sintp != NULL) {
        if (sintp->from > sintp->to) {
          swp = sintp->from;
          sintp->from = sintp->to;
          sintp->to = swp;
        }
        if (sintp->strand == Seq_strand_both) {
          sintp->strand = Seq_strand_plus;
        } else if (sintp->strand == Seq_strand_both_rev) {
          sintp->strand = Seq_strand_minus;
        }
      }
    } else if (loc->choice == SEQLOC_PNT) {
      spp = (SeqPntPtr) loc->data.ptrvalue;
      if (spp != NULL) {
        if (spp->strand == Seq_strand_both) {
          spp->strand = Seq_strand_plus;
        } else if (spp->strand == Seq_strand_both_rev) {
          spp->strand = Seq_strand_minus;
        }
      }
    }
    loc = SeqLocFindNext (slp, loc);
  }

  if (slp->choice == SEQLOC_PACKED_INT) {
    loc = (SeqLocPtr) slp->data.ptrvalue;
    if (loc == NULL || loc->next != NULL) return;
    /* here seqloc_packed_int points to a single location element, so no need for seqloc_packed_int parent */
    slp->choice = loc->choice;
    slp->data.ptrvalue = (Pointer) loc->data.ptrvalue;
    MemFree (loc);
    return;
  }

  if (slp->choice != SEQLOC_MIX) return;
  loc = (SeqLocPtr) slp->data.ptrvalue;
  if (loc == NULL) return;

  if (loc->next != NULL) {
    /* check for null NULL at beginning */
    if (loc->choice == SEQLOC_NULL) {
      slp->data.ptrvalue = (Pointer) loc->next;
      loc->next = NULL;
      ValNodeFree (loc);
    }
    /* check for null NULL at end */
    loc = (SeqLocPtr) slp->data.ptrvalue;
    last = NULL;
    while (loc->next != NULL) {
      last = loc;
      loc = loc->next;
    }
    if (loc->choice == SEQLOC_NULL && last != NULL) {
      last->next = NULL;
      ValNodeFree (loc);
    }
  }

  loc = (SeqLocPtr) slp->data.ptrvalue;
  if (loc == NULL) return;

  if (loc->next == NULL) {
    /* here seqloc_mix points to a single location element, so no need for seqloc_mix parent */
    slp->choice = loc->choice;
    slp->data.ptrvalue = (Pointer) loc->data.ptrvalue;
    MemFree (loc);
    return;
  }

  /* check for nested seqloc_mix, remove nesting */
  curr = loc;
  last = NULL;
  while (curr != NULL) {
    next = curr->next;
    if (curr->choice == SEQLOC_MIX) {
      head = (SeqLocPtr) curr->data.ptrvalue;
      if (head != NULL) {
        tail = head;
        while (tail->next != NULL) {
          tail = tail->next;
        }
        if (last != NULL) {
          last->next = head;
        }
        tail->next = curr->next;
        curr->next = NULL;
        curr = MemFree (curr);
      }
    } else {
      last = curr;
    }
    curr = next;
  }

  NormalizeNullsBetween (slp);

  /*
  FixWrongFuzzOnPlusStrand (slp);
  FixWrongFuzzOnMinusStrand (slp);
  */
}

typedef struct cbloc {
  CodeBreakPtr  cbp;
  Int4          pos;
} CbLoc, PNTR CbLocPtr;

static int LIBCALLBACK SortByCodeBreakLoc (VoidPtr ptr1, VoidPtr ptr2)

{
  CbLocPtr  clp1;
  CbLocPtr  clp2;

  clp1 = (CbLocPtr) ptr1;
  clp2 = (CbLocPtr) ptr2;
  if (clp1 == NULL || clp2 == NULL) return 0;
  if (clp1->pos < clp2->pos) {
    return -1;
  } else if (clp1->pos > clp2->pos) {
    return 1;
  }
  return 0;
}

static CodeBreakPtr SortCodeBreaks (SeqFeatPtr sfp, CodeBreakPtr list)

{
  BioseqPtr     bsp;
  CodeBreakPtr  cbp;
  CbLocPtr      head;
  size_t        count, i;
  Boolean       out_of_order = FALSE;
  Int4          pos;
  SeqLocPtr     slp;

  if (sfp == NULL || list == NULL) return list;
  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return list;

  for (cbp = list, count = 0; cbp != NULL; cbp = cbp->next, count++) continue;
  if (count < 2) return list;

  head = (CbLocPtr) MemNew (sizeof (CbLoc) * (count + 1));
  if (head == NULL) return list;

  for (cbp = list, i = 0; cbp != NULL && i < count; i++) {
    head [i].cbp = cbp;
    slp = dnaLoc_to_aaLoc (sfp, cbp->loc, TRUE, NULL, TRUE);
    head [i].pos = GetOffsetInBioseq (slp, bsp, SEQLOC_START) + 1;
    SeqLocFree (slp);
    cbp = cbp->next;
  }

  pos = head [0].pos;
  for (i = 1; i < count; i++) {
    if (head [i].pos < pos) {
      out_of_order = TRUE;
    }
    pos = head [i].pos;
  }

  if (out_of_order) {
    StableMergeSort (head, count, sizeof (CbLoc), SortByCodeBreakLoc);

    for (i = 0; i < count; i++) {
      cbp = head [i].cbp;
      cbp->next = head [i + 1].cbp;
    }

    list = head [0].cbp;
  }

  MemFree (head);

  return list;
}

static void CleanupDuplicatedCodeBreaks (CodeBreakPtr PNTR prevcbp)

{
  CodeBreakPtr  cbp;
  CodeBreakPtr  last = NULL;
  CodeBreakPtr  next;
  Boolean       unlink;

  if (prevcbp == NULL) return;
  cbp = *prevcbp;
  while (cbp != NULL) {
    next = cbp->next;
    unlink = FALSE;
    if (last != NULL) {
      if (SeqLocCompare (cbp->loc, last->loc) == SLC_A_EQ_B &&
          cbp->aa.choice == last->aa.choice &&
          cbp->aa.value.intvalue == last->aa.value.intvalue) {
        unlink = TRUE;
      }
    } else {
      last = cbp;
    }
    if (unlink) {
      *prevcbp = cbp->next;
      cbp->next = NULL;
      CodeBreakFree (cbp);
    } else {
      last = cbp;
      prevcbp = (CodeBreakPtr PNTR) &(cbp->next);
    }
    cbp = next;
  }
}

//LCOV_EXCL_START
CharPtr ncrnaClassList[] = {
"antisense_RNA",
"autocatalytically_spliced_intron",
"hammerhead_ribozyme",
"ribozyme",
"RNase_P_RNA",
"RNase_MRP_RNA",
"telomerase_RNA",
"guide_RNA",
"rasiRNA",
"scRNA",
"siRNA",
"miRNA",
"piRNA",
"snoRNA",
"snRNA",
"SRP_RNA",
"vault_RNA",
"Y_RNA",
"lncRNA",
"other",
NULL};

Int4 NcrnaOTHER = sizeof (ncrnaClassList) / sizeof (CharPtr) - 1;

extern Boolean IsStringInNcRNAClassList (CharPtr str)
{
  CharPtr PNTR p;

  if (StringHasNoText (str)) return FALSE;
  for (p = ncrnaClassList; *p != NULL; p++)
  {
    if (StringICmp (str, *p) == 0)
    {
      return TRUE;
    }
  }
  return FALSE;
}


CharPtr regulatoryClassList[] = {
"attenuator",
"CAAT_signal",
"DNase_I_hypersensitive_site",
"enhancer_blocking_element",
"enhancer",
"GC_signal",
"imprinting_control_region",
"insulator",
"locus_control_region",
"matrix_attachment_region",
"minus_10_signal",
"minus_35_signal",
"polyA_signal_sequence",
"promoter",
"recoding_stimulatory_region",
"replication_regulatory_region",
"response_element",
"ribosome_binding_site",
"riboswitch",
"silencer",
"TATA_box",
"terminator",
"transcriptional_cis_regulatory_region",
"other",
NULL};

Int4 RegulatoryOTHER = sizeof (regulatoryClassList) / sizeof (CharPtr) - 1;

extern Boolean IsStringInRegulatoryClassList (CharPtr str)

{
  CharPtr PNTR p;

  if (StringHasNoText (str)) return FALSE;
  for (p = regulatoryClassList; *p != NULL; p++)
  {
    if (StringICmp (str, *p) == 0)
    {
      return TRUE;
    }
  }
  return FALSE;
}

CharPtr recombinationClassList[] = {
"chromosome_breakpoint",
"meiotic_recombination",
"mitotic_recombination",
"non_allelic_homologous_recombination",
"other",
NULL};

extern Boolean IsStringInRecombinationClassList (CharPtr str)

{
  CharPtr PNTR p;

  if (StringHasNoText (str)) return FALSE;
  for (p = recombinationClassList; *p != NULL; p++)
  {
    if (StringICmp (str, *p) == 0)
    {
      return TRUE;
    }
  }
  return FALSE;
}
//LCOV_EXCL_STOP

static void AddNonCopiedQual (SeqFeatPtr sfp, CharPtr qual, CharPtr class_val)
{
  GBQualPtr gbq;

  if (sfp == NULL || StringHasNoText (qual) || StringHasNoText (class_val)) 
  {
    return;
  }
  gbq = sfp->qual;
  while (gbq != NULL 
          && (StringCmp (gbq->qual, qual) != 0
              || StringCmp (gbq->val, class_val) != 0))
  {
    gbq = gbq->next;
  }
  if (gbq == NULL)
  {
    gbq = GBQualNew ();
    gbq->qual = StringSave (qual);
    gbq->val = StringSave (class_val);
    gbq->next = sfp->qual;
    sfp->qual = gbq;
  }

}


static CharPtr GetMiRNAProduct (CharPtr str)
{
  Int4    len;
  CharPtr product = NULL;

  if (StringHasNoText (str)) return NULL;
  if (StringNCmp (str, "miRNA ", 6) == 0)
  {
    product = StringSave (str + 6);
  }
  else if (StringNCmp (str, "microRNA ", 9) == 0)
  {
    product = StringSave (str + 9);
  }
  else
  {
    len = StringLen (str);
    if (len > 6 && StringCmp (str + len - 6, " miRNA") == 0
        && (len < 15 || StringCmp (str + len - 15, "precursor miRNA") != 0))
    {
      product = (CharPtr) MemNew (sizeof (Char) * (len - 5));
      StringNCpy (product, str, len - 6);
      product[len - 6] = 0;
    }
    else if (len > 9 && StringCmp (str + len - 9, " microRNA") == 0
             && (len < 18 || StringCmp (str + len - 18, "precursor microRNA") != 0))
    {
      product = (CharPtr) MemNew (sizeof (Char) * (len - 8));
      StringNCpy (product, str, len - 9);
      product[len - 9] = 0;
    }
  }
  return product;
}


static Boolean ConvertToNcRNA (SeqFeatPtr sfp)
{
  GBQualPtr gbq;
  RnaRefPtr rrp;
  Boolean was_converted = FALSE;
  CharPtr miRNAproduct = NULL;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL)
  {
    return FALSE;
  }
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp->type == 5 || rrp->type == 6 || rrp->type == 7)
  {
    if (rrp->type == 5)
    {
      AddNonCopiedQual (sfp, "ncRNA_class", "snRNA");
    }
    else if (rrp->type == 6)
    {
      AddNonCopiedQual (sfp, "ncRNA_class", "scRNA");
    }
    else if (rrp->type == 7)
    {
      AddNonCopiedQual (sfp, "ncRNA_class", "snoRNA");
    }
    if (rrp->ext.choice == 1)
    {
      AddNonCopiedQual (sfp, "product", rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    }
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave ("ncRNA");
    rrp->type = 255;
    was_converted = TRUE;
  }
  else if (rrp->type == 255 && rrp->ext.choice == 1)
  {
    if (IsStringInNcRNAClassList (rrp->ext.value.ptrvalue)) 
    {
      AddNonCopiedQual (sfp, "ncRNA_class", rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave ("ncRNA");
      was_converted = TRUE;
    }
    else if ((miRNAproduct = GetMiRNAProduct (rrp->ext.value.ptrvalue)) != NULL)
    {
      AddNonCopiedQual (sfp, "ncRNA_class", "miRNA");
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave ("ncRNA");
      AddNonCopiedQual (sfp, "product", miRNAproduct);
      miRNAproduct = MemFree (miRNAproduct);
      was_converted = TRUE;
    }
    else if (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
             && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0
             && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") != 0)
    {
      AddNonCopiedQual (sfp, "product", rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave ("misc_RNA");
      was_converted = TRUE;
    }
  }
  if (rrp->type == 255 && rrp->ext.choice == 0) {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave ("misc_RNA");
  }
  if (rrp->type == 255 && rrp->ext.choice == 1 &&
      StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0) {
    for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
      if (StringCmp (gbq->qual, "ncRNA_class") == 0) {
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        rrp->ext.value.ptrvalue = StringSave ("ncRNA");
        was_converted = TRUE;
      } else if (StringCmp (gbq->qual, "tag_peptide") == 0) {
        rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
        rrp->ext.value.ptrvalue = StringSave ("tmRNA");
        was_converted = TRUE;
      }
    }
  }
  return was_converted;
}

static void ModernizeFeatureStrings (SeqFeatPtr sfp, Boolean isEmblOrDdbj)

{
  CharPtr      desc;
  GBQualPtr    gbq;
  CharPtr      name;
  ProtRefPtr   prp;
  RnaRefPtr    rrp;
  CharPtr      str;
  ValNodePtr   vnp;

  if (sfp == NULL) return;

  /* skip feature types that do not use data.value.ptrvalue */
  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
    case SEQFEAT_BOND:
    case SEQFEAT_SITE:
    case SEQFEAT_PSEC_STR:
      return;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      desc = prp->desc;
      if (! isEmblOrDdbj) {
        CleanVisStringList (&(prp->name));
        break;
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringICmp (str, "RbcL") == 0 || StringICmp (str, "rubisco large subunit") == 0) {
          vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit");
          str = MemFree (str);
          if (StringICmp (desc, "RbcL") == 0 || StringICmp (desc, "rubisco large subunit") == 0) {
            prp->desc = MemFree (prp->desc);
          }
        } else if (StringICmp (str, "RbcS") == 0 || StringICmp (str, "rubisco small subunit") == 0) {
          vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit");
          str = MemFree (str);
          if (StringICmp (desc, "RbcS") == 0 || StringICmp (desc, "rubisco small subunit") == 0) {
            prp->desc = MemFree (prp->desc);
          }
        /*
        } else if (StringCmp (desc, str) == 0) {
          prp->desc = MemFree (prp->desc);
        */
        }
        if (StringStr (str, "ribulose") != NULL &&
            StringStr (str, "bisphosphate") != NULL &&
            StringStr (str, "methyltransferase") == NULL &&
            StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit") != 0 &&
            StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit") != 0) {
          if (StringICmp (str, "ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose 1,5-bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose-bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase, large subunit") == 0 ||
              StringICmp (str, "large subunit of ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose bisphosphate carboxylase large chain") == 0 ||
              StringICmp (str, "ribulose 1,5-bisphosphate carboxylase-oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose bisphosphate carboxylase oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose 1,5 bisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase, large subunit") == 0 ||
              StringICmp (str, "large subunit of ribulose-1,5-bisphosphate carboxylase/oxgenase") == 0 ||
              StringICmp (str, "ribulose bisphosphate carboxylase/oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose-1,5-bisphosphate carboxylase oxygenase, large subunit") == 0 ||
              StringICmp (str, "ribulose 5-bisphosphate carboxylase, large subunit") == 0 ||
              StringICmp (str, "ribulosebisphosphate carboxylase large subunit") == 0 ||
              StringICmp (str, "ribulose bisphosphate large subunit") == 0 ||
              StringICmp (str, "ribulose 1,5 bisphosphate carboxylase/oxygenase large subunit") == 0 ||
              StringICmp (str, "ribulose 1,5-bisphosphate carboxylase/oxygenase large chain") == 0 ||
              StringICmp (str, "large subunit ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0 ||
              StringICmp (str, "ribulose-bisphosphate carboxylase, large subunit") == 0 ||
              StringICmp (str, "ribulose-1, 5-bisphosphate carboxylase/oxygenase large-subunit") == 0) {
            vnp->data.ptrvalue = StringSave ("ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit");
            str = MemFree (str);
          }
        }
      }
      CleanVisStringList (&(prp->name));
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            if (StringICmp (name, "its1") == 0 || StringICmp (name, "its 1") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 1");
            } else if (StringICmp (name, "its2") == 0 || StringICmp (name, "its 2") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 2");
            } else if (StringICmp (name, "its3") == 0 || StringICmp (name, "its 3") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 3");
            } else if (StringICmp (name, "Ribosomal DNA internal transcribed spacer 1") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 1");
            } else if (StringICmp (name, "Ribosomal DNA internal transcribed spacer 2") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 2");
            } else if (StringICmp (name, "Ribosomal DNA internal transcribed spacer 3") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 3");
            } else if (StringICmp (name, "internal transcribed spacer 1 (ITS1)") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 1");
            } else if (StringICmp (name, "internal transcribed spacer 2 (ITS2)") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 2");
            } else if (StringICmp (name, "internal transcribed spacer 3 (ITS3)") == 0) {
              gbq->val = MemFree (gbq->val);
              gbq->val = StringSave ("internal transcribed spacer 3");
            }
          }
        }
      }
      break;
    default:
      break;
  }
}

static Boolean IsFeatureCommentRedundant (SeqFeatPtr sfp)

{
  Uint1            aa;
  Choice           cbaa;
  CodeBreakPtr     cbp;
  CharPtr          comment;
  CdRegionPtr      crp;
  SeqFeatPtr       feat;
  Uint1            from;
  GBQualPtr        gbq;
  GeneRefPtr       grp;
  CharPtr          name;
  BioseqPtr        prod;
  ProtRefPtr       prp;
  Uint1            residue;
  RNAGenPtr        rgp;
  RNAQualPtr       rqp;
  RnaRefPtr        rrp;
  SeqAnnotPtr      sap;
  SeqCodeTablePtr  sctp;
  Uint1            seqcode;
  SeqIdPtr         sip;
  SeqMapTablePtr   smtp;
  CharPtr          str;
  tRNAPtr          trp;
  ValNodePtr       vnp;

  if (sfp == NULL) return FALSE;
  comment = sfp->comment;
  if (StringHasNoText (comment)) return FALSE;

  if (sfp->excpt && StringDoesHaveText (sfp->except_text)) {
    if (StringCmp (comment, sfp->except_text) == 0) return TRUE;
  }

  /* skip feature types that do not use data.value.ptrvalue */
  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
    case SEQFEAT_BOND:
    case SEQFEAT_SITE:
    case SEQFEAT_PSEC_STR:
      return FALSE;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return FALSE;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE:
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      /*
      if (StringCmp (comment, grp->locus) == 0) return TRUE;
      if (StringCmp (comment, grp->desc) == 0) return TRUE;
      */
      if (StringCmp (comment, grp->locus_tag) == 0) return TRUE;
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringCmp (comment, str) == 0) return TRUE;
      }
      break;
    case SEQFEAT_CDREGION:
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
        seqcode = 0;
        sctp = NULL;
        cbaa = cbp->aa;
        switch (cbaa.choice) {
          case 1 :
            seqcode = Seq_code_ncbieaa;
            break;
          case 2 :
            seqcode = Seq_code_ncbi8aa;
            break;
          case 3 :
            seqcode = Seq_code_ncbistdaa;
            break;
          default :
            break;
        }
        if (seqcode != 0) {
          sctp = SeqCodeTableFind (seqcode);
          if (sctp != NULL) {
            residue = cbaa.value.intvalue;
            if (residue != 42) {
              if (seqcode != Seq_code_ncbieaa) {
                smtp = SeqMapTableFind (seqcode, Seq_code_ncbieaa);
                residue = SeqMapTableConvert (smtp, residue);
              }
              if (residue == 'U') {
                if (StringCmp (comment, "selenocysteine") == 0) return TRUE;
              } else if (residue == 'O') {
                if (StringCmp (comment, "pyrrolysine") == 0) return TRUE;
              }
            }
          }
        }
      }
      if (sfp->product != NULL) {
        sip = SeqLocId (sfp->product);
        if (sip != NULL) {
          prod = BioseqFind (sip);
          if (prod != NULL) {
            for (sap = prod->annot; sap != NULL; sap = sap->next) {
              if (sap->type != 1) continue;
              for (feat = (SeqFeatPtr) sap->data; feat != NULL; feat = feat->next) {
                if (feat->data.choice != SEQFEAT_PROT) continue;
                prp = (ProtRefPtr) feat->data.value.ptrvalue;
                if (prp == NULL) continue;
                for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
                  str = (CharPtr) vnp->data.ptrvalue;
                  if (StringHasNoText (str)) continue;
                  if (StringCmp (comment, str) == 0) return TRUE;
                }
              }
            }
          }
        }
      }
      break;
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringCmp (comment, str) == 0) return TRUE;
      }
      if (StringDoesHaveText (prp->desc)) {
        if (StringCmp (comment, prp->desc) == 0) return TRUE;
      }
      for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (StringCmp (comment, str) == 0) return TRUE;
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            /*
            if (StringICmp (name, "internal transcribed spacer 1") == 0) {
              if (StringICmp (comment, "its1") == 0 || StringICmp (comment, "its 1") == 0) return TRUE;
            } else if (StringICmp (name, "internal transcribed spacer 2") == 0) {
              if (StringICmp (comment, "its2") == 0 || StringICmp (comment, "its 2") == 0) return TRUE;
            } else if (StringICmp (name, "internal transcribed spacer 3") == 0) {
              if (StringICmp (comment, "its3") == 0 || StringICmp (comment, "its 3") == 0) return TRUE;
            }
            */
          }
        }
      } else if (rrp->type == 3 && rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL) {
          aa = 0;
          if (trp->aatype == 2) {
            aa = trp->aa;
          } else {
            from = 0;
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
            seqcode = Seq_code_ncbieaa;
            smtp = SeqMapTableFind (seqcode, from);
            if (smtp != NULL) {
              aa = SeqMapTableConvert (smtp, trp->aa);
              if (aa == 255 && from == Seq_code_iupacaa) {
                if (trp->aa == 'U') {
                  aa = 'U';
                } else if (trp->aa == 'O') {
                  aa = 'O';
                }
              }
            }
          }
          if (aa > 0 && aa != 255) {
            if (StringNCmp (comment, "aa: ", 4) == 0) {
              comment += 4;
            }
            residue = FindTrnaAA3 (comment);
            if (residue == aa) {
              if (aa == 'M' && StringICmp ("fMet", comment) == 0) return FALSE;
              if (aa == 'M' && StringICmp ("iMet", comment) == 0) return FALSE;
              return TRUE;
            }
            residue = FindTrnaAA (comment);
            if (residue == aa) return TRUE;
          }
        }
      } else if (rrp->ext.choice == 3) {
        rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
        if (rgp != NULL) {
          if (StringCmp (comment, rgp->product) == 0) return TRUE;
          if (StringCmp (comment, rgp->_class) == 0) return TRUE;
          for (rqp = rgp->quals; rqp != NULL; rqp = rqp->next) {
            if (StringCmp (comment, rqp->val) == 0) return TRUE;
          }
        }
      }
      break;
    default:
      break;
  }

  return FALSE;
}


static CharPtr ExtractSatelliteFromComment (CharPtr comment)
{
  CharPtr satellite_type = NULL, satellite_start = NULL;
  CharPtr satellite_qual = NULL;
  Int4    satellite_len, i;

  if (StringHasNoText (comment)) {
    return NULL;
  }

  if (StringNCmp (comment, "microsatellite", 14) == 0) { 
    satellite_type = "microsatellite";
    satellite_start = comment;
  } else if (StringNCmp (comment, "minisatellite", 13) == 0) {
    satellite_type = "minisatellite";
    satellite_start = comment;
  } else if (StringNCmp (comment, "satellite", 9) == 0) {
    satellite_type = "satellite";
    satellite_start = comment;
  }

  if (satellite_start == NULL) {
    return NULL;
  }

  satellite_len = StringLen (satellite_type);
  if (comment[satellite_len] == '\0') {
    satellite_qual = StringSave (satellite_type);
    *comment = 0;
  } else if (comment[satellite_len] == ';') {
    satellite_qual = StringSave (satellite_type);
    for (i = 0; i <= satellite_len; i++) {
      comment [i] = ' ';
    }
    TrimSpacesAroundString (comment);
  }
  if (comment != NULL && comment [0] == '~' && comment [1] != '~') {
    comment [0] = ' ';
    TrimSpacesAroundString (comment);
  }

  return satellite_qual;
}

static void DoModernizeRNAFields (SeqFeatPtr sfp)

{
  RNAQualSetPtr       nextrqp;
  RNAQualSetPtr PNTR  prevrqp;
  RNAGenPtr           rgp;
  RNAQualSetPtr       rqp;
  RnaRefPtr           rrp;
  CharPtr             str;
  Boolean             unlink;
  Int2                i;
  size_t              len;
  CharPtr             ncclass;
  CharPtr             product;
  CharPtr             tmp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;

  ModernizeRNAFields (sfp);
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return;

  if (rrp->ext.choice == 1 && rrp->type == 10) {
    str = rrp->ext.value.ptrvalue;
    if (StringHasNoText (str)) return;

    rgp = (RNAGenPtr) MemNew (sizeof (RNAGen));
    if (rgp == NULL) return;
    rrp->ext.choice = 3;
    rrp->ext.value.ptrvalue = (Pointer) rgp;
    rgp->product = str;
  }

  if (rrp->ext.choice != 3) return;

  rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
  if (rgp == NULL) return;

  rqp = rgp->quals;
  prevrqp = (RNAQualSetPtr PNTR) &(rgp->quals);
  while (rqp != NULL) {
    nextrqp = rqp->next;
    unlink = FALSE;
    if (StringHasNoText (rqp->qual) || StringHasNoText (rqp->val)) {
      unlink = TRUE;
    }
    if (unlink) {
      *(prevrqp) = rqp->next;
      rqp->next = NULL;
      RNAQualFree (rqp);
    } else {
      prevrqp = (RNAQualSetPtr PNTR) &(rqp->next);
    }
    rqp = nextrqp;
  }

  if (rrp->type == 10 && StringDoesHaveText (rgp->product) && rgp->_class == NULL) {
    ncclass = rgp->product;
    for (i = 0; ncrnaClassList [i] != NULL; i++) {
      str = ncrnaClassList [i];
      if (StringHasNoText (str)) continue;
      len = StringLen (str);
      if (len < 1) continue;
      if (StringNICmp (ncclass, str, len) != 0) continue;
      if (ncclass [len] != ' ') continue;
      tmp = ncclass + len + 1;
      if (StringHasNoText (tmp)) continue;
      ncclass [len] = '\0';
      rgp->_class = StringSave (ncclass);
      product = StringSave (tmp);
      rgp->product = MemFree (rgp->product);
      rgp->product = product;
      TrimSpacesAroundString (rgp->_class);
      TrimSpacesAroundString (rgp->product);
      rrp->type = 8;
      sfp->idx.subtype = FEATDEF_ncRNA;
    }
  }

  if (rgp->quals != NULL) return;

  if (rrp->type == 2 || rrp->type == 4) {
    if (StringDoesHaveText (rgp->product) && StringHasNoText (rgp->_class)) {
      str = StringSave (rgp->product);
      rrp->ext.choice = 1;
      rrp->ext.value.ptrvalue = (Pointer) str;
      RNAGenFree (rgp);
      return;
    }
  }

  if (StringDoesHaveText (rgp->_class) || StringDoesHaveText (rgp->product)) return;

  rrp->ext.value.ptrvalue = NULL;
  rrp->ext.choice = 0;
  RNAGenFree (rgp);
}


static void FixncRNAClass (SeqFeatPtr sfp)
{
  RnaRefPtr rrp;
  RNAGenPtr rgp;

  if (sfp == NULL || sfp->idx.subtype != FEATDEF_ncRNA
      || (rrp = (RnaRefPtr) sfp->data.value.ptrvalue) == NULL
      || rrp->ext.choice != 3
      || (rgp = (RNAGenPtr) rrp->ext.value.ptrvalue) == NULL)
  {
    return;
  }

  if (StringICmp (rgp->_class, "antisense") == 0) {
    rgp->_class = MemFree (rgp->_class);
    rgp->_class = StringSave ("antisense_RNA");
  }
}


static void MoveBioSourceFeatureNoteToSubSourceNote (SeqFeatPtr sfp)
{
  ValNode vn;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || StringHasNoText (sfp->comment)) {
    return;
  }

  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = SourceQualChoice_textqual;
  vn.data.intvalue = Source_qual_subsource_note;

  SetSourceQualInBioSource (sfp->data.value.ptrvalue, &vn, NULL, sfp->comment, ExistingTextOption_append_semi);
  sfp->comment = MemFree (sfp->comment);
}


NLM_EXTERN void ConsolidateOneLikeSubSourceModifier (
  SubSourcePtr match_to,
  Boolean use_semicolon
)
{
  SubSourcePtr prev, index;
  Int4         len, num_matches;
  CharPtr      new_value;

  if (match_to == NULL) return;
  len = StringLen (match_to->name) + 1;
  num_matches = 0;
  prev = match_to;
  index = match_to->next;
  while (index != NULL)
  {
    if (index->subtype == match_to->subtype && index->name != NULL)
    {
      len += StringLen (index->name) + 2;
      num_matches++;
    }
    index = index->next;
  }
  if (num_matches == 0) return;

  new_value = MemNew (len * sizeof (char));
  if (new_value == NULL) return;

  StringCpy (new_value, match_to->name);
  index = match_to->next;
  while (index != NULL)
  {
    if (index->subtype == match_to->subtype && index->name != NULL)
    {
      if (use_semicolon)
      {
        StringCat (new_value, "; ");
      }
      else
      {
        StringCat (new_value, " ");
      }
      StringCat (new_value, index->name);
      prev->next = index->next;
      index->next = NULL;
      SubSourceFree (index);
      index = prev;
    }
    prev = index;
    index = index->next;
  }
  MemFree (match_to->name);
  match_to->name = new_value; 
}
  

NLM_EXTERN void ConsolidateOneLikeOrganismModifier (
  OrgModPtr match_to,
  Boolean use_semicolon
)
{
  OrgModPtr prev, index;
  Int4      len, num_matches;
  CharPtr   new_value;

  if (match_to == NULL) return;
  len = StringLen (match_to->subname) + 1;
  num_matches = 0;
  prev = match_to;
  index = match_to->next;
  while (index != NULL)
  {
    if (index->subtype == match_to->subtype && index->subname != NULL)
    {
      len += StringLen (index->subname) + 2;
      num_matches++;
    }
    index = index->next;
  }
  if (num_matches == 0) return;

  new_value = MemNew (len * sizeof (char));
  if (new_value == NULL) return;

  StringCpy (new_value, match_to->subname);
  index = match_to->next;
  while (index != NULL)
  {
    if (index->subtype == match_to->subtype && index->subname != NULL)
    {
      if (use_semicolon)
      {
        StringCat (new_value, "; ");
      }
      else
      {
        StringCat (new_value, " ");
      }
      StringCat (new_value, index->subname);
      prev->next = index->next;
      index->next = NULL;
      OrgModFree (index);
      index = prev;
    }
    prev = index;
    index = index->next;
  }
  MemFree (match_to->subname);
  match_to->subname = new_value; 
}
  
typedef struct reg_feat {
  CharPtr  feat_key;
  CharPtr  reg_class;
} RegFeatData, PNTR RegFeatPtr;

static RegFeatData reg_feat_keys [] = {
  { "enhancer",     "enhancer"              },
  { "promoter",     "promoter"              },
  { "CAAT_signal",  "CAAT_signal"           },
  { "TATA_signal",  "TATA_box"              },
  { "-35_signal",   "minus_35_signal"       },
  { "-10_signal",   "minus_10_signal"       },
  { "GC_signal",    "GC_signal"             },
  { "RBS",          "ribosome_binding_site" },
  { "polyA_signal", "polyA_signal_sequence" },
  { "attenuator",   "attenuator"            },
  { "terminator",   "terminator"            },
  { "misc_signal",  "other"                 },
  { NULL,           NULL                    }
};

NLM_EXTERN void ConsolidateBioSourceNotes (BioSourcePtr biop)
{
  SubSourcePtr ssp, note_ssp;
  OrgModPtr    mod, note_mod;

  if (biop == NULL) return;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next)
  {
    if (ssp->subtype == 255 && ssp->name != NULL)
    {
      ConsolidateOneLikeSubSourceModifier (ssp, TRUE);
      note_ssp = ssp;
    }
  }
    
  if (biop->org == NULL || biop->org->orgname == NULL) return;
  for (mod = biop->org->orgname->mod; mod != NULL; mod = mod->next)
  {
    if (mod->subtype == 255 && mod->subname != NULL)
    {
      ConsolidateOneLikeOrganismModifier (mod, TRUE);
      note_mod = mod;
    }
  }
}


NLM_EXTERN void CleanUpSeqFeat (
  SeqFeatPtr sfp,
  Boolean isEmblOrDdbj,
  Boolean isJscan,
  Boolean stripSerial,
  Boolean modernizeFeats,
  ValNodePtr PNTR publist
)

{
  BioseqPtr     bsp;
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  GBQualPtr     gbq;
  Boolean       emptyRNA;
  IntFuzzPtr    fuzz;
  GeneRefPtr    grp;
  Boolean       hasGibbsq;
  Boolean       hasNulls;
  SeqIdPtr      id;
  ImpFeatPtr    ifp;
  Int2          j;
  MolInfoPtr    mip;
  CharPtr       name;
  CharPtr       note;
  Boolean       partial5;
  Boolean       partial3;
  SeqPntPtr     pntp;
  Uint1         processed;
  ProtRefPtr    prp;
  ValNodePtr    psp;
  RNAGenPtr     rgp;
  RNAQualPtr    rqp;
  RnaRefPtr     rrp;
  Uint1         rrptype;
  CharPtr       satellite_type;
  SeqDescrPtr   sdp;
  SeqIntPtr     sintp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  CharPtr       str;
  Uint1         strand;
  Boolean       sync_mol_info;
  tRNAPtr       trp;
  SeqFeatXrefPtr  xref, next, PNTR prevlink;

  if (sfp == NULL) return;
  crp = NULL;
  if (sfp->data.choice == SEQFEAT_IMP) {
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (ifp != NULL) {
      if (ifp->loc != NULL) {
        str = StringStr (ifp->loc, "replace");
        if (str != NULL) {
          AddReplaceQual (sfp, str);
          ifp->loc = MemFree (ifp->loc);
        }
      }
      if (StringCmp (ifp->key, "CDS") == 0) {
        if (! isEmblOrDdbj) {
          sfp->data.value.ptrvalue = ImpFeatFree (ifp);
          sfp->data.choice = SEQFEAT_CDREGION;
          crp = CdRegionNew ();
          sfp->data.value.ptrvalue = crp;
          sfp->idx.subtype = FEATDEF_CDS;
        }
      } else if (StringCmp (ifp->key, "allele") == 0 ||
                 StringCmp (ifp->key, "mutation") == 0) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("variation");
        sfp->idx.subtype = FEATDEF_variation;
      } else if (StringCmp (ifp->key, "Import") == 0 ||
                 StringCmp (ifp->key, "virion") == 0) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("misc_feature");
        sfp->idx.subtype = FEATDEF_misc_feature;
      } else if (StringCmp (ifp->key, "repeat_unit") == 0 ) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("repeat_region");
        sfp->idx.subtype = FEATDEF_repeat_region;
      } else if (StringCmp (ifp->key, "misc_bind") == 0) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("misc_binding");
        sfp->idx.subtype = FEATDEF_misc_binding;
      } else if (StringCmp (ifp->key, "satellite") == 0 && (! isEmblOrDdbj)) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("repeat_region");
        sfp->idx.subtype = FEATDEF_repeat_region;
        gbq = GBQualNew ();
        if (gbq != NULL) {
          gbq->qual = StringSave ("satellite");
          gbq->val = ExtractSatelliteFromComment (sfp->comment);
          if (gbq->val == NULL) {
            gbq->val = StringSave ("satellite");
          }
          gbq->next = sfp->qual;
          sfp->qual = gbq;
        }
      } else if (StringCmp (ifp->key, "LTR") == 0) {
        ifp->key = MemFree (ifp->key);
        ifp->key = StringSave ("repeat_region");
        sfp->idx.subtype = FEATDEF_repeat_region;
        gbq = GBQualNew ();
        if (gbq != NULL) {
          gbq->qual = StringSave ("rpt_type");
          gbq->val = StringSave ("long_terminal_repeat");
          gbq->next = sfp->qual;
          sfp->qual = gbq;
        }
      } else if (StringHasNoText (ifp->loc)) {
        rrptype = 0;
        if (StringCmp (ifp->key, "precursor_RNA") == 0) {
          rrptype = 1;
        } else if (StringCmp (ifp->key, "mRNA") == 0) {
          rrptype = 2;
        } else if (StringCmp (ifp->key, "tRNA") == 0) {
          rrptype = 3;
        } else if (StringCmp (ifp->key, "rRNA") == 0) {
          rrptype = 4;
        } else if (StringCmp (ifp->key, "snRNA") == 0) {
          rrptype = 5;
        } else if (StringCmp (ifp->key, "scRNA") == 0) {
          rrptype = 6;
        } else if (StringCmp (ifp->key, "snoRNA") == 0) {
          rrptype = 7;
        } else if (StringCmp (ifp->key, "misc_RNA") == 0) {
          rrptype = 255;
        }
        if (rrptype != 0) {
          sfp->data.value.ptrvalue = ImpFeatFree (ifp);
          sfp->data.choice = SEQFEAT_RNA;
          rrp = RnaRefNew ();
          sfp->data.value.ptrvalue = rrp;
          rrp->type = rrptype;
          sfp->idx.subtype = FindFeatDefType (sfp);
        } else {
          processed = 0;
          if (StringCmp (ifp->key, "proprotein") == 0 || StringCmp (ifp->key, "preprotein") == 0) {
            processed = 1;
          } else if (StringCmp (ifp->key, "mat_peptide") == 0) {
            processed = 2;
          } else if (StringCmp (ifp->key, "sig_peptide") == 0) {
            processed = 3;
          } else if (StringCmp (ifp->key, "transit_peptide") == 0) {
            processed = 4;
          } else if (StringCmp (ifp->key, "propeptide") == 0 || StringCmp (ifp->key, "pro_peptide") == 0) {
            processed = 5;
          }
          if (processed != 0 || StringCmp (ifp->key, "Protein") == 0) {
            bsp = BioseqFind (SeqLocId (sfp->location));
            if (bsp != NULL && ISA_aa (bsp->mol)) {
              sfp->data.value.ptrvalue = ImpFeatFree (ifp);
              sfp->data.choice = SEQFEAT_PROT;
              prp = ProtRefNew ();
              sfp->data.value.ptrvalue = prp;
              prp->processed = processed;
              sfp->idx.subtype = FindFeatDefType (sfp);
            }
          }
        }
      }
      if (sfp->data.choice == SEQFEAT_IMP && StringCmp (ifp->key, "repeat_region") == 0 && (! isEmblOrDdbj)) {
        satellite_type = ExtractSatelliteFromComment (sfp->comment);
        if (satellite_type != NULL) {
          gbq = GBQualNew ();
          if (gbq != NULL) {
            gbq->qual = StringSave ("satellite");
            gbq->val = satellite_type;
            gbq->next = sfp->qual;
            sfp->qual = gbq;
          }
        }
      }
    }
  }
  if (crp != NULL && crp->frame == 0 && (! sfp->pseudo)) {
    crp->frame = GetFrameFromLoc (sfp->location);
  }
  if (sfp->data.choice == SEQFEAT_IMP) {
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (ifp != NULL) {
      for (j = 0; reg_feat_keys [j].feat_key != NULL; j++) {
        if (StringICmp (ifp->key, reg_feat_keys [j].feat_key) == 0) {
          ifp->key = MemFree (ifp->key);
          ifp->key = StringSave ("regulatory");
          sfp->idx.subtype = FEATDEF_regulatory;
          gbq = GBQualNew ();
          if (gbq != NULL) {
            gbq->qual = StringSave ("regulatory_class");
            gbq->val = StringSave (reg_feat_keys [j].reg_class);
            gbq->next = sfp->qual;
            sfp->qual = gbq;
          }
          break;
        }
      }
    }
  }
  if (sfp->data.choice == SEQFEAT_IMP) {
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (ifp != NULL && StringCmp (ifp->key, "regulatory") == 0) {
      note = NULL;
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "regulatory_class") != 0) continue;
        str = StringChr (gbq->val, ':');
        if (str == NULL) continue;
        if (StringNCmp (gbq->val, "other:", 6) == 0) continue;
        *str = '\0';
        str++;
        TrimSpacesAroundString (str);
        if (StringHasNoText (str)) continue;
        note = str;
      }
      if (StringDoesHaveText (note)) {
        gbq = GBQualNew ();
        if (gbq != NULL) {
          gbq->qual = StringSave ("note");
          gbq->val = StringSave (note);
          gbq->next = sfp->qual;
          sfp->qual = gbq;
        }
      }
    }
  }
  if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL) {
      if (rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringHasNoText (name)) {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.choice = 0;
        }
      } else if (rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL) {
          if (trp->aatype == 0 && trp->aa == 0 && trp->anticodon == NULL) {
            emptyRNA = TRUE;
            for (j = 0; j < 6; j++) {
              if (trp->codon [j] != 255) {
                emptyRNA = FALSE;
              }
            }
            if (emptyRNA) {
              rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
              rrp->ext.choice = 0;
            }
          }
        }
      } else if (rrp->ext.choice == 3) {
        rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
        if (rgp != NULL) {
          if (StringHasNoText (rgp->_class) && StringHasNoText (rgp->product)) {
            emptyRNA = TRUE;
            for (rqp = rgp->quals; rqp != NULL; rqp = rqp->next) {
              if (StringDoesHaveText (rqp->qual) && StringDoesHaveText (rqp->val)) {
                emptyRNA = FALSE;
              }
            } 
            if (emptyRNA) {
              rrp->ext.value.ptrvalue = RNAGenFree (rrp->ext.value.ptrvalue);
              rrp->ext.choice = 0;
            }
          }
        }
      }
    }
  }
  ModernizeFeatureGBQuals (sfp);
  sfp->qual = SortFeatureGBQuals (sfp->qual);
  CleanupDuplicateGBQuals (&(sfp->qual));
  CleanupFeatureGBQuals (sfp, isEmblOrDdbj);
  sfp->qual = SortIllegalGBQuals (sfp->qual);
  CleanupFeatureStrings (sfp, isJscan, isEmblOrDdbj, stripSerial, modernizeFeats, publist);
  FixOldDbxrefs (sfp->dbxref, isEmblOrDdbj);
  FixNumericDbxrefs (sfp->dbxref);
  sfp->dbxref = ValNodeSort (sfp->dbxref, SortDbxref);
  CleanupDuplicateDbxrefs (&(sfp->dbxref));
  CleanupObsoleteDbxrefs (&(sfp->dbxref));
  CleanupGoDbxrefs (sfp->dbxref);
  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    psp->data.ptrvalue = ValNodeSort ((ValNodePtr) psp->data.ptrvalue, SortCits);
    CleanupDuplicateCits ((ValNodePtr PNTR) &(psp->data.ptrvalue));
  }
  CleanUpSeqLoc (sfp->location);
  strand = SeqLocStrand (sfp->location);
  id = SeqLocId (sfp->location);
  if (sfp->data.choice == SEQFEAT_GENE) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL) {
      if (grp->pseudo) {
        sfp->pseudo = TRUE;
        grp->pseudo = FALSE;
      }
    }
  }
  if (sfp->data.choice == SEQFEAT_CDREGION) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      crp->code_break = SortCodeBreaks (sfp, crp->code_break);
      CleanupDuplicatedCodeBreaks (&(crp->code_break));
      for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
        CleanUpSeqLoc (cbp->loc);
        if (strand == Seq_strand_minus && id != NULL) {
          slp = cbp->loc;
          if (slp != NULL && slp->choice == SEQLOC_INT) {
            sip = SeqLocId (slp);
            if (sip != NULL && SeqIdComp (id, sip) == SIC_YES) {
              sintp = (SeqIntPtr) slp->data.ptrvalue;
              if (sintp != NULL) {
                sintp->strand = Seq_strand_minus;
              }
            }
          }
        }
      }
    }
  }
  if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    if (rrp != NULL) {
      if (rrp->pseudo) {
        sfp->pseudo = TRUE;
        rrp->pseudo = FALSE;
      }
    }
    if (rrp != NULL && rrp->ext.choice == 2) {
      trp = (tRNAPtr) rrp->ext.value.ptrvalue;
      if (trp != NULL && trp->anticodon != NULL) {
        CleanUpSeqLoc (trp->anticodon);
        if (strand == Seq_strand_minus && id != NULL) {
          slp = trp->anticodon;
          if (slp != NULL && slp->choice == SEQLOC_INT) {
            sip = SeqLocId (slp);
            if (sip != NULL && SeqIdComp (id, sip) == SIC_YES) {
              sintp = (SeqIntPtr) slp->data.ptrvalue;
              if (sintp != NULL) {
                sintp->strand = Seq_strand_minus;
              }
            }
          }
        }
      }
    }
    if (ConvertToNcRNA (sfp)) {
      sfp->idx.subtype = FindFeatDefType (sfp);
    }
    if (sfp->idx.subtype == FEATDEF_ncRNA) {
      FixncRNAClass (sfp);
    }       
  }
  if (sfp->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp != NULL && sfp->partial) {
      CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
      if (! partial5 && ! partial3) {
        bsp = BioseqFind (SeqLocId (sfp->location));
        if (bsp != NULL && ISA_aa (bsp->mol)) {
          hasGibbsq = FALSE;
          for (sip = bsp->id; sip != NULL; sip = sip->next) {
            if (sip->choice == SEQID_GIBBSQ) {
              hasGibbsq = TRUE;
            }
          }
          if (hasGibbsq) {
            sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_title, NULL);
            if (sdp != NULL && sdp->choice == Seq_descr_title) {
              str = (CharPtr) sdp->data.ptrvalue;
              if (StringDoesHaveText (str)) {
                sync_mol_info = FALSE;
                if (StringStr (str, "{N-terminal}") != NULL) {
                  partial3 = TRUE;
                  sync_mol_info = TRUE;
                } else if (StringStr (str, "{C-terminal}") != NULL) {
                  partial5 = TRUE;
                  sync_mol_info = TRUE;
                }
                if (sync_mol_info) {
                  SetSeqLocPartial (sfp->location, partial5, partial3);
                  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
                  if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
                    mip = (MolInfoPtr) sdp->data.ptrvalue;
                    if (mip != NULL) {
                      if (partial5 && partial3) {
                        mip->completeness = 5;
                      } else if (partial5) {
                        mip->completeness = 3;
                      } else if (partial3) {
                        mip->completeness = 4;
                      } else if (sfp->partial) {
                        mip->completeness = 2;
                      } else {
                        mip->completeness = 0;
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
  if (sfp->data.choice == SEQFEAT_REGION ||
             sfp->data.choice == SEQFEAT_SITE ||
             sfp->data.choice == SEQFEAT_BOND ||
             sfp->data.choice == SEQFEAT_PROT) {
    bsp = BioseqFind (SeqLocId (sfp->location));
    if (bsp != NULL && ISA_aa (bsp->mol)) {
      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        if (slp->choice == SEQLOC_INT) {
          sintp = (SeqIntPtr) slp->data.ptrvalue;
          if (sintp != NULL) {
            if (sintp->strand != Seq_strand_unknown) {
              sintp->strand = Seq_strand_unknown;
            }
          }
        } else if (slp->choice == SEQLOC_PNT) {
          pntp = (SeqPntPtr) slp->data.ptrvalue;
          if (pntp->strand != Seq_strand_unknown) {
            pntp->strand = Seq_strand_unknown;
          }
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }
    }
  }
  if (sfp->data.choice == SEQFEAT_BIOSRC) {
    /* combine multiple orgmod or subsource note qualifiers */
    ConsolidateBioSourceNotes(sfp->data.value.ptrvalue);
    /* if a BioSource feature has a comment, move the comment to
     * a subsource note.
     */
    MoveBioSourceFeatureNoteToSubSourceNote(sfp);
  }

  ModernizeFeatureStrings (sfp, isEmblOrDdbj);

  if (sfp->data.choice == SEQFEAT_GENE) {
    if (modernizeFeats) {
      ModernizeGeneFields (sfp);
    }
  }

  if (sfp->data.choice == SEQFEAT_RNA) {
    if (modernizeFeats) {
      DoModernizeRNAFields (sfp);
    }
  }

  if (IsFeatureCommentRedundant (sfp)) {
    sfp->comment = MemFree (sfp->comment);
  }

  /* sort and unique gbquals again after recent processing */
  sfp->qual = SortFeatureGBQuals (sfp->qual);
  CleanupDuplicateGBQuals (&(sfp->qual));
  sfp->qual = SortIllegalGBQuals (sfp->qual);

  /* normalize Seq-point fuzz tl to tr and decrement position */
  slp = SeqLocFindNext (sfp->location, NULL);
  for (slp = SeqLocFindNext (sfp->location, NULL);
       slp != NULL;
       slp = SeqLocFindNext (sfp->location, slp)) {
    if (slp->choice != SEQLOC_PNT) continue;
    pntp = (SeqPntPtr) slp->data.ptrvalue;
    if (pntp == NULL) continue;
    fuzz = pntp->fuzz;
    if (fuzz == NULL) continue;
    if (fuzz->choice == 4 /* lim */ && fuzz->a == 4 /* tl */ && pntp->point > 0) {
      (pntp->point)--;
      fuzz->a = 3; /* tr */
    }
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  hasNulls = LocationHasNullsBetween (sfp->location);
  sfp->partial = (sfp->partial || partial5 || partial3 || (hasNulls && ! isEmblOrDdbj));

  prevlink = (SeqFeatXrefPtr PNTR) &(sfp->xref);
  xref = sfp->xref;
  while (xref != NULL) {
    next = xref->next;

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


static void CleanUpSeqGraph (SeqGraphPtr sgp)

{
  if (sgp == NULL) return;
  if (sgp->loc != NULL) {
    CleanUpSeqLoc (sgp->loc);
  }
}

static void RemoveZeroLengthSeqLits (BioseqPtr bsp)
{
  DeltaSeqPtr dsp, prev = NULL, dsp_next;
  SeqLitPtr slip;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) {
    return;
  }

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp_next) {
    dsp_next = dsp->next;
    if (dsp->choice == 2 && (slip = (SeqLitPtr) (dsp->data.ptrvalue)) != NULL 
        && slip->length == 0 && slip->seq_data_type == 1
        && slip->seq_data != NULL) {
      if (prev == NULL) {
        bsp->seq_ext = dsp->next;
      } else {
        prev->next = dsp->next;
      }
      dsp->next = NULL;
      dsp = DeltaSeqFree (dsp);
    } else {
      prev = dsp;
    }
  }
}

/*
static Boolean CleanUpObjId (ObjectIdPtr oip)

{
  size_t   len;
  CharPtr  ptr;
  Boolean  rval = FALSE;
  long     val;

  if (oip == NULL) return FALSE;
  if (StringDoesHaveText (oip->str)) {
    if (isspace (oip->str[0]) || isspace (oip->str[StringLen (oip->str) - 1])) {
      TrimSpacesAroundString (oip->str);
      rval = TRUE;
    }
  }
  ptr = oip->str;
  if (ptr != NULL && *ptr != '0' && StringIsAllDigits(ptr)) {
    len = StringLen (ptr);
    if (len < 10 || (len == 10 && StringCmp (ptr, "2147483647") <= 0)) {
      if (sscanf (oip->str, "%ld", &val) == 1) {
        oip->id = (Int4) val;
        oip->str = MemFree (oip->str);
        rval = TRUE;
      }
    }
  }
  return rval;
}

static Boolean CleanUpSeqIdText (SeqIdPtr sip)
{
  DbtagPtr     dbt;
  ObjectIdPtr  oip;
  Boolean      rval = FALSE;

  if (sip == NULL) return FALSE;
  if (sip->choice == SEQID_LOCAL) {
    oip = (ObjectIdPtr) sip->data.ptrvalue;
    if (oip != NULL) {
      if (CleanUpObjId (oip)) {
        rval = TRUE;
      }
    }
  } else if (sip->choice == SEQID_GENERAL) {
    dbt = (DbtagPtr) sip->data.ptrvalue;
    if (dbt != NULL) {
      oip = dbt->tag;
      if (oip != NULL) {
        if (CleanUpObjId (oip)) {
          rval = TRUE;
        }
      }
    }
  }
  return rval;
}
*/


static Boolean CleanUpSeqIdText (SeqIdPtr sip)
{
  ObjectIdPtr  oip;
  Boolean      rval = FALSE;

  if (sip == NULL) return FALSE;
  if (sip->choice == SEQID_LOCAL) {
    oip = (ObjectIdPtr) sip->data.ptrvalue;
    if (oip != NULL) {
      if (StringDoesHaveText (oip->str)) {
        if (isspace (oip->str[0]) || isspace (oip->str[StringLen (oip->str) - 1])) {
       TrimSpacesAroundString (oip->str);
       rval = TRUE;
        }
      }
    }
  }
  return rval;
}

static void CleanUpSeqId (
  SeqIdPtr sip,
  Pointer userdata
)

{
  CleanUpSeqIdText (sip);
}

static void CleanSeqIdInBioseq (BioseqPtr bsp, Pointer userdata)

{
  SeqIdPtr sip;
  Boolean  need_reindex = FALSE;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (CleanUpSeqIdText (sip)) {
      need_reindex = TRUE;
    }
  }
  if (need_reindex) {
    SeqMgrReplaceInBioseqIndex (bsp);
  }
}

static void CleanSeqIdInSeqFeat (SeqFeatPtr sfp, Pointer userdata)

{
  VisitSeqIdsInSeqFeat (sfp, NULL, CleanUpSeqId);
}

static void CleanSeqIdInSeqAlign (SeqAlignPtr sap, Pointer userdata)

{
  VisitSeqIdsInSeqAlign (sap, NULL, CleanUpSeqId);
}

static void CleanSeqIdInSeqGraph (SeqGraphPtr sgp, Pointer userdata)

{
  VisitSeqIdsInSeqGraph (sgp, NULL, CleanUpSeqId);
}

static void CleanSeqIdInSeqAnnot (SeqAnnotPtr annot, Pointer userdata)

{
  VisitSeqIdsInSeqAnnot (annot, NULL, CleanUpSeqId);
}

typedef struct npcounts {
  Int4     nucs;
  Int4     prots;
  Boolean  make_genbank;
} NPCounts, PNTR NPCountsPtr;

static void CountNucsAndProts (BioseqPtr bsp, Pointer userdata)

{
  NPCountsPtr  ncp;

  if (bsp == NULL) return;
  ncp = (NPCountsPtr) userdata;
  if (ncp == NULL) return;

  if (ISA_na (bsp->mol)) {
    (ncp->nucs)++;
  } else if (ISA_aa (bsp->mol)) {
    (ncp->prots)++;
  }
}

static void CheckInnerSets (BioseqSetPtr bssp, Pointer userdata)

{
  NPCountsPtr  ncp;

  if (bssp == NULL) return;
  ncp = (NPCountsPtr) userdata;
  if (ncp == NULL) return;

  if (bssp->_class == BioseqseqSet_class_segset || bssp->_class == BioseqseqSet_class_parts) return;
  ncp->make_genbank = TRUE;
}

static void FixBadSetClass (BioseqSetPtr bssp, Pointer userdata)

{
  NPCounts  nc;

  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_not_set && bssp->_class != BioseqseqSet_class_other) return;

  MemSet ((Pointer) &nc, 0, sizeof (NPCounts));
  VisitSequencesInSet (bssp, (Pointer) &nc, VISIT_MAINS, CountNucsAndProts);
  VisitSetsInSet (bssp, (Pointer) &nc, CheckInnerSets);
  if (nc.nucs == 1 && nc.prots > 0 && (! nc.make_genbank)) {
    bssp->_class = BioseqseqSet_class_nuc_prot;
  } else {
    bssp->_class = BioseqseqSet_class_genbank;
  }
}

static void RemoveDuplicateSeqIds (BioseqPtr bsp)

{
  SeqIdPtr sip, sip_cmp, sip_prev, sip_next;

  if (bsp == NULL) {
    return;
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    sip_prev = sip;
    for (sip_cmp = sip->next; sip_cmp != NULL; sip_cmp = sip_next) {
      sip_next = sip_cmp->next;
      if (SeqIdComp (sip, sip_cmp) == SIC_YES) {
        sip_prev->next = sip_cmp->next;
        sip_cmp->next = NULL;
        sip_cmp = SeqIdFree (sip_cmp);
      } else {
        sip_prev = sip_cmp;
      }
    }
  }
}


static void BasicSeqEntryCleanupInternal (
  SeqEntryPtr sep,
  ValNodePtr PNTR publist,
  Boolean isEmblOrDdbj,
  Boolean isJscan,
  Boolean stripSerial
)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqDescrPtr   desc;
  Char          div [10];
  GBBlockPtr    gbp;
  MolInfoPtr    mip;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  SeqAnnotPtr   sap = NULL;
  ValNodePtr    sdp = NULL;
  SeqFeatPtr    sfp;
  SeqGraphPtr   sgp;
  SeqEntryPtr   tmp;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    /* remove duplicate SeqIds on the same Bioseq */
    RemoveDuplicateSeqIds (bsp);

    /* repair damaged delta sequences */
    RemoveZeroLengthSeqLits (bsp);

    sap = bsp->annot;
    sdp = bsp->descr;
    desc = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
    if (desc != NULL && desc->choice == Seq_descr_molinfo) {
      mip = (MolInfoPtr) desc->data.ptrvalue;
      if (mip != NULL) {
        /* repair if bsp.mol is not-set */
        if (bsp->mol == 0) {
          switch (mip->biomol) {
            case MOLECULE_TYPE_GENOMIC :
              bsp->mol = Seq_mol_na;
              break;
            case MOLECULE_TYPE_PRE_MRNA :
            case MOLECULE_TYPE_MRNA :
            case MOLECULE_TYPE_RRNA :
            case MOLECULE_TYPE_TRNA :
            case MOLECULE_TYPE_SNRNA :
            case MOLECULE_TYPE_SCRNA :
            case MOLECULE_TYPE_CRNA :
            case MOLECULE_TYPE_SNORNA :
            case MOLECULE_TYPE_TRANSCRIBED_RNA :
            case MOLECULE_TYPE_NCRNA :
            case MOLECULE_TYPE_TMRNA :
              bsp->mol = Seq_mol_rna;
              break;
            case MOLECULE_TYPE_PEPTIDE :
              bsp->mol = Seq_mol_aa;
              break;
            case MOLECULE_TYPE_OTHER_GENETIC_MATERIAL :
              bsp->mol = Seq_mol_other;
              break;
            case MOLECULE_TYPE_GENOMIC_MRNA_MIX :
              bsp->mol = Seq_mol_na;
              break;
            default :
              break;
          }
        } else if (bsp->mol != Seq_mol_rna 
                   && (mip->biomol == MOLECULE_TYPE_CRNA || mip->biomol == MOLECULE_TYPE_MRNA)) {
          bsp->mol = Seq_mol_rna;
        }
      }
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      BasicSeqEntryCleanupInternal (tmp, publist, isEmblOrDdbj, isJscan, stripSerial);
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
        CleanUpSeqFeat (sfp, isEmblOrDdbj, isJscan, stripSerial, TRUE, publist);
        sfp = sfp->next;
      }
    } else if (sap->type == 3) {
      sgp = (SeqGraphPtr) sap->data;
      while (sgp != NULL) {
        CleanUpSeqGraph (sgp);
        sgp = sgp->next;
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
        if (biop != NULL) {
          orp = biop->org;
        }
        break;
      default :
        break;
    }
    CleanupDescriptorStrings (sdp, stripSerial, TRUE, publist, isEmblOrDdbj);
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

static void ReplaceCitOnFeat (CitGenPtr cgp, ValNodePtr publist)

{
  ValNodePtr  nxt;
  ValNodePtr  vnp;

  for (vnp = publist; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != 1) continue;
    if (StringCmp (cgp->cit, (CharPtr) vnp->data.ptrvalue) == 0) {
      nxt = vnp->next;
      if (nxt != NULL && nxt->choice == 2) {
        cgp->cit = MemFree (cgp->cit);
        cgp->cit = StringSaveNoNull ((CharPtr) nxt->data.ptrvalue);
        if (cgp->cit != NULL) {
          if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
            cgp->cit [0] = 'U';
          }
        }
      }
      return;
    }
  }
}

static void ChangeCitsOnFeats (SeqFeatPtr sfp, Pointer userdata)

{
  CitGenPtr   cgp;
  ValNodePtr  ppr;
  ValNodePtr  psp;
  ValNodePtr  vnp;

  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
      vnp = NULL;
      if (ppr->choice == PUB_Gen) {
        vnp = ppr;
      } else if (ppr->choice == PUB_Equiv) {
        for (vnp = (ValNodePtr) ppr->data.ptrvalue;
             vnp != NULL && vnp->choice != PUB_Gen;
             vnp = vnp->next) continue;
      }
      if (vnp != NULL && vnp->choice == PUB_Gen) {
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL && (! StringHasNoText (cgp->cit))) {
          ReplaceCitOnFeat (cgp, (ValNodePtr) userdata);
        }
      }
    }
  }
}

static Int4 GetPmidForMuid (ValNodePtr pairlist, Int4 muid)

{
  ValNodePtr  vnp;

  vnp = pairlist;
  while (vnp != NULL) {
    if (muid == vnp->data.intvalue) {
      vnp = vnp->next;
      if (vnp == NULL) return 0;
      return vnp->data.intvalue;
    } else {
      vnp = vnp->next;
      if (vnp == NULL) return 0;
      vnp = vnp->next;
    }
  }

  return 0;
}

static void ChangeFeatCitsToPmid (SeqFeatPtr sfp, Pointer userdata)

{
  Int4        muid = 0;
  Int4        pmid = 0;
  ValNodePtr  ppr;
  ValNodePtr  psp;
  ValNodePtr  vnp;

  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
      vnp = NULL;
      if (ppr->choice == PUB_Muid) {
        vnp = ppr;
      } else if (ppr->choice == PUB_Equiv) {
        for (vnp = (ValNodePtr) ppr->data.ptrvalue;
             vnp != NULL && vnp->choice != PUB_Muid;
             vnp = vnp->next) continue;
      }
      if (vnp != NULL && vnp->choice == PUB_Muid) {
        muid = vnp->data.intvalue;
        if (muid != 0) {
          pmid = GetPmidForMuid ((ValNodePtr) userdata, muid);
          if (pmid != 0) {
            vnp->choice = PUB_PMid;
            vnp->data.intvalue = pmid;
          }
        }
      }
    }
  }
}

static void GetMuidPmidPairs (PubdescPtr pdp, Pointer userdata)

{
  Int4        muid = 0;
  Int4        pmid = 0;
  ValNodePtr  vnp;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Muid :
        muid = vnp->data.intvalue;
        break;
      case PUB_PMid :
        pmid = vnp->data.intvalue;
        break;
      default :
        break;
    }
  }
  if (muid == 0 || pmid == 0) return;
  ValNodeAddInt ((ValNodePtr PNTR) userdata, 0, muid);
  ValNodeAddInt ((ValNodePtr PNTR) userdata, 0, pmid);
}

static void FlattenPubSet (ValNodePtr PNTR prev)

{
  ValNodePtr  next;
  ValNodePtr  ppr;
  ValNodePtr  vnp;

  if (prev == NULL || *prev == NULL) return;
  ppr = *prev;
  while (ppr != NULL) {
    next = ppr->next;

    if (ppr->choice == PUB_Equiv) {
      vnp = (ValNodePtr) ppr->data.ptrvalue;
      if (vnp != NULL && vnp->next == NULL) {
        ppr->choice = vnp->choice;
        switch (vnp->choice) {
          case PUB_Muid :
          case PUB_PMid :
            ppr->data.intvalue = vnp->data.intvalue;
            break;
          default :
            ppr->data.ptrvalue = vnp->data.ptrvalue;
            break;
        }
        ValNodeFree (vnp);
      }
    }

    ppr = next;
  }
} 

static void FlattenDupInPubSet (ValNodePtr PNTR prev)

{
  ValNodePtr  next;
  ValNodePtr  nxt;
  ValNodePtr  ppr;
  ValNodePtr  vnp;

  if (prev == NULL || *prev == NULL) return;
  ppr = *prev;
  while (ppr != NULL) {
    next = ppr->next;

    if (ppr->choice == PUB_Equiv) {
      vnp = (ValNodePtr) ppr->data.ptrvalue;
      if (vnp != NULL) {
        nxt = vnp->next;
        if (nxt != NULL && nxt->next == NULL && vnp->choice == nxt->choice) {
          switch (vnp->choice) {
            case PUB_Muid :
            case PUB_PMid :
              if (vnp->data.intvalue == nxt->data.intvalue) {
                vnp->next = ValNodeFree (nxt);
              }
              break;
            default :
              break;
          }
        }
      }
    }

    ppr = next;
  }
} 

static void FlattenPubdesc (PubdescPtr pdp, Pointer userdata)

{
  FlattenPubSet (&(pdp->pub));
}

static void FlattenSfpCit (SeqFeatPtr sfp, Pointer userdata)

{
  ValNodePtr  psp;

  psp = sfp->cit;
  if (psp == NULL) return;
  FlattenDupInPubSet ((ValNodePtr PNTR) &(psp->data.ptrvalue));
  FlattenPubSet ((ValNodePtr PNTR) &(psp->data.ptrvalue));
}

typedef struct fastnode {
  ValNodePtr  head;
  ValNodePtr  tail;
} FastNode, PNTR FastNodePtr;

static void GetCitGenLabels (PubdescPtr pdp, Pointer userdata)

{
  Char             buf [121];
  CitGenPtr        cgp;
  FastNodePtr      labellist;
  ValNodePtr       tmp;
  ValNodePtr       vnp;
 
  if (pdp == NULL) return;
  labellist = (FastNodePtr) userdata;
  if (labellist == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != PUB_Gen) continue;
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp == NULL) continue;
    if (cgp->cit == NULL && cgp->journal == NULL &&
        cgp->date == NULL && cgp->serial_number) continue;
    PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE);
    tmp = ValNodeCopyStr (&(labellist->tail), 0, buf);
    if (labellist->head == NULL) {
      labellist->head = tmp;
    }
    labellist->tail = tmp;
  }
}

static void ReplaceShortCitGenOnFeat (CitGenPtr cgp, ValNodePtr labellist)

{
  Char        buf [128];
  Char        ch;
  size_t      len1;
  size_t      len2;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;

  for (vnp = labellist; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    len1 = StringLen (cgp->cit);
    if (len1 < 2 || len1 > 120) continue;
    StringCpy (buf, cgp->cit);
    ptr = StringStr (buf, "Unpublished");
    if (ptr != NULL) {
      ptr += 11;
      *ptr = '\0';
      tmp = StringStr (cgp->cit, "Unpublished");
      if (tmp != NULL) {
        tmp += 11;
        ch = *tmp;
        while (ch == ' ') {
          tmp++;
          ch = *tmp;
        }
        StringCat (buf, tmp);
      }
    }
    len1 = StringLen (buf);
    if (buf [len1 - 1] != '>') continue;
    len1--;
    len2 = StringLen (str);
    if (len1 >= len2) continue;
    if (StringNCmp (str, buf, len1) == 0) {
      cgp->cit = MemFree (cgp->cit);
      cgp->cit = StringSaveNoNull (str);
      if (cgp->cit != NULL) {
        if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
          cgp->cit [0] = 'U';
        }
      }
      return;
    }
  }
}

static void UpdateShortFeatCits (SeqFeatPtr sfp, Pointer userdata)

{
  CitGenPtr   cgp;
  ValNodePtr  ppr;
  ValNodePtr  psp;
  ValNodePtr  vnp;

  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    for (ppr = (ValNodePtr) psp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
      vnp = NULL;
      if (ppr->choice == PUB_Gen) {
        vnp = ppr;
      } else if (ppr->choice == PUB_Equiv) {
        for (vnp = (ValNodePtr) ppr->data.ptrvalue;
             vnp != NULL && vnp->choice != PUB_Gen;
             vnp = vnp->next) continue;
      }
      if (vnp != NULL && vnp->choice == PUB_Gen) {
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL && (! StringHasNoText (cgp->cit))) {
          ReplaceShortCitGenOnFeat (cgp, (ValNodePtr) userdata);
        }
      }
    }
  }
}

//LCOV_EXCL_START
NLM_EXTERN void BasicSeqAnnotCleanup (SeqAnnotPtr sap)

{
  SeqFeatPtr   sfp;
  SeqGraphPtr  sgp;

  if (sap == NULL) return;

  VisitSeqIdsInSeqAnnot (sap, NULL, CleanUpSeqId);

  if (sap->type == 1) {
    sfp = (SeqFeatPtr) sap->data;
    while (sfp != NULL) {
      CleanUpSeqFeat (sfp, FALSE, FALSE, TRUE, TRUE, NULL);
      sfp = sfp->next;
    }
  } else if (sap->type == 3) {
    sgp = (SeqGraphPtr) sap->data;
    while (sgp != NULL) {
      CleanUpSeqGraph (sgp);
      sgp = sgp->next;
    }
  }
}
//LCOV_EXCL_STOP

/*
static CharPtr proteinOrganellePrefix [] = {
  NULL,
  NULL,
  "chloroplast",
  "chromoplast",
  "kinetoplast",
  "mitochondrion",
  "plastid",
  "macronuclear",
  "extrachromosomal",
  "plasmid",
  NULL,
  NULL,
  "cyanelle",
  "proviral",
  "virus",
  "nucleomorph",
  "apicoplast",
  "leucoplast",
  "protoplast",
  "endogenous virus",
  "hydrogenosome",
  "chromosome",
  "chromatophore"
};
*/

static CharPtr proteinOrganellePrefix [] = {
  NULL,
  NULL,
  "chloroplast",
  "chromoplast",
  "kinetoplast",
  "mitochondrion",
  "plastid",
  "macronuclear",
  NULL,
  "plasmid",
  NULL,
  NULL,
  "cyanelle",
  NULL,
  NULL,
  "nucleomorph",
  "apicoplast",
  "leucoplast",
  "protoplast",
  "endogenous virus",
  "hydrogenosome",
  NULL,
  "chromatophore"
};

static CharPtr TitleEndsInOrganism (
  CharPtr title,
  CharPtr organism,
  CharPtr organelle,
  CharPtr PNTR onlp,
  BoolPtr case_diffp
)

{
  int      genome;
  size_t   len1, len2, len3;
  CharPtr  onl, ptr, tmp;

  if (onlp != NULL) {
    *onlp = NULL;
  }
  if (case_diffp != NULL) {
    *case_diffp = FALSE;
  }
  if (StringHasNoText (title) || StringHasNoText (organism)) return NULL;
  len1 = StringLen (title);
  len2 = StringLen (organism);
  if (len2 + 4 > len1) return NULL;

  tmp = title + len1 - len2 - 3;
  if (tmp [0] != ' ' || tmp [1] != '[' || tmp [len2 + 2] != ']') return NULL;
  if (StringNICmp (tmp + 2, organism, len2) != 0) return NULL;
  if (StringNCmp (tmp + 2, organism, len2) != 0 && case_diffp != NULL) {
    *case_diffp = TRUE;
  }

  if (onlp != NULL) {
    len3 = len1 - len2 - 3;
    for (genome = GENOME_chloroplast; genome <= GENOME_chromatophore; genome++) {
      ptr = proteinOrganellePrefix [genome];
      if (ptr == NULL) continue;
      len2 = StringLen (ptr);
      if (len2 + 4 >= len3) continue;
      onl = title + len3 - len2 - 3;
      if (onl [0] != ' ' || onl [1] != '(' || onl [len2 + 2] != ')') continue;
      if (StringNICmp (onl + 2, ptr, len2) != 0) continue;
      *onlp = onl;
      break;
    }
  }

  return tmp;
}

static void RemoveOrgFromEndOfProtein (SeqFeatPtr sfp, Pointer userdata)

{
  CharPtr     cp;
  size_t      len;
  ProtRefPtr  prp;
  CharPtr     str;
  CharPtr     taxname;
  ValNodePtr  vnp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;

  taxname = (CharPtr) userdata;
  if (StringHasNoText (taxname)) return;

  for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    len = StringLen (str);
    if (len < 5) continue;
    if (str [len - 1] != ']') continue;
    cp = StringRChr (str, '[');
    if (cp == NULL) continue;
    if (StringNCmp (cp, "[NAD", 4) == 0) continue;
    len = StringLen (taxname);
    if (StringLen (cp) != len + 2) continue;
    if (StringNICmp (cp + 1, taxname, len - 1) != 0) continue;
    *cp = '\0';
    TrimSpacesAroundString (str);
  }
}

static void AddPartialToProteinTitle (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr             binomial = NULL;
  BioSourcePtr        biop;
  BinomialOrgNamePtr  bonp;
  Boolean             case_difference = FALSE;
  CharPtr             first_super_kingdom = NULL;
  int                 genome = 0;
  CharPtr             genus = NULL;
  Boolean             is_cross_kingdom = FALSE;
  Boolean             is_wp = FALSE;
  size_t              len;
  MolInfoPtr          mip;
  Int2                num_super_kingdom = 0;
  CharPtr             oldname = NULL;
  OrgModPtr           omp;
  OrgNamePtr          onp;
  CharPtr             organelle = NULL;
  OrgRefPtr           orp;
  Boolean             partial = FALSE;
  CharPtr             penult = NULL;
  CharPtr             ptr;
  SeqDescrPtr         sdp;
  CharPtr             second_super_kingdom = NULL;
  SeqIdPtr            sip;
  CharPtr             species = NULL;
  CharPtr             str;
  CharPtr             suffix = NULL;
  Boolean             super_kingdoms_different = FALSE;
  CharPtr             taxname = NULL;
  TaxElementPtr       tep;
  CharPtr             title;
  CharPtr             tmp;
  TextSeqIdPtr        tsip;
  SeqDescrPtr         ttl = NULL;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_SWISSPROT) return;
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNICmp (tsip->accession, "WP_", 3) == 0) {
          is_wp = TRUE;
        }
      }
    }
  }

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && mip->completeness > 1 && mip->completeness < 6) {
      partial = TRUE;
    }
  }

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      genome = biop->genome;
      if (genome >= GENOME_chloroplast && genome <= GENOME_chromatophore) {
        organelle = proteinOrganellePrefix [genome];
      }
      orp = biop->org;
      if (orp != NULL) {
        taxname = orp->taxname;
        /*
        if (StringNICmp (organelle, taxname, StringLen (organelle)) == 0) {
          organelle = NULL;
        }
        */
        onp = orp->orgname;
        if (onp != NULL) {
          if (onp->choice == 1) {
            bonp = (BinomialOrgNamePtr) onp->data;
            if (bonp != NULL) {
              genus = bonp->genus;
              species = bonp->species;
            }
          }
          for (omp = onp->mod; omp != NULL; omp = omp->next) {
            if (omp->subtype == ORGMOD_old_name) {
              oldname = omp->subname;
            }
          }
        }
      }
    }
  }

  VisitFeaturesOnBsp (bsp, (Pointer) taxname, RemoveOrgFromEndOfProtein);

  ttl = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (ttl == NULL || ttl->choice != Seq_descr_title) return;
  str = (CharPtr) ttl->data.ptrvalue;
  if (StringHasNoText (str)) return;

  if (is_wp) {
    for (sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
         sdp != NULL;
         sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, sdp)) {
      if (sdp->choice != Seq_descr_source) continue;
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop == NULL) continue;
      orp = biop->org;
      if (orp == NULL) continue;
      onp = orp->orgname;
      if (onp == NULL) continue;
      if (onp->choice != 5) continue;
      for (tep = (TaxElementPtr) onp->data; tep != NULL; tep = tep->next) {
        if (tep->fixed_level == 0 && StringICmp (tep->level, "superkingdom") == 0) {
          num_super_kingdom++;
          if (first_super_kingdom == NULL) {
            first_super_kingdom = tep->name;
          } else if (StringICmp (first_super_kingdom, tep->name) != 0) {
            second_super_kingdom = tep->name;
            super_kingdoms_different = TRUE;
          }
          if (num_super_kingdom > 1 && super_kingdoms_different) {
            is_cross_kingdom = TRUE;
          }
        }
      }
    }
  }

  /* search for partial, must be just before parenthesized organelle or bracketed organism */
  tmp = StringSearch (str, ", partial [");
  if (tmp == NULL) {
    tmp = StringSearch (str, ", partial (");
  }

  /* find oldname or taxname in brackets at end of protein title */
  if (oldname != NULL && taxname != NULL) {
    suffix = TitleEndsInOrganism (str, oldname, organelle, &penult, &case_difference);
  }
  if (suffix == NULL && taxname != NULL) {
    suffix = TitleEndsInOrganism (str, taxname, organelle, &penult, &case_difference);
    if (suffix == NULL && StringDoesHaveText (genus) && StringDoesHaveText (species)) {
      len = StringLen (genus) + StringLen (species) + 5;
      binomial = (CharPtr) MemNew (len);
      if (binomial != NULL) {
        StringCpy (binomial, genus);
        StringCat (binomial, " ");
        StringCat (binomial, species);
        suffix = TitleEndsInOrganism (str, binomial, organelle, &penult, &case_difference);
      }
    }
    if (suffix == NULL && is_cross_kingdom) {
      ptr = StringStr (str, "][");
      if (ptr != NULL) {
        *(ptr + 1) = '\0';
        suffix = TitleEndsInOrganism (str, taxname, organelle, &penult, &case_difference);
      }
    } else {
      if (organelle == NULL && penult != NULL) {
      } else if (organelle != NULL && penult == NULL) {
      } else if (StringCmp (organelle, penult) != 0) {
      } else if (binomial != NULL) {
      } else if (case_difference) {
      } else {
        /* bail if no need to change partial text (organelle) [organism name] */
        if (partial) {
          if (tmp != NULL) return;
        } else {
          if (tmp == NULL) return;
        }
      }
    }
  }

  binomial = MemFree (binomial);

  /* do not change unless [genus species] was at the end */
  if (suffix == NULL) return;

  /* truncate bracketed info from end of title, will replace with current taxname */
  *suffix = '\0';
  suffix = taxname;

  /* truncate parenthesized info from just before bracketed taxname, will replace with current organelle */
  if (penult != NULL) {
    *penult = '\0';
  }

  /* if ", partial [/(" was indeed just before the [genus species] or (organelle), it will now be ", partial" */
  if (! partial && tmp != NULL && StringCmp (tmp, ", partial") == 0) {
    *tmp = '\0';
  }
  TrimSpacesAroundString (str);

  len = StringLen (str) + StringLen (organelle) + StringLen (suffix) + StringLen (first_super_kingdom) + StringLen (second_super_kingdom) + 20;
  title = MemNew (sizeof (Char) * len);
  if (title == NULL) return;

  StringCpy (title, str);
  if (partial && tmp == NULL) {
    StringCat (title, ", partial");
  }
  if (organelle != NULL) {
    StringCat (title, " (");
    StringCat (title, organelle);
    StringCat (title, ")");
  }
  if (is_cross_kingdom && StringDoesHaveText (first_super_kingdom) && StringDoesHaveText (second_super_kingdom)) {
    StringCat (title, " [");
    StringCat (title, first_super_kingdom);
    StringCat (title, "][");
    StringCat (title, second_super_kingdom);
    StringCat (title, "]");
  } else if (suffix != NULL) {
    StringCat (title, " [");
    StringCat (title, suffix);
    StringCat (title, "]");
  }
  MemFree (str);
  ttl->data.ptrvalue = title;
}

//LCOV_EXCL_START
NLM_EXTERN void CleanUpProteinTitles (SeqEntryPtr sep)

{
  if (sep == NULL) return;
  VisitBioseqsInSep (sep, NULL, AddPartialToProteinTitle);
}
//LCOV_EXCL_STOP

static void BasicSeqEntryCleanupEx (SeqEntryPtr sep, Boolean resync)

{
  AuthorPtr       ap;
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  Uint2           entityID;
  Boolean         isEmblOrDdbj = FALSE;
  Boolean         isJscan = FALSE;
  FastNode        labelnode;
  ValNodePtr      pairlist = NULL;
  ValNodePtr      publist = NULL;
  SeqEntryPtr     oldscope;
  ObjMgrDataPtr   omdp;
  SubmitBlockPtr  sbp;
  SeqSubmitPtr    ssp;
  Boolean         stripSerial = TRUE;

  if (sep == NULL) return;

  /* InGpsGenomic needs idx fields assigned */

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  AssignIDsInEntityEx (entityID, 0, NULL, NULL);

  /* HandleXrefOnCDS call to GetBestProteinFeatureUnindexed now scoped within record */

  oldscope = SeqEntrySetScope (sep);

  /* clean up spaces in local IDs */

  VisitBioseqsInSep (sep, NULL, CleanSeqIdInBioseq);
  VisitFeaturesInSep (sep, NULL, CleanSeqIdInSeqFeat);
  VisitAlignmentsInSep (sep, NULL, CleanSeqIdInSeqAlign);
  VisitGraphsInSep (sep, NULL, CleanSeqIdInSeqGraph);
  VisitAnnotsInSep (sep, NULL, CleanSeqIdInSeqAnnot);

  /* Fix Bioseq-sets with class 0 */

  VisitSetsInSep (sep, NULL, FixBadSetClass);

  /* removed unnecessarily nested Pub-equivs */

  VisitPubdescsInSep (sep, NULL, FlattenPubdesc);
  VisitFeaturesInSep (sep, NULL, FlattenSfpCit);

  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtID);
  SeqEntryExplore (sep, (Pointer) &isEmblOrDdbj, CheckForEmblDdbjID);
  SeqEntryExplore (sep, (Pointer) &isJscan, CheckForJournalScanID);
#ifdef SUPPRESS_STRIP_SERIAL_DIFFERENCES
  stripSerial = FALSE;
#endif

  BasicSeqEntryCleanupInternal (sep, &publist, isEmblOrDdbj, isJscan, stripSerial);
  if (publist != NULL) {
    VisitFeaturesInSep (sep, (Pointer) publist, ChangeCitsOnFeats);
  }
  ValNodeFreeData (publist);

  /* now get muid/pmid pairs, update sfp->cits to pmids */

  VisitPubdescsInSep (sep, (Pointer) &pairlist, GetMuidPmidPairs);
  if (pairlist != NULL) {
    VisitFeaturesInSep (sep, (Pointer) pairlist, ChangeFeatCitsToPmid);
  }
  ValNodeFree (pairlist);

  labelnode.head = NULL;
  labelnode.tail = NULL;
  VisitPubdescsInSep (sep, (Pointer) &labelnode, GetCitGenLabels);
  if (labelnode.head != NULL) {
    VisitFeaturesInSep (sep, (Pointer) labelnode.head, UpdateShortFeatCits);
  }
  ValNodeFreeData (labelnode.head);

  SeqEntrySetScope (oldscope);

  /* also normalize authors on submit block citation */

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL && ssp->datatype == 1) {
      sbp = ssp->sub;
      if (sbp != NULL) {
        csp = sbp->cit;
        if (csp != NULL) {
          NormalizeAuthors (csp->authors, TRUE);
        }
        cip = sbp->contact;
        if (cip != NULL) {
          ap = cip->contact;
          if (ap != NULL) {
            ap->affil = CleanAffil (ap->affil);
          }
        }
      }
    }
  }

  if (resync) {
    ResynchCodingRegionPartials (sep);
    ResynchMessengerRNAPartials (sep);
    ResynchProteinPartials (sep);
  }

  /*
  dynamically add missing partial to already instantiated protein
  titles, in between main title and bracketed organism name
  */

  VisitBioseqsInSep (sep, NULL, AddPartialToProteinTitle);
}

NLM_EXTERN void BasicSeqEntryCleanup (SeqEntryPtr sep)

{
  BasicSeqEntryCleanupEx (sep, FALSE);
}

//LCOV_EXCL_START
NLM_EXTERN void AdvancedSeqEntryCleanup (SeqEntryPtr sep)

{
  BasicSeqEntryCleanupEx (sep, TRUE);
}
//LCOV_EXCL_STOP

typedef struct bsecsmfedata {
  Int4  max;
  Int4  num_at_max;
} BsecSmfeData, PNTR BsecSmfePtr;

static Boolean LIBCALLBACK BsecSMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)


{
  BsecSmfePtr  bsp;
  Int4         len;

  if (sfp == NULL || context == NULL) return TRUE;
  bsp = context->userdata;
  if (bsp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < bsp->max) {
    bsp->max = len;
    bsp->num_at_max = 1;
  } else if (len == bsp->max) {
    (bsp->num_at_max)++;
  }

  return TRUE;
}

NLM_EXTERN void RemoveUnnecessaryGeneXrefs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BsecSmfeData         bsd;
  SeqFeatPtr           cds;
  Int2                 count;
  SeqFeatXrefPtr       curr, next;
  SeqMgrFeatContext    fcontext;
  GeneRefPtr           grp, grpx;
  SeqFeatXrefPtr PNTR  last;
  BioseqPtr            prd;
  SeqFeatPtr           sfpx;
  CharPtr              syn1, syn2;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;

  grpx = NULL;
  sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (sfpx != NULL) {
    if (sfpx->data.choice != SEQFEAT_GENE) return;
    grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  } else {
    prd = BioseqFindFromSeqLoc (sfp->location);
    if (prd != NULL && ISA_aa (prd->mol)) {
      cds = SeqMgrGetCDSgivenProduct (prd, NULL);
      if (cds != NULL) {
        grpx = SeqMgrGetGeneXref (cds);
        if (grpx == NULL) {
          sfpx = SeqMgrGetOverlappingGene (cds->location, &fcontext);
          if (sfpx != NULL && sfpx->data.choice == SEQFEAT_GENE) {
            grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
          }
        }
      }
    }
  }
  if (grpx == NULL || SeqMgrGeneIsSuppressed (grp)) return;

  if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
    if (StringICmp (grp->locus_tag, grpx->locus_tag) != 0) return;
  } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
    if (StringICmp (grp->locus, grpx->locus) != 0) return;
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
      if (StringICmp (syn1, syn2) != 0) return;
    }
  }

  MemSet ((Pointer) &bsd, 0, sizeof (BsecSmfeData));
  bsd.max = INT4_MAX;
  bsd.num_at_max = 0;
  count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE,
                                           NULL, 0, LOCATION_SUBSET,
                                           (Pointer) &bsd, BsecSMFEProc);

  if (bsd.num_at_max < 2) {
    last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
    curr = sfp->xref;
    while (curr != NULL) {
      next = curr->next;
      if (curr->data.choice == SEQFEAT_GENE) {
        *last = next;
        curr->next = NULL;
        SeqFeatXrefFree (curr);
      } else {
        last = &(curr->next);
      }
      curr = next;
    }
  }
}

//LCOV_EXCL_START
static void SortSeqFeatFields (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CdRegionPtr  crp;
  ValNodePtr   psp;

  if (sfp == NULL) return;

  sfp->qual = SortFeatureGBQuals (sfp->qual);

  sfp->qual = SortIllegalGBQuals (sfp->qual);

  sfp->dbxref = ValNodeSort (sfp->dbxref, SortDbxref);

  psp = sfp->cit;
  if (psp != NULL && psp->data.ptrvalue) {
    psp->data.ptrvalue = ValNodeSort ((ValNodePtr) psp->data.ptrvalue, SortCits);
  }

  if (sfp->data.choice == SEQFEAT_CDREGION) {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (crp != NULL) {
      crp->code_break = SortCodeBreaks (sfp, crp->code_break);
    }
  }
}

static void SortBioSourceFields (
  BioSourcePtr biop,
  Pointer userdata
)

{
  OrgNamePtr  onp;
  OrgRefPtr   orp;

  if (biop == NULL) return;

  orp = biop->org;
  if (orp != NULL) {
    orp->db = ValNodeSort (orp->db, SortDbxref);

    orp->syn = ValNodeSort (orp->syn, SortVnpByString);
    orp->syn = UniqueValNode (orp->syn);

    for (onp = orp->orgname; onp != NULL; onp = onp->next) {
      onp->mod = SortOrgModList (onp->mod);
    }
  }

  biop->subtype = SortSubSourceList (biop->subtype);
}

NLM_EXTERN void SortSeqEntryQualifiers (
  SeqEntryPtr sep
)

{
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, SortSeqFeatFields);
  VisitBioSourcesInSep (sep, NULL, SortBioSourceFields);
}
//LCOV_EXCL_STOP

/* end BasicSeqEntryCleanup section */

NLM_EXTERN void CDSPartialsFromTranslation (SeqFeatPtr sfp, Pointer userdata)

{
  Int4          i;
  Int4          len;
  ByteStorePtr  newprot;
  Boolean       partial5 = FALSE;
  Boolean       partial3 = TRUE;
  CharPtr       protseq;
  Int2          residue;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_CDREGION) return;

  newprot = ProteinFromCdRegionExEx (sfp, TRUE, FALSE, NULL, FALSE);
  if (newprot == NULL) return;

  protseq = BSMerge (newprot, NULL);
  if (protseq != NULL) {
    len = StringLen (protseq);

    for (i = 0; i < len; i++) {
      residue = protseq [i];
      if (i == 0 && residue == '-') {
        partial5 = TRUE;
      }
      if (i == len - 1 && residue == '*') {
        partial3 = FALSE;
      }
    }

    MemFree (protseq);

    SetSeqLocPartial (sfp->location, partial5, partial3);
    sfp->partial = (Boolean) (partial5 || partial3);
  }

  BSFree (newprot);
}

NLM_EXTERN void CodingRegionPartialsFromTranslation (SeqEntryPtr sep)

{
  VisitFeaturesInSep (sep, NULL, CDSPartialsFromTranslation);
}

NLM_EXTERN void ImposeGenePartials (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext, gcontext;
  SeqFeatPtr         feat, longest = NULL;
  Int4               len, min = INT4_MAX;
  Boolean            new_partial, partial5, partial3;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_GENE) return;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;

  if (SeqMgrGetDesiredFeature (0, bsp, 0, 0, sfp, &gcontext) != sfp) return;

  feat = SeqMgrGetDesiredFeature (0, bsp, 0, gcontext.index + 1, NULL, &fcontext);
  while (feat != NULL && gcontext.right >= fcontext.left) {
    len = TestFeatOverlap(feat, sfp, CONTAINED_WITHIN);
    if (len >= 0) {
      if (len < min) {
        min = len;
        longest = feat;
      }
    }
    feat = SeqMgrGetNextFeature (bsp, feat, 0, 0, &fcontext);
  }

  if (longest != NULL) {
    CheckSeqLocForPartial (longest->location, &partial5, &partial3);
    new_partial = (Boolean) (longest->partial || partial5 || partial3);
    SetSeqLocPartial (sfp->location, partial5, partial3);
    sfp->partial = new_partial;
  }
}

NLM_EXTERN void ImposeCDSPartials (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr  mrna;
  Boolean     new_partial, partial5, partial3;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_CDREGION) return;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  new_partial = (Boolean) (sfp->partial || partial5 || partial3);
  if (new_partial != sfp->partial) {
    sfp->partial = new_partial;
  }

  mrna = GetmRNAforCDS (sfp);
  if (mrna != NULL) {
    SetSeqLocPartial (mrna->location, partial5, partial3);
    mrna->partial = new_partial;
  }
}

NLM_EXTERN void ImposeCodingRegionPartials (SeqEntryPtr sep)

{
  VisitFeaturesInSep (sep, NULL, ImposeCDSPartials);
  VisitFeaturesInSep (sep, NULL, ImposeGenePartials);
}

NLM_EXTERN void ResynchCDSPartials (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr   bestprot;
  BioseqPtr    bsp;
  MolInfoPtr   mip;
  Boolean      partial5;
  Boolean      partial3;
  ProtRefPtr   prp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  ValNodePtr   vnp;
  /* variables for logging */
  LogInfoPtr    lip;
  CharPtr orig_loc = NULL, new_loc;
  Char    id_buf[100];
  Boolean new_partial;

  if (sfp->data.choice != SEQFEAT_CDREGION) return;
  lip = (LogInfoPtr) userdata;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  new_partial = (Boolean) (sfp->partial || partial5 || partial3);
  if (new_partial != sfp->partial) {
    sfp->partial = new_partial;
    if (lip != NULL) {
      lip->data_in_log = TRUE;
      if (lip->fp != NULL) {
        fprintf (lip->fp, "Changed partial flag for coding region\n");
      }
    }
  }

  /*
  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return;
  */
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp == NULL || !ISA_aa (bsp->mol) || bsp->repr != Seq_repr_raw) return;

  bestprot = SeqMgrGetBestProteinFeature (bsp, NULL);
  if (bestprot == NULL) {
    bestprot = GetBestProteinFeatureUnindexed (sfp->product);
  }

  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) return;

  /* only synchronize and extend best if unprocessed or preprotein, not mature/signal/transit peptide */
  if (bestprot != NULL && bestprot->location != NULL) {
    prp = (ProtRefPtr) bestprot->data.value.ptrvalue;
    slp = bestprot->location;
    if (prp != NULL && prp->processed < 2 && (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE)) {
      
      if (lip != NULL) {
        orig_loc = SeqLocPrintUseBestID (bestprot->location);
      }
      slp = NULL;
      sip = SeqLocId (bestprot->location);
      if (sip != NULL) {
        slp = WholeIntervalFromSeqId (sip);
      }
      if (slp == NULL) {
        slp = CreateWholeInterval (sep);
      }
      SetSeqLocPartial (slp, partial5, partial3);
      if (slp != NULL 
          && (!AsnIoMemComp (slp, bestprot->location, (AsnWriteFunc) SeqLocAsnWrite) || bestprot->partial != sfp->partial)) {
        bestprot->location = SeqLocFree (bestprot->location);
        bestprot->location = slp;
    
        bestprot->partial = sfp->partial;
        if (lip != NULL) {
          new_loc = SeqLocPrintUseBestID (bestprot->location);
          lip->data_in_log = TRUE;
          if (lip->fp != NULL) {
            fprintf (lip->fp, "Synchronized coding region partials for protein feature location at %s\n", orig_loc/*, new_loc*/);
          }
          new_loc = MemFree (new_loc);
        }
      } else {
        slp = SeqLocFree (slp);
      }
      orig_loc = MemFree (orig_loc);
    }
  }

  vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
  id_buf[0] = 0;
  if (vnp == NULL) {
    vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
    if (vnp != NULL) {
      mip = MolInfoNew ();
      vnp->data.ptrvalue = (Pointer) mip;
      if (mip != NULL) {
        mip->biomol = 8; /* peptide */
        mip->tech = 13; /* concept-trans-author */
        if (lip != NULL) {
          if (lip->fp != NULL) {
            SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
            fprintf (lip->fp, "Added MolInfo descriptor for %s\n", id_buf);
          }
          lip->data_in_log = TRUE;
        }
      }
    }
  }

  if (vnp != NULL && (mip = (MolInfoPtr) vnp->data.ptrvalue) != NULL) {
    if (partial5 && partial3) {
      if (mip->completeness != 5) {
        mip->completeness = 5;
        if (lip != NULL) {
          if (lip->fp != NULL) {
            if (id_buf[0] == 0) {
              SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
            }
            fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
            lip->data_in_log = TRUE;
          }
        }
      }
    } else if (partial5) {
      if (mip->completeness != 3) {
        mip->completeness = 3;
        if (lip != NULL) {
          if (lip->fp != NULL) {
            if (id_buf[0] == 0) {
              SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
            }
            fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
          }
          lip->data_in_log = TRUE;
        }
      }
    } else if (partial3) {
      if (mip->completeness != 4) {
        mip->completeness = 4;
        if (lip != NULL) {
          if (lip->fp != NULL) {
            if (id_buf[0] == 0) {
              SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
            }
            fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
          }
          lip->data_in_log = TRUE;
        }
      }
    } else if (sfp->partial) {
      if (mip->completeness != 2) {
        mip->completeness = 2;
        if (lip != NULL) {
          if (lip->fp != NULL) {
            if (id_buf[0] == 0) {
              SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
            }
            fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
          }
          lip->data_in_log = TRUE;
        }
      }
    } else {
      if (mip->completeness != 0 && mip->completeness != 1) {
        mip->completeness = 0;
        if (lip != NULL) {
          if (lip->fp != NULL) {
            if (id_buf[0] == 0) {
              SeqIdWrite (SeqIdFindBest (bsp->id, SEQID_GENBANK), id_buf, PRINTID_FASTA_SHORT, sizeof (id_buf) - 1);
            }
            fprintf (lip->fp, "Adjusted completeness for MolInfo descriptor on %s\n", id_buf);
          }
          lip->data_in_log = TRUE;
        }
      }
    }
  }
}


NLM_EXTERN Boolean ResynchCodingRegionPartialsEx (SeqEntryPtr sep, FILE *log_fp)

{
  LogInfoData lid;
  MemSet (&lid, 0, sizeof (LogInfoData));
  lid.fp = log_fp;
  VisitFeaturesInSep (sep, &lid, ResynchCDSPartials);
  return lid.data_in_log;
}

NLM_EXTERN void ResynchCodingRegionPartials (SeqEntryPtr sep)

{
  ResynchCodingRegionPartialsEx (sep, NULL);
}


NLM_EXTERN void ResynchMRNAPartials (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr    bsp;
  MolInfoPtr   mip;
  Boolean      partial5;
  Boolean      partial3;
  RnaRefPtr    rrp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  ValNodePtr   vnp;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL || rrp->type != 2) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sfp->partial = (Boolean) (sfp->partial || partial5 || partial3);
  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return;
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp != NULL && ISA_na (bsp->mol) && bsp->repr == Seq_repr_raw) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp == NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
      if (vnp != NULL) {
        mip = MolInfoNew ();
        vnp->data.ptrvalue = (Pointer) mip;
        if (mip != NULL) {
          mip->biomol = 3; /* mRNA */
          mip->tech = 1; /* standard */
        }
      }
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip != NULL) {
        if (partial5 && partial3) {
          mip->completeness = 5;
        } else if (partial5) {
          mip->completeness = 3;
        } else if (partial3) {
          mip->completeness = 4;
        } else if (sfp->partial) {
          mip->completeness = 2;
        } else {
          mip->completeness = 0;
        }
      }
    }
  }
}

NLM_EXTERN void ResynchMessengerRNAPartials (SeqEntryPtr sep)

{
  VisitFeaturesInSep (sep, NULL, ResynchMRNAPartials);
}

NLM_EXTERN void ResynchPeptidePartials (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr   bestprot;
  BioseqPtr    bsp;
  MolInfoPtr   mip;
  Boolean      partial5;
  Boolean      partial3;
  ProtRefPtr   prp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  ValNodePtr   vnp;

  if (sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;
  if (prp->processed < 1 || prp->processed > 5) return;
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  sfp->partial = (Boolean) (sfp->partial || partial5 || partial3);
  /*
  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return;
  */
  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp != NULL && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_raw) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep == NULL) return;
    bestprot = SeqMgrGetBestProteinFeature (bsp, NULL);
    if (bestprot == NULL) {
      bestprot = GetBestProteinFeatureUnindexed (sfp->product);
    }
    if (bestprot != NULL && bestprot->location != NULL) {
      /* only synchronize and extend best if unprocessed or preprotein, not mature/signal/transit peptide */
      prp = (ProtRefPtr) bestprot->data.value.ptrvalue;
      slp = bestprot->location;
      if (prp != NULL && prp->processed < 2 && (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE)) {
        slp = NULL;
        sip = SeqLocId (bestprot->location);
        if (sip != NULL) {
          slp = WholeIntervalFromSeqId (sip);
        }
        if (slp == NULL) {
          slp = CreateWholeInterval (sep);
        }
        if (slp != NULL) {
          bestprot->location = SeqLocFree (bestprot->location);
          bestprot->location = slp;
        }
        SetSeqLocPartial (bestprot->location, partial5, partial3);
        bestprot->partial = sfp->partial;
      }
    }
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    if (vnp == NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
      if (vnp != NULL) {
        mip = MolInfoNew ();
        vnp->data.ptrvalue = (Pointer) mip;
        if (mip != NULL) {
          mip->biomol = 8;
          mip->tech = 13;
        }
      }
    }
    if (vnp != NULL) {
      mip = (MolInfoPtr) vnp->data.ptrvalue;
      if (mip != NULL) {
        if (partial5 && partial3) {
          mip->completeness = 5;
        } else if (partial5) {
          mip->completeness = 3;
        } else if (partial3) {
          mip->completeness = 4;
        } else if (sfp->partial) {
          mip->completeness = 2;
        } else {
          mip->completeness = 0;
        }
      }
    }
  }
}

NLM_EXTERN void ResynchProteinPartials (SeqEntryPtr sep)

{
  VisitFeaturesInSep (sep, NULL, ResynchPeptidePartials);
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
      case SEQID_OTHER :
      case SEQID_TPG:
      case SEQID_TPE:
      case SEQID_TPD:
      case SEQID_GPIPE:
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

//LCOV_EXCL_START
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
//LCOV_EXCL_STOP

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

/*
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
*/

static Boolean CheckDataPath (CharPtr dirname, CharPtr subdir)

{
  if (FileExists (dirname, subdir, "seqcode.val")) return TRUE;
  return (Boolean) (FileExists (dirname, subdir, "objprt.prt"));
}

static Boolean CheckErrMsgPath (CharPtr dirname, CharPtr subdir)

{
  return (Boolean) (FileExists (dirname, subdir, "valid.msg"));
}

//LCOV_EXCL_START
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
  Boolean  dataFound;
  Char     path [PATH_MAX];
  Char     appPath[PATH_MAX];
  CharPtr  ptr;

  ProgramPath (appPath, sizeof (appPath));
  StrCpy(path, appPath);
  /* data a sibling of our application? */
  ptr = StringRChr (path, DIRDELIMCHR);
  if (ptr != NULL) {
    ptr++;
    *ptr = '\0';
  }
  dataFound = CheckDataPath (path, "data");
  if (! (dataFound)) {
  /* data an uncle of our application? */
    if (ptr != NULL) {
      ptr--;
      *ptr = '\0';
      ptr = StringRChr (path, DIRDELIMCHR);
      if (ptr != NULL) {
        ptr++;
        *ptr = '\0';
      }
      dataFound = CheckDataPath (path, "data");
    }
  }
#ifdef OS_UNIX_DARWIN
  if (! (dataFound) && IsApplicationPackage (appPath)) {
      /* is data inside our application within Contents/Resources? */
      StrCpy (path, appPath);
      FileBuildPath (path, "Contents", NULL);
      FileBuildPath (path, "Resources", NULL);
      dataFound = CheckDataPath (path, "data");
      if (! dataFound) {
        StrCpy (path, appPath);
        ptr = StringStr (path, "/ncbi/build/");
        if (ptr != NULL) {
          /* see if running under older Xcode 3 build environment */
          ptr [5] = '\0';
          dataFound = CheckDataPath (path, "data");
        }
      }
      if (! dataFound) {
        StrCpy (path, appPath);
        ptr = StringStr (path, "/ncbi/make/");
        if (ptr != NULL) {
          /* see if running under newer Xcode 3 build environment */
          ptr [5] = '\0';
          dataFound = CheckDataPath (path, "data");
        }
      }
      if (! dataFound) {
          StrCpy (path, appPath);
          ptr = StringStr (path, "/Library/Developer/");
          if (ptr != NULL) {
              /* see if running under Xcode 4 build environment */
              ptr [19] = '\0';
              dataFound = CheckDataPath (path, "data");
          }
      }
  }
#endif
  if (dataFound) {
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
        if (ISA_na (bsp->mol)) {
          sip->strand = Seq_strand_plus;
        }
        sip->id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
      }
    }
  }
  return slp;
}
//LCOV_EXCL_STOP


NLM_EXTERN SeqLocPtr WholeIntervalFromSeqId (SeqIdPtr sip)

{
  BioseqPtr  bsp;
  SeqIntPtr  sintp;
  SeqLocPtr  slp;

  if (sip == NULL) return NULL;
  bsp = BioseqFindCore (sip);
  if (bsp == NULL) return NULL;
  slp = ValNodeNew (NULL);
  if (slp == NULL) return NULL;
  sintp = SeqIntNew ();
  if (sintp == NULL) return NULL;
  slp->choice = SEQLOC_INT;
  slp->data.ptrvalue = (Pointer) sintp;
  sintp->from = 0;
  sintp->to = bsp->length - 1;
  if (ISA_na (bsp->mol)) {
    sintp->strand = Seq_strand_plus;
  }
  sintp->id = SeqIdStripLocus (SeqIdDup (sip));
  return slp;
}

//LCOV_EXCL_START
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
//LCOV_EXCL_STOP

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

NLM_EXTERN void NormalizeNullsBetween (SeqLocPtr location)

{
  SeqLocPtr  next, tmp, vnp;

  if (location == NULL) return;
  if (! LocationHasNullsBetween (location)) return;

  if (location->choice != SEQLOC_MIX) return;
  vnp = (ValNodePtr) location->data.ptrvalue;
  if (vnp == NULL) return;

  while (vnp != NULL && vnp->next != NULL) {
    next = vnp->next;
    if (vnp->choice != SEQLOC_NULL && next->choice != SEQLOC_NULL) {
      tmp = ValNodeNew (NULL);
      if (tmp != NULL) {
        tmp->choice = SEQLOC_NULL;
        tmp->next = vnp->next;
        vnp->next = tmp;
      }
    }
    vnp = next;
  }
}

NLM_EXTERN Uint1 FindFeatFromFeatDefType (Uint2 subtype)

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
      if (subtype == FEATDEF_snoRNA) {
        return SEQFEAT_RNA;
      }
      if (subtype >= FEATDEF_ncRNA && subtype <= FEATDEF_tmRNA) {
        return SEQFEAT_RNA;
      }
      if (subtype >= FEATDEF_preprotein && subtype <= FEATDEF_transit_peptide_aa) {
        return SEQFEAT_PROT;
      }
      if (subtype >= FEATDEF_IMP && subtype <= FEATDEF_site_ref) {
        return SEQFEAT_IMP;
      }
      if (subtype >= FEATDEF_gap && subtype <= FEATDEF_oriT) {
        return SEQFEAT_IMP;
      }
      if (subtype >= FEATDEF_mobile_element && subtype <= FEATDEF_propeptide) {
        return SEQFEAT_IMP;
      }
      if (subtype == FEATDEF_propeptide_aa) {
        return SEQFEAT_PROT;
      }
  }
  return 0;
}

//LCOV_EXCL_START
NLM_EXTERN SeqIdPtr MakeSeqID(CharPtr str)

{
  CharPtr   buf;
  Int4      len;
  SeqIdPtr  sip;

  sip = NULL;
  if (str != NULL && *str != '\0') {
    if (StringChr (str, '|') != NULL) {
      sip = SeqIdParse (str);
    } else {
      len = StringLen (str) + 5;
      buf = (CharPtr) MemNew (sizeof (Char) * len);
      sprintf (buf, "lcl|%s", str);
      sip = SeqIdParse (buf);
      buf = MemFree (buf);
    }
  }
  return sip;
}

NLM_EXTERN SeqIdPtr MakeUniqueSeqID (CharPtr prefix)

{
    Char buf[60];
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
    if (len > 0 && len < 52) {
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

NLM_EXTERN SeqIdPtr SeqIdFindWorst (SeqIdPtr sip)

{
  Uint1  order [NUM_SEQID];

  SeqIdBestRank (order, NUM_SEQID);
  order [SEQID_LOCAL] = 10;
  order [SEQID_GENBANK] = 5;
  order [SEQID_EMBL] = 5;
  order [SEQID_PIR] = 5;
  order [SEQID_SWISSPROT] = 5;
  order [SEQID_DDBJ] = 5;
  order [SEQID_PRF] = 5;
  order [SEQID_PDB] = 5;
  order [SEQID_TPG] = 5;
  order [SEQID_TPE] = 5;
  order [SEQID_TPD] = 5;
  order [SEQID_GPIPE] = 9;
  order [SEQID_NAMED_ANNOT_TRACK] = 9;
  order [SEQID_PATENT] = 10;
  order [SEQID_OTHER] = 8;
  order [SEQID_GENERAL] = 15;
  order [SEQID_GIBBSQ] = 15;
  order [SEQID_GIBBMT] = 15;
  order [SEQID_GIIM] = 20;
  order [SEQID_GI] = 20;
  return SeqIdSelect (sip, order, NUM_SEQID);
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
          vnp = SeqDescrNew (descr);
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
            vnp = SeqDescrNew (descr);
            if (descr == NULL) {
              bssp->descr = vnp;
            }
            seqentry = NULL;
          } else if ((_class >= 5 && _class <= 8) || _class == 11 /* || _class == BioseqseqSet_class_gen_prod_set */) {
            seqentry = bssp->seq_set;
          } else {
            descr = bssp->descr;
            vnp = SeqDescrNew (descr);
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


/* common functions to scan binary ASN.1 file of entire release as Bioseq-set */

static Int4 VisitSeqIdList (SeqIdPtr sip, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4  index = 0;

  while (sip != NULL) {
    if (callback != NULL) {
      callback (sip, userdata);
    }
    index++;
    sip = sip->next;
  }
  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqLoc (SeqLocPtr slp, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4           index = 0;
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;

  if (slp == NULL) return index;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        index += VisitSeqIdList (sip, userdata, callback);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          index += VisitSeqIdList (sip, userdata, callback);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          index += VisitSeqIdList (sip, userdata, callback);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          index += VisitSeqIdList (sip, userdata, callback);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          index += VisitSeqIdsInSeqLoc (loc, userdata, callback);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            index += VisitSeqIdList (sip, userdata, callback);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            index += VisitSeqIdList (sip, userdata, callback);
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

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInBioseq (BioseqPtr bsp, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;

  if (bsp->id != NULL) {
    index += VisitSeqIdList (bsp->id, userdata, callback);
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqFeat (SeqFeatPtr sfp, Pointer userdata, VisitSeqIdFunc callback)

{
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  Int4          index = 0;
  RnaRefPtr     rrp;
  tRNAPtr       trp;

  if (sfp == NULL) return index;

  index += VisitSeqIdsInSeqLoc (sfp->location, userdata, callback);
  index += VisitSeqIdsInSeqLoc (sfp->product, userdata, callback);

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          index += VisitSeqIdsInSeqLoc (cbp->loc, userdata, callback);
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 2) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp != NULL && trp->anticodon != NULL) {
          index += VisitSeqIdsInSeqLoc (trp->anticodon, userdata, callback);
        }
      }
      break;
    default :
      break;
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqAlign (SeqAlignPtr sap, Pointer userdata, VisitSeqIdFunc callback)

{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  Int4          index = 0;
  SeqIdPtr      sip;
  SeqLocPtr     slp = NULL;
  StdSegPtr     ssp;

  if (sap == NULL) return index;

  if (sap->bounds != NULL) {
    sip = SeqLocId (sap->bounds);
    index += VisitSeqIdList (sip, userdata, callback);
  }

  if (sap->segs == NULL) return index;

  switch (sap->segtype) {
    case SAS_DENDIAG :
      ddp = (DenseDiagPtr) sap->segs;
      if (ddp != NULL) {
        for (sip = ddp->id; sip != NULL; sip = sip->next) {
          index += VisitSeqIdList (sip, userdata, callback);
        }
      }
      break;
    case SAS_DENSEG :
      dsp = (DenseSegPtr) sap->segs;
      if (dsp != NULL) {
        for (sip = dsp->ids; sip != NULL; sip = sip->next) {
          index += VisitSeqIdList (sip, userdata, callback);
        }
      }
      break;
    case SAS_STD :
      ssp = (StdSegPtr) sap->segs;
      for (slp = ssp->loc; slp != NULL; slp = slp->next) {
        sip = SeqLocId (slp);
        index += VisitSeqIdList (sip, userdata, callback);
      }
      break;
    case SAS_DISC :
      /* recursive */
      for (sap = (SeqAlignPtr) sap->segs; sap != NULL; sap = sap->next) {
        index += VisitSeqIdsInSeqAlign (sap, userdata, callback);
      }
      break;
    default :
      break;
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqGraph (SeqGraphPtr sgp, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4      index = 0;
  SeqIdPtr  sip;

  if (sgp == NULL) return index;

  if (sgp->loc != NULL) {
    sip = SeqLocId (sgp->loc);
    index += VisitSeqIdList (sip, userdata, callback);
  }

  return index;
}

NLM_EXTERN Int4 VisitSeqIdsInSeqAnnot (SeqAnnotPtr annot, Pointer userdata, VisitSeqIdFunc callback)

{
  Int4         index = 0;
  SeqAlignPtr  sap;
  SeqFeatPtr   sfp;
  SeqGraphPtr  sgp;

  if (annot == NULL || annot->data == NULL) return index;

  switch (annot->type) {

    case 1 :
      for (sfp = (SeqFeatPtr) annot->data; sfp != NULL; sfp = sfp->next) {
        index += VisitSeqIdsInSeqFeat (sfp, userdata, callback);
      }
      break;

    case 2 :
      for (sap = (SeqAlignPtr) annot->data; sap != NULL; sap = sap->next) {
        index += VisitSeqIdsInSeqAlign (sap, userdata, callback);
      }
      break;

    case 3 :
      for (sgp = (SeqGraphPtr) annot->data; sgp != NULL; sgp = sgp->next) {
        index += VisitSeqIdsInSeqGraph (sgp, userdata, callback);
      }
      break;

    default :
      break;
  }

  return index;
}

NLM_EXTERN Int4 VisitUserFieldsInUfp (UserFieldPtr ufp, Pointer userdata, VisitUserFieldsFunc callback)

{
  UserFieldPtr  curr;
  Int4          index = 0;
  Boolean       nested = FALSE;

  if (ufp == NULL) return index;
  if (ufp->choice == 11) {
    for (curr = (UserFieldPtr) ufp->data.ptrvalue; curr != NULL; curr = curr->next) {
      index += VisitUserFieldsInUfp (curr, userdata,callback);
      nested = TRUE;
    }
  }
  if (! nested) {
    if (callback != NULL) {
      callback (ufp, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitUserFieldsInUop (UserObjectPtr uop, Pointer userdata, VisitUserFieldsFunc callback)

{
  Int4          index = 0;
  UserFieldPtr  ufp;

  if (uop == NULL) return index;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (callback != NULL) {
      callback (ufp, userdata);
    }
    index++;
  }
  return index;
}

/* Visits only unnested nodes */
NLM_EXTERN Int4 VisitUserObjectsInUop (UserObjectPtr uop, Pointer userdata, VisitUserObjectFunc callback)

{
  Int4           index = 0;
  Boolean        nested = FALSE;
  UserObjectPtr  obj;
  UserFieldPtr   ufp;

  if (uop == NULL) return index;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->choice == 6) {
      obj = (UserObjectPtr) ufp->data.ptrvalue;
      index += VisitUserObjectsInUop (obj, userdata, callback);
      nested = TRUE;
    } else if (ufp->choice == 12) {
      for (obj = (UserObjectPtr) ufp->data.ptrvalue; obj != NULL; obj = obj->next) {
        index += VisitUserObjectsInUop (obj, userdata, callback);
      }
      nested = TRUE;
    }
  }
  if (! nested) {
    if (callback != NULL) {
      callback (uop, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitAllUserObjectsInUop (UserObjectPtr uop, Pointer userdata, VisitUserObjectFunc callback)

{
  Int4           index = 0;
  UserObjectPtr  obj;
  UserFieldPtr   ufp;

  if (uop == NULL) return index;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->choice == 6) {
      obj = (UserObjectPtr) ufp->data.ptrvalue;
      index += VisitAllUserObjectsInUop (obj, userdata, callback);
    } else if (ufp->choice == 12) {
      for (obj = (UserObjectPtr) ufp->data.ptrvalue; obj != NULL; obj = obj->next) {
        index += VisitAllUserObjectsInUop (obj, userdata, callback);
      }
    }
  }
  if (callback != NULL) {
    callback (uop, userdata);
  }
  index++;
  return index;
}
//LCOV_EXCL_STOP

typedef struct uopdata {
  UserObjectPtr  rsult;
  CharPtr        tag;
} UopData, PNTR UopDataPtr;

static void FindUopProc (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;
  UopDataPtr   udp;

  if (uop == NULL || userdata == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  udp = (UopDataPtr) userdata;
  if (StringICmp (oip->str, udp->tag) != 0) return;
  udp->rsult = uop;
}

NLM_EXTERN UserObjectPtr FindUopByTag (UserObjectPtr top, CharPtr tag)

{
  UopData  ud;

  if (top == NULL || StringHasNoText (tag)) return NULL;
  ud.rsult = NULL;
  ud.tag = tag;
  VisitUserObjectsInUop (top, (Pointer) &ud, FindUopProc);
  return ud.rsult;
}

//LCOV_EXCL_START
NLM_EXTERN UserObjectPtr CombineUserObjects (UserObjectPtr origuop, UserObjectPtr newuop)

{
  UserFieldPtr   prev = NULL;
  ObjectIdPtr    oip;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;

  if (newuop == NULL) return origuop;
  if (origuop == NULL) return newuop;

  /* adding to an object that already chaperones at least two user objects */

  oip = origuop->type;
  if (oip != NULL && StringICmp (oip->str, "CombinedFeatureUserObjects") == 0) {

    for (ufp = origuop->data; ufp != NULL; ufp = ufp->next) {
      prev = ufp;
    }

    ufp = UserFieldNew ();
    oip = ObjectIdNew ();
    oip->id = 0;
    ufp->label = oip;
    ufp->choice = 6; /* user object */
    ufp->data.ptrvalue = (Pointer) newuop;

    /* link new set at end of list */

    if (prev != NULL) {
      prev->next = ufp;
    } else {
      origuop->data = ufp;
    }
    return origuop;
  }

  /* creating a new chaperone, link in first two user objects */

  uop = UserObjectNew ();
  oip = ObjectIdNew ();
  oip->str = StringSave ("CombinedFeatureUserObjects");
  uop->type = oip;

  ufp = UserFieldNew ();
  oip = ObjectIdNew ();
  oip->id = 0;
  ufp->label = oip;
  ufp->choice = 6; /* user object */
  ufp->data.ptrvalue = (Pointer) origuop;
  uop->data = ufp;
  prev = ufp;

  ufp = UserFieldNew ();
  oip = ObjectIdNew ();
  oip->id = 0;
  ufp->label = oip;
  ufp->choice = 6; /* user object */
  ufp->data.ptrvalue = (Pointer) newuop;
  prev->next = ufp;

  return uop;
}


static Int4 VisitDescriptorsProc (SeqDescrPtr descr, Pointer userdata, VisitDescriptorsFunc callback)

{
  Int4         index = 0;
  SeqDescrPtr  sdp;

  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (callback != NULL) {
      callback (sdp, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsOnBsp (BioseqPtr bsp, Pointer userdata, VisitDescriptorsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitDescriptorsProc (bsp->descr, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitDescriptorsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitDescriptorsProc (bssp->descr, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsInSet (BioseqSetPtr bssp, Pointer userdata, VisitDescriptorsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitDescriptorsProc (bssp->descr, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitDescriptorsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsOnSep (SeqEntryPtr sep, Pointer userdata, VisitDescriptorsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitDescriptorsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitDescriptorsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitDescriptorsInSep (SeqEntryPtr sep, Pointer userdata, VisitDescriptorsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitDescriptorsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitDescriptorsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitFeaturesProc (SeqAnnotPtr annot, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4         index = 0;
  SeqAnnotPtr  sap;
  SeqFeatPtr   sfp;

  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 1) continue;
    for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
      if (callback != NULL) {
        callback (sfp, userdata);
      }
      index++;
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitFeaturesOnSap (SeqAnnotPtr sap, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4        index = 0;
  SeqFeatPtr  sfp;

  if (sap == NULL) return index;
  if (sap->type != 1) return index;
  for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
    if (callback != NULL) {
      callback (sfp, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitFeaturesOnBsp (BioseqPtr bsp, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitFeaturesProc (bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitFeaturesOnSet (BioseqSetPtr bssp, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitFeaturesProc (bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitFeaturesInSet (BioseqSetPtr bssp, Pointer userdata, VisitFeaturesFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitFeaturesProc (bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitFeaturesInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitFeaturesOnSep (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitFeaturesOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitFeaturesOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitFeaturesInSep (SeqEntryPtr sep, Pointer userdata, VisitFeaturesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitFeaturesOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitFeaturesInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitAlignmentsOnDisc (Pointer segs, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4         index = 0;
  SeqAlignPtr  salp;

  for (salp = (SeqAlignPtr) segs; salp != NULL; salp = salp->next) {
    if (callback != NULL) {
      callback (salp, userdata);
    }
    index++;
    if (salp->segtype == SAS_DISC) {
      index += VisitAlignmentsOnDisc (salp->segs, userdata, callback);
    }
  }
  return index;
}

static Int4 VisitAlignmentsProc (SeqAnnotPtr annot, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4         index = 0;
  SeqAlignPtr  salp;
  SeqAnnotPtr  sap;

  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 2) continue;
    for (salp = (SeqAlignPtr) sap->data; salp != NULL; salp = salp->next) {
      if (callback != NULL) {
        callback (salp, userdata);
      }
      index++;
      if (salp->segtype == SAS_DISC) {
        index += VisitAlignmentsOnDisc (salp->segs, userdata, callback);
      }
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsOnSap (SeqAnnotPtr sap, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4         index = 0;
  SeqAlignPtr  salp;

  if (sap == NULL) return index;
  if (sap->type != 2) return index;
  for (salp = (SeqAlignPtr) sap->data; salp != NULL; salp = salp->next) {
    if (callback != NULL) {
      callback (salp, userdata);
    }
    index++;
    if (salp->segtype == SAS_DISC) {
      index += VisitAlignmentsOnDisc (salp->segs, userdata, callback);
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsOnBsp (BioseqPtr bsp, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitAlignmentsProc (bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitAlignmentsProc (bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsInSet (BioseqSetPtr bssp, Pointer userdata, VisitAlignmentsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitAlignmentsProc (bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitAlignmentsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsOnSep (SeqEntryPtr sep, Pointer userdata, VisitAlignmentsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitAlignmentsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitAlignmentsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitAlignmentsInSep (SeqEntryPtr sep, Pointer userdata, VisitAlignmentsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitAlignmentsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitAlignmentsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitGraphsProc (SeqAnnotPtr annot, Pointer userdata, VisitGraphsFunc callback)

{
  Int4         index = 0;
  SeqAnnotPtr  sap;
  SeqGraphPtr  sgp;

  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 3) continue;
    for (sgp = (SeqGraphPtr) sap->data; sgp != NULL; sgp = sgp->next) {
      if (callback != NULL) {
        callback (sgp, userdata);
      }
      index++;
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitGraphsOnSap (SeqAnnotPtr sap, Pointer userdata, VisitGraphsFunc callback)

{
  Int4         index = 0;
  SeqGraphPtr  sgp;

  if (sap == NULL) return index;
  if (sap->type != 3) return index;
  for (sgp = (SeqGraphPtr) sap->data; sgp != NULL; sgp = sgp->next) {
    if (callback != NULL) {
      callback (sgp, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitGraphsOnBsp (BioseqPtr bsp, Pointer userdata, VisitGraphsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitGraphsProc (bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitGraphsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitGraphsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitGraphsProc (bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitGraphsInSet (BioseqSetPtr bssp, Pointer userdata, VisitGraphsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitGraphsProc (bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitGraphsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitGraphsOnSep (SeqEntryPtr sep, Pointer userdata, VisitGraphsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitGraphsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitGraphsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitGraphsInSep (SeqEntryPtr sep, Pointer userdata, VisitGraphsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitGraphsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitGraphsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitAnnotsProc (SeqAnnotPtr annot, Pointer userdata, VisitAnnotsFunc callback)

{
  Int4         index = 0;
  SeqAnnotPtr  sap;

  for (sap = annot; sap != NULL; sap = sap->next) {
    if (callback != NULL) {
      callback (sap, userdata);
    }
    index++;
  }
  return index;
}

NLM_EXTERN Int4 VisitAnnotsOnBsp (BioseqPtr bsp, Pointer userdata, VisitAnnotsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitAnnotsProc (bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitAnnotsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitAnnotsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitAnnotsProc (bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitAnnotsInSet (BioseqSetPtr bssp, Pointer userdata, VisitAnnotsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitAnnotsProc (bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitAnnotsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitAnnotsOnSep (SeqEntryPtr sep, Pointer userdata, VisitAnnotsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitAnnotsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitAnnotsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitAnnotsInSep (SeqEntryPtr sep, Pointer userdata, VisitAnnotsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitAnnotsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitAnnotsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitAuthorsProc (AuthListPtr alp, Pointer userdata, VisitAuthorFunc callback)

{
  AuthorPtr    ap;
  Int4         index = 0;
  ValNodePtr   names;
  NameStdPtr   nsp;
  PersonIdPtr  pid;

  if (alp == NULL || alp->choice != 1) return index;

  for (names = alp->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (ap == NULL) continue;
    pid = ap->name;
    if (pid == NULL || pid->choice != 2) continue;
    nsp = pid->data;
    if (nsp == NULL) continue;
    if (callback != NULL) {
      callback (nsp, userdata);
    }
    index++;
  }

  return index;
}

NLM_EXTERN Int4 VisitAuthorsInPub (PubdescPtr pdp, Pointer userdata, VisitAuthorFunc callback)

{
  CitArtPtr   cap;
  CitBookPtr  cbp;
  CitGenPtr   cgp;
  CitPatPtr   cpp;
  CitSubPtr   csp;
  Int4        index = 0;
  ValNodePtr  vnp;

  if (pdp == NULL) return index;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_PMid || vnp->choice == PUB_Muid) continue;
    if (vnp->data.ptrvalue == NULL) continue;
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (cgp->authors, userdata, callback);
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (csp->authors, userdata, callback);
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (cap->authors, userdata, callback);
        if (cap->from == 2 || cap->from == 3) {
          cbp = (CitBookPtr) cap->fromptr;
          if (cbp != NULL) {
            index += VisitAuthorsProc (cbp->authors, userdata, callback);
          }
        }
        break;
      case PUB_Book :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (cbp->authors, userdata, callback);
        break;
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp->othertype == 2 && cbp->let_type == 3) {
          index += VisitAuthorsProc (cbp->authors, userdata, callback);
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        index += VisitAuthorsProc (cpp->authors, userdata, callback);
        index += VisitAuthorsProc (cpp->applicants, userdata, callback);
        index += VisitAuthorsProc (cpp->assignees, userdata, callback);
        break;
      default :
        break;
    }
  }

  return index;
}


static Int4 VisitPubdescsProc (SeqDescrPtr descr, SeqAnnotPtr annot, Pointer userdata, VisitPubdescsFunc callback)

{
  Int4         index = 0;
  PubdescPtr   pdp;
  SeqAnnotPtr  sap;
  SeqDescrPtr  sdp;
  SeqFeatPtr   sfp;

  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      if (pdp != NULL) {
        if (callback != NULL) {
          callback (pdp, userdata);
        }
        index++;
      }
    }
  }
  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 1) continue;
    for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
      if (sfp->data.choice == SEQFEAT_PUB) {
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        if (pdp != NULL) {
          if (callback != NULL) {
            callback (pdp, userdata);
          }
          index++;
        }
      }
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitPubdescsOnBsp (BioseqPtr bsp, Pointer userdata, VisitPubdescsFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitPubdescsProc (bsp->descr, bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitPubdescsOnSet (BioseqSetPtr bssp, Pointer userdata, VisitPubdescsFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitPubdescsProc (bssp->descr, bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitPubdescsInSet (BioseqSetPtr bssp, Pointer userdata, VisitPubdescsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitPubdescsProc (bssp->descr, bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitPubdescsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitPubdescsOnSep (SeqEntryPtr sep, Pointer userdata, VisitPubdescsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitPubdescsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitPubdescsOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitPubdescsInSep (SeqEntryPtr sep, Pointer userdata, VisitPubdescsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitPubdescsOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitPubdescsInSet (bssp, userdata, callback);
  }
  return index;
}


static Int4 VisitBioSourcesProc (SeqDescrPtr descr, SeqAnnotPtr annot, Pointer userdata, VisitBioSourcesFunc callback)

{
  BioSourcePtr  biop;
  Int4          index = 0;
  SeqAnnotPtr   sap;
  SeqDescrPtr   sdp;
  SeqFeatPtr    sfp;

  for (sdp = descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        if (callback != NULL) {
          callback (biop, userdata);
        }
        index++;
      }
    }
  }
  for (sap = annot; sap != NULL; sap = sap->next) {
    if (sap->type != 1) continue;
    for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
      if (sfp->data.choice == SEQFEAT_BIOSRC) {
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
        if (biop != NULL) {
          if (callback != NULL) {
            callback (biop, userdata);
          }
          index++;
        }
      }
    }
  }
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesOnBsp (BioseqPtr bsp, Pointer userdata, VisitBioSourcesFunc callback)

{
  Int4  index = 0;

  if (bsp == NULL) return index;
  index += VisitBioSourcesProc (bsp->descr, bsp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesOnSet (BioseqSetPtr bssp, Pointer userdata, VisitBioSourcesFunc callback)

{
  Int4  index = 0;

  if (bssp == NULL) return index;
  index += VisitBioSourcesProc (bssp->descr, bssp->annot, userdata, callback);
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesInSet (BioseqSetPtr bssp, Pointer userdata, VisitBioSourcesFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  index += VisitBioSourcesProc (bssp->descr, bssp->annot, userdata, callback);
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitBioSourcesInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesOnSep (SeqEntryPtr sep, Pointer userdata, VisitBioSourcesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitBioSourcesOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitBioSourcesOnSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitBioSourcesInSep (SeqEntryPtr sep, Pointer userdata, VisitBioSourcesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    index += VisitBioSourcesOnBsp (bsp, userdata, callback);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitBioSourcesInSet (bssp, userdata, callback);
  }
  return index;
}


NLM_EXTERN Int4 VisitBioseqsInSet (BioseqSetPtr bssp, Pointer userdata, VisitBioseqsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitBioseqsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitBioseqsInSep (SeqEntryPtr sep, Pointer userdata, VisitBioseqsFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (callback != NULL) {
      callback (bsp, userdata);
    }
    index++;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitBioseqsInSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitSequencesInSet (BioseqSetPtr bssp, Pointer userdata, Int2 filter, VisitSequencesFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  if (bssp->_class == BioseqseqSet_class_parts) {
    if (filter != VISIT_PARTS) return index;
    filter = VISIT_MAINS;
  }
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitSequencesInSep (tmp, userdata, filter, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitSequencesInSep (SeqEntryPtr sep, Pointer userdata, Int2 filter, VisitSequencesFunc callback)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (filter == VISIT_MAINS ||
        (filter == VISIT_NUCS && ISA_na (bsp->mol)) ||
        (filter == VISIT_PROTS && ISA_aa (bsp->mol))) {
      if (callback != NULL) {
        callback (bsp, userdata);
      }
      index++;
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitSequencesInSet (bssp, userdata, filter, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitSetsInSet (BioseqSetPtr bssp, Pointer userdata, VisitSetsFunc callback)

{
  Int4         index = 0;
  SeqEntryPtr  tmp;

  if (bssp == NULL) return index;
  if (callback != NULL) {
    callback (bssp, userdata);
  }
  index++;
  for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
    index += VisitSetsInSep (tmp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitSetsInSep (SeqEntryPtr sep, Pointer userdata, VisitSetsFunc callback)

{
  BioseqSetPtr  bssp;
  Int4          index = 0;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    index += VisitSetsInSet (bssp, userdata, callback);
  }
  return index;
}

NLM_EXTERN Int4 VisitElementsInSep (SeqEntryPtr sep, Pointer userdata, VisitElementsFunc callback)

{
  BioseqSetPtr  bssp;
  Int4          index = 0;
  SeqEntryPtr   tmp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return index;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return index;
    if (bssp->_class == 7 ||
        (bssp->_class >= 13 && bssp->_class <= 16) ||
        bssp->_class == BioseqseqSet_class_wgs_set ||
        bssp->_class == BioseqseqSet_class_gen_prod_set ||
        bssp->_class == BioseqseqSet_class_small_genome_set) {
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        index += VisitElementsInSep (tmp, userdata, callback);
      }
      return index;
    }
  }
  if (callback != NULL) {
    callback (sep, userdata);
  }
  index++;
  return index;
}

NLM_EXTERN Boolean IsPopPhyEtcSet (Uint1 _class)

{
  if (_class == BioseqseqSet_class_mut_set ||
      _class == BioseqseqSet_class_pop_set ||
      _class == BioseqseqSet_class_phy_set ||
      _class == BioseqseqSet_class_eco_set ||
      _class == BioseqseqSet_class_wgs_set ||
      _class == BioseqseqSet_class_small_genome_set) return TRUE;
  return FALSE;
}


NLM_EXTERN void CleanupStringsForOneDescriptor (SeqDescPtr sdp, SeqEntryPtr sep)
{
  Boolean stripSerial = FALSE;
  Boolean isEmblOrDdbj = FALSE;

  if (sdp == NULL) {
    return;
  }
  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtID);
  SeqEntryExplore (sep, (Pointer) &isEmblOrDdbj, CheckForEmblDdbjID);

  if (sdp->choice == Seq_descr_pub) {
    FlattenPubdesc (sdp->data.ptrvalue, NULL);
  }

  CleanupDescriptorStrings (sdp, stripSerial, TRUE, NULL, isEmblOrDdbj);
}


NLM_EXTERN void CleanupOneSeqFeat (SeqFeatPtr sfp)
{
  Boolean           isEmblOrDdbj = FALSE;
  Boolean           isJscan = FALSE;
  Boolean           stripSerial = TRUE;
  ValNodePtr        publist = NULL;
  SeqEntryPtr       sep;

  if (sfp->idx.entityID == 0) {
    return;
  }
  sep = GetTopSeqEntryForEntityID (sfp->idx.entityID);

  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtID);
  SeqEntryExplore (sep, (Pointer) &isEmblOrDdbj, CheckForEmblDdbjID);
  SeqEntryExplore (sep, (Pointer) &isJscan, CheckForJournalScanID);
  FlattenSfpCit (sfp, NULL);
  CleanUpSeqFeat (sfp, isEmblOrDdbj, isJscan, stripSerial, TRUE, &publist);

  if (publist != NULL) {
   ChangeCitsOnFeats (sfp, publist);
  }
  ValNodeFreeData (publist);
}
//LCOV_EXCL_STOP

NLM_EXTERN void RemoveFeatureLink (SeqFeatPtr sfp1, SeqFeatPtr sfp2)
{
  SeqFeatXrefPtr  xref, next, PNTR prevlink;
  ObjectIdPtr     oip;
  SeqFeatPtr      link_sfp;
  Char            buf [32];
  CharPtr         str = NULL;

  if (sfp1 == NULL) return;

  prevlink = (SeqFeatXrefPtr PNTR) &(sfp1->xref);
  xref = sfp1->xref;
  while (xref != NULL) {
    next = xref->next;
    link_sfp = NULL;

    if (xref->id.choice == 3) {
      oip = (ObjectIdPtr) xref->id.value.ptrvalue;
      if (oip != NULL) {
        if (StringDoesHaveText (oip->str)) {
          str = oip->str;
        } else {
          sprintf (buf, "%ld", (long) oip->id);
          str = buf;
        }
        link_sfp = SeqMgrGetFeatureByFeatID (sfp1->idx.entityID, NULL, str, NULL, NULL);
      }
    }
    if (link_sfp == sfp2) {
      *prevlink = xref->next;
      xref->next = NULL;
      MemFree (xref);
    } else {
      prevlink = (SeqFeatXrefPtr PNTR) &(xref->next);
    }

    xref = next;
  }
}


NLM_EXTERN void LinkTwoFeatures (SeqFeatPtr dst, SeqFeatPtr sfp)

{
  ChoicePtr       cp;
  ObjectIdPtr     oip;
  SeqFeatXrefPtr  xref, prev_xref, next_xref;
  SeqFeatPtr      old_match;

  if (dst == NULL || sfp == NULL) return;

  cp = &(dst->id);
  if (cp == NULL) return;
  if (cp->choice == 3) {
    /* don't create a duplicate xref, remove links to other features */
    xref = sfp->xref;
    prev_xref = NULL;
    while (xref != NULL) {
      next_xref = xref->next;
      if (xref->id.choice == 3 && xref->id.value.ptrvalue != NULL) {
        if (ObjectIdMatch (cp->value.ptrvalue, xref->id.value.ptrvalue)) {
          /* already have this xref */
          return;
        } else {
          old_match = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
          RemoveFeatureLink (sfp, old_match);
          RemoveFeatureLink (old_match, sfp);
        }
      } else {
        prev_xref = xref;
      }
      xref = next_xref;
    }

    oip = (ObjectIdPtr) cp->value.ptrvalue;
    if (oip != NULL) {
      oip = AsnIoMemCopy (oip, (AsnReadFunc) ObjectIdAsnRead,
                          (AsnWriteFunc) ObjectIdAsnWrite);
      if (oip != NULL) {
        xref = SeqFeatXrefNew ();
        if (xref != NULL) {
          xref->id.choice = 3;
          xref->id.value.ptrvalue = (Pointer) oip;
          xref->next = sfp->xref;
          sfp->xref = xref;
        }
      }
    }
  }
}

/* basic cleanup code from sqnutil3.c */

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

extern void ExtendSingleGeneOnMRNA (BioseqPtr bsp, Pointer userdata)

{
  MolInfoPtr        mip;
  SeqDescrPtr       sdp;
  Boolean           is_mrna = FALSE, is_master_seq = FALSE, has_nulls = FALSE;
  SeqFeatPtr        gene = NULL;
  SeqFeatPtr        sfp;
  SeqMgrFeatContext context;
  Int4              num_cds = 0;
  Int4              num_mrna = 0;
  SeqIdPtr          sip;
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
      /* skip this sequence if it has more than one coding region */
      if (num_cds > 1 && !is_master_seq) {
        return;
      }
    } else if (sfp->idx.subtype == FEATDEF_mRNA) {
      num_mrna++;
      /* skip this sequence if it has more than one mRNA */
      if (num_mrna > 1) return;
    }
  }

  if (gene != NULL && gene->location != NULL) {
    slp = gene->location;
    if (slp->choice != SEQLOC_INT) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        /* skip this sequence if it is multi-interval and EMBL or DDBJ */
        if (sip->choice == SEQID_EMBL || sip->choice == SEQID_DDBJ) return;
      }
    }
  }

  if (gene != NULL && BioseqFindFromSeqLoc (gene->location) == bsp) {
    CheckSeqLocForPartial (gene->location, &partial5, &partial3);
    has_nulls = LocationHasNullsBetween (gene->location);
    /* gene should cover entire length of sequence */
    slp = SeqLocIntNew (0, bsp->length - 1, SeqLocStrand (gene->location), SeqIdFindBest (bsp->id, 0));
    SetSeqLocPartial (slp, partial5, partial3);
    gene->location = SeqLocFree (gene->location);
    gene->location = slp;
    if (is_master_seq) {
      MergeFeatureIntervalsToParts (gene, has_nulls);
    }
  }
}

//LCOV_EXCL_START
static DbtagPtr DbtagParse (
  CharPtr str
)

{
  Boolean      all_digits = TRUE;
  Char         ch;
  DbtagPtr     dbt;
  long         num;
  Int2         num_digits = 0;
  ObjectIdPtr  oip;
  CharPtr      ptr;
  CharPtr      tmp;

  if (StringHasNoText (str)) return NULL;
  ptr = StringChr (str, ':');
  if (ptr == NULL) return NULL;

  dbt = DbtagNew ();
  oip = ObjectIdNew ();
  if (dbt == NULL || oip == NULL) return NULL;

  if (ptr != NULL) {
    *ptr = '\0';
    ptr++;
  }

  dbt->db = StringSave (str);
  dbt->tag = oip;

  tmp = ptr;
  ch = *tmp;
  while (ch != '\0') {
    if (IS_DIGIT (ch)) {
      num_digits++;
    } else {
      all_digits = FALSE;
    }
    tmp++;
    ch = *tmp;
  }

  if (all_digits && *ptr != '0') {
    if (num_digits < 10 || (num_digits == 10 && StringCmp (ptr, "2147483647") <= 0)) {
      sscanf (ptr, "%ld", &num);
      oip->id = (Int4) num;
      return dbt;
    }
  }

  oip->str = StringSave (ptr);

  return dbt;
}
//LCOV_EXCL_STOP

static void GetNomenclatureUOP (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr         oip;
  UserObjectPtr PNTR  uopp;

  if (uop == NULL || userdata == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "OfficialNomenclature") != 0) return;
  uopp = (UserObjectPtr PNTR) userdata;
  *uopp = uop;
}


//LCOV_EXCL_START
NLM_EXTERN void ModernizeGeneFields (
  SeqFeatPtr sfp
)

{
  GeneNomenclaturePtr  gnp;
  GeneRefPtr           grp;
  ObjectIdPtr          oip;
  CharPtr              str;
  CharPtr              symbol = NULL, name = NULL, source = NULL;
  Uint2                status = 0;
  UserFieldPtr         ufp;
  UserObjectPtr        uop = NULL;
  UserObjectPtr        curr, next;
  UserObjectPtr PNTR   prev;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_GENE) return;

  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return;

  if (grp->formal_name != NULL) return;

  if (sfp->ext == NULL) return;
  VisitUserObjectsInUop (sfp->ext, (Pointer) &uop, GetNomenclatureUOP);
  if (uop == NULL) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL || oip->str == NULL) continue;
    if (StringICmp (oip->str, "Symbol") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (str != NULL) {
          symbol = str;
        }
      }
    } else if (StringICmp (oip->str, "Name") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (str != NULL) {
          name = str;
        }
      }
    } else if (StringICmp (oip->str, "DataSource") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (str != NULL) {
          source = str;
        }
      }
    } else if (StringICmp (oip->str, "Status") == 0) {
      if (ufp->choice == 1) {
        str = (CharPtr) ufp->data.ptrvalue;
        if (str != NULL) {
          if (StringICmp (str, "Official") == 0) {
            status = 1;
          } else if (StringICmp (str, "Interim") == 0) {
            status = 2;
          }
        }
      }
    }
  }
  if (symbol == NULL && name == NULL && source == NULL && status == 0) return;

  gnp = GeneNomenclatureNew ();
  if (gnp == NULL) return;

  gnp->status = status;
  gnp->symbol = StringSaveNoNull (symbol);
  gnp->name = StringSaveNoNull (name);
  gnp->source = DbtagParse (source);

  grp->formal_name = gnp;

  prev = (UserObjectPtr PNTR) &(sfp->ext);
  curr = sfp->ext;
  while (curr != NULL) {
    next = curr->next;
    if (uop == curr) {
      *(prev) = curr->next;
      curr->next = NULL;
      UserObjectFree (curr);
    } else {
      prev = (UserObjectPtr PNTR) &(curr->next);
    }
    curr = next;
  }
}
//LCOV_EXCL_STOP


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

static ValNodePtr ParsePCRColonString (
  CharPtr strs
)

{
  ValNodePtr  head = NULL;
  size_t      len;
  CharPtr     ptr, str, tmp;

  if (StringHasNoText (strs)) return NULL;

  tmp = StringSave (strs);
  str = tmp;
  len = StringLen (str);
  if (len > 1 && StringChr (str, ':') != NULL) {
    while (StringDoesHaveText (str)) {
      ptr = StringChr (str, ':');
      if (ptr != NULL) {
        *ptr = '\0';
        ptr++;
      }
      TrimSpacesAroundString (str);
      ValNodeCopyStr (&head, 0, str);
      str = ptr;
    }
  } else {
    ValNodeCopyStr (&head, 0, str);
  }

  MemFree (tmp);
  return head;
}

//LCOV_EXCL_START
static CharPtr FusePrimerNames(
  CharPtr first,
  CharPtr second
)

{
  size_t   len;
  CharPtr  str;

  if (first == NULL) return second;
  if (second == NULL) return first;

  len = StringLen (first) + StringLen (second) + 5;
  str = MemNew (len);
  if (str == NULL) return NULL;

  StringCpy (str, first);
  StringCat (str, ":");
  StringCat (str, second);

  return str;
}

static PCRPrimerPtr ModernizePCRPrimerHalf (
  CharPtr seq,
  CharPtr name
)

{
  CharPtr       curr_name = NULL, curr_seq = NULL, fused_name;
  PCRPrimerPtr  curr_primer = NULL, last_primer = NULL, primer_set = NULL;
  ValNodePtr    name_list, seq_list, name_vnp, seq_vnp;

  seq_list = ParsePCRColonString (seq);
  name_list = ParsePCRColonString (name);

  seq_vnp = seq_list;
  name_vnp = name_list;

  while (seq_vnp != NULL /* || name_vnp != NULL */) {
    if (seq_vnp != NULL) {
      curr_seq = (CharPtr) seq_vnp->data.ptrvalue;
      seq_vnp = seq_vnp->next;
    }
    if (name_vnp != NULL) {
      curr_name = (CharPtr) name_vnp->data.ptrvalue;
      name_vnp = name_vnp->next;
    } else {
      curr_name = NULL;
    }

    curr_primer = (PCRPrimerPtr) MemNew (sizeof (PCRPrimer));
    if (curr_primer != NULL) {
      curr_primer->seq = StringSaveNoNull (curr_seq);
      curr_primer->name = StringSaveNoNull (curr_name);

      if (primer_set == NULL) {
        primer_set = curr_primer;
      }
      if (last_primer != NULL) {
        last_primer->next = curr_primer;
      }
      last_primer = curr_primer;
    }
  }

  while (name_vnp != NULL && last_primer != NULL) {
    curr_name = (CharPtr) name_vnp->data.ptrvalue;
    fused_name = FusePrimerNames (last_primer->name, curr_name);
    MemFree (last_primer->name);
    last_primer->name = StringSaveNoNull (fused_name);
    name_vnp = name_vnp->next;
  }

  while (name_vnp != NULL && last_primer == NULL) {
    curr_name = (CharPtr) name_vnp->data.ptrvalue;
    curr_primer = (PCRPrimerPtr) MemNew (sizeof (PCRPrimer));
    if (curr_primer != NULL) {
      curr_primer->name = StringSaveNoNull (curr_name);

      if (primer_set == NULL) {
        primer_set = curr_primer;
      }
      if (last_primer != NULL) {
        last_primer->next = curr_primer;
      }
      last_primer = curr_primer;
    }
    name_vnp = name_vnp->next;
  }

  ValNodeFreeData (seq_list);
  ValNodeFreeData (name_list);

  return primer_set;
}

NLM_EXTERN void ModernizePCRPrimers (
  BioSourcePtr biop
)

{
  PCRReactionSetPtr  curr_reaction, last_reaction = NULL, reaction_set = NULL;
  PCRPrimerPtr       forward, reverse;
  PcrSetPtr          psp;
  ValNodePtr         pset, vnp;
  SubSourcePtr       nextssp;
  SubSourcePtr PNTR  prevssp;
  SubSourcePtr       ssp;
  Boolean            unlink;

  if (biop == NULL) return;
  /* if (biop->pcr_primers != NULL) return; */

  pset = ParsePCRSet (biop);
  if (pset == NULL) return;

  for (vnp = pset; vnp != NULL; vnp = vnp->next) {
    psp = (PcrSetPtr) vnp->data.ptrvalue;
    if (psp == NULL) continue;

    forward = ModernizePCRPrimerHalf (psp->fwd_seq, psp->fwd_name);
    reverse = ModernizePCRPrimerHalf (psp->rev_seq, psp->rev_name);

    if (forward != NULL || reverse != NULL) {

      curr_reaction = (PCRReactionSetPtr) MemNew (sizeof (PCRReactionSet));
      if (curr_reaction != NULL) {
        curr_reaction->forward = forward;
        curr_reaction->reverse = reverse;

        if (reaction_set == NULL) {
          reaction_set = curr_reaction;
        }
        if (last_reaction != NULL) {
          last_reaction->next = curr_reaction;
        }
        last_reaction = curr_reaction;
      }
    }
  }

  FreePCRSet (pset);

  if (reaction_set != NULL) {
    if (last_reaction != NULL) {
      /* merge with existing structured pcr_primers */
      last_reaction->next = biop->pcr_primers;
    }
    biop->pcr_primers = reaction_set;

    ssp = biop->subtype;
    prevssp = (SubSourcePtr PNTR) &(biop->subtype);
    while (ssp != NULL) {
      nextssp = ssp->next;
      unlink= FALSE;

      if (ssp->subtype == SUBSRC_fwd_primer_seq ||
          ssp->subtype == SUBSRC_rev_primer_seq ||
          ssp->subtype == SUBSRC_fwd_primer_name ||
          ssp->subtype == SUBSRC_rev_primer_name) {
        unlink = TRUE;
      }

      if (unlink) {
        *prevssp = ssp->next;
        ssp->next = NULL;
        SubSourceFree (ssp);
      } else {
        prevssp = (SubSourcePtr PNTR) &(ssp->next);
      }
      ssp = nextssp;
    }
  }
}
//LCOV_EXCL_STOP

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




