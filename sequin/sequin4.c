/*   sequin4.c
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
* File Name:  sequin4.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   6/28/96
*
* $Revision: 6.93 $
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

#include "sequin.h"
#include <ncbilang.h>
#include <seqport.h>
#include <gather.h>
#include <objall.h>
#include <objcode.h>
#include <utilpub.h>
#include <vibrant.h>
#include <document.h>
#include <toasn3.h>
#include <asn2ffp.h>
#include <salfiles.h>
#include <salsap.h>
#include <saledit.h>
#include <subutil.h>
#include <pobutil.h>
#include <tfuns.h>
#include <edutil.h>
#include <biosrc.h>
#include <seqgrphx.h>
#include <seqmtrx.h>
#include <bspview.h>
#include <vsmpriv.h>
#include <explore.h>

#define REGISTER_UPDATESEGSET ObjMgrProcLoadEx (OMPROC_FILTER,"Update Segmented Set","UpdateSegSet",0,0,0,0,NULL,UpdateSegSet,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_ADJUSTMULTISEGSEQ ObjMgrProcLoadEx (OMPROC_FILTER,"Adjust SegSeq Length","AdjustSegLength",0,0,0,0,NULL,AdjustSegSeqLength,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_UNDOSEGSET ObjMgrProcLoadEx (OMPROC_FILTER,"Undo Segmented Set","UndoSegSet",0,0,0,0,NULL,UndoSegSet,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_REPACKAGE_PARTS ObjMgrProcLoadEx (OMPROC_FILTER,"Repackage Segmented Parts","RepackageParts",0,0,0,0,NULL,PackagePartsInPartsSet,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_NORMALIZE_NUCPROT ObjMgrProcLoadEx (OMPROC_FILTER,"Normalize Nuc-Prot","NormalizeNucProts",0,0,0,0,NULL,NormalizeNucProts,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_REMOVE_EXTRANEOUS ObjMgrProcLoadEx (OMPROC_FILTER,"Remove Extraneous Sets","RemoveExtraneousSets",0,0,0,0,NULL,RemoveExtraneousSets,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_REMOVE_MESSEDUP ObjMgrProcLoadEx (OMPROC_FILTER,"Repair Messed Up Sets","RepairMessedUpSets",0,0,0,0,NULL,RepairMessedUpRecord,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_MAKESEQALIGN ObjMgrProcLoadEx (OMPROC_FILTER,"Make SeqAlign","CreateSeqAlign",0,0,0,0,NULL,GenerateSeqAlignFromSeqEntry,PROC_PRIORITY_DEFAULT, "Alignment")

#define REGISTER_MAKESEQALIGNP ObjMgrProcLoadEx (OMPROC_FILTER,"Make Protein SeqAlign","CreateSeqAlignProt",0,0,0,0,NULL,GenerateSeqAlignFromSeqEntryProt,PROC_PRIORITY_DEFAULT, "Alignment")

#define REGISTER_NORMSEQALIGN ObjMgrProcLoadEx (OMPROC_FILTER,"Validate SeqAlign","ValidateSeqAlign",0,0,0,0,NULL,ValidateSeqAlignFromSeqEntry,PROC_PRIORITY_DEFAULT, "Alignment")

#define REGISTER_GROUP_EXPLODE ObjMgrProcLoadEx (OMPROC_FILTER, "Explode a group", "GroupExplode", OBJ_SEQFEAT, 0, OBJ_SEQFEAT, 0, NULL, GroupExplodeFunc, PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_INTERVAL_COMBINE ObjMgrProcLoadEx (OMPROC_FILTER, "Combine feature intervals", "IntervalCombine", OBJ_SEQFEAT, 0, OBJ_SEQFEAT, 0, NULL, IntervalCombineFunc, PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_MRNA_FROM_CDS ObjMgrProcLoadEx (OMPROC_FILTER, "mRNA from CDS", "mRNAfromCDS", OBJ_SEQFEAT, FEATDEF_CDS, OBJ_SEQFEAT, FEATDEF_CDS, NULL, MRnaFromCdsFunc, PROC_PRIORITY_DEFAULT, "Utilities")

#define REGISTER_SPLIT_BIOSEQ ObjMgrProcLoadEx (OMPROC_FILTER,"Split Bioseq Into Segments","SplitBioseqIntoSegments",0,0,0,0,NULL,SplitIntoSegmentedBioseq,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_MAP_TO_PROT ObjMgrProcLoadEx (OMPROC_FILTER,"Map to Prot","MapToProt",OBJ_SEQFEAT,0,OBJ_SEQFEAT,0,NULL,MapToProtFunc,PROC_PRIORITY_DEFAULT, "Utilities")

#define REGISTER_MAP_TO_NUC ObjMgrProcLoadEx (OMPROC_FILTER,"Map to Nuc","MapToNuc",OBJ_SEQFEAT,0,OBJ_SEQFEAT,0,NULL,MapToNucFunc,PROC_PRIORITY_DEFAULT, "Utilities")

#define REGISTER_BIOSEQ_ORF ObjMgrProcLoadEx (OMPROC_FILTER, "ORF Finder", "OrfFinder", OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0, NULL, OrfFindFunc, PROC_PRIORITY_DEFAULT, "Utilities")

#define REGISTER_BIOSEQ_REVCOMP_WITHFEAT ObjMgrProcLoadEx (OMPROC_FILTER, "Bioseq and Features RevComp", "BioseqFeatsRevComp", OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0, NULL, RevCompFuncFeat, PROC_PRIORITY_DEFAULT, "Utilities")
#define REGISTER_BIOSEQ_REVCOMP_NOTFEAT ObjMgrProcLoadEx (OMPROC_FILTER, "Bioseq only RevComp", "BioseqOnlyRevComp", OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0, NULL, RevCompFunc, PROC_PRIORITY_DEFAULT, "Utilities")
#define REGISTER_BIOSEQ_REVERSE ObjMgrProcLoadEx (OMPROC_FILTER, "Bioseq Reverse", "BioseqReverse", OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0, NULL, RevFunc, PROC_PRIORITY_DEFAULT, "Utilities")
#define REGISTER_BIOSEQ_COMPLEMENT ObjMgrProcLoadEx (OMPROC_FILTER, "Bioseq Complement", "BioseqComplement", OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0, NULL, CompFunc, PROC_PRIORITY_DEFAULT, "Utilities")

#define REGISTER_FIXUP_RBS ObjMgrProcLoadEx (OMPROC_FILTER,"Fixup RBS","FixupRBS",0,0,0,0,NULL,FixupRBS,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_BSP_INDEX ObjMgrProcLoadEx (OMPROC_FILTER,"Bioseq Index","MakeBioseqIndex",0,0,0,0,NULL,DoBioseqIndexing,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_TRIM_GENES ObjMgrProcLoadEx (OMPROC_FILTER,"Trim Genes","TrimGenes",0,0,0,0,NULL,TrimGenes,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_OPENALED ObjMgrProcLoadEx (OMPROC_FILTER,"Open Align Editor 1","Open protein BelVu", 0,0,0,0,NULL,LaunchAlignEditorFromDesktop, PROC_PRIORITY_DEFAULT, "Alignment")

#define REGISTER_OPENALED2 ObjMgrProcLoadEx (OMPROC_FILTER,"Open Align Editor 2","Open protein PHYLIP", 0,0,0,0,NULL,LaunchAlignEditorFromDesktop2, PROC_PRIORITY_DEFAULT, "Alignment")

#define REGISTER_OPENALED3 ObjMgrProcLoadEx (OMPROC_FILTER,"Open Align Editor 3","Open protein FASTA+GAP", 0,0,0,0,NULL,LaunchAlignEditorFromDesktop3, PROC_PRIORITY_DEFAULT, "Alignment")

#define REGISTER_OPENALED4 ObjMgrProcLoadEx (OMPROC_FILTER,"Open Align Editor 4","Open protein FASTA", 0,0,0,0,NULL,LaunchAlignEditorFromDesktop4, PROC_PRIORITY_DEFAULT, "Alignment")

#define REGISTER_MAKEEXONINTRON ObjMgrProcLoadEx (OMPROC_FILTER,"Make Exons and Introns","MakeExonIntron",OBJ_SEQFEAT,0,OBJ_SEQFEAT,0,NULL,MakeExonIntron,PROC_PRIORITY_DEFAULT, "Misc")

#define REGISTER_PROT_IDS_TO_GENE_SYN ObjMgrProcLoadEx (OMPROC_FILTER,"Protein SeqID to Gene Synonym","ProtLocalIDtoGeneSyn",0,0,0,0,NULL,ProtLocalIDtoGeneSyn,PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_DESKTOP_REPORT ObjMgrProcLoadEx (OMPROC_FILTER, "Desktop Report", "DesktopReport", 0, 0, 0, 0, NULL, DesktopReportFunc, PROC_PRIORITY_DEFAULT, "Indexer")

#define REGISTER_SEQUIN_PROT_TITLES ObjMgrProcLoadEx (OMPROC_FILTER,"Sequin Style Protein Titles","SequinStyleProteinTitles",0,0,0,0,NULL,MakeSequinProteinTitles,PROC_PRIORITY_DEFAULT, "Misc")
static Int2 LIBCALLBACK MakeSequinProteinTitles (Pointer data);

#define REGISTER_SEQUIN_NUC_TITLES ObjMgrProcLoadEx (OMPROC_FILTER,"Sequin Style Nucleotide Titles","SequinStyleNucleotideTitles",0,0,0,0,NULL,MakeSequinNucleotideTitles,PROC_PRIORITY_DEFAULT, "Misc")
static Int2 LIBCALLBACK MakeSequinNucleotideTitles (Pointer data);

#define REGISTER_SEQUIN_FEAT_TABLE ObjMgrProcLoadEx (OMPROC_FILTER,"Sequin Style Feature Table","SequinStyleFeatureTable",0,0,0,0,NULL,MakeSequinFeatureTable,PROC_PRIORITY_DEFAULT, "Misc")
static Int2 LIBCALLBACK MakeSequinFeatureTable (Pointer data);

#define REGISTER_SEQUIN_DUMP_CONTIG ObjMgrProcLoadEx (OMPROC_FILTER,"Make Contig Build Table","MakeContigBuildTable",0,0,0,0,NULL,MakeContigBuildTable,PROC_PRIORITY_DEFAULT, "Misc")
static Int2 LIBCALLBACK MakeContigBuildTable (Pointer data);

#define REGISTER_SEQUIN_GI_TO_ACCN ObjMgrProcLoadEx (OMPROC_FILTER,"Convert Align GI to Accession","ConvertAlignGisToAccn",0,0,0,0,NULL,AlignGiToAccnProc,PROC_PRIORITY_DEFAULT, "Misc")
static Int2 LIBCALLBACK AlignGiToAccnProc (Pointer data);

#define REGISTER_SEQUIN_CACHE_ACCN ObjMgrProcLoadEx (OMPROC_FILTER,"Cache Accessions to Disk","CacheAccnsToDisk",0,0,0,0,NULL,CacheAccnsToDisk,PROC_PRIORITY_DEFAULT, "Misc")
static Int2 LIBCALLBACK CacheAccnsToDisk (Pointer data);

#define REGISTER_REFGENEUSER_DESC_EDIT ObjMgrProcLoad(OMPROC_EDIT,"Edit RefGene UserTrack Desc","RefGene Tracking",OBJ_SEQDESC,Seq_descr_user,OBJ_SEQDESC,Seq_descr_user,NULL,RefGeneUserGenFunc,PROC_PRIORITY_DEFAULT)
extern Int2 LIBCALLBACK RefGeneUserGenFunc (Pointer data);

static void AddBspToSegSet (BioseqPtr segseq, BioseqPtr bsp)

{
  SeqIdPtr   sip;
  SeqLocPtr  slp;

  if (segseq == NULL) return;
  slp = ValNodeNew ((ValNodePtr) segseq->seq_ext);
  if (slp == NULL) return;
  if (segseq->seq_ext == NULL) {
    segseq->seq_ext = (Pointer) slp;
  }
  if (bsp != NULL && bsp->length > 0) {
    segseq->length += bsp->length;
    slp->choice = SEQLOC_WHOLE;
    sip = SeqIdFindBest (bsp->id, 0);
    slp->data.ptrvalue = (Pointer) SeqIdStripLocus (SeqIdDup (sip));
  } else {
    slp->choice = SEQLOC_NULL;
  }
}

static void RemoveMolInfoDescriptors (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     nextsdp;
  Pointer PNTR   prevsdp;
  ValNodePtr     sdp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_molinfo) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

static void MoveSegSetMolInfo (BioseqPtr segseq, BioseqSetPtr parts, BioseqSetPtr segset)

{
  MolInfoPtr   first;
  MolInfoPtr   mip;
  SeqEntryPtr  partssep;
  ValNodePtr   sdp;
  SeqEntryPtr  sep;
  SeqEntryPtr  tmp;

  if (segseq == NULL || parts == NULL || parts->seq_set == NULL || segset == NULL) return;
  sep = SeqMgrGetSeqEntryForData (segset);
  if (sep == NULL) return;
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
  if (sdp != NULL) return;
  first = NULL;
  partssep = SeqMgrGetSeqEntryForData (parts);
  if (partssep != NULL) {
    sdp = SeqEntryGetSeqDescr (partssep, Seq_descr_molinfo, NULL);
    if (sdp != NULL) {
      first = (MolInfoPtr) sdp->data.ptrvalue;
    }
  }
  for (tmp = parts->seq_set; tmp != NULL; tmp = tmp->next) {
    sdp = SeqEntryGetSeqDescr (tmp, Seq_descr_molinfo, NULL);
    if (sdp != NULL) {
      if (first == NULL) {
        first = (MolInfoPtr) sdp->data.ptrvalue;
      } else {
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (first != NULL && mip != NULL) {
          if (mip->biomol != first->biomol) return;
        }
      }
    }
  }
  if (first == NULL) return;
  mip = MolInfoNew ();
  if (mip == NULL) return;
  sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
  if (sdp == NULL) return;
  sdp->data.ptrvalue = (Pointer) mip;
  mip->biomol = first->biomol;
  mip->tech = first->tech;
  mip->completeness = first->completeness;
  mip->techexp = StringSaveNoNull (first->techexp);
  if (partssep != NULL) {
    SeqEntryExplore (partssep, NULL, RemoveMolInfoDescriptors);
  }
  for (tmp = parts->seq_set; tmp != NULL; tmp = tmp->next) {
    SeqEntryExplore (tmp, NULL, RemoveMolInfoDescriptors);
  }
}

static void DoUpdateSegSet (BioseqPtr segseq, BioseqSetPtr parts)

{
  MsgAnswer    ans;
  BioseqPtr    bsp;
  Boolean      notFirst;
  Boolean      nullsBetween;
  SeqFeatPtr   sfp;
  SeqFeatPtr   sfpnext;
  SeqLocPtr    slp;
  SeqEntryPtr  tmp;

  if (segseq == NULL || parts == NULL || parts->seq_set == NULL) return;
  tmp = parts->seq_set;
  notFirst = FALSE;
  nullsBetween = FALSE;
  segseq->length = 0;
  switch (segseq->seq_ext_type) {
    case 1:    /* seg-ext */
      slp = (ValNodePtr) segseq->seq_ext;
      while (slp != NULL) {
        if (slp->choice == SEQLOC_NULL) {
          nullsBetween = TRUE;
        }
        slp = slp->next;
      }
      SeqLocSetFree ((ValNodePtr) segseq->seq_ext);
      break;
    case 2:   /* reference */
      SeqLocFree ((ValNodePtr) segseq->seq_ext);
      break;
    case 3:   /* map */
      sfp = (SeqFeatPtr) segseq->seq_ext;
      while (sfp != NULL) {
        sfpnext = sfp->next;
        SeqFeatFree (sfp);
        sfp = sfpnext;
      }
      break;
    default:
      break;
  }
  segseq->seq_ext = NULL;
  if (nullsBetween) {
    ans = Message (MSG_YN, "Intersperse intervals with gaps (currently has gaps)?");
  } else {
    ans = Message (MSG_YN, "Intersperse intervals with gaps (currently no gaps)?");
  }
  nullsBetween = (Boolean) (ans == ANS_YES);
  while (tmp != NULL) {
    if (nullsBetween && notFirst) {
      AddBspToSegSet (segseq, NULL);
    }
    bsp = (BioseqPtr) tmp->data.ptrvalue;
    if (bsp != NULL) {
      AddBspToSegSet (segseq, bsp);
    }
    notFirst = TRUE;
    tmp = tmp->next;
  }
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

static Int2 DoOneSegFixup (SeqEntryPtr sep)

{
  BioseqSetPtr    bssp;
  Uint1           choice;
  Int2            count;
  SeqEntryPtr     insert;
  SeqEntryPtr     next;
  ObjMgrDataPtr   omdptop;
  ObjMgrData      omdata;
  Uint2           parenttype;
  Pointer         parentptr;
  UpdateSegStruc  uss;
  SeqEntryPtr     tmp;

  if (sep == NULL) return 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      choice = 0;
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        if (choice == 0) {
          choice = tmp->choice;
        } else if (choice != tmp->choice) {
          return 0;
        }
      }
      count = 0;
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        count += DoOneSegFixup (sep);
      }
      return count;
    }
  }
  if (IS_Bioseq (sep) && sep->next != NULL) {
    count = 0;
    SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
    GetSeqEntryParent (sep, &parentptr, &parenttype);
    insert = sep->next;
    sep->next = NULL;
    while (insert != NULL) {
      next = insert->next;
      insert->next = NULL;
      AddSeqEntryToSeqEntry (sep, insert, FALSE);
      count++;
      insert = next;
    }
    SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
    uss.segseq = NULL;
    uss.parts = NULL;
    uss.segset = NULL;
    UpdateSegExplore (sep, (Pointer) &uss, FindSegSetComponentsCallback);
    if (uss.segseq != NULL && uss.parts != NULL && uss.segset != NULL) {
      DoUpdateSegSet (uss.segseq, uss.parts);
      MoveSegSetMolInfo (uss.segseq, uss.parts, uss.segset);
    }
    return count;
  }
  uss.segseq = NULL;
  uss.parts = NULL;
  uss.segset = NULL;
  UpdateSegExplore (sep, (Pointer) &uss, FindSegSetComponentsCallback);
  if (uss.segseq != NULL && uss.parts != NULL && uss.segset != NULL) {
    DoUpdateSegSet (uss.segseq, uss.parts);
    MoveSegSetMolInfo (uss.segseq, uss.parts, uss.segset);
    return 1;
  }
  return 0;
}

static Int2 LIBCALLBACK UpdateSegSet (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  DoOneSegFixup (sep);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static void DoAdjustSegSeqLength (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr   bsp;
  ValNode     head;
  Int4        len;
  Int4        len2;
  ValNodePtr  vnp;

  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL || bsp->repr != Seq_repr_seg) return;
  head.choice = SEQLOC_MIX;
  head.data.ptrvalue = bsp->seq_ext;
  head.next = NULL;
  vnp = NULL;
  len = 0;
  while ((vnp = SeqLocFindNext (&head, vnp)) != NULL) {
    len2 = SeqLocLen (vnp);
    if (len2 > 0) {
      len += len2;
    }
  }
  bsp->length = len;
}

static Int2 LIBCALLBACK AdjustSegSeqLength (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  SeqEntryExplore (sep, NULL, DoAdjustSegSeqLength);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static void DoSegSetUndo (BioseqPtr segseq, BioseqSetPtr parts, BioseqSetPtr segset)

{
  SeqAnnotPtr  annot;
  ValNodePtr   descr;
  SeqAnnotPtr  sap;
  SeqEntryPtr  tmp;

  if (segseq == NULL || parts == NULL || parts->seq_set == NULL || segset == NULL) return;
  segset->_class = 7;
  parts->_class = 14;
  annot = segseq->annot;
  segseq->annot = NULL;
  if (segset->annot == NULL) {
    segset->annot = annot;
    annot = NULL;
  } else {
    sap = segset->annot;
    while (sap->next != NULL) {
      sap = sap->next;
    }
    sap->next = annot;
  }
  descr = segseq->descr;
  segseq->descr = NULL;
  ValNodeLink (&(segset->descr), descr);
  tmp = segset->seq_set;
  if (tmp != NULL && IS_Bioseq (tmp)) {
    segset->seq_set = tmp->next;
    tmp->next = NULL;
    SeqEntryFree (tmp);
  }
}

static Int2 DoOneSegUndo (SeqEntryPtr sep)

{
  BioseqSetPtr    bssp;
  Int2            count = 0;
  UpdateSegStruc  uss;

  if (sep == NULL) return 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        count += DoOneSegUndo (sep);
      }
      return count;
    }
  }
  uss.segseq = NULL;
  uss.parts = NULL;
  uss.segset = NULL;
  UpdateSegExplore (sep, (Pointer) &uss, FindSegSetComponentsCallback);
  if (uss.segseq != NULL && uss.parts != NULL && uss.segset != NULL) {
    DoSegSetUndo (uss.segseq, uss.parts, uss.segset);
    return 1;
  }
  return 0;
}

static SeqEntryPtr BioseqExtract (SeqEntryPtr top, BioseqPtr bsp)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   next;
  SeqEntryPtr   PNTR prev;
  SeqEntryPtr   rsult = NULL;
  SeqEntryPtr   sep;

  if (top == NULL || bsp == NULL) return NULL;
  if (! IS_Bioseq_set (top)) return NULL;
  bssp = (BioseqSetPtr) top->data.ptrvalue;
  if (bssp == NULL) return NULL;
  prev = &(bssp->seq_set);
  sep = bssp->seq_set;
  while (sep != NULL) {
    next = sep->next;
    if (IS_Bioseq_set (sep)) {
      rsult = BioseqExtract (sep, bsp);
      if (rsult != NULL) {
        return rsult;
      }
      prev = &(sep->next);
      sep = next;
    } else if (IS_Bioseq (sep) && sep->data.ptrvalue == (Pointer) bsp) {
      *(prev) = next;
      sep->next = NULL;
      return sep;
    } else {
      prev = &(sep->next);
      sep = next;
    }
  }
  return NULL;
}

static void DoRepairPartsSet (SeqEntryPtr sep)

{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  ValNode         head;
  BioseqSetPtr    parts;
  BioseqPtr       segseq;
  SeqLocPtr       slp;
  SeqEntryPtr     tmp;
  UpdateSegStruc  uss;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) {
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        DoRepairPartsSet (tmp);
      }
      return;
    }
  }
  uss.segseq = NULL;
  uss.parts = NULL;
  uss.segset = NULL;
  UpdateSegExplore (sep, (Pointer) &uss, FindSegSetComponentsCallback);
  if (uss.segseq == NULL || uss.parts == NULL || uss.segset == NULL) return;
  segseq = uss.segseq;
  parts = uss.parts;
  if (segseq->repr != Seq_repr_seg ||
      segseq->seq_ext_type != 1 ||
      segseq->seq_ext == NULL) return;
  head.choice = SEQLOC_MIX;
  head.data.ptrvalue = segseq->seq_ext;
  head.next = NULL;
  slp = NULL;
  while ((slp = SeqLocFindNext (&head, slp)) != NULL) {
    bsp = BioseqFind (SeqLocId (slp));
    if (bsp != NULL) {
      tmp = BioseqExtract (sep, bsp);
      if (tmp != NULL) {
        ValNodeLink (&parts->seq_set, tmp);
      }
    }
  }
}

static Int2 LIBCALLBACK PackagePartsInPartsSet (Pointer data)

{
  OMProcControlPtr  ompcp;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             parenttype;
  Pointer           parentptr;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);
  DoRepairPartsSet (sep);
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static void RemoveDupGenBankSets (SeqEntryPtr sep)

{
  SeqAnnotPtr   annot;
  BioseqSetPtr  bssp;
  ValNodePtr    descr;
  SeqAnnotPtr   sap;
  SeqEntryPtr   tmp;
  BioseqSetPtr  tmpbssp;

  if (sep == NULL || ! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL || bssp->_class != 7) return;
  tmp = bssp->seq_set;
  if (tmp == NULL || ! IS_Bioseq_set (tmp) || tmp->next != NULL) return;
  tmpbssp = (BioseqSetPtr) tmp->data.ptrvalue;
  if (tmpbssp == NULL || tmpbssp->_class != 7 || tmpbssp->seq_set == NULL) return;
  annot = tmpbssp->annot;
  tmpbssp->annot = NULL;
  if (bssp->annot == NULL) {
    tmpbssp->annot = annot;
    annot = NULL;
  } else {
    sap = bssp->annot;
    while (sap->next != NULL) {
      sap = sap->next;
    }
    sap->next = annot;
  }
  descr = tmpbssp->descr;
  tmpbssp->descr = NULL;
  ValNodeLink (&(bssp->descr), descr);
  bssp->seq_set = tmpbssp->seq_set;
  tmpbssp->seq_set = NULL;
  SeqEntryFree (tmp);
}

static Int2 LIBCALLBACK UndoSegSet (Pointer data)

{
  Int2              count;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  OMProcControlPtr  ompcp;
  Uint2             parenttype;
  Pointer           parentptr;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      break;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);
  count = DoOneSegUndo (sep);
  RemoveDupGenBankSets (sep);
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  if (count > 0) {
    PropagateFromGenBankBioseqSet (sep, FALSE);
    ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  }
  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK GenerateSeqAlignFromSeqEntry (Pointer data)

{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  SeqAnnotPtr       curr;
  ObjectIdPtr       oip;
  OMProcControlPtr  ompcp;
  SeqAnnotPtr       sap;
  SeqAnnotPtr PNTR  sapp;
  SeqEntryPtr       sep;
  UserFieldPtr      ufp;
  UserObjectPtr     uop;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = SeqMgrGetSeqEntryForData (ompcp->input_data);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  sap = SeqEntryToSeqAlign (sep, Seq_mol_na);
  if (sap != NULL && sap->type == 2) {

    oip = ObjectIdNew ();
    oip->str = StringSave ("Hist Seqalign");
    uop = UserObjectNew ();
    uop->type = oip;

    oip = ObjectIdNew();
    oip->str = StringSave ("Hist Seqalign");
    ufp = UserFieldNew ();
    ufp->choice = 4;
    ufp->data.boolvalue = TRUE;
    ufp->label = oip;

    uop->data = ufp;
	ValNodeAddPointer (&(sap->desc), Annot_descr_user, (Pointer) uop);

    sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
    if (sep != NULL && sep->data.ptrvalue != NULL) {
      sapp = NULL;
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
      ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
    }
  }
  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK GenerateSeqAlignFromSeqEntryProt (Pointer data)

{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  SeqAnnotPtr       curr;
  ObjectIdPtr       oip;
  OMProcControlPtr  ompcp;
  SeqAnnotPtr       sap;
  SeqAnnotPtr PNTR  sapp;
  SeqEntryPtr       sep;
  UserFieldPtr      ufp;
  UserObjectPtr     uop;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = SeqMgrGetSeqEntryForData (ompcp->input_data);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  sap = SeqEntryToSeqAlign (sep, Seq_mol_aa);
  if (sap != NULL && sap->type == 2) {

    oip = ObjectIdNew ();
    oip->str = StringSave ("Hist Seqalign");
    uop = UserObjectNew ();
    uop->type = oip;

    oip = ObjectIdNew();
    oip->str = StringSave ("Hist Seqalign");
    ufp = UserFieldNew ();
    ufp->choice = 4;
    ufp->data.boolvalue = TRUE;
    ufp->label = oip;

    uop->data = ufp;
	ValNodeAddPointer (&(sap->desc), Annot_descr_user, (Pointer) uop);

    sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
    if (sep != NULL && sep->data.ptrvalue != NULL) {
      sapp = NULL;
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
      ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
    }
  }
  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK ValidateSeqAlignFromSeqEntry (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = SeqMgrGetSeqEntryForData (ompcp->input_data);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  ValidateSeqAlignInSeqEntry (sep, FALSE, TRUE);
  return OM_MSG_RET_DONE;
}

static Boolean RawSeqLaunchFunc (GatherContextPtr gcp)

{
  BioseqPtr  bsp;
  Int2       handled;

  if (gcp == NULL) return TRUE;
  bsp = (BioseqPtr) gcp->userdata;
  if (bsp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ) {
    if (bsp == (BioseqPtr) gcp->thisitem) {
      WatchCursor ();
      handled = GatherProcLaunch (OMPROC_EDIT, FALSE, gcp->entityID, gcp->itemID,
                                  OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
      ArrowCursor ();
      if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
        /*
        Message (MSG_ERROR, "Unable to launch editor on sequence.");
        */
      }
      return FALSE;
    }
  }
  return TRUE;
}

extern Int2 LIBCALLBACK BioseqSegEditFunc (Pointer data)

{
  BioseqPtr         bsp;
  GatherScope       gs;
  SeqIdPtr          sip;
  SeqLocPtr         slp = NULL;
  OMProcControlPtr  ompcp;

  ompcp = (OMProcControlPtr) data;
  slp = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ_SEG :
      slp = (SeqLocPtr) ompcp->input_data;
      break;
   case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (slp == NULL) return OM_MSG_RET_ERROR;
  sip = SeqLocId (slp);
  if (sip == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFind (sip);
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  GatherEntity (ompcp->input_entityID, (Pointer) bsp, RawSeqLaunchFunc, &gs);

  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK CreateLocMsgFunc (OMMsgStructPtr ommsp)
{
  return OM_MSG_RET_OK;
}

static SeqLocPtr SeqLocCopyOne (SeqLocPtr slp)
{
  SeqLocPtr   slpnew, slptemp;

  slptemp = slp->next;
  slp->next = NULL;
  slpnew = AsnIoMemCopy ((Pointer) slp, (AsnReadFunc) SeqLocAsnRead,
                         (AsnWriteFunc) SeqLocAsnWrite);
  slp->next = slptemp;
  return slpnew;
}

static SeqFeatPtr SeqFeatCopy (SeqFeatPtr sfp)
{
  SeqFeatPtr  sfpnew;

  sfpnew = AsnIoMemCopy ((Pointer) sfp, (AsnReadFunc) SeqFeatAsnRead,
                         (AsnWriteFunc) SeqFeatAsnWrite);
  return sfpnew;
}

static Boolean ExplodeGroup (SeqEntryPtr sep, SeqFeatPtr sfp)

{
  SeqFeatPtr     sfpnew, sfpold, sfplast;
  ImpFeatPtr     ifp;
  SeqLocPtr      slphead, slp;
  GBQualPtr      gbq;
  Int2           count;
  Char           str [16];

  if (sfp == NULL || sfp->location == NULL) return FALSE;

/* save the seqloc (chain) */
  slphead = sfp->location;
  slp = SeqLocFindNext (slphead, NULL);
  if (slp == NULL) return FALSE;

/* trash the loc info in the impfeat */
  if (sfp->data.choice == SEQFEAT_IMP) {
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (ifp != NULL) {
      ifp->loc = MemFree (ifp->loc);
    }
  }

/* orig sfp is copied then orig sfp data is replaced */
  sfplast = sfp->next;
  sfpnew = SeqFeatCopy (sfp);

/* if exon, increment /number qualifier */
  gbq = NULL;
  count = 0;
  if (sfpnew != NULL && sfpnew->data.choice == SEQFEAT_IMP) {
    ifp = (ImpFeatPtr) sfpnew->data.value.ptrvalue;
    if (ifp != NULL) {
      if (StringICmp (ifp->key, "exon") == 0) {
        gbq = sfpnew->qual;
        while (gbq != NULL && StringICmp (gbq->qual, "number") != 0) {
          gbq = gbq->next;
        }
        if (gbq != NULL) {
          if (! StrToInt (gbq->val, &count)) {
            gbq = NULL;
          }
        }
      }
    }
  }

  sfp->location = SeqLocCopyOne (slp);
  sfpold = sfp;
  slp = SeqLocFindNext (slphead, slp);

/* clone as many more as required */
  while (slp != NULL) {
    if (slp->choice != SEQLOC_NULL) {
      if (gbq != NULL) {
        gbq->val = MemFree (gbq->val);
        count++;
        sprintf (str, "%d", (int) count);
        gbq->val = StringSave (str);
      }
      sfp = SeqFeatCopy (sfpnew);
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = SeqLocCopyOne (slp);
      sfp->partial = CheckSeqLocForPartial (sfp->location, NULL, NULL);
      sfpold->next = sfp;
      sfpold = sfp;
    }
    slp = SeqLocFindNext (slphead, slp);
  }
  sfpold->next = sfplast;
  sfpnew = SeqFeatFree (sfpnew);
  slphead = SeqLocFree (slphead);
  return TRUE;
}

static Int2 LIBCALLBACK GroupExplodeFunc (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqFeatPtr        sfp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_itemtype == 0)
    return OM_MSG_RET_ERROR;

  switch (ompcp->input_itemtype)
  {
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) ompcp->input_data;
      break;
    default:
      return OM_MSG_RET_ERROR;
  }

  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);

  if (ExplodeGroup (sep, sfp))
  {
    ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, ompcp->input_itemID,
                   ompcp->input_itemtype);
    return OM_MSG_RET_DONE;
  }
  else
    return OM_MSG_RET_ERROR;
}

static Boolean MakeExonsAndIntronsFromFeature (SeqEntryPtr sep, BioseqPtr bsp,
                                               SeqLocPtr location, SeqFeatPtr putafterhere)

{
  SeqFeatPtr  curr;
  Boolean     first;
  Int2        fuzz_from;
  Int2        fuzz_to;
  ImpFeatPtr  ifp;
  Int4        last;
  SeqLocPtr   next;
  SeqFeatPtr  putbeforehere;
  SeqFeatPtr  sfp;
  SeqLocPtr   slp;
  Int4        start;
  Int4        stop;
  Uint1       strand;
  Int4        tmp;

  if (sep == NULL || bsp == NULL || location == NULL || putafterhere == NULL) return FALSE;
  putbeforehere = putafterhere->next;
  curr = putafterhere;
  slp = SeqLocFindNext (location, NULL);
  if (slp == NULL) return FALSE;
  first = TRUE;
  last = 0;
  while (slp != NULL) {
    next = SeqLocFindNext (location, slp);
    if (slp->choice != SEQLOC_NULL) {
      start = GetOffsetInBioseq (slp, bsp, SEQLOC_START);
      stop = GetOffsetInBioseq (slp, bsp, SEQLOC_STOP);
      strand = SeqLocStrand (slp);
      if (strand > Seq_strand_both_rev && strand != Seq_strand_other) {
        strand = Seq_strand_unknown;
      }
      fuzz_from = -1;
      fuzz_to = -1;
      if (start > stop) {
        tmp = start;
        start = stop;
        stop = tmp;
      }
      if (! first) {
        sfp = SeqFeatNew ();
        if (sfp != NULL) {
          sfp->data.choice = SEQFEAT_IMP;
          if (strand == Seq_strand_minus) {
            AddIntToSeqFeat (sfp, stop + 1, last - 1, bsp,
                             fuzz_from, fuzz_to, strand);
          } else {
            AddIntToSeqFeat (sfp, last + 1, start - 1, bsp,
                             fuzz_from, fuzz_to, strand);
          }
          ifp = ImpFeatNew ();
          if (ifp != NULL) {
            sfp->data.value.ptrvalue = (Pointer) ifp;
            ifp->key = StringSave ("intron");
          }
          curr->next = sfp;
          curr = sfp;
        }
      }
      first = FALSE;
      if (strand == Seq_strand_minus) {
        last = start;
      } else {
        last = stop;
      }
      sfp = SeqFeatNew ();
      if (sfp != NULL) {
        sfp->data.choice = SEQFEAT_IMP;
        AddIntToSeqFeat (sfp, start, stop, bsp,
                         fuzz_from, fuzz_to, strand);
        ifp = ImpFeatNew ();
        if (ifp != NULL) {
          sfp->data.value.ptrvalue = (Pointer) ifp;
          ifp->key = StringSave ("exon");
        }
        curr->next = sfp;
        curr = sfp;
      }
    }
    slp = next;
  }
  curr->next = putbeforehere;
  return TRUE;
}

static Int2 LIBCALLBACK MakeExonIntron (Pointer data)

{
  BioseqPtr         nbsp;
  SeqEntryPtr       nsep;
  OMProcControlPtr  ompcp;
  SeqFeatPtr        sfp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_itemtype == 0 || ompcp->input_data == NULL)
    return OM_MSG_RET_ERROR;

  switch (ompcp->input_itemtype)
  {
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp->data.choice != SEQFEAT_CDREGION && sfp->data.choice != SEQFEAT_RNA) {
        return OM_MSG_RET_ERROR;
      }
      break;
    default:
      return OM_MSG_RET_ERROR;
  }

  sep = GetBestTopParentForItemID (ompcp->input_entityID,
                                   ompcp->input_itemID,
                                   ompcp->input_itemtype);
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL || nsep->choice != 1) return OM_MSG_RET_ERROR;
  nbsp = (BioseqPtr) nsep->data.ptrvalue;
  if (nbsp == NULL) return OM_MSG_RET_ERROR;

  if (MakeExonsAndIntronsFromFeature (sep, nbsp, sfp->location, sfp))
  {
    ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, ompcp->input_itemID,
                   ompcp->input_itemtype);
    return OM_MSG_RET_DONE;
  }
  else
    return OM_MSG_RET_ERROR;
}

static Int2 LIBCALLBACK ProtLocalIDtoGeneSyn (Pointer data)

{
  BioseqPtr          bsp = NULL;
  Char               buf [41];
  SeqMgrFeatContext  ccontext;
  SeqFeatPtr         cds;
  SeqMgrFeatContext  gcontext;
  SeqFeatPtr         gene;
  GeneRefPtr         grp;
  OMProcControlPtr   ompcp;
  BioseqPtr          pbsp;
  SeqEntryPtr        sep;
  SeqIdPtr           sip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_itemtype == 0 || ompcp->input_data == NULL)
    return OM_MSG_RET_ERROR;

  switch (ompcp->input_itemtype)
  {
    case OBJ_BIOSEQ:
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    default:
      return OM_MSG_RET_ERROR;
  }

  if (bsp == NULL) {
    return OM_MSG_RET_ERROR;
  }
  sep = SeqMgrGetSeqEntryForData (bsp);
  if (sep == NULL) {
    return OM_MSG_RET_ERROR;
  }

  cds = NULL;
  while ((cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, 0, &ccontext)) != NULL) {
    sip = SeqLocId (cds->product);
    if (sip != NULL) {
      pbsp = BioseqFind (sip);
      if (pbsp != NULL) {
        sip = SeqIdFindBest (pbsp->id, SEQID_LOCAL);
        if (sip != NULL) {
          SeqIdWrite (sip, buf, PRINTID_REPORT, sizeof (buf) - 1);
          grp = SeqMgrGetGeneXref (cds);
          if (grp != NULL && SeqMgrGeneIsSuppressed (grp)) {
            continue;
          }
          if (grp == NULL) {
            gene = SeqMgrGetOverlappingGene (cds->location, &gcontext);
            if (gene == NULL) {
              gene = CreateNewFeature (sep, NULL, SEQFEAT_GENE, NULL);
              if (gene != NULL) {
                grp = GeneRefNew ();
                gene->data.value.ptrvalue = (Pointer) grp;
                gene->location = SeqLocFree (gene->location);
                gene->location = AsnIoMemCopy ((Pointer) cds->location,
                                               (AsnReadFunc) SeqLocAsnRead,
                                               (AsnWriteFunc) SeqLocAsnWrite);
              }
            }
            if (gene != NULL) {
              grp = (GeneRefPtr) gene->data.value.ptrvalue;
            }
          }
          if (grp != NULL) {
            ValNodeCopyStr (&(grp->syn), 0, buf);
          }
        }
      }
    }
  }

  BasicSeqEntryCleanup (sep);
  SeriousSeqEntryCleanup (sep, NULL, NULL);

  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, ompcp->input_itemID,
                 ompcp->input_itemtype);
  return OM_MSG_RET_DONE;
}

#define BLACK     0
#define RED       4
#define GREEN     2
#define BLUE      1
#define CYAN      3
#define MAGENTA   5
#define YELLOW    6
#define WHITE    15
#define GRAY      8
#define LTGRAY    7

#define DKCYAN   21
#define DKGREEN  22
#define DKBLUE  23

#define ORF_LENGTH 10


typedef struct orfviewform
{
	FORM_MESSAGE_BLOCK
	IcoN 		icon;
	WindoW 		w;
	DoC 		doc;
	BioseqPtr	bsp;
	Int2		gcode;
	ValNodePtr	orfs;
	Int2 		frame, strand;
	Int4 		from, to;
	double 		dx, dy;
	RecT		mi;
	Boolean		orf_only;
	SeqLocPtr	select_orf;
	Uint2 		bsp_entityID;
	Uint2 		bsp_itemID;
	Uint2 		len;		/* minimum length of the ORF shown */
	Boolean     standAlone;
	Boolean     alt_start;
} OrfViewForm, PNTR OrfViewFormPtr;

static RecT mi0 = {  24,  56, 460, 200 };

Uint1 AAForCodon (Uint1Ptr codon, CharPtr codes); /* in seqport.c */

static void dkCyan (void)
{
	SelectColor(0, 203, 196);
}

static void dkBlue (void)
{
	SelectColor(0, 0, 196);
}

static void dkGreen (void)
{
	SelectColor(0, 203, 0);
}

static void SetIntColor (int color)
{
	switch (color) {
		case BLACK: Black(); break;
		case RED: Red(); break;
		case GREEN: Green(); break;
		case BLUE: Blue(); break;
		case CYAN: Cyan(); break;
		case MAGENTA: Magenta(); break;
		case YELLOW: Yellow(); break;
		case WHITE: White(); break;
		case GRAY: Gray(); break;
		case LTGRAY: LtGray(); break;
		case DKCYAN: dkCyan(); break;
		case DKBLUE: dkBlue(); break;
		case DKGREEN: dkGreen(); break;
	}
}

static void Frame2Rect (RecT PNTR a, int color1, int color2)
{
	MoveTo(a->left, a->top);
	SetIntColor(color1);
	LineTo(a->right, a->top);
	LineTo(a->right, a->bottom);
	SetIntColor(color2);
	LineTo(a->left, a->bottom);
	LineTo(a->left, a->top);
}

static void Frame3d (RecT PNTR a)
{
	int i, in = 3;
	
	SetIntColor(LTGRAY);
	PaintRect(a);
	for (i=0; i < in; i++) {
		if (i > 0) {
			InsetRect(a, 1, 1);
		}
		Frame2Rect(a, WHITE, GRAY);
	}
	for (i=0; i < in; i++) {
		InsetRect(a, 1, 1);
	}
	for (i=0; i < in; i++) {
		InsetRect(a, 1, 1);
		Frame2Rect(a, GRAY, WHITE);
	}
	SelectColor(0, 203, 196);
	PaintRect(a);
}

static void DrawColorLine (int x0, int y0, int x1, int y1, int color)
{
	SetIntColor(color);
	MoveTo(x0, y0);
	LineTo(x1, y1);
}

static void Rect3d (RecT PNTR a, int color)
{
	Frame2Rect(a, GRAY, WHITE);
	InsetRect(a, 1, 1);
	Frame2Rect(a, GRAY, WHITE);
	InsetRect(a, 1, 1);
	if (color != -1) {
		SetIntColor(color);
		PaintRect(a);
	}
}



static void QuitProc (ButtoN b)

{
  QuitProgram ();
}

static void CloseProc (ButtoN b)
{ 
	WindoW   w;

    w = ParentWindow (b);
    Remove (w);
}

static void draw_rect(SeqPortPtr spp, Int2 ir, RecT PNTR frect, CharPtr vals, CharPtr codes, Boolean paint, OrfViewFormPtr ovp, Int2 strand)
{
	Int4 pos;
	Uint1 codon[3], aa;
	Int4 len;
	double x;
	RecT r;
	
    	SeqPortSeek(spp, ir, SEEK_SET);
		frect->bottom = frect->top + 15;
		len = spp->totlen;
		if (paint) {
			dkGreen();
		} else {
			dkCyan();
		}
		dkCyan();
		PaintRect(frect);
		Black();
		FrameRect(frect);
		for (pos=0; pos < len-2; pos += 3) {
			codon[0] = SeqPortGetResidue(spp);
			codon[1] = SeqPortGetResidue(spp);
			codon[2] = SeqPortGetResidue(spp);
			aa = AAForCodon(codon, codes);
			if (aa == '*') {
				if (strand == Seq_strand_plus) {
					x = frect->left + (pos+ir)*ovp->dx;
				} else {
					x = frect->left + (len+2-(pos+ir))*ovp->dx;
				}
				DrawColorLine(x, frect->top+1, x, frect->bottom-1, RED);
				DrawColorLine(x+1, frect->top+1, x+1, frect->bottom-1, RED);
			}
			if (ovp->alt_start == TRUE) {
				aa = AAForCodon(codon, vals);
			}
			if (aa == 'M') {
				if (strand == Seq_strand_plus) {
					x = frect->left + (pos+ir)*ovp->dx;
				} else {
					x = frect->left + (len+2-(pos+ir))*ovp->dx;
				}
				DrawColorLine(x, frect->top+1, x, frect->bottom-1, WHITE);
				DrawColorLine(x+1, frect->top+1, x+1, frect->bottom-1, WHITE);
			}
		}
		if (paint) {
			r.top = frect->top + 1;
			r.bottom = frect->bottom - 1;
			r.left = frect->left + ovp->from*ovp->dx + 2;
			r.right = frect->left + ovp->to*ovp->dx - 2;
			Magenta();
			PaintRect(&r);
		}
		frect->top = frect->bottom + 4;
}

static void draw_frame(Int2 ir, RecT PNTR frect, Boolean paint, OrfViewFormPtr ovp, Int2 strand)
{
	RecT r;  /* for ORF */
	ValNodePtr vnp;
	SeqLocPtr slp;
	SeqIntPtr sip;

	frect->bottom = frect->top + 15;
	LtGray();
	PaintRect(frect);
	Black();
	FrameRect(frect);
	for (vnp = ovp->orfs; vnp; vnp=vnp->next) {
		if (vnp->choice == ir) {
			slp = vnp->data.ptrvalue;
			if (slp == NULL) {
				continue;
			}
			sip = slp->data.ptrvalue;
			if (sip == NULL) {
				continue;
			}
			if (sip->strand == strand) {
				r.top = frect->top + 1;
				r.bottom = frect->bottom - 1;
				r.left = frect->left + sip->from*ovp->dx;
				r.right = frect->left + sip->to*ovp->dx;
				dkCyan();
				PaintRect(&r);
				Black();
				r.top = frect->top;
				r.bottom = frect->bottom;
				FrameRect(&r);
			}
			if (paint) {
				r.top = frect->top + 1;
				r.bottom = frect->bottom - 1;
				r.left = frect->left + ovp->from*ovp->dx + 2;
				r.right = frect->left + ovp->to*ovp->dx - 2;
				Magenta();
				PaintRect(&r);
			}
		}
	}
	frect->top = frect->bottom + 4;
	return;
}

static void draw_strands(OrfViewFormPtr ovp)
{
	Int2 ir, gcode;
	RecT frect;
	SeqPortPtr spp;
	GeneticCodePtr gcp;
	CharPtr vals, codes;
	ValNodePtr vnp;
	Boolean paint;
	RecT PNTR r;
	
	r = &(ovp->mi);
	gcode = ovp->gcode;
	gcp = GeneticCodeFind(gcode, NULL);   /* use universal */
	vals = NULL;
	codes = NULL;
	for (vnp = (ValNodePtr)gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next)
	{
		if (vnp->choice == 6)   /* sncbieaa */
			vals = (CharPtr)vnp->data.ptrvalue;
		else if (vnp->choice == 3)  /* ncbieaa */
			codes = (CharPtr)vnp->data.ptrvalue;
	}
	if (vals == NULL) {
		vals = codes;
	}
	frect.left = r->left + 11;
	frect.right = r->right - 11;
	frect.top = r->top + 4;
	ovp->dx = (frect.right - frect.left - 1.) / (ovp->bsp)->length;
	ovp->dy = (r->bottom - r->top - 1.) / 6.;
	spp = SeqPortNew(ovp->bsp, 0, -1, Seq_strand_plus, Seq_code_ncbi4na);
	for (ir=0; ir < 3; ir++) {
		if (ovp->from == 0 && ovp->to == 0) {
			paint = FALSE;
		} else if (ovp->strand == Seq_strand_plus && ovp->frame == ir) {
			paint = TRUE;
		} else {
			paint = FALSE;
		}
		if (ovp->orf_only) {
			draw_frame(ir, &frect, paint, ovp, Seq_strand_plus);
		} else {
			draw_rect(spp, ir, &frect, vals, codes, paint, ovp, Seq_strand_plus);
		}
	}
	SeqPortFree(spp);
	spp = SeqPortNew(ovp->bsp, 0, -1, Seq_strand_minus, Seq_code_ncbi4na);
	frect.top += 7;
	for (ir=0; ir < 3; ir++) {
		if (ovp->from == 0 && ovp->to == 0) {
			paint = FALSE;
		} else if (ovp->strand == Seq_strand_minus && ovp->frame == ir) {
			paint = TRUE;
		} else {
			paint = FALSE;
		}
		if (ovp->orf_only) {
			draw_frame(ir, &frect, paint, ovp, Seq_strand_minus);
		} else {
			draw_rect(spp, ir, &frect, vals, codes, paint, ovp, Seq_strand_minus);
		}
	}
	
	SeqPortFree(spp);
}

static void DrawIcon(IcoN ic0)
{
	RecT r;
	OrfViewFormPtr ovp;

	ovp = (OrfViewFormPtr) GetObjectExtra (ic0);
	ObjectRect(ovp->icon, &r);
	ovp->mi.left = mi0.left + r.left;
	ovp->mi.right = mi0.right + r.left;
	ovp->mi.top = mi0.top + r.top;
	ovp->mi.bottom = mi0.bottom + r.top;
	Frame3d(&r);
	Rect3d(&(ovp->mi), LTGRAY);
	draw_strands(ovp);
}

static void notify( DoC d, Int2 item, Int2 raw, Int2 col, Boolean event)
{
	ValNodePtr vnp;
	SeqLocPtr slp, slptmp;
	SeqIntPtr sip;
	Int2 i;
	Int2	itemOld1;
	Int2	itemOld2;
	Boolean status;
	OrfViewFormPtr ovp;
	Int2 top, bottom;
	BaR sb;
	Int2 startsAt;
	
	if( item == 0) {
		return;
	}
	ovp = (OrfViewFormPtr) GetObjectExtra (d);
	for (vnp = ovp->orfs, i = 1; i < item && vnp; i++, vnp=vnp->next) continue;
	if (vnp == NULL) {
		return;
	}
	if (ItemIsVisible (d, item, &top, &bottom, NULL) == FALSE) {
      GetItemParams (d, item, &startsAt, NULL, NULL, NULL, NULL);
      sb = GetSlateVScrollBar ((SlatE) d);
      CorrectBarValue (sb, startsAt);
	}
    GetDocHighlight(d, &itemOld1, &itemOld2);
	SetDocHighlight(d, item, item);
	UpdateDocument(d, itemOld1, itemOld2);
	UpdateDocument(d, item, item);

	ovp->frame = vnp->choice;
	slp = vnp->data.ptrvalue;
	sip = slp->data.ptrvalue;
	ovp->strand = sip->strand;
	ovp->from = sip->from;
	ovp->to = sip->to;
	DrawIcon(ovp->icon);
	status = GetStatus(ovp->icon);
	SetStatus(ovp->icon, !status);

    if (! ovp->standAlone) {
		slptmp = AsnIoMemCopy(slp, (AsnReadFunc) SeqLocAsnRead, 
				(AsnWriteFunc) SeqLocAsnWrite);
		ovp->select_orf = slptmp;
		ObjMgrSelect (ovp->bsp_entityID, ovp->bsp_itemID, OBJ_BIOSEQ, 
				OM_REGION_SEQLOC, slptmp);
	}

	Update();
}

static void myrelease(IcoN ic, PoinT pt)
{
	RecT r;
	Int2 ir;
	Int4 pos;
	ValNodePtr vnp;
	SeqLocPtr slp;
	SeqIntPtr sip;
	Boolean status;
	OrfViewFormPtr ovp;
	Int2	item;
	Int2 top, bottom;
	BaR sb;
	Int2 startsAt;
	DoC doc;
	
	ovp = GetObjectExtra((IcoN) ic);
	ObjectRect(ic, &r);
	if (!PtInRect(pt, &(ovp->mi))) {
		return;
	}
	doc = ovp->doc;
	r.left = ovp->mi.left + 11;
	r.right = ovp->mi.right - 11;
	r.top = ovp->mi.top + 4;
	ovp->dx = (r.right - r.left - 1.) / (ovp->bsp)->length;
	ovp->dy = (ovp->mi.bottom - ovp->mi.top - 1.) / 6.;
	for (ir=0; ir < 3; ir++) {
		r.bottom = r.top + 15;
		if (pt.y < r.bottom && pt.y > r.top) {
			break;
		}
		r.top = r.bottom + 4;
	}
	if (ir < 3) {
		pos = (pt.x - r.left)/ovp->dx - ir;
		for (vnp=ovp->orfs, item=1; vnp; vnp=vnp->next, item++) {
			slp = vnp->data.ptrvalue;
			sip = slp->data.ptrvalue;
			if (vnp->choice != ir || sip->strand != Seq_strand_plus) {
				continue;
			}
			if (pos < sip->to && pos > sip->from) {
				break;
			}
		}
		if (vnp == NULL) {
			Beep();
			return;
		}
		ovp->frame = ir;
		ovp->strand = sip->strand;
		ovp->from = sip->from;
		ovp->to = sip->to;
		ovp->select_orf = AsnIoMemCopy(slp, (AsnReadFunc) SeqLocAsnRead, 
		(AsnWriteFunc) SeqLocAsnWrite);
		status = GetStatus(ic);
		SetStatus(ic, !status);
		Update();
		Select(doc);
		SetDocHighlight(doc, item, item);
		if (ItemIsVisible (doc, item, &top, &bottom, NULL) == FALSE) {
      	GetItemParams (doc, item, &startsAt, NULL, NULL, NULL, NULL);
      	sb = GetSlateVScrollBar ((SlatE) doc);
      	CorrectBarValue (sb, startsAt);
		}
		UpdateDocument(doc, 0, 0);

		if (! ovp->standAlone) {
			if (ovp->select_orf != NULL) {
				ObjMgrSelect (ovp->bsp_entityID, ovp->bsp_itemID, OBJ_BIOSEQ, 
					OM_REGION_SEQLOC, ovp->select_orf);
			}
		}
		Update();
		return;
	}
	r.top += 7;
	for (ir=0; ir < 3; ir++) {
		r.bottom = r.top + 15;
		if (pt.y < r.bottom && pt.y > r.top) {
			break;
		}
		r.top = r.bottom + 4;
	}
	if (ir < 3) {
		pos = (pt.x - r.left)/ovp->dx - ir;
		for (vnp=ovp->orfs, item=1; vnp; vnp=vnp->next, item++) {
			slp = vnp->data.ptrvalue;
			sip = slp->data.ptrvalue;
			if (vnp->choice != ir  || sip->strand != Seq_strand_minus) {
				continue;
			}
			if (pos < sip->to && pos > sip->from) {
				break;
			}
		}
		if (vnp == NULL) {
			Beep();
			return;
		}
		ovp->frame = ir;
		ovp->strand = sip->strand;
		ovp->from = sip->from;
		ovp->to = sip->to;
		ovp->select_orf = AsnIoMemCopy(slp, (AsnReadFunc) SeqLocAsnRead, 
					(AsnWriteFunc) SeqLocAsnWrite);
		status = GetStatus(ic);
		SetStatus(ic, !status);
		Update();
		Select(doc);
		SetDocHighlight(doc, item, item);
		if (ItemIsVisible (doc, item, &top, &bottom, NULL) == FALSE) {
      	GetItemParams (doc, item, &startsAt, NULL, NULL, NULL, NULL);
      	sb = GetSlateVScrollBar ((SlatE) doc);
      	CorrectBarValue (sb, startsAt);
		}
		UpdateDocument(doc, 0, 0); 
	}
	if (! ovp->standAlone) {
		if (ovp->select_orf != NULL) {
			ObjMgrSelect (ovp->bsp_entityID, ovp->bsp_itemID, OBJ_BIOSEQ, 
				OM_REGION_SEQLOC, ovp->select_orf);
		}
	}
	Update();
	return;
}

/* show all ORF from List without starts and stops */
static void ORFProc(ButtoN b)
{
	OrfViewFormPtr ovp;
	Boolean status;
	
	if ((ovp = (OrfViewFormPtr) GetObjectExtra (b)) == NULL) {
		return;
	}
	ovp->orf_only = !ovp->orf_only;
	draw_strands(ovp);
	status = GetStatus(ovp->icon);
	SetStatus(ovp->icon, !status);
	Update();
}
static void AddOrfListToDoc(DoC doc, ValNodePtr list)
{
	ParPtr	par;
	ValNodePtr vnp;
	SeqLocPtr slp, tmp;
	SeqIntPtr sip;
	Boolean minus;
	Int2 l, i;
	CharPtr buf, str;
	OrfViewFormPtr ovp;

	Reset(doc);
	ovp = (OrfViewFormPtr) GetObjectExtra(doc);
	par = (ParPtr) MemNew(sizeof (ParData));
  	par->openSpace = FALSE;
  	par->keepWithNext = 1;
  	par->keepTogether = 1;
	par->newPage  = 1;
  	par->tabStops     = 1;
	for (vnp = list, i=0; vnp; vnp=vnp->next, i++) {
		minus = FALSE;
		tmp = (SeqLocPtr) vnp->data.ptrvalue;
		slp = AsnIoMemCopy ((Pointer) tmp, (AsnReadFunc) SeqLocAsnRead,
                                        (AsnWriteFunc) SeqLocAsnWrite);
		if (slp == NULL) {
			continue;
		}                                        
        slp->next = NULL;
		sip = (SeqIntPtr) slp->data.ptrvalue;
		if (sip->strand == Seq_strand_minus) {
			sip->strand = Seq_strand_plus;
			minus = TRUE;
		}
		str = FlatLoc(ovp->bsp, slp);
		MemFree(slp);
		l = StringLen(str);
		if (minus == TRUE) {
			buf = MemNew(l+4);
			sprintf(buf, "c(%s)", str);
			AppendText(doc, buf, par, NULL, NULL); 
			MemFree(buf);
			sip->strand = Seq_strand_minus;
		} else {
			AppendText(doc, str, par, NULL, NULL); 
		}		
	}
	UpdateDocument(doc, 0, 0);
}

static void AltProc(ButtoN b)
{
	OrfViewFormPtr ovp;
	Boolean status;
	
	if ((ovp = (OrfViewFormPtr) GetObjectExtra (b)) == NULL) {
		return;
	}
	ovp->orfs =  GetAltList(ovp->orfs);
	AddOrfListToDoc(ovp->doc, ovp->orfs);
	ovp->from = 0;
	ovp->to = 0;
	draw_strands(ovp); 
	status = GetStatus(ovp->icon);
	SetStatus(ovp->icon, !status);
	if (ovp->alt_start) {
		SetTitle(b, "Alternative Initiation Codon");
		ovp->alt_start = FALSE;
	} else {
		SetTitle(b, "Standard Initiation Codon");
		ovp->alt_start = TRUE;
	}
	Update();
}

static void ChangeAAcutoff(PopuP p)
{
	OrfViewFormPtr ovp;
	Int2 i;
	Boolean status;

	ovp = (OrfViewFormPtr) GetObjectExtra(p);
	i = GetValue((PopuP) p);
	switch (i) {
		case 1:
			ovp->len = 10;
			break;
		case 2:
			ovp->len = 50;
			break;
		case 3:
			ovp->len = 100;
			break;
		default:
			ovp->len = 3;
			break;
	}
	ovp->orfs =  GetOrfList(ovp->bsp, ovp->len);
	AddOrfListToDoc(ovp->doc, ovp->orfs);
	draw_strands(ovp);
	status = GetStatus(ovp->icon);
	SetStatus(ovp->icon, !status);
	Update();
}

static Int2 LIBCALLBACK OrfViewerMsgFunc (OMMsgStructPtr ommsp)

{
  ObjMgrDataPtr   omdp;
  OMUserDataPtr   omudp;
  OrfViewFormPtr  ovp;

  omudp = (OMUserDataPtr)(ommsp->omuserdata);
  if (omudp == NULL) return OM_MSG_RET_ERROR;
  ovp = (OrfViewFormPtr) omudp->userdata.ptrvalue;
  if (ovp == NULL) return OM_MSG_RET_ERROR;
  switch (ommsp->message) {
    case OM_MSG_DEL:
      omdp = ObjMgrGetData (ommsp->entityID);
      if (omdp != NULL) {
        if (ObjMgrWholeEntity (omdp, ommsp->itemID, ommsp->itemtype)) {
          if (ovp != NULL) {
            Remove (ovp->form);
          }
          return OM_MSG_RET_OK;
        }
      }
      break;
    default :
      break;
  }
  return OM_MSG_RET_OK;
}

static void CleanupOrfViewer (GraphiC g, VoidPtr data)

{
  OrfViewFormPtr  ovp;

  ovp = (OrfViewFormPtr) data;
  if (ovp != NULL && ovp->input_entityID > 0) {
    ObjMgrFreeUserData (ovp->input_entityID, ovp->procid, ovp->proctype, ovp->userkey);
  }
  StdCleanupFormProc (g, data);
}

extern void LaunchOrfViewer (BioseqPtr bsp, Uint2 entityID, Uint2 itemID, Boolean standAlone)
{
	ButtoN qu, b1, b2;
	GrouP g;
	Int2 i, l;
	ValNodePtr vnp;
	SeqLocPtr slp, tmp;
	SeqIntPtr sip;
	CharPtr str, buf;
	DoC doc;
	OrfViewFormPtr ovp;
	WindoW w;
	IcoN ic;
	ParPtr	par;
	Boolean minus;
	PopuP pu;
	ObjMgrPtr omp;
	ObjMgrProcPtr ompp;
	OMUserDataPtr omudp;

	WatchCursor ();
	Update ();
	ovp = (OrfViewFormPtr) MemNew (sizeof (OrfViewForm));
	ovp->bsp = bsp;
	ovp->select_orf = NULL;
	ovp->input_entityID = entityID;
	ovp->bsp_entityID = entityID;
	ovp->bsp_itemID = itemID;
	ovp->orfs =  GetOrfList(bsp, ORF_LENGTH);
	ovp->orf_only = TRUE;
	ovp->gcode = 1;
	ovp->standAlone = standAlone;
	ovp->alt_start = FALSE;
	w = FixedWindow(-50, -33, -10, -10, "Orf Finder", NULL);
	SetObjectExtra (w, ovp, CleanupOrfViewer);
	ovp->form = (ForM) w;
    g = HiddenGroup (w, 5, 0, NULL);
    SetGroupSpacing (g, 6, 0);
	if (ovp->standAlone) {
		qu = PushButton(g, "Quit", QuitProc);
	} else {
		qu = PushButton(g, "Close", CloseProc);
	}
	b1 = PushButton(g, "ORF", ORFProc);
	SetObjectExtra (b1, ovp, NULL);
	b2 = PushButton(g, "Alternative Initiation Codon", AltProc);
	SetObjectExtra (b2, ovp, NULL);
	StaticPrompt(g, "ORF length", 0, 16, systemFont, '1');
	pu = PopupList(g, TRUE, ChangeAAcutoff );
	PopupItem(pu, "10");
	PopupItem(pu, "50");
	PopupItem(pu, "100");
	SetValue(pu, 1);
	SetObjectExtra (pu, ovp, NULL);
	
    g = HiddenGroup (w, 2, 0, NULL);
	ic = IconButton(g, 500, 240, 
		DrawIcon, NULL, NULL, NULL, NULL, myrelease);
	if (ovp->select_orf != NULL) {
	}
	ovp->icon = ic;
    SetObjectExtra (ic, ovp, NULL);
	doc = DocumentPanel(g, 10 * stdCharWidth, 240);
    SetObjectExtra (doc, ovp, NULL);
/****/
	par = (ParPtr) MemNew(sizeof (ParData));
  	par->openSpace = FALSE;
  	par->keepWithNext = 1;
  	par->keepTogether = 1;
	par->newPage  = 1;
  	par->tabStops     = 1;
	for (vnp = ovp->orfs, i=0; vnp; vnp=vnp->next, i++) {
		minus = FALSE;
		tmp = (SeqLocPtr) vnp->data.ptrvalue;
		slp = AsnIoMemCopy ((Pointer) tmp, (AsnReadFunc) SeqLocAsnRead,
                                        (AsnWriteFunc) SeqLocAsnWrite);
        if (slp == NULL) {
        	continue;
        }
        slp->next = NULL;
		sip = (SeqIntPtr) slp->data.ptrvalue;
		if (sip->strand == Seq_strand_minus) {
			sip->strand = Seq_strand_plus;
			minus = TRUE;
		}
		str = FlatLoc(ovp->bsp, slp);
		MemFree(slp);
		l = StringLen(str);
		if (minus == TRUE) {
			buf = MemNew(l+4);
			sprintf(buf, "c(%s)", str);
			AppendText(doc, buf, par, NULL, NULL); 
			MemFree(buf);
			sip->strand = Seq_strand_minus;
		} else {
			AppendText(doc, str, par, NULL, NULL); 
		}		
	}
/****/
    ovp->doc = doc;
	SetDocNotify(doc, notify);
	/*ovp = (OrfViewFormPtr) GetObjectExtra(ic);*/
    omp = ObjMgrGet ();
    if (omp != NULL) {
		ompp = ObjMgrProcFind (omp, 0, "ORF Finder", OMPROC_FILTER);
		if (ompp != NULL) {
			ovp->procid = ompp->procid;
			ovp->proctype = OMPROC_FILTER;
			ovp->userkey = OMGetNextUserKey ();
			omudp = ObjMgrAddUserData (ovp->input_entityID, ompp->procid,
	                           OMPROC_FILTER, ovp->userkey);
			if (omudp != NULL) {
				omudp->userdata.ptrvalue = (Pointer) ovp;
				omudp->messagefunc = OrfViewerMsgFunc;
			}
		}
	}
	RealizeWindow (w);
	Show(w);
	ArrowCursor ();
	Update ();
	if (ovp->standAlone) {
		ProcessEvents();
	} else {
		return;
	}
}

static Int2 LIBCALLBACK OrfFindFunc (Pointer data)

{
  BioseqPtr         bsp;
  OMProcControlPtr  ompcp;


	ompcp = (OMProcControlPtr) data;   /* always do this cast */

	if (ompcp == NULL || ompcp->input_itemtype == 0)
		return OM_MSG_RET_ERROR;

	switch (ompcp->input_itemtype)
	{
		case OBJ_BIOSEQ:
			bsp = (BioseqPtr) ompcp->input_data;
		break;
		default:
			return OM_MSG_RET_ERROR;
	}

	LaunchOrfViewer (bsp, ompcp->input_entityID, ompcp->input_itemID, FALSE);
	return OM_MSG_RET_DONE;
}

typedef struct mergedata {
  SeqLocPtr     slp;
} MergeData, PNTR MergeDataPtr;

static Boolean AddToSeqLoc (GatherContextPtr gcp)

{
  BioseqPtr     bsp;
  MergeDataPtr  mdp;
  SeqFeatPtr    sfp;
  SeqLocPtr     slp;

  mdp = (MergeDataPtr) gcp->userdata;
  if (mdp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->location == NULL) return TRUE;
  bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
  if (bsp == NULL) return TRUE;
  slp = SeqLocMerge (bsp, sfp->location, mdp->slp, FALSE, TRUE, FALSE);
  mdp->slp = SeqLocFree (mdp->slp);
  mdp->slp = slp;
  return TRUE;
}

static SeqLocPtr MergeSelectedFeatureIntervals (void)

{
  MergeData     md;
  SelStructPtr  sel;

  md.slp = NULL;
  for (sel = ObjMgrGetSelected (); sel != NULL; sel = sel->next) {
    GatherItem (sel->entityID, sel->itemID, sel->itemtype,
                (Pointer) &md, AddToSeqLoc);
  }
  return md.slp;
}

static Int2 LIBCALLBACK IntervalCombineFunc (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqLocPtr         slp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  slp = MergeSelectedFeatureIntervals ();
  if (slp == NULL) return OM_MSG_RET_ERROR;
  ObjMgrRegister (OBJ_SEQLOC, (Pointer) slp);
  return OM_MSG_RET_DONE;
}

static CharPtr objmgrtypestrs [] = {
  "OBJ_ALL", "OBJ_SEQENTRY", "OBJ_BIOSEQ", "OBJ_BIOSEQSET", "OBJ_SEQDESC",
  "OBJ_SEQANNOT", "OBJ_ANNOTDESC", "OBJ_SEQFEAT", "OBJ_SEQALIGN", "OBJ_SEQGRAPH",
  "OBJ_SEQSUB", "OBJ_SUBMIT_BLOCK", "OBJ_SEQSUB_CONTACT", "13", "OBJ_BIOSEQ_MAPFEAT",
  "OBJ_BIOSEQ_SEG", "OBJ_SEQHIST", "OBJ_SEQHIST_ALIGN", "OBJ_BIOSEQ_DELTA", "19",
  "OBJ_PUB", "OBJ_SEQFEAT_CIT", "OBJ_SEQSUB_CIT", "OBJ_MEDLINE_ENTRY", "OBJ_PUB_SET",
  "OBJ_SEQLOC", "OBJ_SEQID", "OBJ_SEQCODE", "OBJ_SEQCODE_SET", "OBJ_GENETIC_CODE",
  "OBJ_GENETIC_CODE_SET", "OBJ_TEXT_REPORT", "OBJ_FASTA", "OBJ_VIBRANT_PICTURE", "OBJ_PROJECT"
};

static CharPtr temploadstrs [] = {
  "TL_NOT_TEMP", "TL_LOADED", "TL_CACHED"
};

static CharPtr proctypestrs [] = {
  "0", "OMPROC_OPEN", "OMPROC_DELETE", "OMPROC_VIEW", "OMPROC_EDIT",
  "OMPROC_SAVE", "OMPROC_CUT", "OMPROC_COPY", "OMPROC_PASTE", "OMPROC_ANALYZE",
  "OMPROC_FIND", "OMPROC_REPLACE", "OMPROC_FILTER", "OMPROC_FETCH",
};

static void PrintABool (FILE *fp, CharPtr str, Boolean val)

{
  if (val) {
    fprintf (fp, "%s TRUE\n", str);
  } else {
    fprintf (fp, "%s FALSE\n", str);
  }
}

Int2 LIBCALLBACK VSMPictMsgFunc PROTO((OMMsgStructPtr ommsp));

static void ReportOnEntity (ObjMgrDataPtr omdp, ObjMgrPtr omp, Boolean selected, Uint2 itemID,
                            Uint2 itemtype, Int2 index, FILE *fp)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Char           buf [50];
  OMUserDataPtr  omudp;
  VSMPictPtr     vsmpp;

  if (omdp == NULL || fp == NULL) return;
  if (selected) {
    fprintf (fp, "Data Element\n\n");
    fprintf (fp, "  EntityID %d selected\n", (int) omdp->EntityID);
    fprintf (fp, "  ItemID %d, Itemtype %d\n", (int) itemID, (int) itemtype);
  } else if (omdp->parentptr == NULL) {
    fprintf (fp, "Top Data Element %d\n\n", (int) index);
    fprintf (fp, "  EntityID %d\n", (int) omdp->EntityID);
  } else {
    fprintf (fp, "Inner Data Element %d\n\n", (int) index);
    fprintf (fp, "  EntityID %d\n", (int) omdp->EntityID);
  }
  if (omdp->datatype < OBJ_MAX) {
    fprintf (fp, "  Datatype %s", objmgrtypestrs [omdp->datatype]);
    if (omdp->datatype == OBJ_BIOSEQ) {
      bsp = (BioseqPtr) omdp->dataptr;
      if (bsp != NULL) {
        SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
        fprintf (fp, " %s, length %ld", buf, (long) bsp->length);
      }
    } else if (omdp->datatype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) omdp->dataptr;
      if (bssp != NULL) {
        fprintf (fp, " class %d", (int) bssp->_class);
      }
    }
    fprintf (fp, "\n");
  } else {
    fprintf (fp, "  Unregistered datatype %d\n", (int) omdp->datatype);
  }
  fprintf (fp, "  Lockcnt %d\n", (int) omdp->lockcnt);
  if (omdp->tempload < 3) {
    fprintf (fp, "  Tempload %s\n", temploadstrs [omdp->tempload]);
  } else {
    fprintf (fp, "  Unrecognized tempload %d\n", (int) omdp->tempload);
  }
  PrintABool (fp, "  Clipboard", omdp->clipboard);
  PrintABool (fp, "  Dirty", omdp->dirty);
  PrintABool (fp, "  Being_freed", omdp->being_freed);
  PrintABool (fp, "  Free", omdp->free);
  fprintf (fp, "\n");
  for (omudp = omdp->userdata; omudp != NULL; omudp = omudp->next) {
    if (omudp->proctype <= OMPROC_MAX) {
      fprintf (fp, "    Proctype %s\n", proctypestrs [omudp->proctype]);
    } else {
      fprintf (fp, "    Unrecognized proctype %d\n", (int) omudp->proctype);
    }
    fprintf (fp, "    Procid %d\n", (int) omudp->procid);
    fprintf (fp, "    Userkey %d\n", (int) omudp->userkey);
    if (omudp->messagefunc == VSMPictMsgFunc) {
      vsmpp = (VSMPictPtr) omudp->userdata.ptrvalue;
      if (vsmpp != NULL) {
        if (vsmpp->s != NULL) {
          fprintf (fp, "    VSMPictPtr segment not NULL\n");
        } else {
          fprintf (fp, "    VSMPictPtr segment is NULL\n");
        }
      }
    }
    fprintf (fp, "\n");
  }
}

static Int2 LIBCALLBACK DesktopReportFunc (Pointer data)

{
  FILE           *fp;
  Int2           j;
  Int2           num;
  ObjMgrPtr      omp;
  ObjMgrDataPtr  omdp;
  ObjMgrDataPtr  PNTR omdpp;
  Char           path [PATH_MAX];
  SelStructPtr   sel;

  omp = ObjMgrGet ();
  if (omp == NULL) return OM_MSG_RET_DONE;
  TmpNam (path);
  fp = FileOpen (path, "w");
  fprintf (fp, "Object Manager\n\n");
  fprintf (fp, "  HighestEntityID %d\n", (int) omp->HighestEntityID);
  fprintf (fp, "  Totobj %d\n", (int) omp->totobj);
  fprintf (fp, "  Currobj %d\n", (int) omp->currobj);
  fprintf (fp, "  Maxtemp %d\n", (int) omp->maxtemp);
  fprintf (fp, "  Tempcnt %d\n", (int) omp->tempcnt);
  fprintf (fp, "  Hold %d\n", (int) omp->hold);
  PrintABool (fp, "  Reaping", omp->reaping);
  PrintABool (fp, "  Is_write_locked", omp->is_write_locked);
  fprintf (fp, "\n");
  sel = ObjMgrGetSelected ();
  if (sel != NULL) {
    omdp = ObjMgrGetData (sel->entityID);
    ReportOnEntity (omdp, omp, TRUE, sel->itemID, sel->itemtype, 0, fp);
  } else {
    num = omp->currobj;
    for (j = 0, omdpp = omp->datalist; j < num && omdpp != NULL; j++, omdpp++) {
      omdp = *omdpp;
      if (omdp->parentptr == NULL) {
        ReportOnEntity (omdp, omp, FALSE, 0, 0, j + 1, fp);
      }
    }
    for (j = 0, omdpp = omp->datalist; j < num && omdpp != NULL; j++, omdpp++) {
      omdp = *omdpp;
      if (omdp->parentptr != NULL) {
        ReportOnEntity (omdp, omp, FALSE, 0, 0, j + 1, fp);
      }
    }
  }
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Object Manager Report");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}

static void ConvertGiToAccn (SeqIdPtr sip)

{
  Int4      gi;
  SeqIdPtr  newsip;

  if (sip == NULL) return;
  if (sip->choice != SEQID_GI) return;
  gi = sip->data.intvalue;
  newsip = GetSeqIdForGI (gi);
  if (newsip == NULL) return;
  if (newsip->choice == SEQID_GIBBSQ ||
      newsip->choice == SEQID_GIBBMT ||
      newsip->choice == SEQID_GI) {
    SeqIdFree (newsip);
    return;
  }
  sip->choice = newsip->choice;
  sip->data.ptrvalue = newsip->data.ptrvalue;
  newsip->choice = SEQID_NOT_SET;
  newsip->data.ptrvalue = NULL;
  SeqIdFree (newsip);
}

static Boolean GiToAccnAlignCallback (GatherContextPtr gcp)

{
  SeqAlignPtr   align;
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  SeqIdPtr      sip;
  StdSegPtr     ssp;
  SeqLocPtr     tloc;

  if (gcp == NULL) return TRUE;
  switch (gcp->thistype) {
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      align = (SeqAlignPtr) gcp->thisitem;
      sip = NULL;
      if (align->segtype == 1) {
        ddp = (DenseDiagPtr) align->segs;
        if (ddp != NULL) {
          for (sip = ddp->id; sip != NULL; sip = sip->next) {
            ConvertGiToAccn (sip);
          }
        }
      } else if (align->segtype == 2) {
        dsp = (DenseSegPtr) align->segs;
        if (dsp != NULL) {
          for (sip = dsp->ids; sip != NULL; sip = sip->next) {
            ConvertGiToAccn (sip);
          }
        }
      } else if (align->segtype == 3) {
        ssp = (StdSegPtr) align->segs;
        if (ssp != NULL) {
          for (tloc = ssp->loc; tloc != NULL; tloc = tloc->next) {
            sip = SeqLocId (tloc);
            ConvertGiToAccn (sip);
          }
        }
      }
      break;
    default :
      break;
  }
  return TRUE;
}

static Int2 LIBCALLBACK AlignGiToAccnProc (Pointer data)

{
  MsgAnswer         ans;
  GatherScope       gs;
  OMProcControlPtr  ompcp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ans = Message (MSG_OKC, "Are you sure you want to convert alignment GIs to accessions?");
  if (ans == ANS_CANCEL) return OM_MSG_RET_DONE;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQALIGN] = FALSE;
  gs.ignore[OBJ_SEQHIST_ALIGN] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (ompcp->input_entityID, NULL, GiToAccnAlignCallback, &gs);
  return OM_MSG_RET_DONE;
}

static void CacheByAccn (SeqIdPtr sip, CharPtr directory)

{
  AsnIoPtr     aip;
  BioseqPtr    bsp;
  Char         buf [45];
  Uint2        entityID;
  Int4         gi;
  SeqIdPtr     newsip = NULL;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;

  if (sip == NULL || sip->choice == SEQID_LOCAL) return;
  if (sip->choice == SEQID_GI) {
    gi = sip->data.intvalue;
    newsip = GetSeqIdForGI (gi);
    sip = newsip;
  }
  if (sip == NULL) return;

  SeqIdWrite(sip, buf, PRINTID_TEXTID_ACCESSION, 40);
  StringCat (buf, ".ent");
  StringCpy (path, directory);
  FileBuildPath (path, NULL, buf);
  if (FileLength (path) > 0) return;

  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    entityID = ObjMgrGetEntityIDForChoice (sep);
    sep = GetTopSeqEntryForEntityID (entityID);
    if (sep != NULL) {
      aip = AsnIoOpen (path, "w");
      if (aip != NULL) {
        SeqEntryAsnWrite (sep, aip, NULL);
        AsnIoClose (aip);
      }
    }
  }
  BioseqUnlock (bsp);

  if (newsip != NULL) {
    SeqIdFree (newsip);
  }
}

static Boolean CacheAccnsCallback (GatherContextPtr gcp)

{
  SeqAlignPtr   align;
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  CharPtr       path;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  StdSegPtr     ssp;
  SeqLocPtr     tloc;

  if (gcp == NULL) return TRUE;
  path = (CharPtr) gcp->userdata;
  switch (gcp->thistype) {
    case OBJ_BIOSEQ_SEG :
      slp = (SeqLocPtr) gcp->thisitem;
      if (slp != NULL) {
        sip = SeqLocId (slp);
        CacheByAccn (sip, path);
      }
      break;
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      align = (SeqAlignPtr) gcp->thisitem;
      sip = NULL;
      if (align->segtype == 1) {
        ddp = (DenseDiagPtr) align->segs;
        if (ddp != NULL) {
          for (sip = ddp->id; sip != NULL; sip = sip->next) {
            CacheByAccn (sip, path);
          }
        }
      } else if (align->segtype == 2) {
        dsp = (DenseSegPtr) align->segs;
        if (dsp != NULL) {
          for (sip = dsp->ids; sip != NULL; sip = sip->next) {
            CacheByAccn (sip, path);
          }
        }
      } else if (align->segtype == 3) {
        ssp = (StdSegPtr) align->segs;
        if (ssp != NULL) {
          for (tloc = ssp->loc; tloc != NULL; tloc = tloc->next) {
            sip = SeqLocId (tloc);
            CacheByAccn (sip, path);
          }
        }
      }
      break;
    default :
      break;
  }
  return TRUE;
}

static Int2 LIBCALLBACK CacheAccnsToDisk (Pointer data)

{
  GatherScope       gs;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  CharPtr           ptr;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  if (! GetOutputFileName (path, sizeof (path), NULL)) return OM_MSG_RET_DONE;
  ptr = StringRChr (path, DIRDELIMCHR);
  if (ptr != NULL) {
    *ptr = '\0';
  }
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQALIGN] = FALSE;
  gs.ignore[OBJ_SEQHIST_ALIGN] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (ompcp->input_entityID, (Pointer) path, CacheAccnsCallback, &gs);
  return OM_MSG_RET_DONE;
}

static Boolean AddFeatToVnp (GatherContextPtr gcp)

{
  SeqFeatPtr       sfp;
  ValNodePtr       vnp;
  ValNodePtr PNTR  vnpp;

  vnpp = (ValNodePtr PNTR) gcp->userdata;
  if (vnpp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || sfp->location == NULL) return TRUE;
  vnp = ValNodeAdd (vnpp);
  if (vnp == NULL) return TRUE;
  vnp->data.ptrvalue = (Pointer) sfp;
  return TRUE;
}

typedef struct protfindlist {
  SeqLocPtr   slp;
  ProtRefPtr  prp;
  Int4        min;
} ProtFindList, PNTR ProtFindPtr;

static Boolean ProtFindFunc (GatherContextPtr gcp)

{
  Int4         diff;
  ProtFindPtr  pfp;
  SeqFeatPtr   sfp;

  if (gcp == NULL) return TRUE;

  pfp = (ProtFindPtr) gcp->userdata;
  if (pfp == NULL) return TRUE;

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT && sfp->data.value.ptrvalue != NULL) {
      diff = SeqLocAinB (pfp->slp, sfp->location);
      if (diff >= 0) {
        if (diff < pfp->min) {
          pfp->min = diff;
          pfp->prp = (ProtRefPtr) sfp->data.value.ptrvalue;
        }
      }
    }
  }

  return TRUE;
}

static ProtRefPtr FindBestProtRef (Uint2 entityID, SeqFeatPtr cds)

{
  GatherScope   gs;
  ProtFindList  pfl;

  if (entityID == 0 || cds == NULL || cds->product == NULL) return NULL;
  pfl.slp = cds->product;
  pfl.prp = NULL;
  pfl.min = INT4_MAX;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (entityID, (Pointer) &pfl, ProtFindFunc, &gs);
  return pfl.prp;
}

extern void MRnaFromCdsProc (Uint2 entityID)

{
  BioseqPtr     bsp;
  SeqFeatPtr    cds;
  ValNodePtr    head;
  Char          mRnaName [128];
  ValNodePtr    name;
  Boolean       noLeft;
  Boolean       noRight;
  ProtRefPtr    prp;
  SeqFeatPtr    rna;
  RnaRefPtr     rrp;
  SelStructPtr  sel;
  SeqEntryPtr   sep;
  ValNodePtr    vnp;

  if (entityID == 0) return;
  head = NULL;
  for (sel = ObjMgrGetSelected (); sel != NULL; sel = sel->next) {
    GatherItem (sel->entityID, sel->itemID, sel->itemtype,
                (Pointer) &head, AddFeatToVnp);
  }
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    cds = (SeqFeatPtr) vnp->data.ptrvalue;
    if (cds != NULL) {
      rrp = RnaRefNew ();
      if (rrp != NULL) {
        rrp->type = 2;
        mRnaName [0] = '\0';
        prp = FindBestProtRef (entityID, cds);
        if (prp != NULL) {
          name = prp->name;
          if (name != NULL) {
            StringNCpy_0 (mRnaName, (CharPtr) name->data.ptrvalue, sizeof (mRnaName));
          }
          if (StringHasNoText (mRnaName)) {
            StringNCpy_0 (mRnaName, prp->desc, sizeof (mRnaName));
          }
        }
        if (! StringHasNoText (mRnaName)) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave (mRnaName);
        }
        rna = SeqFeatNew ();
        if (rna != NULL) {
          rna->data.choice = SEQFEAT_RNA;
          rna->data.value.ptrvalue = (Pointer) rrp;
          rna->location = AsnIoMemCopy ((Pointer) cds->location,
                                        (AsnReadFunc) SeqLocAsnRead,
                                        (AsnWriteFunc) SeqLocAsnWrite);
          CheckSeqLocForPartial (rna->location, &noLeft, &noRight);
          rna->partial = (rna->partial || noLeft || noRight);
          bsp = GetBioseqGivenSeqLoc (rna->location, entityID);
          if (bsp != NULL) {
            sep = SeqMgrGetSeqEntryForData (bsp);
            if (sep != NULL) {
              CreateNewFeature (sep, NULL, SEQFEAT_RNA, rna);
            } else {
              rna->next = cds->next;
              cds->next = rna;
            }
          } else {
            rna->next = cds->next;
            cds->next = rna;
          }
        }
      }
    }
  }
  ValNodeFree (head);
  ObjMgrSetDirtyFlag (entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
}

static Int2 LIBCALLBACK MRnaFromCdsFunc (Pointer data)

{
  OMProcControlPtr  ompcp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  MRnaFromCdsProc (ompcp->input_entityID);
  return OM_MSG_RET_DONE;
}

typedef struct cdsfindlist {
  SeqLocPtr   loc;
  SeqLocPtr   prod;
  SeqFeatPtr  cds;
  Int4        min;
} CdsFindList, PNTR CdsFindPtr;

static Boolean CdsFindFunc (GatherContextPtr gcp)

{
  CdsFindPtr      cfp;
  Int4            diff;
  SeqFeatPtr      sfp;

  if (gcp == NULL) return TRUE;

  cfp = (CdsFindPtr) gcp->userdata;
  if (cfp == NULL) return TRUE;

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION && sfp->data.value.ptrvalue != NULL) {
      if (cfp->loc != NULL) {
        diff = SeqLocAinB (cfp->loc, sfp->location);
        if (diff >= 0) {
          if (diff < cfp->min) {
            cfp->min = diff;
            cfp->cds = sfp;
          }
        }
      } else if (cfp->prod != NULL) {
        diff = SeqLocAinB (cfp->prod, sfp->product);
        if (diff >= 0) {
          if (diff < cfp->min) {
            cfp->min = diff;
            cfp->cds = sfp;
          }
        }
      }
    }
  }

  return TRUE;
}

extern SeqFeatPtr FindBestCds (Uint2 entityID, SeqLocPtr loc, SeqLocPtr prod, SeqEntryPtr scope)

{
  CdsFindList  cfl;
  GatherScope  gs;

  if (entityID == 0) return NULL;
  cfl.loc = loc;
  cfl.prod = prod;
  cfl.cds = NULL;
  cfl.min = INT4_MAX;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.scope = scope;
  GatherEntity (entityID, (Pointer) &cfl, CdsFindFunc, &gs);
  return cfl.cds;
}

static Int2 LIBCALLBACK MapToProtFunc (Pointer data)

{
  SeqFeatPtr        cds;
  OMProcControlPtr  ompcp;
  SeqEntryPtr       scope;
  SeqFeatPtr        sfp;
  SeqLocPtr         slp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_itemtype == 0)
    return OM_MSG_RET_ERROR;

  switch (ompcp->input_itemtype)
  {
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) ompcp->input_data;
      break;
    default:
      return OM_MSG_RET_ERROR;
  }

  scope = GetBestTopParentForItemID (ompcp->input_entityID, ompcp->input_itemID,
                                     ompcp->input_itemtype);
  cds = FindBestCds (ompcp->input_entityID, sfp->location, NULL, scope);
  if (cds != NULL) {
    slp = dnaLoc_to_aaLoc (cds, sfp->location, TRUE, NULL, FALSE);
    if (slp != NULL) {
      ObjMgrRegister (OBJ_SEQLOC, (Pointer) slp);
    }
    return OM_MSG_RET_DONE;
  } else {
    Message (MSG_ERROR, "Unable to find CDS for peptide feature");
    return OM_MSG_RET_ERROR;
  }
}

static Int2 LIBCALLBACK MapToNucFunc (Pointer data)

{
  SeqFeatPtr        cds;
  OMProcControlPtr  ompcp;
  SeqEntryPtr       scope;
  SeqFeatPtr        sfp;
  SeqLocPtr         slp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_itemtype == 0)
    return OM_MSG_RET_ERROR;

  switch (ompcp->input_itemtype)
  {
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) ompcp->input_data;
      break;
    default:
      return OM_MSG_RET_ERROR;
  }

  scope = GetBestTopParentForItemID (ompcp->input_entityID, ompcp->input_itemID,
                                     ompcp->input_itemtype);
  cds = FindBestCds (ompcp->input_entityID, NULL, sfp->location, scope);
  if (cds != NULL) {
    slp = aaLoc_to_dnaLoc (cds, sfp->location);
    if (slp != NULL) {
      ObjMgrRegister (OBJ_SEQLOC, (Pointer) slp);
    }
    return OM_MSG_RET_DONE;
  } else {
    Message (MSG_ERROR, "Unable to find CDS for peptide feature");
    return OM_MSG_RET_ERROR;
  }
}

static void RevCompFeats (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  bsp = (BioseqPtr) mydata;
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        sip = SeqLocId (sfp->location);
        if (sip != NULL) {
          if (SeqIdIn (sip, bsp->id)) {
            slp = SeqLocCopyRegion (sip, sfp->location, bsp, 0,
                                    bsp->length - 1, Seq_strand_minus, FALSE);
            sfp->location = SeqLocFree (sfp->location);
            sfp->location = slp;
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

static void AddBspToVnpCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr        bsp;
  ValNodePtr PNTR  vnpp;

  vnpp = (ValNodePtr PNTR) mydata;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL && ISA_na (bsp->mol)) {
      ValNodeAddPointer (vnpp, 0, (Pointer) bsp);
    }
  }
}

static Boolean AddBspToVnp (GatherContextPtr gcp)

{
  BioseqPtr        bsp;
  BioseqSetPtr     bssp;
  SeqEntryPtr      sep;
  ValNodePtr PNTR  vnpp;

  vnpp = (ValNodePtr PNTR) gcp->userdata;
  if (vnpp == NULL) return TRUE;
  bsp = NULL;
  if (gcp->thistype == OBJ_BIOSEQ) {
    bsp = (BioseqPtr) gcp->thisitem;
  } else if (gcp->thistype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) gcp->thisitem;
    if (bssp != NULL) {
      sep = SeqMgrGetSeqEntryForData (bssp);
      if (sep != NULL) {
        SeqEntryExplore (sep, (Pointer) vnpp, AddBspToVnpCallback);
      }
    }
    return TRUE;
  } else return TRUE;
  if (bsp == NULL) return TRUE;
  ValNodeAddPointer (vnpp, 0, (Pointer) bsp);
  return TRUE;
}

typedef Boolean (LIBCALL *BioseqFunc) (BioseqPtr);

static void ProcessMultipleBioseqFunctions (OMProcControlPtr ompcp, BioseqFunc func,
                                            Boolean revCompFeats)

{
  BioseqPtr     bsp;
  ValNodePtr    head;
  SelStructPtr  sel;
  SeqEntryPtr   sep;
  ValNodePtr    vnp;

  if (ompcp == NULL || ompcp->input_entityID == 0 || func == NULL) return;
  head = NULL;
  for (sel = ObjMgrGetSelected (); sel != NULL; sel = sel->next) {
    GatherItem (sel->entityID, sel->itemID, sel->itemtype,
                (Pointer) &head, AddBspToVnp);
  }
  if (head == NULL) return;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    func (bsp);
    if (revCompFeats) {
      if (bsp->repr == Seq_repr_raw || bsp->repr == Seq_repr_const) {
        sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
        if (sep != NULL) {
          SeqEntryExplore (sep, (Pointer) bsp, RevCompFeats);
        }
      }
    }
  }
  ValNodeFree (head);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
}

static Int2 LIBCALLBACK RevCompFuncFeat (Pointer data)

{
  /*
  BioseqPtr         bsp;
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    default :
      return OM_MSG_RET_ERROR;
      break;
  }
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  BioseqRevComp (bsp);
  if (bsp->repr == Seq_repr_raw || bsp->repr == Seq_repr_const) {
    sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
    if (sep != NULL) {
      SeqEntryExplore (sep, (Pointer) bsp, RevCompFeats);
    }
  }
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  */
  ProcessMultipleBioseqFunctions ((OMProcControlPtr) data, BioseqRevComp, TRUE);
  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK RevCompFunc (Pointer data)

{
  /*
  BioseqPtr         bsp;
  OMProcControlPtr  ompcp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    default :
      return OM_MSG_RET_ERROR;
      break;
  }
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  BioseqRevComp (bsp);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  */
  ProcessMultipleBioseqFunctions ((OMProcControlPtr) data, BioseqRevComp, FALSE);
  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK RevFunc (Pointer data)

{
  /*
  BioseqPtr         bsp;
  OMProcControlPtr  ompcp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    default :
      return OM_MSG_RET_ERROR;
      break;
  }
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  BioseqReverse (bsp);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  */
  ProcessMultipleBioseqFunctions ((OMProcControlPtr) data, BioseqReverse, FALSE);
  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK CompFunc (Pointer data)

{
  /*
  BioseqPtr         bsp;
  OMProcControlPtr  ompcp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    default :
      return OM_MSG_RET_ERROR;
      break;
  }
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  BioseqComplement (bsp);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  */
  ProcessMultipleBioseqFunctions ((OMProcControlPtr) data, BioseqComplement, FALSE);
  return OM_MSG_RET_DONE;
}

typedef struct orphandata {
  ValNodePtr  bspfromcds;
  ValNodePtr  bspinentry;
} OrphanData, PNTR OrphanDataPtr;

static void GetCdsProducts (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  OrphanDataPtr  odp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  SeqIdPtr       sip;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  odp = (OrphanDataPtr) mydata;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    if (ISA_aa (bsp->mol)) {
      ValNodeAddPointer (&(odp->bspinentry), 0, (Pointer) bsp);
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_CDREGION) {
          if (sfp->product != NULL) {
            sip = SeqLocId (sfp->product);
            if (sip != NULL) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) {
                ValNodeAddPointer (&(odp->bspfromcds), 0, (Pointer) bsp);
              }
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

static int LIBCALLBACK SortByVnpDataPtrvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->data.ptrvalue > vnp2->data.ptrvalue) {
        return 1;
      } else if (vnp1->data.ptrvalue < vnp2->data.ptrvalue) {
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

static void RemoveProteinBioseq (SeqEntryPtr sep, BioseqPtr bsp)

{
  BioseqSetPtr      bssp;
  SeqEntryPtr       next;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  Uint2             parenttype;
  Pointer           parentptr;
  SeqEntryPtr PNTR  prev;
  SeqEntryPtr       top;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    top = sep;
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL) {
      prev = &(bssp->seq_set);
      sep = bssp->seq_set;
      while (sep != NULL) {
        next = sep->next;
        if (IS_Bioseq_set (sep)) {
          RemoveProteinBioseq (sep, bsp);
          sep = next;
        } else if (IS_Bioseq (sep) && sep->data.ptrvalue == (Pointer) bsp) {
          *(prev) = next;
          sep->next = NULL;
          SaveSeqEntryObjMgrData (top, &omdptop, &omdata);
          GetSeqEntryParent (top, &parentptr, &parenttype);
          SeqEntryFree (sep);
          SeqMgrLinkSeqEntry (top, parenttype, parentptr);
          RestoreSeqEntryObjMgrData (top, omdptop, &omdata);
          sep = next;
        } else {
          prev = &(sep->next);
          sep = next;
        }
      }
    }
  }
}

static void RemoveOrphanProteins (Uint2 entityID, SeqEntryPtr sep)

{
  Int2        doit;
  OrphanData  od;
  ValNodePtr  vnp;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (sep == NULL || entityID == 0) return;
  od.bspfromcds = NULL;
  od.bspinentry = NULL;
  SeqEntryExplore (sep, (Pointer) &od, GetCdsProducts);
  if (od.bspinentry != NULL && od.bspfromcds != NULL) {
    od.bspinentry = SortValNode (od.bspinentry, SortByVnpDataPtrvalue);
    od.bspfromcds = SortValNode (od.bspfromcds, SortByVnpDataPtrvalue);
    vnp1 = od.bspfromcds;
    vnp2 = od.bspinentry;
    while (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->data.ptrvalue < vnp2->data.ptrvalue) {
        vnp1 = vnp1->next;
      } else if (vnp2->data.ptrvalue < vnp1->data.ptrvalue) {
        vnp2 = vnp2->next;
      } else {
        vnp2->data.ptrvalue = NULL;
        vnp1 = vnp1->next;
        vnp2 = vnp2->next;
      }
    }
    doit = 0;
    for (vnp = od.bspinentry; vnp != NULL; vnp = vnp->next) {
      if (vnp->data.ptrvalue != NULL) {
        doit++;
      }
    }
    if (doit > 0) {
      if (Message (MSG_YN, "Do you want to remove %d orphaned proteins?", (int) doit) == ANS_YES) {
        vnp = od.bspinentry;
        while (vnp != NULL) {
          if (vnp->data.ptrvalue != NULL) {
            RemoveProteinBioseq (sep, (BioseqPtr) vnp->data.ptrvalue);
          }
          vnp = vnp->next;
        }
      }
    }
  }
  ValNodeFree (od.bspinentry);
  ValNodeFree (od.bspfromcds);
}

static Int2 LIBCALLBACK NormalizeNucProts (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  RemoveOrphanProteins (ompcp->input_entityID, sep);
  RenormalizeNucProtSets (sep, TRUE);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static void RemovePopPhyMutSet (SeqEntryPtr sep, Boolean forceGBClass, Boolean forceAll)

{
  SeqAnnotPtr   annot;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    descr;
  SeqAnnotPtr   sap;
  SeqEntryPtr   seqentry;

  if (sep == NULL || IS_Bioseq (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (forceAll) {
  } else if (forceGBClass) {
    if ((bssp->_class < 13 || bssp->_class > 15) && bssp->_class != 7) return;
  } else {
    if (bssp->_class < 13 || bssp->_class > 15) return;
  }
  seqentry = bssp->seq_set;
  if (seqentry == NULL) return;
  if (! forceAll) {
    if (! forceGBClass) {
      if (seqentry->next != NULL) return;
    }
  }
  descr = bssp->descr;
  bssp->descr = NULL;
  annot = bssp->annot;
  bssp->annot = NULL;
  sep->choice = seqentry->choice;
  sep->data.ptrvalue = seqentry->data.ptrvalue;
  sep->next = seqentry->next;
  seqentry->data.ptrvalue = NULL;
  seqentry->next = NULL;
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
}

static void InsetNewSet (BioseqSetPtr bssp, Uint1 _class)

{
  BioseqSetPtr  popphymut;
  SeqEntryPtr   tmp;

  if (bssp == NULL) return;
  popphymut = BioseqSetNew ();
  if (popphymut == NULL) return;
  popphymut->_class = _class;
  popphymut->seq_set = bssp->seq_set;
  tmp = SeqEntryNew ();
  if (tmp == NULL) return;
  tmp->choice = 2;
  tmp->data.ptrvalue = (Pointer) popphymut;
  bssp->seq_set = tmp;
}

static Int2 LIBCALLBACK RemoveExtraneousSets (Pointer data)

{
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  Uint1             _class;
  ValNodePtr        descr1 = NULL;
  ValNodePtr        descr2 = NULL;
  BioseqSetPtr      first = NULL;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  OMProcControlPtr  ompcp;
  Uint2             parenttype;
  Pointer           parentptr;
  SeqEntryPtr       sep;
  SeqEntryPtr       tmp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == 7) {
      descr1 = bssp->descr;
      bssp->descr = NULL;
      _class = 0;
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        if (IS_Bioseq (tmp)) return OM_MSG_RET_DONE;
        if (_class == 0 && IS_Bioseq_set (tmp)) {
          first = (BioseqSetPtr) tmp->data.ptrvalue;
          if (first != NULL) {
            _class = first->_class;
          }
        }
      }
      if (first != NULL && bssp->seq_set != NULL && bssp->seq_set->next == NULL) {
        descr2 = first->descr;
        first->descr = NULL;
      }
      SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
      GetSeqEntryParent (sep, &parentptr, &parenttype);
      if (_class >= 13 && _class <= 15) {
        for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
          RemovePopPhyMutSet (tmp, FALSE, FALSE);
        }
        InsetNewSet (bssp, _class);
      }
      if (_class == 7) {
        tmp = bssp->seq_set;
        if (tmp != NULL && IS_Bioseq_set (tmp)) {
          RemovePopPhyMutSet (tmp, TRUE, FALSE);
        }
      }
      if (IS_Bioseq (sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        ValNodeLink (&(bsp->descr), descr1);
        ValNodeLink (&(bsp->descr), descr2);
      } else if (IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        ValNodeLink (&(bssp->descr), descr1);
        ValNodeLink (&(bssp->descr), descr2);
      }
      SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
      RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
      ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
    }
  }
  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK RepairMessedUpRecord (Pointer data)

{
  BioseqSetPtr      bssp;
  ValNodePtr        descr1 = NULL;
  ValNodePtr        descr2 = NULL;
  BioseqSetPtr      first = NULL;
  SeqEntryPtr       last;
  SeqEntryPtr       next;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  OMProcControlPtr  ompcp;
  Uint2             parenttype;
  Pointer           parentptr;
  SeqEntryPtr       sep;
  SeqEntryPtr       tmp;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  if (! IS_Bioseq_set (sep)) return OM_MSG_RET_ERROR;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL || bssp->_class != 7) return OM_MSG_RET_ERROR;
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);
  tmp = bssp->seq_set;
  while (tmp != NULL) {
    next = tmp->next;
    tmp->next = NULL;
    if (IS_Bioseq_set (tmp)) {
      RemovePopPhyMutSet (tmp, FALSE, TRUE);
    }
    last = bssp->seq_set;
    while (last->next != NULL) {
      last = last->next;
    }
    last->next = next;
    tmp = next;
  }
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

typedef struct featsints {
  SeqFeatPtr  sfp;
  Int4        leftend;
  Int4        rightend;
  Uint1       strand;
} FeatsIntsData, PNTR FeatsIntsPtr;

static void CommonGatherRBSorGene (ValNodePtr PNTR vnpp, SeqFeatPtr sfp, Uint2 entityID)

{
  BioseqPtr     bsp;
  FeatsIntsPtr  fip;
  GeneRefPtr    grp;
  Int4          start;
  Int4          stop;
  Char          str [128];
  ValNodePtr    vnp;

  if (vnpp == NULL || sfp == NULL || entityID == 0) return;
  bsp = GetBioseqGivenSeqLoc (sfp->location, entityID);
  if (bsp == NULL) return;
  fip = MemNew (sizeof (FeatsIntsData));
  if (fip == NULL) return;
  fip->sfp = sfp;
  fip->leftend = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_LEFT_END);
  fip->rightend = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_RIGHT_END);
  fip->strand = SeqLocStrand (sfp->location);
  if (sfp->data.choice == SEQFEAT_GENE) {
    start = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_START);
    stop = GetOffsetInBioseq (sfp->location, bsp, SEQLOC_STOP);
    if ((start > stop && fip->strand == Seq_strand_plus) ||
        (start < stop && fip->strand == Seq_strand_minus)) {
      str [0] = '\0';
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        StringNCpy_0 (str, grp->locus, sizeof (str));
        if (StringHasNoText (str)) {
          vnp = grp->syn;
          if (vnp != NULL) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
          }
        }
        if (StringHasNoText (str)) {
          StringNCpy_0 (str, grp->desc, sizeof (str));
        }
        if (StringHasNoText (str)) {
          StringNCpy_0 (str, "?", sizeof (str));
        }
      }
      Message (MSG_POST, "Ignoring gene '%s' going across origin, please fix manually", str);
      return;
    }
  }
  ValNodeAddPointer (vnpp, 0, (Pointer) fip);
}

static Boolean GetRBSGatherFunc (GatherContextPtr gcp)

{
  SeqFeatPtr  sfp;
  Uint1       type;
  ValNodePtr  PNTR vnpp;

  if (gcp == NULL || gcp->thisitem == NULL || gcp->userdata == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  switch (sfp->data.choice) {
    case SEQFEAT_IMP :
      type = FindFeatDefType (sfp);
      if (type != FEATDEF_RBS) return TRUE;
      vnpp = (ValNodePtr PNTR) gcp->userdata;
      CommonGatherRBSorGene (vnpp, sfp, gcp->entityID);
      break;
    default :
      break;
  }
  return TRUE;
}

static Boolean GetGeneGatherFunc (GatherContextPtr gcp)

{
  SeqFeatPtr  sfp;
  ValNodePtr  PNTR vnpp;

  if (gcp == NULL || gcp->thisitem == NULL || gcp->userdata == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      vnpp = (ValNodePtr PNTR) gcp->userdata;
      CommonGatherRBSorGene (vnpp, sfp, gcp->entityID);
      break;
    default :
      break;
  }
  return TRUE;
}

static Boolean GetNonGeneGatherFunc (GatherContextPtr gcp)

{
  SeqFeatPtr  sfp;
  ValNodePtr  PNTR vnpp;

  if (gcp == NULL || gcp->thisitem == NULL || gcp->userdata == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
    case SEQFEAT_PROT :
      break;
    default :
      vnpp = (ValNodePtr PNTR) gcp->userdata;
      CommonGatherRBSorGene (vnpp, sfp, gcp->entityID);
      break;
  }
  return TRUE;
}

static int LIBCALLBACK SortFeatByLocation (VoidPtr ptr1, VoidPtr ptr2)

{
  FeatsIntsPtr  fip1;
  FeatsIntsPtr  fip2;
  ValNodePtr    vnp1;
  ValNodePtr    vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      fip1 = (FeatsIntsPtr) vnp1->data.ptrvalue;
      fip2 = (FeatsIntsPtr) vnp2->data.ptrvalue;
      if (fip1 != NULL && fip2 != NULL) {
        if (fip1->leftend > fip2->leftend) {
          return 1;
        } else if (fip1->leftend < fip2->leftend) {
          return -1;
        } else if (fip1->rightend > fip2->rightend) {
          return 1;
        } else if (fip1->rightend < fip2->rightend) {
          return -1;
        } else {
          return 0;
        }
      }
    }
  }
  return 0;
}

static void ExtendGeneByRbs (SeqFeatPtr sfpr, SeqFeatPtr sfpg, Uint2 entityID)

{
  BioseqPtr       bsp;
  SeqFeatXrefPtr  curr;
  SeqFeatXrefPtr  PNTR last;
  SeqFeatXrefPtr  next;
  SeqLocPtr       slp;

  if (sfpr == NULL || sfpg == NULL || sfpg->data.choice != SEQFEAT_GENE) return;
  bsp = GetBioseqGivenSeqLoc (sfpg->location, entityID);
  if (bsp == NULL) return;
  slp = SeqLocMerge (bsp, sfpg->location, sfpr->location, TRUE, TRUE, FALSE);
  sfpg->location = SeqLocFree (sfpg->location);
  sfpg->location = slp;
  last = (SeqFeatXrefPtr PNTR) &(sfpr->xref);
  curr = sfpr->xref;
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

static void MakeGeneXref (SeqFeatPtr sfpr, SeqFeatPtr sfpg)

{
  SeqFeatXrefPtr  curr;
  GeneRefPtr      grp;
  SeqFeatXrefPtr  PNTR last;
  CharPtr         locus;
  SeqFeatXrefPtr  next;
  SeqFeatXrefPtr  xref;

  if (sfpr == NULL || sfpg == NULL || sfpg->data.choice != SEQFEAT_GENE) return;
  grp = (GeneRefPtr) sfpg->data.value.ptrvalue;
  if (grp == NULL) return;
  locus = grp->locus;
  if (StringHasNoText (locus)) return;
  last = (SeqFeatXrefPtr PNTR) &(sfpr->xref);
  curr = sfpr->xref;
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
  grp = GeneRefNew ();
  if (grp == NULL) return;
  grp->locus = StringSave (locus);
  TrimSpacesAroundString (grp->locus);
  xref = SeqFeatXrefNew ();
  sfpr->xref = xref;
  if (xref != NULL) {
    xref->data.choice = SEQFEAT_GENE;
    xref->data.value.ptrvalue = (Pointer) grp;
  }
}

static Int2 FindNextGeneFeat (ValNodePtr PNTR genelist, Int4 max, FeatsIntsPtr fipr)

{
  Int4          compare;
  FeatsIntsPtr  fipg;
  Int4          left;
  Int4          mid;
  Int4          right;
  ValNodePtr    vnpg;

  if (genelist == NULL || fipr == NULL || max < 1) return 0;
  mid = 0;
  left = 1;
  right = max;
  while (left <= right) {
    mid = (left + right) / 2;
    vnpg = genelist [mid - 1];
    fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
    if (fipr->strand == Seq_strand_minus) {
      compare = fipg->rightend - fipr->leftend;
    } else {
      compare = fipr->rightend - fipg->leftend;
    }
    if (compare <= 0) {
      right = mid - 1;
    }
    if (compare >= 0) {
      left = mid + 1;
    }
  }
  if (left <= right + 1) {
    if (fipr->strand == Seq_strand_minus) {
      mid += 2;
      mid = MIN (mid, max);
      vnpg = genelist [mid - 1];
      fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
      while (fipg->rightend > fipr->leftend && mid > 1) {
        mid--;
        vnpg = genelist [mid - 1];
        fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
      }
    } else {
      mid -= 2;
      mid = MAX (mid, 1);
      vnpg = genelist [mid - 1];
      fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
      while (fipr->rightend > fipg->leftend && mid < max) {
        mid++;
        vnpg = genelist [mid - 1];
        fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
      }
    }
  }
  return mid;
}

static Boolean  rbsDirty;
static Boolean  rbsChoiceAsked;
static Boolean  rbsExtend;

static void FixupRBSGenes (Uint2 entityID, SeqEntryPtr sep)

{
  MsgAnswer     ans;
  BioseqSetPtr  bssp;
  FeatsIntsPtr  fipg;
  FeatsIntsPtr  fipr;
  ValNodePtr    PNTR genelist;
  ValNodePtr    genes;
  GatherScope   gs;
  Int2          i;
  Int4          max;
  ValNodePtr    rbss;
  ValNodePtr    tmp;
  ValNodePtr    vnpg;
  ValNodePtr    vnpr;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        FixupRBSGenes (entityID, sep);
      }
      return;
    }
  }

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  rbss = NULL;
  GatherEntity (entityID, (Pointer) (&rbss), GetRBSGatherFunc, &gs);
  rbss = SortValNode (rbss, SortFeatByLocation);
  genes = NULL;
  GatherEntity (entityID, (Pointer) (&genes), GetGeneGatherFunc, &gs);
  genes = SortValNode (genes, SortFeatByLocation);

  for (vnpg = genes, i = 0; vnpg != NULL; vnpg = vnpg->next) {
    i++;
  }
  genelist = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * (i + 1));
  if (genelist != NULL) {
    for (tmp = genes, i = 0; tmp != NULL; tmp = tmp->next) {
      genelist [i] = tmp;
      i++;
    }
    max = (Int4) i;

    for (vnpr = rbss; vnpr != NULL; vnpr = vnpr->next) {
      fipr = (FeatsIntsPtr) vnpr->data.ptrvalue;
      i = FindNextGeneFeat (genelist, max, fipr);
      if (i > 0) {
        vnpg = genelist [i - 1];
        fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
        rbsDirty = TRUE;
        if (! rbsChoiceAsked) {
          rbsChoiceAsked = TRUE;
          ans = Message (MSG_YN, "Extend Gene Feature?");
          rbsExtend = (Boolean) (ans == ANS_YES);
        }
        if (rbsExtend) {
          ExtendGeneByRbs (fipr->sfp, fipg->sfp, entityID);
        } else {
          MakeGeneXref (fipr->sfp, fipg->sfp);
        }
      }
    }
  }

  MemFree (genelist);
  ValNodeFreeData (rbss);
  ValNodeFreeData (genes);
}

static Int2 LIBCALLBACK FixupRBS (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_entityID == 0) {
    Message (MSG_ERROR, "Please select a Bioseq or BioseqSet");
    return OM_MSG_RET_ERROR;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) {
    Message (MSG_ERROR, "Please select a Bioseq or BioseqSet");
    return OM_MSG_RET_ERROR;
  }
  rbsDirty = FALSE;
  rbsChoiceAsked = FALSE;
  rbsExtend = FALSE;
  WatchCursor ();
  Update ();
  FixupRBSGenes (ompcp->input_entityID, sep);
  ArrowCursor ();
  Update ();
  if (rbsDirty) {
    ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  }
  return OM_MSG_RET_DONE;
}

static Boolean LIBCALLBACK GetSeqs (BioseqPtr bsp, SeqMgrBioseqContextPtr context)

{
  FILE  *fp;
  Char  str [256];

  fp = (FILE *) context->userdata;
  if (BioseqLabel (bsp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
    fprintf (fp, "  Bioseq item %d %s\n", (int) context->itemID, str);
  }
  return TRUE;
}

static Boolean LIBCALLBACK GetSegs (SeqLocPtr slp, SeqMgrSegmentContextPtr context)

{
  FILE  *fp;
  Char  str [256];

  fp = (FILE *) context->userdata;
  if (SeqLocLabel (slp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
    fprintf (fp, "  SeqLoc item %d %s cumulative %ld\n", (int) context->itemID,
             str, (long) context->cumOffset);
  }
  return TRUE;
}

static Boolean LIBCALLBACK GetDescs (ValNodePtr sdp, SeqMgrDescContextPtr context)

{
  FILE  *fp;
  Char  str [256];

  fp = (FILE *) context->userdata;
  if (SeqDescLabel (sdp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
    fprintf (fp, "  Descriptor item %d index %d %s\n",
            (int) context->itemID, (int) context->index, str);
  }
  return TRUE;
}

static Boolean LIBCALLBACK GetFeats (SeqFeatPtr sfp, SeqMgrFeatContextPtr context)

{
  FILE  *fp;
  Char  str [256];

  fp = (FILE *) context->userdata;
  if (FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
    fprintf (fp, "  Feature item %d index %d %s\n",
             (int) context->itemID, (int) context->index, str);
  }
  return TRUE;
}

static void DoBioseqReport (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr          bsp;
  SeqMgrDescContext  desccontext;
  SeqMgrFeatContext  featcontext;
  FILE               *fp;
  SeqFeatPtr         gene;
  Int2               i;
  Int4Ptr            ivals;
  Int2               numivals;
  ValNodePtr         sdp;
  SeqFeatPtr         sfp;
  Char               str [256];

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;

  fp = (FILE*) mydata;
  if (fp == NULL) return;

  if (BioseqLabel (bsp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
    fprintf (fp, "Bioseq %s\n\n", str);
  }

  if (bsp->repr == Seq_repr_seg) {
    fprintf (fp, " Exploring segments\n\n");
    SeqMgrExploreSegments (bsp, (Pointer) fp, GetSegs);
    fprintf (fp, "\n");
  }

  if (ISA_aa (bsp->mol)) {
    sfp = SeqMgrGetBestProteinFeature (bsp, NULL);
    if (FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
      fprintf (fp, "  Best protein %s", str);
    }
    sfp = SeqMgrGetCDSgivenProduct (bsp, NULL);
    if (FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
      fprintf (fp, ", best cds %s", str);
    }
    fprintf (fp, "\n\n");
  }


  fprintf (fp, " Exploring descriptors\n\n");
  SeqMgrExploreDescriptors (bsp, (Pointer) fp, GetDescs, NULL);
  fprintf (fp, "\n");

  fprintf (fp, " Exploring features\n\n");
  SeqMgrExploreFeatures (bsp, (Pointer) fp, GetFeats, NULL, NULL, NULL);
  fprintf (fp, "\n");


  fprintf (fp, " Collecting descriptors\n\n");
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &desccontext);
  while (sdp != NULL) {
    if (SeqDescLabel (sdp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
      fprintf (fp, "  Descriptor item %d index %d %s\n",
               (int) desccontext.itemID, (int) desccontext.index, str);
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &desccontext);
  }
  fprintf (fp, "\n");

  fprintf (fp, " Collecting features\n\n");
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &featcontext);
  while (sfp != NULL) {
    if (FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
      fprintf (fp, "  Feature item %d index %d %s",
               (int) featcontext.itemID, (int) featcontext.index, str);
    }
    if (featcontext.seqfeattype != SEQFEAT_GENE) {
      gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
      if (gene != NULL) {
        if (FeatDefLabel (gene, str, sizeof (str) - 1, OM_LABEL_BOTH)) {
          fprintf (fp, ", gene %s", str);
        }
      }
    }
    fprintf (fp, "\n");
    ivals = featcontext.ivals;
    numivals = featcontext.numivals * 2;
    if (ivals != NULL && numivals > 0) {
      fprintf (fp, "    (");
      for (i = 0; i < numivals; i += 2) {
        if (i > 0) {
          fprintf (fp, ", ");
        }
        fprintf (fp, "%ld-%ld", (long) ivals [i], (long) ivals [i + 1]);
      }
      fprintf (fp, ")\n");
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &featcontext);
  }
  fprintf (fp, "\n");
}

static Int2 LIBCALLBACK DoBioseqIndexing (Pointer data)

{
  BioseqPtr         bsp;
  Uint2             entityID;
  FILE              *fp;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_entityID == 0) {
    Message (MSG_ERROR, "Please select a Bioseq");
    return OM_MSG_RET_ERROR;
  }
  bsp = (BioseqPtr) ompcp->input_data;
  if (bsp == NULL) {
    Message (MSG_ERROR, "Please select a Bioseq");
    return OM_MSG_RET_ERROR;
  }
  entityID = SeqMgrIndexFeatures (ompcp->input_entityID, NULL);
  if (entityID == 0) {
    Message (MSG_ERROR, "SeqMgrReindexBioseqExtraData failed");
    return OM_MSG_RET_ERROR;
  }
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_DONE;
  TmpNam (path);
  fp = FileOpen (path, "w");

  fprintf (fp, "Exploring bioseqs\n\n");
  SeqMgrExploreBioseqs (entityID, 0, (Pointer) fp, GetSeqs, TRUE, TRUE, TRUE);
  fprintf (fp, "\n");

  SeqEntryExplore (sep, (Pointer) fp, DoBioseqReport);
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Bioseq Index Report");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}

static Int2 FindBestGeneFeat (ValNodePtr PNTR genelist, Int4 max, FeatsIntsPtr fipr)

{
  Int4          compare;
  FeatsIntsPtr  fipg;
  Int4          left;
  Int4          mid;
  Int4          right;
  ValNodePtr    vnpg;

  if (genelist == NULL || fipr == NULL || max < 1) return 0;
  mid = 0;
  left = 1;
  right = max;
  while (left <= right) {
    mid = (left + right) / 2;
    vnpg = genelist [mid - 1];
    fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
    if (fipg->rightend < fipr->leftend) {
      compare =  1;
    } else if (fipr->rightend < fipg->leftend) {
      compare =  -1;
    } else {
      compare =  0;
    }
    if (compare <= 0) {
      right = mid - 1;
    }
    if (compare >= 0) {
      left = mid + 1;
    }
  }
  if (left > right + 1) {
    if (fipg->strand == fipr->strand) {
      if (SeqLocAinB (fipr->sfp->location, fipg->sfp->location)) {
        return mid;
      }
    }
  }
  return 0;
}

static void DoTrimGenesGenes (Uint2 entityID, SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  ValNodePtr    PNTR featlist;
  FeatsIntsPtr  fipf;
  FeatsIntsPtr  fipg;
  FeatsIntsPtr  fipr;
  ValNodePtr    PNTR genelist;
  ValNodePtr    genes;
  GatherScope   gs;
  Int2          i;
  Int4          max;
  ValNodePtr    rbss;
  ValNodePtr    tmp;
  ValNodePtr    vnpf;
  ValNodePtr    vnpg;
  ValNodePtr    vnpr;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        DoTrimGenesGenes (entityID, sep);
      }
      return;
    }
  }

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  rbss = NULL;
  GatherEntity (entityID, (Pointer) (&rbss), GetNonGeneGatherFunc, &gs);
  rbss = SortValNode (rbss, SortFeatByLocation);
  genes = NULL;
  GatherEntity (entityID, (Pointer) (&genes), GetGeneGatherFunc, &gs);
  genes = SortValNode (genes, SortFeatByLocation);

  for (vnpg = genes, i = 0; vnpg != NULL; vnpg = vnpg->next) {
    i++;
  }
  genelist = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * (i + 1));
  featlist = (ValNodePtr PNTR) MemNew (sizeof (ValNodePtr) * (i + 1));
  if (genelist != NULL && featlist != NULL) {
    for (tmp = genes, i = 0; tmp != NULL; tmp = tmp->next) {
      genelist [i] = tmp;
      i++;
    }
    max = (Int4) i;

    for (vnpr = rbss; vnpr != NULL; vnpr = vnpr->next) {
      fipr = (FeatsIntsPtr) vnpr->data.ptrvalue;
      /*
      for (i = 0; i < (Int2) max; i++) {
        vnpg = genelist [i];
        fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
        if (SeqLocAinB (fipr->sfp->location, fipg->sfp->location)) {
          ValNodeAddPointer (&(featlist [i]), 0, (Pointer) fipr);
        }
      }
      */
      i = FindBestGeneFeat (genelist, max, fipr);
      if (i > 0) {
        ValNodeAddPointer (&(featlist [i - 1]), 0, (Pointer) fipr);
      }
    }

    for (i = 0; i < (Int2) max; i++) {
      vnpf = featlist [i];
      vnpg = genelist [i];
      if (vnpf != NULL && vnpg != NULL) {
        fipg = (FeatsIntsPtr) vnpg->data.ptrvalue;
        fipf = (FeatsIntsPtr) vnpf->data.ptrvalue;
        fipg->sfp->location = SeqLocFree (fipg->sfp->location);
        fipg->sfp->location = AsnIoMemCopy ((Pointer) (fipf->sfp->location),
                                            (AsnReadFunc) SeqLocAsnRead,
                                            (AsnWriteFunc) SeqLocAsnWrite);
        for (tmp = vnpf->next; tmp != NULL; tmp = tmp->next) {
          fipf = (FeatsIntsPtr) tmp->data.ptrvalue;
          ExtendGeneByRbs (fipf->sfp, fipg->sfp, entityID);
        }
      }
    }

    for (i = 0; i < (Int2) max; i++) {
      ValNodeFree (featlist [i]);
    }
  }

  MemFree (genelist);
  MemFree (featlist);
  ValNodeFreeData (rbss);
  ValNodeFreeData (genes);
}

static Int2 LIBCALLBACK TrimGenes (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->input_entityID == 0) {
    Message (MSG_ERROR, "Please select a Bioseq or BioseqSet");
    return OM_MSG_RET_ERROR;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) {
    Message (MSG_ERROR, "Please select a Bioseq or BioseqSet");
    return OM_MSG_RET_ERROR;
  }
  WatchCursor ();
  Update ();
  DoTrimGenesGenes (ompcp->input_entityID, sep);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static Int2 LIBCALLBACK SplitIntoSegmentedBioseq (Pointer data)

{
  BioseqPtr         bsp, bsp2;
  BioseqSetPtr      bssp;
  SeqIdPtr          id, id1, id2;
  Int4              left;
  ObjMgrDataPtr     omdptop;
  ObjMgrData        omdata;
  OMProcControlPtr  ompcp;
  Uint2             parenttype;
  Pointer           parentptr;
  BioseqPtr         part1, part2;
  Int4              right;
  SelStructPtr      sel;
  SeqEntryPtr       sep, sep1, sep2;
  SeqIdPtr          sip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  sel = ObjMgrGetSelected ();
  if (sel == NULL || sel->next != NULL) return OM_MSG_RET_ERROR;
  if (sel->entityID != ompcp->input_entityID) return OM_MSG_RET_ERROR;
  if (sel->itemtype != OBJ_BIOSEQ || sel->region == NULL) return OM_MSG_RET_ERROR;

  bsp = GetBioseqGivenSeqLoc (sel->region, ompcp->input_entityID);
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  sep = SeqMgrGetSeqEntryForData (bsp);
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  left = SeqLocStart (sel->region);
  right = SeqLocStop (sel->region);
  if (right != left + 1) {
    Message (MSG_OK, "Please select two adjacent bases");
    return OM_MSG_RET_ERROR;
  }
  id = SeqIdFindBest (bsp->id, 0);

  id1 = MakeNewProteinSeqId (NULL, id);
  part1 = BioseqCopy (id1, id, 0, right - 1, Seq_strand_plus, FALSE);
  sep1 = SeqEntryNew ();
  if (sep1 == NULL) return OM_MSG_RET_ERROR;
  sep1->choice = 1;
  sep1->data.ptrvalue = part1;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) part1, sep1);

  id2 = MakeNewProteinSeqId (NULL, id);
  part2 = BioseqCopy (id2, id, right, bsp->length - 1, Seq_strand_plus, FALSE);
  sep2 = SeqEntryNew ();
  if (sep2 == NULL) return OM_MSG_RET_ERROR;
  sep2->choice = 1;
  sep2->data.ptrvalue = part2;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) part2, sep2);

  AddSeqEntryToSeqEntry (sep1, sep2, TRUE);
  sep2 = FindNucSeqEntry (sep1);
  if (sep2 == NULL || sep2->choice != 1 || sep2->data.ptrvalue == NULL) return OM_MSG_RET_ERROR;
  bsp2 = (BioseqPtr) sep2->data.ptrvalue;

  sip = bsp->id;
  bsp->id = bsp2->id;
  bsp2->id = sip;

  if (sep1 == NULL || sep1->choice != 2 || sep1->data.ptrvalue == NULL) return OM_MSG_RET_ERROR;
  bssp = (BioseqSetPtr) sep1->data.ptrvalue;

  bssp->descr = bsp->descr;
  bsp->descr = NULL;
  bssp->annot = bsp->annot;
  bsp->annot = NULL;

  sep = SeqMgrGetSeqEntryForData (bsp);
  sep->choice = sep1->choice;
  sep->data.ptrvalue = sep1->data.ptrvalue;

  SeqMgrReplaceInBioseqIndex (bsp2);
  SeqMgrReplaceInBioseqIndex (bsp);
  BioseqFree (bsp);

  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);

  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

extern void SetupSequinFilters (void)

{
  Char  str [32];

  if (extraServices) {
    REGISTER_SEQUIN_PROT_TITLES;
    REGISTER_SEQUIN_NUC_TITLES;
    REGISTER_SEQUIN_FEAT_TABLE;
  }

  if (genomeCenter != NULL || indexerVersion) {
    REGISTER_SEQUIN_DUMP_CONTIG;
  }

  if (indexerVersion) {
    REGISTER_SEQUIN_CACHE_ACCN;
    REGISTER_SEQUIN_GI_TO_ACCN;
    REGISTER_REFGENEUSER_DESC_EDIT;
    REGISTER_PROT_IDS_TO_GENE_SYN;
  }

  if (indexerVersion) {
    if (GetAppParam ("SEQUIN", "SETTINGS", "ALLOWAUTOINTRONEXON", NULL, str, sizeof (str))) {
      if (StringICmp (str, "TRUE") == 0) {
        REGISTER_MAKEEXONINTRON;
      }
    }
  }

  REGISTER_BIOSEQ_COMPLEMENT;
  REGISTER_BIOSEQ_REVERSE;
  REGISTER_BIOSEQ_REVCOMP_WITHFEAT;
  REGISTER_BIOSEQ_REVCOMP_NOTFEAT;

  REGISTER_ALUFEAT;
  REGISTER_BIOSEQ_ORF;
  REGISTER_MRNA_FROM_CDS;
  REGISTER_MAP_TO_PROT;
  REGISTER_MAP_TO_NUC;

  REGISTER_MAKESEQALIGN;
  REGISTER_MAKESEQALIGNP;
  REGISTER_NORMSEQALIGN;
  REGISTER_OPENALED;
  REGISTER_OPENALED2;
  REGISTER_OPENALED3;
  REGISTER_OPENALED4;

  if (extraServices) {
    REGISTER_TRIM_GENES;
    REGISTER_FIXUP_RBS;
    REGISTER_INTERVAL_COMBINE;
    REGISTER_REPACKAGE_PARTS;
    REGISTER_UNDOSEGSET;
    REGISTER_ADJUSTMULTISEGSEQ;
    REGISTER_UPDATESEGSET;
  }

  if (indexerVersion) {
    REGISTER_NORMALIZE_NUCPROT;
    REGISTER_REMOVE_EXTRANEOUS;
    REGISTER_REMOVE_MESSEDUP;
    REGISTER_GROUP_EXPLODE;
    REGISTER_SPLIT_BIOSEQ;
    REGISTER_DESKTOP_REPORT;
#ifdef WIN_MAC
    REGISTER_BSP_INDEX;
#endif
  }

  REGISTER_GROUP_FILTERGC;
  REGISTER_GROUP_DUST;
  REGISTER_GROUP_PCC;
  REGISTER_GROUP_HOPP;
  REGISTER_GROUP_KYTE;
  REGISTER_GROUP_MATRIX;
}

extern CharPtr MergeValNodeStrings (ValNodePtr list, Boolean useReturn)

{
  size_t      len;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;


  ptr = NULL;
  if (list != NULL) {
    vnp = list;
    len = 0;
    while (vnp != NULL) {
      if (vnp->data.ptrvalue != NULL) {
        len += StringLen ((CharPtr) vnp->data.ptrvalue) + 1;
      }
      vnp = vnp->next;
    }
    if (len > 0) {
      ptr = MemNew (sizeof (Char) * (len + 2));
      if (ptr != NULL) {
        vnp = list;
        tmp = NULL;
        while (vnp != NULL) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
            if (tmp == NULL) {
              tmp = ptr;
            } else if (useReturn) {
              tmp = StringMove (tmp, "\n");
            } else if (IsJapanese () && (tmp - ptr > 2) &&
            		IsMBLetter (tmp - 2) && IsMBLetter (str)) {
              /* no space required between two Japanese letters. */
              tmp = tmp;
            } else if (str [0] != ',' && str [0] != ';' && str [0] != ':') {
              tmp = StringMove (tmp, " ");
            } else {
              tmp = StringMove (tmp, " ");
            }
            tmp = StringMove (tmp, str);
          }
          vnp = vnp->next;
        }
      }
    }
  }
  return ptr;
}

typedef struct deffeats {
  SeqFeatPtr  sfp;
  SeqFeatPtr  gene;
  SeqFeatPtr  prot;
  CharPtr     genename;
  CharPtr     protname;
  Boolean     alreadyTrimmed;
  Uint2       entityID;
  Uint2       itemID;
  Uint2       subtype;
  Boolean     lastInString;
  Boolean     lastInGroup;
  Boolean     lastInType;
  Boolean     lastInPenultimate;
  Boolean     pseudo;
  Boolean     ignore;
  Int2        altSplices;
  Int2        numUnknown;
} DefFeatsData, PNTR DefFeatsPtr;

static Boolean GetMolBioFeatsGatherFunc (GatherContextPtr gcp, Boolean getGene, Boolean getMrna)

{
  DefFeatsPtr  dfp;
  RnaRefPtr    rrp;
  SeqFeatPtr   sfp;
  CharPtr      str;
  Uint1        type;
  ValNodePtr   PNTR vnpp;

  if (gcp == NULL || gcp->thisitem == NULL || gcp->userdata == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  vnpp = (ValNodePtr PNTR) gcp->userdata;
  sfp = (SeqFeatPtr) gcp->thisitem;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      if (getGene) {
        dfp = MemNew (sizeof (DefFeatsData));
        if (dfp == NULL) return TRUE;
        dfp->entityID = gcp->entityID;
        dfp->itemID = gcp->itemID;
        dfp->sfp = sfp;
        dfp->subtype = FEATDEF_GENE;
        ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      }
      break;
    case SEQFEAT_CDREGION :
      dfp = MemNew (sizeof (DefFeatsData));
      if (dfp == NULL) return TRUE;
      dfp->entityID = gcp->entityID;
      dfp->itemID = gcp->itemID;
      dfp->sfp = sfp;
      dfp->subtype = FEATDEF_CDS;
      dfp->altSplices = 1;
      ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp == NULL) return TRUE;
      switch (rrp->type) {
        case 2 :
          if (getMrna) {
            dfp = MemNew (sizeof (DefFeatsData));
            if (dfp == NULL) return TRUE;
            dfp->entityID = gcp->entityID;
            dfp->itemID = gcp->itemID;
            dfp->sfp = sfp;
            dfp->subtype = FEATDEF_mRNA;
            ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          }
          break;
        case 3 :
          dfp = MemNew (sizeof (DefFeatsData));
          if (dfp == NULL) return TRUE;
          dfp->entityID = gcp->entityID;
          dfp->itemID = gcp->itemID;
          dfp->sfp = sfp;
          dfp->subtype = FEATDEF_tRNA;
          ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          break;
        case 4 :
          dfp = MemNew (sizeof (DefFeatsData));
          if (dfp == NULL) return TRUE;
          dfp->entityID = gcp->entityID;
          dfp->itemID = gcp->itemID;
          dfp->sfp = sfp;
          dfp->subtype = FEATDEF_rRNA;
          ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          break;
        case 255 :
          if (rrp->ext.choice == 1) {
            str = (CharPtr) rrp->ext.value.ptrvalue;
            if (StringICmp (str, "internal transcribed spacer 1") == 0 ||
                StringICmp (str, "internal transcribed spacer 2") == 0 ||
                StringICmp (str, "internal transcribed spacer 3") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS1") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS2") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS3") == 0 ||
                StringICmp (str, "ITS1") == 0 ||
                StringICmp (str, "ITS2") == 0 ||
                StringICmp (str, "ITS3") == 0) {
              dfp = MemNew (sizeof (DefFeatsData));
              if (dfp == NULL) return TRUE;
              dfp->entityID = gcp->entityID;
              dfp->itemID = gcp->itemID;
              dfp->sfp = sfp;
              dfp->subtype = FEATDEF_otherRNA;
              ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
            }
          }
          break;
        default :
          break;
      }
      break;
    case SEQFEAT_IMP :
      type = FindFeatDefType (sfp);
      if (type == FEATDEF_LTR || type == FEATDEF_exon) {
        dfp = MemNew (sizeof (DefFeatsData));
        if (dfp == NULL) return TRUE;
        dfp->entityID = gcp->entityID;
        dfp->itemID = gcp->itemID;
        dfp->sfp = sfp;
        dfp->subtype = type;
        ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      }
      break;
    default :
      break;
  }
  return TRUE;
}

static Boolean GetCDStRNArRNAGatherFunc (GatherContextPtr gcp)

{
  return GetMolBioFeatsGatherFunc (gcp, FALSE, FALSE);
}

static Boolean GetGeneCDStRNArRNAmRNAGatherFunc (GatherContextPtr gcp)

{
  return GetMolBioFeatsGatherFunc (gcp, TRUE, TRUE);
}

static void LabelAModifier (CharPtr str, CharPtr text, Boolean labelMods)

{
  Char     ch;
  CharPtr  ptr;

  if (str == NULL || text == NULL) return;
  if (StringHasNoText (text)) {
    str [0] = '\0';
    text [0] = '\0';
  } else {
    if (labelMods) {
      ptr = str;
      while (*ptr != '\0') {
        ch = *ptr;
        *ptr = TO_LOWER (ch);
        ptr++;
      }
      StringCat (str, " ");
    } else {
      str [0] = '\0';
    }
    StringCat (str, text);
  }
}

static int LIBCALLBACK SortCDStRNArRNAByLocation (VoidPtr ptr1, VoidPtr ptr2)

{
  BioseqPtr    bsp1;
  BioseqPtr    bsp2;
  DefFeatsPtr  dfp1;
  DefFeatsPtr  dfp2;
  Int4         leftend1;
  Int4         leftend2;
  Int4         rightend1;
  Int4         rightend2;
  SeqFeatPtr   sfp1;
  SeqFeatPtr   sfp2;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      dfp1 = (DefFeatsPtr) vnp1->data.ptrvalue;
      dfp2 = (DefFeatsPtr) vnp2->data.ptrvalue;
      if (dfp1 != NULL && dfp2 != NULL) {
        sfp1 = dfp1->sfp;
        sfp2 = dfp2->sfp;
        if (sfp1 != NULL && sfp2 != NULL) {
          bsp1 = GetBioseqGivenSeqLoc (sfp1->location, dfp1->entityID);
          bsp2 = GetBioseqGivenSeqLoc (sfp2->location, dfp2->entityID);
          if (bsp1 != NULL && bsp2 != NULL) {
            leftend1 = GetOffsetInBioseq (sfp1->location, bsp1, SEQLOC_LEFT_END);
            leftend2 = GetOffsetInBioseq (sfp2->location, bsp2, SEQLOC_LEFT_END);
            rightend1 = GetOffsetInBioseq (sfp1->location, bsp1, SEQLOC_RIGHT_END);
            rightend2 = GetOffsetInBioseq (sfp2->location, bsp2, SEQLOC_RIGHT_END);
            if (leftend1 > leftend2) {
              return 1;
            } else if (leftend1 < leftend2) {
              return -1;
            } else if (rightend1 > rightend2) {
              return 1;
            } else if (rightend1 < rightend2) {
              return -1;
            } else {
              return 0;
            }
          }
        }
      }
    }
  }
  return 0;
}

static int LIBCALLBACK SortCDSAfterExons (VoidPtr ptr1, VoidPtr ptr2)

{
  DefFeatsPtr  dfp1;
  DefFeatsPtr  dfp2;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      dfp1 = (DefFeatsPtr) vnp1->data.ptrvalue;
      dfp2 = (DefFeatsPtr) vnp2->data.ptrvalue;
      if (dfp1 != NULL && dfp2 != NULL) {
        if (dfp1->subtype == FEATDEF_CDS && dfp2->subtype == FEATDEF_exon) {
          return 1;
        } else if (dfp1->subtype == FEATDEF_exon && dfp2->subtype == FEATDEF_CDS) {
          return -1;
        } else {
          /* return 0; */
          return SortCDStRNArRNAByLocation (ptr1, ptr2);
        }
      }
    }
  }
  return 0;
}

extern EnumFieldAssoc  orgmod_subtype_alist [];
extern EnumFieldAssoc  subsource_subtype_alist [];

static Int2  orgmod_rank [24];
static Int2  subsource_rank [24];

static Boolean StrainAlreadyInParentheses (CharPtr taxname, CharPtr strain)

{
  size_t   len;
  CharPtr  ptr;

  ptr = StringChr (taxname, '(');
  if (ptr == NULL) return FALSE;
  ptr++;
  len = StringLen (strain);
  if (StringNCmp (taxname, strain, len) != 0) return FALSE;
  ptr += len;
  if (*ptr != ')') return FALSE;
  return TRUE;
}

static void AddOrgModsToDef (ValNodePtr PNTR stringsPtr, BioSourcePtr biop, Boolean labelMods)

{
  EnumFieldAssocPtr  ap;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  CharPtr            ptr;
  SubSourcePtr       ssp;
  Char               str [128];
  Char               text [64];

  if (stringsPtr != NULL && biop != NULL) {
    orp = biop->org;
    if (orp != NULL) {
      onp = orp->orgname;
      if (onp != NULL) {
        mod = onp->mod;
        while (mod != NULL) {
          if (mod->subtype < 24 && orgmod_rank [mod->subtype] > 0) {
            if (mod->subtype == 2 && StrainAlreadyInParentheses (orp->taxname, mod->subname)) {
              /* do not add strain if already parenthetical in organism name */
            } else {
              text [0] = '\0';
              str [0] = '\0';
              StringNCpy_0 (text, mod->subname, sizeof (text));
              for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
                if (ap->value == mod->subtype) {
                  StringNCpy_0 (str, ap->name, sizeof (str));
                }
              }
              LabelAModifier (str, text, labelMods);
              if (! StringHasNoText (str)) {
                ValNodeCopyStr (stringsPtr, orgmod_rank [mod->subtype], str);
              }
            }
          }
          mod = mod->next;
        }
      }
    }

    ssp = biop->subtype;
    while (ssp != NULL) {
      if (ssp->subtype < 24 && subsource_rank [ssp->subtype] > 0) {
        text [0] = '\0';
        str [0] = '\0';
        StringNCpy_0 (text, ssp->name, sizeof (text));
        for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
          if (ap->value == ssp->subtype) {
            StringNCpy_0 (str, ap->name, sizeof (str));
            ptr = StringStr (str, "-name");
            if (ptr != NULL) {
              *ptr = '\0';
            }
          }
        }
        LabelAModifier (str, text, labelMods);
        if (! StringHasNoText (str)) {
          ValNodeCopyStr (stringsPtr, subsource_rank [ssp->subtype], str);
        }
      }
      ssp = ssp->next;
    }
  }
}

static void FinishAutoDefProc (Uint2 entityID, SeqEntryPtr sep,
                               ValNodePtr head, BioseqPtr target,
                               SeqEntryPtr nsep, MolInfoPtr mip,
                               ValNodePtr strings, ValNodePtr xtras)

{
  Int2          count;
  DefFeatsPtr   dfp;
  GBQualPtr     gbq;
  DefFeatsPtr   nextdfp;
  CharPtr       ptr;
  RnaRefPtr     rrp;
  SeqFeatPtr    sfp;
  Char          str [380];
  tRNAPtr       trna;
  ValNodePtr    ttl;
  Char          text [64];
  ValNodePtr    vnp;

  vnp = head;
  count = 0;
  while (vnp != NULL) {
    str [0] = '\0';
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL) {
      sfp = dfp->sfp;
      if (sfp != NULL || dfp->numUnknown > 0) {
        count++;
        /* FindGeneAndProtForCDS (entityID, sfp, &(dfp->gene), &(dfp->prot)); */
        /* StringCpy (str, "unknown"); */
        text [0] = '\0';
        if (dfp->subtype == FEATDEF_CDS && dfp->prot != NULL) {
          if (dfp->protname != NULL) {
            if (StringICmp (dfp->protname, "unknown") == 0 && (! StringHasNoText (dfp->genename))) {
              StringNCpy_0 (text, dfp->genename, sizeof (text));
              StringCat (str, "(");
              StringCat (str, text);
              StringCat (str, ")");
            } else {
              StringNCpy_0 (str, dfp->protname, sizeof (str) - 100);
              if (dfp->genename != NULL) {
                StringNCpy_0 (text, dfp->genename, sizeof (text));
                if (! StringHasNoText (text)) {
                  StringCat (str, " (");
                  StringCat (str, text);
                  StringCat (str, ")");
                }
              }
            }
            if (dfp->lastInGroup || dfp->lastInType) {
              if (mip != NULL) {
                ptr = StringStr (str, "precursor");
                if (ptr != NULL && StringICmp (ptr, "precursor") == 0) {
                  StringCat (str, ",");
                }
                if (mip->biomol == MOLECULE_TYPE_MRNA) {
                  StringCat (str, " mRNA");
                } else {
                  StringCat (str, " gene");
                }
              }
              if (count > 1) {
                StringCat (str, "s");
              }
              if (dfp->altSplices > 1) {
                StringCat (str, ", alternative splice products");
              }
              if (dfp->sfp->partial) {
                StringCat (str, ", partial cds");
              } else {
                StringCat (str, ", complete cds");
              }
            }
          }
        } else if (dfp->subtype == FEATDEF_CDS && dfp->pseudo) {
          if (dfp->genename != NULL) {
            StringNCpy_0 (str, dfp->genename, sizeof (str) - 50);
          }
          if (dfp->lastInGroup || dfp->lastInType) {
            StringCat (str, " pseudogene");
            if (count > 1) {
              StringCat (str, "s");
            }
            if (dfp->sfp->partial) {
              StringCat (str, ", partial sequence");
            } else {
              StringCat (str, ", complete sequence");
            }
          }
        } else if (dfp->subtype == FEATDEF_CDS && dfp->prot == NULL && dfp->genename != NULL) {
          StringNCpy_0 (str, dfp->genename, sizeof (str) - 50);
          if (dfp->lastInGroup || dfp->lastInType) {
            if (mip != NULL) {
              ptr = StringStr (str, "precursor");
              if (ptr != NULL && StringICmp (ptr, "precursor") == 0) {
                StringCat (str, ",");
              }
              if (mip->biomol == MOLECULE_TYPE_MRNA) {
                StringCat (str, " mRNA");
              } else {
                StringCat (str, " gene");
              }
            }
            if (count > 1) {
              StringCat (str, "s");
            }
            if (dfp->sfp->partial) {
              StringCat (str, ", partial cds");
            } else {
              StringCat (str, ", complete cds");
            }
          }
        } else if (dfp->subtype == FEATDEF_rRNA || dfp->subtype == FEATDEF_otherRNA) {
          rrp = (RnaRefPtr) dfp->sfp->data.value.ptrvalue;
          if (rrp != NULL) {
            if (rrp->ext.choice == 1) {
              StringNCpy_0 (str, (CharPtr) rrp->ext.value.ptrvalue, sizeof (str) - 50);
              if (dfp->subtype == FEATDEF_rRNA) {
                ptr = StringStr (str, " rRNA");
                if (ptr != NULL) {
                  *ptr = '\0';
                }
                ptr = StringStr (str, " ribosomal RNA");
                if (ptr != NULL) {
                  *ptr = '\0';
                }
                if (! StringHasNoText (str)) {
                  StringCat (str, " ribosomal RNA");
                }
              } else if (dfp->subtype == FEATDEF_otherRNA) {
              }
              if (dfp->genename != NULL) {
                StringNCpy_0 (text, dfp->genename, sizeof (text));
                if (! StringHasNoText (text)) {
                  StringCat (str, " (");
                  StringCat (str, text);
                  StringCat (str, ")");
                }
              }
              if (dfp->lastInString || dfp->lastInGroup || dfp->lastInType) {
                if (dfp->subtype == FEATDEF_rRNA) {
                  StringCat (str, " gene");
                  if (count > 1) {
                    StringCat (str, "s");
                  }
                }
              }
              if (dfp->lastInGroup || dfp->lastInType) {
                if (dfp->sfp->partial) {
                  StringCat (str, ", partial sequence");
                } else {
                  StringCat (str, ", complete sequence");
                }
              }
            }
          }
        } else if (dfp->subtype == FEATDEF_tRNA) {
          rrp = (RnaRefPtr) dfp->sfp->data.value.ptrvalue;
          if (rrp != NULL) {
            if (rrp->ext.choice == 2) {
              trna = rrp->ext.value.ptrvalue;
              if (trna != NULL) {
                if (FeatDefLabel (dfp->sfp, str, sizeof (str) - 2, OM_LABEL_CONTENT) > 0) {
                  if (dfp->genename != NULL) {
                    StringNCpy_0 (text, dfp->genename, sizeof (text));
                    if (! StringHasNoText (text)) {
                      StringCat (str, " (");
                      StringCat (str, text);
                      StringCat (str, ")");
                    }
                  }
                  if (dfp->lastInGroup || dfp->lastInType) {
                    StringCat (str, " gene");
                    if (count > 1) {
                      StringCat (str, "s");
                    }
                    if (dfp->sfp->partial) {
                      StringCat (str, ", partial sequence");
                    } else {
                      StringCat (str, ", complete sequence");
                    }
                  }
                }
              }
            }
          }
        } else if (dfp->subtype == FEATDEF_LTR) {
          if (! StringHasNoText (sfp->comment)) {
            StringNCpy_0 (str, sfp->comment, sizeof (str));
            ptr = StringStr (str, " long terminal repeat");
            if (ptr != NULL) {
              *ptr = '\0';
            }
            ptr = StringStr (str, " long terminal repeat");
            if (ptr != NULL) {
              *ptr = '\0';
            }
            if (! StringHasNoText (str)) {
              StringCat (str, " long terminal repeat");
            }
          } else {
            /* StringCpy (str, "uncharacterized"); */
            StringCpy (str, "long terminal repeat");
          }
          if (dfp->lastInGroup || dfp->lastInType) {
            if (dfp->sfp->partial) {
              StringCat (str, ", partial sequence");
            } else {
              StringCat (str, ", complete sequence");
            }
          }
        } else if (dfp->subtype == FEATDEF_exon && target != NULL) {
          if (dfp->genename != NULL) {
            StringNCpy_0 (str, dfp->genename, sizeof (text));
            if (! StringHasNoText (str)) {
              StringCat (str, " exon");
              gbq = sfp->qual;
              while (gbq != NULL) {
                if (StringICmp (gbq->qual, "number") == 0) {
                  StringCat (str, " ");
                  StringCat (str, gbq->val);
                }
                gbq = gbq->next;
              }
            }
          } else {
            StringCpy (str, "uncharacterized exon");
          }
          if (dfp->lastInGroup || dfp->lastInType) {
            if (dfp->sfp->partial) {
              StringCat (str, ", partial sequence");
            } else {
              StringCat (str, ", complete sequence");
            }
          }
        } else if (dfp->numUnknown > 0) {
          if (mip != NULL && mip->biomol == MOLECULE_TYPE_MRNA) {
            StringCat (str, "unknown mRNA");
          } else {
            StringCat (str, "unknown gene");
          }
          if (dfp->numUnknown > 1) {
            StringCat (str, "s");
          }
        }
      }
    }
    vnp = vnp->next;
    if (! StringHasNoText (str)) {
      if (vnp == NULL) {
        StringCat (str, ".");
      } else if (vnp->next == NULL) {
        nextdfp = (DefFeatsPtr) vnp->data.ptrvalue;
        if (dfp->lastInPenultimate) {
          if ((dfp->subtype == FEATDEF_rRNA && nextdfp->subtype == FEATDEF_otherRNA &&
               dfp->sfp->partial == nextdfp->sfp->partial) ||
              (dfp->subtype == FEATDEF_otherRNA && nextdfp->subtype == FEATDEF_rRNA &&
               dfp->sfp->partial == nextdfp->sfp->partial)) {
            StringCat (str, ", and");
          } else {
            StringCat (str, "; and");
          }
        } else if (dfp->lastInType || dfp->lastInGroup) {
          if (count > 1) {
            StringCat (str, ", and");
          } else {
            StringCat (str, " and");
          }
        } else if (nextdfp->lastInType || nextdfp->lastInGroup) {
          if (count > 1) {
            StringCat (str, ", and");
          } else {
            StringCat (str, " and");
          }
        } else {
          if (count > 1) {
            StringCat (str, ", and");
          } else {
            StringCat (str, " and");
          }
        }
      } else {
        nextdfp = (DefFeatsPtr) vnp->data.ptrvalue;
        if (dfp->lastInPenultimate) {
          if ((dfp->subtype == FEATDEF_rRNA && nextdfp->subtype == FEATDEF_otherRNA &&
               dfp->sfp->partial == nextdfp->sfp->partial) ||
              (dfp->subtype == FEATDEF_otherRNA && nextdfp->subtype == FEATDEF_rRNA &&
               dfp->sfp->partial == nextdfp->sfp->partial)) {
            StringCat (str, ", and");
          } else {
            StringCat (str, "; and");
          }
        } else if (dfp->lastInType || dfp->lastInGroup) {
          StringCat (str, ";");
        } else if (nextdfp->lastInType || nextdfp->lastInGroup) {
          if (count > 1) {
            StringCat (str, ", and");
          } else {
            StringCat (str, " and");
          }
        } else {
          StringCat (str, ",");
        }
      }
      ValNodeCopyStr (&strings, 0, str);
    }
    if (dfp->lastInString || dfp->lastInGroup || dfp->lastInType) {
      count = 0;
    }
  }

  ptr = MergeValNodeStrings (strings, FALSE);
  if (nsep != NULL) {
    ttl = SeqEntryGetSeqDescr (nsep, Seq_descr_title, NULL);
    if (ttl == NULL) {
      ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
    }
    if (ttl == NULL) {
      ttl = CreateNewDescriptor (nsep, Seq_descr_title);
    }
    if (ttl != NULL) {
      MemFree (ttl->data.ptrvalue);
      ttl->data.ptrvalue = ptr;
      ptr = NULL;
    }
  }
  MemFree (ptr);
  ValNodeFreeData (strings);
  ValNodeFreeData (head);
  ValNodeFreeData (xtras);
}

static void CombineProteinNames (DefFeatsPtr dfp1, DefFeatsPtr dfp2)

{
  Int2     i, j, lastspace, lastdash, lastcomma;
  size_t   len1, len2;
  CharPtr  str1 = NULL, str2 = NULL;

  if (dfp1 == NULL || dfp2 == NULL) return;
  if (dfp1->protname == NULL && dfp2->protname == NULL) return;
  len1 = StringLen (dfp1->protname);
  len2 = StringLen (dfp2->protname);
  if (len1 < 1 || len2 < 1) return;
  i = 0;
  j = 0;
  lastspace = 0;
  lastdash = 0;
  lastcomma = 0;
  while (i < len1 && j < len2 && dfp1->protname [i] == dfp2->protname [j]) {
    if (dfp1->protname [i] == ' ') {
      lastspace = i;
    }
    if (dfp1->protname [i] == '-') {
      lastdash = i;
    }
    if (dfp1->protname [i] == ',') {
      lastcomma = i;
    }
    i++;
    j++;
  }
  str1 = StringSave (dfp1->protname);
  if (str1 != NULL) {
    str1 [i] = '\0';
    if (! dfp1->alreadyTrimmed) {
      if (lastcomma > 0) {
        str1 [lastcomma] = '\0';
        dfp1->alreadyTrimmed = TRUE;
      } else if (lastdash > 0) {
        str1 [lastdash] = '\0';
        dfp1->alreadyTrimmed = TRUE;
      } else if (lastspace > 0) {
        str1 [lastspace] = '\0';
        dfp1->alreadyTrimmed = TRUE;
      }
    }
  }

  i = len1;
  j = len2;
  while (i > 0 && j > 0 && dfp1->protname [i - 1] == dfp2->protname [j - 1]) {
    i--;
    j--;
  }
  str2 = StringSave (dfp1->protname + i);
  TrimSpacesAroundString (str1);
  TrimSpacesAroundString (str2);
  len1 = StringLen (str1);
  len2 = StringLen (str2);
  if (len1 > len2) {
    dfp1->protname = str1;
    MemFree (str2);
  } else {
    dfp1->protname = str2;
    MemFree (str1);
  }
}

static Boolean AreAltSpliceGenes (DefFeatsPtr dfp1, DefFeatsPtr dfp2, ValNodePtr PNTR xtras)

{
  Int2        comp;
  SeqFeatPtr  sfp1, sfp2;

  if (dfp1 == NULL || dfp2 == NULL) return FALSE;
  if (dfp1->prot == NULL || dfp2->prot == NULL) return FALSE;
  if (dfp1->pseudo || dfp2->pseudo) return FALSE;
  sfp1 = dfp1->sfp;
  sfp2 = dfp2->sfp;
  if (sfp1 == NULL || sfp2 == NULL) return FALSE;
  if (sfp1->partial != sfp2->partial) return FALSE;
  if (SeqLocStrand (sfp1->location) != SeqLocStrand (sfp2->location)) return FALSE;
  comp = SeqLocCompare (sfp1->location, sfp2->location);
  if (comp == SLC_NO_MATCH) return FALSE;
  if (dfp1->genename == NULL || dfp2->genename == NULL) return FALSE;
  if (StringICmp (dfp1->genename, dfp2->genename) != 0) return FALSE;
  CombineProteinNames (dfp1, dfp2);
  return TRUE;
}

static ValNodePtr MergeAltSpliceCDSs (ValNodePtr head)

{
  DefFeatsPtr  dfp1, dfp2;
  ValNodePtr   nextvnp;
  ValNodePtr   PNTR prevvnp;
  ValNodePtr   vnp1, vnp2;
  ValNodePtr   xtras;

  xtras = NULL;
  vnp1 = head;
  while (vnp1 != NULL) {
    dfp1 = (DefFeatsPtr) vnp1->data.ptrvalue;
    if (dfp1 != NULL && dfp1->subtype == FEATDEF_CDS) {
      vnp2 = vnp1->next;
      prevvnp = &(vnp1->next);
      while (vnp2 != NULL) {
        nextvnp = vnp2->next;
        dfp2 = (DefFeatsPtr) vnp2->data.ptrvalue;
        if (dfp2 != NULL && dfp2->subtype == FEATDEF_CDS) {
          if (AreAltSpliceGenes (dfp1, dfp2, &xtras)) {
            (dfp1->altSplices)++;
            *prevvnp = vnp2->next;
            vnp2->next = NULL;
            ValNodeFreeData (vnp2);
          } else {
            prevvnp = (ValNodePtr PNTR) &(vnp2->next);
          }
        } else {
          prevvnp = (ValNodePtr PNTR) &(vnp2->next);
        }
        vnp2 = nextvnp;
      }
    }
    vnp1 = vnp1->next;
  }
  return xtras;
}

static void AutoDefProc (Uint2 entityID, SeqEntryPtr sep, Boolean addMods,
                         Boolean labelMods, Int2 maxMods, Boolean leaveInParen,
                         BioseqPtr target, BioseqPtr seg, ValNodePtr nonUniqueOrgs)

{
  BioseqContextPtr  bcp;
  BioSourcePtr    biop;
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  Boolean         change;
  CdRegionPtr     crp;
  DefFeatsPtr     dfp;
  GBQualPtr       gbq;
  Int2            group;
  GeneRefPtr      grp;
  GatherScope     gs;
  ValNodePtr      head;
  SeqLocPtr       lastslp;
  MolInfoPtr      mip;
  DefFeatsPtr     nextdfp;
  SeqLocPtr       nextslp;
  ValNodePtr      nextvnp;
  SeqEntryPtr     nsep;
  Int2            numUnknown;
  BioseqPtr       part;
  DefFeatsPtr     penult;
  ValNodePtr      PNTR prevvnp;
  ProtRefPtr      prp;
  CharPtr         ptr;
  ValNodePtr      sdp;
  SeqFeatPtr      sfp;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  Char            str [128];
  ValNodePtr      ttl;
  ValNodePtr      strings;
  ValNode         vn;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;
  ValNodePtr      xtras = NULL;

  if (sep == NULL) return;
  if (target == NULL && seg == NULL) {
    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp == NULL) return;
      if (bssp->_class == 7 || bssp->_class == 13 ||
          bssp->_class == 14 || bssp->_class == 15) {
        for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
          AutoDefProc (entityID, sep, addMods, labelMods, maxMods, leaveInParen, NULL, NULL, nonUniqueOrgs);
        }
        return;
      }
    }

    nsep = FindNucSeqEntry (sep);
    if (nsep != NULL) {
      bsp = (BioseqPtr) nsep->data.ptrvalue;
      if (bsp != NULL && bsp->repr == Seq_repr_seg && bsp->seq_ext != NULL && bsp->seq_ext_type == 1) {
        vn.choice = SEQLOC_MIX;
        vn.next = NULL;
        vn.data.ptrvalue = bsp->seq_ext;
        slp = SeqLocFindNext (&vn, NULL);
        while (slp != NULL) {
          nextslp = SeqLocFindNext (&vn, slp);
          sip = SeqLocId (slp);
          if (sip != NULL) {
            part = BioseqFind (sip);
            if (part != NULL) {
              AutoDefProc (entityID, sep, addMods, labelMods, maxMods, leaveInParen, part, NULL, nonUniqueOrgs);
            }
          }
          slp = nextslp;
        }
        AutoDefProc (entityID, sep, addMods, labelMods, maxMods, leaveInParen, NULL, bsp, nonUniqueOrgs);
        return;
      }
    }
  } else if (target != NULL) {
    nsep = SeqMgrGetSeqEntryForData (target);
  } else {
    nsep = SeqMgrGetSeqEntryForData (seg);
  }

  biop = NULL;
  strings = NULL;
  str [0] = '\0';

  SeqEntryToBioSource (sep, NULL, str, sizeof (str) - 1, &biop);
  if (! leaveInParen) {
    ptr = StringStr (str, "(");
    if (ptr != NULL) {
      *ptr = '\0';
    }
  }
  TrimSpacesAroundString (str);
  if (StringICmp (str, "Human immunodeficiency virus type 1") == 0) {
    StringCpy (str, "HIV-1");
  } else if (StringICmp (str, "Human immunodeficiency virus type 2") == 0) {
    StringCpy (str, "HIV-2");
  }
  str [0] = TO_UPPER (str [0]);
  ValNodeCopyStr (&strings, 0, str);

  mip = NULL;
  if (nsep != NULL) {
    bsp = (BioseqPtr) nsep->data.ptrvalue;
    bcp = BioseqContextNew (bsp);
    sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
    BioseqContextFree (bcp);
    if (sdp != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        if (mip->tech == MI_TECH_htgs_1 ||
            mip->tech == MI_TECH_htgs_2 ||
            mip->tech == MI_TECH_est ||
            mip->tech == MI_TECH_sts ||
            mip->tech == MI_TECH_survey) {
          ttl = ValNodeExtract (&(bsp->descr), Seq_descr_title);
          if (ttl != NULL) {
            ttl = ValNodeFreeData (ttl);
          }
          return;
        }
      }
    }
  }

  /*
  nsep = FindNucSeqEntry (sep);
  if (nsep != NULL) {
    sdp = SeqEntryGetSeqDescr (nsep, Seq_descr_molinfo, NULL);
    if (sdp == NULL) {
      sdp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    }
    if (sdp != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
    }
  }
  */

  if (nonUniqueOrgs != NULL) {
    for (vnp = nonUniqueOrgs; vnp != NULL; vnp = vnp->next) {
      if (StringICmp ((CharPtr) vnp->data.ptrvalue, str) == 0) {
        if (vnp->choice == 1) {
          addMods = FALSE; /* if only one organism in record, already unique defline */
        }
      }
    }
  }
  if (addMods) {
    AddOrgModsToDef (&strings, biop, labelMods);
    strings = SortValNode (strings, SortByVnpChoice);
    vnp = strings;
    while (vnp != NULL && maxMods > 0) {
      maxMods--;
      vnp = vnp->next;
    }
    if (vnp != NULL) {
      vnp->next = ValNodeFreeData (vnp->next);
    }
  }

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  gs.target = NULL;
  if (target != NULL) {
    slp = ValNodeNew (NULL);
    slp->choice = SEQLOC_WHOLE;
    sip = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (target->id, 0)));
    slp->data.ptrvalue = sip;
    gs.target = slp;
  }
  head = NULL;
  GatherEntity (entityID, (Pointer) (&head), GetCDStRNArRNAGatherFunc, &gs);
  gs.target = SeqLocFree (gs.target);
  head = SortValNode (head, SortCDStRNArRNAByLocation);

  if (head == NULL && mip != NULL && mip->tech == MI_TECH_survey) {
    ValNodeCopyStr (&strings, 0, ", genome survey sequence.");
  }

  numUnknown = 0;
  vnp = head;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL) {
      sfp = dfp->sfp;
      if (sfp != NULL) {
        FindGeneAndProtForCDS (entityID, sfp, &(dfp->gene), &(dfp->prot));
        dfp->pseudo = FALSE;
        dfp->genename = NULL;
        grp = NULL;
        if (dfp->gene != NULL) {
          grp = (GeneRefPtr) dfp->gene->data.value.ptrvalue;
          if (grp != NULL && grp->pseudo) {
            dfp->pseudo = TRUE;
          }
        }
        xref = sfp->xref;
        while (xref != NULL && xref->data.choice != SEQFEAT_GENE) {
          xref = xref->next;
        }
        if (xref != NULL) {
          grp = (GeneRefPtr) xref->data.value.ptrvalue;
        }
        if (grp != NULL) {
          dfp->genename = (CharPtr) grp->locus;
          if (grp->pseudo) {
            dfp->pseudo = TRUE;
          }
        }
        if (dfp->subtype == FEATDEF_CDS) {
          if (target != NULL) {
            lastslp = NULL;
            slp = SeqLocFindNext (sfp->location, NULL);
            while (slp != NULL) {
              lastslp = slp;
              slp = SeqLocFindNext (sfp->location, slp);
            }
            if (lastslp != NULL && nsep != NULL) {
              bsp = (BioseqPtr) nsep->data.ptrvalue;
              if (GetBioseqGivenSeqLoc (lastslp, entityID) != bsp) {
                dfp->ignore = TRUE;
              }
            }
          }
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
            if (crp->orf) {
              dfp->ignore = TRUE;
            }
          }
          gbq = sfp->qual;
          while (gbq != NULL) {
            if (StringICmp (gbq->qual, "pseudo") == 0) {
              dfp->pseudo = TRUE;
            }
            gbq = gbq->next;
          }
          if (dfp->pseudo) {
          } else if (dfp->prot == NULL) {
            if (dfp->gene == NULL) {
              dfp->ignore = TRUE;
            }
          } else {
            prp = (ProtRefPtr) dfp->prot->data.value.ptrvalue;
            if (prp != NULL) {
              if (prp->name == NULL || StringHasNoText ((CharPtr) prp->name->data.ptrvalue)) {
                if (prp->desc == NULL || StringHasNoText (prp->desc)) {
                  dfp->ignore = TRUE;
                } else {
                  dfp->protname = prp->desc;
                }
              } else {
                dfp->protname = (CharPtr) prp->name->data.ptrvalue;
              }
              if (! dfp->ignore) {
                if (StringICmp (dfp->protname, "unknown") == 0) {
                  if (StringHasNoText (dfp->genename)) {
                    numUnknown++;
                    dfp->ignore = TRUE;
                  }
                }
              }
            }
          }
        } else if (dfp->subtype == FEATDEF_exon && target == NULL) {
          dfp->ignore = TRUE;
        }
      }
    }
    vnp = vnp->next;
  }

  vnp = head;
  prevvnp = &head;
  while (vnp != NULL) {
    nextvnp = vnp->next;
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp->ignore) {
      *prevvnp = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      prevvnp = (ValNodePtr PNTR) &(vnp->next);
    }
    vnp = nextvnp;
  }

  xtras = MergeAltSpliceCDSs (head);

  if (target != NULL) {
    head = SortValNode (head, SortCDSAfterExons);
  }

  if (numUnknown > 0) {
    dfp = (DefFeatsPtr) MemNew (sizeof (DefFeatsData));
    if (dfp != NULL) {
      dfp->entityID = entityID;
      dfp->subtype = 0;
      dfp->numUnknown = numUnknown;
      ValNodeAddPointer (&head, 0, (Pointer) dfp);
    }
  }

  vnp = head;
  group = 0;
  penult = NULL;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    vnp->choice = (Uint1) group;
    vnp = vnp->next;
    if (vnp != NULL) {
      change = FALSE;
      nextdfp = (DefFeatsPtr) vnp->data.ptrvalue;
      if ((dfp->subtype == FEATDEF_rRNA && nextdfp->subtype == FEATDEF_otherRNA) ||
          (dfp->subtype == FEATDEF_otherRNA && nextdfp->subtype == FEATDEF_rRNA)) {
        dfp->lastInString = TRUE;
        if (dfp->sfp->partial != nextdfp->sfp->partial) {
          dfp->lastInGroup = TRUE;
        }
        change = TRUE;
      } else if (dfp->subtype != nextdfp->subtype) {
        dfp->lastInType = TRUE;
        change = TRUE;
      } else if (dfp->sfp->partial != nextdfp->sfp->partial) {
        dfp->lastInGroup = TRUE;
        change = TRUE;
      } else if (dfp->pseudo != nextdfp->pseudo) {
        dfp->lastInGroup = TRUE;
        change = TRUE;
      } else if (dfp->altSplices > 1 || nextdfp->altSplices > 1) {
        dfp->lastInGroup = TRUE;
        change = TRUE;
      }
      if (change) {
        group++;
        penult = dfp;
      }
    } else {
      dfp->lastInString = TRUE;
      dfp->lastInGroup = TRUE;
      dfp->lastInType = TRUE;
      group++;
    }
  }
  if (penult != NULL) {
    penult->lastInPenultimate = TRUE;
  }

  FinishAutoDefProc (entityID, sep, head, target, nsep, mip, strings, xtras);
}

static CharPtr sourceModRankList [] = {
  "Strain", "Isolate", "Clone", "Type", "Cultivar", "Haplotype",
  "Substrain", "Subclone", "Subtype", "Serotype", "Serogroup", "Serovar",
  "Variety", "Pathovar", "Chemovar", "Biovar", "Biotype", "Group", "Subgroup",
  "Cell-line", "Cell-type", "Tissue-type", "Clone-lib", "Tissue-lib", "Dev-stage",
  "Lab-host", "Pop-variant", "Frequency", "Germline", "Rearranged",
  "Chromosome", "Map", "Genotype", "Sex", "Plasmid-name", "Transposon-name",
  "Ins-seq-name", "Plastid-name", "Country", "Old Name", "Common", "Acronym",
  "Dosage", "Natural-host", "Sub-species", "Specimen-voucher",
  NULL
};

typedef struct deflineform {
  FORM_MESSAGE_BLOCK
  SeqEntryPtr    sep;
  ButtoN         addLabels;
  GrouP          customGrp;
  GrouP          sourceListGrp;
  ButtoN         PNTR sourceBoxList;
  PopuP          modLimit;
  ButtoN         onlyModifyTargeted;
  BioseqPtr      target;
  ButtoN         leaveInParentheses;
  Boolean        smartMods;
} DeflineForm, PNTR DeflineFormPtr;

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

static Boolean GatherOrgnamesFunc (GatherContextPtr gcp)

{
  BioSourcePtr     biop;
  OrgRefPtr        orp;
  ValNodePtr       sdp;
  SeqFeatPtr       sfp;
  CharPtr          str;
  ValNodePtr PNTR  vnpp;

  if (gcp == NULL || gcp->thisitem == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT  && gcp->thistype != OBJ_SEQDESC) return TRUE;
  vnpp = (ValNodePtr PNTR) gcp->userdata;
  if (vnpp == NULL) return TRUE;

  orp = NULL;
  biop = NULL;
  switch (gcp->thistype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      switch (sfp->data.choice) {
        case SEQFEAT_ORG :
          orp = (OrgRefPtr) sfp->data.value.ptrvalue;
          break;
        case SEQFEAT_BIOSRC :
          biop = (BioSourcePtr) sfp->data.value.ptrvalue;
          break;
        default :
          break;
      }
      break;
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) gcp->thisitem;
      switch (sdp->choice) {
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
  }
  if (orp != NULL) {
    str = StringSave (orp->taxname);
    if (str == NULL) return TRUE;
    TrimSpacesAroundString (str);
    if (StringICmp (str, "Human immunodeficiency virus type 1") == 0) {
      StringCpy (str, "HIV-1");
    } else if (StringICmp (str, "Human immunodeficiency virus type 2") == 0) {
      StringCpy (str, "HIV-2");
    }
    str [0] = TO_UPPER (str [0]);
    ValNodeAddStr (vnpp, 1, str);
  }
  return TRUE;
}

static void DefLineModFormAcceptProc (ButtoN b)

{
  EnumFieldAssocPtr  ap;
  Int2               count;
  DeflineFormPtr     dfp;
  GatherScope        gs;
  Boolean            labelMods;
  Boolean            leaveInParen;
  Int2               maxMods;
  ValNodePtr         nextvnp;
  ValNodePtr         nonUniqueOrgs;
  Int2               val;
  ValNodePtr         vnp;

  dfp = (DeflineFormPtr) GetObjectExtra (b);
  if (dfp == NULL) return;
  WatchCursor ();
  Hide (dfp->form);
  Update ();
  labelMods = GetStatus (dfp->addLabels);

  if (dfp->sourceBoxList != NULL && GetValue (dfp->customGrp) == 2) {
    count = 0;
    while (sourceModRankList [count] != NULL && dfp->sourceBoxList [count] != NULL) {
      if (! GetStatus (dfp->sourceBoxList [count])) {
        for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
          if (StringICmp (ap->name, sourceModRankList [count]) == 0 &&
              ap->value > 0 && ap->value < 24) {
            orgmod_rank [ap->value] = 0;
          }
        }
        for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
          if (StringICmp (ap->name, sourceModRankList [count]) == 0 &&
              ap->value > 0 && ap->value < 24) {
            subsource_rank [ap->value] = 0;
          }
        }
      }
      count++;
    }
  }

  maxMods = INT2_MAX;
  val = GetValue (dfp->modLimit);
  if (val > 1) {
    maxMods = val - 1;
  }

  if (dfp->onlyModifyTargeted != NULL && GetStatus (dfp->onlyModifyTargeted)) {
    dfp->sep = GetBestTopParentForData (dfp->input_entityID, dfp->target);
  }

  leaveInParen = FALSE;
  if (dfp->leaveInParentheses != NULL && GetStatus (dfp->leaveInParentheses)) {
    leaveInParen = TRUE;
  }

  nonUniqueOrgs = NULL;
  if (dfp->smartMods) {
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    gs.ignore[OBJ_SEQDESC] = FALSE;
    GatherEntity (dfp->input_entityID, (Pointer) (&nonUniqueOrgs), GatherOrgnamesFunc, &gs);
    nonUniqueOrgs = SortValNode (nonUniqueOrgs, SortByName);
    vnp = nonUniqueOrgs;
    while (vnp != NULL) {
      nextvnp = vnp->next;
      if (nextvnp != NULL && StringICmp ((CharPtr) vnp->data.ptrvalue, (CharPtr) nextvnp->data.ptrvalue) == 0) {
        vnp->next = nextvnp->next;
        nextvnp->next = NULL;
        ValNodeFreeData (nextvnp);
        nextvnp = vnp;
        (vnp->choice)++;
      }
      vnp = nextvnp;
    }
  }

  AutoDefProc (dfp->input_entityID, dfp->sep, TRUE, labelMods, maxMods, leaveInParen, NULL, NULL, nonUniqueOrgs);
  ValNodeFreeData (nonUniqueOrgs);
  ArrowCursor ();
  Remove (dfp->form);
  Update ();
  ObjMgrSetDirtyFlag (dfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dfp->input_entityID, 0, 0);
}

static void ChangeCustomGrp (GrouP g)

{
  DeflineFormPtr  dfp;

  dfp = (DeflineFormPtr) GetObjectExtra (g);
  if (dfp == NULL) return;
  if (GetValue (g) == 1) {
    SafeDisable (dfp->sourceListGrp);
  } else {
    SafeEnable (dfp->sourceListGrp);
  }
}

static void CleanupDeflineForm (GraphiC g, VoidPtr data)

{
  DeflineFormPtr  dfp;

  dfp = (DeflineFormPtr) data;
  if (dfp != NULL) {
    MemFree (dfp->sourceBoxList);
  }
  StdCleanupFormProc (g, data);
}

static void DeflineMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

static ForM CreateDefLineModForm (Uint2 entityID, SeqEntryPtr sep, BioseqPtr target,
                                  Boolean smartMods, Int2 maxCount)

{
  ButtoN          b;
  GrouP           c;
  Int2            count;
  DeflineFormPtr  dfp;
  GrouP           h;
  Int2            i;
  GrouP           p;
  GrouP           q;
  Char            str [16];
  WindoW          w;

  w = NULL;
  dfp = MemNew (sizeof (DeflineForm));
  if (dfp != NULL) {
    w = FixedWindow (-50, -33, -10, -10, "Defline Modifier Customization", StdCloseWindowProc);
    SetObjectExtra (w, dfp, CleanupDeflineForm);
    dfp->form = (ForM) w;
    dfp->formmessage = DeflineMessageProc;
    dfp->input_entityID = entityID;
    dfp->sep = sep;
    dfp->target = target;
    dfp->smartMods = smartMods;

    h = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (h, 10, 10);

    dfp->addLabels = CheckBox (h, "Use Labels (e.g., 'strain BALB/c')", NULL);
    if (smartMods) {
      SetStatus (dfp->addLabels, TRUE);
    }

    p = NormalGroup (h, -1, 0, "Modifier Classes", programFont, ChangeCustomGrp);
    SetGroupSpacing (p, 10, 10);

    dfp->customGrp = HiddenGroup (p, 2, 0, ChangeCustomGrp);
    SetObjectExtra (dfp->customGrp, dfp, NULL);
    RadioButton (dfp->customGrp, "All");
    RadioButton (dfp->customGrp, "Custom");
    SetValue (dfp->customGrp, 1);

    dfp->sourceListGrp = HiddenGroup (p, 3, 0, NULL);
    count = 0;
    while (sourceModRankList [count] != NULL) {
      count++;
    }
    dfp->sourceBoxList = MemNew (sizeof (ButtoN) * (count + 3));
    count = 0;
    if (dfp->sourceBoxList != NULL) {
      while (sourceModRankList [count] != NULL && count < maxCount) {
        dfp->sourceBoxList [count] = CheckBox (dfp->sourceListGrp, sourceModRankList [count], NULL);
        count++;
      }
    }
    Disable (dfp->sourceListGrp);

    q = HiddenGroup (h, 4, 0, NULL);
    StaticPrompt (q, "Maximum modifiers per line", 0, popupMenuHeight, programFont, 'l');
    dfp->modLimit = PopupList (q, TRUE, NULL);
    PopupItem (dfp->modLimit, "no limit");
    for (i = 1; i <= count; i++) {
      sprintf (str, "%d", (int) i);
      PopupItem (dfp->modLimit, str);
    }
    if (smartMods) {
      SetValue (dfp->modLimit, 2);
    } else {
      SetValue (dfp->modLimit, 1);
    }

    dfp->leaveInParentheses = CheckBox (w, "Leave in parenthetical organism info", NULL);
    SetStatus (dfp->leaveInParentheses, TRUE);

    dfp->onlyModifyTargeted = NULL;
    if (target != NULL) {
      dfp->onlyModifyTargeted = CheckBox (w, "Only modify targeted record", NULL);
    }

    c = HiddenGroup (w, 4, 0, NULL);
    b = PushButton (c, "Accept", DefLineModFormAcceptProc);
    SetObjectExtra (b, dfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);

    AlignObjects (ALIGN_CENTER, (HANDLE) dfp->addLabels, (HANDLE) dfp->customGrp,
                  (HANDLE) dfp->sourceListGrp, (HANDLE) q,
                  (HANDLE) dfp->leaveInParentheses, (HANDLE) c,
                  (HANDLE) dfp->onlyModifyTargeted, NULL);

    RealizeWindow (w);
  }
  if (smartMods) {
    DefLineModFormAcceptProc (b);
    return NULL;
  }
  return (ForM) w;
}

extern void GenerateAutomaticDefLinesCommon (IteM i, Boolean addMods, Boolean smartMods, ButtoN b)

{
  EnumFieldAssocPtr  ap;
  BaseFormPtr        bfp;
  Int2               count;
  ForM               f;
  Int2               maxCount;
  SeqEntryPtr        sep;
  BioseqPtr          target;

  if (b != NULL) {
    bfp = GetObjectExtra (b);
  } else {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
  }
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  MemSet ((Pointer) (orgmod_rank), (int)(0), sizeof(orgmod_rank));
  MemSet ((Pointer) (subsource_rank), (int)(0), sizeof(subsource_rank));

  count = 0;
  while (sourceModRankList [count] != NULL) {
    for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
      if (StringICmp (ap->name, sourceModRankList [count]) == 0 &&
          ap->value > 0 && ap->value < 24) {
        orgmod_rank [ap->value] = count + 1;
      }
    }
    for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
      if (StringICmp (ap->name, sourceModRankList [count]) == 0 &&
          ap->value > 0 && ap->value < 24) {
        subsource_rank [ap->value] = count + 1;
      }
    }
    count++;
  }

  if (smartMods) {
    if (CountSeqEntryComponents (sep) == 1) {
      smartMods = FALSE;
    }
  }

  if (addMods || smartMods) {
    target =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
    if (addMods) {
      maxCount = INT2_MAX;
    } else {
      maxCount = 6;
    }
    f = CreateDefLineModForm (bfp->input_entityID, sep, target, smartMods, maxCount);
    if (addMods && f != NULL) {
      Show (f);
      Select (f);
    } else if (bfp->activate != NULL) {
      bfp->activate ((WindoW) bfp->form);
    }
    return;
  }

  WatchCursor ();
  Update ();
  AutoDefProc (bfp->input_entityID, sep, FALSE, FALSE, INT2_MAX, TRUE, NULL, NULL, NULL);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

void GenerateAutoDefLinesNoMods (IteM i)

{
  GenerateAutomaticDefLinesCommon (i, FALSE, FALSE, NULL);
}

void GenerateAutoDefLinesWithMods (IteM i)

{
  GenerateAutomaticDefLinesCommon (i, TRUE, FALSE, NULL);
}

void GenerateAutoDefLinesSmartMods (IteM i)

{
  GenerateAutomaticDefLinesCommon (i, FALSE, TRUE, NULL);
}

static void StringToLower (CharPtr str)

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

static void MakeNucleotideTitlesInSequinStyle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  EnumFieldAssocPtr  ap;
  BioseqContextPtr   bcp;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SubSourcePtr       ssp;
  CharPtr            str;
  Char               text [256];
  ValNodePtr         ttl;
  ValNodePtr         vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  bcp = BioseqContextNew (bsp);
  sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_source, NULL, NULL);
  BioseqContextFree (bcp);
  if (sdp == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  if (bsp->descr != NULL) {
    vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
    vnp = ValNodeFreeData (vnp);
  }
  str = MemNew (2000);

  orp = biop->org;
  if (orp != NULL) {
    StringCpy (text, "[org=");
    StringCat (text, orp->taxname);
    StringCat (text, "] ");
    StringCat (str, text);
  }

  ssp = biop->subtype;
  while (ssp != NULL) {
    for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
      if (ssp->subtype == ap->value) {
        StringCpy (text, "[");
        if (ap->value != 255) {
          StringCat (text, ap->name);
        } else {
          StringCat (text, "subsource");
        }
        StringToLower (text);
        StringCat (text, "=");
        StringCat (text, ssp->name);
        StringCat (text, "] ");
        StringCat (str, text);
      }
    }
    ssp = ssp->next;
  }
  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      mod = onp->mod;
      while (mod != NULL) {
        for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
          if (mod->subtype == ap->value) {
            StringCpy (text, "[");
            if (ap->value != 255) {
              StringCat (text, ap->name);
            } else {
              StringCat (text, "note");
            }
            StringToLower (text);
            StringCat (text, "=");
            StringCat (text, mod->subname);
            StringCat (text, "] ");
            StringCat (str, text);
          }
        }
        mod = mod->next;
      }
    }
  }

  TrimSpacesAroundString (str);
  if (! StringHasNoText (str)) {
    ttl = CreateNewDescriptor (sep, Seq_descr_title);
    if (ttl != NULL) {
      ttl->data.ptrvalue = StringSave (str);
    }
  }
  MemFree (str);
}

static Int2 LIBCALLBACK MakeSequinNucleotideTitles (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  SeqEntryExplore (sep, NULL, MakeNucleotideTitlesInSequinStyle);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static void MakeProteinTitlesInSequinStyle (Uint2 entityID, SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  DefFeatsPtr   dfp;
  GeneRefPtr    grp;
  GatherScope   gs;
  ValNodePtr    head;
  SeqEntryPtr   nsep;
  ProtRefPtr    prp;
  SeqEntryPtr   psep;
  Char          str [256];
  Char          text [256];
  ValNodePtr    ttl;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        MakeProteinTitlesInSequinStyle (entityID, sep);
      }
      return;
    }
  }
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  gs.target = NULL;
  head = NULL;
  GatherEntity (entityID, (Pointer) (&head), GetCDStRNArRNAGatherFunc, &gs);
  /* head = SortValNode (head, SortCDStRNArRNAByLocation); */
  if (head == NULL) return;

  vnp = head;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL && dfp->sfp != NULL && dfp->subtype == FEATDEF_CDS) {
      FindGeneAndProtForCDS (entityID, dfp->sfp, &(dfp->gene), &(dfp->prot));
      bsp = GetBioseqGivenSeqLoc (dfp->sfp->product, entityID);
      if (bsp != NULL) {
        str [0] = '\0';
        if (dfp->gene != NULL) {
          grp = (GeneRefPtr) dfp->gene->data.value.ptrvalue;
          if (grp != NULL) {
            StringNCpy_0 (text, (CharPtr) grp->locus, sizeof (text));
            if (! StringHasNoText (text)) {
              StringCat (str, "[gene=");
              StringCat (str, text);
              StringCat (str, "]");
            }
          }
        }
        if (dfp->prot != NULL) {
          prp = (ProtRefPtr) dfp->prot->data.value.ptrvalue;
          if (prp != NULL && prp->name != NULL) {
            StringNCpy_0 (text, (CharPtr) prp->name->data.ptrvalue, sizeof (text));
            if (! StringHasNoText (text)) {
              if (str [0] != '\0') {
                StringCat (str, " ");
              }
              StringCat (str, "[prot=");
              StringCat (str, text);
              StringCat (str, "]");
            }
          }
        }
        if (! StringHasNoText (str)) {
          psep = SeqMgrGetSeqEntryForData (bsp);
          if (psep != NULL) {
            ttl = CreateNewDescriptor (psep, Seq_descr_title);
            if (ttl != NULL) {
              ttl->data.ptrvalue = StringSave (str);
            }
          }
        }
      }
    }
    vnp = vnp->next;
  }
  ValNodeFreeData (head);
}

static Int2 LIBCALLBACK MakeSequinProteinTitles (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  MakeProteinTitlesInSequinStyle (ompcp->input_entityID, sep);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static void PrintFtableIntervals (FILE *fp, SeqLocPtr location, Uint2 entityID, CharPtr label)

{
  BioseqPtr  bsp;
  SeqLocPtr  slp;
  Int4       start;
  Int4       stop;

  if (fp == NULL || location == NULL || StringHasNoText (label)) return;
  bsp = GetBioseqGivenSeqLoc (location, entityID);
  if (bsp == NULL) return;
  slp = SeqLocFindNext (location, NULL);
  if (slp == NULL) return;
  start = GetOffsetInBioseq (slp, bsp, SEQLOC_START) + 1;
  stop = GetOffsetInBioseq (slp, bsp, SEQLOC_STOP) + 1;
  fprintf (fp, "%ld\t%ld\t%s\n", (long) start, (long) stop, label);
  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    start = GetOffsetInBioseq (slp, bsp, SEQLOC_START) + 1;
    stop = GetOffsetInBioseq (slp, bsp, SEQLOC_STOP) + 1;
    fprintf (fp, "%ld\t%ld\n", (long) start, (long) stop);
  }
}

static void MakeFeatureInSequinStyle (Uint2 entityID, SeqEntryPtr sep, FILE *fp)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  DefFeatsPtr   dfp;
  GeneRefPtr    grp;
  GatherScope   gs;
  ValNodePtr    head;
  SeqEntryPtr   nsep;
  ProtRefPtr    prp;
  RnaRefPtr     rrp;
  SeqIdPtr      sip;
  Char          str [256];
  ValNodePtr    syn;
  ValNodePtr    vnp;

  if (sep == NULL || fp == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        MakeFeatureInSequinStyle (entityID, sep, fp);
      }
      return;
    }
  }
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;
  bsp = nsep->data.ptrvalue;
  if (bsp == NULL) return;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  gs.target = NULL;
  head = NULL;
  GatherEntity (entityID, (Pointer) (&head), GetGeneCDStRNArRNAmRNAGatherFunc, &gs);
  head = SortValNode (head, SortCDStRNArRNAByLocation);
  if (head == NULL) return;

  sip = SeqIdFindBest (bsp->id, 0);
  SeqIdWrite (sip, str, PRINTID_FASTA_LONG, sizeof (str) - 1);
  if (! StringHasNoText (str)) {
    fprintf (fp, ">Feature %s\n", str);
  }

  vnp = head;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL && dfp->sfp != NULL) {
      str [0] = '\0';
      switch (dfp->subtype) {
        case FEATDEF_GENE :
          PrintFtableIntervals (fp, dfp->sfp->location, entityID, "gene");
          grp = (GeneRefPtr) dfp->sfp->data.value.ptrvalue;
          if (grp != NULL) {
            StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tgene\t%s\n", str);
            }
            for (syn = grp->syn; syn != NULL; syn = syn->next) {
              StringNCpy_0 (str, (CharPtr) syn->data.ptrvalue, sizeof (str));
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tgene\t%s\n", str);
              }
            }
          }
          break;
        case FEATDEF_CDS :
          PrintFtableIntervals (fp, dfp->sfp->location, entityID, "CDS");
          FindGeneAndProtForCDS (entityID, dfp->sfp, &(dfp->gene), &(dfp->prot));
          if (dfp->gene != NULL) {
            grp = (GeneRefPtr) dfp->gene->data.value.ptrvalue;
            if (grp != NULL) {
              StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tgene\t%s\n", str);
              }
            }
          }
          if (dfp->prot != NULL) {
            prp = (ProtRefPtr) dfp->prot->data.value.ptrvalue;
            if (prp != NULL) {
              if (prp->name != NULL) {
                StringNCpy_0 (str, (CharPtr) prp->name->data.ptrvalue, sizeof (str));
              } else if (prp->desc != NULL) {
                StringNCpy_0 (str, prp->desc, sizeof (str));
              }
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tproduct\t%s\n", str);
              }
            }
          }
          break;
        case FEATDEF_mRNA :
          PrintFtableIntervals (fp, dfp->sfp->location, entityID, "mRNA");
          FindGeneAndProtForCDS (entityID, dfp->sfp, &(dfp->gene), NULL);
          rrp = (RnaRefPtr) dfp->sfp->data.value.ptrvalue;
          if (rrp != NULL) {
            if (rrp->ext.choice == 1) {
              StringNCpy_0 (str, (CharPtr) rrp->ext.value.ptrvalue, sizeof (str));
            }
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tproduct\t%s\n", str);
            }
          }
          break;
      }
    }
    vnp = vnp->next;
  }
  ValNodeFreeData (head);
}

static Boolean LIBCALLBACK SequinFTableFeature (SeqFeatPtr sfp, SeqMgrFeatContextPtr context)

{
  DbtagPtr     dbt;
  FILE         *fp;
  GBQualPtr    gbq;
  GeneRefPtr   grp;
  CharPtr      label;
  ObjectIdPtr  oip;
  BioseqPtr    prod;
  SeqFeatPtr   prot;
  ProtRefPtr   prp;
  RnaRefPtr    rrp;
  Char         str [256];
  tRNAPtr      trna;
  ValNodePtr   vnp;

  if (sfp == NULL) return TRUE;
  fp = (FILE *) context->userdata;
  label = (CharPtr) FeatDefTypeLabel (sfp);
  if (StringCmp (label, "Gene") == 0) {
    label = "gene";
  }
  PrintFtableIntervals (fp, sfp->location, context->entityID, label);
  switch (context->seqfeattype) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
        if (! StringHasNoText (str)) {
          fprintf (fp, "\t\t\tgene\t%s\n", str);
        }
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
          if (! StringHasNoText (str)) {
            fprintf (fp, "\t\t\tgene\t%s\n", str);
          }
        }
      }
      break;
    case SEQFEAT_CDREGION :
      prod = BioseqFind (SeqLocId (sfp->product));
      prot = SeqMgrGetBestProteinFeature (prod, NULL);
      if (prot != NULL) {
        prp = (ProtRefPtr) prot->data.value.ptrvalue;
        if (prp != NULL) {
          if (prp->name != NULL) {
            for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
              StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tproduct\t%s\n", str);
              }
            }
          } else if (prp->desc != NULL) {
            StringNCpy_0 (str, prp->desc, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tproduct\t%s\n", str);
            }
          }
          for (vnp = prp->activity; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tfunction\t%s\n", str);
            }
          }
          for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tEC_number\t%s\n", str);
            }
          }
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL) {
        switch (rrp->ext.choice) {
          case 1 :
            StringNCpy_0 (str, (CharPtr) rrp->ext.value.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tproduct\t%s\n", str);
            }
            break;
          case 2 :
            trna = rrp->ext.value.ptrvalue;
            if (trna != NULL) {
              FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_CONTENT);
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tproduct\t%s\n", str);
              }
            }
            break;
          default :
            break;
        }
      }
      break;
    default :
      break;
  }
  if (! StringHasNoText (sfp->comment)) {
    fprintf (fp, "\t\t\tnote\t%s\n", sfp->comment);
  }
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (! StringHasNoText (gbq->qual)) {
      if (! StringHasNoText (gbq->val)) {
        fprintf (fp, "\t\t\t%s\t%s\n", gbq->qual, gbq->val);
      }
    }
  }
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (! StringHasNoText (dbt->db)) {
        oip = dbt->tag;
        if (oip->str != NULL && (! StringHasNoText (oip->str))) {
          fprintf (fp, "\t\t\tdb_xref\t%s:%s\n", dbt->db, oip->str);
        } else {
          fprintf (fp, "\t\t\tdb_xref\t%s:%ld\n", dbt->db, (long) oip->id);
        }
      }
    }
  }
  return TRUE;
}

static Boolean LIBCALLBACK SequinFTableBioseq (BioseqPtr bsp, SeqMgrBioseqContextPtr context)

{
  FILE      *fp;
  SeqIdPtr  sip;
  Char      str [42];

  if (bsp == NULL) return TRUE;
  fp = (FILE *) context->userdata;
  sip = SeqIdFindBest (bsp->id, 0);
  SeqIdWrite (sip, str, PRINTID_FASTA_LONG, sizeof (str) - 1);
  if (! StringHasNoText (str)) {
    fprintf (fp, ">Feature %s\n", str);
  }
  return SeqMgrExploreFeatures (bsp, (Pointer) fp, SequinFTableFeature, NULL, NULL, NULL);
}

static Int2 LIBCALLBACK MakeSequinFeatureTable (Pointer data)

{
  FILE              *fp;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return OM_MSG_RET_ERROR;
  if (SeqMgrFeaturesAreIndexed (ompcp->input_entityID) != 0) {
    SeqMgrExploreBioseqs (ompcp->input_entityID, NULL, (Pointer) fp, SequinFTableBioseq, TRUE, FALSE, FALSE);
  } else {
    MakeFeatureInSequinStyle (ompcp->input_entityID, sep, fp);
  }
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Gene - CDS Feature Table");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}

/* the following two functions are modified from PrintGenome in wprint.c */

static void SqLitPrintGenome(SeqLitPtr slp, FILE *fp)
{
	static Char		val[166];

	if (slp->seq_data != NULL)         /* not a gap */
	{
		if (slp->length == 0)  /* unknown length */
		{
			fprintf(fp, "gap\t0\t0\t0\n");
		} else {
/* don't know what to do here */
		}
	} else {                  /* gap length was set */
			fprintf(fp, "gap\t0\t0\t%ld\n", (long) slp->length);
	}
}

static void SqLocPrintGenome(SeqLocPtr slp_head, FILE *fp)
{
	SeqLocPtr	slp;
	static Char		buf[11];
	SeqIdPtr	sid, newid;
	Int4 		start, stop;
		
	for (slp = slp_head; slp; slp = slp->next) {
		sid = SeqLocId(slp);
		if (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE) {
			start = SeqLocStart(slp);
			stop = SeqLocStop(slp);
		} else if (slp->choice == SEQLOC_NULL){
			fprintf(fp, "gap\t0\t0\t0\n");
			continue;
		} else {
			continue;
		}
		if (sid->choice == SEQID_GI) {
			newid = GetSeqIdForGI(sid->data.intvalue);
		} else if (sid->choice == SEQID_GENERAL) {
			newid = sid;
		} else {
			newid = sid;
		}
		SeqIdWrite(newid, buf, PRINTID_TEXTID_ACCESSION, 10);
		fprintf(fp, "%s\t%ld\t%ld\n", buf, (long) start+1, (long) stop+1);
	}
}

static Boolean DeltaLitOnly (BioseqPtr bsp)

{
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return FALSE;
  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) return FALSE;
  }
  return TRUE;
}

static Int2 LIBCALLBACK MakeContigBuildTable (Pointer data)

{
  BioseqPtr         bsp;
  DeltaSeqPtr       dsp;
  FILE              *fp;
  SeqLitPtr 	    litp;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  bsp = NULL;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  sep = GetBestTopParentForData (ompcp->input_entityID, bsp);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  if (IsAGenomeRecord (sep) ||
      IsSegmentedBioseqWithoutParts (sep)) {
  } else if (IsADeltaBioseq (sep) && (! DeltaLitOnly (bsp))) {
  } else return OM_MSG_RET_ERROR;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return OM_MSG_RET_ERROR;
/* the following code is modified from PrintGenome in wprint.c */
	if (bsp->seq_ext_type == 1) {
		SqLocPrintGenome((SeqLocPtr) bsp->seq_ext, fp);
	} else if (bsp->seq_ext_type == 4) {
		for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp; dsp=dsp->next) {
			if (dsp->choice == 1) {  /* SeqLoc */
				SqLocPrintGenome((SeqLocPtr)(dsp->data.ptrvalue), fp);
			} else {
				litp = (SeqLitPtr)(dsp->data.ptrvalue);
				if (litp == NULL) continue;
				SqLitPrintGenome(litp, fp);
			}
		}
	}
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Contig Build Table");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}
