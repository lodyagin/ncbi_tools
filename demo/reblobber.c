/*   reblobber.c
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
* File Name:  reblobber.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   Dec. 2, 2011
*
* $Revision: 1.7 $
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

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <sequtil.h>
#include <salutil.h>
#include <edutil.h>
#include <seqport.h>
#include <gather.h>
#include <sqnutils.h>
#include <subutil.h>
#include <toasn3.h>
#include <valid.h>
#include <asn2gnbk.h>
#include <explore.h>
#include <tofasta.h>
#include <simple.h>
#include <suggslp.h>
#include <toporg.h>
#include <aliparse.h>
#include <util/creaders/alnread.h>
#include <pmfapi.h>
#include <tax3api.h>
#ifdef INTERNAL_NCBI_TBL2ASN
#include <accpubseq.h>
#endif
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

#define REBLOBBER_APP_VER "1.0"

CharPtr REBLOBBER_APPLICATION = REBLOBBER_APP_VER;

typedef struct cleanupargs {
  Boolean collection_dates;
  Boolean collection_dates_month_first;
  Boolean add_notes_to_overlapping_cds_without_abc;
  Boolean extend_partial_features_to_gaps_or_ends;
  Boolean add_exception_to_nonextendable_partials;
  Boolean add_exception_to_short_introns;
  FILE *  cleanup_log;
} CleanupArgsData, PNTR CleanupArgsPtr;

typedef struct tblargs {
  Boolean     raw2delt;
  Int2        r2dmin;
  Boolean     r2dunk100;
  Boolean     fastaset;
  Int2        whichclass;
  Boolean     deltaset;
  Boolean     alignset;
  Boolean     gapped;
  Boolean     phrapace;
  Boolean     ftable;
  Boolean     genprodset;
  Boolean     delaygenprodset;
  Boolean     linkbyoverlap;
  Boolean     linkbyproduct;
  Boolean     implicitgaps;
  Boolean     forcelocalid;
  Boolean     gpstonps;
  Boolean     gnltonote;
  Boolean     removeunnecxref;
  Boolean     dotaxlookup;
  Boolean     dopublookup;
  CharPtr     accn;
  CharPtr     center;
  Int4        project_version;
  CharPtr     organism;
  CharPtr     srcquals;
  CharPtr     comment;
  CharPtr     commentFile;
  CharPtr     tableFile;
  Boolean     findorf;
  Boolean     runonorf;
  Boolean     altstart;
  Boolean     conflict;
  Boolean     validate;
  Boolean     relaxed;
  Boolean     validate_barcode;
  Boolean     flatfile;
  Boolean     genereport;
  Boolean     seqidfromfile;
  Boolean     smartfeats;
  Boolean     refSeqTitles;
  Boolean     smarttitle;
  Boolean     logtoterminal;
  CharPtr     aln_beginning_gap;
  CharPtr     aln_end_gap;
  CharPtr     aln_middle_gap;
  CharPtr     aln_missing;
  CharPtr     aln_match;
  Boolean     aln_is_protein;
  Boolean     save_bioseq_set;
  Boolean     auto_def;
  Boolean     apply_cmt_to_all;
  Boolean     adjust_mrna_for_cds_stop_codon;
  Boolean     seq_fetch_failure;

  GlobalDiscrepReportPtr global_report;

  CleanupArgsData cleanup_args;
} TblArgs, PNTR TblArgsPtr;


static void RetryRead (CharPtr filename)
{
  AsnIoPtr aip;
  SeqEntryPtr orig_sep;
  FILE *fp;

  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open %s", filename);
  }
  FileClose (fp);
  
  /* for debugging purposes, try to read it again */
  aip = AsnIoOpen (filename, "r");
  orig_sep = SeqEntryAsnRead (aip, NULL);
  AsnIoClose (aip);
}


/* source information for several common organisms sequenced by genome centers */


typedef struct gcmdata {
  SeqFeatPtr  gene;
  SeqFeatPtr  feat;
  CharPtr     label;
} GmcData, PNTR GmcDataPtr;

static int LIBCALLBACK SortByGenePtr (
  VoidPtr vp1,
  VoidPtr vp2
)

{
  GmcDataPtr gdp1, gdp2;

  if (vp1 == NULL || vp2 == NULL) return 0;
  gdp1 = (GmcDataPtr) vp1;
  gdp2 = (GmcDataPtr) vp2;
  if (gdp1 == NULL || gdp2 == NULL) return 0;

  if (gdp1->gene > gdp2->gene) return -1;
  if (gdp1->gene < gdp2->gene) return 1;

  if (gdp1->feat > gdp2->feat) return -1;
  if (gdp1->feat < gdp2->feat) return 1;

  return 0;
}



static void EnhanceOneCDS (
  SeqFeatPtr sfp,
  Boolean alt_splice
)

{
  DbtagPtr        dbt;
  GBQualPtr       gbq;
  Char            id [64];
  SeqIdPtr        ids, sip;
  size_t          len;
  CharPtr         name, nwstr, ptr, str;
  ObjectIdPtr     oip;
  ProtRefPtr      prp;
  Char            tmp [256];
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;

  name = NULL;
  vnp = NULL;
  prp = NULL;

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->data.choice == SEQFEAT_PROT) {
      prp = (ProtRefPtr) xref->data.value.ptrvalue;
    }
  }

  id [0] = '\0';
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "protein_id") == 0) {
      StringNCpy_0 (id, gbq->val, sizeof (id));
    }
  }
  if (StringDoesHaveText (id) && StringChr (id, '|') != NULL) {
    str = NULL;
    ids = SeqIdParse (id);
    for (sip = ids; sip != NULL; sip = sip->next) {
      if (sip->choice != SEQID_GENERAL) continue;
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt == NULL) continue;
      if (IsSkippableDbtag (dbt)) continue;
      oip = dbt->tag;
      if (oip == NULL) continue;
      str = oip->str;
    }

    if (StringDoesHaveText (str)) {
      if (prp != NULL && prp->name != NULL) {
        vnp = prp->name;
        name = (CharPtr) vnp->data.ptrvalue;
      }
      if (StringDoesHaveText (name) && vnp != NULL) {
        if (alt_splice) {
          ptr = StringChr (str, '-');
          if (ptr != NULL && StringLen (ptr) == 3) {
            ptr++;
            ptr++;
            sprintf (tmp, "%s, isoform %s", str, ptr);
            len = StringLen (name) + StringLen (", ") + StringLen (tmp);
            nwstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
            if (nwstr != NULL) {
              StringCpy (nwstr, name);
              /*
              StringCat (nwstr, ", ");
              */
              StringCat (nwstr, " ");
              StringCat (nwstr, tmp);
              vnp->data.ptrvalue = (Pointer) nwstr;
              MemFree (name);
            }
          } else {
            AddQualifierToFeature (sfp, "product", str);
          }
        } else {
          len = StringLen (name) + StringLen (", ") + StringLen (str);
          nwstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
          if (nwstr != NULL) {
            StringCpy (nwstr, name);
            /*
            StringCat (nwstr, ", ");
            */
            StringCat (nwstr, " ");
            StringCat (nwstr, str);
            vnp->data.ptrvalue = (Pointer) nwstr;
            MemFree (name);
          }
        }
      } else {
        if (alt_splice) {
          ptr = StringChr (str, '-');
          if (ptr != NULL && StringLen (ptr) == 3) {
            ptr++;
            ptr++;
            sprintf (tmp, "%s, isoform %s", str, ptr);
            AddQualifierToFeature (sfp, "product", tmp);
          } else {
            AddQualifierToFeature (sfp, "product", str);
          }
        } else {
          AddQualifierToFeature (sfp, "product", str);
        }
      }
    }

    SeqIdSetFree (ids);
  }
}

static void EnhanceOneRna (
  SeqFeatPtr sfp,
  Boolean alt_splice
)

{
  DbtagPtr     dbt;
  GBQualPtr    gbq, nm_gbq;
  Char         id [64];
  SeqIdPtr     ids, sip;
  size_t       len;
  CharPtr      name, nwstr, ptr, str;
  ObjectIdPtr  oip;
  RnaRefPtr    rrp;
  Char         tmp [256];

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;

  name = NULL;
  nm_gbq = NULL;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp != NULL && rrp->ext.choice == 1) {
    switch (rrp->type) {
      case 1 :  /* precurrsor_RNA */
      case 2 :  /* mRNA */
      case 4 :  /* rRNA */
        name = rrp->ext.value.ptrvalue;
        break;
      case 255 :  /* misc_RNA, ncRNA, tmRNA */
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringICmp (gbq->qual, "product") == 0) {
            nm_gbq = gbq;
            name = gbq->val;
          }
        }
        break;
      case 3:  /* tRNA */
        return;
      default :
        break;
    }
  }

  id [0] = '\0';
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "transcript_id") == 0) {
      StringNCpy_0 (id, gbq->val, sizeof (id));
    }
  }
  if (StringDoesHaveText (id) && StringChr (id, '|') != NULL) {
    str = NULL;
    ids = SeqIdParse (id);
    for (sip = ids; sip != NULL; sip = sip->next) {
      if (sip->choice != SEQID_GENERAL) continue;
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt == NULL) continue;
      if (IsSkippableDbtag(dbt)) continue;
      oip = dbt->tag;
      if (oip == NULL) continue;
      str = oip->str;
    }

    if (StringDoesHaveText (str)) {
      if (StringDoesHaveText (name) && StringCmp (str, name) != 0) {
        if (alt_splice) {
          ptr = StringChr (str, '-');
          if (ptr != NULL && StringLen (ptr) == 3) {
            ptr++;
            ptr++;
            sprintf (tmp, "%s, transcript variant %s", str, ptr);
            len = StringLen (name) + StringLen (", ") + StringLen (tmp);
            nwstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
            if (nwstr != NULL) {
              StringCpy (nwstr, name);
              /*
              StringCat (nwstr, ", ");
              */
              StringCat (nwstr, " ");
              StringCat (nwstr, tmp);
              if (nm_gbq != NULL) {
                nm_gbq->val = (Pointer) nwstr;
              } else if (rrp != NULL) {
                rrp->ext.value.ptrvalue = (Pointer) nwstr;
              }
              MemFree (name);
            }
          } else {
            AddQualifierToFeature (sfp, "product", str);
          }
        } else {
          len = StringLen (name) + StringLen (", ") + StringLen (str);
          nwstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
          if (nwstr != NULL) {
            StringCpy (nwstr, name);
            /*
            StringCat (nwstr, ", ");
            */
            StringCat (nwstr, " ");
            StringCat (nwstr, str);
            if (nm_gbq != NULL) {
              nm_gbq->val = (Pointer) nwstr;
            } else if (rrp != NULL) {
              rrp->ext.value.ptrvalue = (Pointer) nwstr;
            }
            MemFree (name);
          }
        }
      } else {
        if (alt_splice) {
          ptr = StringChr (str, '-');
          if (ptr != NULL && StringLen (ptr) == 3) {
            ptr++;
            ptr++;
            sprintf (tmp, "%s, transcript variant %s", str, ptr);
            AddQualifierToFeature (sfp, "product", tmp);
          } else {
            AddQualifierToFeature (sfp, "product", str);
          }
        } else {
          AddQualifierToFeature (sfp, "product", str);
        }
      }
    }

    SeqIdSetFree (ids);
  }
}

static void EnhanceFeatureAnnotation (
  SeqFeatPtr features,
  BioseqPtr bsp
)

{
  GmcDataPtr  gdp, head;
  GeneRefPtr  grp;
  Int2        i, j, k, numgene, numcds, numrna;
  SeqFeatPtr  sfp;

  if (features == NULL || bsp == NULL) return;

  numgene = 0;
  numcds = 0;
  numrna = 0;

  for (sfp = features; sfp != NULL; sfp = sfp->next) {
    switch (sfp->data.choice) {
      case SEQFEAT_GENE :
        numgene++;
        break;
      case SEQFEAT_CDREGION :
        numcds++;
        break;
      case SEQFEAT_RNA :
        numrna++;
        break;
      default :
        break;
    }
  }

  if (numgene == 0) return;

  if (numcds > 0) {
    head = (GmcDataPtr) MemNew (sizeof (GmcData) * (size_t) (numcds + 1));
    if (head != NULL) {
      gdp = head;
      for (sfp = features; sfp != NULL; sfp = sfp->next) {
        if (sfp->idx.subtype == FEATDEF_CDS) {
          gdp->feat = sfp;
          grp = SeqMgrGetGeneXref (sfp);
          if (grp == NULL || (! SeqMgrGeneIsSuppressed (grp))) {
            gdp->gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          }
          gdp++;
        }
      }
      HeapSort (head, (size_t) numcds, sizeof (GmcData), SortByGenePtr);
      for (i = 0; i < numcds; i += j) {
        sfp = head [i].gene;
        for (j = 1; i + j < numcds && sfp == head [i + j].gene; j++) continue;
        if (j == 1) {
          /* no alt splicing */
          EnhanceOneCDS (head [i].feat, FALSE);
        } else {
          /* is alt splicing */
          for (k = 0; k < j; k++) {
            EnhanceOneCDS (head [i + k].feat, TRUE);
          }
        }
      }
    }
    MemFree (head);
  }

  if (numrna > 0) {
    head = (GmcDataPtr) MemNew (sizeof (GmcData) * (size_t) (numrna + 1));
    if (head != NULL) {
      gdp = head;
      for (sfp = features; sfp != NULL; sfp = sfp->next) {
        if (sfp->data.choice == SEQFEAT_RNA) {
          gdp->feat = sfp;
          grp = SeqMgrGetGeneXref (sfp);
          if (grp == NULL || (! SeqMgrGeneIsSuppressed (grp))) {
            gdp->gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          }
          gdp++;
        }
      }
      HeapSort (head, (size_t) numrna, sizeof (GmcData), SortByGenePtr);
      for (i = 0; i < numrna; i += j) {
        sfp = head [i].gene;
        for (j = 1; i + j < numrna && sfp == head [i + j].gene; j++) continue;
        if (j == 1) {
          /* no alt splicing */
          EnhanceOneRna (head [i].feat, FALSE);
        } else {
          /* is alt splicing */
          for (k = 0; k < j; k++) {
            EnhanceOneRna (head [i + k].feat, TRUE);
          }
        }
      }
    }
    MemFree (head);
  }
}

static BioseqPtr AttachSeqAnnotEntity (
  Uint2 entityID,
  SeqAnnotPtr sap,
  TblArgsPtr tbl
)

{
  SeqAnnotPtr  anp;
  BioseqPtr    bsp;
  Char         buf [80];
  Int2         genCode;
  SeqEntryPtr  oldscope;
  SeqEntryPtr  sep;
  SeqFeatPtr   sfp = NULL;
  SeqIdPtr     sip;
  SeqLocPtr    slp;

  if (sap == NULL || tbl == NULL) return NULL;

  bsp = GetBioseqReferencedByAnnot (sap, entityID);
  if (bsp == NULL) {
    oldscope = SeqEntrySetScope (NULL);
    if (oldscope != NULL) {
      bsp = GetBioseqReferencedByAnnot (sap, entityID);
      SeqEntrySetScope (oldscope);
    }
  }
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    entityID = ObjMgrGetEntityIDForChoice (sep);
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      genCode = GetGenCodeForBsp (bsp);
      SetEmptyGeneticCodes (sap, genCode);
    }
    if (bsp->annot == NULL) {
      bsp->annot = sap;
    } else {
      anp = bsp->annot;
      while (anp->next != NULL) {
        anp = anp->next;
      }
      anp->next = sap;
    }
    if (sfp != NULL) {
      if (tbl->smartfeats) {

        /* indexing needed to find mRNA and CDS within each gene */

        SeqMgrIndexFeatures (entityID, NULL);

        EnhanceFeatureAnnotation (sfp, bsp);
      }

      PromoteXrefsExEx (sfp, bsp, entityID, TRUE, FALSE, tbl->genprodset, tbl->forcelocalid, &(tbl->seq_fetch_failure));
      sep = GetTopSeqEntryForEntityID (entityID);
    }
  } else {
    buf [0] = '\0';
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      if (sfp != NULL && sfp->location != NULL) {
        slp = SeqLocFindNext (sfp->location, NULL);
        if (slp != NULL) {
          sip = SeqLocId (slp);
          if (sip != NULL) {
            SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
          }
        }
      }
    }
//    Message (MSG_POSTERR, "Feature table identifiers %s do not match record", buf);
  }
  sep = GetTopSeqEntryForEntityID (entityID);
  return bsp;
}




static void ProcessOneAnnot (
  SeqAnnotPtr sap,
  Uint2 entityID,
  TblArgsPtr tbl
)

{
  BioseqPtr   bsp;
  GBQualPtr   gbq;
  Int2        genCode;
  SeqFeatPtr  sfp;
  SeqEntryPtr sep;

  if (sap == NULL || tbl == NULL) return;

  if (tbl->delaygenprodset) {
    if (sap->type == 1) {
      for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
        for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
          if (StringICmp (gbq->qual, "protein_id") == 0 && sfp->data.choice == SEQFEAT_RNA) {
            gbq->qual = MemFree (gbq->qual);
            gbq->qual = StringSave ("orig_protein_id");
          }
          if (StringICmp (gbq->qual, "transcript_id") == 0) {
            gbq->qual = MemFree (gbq->qual);
            gbq->qual = StringSave ("orig_transcript_id");
          }
        }
      }
    }
  }

  bsp = AttachSeqAnnotEntity (entityID, sap, tbl);
  if (bsp == NULL) return;

  sep = GetTopSeqEntryForEntityID (entityID);

  /* correct all idx parent pointers */

  AssignIDsInEntity (entityID, 0, NULL);

  genCode = GetGenCodeForBsp (bsp);

  /* coercion of SeqIds to accession moved to ProcessOneRecord->MakeAccessionID */

  /* for parsed in features or best ORF, promote CDS products to protein bioseq */

  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 1) {
      SetEmptyGeneticCodes (sap, genCode);
      sfp = (SeqFeatPtr) sap->data;
      PromoteXrefsExEx (sfp, bsp, entityID, TRUE, FALSE, tbl->genprodset, tbl->forcelocalid, &(tbl->seq_fetch_failure));
    }
  }
  sep = GetTopSeqEntryForEntityID (entityID);
  move_cds_ex (sep, FALSE);
}





typedef struct raw2deltdata {
  Uint2           entityID;
  SubmitBlockPtr  sbp;
  BioSourcePtr    src;
  TblArgsPtr      tbl;
  MolInfoPtr      template_molinfo;
} Raw2DeltData, PNTR Raw2DeltPtr;




typedef struct resqseqgph {
  Int2         index;
  SeqGraphPtr  sgp;
} ResqSeqgph, PNTR ResqSeqgphPtr;





typedef struct reqcontig {
  Int2  index;
  Char  str [41];
} ResqContig, PNTR ResqContigPtr;

#define MAX_FIELDS  8





static Boolean DoSequenceLengthsMatch (
  TAlignmentFilePtr afp
)

{
  int    seq_index;
  Int4   seq_len;

  if (afp == NULL || afp->sequences == NULL || afp->num_sequences == 0) {
    return TRUE;
  }
  seq_len = StringLen (afp->sequences[0]);
  for (seq_index = 1; seq_index < afp->num_sequences; seq_index++) {
    if (StringLen (afp->sequences[seq_index]) != seq_len) {
      return FALSE;
    }
  }
  return TRUE;
}



#ifdef INTERNAL_NCBI_ASNDISC
const PerformDiscrepancyTest taxlookup = CheckTaxNamesAgainstTaxDatabase;
#else
const PerformDiscrepancyTest taxlookup = NULL;
#endif




typedef struct reblobseq {
  SeqIdPtr  sip;
  CharPtr   orig_file;
  CharPtr   feature_file;
  CharPtr   src_file;
  CharPtr   quality_file;
  BioseqPtr bsp;
} ReblobSeqData, PNTR ReblobSeqPtr;


static ReblobSeqPtr ReblobSeqNew (ValNodePtr token_list)
{
  ReblobSeqPtr rs;
  Char         id_buf[PATH_MAX];

  rs = (ReblobSeqPtr) MemNew (sizeof (ReblobSeqData));
  if (token_list != NULL) {
    sprintf (id_buf, "gb|%s", (CharPtr) token_list->data.ptrvalue);
    rs->sip = MakeSeqID (id_buf);
    if (token_list->next != NULL) {
      rs->orig_file = token_list->next->data.ptrvalue;
      token_list->next->data.ptrvalue = NULL;
      if (token_list->next->next != NULL) {
        rs->feature_file = token_list->next->next->data.ptrvalue;
        token_list->next->next->data.ptrvalue = NULL;
        if (token_list->next->next->next != NULL) {
          rs->src_file = token_list->next->next->next->data.ptrvalue;
          token_list->next->next->next->data.ptrvalue = NULL;
          if (token_list->next->next->next->next != NULL) {
            rs->quality_file = token_list->next->next->next->next->data.ptrvalue;
            token_list->next->next->next->next->data.ptrvalue = NULL;
          }
        }
      }
    }
  }
  return rs;
}


static ReblobSeqPtr ReblobSeqFree (ReblobSeqPtr rs)
{
  if (rs != NULL) {
    rs->sip = SeqIdFree (rs->sip);
    rs->orig_file = MemFree (rs->orig_file);
    rs->feature_file = MemFree (rs->feature_file);
    rs->src_file = MemFree (rs->src_file);
    rs->quality_file = MemFree (rs->quality_file);
    rs = MemFree (rs);
  }
  return rs;
}


static ValNodePtr ReblobSeqListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = ReblobSeqFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


typedef struct reblobgroup {
  CharPtr filename;
  ValNodePtr seq_list;
} ReblobGroupData, PNTR ReblobGroupPtr;


static ReblobGroupPtr ReblobGroupNew (CharPtr filename, ValNodePtr seq_list)
{
  ReblobGroupPtr rg;

  rg = (ReblobGroupPtr) MemNew (sizeof (ReblobGroupData));
  rg->filename = StringSave (filename);
  rg->seq_list = seq_list;
  return rg;
}


static ReblobGroupPtr ReblobGroupFree (ReblobGroupPtr rg)
{
  if (rg != NULL) {
    rg->filename = MemFree (rg->filename);
    rg->seq_list = ReblobSeqListFree (rg->seq_list);
    rg = MemFree (rg);
  }
  return rg;
}


static ValNodePtr ReblobGroupListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;

  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = ReblobGroupFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


const char *nucleotide_alphabet = "ABCDGHKMRSTUVWYabcdghkmrstuvwy";
const char *protein_alphabet = "ABCDEFGHIKLMPQRSTUVWXYZabcdefghiklmpqrstuvwxyz";

extern TSequenceInfoPtr GetDefaultSequenceInfo (void)
{
  TSequenceInfoPtr sequence_info = SequenceInfoNew();

  sequence_info->missing = MemFree (sequence_info->missing);
  sequence_info->missing = StringSave ("?Nn");

  sequence_info->beginning_gap = MemFree (sequence_info->beginning_gap);
  sequence_info->beginning_gap = StringSave ("-.Nn?");

  sequence_info->middle_gap = MemFree (sequence_info->middle_gap);
  sequence_info->middle_gap = StringSave ("-.");

  sequence_info->end_gap = MemFree (sequence_info->end_gap);
  sequence_info->end_gap = StringSave ("-.Nn?");

  sequence_info->match = MemFree (sequence_info->match);
  sequence_info->match = StringSave (":");

  sequence_info->alphabet = nucleotide_alphabet;

  return sequence_info;
}

static void MakeFeatureTableFileIdsGenbank (SeqAnnotPtr sap);

static Boolean AddOneFeatureTable (CharPtr filename, SeqEntryPtr sep, Uint2 entityID)
{
  FILE *fp;
  Pointer      dataptr;
  Uint2        datatype;
  Boolean      failure = FALSE;
  Int4         linenum = 0;
  SeqAnnotPtr  sap;
  TblArgs      tbl_args;
  SeqEntryPtr  orig_scope;


  if (StringHasNoText (filename)) {
    return FALSE;
  }
  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open feature table file %s\n", filename);
    return FALSE;
  }
  orig_scope = SeqEntrySetScope (sep);

  MemSet (&tbl_args, 0, sizeof (TblArgs));

  while ((! failure) && (dataptr = ReadFeatureTableFile (fp, &datatype, NULL, &linenum, &failure, TRUE)) != NULL) {
    if (datatype == OBJ_SEQANNOT) {
      sap = (SeqAnnotPtr) dataptr;
      MakeFeatureTableFileIdsGenbank (sap);
      ProcessOneAnnot (sap, entityID, &tbl_args);

    } else {
      ObjMgrFree (datatype, dataptr);
    }
  }
  FileClose (fp);

  SeqEntrySetScope (orig_scope);

  if (failure) {
    Message (MSG_POSTERR, "Bad feature table at line %ld of file %s", (long) linenum, filename);
    return FALSE;
  } else {
    return TRUE;
  }
}
  

static Boolean AddFeaturesToGroup (ReblobGroupPtr rg, SeqEntryPtr sep, Uint2 entityID, CharPtr file_filter)
{
  ValNodePtr vnp;
  ReblobSeqPtr rs;
  ValNodeBlock file_list;
  Boolean rval = TRUE;

  if (rg == NULL) {
    return FALSE;
  }

  InitValNodeBlock (&file_list, NULL);
 
  for (vnp = rg->seq_list; vnp != NULL; vnp = vnp->next) {
    rs = vnp->data.ptrvalue;
    if (rs != NULL && !StringHasNoText (rs->feature_file)
        && (StringHasNoText (file_filter) || StringCmp (file_filter, rs->orig_file) == 0)) {
      ValNodeAddPointerToEnd (&file_list, 0, rs->feature_file);
    }
  }

  file_list.head = ValNodeSort (file_list.head, SortVnpByString);
  ValNodeUnique (&file_list.head, SortVnpByString, ValNodeFree);
  
  for (vnp = file_list.head; vnp != NULL; vnp = vnp->next) {
    AddOneFeatureTable (vnp->data.ptrvalue, sep, entityID);
  }
  file_list.head = ValNodeFree (file_list.head);
  return rval; 
}


static ValNodePtr 
AutomatchQualsFromTableHeader
(ValNodePtr header,
 Boolean force_first_accession,
 Boolean allow_no_match)
{
  ValNodePtr val;
  ValNodeBlock columns;
  TabColumnConfigPtr t;

  InitValNodeBlock (&columns, NULL);

  for (val = header;
       val != NULL;
       val = val->next) {
    t = TabColumnConfigNew ();
    if (force_first_accession && columns.head == NULL) {
      t->match_type = MatchTypeNew();
      t->match_type->choice = eTableMatchNucID;
      t->match_type->match_location = String_location_equals;
    } else {
      t->field = FieldTypeFromString (val->data.ptrvalue);
      if (t->field == NULL) {
        Message (MSG_POSTERR, "Unrecognized source qualifier: %s", val->data.ptrvalue);
        t = TabColumnConfigFree (t);
        if (!allow_no_match) {
          columns.head = TabColumnConfigListFree (columns.head);
          return NULL; 
        }
      } else {
        t->existing_text = ExistingTextOption_replace_old;
      }
    }
    ValNodeAddPointerToEnd (&columns, 0, t);
  }
  return columns.head;
}


static void PostErrorList (ValNodePtr vnp, CharPtr filename)
{
  while (vnp != NULL) {
    Message (MSG_POSTERR, "%s in %s", (CharPtr) vnp->data.ptrvalue, filename);
    vnp = vnp->next;
  }
}


static void SpecialFixSrcMods (ValNodePtr list)
{
  while (list != NULL) {
    if (StringICmp (list->data.ptrvalue, "note") == 0) {
      list->data.ptrvalue = MemFree (list->data.ptrvalue);
      list->data.ptrvalue = StringSave ("note-subsrc");
    }
    list = list->next;
  }
}


static Boolean AddOneSrcTable (CharPtr filename, SeqEntryPtr sep, Uint2 entityID)
{
  ValNodePtr err_list;
  ValNodePtr columns, obj_table;
  ValNodePtr table;
  FILE *fp;

  if (StringHasNoText (filename)) {
    return FALSE;
  }
  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    Message (MSG_POSTERR, "Unable to open source file %s", filename);
    return FALSE;
  }
  table = ReadTabTableFromFile (fp);
  FileClose (fp);
  err_list = AutoReplaceSpecialCharactersInTabTable (table);
  PostErrorList(err_list, filename); 
  err_list = ValNodeFreeData(err_list);

  if (table == NULL || table->next == NULL || table->data.ptrvalue == NULL) {
    Message (MSG_POSTERR, "Unable to read table from file %s", filename);
    table = FreeTabTable (table);
    return FALSE;
  }
  SpecialFixSrcMods (table->data.ptrvalue);
  columns = AutomatchQualsFromTableHeader(table->data.ptrvalue, TRUE, FALSE);
  if (columns == NULL) {
    Message (MSG_POSTERR, "Unable to read table from file %s", filename);
    table = FreeTabTable (table);
    return FALSE;
  }

  err_list = ValidateTabTableValues (table->next, columns);
  if (err_list != NULL) {
    PostErrorList (err_list, filename);
    err_list = ValNodeFreeData (err_list);
    table = FreeTabTable (table);
    return FALSE;
  }

  obj_table = GetObjectTableForTabTable (sep, table, columns, &err_list);

  err_list = CheckObjTableForRowsThatApplyToTheSameDestination (obj_table);
  if (err_list != NULL) {
    PostErrorList (err_list, filename);
    Message (MSG_POSTERR, "For one or more columns in %s, two or more rows in the table apply to the same object.  Cannot apply table.", filename);
    err_list = ValNodeFreeData (err_list);
    obj_table = FreeObjectTableForTabTable (obj_table);
    DeleteMarkedObjects (entityID, 0, NULL);
    return FALSE;
  }

  err_list = ApplyTableValuesToObjectTable (sep, table, columns, obj_table);
  /* PostErrorList (err_list, filename); */

  err_list = ValNodeFreeData (err_list);
  obj_table = FreeObjectTableForTabTable (obj_table);
  DeleteMarkedObjects (entityID, 0, NULL);
  table = FreeTabTable(table);
  return TRUE;

}


static Boolean ApplySrcTableToGroup (ReblobGroupPtr rg, SeqEntryPtr sep, Uint2 entityID, CharPtr file_filter)
{
  ValNodePtr vnp;
  ReblobSeqPtr rs;
  ValNodeBlock file_list;
  Boolean rval = TRUE;

  if (rg == NULL) {
    return FALSE;
  }

  InitValNodeBlock (&file_list, NULL);
 
  for (vnp = rg->seq_list; vnp != NULL; vnp = vnp->next) {
    rs = vnp->data.ptrvalue;
    if (rs != NULL && !StringHasNoText (rs->src_file)
        && (StringHasNoText (file_filter) || StringCmp (file_filter, rs->orig_file) == 0)) {
      ValNodeAddPointerToEnd (&file_list, 0, rs->src_file);
    }
  }

  file_list.head = ValNodeSort (file_list.head, SortVnpByString);
  ValNodeUnique (&file_list.head, SortVnpByString, ValNodeFree);
  
  for (vnp = file_list.head; vnp != NULL; vnp = vnp->next) {
    rval &= AddOneSrcTable (vnp->data.ptrvalue, sep, entityID);
  }
  file_list.head = ValNodeFree (file_list.head);
  return rval; 
  
}


static SeqAnnotPtr NewGraphSeqAnnot (
  CharPtr name,
  SeqGraphPtr sgp
)

{
  SeqAnnotPtr  sap = NULL;

  if (sgp == NULL) return NULL;
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;

  if (StringDoesHaveText (name)) {
    SeqDescrAddPointer (&(sap->desc), Annot_descr_name, StringSave (name));
  }
  sap->type = 3;
  sap->data = (Pointer) sgp;

  return sap;
}


static SeqIdPtr MakeSeqGraphFileId (CharPtr str)
{
  SeqIdPtr     sip_new;
  TextSeqIdPtr tsip;
  CharPtr      cp;

  tsip = TextSeqIdNew ();
  tsip->accession = StringSave (str);
  if ((cp = StringChr(tsip->accession, '.')) != NULL) {
    tsip->version = atoi (cp + 1);
    *cp = 0;
  }
  sip_new = ValNodeNew (NULL);
  sip_new->choice = SEQID_GENBANK;
  sip_new->data.ptrvalue = tsip;

  return sip_new;
}


static Boolean ApplyOneQualityFile (CharPtr filename, SeqEntryPtr sep, Uint2 entityID)
{
  FILE              *fp;
  FileCache          fc;
  Boolean            goOn = TRUE;
  Char               buf [1024];
  Boolean            nonewline;
  CharPtr            str, ptr;
  SeqIdPtr           sip;
  BioseqPtr          bsp;
  SeqGraphPtr        sgp, lastsgp;
  SeqAnnotPtr        sap;

  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    return FALSE;
  }

  FileCacheSetup (&fc, fp);

  goOn = TRUE;
  while (goOn) {
    str = FileCacheReadLine (&fc, buf, sizeof (buf), &nonewline);
    if (str == NULL) {
      goOn = FALSE;
    } else if (StringDoesHaveText (str)) {
      if (str [0] == '>') {
        ptr = StringChr (str, ' ');
        if (ptr == NULL) {
          ptr = StringChr (str, '\t');
        }
        if (ptr != NULL) {
          *ptr = '\0';
        }
        sip = MakeSeqGraphFileId (str + 1);
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          sgp = ReadPhrapQualityFC (&fc, bsp);
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
                sap->next = NewGraphSeqAnnot ("Phrap Graph", sgp);
              } else {
                bsp->annot = NewGraphSeqAnnot ("Phrap Graph", sgp);
              }
            }
          }
        }
        SeqIdFree (sip);
      }
    }
  }
  FileClose (fp);
  return TRUE;
}


static Boolean ApplyQualityScoresToGroup (ReblobGroupPtr rg, SeqEntryPtr sep, Uint2 entityID, CharPtr file_filter)
{
  ValNodePtr vnp;
  ReblobSeqPtr rs;
  ValNodeBlock file_list;
  Boolean rval = TRUE;

  if (rg == NULL) {
    return FALSE;
  }

  InitValNodeBlock (&file_list, NULL);
 
  for (vnp = rg->seq_list; vnp != NULL; vnp = vnp->next) {
    rs = vnp->data.ptrvalue;
    if (rs != NULL && !StringHasNoText (rs->quality_file)
        && (StringHasNoText (file_filter) || StringCmp (file_filter, rs->orig_file) == 0)) {
      ValNodeAddPointerToEnd (&file_list, 0, rs->quality_file);
    }
  }

  file_list.head = ValNodeSort (file_list.head, SortVnpByString);
  ValNodeUnique (&file_list.head, SortVnpByString, ValNodeFree);
  
  for (vnp = file_list.head; vnp != NULL; vnp = vnp->next) {
    rval &= ApplyOneQualityFile (vnp->data.ptrvalue, sep, entityID);
  }
  file_list.head = ValNodeFree (file_list.head);
  return rval; 
}


static void CopyDescriptors (BioseqPtr orig_bsp, BioseqPtr new_bsp)
{
  SeqMgrDescContext context;
  SeqDescPtr sdp;

  /* only copy title descriptors if they are on the bioseq itself */
  for (sdp = SeqMgrGetNextDescriptor (orig_bsp, NULL, 0, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (orig_bsp, sdp, 0, &context)) {
    if (sdp->choice != Seq_descr_title) {
      ValNodeLink (&(new_bsp->descr),
                   AsnIoMemCopy ((Pointer) sdp,
                                 (AsnReadFunc) SeqDescAsnRead,
                                 (AsnWriteFunc) SeqDescAsnWrite));
    }
  } 

  /* copy titles on bioseq */
  for (sdp = orig_bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice == Seq_descr_title) {
      ValNodeLink (&(new_bsp->descr),
                   AsnIoMemCopy ((Pointer) sdp,
                                 (AsnReadFunc) SeqDescAsnRead,
                                 (AsnWriteFunc) SeqDescAsnWrite));
    }
  } 
      
}


static void 
CopyFeaturesAndProteins 
(Uint2       entityID,
 BioseqPtr   orig_bsp,
 BioseqPtr   new_bsp,
 SeqEntryPtr orig_scope,
 SeqEntryPtr new_scope)
{
  SeqMgrFeatContext context;
  SeqFeatPtr sfp, sfp_cpy;
  BioseqPtr prod_bsp, prod_cpy;
  SeqEntryPtr sep, sep_prod;

  sep = GetBestTopParentForData (entityID, new_bsp);
  SeqEntrySetScope (orig_scope);
  for (sfp = SeqMgrGetNextFeature (orig_bsp, NULL, 0, 0, &context);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (orig_bsp, sfp, 0, 0, &context)) {
    if (sfp->product != NULL && (prod_bsp = BioseqFindFromSeqLoc (sfp->product)) != NULL) {
      prod_cpy = AsnIoMemCopy (prod_bsp, (AsnReadFunc) BioseqAsnRead, (AsnWriteFunc) BioseqAsnWrite);
      sep_prod = SeqEntryNew ();
      sep_prod->choice = 1;
      sep_prod->data.ptrvalue = prod_cpy;
      AddSeqEntryToSeqEntry (sep, sep_prod, TRUE);
      sep = GetBestTopParentForData (entityID, new_bsp);
    }
    sfp_cpy = AsnIoMemCopy (sfp, (AsnReadFunc) SeqFeatAsnRead, (AsnWriteFunc) SeqFeatAsnWrite);
    CreateNewFeature (SeqMgrGetSeqEntryForData (new_bsp), NULL, sfp->data.choice, sfp_cpy);
  }

}


static void CopyExtraIds (BioseqPtr new_bsp, BioseqPtr orig_bsp)
{
  SeqIdPtr sip_new, sip_orig, sip_last;
  Boolean  found;
  TextSeqIdPtr orig_acc, new_acc;

  for (sip_orig = orig_bsp->id; sip_orig != NULL; sip_orig = sip_orig->next) {
    found = FALSE;
    sip_last = NULL;
    for (sip_new = new_bsp->id; sip_new != NULL && !found; sip_new = sip_new->next) {
      if (SeqIdComp (sip_new, sip_orig) == SIC_YES) {
        found = TRUE;
      }
      sip_last = sip_new;
    }
    if (!found) {
      sip_new = AsnIoMemCopy (sip_orig, (AsnReadFunc) SeqIdAsnRead, (AsnWriteFunc) SeqIdAsnWrite);
      if (sip_last == NULL) {
        new_bsp->id = sip_new;
      } else {
        sip_last->next = sip_new;
      }
    } else if (sip_orig->choice == SEQID_GENBANK && sip_last->choice == SEQID_GENBANK) {
      orig_acc = (TextSeqIdPtr) sip_orig->data.ptrvalue;
      new_acc = (TextSeqIdPtr) sip_last->data.ptrvalue;
      if (orig_acc != NULL && new_acc != NULL) {
        new_acc->version = orig_acc->version;
      }
    }
  }
}


static Boolean CopySetType (SeqEntryPtr new_sep, SeqEntryPtr orig_sep)
{
  BioseqSetPtr orig_bssp, new_bssp;

  if (orig_sep != NULL && IS_Bioseq_set (orig_sep) && (orig_bssp = orig_sep->data.ptrvalue) != NULL
     && new_sep != NULL && IS_Bioseq_set (new_sep) && (new_bssp = new_sep->data.ptrvalue) != NULL
     && orig_bssp->_class != BioseqseqSet_class_nuc_prot) {
    new_bssp->_class = orig_bssp->_class;
    return TRUE;
  } else {
    return FALSE;
  }
  
}


static Boolean CopySetTitle (SeqEntryPtr new_sep, SeqEntryPtr orig_sep)
{
  BioseqSetPtr orig_bssp, new_bssp;
  SeqDescPtr sdp;
  Boolean rval = FALSE;

  if (orig_sep != NULL && IS_Bioseq_set (orig_sep) && (orig_bssp = orig_sep->data.ptrvalue) != NULL
     && new_sep != NULL && IS_Bioseq_set (new_sep) && (new_bssp = new_sep->data.ptrvalue) != NULL) {
    for (sdp = orig_bssp->descr; sdp != NULL && sdp->choice != Seq_descr_title; sdp = sdp->next) {
    }
    if (sdp != NULL) {
      ValNodeLink (&(new_bssp->descr),
                   AsnIoMemCopy ((Pointer) sdp,
                                 (AsnReadFunc) SeqDescAsnRead,
                                 (AsnWriteFunc) SeqDescAsnWrite));
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean CopyQualityGraphs (BioseqPtr bsp_new, BioseqPtr bsp_orig)
{
  SeqAnnotPtr sap, sap_last, sap_new;

  if (bsp_new == NULL || bsp_orig == NULL || CompareSequences (bsp_new, bsp_orig, 0) != 0) {
    return FALSE;
  }

  sap_last = bsp_new->annot;
  while (sap_last != NULL && sap_last->next != NULL) {
    sap_last = sap_last->next;
  }

  for (sap = bsp_orig->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 3) {
      sap_new = (SeqAnnotPtr) AsnIoMemCopy (sap, (AsnReadFunc) SeqAnnotAsnRead, (AsnWriteFunc) SeqAnnotAsnWrite);
      if (sap_last == NULL) {
        bsp_new->annot = sap_new; 
      } else {
        sap_last->next = sap_new;
      }
      sap_last = sap_new;
    }
  }
  return TRUE;
}


/* NOTES:
 * if the orig is in a nuc-prot set, need to add all the proteins from the original nuc-prot-set to
 * the new file.
 * also, should look for descriptors on the nuc-prot set and copy those as well.
 */
static Boolean CopyDataFromOriginals (ReblobGroupPtr rg, SeqEntryPtr sep, Uint2 entityID)
{
  ValNodePtr vnp, vnp_r;
  ReblobSeqPtr rs;
  ValNodeBlock file_list;
  Boolean rval = TRUE;
  SeqEntryPtr orig_sep;
  AsnIoPtr       aip;
  BioseqPtr      orig_bsp;
  CharPtr        filename;
  Uint2          orig_entityID;
  Char           buf[PATH_MAX];

  if (rg == NULL) {
    return FALSE;
  }

  InitValNodeBlock (&file_list, NULL);
 
  SeqEntrySetScope (sep);
  SeqMgrIndexFeatures (entityID, NULL);
  for (vnp = rg->seq_list; vnp != NULL; vnp = vnp->next) {
    rs = vnp->data.ptrvalue;
    if (rs != NULL && !StringHasNoText (rs->orig_file)) {
      ValNodeAddPointerToEnd (&file_list, 0, rs->orig_file);
      rs->bsp = BioseqFind (rs->sip);
    }
  }

  file_list.head = ValNodeSort (file_list.head, SortVnpByString);
  ValNodeUnique (&file_list.head, SortVnpByString, ValNodeFree);
  
  for (vnp = file_list.head; vnp != NULL; vnp = vnp->next) {
    filename = vnp->data.ptrvalue;
    aip = AsnIoOpen (filename, "r");
    orig_sep = SeqEntryAsnRead (aip, NULL);
    AsnIoClose (aip);
    if (orig_sep == NULL) {
      Message (MSG_ERROR, "Unable to read Seq-entry from %s", filename);
      continue;
    }
    SeqEntrySetScope (orig_sep);
    orig_entityID = ObjMgrGetEntityIDForChoice (orig_sep);
    AssignIDsInEntity (orig_entityID, 0, NULL);
    SeqMgrIndexFeatures (orig_entityID, NULL);
    CopySetType (sep, orig_sep);
    CopySetTitle (sep, orig_sep);
    for (vnp_r = rg->seq_list; vnp_r != NULL; vnp_r = vnp_r->next) {
      rs = (ReblobSeqPtr) vnp_r->data.ptrvalue;
      if (StringCmp (rs->orig_file, filename) == 0) {
        /* find original Bioseq */
        orig_bsp = BioseqFind (rs->sip);
        if (orig_bsp != NULL && rs->bsp != NULL && orig_bsp != rs->bsp) {
          /* copy over descriptors */
          CopyDescriptors (orig_bsp, rs->bsp);
          
          /* copy over features if not retrieving from feature file */
          if (StringHasNoText (rs->feature_file)) {
            CopyFeaturesAndProteins (entityID, orig_bsp, rs->bsp, orig_sep, sep);
          }
          /* additional Seq-ids */
          CopyExtraIds (rs->bsp, orig_bsp);
          
          /* copy mol */
          rs->bsp->mol = orig_bsp->mol;

#if 0
          /* do not copy quality graphs */
          /* quality graphs */
          CopyQualityGraphs (rs->bsp, orig_bsp);
#endif

          /* anything else? */
          
        } else {
          SeqIdWrite (rs->sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
          Message (MSG_POSTERR, "Unable to find original sequence for %s", buf);
        }
      }
    }
    SeqEntrySetScope (NULL);
    SeqEntryFree (orig_sep);
  }
  SeqEntrySetScope (sep);
  return rval;
}


static void MarkFeaturesForDeletion (SeqFeatPtr sfp, Pointer data)
{
  if (sfp != NULL) {
    sfp->idx.deleteme = TRUE;
  }
}


static Boolean StripFeaturesFromSelectedSequences (ReblobGroupPtr rg, Uint2 orig_entityID, CharPtr filename)
{
  ValNodePtr   vnp_r;
  ReblobSeqPtr rs;
  BioseqPtr    orig_bsp;
  SeqEntryPtr  seq_sep;
  Boolean      rval = TRUE;
  Char         buf[PATH_MAX];

  if (rg == NULL) {
    return FALSE;
  }

  for (vnp_r = rg->seq_list; vnp_r != NULL; vnp_r = vnp_r->next) {
    rs = (ReblobSeqPtr) vnp_r->data.ptrvalue;
    if (StringCmp (rs->orig_file, filename) == 0) {
      /* find original Bioseq */
      orig_bsp = BioseqFind (rs->sip);
      seq_sep = GetBestTopParentForData (orig_entityID, orig_bsp);
      if (seq_sep == NULL) {
        SeqIdWrite (rs->sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
        Message (MSG_POSTERR, "Unable to find original sequence for %s", buf);
        rval = FALSE;
      } else {
        /* strip features if retrieving from feature file */
        if (!StringHasNoText (rs->feature_file)) {
          VisitFeaturesInSep (seq_sep, NULL, MarkFeaturesForDeletion);
        }
      }
    }
  }
  DeleteMarkedObjects (orig_entityID, 0, NULL);
  return rval;
}


static void MakeAlignmentFileIdsGenbank (TAlignmentFilePtr afp)
{
  int i;
  CharPtr fmt = "gb|%s";
  CharPtr tmp;


  if (afp == NULL) {
    return;
  }
  for (i = 0; i < afp->num_sequences; i++) {
    if (StringNICmp(afp->ids[i], "gb|", 3) != 0) {
      tmp = (CharPtr) MemNew (sizeof (Char) * (StringLen (afp->ids[i]) + StringLen (fmt)));
      sprintf (tmp, fmt, afp->ids[i]);
      afp->ids[i] = MemFree (afp->ids[i]);
      afp->ids[i] = tmp;
    }
  }
}


static void MakeFeatureTableFileIdsGenbank (SeqAnnotPtr sap)
{
  SeqFeatPtr   sfp;
  SeqIdPtr     sip_orig, sip_new;
  ObjectIdPtr  oip;
  TextSeqIdPtr tsip;
  Boolean      found = FALSE;
  CharPtr      cp;

  if (sap == NULL || sap->type != 1) {
    return;
  }
  for (sfp = sap->data; sfp != NULL; sfp = sfp->next) {
    sip_orig = SeqLocId (sfp->location);
    if (sip_orig != NULL && sip_orig->choice == 1) {
      oip = (ObjectIdPtr) sip_orig->data.ptrvalue;
      if (oip != NULL && oip->str != NULL) { 
        found = TRUE;
        break;
      }
    }
  }
  if (found) {
    tsip = TextSeqIdNew ();
    tsip->accession = StringSave (oip->str);
    if ((cp = StringChr(tsip->accession, '.')) != NULL) {
      tsip->version = atoi (cp + 1);
      *cp = 0;
    }
    sip_new = ValNodeNew (NULL);
    sip_new->choice = SEQID_GENBANK;
    sip_new->data.ptrvalue = tsip;

    sip_orig = SeqIdDup (sip_orig);
    for (sfp = sap->data; sfp != NULL; sfp = sfp->next) {
      ReplaceSeqIdWithSeqIdInFeat (sip_orig, sip_new, sfp);
    }
    sip_orig = SeqIdFree (sip_orig);
    sip_new = SeqIdFree (sip_new);
  }
}


static void RemoveAlignmentTitlesAndOrganisms (TAlignmentFilePtr afp)
{
  int i;

  if (afp == NULL || afp->deflines == NULL || afp->num_deflines == 0) {
    return;
  }
  for (i = 0; i < afp->num_deflines; i++) {
    free (afp->deflines[i]); 
  }
  free (afp->deflines);
  afp->deflines = NULL;
  afp->num_deflines = 0;
  for (i = 0; i < afp->num_organisms; i++) {
    free (afp->organisms [i]);
  }  
  free (afp->organisms);
  afp->organisms = NULL;
  afp->num_organisms = 0;
}


static Boolean EndsWithAsn (CharPtr filename)
{
  Int4 len;

  if (StringHasNoText (filename)) {
    return FALSE;
  }
  len = StringLen (filename);
  if (StringCmp (filename + len - 4, ".asn") == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static void RemoveSingleItemAlignments (SeqAnnotPtr sap, Pointer data)
{
  SeqAlignPtr salp;

  if (sap != NULL && sap->type == 2 && (salp = sap->data) != NULL && salp->dim == 1) {
    sap->idx.deleteme = TRUE;
  }
}


static SeqEntryPtr ProcessGroup (ReblobGroupPtr rg)
{
  TErrorInfoPtr     error_list = NULL;
  ReadBufferData    rbd;
  TAlignmentFilePtr afp;
  SeqEntryPtr       sep = NULL;
  TSequenceInfoPtr  sequence_info;
  Uint2             entityID;
  AsnIoPtr          aip;

  if (rg == NULL) {
    return NULL;
  }

  if (EndsWithAsn(rg->filename)) {
    /* not an alignment, no reblobbing, just do features and source table */
    aip = AsnIoOpen (rg->filename, "r");
    sep = SeqEntryAsnRead (aip, NULL);
    AsnIoClose (aip);
    if (sep == NULL) {
      Message (MSG_ERROR, "Unable to read Seq-entry from file %s\n", rg->filename);
    }
    SeqEntrySetScope (sep);
    entityID = ObjMgrGetEntityIDForChoice (sep);
    AssignIDsInEntity (entityID, 0, NULL);
    SeqMgrIndexFeatures (entityID, NULL);
    StripFeaturesFromSelectedSequences (rg, entityID, rg->filename);
    AddFeaturesToGroup (rg, sep, entityID, rg->filename);
    sep = GetTopSeqEntryForEntityID (entityID);
    ApplySrcTableToGroup (rg, sep, entityID, rg->filename);
  } else {
    rbd.fp = FileOpen (rg->filename, "r");
    if (rbd.fp == NULL) {
      Message (MSG_ERROR, "Unable to open alignment file %s", rg->filename);
      return NULL;
    }
    rbd.current_data = NULL;
    sequence_info = GetDefaultSequenceInfo();
    afp = ReadAlignmentFile ( AbstractReadFunction,
                              (Pointer) &rbd,
                              AbstractReportError,
                              (Pointer) &error_list,
                              sequence_info);
    FileClose (rbd.fp);
    ErrorInfoFree (error_list);
    if (afp == NULL) {
      Message (MSG_ERROR, "Unable to read alignment from %s", rg->filename);
    } else if (! DoSequenceLengthsMatch (afp)) {
      Message (MSG_ERROR, "Your alignment in %s is incorrectly formatted.  Sequence plus gaps should be the same length for all sequences", rg->filename);
    } else {
      MakeAlignmentFileIdsGenbank (afp);
      RemoveAlignmentTitlesAndOrganisms (afp);
      sep = MakeSequinDataFromAlignment (afp, Seq_mol_na);
      entityID = ObjMgrGetEntityIDForChoice (sep);
      /* now - copy extras from original, add in new feature annotation */
      CopyDataFromOriginals (rg, sep, entityID);
      AddFeaturesToGroup (rg, sep, entityID, NULL);
      sep = GetTopSeqEntryForEntityID (entityID);
      ApplySrcTableToGroup (rg, sep, entityID, NULL);
      /* add quality scores */
      ApplyQualityScoresToGroup (rg, sep, entityID, NULL);
      /* collect pubs on set level */
      MovePopPhyMutPubs(sep);
      VisitAnnotsInSep (sep, NULL, RemoveSingleItemAlignments);
      DeleteMarkedObjects (entityID, 0, NULL);
    }

    SequenceInfoFree (sequence_info);
    AlignmentFileFree (afp);
  }

  return sep;
}


static Boolean NotTitle (SeqDescPtr sdp, Pointer extradata)
{

  if (sdp == NULL 
      || sdp->choice == Seq_descr_title) {
    return FALSE;
  } else {
    return TRUE;
  }
}


static SeqFeatPtr GetCDS (BioseqPtr bsp)
{
  SeqFeatPtr sfp = NULL, tmp_sfp;
  SeqAnnotPtr sap;
  BioseqSetPtr bssp;

  if (bsp == NULL) {
    return NULL;
  }
  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 1) {
      tmp_sfp = sap->data;
      while (tmp_sfp != NULL) {
        if (tmp_sfp->data.choice == SEQFEAT_CDREGION) {
          if (sfp == NULL) {
            sfp = tmp_sfp;
          } else {
            /* found second coding region */
            return NULL;
          }
        }
        tmp_sfp = tmp_sfp->next;
      }
    }
  }
  if (bsp->idx.parenttype == OBJ_BIOSEQSET && (bssp = (BioseqSetPtr) bsp->idx.parentptr) != NULL) {
    for (sap = bssp->annot; sap != NULL; sap = sap->next) {
      if (sap->type == 1) {
        tmp_sfp = sap->data;
        while (tmp_sfp != NULL) {
          if (tmp_sfp->data.choice == SEQFEAT_CDREGION) {
            if (sfp == NULL) {
              sfp = tmp_sfp;
            } else {
              /* found second coding region */
              return NULL;
            }
          }
          tmp_sfp = tmp_sfp->next;
        }
      }
    }
  }
  return sfp;
}


static Boolean CopyProteinId (SeqEntryPtr orig_sep, SeqEntryPtr new_sep, BioseqPtr prot_bsp)
{
  SeqEntryPtr orig_scope;
  SeqFeatPtr  orig_cds, new_cds;
  BioseqPtr   orig_nuc, new_nuc, new_prot;
  Boolean           rval = FALSE;
  SeqIdPtr          sip, sip_tmp;
  TextSeqIdPtr      tsip;

  orig_scope = SeqEntrySetScope (orig_sep);
  orig_cds = SeqMgrGetCDSgivenProduct (prot_bsp, NULL);
  if (orig_cds != NULL) {
    orig_nuc = BioseqFindFromSeqLoc (orig_cds->location);
    SeqEntrySetScope (new_sep);
    new_nuc = NULL;
    for (sip = orig_nuc->id; sip != NULL && new_nuc == NULL; sip = sip->next) {
      new_nuc = BioseqFind (sip);
      if (new_nuc == NULL && sip->choice == SEQID_GENBANK) {
        sip_tmp = SeqIdDup (sip);
        tsip = (TextSeqIdPtr) sip_tmp->data.ptrvalue;
        tsip->version = 0;
        new_nuc = BioseqFind (sip_tmp);
        sip_tmp = SeqIdFree (sip_tmp);
      }
    }
    if (new_nuc != NULL) {
      new_cds = GetCDS (new_nuc);
      if (new_cds != NULL) {
        new_prot = BioseqFindFromSeqLoc (new_cds->product);
        if (new_prot != NULL) {
          sip = SeqIdFindBest (prot_bsp->id, SEQID_GENBANK);
          if (sip != NULL) {
            sip = SeqIdDup (sip);
            BioseqReplaceID (new_prot, sip);
            sip = SeqIdFree (sip);
            rval = TRUE;
          }
        }
      }
    }
  }    

  SeqEntrySetScope (orig_scope);
  return rval;
}


static Boolean ProcessSingletons (ReblobGroupPtr rg, SeqEntryPtr top_sep)
{
  ValNodePtr     vnp, vnp_r;
  ReblobSeqPtr   rs;
  ValNodeBlock   file_list;
  Boolean        rval = TRUE;
  SeqEntryPtr    orig_sep, seq_sep;
  AsnIoPtr       aip;
  BioseqPtr      orig_bsp;
  CharPtr        filename;
  Uint2          orig_entityID;

  if (rg == NULL) {
    return FALSE;
  }

  InitValNodeBlock (&file_list, NULL);

  SeqEntrySetScope (top_sep);
 
  for (vnp = rg->seq_list; vnp != NULL; vnp = vnp->next) {
    rs = vnp->data.ptrvalue;
    if (rs != NULL && !StringHasNoText (rs->orig_file)) {
      ValNodeAddPointerToEnd (&file_list, 0, rs->orig_file);
      rs->bsp = BioseqFind (rs->sip);
    }
  }

  file_list.head = ValNodeSort (file_list.head, SortVnpByString);
  ValNodeUnique (&file_list.head, SortVnpByString, ValNodeFree);

  for (vnp = file_list.head; vnp != NULL; vnp = vnp->next) {
    filename = vnp->data.ptrvalue;
    aip = AsnIoOpen (filename, "r");
    orig_sep = SeqEntryAsnRead (aip, NULL);
    AsnIoClose (aip);
    if (orig_sep == NULL) {
      Message (MSG_ERROR, "Unable to read from file %s", filename);
      continue;
    }
    SeqEntrySetScope (orig_sep);
    orig_entityID = ObjMgrGetEntityIDForChoice (orig_sep);
    AssignIDsInEntity (orig_entityID, 0, NULL);
    SeqMgrIndexFeatures (orig_entityID, NULL);
    /* push descriptors so they will be attached to the singleton sequences */
    if (IS_Bioseq_set(orig_sep)) {
      PropagateSomeDescriptors(orig_sep, NotTitle, NULL);
    }

    rval &= StripFeaturesFromSelectedSequences (rg, orig_entityID, filename);

    /* add features */
    rval &= AddFeaturesToGroup (rg, orig_sep, orig_entityID, filename);

    /* add source table information */
    rval &= ApplySrcTableToGroup (rg, orig_sep, orig_entityID, filename);

    /* add quality scores */
    rval &= ApplyQualityScoresToGroup (rg, orig_sep, orig_entityID, filename);

    /* now make copies and add to top_sep */
    for (vnp_r = rg->seq_list; vnp_r != NULL; vnp_r = vnp_r->next) {
      rs = (ReblobSeqPtr) vnp_r->data.ptrvalue;
      if (StringCmp (rs->orig_file, filename) == 0) {
        /* find original Bioseq */
        orig_bsp = BioseqFind (rs->sip);
        /* don't add if protein, will have been added with nucleotide already */
        if (orig_bsp != NULL) {
          if (ISA_aa (orig_bsp->mol)) {
            /* if there's a replacement protein for this one, copy the accession */
            CopyProteinId (orig_sep, top_sep, orig_bsp);
          } else {
            seq_sep = GetBestTopParentForData (orig_entityID, orig_bsp);
            if (seq_sep != NULL) {
              seq_sep = AsnIoMemCopy (seq_sep, (AsnReadFunc) SeqEntryAsnRead, (AsnWriteFunc) SeqEntryAsnWrite);
              PromoteAllToWorstID (seq_sep);
              AddSeqEntryToSeqEntry (top_sep, seq_sep, TRUE);
            }
          }
        }
      }
    }
    
    SeqEntrySetScope (NULL);
    SeqEntryFree (orig_sep);
  }
  return rval;
}


static SeqEntryPtr ProcessGroupList (ValNodePtr group_list)
{
  ValNodePtr vnp;
  SeqEntryPtr top_sep, sep;
  BioseqSetPtr bssp;
  ReblobGroupPtr rg;
 
  bssp = BioseqSetNew ();
  bssp->_class = BioseqseqSet_class_genbank;
  top_sep = SeqEntryNew();
  top_sep->choice = 2;
  top_sep->data.ptrvalue = bssp;

  for (vnp = group_list; vnp != NULL; vnp = vnp->next) {
    rg = (ReblobGroupPtr) vnp->data.ptrvalue;
    if (StringHasNoText (rg->filename)) {
      /* process as singletons */
      ProcessSingletons (rg, top_sep);
    } else {
      sep = ProcessGroup (vnp->data.ptrvalue);
      if (sep != NULL) {
        AddSeqEntryToSeqEntry (top_sep, sep, TRUE);
      }
    }
  }
  
  return top_sep;
}


static Int4 ReadMapFile (CharPtr filename, CharPtr result_dir) 
{
  FILE *fp;
  ReadBufferData rbd;
  CharPtr        line;
  Boolean        collecting = FALSE;
  CharPtr        grp_file = NULL;
  ValNodeBlock   group_list;
  ValNodeBlock   seq_list;
  ValNodePtr     tokens;
  SeqEntryPtr    sep;
  AsnIoPtr       aip;
  Char           output_file[PATH_MAX];
  Int4           cluster_number = 1;

  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    Message (MSG_FATAL, "Unable to open map file %s", filename);
    return 1;
  }

  rbd.fp = fp;
  rbd.current_data = NULL;

  line = AbstractReadFunction (&rbd);
  while (line != NULL) {
    if (StringNICmp (line, "START", 5) == 0) {
      /* start collecting group */
      collecting = TRUE;
      InitValNodeBlock (&group_list, NULL);
      InitValNodeBlock (&seq_list, NULL);
      grp_file = MemFree (grp_file);
    } else if (StringNICmp (line, "STOP", 4) == 0) {
      /* finish last group */
      ValNodeAddPointerToEnd (&group_list, 0, ReblobGroupNew (grp_file, seq_list.head));

      /* process group list */
      sep = ProcessGroupList (group_list.head);
      if (sep == NULL) {
        Message (MSG_ERROR, "Unable to process cluster %d", cluster_number);
      } else {
        sprintf (output_file, "%s/cluster_%d.sqn", result_dir, cluster_number);
        aip = AsnIoOpen (output_file, "w");
        SeqEntryAsnWrite (sep, aip, NULL);
        AsnIoClose (aip);
        sep = SeqEntryFree (sep);
      }

      /* cleanup */
      cluster_number++;
      group_list.head = ReblobGroupListFree (group_list.head);
      collecting = FALSE;
      grp_file = MemFree (grp_file);
    } else if (collecting) {
      tokens = ReadOneColumnList (line);
      if (tokens != NULL) {
        if (grp_file == NULL) {
          grp_file = StringSave (tokens->data.ptrvalue);
        } else if (StringCmp (grp_file, tokens->data.ptrvalue) != 0) {
          ValNodeAddPointerToEnd (&group_list, 0, ReblobGroupNew (grp_file, seq_list.head));
          InitValNodeBlock (&seq_list, NULL);
          grp_file = MemFree (grp_file);
          grp_file = StringSave (tokens->data.ptrvalue);
        }  
        if (tokens->next != NULL) {
          ValNodeAddPointerToEnd (&seq_list, 0, ReblobSeqNew (tokens->next));
        }
      }
      tokens = ValNodeFreeData (tokens);
    }    
    line = AbstractReadFunction (&rbd);
  }
 
  FileClose(fp);
  return 0;
}


/* Args structure contains command-line arguments */

typedef enum {
  m_argMapFile = 0,
  r_argOutputPath,
} Arguments;


Args myargs [] = {
  {"Input Map File", NULL, NULL, NULL,
    TRUE, 'm', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char            app [64];
  CharPtr         base;
  CharPtr         results;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrSetMessageLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  /* finish resolving internal connections in ASN.1 parse tables */

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }


  /* process command line arguments */

  sprintf (app, "reblobber %s", REBLOBBER_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  SetAppProperty ("NcbiTbl2Asn", (void *) 1024);

  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  base = (CharPtr) myargs [m_argMapFile].strvalue;
  if (StringHasNoText (results) || StringHasNoText (base)) {
    Message (MSG_FATAL, "Must supply map file and results directory.");
    return 1;
  }
   
  return ReadMapFile (base, results);

}

