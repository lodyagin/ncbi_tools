/*   cleanasn.c
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
* File Name:  cleanasn.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   10/19/99
*
* $Revision: 6.395 $
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
#include <objfdef.h>
#include <objsub.h>
#include <sequtil.h>
#include <gather.h>
#include <sqnutils.h>
#include <explore.h>
#include <tofasta.h>
#include <toasn3.h>
#include <toporg.h>
#include <subutil.h>
#include <asn2gnbk.h>
#include <pmfapi.h>
#include <tax3api.h>
#include <asn2gnbi.h>
#include <ent2api.h>
#include <gbftdef.h>
#include <valid.h>
#include <utilpub.h>
#include <actutils.h>
#include <algo/blast/api/twoseq_api.h>
#ifdef INTERNAL_NCBI_CLEANASN
#include <accpubseq.h>
#endif
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>
#include <findrepl.h>
#include <connect/ncbi_gnutls.h>

#define CLEANASN_APP_VER "17.9"

CharPtr CLEANASN_APPLICATION = CLEANASN_APP_VER;

typedef struct sums {
  Int4          nucs;
  Int4          prts;
  Int4          recs;
} SumData, PNTR SumDataPtr;

typedef struct dbsums {
  SumData      genbank;
  SumData      embl;
  SumData      ddbj;
  SumData      refseq;
  SumData      other;
} DbSumData, PNTR DbSumPtr;

typedef struct counts {
  Int4          auth;
  Int4          bsec;
  Int4          clnr;
  Int4          gbbk;
  Int4          modr;
  Int4          move;
  Int4          norm;
  Int4          nucs;
  Int4          okay;
  Int4          othr;
  Int4          pack;
  Int4          prot;
  Int4          prts;
  Int4          publ;
  Int4          recs;
  Int4          sloc;
  Int4          sort;
  Int4          ssec;
  Int4          sync;
  Int4          titl;
} CountData, PNTR CountDataPtr;

typedef struct cleanflags {
  CharPtr       filename;
  Char          buf [64];
  Int4          gi;
  Int2          year;
  Boolean       stripSerial;
  Boolean       isRefSeq;
  Boolean       isEmblDdbj;
  Boolean       isSmallGenome;
  Boolean       batch;
  Boolean       binary;
  Boolean       compressed;
  Int2          type;
  CharPtr       results;
  CharPtr       outfile;
  CharPtr       firstfile;
  CharPtr       lastfile;
  Boolean       foundfirst;
  Boolean       foundlast;
  CharPtr       sourcedb;
  CharPtr       branch;
  CharPtr       report;
  CharPtr       selective;
  ModType       ffmode;
  CharPtr       ffdiff;
  CharPtr       asn2flat;
  CharPtr       asnval;
  CharPtr       clean;
  CharPtr       modernize;
  CharPtr       link;
  CharPtr       feat;
  CharPtr       desc;
  CharPtr       conv;
  CharPtr       zapuser;
  CharPtr       mods;
  CharPtr       prot;
  CharPtr       valid;
  ValNodePtr    action_list;
  Boolean       taxon;
  CharPtr       pub;
  Int4          unpubcount;
  ValNodePtr    unpublist;
  ValNodePtr    lastunpub;
  CountData     rawcounts;
  CountData     cumcounts;
  DbSumData     dbsums;
  AsnModulePtr  amp;
  AsnTypePtr    atp_bss;
  AsnTypePtr    atp_bsss;
  AsnTypePtr    atp_se;
  AsnTypePtr    atp_bsc;
  AsnTypePtr    bssp_atp;
  AsnTypePtr    atp_sbp;
  AsnTypePtr    atp_ssp;
  BioseqSet     bss;
  FILE          *logfp;
} CleanFlagData, PNTR CleanFlagPtr;

#ifdef INTERNAL_NCBI_CLEANASN
static CharPtr smartfetchproc = "SmartBioseqFetch";

static CharPtr smartfetchcmd = NULL;

extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", smartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}

static Int2 LIBCALLBACK SmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean SmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, smartfetchproc, smartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  SmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr tpasmartfetchproc = "TPASmartBioseqFetch";

static CharPtr tpasmartfetchcmd = NULL;

extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}

static Int2 LIBCALLBACK TPASmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_TPG) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean TPASmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, tpasmartfetchproc, tpasmartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  TPASmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr hupfetchproc = "HUPBioseqFetch";

static CharPtr hupfetchcmd = NULL;

extern Pointer ReadFromHUP (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromHUP (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (hupfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "HUP", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      hupfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (hupfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", hupfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", hupfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}

static Int2 LIBCALLBACK HUPBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (hupfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "HUP", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      hupfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (hupfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", hupfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", hupfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean HUPFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, hupfetchproc, hupfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  HUPBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}
#endif

//LCOV_EXCL_START
static void RemoveFeatUser (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL) return;
  if (sfp->ext != NULL) {
    sfp->ext = UserObjectFree (sfp->ext);
  }
}

static void RemoveFeatDbxref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  DbtagPtr    dbt;
  ValNodePtr  next, vnp;

  if (sfp == NULL) return;
  for (vnp = sfp->dbxref; vnp != NULL; vnp = next) {
    next = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    DbtagFree (dbt);
    MemFree (vnp);
  }
  sfp->dbxref = NULL;
}

static void RemoveFeatEvInf (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  GBQualPtr       gbq;
  GBQualPtr       nextqual;
  GBQualPtr PNTR  prevqual;
  Boolean         unlink;

  if (sfp == NULL) return;

  gbq = sfp->qual;
  prevqual = (GBQualPtr PNTR) &(sfp->qual);
  while (gbq != NULL) {
    nextqual = gbq->next;
    unlink = FALSE;
    if (StringICmp (gbq->qual, "experiment") == 0) {
      unlink = TRUE;
    } else if (StringICmp (gbq->qual, "inference") == 0) {
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

static Boolean GeneSpansOrigin (SeqMgrFeatContextPtr context, Int4 bsplength)

{
  Int4Ptr  ivals;

  if (context == NULL || bsplength < 1) return FALSE;
  ivals = context->ivals;
  if (ivals == NULL || context->numivals != 2) return FALSE;
  if (context->strand == Seq_strand_minus) {
    if (ivals [1] == 0 && ivals [2] == bsplength - 1) return TRUE;
  } else {
    if (ivals [2] == 0 && ivals [1] == bsplength - 1) return TRUE;
  }
  return FALSE;
}

static Boolean LIBCALLBACK FuseMultiIntervalGenes (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  BioseqPtr  bsp;
  SeqLocPtr  gslp;
  Boolean    noLeft;
  Boolean    noRight;

  if (sfp == NULL || sfp->location == NULL || context == NULL) return TRUE;

  if (sfp->data.choice != SEQFEAT_GENE) return TRUE;
  if (context->numivals < 2) return TRUE;
  bsp = context->bsp;
  if (bsp == NULL) return TRUE;

  if (sfp->excpt) {
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) return TRUE;
  }

  if (bsp->topology == 2) {
    if (context->numivals == 2 && GeneSpansOrigin (context, bsp->length)) return TRUE;
  }

  if (SeqLocId (sfp->location) == NULL) return TRUE;
  if (SeqLocMixedStrands (bsp, sfp->location)) return TRUE;

  gslp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
  if (gslp != NULL) {
    CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = gslp;
    FreeAllFuzz (sfp->location);
    SetSeqLocPartial (sfp->location, noLeft, noRight);
  }
  if (SeqLocPartialCheck (sfp->location) != SLP_COMPLETE) {
    sfp->partial = TRUE;
  }

  return TRUE;
}

static Boolean HasAdjacentIntervals (SeqLocPtr slp)

{
  Boolean    adjacent = FALSE;
  SeqLocPtr  tmp = NULL;
  SeqIntPtr  sip1 = NULL, sip2 = NULL, prevsip = NULL;
  Uint1      strand2 = 0;

  if (slp == NULL) return FALSE;

  while ((tmp = SeqLocFindNext (slp, tmp)) != NULL) {
    if (tmp->choice != SEQLOC_NULL) {
      if (tmp->choice != SEQLOC_INT) return FALSE;
      sip1 = prevsip;
      sip2 = (SeqIntPtr) (tmp->data.ptrvalue);
      strand2 = sip2->strand;

      if (sip1 != NULL && SeqIntCheck (sip2)) {
        if (SeqIdForSameBioseq (sip1->id, sip2->id)) {
          if (strand2 == Seq_strand_minus) {
            if (sip2->to + 1 == sip1->from) {
              adjacent = TRUE;
            }
          } else {
            if (sip1->to + 1 == sip2->from) {
              adjacent = TRUE;
            }
          }
        }
      }
      prevsip = sip2;
    }
  }

  return adjacent;
}

static Boolean LIBCALLBACK FuseAdjacentIntervalImpFeats (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  BioseqPtr  bsp;
  SeqLocPtr  gslp;
  Boolean    noLeft;
  Boolean    noRight;

  if (sfp == NULL || sfp->location == NULL || context == NULL) return TRUE;

  if (sfp->data.choice != SEQFEAT_IMP) return TRUE;
  if (context->numivals < 2) return TRUE;
  bsp = context->bsp;
  if (bsp == NULL) return TRUE;

  if (sfp->excpt) {
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) return TRUE;
  }

  if (bsp->topology == 2) {
    if (context->numivals == 2 && GeneSpansOrigin (context, bsp->length)) return TRUE;
  }

  if (SeqLocId (sfp->location) == NULL) return TRUE;
  if (SeqLocMixedStrands (bsp, sfp->location)) return TRUE;
  if (! HasAdjacentIntervals (sfp->location)) return TRUE;

  gslp = SeqLocMergeEx (bsp, sfp->location, NULL, FALSE, TRUE, TRUE, FALSE);
  if (gslp != NULL) {
    CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
    sfp->location = SeqLocFree (sfp->location);
    sfp->location = gslp;
    FreeAllFuzz (sfp->location);
    SetSeqLocPartial (sfp->location, noLeft, noRight);
  }
  if (SeqLocPartialCheck (sfp->location) != SLP_COMPLETE) {
    sfp->partial = TRUE;
  }

  return TRUE;
}
//LCOV_EXCL_STOP

static void IsSmallGenomeSet (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  BoolPtr  bp;

  if (bssp == NULL || userdata == NULL) return;
  if (bssp->_class != BioseqseqSet_class_small_genome_set) return;
  bp = (BoolPtr) userdata;
  *bp = TRUE;
}

//LCOV_EXCL_START
static void MultiIntervalGenesCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  if (bsp->repr != Seq_repr_delta) return;

  SeqMgrExploreFeatures (bsp, NULL, FuseMultiIntervalGenes, NULL, NULL, NULL);
}

static void AdjacentIntervalImpFeatCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  SeqMgrExploreFeatures (bsp, NULL, FuseAdjacentIntervalImpFeats, NULL, NULL, NULL);
}

static Boolean StringsDoMatch (
  CharPtr str1,
  CharPtr str2,
  Boolean case_sensitive
)

{
  if (case_sensitive) {
    if (StringCmp (str1, str2) == 0) {
      return TRUE;
    }
  } else {
    if (StringICmp (str1, str2) == 0) {
      return TRUE;
    }
  }

  return FALSE;
}

static Boolean LocationPartialsDiffer (
  SeqLocPtr slp1,
  SeqLocPtr slp2
)

{
  Boolean  partial5_1, partial3_1, partial1;
  Boolean  partial5_2, partial3_2, partial2;

  if (slp1 == NULL || slp2 == NULL) return TRUE;

  partial1 = CheckSeqLocForPartial (slp1, &partial5_1, &partial3_1);
  partial2 = CheckSeqLocForPartial (slp2, &partial5_2, &partial3_2);

  if (partial1 != partial2) return TRUE;
  if (partial5_1 != partial5_2) return TRUE;
  if (partial3_1 != partial3_2) return TRUE;

  return FALSE;
}

static Boolean LocationsDoMatch (
  SeqLocPtr slp1,
  SeqLocPtr slp2,
  Boolean allow_different_sequences,
  Boolean strict_partial
)

{
  SeqLocPtr  slp_tmp1, slp_tmp2;

  if (slp1 == NULL || slp2 == NULL) return FALSE;

  if (strict_partial && LocationPartialsDiffer (slp1, slp2)) return FALSE;

  if (allow_different_sequences) {
    for (slp_tmp1 = SeqLocFindNext (slp1, NULL), slp_tmp2 = SeqLocFindNext (slp2, NULL);
         slp_tmp1 != NULL && slp_tmp2 != NULL;
         slp_tmp1 = SeqLocFindNext (slp1, slp_tmp1), slp_tmp2 = SeqLocFindNext (slp2, slp_tmp2)) {
      if (SeqLocStart (slp_tmp1) != SeqLocStart (slp_tmp2)) return FALSE;
      if (SeqLocStop (slp_tmp1) != SeqLocStop (slp_tmp2)) return FALSE;
      if (strict_partial && LocationPartialsDiffer (slp_tmp1, slp_tmp2)) return FALSE;
    }

    return TRUE;
  }
  
  if (SeqLocCompare (slp1, slp2) != SLC_A_EQ_B) return FALSE;

  return TRUE;
}

static CharPtr GetFeatureProductName (
  SeqFeatPtr sfp
)

{
  Uint1           aacode = 0;
  GBQualPtr       gbp;
  GeneRefPtr      grp;
  BioseqPtr       pbsp;
  ProtRefPtr      prp = NULL;
  SeqFeatPtr      psfp;
  RNAGenPtr       rgp;
  RnaRefPtr       rrp;
  SeqCodeTablePtr sctp;
  SeqFeatXrefPtr  sfxp;
  SeqMapTablePtr  smtp;
  CharPtr         str;
  tRNAPtr         trp;
  ValNodePtr      vnp;

  if (sfp == NULL) return NULL;

  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
      if (sfp->comment != NULL) return sfp->comment;
      return NULL;
    case SEQFEAT_BOND:
      str = AsnEnumStr ("SeqFeatData.bond", (Int2) sfp->data.value.intvalue);
      if (str != NULL) return str;
      return NULL;
    case SEQFEAT_SITE:
      str = AsnEnumStr ("SeqFeatData.site", (Int2) sfp->data.value.intvalue);
      if (str != NULL) return str;
      return NULL;
    case SEQFEAT_PSEC_STR:
      return NULL;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return NULL;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp->locus != NULL) return grp->locus;
      if (grp->desc != NULL) return grp->desc;
      if (grp->locus_tag != NULL) return grp->locus_tag;
      if (grp->syn != NULL) {
        vnp = grp->syn;
        str = (CharPtr) vnp->data.ptrvalue;
        if (str != NULL) return str;
      }
      break;
    case SEQFEAT_CDREGION :
      for (sfxp = sfp->xref; sfxp != NULL; sfxp = sfxp->next) {
        if (sfxp->data.choice != SEQFEAT_PROT) continue;
        prp = (ProtRefPtr) sfxp->data.value.ptrvalue;
        if (prp != NULL) break;
      }
      if (prp == NULL && sfp->product != NULL) {
        pbsp = BioseqFindFromSeqLoc (sfp->product);
        if (pbsp != NULL) {
          psfp = SeqMgrGetBestProteinFeature (pbsp, NULL);
          if (psfp != NULL) {
            prp = (ProtRefPtr) psfp->data.value.ptrvalue;
          }
        }
      }
      if (prp == NULL) return NULL;
      if (prp->name != NULL) {
        vnp = prp->name;
        str = (CharPtr) vnp->data.ptrvalue;
        if (str != NULL) return str;
      }
      if (prp->desc != NULL) return prp->desc;
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      if (prp->name != NULL) {
        vnp = prp->name;
        str = (CharPtr) vnp->data.ptrvalue;
        if (str != NULL) return str;
      }
      if (prp->desc != NULL) return prp->desc;
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      switch (rrp->ext.choice) {
        case 1 :
          str = (CharPtr) rrp->ext.value.ptrvalue;
          if (str == NULL) return NULL;
          if (StringCmp (str, "ncRNA") == 0 ||
                         StringCmp (str, "tmRNA") == 0 ||
                         StringCmp (str, "misc_RNA") == 0) {
            for (gbp = sfp->qual; gbp != NULL; gbp = gbp->next) {
              if (StringICmp ("product", gbp->qual) != 0) continue;
              if (gbp->val != NULL) return gbp->val;
            }
          }
          if (str != NULL) return str;
          break;
        case 2 :
          trp = (tRNAPtr) rrp->ext.value.ptrvalue;
          if (trp == NULL) return NULL;
          switch (trp->aatype) {
            case 1:
              aacode = Seq_code_iupacaa;
              break;
            case 2:
              aacode = Seq_code_ncbieaa;
              break;
            case 3:
              aacode = Seq_code_ncbi8aa;
              break;
            case 4:
              aacode = Seq_code_ncbistdaa;
              break;
            default:
              aacode = 0;
              break;
          }
          if (aacode == 0) return NULL;
          sctp = SeqCodeTableFind (Seq_code_iupacaa3);
          if (sctp == NULL) return NULL;
          if (aacode != Seq_code_ncbistdaa) {
            smtp = SeqMapTableFind (Seq_code_ncbistdaa, aacode);
            if (smtp == NULL) return NULL;
            aacode = SeqMapTableConvert (smtp, trp->aa);
          } else {
            aacode = trp->aa;
          }
          if (aacode == 255) {
            if (trp->aatype == Seq_code_iupacaa || trp->aatype == Seq_code_ncbieaa) {
              if (trp->aa == 74) {
                aacode = 27; /* Xle */
              } else if (trp->aa == 79) {
                aacode = 26; /* Pyl */
              }
            }
          }
          if (aacode == 255) return NULL;
          str = sctp->symbols [aacode - sctp->start_at];
          if (str != NULL) return str;
          break;
        case 3 :
          rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
          if (rgp == NULL) return NULL;
          if (rgp->product != NULL) return rgp->product;
          if (rgp->_class != NULL) return rgp->_class;
          break;
        default :
          break;
      }
      break;
    case SEQFEAT_REGION :
      str = (CharPtr) sfp->data.value.ptrvalue;
      if (str == NULL) return NULL;
      if (sfp->comment != NULL) {
        if (StringICmp (str, "Domain") == 0) return sfp->comment;
        if (StringICmp (str, "Variant") == 0) return sfp->comment;
      }
      return str;
    default :
      break;
  }

  return NULL;
}

static Boolean ProductNamesMatch (
  SeqFeatPtr sfp1,
  SeqFeatPtr sfp2,
  Boolean case_sensitive
)

{
  CharPtr  str1, str2;

  if (sfp1 == NULL || sfp2 == NULL) return FALSE;
  if (sfp1->data.choice != sfp2->data.choice) return FALSE;
  if (sfp1->idx.subtype != sfp2->idx.subtype) return FALSE;

  str1 = GetFeatureProductName (sfp1);
  str2 = GetFeatureProductName (sfp2);
  if (str1 != NULL || str2 != NULL) {
    if (str1 == NULL || str2 == NULL) return FALSE;
    if (! StringsDoMatch (str1, str2, case_sensitive)) {
      return FALSE;
    }
  }

  return TRUE;
}

static Boolean CdRegionsAreCompatible (
  CdRegionPtr crp1,
  CdRegionPtr crp2
)

{
  if (crp1 == NULL || crp2 == NULL) return FALSE;
  if (crp1->orf != crp2->orf) return FALSE;
  if (crp1->conflict != crp2->conflict) return FALSE;
  if (crp1->gaps != crp2->gaps) return FALSE;
  if (crp1->mismatch != crp2->mismatch) return FALSE;
  if (crp1->stops != crp2->stops) return FALSE;

  if (crp1->frame != crp2->frame) {
    if (crp1->frame > 1 || crp2->frame > 1) return FALSE;
  }

  if (crp1->genetic_code != NULL || crp2->genetic_code != NULL) {
    if (crp1->genetic_code == NULL || crp2->genetic_code == NULL) return FALSE;
    if (! AsnIoMemComp (crp1->genetic_code, crp2->genetic_code, (AsnWriteFunc) GeneticCodeAsnWrite)) return FALSE;
  }

  if (crp1->code_break != NULL || crp2->code_break != NULL) {
    if (crp1->code_break == NULL || crp2->code_break == NULL) return FALSE;
    if (! AsnIoMemComp (crp1->code_break, crp2->code_break, (AsnWriteFunc) CodeBreakAsnWrite)) return FALSE;
  }

  return TRUE;
}

static Boolean CdProductsAreCompatible (
  SeqLocPtr slp1,
  SeqLocPtr slp2,
  Boolean case_sensitive,
  Boolean strict_partial
)

{
  BioseqPtr          bsp1, bsp2;
  SeqMgrFeatContext  fcontext1, fcontext2;
  ProtRefPtr         prp1, prp2;
  SeqFeatPtr         sfp1, sfp2;

  if (slp1 == NULL || slp2 == NULL) return FALSE;

  bsp1 = BioseqFindFromSeqLoc (slp1);
  bsp2 = BioseqFindFromSeqLoc (slp2);
  if (bsp1 == NULL || bsp2 == NULL) return FALSE;

  sfp1 = SeqMgrGetNextFeature (bsp1, NULL, SEQFEAT_PROT, 0, &fcontext1);
  sfp2 = SeqMgrGetNextFeature (bsp2, NULL, SEQFEAT_PROT, 0, &fcontext2);

  while (sfp1 != NULL && sfp2 != NULL) {
    if (! LocationsDoMatch (sfp1->location, sfp2->location, TRUE, strict_partial)) return FALSE;
    if (! ProductNamesMatch (sfp1, sfp2, case_sensitive)) return FALSE;

    prp1 = (ProtRefPtr) sfp1->data.value.ptrvalue;
    prp2 = (ProtRefPtr) sfp2->data.value.ptrvalue;
    if (prp1 == NULL || prp2 == NULL) return FALSE;
    if (prp1->processed != prp2->processed) return FALSE;

    sfp1 = SeqMgrGetNextFeature (bsp1, sfp1, SEQFEAT_PROT, 0, &fcontext1);
    sfp2 = SeqMgrGetNextFeature (bsp2, sfp2, SEQFEAT_PROT, 0, &fcontext2);
  }

  if (sfp1 != NULL || sfp2 != NULL) return FALSE;

  return TRUE;
}

static Boolean IsConservativelyFuseable (
  SeqFeatPtr sfp1,
  SeqFeatPtr sfp2,
  Boolean case_sensitive,
  Boolean strict_partial
)

{
  CdRegionPtr  crp1, crp2;

  if (sfp1 == NULL || sfp2 == NULL) return FALSE;
  if (sfp1->data.choice != sfp2->data.choice) return FALSE;
  if (sfp1->idx.subtype != sfp2->idx.subtype) return FALSE;

  if (sfp1->pseudo != sfp2->pseudo) return FALSE;
  if (sfp1->exp_ev != sfp2->exp_ev) return FALSE;

  if (sfp1->excpt != sfp2->excpt) return FALSE;
  if (! StringsDoMatch (sfp1->except_text, sfp2->except_text, case_sensitive)) return FALSE;

  if (strict_partial && sfp1->partial != sfp2->partial) return FALSE;
  if (! LocationsDoMatch (sfp1->location, sfp2->location, FALSE, strict_partial)) return FALSE;

  if (! ProductNamesMatch (sfp1, sfp2, case_sensitive)) {
    return FALSE;
  }

  if (sfp1->data.choice == SEQFEAT_CDREGION) {
    crp1 = (CdRegionPtr) sfp1->data.value.ptrvalue;
    crp2 = (CdRegionPtr) sfp2->data.value.ptrvalue;
    if (! CdRegionsAreCompatible (crp1, crp2)) return FALSE;
    if (! CdProductsAreCompatible (sfp1->product, sfp2->product, case_sensitive, strict_partial)) return FALSE;
  }

  return TRUE;
}

// Used in FuseCommonFeatureFields (unused)
static void CombineTexts (
  CharPtr PNTR txtptr,
  CharPtr PNTR oldtxtptr
)

{
  size_t   len;
  CharPtr  str;

  if (txtptr == NULL || oldtxtptr == NULL) return;

  if (*txtptr == NULL) {

    *txtptr = *oldtxtptr;
    *oldtxtptr = NULL;

  } else if (*oldtxtptr != NULL && StringICmp (*txtptr, *oldtxtptr) != 0) {

    len = StringLen (*txtptr) + StringLen (*oldtxtptr) + 5;
    str = MemNew (sizeof (Char) * len);
    StringCpy (str, *txtptr);
    StringCat (str, "; ");
    StringCat (str, *oldtxtptr);
    *txtptr = MemFree (*txtptr);
    *txtptr = str;
  }
}

// used in FuseFeatureQualifiers (unused)
static void FuseCommonFeatureFields (
  SeqFeatPtr sfp,
  SeqFeatPtr oldsfp
)

{
  GBQualPtr       lastgbq;
  SeqFeatXrefPtr  lastxref;
  Boolean         oldpartial5;
  Boolean         oldpartial3;
  Boolean         partial5;
  Boolean         partial3;

  if (sfp == NULL || oldsfp == NULL) return;

  CombineTexts (&(sfp->comment), &(oldsfp->comment));
  CombineTexts (&(sfp->title), &(oldsfp->title));
  CombineTexts (&(sfp->except_text), &(oldsfp->except_text));

  if (sfp->qual == NULL) {
    sfp->qual = oldsfp->qual;
    oldsfp->qual = NULL;
  } else if (oldsfp->qual != NULL) {
    for (lastgbq = sfp->qual; lastgbq->next != NULL; lastgbq = lastgbq->next) continue;
    lastgbq->next = oldsfp->qual;
    oldsfp->qual = NULL;
  }

  ValNodeLink (&(sfp->dbxref), oldsfp->dbxref);
  oldsfp->dbxref = NULL;

  ValNodeLink (&(sfp->cit), oldsfp->cit);
  oldsfp->cit = NULL;

  if (sfp->xref == NULL) {
    sfp->xref = oldsfp->xref;
    oldsfp->xref = NULL;
  } else if (oldsfp->xref != NULL) {
    for (lastxref = sfp->xref; lastxref->next != NULL; lastxref = lastxref->next) continue;
    lastxref->next = oldsfp->xref;
    oldsfp->xref = NULL;
  }

  if (sfp->ext == NULL) {
    sfp->ext = oldsfp->ext;
    oldsfp->ext = NULL;
  } else if (oldsfp->ext != NULL) {
    sfp->ext = CombineUserObjects (sfp->ext, oldsfp->ext);
    oldsfp->ext = NULL;
  }

  CheckSeqLocForPartial (oldsfp->location, &oldpartial5, &oldpartial3);
  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);

  SetSeqLocPartial (sfp->location, oldpartial5 | partial5, oldpartial3 | partial3);

  sfp->partial |= oldsfp->partial;
  sfp->excpt |= oldsfp->excpt;
  sfp->pseudo |= oldsfp->pseudo;
}

// used in FuseDuplicateFeats (unused)
static void FuseFeatureQualifiers (
  SeqFeatPtr sfp,
  SeqFeatPtr oldsfp
)

{
  GeneRefPtr  grp, oldgrp;
  BioseqPtr   prod, oldprod;
  SeqFeatPtr  prot, oldprot;
  ProtRefPtr  prp, oldprp;
  RnaRefPtr   rrp, oldrrp;
  SeqIdPtr    sip;

  if (sfp == NULL || oldsfp == NULL) return;

  /* merge common fields */

  FuseCommonFeatureFields (sfp, oldsfp);

  /* now deal with type-specific data */

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      oldgrp = (GeneRefPtr) oldsfp->data.value.ptrvalue;
      if (grp == NULL || oldgrp == NULL) return;
      CombineTexts (&(grp->locus), &(oldgrp->locus));
      CombineTexts (&(grp->allele), &(oldgrp->allele));
      CombineTexts (&(grp->desc), &(oldgrp->desc));
      CombineTexts (&(grp->maploc), &(oldgrp->maploc));
      CombineTexts (&(grp->locus_tag), &(oldgrp->locus_tag));
      grp->pseudo |= oldgrp->pseudo;
      ValNodeLink (&(grp->db), oldgrp->db);
      oldgrp->db = NULL;
      ValNodeLink (&(grp->syn), oldgrp->syn);
      oldgrp->syn = NULL;
      break;
    case SEQFEAT_CDREGION :
      sip = SeqLocId (sfp->product);
      prod = BioseqFind (sip);
      sip = SeqLocId (oldsfp->product);
      oldprod = BioseqFind (sip);
      if (prod == NULL || oldprod == NULL) return;
      prot = SeqMgrGetBestProteinFeature (prod, NULL);
      oldprot = SeqMgrGetBestProteinFeature (oldprod, NULL);
      if (prot == NULL || oldprot == NULL) return;
      FuseCommonFeatureFields (prot, oldprot);
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
      oldprp = (ProtRefPtr) oldprot->data.value.ptrvalue;
      if (prp == NULL || oldprp == NULL) return;
      ValNodeLink (&(prp->name), oldprp->name);
      oldprp->name = NULL;
      ValNodeLink (&(prp->ec), oldprp->ec);
      oldprp->ec = NULL;
      ValNodeLink (&(prp->activity), oldprp->activity);
      oldprp->activity = NULL;
      ValNodeLink (&(prp->db), oldprp->db);
      oldprp->db = NULL;
      CombineTexts (&(prp->desc), &(oldprp->desc));
      /* mark protein product for deletion */
      /*
      oldprot->idx.deleteme = TRUE;
      */
      oldprod->idx.deleteme = TRUE;
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      oldrrp = (RnaRefPtr) oldsfp->data.value.ptrvalue;
      if (rrp == NULL || oldrrp == NULL) return;
      if (rrp->ext.choice == 1 && oldrrp->ext.choice == 1) {
        CombineTexts ((CharPtr PNTR) &(rrp->ext.value.ptrvalue), (CharPtr PNTR) &(oldrrp->ext.value.ptrvalue));
      }
      break;
    case SEQFEAT_REGION :
    case SEQFEAT_COMMENT :
      if (sfp->data.value.ptrvalue == NULL || oldsfp->data.value.ptrvalue == NULL) return;
      CombineTexts ((CharPtr PNTR) &(sfp->data.value.ptrvalue), (CharPtr PNTR) &(oldsfp->data.value.ptrvalue));
      break;
    default :
      break;
  }
}

// used for -F f (unused option)
static void FuseDuplicateFeats (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         sfp, lastsfp = NULL;

  if (bsp == NULL) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  while (sfp != NULL) {
    if (lastsfp != NULL && IsConservativelyFuseable (sfp, lastsfp, FALSE, FALSE)) {
      /* fuse qualifiers */
      FuseFeatureQualifiers (lastsfp, sfp);
      /* then mark for deletion */
      sfp->idx.deleteme = TRUE;
    } else {
      lastsfp = sfp;
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
  }
}
//LCOV_EXCL_STOP

static void SetBestFrameByLoc (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CdRegionPtr  crp;
  Uint1        new_frame = 0, i;
  ByteStorePtr bs;
  Int4         lens [3];
  Int4         max;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  
  crp = sfp->data.value.ptrvalue;
  if (crp == NULL) return;

  max = 0;
  for (i = 1; i <= 3; i++) {
    crp->frame = i;
    bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
    lens[i - 1] = BSLen (bs);
    BSFree (bs);
    if (lens[i - 1] > max) {
      max = lens[i - 1];
      new_frame = i;
    }
  }
  crp->frame = new_frame;
}

static void RetranslateCDS (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;

  RetranslateOneCDS (sfp, sfp->idx.entityID, TRUE, TRUE);
}

// Used in AdjustForStopCodon
static CharPtr GetProtString (SeqFeatPtr sfp)

{
  ByteStorePtr bs;
  CharPtr      prot = NULL;

  if (sfp == NULL || sfp->location == NULL) return NULL;

  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (bs != NULL) {
    prot = BSMerge (bs, NULL);
    bs = BSFree (bs);
  }

  return prot;
}


static SeqLocPtr TruncateLoc (SeqLocPtr head, Int4 len)

{
  SeqLocPtr slp = NULL;
  SeqIdPtr sip;
  Int4 from = 0, to = 0;
  Boolean changed = FALSE;
  Int4       loc_len = 0;
  Int4       cum_len = 0;
  SeqLocPtr  del_slp;
  ValNodePtr del_slp_list = NULL, vnp;
  Uint1      strand;

  if (head == NULL || len < 1) return head;

  slp = SeqLocFindNext (head, slp);
  while (slp != NULL) {
    if (cum_len >= len) {
      del_slp = AsnIoMemCopy (slp,
                              (AsnReadFunc) SeqLocAsnRead,
                              (AsnWriteFunc) SeqLocAsnWrite);
      ValNodeAddPointer (&del_slp_list, 0, del_slp);
    } else {
      sip = SeqLocId(slp);
      from = SeqLocStart(slp);
      to = SeqLocStop(slp);
      loc_len = (to - from + 1);
      cum_len += loc_len;
      if (cum_len > len) {
        strand = SeqLocStrand (slp);
        /* remove part of this location */
        if (strand == Seq_strand_minus) {
          to = cum_len - len + from - 1;
        } else {
          from = to - (cum_len - len) + 1;
        }
        del_slp = SeqLocIntNew (from, to, strand, sip);
        ValNodeAddPointer (&del_slp_list, 0, del_slp);
      }
    }
    slp = SeqLocFindNext (head, slp);
  }
  
  for (vnp = del_slp_list; vnp != NULL; vnp = vnp->next) {
    slp = (SeqLocPtr) vnp->data.ptrvalue;
    if (slp == NULL) continue;
    sip = SeqLocId(slp);
    from = SeqLocStart(slp);
    to = SeqLocStop(slp);
    head = SeqLocDelete(head, sip, from, to, FALSE, &changed);    
  }
  
  return head;
}

// Used in AdjustForStopCodon (active)
static Boolean LengthenLoc (SeqLocPtr head, Int4 len)

{
  SeqLocPtr slp = NULL, last_slp = NULL;
  SeqIntPtr sintp;
  SeqPntPtr spp;
  Int4        orig_len, len_diff;
  Boolean     rval = FALSE;
  BioseqPtr   bsp;

  if (head == NULL || len < 1) return FALSE;

  orig_len = SeqLocLen (head);
  len_diff = len - orig_len;
  if (len_diff <= 0) return FALSE;

  slp = SeqLocFindNext (head, slp);
  while (slp != NULL) {
    last_slp = slp;
    slp = SeqLocFindNext (head, slp);
  }

  if (last_slp != NULL && last_slp->choice == SEQLOC_INT) {
    sintp = (SeqIntPtr) last_slp->data.ptrvalue;
    if (sintp->strand == Seq_strand_minus) {
      if (sintp->from - len_diff >= 0) {
        sintp->from -= len_diff;
        rval = TRUE;
      }
    } else {
      bsp = BioseqFind (sintp->id);
      if (bsp != NULL && sintp->to + len_diff < bsp->length) {
        sintp->to += len_diff;
        rval = TRUE;
      }
    }      
  } else if (last_slp != NULL && last_slp->choice == SEQLOC_PNT) {
    spp = (SeqPntPtr) last_slp->data.ptrvalue;
    rval = TRUE;
    if (spp->strand == Seq_strand_minus) {
      if (spp->point - len_diff < 0) {
        rval = FALSE;
      }
    } else {
      bsp = BioseqFind (spp->id);
      if (bsp == NULL || spp->point + len_diff >= bsp->length) {
        rval = FALSE;
      }
    }
    if (rval) {
      sintp = SeqIntNew ();
      sintp->id = SeqIdDup (spp->id);
      sintp->strand = spp->strand;
      sintp->if_from = AsnIoMemCopy (spp->fuzz, (AsnReadFunc) IntFuzzAsnRead, (AsnWriteFunc) IntFuzzAsnWrite);
      sintp->if_from = AsnIoMemCopy (spp->fuzz, (AsnReadFunc) IntFuzzAsnRead, (AsnWriteFunc) IntFuzzAsnWrite);
      if (sintp->strand == Seq_strand_minus) {
        sintp->from = spp->point - len_diff;
        sintp->to = spp->point;
      } else {
        sintp->from = spp->point;
        sintp->to = spp->point + len_diff;
      }
      spp = SeqPntFree (spp);
      last_slp->data.ptrvalue = sintp;
      last_slp->choice = SEQLOC_INT;
    }
  }
  
  return rval;
}

static void CopyCDSproductToMrna (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  SeqFeatPtr  bestprot;
  SeqFeatPtr  mrna;
  CharPtr     mRNAname = NULL;
  ValNodePtr  name;
  ProtRefPtr  prp;
  RnaRefPtr   rrp;
  CharPtr     str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || sfp->location == NULL) return;

  mrna = GetmRNAforCDS (sfp);
  if (mrna == NULL) return;

  bestprot = FindBestProtein (sfp->idx.entityID, sfp->product);
  if (bestprot == NULL) return;
  prp = bestprot->data.value.ptrvalue;
  if (prp == NULL) return;

  name = prp->name;
  if (name != NULL && !StringHasNoText (name->data.ptrvalue)) {
    mRNAname = StringSave (name->data.ptrvalue);
  } else if (!StringHasNoText (prp->desc)) {
    mRNAname = StringSave (prp->desc);
  }
  if (mRNAname == NULL) return;

  rrp = (RnaRefPtr) mrna->data.value.ptrvalue;
  if (rrp == NULL) {
    rrp = RnaRefNew ();
    mrna->data.value.ptrvalue = (Pointer) rrp;
  }

  if (rrp != NULL) {
    rrp->type = 2;
    rrp->ext.choice = 1;
    str = rrp->ext.value.ptrvalue;
    if (StringDoesHaveText (str)) return;
    rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    rrp->ext.value.ptrvalue = mRNAname;
  }

}

static void AdjustForStopCodon (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  SeqFeatPtr   gene = NULL, mrna = NULL;
  Boolean      mrna_same = FALSE;
  SeqLocPtr    slp, new_slp = NULL;
  SeqLocPtr    orig_slp;
  BioseqPtr    nucBsp;
  Boolean      partial5, partial3;
  Uint1        strand;
  CharPtr      prot, cp;
  Int4         desired_cds_len, loc_len, max_len;
  CdRegionPtr  crp;
  Uint1        frame = 0;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || sfp->location == NULL) return;

  orig_slp = sfp->location;
  if (orig_slp == NULL) {
    return;
  }

  CheckSeqLocForPartial (orig_slp, &partial5, &partial3);

  nucBsp = BioseqFindFromSeqLoc (orig_slp);
  if (nucBsp == NULL) {
    return;
  }

  gene = GetGeneForFeature (sfp);
  mrna = GetmRNAforCDS (sfp);

  if (mrna != NULL && SeqLocCompareEx (sfp->location, mrna->location, TRUE) == SLC_A_EQ_B) {
    mrna_same = TRUE;
  }

  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp != NULL) {
    frame = crp->frame;
    if (frame > 0) {
      frame--;
    }
  }

  prot = GetProtString (sfp);
  if (StringHasNoText (prot)) {
    /*
    Message (MSG_ERROR, "Translation failed");
    */
    prot = MemFree (prot);
    return;
  } 

  cp = StringChr (prot, '*');

  if (cp == NULL) {
    /* lengthen coding region, look for stop codon */
    strand = SeqLocStrand (orig_slp);
    loc_len = SeqLocLen (orig_slp);
    if (strand == Seq_strand_minus) {
      max_len = loc_len + SeqLocStart (orig_slp);
    } else {
      max_len = loc_len + nucBsp->length - SeqLocStop (orig_slp) - 1;
    }
    
    if (LengthenLoc (orig_slp, max_len)) {
      prot = MemFree (prot);
      prot = GetProtString (sfp);
      cp = StringChr (prot, '*');
      if (cp == NULL) {
        /*
        if (ANS_YES == Message (MSG_YN, "No stop codon found - extend to full length of sequence (and make partial)?")) {
          SetSeqLocPartial (orig_slp, partial5, TRUE);
          RetranslateOneCDS (sfp, sfp->idx.entityID, TRUE, TRUE);
        }
        */
        prot = MemFree (prot);
        if (mrna_same) {
          if (mrna != NULL) {
            slp = SeqLocMerge (nucBsp, mrna->location, sfp->location, FALSE, TRUE, FALSE);
            if (slp != NULL) {
              mrna->location = SeqLocFree (mrna->location);
              mrna->location = slp;
            }
          }
          if (gene != NULL) {
            slp = SeqLocMerge (nucBsp, gene->location, sfp->location, TRUE, FALSE, FALSE);
            if (slp != NULL) {
              gene->location = SeqLocFree (gene->location);
              gene->location = slp;
            }
          }
        }
        return;
      }
    } else {
      /*
      Message (MSG_ERROR, "Unable to lengthen location!");
      */
      prot = MemFree (prot);
      return;
    }
  }

  desired_cds_len = ((cp - prot) + 1) * 3;
  if (frame > 0) {
    desired_cds_len += frame;
  }
    
  loc_len = SeqLocLen (orig_slp);
  if (desired_cds_len < loc_len) {
    /* truncate to correct length */
    new_slp = TruncateLoc (orig_slp, desired_cds_len);
    if (new_slp != NULL) {
      sfp->location = new_slp;
      orig_slp = new_slp;
    }
  }
  SetSeqLocPartial (orig_slp, partial5, FALSE);

  RetranslateOneCDS (sfp, sfp->idx.entityID, TRUE, TRUE);

  prot = MemFree (prot);

  if (mrna_same) {
    if (mrna != NULL) {
      slp = SeqLocMerge (nucBsp, mrna->location, sfp->location, FALSE, TRUE, FALSE);
      if (slp != NULL) {
        mrna->location = SeqLocFree (mrna->location);
        mrna->location = slp;
      }
    }
    if (gene != NULL) {
      slp = SeqLocMerge (nucBsp, gene->location, sfp->location, TRUE, FALSE, FALSE);
      if (slp != NULL) {
        gene->location = SeqLocFree (gene->location);
        gene->location = slp;
      }
    }
  }
}

static void AddMissingProtRef (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  GBQualPtr   gbq;
  CharPtr     newname;
  BioseqPtr   pbsp;
  ProtRefPtr  prp;
  SeqFeatPtr  psfp;
  CharPtr     str;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || sfp->product == NULL) return;
  newname = userdata;
  if (StringHasNoText (newname)) return;

  pbsp = BioseqFindFromSeqLoc (sfp->product);
  if (pbsp == NULL) return;
  psfp = SeqMgrGetBestProteinFeature (pbsp, NULL);
  if (psfp != NULL) return;

  psfp = CreateNewFeatureOnBioseq (pbsp, SEQFEAT_PROT, NULL);
  if (psfp == NULL) return;

  prp = ProtRefNew ();
  if (prp == NULL) return;
  psfp->data.value.ptrvalue = (Pointer) prp;

  if (StringCmp (newname, "unnamed protein product") == 0) {
    str = StringSave ("unnamed protein product");
    if (str == NULL) return;
    ValNodeAddPointer (&prp->name, 0, str);
    return;
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, newname) != 0) continue;
    if (StringHasNoText (gbq->val)) continue;
    str = StringSave (gbq->val);
    if (str == NULL) return;
    ValNodeAddPointer (&prp->name, 0, str);
    return;
  }

  /*
  str = StringSave ("unnamed protein product");
  */
  str = StringSave (newname);
  if (str == NULL) return;
  ValNodeAddPointer (&prp->name, 0, str);
}

// Used for -N f
static void FixFeatureIDs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Boolean         hasxref;
  SeqFeatPtr      matchsfp;
  SeqFeatXrefPtr  xref, matchxref;

  if (sfp == NULL) return;

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->id.choice == 0) continue;
    matchsfp = SeqMgrGetFeatureByFeatID (sfp->idx.entityID, NULL, NULL, xref, NULL);
    if (matchsfp == NULL) continue;
    hasxref = FALSE;
    for (matchxref = matchsfp->xref; matchxref != NULL; matchxref = matchxref->next) {
      if (matchxref->id.choice == 0) continue;
      hasxref = TRUE;
    }
    if (hasxref) continue;
    LinkTwoFeatures (sfp, matchsfp);
  }
}

static void MarkTitles (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  if (sdp->extended == 0) return;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}

static void MarkNucTitles (
  BioseqPtr bsp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;
  SeqDescrPtr    sdp;

  if (bsp == NULL) return;
  for (sdp = bsp->descr; sdp != NULL; sdp = sdp->next) {
    if (sdp->choice != Seq_descr_title) continue;
    if (sdp->extended == 0) continue;
    ovn = (ObjValNodePtr) sdp;
    ovn->idx.deleteme = TRUE;
  }
}

//LCOV_EXCL_START
static void MarkComments (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;

  if (sdp == NULL || sdp->choice != Seq_descr_comment) return;
  if (sdp->extended == 0) return;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}
//LCOV_EXCL_STOP

static void MarkUserObjects (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;
  CharPtr        str;
  UserObjectPtr  uop;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  if (sdp->extended == 0) return;

  str = (CharPtr) userdata;
  if (StringHasNoText (str)) return;

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL || uop->type == NULL) return;
  if (StringCmp (uop->type->str, str) != 0) return;

  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}

//LCOV_EXCL_START
// Only for SegSets
static void FindSegSeq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioseqPtr PNTR bspp;

  if (bsp == NULL || userdata == NULL) return;
  if (bsp->repr != Seq_repr_seg) return;

  bspp = (BioseqPtr PNTR) userdata;
  if (*bspp != NULL) return;

  *bspp = bsp;
}

// Only for SegSets
static void FindMitoChloroTitle (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  Char     ch;
  Int2     PNTR i2p;
  CharPtr  ptr, str;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;

  str = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (str)) return;

  if (userdata == NULL) return;
  i2p = (Int2 PNTR) userdata;

  ptr = StringStr (str, "; nuclear gene");
  if (ptr == NULL) return;
  ptr += 14;
  ch = *ptr;
  if (ch == 's') {
    ptr++;
    ch = *ptr;
  }
  while (ch == ' ') {
    ptr++;
    ch = *ptr;
  }

  if (StringNCmp (ptr, "for mitochondrial", 17) == 0) {
    *i2p = 1;
  } else if (StringNCmp (ptr, "for chloroplast", 15) == 0) {
    *i2p = 2;
  }
}

// Only for SegSets
static void CountSegSetFeatures (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Int4Ptr  featcountsp;
  Uint1    subtype;

  if (sfp == NULL || userdata == NULL) return;
  featcountsp = (Int4Ptr) userdata;

  subtype = sfp->idx.subtype;
  if (subtype > 0 && subtype < FEATDEF_MAX) {
    (featcountsp [subtype])++;
  }
}


static const Int4 s_auto_def_id_preferred_quals[] = {
 DEFLINE_POS_Strain,
 DEFLINE_POS_Clone,
 DEFLINE_POS_Isolate,
 DEFLINE_POS_Cultivar,
 DEFLINE_POS_Haplotype,
 DEFLINE_POS_Ecotype,
 DEFLINE_POS_Chromosome,
 DEFLINE_POS_Breed,
 DEFLINE_POS_Pop_variant,
 DEFLINE_POS_Segment,
 DEFLINE_POS_Specimen_voucher,
 DEFLINE_POS_Culture_collection
};

static const Int4 k_num_auto_def_id_preferred_quals = sizeof (s_auto_def_id_preferred_quals) / sizeof (Int4);

/* This function generates the modifiers for "AutoDefID" */
// For Autodef option, no longer used
static void Tbl_SetAutoDefIDModifiers (ModifierItemLocalPtr modList)
{
 Int4 index;
 Boolean added_required = FALSE;

 if (modList == NULL) return;
 /* first look for first modifier in list that is present on all sources */
 for (index = 0; index < k_num_auto_def_id_preferred_quals && !added_required; index++) {
   if (modList[s_auto_def_id_preferred_quals[index]].all_present) {
     modList[s_auto_def_id_preferred_quals[index]].required = TRUE;
     added_required = TRUE;
   }
 }
 /* if not found, then look for first modifier in list that is present on any sources */
 for (index = 0; index < k_num_auto_def_id_preferred_quals && !added_required; index++) {
   if (modList[s_auto_def_id_preferred_quals[index]].any_present) {
     modList[s_auto_def_id_preferred_quals[index]].required = TRUE;
     added_required = TRUE;
   }
 }
}

static void OnlyCommentFeats (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BoolPtr  notCommP;

  if (sfp == NULL) return;
  if (sfp->data.choice == SEQFEAT_COMMENT) return;
  if (sfp->idx.subtype == FEATDEF_repeat_region) return;
  notCommP = (BoolPtr) userdata;
  if (notCommP == NULL) return;
  *notCommP = TRUE;
}

// Use macro to add definition lines to a record instead
static void DoAutoDef (
  SeqEntryPtr sep,
  Uint2 entityID,
  Boolean specialSegSeq,
  Boolean autoDefID
)

{
  ValNodePtr                    defline_clauses = NULL;
  DeflineFeatureRequestList     feature_requests;
  Int4                          index;
  Int2                          mitochloroflag = DEFAULT_ORGANELLE_CLAUSE;
  ValNodePtr                    modifier_indices = NULL;
  ModifierItemLocalPtr          modList;
  OrganismDescriptionModifiers  odmp;
  SeqEntryPtr                   oldscope;
  Boolean                       notOnlyComm = FALSE;
  Int4                          featcounts [FEATDEF_MAX];

  if (sep == NULL) return;
  if (entityID < 1) return;

  VisitFeaturesInSep (sep, (Pointer) (&notOnlyComm), OnlyCommentFeats);
  if (! notOnlyComm) return;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;

  InitFeatureRequests (&feature_requests);

  if (specialSegSeq) {
// Only for SegSets
    VisitDescriptorsInSep (sep, (Pointer) &mitochloroflag, FindMitoChloroTitle);

    MemSet ((Pointer) featcounts, 0, sizeof (featcounts));
    VisitFeaturesInSep (sep, (Pointer) featcounts, CountSegSetFeatures);
    if (featcounts [FEATDEF_CDS] > 0) {
      ValNodeAddPointer (&(feature_requests.suppressed_feature_list), FEATDEF_exon, NULL);
      ValNodeAddPointer (&(feature_requests.suppressed_feature_list), FEATDEF_intron, NULL);
    }
  }

  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);

  if (autoDefID) {
    Tbl_SetAutoDefIDModifiers (modList);
  }

  InitOrganismDescriptionModifiers (&odmp, sep);
  odmp.use_modifiers = TRUE;

  RemoveNucProtSetTitles (sep);
  oldscope = SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, entityID, &feature_requests,
                              mitochloroflag, FALSE, FALSE,
                              &defline_clauses);

  if (AreFeatureClausesUnique (defline_clauses)) {
    modifier_indices = GetModifierIndicesFromModList (modList);
  } else {
    modifier_indices = FindBestModifiers (sep, modList);
  }

  AutoDefForSeqEntry(sep, entityID, &odmp, modList, modifier_indices, &feature_requests, mitochloroflag, FALSE, FALSE);

  DefLineFeatClauseListFree (defline_clauses);
  if (modList != NULL) {
    for (index = 0; index < NumDefLineModifiers (); index++) {
      ValNodeFree (modList [index].values_seen);
    }
    MemFree (modList);
  }
  modifier_indices = ValNodeFree (modifier_indices);

  SeqEntrySetScope (oldscope);
}

// Regenerate Definition Line using options object
static void RefreshDefLineCallback(BioseqPtr bsp, Pointer data)
{
    RegenerateAutoDef(bsp);
}
static void RegenerateDeflines(SeqEntryPtr sep)
{
    VisitBioseqsInSep(sep, NULL, RefreshDefLineCallback);
    NormalizeDescriptorOrder(sep);
}

// Used of -C c (unused)
static void CompressSeqs (BioseqPtr bsp, Pointer userdata)

{
  if (ISA_na (bsp->mol)) {
    BioseqConvert (bsp, Seq_code_ncbi2na);
  } else if (ISA_aa (bsp->mol)) {
    BioseqConvert (bsp, Seq_code_ncbieaa);
  }
}

// Used for -C d (unused)
static void ExpandSeqs (BioseqPtr bsp, Pointer userdata)

{
  if (ISA_na (bsp->mol)) {
    BioseqConvert (bsp, Seq_code_iupacna);
  } else if (ISA_aa (bsp->mol)) {
    BioseqConvert (bsp, Seq_code_ncbieaa);
  }
}
//LCOV_EXCL_STOP

// Used for -P j
static CharPtr FindISOJTA (CharPtr title)

{
  CharPtr       count = NULL, jids = NULL, iso = NULL, str = NULL;
  Boolean       is_issn = FALSE;
  ValNodePtr    head = NULL, tail = NULL;
  Int2          num;
  XmlObjPtr     nxt, sub, tmp, xop;

  if (StringHasNoText (title)) return NULL;

  ConvertToTermListForm (title);

  if (LooksLikeISSN (title)) {
    is_issn = TRUE;
    if (title [4] == '+' || title [4] == ' ') {
      title [4] = '-';
    }
    if (title [8] == 'x') {
      title [8] = 'X';
    }
    str = QUERY_UrlSynchronousQuery ("eutils.ncbi.nlm.nih.gov", 80,
                                     "/entrez/eutils/esearch.fcgi",
                                     "db=nlmcatalog&retmax=200&term=",
                                     title, "%5Bissn%5D", NULL);
  }

  if (str == NULL) {
    str = QUERY_UrlSynchronousQuery ("eutils.ncbi.nlm.nih.gov", 80,
                                     "/entrez/eutils/esearch.fcgi",
                                     "db=nlmcatalog&retmax=200&term=",
                                     title, "%5Bmulti%5D+AND+ncbijournals%5Bsb%5D", NULL);
  }

  if (str == NULL) return NULL;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return NULL;

  head = NULL;
  tail = NULL;
  count = NULL;

  for (tmp = xop; tmp != NULL; tmp = tmp->successor) {
    if (XmlPathSuffixIs (tmp, "/Id")) {
      ValNodeCopyStrEx (&head, &tail, 0, tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/eSearchResult/Count")) {
      count = StringSave (tmp->contents);
    }
  }

  FreeXmlObject (xop);

  if (StringCmp (count, "0") == 0 || head == NULL) {

    count = MemFree (count);
    ValNodeFreeData (head);

    /*
    if [multi] failed, try [jour]
    Microbiology (Reading, Engl.)
    is indexed in [multi] as
    microbiology reading, england
    and in [jour] as
    microbiology reading, engl
    microbiology reading, england
    */

    str = QUERY_UrlSynchronousQuery ("eutils.ncbi.nlm.nih.gov", 80,
                                     "/entrez/eutils/esearch.fcgi",
                                     "db=nlmcatalog&retmax=200&term=",
                                     title, "%5Bjour%5D", NULL);
    if (str == NULL) return NULL;

    xop = ParseXmlString (str);
    MemFree (str);
    if (xop == NULL) return NULL;

    head = NULL;
    tail = NULL;
    count = NULL;

    for (tmp = xop; tmp != NULL; tmp = tmp->successor) {
      if (XmlPathSuffixIs (tmp, "/Id")) {
        ValNodeCopyStrEx (&head, &tail, 0, tmp->contents);
      } else if (XmlPathSuffixIs (tmp, "/eSearchResult/Count")) {
        count = StringSave (tmp->contents);
      }
    }

    FreeXmlObject (xop);
  }

  MemFree (count);

  if (head != NULL) {
    jids = ValNodeMergeStrsEx (head, ",");
    ValNodeFreeData (head);
  }

  if (jids == NULL) return NULL;

  str = QUERY_UrlSynchronousQuery ("eutils.ncbi.nlm.nih.gov", 80,
                                   "/entrez/eutils/esummary.fcgi",
                                   "db=nlmcatalog&retmax=200&version=2.0&id=",
                                   jids, NULL, NULL);
  MemFree (jids);
  if (str == NULL) return NULL;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return NULL;

  iso = NULL;
  num = 0;

  for (tmp = xop; tmp != NULL; tmp = nxt) {
    nxt = tmp->successor;
    if (XmlPathSuffixIs (tmp, "/DocumentSummary")) {
      for (sub = tmp->successor; sub != NULL && sub->level > tmp->level; sub = sub->successor) {
        if (XmlPathSuffixIs (sub, "/ISOAbbreviation")) {
          if (num == 0) {
            iso = StringSave (sub->contents);
          }
          num++;
        }
      }
      nxt = sub;
    }
  }

  FreeXmlObject (xop);

  if (num > 1) {
    iso = MemFree (iso);
  }

  return iso;
}

// Used for -P j
static void LookupISOJTA (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CitArtPtr   cap;
  CitGenPtr   cgp;
  CitJourPtr  cjp;
  CharPtr     iso, str;
  ValNodePtr  jta = NULL, ttl, vnp;

  if (pdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp == NULL) continue;
        for (ttl = cgp->journal; ttl != NULL; ttl = ttl->next) {
          if (ttl->choice == Cit_title_iso_jta) return;
          if (ttl->choice == Cit_title_jta) {
            jta = ttl;
          }
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap == NULL) return;
        if (cap->from != 1) continue;
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp == NULL) continue;
        for (ttl = cjp->title; ttl != NULL; ttl = ttl->next) {
          if (ttl->choice == Cit_title_iso_jta) return;
          if (ttl->choice == Cit_title_jta) {
            jta = ttl;
          }
        }
        break;
      default :
        break;
    }
  }

  if (jta == NULL) return;

  str = (CharPtr) jta->data.ptrvalue;
  if (StringHasNoText (str)) return;

  iso = FindISOJTA (str);
  if (iso == NULL) return;

  jta->data.ptrvalue = MemFree (jta->data.ptrvalue);
  jta->data.ptrvalue = StringSave (iso);
  jta->choice = Cit_title_iso_jta;
}

// used for -P u
static void LookupPubdesc (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CitArtPtr        cap;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pep;
  Int4             pmid = 0;
  ValNodePtr       vnp;

  if (pdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Muid :
        /* ignore obsolete muids */
        break;
      case PUB_PMid :
        pmid = vnp->data.intvalue;
        break;
      default :
        /* return on real pub */
        return;
        break;
    }
  }

  if (pmid == 0) return;

  pep = GetPubMedForUid (pmid);
  if (pep == NULL) return;
  mep = (MedlineEntryPtr) pep->medent;
  if (mep != NULL && mep->cit != NULL) {
    cap = AsnIoMemCopy ((Pointer) mep->cit,
                        (AsnReadFunc) CitArtAsnRead,
                        (AsnWriteFunc) CitArtAsnWrite);
    ValNodeAddPointer (&(pdp->pub), PUB_Article, (Pointer) cap);
  }

  PubmedEntryFree (pep);
}

// used for -P #
static AuthListPtr AppendConsortiumToAuthList (
  AuthListPtr alp,
  CharPtr consortium
)

{
  AuthorPtr    ap;
  ValNodePtr   names;
  PersonIdPtr  pid;

  if (StringHasNoText (consortium)) return alp;
  if (alp == NULL) {
    alp = AuthListNew ();
    alp->choice = 1;
  }
  pid = PersonIdNew ();
  if (pid == NULL) return NULL;
  pid->choice = 5;
  pid->data = StringSave (consortium);
  ap = AuthorNew ();
  if (ap == NULL) return NULL;
  ap->name = pid;
  names = ValNodeAdd (&(alp->names));
  names->choice = 1;
  names->data.ptrvalue = ap;
  return alp;
}

// used for -P #
static void ReplaceUnpub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CharPtr          authors = NULL, consortium = NULL;
  CitArtPtr        cap;
  CitGenPtr        cgp = NULL;
  Int2             count = 0;
  Boolean          hasUnpublished = FALSE;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pep;
  Int4             pmid;
  Int4Ptr          pmP;
  ValNodePtr       vnp = NULL;

  if (pdp == NULL || userdata == NULL) return;
  pmP = (Int4Ptr) userdata;
  pmid = *pmP;
  if (pmid == 0) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    count++;
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            hasUnpublished = TRUE;
          }
        }
      }
    } else if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      return;
    } else if (vnp->choice == PUB_Article || vnp->choice == PUB_Book || vnp->choice == PUB_Man) {
      return;
    }
  }

  if (! hasUnpublished) return;

  /* assume only a cit-gen for now */
  if (count != 1 || cgp == NULL) return;

  authors = GetAuthorsString (GENBANK_FMT, cgp->authors, &consortium, NULL, NULL);
  MemFree (authors);

  pep = GetPubMedForUid (pmid);
  if (pep == NULL) return;
  mep = (MedlineEntryPtr) pep->medent;
  if (mep != NULL && mep->cit != NULL) {
    /* remove cit-gen */
    CitGenFree (cgp);
    ValNodeFree (pdp->pub);
    pdp->pub = NULL;

    ValNodeAddInt (&(pdp->pub), PUB_PMid, pmid);
    cap = AsnIoMemCopy ((Pointer) mep->cit,
                        (AsnReadFunc) CitArtAsnRead,
                        (AsnWriteFunc) CitArtAsnWrite);
    ValNodeAddPointer (&(pdp->pub), PUB_Article, (Pointer) cap);
    if (consortium != NULL) {
      AppendConsortiumToAuthList (cap->authors, consortium);
    }
  }

  MemFree (consortium);

  PubmedEntryFree (pep);
}

//LCOV_EXCL_START
// unused (selective and reports)
static void CleanupLocation (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr  bsp;
  SeqIntPtr  sintp;
  SeqLocPtr  slp;

  if (sfp == NULL || sfp->location == NULL) return;

  CleanUpSeqLoc (sfp->location);

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
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }
    }
  }
}

static void CleanupMostRNAs(
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  RnaRefPtr  rrp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL || rrp->type == 255) return;

  CleanUpSeqFeat (sfp, FALSE, FALSE, TRUE, FALSE, NULL);
}

static void CleanupRemainingRNAs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;

  CleanUpSeqFeat (sfp, FALSE, FALSE, TRUE, FALSE, NULL);
}

static void CleanupPubAuthors (
  PubdescPtr pdp,
  Pointer userdata
)

{
  if (pdp == NULL) return;

  CleanUpPubdescAuthors (pdp);
}

static void CleanupPubBody (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  CleanUpPubdescBody (pdp, cfp->stripSerial);
}

static void ModGenes (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeGeneFields (sfp);
}

static void ModRNAs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeRNAFields (sfp);
}

static void ModPCRs (
  BioSourcePtr biop,
  Pointer userdata
)

{
  ModernizePCRPrimers (biop);
}

static void StrucCommUop (
  UserObjectPtr uop,
  Pointer userdata
)

{
  CleanStructuredComment (uop);
}

static void CleanupStrucComment (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  UserObjectPtr  uop;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;

  VisitUserObjectsInUop (uop, NULL, StrucCommUop);
}


// Unused option
static void FixBspMol (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqDescrPtr  desc;
  MolInfoPtr   mip;

  if (bsp == NULL) return;

  desc = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (desc == NULL || desc->choice != Seq_descr_molinfo) return;
  mip = (MolInfoPtr) desc->data.ptrvalue;
  if (mip == NULL) return;
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

// Unused option
static ByteStorePtr Se2Bs (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("wb", bs);
  if (aibp == NULL || aibp->aip == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  return bs;
}
//LCOV_EXCL_STOP

static ByteStorePtr Se2BsX (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("wb", bs);
  if (aibp == NULL || aibp->aip == NULL) return NULL;

  aibp->aip->asn_no_newline = TRUE;
  aibp->aip->asn_alt_struct = TRUE;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  return bs;
}

/*
static CharPtr Se2Str (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;
  CharPtr       str;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  str = BSMerge (bs, NULL);
  BSFree (bs);

  return str;
}
*/

typedef struct chgdata {
  Boolean       isRefSeq;
  Boolean       sgml;
  Boolean       cdscodon;
  Boolean       conflict;
  Boolean       bad_fuzz;
  Boolean       bad_part;
  Boolean       rubisco;
  Boolean       rbc;
  Boolean       its;
  Boolean       rnaother;
  Boolean       trnanote;
  Boolean       oldbiomol;
  Boolean       oldgbqual;
  Boolean       badDbxref;
  Boolean       refDbxref;
  Boolean       srcDbxref;
  Boolean       capDbxref;
  Boolean       oldDbxref;
  Boolean       privDbxref;
  Boolean       multDbxref;
  Boolean       rareDbxref;
  Boolean       gbDbxref;
  Boolean       badOrg;
  Boolean       rpt_unit_seq;
  Boolean       hasUnpublished;
  Boolean       hasPublished;
  Boolean       protNucSet;
  Boolean       noPopsetTitle;
  Int4          protdesc;
  Int4          sfpnote;
  Int4          gbsource;
  Int4          cdsconf;
} ChangeData, PNTR ChangeDataPtr;

//LCOV_EXCL_START
// Used in ScoreFeature
static Boolean IsRubisco (
  CharPtr name
)

{
  return (StringICmp (name, "rubisco large subunit") == 0 ||
          StringICmp (name, "rubisco small subunit") == 0);
}

// Used in ScoreFeature
static Boolean IsRbc (
  CharPtr name
)

{
  return (StringICmp (name, "RbcL") == 0 ||
          StringICmp (name, "RbcS") == 0);
}

// Used in ScoreFeature
static Boolean IsITS (
  CharPtr name
)

{
  return (StringICmp (name, "its1") == 0 ||
          StringICmp (name, "its 1") == 0 ||
          StringICmp (name, "its2") == 0 ||
          StringICmp (name, "its 2") == 0 ||
          StringICmp (name, "its3") == 0 ||
          StringICmp (name, "its 3") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 1") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 2") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 3") == 0 ||
          StringICmp (name, "internal transcribed spacer 1 (ITS1)") == 0 ||
          StringICmp (name, "internal transcribed spacer 2 (ITS2)") == 0 ||
          StringICmp (name, "internal transcribed spacer 3 (ITS3)") == 0);
}

// Used in ScoreFeature (unused)
static Boolean HasSgml (
  CharPtr str
)

{
  Int2  ascii_len;
  Char  buf [1024];

  if (StringHasNoText (str)) return FALSE;

  ascii_len = Sgml2AsciiLen (str);
  if (ascii_len + 2 > sizeof (buf)) return FALSE;

  Sgml2Ascii (str, buf, ascii_len + 1);
  if (StringCmp (str, buf) != 0) {
    return TRUE;
  }

  return FALSE;
}

static void LookForBadDbxref (
  ValNodePtr list,
  ChangeDataPtr cdp,
  Boolean isSource
)

{
  Boolean      cap;
  DbtagPtr     dp;
  CharPtr      good;
  ObjectIdPtr  oip;
  Boolean      ref;
  Boolean      src;
  CharPtr      str;
  ValNodePtr   vnp;

  if (list == NULL || cdp == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    dp = (DbtagPtr) vnp->data.ptrvalue;
    if (dp != NULL && StringDoesHaveText (dp->db)) {

      oip = dp->tag;
      if (oip != NULL && StringDoesHaveText (oip->str)) {
        if (StringChr (oip->str, ':') != NULL) {
          cdp->multDbxref = TRUE;
        }
      }

      str = dp->db;
      if (StringICmp (str, "PID") == 0 ||
          StringICmp (str, "PIDg") == 0 ||
          StringICmp (str, "PIDd") == 0 ||
          StringICmp (str, "PIDe") == 0 ||
          StringICmp (str, "NID") == 0 ||
          StringICmp (str, "GI") == 0) {
        cdp->privDbxref = TRUE;
        continue;
      }
      if (StringICmp (str, "SWISS-PROT") == 0 ||
          StringICmp (str, "SWISSPROT") == 0 ||
          StringICmp (str, "SPTREMBL") == 0 ||
          StringICmp (str, "SUBTILIS") == 0 ||
          StringICmp (str, "MGD") == 0 ||
          StringCmp (str, "cdd") == 0 ||
          StringICmp (str, "TrEMBL") == 0 ||
          StringICmp (str, "LocusID") == 0 ||
          StringICmp (str, "MaizeDB") == 0 ||
          StringICmp (str, "UniProt/Swiss-Prot") == 0 ||
          StringICmp (str, "UniProt/TrEMBL") == 0 ||
          StringICmp (str, "Genew") == 0 ||
          StringICmp (str, "GENEDB") == 0 ||
          StringICmp (str, "GreengenesID") == 0 ||
          StringICmp (str, "HMPID") == 0 ||
          StringICmp (str, "IFO") == 0 ||
          StringICmp (str, "BHB") == 0 ||
          StringICmp (str, "BioHealthBase") == 0) {
        cdp->oldDbxref = TRUE;
        continue;
      }
      if (StringICmp (str, "ATCC(dna)") == 0 ||
          StringICmp (str, "ATCC(in host)") == 0 ||
          StringICmp (str, "BDGP_EST") == 0 ||
          StringICmp (str, "BDGP_INS") == 0 ||
          StringICmp (str, "CGNC") == 0 ||
          StringICmp (str, "CloneID") == 0 ||
          StringICmp (str, "ENSEMBL") == 0 ||
          StringICmp (str, "ESTLIB") == 0 ||
          StringICmp (str, "GDB") == 0 ||
          /*
          StringICmp (str, "GOA") == 0 ||
          */
          StringICmp (str, "IMGT/HLA") == 0 ||
          StringICmp (str, "PIR") == 0 ||
          StringICmp (str, "PSEUDO") == 0 ||
          StringICmp (str, "RZPD") == 0 ||
          StringICmp (str, "SoyBase") == 0 ||
          StringICmp (str, "UNILIB") == 0) {
        cdp->rareDbxref = TRUE;
        continue;
      }
      if (StringICmp (str, "genbank") == 0) {
        cdp->gbDbxref = TRUE;
        continue;
      }
      if (StringICmp (str, "MGD") == 0 || StringICmp (str, "MGI") == 0) {
        oip = dp->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "MGI:", 4) == 0 || StringNICmp (str, "MGD:", 4) == 0) {
            cdp->oldDbxref = TRUE;
            continue;
          }
        }
      } else if (StringICmp (str, "HPRD") == 0) {
        oip = dp->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "HPRD_", 5) == 0) {
            cdp->oldDbxref = TRUE;
            continue;
          }
        }
      }

      if (isSource && StringCmp (str, "taxon") == 0) continue;

      if (DbxrefIsValid (str, &ref, &src, &cap, &good)) {
        if (ref && (! cdp->isRefSeq)) {
          cdp->refDbxref = TRUE;
        }
        if (isSource && (! src)) {
          cdp->srcDbxref = TRUE;
        }
        if (cap) {
          cdp->capDbxref = TRUE;
        }
      } else {
        cdp->badDbxref = TRUE;
      }
    }
  }
}

// Used in CheckForChanges (unused)
static void ScoreFeature (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  ChangeDataPtr  cdp;
  Char           ch;
  CharPtr        comment;
  CdRegionPtr    crp;
  CharPtr        desc;
  GBQualPtr      gbq;
  GeneRefPtr     grp;
  ImpFeatPtr     ifp;
  CharPtr        name;
  OrgRefPtr      orp;
  ProtRefPtr     prp;
  CharPtr        ptr;
  Uint1          residue;
  RnaRefPtr      rrp;
  CharPtr        str;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
   cdp = (ChangeDataPtr) userdata;
   if (cdp == NULL) return;

  comment = sfp->comment;
  if (StringDoesHaveText (comment)) {
    (cdp->sfpnote)++;
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringCmp (gbq->qual, "partial") == 0 ||
        StringCmp (gbq->qual, "evidence") == 0 ||
        StringCmp (gbq->qual, "exception") == 0 ||
        StringCmp (gbq->qual, "note") == 0 ||
        StringCmp (gbq->qual, "notes") == 0 ||
        StringCmp (gbq->qual, "comment") == 0 ||
        StringCmp (gbq->qual, "db_xref") == 0 ||
        StringCmp (gbq->qual, "gdb_xref") == 0 ||
        StringCmp (gbq->qual, "rpt_unit") == 0 ||
        StringCmp (gbq->qual, "pseudo") == 0 ||
        StringCmp (gbq->qual, "gene") == 0 ||
        StringCmp (gbq->qual, "codon_start") == 0 ||
        StringCmp (gbq->qual, "transposon") == 0 ||
        StringCmp (gbq->qual, "insertion_seq") == 0) {
      cdp->oldgbqual = TRUE;
    } else if (StringICmp (gbq->qual, "rpt_unit_seq") == 0) {
      if (StringHasNoText (gbq->val)) continue;
      ptr = gbq->val;
      ch = *ptr;
      while (ch != '\0') {
        if (IS_UPPER (ch)) {
          cdp->rpt_unit_seq = TRUE;
        }
        ptr++;
        ch = *ptr;
      }
    }
  }

  LookForBadDbxref (sfp->dbxref, cdp, FALSE);

  if (SeqLocPartialCheckEx (sfp->location, FALSE) != SLP_COMPLETE) {
    if (! sfp->partial) {
      cdp->bad_fuzz = TRUE;
    }
  } else {
    if (sfp->partial) {
      cdp->bad_part = TRUE;
    }
  }

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
    case SEQFEAT_GENE:
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (HasSgml (grp->locus)) {
        cdp->sgml = TRUE;
      }
      if (HasSgml (grp->desc)) {
        cdp->sgml = TRUE;
      }
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (HasSgml (str)) {
          cdp->sgml = TRUE;
        }
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "map") == 0 ||
            StringCmp (gbq->qual, "allele") == 0 ||
            StringCmp (gbq->qual, "locus_tag") == 0 ||
            StringCmp (gbq->qual, "old_locus_tag") == 0) {
          cdp->oldgbqual = TRUE;
        }
      }
      LookForBadDbxref (grp->db, cdp, FALSE);
      break;
    case SEQFEAT_CDREGION:
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp->conflict) {
        (cdp->cdsconf)++;
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "codon") != 0) continue;
        if (StringHasNoText (gbq->val)) continue;
        cdp->cdscodon = TRUE;
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "product") == 0 ||
            StringCmp (gbq->qual, "function") == 0 ||
            StringCmp (gbq->qual, "EC_number") == 0 ||
            StringCmp (gbq->qual, "prot_note") == 0) {
          cdp->oldgbqual = TRUE;
        }
      }
      break;
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      desc = prp->desc;
      if (StringDoesHaveText (desc)) {
        (cdp->protdesc)++;
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (IsRubisco (str)) {
          cdp->rubisco = TRUE;
        }
        if (IsRbc (str)) {
          cdp->rbc = TRUE;
        }
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "product") == 0 ||
            StringCmp (gbq->qual, "function") == 0 ||
            StringCmp (gbq->qual, "EC_number") == 0 ||
            StringCmp (gbq->qual, "label") == 0 ||
            StringCmp (gbq->qual, "allele") == 0) {
          cdp->oldgbqual = TRUE;
        }
        if (StringCmp (gbq->qual, "standard_name") == 0 && prp->name == NULL) {
          cdp->oldgbqual = TRUE;
        }
      }
      LookForBadDbxref (prp->db, cdp, FALSE);
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
            if (IsITS (name)) {
              cdp->its = TRUE;
            }
          }
        } else if (StringCmp (name, "ncRNA") == 0 || StringCmp (name, "tmRNA") == 0) {
        } else {
          cdp->rnaother = TRUE;
          if (IsITS (name)) {
            cdp->its = TRUE;
          }
        }
      } else if (rrp->type == 3 && rrp->ext.choice == 2) {
        if (StringDoesHaveText (comment)) {
          if (StringNCmp (comment, "aa: ", 4) == 0) {
            comment += 4;
          }
          residue = FindTrnaAA3 (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
          }
          residue = FindTrnaAA (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
          }
        }
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "product") == 0 ||
            StringCmp (gbq->qual, "ncRNA_class") == 0 ||
            StringCmp (gbq->qual, "tag_peptide") == 0) {
          cdp->oldgbqual = TRUE;
        }
      }
      break;
    case SEQFEAT_IMP :
      ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (StringICmp (ifp->key, "conflict") == 0) {
        cdp->conflict = TRUE;
      }
      break;
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      LookForBadDbxref (orp->db, cdp, TRUE);
      cdp->badOrg = TRUE;
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      orp = biop->org;
      if (orp != NULL) {
        LookForBadDbxref (orp->db, cdp, TRUE);
      }
    default:
      break;
  }
}

// Used in CheckForChanges (unused)
static void ScoreDescriptor (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  ChangeDataPtr  cdp;
  GBBlockPtr     gbp;
  MolInfoPtr     mip;
  OrgRefPtr      orp;

  if (sdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL) {
        if (StringDoesHaveText (gbp->source)) {
          (cdp->gbsource)++;
        }
      }
      break;
    case Seq_descr_molinfo :
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        switch (mip->biomol) {
          case MOLECULE_TYPE_SNRNA:
          case MOLECULE_TYPE_SCRNA:
          case MOLECULE_TYPE_SNORNA:
            cdp->oldbiomol = TRUE;
            break;
          default :
            break;
        }
      }
      break;
    case Seq_descr_org :
      orp = (OrgRefPtr) sdp->data.ptrvalue;
      if (orp != NULL) {
        LookForBadDbxref (orp->db, cdp, TRUE);
      }
      cdp->badOrg = TRUE;
      break;
    case Seq_descr_source :
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          LookForBadDbxref (orp->db, cdp, TRUE);
        }
      }
      break;
    default :
      break;
  }
}


// Used in CheckForChanges (unused)
static void CheckForUnpubPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  CitGenPtr      cgp;
  ValNodePtr     vnp;

  if (pdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            cdp->hasUnpublished = TRUE;
          }
        }
      }
    } else if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      cdp->hasPublished = TRUE;
    } else if (vnp->choice == PUB_Article || vnp->choice == PUB_Book || vnp->choice == PUB_Man) {
      cdp->hasPublished = TRUE;
    }
  }
}


// Used in CheckForChanges (unused)
static void CheckForSetTitle (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  BioseqPtr      bsp;
  ChangeDataPtr  cdp;
  SeqDescrPtr    sdp;
  SeqEntryPtr    sep;
  CharPtr        title;

  if (bssp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  if (bssp->_class == BioseqseqSet_class_nuc_prot) {
    sep = bssp->seq_set;
    if (sep != NULL && IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp != NULL) {
        if (ISA_aa (bsp->mol)) {
          cdp->protNucSet = TRUE;
        }
      }
    }
  } else if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) {
    for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
      if (sdp->choice != Seq_descr_title) continue;
      title = (CharPtr) sdp->data.ptrvalue;
      if (StringDoesHaveText (title)) continue;
      cdp->noPopsetTitle = TRUE;
    }
  }
}

// Used in DoASNReport (unused)
static void CheckForChanges (
  SeqEntryPtr sep,
  ChangeDataPtr cdp
)

{
  if (sep == NULL || cdp == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) cdp, ScoreFeature);
  VisitDescriptorsInSep (sep, (Pointer) cdp, ScoreDescriptor);
  VisitPubdescsInSep (sep, (Pointer) cdp, CheckForUnpubPub);
  VisitSetsInSep (sep, (Pointer) cdp, CheckForSetTitle);
}
//LCOV_EXCL_STOP

static void StripBadProtTitles (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr            buf;
  size_t             buflen = 1001;
  ObjValNodePtr      ovp;
  SeqIdPtr           sip;
  CharPtr            title;
  ValNodePtr         vnp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) return;
  }

  vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (vnp == NULL) return;
  title = (CharPtr) vnp->data.ptrvalue;
  if (StringHasNoText (title)) return;

  buf = MemNew (sizeof (Char) * (buflen + 1));
  if (buf == NULL) return;

  if (NewCreateDefLineBuf (NULL, bsp, buf, buflen, TRUE, FALSE)) {
    if (StringICmp (buf, title) != 0) {
      if (vnp->extended != 0) {
        ovp = (ObjValNodePtr) vnp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }

  MemFree (buf);
}

static void BadProtTitleProc (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;
  VisitBioseqsInSep (sep, NULL, StripBadProtTitles);
}

//LCOV_EXCL_START
static void BadNCTitleProc (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqPtr      bsp;
  Boolean        is_nc;
  ObjValNodePtr  ovp;
  SeqIdPtr       sip;
  TextSeqIdPtr   tsip;
  ValNodePtr     vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;

  is_nc = FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
          is_nc = TRUE;
        }
      }
    }
  }
  if (! is_nc) return;

  vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (vnp == NULL) return;
  if (vnp->extended != 0) {
    ovp = (ObjValNodePtr) vnp;
    ovp->idx.deleteme = TRUE;
  }
}

static void BadNMTitleProc (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqPtr      bsp;
  Boolean        is_nm;
  ObjValNodePtr  ovp;
  SeqIdPtr       sip;
  TextSeqIdPtr   tsip;
  ValNodePtr     vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;

  is_nm = FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNICmp (tsip->accession, "NM_", 3) == 0) {
          is_nm = TRUE;
        } else if (StringNICmp (tsip->accession, "XM_", 3) == 0) {
          is_nm = TRUE;
        }
      }
    }
  }
  if (! is_nm) return;

  vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (vnp == NULL) return;
  if (vnp->extended != 0) {
    ovp = (ObjValNodePtr) vnp;
    ovp->idx.deleteme = TRUE;
  }
}

typedef struct xmfeatdata {
  SeqFeatPtr  gene;
  SeqFeatPtr  cds;
  Int2        numgenes;
  Int2        numcds;
  Int2        numprots;
} XmFeatData, PNTR XmFeatPtr;

static void FindXMFeats (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  XmFeatPtr  xfp;

  if (sfp == NULL) return;
  xfp = (XmFeatPtr) userdata;
  if (xfp == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      xfp->gene = sfp;
      (xfp->numgenes)++;
      break;
    case SEQFEAT_CDREGION :
      xfp->cds = sfp;
      (xfp->numcds++);
      break;
    case SEQFEAT_PROT :
      (xfp->numprots)++;
      break;
    default :
      break;
  }
}


static CharPtr GetNmTaxname (BioseqPtr bsp)

{
  BioSourcePtr  biop;
  OrgRefPtr     orp;
  SeqDescrPtr   sdp;

  if (bsp == NULL) return NULL;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (sdp == NULL) return NULL;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return NULL;
  orp = biop->org;
  if (orp == NULL) return NULL;
  if (StringHasNoText (orp->taxname)) return NULL;

  return orp->taxname;
}

static CharPtr CreateSpecialXmTitle (BioseqPtr bsp)

{
  Char         buf [512];
  CharPtr      cds = NULL, gene = NULL, result = NULL, taxname = NULL;
  Uint2        entityID;
  size_t       len;
  SeqEntryPtr  sep;
  XmFeatData   xfd;

  if (bsp == NULL) return NULL;

  taxname = GetNmTaxname (bsp);
  if (StringHasNoText (taxname)) return NULL;

  MemSet ((Pointer) &xfd, 0, sizeof (XmFeatData));

  entityID = ObjMgrGetEntityIDForPointer (bsp);
  sep = GetBestTopParentForDataEx (entityID, bsp, TRUE);

  VisitFeaturesInSep (sep, (Pointer) &xfd, FindXMFeats);
  if (xfd.numgenes != 1 || xfd.numcds != 1 || xfd.numprots < 1) return NULL;

  FeatDefLabel (xfd.gene, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);
  gene = StringSaveNoNull (buf);

  FeatDefLabel (xfd.cds, buf, sizeof (buf) - 1, OM_LABEL_CONTENT);
  cds = StringSaveNoNull (buf);

  len = StringLen (taxname) + StringLen (cds) +
        StringLen (gene) + StringLen ("  (), partial mRNA") + 10;

  result = (CharPtr) MemNew (sizeof (Char) * len);

  if (result != NULL) {
    if (xfd.cds != NULL && xfd.cds->partial) {
      sprintf (result, "%s %s partial mRNA", taxname, cds);
    } else {
      sprintf (result, "%s %s mRNA", taxname, cds);
    }
  }

  MemFree (gene);
  MemFree (cds);

  return result;
}

static Boolean AddSpecialXmTitles(GatherObjectPtr gop)
{
  BioseqPtr     bsp;
  Boolean       is_nm;
  SeqIdPtr      sip;
  CharPtr       str;
  TextSeqIdPtr  tsip;

  if (gop == NULL || gop->itemtype != OBJ_BIOSEQ) return TRUE;
  bsp = (BioseqPtr) gop->dataptr;
  if (bsp == NULL) return TRUE;
  is_nm = FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringNICmp (tsip->accession, "NM_", 3) == 0) {
          is_nm = TRUE;
        } else if (StringNICmp (tsip->accession, "XM_", 3) == 0) {
          is_nm = TRUE;
        }
      }
    }
  }
  if (! is_nm) return TRUE;
  str = CreateSpecialXmTitle (bsp);
  if (str == NULL) return TRUE;
  SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) str);
  return TRUE;
}

static void InstantiateSpecialXMTitles (
  Uint2 entityID,
  Pointer ptr
)

{
  Boolean  objMgrFilt [OBJ_MAX];

  if (entityID == 0) {
    entityID = ObjMgrGetEntityIDForPointer (ptr);
  }
  if (entityID == 0) return;

  AssignIDsInEntity (entityID, 0, NULL);
  MemSet ((Pointer) objMgrFilt, FALSE, sizeof (objMgrFilt));
  objMgrFilt [OBJ_BIOSEQ] = TRUE;
  GatherObjectsInEntity (entityID, 0, NULL, AddSpecialXmTitles, NULL, objMgrFilt);
}

static void InstantiateGpipeTitles (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr  title;

  if (bsp == NULL) return;
  if (ISA_aa (bsp->mol)) return;

  title = NewCreateDefLineEx (NULL, bsp, FALSE, FALSE, TRUE);
  if (title == NULL) return;

  SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, (Pointer) title);
}

static void BSSaveToFile (
  ByteStorePtr bs,
  CharPtr path
)

{
  Byte  buf [256];
  Int4  count;
  FILE  *fp;

  if (bs == NULL || StringHasNoText (path)) return;

  fp = FileOpen (path, "w");
  if (fp != NULL) {
    Nlm_BSSeek (bs, 0, SEEK_SET);
    count = BSRead (bs, buf, sizeof (buf));
    while (count > 0) {
      FileWrite (buf, count, 1, fp);
      count = BSRead (bs, buf, sizeof (buf));
    }
    FileClose (fp);
  }
}

typedef struct diffblock {
  ValNodePtr  head;
  ValNodePtr  tail;
} DiffBlock, PNTR DiffBlockPtr;

// used by ReportAsnDiffs (unused)
static void RecordDiffBlock (
  DiffBlockPtr dbp,
  CharPtr str
)

{
  ValNodePtr  vnp;

  if (dbp == NULL || StringHasNoText (str)) return;

  vnp = ValNodeCopyStr (&(dbp->tail), 0, str);
  if (dbp->head == NULL) {
    dbp->head = vnp;
  }
  dbp->tail = vnp;
}

static void WriteDiffBlock (
  DiffBlockPtr dbp,
  FILE *fp
)

{
  Char        ch;
  Int2        idx;
  Int2        margin = INT2_MAX;
  CharPtr     ptr;
  Int2        spaces;
  CharPtr     str;
  ValNodePtr  vnp;

  if (dbp == NULL || dbp->head == NULL || fp == NULL) return;

  for (vnp = dbp->head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    ch = str [0];
    if (ch == '<' || ch == '>') {
      ptr = str + 1;
      ch = *ptr;
      spaces = 0;
      while (ch == ' ') {
        spaces++;
        ptr++;
        ch = *ptr;
      }
      if (spaces < margin) {
        margin = spaces;
      }
    }
  }

  if (margin > 80) {
    margin = 80;
  }

  for (vnp = dbp->head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    ch = str [0];
    if (ch == '<' || ch == '>') {
      ptr = str + 1;
      ch = *ptr;
      idx = 0;
      while (idx < margin && ch == ' ') {
        idx++;
        ptr++;
        ch = *ptr;
      }
      fprintf (fp, "%c %s\n", str [0], ptr);
    } else if (ch == '-') {
      fprintf (fp, "---\n");
    } else if (ch == '=') {
      fprintf (fp, "===\n");
    }
  }

  fprintf (fp, "\n");
  fflush (fp);
}

static void ResetDiffBlock (
  DiffBlockPtr dbp
)

{
  if (dbp == NULL) return;

  dbp->head = ValNodeFreeData (dbp->head);
  dbp->tail = NULL;
}

// Unused option
static void ReportAsnDiffs (
  FILE *logfp,
  CharPtr id,
  ByteStorePtr bs1,
  ByteStorePtr bs2)

{
#ifdef OS_UNIX
  Char       ch;
  Char       cmmd [512];
  DiffBlock  db;
  int        diff;
  FileCache  fc;
  FILE       *fp;
  Char       line [512];
  Char       path1 [PATH_MAX];
  Char       path2 [PATH_MAX];
  Char       path3 [PATH_MAX];
  CharPtr    str;

  if (logfp == NULL || StringHasNoText (id)) return;
  if (bs1 == NULL || bs2 == NULL) return;

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);

  BSSaveToFile (bs1, path1);
  BSSaveToFile (bs2, path2);

  db.head = NULL;
  db.tail = NULL;

  sprintf (cmmd, "diff -b -h %s %s > %s", path1, path2, path3);
  diff = system (cmmd);

  if (diff > 0) {
    fp = FileOpen (path3, "r");
    if (fp != NULL) {
      fprintf (logfp, "\n\n%s\n\n", id);
      fflush (logfp);
      if (FileCacheSetup (&fc, fp)) {
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
        while (str != NULL) {
          ch = line [0];
          if (ch == '<' || ch == '>') {
            RecordDiffBlock (&db, line);
          } else if (ch == '-') {
            RecordDiffBlock (&db, "---");
          } else if (IS_DIGIT (ch)) {
            WriteDiffBlock (&db, logfp);
            ResetDiffBlock (&db);
            RecordDiffBlock (&db, "===");
          } else if (StringHasNoText (str)) {
            WriteDiffBlock (&db, logfp);
            ResetDiffBlock (&db);
          }
          str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
        }
        WriteDiffBlock (&db, logfp);
        ResetDiffBlock (&db);
      }
      fprintf (logfp, "//\n\n");
      FileClose (fp);
    }
  }

  sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
  system (cmmd);
#endif
}

static void LogPopSetTitleResults (FILE *fp, PopSetRetroStatPtr stat, SeqEntryPtr sep)
{
  ValNode vn;
  CharPtr txt;

  if (fp == NULL || stat == NULL || sep == NULL || !IS_Bioseq_set (sep)) {
    return;
  }
  /* if test wasn't run, return */
  if (stat->feature_clause == 0 && stat->common_title == 0 && stat->uncalculatable == 0) {
    return;
  }

  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = OBJ_BIOSEQSET;
  vn.data.ptrvalue = sep->data.ptrvalue;
  txt = GetDiscrepancyItemText (&vn);
  fprintf (fp, "AutoDefPopSetResults: %s", txt);
  txt = MemFree (txt);
  fprintf (fp, "\tFeature Clause: %d\n\tCommon Title: %d\n\tUncalculatable: %d\n",
           stat->feature_clause, stat->common_title, stat->uncalculatable);
}

static void DoAsnDiffReport(
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  ByteStorePtr         bs = NULL, tmp = NULL;
  Uint2                entityID;
  Boolean              okay = FALSE;
  PopSetRetroStatData  stat;

  if (sep == NULL || cfp == NULL) return;

  MemSet ((Pointer) &stat, 0, sizeof (PopSetRetroStatData));

  RemoveAllNcbiCleanupUserObjects (sep);

  entityID = ObjMgrGetEntityIDForChoice (sep);

  /* Capital letters avoid unwanted diffs on issues already fixed in ID */

  if (StringChr (cfp->selective, 'S') != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
  }
  if (StringChr (cfp->selective, 'B') != NULL) {
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->selective, 'A') != NULL) {
    VisitPubdescsInSep (sep, NULL, CleanupPubAuthors);
  }
  if (StringChr (cfp->selective, 'P') != NULL) {
    VisitPubdescsInSep (sep, (Pointer) cfp, CleanupPubBody);
  }
  if (StringChr (cfp->selective, 'L') != NULL) {
    VisitFeaturesInSep (sep, NULL, CleanupLocation);
  }
  if (StringChr (cfp->selective, 'R') != NULL) {
    VisitFeaturesInSep (sep, NULL, CleanupMostRNAs);
    VisitFeaturesInSep (sep, NULL, CleanupRemainingRNAs);
    VisitFeaturesInSep (sep, NULL, ModRNAs);
  }
  if (StringChr (cfp->selective, 'Q') != NULL) {
    SortSeqEntryQualifiers (sep);
  }
  if (StringChr (cfp->selective, 'G') != NULL) {
    EntryChangeGBSource (sep);
    EntryCheckGBBlock (sep);
  }
  if (StringChr (cfp->selective, 'K') != NULL) {
    MoveFeatsFromPartsSet (sep);
    move_cds_ex (sep, TRUE);
  }
  if (StringChr (cfp->selective, 'M') != NULL) {
    SeqEntryPubsAsn4 (sep, cfp->isEmblDdbj);
  }
  if (StringChr (cfp->selective, 'O') != NULL) {
    SeqEntryPubsAsn4Ex (sep, cfp->isEmblDdbj, FALSE);
  }
  if (StringChr (cfp->selective, 'D') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID, FALSE, FALSE);
  }
  if (StringChr (cfp->selective, 'E') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    PopSetAutoDefRetro (sep, &stat);
    LogPopSetTitleResults (cfp->logfp, &stat, sep);
  }

  NormalizeDescriptorOrder (sep);

  /* Look for change in single issue */

  bs = Se2BsX (sep);
  if (StringHasNoText (cfp->selective)) {
    okay = TRUE;
  } else if (StringChr (cfp->selective, 's') != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'b') != NULL) {
    BasicSeqEntryCleanup (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'a') != NULL) {
    VisitPubdescsInSep (sep, NULL, CleanupPubAuthors);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'p') != NULL) {
    VisitPubdescsInSep (sep, (Pointer) cfp, CleanupPubBody);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'l') != NULL) {
    VisitFeaturesInSep (sep, NULL, CleanupLocation);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'r') != NULL) {
    VisitFeaturesInSep (sep, NULL, CleanupMostRNAs);
    VisitFeaturesInSep (sep, NULL, CleanupRemainingRNAs);
    VisitFeaturesInSep (sep, NULL, ModRNAs);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'q') != NULL) {
    SortSeqEntryQualifiers (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'g') != NULL) {
    EntryChangeGBSource (sep);
    EntryCheckGBBlock (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'k') != NULL) {
    MoveFeatsFromPartsSet (sep);
    move_cds_ex (sep, TRUE);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'm') != NULL) {
    SeqEntryPubsAsn4 (sep, cfp->isEmblDdbj);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'o') != NULL) {
    SeqEntryPubsAsn4Ex (sep, cfp->isEmblDdbj, FALSE);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'd') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID, FALSE, FALSE);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'e') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    PopSetAutoDefRetro (sep, &stat);
    LogPopSetTitleResults (cfp->logfp, &stat, sep);
    NormalizeDescriptorOrder (sep);
    if (stat.title_added) {
      okay = TRUE;
    }
  } else {
    okay = TRUE;
  }

  /* Report incremental diff */

  if (okay) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      if (cfp->logfp != NULL) {
        ReportAsnDiffs (cfp->logfp, cfp->buf, bs, tmp);
      }
    }
    tmp = BSFree (tmp);
  }

  BSFree (bs);
}
//LCOV_EXCL_STOP

static void GetPopPhyMutEcoType (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CharPtr PNTR  set_type;

  if (bssp == NULL || userdata == NULL) return;
  set_type = (CharPtr PNTR) userdata;

  if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) {
    switch (bssp->_class) {
      case BioseqseqSet_class_mut_set :
        *set_type = "MUT";
        break;
      case BioseqseqSet_class_pop_set :
        *set_type = "POP";
        break;
      case BioseqseqSet_class_phy_set :
        *set_type = "PHY";
        break;
      case BioseqseqSet_class_eco_set :
        *set_type = "ECO";
        break;
      default :
        break;
    }
  }
}

/*
static void ChromosomeScanCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioSourcePtr       biop = NULL;
  CleanFlagPtr       cfp;
  Boolean            complete = FALSE;
  SeqMgrDescContext  dcontext;
  CharPtr            genome = NULL;
  MolInfoPtr         mip = NULL;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;

  if (bsp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }
  if (biop == NULL) return;
  if (biop->genome == GENOME_chromosome) {
    genome = "CHROM";
  } else if (biop->genome == GENOME_genomic) {
    genome = "GENOM";
  } else if (biop->genome == GENOME_unknown) {
    genome = "UNKWN";
  } else {
    return;
  }
  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;
  if (StringNICmp (onp->lineage, "Eukaryota; ", 11) != 0) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
  }
  if (mip != NULL && mip->completeness == 1) {
    complete = TRUE;
  }

  if (cfp->logfp != NULL) {
    if (complete) {
      fprintf (cfp->logfp, "COMP %s %s\n", genome, cfp->buf);
    } else {
      fprintf (cfp->logfp, "PART %s %s\n", genome, cfp->buf);
    }
    fflush (cfp->logfp);
  }
}
*/

/*
static void PopPhyMutEcoScanCallback (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  SeqDescrPtr   sdp;
  CharPtr       set_type = NULL;
  CharPtr       title = NULL;

  if (bssp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) {
    GetPopPhyMutEcoType (bssp, (Pointer) &set_type);
    switch (bssp->_class) {
      case BioseqseqSet_class_mut_set :
        set_type = "MUT";
        break;
      case BioseqseqSet_class_pop_set :
        set_type = "POP";
        break;
      case BioseqseqSet_class_phy_set :
        set_type = "PHY";
        break;
      case BioseqseqSet_class_eco_set :
        set_type = "ECO";
        break;
      default :
        break;
    }
    for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
      if (sdp->choice != Seq_descr_title) continue;
      title = (CharPtr) sdp->data.ptrvalue;
    }
    if (StringHasNoText (set_type)) {
      set_type = "UNK";
    }
    if (StringDoesHaveText (title)) {
      fprintf (cfp->logfp, "%s set TITLE %s\n", set_type, cfp->buf);
    } else {
      fprintf (cfp->logfp, "%s set ANNON %s\n", set_type, cfp->buf);
    }
  }
}
*/

/*
static void PseudoScanCallback (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       str;

  if (sfp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (! sfp->pseudo) return;

  str = FindKeyFromFeatDefType (sfp->idx.subtype, FALSE);
  if (StringHasNoText (str)) {
    str = "?";
  }

  fprintf (cfp->logfp, "%s Pseudo %s\n", cfp->buf, str);
}
*/

/*
static void PubScanCallback (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CitArtPtr     cap;
  CleanFlagPtr  cfp;
  CitJourPtr    cjp;
  ImprintPtr    imp;
  ValNodePtr    vnp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
      if (cap != NULL && cap->from == 1) {
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp != NULL) {
          imp = cjp->imp;
          if (imp != NULL) {
            if (imp->part_sup != NULL && StringHasNoText (imp->part_sup)) {
              fprintf (cfp->logfp, "EMPTY PART_SUP %s\n", cfp->buf);
            }
            if (StringDoesHaveText (imp->issue)) {
              if (StringNICmp (imp->issue, "PT ", 3) == 0) {
                fprintf (cfp->logfp, "IMP_PT %s\n", cfp->buf);
              }
              if (StringDoesHaveText (imp->part_sup)) {
                fprintf (cfp->logfp, "SUPN %s\n", cfp->buf);
              }
              if (StringDoesHaveText (imp->part_supi)) {
                fprintf (cfp->logfp, "SUPI %s\n", cfp->buf);
              }
            }
          }
        }
      }
    }
  }
}
*/

/*
static void CheckForParts (BioseqSetPtr bssp, Pointer userdata)
{
  BoolPtr  hasPartsP;

  if (bssp->_class == BioseqseqSet_class_parts) {
    hasPartsP = (BoolPtr) userdata;
    *hasPartsP = TRUE;
  }
}

static void SegmentedFeatCallback (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr     bsp;
  CleanFlagPtr  cfp;
  Boolean       on_seg = FALSE;
  SeqIdPtr      sip;
  SeqLocPtr     slp;

  if (sfp == NULL || sfp->location == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;

  if (bsp->repr == Seq_repr_seg) {
    slp = SeqLocFindNext (sfp->location, NULL);
    while (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        if (SeqIdIn (sip, bsp->id)) {
          on_seg = TRUE;
        }
      }
      slp = SeqLocFindNext (sfp->location, slp);
    }
    if (on_seg) {
      fprintf (cfp->logfp, "SEGSFP %s\n", cfp->buf);
    }
  }
}
*/

/*
static void TildeScanCallback (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;

  if (sfp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (sfp->comment == NULL) return;
  if (StringChr (sfp->comment, '~') == NULL) return;

  fprintf (cfp->logfp, "%s Tilde '%s'\n", cfp->buf, sfp->comment);
}

static void CitGenCallback (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CitGenPtr     cgp;
  ValNodePtr    vnp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != PUB_Gen) continue;
    cgp = (CitGenPtr) vnp->data.ptrvalue;
    if (cgp == NULL) continue;
    if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
      if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
        continue;
      }
    }
    if (cgp->journal == NULL && StringNICmp (cgp->cit, "unpublished", 11) == 0) {
      if (cgp->date != NULL) {
        fprintf (cfp->logfp, "%s No date Cit-gen\n", cfp->buf);
      } else {
        fprintf (cfp->logfp, "%s Unpub Cit-gen\n", cfp->buf);
      }
    }
  }
}
*/

/*
static Boolean FeatIsFullLength (
  BioseqPtr bsp,
  SeqLocPtr slp
)

{
  SeqIntPtr  sintp;

  if (bsp == NULL || slp == NULL) return FALSE;

  if (slp->choice == SEQLOC_WHOLE) return TRUE;
  if (slp->choice == SEQLOC_INT) {
    sintp = (SeqIntPtr) slp->data.ptrvalue;
    if (sintp == NULL) return FALSE;
    if (sintp->from == 0 && sintp->to == bsp->length - 1) return TRUE;
  }

  return FALSE;
}

static void FullMatProtCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  ProtRefPtr    prp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type != 1) continue;
    for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
      if (sfp->data.choice != SEQFEAT_PROT) continue;
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      if (prp == NULL) continue;
      if (prp->processed == 0 && FeatIsFullLength (bsp, sfp->location)) return;
    }
  }

  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type != 1) continue;
    for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
      if (sfp->data.choice != SEQFEAT_PROT) continue;
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      if (prp == NULL) continue;
      if (prp->processed == 0) continue;
      if (! FeatIsFullLength (bsp, sfp->location)) continue;
      fprintf (cfp->logfp, "%s FULL_LENGTH processed %d\n", cfp->buf, (int) prp->processed);
    }
  }
}
*/

/*
static void OctoberReleaseCallback (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  GBQualPtr     gbq;
  ImpFeatPtr    ifp;
  CharPtr       key = NULL;

  if (sfp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (sfp->data.choice == SEQFEAT_IMP) {
    ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (ifp != NULL) {
      key = ifp->key;
      if (key != NULL && StringICmp (key, "conflict") == 0) {
        fprintf (cfp->logfp, "CFT %s\n", cfp->buf);
      }
    }
  }

  if (StringDoesHaveText (sfp->except_text)) {
    if (StringStr (sfp->except_text, "artificial location") != NULL) {
      if (StringICmp (sfp->except_text, "artificial location") == 0) {
        fprintf (cfp->logfp, "ALC %s\n", cfp->buf);
      } else {
        fprintf (cfp->logfp, "EXC %s\n", cfp->buf);
      }
    }
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringCmp (gbq->qual, "label") == 0) {
      fprintf (cfp->logfp, "LBL %s\n", cfp->buf);
    }
    if (StringCmp (gbq->qual, "mobile_element") == 0) {
      if (key != NULL && StringICmp (key, "repeat_region") == 0) {
        fprintf (cfp->logfp, "MOB %s\n", cfp->buf);
      } else {
        fprintf (cfp->logfp, "BOM %s\n", cfp->buf);
      }
    }
  }
}
*/

/*
static void SplitAnticodonCallback (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  SeqLocPtr     anticodon;
  CleanFlagPtr  cfp;
  RnaRefPtr     rrp;
  tRNAPtr       trna;

  if (sfp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->idx.subtype != FEATDEF_tRNA) return;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return;
  if (rrp->type != 3 || rrp->ext.choice != 2) return;
  trna = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trna == NULL) return;

  anticodon = trna->anticodon;
  if (anticodon == NULL) return;
  if (anticodon->choice != SEQLOC_INT) {
    fprintf (cfp->logfp, "ANTICODON %s\n", cfp->buf);
  }
}
*/

/*
static Boolean GeneSpansOrigin (SeqMgrFeatContextPtr context, Int4 bsplength)

{
  Int4Ptr  ivals;

  if (context == NULL || bsplength < 1) return FALSE;
  ivals = context->ivals;
  if (ivals == NULL || context->numivals != 2) return FALSE;
  if (context->strand == Seq_strand_minus) {
    if (ivals [1] == 0 && ivals [2] == bsplength - 1) return TRUE;
  } else {
    if (ivals [2] == 0 && ivals [1] == bsplength - 1) return TRUE;
  }
  return FALSE;
}

static Boolean LIBCALLBACK SplitGeneCallback (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  BioseqPtr     bsp;
  CleanFlagPtr  cfp;
  SeqLocPtr     slp;

  if (sfp == NULL || context == NULL) return TRUE;
  cfp = (CleanFlagPtr) context->userdata;
  if (cfp->logfp == NULL) return TRUE;

  if (sfp->data.choice != SEQFEAT_GENE) return TRUE;
  if (context->numivals < 2) return TRUE;
  bsp = context->bsp;
  if (bsp == NULL) return TRUE;

  if (sfp->excpt) {
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) return TRUE;
  }

  if (bsp->topology == 2) {
    if (context->numivals == 2 && GeneSpansOrigin (context, bsp->length)) return TRUE;
  }

  slp = sfp->location;
  if (slp == NULL) return TRUE;

  if (slp->choice != SEQLOC_INT) {
    fprintf (cfp->logfp, "SPLIT_GENE %s\n", cfp->buf);
  }

  return TRUE;
}

static void SplitGeneBspCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  SeqMgrExploreFeatures (bsp, userdata, SplitGeneCallback, NULL, NULL, NULL);
}
*/

/*
static void MarkComments (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;

  if (sdp == NULL || sdp->choice != Seq_descr_comment) return;
  if (sdp->extended == 0) return;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}

static void RemoveSfpComment (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL || sfp->comment == NULL) return;
  sfp->comment = MemFree (sfp->comment);
}
*/

/*
static void ReportSeqLocMixMix (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  SeqLocPtr     loc;
  SeqLocPtr     slp;

  if (sfp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;
  slp = sfp->location;
  if (slp == NULL || slp->choice != SEQLOC_MIX) return;
  for (loc = (SeqLocPtr) slp->data.ptrvalue; loc != NULL; loc = loc->next) {
    if (loc->choice == SEQLOC_MIX) {
      fprintf (cfp->logfp, "NESTED_SEQLOC_MIX %s\n", cfp->buf);
    }
  }
}
*/

/*
static void ReportSeqLocMixGene (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  SeqLocPtr     slp;

  if (sfp == NULL || userdata == NULL) return;
  if (sfp->data.choice != SEQFEAT_GENE) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;
  slp = sfp->location;
  if (slp == NULL) return;
  if (slp->choice == SEQLOC_MIX) {
    fprintf (cfp->logfp, "MULTI_GENE %s\n", cfp->buf);
  }
}
*/

/*
static void PlasmidScanCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioSourcePtr       biop = NULL;
  CleanFlagPtr       cfp;
  SeqMgrDescContext  dcontext;
  CharPtr            genome = NULL;
  Boolean            okay = FALSE;
  SeqDescrPtr        sdp;
  SubSourcePtr       ssp;

  if (bsp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }
  if (biop == NULL) return;
  if (biop->genome == GENOME_genomic) {
    genome = "GENOME";
  } else if (biop->genome == GENOME_unknown) {
    genome = "ABSENT";
  } else if (biop->genome == GENOME_mitochondrion) {
    genome = "MITOCH";
  } else {
    return;
  }

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_plasmid_name) {
      okay = TRUE;
    }
  }
  if (! okay) return;

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s %s\n", genome, cfp->buf);
    fflush (cfp->logfp);
  }
}
*/

/*
static void SegmentedBioseqCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CleanFlagPtr       cfp;
  SeqMgrDescContext  dcontext;
  SeqDescrPtr        sdp;
  CharPtr            ttl = NULL;

  if (bsp == NULL || userdata == NULL) return;
  if (bsp->repr != Seq_repr_seg) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
  if (sdp != NULL) {
    ttl = (CharPtr) sdp->data.ptrvalue;
  }

  if (ttl != NULL) {
    fprintf (cfp->logfp, "SEG_TITL %s\n", cfp->buf);
    fflush (cfp->logfp);
  } else {
    fprintf (cfp->logfp, "SEG_ANON %s\n", cfp->buf);
    fflush (cfp->logfp);
  }
}
*/

static void RescuePartsTitle (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr     str;
  ValNodePtr  vnp = NULL;

  if (bsp == NULL) return;

  for (vnp = bsp->descr; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Seq_descr_title) break;
  }

  if (vnp != NULL) return;

  str = (CharPtr) userdata;
  SeqDescrAddPointer (&(bsp->descr), Seq_descr_title, StringSave (str));
}

static void NormalizePartsSetTitle (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_parts) return;

  VisitBioseqsInSet (bssp, userdata, RescuePartsTitle);
}

static void NormalizeSegSetTitle (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  BioseqPtr      bsp;
  ObjValNodePtr  ovp;
  SeqEntryPtr    sep;
  CharPtr        str;
  ValNodePtr     vnp1 = NULL, vnp2 = NULL;

  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_segset) return;
  sep = bssp->seq_set;
  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;

  for (vnp1 = bssp->descr; vnp1 != NULL; vnp1 = vnp1->next) {
    if (vnp1->choice == Seq_descr_title) break;
  }
  for (vnp2 = bsp->descr; vnp2 != NULL; vnp2 = vnp2->next) {
    if (vnp2->choice == Seq_descr_title) break;
  }

  if (vnp1 == NULL) return;
  str = (CharPtr) vnp1->data.ptrvalue;

  if (StringDoesHaveText (str)) {
    VisitSetsInSep (sep, (Pointer) str, NormalizePartsSetTitle);
  }

  if (vnp2 == NULL) {
    str = (CharPtr) vnp1->data.ptrvalue;
    if (StringDoesHaveText (str)) {
      SeqDescrAddPointer (&bsp->descr, Seq_descr_title, StringSave (str));
    }
    if (vnp1->extended != 0) {
      ovp = (ObjValNodePtr) vnp1;
      ovp->idx.deleteme = TRUE;
    }
    return;
  }
  if (vnp1->extended != 0) {
    ovp = (ObjValNodePtr) vnp1;
    ovp->idx.deleteme = TRUE;
  }
}

static void CountParts (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  Int4          count;

  if (bssp == NULL || userdata == NULL) return;
  if (bssp->_class != BioseqseqSet_class_parts) return;
  cfp = (CleanFlagPtr) userdata;

  count = VisitBioseqsInSet (bssp, NULL, NULL);
  fprintf (cfp->logfp, "%ld %s", (long) count, cfp->buf);
}

//LCOV_EXCL_START
// Only for SegSets
static void RemoveSegSeqTitle (
  BioseqPtr bsp,
  Pointer userdata
)

{
  ObjValNodePtr  ovp;
  ValNodePtr     vnp = NULL;

  if (bsp == NULL) return;
  if (bsp->repr != Seq_repr_seg) return;

  for (vnp = bsp->descr; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Seq_descr_title) break;
  }

  if (vnp == NULL) return;
  if (vnp->extended != 0) {
    ovp = (ObjValNodePtr) vnp;
    ovp->idx.deleteme = TRUE;
  }
}

static void RemoveSegSetTitle (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CharPtr PNTR   saveP;
  ObjValNodePtr  ovp;
  ValNodePtr     vnp = NULL;

  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_segset) return;

  for (vnp = bssp->descr; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Seq_descr_title) break;
  }

  if (vnp == NULL) return;
  if (vnp->extended != 0) {
    ovp = (ObjValNodePtr) vnp;
    ovp->idx.deleteme = TRUE;
    if (userdata != NULL) {
      saveP = (CharPtr PNTR) userdata;
      *saveP = StringSave ((CharPtr) vnp->data.ptrvalue);
    }
  }
}

static void RemoveAllTitles (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ObjValNodePtr  ovn;
  SeqDescrPtr    sdp = NULL;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
  } else return;

  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_title && sdp->extended != 0) {
      ovn = (ObjValNodePtr) sdp;
      ovn->idx.deleteme = TRUE;
    }
    sdp = sdp->next;
  }
}

static void SaveAllTitles (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ValNodePtr PNTR titleP;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  titleP = (ValNodePtr PNTR) userdata;
  if (titleP == NULL) return;

  ValNodeCopyStr (titleP, Seq_descr_title, (CharPtr) sdp->data.ptrvalue);
}

//LCOV_EXCL_STOP

static void MarkPubFeats (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_PUB) return;
  sfp->idx.deleteme = TRUE;
}

static void MarkPubDescs (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;

  if (sdp == NULL || sdp->choice != Seq_descr_pub) return;
  if (sdp->extended == 0) return;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}

/*
static void RemovePubSerial (
  PubdescPtr pdp,
  Pointer userdata
)

{
  if (pdp == NULL) return;

  CleanUpPubdescBody (pdp, TRUE);
}
*/

static void RemovePubFigNumName (
  PubdescPtr pdp,
  Pointer userdata
)

{
  if (pdp == NULL) return;

  pdp->name = MemFree (pdp->name);
  pdp->fig = MemFree (pdp->fig);
  pdp->maploc = MemFree (pdp->maploc);
  pdp->num = NumberingFree (pdp->num);
}

static void RemovePubRemark (
  PubdescPtr pdp,
  Pointer userdata
)

{
  if (pdp == NULL) return;

  pdp->comment = MemFree (pdp->comment);
}

static void MergePubFeats (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqMgrFeatContext  context1, context2;
  PubdescPtr         pdp1, pdp2;
  SeqFeatPtr         sfp1, sfp2;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  for (sfp1 = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PUB, 0, &context1);
       sfp1 != NULL;
       sfp1 = SeqMgrGetNextFeature (bsp, sfp1, SEQFEAT_PUB, 0, &context1)) {
   if (sfp1->data.choice != SEQFEAT_PUB) continue;
   pdp1 = (PubdescPtr) sfp1->data.value.ptrvalue;
    if (pdp1 == NULL) continue;
    for (sfp2 = SeqMgrGetDesiredFeature (0, bsp, 0, context1.index + 1, NULL, &context2);
         sfp2 != NULL;
         sfp2 = SeqMgrGetNextFeature (bsp, sfp2, SEQFEAT_PUB, 0, &context2)) {
      if (sfp2->data.choice != SEQFEAT_PUB) continue;
      pdp2 = (PubdescPtr) sfp2->data.value.ptrvalue;
      if (pdp2 == NULL) continue;
      if (IsIdenticalPublication (pdp1, pdp2)) {
        sfp1->location = SeqLocMergeEx (bsp, sfp1->location, sfp2->location, FALSE, TRUE, TRUE, FALSE);
        sfp2->idx.deleteme = TRUE;
      }
    }
  }
}

static void PubDescTrumpsFeat (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  PubdescPtr         pdp1, pdp2;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;

  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PUB, 0, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PUB, 0, &fcontext)) {
    pdp1 = (PubdescPtr) sfp->data.value.ptrvalue;
    if (pdp1 == NULL) continue;

    for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
         sdp != NULL;
         sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext)) {
      pdp2 = (PubdescPtr) sdp->data.ptrvalue;
      if (pdp2 == NULL) continue;
      if (IsIdenticalPublication (pdp1, pdp2)) {
        sfp->idx.deleteme = TRUE;
      }
    }
  }
}

/*
static void MarkExonsAndIntrons (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL) return;

  if (sfp->idx.subtype == FEATDEF_exon || sfp->idx.subtype == FEATDEF_intron) {
    sfp->idx.deleteme = TRUE;
  }
}
*/

static void PrintBioseqTitle (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       ttl;

  if (bsp == NULL) return;
  cfp = (CleanFlagPtr) userdata;

  ttl = BioseqGetTitle (bsp);
  if (ttl == NULL) return;

  fprintf (cfp->logfp, "PART %s\n", ttl);
}

static void PrintPartsTitle (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  if (bssp == NULL) return;

  if (bssp->_class != BioseqseqSet_class_parts) return;
  VisitBioseqsInSet (bssp, userdata, PrintBioseqTitle);

}

//LCOV_EXCL_START
// Only for SegSets
static void PrintSegSeqTitleBefore (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       ttl;

  if (bsp == NULL) return;
  if (bsp->repr != Seq_repr_seg) return;

  cfp = (CleanFlagPtr) userdata;

  ttl = BioseqGetTitle (bsp);
  if (ttl == NULL) return;

  fprintf (cfp->logfp, "ORIG %s\n", ttl);
}

// Only for SegSets
static void PrintSegSeqTitleAfter (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       ttl;

  if (bsp == NULL) return;
  if (bsp->repr != Seq_repr_seg) return;

  cfp = (CleanFlagPtr) userdata;

  ttl = BioseqGetTitle (bsp);
  if (ttl == NULL) return;

  fprintf (cfp->logfp, "MSTR %s\n", ttl);
}
//LCOV_EXCL_STOP

static void FindOneCds (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  SeqFeatPtr PNTR sfpp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION || userdata == NULL) return;
  sfpp = (SeqFeatPtr PNTR) userdata;
  *sfpp = sfp;
}

/*
typedef struct partstitles {
  CharPtr  title;
  Boolean  different;
} PartsTitles, PNTR PartsTitlesPtr;

static void PartsTitlesCompare (
  BioseqPtr bsp,
  Pointer userdata
)

{
  PartsTitlesPtr  ptp;
  CharPtr         ttl;

  if (bsp == NULL || userdata == NULL) return;
  ptp = (PartsTitlesPtr) userdata;

  ttl = BioseqGetTitle (bsp);
  if (ttl == NULL) return;

  if (ptp->title == NULL) {
    ptp->title = ttl;
    return;
  }

  if (StringCmp (ttl, ptp->title) != 0) {
    ptp->different = TRUE;
  }
}

static void PartsTitlesFindPartsSet (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_parts) return;

  VisitBioseqsInSet (bssp, userdata, PartsTitlesCompare);
}

static Boolean PartsTitlesIdentical (
  SeqEntryPtr sep
)

{
  PartsTitles  pt;

  if (sep == NULL) return FALSE;

  pt.title = NULL;
  pt.different = FALSE;

  VisitSetsInSep (sep, (Pointer) &pt, PartsTitlesFindPartsSet);

  if (pt.different) return FALSE;

  return TRUE;
}
*/

/*
typedef struct fgbdata {
  CleanFlagPtr  cfp;
  GBBlockPtr    gbp;
  CharPtr       source;
  CharPtr       origin;
  Boolean       mixedsources;
  Boolean       mixedorigins;
} FGBData, PNTR FBGPtr;

static void FindGenBankDiffs (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  FBGPtr        fbgp;
  CleanFlagPtr  cfp;
  GBBlockPtr    gbp;

  if (sdp == NULL || sdp->choice != Seq_descr_genbank) return;
  fbgp = (FBGPtr) userdata;
  if (fbgp == NULL) return;
  cfp = (CleanFlagPtr) fbgp->cfp;
  if (cfp == NULL || cfp->logfp == NULL) return;

  gbp = (GBBlockPtr) sdp->data.ptrvalue;
  if (gbp == NULL) return;

  if (fbgp->gbp == NULL) {
    fbgp->gbp = gbp;
    fbgp->source = gbp->source;
    fbgp->origin = gbp->origin;
    return;
  }

  if (StringDoesHaveText (gbp->source) && StringICmp (gbp->source, fbgp->source) != 0) {
    fprintf (cfp->logfp, "SOURCE  %s  '%s' vs. '%s'\n", cfp->buf, gbp->source, fbgp->source);
  }
  if (StringDoesHaveText (gbp->origin) && StringICmp (gbp->origin, fbgp->origin) != 0) {
    fprintf (cfp->logfp, "ORIGIN  %s  '%s'vs. '%s'\n", cfp->buf, gbp->origin, fbgp->origin);
  }
}

static void FindGenBankBlocks (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  GBBlockPtr    gbp;

  if (sdp == NULL || sdp->choice != Seq_descr_genbank) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  gbp = (GBBlockPtr) sdp->data.ptrvalue;
  if (gbp == NULL) return;

  if (StringDoesHaveText (gbp->source)) {
    fprintf (cfp->logfp, "SOURCE  %s  '%s'\n", cfp->buf, gbp->source);
  }
  if (StringDoesHaveText (gbp->origin)) {
    fprintf (cfp->logfp, "ORIGIN  %s  '%s'\n", cfp->buf, gbp->origin);
  }
}
*/

/*
typedef struct udxprotdata {
  SeqIdPtr    bspid;
  Int4        longest;
  Uint1       processed;
  SeqFeatPtr  sfp;
} UdxProtData, PNTR UdxProtPtr;

static void my_GetLongestProtFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Int4        len;
  ProtRefPtr  prp;
  SeqIdPtr    sip;
  UdxProtPtr  upp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;

  upp = (UdxProtPtr) userdata;
  if (upp == NULL) return;

  sip = SeqLocId (sfp->location);
  if (sip == NULL) return;

  if (! SeqIdIn (sip, upp->bspid)) return;
  len = SeqLocLen (sfp->location);
  if (len == -1) return;

  if (len > upp->longest) {
    upp->sfp = sfp;
    upp->longest = len;
    upp->processed = prp->processed;
  } else if (len == upp->longest) {
    if (prp->processed < upp->processed) {
      upp->sfp = sfp;
      upp->longest = len;
      upp->processed = prp->processed;
    }
  }
}

static SeqFeatPtr my_GetLongestProteinUnindexed (
  BioseqPtr bsp
)

{
  BioseqSetPtr  bssp = NULL;
  UdxProtData   ufd;

  if (bsp == NULL) return NULL;

  MemSet ((Pointer) &ufd, 0, sizeof (UdxProtData));
  ufd.bspid = bsp->id;
  ufd.longest = 0;
  ufd.sfp = NULL;

  VisitFeaturesOnBsp (bsp, (Pointer) &ufd, my_GetLongestProtFeat);

  if (ufd.sfp != NULL && ufd.longest == bsp->length) return ufd.sfp;

  if (bsp->idx.parenttype == OBJ_BIOSEQSET) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
  }

  if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) {
    VisitFeaturesOnSet (bssp, (Pointer) &ufd, my_GetLongestProtFeat);

    if (bssp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) bssp->idx.parentptr;
    }
  }

  if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset) {
    VisitFeaturesOnSet (bssp, (Pointer) &ufd, my_GetLongestProtFeat);
  }

  return ufd.sfp;
}

static void InconsistentPartialsCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqFeatPtr   cds;
  CleanFlagPtr cfp;
  MolInfoPtr   mip;
  SeqDescrPtr  sdp;
  SeqFeatPtr   sfp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  sdp = BioseqGetSeqDescr (bsp, Seq_descr_molinfo, NULL);
  if (sdp == NULL || sdp->choice != Seq_descr_molinfo) {
    fprintf (cfp->logfp, "MISS  %s\n", cfp->buf);
    return;
  }
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;

  sfp = my_GetLongestProteinUnindexed (bsp);
  if (sfp == NULL) return;

  cds = SeqMgrGetCDSgivenProduct (bsp, NULL);

  if (sfp->partial) {
    if (mip->completeness < 2 || mip->completeness > 5) {
      if (cds == NULL || cds->partial) {
        fprintf (cfp->logfp, "BMOL  %s\n", cfp->buf);
      } else {
        fprintf (cfp->logfp, "BPRT  %s\n", cfp->buf);
      }
    }
  } else {
    if (mip->completeness > 1 && mip->completeness < 6) {
      if (cds == NULL || ! cds->partial) {
        fprintf (cfp->logfp, "BCDS  %s\n", cfp->buf);
      } else {
        fprintf (cfp->logfp, "ENDS  %s\n", cfp->buf);
      }
    }
  }
}
*/

/*
static void DoSegSetFix (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  BioseqSetPtr  bssp;
  Uint2         entityID;
  Int4          featcounts [FEATDEF_MAX];
  SeqFeatPtr    cds = NULL;

  if (sep == NULL || cfp == NULL || cfp->logfp == NULL) return;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) return;
    if (bssp->_class == BioseqseqSet_class_wgs_set) return;
  }

  RemoveNucProtSetTitles (sep);
  entityID = ObjMgrGetEntityIDForChoice (sep);
  AssignIDsInEntity (entityID, 0, NULL);
  VisitSetsInSep (sep, NULL, NormalizeSegSetTitle);
  DeleteMarkedObjects (entityID, 0, NULL);
  NormalizeDescriptorOrder (sep);

  MemSet ((Pointer) featcounts, 0, sizeof (featcounts));
  VisitFeaturesInSep (sep, (Pointer) featcounts, CountSegSetFeatures);

  if (featcounts [FEATDEF_GENE] > 1) return;
  if (featcounts [FEATDEF_CDS] > 1) return;

  VisitFeaturesInSep (sep, (Pointer) &cds, FindOneCds);
  if (cds == NULL) return;

  if (! cds->partial) return;

  if (cds->partial) {
    fprintf (cfp->logfp, "PARTIAL  %s", cfp->buf);
  } else {
    fprintf (cfp->logfp, "COMPLETE %s", cfp->buf);
  }

  if (featcounts [FEATDEF_GENE] < 2) {
    if (featcounts [FEATDEF_CDS] < 2) {

      if (featcounts [FEATDEF_GENE] == 0) {
        fprintf (cfp->logfp, ", 0 Genes");
      } else if (featcounts [FEATDEF_GENE] == 1) {
        fprintf (cfp->logfp, ", 1 Gene");
      }

      if (featcounts [FEATDEF_CDS] == 0) {
        fprintf (cfp->logfp, ", 0 CDSs");
      } else if (featcounts [FEATDEF_CDS] == 1) {
        fprintf (cfp->logfp, ", 1 CDS");
      }

      if (featcounts [FEATDEF_exon] > 0) {
        fprintf (cfp->logfp, ", %ld exons", (long) featcounts [FEATDEF_exon]);
      }

      if (featcounts [FEATDEF_intron] > 0) {
        fprintf (cfp->logfp, ", %ld introns", (long) featcounts [FEATDEF_intron]);
      }
    }
  }

  fprintf (cfp->logfp, "\n");

  VisitSetsInSep (sep, (Pointer) cfp, PrintPartsTitle);
  VisitBioseqsInSep (sep, (Pointer) cfp, PrintSegSeqTitleBefore);

  VisitBioseqsInSep (sep, NULL, RemoveSegSeqTitle);
  VisitSetsInSep (sep, NULL, RemoveSegSetTitle);

  DeleteMarkedObjects (entityID, 0, NULL);
  SeqMgrIndexFeatures (entityID, 0);
  DoAutoDef (sep, entityID, TRUE, FALSE);
  NormalizeDescriptorOrder (sep);
  VisitBioseqsInSep (sep, (Pointer) cfp, PrintSegSeqTitleAfter);

  fprintf (cfp->logfp, "\n\n\n");
  fflush (cfp->logfp);
}
*/

/*
static Boolean XmlUnicodeEndOkay (
  Char ch
)

{
  if (ch == ';' || ch == ' ' || ch == '\0') return TRUE;

  return FALSE;
}

static void CheckForXml (
  CharPtr str,
  CleanFlagPtr cfp
)

{
  Char     buf [64];
  Char     ch;
  Int2     i;
  CharPtr  ptr;

  if (StringHasNoText (str)) return;
  if (cfp == NULL) return;

  ch = *str;
  while (ch != '\0') {
    if (ch == '&') {
      if (StringNICmp (str, "&amp", 4) == 0) {
        if (XmlUnicodeEndOkay (str [4])) {
          fprintf (cfp->logfp, "AMP %s\n", cfp->buf);
        }
      } else if (StringNICmp (str, "&apos", 5) == 0) {
        if (XmlUnicodeEndOkay (str [5])) {
          fprintf (cfp->logfp, "POS %s\n", cfp->buf);
        }
      } else if (StringNICmp (str, "&gt", 3) == 0) {
        if (XmlUnicodeEndOkay (str [3])) {
          fprintf (cfp->logfp, "GTR %s\n", cfp->buf);
        }
      } else if (StringNICmp (str, "&lt", 3) == 0) {
        if (XmlUnicodeEndOkay (str [3])) {
         fprintf (cfp->logfp, "LSS %s\n", cfp->buf);
        }
       } else if (StringNICmp (str, "&quot", 5) == 0) {
        if (XmlUnicodeEndOkay (str [5])) {
          fprintf (cfp->logfp, "QUT  %s\n", cfp->buf);
        }
      } else if (StringNICmp (str, "&#", 2) == 0) {
        ptr = str + 2;
        ch = *ptr;
        i = 0;
        while (IS_DIGIT (ch) && i < 60) {
          buf [i] = ch;
          i++;
          ptr++;
          ch = *ptr;
        }
        buf [i] = '\0';
        if (XmlUnicodeEndOkay (ch)) {
          fprintf (cfp->logfp, "NUM %s %s\n", buf, cfp->buf);
        } else {
          fprintf (cfp->logfp, "BAD %s %s\n", buf, cfp->buf);
        }
      }
    }
    str++;
    ch = *str;
  }
}

static void FindXmlInFeatComments (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  GeneRefPtr    grp;
  ProtRefPtr    prp;
  RnaRefPtr     rrp;
  CharPtr       str;
  ValNodePtr    vnp;

  if (sfp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  CheckForXml (sfp->comment, cfp);
  CheckForXml (sfp->title, cfp);

  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      CheckForXml (grp->locus, cfp);
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        CheckForXml (str, cfp);
      }
      break;
    case SEQFEAT_PROT :
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        CheckForXml (str, cfp);
      }
      CheckForXml (prp->desc, cfp);
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->ext.choice == 1) {
        CheckForXml ((CharPtr) rrp->ext.value.ptrvalue, cfp);
      }
      break;
    default :
      break;
  }
}

static void FindXmlInCommentDescs (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       str;

  if (sdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_title :
      str = (CharPtr) sdp->data.ptrvalue;
      CheckForXml (str, cfp);
      break;
    case Seq_descr_comment :
      str = (CharPtr) sdp->data.ptrvalue;
      CheckForXml (str, cfp);
      break;
    default :
      break;
  }
}
*/

/*
static void FindCDS (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr     bsp;
  CleanFlagPtr  cfp;
  CharPtr       dir;
  MolInfoPtr    mip;
  SeqDescrPtr   sdp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp == NULL) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;

  if (mip->biomol == MOLECULE_TYPE_RRNA ||
      mip->biomol == MOLECULE_TYPE_TRNA ||
      mip->biomol == MOLECULE_TYPE_SNRNA ||
      mip->biomol == MOLECULE_TYPE_SCRNA ||
      mip->biomol == MOLECULE_TYPE_SNORNA ||
      mip->biomol == MOLECULE_TYPE_NCRNA) {
    dir = StringRChr (cfp->filename, DIRDELIMCHR);
    if (dir != NULL) {
      dir++;
    }
    fprintf (cfp->logfp, "RNA %s %s\n", dir, cfp->buf);
    fflush (cfp->logfp);
  } else if (mip->biomol == MOLECULE_TYPE_MRNA && SeqLocStrand (sfp->location) == Seq_strand_minus) {
    dir = StringRChr (cfp->filename, DIRDELIMCHR);
    if (dir != NULL) {
      dir++;
    }
    fprintf (cfp->logfp, "NEG %s %s\n", dir, cfp->buf);
    fflush (cfp->logfp);
  }
}
*/

/*
static void FindCodon (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  GBQualPtr     gbq;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringCmp (gbq->qual, "codon") != 0) continue;
    if (StringHasNoText (gbq->val)) continue;
    fprintf (cfp->logfp, "codon %s\n", cfp->buf);
    fflush (cfp->logfp);
  }
}
*/

/*
static void FindArrow (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       str;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  str = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (str)) return;
  if (StringChr (str, '>') != NULL) {
    fprintf (cfp->logfp, "ARROW %s\n", cfp->buf);
  }
}
*/

/*
static void FindHGNC (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr         cfp;
  DbtagPtr             dbt;
  GeneNomenclaturePtr  gnp;
  GeneRefPtr           grp;
  ValNodePtr           vnp;

  if (sfp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL) continue;
    if (StringCmp (dbt->db, "HGNC") == 0) {
      fprintf (cfp->logfp, "HGNC %s\n", cfp->buf);
    } else if (StringCmp (dbt->db, "MGI") == 0) {
      fprintf (cfp->logfp, "MGI %s\n", cfp->buf);
    }
  }

  if (sfp->data.choice == SEQFEAT_GENE) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL) {
      gnp = grp->formal_name;
      if (gnp != NULL) {
        dbt = gnp->source;
        if (dbt != NULL) {
          if (StringCmp (dbt->db, "HGNC") == 0) {
            fprintf (cfp->logfp, "HGNC %s\n", cfp->buf);
          } else if (StringCmp (dbt->db, "MGI") == 0) {
            fprintf (cfp->logfp, "MGI %s\n", cfp->buf);
          }
        }
      }
    }
  }
}

static CharPtr lim_types [] = {
  "un",
  "gt",
  "lt",
  "tr",
  "tl",
  "cr",
  "ot",
  "bd"
};

static void FindPointFuzz (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  IntFuzzPtr    ifp;
  Int4          lim;
  SeqLocPtr     slp;
  SeqPntPtr     spp;

  if (sfp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (slp = SeqLocFindNext (sfp->location, NULL);
       slp != NULL;
       slp = SeqLocFindNext (sfp->location, slp)) {
    if (slp->choice != SEQLOC_PNT) continue;
    spp = (SeqPntPtr) slp->data.ptrvalue;
    if (spp == NULL) continue;
    ifp = spp->fuzz;
    if (ifp == NULL) continue;
    if (ifp->choice != 4) continue;
    lim = ifp->a;
    if (lim == 255) {
      lim = 6;
    } else if (lim < 0 || lim > 5) {
      lim = 7;
    }
    fprintf (cfp->logfp, "LIM %s %s\n", lim_types [(int) lim], cfp->buf);
  }
}
*/

static CharPtr stopWords [] = {
  "a",
  "about",
  "again",
  "all",
  "almost",
  "also",
  "although",
  "always",
  "among",
  "an",
  "and",
  "another",
  "any",
  "are",
  "as",
  "at",
  "be",
  "because",
  "been",
  "before",
  "being",
  "between",
  "both",
  "but",
  "by",
  "can",
  "could",
  "did",
  "do",
  "does",
  "done",
  "due",
  "during",
  "each",
  "either",
  "enough",
  "especially",
  "etc",
  "for",
  "found",
  "from",
  "further",
  "had",
  "has",
  "have",
  "having",
  "here",
  "how",
  "however",
  "i",
  "if",
  "in",
  "into",
  "is",
  "it",
  "its",
  "itself",
  "just",
  "kg",
  "km",
  "made",
  "mainly",
  "make",
  "may",
  "mg",
  "might",
  "ml",
  "mm",
  "most",
  "mostly",
  "must",
  "nearly",
  "neither",
  "no",
  "nor",
  "obtained",
  "of",
  "often",
  "on",
  "our",
  "overall",
  "perhaps",
  "pmid",
  "quite",
  "rather",
  "really",
  "regarding",
  "seem",
  "seen",
  "several",
  "should",
  "show",
  "showed",
  "shown",
  "shows",
  "significantly",
  "since",
  "so",
  "some",
  "such",
  "than",
  "that",
  "the",
  "their",
  "theirs",
  "them",
  "then",
  "there",
  "therefore",
  "these",
  "they",
  "this",
  "those",
  "through",
  "thus",
  "to",
  "upon",
  "use",
  "used",
  "using",
  "various",
  "very",
  "was",
  "we",
  "were",
  "what",
  "when",
  "which",
  "while",
  "with",
  "within",
  "without",
  "would",
  NULL
};

static Boolean IsStopWord (
  CharPtr str
)

{
  Int2  i;

  if (StringHasNoText (str)) return FALSE;

  for (i = 0; stopWords [i] != NULL; i++) {
    if (StringICmp (str, stopWords [i]) == 0) return TRUE;
  }

  return FALSE;
}

static ValNodePtr GetTitleWords (
  CharPtr title
)

{
  Char        ch;
  Boolean     goOn = TRUE;
  ValNodePtr  head = NULL;
  ValNodePtr  last = NULL;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;

  if (StringHasNoText (title)) return NULL;

  tmp = StringSave (title);
  if (tmp == NULL) return NULL;

  ptr = tmp;
  ch = *ptr;
  if (ch == '\0') {
    goOn = FALSE;
  }
  while (goOn) {
    while (ch != '\0' && (! IS_ALPHANUM (ch))) {
      ptr++;
      ch = *ptr;
    }
    str = ptr;
    while (ch != '\0' && IS_ALPHANUM (ch)) {
      ptr++;
      ch = *ptr;
    }
    if (ch == '\0') {
      goOn = FALSE;
    }
    *ptr = '\0';
    ptr++;
    ch = *ptr;
    TrimSpacesAroundString (str);
#if 0
    if (! IsStopWord (str)) {
      vnp = ValNodeCopyStr (&last, 0, str);
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;
    }
#endif
    vnp = ValNodeCopyStr (&last, 0, str);
    if (head == NULL) {
      head = vnp;
    }
    last = vnp;
  }

  MemFree (tmp);

  return head;
}

/*
static Boolean TitleEndsWithOrganism (
  CharPtr title,
  CharPtr organism
)

{
  Char     ch;
  size_t   len1, len2;

  if (StringHasNoText (title)) return FALSE;
  if (StringHasNoText (organism)) return FALSE;

  len1 = StringLen (title);
  len2 = StringLen (organism);
  if (len2 + 4 >= len1) return FALSE;

  ch = title [len1 - 1];
  if (ch != ']') return FALSE;
  ch = title [len1 - len2 - 2];
  if (ch != '[') return FALSE;

  return TRUE;
}

static Boolean TitleAndOrganismWordMatch (
  CharPtr title,
  CharPtr organism
)

{
  Boolean     match = FALSE;
  ValNodePtr  organismwords, titlewords;
  CharPtr     orgstr, ttlstr;
  ValNodePtr  vnp1, vnp2;

  if (StringHasNoText (title)) return FALSE;
  if (StringHasNoText (organism)) return FALSE;

  titlewords = GetTitleWords (title);
  organismwords = GetTitleWords (organism);


  for (vnp1 = titlewords; vnp1 != NULL; vnp1 = vnp1->next) {
    ttlstr = (CharPtr) vnp1->data.ptrvalue;
    if (StringHasNoText (ttlstr)) continue;
    if (StringLen (ttlstr) < 5) continue;
    for (vnp2 = organismwords; vnp2 != NULL; vnp2 = vnp2->next) {
      orgstr = (CharPtr) vnp2->data.ptrvalue;
      if (StringHasNoText (orgstr)) continue;
      if (StringLen (orgstr) < 5) continue;
      if (StringICmp (ttlstr, orgstr) != 0) continue;
      match = TRUE;
    }
  }

  ValNodeFreeData (titlewords);
  ValNodeFreeData (organismwords);

  return match;
}

static void CheckProtTitle (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  Char           buf [128];
  CleanFlagPtr   cfp;
  CharPtr        oldname = NULL;
  OrgModPtr      omp;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  SeqDescrPtr    sdp;
  CharPtr        str;
  CharPtr        taxname = NULL;
  CharPtr        title;
  SeqDescrPtr    ttl;

  if (bsp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  if (! ISA_aa (bsp->mol)) return;

  ttl = GetNextDescriptorUnindexed (bsp, Seq_descr_title, NULL);
  if (ttl == NULL || ttl->choice != Seq_descr_title) return;
  str = (CharPtr) ttl->data.ptrvalue;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        taxname = orp->taxname;
        onp = orp->orgname;
        if (onp != NULL) {
          for (omp = onp->mod; omp != NULL; omp = omp->next) {
            if (omp->subtype == ORGMOD_old_name) {
              oldname = omp->subname;
            }
          }
        }
      }
    }
  }

  if (TitleEndsWithOrganism (str, oldname)) return;
  if (TitleEndsWithOrganism (str, taxname)) return;

  if (TitleAndOrganismWordMatch (str, oldname)) return;
  if (TitleAndOrganismWordMatch (str, taxname)) return;

  SeqIdWrite (bsp->id, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);

#if 0
  fprintf (cfp->logfp, "ORG\t%s\t%s\t%s\t%s\n", cfp->buf, buf, str, taxname);
#endif

  title = NewCreateDefLineEx (NULL, bsp, TRUE, FALSE, FALSE);
  if (title == NULL) return;

  fprintf (cfp->logfp, "%s\t%s\n", buf, title);

  MemFree (title);
}
*/

/*
static void FindTerminalInProtTitle (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr      bsp;
  CleanFlagPtr   cfp;
  Boolean        hasGibbsq;
  Boolean        partial5, partial3;
  ProtRefPtr     prp;
  SeqIdPtr       sip;
  CharPtr        str;
  SeqDescrPtr    ttl;

  if (sfp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  if (! sfp->partial) return;

  if (sfp->data.choice != SEQFEAT_PROT) return;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  if (partial5  || partial3) return;

  bsp = BioseqFind (SeqLocId (sfp->location));
  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;

  hasGibbsq = FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GIBBSQ) {
      hasGibbsq = TRUE;
    }
  }
  if (! hasGibbsq) return;

  ttl = GetNextDescriptorUnindexed (bsp, Seq_descr_title, NULL);
  if (ttl == NULL || ttl->choice != Seq_descr_title) return;
  str = (CharPtr) ttl->data.ptrvalue;
  if (StringHasNoText (str)) return;

  if (StringStr (str, "{N-terminal}") == NULL && StringStr (str, "{C-terminal}") == NULL) return;

  fprintf (cfp->logfp, "TRM\t%s\t%s\n", cfp->buf, str);
}
*/

/*
static CharPtr findrepstrs [] = {
  " ,", ",,", NULL
};

static void FindBadPunctuation (Uint2 entityID, Uint4 itemID, Uint2 itemtype, Pointer userdata)

{
  CleanFlagPtr  cfp;

  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  fprintf (cfp->logfp, "COM\t%s\n", cfp->buf);
}
*/


/*
typedef struct dummysmfedata {
  Int4        max;
  Int4        num_at_max;
  Int4        num_trans_spliced;
  Boolean     equivalent_genes;
  GeneRefPtr  grp_at_max;
} DummySmfeData, PNTR DummySmfePtr;

static Boolean LIBCALLBACK DummySMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)


{
  DummySmfePtr  dsp;
  GeneRefPtr    grp, grpx;
  Int4          len;
  Boolean       redundantgenexref = FALSE;
  CharPtr       syn1, syn2;

  if (sfp == NULL || context == NULL) return TRUE;
  dsp = context->userdata;
  if (dsp == NULL) return TRUE;
  if (sfp->data.choice != SEQFEAT_GENE) return TRUE;
  grp = (GeneRefPtr) sfp->data.value.ptrvalue;
  if (grp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < dsp->max) {
    dsp->max = len;
    dsp->num_at_max = 1;
    dsp->num_trans_spliced = 0;
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) {
      (dsp->num_trans_spliced)++;
    }
    dsp->equivalent_genes = FALSE;
    dsp->grp_at_max = grp;
  } else if (len == dsp->max) {
    (dsp->num_at_max)++;
    if (StringISearch (sfp->except_text, "trans-splicing") != NULL) {
      (dsp->num_trans_spliced)++;
    }
    grpx = dsp->grp_at_max;
    if (grpx != NULL) {
      redundantgenexref = FALSE;
      if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
        if (StringICmp (grp->locus_tag, grpx->locus_tag) == 0) {
          redundantgenexref = TRUE;
        }
      } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
        if (StringICmp (grp->locus, grpx->locus) == 0) {
          redundantgenexref = TRUE;
        }
      } else if (grp->syn != NULL && grpx->syn != NULL) {
        syn1 = (CharPtr) grp->syn->data.ptrvalue;
        syn2 = (CharPtr) grpx->syn->data.ptrvalue;
        if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
          if (StringICmp (syn1, syn2) == 0) {
            redundantgenexref = TRUE;
          }
        }
      }
    }
    if (redundantgenexref) {
      dsp->equivalent_genes = TRUE;
    }
  }

  return TRUE;
}

static void FindMissingGeneXref (SeqFeatPtr sfp, Pointer userdata)

{
  CleanFlagPtr   cfp;
  Int2           count;
  DummySmfeData  dsd;
  GeneRefPtr     grp;

  if (sfp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  if (sfp->idx.subtype == FEATDEF_GENE ||
      sfp->idx.subtype == FEATDEF_operon ||
      sfp->idx.subtype == FEATDEF_gap ||
      sfp->idx.subtype == FEATDEF_repeat_region ||
      sfp->idx.subtype == FEATDEF_mobile_element) return;

  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL) return;

  MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
  dsd.max = INT4_MAX;
  dsd.num_at_max = 0;
  dsd.num_trans_spliced = 0;
  dsd.equivalent_genes = FALSE;
  dsd.grp_at_max = NULL;
  count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE, NULL, 0,
                                           LOCATION_SUBSET, (Pointer) &dsd, DummySMFEProc);

  if (dsd.num_at_max <= 1) return;
  if (cfp->isSmallGenome && dsd.num_at_max == dsd.num_trans_spliced) return;

  if (dsd.equivalent_genes) {
    fprintf (cfp->logfp, "EQUV\t%s\n", cfp->buf);
  } else {
    fprintf (cfp->logfp, "XREF\t%s\n", cfp->buf);
  }
}

static void FindUnsyncedCDS (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr     bsp;
  CleanFlagPtr  cfp;
  MolInfoPtr    mip;
  Boolean       partial = FALSE;
  SeqDescrPtr   sdp;
  SeqIdPtr      sip;

  if (sfp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  if (sfp->data.choice != SEQFEAT_CDREGION) return;

  sip = SeqLocId (sfp->product);
  if (sip == NULL) return;
  bsp = BioseqFind (sip);
  if (bsp == NULL || !ISA_aa (bsp->mol)) return;

  sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (sdp == NULL || sdp->choice != Seq_descr_molinfo) return;
  mip = (MolInfoPtr) sdp->data.ptrvalue;
  if (mip == NULL) return;

  if (mip != NULL && mip->completeness > 1 && mip->completeness < 6) {
    partial = TRUE;
  }

  if (sfp->partial == partial) return;

  fprintf (cfp->logfp, "ASEC\t%s\n", cfp->buf);
}

static void FindNorCTerminal (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       str;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  str = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (str)) return;
  if (StringStr (str, "{N-terminal}") != NULL || StringStr (str, "{C-terminal}") != NULL) {
    fprintf (cfp->logfp, "TERM\t%s\n", cfp->buf);
  }
}

static void FindStrandOther (SeqFeatPtr sfp, Pointer userdata)

{
  CleanFlagPtr  cfp;
  SeqLocPtr     slp;

  if (sfp == NULL || sfp->location == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  slp = SeqLocFindNext (sfp->location, NULL);
  while (slp != NULL) {
    if (SeqLocStrand (slp) == Seq_strand_other) {
      fprintf (cfp->logfp, "OTHR\t%s\n", cfp->buf);
    }
    slp = SeqLocFindNext (sfp->location, slp);
  }
}

static Boolean FoundSourceTilde (CharPtr str)

{
  if (StringHasNoText (str)) return FALSE;


  return FALSE;
}

static void FindSourceTildes (BioSourcePtr biop, Pointer userdata)

{
  CleanFlagPtr  cfp;
  OrgModPtr     omp;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  SubSourcePtr  ssp;

  if (biop == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      for (omp = onp->mod; omp != NULL; omp = omp->next) {
        if (omp->subtype == ORGMOD_other && StringDoesHaveText (omp->subname)) {
          if (FoundSourceTilde (omp->subname)) {
            fprintf (cfp->logfp, "TILD\t%s\n", cfp->buf);
          }
        }
      }
    }
  }

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_other && StringDoesHaveText (ssp->name)) {
      if (FoundSourceTilde (ssp->name)) {
        fprintf (cfp->logfp, "TILD\t%s\n", cfp->buf);
      }
    }
  }
}

static CharPtr findtildestrs [] = {
  "were ~25 cm in height (~3 weeks)",
  "generally ~3 weeks",
  "sequencing (~4 96-well plates)",
  "size distribution (~2 kb)",
  "sequencing (~3 96-well plates)",
  "vector. 1~2 ul of ligated",
  "different strains~of mice",
  "oligo-dT-NotI primer~(5'-biotin",
  "sizes of 200~800 bp were purified",
  "Tween 20 (~50 ml per tree)",
  "the SMART approach (~http://www.evrogen.com",
  "the morning (~10 am) with",
  "(host) sequences (~10%)",
  "Plasmid; ~The F03-1270",
  "using STS-PCR~from Eb",
  "specific to~the Eb",
  "side of insert); , M.F., Lennon",
  "Uni-ZAP XR vector. 1~2 ul of",
  "from diploid~Secale montanum",
  "homology with~U43516,",
  "plasmid; ~Assembled EST",
  "databases.~Different cDNA",
  "enzyme PstI.~DH5-alpha",
  "as they~were prepared",
  "loci in~the genome",
  "P{CaSpeR}Cp1~50C (FBti0004219)",
  "seedlings with 2~4 leaves",
  "tween 20 (~50mLs per tree)",
  "nnotated by GenomeRefine~~",
  "based on SOLiD3 (Applied Biosystems)~~",
  "Biological resourse center, NITE (NRBC)~~",
  "developmental01.html~~",
  "http://bionano.toyo.ac.jp/~~",
  "http://dictycdb1.biol.tsukuba.ac.jp/acytodb/~~",
  "http://egg.umh.es~~",
  "http://www.aist.go.jp/~~",
  "http://www.bio.nite.go.jp/~~DOGAN ; Database",
  "http://www.bio.nite.go.jp/ngac/e/~~",
  "http://www.brs.kyushu-u.ac.jp/~fcmic/~~",
  "http://www.miyazaki-u.ac.jp/ir/english/index.html~~",
  "URL:http://www.bio.nite.go.jp/~~",
  "RAST version 2.0 (http://rast.nmpdr.org/)~~",
  "URL:http://www.tmd.ac.jp/grad/bac/database.html~~",
  NULL
};

static void FindTildeProc (Uint2 entityID, Uint4 itemID, Uint2 itemtype, Pointer userdata)

{
  CleanFlagPtr cfp;

  cfp = (CleanFlagPtr) userdata;
  fprintf (cfp->logfp, "TILD\t%s\n", cfp->buf);
}
*/

/*
typedef struct genexrefskew {
  Int4  gene_locus;
  Int4  gene_locus_tag;
  Int4  xref_locus;
  Int4  xref_locus_tag;
} GeneSkew, PNTR GeneSkewPtr;

static void FindGeneXrefMismatch (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  GeneRefPtr      grp;
  GeneSkewPtr     gsp;
  SeqFeatXrefPtr  xref;

  if (sfp == NULL) return;
  gsp = (GeneSkewPtr) userdata;
  if (gsp == NULL) return;

  if (sfp->data.choice == SEQFEAT_GENE) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL) {
      if (StringDoesHaveText (grp->locus)) {
        (gsp->gene_locus)++;
      }
      if (StringDoesHaveText (grp->locus_tag)) {
        (gsp->gene_locus_tag)++;
      }
    }
  }

  for (xref = sfp->xref; xref != NULL; xref = xref->next) {
    if (xref->data.choice != SEQFEAT_GENE) continue;
    grp = (GeneRefPtr) xref->data.value.ptrvalue;
    if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) continue;
      if (StringDoesHaveText (grp->locus)) {
        (gsp->xref_locus)++;
      }
      if (StringDoesHaveText (grp->locus_tag)) {
        (gsp->xref_locus_tag)++;
      }
  }
}

static void FindBadAltitude (BioSourcePtr biop, Pointer userdata)

{
  CleanFlagPtr  cfp;
  size_t        len;
  SubSourcePtr  ssp;
  CharPtr       str;

  if (biop == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype != SUBSRC_altitude) continue;
    str = ssp->name;
    if (StringHasNoText (str)) continue;
    len = StringLen (str);
    if (len < 3) continue;
    if (str [len - 1] == '.') {
      fprintf (cfp->logfp, "PERI\t%s\n", cfp->buf);
    }
    if (StringStr (str, "f") != NULL) {
      fprintf (cfp->logfp, "FEET\t%s\n", cfp->buf);
    }
  }
}
*/

/*
static void AuthStdMlCallback (
  PubdescPtr pdp,
  Pointer userdata
)

{
  AuthListPtr   alp;
  AuthorPtr     ap;
  CleanFlagPtr  cfp;
  ValNodePtr    names;
  PersonIdPtr   pid;
  PubPtr        pub;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (pub = pdp->pub; pub != NULL; pub = pub->next) {
    alp = GetAuthorListForPub (pub);
    if (alp == NULL || alp->choice != 1) continue;

    for (names = alp->names; names != NULL; names = names->next) {
      ap = names->data.ptrvalue;
      if (ap == NULL) continue;
      pid = ap->name;
      if (pid == NULL) continue;
      if (pid->choice == 3) {
        fprintf (cfp->logfp, "MDLN\t%s\n", cfp->buf);
      } else if (pid->choice == 4) {
        fprintf (cfp->logfp, "STNG\t%s\n", cfp->buf);
      }
    }
  }
}
*/

#define IS_INSD 1
#define IS_REFSEQ 2

static Int2 RefSeqOrINSDProc (
  SeqIdPtr sip
)

{
  while (sip != NULL) {
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
        return IS_INSD;
      case SEQID_OTHER :
        return IS_REFSEQ;
        break;
      default:
        break;
    }

    sip = sip->next;
  }

  return 0;
}

static Int2 RefSeqOrINSD (
  SeqIdPtr sip
)

{
  Int4      gi;
  Int2      res = 0;
  SeqIdPtr  sipset;

  res = RefSeqOrINSDProc (sip);
  if (res != 0) {
    return res;
  }

  while (sip != NULL) {
    if (sip->choice == SEQID_GI) {
      gi = (Int4) sip->data.intvalue;
      if (gi > 0) {
        sipset = GetSeqIdSetForGI (gi);
        if (sipset != NULL) {
          res = RefSeqOrINSD (sipset);
          SeqIdSetFree (sipset);
          return res;
        }
      }
    }
    sip = sip->next;
  }

  return res;
}

static void FarFeatScan (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CleanFlagPtr       cfp;
  Int2               child;
  SeqMgrFeatContext  context;
  SeqMgrDescContext  dcontext;
  DeltaSeqPtr        dsp;
  Int4               gi;
  Boolean            speciallocal = FALSE;
  Int2               master;
  ObjectIdPtr        oip;
  SeqDescrPtr        sdp;
  DeltaSeqPtr        sdsp;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  SeqIdPtr           sipset;
  SeqLocPtr          slp;
  BioseqPtr          sub;
  CharPtr            str;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;

  if (bsp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  if (bsp->repr != Seq_repr_delta) return;
  if (ISA_aa (bsp->mol)) return;

  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &context);
  while (sfp != NULL) {
    if (sfp->idx.subtype != FEATDEF_centromere &&
        sfp->idx.subtype != FEATDEF_telomere &&
        sfp->idx.subtype != FEATDEF_REGION &&
        sfp->idx.subtype != FEATDEF_rep_origin) return;
    speciallocal = TRUE;
    sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &context);
  }

  if (speciallocal) {
    fprintf (cfp->logfp, "LOCL\t%s\n", cfp->buf);
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL) {
        if (StringCmp (oip->str, "FeatureFetchPolicy") == 0) {
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL || ufp->data.ptrvalue == NULL) continue;
            if (StringCmp (oip->str, "Policy") == 0) {
              str = (CharPtr) ufp->data.ptrvalue;
              if (StringICmp (str, "OnlyNearFeatures") == 0) {
                fprintf (cfp->logfp, "ONLY\t%s\n", cfp->buf);
              } else if (StringICmp (str, "AllowFarFeatures") == 0) {
              }
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  master = RefSeqOrINSD (bsp->id);

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice != 1) continue;
    slp = (SeqLocPtr) dsp->data.ptrvalue;
    if (slp == NULL || slp->choice == SEQLOC_NULL) continue;
    sip = SeqLocId (slp);
    if (sip == NULL) continue;

    child = 0;
    if (sip->choice == SEQID_GI) {
      gi = (Int4) sip->data.intvalue;
      if (gi > 0) {
        sipset = GetSeqIdSetForGI (gi);
        if (sipset != NULL) {
          child = RefSeqOrINSD (sipset);
          SeqIdSetFree (sipset);
        }
      }
    } else {
      child = RefSeqOrINSD (sip);
    }
    if (child == master) continue;

    sub = BioseqLockById (sip);
    if (sub != NULL) {
      if (VisitFeaturesOnBsp (sub, NULL, NULL) == 0) {
        if (sub->repr == Seq_repr_delta) {
          for (sdsp = (DeltaSeqPtr) sub->seq_ext; sdsp != NULL; sdsp = sdsp->next) {
            if (sdsp->choice != 1) continue;
            fprintf (cfp->logfp, "DLTA\t%s\n", cfp->buf);
            break;
          }
        } else {
          fprintf (cfp->logfp, "TRPL\t%s\n", cfp->buf);
        }
      }
      BioseqUnlock (sub);
    }
  }  
}

static void DoASNScan (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  Uint2 entityID;

  if (sep == NULL || cfp == NULL || cfp->logfp == NULL) return;

  entityID = ObjMgrGetEntityIDForChoice (sep);
  SeqMgrIndexFeatures (entityID, 0);
  VisitBioseqsInSep (sep, (Pointer) cfp, FarFeatScan);

  /*
  VisitPubdescsInSep (sep, (Pointer) cfp, AuthStdMlCallback);
  */

  /*
  GeneSkew  gs;

  MemSet ((Pointer) &gs, 0, sizeof (GeneSkew));
  VisitFeaturesInSep (sep, (Pointer) &gs, FindGeneXrefMismatch);
  if (gs.gene_locus == 0 && gs.gene_locus_tag > 0) {
    if (gs.xref_locus > 0 && gs.xref_locus_tag == 0) {
      fprintf (cfp->logfp, "XREF\t%s\n", cfp->buf);
    }
  } else if (gs.gene_locus > 0 && gs.gene_locus_tag == 0) {
    if (gs.xref_locus == 0 && gs.xref_locus_tag > 0) {
      fprintf (cfp->logfp, "XREF\t%s\n", cfp->buf);
    }
  }
  VisitBioSourcesInSep (sep, (Pointer) cfp, FindBadAltitude);
  */

  /*
  Uint2 entityID;

  entityID = ObjMgrGetEntityIDForChoice (sep);

  FindStringsInEntity (entityID, findtildestrs, FALSE, FALSE, FALSE, UPDATE_NEVER,
                       NULL, NULL, NULL, TRUE, FindTildeProc, (Pointer) cfp);

  SeqMgrIndexFeatures (entityID, NULL);

  VisitFeaturesInSep (sep, (Pointer) cfp, FindUnsyncedCDS);
  VisitFeaturesInSep (sep, (Pointer) cfp, FindMissingGeneXref);
  VisitDescriptorsInSep (sep, (Pointer) cfp, FindNorCTerminal);
  VisitFeaturesInSep (sep, (Pointer) cfp, FindStrandOther);
  VisitBioSourcesInSep (sep, (Pointer) cfp, FindSourceTildes);
  */

  /*
  entityID = ObjMgrGetEntityIDForChoice (sep);
  FindStringsInEntity (entityID, findrepstrs, FALSE, FALSE, FALSE, UPDATE_NEVER,
                       NULL, NULL, NULL, TRUE, FindBadPunctuation, (Pointer) cfp);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, FindTerminalInProtTitle);
  */

  /*
  VisitBioseqsInSep (sep, (Pointer) cfp, CheckProtTitle);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, FindHGNC);
  VisitFeaturesInSep (sep, (Pointer) cfp, FindPointFuzz);
  */

  /*
  VisitDescriptorsInSep (sep, (Pointer) cfp, FindArrow);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, FindCodon);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, FindCDS);
  */

  /*
  ByteStorePtr  bs = NULL, tmp = NULL;

  if (sep == NULL || cfp == NULL || cfp->logfp == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) cfp, FindXmlInFeatComments);
  VisitDescriptorsInSep (sep, (Pointer) cfp, FindXmlInCommentDescs);

  bs = Se2Bs (sep);

  UpdateReplacedECNumbers (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    fprintf (cfp->logfp, "REP %s\n", cfp->buf);
    fflush (cfp->logfp);
  }
  BSFree (bs);
  bs = tmp;

  DeleteBadECNumbers (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    fprintf (cfp->logfp, "DEL %s\n", cfp->buf);
    fflush (cfp->logfp);
  }
  BSFree (bs);
  bs = tmp;

  BSFree (bs);
  */

  /*
  ByteStorePtr  bs = NULL, tmp = NULL;

  if (sep == NULL || cfp == NULL || cfp->logfp == NULL) return;
  */

  /*
  VisitBioseqsInSep (sep, (Pointer) cfp, ChromosomeScanCallback);
  */

  /*
  VisitSetsInSep (sep, (Pointer) cfp, PopPhyMutEcoScanCallback);
  */

  /*
  SeriousSeqEntryCleanup (sep, NULL, NULL);
  RemoveAllNcbiCleanupUserObjects (sep);
  NormalizeDescriptorOrder (sep);

  bs = Se2Bs (sep);

  SeqEntryPubsAsn4Ex (sep, cfp->isEmblDdbj, TRUE);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    fprintf (cfp->logfp, "%s\n", cfp->buf);
    fflush (cfp->logfp);
  }
  BSFree (bs);
  bs = tmp;

  BSFree (bs);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, PseudoScanCallback);
  */

  /*
  VisitPubdescsInSep (sep, (Pointer) cfp, PubScanCallback);
  */

  /*
  Boolean  hasParts = FALSE;

  VisitSetsInSep (sep, (Pointer) &hasParts, CheckForParts);
  if (hasParts) {
    VisitFeaturesInSep (sep, (Pointer) cfp, SegmentedFeatCallback);
  }
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, TildeScanCallback);
  VisitPubdescsInSep (sep, (Pointer) cfp, CitGenCallback);
  */

  /*
  VisitBioseqsInSep (sep, (Pointer) cfp, FullMatProtCallback);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, OctoberReleaseCallback);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, SplitAnticodonCallback);
  */

  /*
  Uint2 entityID = ObjMgrGetEntityIDForChoice (sep);
  SeqMgrIndexFeatures (entityID, NULL);
  VisitBioseqsInSep (sep, (Pointer) cfp, SplitGeneBspCallback);
  */

  /*
  ByteStorePtr bs = NULL, tmp = NULL;
  Uint2        entityID;

  if (sep == NULL || cfp == NULL || cfp->logfp == NULL) return;

  entityID = ObjMgrGetEntityIDForChoice (sep);
  VisitDescriptorsInSep (sep, NULL, MarkComments);
  DeleteMarkedObjects (entityID, 0, NULL);

  VisitFeaturesInSep (sep, NULL, RemoveSfpComment);

  SeriousSeqEntryCleanup (sep, NULL, NULL);
  RemoveAllNcbiCleanupUserObjects (sep);
  NormalizeDescriptorOrder (sep);

  bs = Se2Bs (sep);
  FindReplaceInEntity (entityID, "  ", " ", FALSE, FALSE, TRUE, FALSE, 0, NULL, NULL, NULL, FALSE, NULL, NULL);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    fprintf (cfp->logfp, "SPACES %s\n", cfp->buf);
  }
  BSFree (bs);
  BSFree (tmp);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, ReportSeqLocMixGene);
  */

  /*
  VisitBioseqsInSep (sep, (Pointer) cfp, PlasmidScanCallback);
  */

  /*
  VisitBioseqsInSep (sep, (Pointer) cfp, SegmentedBioseqCallback);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, ReportSeqLocMixMix);
  */

  /*
  Uint2 entityID = ObjMgrGetEntityIDForChoice (sep);
  SeqMgrIndexFeatures (entityID, NULL);
  VisitBioseqsInSep (sep, (Pointer) cfp, InconsistentPartialsCallback);
  */

  /*
  DoSegSetFix (sep, cfp);
  */
}

static Boolean LatLonInRange (
  FloatHi lat,
  FloatHi lon
)

{
  if (lat < -90.0001 || lat > 90.0001) return FALSE;
  if (lon < -180.0001 || lon > 180.0001) return FALSE;

  return TRUE;
}

static Boolean CountryIsClosestToLatLon (
  CharPtr country,
  FloatHi lat,
  FloatHi lon,
  FloatHi range,
  FloatHi PNTR distanceP
)

{
  CharPtr  guess;

  if (StringHasNoText (country)) return FALSE;

  guess = CountryClosestToLatLon (lat, lon, range, distanceP);

  if (StringCmp (guess, country) == 0) return TRUE;

  return FALSE;
}

//LCOV_EXCL_START
// Unused option
static void LatLonScanCallback (
  BioSourcePtr biop,
  Pointer userdata
)

{
  Boolean       format_ok = FALSE, lat_in_range = FALSE, lon_in_range = FALSE, precision_ok = FALSE;
  Char          buf [512];
  CleanFlagPtr  cfp;
  CharPtr       countryname = NULL;
  FloatHi       distance = 0.0;
  CharPtr       guess;
  FloatHi       lat = 0.0;
  FloatHi       lon = 0.0;
  CharPtr       lat_lon = NULL;
  CharPtr       ptr;
  SubSourcePtr  ssp;
  Boolean       strict = TRUE;
  Char          tmp [128];

  if (biop == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    if (ssp->subtype == SUBSRC_lat_lon) {
      lat_lon = ssp->name;
    } else if (ssp->subtype == SUBSRC_country) {
      countryname = ssp->name;
    }
  }

  if (StringHasNoText (countryname)) return;

  /*
  fprintf (cfp->logfp, "CTRY\t%s\t%s\n", countryname, cfp->buf);
  */

  if (StringHasNoText (lat_lon)) return;

  IsCorrectLatLonFormat (lat_lon, &format_ok, &precision_ok, &lat_in_range, &lon_in_range);
  if (! format_ok) {
    /* may have comma and then altitude, so just get lat_lon component */
    StringNCpy_0 (tmp, lat_lon, sizeof (tmp));
    ptr = StringChr (tmp, ',');
    if (ptr != NULL) {
      *ptr = '\0';
      lat_lon = tmp;
      IsCorrectLatLonFormat (tmp, &format_ok, &precision_ok, &lat_in_range, &lon_in_range);
    }
  }

  if (format_ok && ParseLatLon (lat_lon, &lat, &lon)) {
    StringNCpy_0 (buf, countryname, sizeof (buf));
    ptr = StringChr (buf, ':');
    if (ptr != NULL) {
      *ptr = '\0';
      strict = FALSE;
    }
    if (CountryIsInLatLonList (buf)) {
      if (CountryContainsLatLon (buf, lat, lon)) {
        /* match */
        if (! strict) {
          StringNCpy_0 (buf, countryname, sizeof (buf));
          ptr = StringChr (buf, ',');
          if (ptr != NULL) {
            *ptr = '\0';
          }
          ptr = StringChr (buf, ';');
          if (ptr != NULL) {
            *ptr = '\0';
          }
          if (CountryIsInLatLonList (buf)) {
            if (CountryContainsLatLon (buf, lat, lon)) {
              fprintf (cfp->logfp, "%s\t%s\t%s\tCountry and subregion match\n", buf, lat_lon, cfp->buf);
            } else {
              if (! StringContainsBodyOfWater (countryname)) {
                guess = LookupCountryByLatLon (lat, lon);
                if (StringDoesHaveText (guess)) {
                  if (StringICmp (guess, buf) != 0) {
                    if (CountryContainsLatLon (guess, lat, lon) && CountryContainsLatLon (buf, lat, lon)) {
                      fprintf (cfp->logfp, "%s\t%s\t%s\tSubregion ambiguous with '%s'\n", buf, lat_lon, cfp->buf, guess);
                    } else if (CountryIsClosestToLatLon (buf, lat, lon, 1.0, &distance)) {
                      fprintf (cfp->logfp, "%s\t%s\t%s\tSelf subregion is closest\n", buf, lat_lon, cfp->buf);
                    } else {
                      fprintf (cfp->logfp, "%s\t%s\t%s\tSubregion might instead be '%s'\n", buf, lat_lon, cfp->buf, guess);
                    }
                  } else {
                    fprintf (cfp->logfp, "%s\t%s\t%s\tConfused subregion guess is '%s'\n", buf, lat_lon, cfp->buf, guess);
                  }
                } else {
                  if (CountryIsNearLatLon (buf, lat, lon, 1.0, &distance)) {
                    fprintf (cfp->logfp, "%s\t%s\t%s\tIs nearby but not inside\n", buf, lat_lon, cfp->buf);
                  } else {
                    fprintf (cfp->logfp, "%s\t%s\t%s\tDoes not map to subregion\n", buf, lat_lon, cfp->buf);
                  }
                }
              }
            }
          } else {
            fprintf (cfp->logfp, "%s\t%s\t%s\tCountry matches, unrecognized subregion\n", buf, lat_lon, cfp->buf);
          }
        } else {
          fprintf (cfp->logfp, "%s\t%s\t%s\tStrict match\n", buf, lat_lon, cfp->buf);
        }
      } else if (CountryIsClosestToLatLon (buf, lat, lon, 1.0, &distance)) {
        fprintf (cfp->logfp, "%s\t%s\t%s\tSelf country is closest\n", buf, lat_lon, cfp->buf);
      } else if (StringCmp (buf, "Antarctica") == 0 && lat < -60.0) {
        /* ignore Antarctica below 60 degrees S */
      } else if (LatLonInRange (-lat, lon) && CountryContainsLatLon (buf, -lat, lon)) {
        if (lat < 0.0) {
          fprintf (cfp->logfp, "%s\t%s\t%s\tLatitude should be set to N (northern hemisphere)\n", buf, lat_lon, cfp->buf);
        } else {
          fprintf (cfp->logfp, "%s\t%s\t%s\tLatitude should be set to S (southern hemisphere)\n", buf, lat_lon, cfp->buf);
        }
      } else if (LatLonInRange (lat, -lon) && CountryContainsLatLon (buf, lat, -lon)) {
        if (lon < 0.0) {
          fprintf (cfp->logfp, "%s\t%s\t%s\tLongitude should be set to E (eastern hemisphere)\n", buf, lat_lon, cfp->buf);
        } else {
          fprintf (cfp->logfp, "%s\t%s\t%s\tLongitude should be set to W (western hemisphere)\n", buf, lat_lon, cfp->buf);
        }
      } else if (LatLonInRange (lon, lat) && CountryContainsLatLon (buf, lon, lat)) {
        fprintf (cfp->logfp, "%s\t%s\t%s\tLatitude and longitude values appear to be exchanged\n", buf, lat_lon, cfp->buf);
      } else {
        if (! StringContainsBodyOfWater (countryname)) {
          guess = LookupCountryByLatLon (lat, lon);
          if (StringDoesHaveText (guess)) {
            if (StringICmp (guess, buf) != 0) {
              if (CountryContainsLatLon (guess, lat, lon) && CountryContainsLatLon (buf, lat, lon)) {
                fprintf (cfp->logfp, "%s\t%s\t%s\tCountry ambiguous with '%s'\n", buf, lat_lon, cfp->buf, guess);
              } else {
                fprintf (cfp->logfp, "%s\t%s\t%s\tCountry might instead be '%s'\n", buf, lat_lon, cfp->buf, guess);
              }
            } else {
              fprintf (cfp->logfp, "%s\t%s\t%s\tConfused country guess is '%s'\n", buf, lat_lon, cfp->buf, guess);
            }
          } else {
            if (CountryIsNearLatLon (buf, lat, lon, 1.0, &distance)) {
              fprintf (cfp->logfp, "%s\t%s\t%s\tNearby but not inside\n", buf, lat_lon, cfp->buf);
            } else {
              guess = CountryClosestToLatLon (lat, lon, 5.0, &distance);
              if (StringDoesHaveText (guess)) {
                fprintf (cfp->logfp, "%s\t%s\t%s\tClosest guess is '%s' at distance %f\n",
                         cfp->buf, lat_lon, buf, guess, (float) distance);
              } else {
                fprintf (cfp->logfp, "%s\t%s\t%s\tDoes not map to country\n", buf, lat_lon, cfp->buf);
              }
            }
          }
        }
      }
    } else {
      /*
      fprintf (cfp->logfp, "%s\t%s\t%s\tUnapproved country\n", buf, lat_lon, cfp->buf);
      */
    }

    fflush (cfp->logfp);
  }
}

// Unused option
static void DoLatLonReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  if (sep == NULL || cfp == NULL || cfp->logfp == NULL) return;

  VisitBioSourcesInSep (sep, (Pointer) cfp, LatLonScanCallback);
}

// Unused option
static void DoASNReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp,
  Boolean dossec,
  Boolean quick,
  Boolean popphymutdef
)

{
  Boolean              auth = FALSE, bsec = FALSE, clnr = FALSE, gbbk = FALSE,
                       modg = FALSE, modp = FALSE, modr = FALSE, move = FALSE,
                       norm = FALSE, orgm = FALSE, othr = FALSE, pack = FALSE,
                       prot = FALSE, publ = FALSE, ssec = FALSE, sloc = FALSE,
                       sort = FALSE, stcc = FALSE, sync = FALSE, titl = FALSE,
                       chkseg = FALSE, chknps = FALSE, dnarna = FALSE, ecnm = FALSE,
                       ncbiusrobj = FALSE, popphymuttitle = FALSE;
  ByteStorePtr         bs = NULL, tmp = NULL;
  ChangeData           cdbefore, cdafter;
  Int2                 chk_nuc_prot_val = 0;
  Uint2                entityID;
  CharPtr              set_type = NULL;
  PopSetRetroStatData  stat;

  if (sep == NULL || cfp == NULL) return;

  if (FindNcbiCleanupUserObject (sep) != NULL) {
    ncbiusrobj = TRUE;
  }

  MemSet ((Pointer) &stat, 0, sizeof (PopSetRetroStatData));

  RemoveAllNcbiCleanupUserObjects (sep);

  if (popphymutdef) {
    VisitSetsInSep (sep, (Pointer) &set_type, GetPopPhyMutEcoType);
    if (StringHasNoText (set_type)) return;

    entityID = ObjMgrGetEntityIDForChoice (sep);

    SeqMgrIndexFeatures (entityID, 0);
    PopSetAutoDefRetro (sep, &stat);
    LogPopSetTitleResults (cfp->logfp, &stat, sep);

    if (stat.title_added) {
      popphymuttitle = TRUE;
    }

    if (popphymuttitle) {
      if (cfp->logfp != NULL) {
        if (StringHasNoText (set_type)) {
          set_type = "UNK";
        }
        fprintf (cfp->logfp, "%s_SET %s\n", set_type, cfp->buf);
        fflush (cfp->logfp);
      }
    }

    return;
  }

  if (quick) {
    bs = Se2Bs (sep);

    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      norm = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    ResynchCodingRegionPartials (sep);
    ResynchMessengerRNAPartials (sep);
    ResynchProteinPartials (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      sync = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    CleanUpProteinTitles (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      prot = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    BasicSeqEntryCleanup (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      bsec = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    entityID = ObjMgrGetEntityIDForChoice (sep);
    SeqMgrIndexFeatures (entityID, NULL);
    SeqEntryExplore (sep, NULL, BadProtTitleProc);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateProteinTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      titl = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ssec = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    BSFree (bs);

    if (ssec) {
      (cfp->rawcounts.ssec)++;
      (cfp->cumcounts.ssec)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "SSEC %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    } else if (titl) {
      (cfp->rawcounts.titl)++;
      (cfp->cumcounts.titl)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "TITL %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    } else if (bsec) {
      (cfp->rawcounts.bsec)++;
      (cfp->cumcounts.bsec)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "BSEC %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    } else if (prot) {
      (cfp->rawcounts.prot)++;
      (cfp->cumcounts.prot)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "PROT %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    } else if (sync) {
      (cfp->rawcounts.sync)++;
      (cfp->cumcounts.sync)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "SYNC %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    } else if (norm) {
      (cfp->rawcounts.norm)++;
      (cfp->cumcounts.norm)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "NORM %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    } else {
      (cfp->rawcounts.okay)++;
      (cfp->cumcounts.okay)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "OKAY %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    }

    return;
  }

  MemSet ((Pointer) &cdbefore, 0, sizeof (ChangeData));
  MemSet ((Pointer) &cdafter, 0, sizeof (ChangeData));

  cdbefore.isRefSeq = cfp->isRefSeq;
  cdafter.isRefSeq = cfp->isRefSeq;

  CheckForChanges (sep, &cdbefore);

  bs = Se2Bs (sep);

  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    norm = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, CleanupLocation);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    sloc = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, CleanupMostRNAs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    clnr = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, CleanupRemainingRNAs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    othr = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, ModGenes);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    modg = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, ModRNAs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    modr = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitBioSourcesInSep (sep, NULL, ModPCRs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    modp = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, CleanupSubSourceOrgModOtherFeat);
  VisitDescriptorsInSep (sep, NULL, CleanupSubSourceOrgModOtherDesc);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    orgm = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitDescriptorsInSep (sep, NULL, CleanupStrucComment);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    stcc = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitPubdescsInSep (sep, NULL, CleanupPubAuthors);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    auth = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitPubdescsInSep (sep, (Pointer) cfp, CleanupPubBody);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    publ = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  SortSeqEntryQualifiers (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    sort = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitBioseqsInSep (sep, NULL, FixBspMol);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    dnarna = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  ResynchCodingRegionPartials (sep);
  ResynchMessengerRNAPartials (sep);
  ResynchProteinPartials (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    sync = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  CleanUpProteinTitles (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    prot = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  BasicSeqEntryCleanup (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    bsec = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  EntryChangeGBSource (sep);
  EntryCheckGBBlock (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    gbbk = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  entityID = ObjMgrGetEntityIDForChoice (sep);
  SeqMgrIndexFeatures (entityID, NULL);
  SeqEntryExplore (sep, NULL, BadProtTitleProc);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
  SeqMgrIndexFeatures (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);
  SeqMgrClearFeatureIndexes (entityID, NULL);
  BasicSeqEntryCleanup (sep);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    titl = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  UpdateReplacedECNumbers (sep);
  DeleteBadECNumbers (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    ecnm = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  MoveFeatsFromPartsSet (sep);
  move_cds_ex (sep, TRUE);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    pack = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  SeqEntryExplore(sep, NULL, ChkSegset);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    chkseg = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  SeqEntryExplore(sep, (Pointer) &chk_nuc_prot_val, ChkNucProt);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    chknps = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  SeqEntryPubsAsn4 (sep, cfp->isEmblDdbj);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    move = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  if (dossec) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ssec = TRUE;
    }
    BSFree (bs);
    bs = tmp;
  }

  BSFree (bs);

  CheckForChanges (sep, &cdafter);

  if (ssec || chkseg || chknps) {
    (cfp->rawcounts.ssec)++;
    (cfp->cumcounts.ssec)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SSEC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (move) {
    (cfp->rawcounts.move)++;
    (cfp->cumcounts.move)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MOVE %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (pack) {
    (cfp->rawcounts.pack)++;
    (cfp->cumcounts.pack)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PACK %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (titl) {
    (cfp->rawcounts.titl)++;
    (cfp->cumcounts.titl)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "TITL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (gbbk) {
    (cfp->rawcounts.gbbk)++;
    (cfp->cumcounts.gbbk)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GBBK %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (bsec || dnarna || orgm) {
    (cfp->rawcounts.bsec)++;
    (cfp->cumcounts.bsec)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "BSEC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (prot) {
    (cfp->rawcounts.prot)++;
    (cfp->cumcounts.prot)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PROT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (sync) {
    (cfp->rawcounts.sync)++;
    (cfp->cumcounts.sync)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SYNC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (sort) {
    (cfp->rawcounts.sort)++;
    (cfp->cumcounts.sort)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SORT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (sloc) {
    (cfp->rawcounts.sloc)++;
    (cfp->cumcounts.sloc)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SLOC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (clnr) {
    (cfp->rawcounts.clnr)++;
    (cfp->cumcounts.clnr)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CLNR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (othr) {
    (cfp->rawcounts.othr)++;
    (cfp->cumcounts.othr)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "OTHR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (modr || modg || modp) {
    (cfp->rawcounts.modr)++;
    (cfp->cumcounts.modr)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MODR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (publ) {
    (cfp->rawcounts.publ)++;
    (cfp->cumcounts.publ)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PUBL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (auth) {
    (cfp->rawcounts.auth)++;
    (cfp->cumcounts.auth)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "AUTH %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (norm) {
    (cfp->rawcounts.norm)++;
    (cfp->cumcounts.norm)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NORM %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else {
    (cfp->rawcounts.okay)++;
    (cfp->cumcounts.okay)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "OKAY %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (cdbefore.protNucSet) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PNS %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.noPopsetTitle) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PST %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.oldgbqual) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GBQ %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.sgml) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SGM %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.cdscodon) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CDN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.conflict) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CFT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.bad_part) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PAR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.bad_fuzz) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "FUZ %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rubisco) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RUB %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rbc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RBC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.its) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ITS %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rnaother) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RNA %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.trnanote) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "TRN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.oldbiomol) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MOL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.badOrg) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ORG %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rpt_unit_seq) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RUS %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (cdbefore.badDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "BDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.refDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "FDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.srcDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.capDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.privDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.oldDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ODX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.multDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rareDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.gbDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdafter.hasUnpublished && ! cdafter.hasPublished) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "UNP %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (sort) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SRT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (sloc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SLC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (clnr) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RCN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (othr) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RNO %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (modg) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GEN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (modr) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NCR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (modp) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PCR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (orgm) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SCN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (stcc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "STC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (publ) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PBC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (auth) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ATH %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (ecnm) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ECN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (pack) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PKG %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (move) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MVP %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (sync) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SNC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (prot) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PTL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (titl) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "TTL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (chkseg) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SEG %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (chknps) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NSS%d %s\n", (int) chk_nuc_prot_val, cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (dnarna) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MRN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (cdbefore.protdesc != cdafter.protdesc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PRT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.sfpnote != cdafter.sfpnote) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "COM %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.gbsource != cdafter.gbsource) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SRC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.cdsconf != cdafter.cdsconf) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CNF %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (ncbiusrobj) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "USR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
}

static CharPtr ffmod [] = {
  "",
  "release",
  "entrez",
  "gbench",
  "dump",
  NULL
};

static void DoGBFFReport (
  SeqEntryPtr sep,
  Uint2 entityID,
  CleanFlagPtr cfp,
  Int2 batch
)

{
#ifdef OS_UNIX
  AsnIoPtr      aip;
  Char          arguments [128];
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Char          ch;
  Char          cmmd [512];
  int           diff;
  FileCache     fc;
  FILE          *fp;
  SeqEntryPtr   fsep;
  Char          line [512];
  FILE          *ofp;
  Char          path1 [PATH_MAX];
  Char          path2 [PATH_MAX];
  Char          path3 [PATH_MAX];
  CharPtr       rep = "reports";
  SeqIdPtr      sip;
  CharPtr       str;
  SeqEntryPtr   tmp;

  if (sep == NULL || cfp == NULL) return;

  fsep = FindNthBioseq (sep, 1);
  if (fsep != NULL && fsep->choice == 1) {
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp != NULL) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_GENBANK :
            rep = "gbreports";
            break;
          case SEQID_EMBL :
            rep = "ebreports";
            break;
          case SEQID_DDBJ :
            rep = "djreports";
            break;
          case SEQID_OTHER :
            rep = "rfreports";
            break;
          default :
            break;
        }
      }
    }
  }

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s\n", cfp->buf);
    fflush (cfp->logfp);
  }

  if (batch == 1) {

    TmpNam (path1);
    TmpNam (path2);

    fp = FileOpen (path1, "w");
    if (fp != NULL) {
      SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, NORMAL_STYLE, 0, 0, 0, NULL, fp);
    }
    FileClose (fp);
    SeriousSeqEntryCleanupBulk (sep);
    fp = FileOpen (path2, "w");
    if (fp != NULL) {
      SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, NORMAL_STYLE, 0, 0, 0, NULL, fp);
    }
    FileClose (fp);

    sprintf (cmmd, "%s -o %s -n %s -d %s", cfp->ffdiff, path1, path2, rep);
    system (cmmd);

    sprintf (cmmd, "rm %s; rm %s", path1, path2);
    system (cmmd);

  } else if (batch == 2) {

    TmpNam (path1);
    TmpNam (path2);
    TmpNam (path3);

    aip = AsnIoOpen (path3, "w");
    if (aip == NULL) return;

    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);

    fp = FileOpen (path1, "w");
    if (fp != NULL) {
      SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, MASTER_STYLE, 0, 0, 0, NULL, fp);
    }
    FileClose (fp);

    arguments [0] = '\0';
    sprintf (arguments,
             "-format genbank -mode %s -style master -view nuc -gbload -cleanup -resolve-all -no-external -depth 1",
             ffmod [(int) cfp->ffmode]);

    sprintf (cmmd, "%s %s -i %s -o %s", cfp->asn2flat, arguments, path3, path2);
    system (cmmd);

    sprintf (cmmd, "diff -h %s %s > %s", path1, path2, path3);
    diff = system (cmmd);

    if (diff > 0) {
      fp = FileOpen (path3, "r");
      ofp = FileOpen (rep, "a");
      if (fp != NULL && ofp != NULL) {
        fprintf (ofp, "\n\n%s\n", cfp->buf);
        fflush (ofp);
        if (FileCacheSetup (&fc, fp)) {
          str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
          while (str != NULL) {
            ch = line [0];
            if (ch == '<' || ch == '>' || ch == '-') {
              fprintf (ofp, "%s\n", line);
            } else if (IS_DIGIT (ch)) {
              fprintf (ofp, "\n%s\n", "===");
            }
            str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
          }
        }
      }
      FileClose (ofp);
      FileClose (fp);
    }

    sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
    system (cmmd);

  } else if (batch == 3) {

    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp == NULL) return;
      /* skip if no parent nuc-prot set, code will handle this case after other records are retroed */
      if (bssp->_class == BioseqseqSet_class_segset) return;
      /* skip pop/phy/mut/eco/wgs sets for now */
      if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) return;
      if (bssp->_class == BioseqseqSet_class_wgs_set) return;
    }

    TmpNam (path1);
    TmpNam (path2);
    TmpNam (path3);

    SegSeqNullToVirtual (sep);

    fp = FileOpen (path1, "w");
    if (fp != NULL) {
      SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, MASTER_STYLE, 0, 0, 0, NULL, fp);
    }
    FileClose (fp);

    VisitBioseqsInSep (sep, NULL, RemoveSegSeqTitle);
    VisitSetsInSep (sep, NULL, RemoveSegSetTitle);
    /*
    VisitPubdescsInSep (sep, NULL, RemovePubSerial);
    */
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID, TRUE, FALSE);
    ConvertSegSetToDeltaSeq (sep);

    aip = AsnIoOpen (path3, "w");
    if (aip == NULL) return;

    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);

    aip = AsnIoOpen (path3, "r");
    if (aip == NULL) return;

    tmp = SeqEntryAsnRead (aip, NULL);
    AsnIoClose (aip);

    if (tmp != NULL) {
      fp = FileOpen (path2, "w");
      if (fp != NULL) {
        SeqEntryToGnbk (tmp, NULL, GENBANK_FMT, cfp->ffmode, MASTER_STYLE, 0, 0, 0, NULL, fp);
      }
      FileClose (fp);
      SeqEntryFree (tmp);
    }

    sprintf (cmmd, "%s -o %s -n %s -d %s", cfp->ffdiff, path1, path2, rep);
    system (cmmd);

    sprintf (cmmd, "rm %s; rm %s", path1, path2);
    system (cmmd);

  }
#endif
}

static void DoValidatorReport (
  SeqEntryPtr sep,
  FILE *logfp,
  CharPtr id,
  CharPtr asnval
)

{
#ifdef OS_UNIX
  AsnIoPtr   aip;
  Char       ch;
  Char       cmmd [512];
  int        diff;
  FileCache  fc;
  FILE       *fp;
  Char       line [512];
  Char       path1 [PATH_MAX];
  Char       path2 [PATH_MAX];
  Char       path3 [PATH_MAX];
  Char       path4 [PATH_MAX];
  Char       path5 [PATH_MAX];
  CharPtr    str;

  if (sep == NULL || logfp == NULL) return;
  if (StringHasNoText (id) || StringHasNoText (asnval)) return;

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);
  TmpNam (path4);
  TmpNam (path5);

  RemoveAllNcbiCleanupUserObjects (sep);

  aip = AsnIoOpen (path3, "w");
  if (aip == NULL) return;

  SeqEntryAsnWrite (sep, aip, NULL);
  AsnIoClose (aip);

  SeriousSeqEntryCleanup (sep, NULL, NULL);
  RemoveAllNcbiCleanupUserObjects (sep);

  aip = AsnIoOpen (path4, "w");
  if (aip == NULL) return;

  SeqEntryAsnWrite (sep, aip, NULL);
  AsnIoClose (aip);

  sprintf (cmmd, "%s -i %s -o stdout -Q 1 -r -l | sort > %s", asnval, path3, path1);
  system (cmmd);

  sprintf (cmmd, "%s -i %s -o stdout -Q 1 -r -l | sort > %s", asnval, path4, path2);
  system (cmmd);

  sprintf (cmmd, "diff -h %s %s > %s", path1, path2, path5);
  diff = system (cmmd);

  if (diff > 0) {
    fp = FileOpen (path5, "r");
    if (fp != NULL) {
      fprintf (logfp, "\n\n%s\n", id);
      fflush (logfp);
      if (FileCacheSetup (&fc, fp)) {
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
        while (str != NULL) {
          ch = line [0];
          if (ch == '<' || ch == '>' || ch == '-') {
            fprintf (logfp, "%s\n", line);
          } else if (IS_DIGIT (ch)) {
            fprintf (logfp, "\n%s\n", "===");
          }
          str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
        }
      }
      fflush (logfp);
    }
    FileClose (fp);
  }

  sprintf (cmmd, "rm %s; rm %s; rm %s; rm %s; rm %s", path1, path2, path3, path4, path5);
  system (cmmd);
#endif
}

static void DoModernizeReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  ByteStorePtr  bs = NULL, tmp = NULL;

  bs = Se2Bs (sep);

  VisitFeaturesInSep (sep, NULL, ModGenes);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GEN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, ModRNAs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NCR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  BSFree (bs);
  bs = tmp;

  VisitBioSourcesInSep (sep, NULL, ModPCRs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PCR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  BSFree (bs);
  bs = tmp;

  BSFree (bs);
}

// Used in DoAnticodonReport and SeqToAnticodon (both unused)
static void ExtractString (
  CharPtr str,
  Char PNTR buf,
  Int2 from,
  Int2 to
)

{
  Int2    i;
  Int2    j;
  size_t  len;

  if (buf == NULL) return;
  *buf = '\0';
  if (str == NULL) return;

  if (from < 0 || to < 0) return;
  len = StringLen (str);
  if (from > len || to > len) return;

  j = 0;
  for (i = from; i < to; i++) {
    buf [j] = str [i];
    j++;
  }
  buf [j] = '\0';
}

static void RevcompString (
  CharPtr str
)

{
  Char    ch;
  int     i;
  int     j;
  size_t  len;

  if (str == NULL) return;
  len = (int) StringLen (str);
  if (len < 1) return;

  for (i = 0, j = len - 1; i < j; i++, j--) {
    ch = str [i];
    str [i] = str [j];
    str [j] = ch;
  }

  for (i = 0; i < len; i++) {
    ch = str [i];
    if (ch == 'A') {
      ch = 'T';
    } else if (ch == 'C') {
      ch = 'G';
    } else if (ch == 'G') {
      ch = 'C';
    } else if (ch == 'T') {
      ch = 'A';
    }
    str [i] = ch;
  }
}

// Used in LongestSlidingMatch
static Int2 LongestExactMatch (
  CharPtr str1,
  CharPtr str2,
  Int2 PNTR gcp
)

{
  Char  ch1, ch2;
  Int2  gc = 0;
  int   i;
  Int2  len = 0;

  i = 0;
  ch1 = str1 [i];
  ch2 = str2 [i];
  while (ch1 != '\0' && ch2 != '\0' && ch1 == ch2) {
    len++;
    i++;
    if (ch1 == 'C' || ch1 == 'G') {
      gc++;
    }
    ch1 = str1 [i];
    ch2 = str2 [i];
  }

  *gcp = gc;

  return len;
}

// used in LongestMatch
static Int2 LongestSlidingMatch (
  CharPtr str1,
  CharPtr str2,
  Int2 PNTR gcp
)

{
  Int2    gc;
  Int2    len;
  Int2    max = 0;
  size_t  remains = 0;

  if (StringHasNoText (str1) || StringHasNoText (str2)) return 0;

  remains = StringLen (str1);

  while (remains > 3 && *str1 != '\0') {
    len = LongestExactMatch (str1, str2, &gc);
    if (len > max) {
      max = len;
      *gcp = gc;
    }
    str1++;
    remains--;
  }

  return max;
}

// Used in SeqToAnticodon and DoAnticodonReport
static Int2 LongestMatch (
  CharPtr str1,
  CharPtr str2,
  Int2 PNTR gcp
)

{
  Int2    gc = 0;
  Int2    len;
  Int2    max = 0;
  size_t  remains = 0;

  if (StringHasNoText (str1) || StringHasNoText (str2)) return 0;

  remains = StringLen (str2);

  while (remains > 3 && *str2 != '\0') {
    len = LongestSlidingMatch (str1, str2, &gc);
    if (len > max) {
      max = len;
      *gcp = gc;
    }
    str2++;
    remains--;
  }

  return max;
}

static void DoAnticodonReport (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Char          anticodon [16];
  Char          buf5 [16], buf3 [16], buf4 [16];
  CleanFlagPtr  cfp;
  Char          ch;
  Boolean       do_print;
  GBQualPtr     gbq;
  Int2          gc;
  int           i;
  int           j;
  CharPtr       locstr;
  Int2          max;
  CharPtr       palindrome;
  CharPtr       partial = " ";
  long          pos;
  CharPtr       prod;
  CharPtr       ptr;
  CharPtr       str;

  if (sfp == NULL || sfp->location == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->idx.subtype != FEATDEF_tRNA) return;

  anticodon [0] = '\0';
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next)  {
    if (StringICmp (gbq->qual, "anticodon") != 0) continue;
    if (StringLen (gbq->val) > 3) continue;
    StringCpy (anticodon, gbq->val);
    ptr = anticodon;
    ch = *ptr;
    while (ch != '\0') {
      ch = TO_UPPER (ch);
      if (ch == 'U') {
        ch = 'T';
      }
      *ptr = ch;
      ptr++;
      ch = *ptr;
    }
  }

  if (StringHasNoText (anticodon)) return;

  str = GetSequenceByFeature (sfp);
  if (str == NULL) return;

  if (sfp->partial) {
    partial = "*";
  }

  prod = GetFeatureProductName (sfp);
  if (StringHasNoText (prod)) {
    prod = "Xxx";
  }

  ptr = str;
  while (ptr != NULL) {
    ptr = StringStr (ptr, anticodon);
    if (ptr != NULL) {
      pos = (long) (ptr + 1 - str);
      ExtractString (str, buf5, pos - 10, pos - 3);
      ExtractString (str, buf3, pos + 4, pos + 11);
      if (StringHasNoText (buf5)) {
        StringCpy (buf5, "       ");
      }
      if (StringHasNoText (buf3)) {
        StringCpy (buf3, "       ");
      }
      StringCpy (buf4, "       ");

      palindrome = "-";
      max = 0;
      gc = 0;
      do_print = FALSE;
      if (StringDoesHaveText (buf5) && StringDoesHaveText (buf3)) {
        StringCpy (buf4, buf3);
        RevcompString (buf4);
        max = LongestMatch (buf5, buf4, &gc);
        if (max > 4) {
          palindrome = "+";
          do_print = TRUE;
        } else if (max == 4 && gc > 1) {
          palindrome = "~";
          do_print = TRUE;
        } else {
          do_print = TRUE;
        }
      }

      if (do_print) {
        locstr = SeqLocPrint (sfp->location);
        if (locstr == NULL) {
          locstr = StringSave ("?");
        }
        for (i = 0, j = (int) (ptr - str); i < 3; i++, j++) {
          ch = str [j];
          str [j] = TO_LOWER (ch);
        }
        fprintf (cfp->logfp, "%s (%s) - trna %s (%s) @ %2ld%s %s %d %s .. %s - %s\n",
                 cfp->buf, locstr, prod, anticodon, (long) pos, partial,
                 palindrome, (int) max, buf5, buf3, str);
        for (i = 0, j = (int) (ptr - str); i < 3; i++, j++) {
          ch = str [j];
          str [j] = TO_UPPER (ch);
        }
        MemFree (locstr);
        fflush (cfp->logfp);
      }

      ptr++;
    }
  }

  fprintf (cfp->logfp, "\n");

  MemFree (str);
}

static void DoOverlapReport (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  Char               buf0 [1000], buf1 [1000], buf2 [1000];
  CleanFlagPtr       cfp;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         feat;
  CharPtr            key, locstr1, locstr2;
  SeqIdPtr           sip, siphead;

  if (sfp == NULL || sfp->location == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (sfp->data.choice == SEQFEAT_GENE) return;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  siphead = SeqIdSetDup (bsp->id);
  for (sip = siphead; sip != NULL; sip = sip->next) {
    SeqIdStripLocus (sip);
  }
  SeqIdWrite (siphead, buf0, PRINTID_FASTA_LONG, sizeof (buf0));
  SeqIdSetFree (siphead);

  key = FindKeyFromFeatDefType (sfp->idx.subtype, FALSE);
  if (key == NULL) {
    key = "?";
  }

  locstr1 = SeqLocPrint (sfp->location);
  if (locstr1 == NULL) {
    locstr1 = StringSave ("?");
  }

  buf1 [0] = '\0';
  FeatDefLabel (sfp, buf1, sizeof (buf1) - 1, OM_LABEL_BOTH);

  feat = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_GENE, NULL, 0, NULL, CONTAINED_WITHIN, &fcontext);
  if (feat != NULL) {
    buf2 [0] = '\0';
    FeatDefLabel (feat, buf2, sizeof (buf2) - 1, OM_LABEL_BOTH);
    locstr2 = SeqLocPrint (feat->location);
    if (locstr2 == NULL) {
      locstr2 = StringSave ("?");
    }
    fprintf (cfp->logfp, "%s\t%s\t%s\t%s\t%s\n", buf0, buf1, locstr1, buf2, locstr2);
    MemFree (locstr2);
  }

  if (sfp->data.choice != SEQFEAT_CDREGION) return;

  feat = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, CHECK_INTERVALS, &fcontext);
  if (feat != NULL) {
    buf2 [0] = '\0';
    FeatDefLabel (feat, buf2, sizeof (buf2) - 1, OM_LABEL_BOTH);
    locstr2 = SeqLocPrint (feat->location);
    if (locstr2 == NULL) {
      locstr2 = StringSave ("?");
    }
    fprintf (cfp->logfp, "%s\t%s\t%s\t%s\t%s\n", buf0, buf1, locstr1, buf2, locstr2);
    MemFree (locstr2);
  }

  MemFree (locstr1);
}

static void DoContifOffsetReport (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  DeltaSeqPtr   dsp;
  Char          id [256];
  Int4          offset = 0;
  SeqIdPtr      sip;
  SeqLitPtr     slitp;
  SeqLocPtr     slp;
  /*
  Int4          gi;
  SeqIdPtr      sipset;
  */

  if (bsp == NULL || userdata == NULL) return;
  if (bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4 || bsp->seq_ext == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {

      slp = (SeqLocPtr) dsp->data.ptrvalue;
      if (slp == NULL || slp->choice == SEQLOC_NULL) continue;

      sip = SeqLocId (slp);
      if (sip == NULL) continue;

      SeqIdWrite (sip, id, PRINTID_FASTA_LONG, sizeof (id) - 1);

      /*
      if (sip->choice == SEQID_GI) {
        gi = (Int4) sip->data.intvalue;
        if (gi > 0) {
          sipset = GetSeqIdSetForGI (gi);
          if (sipset != NULL) {
            SeqIdWrite (sipset, id, PRINTID_FASTA_ALL, sizeof (id) - 1);
            SeqIdSetFree (sipset);
          }
        }
      }
      */

      fprintf (cfp->logfp, "%s\t%ld\n", id, (long) offset);

      offset += SeqLocLen (slp);

    } else if (dsp->choice == 2) {

      slitp = (SeqLitPtr) dsp->data.ptrvalue;
      if (slitp == NULL) continue;

      offset += slitp->length;
    }
  }
}
//LCOV_EXCL_STOP

static ValNodePtr GetAuthorMLNameList (
  AuthListPtr alp
)

{
  AuthorPtr    ap;
  Char         buf [512];
  Char         ch;
  Char         chr [4];
  ValNodePtr   head = NULL;
  Char         initials [32];
  ValNodePtr   last = NULL;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  CharPtr      ptr;
  CharPtr      str;
  ValNodePtr   tmp;
  ValNodePtr   vnp;

  if (alp == NULL) return NULL;

  for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
    buf [0] = '\0';
    initials [0] = '\0';
    switch (alp->choice) {
      case 1 :
        ap = (AuthorPtr) vnp->data.ptrvalue;
        if (ap == NULL) continue;
        pid = ap->name;
        if (pid == NULL) continue;
        if (pid->choice == 2) {
          nsp = pid->data;
          if (nsp == NULL) continue;
          str = nsp->names [0];
          if (StringHasNoText (str)) continue;
          StringNCpy_0 (buf, str, sizeof (buf) - 50);
          StringNCpy_0 (initials, nsp->names [4], sizeof (initials));
        } else if (pid->choice == 3) {
          str = (CharPtr) pid->data;
          if (StringHasNoText (str)) continue;
          StringNCpy_0 (buf, str, sizeof (buf) - 50);
        }
        break;
      case 2 :
      case 3 :
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        StringNCpy_0 (buf, str, sizeof (buf));
        ptr = StringChr (buf, ',');
        if (ptr == NULL) {
          ptr = StringChr (buf, ' ');
        }
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
          StringNCpy_0 (initials, ptr, sizeof (initials));
        }
        break;
      default :
        break;
    }
    if (StringHasNoText (buf)) continue;
    if (StringDoesHaveText (initials)) {
      StringCat (buf, " ");
      chr [1] = '\0';
      ptr = initials;
      ch = *ptr;
      while (ch != '\0') {
        if (ch != ' ' && ch != '.' && ch != ',') {
          chr [0] = ch;
          StringCat (buf, chr);
        }
        ptr++;
        ch = *ptr;
      }
    }
    TrimSpacesAroundString (buf);
    tmp = ValNodeCopyStr (&last, 0, buf);
    if (head == NULL) {
      head = tmp;
    }
    last = tmp;
  }

  return head;
}

static ValNodePtr DuplicateStringList (
  ValNodePtr list
)

{
  ValNodePtr  head = NULL;
  ValNodePtr  last = NULL;
  CharPtr     str;
  ValNodePtr  tmp;
  ValNodePtr  vnp;

  if (list == NULL) return NULL;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    tmp = ValNodeCopyStr (&last, 0, str);
    if (head == NULL) {
      head = tmp;
    }
    last = tmp;
  }

  return head;
}

typedef enum {
  FULL_INITIALS,
  NO_HYPHENS,
  TWO_INITIALS,
  ONE_INITIAL,
  NO_INITIALS,
  NO_AUTHORS
} InitialsPolicy;

static void TrimInitials (
  CharPtr auth,
  InitialsPolicy initials
)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (StringHasNoText (auth)) return;

  switch (initials) {
    case FULL_INITIALS :
      break;
    case NO_HYPHENS :
      dst = auth;
      ptr = auth;
      ch = *ptr;
      while (ch != '\0') {
        if (ch != '-') {
          *dst = ch;
          dst++;
        }
        ptr++;
        ch = *ptr;
      }
      *dst = '\0';
      break;
    case TWO_INITIALS :
      ptr = StringRChr (auth, ' ');
      if (ptr != NULL) {
        ptr++;
        ch = *ptr;
        if (IS_ALPHANUM (ch)) {
          ptr++;
          ch = *ptr;
          if (IS_ALPHANUM (ch)) {
            ptr++;
            *ptr = '\0';
          }
        }
      }
      break;
    case ONE_INITIAL :
      ptr = StringRChr (auth, ' ');
      if (ptr != NULL) {
        ptr++;
        ch = *ptr;
        if (IS_ALPHANUM (ch)) {
          ptr++;
          *ptr = '\0';
        }
      }
      break;
    case NO_INITIALS :
      ptr = StringRChr (auth, ' ');
      if (ptr != NULL) {
        *ptr = '\0';
      }
      break;
    default :
      break;
  }
}

static Boolean FirstOrLastAuthorMatches (
  ValNodePtr oldauthors,
  ValNodePtr newauthors
)

{
  Char        buf [128];
  CharPtr     firstauthor = NULL, lastauthor = NULL, str;
  ValNodePtr  vnp;

  if (oldauthors == NULL || newauthors == NULL) return FALSE;

  for (vnp = oldauthors; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (vnp != oldauthors && vnp->next != NULL) continue;
    StringNCpy_0 (buf, str, sizeof (buf));
    TrimInitials (buf, ONE_INITIAL);
    if (firstauthor == NULL) {
    firstauthor = StringSave (buf);
    } else if (lastauthor == NULL) {
      lastauthor = StringSave (buf);
    }
  }

  for (vnp = newauthors; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    StringNCpy_0 (buf, str, sizeof (buf));
    TrimInitials (buf, ONE_INITIAL);
    if (StringICmp (buf, firstauthor) == 0) return TRUE;
    if (StringICmp (buf, lastauthor) == 0) return TRUE;
  }

  return FALSE;
}

static Int4 DoUnpubHydraQuery (
  ValNodePtr authors,
  InitialsPolicy initials,
  Boolean firstLastOnly,
  ValNodePtr titlewords,
  Int2 year,
  Boolean expand,
  Uint4Ptr uidp
)

{
  XmlObjPtr        attr, tmp, xop;
  Char             buf [128];
  CitArtPtr        cap;
  Int4             count = 0, pmid = 0, pmval;
  ValNodePtr       head = NULL, tail = NULL, vnp;
  CharPtr          idstr, query, score, str;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pmep;
  long int         val;

  if (uidp != NULL) {
    *uidp = 0;
  }

  for (vnp = titlewords; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    /*
    if (IsStopWord (str)) continue;
    */
    StringNCpy_0 (buf, str, sizeof (buf));
    ValNodeCopyStrEx (&head, &tail, 0, buf);
  }

  if (initials != NO_AUTHORS) {
    for (vnp = authors; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      if (firstLastOnly) {
        if (vnp != authors && vnp->next != NULL) continue;
      }
      StringNCpy_0 (buf, str, sizeof (buf));
      switch (initials) {
        case NO_HYPHENS :
          TrimInitials (buf, NO_HYPHENS);
          break;
        case TWO_INITIALS :
          TrimInitials (buf, TWO_INITIALS);
          break;
        case ONE_INITIAL :
          TrimInitials (buf, ONE_INITIAL);
          break;
        case NO_INITIALS :
          TrimInitials (buf, NO_INITIALS);
          break;
        default :
          break;
      }
      ValNodeCopyStrEx (&head, &tail, 0, buf);
    }
  }

  if (year > 0) {
    if (expand) {
      sprintf (buf, "(%d OR %d OR %d)", (int) year - 1, (int) year, (int) year + 1);
      ValNodeCopyStrEx (&head, &tail, 0, buf);
    } else {
      sprintf (buf, "%d", (int) year);
      ValNodeCopyStrEx (&head, &tail, 0, buf);
    }
  }

  query = ValNodeMergeStrsEx (head, "+");
  ValNodeFreeData (head);

  if (StringHasNoText (query)) return 0;

  str = QUERY_UrlSynchronousQuery ("intranet.ncbi.nlm.nih.gov", 80,
                                   "/projects/hydra/hydra_search.cgi",
                                   "search=pubmed_search_citation_top_20.1&query=",
                                   query, NULL, NULL);
  if (str == NULL) return 0;

  xop = ParseXmlString (str);
  if (xop == NULL) return 0;

  head = NULL;
  tail = NULL;

  for (tmp = xop; tmp != NULL; tmp = tmp->successor) {
    if (XmlPathSuffixIs (tmp, "/IdList/Id")) {
      if (StringHasNoText (tmp->contents)) continue;
      for (attr = tmp->attributes; attr != NULL; attr = attr->next) {
        if (StringICmp (attr->name, "score") != 0) continue;
        score = attr->contents;
        if (StringHasNoText (score)) continue;
        if (StringChr (score, '-') != NULL) continue;
        if (StringNCmp (score, "1", 1) == 0 || StringNCmp (score, "0.9", 3) == 0 || StringNCmp (score, "0.8", 3) == 0) {
          ValNodeCopyStrEx (&head, &tail, 0, tmp->contents);
        }
      }
    }
  }

  if (head == NULL) return 0;

  count = ValNodeLen (head);

  for (vnp = head; vnp != NULL && pmid == 0; vnp = vnp->next) {
    idstr = (CharPtr) vnp->data.ptrvalue;
    if (StringDoesHaveText (idstr)) {
      if (sscanf (idstr, "%ld", &val) == 1) {
        pmval = (Int4) val;
        if (pmval == 0) continue;
        pmep = PubMedSynchronousQuery (pmval);
        if (pmep != NULL) {
          mep = (MedlineEntryPtr) pmep->medent;
          if (mep != NULL) {
            cap = mep->cit;
            if (cap != NULL && cap->from == 1) {
              if (FirstOrLastAuthorMatches (authors, GetAuthorMLNameList (cap->authors))) {
                pmid = pmval;
              }
            }
          }
          pmep = PubmedEntryFree (pmep);
        }
      }
    }
  }

  ValNodeFreeData (head);

  FreeXmlObject (xop);

  MemFree (str);

  if (uidp != NULL) {
    *uidp = pmid;
  }

  return count;
}

static CharPtr GetBestJournal (
  ValNodePtr journaltitle
)

{
  CharPtr     str;
  ValNodePtr  vnp;

  if (journaltitle == NULL) return NULL;

  for (vnp = journaltitle; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Cit_title_iso_jta) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      return str;
    }
  }

  for (vnp = journaltitle; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Cit_title_name || vnp->choice == Cit_title_jta) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      return str;
    }
  }

  return NULL;
}

typedef struct pubref {
  Int4        unpubcount;
  ValNodePtr  seqids;
  ValNodePtr  authors;
  Boolean     firstLastOnly;
  ValNodePtr  titlewords;
  CharPtr     fulltitle;
  CharPtr     uniquestr;
  CharPtr     journal;
  ImprintPtr  imp;
  Int2        year;
  Uint4       pmid;
} PubRef, PNTR PubRefPtr;

static void PrintPubAuthors (
  CleanFlagPtr cfp,
  PubRefPtr prp
)

{
  CharPtr     prefix = "";
  CharPtr     str;
  ValNodePtr  vnp;

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL) return;

  for (vnp = prp->authors; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    fprintf (cfp->logfp, "%s%s", prefix, str);
    prefix = ", ";
  }
}

static void PrintPubTitle (
  CleanFlagPtr cfp,
  PubRefPtr prp
)

{
  CharPtr     prefix = "";
  CharPtr     str;
  ValNodePtr  vnp;

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL) return;

  if (StringDoesHaveText (prp->fulltitle)) {
    fprintf (cfp->logfp, "%s%s", prefix, prp->fulltitle);
  } else {
    for (vnp = prp->titlewords; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      fprintf (cfp->logfp, "%s%s", prefix, str);
      prefix = " ";
    }
  }
}

static void PrintPubJournal (
  CleanFlagPtr cfp,
  PubRefPtr prp
)

{
  DatePtr     dp = NULL;
  ImprintPtr  imp;
  CharPtr     prefix = "";
  Int2        year;

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL) return;

  if (StringHasNoText (prp->journal) && prp->imp == NULL) {
    fprintf (cfp->logfp, "Unpublished");
    prefix = " ";
    if (prp->year > 0) {
      fprintf (cfp->logfp, "%s[%d]", prefix, (int) prp->year);
      prefix = " ";
    }
    return;
  }

  if (StringDoesHaveText (prp->journal)) {
    fprintf (cfp->logfp, "%s%s", prefix, prp->journal);
    prefix = " ";
  }

  imp = prp->imp;
  if (imp != NULL) {
    dp = imp->date;
    if (dp != NULL && dp->data [0] == 1) {
      year = (Int2) dp->data [1] + 1900;
      fprintf (cfp->logfp, "%s[%d]", prefix, (int) year);
      prefix = " ";
    }
    if (StringDoesHaveText (imp->volume)) {
      fprintf (cfp->logfp, "%s%s", prefix, imp->volume);
      prefix = " ";
    }
    /*
    if (StringDoesHaveText (imp->issue)) {
      fprintf (cfp->logfp, "%s(%s)", prefix, imp->issue);
      prefix = " ";
    }
    */
    if (StringDoesHaveText (imp->pages)) {
      fprintf (cfp->logfp, "%s: %s", prefix, imp->pages);
      prefix = " ";
    }
  }

  if (prp->pmid > 0) {
    fprintf (cfp->logfp, "%s<%ld>", prefix, (long) prp->pmid);
    prefix = " ";
  }
}

typedef enum {
  NO_NAME_MATCH,
  LAST_NAME_MATCH,
  ONE_INIT_MATCH,
  TWO_INIT_MATCH,
  NO_HYPHEN_MATCH,
  FULL_NAME_MATCH
} AuthComp;

static CharPtr authlabel [] = {
  "AUTH_MISMATCH", "LAST_NAMES", "ONE_INIT", "TWO_INITS", "NO_HYPHENS", "FULL_NAMES"
};

static AuthComp AuthorCompare (
  CharPtr auth1,
  CharPtr auth2
)

{
  Char  buf1 [128];
  Char  buf2 [128];

  if (StringHasNoText (auth1) || StringHasNoText (auth2)) return NO_NAME_MATCH;

  StringNCpy_0 (buf1, auth1, sizeof (buf1));
  StringNCpy_0 (buf2, auth2, sizeof (buf2));

  if (StringICmp (buf1, buf2) == 0) return FULL_NAME_MATCH;

  TrimInitials (buf1, NO_HYPHENS);
  TrimInitials (buf2, NO_HYPHENS);

  if (StringICmp (buf1, buf2) == 0) return NO_HYPHEN_MATCH;

  TrimInitials (buf1, TWO_INITIALS);
  TrimInitials (buf2, TWO_INITIALS);

  if (StringICmp (buf1, buf2) == 0) return TWO_INIT_MATCH;

  TrimInitials (buf1, ONE_INITIAL);
  TrimInitials (buf2, ONE_INITIAL);

  if (StringICmp (buf1, buf2) == 0) return ONE_INIT_MATCH;

  TrimInitials (buf1, NO_INITIALS);
  TrimInitials (buf2, NO_INITIALS);

  if (StringICmp (buf1, buf2) == 0) return LAST_NAME_MATCH;

  return NO_NAME_MATCH;
}

static AuthComp AuthorsIdentical (
  ValNodePtr oldauthors,
  ValNodePtr newauthors
)

{
  AuthComp    curr, rsult = FULL_NAME_MATCH;
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (oldauthors == NULL || newauthors == NULL) return NO_NAME_MATCH;

  for (vnp1 = oldauthors, vnp2 = newauthors;
       vnp1 != NULL && vnp2 != NULL;
       vnp1 = vnp1->next, vnp2 = vnp2->next) {
    str1 = (CharPtr) vnp1->data.ptrvalue;
    str2 = (CharPtr) vnp2->data.ptrvalue;
    curr = AuthorCompare (str1, str2);
    if (curr == NO_NAME_MATCH) return NO_NAME_MATCH;
    if (curr < rsult) {
      rsult = curr;
    }
  }

  if (vnp1 != NULL || vnp2 != NULL) return NO_NAME_MATCH;

  return rsult;
}

static AuthComp AuthorInList (
  CharPtr author,
  ValNodePtr newauthors
)

{
  AuthComp    curr, rsult = FULL_NAME_MATCH;
  Boolean     okay = FALSE;
  CharPtr     str;
  ValNodePtr  vnp;

  if (StringHasNoText (author) || newauthors == NULL) return NO_NAME_MATCH;

  for (vnp = newauthors; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    curr = AuthorCompare (author, str);
    if (curr == NO_NAME_MATCH) continue;
    okay = TRUE;
    if (curr < rsult) {
      rsult = curr;
    }
  }

  if (! okay) return NO_NAME_MATCH;

  return rsult;
}

static Boolean WordInList (
  CharPtr word,
  ValNodePtr newtitlewords
)

{
  CharPtr     str;
  ValNodePtr  vnp;

  if (StringHasNoText (word) || newtitlewords == NULL) return NO_NAME_MATCH;

  for (vnp = newtitlewords; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringICmp (word, str) == 0) return TRUE;
  }

  return FALSE;
}

static void PrintComparison (
  CleanFlagPtr cfp,
  PubRefPtr oldprp,
  PubRefPtr newprp
)

{
  AuthComp    authcomp, curr, best = FULL_NAME_MATCH;
  Boolean     bothokay = TRUE, titlsame;
  Int2        matches, newcount, oldcount;
  CharPtr     str, str1, str2;
  ValNodePtr  vnp;

  if (cfp == NULL || cfp->logfp == NULL || oldprp == NULL || newprp == NULL) return;

  authcomp = AuthorsIdentical (oldprp->authors, newprp->authors);
  titlsame = (Boolean) (StringICmp (oldprp->fulltitle, newprp->fulltitle) == 0);

  fprintf (cfp->logfp, "PMID %ld", (long) newprp->pmid);
  fprintf (cfp->logfp, "\t");

  if (oldprp->seqids != NULL) {
    oldprp->seqids = ValNodeSort (oldprp->seqids, SortVnpByString);
    for (vnp = oldprp->seqids; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      fprintf (cfp->logfp, "SEQID %s", str);
      fprintf (cfp->logfp, "\t");
    }
  } else {
    fprintf (cfp->logfp, "%s", cfp->buf);
    fprintf (cfp->logfp, "\t");
  }

  /*
  fprintf (cfp->logfp, "REF_COUNT %ld", (long) oldprp->unpubcount);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "ORIG_NAMES %ld", (long) ValNodeLen (oldprp->authors));
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "ADDL_NAMES %ld", (long) (ValNodeLen (newprp->authors) - ValNodeLen (oldprp->authors)));
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "ORIG_WORDS %ld", (long) ValNodeLen (oldprp->titlewords));
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "ADDL_WORDS %ld", (long) (ValNodeLen (newprp->titlewords) - ValNodeLen (oldprp->titlewords)));
  fprintf (cfp->logfp, "\t");
  */

  if (StringDoesHaveText (oldprp->uniquestr)) {
    fprintf (cfp->logfp, "UNIQ_CIT %s", oldprp->uniquestr);
  } else {
    fprintf (cfp->logfp, "?");
  }
  fprintf (cfp->logfp, "\t");

  if (authcomp != NO_NAME_MATCH) {
    str = authlabel [(int) authcomp];
    if (oldprp->firstLastOnly) {
      fprintf (cfp->logfp, "AUTHORS_FIRST_LAST [%s]", str);
    } else {
      fprintf (cfp->logfp, "AUTHORS_SAME [%s]", str);
    }
    fprintf (cfp->logfp, "\t");
  } else {
    newcount = ValNodeLen (newprp->authors);
    oldcount = ValNodeLen (oldprp->authors);

    matches = 0;
    for (vnp = oldprp->authors; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      curr = AuthorInList (str, newprp->authors);
      if (curr == NO_NAME_MATCH) continue;
      matches++;
      if (curr < best) {
        best = curr;
      }
    }

    str = authlabel [(int) best];
    if (newcount > 0 && matches == oldcount) {
      if (oldcount < 3 && newcount > 4) {
        fprintf (cfp->logfp, "AUTHORS_QUESTIONABLE [%s] %d -> %d ", str, (int) oldcount, (int) newcount);
        bothokay = FALSE;
      } else if (oldcount < newcount) {
        fprintf (cfp->logfp, "AUTHORS_ADDED [%s] %d ", str, (int) (newcount - oldcount));
      } else {
        fprintf (cfp->logfp, "AUTHORS_REORDERED [%s]", str);
      }
    } else if (oldprp->firstLastOnly) {
      fprintf (cfp->logfp, "AUTHORS_REMOVED [%s] %d", str, (int) matches);
      bothokay = FALSE;
    } else {
      fprintf (cfp->logfp, "AUTHORS_CHANGED [%s] %d / %d", str, (int) matches, (int) newcount);
      bothokay = FALSE;
    }
    fprintf (cfp->logfp, "\t");
  }

  if (titlsame) {
    fprintf (cfp->logfp, "TITLE_SAME [IDENTICAL]");
  } else {
    newcount = 0;
    oldcount = 0;
    matches = 0;

    for (vnp = newprp->titlewords; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      /*
      if (IsStopWord (str)) continue;
      */
      newcount++;
    }

    for (vnp = oldprp->titlewords; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      /*
      if (IsStopWord (str)) continue;
      */
      oldcount++;
      if (! WordInList (str, newprp->titlewords)) continue;
      matches++;
    }

    str1 = NULL;
    str2 = NULL;
    vnp = oldprp->titlewords;
    if (vnp != NULL) {
      str1 = (CharPtr) vnp->data.ptrvalue;
      /*
      vnp = vnp->next;
      while (StringDoesHaveText (str1) && IsStopWord (str1) && vnp != NULL) {
        str1 = (CharPtr) vnp->data.ptrvalue;
        vnp = vnp->next;
      }
      */
    }
    vnp = newprp->titlewords;
    if (vnp != NULL) {
      str2 = (CharPtr) vnp->data.ptrvalue;
      /*
      vnp = vnp->next;
      while (StringDoesHaveText (str2) && IsStopWord (str2) && vnp != NULL) {
        str2 = (CharPtr) vnp->data.ptrvalue;
        vnp = vnp->next;
      }
      */
    }

    if (oldcount < 3 && newcount > 4) {
      fprintf (cfp->logfp, "TITLE_QUESTIONABLE %d -> %d", (int) oldcount, (int) newcount);
      bothokay = FALSE;
    } else if (StringDoesHaveText (str1) && StringDoesHaveText (str2) &&
        StringICmp (str1, str2) == 0 && newcount > 0 && matches == newcount) {
      fprintf (cfp->logfp, "TITLE_SAME [SIMILAR] %d", (int) matches);
    } else if (newcount > 0 && matches == newcount) {
      fprintf (cfp->logfp, "TITLE_ALTERED %d", (int) matches);
      bothokay = FALSE;
    } else {
      fprintf (cfp->logfp, "TITLE_DIFFERS %d / %d", (int) matches, (int) newcount);
      bothokay = FALSE;
    }
  }
  fprintf (cfp->logfp, "\t");

  if (bothokay) {
    fprintf (cfp->logfp, "PROBABLE");
  } else {
    fprintf (cfp->logfp, "POSSIBLE");
  }
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "OLD_AUTH ");
  PrintPubAuthors (cfp, oldprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "NEW_AUTH ");
  PrintPubAuthors (cfp, newprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "OLD_TITL ");
  PrintPubTitle (cfp, oldprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "NEW_TITL ");
  PrintPubTitle (cfp, newprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "OLD_JOUR ");
  PrintPubJournal (cfp, oldprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "NEW_JOUR ");
  PrintPubJournal (cfp, newprp);

  fprintf (cfp->logfp, "\n");
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
        new_str--;
      }
    } else {
      str++;
    }
  }
  *new_str = '\0';
}

static void StrStripBrackets (
  CharPtr str
)

{
  size_t  len;

  if (str == NULL) return;

  len = StringLen (str);
  if (len < 2) return;

  if (str [0] == '[') {
    str [0] = ' ';
  }

  if (str [len - 1] == ']') {
    str [len - 1] = ' ';
  }
}

static void PrintPubMedCit (
  CleanFlagPtr cfp,
  Uint4 pmid,
  PubRefPtr oldprp
)

{
  CitArtPtr        cap;
  CitJourPtr       cjp;
  DatePtr          dp;
  ImprintPtr       imp;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pep;
  PubRef           pr;
  CharPtr          str;
  CharPtr          tmp;
  ValNodePtr       vnp;

  if (cfp == NULL || cfp->logfp == NULL || pmid < 1) return;

  MemSet ((Pointer) &pr, 0, sizeof (PubRef));

  pep = PubMedSynchronousQuery (pmid);
  if (pep == NULL) return;

  mep = (MedlineEntryPtr) pep->medent;
  if (mep != NULL && mep->cit != NULL) {
    cap = mep->cit;
    if (cap != NULL) {
      pr.authors = GetAuthorMLNameList (cap->authors);
      for (vnp = cap->title; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == Cit_title_name) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          pr.titlewords = GetTitleWords (str);
          tmp = StringSave (str);
          TrimSpacesAndJunkFromEnds (tmp, TRUE);
          s_RemovePeriodFromEnd (tmp);
          StrStripBrackets (tmp);
          StrStripSpaces (tmp);
          TrimSpacesAroundString (tmp);
          pr.fulltitle = tmp;
        }
      }
      if (cap->from == 1) {
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp != NULL) {
          pr.journal = GetBestJournal (cjp->title);
          imp = cjp->imp;
          pr.imp = imp;
          if (imp != NULL) {
            dp = imp->date;
            if (dp != NULL && dp->data [0] == 1) {
              pr.year = (Int2) dp->data [1] + 1900;
            }
          }
        }
      }
      pr.pmid = pmid;
      if (pr.authors != NULL && pr.titlewords != NULL) {
        PrintComparison (cfp, oldprp, &pr);
      }
      ValNodeFreeData (pr.authors);
      ValNodeFreeData (pr.titlewords);
      MemFree (pr.fulltitle);
    }
  }

  pep = PubmedEntryFree (pep);
}

static void TryEntrezQueries (
  CleanFlagPtr cfp,
  PubRefPtr prp
)

{
  Int4        count;
  Uint4       pmid = 0;
  /*
  CharPtr     str;
  ValNodePtr  vnp;
  */

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL || prp->authors == NULL) return;

  count = DoUnpubHydraQuery (prp->authors, NO_AUTHORS, FALSE, prp->titlewords, 0, FALSE, &pmid);

  /*
  count = DoUnpubHydraQuery (prp->authors, NO_INITIALS, TRUE, prp->titlewords, prp->year, FALSE, &pmid);

  if (count != 1) {
    count = DoUnpubHydraQuery (prp->authors, ONE_INITIAL, FALSE, prp->titlewords, prp->year, TRUE, &pmid);
  }

  if (count > 1) {
    count = DoUnpubHydraQuery (prp->authors, TWO_INITIALS, FALSE, prp->titlewords, prp->year, TRUE, &pmid);
  }

  if (count < 1) {
    count = DoUnpubHydraQuery (prp->authors, ONE_INITIAL, TRUE, prp->titlewords, prp->year, TRUE, &pmid);
    if (count == 1) {
      prp->firstLastOnly = TRUE;
    }
  }
  */

  if (count < 1) {

    /*
    if (prp->seqids != NULL) {
      fprintf (cfp->logfp, "0\t");
      for (vnp = prp->seqids; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        fprintf (cfp->logfp, "SEQID %s", str);
        fprintf (cfp->logfp, "\t");
      }
      fprintf (cfp->logfp, "UNPUB\t");
    } else {
      fprintf (cfp->logfp, "0\t%s\tUNPUB\t", cfp->buf);
    }

    PrintPubAuthors (cfp, prp);
    fprintf (cfp->logfp, "\t");
    PrintPubTitle (cfp, prp);
    fprintf (cfp->logfp, "\t");
    PrintPubJournal (cfp, prp);
    fprintf (cfp->logfp, "\n");
    */

  } else if (count > 1) {

    /*
    if (prp->seqids != NULL) {
      fprintf (cfp->logfp, "0\t");
      for (vnp = prp->seqids; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        fprintf (cfp->logfp, "SEQID %s", str);
        fprintf (cfp->logfp, "\t");
      }
      fprintf (cfp->logfp, "COUNT %ld\t", (long) count);
    } else {
      fprintf (cfp->logfp, "0\t%s\tCOUNT %ld\t", cfp->buf, (long) count);
    }

    PrintPubAuthors (cfp, prp);
    fprintf (cfp->logfp, "\t");
    PrintPubTitle (cfp, prp);
    fprintf (cfp->logfp, "\t");
    PrintPubJournal (cfp, prp);
    fprintf (cfp->logfp, "\n");
    */

  } else {

    PrintPubMedCit (cfp, pmid, prp);
  }

  fflush (cfp->logfp);
}

static void CountUnpubPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CitGenPtr     cgp = NULL;
  Boolean       hasUnpublished = FALSE;
  ValNodePtr    vnp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            hasUnpublished = TRUE;
          }
        }
      }
    } else if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      return;
    } else if (vnp->choice == PUB_Article || vnp->choice == PUB_Book || vnp->choice == PUB_Man) {
      return;
    }
  }

  if (! hasUnpublished) return;
  if (cgp == NULL) return;

  (cfp->unpubcount)++;
}

//LCOV_EXCL_START
static void ProcessUnpubPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  AuthListPtr   authors = NULL;
  Char          buf [521];
  CitArtPtr     cap = NULL;
  CleanFlagPtr  cfp;
  CitGenPtr     cgp = NULL;
  CitJourPtr    cjp;
  DatePtr       dp = NULL;
  Boolean       hasUnpublished = FALSE;
  ImprintPtr    imp;
  CharPtr       journal = NULL;
  PubRef        pr;
  CharPtr       str;
  CharPtr       title = NULL;
  CharPtr       tmp;
  ValNodePtr    vnp, vnp2, vnpcgp = NULL;
  Int2          year = 0;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            hasUnpublished = TRUE;
            vnpcgp = vnp;
            authors = cgp->authors;
            title = cgp->title;
            journal = GetBestJournal (cgp->journal);
            dp = cgp->date;
          }
        }
      }
    } else if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      return;
    } else if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
      if (cap != NULL) {
        if (cap->from == 1) {
          cjp = (CitJourPtr) cap->fromptr;
          if (cjp != NULL) {
            imp = cjp->imp;
            if (imp != NULL) {
              if (imp->prepub == 2 || imp->pubstatus == PUBSTATUS_aheadofprint) {
                hasUnpublished = TRUE;
                vnpcgp = vnp;
                authors = cap->authors;
                for (vnp2 = cap->title; vnp2 != NULL; vnp2 = vnp2->next) {
                  if (vnp2->choice != Cit_title_name) continue;
                  str = (CharPtr) vnp2->data.ptrvalue;
                  if (StringHasNoText (str)) continue;
                  title = str;
                }
                journal = GetBestJournal (cjp->title);
                dp = imp->date;
              }
            }
          }
        }
      }
    } else if (vnp->choice == PUB_Book || vnp->choice == PUB_Man) {
      return;
    }
  }

  if (! hasUnpublished) return;

  MemSet ((Pointer) &pr, 0, sizeof (PubRef));

  pr.unpubcount = cfp->unpubcount;

  pr.seqids = ValNodeCopyStr (NULL, 0, cfp->buf);

  pr.authors = GetAuthorMLNameList (authors);
  pr.titlewords = GetTitleWords (title);
  if (vnpcgp != NULL) {
    if (PubLabelUnique (vnpcgp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
      pr.uniquestr = StringSaveNoNull (buf);
    }
  }

  tmp = StringSave (title);
  TrimSpacesAndJunkFromEnds (tmp, TRUE);
  s_RemovePeriodFromEnd (tmp);
  StrStripBrackets (tmp);
  StrStripSpaces (tmp);
  TrimSpacesAroundString (tmp);
  pr.fulltitle = tmp;

  pr.journal = journal;
  pr.imp = NULL;

  if (dp != NULL && dp->data [0] == 1) {
    year = (Int2) dp->data [1] + 1900;
  }
  if (year == 0) {
    year = cfp->year;
  }
  pr.year = year;
  pr.pmid = 0;

  if (pr.authors != NULL && pr.titlewords != NULL) {
    TryEntrezQueries (cfp, &pr);
  }

  ValNodeFreeData (pr.seqids);
  ValNodeFreeData (pr.authors);
  ValNodeFreeData (pr.titlewords);
  MemFree (pr.fulltitle);
  MemFree (pr.uniquestr);
}

static void DoUnpublishedReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
   if (sep == NULL || cfp == NULL) return;

  cfp->unpubcount = 0;
  cfp->unpublist = NULL;
  cfp->lastunpub = NULL;
  VisitPubdescsInSep (sep, (Pointer) cfp, CountUnpubPub);
  VisitPubdescsInSep (sep, (Pointer) cfp, ProcessUnpubPub);
}
//LCOV_EXCL_STOP

static void CreateUnpubList (
  PubdescPtr pdp,
  Pointer userdata
)

{
  AuthComp      authcomp;
  AuthListPtr   authors = NULL;
  Char          buf [521];
  CitArtPtr     cap = NULL;
  CleanFlagPtr  cfp;
  CitGenPtr     cgp = NULL;
  CitJourPtr    cjp;
  PubRefPtr     curr, prp;
  DatePtr       dp = NULL;
  Boolean       hasUnpublished = FALSE;
  ImprintPtr    imp;
  CharPtr       journal = NULL;
  CharPtr       str, tmp;
  CharPtr       title = NULL;
  Boolean       titlsame, uniqsame;
  ValNodePtr    vnp, vnpx, vnpy, vnpcgp = NULL;
  Int2          year = 0;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            hasUnpublished = TRUE;
            vnpcgp = vnp;
            authors = cgp->authors;
            title = cgp->title;
            journal = GetBestJournal (cgp->journal);
            dp = cgp->date;
          }
        }
      }
    } else if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
      if (cap != NULL) {
        if (cap->from == 1) {
          cjp = (CitJourPtr) cap->fromptr;
          if (cjp != NULL) {
            imp = cjp->imp;
            if (imp != NULL) {
              if (imp->prepub == 2 || imp->pubstatus == PUBSTATUS_aheadofprint) {
                hasUnpublished = TRUE;
                vnpcgp = vnp;
                authors = cap->authors;
                for (vnpy = cap->title; vnpy != NULL; vnpy = vnpy->next) {
                  if (vnpy->choice != Cit_title_name) continue;
                  str = (CharPtr) vnpy->data.ptrvalue;
                  if (StringHasNoText (str)) continue;
                  title = str;
                }
                journal = GetBestJournal (cjp->title);
                dp = imp->date;
              }
            }
          }
        }
      }
    }
  }

  if (! hasUnpublished) return;

  prp = (PubRefPtr) MemNew (sizeof (PubRef));
  if (prp == NULL) return;

  prp->unpubcount = cfp->unpubcount;

  prp->seqids = ValNodeCopyStr (NULL, 0, cfp->buf);

  prp->authors = GetAuthorMLNameList (authors);
  prp->titlewords = GetTitleWords (title);
  if (vnpcgp != NULL) {
    if (PubLabelUnique (vnpcgp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
      prp->uniquestr = StringSaveNoNull (buf);
    }
  }

  tmp = StringSave (title);
  TrimSpacesAndJunkFromEnds (tmp, TRUE);
  s_RemovePeriodFromEnd (tmp);
  StrStripBrackets (tmp);
  StrStripSpaces (tmp);
  TrimSpacesAroundString (tmp);
  prp->fulltitle = tmp;

  prp->journal = journal;
  prp->imp = NULL;

  if (dp != NULL && dp->data [0] == 1) {
    year = (Int2) dp->data [1] + 1900;
  }
  if (year == 0) {
    year = cfp->year;
  }
  prp->year = year;
  prp->pmid = 0;

  for (vnp = cfp->unpublist; vnp != NULL; vnp = vnp->next) {
    curr = (PubRefPtr) vnp->data.ptrvalue;
    if (curr == NULL) continue;

    authcomp = AuthorsIdentical (prp->authors, curr->authors);
    titlsame = (Boolean) (StringICmp (prp->fulltitle, curr->fulltitle) == 0);
    uniqsame = (Boolean) (StringICmp (prp->uniquestr, curr->uniquestr) == 0);

    if (authcomp == FULL_NAME_MATCH && titlsame && uniqsame) {
      for (vnpx = curr->seqids; vnpx != NULL; vnpx = vnpx->next) {
        str = (CharPtr) vnpx->data.ptrvalue;
        if (str == NULL) continue;
        if (StringCmp (str, cfp->buf) == 0) {
          ValNodeFreeData (prp->seqids);
          ValNodeFreeData (prp->authors);
          ValNodeFreeData (prp->titlewords);
          MemFree (prp->fulltitle);
          MemFree (prp->uniquestr);
          MemFree (prp);
          return;
        }
      }

      ValNodeCopyStr (&(curr->seqids), 0, cfp->buf);
      ValNodeFreeData (prp->seqids);
      ValNodeFreeData (prp->authors);
      ValNodeFreeData (prp->titlewords);
      MemFree (prp->fulltitle);
      MemFree (prp->uniquestr);
      MemFree (prp);
      return;
    }
  }

  vnp = ValNodeAddPointer (&(cfp->lastunpub), 0, (Pointer) prp);
  if (cfp->unpublist == NULL) {
    cfp->unpublist = vnp;
  }
  cfp->lastunpub = vnp;
}

// Used in JournalIsPresent
static Int4 GetNumericValue (
  CharPtr str
)

{
  long int  val;

  if (StringHasNoText (str)) return 0;
  if (! StringIsAllDigits (str)) return 0;

  if (sscanf (str, "%ld", &val) != 1) return 0;

  return val;
}


static void CheckPubmedSI (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr       accn = NULL;
  CleanFlagPtr  cfp;
  Int4          pmid;
  SeqIdPtr      sip;
  CharPtr       str;
  TextSeqIdPtr  tsip;
  XmlObjPtr     xop, tmp;

  if (bsp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_TPG:
      case SEQID_TPE:
      case SEQID_TPD:
      case SEQID_OTHER :
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip == NULL) continue;
        accn = tsip->accession;
        break;
      default :
        break;
    }
  }

  if (accn == NULL) return;

  str = DoEsearchQuery ("pubmed", accn, "[SI]", NULL, 0, 10, TRUE);
  if (str == NULL) return;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return;

  for (tmp = xop; tmp != NULL; tmp = tmp->successor) {
    if (XmlPathSuffixIs (tmp, "/IdList/Id")) {
      pmid = GetNumericValue (tmp->contents);
      if (pmid > 0) {
        if (cfp->logfp != NULL) {
          fprintf (cfp->logfp, "PMID %ld\tSEQID %s\n", (long) pmid, cfp->buf);
        }
      }
    }
  }

  FreeXmlObject (xop);
}

static void DoUnpublishedList (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  if (sep == NULL || cfp == NULL) return;

  cfp->unpubcount = 0;
  VisitPubdescsInSep (sep, (Pointer) cfp, CountUnpubPub);
  VisitPubdescsInSep (sep, (Pointer) cfp, CreateUnpubList);
  /*
  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "Trying ESearch Accession [SI]\n\n");
  }
  VisitBioseqsInSep (sep, (Pointer) cfp, CheckPubmedSI);
  */
}

static void FinishUnpublishedList (
  CleanFlagPtr cfp
)

{
  PubRefPtr   prp;
  ValNodePtr  vnp;

  if (cfp == NULL) return;

  for (vnp = cfp->unpublist; vnp != NULL; vnp = vnp->next) {
    prp = (PubRefPtr) vnp->data.ptrvalue;
    if (prp == NULL) continue;
    if (prp->authors != NULL && prp->titlewords != NULL) {
      TryEntrezQueries (cfp, prp);
    }
  }

  for (vnp = cfp->unpublist; vnp != NULL; vnp = vnp->next) {
    prp = (PubRefPtr) vnp->data.ptrvalue;
    if (prp == NULL) continue;
    ValNodeFreeData (prp->seqids);
    ValNodeFreeData (prp->authors);
    ValNodeFreeData (prp->titlewords);
    MemFree (prp->fulltitle);
    MemFree (prp->uniquestr);
    MemFree (prp);
  }
  ValNodeFree (cfp->unpublist);
}

static void CountPublishedPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CitArtPtr     cap = NULL;
  ValNodePtr    vnp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
    } else if (vnp->choice == PUB_PMid) {
      return;
    }
  }

  if (cap == NULL) return;

  (cfp->unpubcount)++;
}

static void ProcessPublishedPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CitArtPtr     cap = NULL;
  CitJourPtr    cjp;
  DatePtr       dp = NULL;
  ImprintPtr    imp;
  PubRef        pr;
  CharPtr       str;
  CharPtr       tmp;
  ValNodePtr    vnp;
  Int2          year = 0;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
    } else if (vnp->choice == PUB_PMid) {
      return;
    }
  }

  if (cap == NULL) return;

  MemSet ((Pointer) &pr, 0, sizeof (PubRef));

  pr.seqids = ValNodeCopyStr (NULL, 0, cfp->buf);

  pr.authors = GetAuthorMLNameList (cap->authors);
  for (vnp = cap->title; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Cit_title_name) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      pr.titlewords = GetTitleWords (str);
      tmp = StringSave (str);
      TrimSpacesAndJunkFromEnds (tmp, TRUE);
      s_RemovePeriodFromEnd (tmp);
      StrStripBrackets (tmp);
      StrStripSpaces (tmp);
      TrimSpacesAroundString (tmp);
      pr.fulltitle = tmp;
    }
  }

  if (cap->from == 1) {
    cjp = (CitJourPtr) cap->fromptr;
    if (cjp != NULL) {
      pr.journal = GetBestJournal (cjp->title);
      imp = cjp->imp;
      pr.imp = imp;
      if (imp != NULL) {
        dp = imp->date;
        if (dp != NULL && dp->data [0] == 1) {
          year = (Int2) dp->data [1] + 1900;
        }
      }
    }
  }

  if (year == 0) {
    year = cfp->year;
  }
  pr.year = year;
  pr.pmid = 0;

  if (pr.authors != NULL && pr.titlewords != NULL) {
    TryEntrezQueries (cfp, &pr);
  }

  ValNodeFreeData (pr.seqids);
  ValNodeFreeData (pr.authors);
  ValNodeFreeData (pr.titlewords);
  MemFree (pr.fulltitle);
}

static void DoPublishedReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  if (sep == NULL || cfp == NULL) return;

  cfp->unpubcount = 0;
  cfp->unpublist = NULL;
  cfp->lastunpub = NULL;
  VisitPubdescsInSep (sep, (Pointer) cfp, CountPublishedPub);
  VisitPubdescsInSep (sep, (Pointer) cfp, ProcessPublishedPub);
}

// Used in ProcessMissingJournal
static Boolean JournalIsPresent (
  CharPtr title
)

{
  Int4       count = 0;
  Boolean    is_issn = FALSE;
  CharPtr    str = NULL;
  XmlObjPtr  tmp, xop;

  if (StringHasNoText (title)) return FALSE;

  ConvertToTermListForm (title);

  if (LooksLikeISSN (title)) {
    is_issn = TRUE;
    if (title [4] == '+' || title [4] == ' ') {
      title [4] = '-';
    }
    if (title [8] == 'x') {
      title [8] = 'X';
    }
    str = DoEsearchQuery ("nlmcatalog", title, "[issn]", NULL, 0, 10, FALSE);
  }

  if (str == NULL) {
    str = DoEsearchQuery ("nlmcatalog", title, "[multi] AND ncbijournals[sb]", NULL, 0, 10, FALSE);
  }

  if (str == NULL) return FALSE;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return FALSE;

  for (tmp = xop; tmp != NULL; tmp = tmp->successor) {
    if (XmlPathSuffixIs (tmp, "/eSearchResult/Count")) {
      count = GetNumericValue (tmp->contents);
    }
  }

  FreeXmlObject (xop);

  if (count > 0) return TRUE;

  str = DoEsearchQuery ("nlmcatalog", title, "[jour]", NULL, 0, 10, FALSE);

  if (str == NULL) return FALSE;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return FALSE;

  for (tmp = xop; tmp != NULL; tmp = tmp->successor) {
    if (XmlPathSuffixIs (tmp, "/eSearchResult/Count")) {
      count = GetNumericValue (tmp->contents);
    }
  }

  FreeXmlObject (xop);

  if (count > 0) return TRUE;

  return FALSE;
}

// Used in DoJournalReport
static void ProcessMissingJournal (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CitArtPtr     cap;
  CitJourPtr    cjp;
  CharPtr       journal;
  ValNodePtr    vnp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice != PUB_Article) continue;
    cap = (CitArtPtr) vnp->data.ptrvalue;
    if (cap == NULL || cap->from != 1) continue;
    cjp = (CitJourPtr) cap->fromptr;
    if (cjp == NULL) continue;
    journal = GetBestJournal (cjp->title);
    if (! JournalIsPresent (journal)) {
      fprintf (cfp->logfp, "%s\t%s\n", cfp->buf, journal);
    }
  }
}

static void DoJournalReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  if (sep == NULL || cfp == NULL) return;

  VisitPubdescsInSep (sep, (Pointer) cfp, ProcessMissingJournal);
}

static void RemoveFeatureCitations (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL || sfp->cit == NULL) return;

  sfp->cit = PubSetFree (sfp->cit);
}

static SeqIdPtr SeqIdFindBestForPromotion (SeqIdPtr sip)

{
  return SeqIdFindBest (sip, 0);
}

static SeqIdPtr SeqIdFindWorstForPromotion (SeqIdPtr sip)

{
  return SeqIdFindWorst (sip);
}

static void PromoteSeqId (SeqIdPtr sip, Boolean alsoCheckLocalAccn, Boolean findWorst)

{
  SeqIdPtr     bestid, newid, oldid;
  BioseqPtr    bsp;
  ObjectIdPtr  oip;
  TextSeqId    tsi;
  SeqId        vn;

  bsp = BioseqFind (sip);
  if (bsp == NULL && alsoCheckLocalAccn && sip->choice == SEQID_LOCAL) {
    oip = (ObjectIdPtr) sip->data.ptrvalue;
    if (oip != NULL && (! StringHasNoText (oip->str))) {
      MemSet ((Pointer) &vn, 0, sizeof (SeqId));
      MemSet ((Pointer) &tsi, 0, sizeof (TextSeqId));
      tsi.accession = oip->str;
      vn.choice = SEQID_GENBANK;
      vn.data.ptrvalue = (Pointer) &tsi;
      bsp = BioseqFind (&vn);
    }
  }
  if (bsp == NULL) return;

  if (findWorst) {
    bestid = SeqIdFindWorstForPromotion (bsp->id);
  } else {
    bestid = SeqIdFindBestForPromotion (bsp->id);
  }
  if (bestid == NULL) return;
  newid = SeqIdDup (bestid);
  if (newid == NULL) return;

  oldid = ValNodeNew (NULL);
  if (oldid == NULL) return;

  MemCopy (oldid, sip, sizeof (ValNode));
  oldid->next = NULL;

  sip->choice = newid->choice;
  sip->data.ptrvalue = newid->data.ptrvalue;

  SeqIdFree (oldid);
  ValNodeFree (newid);

  SeqIdStripLocus (sip);
}

static void PromoteSeqIdList (SeqIdPtr sip, Boolean alsoCheckLocalAccn, Boolean findWorst)

{
  while (sip != NULL) {
    PromoteSeqId (sip, alsoCheckLocalAccn, findWorst);
    sip = sip->next;
  }
}

static void PromoteSeqLocList (SeqLocPtr slp, Boolean alsoCheckLocalAccn, Boolean findWorst)

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
        PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          PromoteSeqLocList (loc, alsoCheckLocalAccn, findWorst);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            PromoteSeqIdList (sip, alsoCheckLocalAccn, findWorst);
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

static Boolean PromoteIDsProc (GatherObjectPtr gop, Boolean findWorst)

{
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  RnaRefPtr     rrp;
  SeqFeatPtr    sfp;
  tRNAPtr       trp;

  if (gop->itemtype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gop->dataptr;
  if (sfp == NULL) return TRUE;

  PromoteSeqLocList (sfp->location, FALSE, findWorst);

  PromoteSeqLocList (sfp->product, FALSE, findWorst);

  switch (sfp->data.choice) {
    case SEQFEAT_CDREGION :
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL && crp->code_break != NULL) {
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          PromoteSeqLocList (cbp->loc, FALSE, findWorst);
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->type == 3 && rrp->ext.choice == 2) {
        trp = rrp->ext.value.ptrvalue;
        if (trp != NULL && trp->anticodon != NULL) {
          PromoteSeqLocList (trp->anticodon, FALSE, findWorst);
        }
      }
      break;
    default :
      break;
  }

  return TRUE;
}

static Boolean PromoteWorstIDsProc (GatherObjectPtr gop)

{
  return PromoteIDsProc (gop, TRUE);
}

static Boolean RemoveGIProc (GatherObjectPtr gop)

{
  BioseqPtr      bsp;
  SeqIdPtr       nextid, sip;
  SeqIdPtr PNTR  previd;

  if (gop->itemtype != OBJ_BIOSEQ) return TRUE;
  bsp = (BioseqPtr) gop->dataptr;
  if (bsp == NULL) return TRUE;

  previd = (SeqIdPtr PNTR) &(bsp->id);
  sip = bsp->id;
  while (sip != NULL) {
    nextid = sip->next;
    if (sip->choice == SEQID_GI) {
      *previd = sip->next;
      sip->next = NULL;
      SeqIdFree (sip);
    } else {
      previd = (SeqIdPtr PNTR) &(sip->next);
    }
    sip = sip->next;
  }

  return TRUE;
}

static void StripLocusFromSeqId (SeqIdPtr sip, Pointer userdata)

{
  TextSeqIdPtr  tsip;

  if (sip == NULL) return;

  switch (sip->choice) {
    case SEQID_GENBANK :
    case SEQID_TPG :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip == NULL) return;
      if (tsip->name == NULL) return;
      tsip->name = MemFree (tsip->name);
      break;
    default :
      break;
  }
}

static void StripLocusFromBsp (BioseqPtr bsp, Pointer userdata)

{
  if (bsp == NULL) return;

  VisitSeqIdsInBioseq (bsp, NULL, StripLocusFromSeqId);
}

// Used for -K i
static void StripVersionsFromSeqId (SeqIdPtr sip, Pointer userdata)

{
  TextSeqIdPtr  tsip;

  if (sip == NULL) return;

  switch (sip->choice) {
    case SEQID_GENBANK :
    case SEQID_TPG :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip == NULL) return;
      tsip->version = INT2_MIN;
      break;
    default :
      break;
  }
}

static void StripVersionsFromBsp (BioseqPtr bsp, Pointer userdata)

{
  if (bsp == NULL) return;

  VisitSeqIdsInBioseq (bsp, NULL, StripVersionsFromSeqId);
}

// Used for -K i
static void StripVersionsFromSfp (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp == NULL) return;

  VisitSeqIdsInSeqFeat (sfp, NULL, StripVersionsFromSeqId);
}

// Used for -K i
static void PromoteToWorstIDProc (SeqEntryPtr sep, Uint2 entityID)

{
  SeqEntryPtr  oldscope;

  if (sep == NULL || entityID < 1) return;

  oldscope = SeqEntrySetScope (sep);

  GatherObjectsInEntity (entityID, 0, NULL, PromoteWorstIDsProc, NULL, NULL);
  GatherObjectsInEntity (entityID, 0, NULL, RemoveGIProc, NULL, NULL);

  VisitBioseqsInSep (sep, NULL, StripVersionsFromBsp);
  VisitBioseqsInSep (sep, NULL, StripLocusFromBsp);

  VisitFeaturesInSep (sep, NULL, StripVersionsFromSfp);

  SeqEntrySetScope (oldscope);
}

// Used for -K r
static void ReassignLocalSIP (SeqIdPtr sip, Pointer userdata)

{
  Int4Ptr      countp;
  ObjectIdPtr  oip;

  if (sip == NULL || sip->choice != SEQID_LOCAL) return;
  countp = (Int4Ptr) userdata;
  if (countp == NULL) return;

  oip = (ObjectIdPtr) sip->data.ptrvalue;
  if (oip == NULL) return;

  if (oip->str != NULL) {
    oip->str = MemFree (oip->str);
  }

  oip->id = *countp;
}

// Used for -K r
static void ReassignLocalSFP (SeqFeatPtr sfp, Pointer userdata)

{
  if (sfp == NULL || userdata == NULL) return;

  VisitSeqIdsInSeqLoc (sfp->location, userdata, ReassignLocalSIP);
  VisitSeqIdsInSeqLoc (sfp->product, userdata, ReassignLocalSIP);
}

// Used for -K r
static void ReassignLocalBSP (BioseqPtr bsp, Pointer userdata)

{
  Int4Ptr  countp;

  if (bsp == NULL || userdata == NULL) return;
  countp = (Int4Ptr) userdata;
  (*countp)++;

  VisitSeqIdsInBioseq (bsp, userdata, ReassignLocalSIP);
  VisitFeaturesOnBsp (bsp, userdata, ReassignLocalSFP);
}

// Used for -K r
static void ReassignLocalIDProc (SeqEntryPtr sep, Uint2 entityID)

{
  Int4         count = 0;
  SeqEntryPtr  oldscope;

  if (sep == NULL || entityID < 1) return;

  oldscope = SeqEntrySetScope (sep);

  VisitBioseqsInSep (sep, (Pointer) &count, ReassignLocalBSP);

  SeqEntrySetScope (oldscope);
}

static void MakeMrnaFromCDS (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr  mrna;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  mrna = AddmRNAForCDS (sfp);
  if (mrna == NULL) return;

  sfp->idx.scratch = (Pointer) mrna;
}

// Used for -C c
static void LinkMrnaToCDS (SeqFeatPtr sfp, Pointer userdata)

{
  SeqFeatPtr  mrna;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return;
  mrna = sfp->idx.scratch;
  if (mrna == NULL) return;
  sfp->idx.scratch = NULL;

  LinkTwoFeatures (sfp, mrna);
  LinkTwoFeatures (mrna, sfp);
}

// Used for -C c
static void CopyCDStomRNA (SeqEntryPtr sep)

{
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, MakeMrnaFromCDS);
  AssignFeatureIDs (sep);
  VisitFeaturesInSep (sep, NULL, LinkMrnaToCDS);
}

typedef struct recipdata {
  CharPtr     orig_trans_id;
  CharPtr     orig_prot_id;
  SeqFeatPtr  mrna;
  SeqFeatPtr  cds;
} RecipData, PNTR RecipDataPtr;

typedef struct reciplist {
  RecipDataPtr PNTR  array;
  Int4               num;
} RecipList, PNTR RecipListPtr;

// Used for -X f
static void CollectRecipIDs (SeqFeatPtr sfp, Pointer userdata)

{
  GBQualPtr        gbq;
  RecipDataPtr     rdp;
  ValNodeBlockPtr  vbp;

  if (sfp == NULL) return;
  if (sfp->idx.subtype != FEATDEF_mRNA) return;
  vbp = (ValNodeBlockPtr) userdata;
  if (vbp == NULL) return;

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "transcript_id") != 0 && StringICmp (gbq->qual, "orig_transcript_id") != 0) continue;
    rdp = (RecipDataPtr) MemNew (sizeof (RecipData));
    if (rdp == NULL) continue;
    rdp->orig_trans_id = StringSave (gbq->val);
    rdp->mrna = sfp;
    ValNodeAddPointerEx (&(vbp->head), &(vbp->tail), 0, (Pointer) rdp);
  }
}

// Used for -X f
static int LIBCALLBACK SortRecipIDs (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  int           cmp;
  RecipDataPtr  rdp1;
  RecipDataPtr  rdp2;
  ValNodePtr    vnp1;
  ValNodePtr    vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  rdp1 = (RecipDataPtr) vnp1->data.ptrvalue;
  rdp2 = (RecipDataPtr) vnp2->data.ptrvalue;
  if (rdp1 == NULL || rdp2 == NULL) return 0;

  cmp = StringICmp (rdp1->orig_trans_id, rdp2->orig_trans_id);
  if (cmp > 0) {
    return 1;
  } else if (cmp < 0) {
    return -1;
  }

  return 0;
}

// Used for -X f
static void FillInRecipIDs (SeqFeatPtr sfp, Pointer userdata)

{
  RecipDataPtr PNTR  array;
  BioseqPtr          bsp;
  DbtagPtr           dbt;
  GBQualPtr          gbq;
  Char               id [128];
  Int2               L, R, mid;
  RecipDataPtr       rdp;
  RecipListPtr       rlp;
  SeqIdPtr           sip;
  CharPtr            str;

  if (sfp == NULL) return;
  if (sfp->idx.subtype != FEATDEF_CDS) return;
  rlp = (RecipListPtr) userdata;
  if (rlp == NULL) return;

  array = rlp->array;
  if (array == NULL) return;

  str = NULL;
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "transcript_id") != 0 && StringICmp (gbq->qual, "orig_transcript_id") != 0) continue;
    str = gbq->val;
  }
  if (StringHasNoText (str)) return;

  L = 0;
  R = rlp->num - 1;

  while (L < R) {
    mid = (L + R) / 2;
    rdp = array [mid];
    if (rdp != NULL && StringICmp (rdp->orig_trans_id, str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  rdp = array [R];
  if (rdp != NULL && StringICmp (rdp->orig_trans_id, str) == 0 && rdp->cds == NULL && rdp->orig_prot_id == NULL) {
    rdp->cds = sfp;
    if (sfp->product != NULL) {
      bsp = BioseqFindFromSeqLoc (sfp->product);
      if (bsp != NULL) {
        id [0] = '\0';
        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          if (sip->choice != SEQID_GENERAL) continue;
          dbt = (DbtagPtr) sip->data.ptrvalue;
          if (dbt == NULL) continue;
          if (IsSkippableDbtag (dbt)) continue;
          SeqIdWrite (sip, id, PRINTID_FASTA_SHORT, sizeof (id) - 1);
        }
        if (StringDoesHaveText (id)) {
          rdp->orig_prot_id = StringSave (id);
        }
      }
    } else {
      str = NULL;
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringICmp (gbq->qual, "protein_id") != 0 && StringICmp (gbq->qual, "orig_protein_id") != 0) continue;
        str = gbq->val;
      }
      if (StringDoesHaveText (str)) {
        rdp->orig_prot_id = StringSave (str);
      }
    }
  }
}

// Used for -X f
static void AddMissingGbQual (SeqFeatPtr sfp, CharPtr qual, CharPtr val)

{
  GBQualPtr  gbq;

  if (sfp == NULL || StringHasNoText (qual) || StringHasNoText (val)) return;

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, qual) == 0) return;
  }

  gbq = GBQualNew ();
  if (gbq == NULL) return;
  gbq->qual = StringSave (qual);
  gbq->val = StringSave (val);

  gbq->next = sfp->qual;
  sfp->qual = gbq;
}

// Used for -X f
static void FixReciprocalProtTransIDs (SeqEntryPtr sep)

{
  RecipDataPtr PNTR  array;
  Int4               i;
  Int4               num;
  RecipDataPtr       rdp;
  RecipList          rl;
  ValNodeBlock       vb;
  ValNodePtr         vnp;

  if (sep == NULL) return;

  MemSet ((Pointer) &vb, 0, sizeof (ValNodeBlock));
  VisitFeaturesInSep (sep, (Pointer) &vb, CollectRecipIDs);
  if (vb.head == NULL) return;

  vb.head = ValNodeSort (vb.head, SortRecipIDs);
  num = ValNodeLen (vb.head);
  array = (RecipDataPtr PNTR) MemNew (sizeof (RecipDataPtr) * (num + 1));
  if (array == NULL) return;

  for (i = 0, vnp = vb.head; i < num && vnp != NULL; i++, vnp = vnp->next) {
    rdp = (RecipDataPtr) vnp->data.ptrvalue;
    array [i] = rdp;
  }

  MemSet ((Pointer) &rl, 0, sizeof (RecipList));
  rl.array = array;
  rl.num = num;

  VisitFeaturesInSep (sep, (Pointer) &rl, FillInRecipIDs);

  for (vnp = vb.head; vnp != NULL; vnp = vnp->next) {
    rdp = (RecipDataPtr) vnp->data.ptrvalue;
    if (rdp == NULL) continue;
    AddMissingGbQual (rdp->mrna, "orig_transcript_id", rdp->orig_trans_id);
    AddMissingGbQual (rdp->mrna, "orig_protein_id", rdp->orig_prot_id);
    AddMissingGbQual (rdp->cds, "orig_transcript_id", rdp->orig_trans_id);
  }

  for (vnp = vb.head; vnp != NULL; vnp = vnp->next) {
    rdp = (RecipDataPtr) vnp->data.ptrvalue;
    if (rdp == NULL) continue;
    MemFree (rdp->orig_trans_id);
    MemFree (rdp->orig_prot_id);
  }

  ValNodeFreeData (vb.head);
}

//LCOV_EXCL_START
static void DoReversePreOrNcRNATransIDs(SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr  bsp;
  Char       buf [256];
  GBQualPtr  gbq;

  if (sfp == NULL || sfp->product == NULL) return;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->idx.subtype != FEATDEF_preRNA && sfp->idx.subtype != FEATDEF_ncRNA) return;

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL) return;

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "transcript_id") == 0) return;
  }

  SeqIdWrite (bsp->id, buf, PRINTID_FASTA_ALL, sizeof (buf) - 1);
  if (StringHasNoText (buf)) return;

  gbq = GBQualNew ();
  if (gbq == NULL) return;

  gbq->qual = StringSave ("transcript_id");
  gbq->val = StringSave (buf);
  gbq->next = sfp->qual;
  sfp->qual = gbq->next;

  sfp->product = SeqLocFree (sfp->product);

  bsp->idx.deleteme = TRUE;
}

static void ReversePreOrNcRNATransIDs(SeqEntryPtr sep)

{
  if (sep == NULL) return;

  VisitFeaturesInSep (sep, NULL, DoReversePreOrNcRNATransIDs);
}

// Unused
static void SeqToAnticodon (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  GBQualPtr     acgbq = NULL;
  Char          anticodon [16];
  Char          buf5 [16], buf3 [16], buf4 [16];
  Char          ch;
  GBQualPtr     gbq;
  Int2          gc;
  Int2          max;
  long          pos;
  CharPtr       ptr;
  RnaRefPtr     rrp;
  SeqIntPtr     sintp;
  SeqLocPtr     slp;
  CharPtr       str;
  tRNAPtr       trna;

  if (sfp == NULL || sfp->location == NULL) return;

  if (sfp->data.choice != SEQFEAT_RNA) return;
  if (sfp->idx.subtype != FEATDEF_tRNA) return;

  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL) return;
  if (rrp->type != 3 || rrp->ext.choice != 2) return;
  trna = (tRNAPtr) rrp->ext.value.ptrvalue;
  if (trna == NULL) return;

  if (trna->anticodon != NULL) return;

  anticodon [0] = '\0';
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next)  {
    if (StringICmp (gbq->qual, "anticodon") != 0) continue;
    if (StringLen (gbq->val) > 3) continue;
    acgbq = gbq;
    StringCpy (anticodon, gbq->val);
    ptr = anticodon;
    ch = *ptr;
    while (ch != '\0') {
      ch = TO_UPPER (ch);
      if (ch == 'U') {
        ch = 'T';
      }
      *ptr = ch;
      ptr++;
      ch = *ptr;
    }
  }

  if (StringHasNoText (anticodon)) return;

  str = GetSequenceByFeature (sfp);
  if (str == NULL) return;

  ptr = str;
  while (ptr != NULL && trna->anticodon == NULL) {
    ptr = StringStr (ptr, anticodon);
    if (ptr != NULL) {
      pos = (long) (ptr + 1 - str);
      ExtractString (str, buf5, pos - 10, pos - 3);
      ExtractString (str, buf3, pos + 4, pos + 11);
      if (StringHasNoText (buf5)) {
        StringCpy (buf5, "       ");
      }
      if (StringHasNoText (buf3)) {
        StringCpy (buf3, "       ");
      }
      StringCpy (buf4, "       ");

      max = 0;
      gc = 0;
      if (StringDoesHaveText (buf5) && StringDoesHaveText (buf3) && pos > 30 && pos < 40) {
        StringCpy (buf4, buf3);
        RevcompString (buf4);
        max = LongestMatch (buf5, buf4, &gc);
        if (max > 4 || (max == 4 && gc > 1)) {
          slp = sfp->location;
          if (slp->choice == SEQLOC_INT) {
            slp = SeqLocCopy (sfp->location);
            if (slp != NULL && slp->choice == SEQLOC_INT) {
              sintp = (SeqIntPtr) slp->data.ptrvalue;
              if (sintp != NULL) {
                sintp->if_from = IntFuzzFree (sintp->if_from);
                sintp->if_to = IntFuzzFree (sintp->if_to);
                if (sintp->strand == Seq_strand_minus) {
                  sintp->to -= pos - 1;
                  sintp->from = sintp->to - 2;
                } else {
                  sintp->from += pos - 1;
                  sintp->to = sintp->from + 2;
                }
                trna->anticodon = slp;
                if (acgbq != NULL) {
                  acgbq->qual = MemFree (acgbq->qual);
                  acgbq->val = MemFree (acgbq->val);
                }
              }
            }
          }
        }
      }

      ptr++;
    }
  }

  MemFree (str);
}

#ifdef OS_UNIX
static SeqEntryPtr CppBasicCleanup (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  AsnIoPtr      aip, aop;
  ByteStorePtr  bs1, bs2;
  Char          cmmd [512];
  SeqEntryPtr   csep, nsep;
  Char          path1 [PATH_MAX];
  Char          path2 [PATH_MAX];
  Char          path3 [PATH_MAX];

  if (sep == NULL || cfp == NULL) return NULL;

  VisitFeaturesInSep (sep, NULL, RemoveFeatureCitations);

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);

  aop = AsnIoOpen (path1, "w");
  SeqEntryAsnWrite (sep, aop, NULL);
  AsnIoClose (aop);

  sprintf (cmmd, "%s -i %s | cleanasn -a e -o %s",
           "~/ncbi_cxx/compilers/xCode/build/bin/Debug/test_basic_cleanup",
           path1, path2);
  system (cmmd);

  sprintf (cmmd, "cleanasn -i %s -o %s -K b",
           path1, path3);
  system (cmmd);

  aip = AsnIoOpen (path3, "r");
  csep = SeqEntryAsnRead (aip, NULL);
  AsnIoClose (aip);

  bs1 = Se2Bs (csep);

  aip = AsnIoOpen (path2, "r");
  nsep = SeqEntryAsnRead (aip, NULL);
  AsnIoClose (aip);

  bs2 = Se2Bs (nsep);

  if (nsep == NULL) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "EMPTY %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (! BSEqual (bs1, bs2)) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "BSEC DIFF %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
    if (cfp->gi > 0) {
      sprintf (cmmd, "echo '' >> ~/Desktop/diffclean.txt");
      system (cmmd);
      sprintf (cmmd, "echo '' >> ~/Desktop/diffclean.txt");
      system (cmmd);
      sprintf (cmmd, "echo '********** gi|%ld **********' >> ~/Desktop/diffclean.txt", (long) cfp->gi);
      system (cmmd);
      sprintf (cmmd, "echo '' >> ~/Desktop/diffclean.txt");
      system (cmmd);
      sprintf (cmmd, "diff %s %s >> ~/Desktop/diffclean.txt", path3, path2);
      system (cmmd);
    }
  }

  BSFree (bs1);
  BSFree (bs2);

  SeqEntryFree (csep);

  sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
  system (cmmd);

  return nsep;
}
#endif
//LCOV_EXCL_STOP

/* now only strips serials for local, general, refseq, and 2+6 genbank ids */
static void CheckForSwissProtIDX (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

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

//LCOV_EXCL_START
static void CheckForNucProt (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  BoolPtr  isNpsP;

  if (bssp == NULL || userdata == NULL) return;
  if (bssp->_class == BioseqseqSet_class_nuc_prot) {
    isNpsP = (BoolPtr) userdata;
    *isNpsP = TRUE;
  }
}
//LCOV_EXCL_STOP

//LCOV_EXCL_START
static void CheckForSegSeq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BoolPtr  isSegP;

  if (bsp == NULL || userdata == NULL) return;
  if (bsp->repr == Seq_repr_seg) {
    isSegP = (BoolPtr) userdata;
    *isSegP = TRUE;
  }
}
//LCOV_EXCL_STOP

//LCOV_EXCL_START
static void CheckForPPMEW (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  BoolPtr  isNpsP;

  if (bssp == NULL || userdata == NULL) return;
  if (bssp->_class == BioseqseqSet_class_mut_set ||
      bssp->_class == BioseqseqSet_class_pop_set ||
      bssp->_class == BioseqseqSet_class_phy_set ||
      bssp->_class == BioseqseqSet_class_eco_set ||
      bssp->_class == BioseqseqSet_class_wgs_set) {
    isNpsP = (BoolPtr) userdata;
    *isNpsP = TRUE;
  }
}
//LCOV_EXCL_STOP

//LCOV_EXCL_START
static void CheckForCDSs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BoolPtr  HasCdsP;

  if (sfp == NULL || userdata == NULL) return;
  if (sfp->data.choice != SEQFEAT_CDREGION) return;
  HasCdsP = (BoolPtr) userdata;
  *HasCdsP = TRUE;
}
//LCOV_EXCL_STOP

static SeqAlignPtr LIBCALLBACK GetSeqAlignEx (BioseqPtr bsp1, BioseqPtr bsp2, Boolean no_extend)
{
   BLAST_SummaryOptions *options = NULL;
   SeqAlignPtr           salp = NULL;
   
   if (bsp1 == NULL || bsp2 == NULL) return NULL;

   BLAST_SummaryOptionsInit(&options);
   if (bsp1->length > 10000 || bsp2->length > 10000 || no_extend)
   {
      options->filter_string = StringSave ("m L");
      if (no_extend) {
        options->word_size = 28;
      } else {
        options->word_size = 20;
        options->cutoff_evalue = act_get_eval (60);
      }
      options->hint = eNone;
   }
   else
   {
    options->filter_string = StringSave ("m F");
   }
   if (ISA_na (bsp1->mol))
   {
     options->program = eBlastn;
   }
   else
   {
     options->program = eBlastp;
   }

   BLAST_TwoSequencesSearch(options, bsp1, bsp2, &salp);
   BLAST_SummaryOptionsFree(options);
   return salp;
}


static SeqAlignPtr LIBCALLBACK GetSeqAlign (BioseqPtr bsp1, BioseqPtr bsp2)
{
  return GetSeqAlignEx (bsp1, bsp2, FALSE);
}


static SeqAlignPtr LIBCALLBACK GetSeqAlignPiece (SeqLocPtr slp1, SeqLocPtr slp2)
{
   BLAST_SummaryOptions *options = NULL;
   SBlastSeqalignArray * seqalign_arr=NULL;
   SeqAlignPtr           salp = NULL;
   BioseqPtr             bsp;
   
   if (slp1 == NULL || slp2 == NULL) return NULL;


   bsp = BioseqFindFromSeqLoc (slp1);
   if (bsp == NULL)
   {
     return NULL;
   }
   
   BLAST_SummaryOptionsInit(&options);
   
   if (ISA_na (bsp->mol))
   {
     options->program = eBlastn;
   }
   else
   {
     options->program = eBlastp;
   }

   options->gapped_calculation = TRUE;
   options->cutoff_evalue = 10;
   options->gap_x_dropoff = 100;
   options->gap_open = 4;
   options->gap_extend = 1;
   options->nucleotide_mismatch = -1; 
   options->word_size = 7;

   BLAST_TwoSeqLocSets(options, slp1, slp2, NULL, &seqalign_arr, NULL, NULL, NULL);
   
   if (seqalign_arr != NULL)
   {
     salp = seqalign_arr->array[0];
     seqalign_arr->array[0] = NULL;
     seqalign_arr = SBlastSeqalignArrayFree(seqalign_arr);
   }
    
   BLAST_SummaryOptionsFree(options);
   return salp;
}

static SeqAlignPtr Cln_GlobalAlign2Seq (BioseqPtr bsp1, BioseqPtr bsp2, BoolPtr revcomp)

{
  return Sqn_GlobalAlign2SeqEx (bsp1, bsp2, revcomp, GetSeqAlign, GetSeqAlignPiece, TRUE);
}


static void LIBCALLBACK ValidCallback (
  ErrSev severity,
  int errcode,
  int subcode,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  CharPtr accession,
  CharPtr seqid,
  CharPtr featureID,
  CharPtr message,
  CharPtr objtype,
  CharPtr label,
  CharPtr context,
  CharPtr location,
  CharPtr product,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  CharPtr            level;
  SeqFeatPtr         sfp;

  if (userdata == NULL) return;
  level = (CharPtr) userdata;

  if (severity == SEV_REJECT && StringChr (level, 'r') == NULL) return;
  if (severity == SEV_ERROR && StringChr (level, 'e') == NULL) return;
  if (severity == SEV_WARNING && StringChr (level, 'w') == NULL) return;
  if (severity == SEV_INFO && StringChr (level, 'i') == NULL) return;

  switch (itemtype) {
    case OBJ_SEQFEAT :
      sfp = SeqMgrGetDesiredFeature (entityID, NULL, itemID, 0, NULL, &fcontext);
      if (sfp == NULL) return;
      sfp->idx.deleteme = TRUE;
      if (sfp->product == NULL) return;
      bsp = BioseqFindFromSeqLoc (sfp->product);
      if (bsp == NULL) return;
      bsp->idx.deleteme = TRUE;
      break;
    default :
      break;
  }
}

static void CDSPartialIfFrame (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CdRegionPtr  crp;
  Boolean      partial5 = FALSE;
  Boolean      partial3 = FALSE;

  if (sfp == NULL) return;
  if (sfp->data.choice != SEQFEAT_CDREGION) return;

  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return;

  if (crp->frame < 2) return;

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
  SetSeqLocPartial (sfp->location, TRUE, partial3);
  sfp->partial = TRUE;
}

static void AddFeatureToBioseq (SeqFeatPtr sfp, BioseqPtr bsp)

{
  SeqFeatPtr   prev;
  SeqAnnotPtr  sap;

  if (sfp == NULL || bsp == NULL) return;
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
  }
  sap = bsp->annot;
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

static void PackageFeatureCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp = NULL;
  BioseqSetPtr   bssp = NULL;
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  BioseqPtr      target;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        target = BioseqFindFromSeqLoc (sfp->location);
        if (target != bsp) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          AddFeatureToBioseq (sfp, target);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static void RemoveLocalIDs (BioseqPtr bsp, Pointer userdata)

{
  SeqIdPtr      nextsip;
  Pointer PNTR  prevsip;
  Boolean       replaced;
  SeqIdPtr      sip;

  if (bsp == NULL) return;
  sip = bsp->id;
  prevsip = (Pointer PNTR) &(bsp->id);
  replaced = FALSE;
  while (sip != NULL) {
    nextsip = sip->next;
    if (sip->choice == SEQID_LOCAL) {
      (*prevsip) = sip->next;
      sip->next = NULL;
      SeqIdFree (sip);
      replaced = TRUE;
    } else {
      prevsip = (Pointer PNTR) &(sip->next);
    }
    sip = nextsip;
  }
  if (replaced) {
    SeqMgrReplaceInBioseqIndex (bsp);
  }
}

static void ZapBadTitles (SeqDescrPtr sdp, Pointer userdata)

{
  ObjValNodePtr  ovp;
  CharPtr        str;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  str = (CharPtr) sdp->data.ptrvalue;
  if (str == NULL) return;
  if (*str != '[' && StringNICmp (str, "hypothetical protein", 20) != 0) return;
  if (sdp->extended == 0) return;
  ovp = (ObjValNodePtr) sdp;
  ovp->idx.deleteme = TRUE;
}

//LCOV_EXCL_START
//only for segsets
static void FindSegSet (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqSetPtr      bssp;
  SeqEntryPtr PNTR  sepp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_segset) return;
  sepp = (SeqEntryPtr PNTR) mydata;
  if (sepp == NULL) return;
  *sepp = sep;
}

//only for segsets
static void MarkPartsSets (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_parts) return;
  bssp->idx.deleteme = TRUE;
}

// only for segsets
static void ChangeSegSets (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  if (bssp->_class != BioseqseqSet_class_segset) return;
  bssp->_class = BioseqseqSet_class_nuc_prot;
}
//LCOV_EXCL_STOP

static SeqEntryPtr MakeVirtualProtein (void)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;

  bsp = BioseqNew ();
  if (bsp == NULL) return NULL;

  bsp->repr = Seq_repr_virtual;
  bsp->mol = Seq_mol_aa;
  bsp->length = 100;
  bsp->id = MakeSeqID ("dummy");

  sep = SeqEntryNew ();
  if (sep == NULL) return NULL;

  sep->choice = 1;
  sep->data.ptrvalue = (Pointer) bsp;

  SeqMgrAddToBioseqIndex (bsp);

  return sep;
}

static void LIBCALLBACK ErrRejCallback (
  ErrSev severity,
  int errcode,
  int subcode,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  CharPtr accession,
  CharPtr seqid,
  CharPtr featureID,
  CharPtr message,
  CharPtr objtype,
  CharPtr label,
  CharPtr context,
  CharPtr location,
  CharPtr product,
  Pointer userdata
)

{
  BoolPtr  hasErrOrRejP;

  if (userdata == NULL) return;
  hasErrOrRejP = (BoolPtr) userdata;

  if (severity == SEV_REJECT || severity == SEV_ERROR) {
    *hasErrOrRejP = TRUE;
  }
}

static Boolean HasErrorsOrReject (SeqEntryPtr sep)

{
  Boolean         hasErrOrRej = FALSE;
  ValidStructPtr  vsp;

  vsp = ValidStructNew ();
  if (vsp == NULL) return FALSE;

  vsp->useSeqMgrIndexes = TRUE;
  vsp->errfunc = ErrRejCallback;
  vsp->userdata = (Pointer) &hasErrOrRej;
  ValidateSeqEntry (sep, vsp);
  ValidStructFree (vsp);

  return hasErrOrRej;
}

static void MarkProteinBioseqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  if (bsp == NULL) return;
  if (ISA_aa (bsp->mol)) {
    bsp->idx.deleteme = TRUE;
  }
}

static void MarkCodingRegions (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL) return;
  if (sfp->data.choice == SEQFEAT_CDREGION) {
    sfp->idx.deleteme = TRUE;
  }
}

static void WholeCommentLoc (SeqFeatPtr sfp, Pointer userdata)

{
  BioseqPtr bsp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_COMMENT) return;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = SeqLocWholeNew (bsp);
}

static void RescueSegSetTitlesInt (SeqEntryPtr sep)

{
  CharPtr        first_title = NULL, newstr, str;
  BioseqPtr      nbsp;
  ObjValNodePtr  ovp;
  SeqEntryPtr    segsep;
  CharPtr        segTitlePrefix = "Original segmented set titles: ";
  SeqFeatPtr     sfp;
  ValNodePtr     titles = NULL, vnp1, vnp2;

  segsep = NULL;
  SeqEntryExplore (sep, (Pointer) &segsep, FindSegSet);
  if (segsep == NULL) {
    return;
  }
  VisitSetsInSep (sep, &first_title, RemoveSegSetTitle);
  VisitDescriptorsInSep (sep, (Pointer) &titles, SaveAllTitles);
  SeqEntryExplore (sep, NULL, RemoveAllTitles);

  if (titles != NULL) {
    nbsp = FindNucBioseq (sep);
    if (nbsp != NULL) {
      ValNodeUnique (&titles, SortVnpByString, ValNodeFreeData);
      if (titles != NULL && first_title == NULL) {
        first_title = StringSave ((CharPtr) titles->data.ptrvalue);
      }
      str = ValNodeMergeStrsEx (titles, "; ");
      if (str != NULL) {
        sfp = CreateNewFeatureOnBioseq (nbsp, SEQFEAT_COMMENT, NULL);
        if (sfp != NULL) {
          newstr = (CharPtr) MemNew (StringLen (segTitlePrefix) + StringLen (str) + 5);
          if (newstr != NULL) {
            StringCpy (newstr, segTitlePrefix);
            StringCat (newstr, str);
            sfp->comment = newstr;
          }
        }
        MemFree (str);
      }
    }
  }
  nbsp = FindNucBioseq (sep);
  if (nbsp != NULL) {
    for (vnp1 = nbsp->descr; vnp1 != NULL; vnp1 = vnp1->next) {
      if (vnp1->choice != Seq_descr_comment) continue;
      for (vnp2 = vnp1->next; vnp2 != NULL; vnp2 = vnp2->next) {
        if (vnp2->choice != Seq_descr_comment) continue;
        if (StringICmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0) continue;
        if (vnp2->extended != 0) {
          ovp = (ObjValNodePtr) vnp2;
          ovp->idx.deleteme = TRUE;
        }
      }
    }
    if (StringDoesHaveText (first_title)) {
      SeqDescrAddPointer (&(nbsp->descr), Seq_descr_title, first_title);
    }
  }
}


static void RescueSegSetTitles (SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   tmp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;

  if ((bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) ||
      bssp->_class == BioseqseqSet_class_wgs_set ||
      bssp->_class == BioseqseqSet_class_small_genome_set) {
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      RescueSegSetTitlesInt (tmp);
    }
  } else {
    RescueSegSetTitlesInt (sep);
  }
}

static void MoveTitlesToName (SeqDescrPtr sdp, Pointer userdata)

{
  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  sdp->choice = Seq_descr_name;
}

static void ChooseTitleOrName (BioseqPtr bsp, Pointer userdata)

{
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  ObjValNodePtr      ovp;
  SeqDescrPtr        sdp;

  if (bsp == NULL) return;
  if (SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_RNA, 0, &fcontext) == NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
    if (sdp != NULL && sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_name, &dcontext);
    if (sdp != NULL) {
      sdp->choice = Seq_descr_title;
    }
  } else {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_name, &dcontext);
    if (sdp != NULL && sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      ovp->idx.deleteme = TRUE;
    }
  }
}

static void PrefixTitleWithOrganism (BioseqPtr bsp, Pointer userdata)

{
  BioSourcePtr       biop;
  SeqMgrDescContext  dcontext;
  size_t             len;
  OrgModPtr          mod;
  CharPtr            old_name = NULL;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  CharPtr            str;
  CharPtr            title;
  CharPtr            taxname;

  if (bsp == NULL) return;
  biop = GetBiopForBsp (bsp);
  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  taxname = orp->taxname;
  if (StringHasNoText (taxname)) return;
  onp = orp->orgname;
  if (onp != NULL) {
    for (mod = onp->mod; mod != NULL; mod = mod->next) {
      if (mod->subtype != ORGMOD_old_name) continue;
      old_name = mod->subname;
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
  if (sdp == NULL) return;
  str = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (str)) return;

  if (StringNICmp (str, taxname, StringLen (taxname)) == 0) return;
  if (old_name != NULL && StringNICmp (str, old_name, StringLen (old_name)) == 0) return;

  len = StringLen (str) + StringLen (taxname) + 5;
  title = (CharPtr) MemNew (sizeof (Char) * len);
  if (title == NULL) return;
  StringCpy (title, taxname);
  StringCat (title, " ");
  StringCat (title, str);

  sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
  sdp->data.ptrvalue = title;
}

static time_t DoCleanup (
  SeqEntryPtr sep,
  Uint2 entityID,
  CleanFlagPtr cfp,
  AsnIoPtr aop,
  AsnTypePtr atp,
  SeqSubmitPtr ssp,
  BoolPtr failed
)

{
  Boolean              all_digits = TRUE, isDdbj = FALSE, isEmbl = FALSE,
                       isGenBank = FALSE, isNcbi = FALSE, isRefSeq = FALSE,
                       isNucProt = FALSE, isPPMEW = FALSE, isSegmented = FALSE,
                       stripSerial = TRUE, hasCDS = FALSE,
                       small_genome_set = FALSE;
  BioseqPtr            bsp, nbsp;
  BioseqSetPtr         bssp, nbssp;
  Char                 ch;
  Int4                 count, nucs, prts, pmid = 0, offset = 0;
  DatePtr              dp;
  SeqEntryPtr          fsep, nsep = NULL, tmp, segsep, ssep;
  ValNodePtr           head = NULL, tail = NULL;
  ObjValNodePtr        ovp;
  CharPtr              first_title = NULL, newstr, ptr, str;
  SumDataPtr           sdp;
  SeqFeatPtr           sfp;
  SeqIdPtr             sip, siphead;
  time_t               starttime, stoptime;
  ValNodePtr           titles = NULL;
  long int             val;
  SeqDescrPtr          vnp, vnp1, vnp2;
  ValidStructPtr       vsp;
  PopSetRetroStatData  stat;
  CharPtr              segTitlePrefix = "Original segmented set titles: ";

  if (sep == NULL || cfp == NULL) return 0;

  MemSet ((Pointer) &stat, 0, sizeof (PopSetRetroStatData));

  AssignIDsInEntityEx (entityID, 0, NULL, NULL);

  starttime = GetSecs ();

  StringCpy (cfp->buf, "");
  cfp->gi = 0;
  cfp->year = 0;
  cfp->isRefSeq = FALSE;
  cfp->isEmblDdbj = FALSE;
  cfp->isSmallGenome = FALSE;

  fsep = FindNthBioseq (sep, 1);
  if (fsep != NULL && fsep->choice == 1) {
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp != NULL) {
      siphead = SeqIdSetDup (bsp->id);
      for (sip = siphead; sip != NULL; sip = sip->next) {
        SeqIdStripLocus (sip);
        if (sip->choice == SEQID_GI) {
          cfp->gi = (Int4) sip->data.intvalue;
        } else if (sip->choice == SEQID_GENBANK || sip->choice == SEQID_TPG) {
          isGenBank = TRUE;
          isNcbi = TRUE;
        } else if (sip->choice == SEQID_EMBL || sip->choice == SEQID_TPE) {
          isEmbl = TRUE;
          cfp->isEmblDdbj = TRUE;
        } else if (sip->choice == SEQID_DDBJ || sip->choice == SEQID_TPD) {
          isDdbj = TRUE;
          cfp->isEmblDdbj = TRUE;
        } else if (sip->choice == SEQID_OTHER) {
          isRefSeq = TRUE;
          isNcbi = TRUE;
          cfp->isRefSeq = TRUE;
        }
      }
      SeqIdWrite (siphead, cfp->buf, PRINTID_FASTA_LONG, sizeof (cfp->buf));
      SeqIdSetFree (siphead);
    }
    vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_update_date, NULL);
    if (vnp == NULL) {
      vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_create_date, NULL);
    }
    if (vnp != NULL) {
      dp = (DatePtr) vnp->data.ptrvalue;
      if (dp != NULL && dp->data [0] == 1) {
        cfp->year = (Int2) dp->data [1] + 1900;
      }
    }
  }

  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtIDX);
  cfp->stripSerial = stripSerial;

  VisitSetsInSep (sep, (Pointer) &small_genome_set, IsSmallGenomeSet);
  cfp->isSmallGenome = small_genome_set;

  if (StringDoesHaveText (cfp->sourcedb)) {
    if (StringChr (cfp->sourcedb, 'g') != NULL) {
      if (! isGenBank) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->sourcedb, 'b') != NULL) {
      if ((! isEmbl) && (! isDdbj)) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->sourcedb, 'e') != NULL) {
      if (! isEmbl) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->sourcedb, 'd') != NULL) {
      if (! isDdbj) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->sourcedb, 'i') != NULL) {
      if ((! isGenBank) && (! isEmbl) && (! isDdbj)) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->sourcedb, 'r') != NULL) {
      if (! isRefSeq) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->sourcedb, 'n') != NULL) {
      if (! isNcbi) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->sourcedb, 'x') != NULL) {
      if (isEmbl || isDdbj) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
  }
  if (StringDoesHaveText (cfp->branch)) {
    if (StringChr (cfp->branch, 'c') != NULL) {
      VisitFeaturesInSep (sep, (Pointer) &hasCDS, CheckForCDSs);
      if (! hasCDS) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->branch, 'd') != NULL) {
      VisitFeaturesInSep (sep, (Pointer) &hasCDS, CheckForCDSs);
      if (hasCDS) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->branch, 'p') != NULL) {
      SeqMgrIndexFeatures (entityID, NULL);
      if (HasErrorsOrReject (sep)) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->branch, 'q') != NULL) {
      SeqMgrIndexFeatures (entityID, NULL);
      if (! HasErrorsOrReject (sep)) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
    }
    if (StringChr (cfp->branch, 'r') != NULL) {
//LCOV_EXCL_START
      VisitSetsInSep (sep, (Pointer) &isPPMEW, CheckForPPMEW);
      if (! isPPMEW) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
//LCOV_EXCL_STOP
    }
    if (StringChr (cfp->branch, 's') != NULL) {
//LCOV_EXCL_START
      VisitSetsInSep(sep, (Pointer)&isPPMEW, CheckForPPMEW);
      if (isPPMEW) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
//LCOV_EXCL_STOP
    }
    if (StringChr (cfp->branch, 't') != NULL) {
//LCOV_EXCL_START
      VisitSetsInSep (sep, (Pointer) &isNucProt, CheckForNucProt);
      if (! isNucProt) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
//LCOV_EXCL_STOP
    }
    if (StringChr (cfp->branch, 'u') != NULL) {
//LCOV_EXCL_START
      VisitSetsInSep(sep, (Pointer)&isNucProt, CheckForNucProt);
      if (isNucProt) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
//LCOV_EXCL_STOP
    }
    if (StringChr (cfp->branch, 'v') != NULL) {
//LCOV_EXCL_START
      VisitBioseqsInSep (sep, (Pointer) &isSegmented, CheckForSegSeq);
      if (! isSegmented) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
//LCOV_EXCL_STOP
    }
    if (StringChr (cfp->branch, 'w') != NULL) {
//LCOV_EXCL_START
      VisitBioseqsInSep(sep, (Pointer)&isSegmented, CheckForSegSeq);
      if (isSegmented) {
        if (failed != NULL) {
          *failed = TRUE;
        }
        return 0;
      }
//LCOV_EXCL_STOP
    }
  }

  nucs = VisitSequencesInSep (sep, NULL, VISIT_NUCS, NULL);
  prts = VisitSequencesInSep (sep, NULL, VISIT_PROTS, NULL);
  cfp->rawcounts.nucs += nucs;
  cfp->rawcounts.prts += prts;
  (cfp->rawcounts.recs)++;
  cfp->cumcounts.nucs += nucs;
  cfp->cumcounts.prts += prts;
  (cfp->cumcounts.recs)++;

  sdp = NULL;
  if (isGenBank) {
    sdp = &(cfp->dbsums.genbank);
  } else if (isEmbl) {
    sdp = &(cfp->dbsums.embl);
  } else if (isDdbj) {
    sdp = &(cfp->dbsums.ddbj);
  } else if (isRefSeq) {
    sdp = &(cfp->dbsums.refseq);
  } else {
    sdp = &(cfp->dbsums.other);
  }
  if (sdp != NULL) {
    sdp->nucs += nucs;
    sdp->prts += prts;
    (sdp->recs)++;
  }

  if (StringChr (cfp->report, 'c') != NULL) {
//LCOV_EXCL_START
    return 0;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'r') != NULL) {
//LCOV_EXCL_START
    DoASNReport(sep, cfp, FALSE, FALSE, FALSE);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 's') != NULL) {
//LCOV_EXCL_START
    DoASNReport(sep, cfp, TRUE, FALSE, FALSE);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'n') != NULL) {
//LCOV_EXCL_START
    DoASNReport (sep, cfp, TRUE, TRUE, FALSE);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'e') != NULL) {
//LCOV_EXCL_START
    DoASNReport (sep, cfp, FALSE, FALSE, TRUE);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'd') != NULL) {
//LCOV_EXCL_START
    DoAsnDiffReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'g') != NULL) {
//LCOV_EXCL_START
    DoGBFFReport (sep, entityID, cfp, 1);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'f') != NULL) {
//LCOV_EXCL_START
    DoGBFFReport (sep, entityID, cfp, 2);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'h') != NULL) {
//LCOV_EXCL_START
    DoGBFFReport (sep, entityID, cfp, 3);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'v') != NULL) {
//LCOV_EXCL_START
    DoValidatorReport (sep, cfp->logfp, cfp->buf, cfp->asnval);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'l') != NULL) {
//LCOV_EXCL_START
    DoLatLonReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'm') != NULL) {
//LCOV_EXCL_START
    DoModernizeReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'o') != NULL) {
//LCOV_EXCL_START
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, (Pointer) cfp, DoOverlapReport);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'u') != NULL) {
    DoUnpublishedList (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'p') != NULL) {
    DoPublishedReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'j') != NULL) {
    DoJournalReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'x') != NULL) {
    DoASNScan (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 't') != NULL) {
//LCOV_EXCL_START
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, (Pointer) cfp, DoAnticodonReport);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->report, 'w') != NULL) {
//LCOV_EXCL_START
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, (Pointer) cfp, DoContifOffsetReport);
    stoptime = GetSecs ();
    return stoptime - starttime;
//LCOV_EXCL_STOP
  }

  if (StringDoesHaveText (cfp->report)) return 0;

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s\n", cfp->buf);
    fflush (cfp->logfp);
  }

  if (cfp->taxon) {
    Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
  }

//LCOV_EXCL_START
  if (StringChr (cfp->pub, 'a') != NULL) {
    VisitFeaturesInSep (sep, (Pointer) cfp, MarkPubFeats);
    VisitDescriptorsInSep (sep, NULL, MarkPubDescs);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
  if (StringChr (cfp->pub, 's') != NULL) {
    ForceStripSerialNumber (sep);
  }
  if (StringChr (cfp->pub, 'f') != NULL) {
    VisitPubdescsInSep (sep, NULL, RemovePubFigNumName);
  }
  if (StringChr (cfp->pub, 'r') != NULL) {
    VisitPubdescsInSep (sep, NULL, RemovePubRemark);
  }
  if (StringChr (cfp->pub, 'm') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, NULL, MergePubFeats);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, NULL, PubDescTrumpsFeat);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }
//LCOV_EXCL_STOP
  if (StringChr (cfp->pub, 'u') != NULL) {
    VisitPubdescsInSep (sep, NULL, LookupPubdesc);
  }
  if (StringChr (cfp->pub, 'j') != NULL) {
    VisitPubdescsInSep (sep, NULL, LookupISOJTA);
  }
  if (cfp->pub != NULL) {
    ptr = cfp->pub;
    ch = *ptr;
    while (ch != '\0') {
      if (! IS_DIGIT (ch)) {
        all_digits = FALSE;
      }
      ptr++;
      ch = *ptr;
    }
    if (all_digits && sscanf (cfp->pub, "%ld", &val) == 1 && val != 0) {
      pmid = (Int4) val;
      VisitPubdescsInSep (sep, (Pointer) &pmid, ReplaceUnpub);
    }
  }

  if (StringChr (cfp->clean, 'b') != NULL) {
    BasicSeqEntryCleanup (sep);
  }
#ifdef OS_UNIX
  if (StringChr (cfp->clean, 'p') != NULL) {
//LCOV_EXCL_START
    nsep = CppBasicCleanup (sep, cfp);
//LCOV_EXCL_STOP
  }
#endif
  if (StringChr (cfp->clean, 'v') != NULL) {
    AdvancedSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->clean, 's') != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
  }
//LCOV_EXCL_START
  // option not used
  if (StringChr (cfp->clean, 'x') != NULL) {
    ExtendedSeqEntryCleanup (sep);
  }
  // GPIPE doesn't use cleanasn
  if (StringChr (cfp->clean, 'g') != NULL) {
    GpipeSeqEntryCleanup (sep);
  }
//LCOV_EXCL_STOP
  if (StringChr (cfp->clean, 'n') != NULL) {
    NormalizeDescriptorOrder (sep);
  }
  if (StringChr (cfp->clean, 'u') != NULL) {
    RemoveAllNcbiCleanupUserObjects (sep);
  }
  if (StringChr (cfp->clean, 'c') != NULL) {
//LCOV_EXCL_START
    CorrectGenCodes (sep, entityID);
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->clean, 'f') != NULL) {
//LCOV_EXCL_START
    CodingRegionPartialsFromTranslation (sep);
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->clean, 'e') != NULL) {
//LCOV_EXCL_START
    ImposeCodingRegionPartials (sep);
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->clean, 'd') != NULL) {
//LCOV_EXCL_START
    ResynchCodingRegionPartials (sep);
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->clean, 'm') != NULL) {
//LCOV_EXCL_START
    ResynchMessengerRNAPartials (sep);
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->clean, 't') != NULL) {
//LCOV_EXCL_START
    ResynchProteinPartials (sep);
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->clean, 'a') != NULL) {
//LCOV_EXCL_START
    AdjustSeqEntryForConsensusSplice (sep);
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->clean, 'i') != NULL) {
    PromoteToWorstIDProc (sep, entityID);
  }
  if (StringChr (cfp->clean, 'r') != NULL) {
    ReassignLocalIDProc (sep, entityID);
  }
  if (StringChr (cfp->clean, 'l') != NULL) {
    VisitBioseqsInSep (sep, NULL, StripLocusFromBsp);
  }

//LCOV_EXCL_START
  if (StringChr (cfp->modernize, 'g') != NULL) {
    VisitFeaturesInSep (sep, NULL, ModGenes);
  }
  if (StringChr (cfp->modernize, 'r') != NULL) {
    VisitFeaturesInSep (sep, NULL, ModRNAs);
  }
  if (StringChr (cfp->modernize, 'p') != NULL) {
    VisitBioSourcesInSep (sep, NULL, ModPCRs);
  }

  if (StringChr (cfp->feat, 'u') != NULL) {
    VisitFeaturesInSep (sep, NULL, RemoveFeatUser);
  }
  if (StringChr (cfp->feat, 'd') != NULL) {
    VisitFeaturesInSep (sep, NULL, RemoveFeatDbxref);
  }
  if (StringChr (cfp->feat, 'e') != NULL) {
    VisitFeaturesInSep (sep, NULL, RemoveFeatEvInf);
  }
  if (StringChr (cfp->feat, 'g') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    if (! small_genome_set) {
      VisitBioseqsInSep (sep, NULL, MultiIntervalGenesCallback);
    }
  }
  if (StringChr (cfp->feat, 'i') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    if (! small_genome_set) {
      VisitBioseqsInSep (sep, NULL, AdjacentIntervalImpFeatCallback);
    }
  }
  if (StringChr (cfp->feat, 'r') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, RemoveUnnecessaryGeneXrefs);
  }
  if (StringChr (cfp->feat, 'f') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, NULL, FuseDuplicateFeats);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
//LCOV_EXCL_STOP

  if (StringChr (cfp->feat, 's') != NULL) {
    SeqEntryExplore (sep, NULL, PackageFeatureCallback);
  }
  if (StringChr (cfp->feat, 'k') != NULL) {
    MoveFeatsFromPartsSet (sep);
    move_cds_ex (sep, TRUE);
  }
  if (StringChr (cfp->feat, 'z') != NULL) {
    BasicSeqEntryCleanup (sep);
    UpdateReplacedECNumbersEx (sep, &head, &tail, FALSE, FALSE);
    DeleteBadECNumbersEx (sep, &head, &tail, FALSE);
    BasicSeqEntryCleanup (sep);
    if (head != NULL) {
      if (cfp->logfp != NULL) {
        for (vnp = head; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          fprintf (cfp->logfp, "%s\n", str);
        }
      }
      ValNodeFreeData (head);
    }
  }
  if (StringChr (cfp->feat, 'b') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, SetBestFrameByLoc);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }
  if (StringChr (cfp->feat, 'x') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, RetranslateCDS);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }
  if (StringChr (cfp->feat, 'a') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, AdjustForStopCodon);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }

//LCOV_EXCL_START
  // unused options
  if (StringChr (cfp->conv, 'c') != NULL) {
    VisitBioseqsInSep (sep, NULL, CompressSeqs);
  }
  if (StringChr (cfp->conv, 'd') != NULL) {
    VisitBioseqsInSep (sep, NULL, ExpandSeqs);
  }
//LCOV_EXCL_STOP
  if (StringChr (cfp->conv, 'v') != NULL) {
//LCOV_EXCL_START
    // Only for SegSets
    SeqMgrIndexFeatures (entityID, 0);
    SegSeqNullToVirtual (sep);
//LCOV_EXCL_STOP
  }
  if (StringChr (cfp->conv, 's') != NULL) {
//LCOV_EXCL_START
    // Only for SegSets
    VisitBioseqsInSep (sep, NULL, RemoveSegSeqTitle);
    VisitSetsInSep (sep, NULL, RemoveSegSetTitle);
    /*
    VisitDescriptorsInSep (sep, (Pointer) &titles, SaveAllTitles);
    VisitBioseqsInSep (sep, NULL, RemoveSegSeqTitle);
    VisitSetsInSep (sep, &first_title, RemoveSegSetTitle);
    */
    VisitBioseqsInSep (sep, NULL, StripLocusFromBsp);
    /*
    VisitPubdescsInSep (sep, NULL, RemovePubSerial);
    */
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID, TRUE, TRUE);
    if (! ConvertSegSetToDeltaSeq (sep)) {
      /* do not output record if it was not converted */
      if (failed != NULL) {
        *failed = TRUE;
      }
    }

    /*
    if (titles != NULL) {
      nbsp = FindNucBioseq (sep);
      if (nbsp != NULL) {
        ValNodeUnique (&titles, SortVnpByString, ValNodeFreeData);
        if (titles != NULL && first_title == NULL) {
          first_title = StringSave ((CharPtr) titles->data.ptrvalue);
        }
        str = ValNodeMergeStrsEx (titles, "; ");
        if (str != NULL) {
          sfp = CreateNewFeatureOnBioseq (nbsp, SEQFEAT_COMMENT, NULL);
          if (sfp != NULL) {
            newstr = (CharPtr) MemNew (StringLen (segTitlePrefix) + StringLen (str) + 5);
            if (newstr != NULL) {
              StringCpy (newstr, segTitlePrefix);
              StringCat (newstr, str);
              sfp->comment = newstr;
            }
          }
        }
      }
    }
    nbsp = FindNucBioseq (sep);
    if (nbsp != NULL) {
      for (vnp1 = nbsp->descr; vnp1 != NULL; vnp1 = vnp1->next) {
        if (vnp1->choice != Seq_descr_comment) continue;
        for (vnp2 = vnp1->next; vnp2 != NULL; vnp2 = vnp2->next) {
          if (vnp2->choice != Seq_descr_comment) continue;
          if (StringICmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0) continue;
          if (vnp2->extended != 0) {
            ovp = (ObjValNodePtr) vnp2;
            ovp->idx.deleteme = TRUE;
          }
        }
      }
      DeleteMarkedObjects (entityID, 0, NULL);
      if (StringDoesHaveText (first_title)) {
        SeqDescrAddPointer (&(nbsp->descr), Seq_descr_title, first_title);
      }
    }
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    */
//LCOV_EXCL_STOP
  }

  if (StringChr (cfp->conv, 't') != NULL) {
//LCOV_EXCL_START
    // Only for SegSets
    SeqMgrClearFeatureIndexes (entityID, NULL);
    VisitBioseqsInSep (sep, NULL, RemoveLocalIDs);
    segsep = NULL;
    SeqEntryExplore (sep, (Pointer) &segsep, FindSegSet);
    if (segsep == NULL) {
      return 0;
    }
    VisitSetsInSep (sep, &first_title, RemoveSegSetTitle);
    DeleteMarkedObjects (entityID, 0, NULL);
    VisitDescriptorsInSep (sep, (Pointer) &titles, SaveAllTitles);
    SeqEntryExplore (sep, NULL, RemoveAllTitles);
    VisitBioseqsInSep (sep, NULL, StripLocusFromBsp);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrIndexFeatures (entityID, 0);
    /*
    DoAutoDef (sep, entityID, TRUE, TRUE);
    */
    ConvertSegSetsToDeltaSequences (sep);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    SeqEntryExplore (sep, NULL, MarkPartsSets);
    SeqEntryExplore (sep, NULL, ChangeSegSets);
    DeleteMarkedObjects (entityID, 0, NULL);
    if (titles != NULL) {
      nbsp = FindNucBioseq (sep);
      if (nbsp != NULL) {
        ValNodeUnique (&titles, SortVnpByString, ValNodeFreeData);
        if (titles != NULL && first_title == NULL) {
          first_title = StringSave ((CharPtr) titles->data.ptrvalue);
        }
        str = ValNodeMergeStrsEx (titles, "; ");
        if (str != NULL) {
          sfp = CreateNewFeatureOnBioseq (nbsp, SEQFEAT_COMMENT, NULL);
          if (sfp != NULL) {
            newstr = (CharPtr) MemNew (StringLen (segTitlePrefix) + StringLen (str) + 5);
            if (newstr != NULL) {
              StringCpy (newstr, segTitlePrefix);
              StringCat (newstr, str);
              sfp->comment = newstr;
            }
          }
        }
      }
    }
    nbsp = FindNucBioseq (sep);
    if (nbsp != NULL) {
      for (vnp1 = nbsp->descr; vnp1 != NULL; vnp1 = vnp1->next) {
        if (vnp1->choice != Seq_descr_comment) continue;
        for (vnp2 = vnp1->next; vnp2 != NULL; vnp2 = vnp2->next) {
          if (vnp2->choice != Seq_descr_comment) continue;
          if (StringICmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0) continue;
          if (vnp2->extended != 0) {
            ovp = (ObjValNodePtr) vnp2;
            ovp->idx.deleteme = TRUE;
          }
        }
      }
      DeleteMarkedObjects (entityID, 0, NULL);
      if (StringDoesHaveText (first_title)) {
        SeqDescrAddPointer (&(nbsp->descr), Seq_descr_title, first_title);
      }
    }
    SeriousSeqEntryCleanup (sep, NULL, NULL);

    /*
    after -C x do separate pass through -X s to regenerate definition lines,
    and -P m to merge publications, will have object manager problem if DoAutoDef
    or merge pubs attempted here
    */
//LCOV_EXCL_STOP
  }

  if (StringChr (cfp->conv, 'u') != NULL) {
//LCOV_EXCL_START
    // Only for SegSets
    SeqMgrClearFeatureIndexes (entityID, NULL);
    VisitBioseqsInSep (sep, NULL, RemoveLocalIDs);
    VisitDescriptorsInSep (sep, NULL, ZapBadTitles);
    DeleteMarkedObjects (entityID, 0, NULL);
    RescueSegSetTitles(sep);
    VisitBioseqsInSep (sep, NULL, StripLocusFromBsp);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrIndexFeatures (entityID, 0);
    if (! ConvertSegSetToDeltaSeqEx (sep, FALSE)) {
      /* do not output record if it was not converted */
      if (failed != NULL) {
        *failed = TRUE;
      }
    }
    SeqMgrClearFeatureIndexes (entityID, NULL);
    VisitBioseqsInSep (sep, NULL, MarkProteinBioseqs);
    VisitFeaturesInSep (sep, NULL, MarkCodingRegions);
    DeleteMarkedObjects (entityID, 0, NULL);
    VisitFeaturesInSep (sep, NULL, WholeCommentLoc);
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
//LCOV_EXCL_STOP
  }

  if (StringChr (cfp->conv, 'g') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, Cln_GlobalAlign2Seq, BioseqToDeltaByGapFeat);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
  if (StringChr (cfp->conv, 'm') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, NULL, BioseqToDeltaMergeGapFeat);
    DeleteMarkedObjects (entityID, 0, NULL);
  }

  if (StringChr (cfp->link, 'o') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    LinkCDSmRNAbyOverlap (sep);
  }
  if (StringChr (cfp->link, 'p') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    LinkCDSmRNAbyProduct (sep);
  }
  if (StringChr (cfp->link, 'l') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    LinkCDSmRNAbyLabelAndLocation (sep);
  }
  if (StringChr (cfp->link, 'r') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    ReassignFeatureIDs (sep);
  }
  if (StringChr (cfp->link, 'm') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp->_class == BioseqseqSet_class_genbank) {
        for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
          count = FindHighestFeatureID (tmp);
          OffsetFeatureIDs (tmp, offset);
          OffsetFeatureIDXrefs (tmp, offset);
          offset += count;
        }
      }
    }
  }
  if (StringChr (cfp->link, 'f') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, FixFeatureIDs);
  }
  if (StringChr (cfp->link, 'c') != NULL) {
    ClearFeatureIDs (sep);
  }

  if (StringChr (cfp->desc, 't') != NULL) {
    VisitDescriptorsInSep (sep, NULL, MarkTitles);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
//LCOV_EXCL_START
  if (StringChr (cfp->desc, 'c') != NULL) {
    VisitDescriptorsInSep (sep, NULL, MarkComments);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
  if (StringChr (cfp->desc, 'n') != NULL) {
    RemoveNucProtSetTitles (sep);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
//LCOV_EXCL_STOP
  if (StringChr (cfp->desc, 'e') != NULL) {
    RemovePopsetTitles (sep);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
//LCOV_EXCL_START
  if (StringChr (cfp->desc, 'm') != NULL) {
    RemoveMRnaTitles (sep);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
  if (StringChr (cfp->desc, 'p') != NULL) {
    RemoveProteinTitles (sep);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
  if (StringChr (cfp->desc, 'a') != NULL) {
    VisitDescriptorsInSep (sep, NULL, MoveTitlesToName);
    NormalizeDescriptorOrder (sep);
  }
  if (StringChr (cfp->desc, 'b') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, NULL, ChooseTitleOrName);
    DeleteMarkedObjects (entityID, 0, NULL);
    NormalizeDescriptorOrder (sep);
  }
  if (StringChr (cfp->desc, 'x') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitBioseqsInSep (sep, NULL, PrefixTitleWithOrganism);
     NormalizeDescriptorOrder (sep);
  }
//LCOV_EXCL_STOP

  if (StringDoesHaveText (cfp->zapuser)) {
    VisitDescriptorsInSep (sep, (Pointer) cfp->zapuser, MarkUserObjects);
    DeleteMarkedObjects (entityID, 0, NULL);
  }

//LCOV_EXCL_START
  // Use macro to add definition lines to a record instead
  if (StringChr (cfp->mods, 'd') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID, FALSE, FALSE);
    NormalizeDescriptorOrder (sep);
  }
  if (StringChr (cfp->mods, 's') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID, FALSE, TRUE);
    NormalizeDescriptorOrder (sep);
  }
//LCOV_EXCL_STOP
  // Regenerate definition lines
  if (StringChr(cfp->mods, 'r') != NULL) {
      SeqMgrIndexFeatures(entityID, 0);
      RegenerateDeflines(sep);
  }
//LCOV_EXCL_START
  if (StringChr (cfp->mods, 'e') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    PopSetAutoDefRetro (sep, &stat);
    LogPopSetTitleResults (cfp->logfp, &stat, sep);
  }
  if (StringChr (cfp->mods, 'n') != NULL) {
    SeqMgrIndexFeatures (entityID, NULL);
    SeqEntryExplore (sep, NULL, BadNCTitleProc);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateNCTitle (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->mods, 'm') != NULL) {
    SeqMgrIndexFeatures (entityID, NULL);
    SeqEntryExplore (sep, NULL, BadNMTitleProc);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateNMTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->mods, 'x') != NULL) {
    SeqMgrIndexFeatures (entityID, NULL);
    SeqEntryExplore (sep, NULL, BadNMTitleProc);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateSpecialXMTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->mods, 'p') != NULL) {
    SeqMgrIndexFeatures (entityID, NULL);
    SeqEntryExplore (sep, NULL, BadProtTitleProc);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateProteinTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->mods, 'g') != NULL) {
    SeqMgrIndexFeatures (entityID, NULL);
    VisitBioseqsInSep (sep, NULL, MarkNucTitles);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    VisitBioseqsInSep (sep, NULL, InstantiateGpipeTitles);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
//LCOV_EXCL_STOP
  if (StringChr (cfp->mods, 'c') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    CopyCDStomRNA (sep);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->mods, 'f') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    FixReciprocalProtTransIDs (sep);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
//LCOV_EXCL_START
  if (StringChr (cfp->mods, 'v') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    ReversePreOrNcRNATransIDs (sep);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->mods, 't') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, SeqToAnticodon);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->mods, 'z') != NULL) {
    // only for segsets
    SeqMgrClearFeatureIndexes (entityID, NULL);
    segsep = NULL;
    SeqEntryExplore (sep, (Pointer) &segsep, FindSegSet);
    if (segsep != NULL) {
      nbssp = BioseqSetNew();
      if (nbssp != NULL) {
        nbssp->_class = 1;
        ssep = SeqEntryNew ();
        if (ssep != NULL) {
          nbssp->seq_set = ssep;
          ssep->choice = 2;
          ssep->data.ptrvalue = segsep->data.ptrvalue;
          segsep->data.ptrvalue = nbssp;
        }
      }
    }
    BasicSeqEntryCleanup (sep);
  }
  //LCOV_EXCL_STOP
  if (StringChr (cfp->mods, 'b') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    CorrectGenCodes (sep, entityID);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, SetBestFrameByLoc);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, RetranslateCDS);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, AdjustForStopCodon);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeriousSeqEntryCleanup (sep, NULL, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, (Pointer) "unnamed protein product", AddMissingProtRef);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }
  if (StringChr (cfp->mods, 'w') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    CorrectGenCodes (sep, entityID);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, SetBestFrameByLoc);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    CodingRegionPartialsFromTranslation (sep);
    VisitFeaturesInSep (sep, NULL, CDSPartialIfFrame);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, RetranslateCDS);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    /*
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, AdjustForStopCodon);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    CodingRegionPartialsFromTranslation (sep);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    */

    SeriousSeqEntryCleanup (sep, NULL, NULL);

    RemoveProteinTitles (sep);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, (Pointer) "hypothetical protein", AddMissingProtRef);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, CopyCDSproductToMrna);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateProteinTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    ImposeCodingRegionPartials (sep);
    ResynchCodingRegionPartials (sep);
    ResynchMessengerRNAPartials (sep);
    ResynchProteinPartials (sep);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateProteinTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    NormalizeDescriptorOrder (sep);
    SeqMgrClearFeatureIndexes (entityID, NULL);
  }

  if (StringDoesHaveText (cfp->prot)) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, (Pointer) cfp->prot, AddMissingProtRef);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, 0);
    RemoveProteinTitles (sep);
    DeleteMarkedObjects (entityID, 0, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);

    SeqMgrIndexFeatures (entityID, NULL);
    SeqEntryExplore (sep, NULL, BadProtTitleProc);
    DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
    SeqMgrIndexFeatures (entityID, NULL);
    InstantiateProteinTitles (entityID, NULL);
    SeqMgrClearFeatureIndexes (entityID, NULL);
    BasicSeqEntryCleanup (sep);

    NormalizeDescriptorOrder (sep);
  }

  if (StringDoesHaveText (cfp->valid)) {
    SeqMgrIndexFeatures (entityID, NULL);
    vsp = ValidStructNew ();
    if (vsp != NULL) {
      vsp->useSeqMgrIndexes = TRUE;
      vsp->errfunc = ValidCallback;
      vsp->userdata = (Pointer) cfp->valid;
      ValidateSeqEntry (sep, vsp);
      ValidStructFree (vsp);
      DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
      SeqMgrIndexFeatures (entityID, NULL);
    }
  }

  if (cfp->action_list != NULL) {
    ApplyMacroToSeqEntry (sep, cfp->action_list);
  }

  /* normalize order again at end */
  if (StringChr (cfp->clean, 'n') != NULL) {
    NormalizeDescriptorOrder (sep);
  }

  stoptime = GetSecs ();

  if (aop != NULL) {
    if (failed != NULL && *failed) {
      /* skip if ConvertSegSetToDeltaSeq failed */
    } else if (ssp != NULL) {
      SeqSubmitAsnWrite (ssp, aop, atp);
    } else if (nsep != NULL) {
      SeqEntryAsnWrite (nsep, aop, atp);
      SeqEntryFree (nsep);
    } else {
      SeqEntryAsnWrite (sep, aop, atp);
    }
  }

  return stoptime - starttime;
}

static void CleanupSingleRecord (
  CharPtr filename,
  CleanFlagPtr cfp
)

{
  AsnIoPtr      aip, aop = NULL;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype, entityID = 0;
  Boolean       failed;
  FILE          *fp;
  Char          path [PATH_MAX];
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqSubmitPtr  ssp = NULL;

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

  cfp->filename = filename;

  if (cfp->type == 1) {
    fp = FileOpen (filename, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", filename);
      return;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (cfp->type >= 2 && cfp->type <= 5) {
    aip = AsnIoOpen (filename, cfp->binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", filename);
      return;
    }

    SeqMgrHoldIndexing (TRUE);
    switch (cfp->type) {
      case 2 :
        dataptr = (Pointer) SeqEntryAsnRead (aip, NULL);
        datatype = OBJ_SEQENTRY;
        break;
      case 3 :
        dataptr = (Pointer) BioseqAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQ;
        break;
      case 4 :
        dataptr = (Pointer) BioseqSetAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQSET;
        break;
      case 5 :
        dataptr = (Pointer) SeqSubmitAsnRead (aip, NULL);
        ssp = (SeqSubmitPtr) dataptr;
        datatype = OBJ_SEQSUB;
        break;
      default :
        break;
    }
    SeqMgrHoldIndexing (FALSE);

    AsnIoClose (aip);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else {
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) cfp->type);
    return;
  }

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", filename);
    return;
  }

  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

    sep = GetTopSeqEntryForEntityID (entityID);

    if (sep == NULL) {
      sep = SeqEntryNew ();
      if (sep != NULL) {
        if (datatype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) dataptr;
          sep->choice = 1;
          sep->data.ptrvalue = bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
        } else if (datatype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) dataptr;
          sep->choice = 2;
          sep->data.ptrvalue = bssp;
          SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
        } else {
          sep = SeqEntryFree (sep);
        }
      }
      sep = GetTopSeqEntryForEntityID (entityID);
    }

    if (sep != NULL) {

      path [0] = '\0';
      if (StringDoesHaveText (cfp->outfile)) {

        StringNCpy_0 (path, cfp->outfile, sizeof (path));

      } else if (StringDoesHaveText (cfp->results)) {

        ptr = StringRChr (filename, DIRDELIMCHR);
        if (ptr != NULL) {
          StringNCpy_0 (path, cfp->results, sizeof (path));
          ptr++;
          FileBuildPath (path, NULL, ptr);
        }
      }

      sep = GetTopSeqEntryForEntityID (entityID);
      if (sep != NULL) {

        if (StringHasNoText (cfp->report) && StringDoesHaveText (path)) {
          aop = AsnIoOpen (path, "w");
        }

        failed = FALSE;
        DoCleanup (sep, entityID, cfp, aop, NULL, ssp, &failed);

        if (aop != NULL) {
          AsnIoFlush (aop);
          AsnIoClose (aop);
          if (failed) {
            FileRemove (path);
          }
        }
      }

      ObjMgrFreeByEntityID (entityID);
    }

  } else {

    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }
}

static void CleanupMultipleRecord (
  CharPtr filename,
  CleanFlagPtr cfp
)

{
  AsnIoPtr     aip, aop = NULL;
  AsnTypePtr   atp = NULL;
  DataVal      av;
  Char         ch;
  Uint2        entityID;
  Boolean      failed;
  FILE         *fp;
  size_t       len;
  Char         longest [64];
  Int4         numrecords;
  Char         path [PATH_MAX];
  CharPtr      ptr;
  SeqEntryPtr  sep;
  time_t       timediff, worsttime;
#ifdef OS_UNIX
  Char         cmmd [512];
  CharPtr      gzcatprog;
  int          ret;
  Boolean      usedPopen = FALSE;
#endif

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

  cfp->filename = filename;

  path [0] = '\0';
  if (StringDoesHaveText (cfp->outfile)) {

    StringNCpy_0 (path, cfp->outfile, sizeof (path));

  } else if (StringDoesHaveText (cfp->results)) {

    ptr = StringRChr (filename, DIRDELIMCHR);
    if (ptr != NULL) {
      StringNCpy_0 (path, cfp->results, sizeof (path));
      ptr++;
      if (cfp->compressed) {
        len = StringLen (ptr);
        if (len > 4 && StringCmp (ptr + len - 3, ".gz") == 0) {
          ptr [len - 3] = '\0';
        }
      }
      FileBuildPath (path, NULL, ptr);
    }
  }
  if (StringHasNoText (cfp->report) && StringHasNoText (path)) return;

#ifndef OS_UNIX
  if (cfp->compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return;
  }
#endif

#ifdef OS_UNIX
  if (cfp->compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, filename);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", filename);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", filename);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return;
        }
      }
    }
    fp = popen (cmmd, /* cfp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (filename, cfp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (filename, cfp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", filename);
    return;
  }

  aip = AsnIoNew (cfp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", filename);
    return;
  }

  if (cfp->logfp != NULL) {
    if (StringChr (cfp->report, 'c') == NULL) {
      fprintf (cfp->logfp, "%s\n\n", filename);
      fflush (cfp->logfp);
    }
  }

  longest [0] = '\0';
  worsttime = 0;
  numrecords = 0;

  if (StringHasNoText (cfp->report)) {
    aop = AsnIoOpen (path, cfp->binary? "wb" : "w");
    if (aop != NULL) {
      if (cfp->type == 4) {
        AsnOpenStruct (aop, cfp->bssp_atp, (Pointer) &(cfp->bss));
        av.intvalue = 7;
        AsnWrite (aop, cfp->atp_bsc, &av);
        AsnOpenStruct (aop, cfp->atp_bsss, (Pointer) &(cfp->bss.seq_set));
      }
    }
  }

  if (cfp->type == 4) {
    atp = cfp->atp_bss;
  } else if (cfp->type == 5) {
    atp = cfp->atp_ssp;
  }

  if (atp == NULL) {
    Message (MSG_POSTERR, "Batch processing type not set properly");
    return;
  }

  while ((atp = AsnReadId (aip, cfp->amp, atp)) != NULL) {
    if (atp == cfp->atp_se) {

      SeqMgrHoldIndexing (TRUE);
      sep = SeqEntryAsnRead (aip, atp);
      SeqMgrHoldIndexing (FALSE);

      if (sep != NULL) {

        entityID = ObjMgrGetEntityIDForChoice (sep);

        failed = FALSE;
        timediff = DoCleanup (sep, entityID, cfp, aop, cfp->atp_se, NULL, &failed);

        if (timediff > worsttime) {
          worsttime = timediff;
          StringCpy (longest, cfp->buf);
          ptr = longest;
          ch = *ptr;
          while (ch != '\0') {
            if (ch == '|') {
              *ptr = ' ';
            }
            ptr++;
            ch = *ptr;
          }
        }
        numrecords++;

        ObjMgrFreeByEntityID (entityID);
      }

    } else if (cfp->type == 4) {

      AsnReadVal (aip, atp, NULL);

    } else {

      AsnReadVal (aip, atp, &av);
      AsnWrite (aop, atp, &av);
      AsnKillValue (atp, &av);
    }
  }

  if (aop != NULL) {
    if (cfp->type == 4) {
      AsnCloseStruct (aop, cfp->atp_bsss, (Pointer) &(cfp->bss.seq_set));
      AsnCloseStruct (aop, cfp->bssp_atp, (Pointer) &(cfp->bss));
    }
    AsnIoClose (aop);
  }

  AsnIoFree (aip, FALSE);

#ifdef OS_UNIX
  if (usedPopen) {
    pclose (fp);
  } else {
    FileClose (fp);
  }
#else
  FileClose (fp);
#endif
  if (cfp->logfp != NULL) {
    if (StringChr (cfp->report, 'c') == NULL) {
      fprintf (cfp->logfp, "\nTotal number of records %ld\n", (long) numrecords);
      if (StringDoesHaveText (longest)) {
        fprintf (cfp->logfp, "Longest processing time %ld seconds on %s\n",
                 (long) worsttime, longest);
      }
      fprintf (cfp->logfp, "Counts ");
      fprintf (cfp->logfp, "- %9ld RECS", (long) cfp->rawcounts.recs);
      fprintf (cfp->logfp, ", %9ld NUCS", (long) cfp->rawcounts.nucs);
      fprintf (cfp->logfp, ", %9ld PRTS", (long) cfp->rawcounts.prts);
      fprintf (cfp->logfp, ", %9ld OKAY", (long) cfp->rawcounts.okay);
      fprintf (cfp->logfp, ", %9ld NORM", (long) cfp->rawcounts.norm);
      fprintf (cfp->logfp, ", %9ld CLNR", (long) cfp->rawcounts.clnr);
      fprintf (cfp->logfp, ", %9ld OTHR", (long) cfp->rawcounts.othr);
      fprintf (cfp->logfp, ", %9ld MODR", (long) cfp->rawcounts.modr);
      fprintf (cfp->logfp, ", %9ld SLOC", (long) cfp->rawcounts.sloc);
      fprintf (cfp->logfp, ", %9ld PUBL", (long) cfp->rawcounts.publ);
      fprintf (cfp->logfp, ", %9ld AUTH", (long) cfp->rawcounts.auth);
      fprintf (cfp->logfp, ", %9ld SORT", (long) cfp->rawcounts.sort);
      fprintf (cfp->logfp, ", %9ld BSEC", (long) cfp->rawcounts.bsec);
      fprintf (cfp->logfp, ", %9ld GBBK", (long) cfp->rawcounts.gbbk);
      fprintf (cfp->logfp, ", %9ld TITL", (long) cfp->rawcounts.titl);
      fprintf (cfp->logfp, ", %9ld PACK", (long) cfp->rawcounts.pack);
      fprintf (cfp->logfp, ", %9ld MOVE", (long) cfp->rawcounts.move);
      fprintf (cfp->logfp, ", %9ld SSEC", (long) cfp->rawcounts.ssec);
      fprintf (cfp->logfp, "\n");
      fflush (cfp->logfp);
    }
  }
}

static void CleanupOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       ptr;
  SumDataPtr    sdp;

  if (StringHasNoText (filename)) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  cfp->filename = NULL;

  MemSet ((Pointer) &(cfp->rawcounts), 0, sizeof (CountData));
  MemSet ((Pointer) &(cfp->dbsums), 0, sizeof (DbSumData));

  if (StringChr (cfp->sourcedb, 'y') != NULL) {
    ptr = StringRChr (filename, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      if (StringStr (ptr, "gbcon") != NULL ||
          StringStr (ptr, "gbest") != NULL ||
          StringStr (ptr, "gbgss") != NULL ||
          StringStr (ptr, "gbhtg") != NULL ||
          StringStr (ptr, "gbpat") != NULL ||
          StringStr (ptr, "gbsts") != NULL) return;
    }
  }

  if (cfp->batch) {
    ptr = StringRChr (filename, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      if (StringDoesHaveText (cfp->firstfile)) {
        if (StringICmp (cfp->firstfile, ptr) == 0) {
          cfp->foundfirst = TRUE;
        }
        if (! cfp->foundfirst) return;
      }

      if (StringDoesHaveText (cfp->lastfile)) {
        if (cfp->foundlast) return;
        if (StringICmp (cfp->lastfile, ptr) == 0) {
          cfp->foundlast = TRUE;
        }
      }
    }

    CleanupMultipleRecord (filename, cfp);
  } else {
    CleanupSingleRecord (filename, cfp);
  }

  if (cfp->logfp != NULL) {
    if (StringChr (cfp->report, 'c') != NULL) {
      ptr = StringRChr (filename, DIRDELIMCHR);
      if (ptr != NULL) {
        ptr++;
        fprintf (cfp->logfp, "%s", ptr);
      }
      sdp = &(cfp->dbsums.genbank);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      sdp = &(cfp->dbsums.embl);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      sdp = &(cfp->dbsums.ddbj);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      sdp = &(cfp->dbsums.refseq);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      sdp = &(cfp->dbsums.other);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) cfp->rawcounts.recs, (long) cfp->rawcounts.nucs, (long) cfp->rawcounts.prts);
      fprintf (cfp->logfp, "\n");
      fflush (cfp->logfp);
    }
  }
}

static Boolean PvtStrToLong (CharPtr str, Int4Ptr longval)

{
  Char     ch;
  Int2     i;
  Int2     len;
  Char     local [64];
  Boolean  nodigits;
  Boolean  rsult;
  long int     val;

  rsult = FALSE;
  if (longval != NULL) {
    *longval = (Int4) 0;
  }
  len = (Int2) StringLen (str);
  if (len != 0) {
    rsult = TRUE;
    nodigits = TRUE;
    for (i = 0; i < len; i++) {
      ch = str [i];
      if (ch == ' ' || ch == '+' || ch == '-') {
      } else if (ch < '0' || ch > '9') {
        rsult = FALSE;
      } else {
        nodigits = FALSE;
      }
    }
    if (nodigits) {
      rsult = FALSE;
    }
    if (rsult && longval != NULL) {
      StringNCpy_0 (local, str, sizeof (local));
      if (sscanf (local, "%ld", &val) == 1) {
        *longval = val;
      }
    }
  }
  return rsult;
}

static Int4 AccessionToGi (CharPtr string)

{
   Int4      gi;
   SeqIdPtr  sip;

   sip = SeqIdFromAccessionDotVersion (string);
   if (sip == NULL) return 0;
   gi = GetGIForSeqId (sip);
   SeqIdFree (sip);
   return gi;
}

static void UnlockSeq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  ObjMgrDataPtr  omdp;

  if (bsp == NULL) return;
  omdp = SeqMgrGetOmdpForBioseq (bsp);
  if (omdp == NULL) return;
  omdp->lockcnt = 0;
}

static void UnlockSet (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  ObjMgrDataPtr  omdp;
  ObjMgrPtr      omp;

  if (bssp == NULL) return;
  omp = ObjMgrWriteLock ();
  omdp = ObjMgrFindByData (omp, bssp); 
  if (omdp != NULL) {
    omdp->lockcnt = 0;
  }
  ObjMgrUnlock ();
}

static void UnlockAll (SeqEntryPtr sep)

{
  if (sep == NULL) return;

  SeqMgrWriteLock ();

  VisitBioseqsInSep (sep, NULL, UnlockSeq);
  VisitSetsInSep (sep, NULL, UnlockSet);

  SeqMgrUnlock ();
}

static void ProcessAccessionList (
  CharPtr accnfile,
  CleanFlagPtr cfp
)

{
  AsnIoPtr     aop = NULL;
  BioseqPtr    bsp;
  Uint2        entityID;
  Boolean      failed;
  FileCache    fc;
  FILE         *fp;
  Char         line [1023];
  ObjMgrPtr    omp;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  ErrSev       sev;
  SeqIdPtr     sip;
  CharPtr      str;
  Char         tmp [64];

  if (StringHasNoText (accnfile) || cfp == NULL) return;

  fp = FileOpen (accnfile, "r");
  if (fp == NULL) return;

  if (FileCacheSetup (&fc, fp)) {
    str = FileCacheGetString (&fc, line, sizeof (line));
    while (str != NULL) {
      TrimSpacesAroundString (str);

      bsp = NULL;
      entityID = 0;
      sep = NULL;
      sip = NULL;

      if (StringIsAllDigits (str)) {
        sprintf (tmp, "gi|%s", str);
        sip = SeqIdParse (tmp);
      } else {
        sip = SeqIdFromAccessionDotVersion (str);
      }

      if (sip != NULL) {
        sev = ErrSetMessageLevel (SEV_MAX);
        bsp = BioseqLockById (sip);
        ErrSetMessageLevel (sev);
        sip = SeqIdFree (sip);
      }

      if (bsp != NULL) {
        entityID = ObjMgrGetEntityIDForPointer (bsp);
        sep = GetTopSeqEntryForEntityID (entityID);

        if (sep != NULL) {

          path [0] = '\0';
          if (StringDoesHaveText (cfp->outfile)) {

            StringNCpy_0 (path, cfp->outfile, sizeof (path));

          } else if (StringDoesHaveText (cfp->results)) {

            StringNCpy_0 (path, cfp->results, sizeof (path));
            FileBuildPath (path, NULL, str);
            StringCat (path, ".asn");
          }

          entityID = ObjMgrGetEntityIDForChoice (sep);
          if (entityID > 0) {

            if (StringHasNoText (cfp->report) && StringDoesHaveText (path)) {
              aop = AsnIoOpen (path, "w");
            }

            failed = FALSE;
            DoCleanup (sep, entityID, cfp, aop, NULL, NULL, &failed);

            if (aop != NULL) {
              AsnIoFlush (aop);
              AsnIoClose (aop);
              if (failed) {
                FileRemove (path);
              }
            }
          }
        }

        /*
        BioseqUnlock (bsp);
        SeqEntryFree (sep);
        */

        UnlockAll (sep);

        ObjMgrFree (OBJ_SEQENTRY, (Pointer) sep);
  
        omp = ObjMgrGet ();
        ObjMgrReapOne (omp);
        SeqMgrClearBioseqIndex ();
        ObjMgrFreeCache (0);
        FreeSeqIdGiCache ();
  
        SeqEntrySetScope (NULL);
      }

      str = FileCacheGetString (&fc, line, sizeof (line));
    }
  }

  FileClose (fp);
}

/* Args structure contains command-line arguments */

typedef enum {
  p_argInputPath = 0,
  r_argOutputPath,
  i_argInputFile,
  o_argOutputFile,
  f_argFilter,
  x_argSuffix,
  j_argFirstFile,
  k_argLastFile,
  d_argSourceDb,
  a_argType,
  b_argBinary,
  c_argCompressed,
  L_argLogFile,
  R_argRemote,
  B_argBranch,
  Q_argReport,
  S_argSelective,
  m_argFfMode,
  q_argFfDiff,
  n_argAsn2Flat,
  v_argAsnVal,
  K_argClean,
  U_argModernize,
  N_argLink,
  F_argFeat,
  D_argDesc,
  C_argConvert,
  P_argPub,
  Z_argUser,
  X_argMods,
  O_argProtName,
  V_argValid,
  T_argTaxon,
  M_argMacro,
  A_argAccnFile,
#ifdef INTERNAL_NCBI_CLEANASN
  H_argHup,
#endif
} Arguments;

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", "stdout", NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".ent", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"First File Name", NULL, NULL, NULL,
    TRUE, 'j', ARG_STRING, 0.0, 0, NULL},
  {"Last File Name", NULL, NULL, NULL,
    TRUE, 'k', ARG_STRING, 0.0, 0, NULL},
  {"Source Database\n"
   "      a Any\n"
   "      g GenBank\n"
   "      e EMBL\n"
   "      d DDBJ\n"
   "      b EMBL or DDBJ\n"
   "      i INSD\n"
   "      r RefSeq\n"
   "      n NCBI\n"
   "      x Exclude EMBL/DDBJ\n"
   "      y Exclude gbcon, gbest, gbgss, gbhtg, gbpat, gbsts\n", "a", NULL, NULL,
    TRUE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"ASN.1 Type\n"
   "      a Any\n"
   "      e Seq-entry\n"
   "      b Bioseq\n"
   "      s Bioseq-set\n"
   "      m Seq-submit\n"
   "      t Batch Bioseq-set\n"
   "      u Batch Seq-submit\n", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log File", NULL, NULL, NULL,
    TRUE, 'L', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Branch\n"
   "      c Has Coding Regions\n"
   "      d No Coding Regions\n"
   "      p Passes Validation\n"
   "      q Validator Errors or Rejects\n"
   "      r Only Pop/Phy/Mut/Eco/WGS Sets\n"
   "      s Exclude Pop/Phy/Mut/Eco/WGS Sets\n"
   "      t Only Nuc-Prot Sets\n"
   "      u Exclude Nuc-Prot Sets\n"
   "      v Only Segmented Sequences\n"
   "      w Exclude Segmented Sequences\n", NULL, NULL, NULL,
    TRUE, 'B', ARG_STRING, 0.0, 0, NULL},
  {"Report\n"
   "      c Record Count\n"
   "      r ASN.1 BSEC Report\n"
   "      s ASN.1 SSEC Report\n"
   "      n NORM vs. SSEC Report\n"
   "      e PopPhyMutEco AutoDef Report\n"
   "      o Overlap Report\n"
   "      l Latitude-Longitude Country Diff\n"
   "      d Log SSEC Differences\n"
   "      g GenBank SSEC Diff\n"
   "      f asn2gb/asn2flat Diff\n"
   "      h Seg-to-Delta GenBank Diff\n"
   "      v Validator SSEC Diff\n"
   "      m Modernize Gene/RNA/PCR\n"
   "      u Unpublished Pub Lookup\n"
   "      p Published Pub Lookup\n"
   "      j Unindexed Journal Report\n"
   "      t tRNA Anticodon Report\n"
   "      w Component Offset Report\n"
   "      x Custom Scan\n", NULL, NULL, NULL,
    TRUE, 'Q', ARG_STRING, 0.0, 0, NULL},
  {"Selective Difference Filter\n"
   "      s SSEC\n"
   "      b BSEC\n"
   "      a Author\n"
   "      p Publication\n"
   "      l Location\n"
   "      r RNA\n"
   "      q Qualifier Sort Order\n"
   "      g Genbank Block\n"
   "      k Package CdRegion or Parts Features\n"
   "      m Move Publication\n"
   "      o Leave Duplicate Bioseq Publication\n"
   "      d Automatic Definition Line\n"
   "      e Pop/Phy/Mut/Eco Set Definition Line\n"
   "      (Capital Letters Skip)\n", NULL, NULL, NULL,
    TRUE, 'S', ARG_STRING, 0.0, 0, NULL},
  {"Flatfile Mode\n"
   "      r Release\n"
   "      e Entrez\n"
   "      s Sequin\n"
   "      d Dump\n", NULL, NULL, NULL,
    TRUE, 'm', ARG_STRING, 0.0, 0, NULL},
  {"ffdiff Executable", "/netopt/genbank/subtool/bin/ffdiff", NULL, NULL,
    TRUE, 'q', ARG_FILE_IN, 0.0, 0, NULL},
  {"asn2flat Executable", "/netopt/ncbi_tools/bin/asn2flat", NULL, NULL,
    TRUE, 'n', ARG_FILE_IN, 0.0, 0, NULL},
  {"asnval Executable", "/netopt/ncbi_tools/bin/asnval", NULL, NULL,
    TRUE, 'v', ARG_FILE_IN, 0.0, 0, NULL},
  {"Cleanup\n"
   "      b BasicSeqEntryCleanup\n"
   "      p C++ BasicCleanup\n"
   "      v AdvancedSeqEntryCleanup\n"
   "      s SeriousSeqEntryCleanup\n"
   "      x ExtendedSeqEntryCleanup\n"
   "      g GpipeSeqEntryCleanup\n"
   "      n Normalize Descriptor Order\n"
   "      u Remove NcbiCleanup User Objects\n"
   "      c Synchronize Genetic Codes\n"
   "      f CDS Partials from Translation\n"
   "      e Impose CDS Partials\n"
   "      d Resynchronize CDS Partials\n"
   "      m Resynchronize mRNA Partials\n"
   "      t Resynchronize Peptide Partials\n"
   "      a Adjust Consensus Splice\n"
   "      i Promote to Worst Seq-ID\n"
   "      r Reassign Local IDs\n"
   "      l Remove Locus\n", NULL, NULL, NULL,
    TRUE, 'K', ARG_STRING, 0.0, 0, NULL},
  {"Modernize\n"
   "      g Gene\n"
   "      r RNA\n"
   "      p PCR Primers\n", NULL, NULL, NULL,
    TRUE, 'U', ARG_STRING, 0.0, 0, NULL},
  {"Link\n"
   "      o Link CDS mRNA by Overlap\n"
   "      p Link CDS mRNA by Product\n"
   "      l Link CDS mRNA by Label and Location\n"
   "      r Reassign Feature IDs\n"
   "      m Merge Colliding Feature ID\n"
   "      f Fix Missing Reciprocal Feature IDs\n"
   "      c Clear Feature IDs\n", NULL, NULL, NULL,
    TRUE, 'N', ARG_STRING, 0.0, 0, NULL},
  {"Feature\n"
   "      u Remove User Object\n"
   "      d Remove db_xref\n"
   "      e Remove /evidence and /inference\n"
   "      g Fuse Multi-Interval Genes\n"
   "      i Fuse Adjacent-Interval Imp Feats\n"
   "      r Remove Redundant Gene xref\n"
   "      f Fuse Duplicate Features\n"
   "      s Package Features on Referenced Bioseq\n"
   "      k Package CdRegion or Parts Features\n"
   "      z Delete or Update EC Numbers\n"
   "      b Set Best CdRegion Reading Frame\n"
   "      x Retranslate CdRegions\n"
   "      a Adjust for Missing Stop Codon\n", NULL, NULL, NULL,
    TRUE, 'F', ARG_STRING, 0.0, 0, NULL},
  {"Descriptor\n"
   "      t Remove Title\n"
   "      c Remove Comment\n"
   "      n Remove Nuc-Prot Set Title\n"
   "      e Remove Pop/Phy/Mut/Eco Set Title\n"
   "      m Remove mRNA Title\n"
   "      p Remove Protein Title\n"
   "      a Title to Name\n"
   "      b AutoDef Title or Name\n"
   "      x Prefixf Title with Organism Name\n", NULL, NULL, NULL,
    TRUE, 'D', ARG_STRING, 0.0, 0, NULL},
  {"Sequence\n"
   "      c Compress\n"
   "      d Decompress\n"
   "      v Virtual Gaps inside Seg Seq\n"
   "      s Convert Seg Set to Delta Seq\n"
   "      t Non-NucProt Seg Set to Delta Seq\n"
   "      u Improved Non-NucProt Seg Set to Delta Seq\n"
   "      g Raw to Delta by Assembly Gap\n"
   "      m Merge Assembly Gap Features\n", NULL, NULL, NULL,
    TRUE, 'C', ARG_STRING, 0.0, 0, NULL},
  {"Publication\n"
   "      a Remove All Publications\n"
   "      s Remove Serial Number\n"
   "      f Remove Figure, Numbering, and Name\n"
   "      r Remove Remark\n"
   "      u Update PMID-only Publication\n"
   "      j Lookup ISO Journal Title Abbreviation\n"
   "      m Merge Identical Publication Features\n"
   "      # Replace Unpublished With PMID\n", NULL, NULL, NULL,
    TRUE, 'P', ARG_STRING, 0.0, 0, NULL},
  {"Remove Indicated User Object", NULL, NULL, NULL,
    TRUE, 'Z', ARG_STRING, 0.0, 0, NULL},
  {"Miscellaneous\n"
   "      d Automatic Definition Line\n"
   "      s Automatic Definition Line with Source Qualifiers\n"
   "      e Pop/Phy/Mut/Eco Set Definition Line\n"
   "      n Instantiate NC Title\n"
   "      m Instantiate NM Titles\n"
   "      x Special XM Titles\n"
   "      p Instantiate Protein Titles\n"
   "      g Gpipe Instantiate Titles\n"
   "      c Create mRNAs for CDSs\n"
   "      f Fix Reciprocal protein_id/transcript_id\n"
   "      v Revert preRNA or ncRNA transcript_id\n"
   "      t Parse Anticodon from Sequence\n"
   "      b Batch Cleanup of Multireader Output\n"
   "      z Wrap SegSet with NucProt Set\n"
   "      w GFF/WGS Genome Cleanup\n", NULL, NULL, NULL,
    TRUE, 'X', ARG_STRING, 0.0, 0, NULL},
  {"Missing Prot-Ref Name", NULL, NULL, NULL,
    TRUE, 'O', ARG_STRING, 0.0, 0, NULL},
  {"Validation\n"
   "      Remove Features by Validator Severity - r Reject, e Error, w Warning, i Info\n", NULL, NULL, NULL,
    TRUE, 'V', ARG_STRING, 0.0, 0, NULL},
  {"Taxonomy Lookup", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Macro File", NULL, NULL, NULL,
    TRUE, 'M', ARG_FILE_IN, 0.0, 0, NULL},
  {"Accession List File", NULL, NULL, NULL,
    TRUE, 'A', ARG_FILE_IN, 0.0, 0, NULL},
#ifdef INTERNAL_NCBI_CLEANASN
  {"Internal Access to HUP", "F", NULL, NULL,
    TRUE, 'H', ARG_BOOLEAN, 0.0, 0, NULL},
#endif
};

Int2 Main (void)

{
  ValNodePtr     action_list, parflat_list, vnp;
  AsnIoPtr       aip;
  Char           app [64], mode, type;
  CleanFlagData  cfd;
  CharPtr        directory, filter, accnfile, infile, logfile, outfile,
                 macro_file, results, str, suffix;
  Boolean        remote;
  time_t         runtime, starttime, stoptime;
#ifdef INTERNAL_NCBI_CLEANASN
  Boolean        hup;
#endif

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  SOCK_SetupSSL(NcbiSetupGnuTls);

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

  parflat_list = Validate_ParFlat_GBFeat ();
  if (parflat_list != NULL) {
    Message (MSG_POSTERR, "Validate_ParFlat_GBFeat warnings");
    for (vnp = parflat_list; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      Message (MSG_POSTERR, "%s", str);
    }
    ValNodeFreeData (parflat_list);
  }

  /* process command line arguments */

  sprintf (app, "cleanasn %s", CLEANASN_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &cfd, 0, sizeof (CleanFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = directory;
  }
  accnfile = (CharPtr) myargs [A_argAccnFile].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  cfd.batch = FALSE;
  cfd.binary = (Boolean) myargs [b_argBinary].intvalue;
  cfd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  cfd.type = 1;

  cfd.foundfirst = FALSE;
  cfd.foundlast = FALSE;
  cfd.sourcedb = myargs [d_argSourceDb].strvalue;
  cfd.branch = myargs [B_argBranch].strvalue;

  str = myargs [a_argType].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    type = str [0];
  } else {
    type = 'a';
  }

  type = TO_LOWER (type);
  switch (type) {
    case 'a' :
      cfd.type = 1;
      break;
    case 'e' :
      cfd.type = 2;
      break;
    case 'b' :
      cfd.type = 3;
      break;
    case 's' :
      cfd.type = 4;
      break;
    case 'm' :
      cfd.type = 5;
      break;
    case 't' :
      cfd.type = 4;
      cfd.batch = TRUE;
      break;
    case 'u' :
      cfd.type = 5;
      cfd.batch = TRUE;
      break;
    default :
      cfd.type = 1;
      break;
  }

  remote = (Boolean) myargs [R_argRemote].intvalue;
  if (StringDoesHaveText (accnfile)) {
    remote = TRUE;
  }
#ifdef INTERNAL_NCBI_CLEANASN
  hup = (Boolean) myargs [H_argHup].intvalue;
#endif

  cfd.report = myargs [Q_argReport].strvalue;
  cfd.selective = myargs [S_argSelective].strvalue;
  cfd.ffdiff = myargs [q_argFfDiff].strvalue;
  cfd.asn2flat = myargs [n_argAsn2Flat].strvalue;
  cfd.asnval = myargs [v_argAsnVal].strvalue;

  str = myargs [m_argFfMode].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    mode = str [0];
  } else {
    mode = 'e';
  }

  mode = TO_LOWER (mode);
  switch (mode) {
    case 'r' :
      cfd.ffmode = RELEASE_MODE;
      break;
    case 'e' :
      cfd.ffmode = ENTREZ_MODE;
      break;
    case 's' :
      cfd.ffmode = SEQUIN_MODE;
      break;
    case 'd' :
      cfd.ffmode = DUMP_MODE;
      break;
    default :
      cfd.ffmode = ENTREZ_MODE;
      break;
  }

  cfd.clean = myargs [K_argClean].strvalue;
  cfd.modernize = myargs [U_argModernize].strvalue;
  cfd.link = myargs [N_argLink].strvalue;
  cfd.feat = myargs [F_argFeat].strvalue;
  cfd.desc = myargs [D_argDesc].strvalue;
  cfd.conv = myargs [C_argConvert].strvalue;
  cfd.zapuser = myargs [Z_argUser].strvalue;
  cfd.mods = myargs [X_argMods].strvalue;
  cfd.prot = myargs [O_argProtName].strvalue;
  cfd.valid = myargs [V_argValid].strvalue;
  cfd.taxon = (Boolean) myargs [T_argTaxon].intvalue;
  cfd.pub = myargs [P_argPub].strvalue;
  cfd.unpubcount = 0;
  cfd.unpublist = NULL;
  cfd.lastunpub = NULL;

  macro_file = myargs [M_argMacro].strvalue;
  if (StringDoesHaveText (macro_file)) {
    aip = AsnIoOpen (macro_file, "r");
    if (aip == NULL) {
      Message (MSG_FATAL, "Unable to open macro file '%s'", macro_file);
      return 1;
    }
    action_list = MacroActionListAsnRead (aip, NULL);
    AsnIoClose (aip);
    if (action_list == NULL) {
      Message (MSG_FATAL, "Unable to read macro file '%s'", macro_file);
    }
    cfd.action_list = action_list;
  }

  cfd.amp = AsnAllModPtr ();
  cfd.atp_ssp = AsnFind ("Seq-submit");
  cfd.atp_sbp = AsnFind ("Seq-submit.sub");
  cfd.atp_bss = AsnFind ("Bioseq-set");
  cfd.atp_bsss = AsnFind ("Bioseq-set.seq-set");
  cfd.atp_se = AsnFind ("Bioseq-set.seq-set.E");
  cfd.atp_bsc = AsnFind ("Bioseq-set.class");
  cfd.bssp_atp = AsnLinkType (NULL, cfd.atp_bss);

  logfile = (CharPtr) myargs [L_argLogFile].strvalue;
  if (StringDoesHaveText (logfile)) {
    cfd.logfp = FileOpen (logfile, "w");
  }

  if (remote) {
#ifdef INTERNAL_NCBI_CLEANASN
    if (! PUBSEQBioseqFetchEnable ("cleanasn", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
    if (hup) {
      SmartFetchEnable ();
      TPASmartFetchEnable ();
      HUPFetchEnable ();
    }
#else
    PubSeqFetchEnableEx (TRUE, TRUE, TRUE, TRUE, TRUE, 0);
#endif
  }

  if (remote || cfd.pub) {
    PubMedFetchEnable ();
  }

  /*
  if (cfd.logfp != NULL && StringChr (cfd.report, 'c') != NULL) {
    fprintf (cfd.logfp, "FILE\t\tGENBANK\t\t\tEMBL\t\t\tDDBJ\t\t\tREFSEQ\t\t\tOTHER\n");
    fprintf (cfd.logfp, "\tREC\tNUC\tPRT\tREC\tNUC\tPRT\tREC\tNUC\tPRT\tREC\tNUC\tPRT\tREC\tNUC\tPRT\n");
    fflush (cfd.logfp);
  }
  */

  starttime = GetSecs ();

  if (StringDoesHaveText (directory)) {
    if (StringHasNoText (cfd.report) && StringCmp (directory, results) == 0) {
      Message (MSG_POSTERR, "-r results path must be different than -p data path");
      if (cfd.logfp != NULL) {
        fprintf (cfd.logfp, "-r results path must be different than -p data path\n");
      }
    } else {

      cfd.firstfile = (CharPtr) myargs [j_argFirstFile].strvalue;
      cfd.lastfile = (CharPtr) myargs [k_argLastFile].strvalue;

      cfd.results = results;

      DirExplore (directory, filter, suffix, FALSE, CleanupOneRecord, (Pointer) &cfd);
    }

  } else if (StringDoesHaveText (accnfile)) {
    if (StringHasNoText (cfd.report) && StringHasNoText (results)) {
      Message (MSG_POSTERR, "-r results path must be set when -A accession list is used");
      if (cfd.logfp != NULL) {
        fprintf (cfd.logfp, "-r results path must be set when -A accession list is used\n");
      }
    } else {

      cfd.results = results;

      ProcessAccessionList (accnfile, (Pointer) &cfd);
    }

  } else if (StringDoesHaveText (infile) && StringDoesHaveText (outfile)) {

    cfd.outfile = outfile;

    CleanupOneRecord (infile, (Pointer) &cfd);
  }

  if (StringChr (cfd.report, 'u') != NULL) {
    if (cfd.logfp != NULL) {
      fprintf (cfd.logfp, "\n\nTrying %ld Entrez Queries\n\n", (long) ValNodeLen (cfd.unpublist));
    }
    FinishUnpublishedList (&cfd);
  }

  stoptime = GetSecs ();
  runtime = stoptime - starttime;

  if (cfd.logfp != NULL) {
    if (StringChr (cfd.report, 'c') == NULL &&
        StringChr (cfd.report, 'j') == NULL &&
        StringChr (cfd.report, 'w') == NULL) {
      fprintf (cfd.logfp, "\nFinished in %ld seconds\n", (long) runtime);
      fprintf (cfd.logfp, "Cumulative counts ");
      fprintf (cfd.logfp, "- %9ld RECS", (long) cfd.cumcounts.recs);
      fprintf (cfd.logfp, ", %9ld NUCS", (long) cfd.cumcounts.nucs);
      fprintf (cfd.logfp, ", %9ld PRTS", (long) cfd.cumcounts.prts);
      fprintf (cfd.logfp, ", %9ld OKAY", (long) cfd.cumcounts.okay);
      fprintf (cfd.logfp, ", %9ld NORM", (long) cfd.cumcounts.norm);
      fprintf (cfd.logfp, ", %9ld CLNR", (long) cfd.cumcounts.clnr);
      fprintf (cfd.logfp, ", %9ld OTHR", (long) cfd.cumcounts.othr);
      fprintf (cfd.logfp, ", %9ld MODR", (long) cfd.cumcounts.modr);
      fprintf (cfd.logfp, ", %9ld SLOC", (long) cfd.cumcounts.sloc);
      fprintf (cfd.logfp, ", %9ld PUBL", (long) cfd.cumcounts.publ);
      fprintf (cfd.logfp, ", %9ld AUTH", (long) cfd.cumcounts.auth);
      fprintf (cfd.logfp, ", %9ld SORT", (long) cfd.cumcounts.sort);
      fprintf (cfd.logfp, ", %9ld BSEC", (long) cfd.cumcounts.bsec);
      fprintf (cfd.logfp, ", %9ld GBBK", (long) cfd.cumcounts.gbbk);
      fprintf (cfd.logfp, ", %9ld TITL", (long) cfd.cumcounts.titl);
      fprintf (cfd.logfp, ", %9ld PACK", (long) cfd.cumcounts.pack);
      fprintf (cfd.logfp, ", %9ld MOVE", (long) cfd.cumcounts.move);
      fprintf (cfd.logfp, ", %9ld SSEC", (long) cfd.cumcounts.ssec);
      fprintf (cfd.logfp, "\n");
      fflush (cfd.logfp);
    }
    FileClose (cfd.logfp);
  }

  if (remote || cfd.pub) {
    PubMedFetchDisable ();
  }

  if (remote) {
#ifdef INTERNAL_NCBI_CLEANASN
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
  }

  return 0;
}

