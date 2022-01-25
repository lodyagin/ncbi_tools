/*   tbl2asn.c
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
* File Name:  tbl2asn.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   5/5/00
*
* $Revision: 6.4 $
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
#include <edutil.h>
#include <seqport.h>
#include <gather.h>
#include <sqnutils.h>
#include <subutil.h>
#include <toasn3.h>
#include <valid.h>
#include <asn2ff.h>
#include <explore.h>

static FILE* OpenOneFile (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix
)

{
  Char  file [FILENAME_MAX], path [PATH_MAX];

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  return FileOpen (path, "r");
}

static void WriteOneFile (
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  SeqEntryPtr sep,
  SubmitBlockPtr sbp
)

{
  AsnIoPtr   aip;
  Char       file [FILENAME_MAX], path [PATH_MAX];
  SeqSubmit  ssb;

  MemSet ((Pointer) &ssb, 0, sizeof (SeqSubmit));
  ssb.sub = sbp;
  ssb.datatype = 1;
  ssb.data = (Pointer) sep;

  StringNCpy_0 (path, results, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  aip = AsnIoOpen (path, "w");
  if (aip == NULL) return;

  if (sbp != NULL) {
    SeqSubmitAsnWrite (&ssb, aip, NULL);
  } else {
    SeqEntryAsnWrite (sep, aip, NULL);
  }

  AsnIoFlush (aip);
  AsnIoClose (aip);
}

static void ValidateOneFile (
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  SeqEntryPtr sep
)

{
  Char            file [FILENAME_MAX], path [PATH_MAX];
  ErrSev          oldErrSev;
  ValidStructPtr  vsp;

  StringNCpy_0 (path, results, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  ErrSetOptFlags (EO_LOGTO_USRFILE);
  ErrSetLogfile (path, ELOG_APPEND | ELOG_NOCREATE);

  vsp = ValidStructNew ();
  if (vsp != NULL) {
    vsp->useSeqMgrIndexes = TRUE;
    vsp->suppressContext = TRUE;
    oldErrSev = ErrSetMessageLevel (SEV_NONE);
    ValidateSeqEntry (sep, vsp);
    ValidStructFree (vsp);
    ErrSetMessageLevel (oldErrSev);
  }

  ErrSetLogfile (NULL, ELOG_APPEND | ELOG_NOCREATE);
  ErrClearOptFlags  (EO_LOGTO_USRFILE);
}

static void FlatfileOneFile (
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  SeqEntryPtr sep
)

{
  Char    file [FILENAME_MAX], path [PATH_MAX];
  FILE    *fp;
  ErrSev  oldErrSev;

  StringNCpy_0 (path, results, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  fp = FileOpen (path, "w");
  if (fp == NULL) return;

  oldErrSev = ErrSetMessageLevel (SEV_MAX);
  SeqEntryToFlat (sep, fp, GENBANK_FMT, SEQUIN_MODE);
  ErrSetMessageLevel (oldErrSev);

  FileClose (fp);
}

/* for full-length cDNAs, allow automatic annotation of largest internal ORF */

typedef struct orfdata {
  Int4     curlen [6], bestlen [6], currstart [6], beststart [6], sublen [6];
  Boolean  inorf [6], altstart;
} OrfData, PNTR OrfDataPtr;

static void LIBCALLBACK LookForOrfs (
  Int4 position,
  Char residue,
  Boolean atgStart,
  Boolean altStart,
  Boolean orfStop,
  Int2 frame,
  Uint1 strand,
  Pointer userdata
)

{
  Int2        idx;
  OrfDataPtr  odp;

  odp = (OrfDataPtr) userdata;
  if (strand == Seq_strand_plus) {

    /* top strand */

    idx = frame;
    if (odp->inorf [idx]) {
      if (orfStop) {
        odp->inorf [idx] = FALSE;
        if (odp->curlen [idx] > odp->bestlen [idx]) {
          odp->bestlen [idx] = odp->curlen [idx];
          odp->beststart [idx] = odp->currstart [idx];
        }
      } else {
        (odp->curlen [idx])++;
      }
    } else if (atgStart || (altStart && odp->altstart)) {
      odp->inorf [idx] = TRUE;
      odp->curlen [idx] = 1;
      odp->currstart [idx] = position - frame;
    }
  } else {

    /* bottom strand */

    idx = frame + 3;
    if (orfStop) {
      odp->curlen [idx] = 0;
      odp->sublen [idx] = 0;
      odp->currstart [idx] = position - frame;
    } else if (atgStart || (altStart && odp->altstart)) {
      (odp->sublen [idx])++;
      odp->curlen [idx] = odp->sublen [idx];
      if (odp->curlen [idx] > odp->bestlen [idx]) {
        odp->bestlen [idx] = odp->curlen [idx];
        odp->beststart [idx] = odp->currstart [idx];
      }
    } else {
      (odp->sublen [idx])++;
    }
  }
}

static void AnnotateBestOrf (
  BioseqPtr bsp,
  Int2 genCode,
  Boolean altstart
  
)

{
  CdRegionPtr     crp;
  Int2            i, best, idx;
  OrfData         od;
  ProtRefPtr      prp;
  SeqFeatPtr      sfp;
  SeqInt          sint;
  TransTablePtr   tbl;
  ValNode         vn;
  SeqFeatXrefPtr  xref;

  if (bsp == NULL) return;
  for (i = 0; i < 6; i++) {
    od.curlen [i] = INT4_MIN;
    od.bestlen [i] = 0;
    od.currstart [i] = 0;
    od.beststart [i] = 0;
    od.sublen [i] = INT4_MIN;
    od.inorf [i] = FALSE;
  }
  od.altstart = altstart;

  /* use simultaneous 6-frame translation finite state machine */

  tbl = TransTableNew (genCode);
  if (tbl != NULL) {
    TransTableProcessBioseq (tbl, LookForOrfs, (Pointer) &od, bsp);
  }
  TransTableFree (tbl);
  best = -1;
  idx = -1;
  for (i = 0; i < 6; i++) {
    if (od.bestlen [i] > best) {
      best = od.bestlen [i];
      idx = i;
    }
  }
  if (idx == -1) return;

  /* make feature location on largest ORF */

  if (idx < 3) {
    MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
    sint.from = od.beststart [idx] + idx;
    sint.to = sint.from + (od.bestlen [idx]) * 3 + 2;
    sint.id = SeqIdFindBest (bsp->id, 0);
    sint.strand = Seq_strand_plus;
    vn.choice = SEQLOC_INT;
    vn.extended = 0;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
  } else {
    MemSet ((Pointer) &sint, 0, sizeof (SeqInt));
    sint.from = od.beststart [idx] + idx - 3;
    sint.to = sint.from + (od.bestlen [idx]) * 3 + 2;
    sint.id = SeqIdFindBest (bsp->id, 0);
    sint.strand = Seq_strand_minus;
    vn.choice = SEQLOC_INT;
    vn.extended = 0;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
  }

  /* make CDS feature with unknown product */

  sfp = CreateNewFeatureOnBioseq (bsp, SEQFEAT_CDREGION, &vn);
  if (sfp == NULL) return;
  crp = CreateNewCdRgn (1, FALSE, genCode);
  if (crp == NULL) return;
  crp->frame = 1;
  sfp->data.value.ptrvalue = (Pointer) crp;

  prp = ProtRefNew ();
  if (prp == NULL) return;
  xref = SeqFeatXrefNew ();
  if (xref == NULL) return;
  xref->data.choice = SEQFEAT_PROT;
  xref->data.value.ptrvalue = (Pointer) prp;
  xref->next = sfp->xref;
  sfp->xref = xref;
  prp->name = ValNodeCopyStr (NULL, 0, "unknown");
}

/* change all feature IDs to entered accession */

static void PromoteSeqId (SeqIdPtr sip, Pointer userdata)

{
  SeqIdPtr  bestid, newid, oldid;

  bestid = (SeqIdPtr) userdata;

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

static void CorrectFeatureSeqIds (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  VisitSeqIdsInSeqLoc (sfp->location, userdata, PromoteSeqId);
}

/* source information for several common organisms sequenced by genome centers */

typedef struct orgstuff {
  CharPtr  taxname;
  CharPtr  lineage;
  CharPtr  division;
  Uint1    gcode;
  Uint1    mgcode;
  Int4     taxID;
} OrgStuff, PNTR OrfStuffPtr;

static OrgStuff commonOrgStuff [] = {
  {
    "Saccharomyces cerevisiae",
    "Eukaryota; Fungi; Ascomycota; Hemiascomycetes; Saccharomycetales; Saccharomycetaceae; Saccharomyces",
    "PLN", 1, 3, 4932
  },
  {
    "Drosophila melanogaster",
    "Eukaryota; Metazoa; Arthropoda; Tracheata; Hexapoda; Insecta; Pterygota; Diptera; Brachycera; Muscomorpha; Ephydroidea; Drosophilidae; Drosophila",
    "INV", 1, 5, 7227
  },
  {
    "Homo sapiens",
    "Eukaryota; Metazoa; Chordata; Vertebrata; Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Homo",
    "PRI", 1, 2, 9606
  },
  {
    "Escherichia coli",
    "Bacteria; Proteobacteria; gamma subdivision; Enterobacteriaceae; Escherichia",
    "BCT", 11, 0, 562
  },
  {
    "Helicobacter pylori",
    "Bacteria; Proteobacteria; epsilon subdivision; Helicobacter group; Helicobacter",
    "BCT", 11, 0, 210
  },
  {
    NULL, NULL, NULL, 0, 0, 0
  }
};

static Boolean HasTaxon (OrgRefPtr orp)

{
  ValNodePtr  db;
  DbtagPtr    dbt;

  if (orp == FALSE) return FALSE;
  for (db = orp->db; db != NULL; db = db->next) {
    dbt = (DbtagPtr) db->data.ptrvalue;
    if (dbt != NULL && dbt->db != NULL &&
        StringICmp (dbt->db, "taxon") == 0) return TRUE;
  }
  return FALSE;
}

static void AddMissingSourceInfo (BioSourcePtr biop)

{
  ValNodePtr   db;
  DbtagPtr     dbt;
  Int2         idx;
  ObjectIdPtr  oip;
  OrgNamePtr   onp;
  OrgRefPtr    orp;

  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;

  /* look for entry of organisms in commonOrgStuff table */

  for (idx = 0; commonOrgStuff [idx].taxname != NULL; idx++) {
    if (StringICmp (orp->taxname, commonOrgStuff [idx].taxname) == 0) {
      if (onp->gcode == 0) {
        onp->gcode = commonOrgStuff [idx].gcode;
      }
      if (onp->mgcode == 0) {
        onp->mgcode = commonOrgStuff [idx].mgcode;
      }
      if (StringHasNoText (onp->div)) {
        onp->div = StringSave (commonOrgStuff [idx].division);
      }
      if (StringHasNoText (onp->lineage)) {
        onp->lineage = StringSave (commonOrgStuff [idx].lineage);
      }
      if (! HasTaxon (orp)) {
        db = ValNodeNew (NULL);
        if (db != NULL) {
          dbt = DbtagNew ();
          if (dbt != NULL) {
            oip = ObjectIdNew ();
            if (oip != NULL) {
              oip->id = commonOrgStuff [idx].taxID;
              dbt->db = StringSave ("taxon");
              dbt->tag = oip;
              db->data.ptrvalue = (Pointer) dbt;
              orp->db = db;
            }
          }
        }
      }
    }
  }
}

static BioseqPtr SqnGetBioseqGivenSeqLoc (SeqLocPtr slp, Uint2 entityID)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    tmp;

  if (slp == NULL) return NULL;
  bsp = NULL;
  sip = SeqLocId (slp);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
  } else {
    tmp = SeqLocFindNext (slp, NULL);
    if (tmp != NULL) {
      sip = SeqLocId (tmp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          sep = SeqMgrGetSeqEntryForData (bsp);
          entityID = ObjMgrGetEntityIDForChoice (sep);
          bsp = GetBioseqGivenSeqLoc (slp, entityID);
        }
      }
    }
  }
  return bsp;
}

static BioseqPtr GetBioseqReferencedByAnnot (SeqAnnotPtr sap, Uint2 entityID)

{
  SeqAlignPtr   align;
  BioseqPtr     bsp;
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  SeqFeatPtr    feat;
  SeqGraphPtr   graph;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  StdSegPtr     ssp;
  SeqLocPtr     tloc;

  if (sap == NULL) return NULL;
  switch (sap->type) {
    case 1 :
      feat = (SeqFeatPtr) sap->data;
      while (feat != NULL) {
        slp = feat->location;
        if (slp != NULL) {
          bsp = SqnGetBioseqGivenSeqLoc (slp, entityID);
          if (bsp != NULL) return bsp;
        }
        feat = feat->next;
      }
      break;
    case 2 :
      align = (SeqAlignPtr) sap->data;
      while (align != NULL) {
        if (align->segtype == 1) {
          ddp = (DenseDiagPtr) align->segs;
          if (ddp != NULL) {
            for (sip = ddp->id; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 2) {
          dsp = (DenseSegPtr) align->segs;
          if (dsp != NULL) {
            for (sip = dsp->ids; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 3) {
          ssp = (StdSegPtr) align->segs;
          if (ssp != NULL && ssp->loc != NULL) {
            for (tloc = ssp->loc; tloc != NULL; tloc = tloc->next) {
              bsp = BioseqFind (SeqLocId (tloc));
              if (bsp != NULL) return bsp;
            }
          }
        }
        align = align->next;
      }
      break;
    case 3 :
      graph = (SeqGraphPtr) sap->data;
      while (graph != NULL) {
        slp = graph->loc;
        if (slp != NULL) {
          bsp = SqnGetBioseqGivenSeqLoc (slp, entityID);
          if (bsp != NULL) return bsp;
        }
        graph = graph->next;
      }
      break;
    default :
      break;
  }
  return NULL;
}

static BioseqPtr AttachSeqAnnotEntity (Uint2 entityID, SeqAnnotPtr sap)

{
  BioseqPtr      bsp;
  Int2           genCode;
  SeqEntryPtr    oldscope;
  OMProcControl  ompc;
  SeqEntryPtr    sep;
  SeqFeatPtr     sfp = NULL;

  if (sap == NULL) return NULL;
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
      sep = GetBestTopParentForData (entityID, bsp);
      genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
      SetEmptyGeneticCodes (sap, genCode);
    }
    MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
    ompc.input_entityID = entityID;
    ompc.input_itemID = GetItemIDGivenPointer (entityID, OBJ_BIOSEQ, (Pointer) bsp);
    ompc.input_itemtype = OBJ_BIOSEQ;
    ompc.output_itemtype = OBJ_SEQANNOT;
    ompc.output_data = (Pointer) sap;
    if (! AttachDataForProc (&ompc, FALSE)) {
      Message (MSG_POSTERR, "AttachSeqAnnotEntity failed");
    } else if (sfp != NULL) {
      PromoteXrefs (sfp, bsp, entityID);
    }
  } else {
    Message (MSG_POSTERR, "Feature table identifiers do not match record");
  }
  return bsp;
}

static void ProcessOneRecord (
  SubmitBlockPtr sbp,
  BioSourcePtr src,
  CharPtr directory,
  CharPtr results,
  CharPtr base,
  CharPtr suffix,
  CharPtr accn,
  CharPtr organism,
  Boolean findorf,
  Boolean altstart,
  Boolean validate,
  Boolean flatfile
)

{
  BioSourcePtr  biop = NULL;
  BioseqPtr     bsp = NULL;
  Pointer       dataptr;
  Uint2         datatype, entityID;
  FILE          *fp;
  Int2          genCode;
  MolInfoPtr    mip;
  SeqAnnotPtr   sap;
  SeqEntryPtr   nsep, sep;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SqnTagPtr     stp;
  CharPtr       ttl;
  ValNodePtr    vnp;

  fp = OpenOneFile (directory, base, suffix);
  if (fp == NULL) return;

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);

  if (dataptr == NULL) return;

  sep = GetTopSeqEntryForEntityID (entityID);
  nsep = FindNucSeqEntry (sep);
  if (nsep != NULL && IS_Bioseq (nsep)) {
    bsp = (BioseqPtr) nsep->data.ptrvalue;
  }
  if (bsp == NULL) {
    ObjMgrFreeByEntityID (entityID);
    return;
  }

  if (bsp->mol == Seq_mol_na) {
    bsp->mol = Seq_mol_dna;
  }

  if (src != NULL) {
    src = AsnIoMemCopy ((Pointer) src,
                        (AsnReadFunc) BioSourceAsnRead,
                        (AsnWriteFunc) BioSourceAsnWrite);
  }

  vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
  if (vnp != NULL) {
    ttl = (CharPtr) vnp->data.ptrvalue;
    if (ttl != NULL) {
      stp = SqnTagParse (ttl);
      if (stp != NULL) {
        biop = ParseTitleIntoBioSource (stp, organism, src);
        ParseTitleIntoBioseq (stp, bsp);
      }
      SqnTagFree (stp);
    }
    ValNodeFreeData (vnp);
  }
  if (biop == NULL) {
    biop = ParseTitleIntoBioSource (NULL, organism, src);
  }
  if (biop != NULL) {
    SeqDescrAddPointer (&(bsp->descr), Seq_descr_source, (Pointer) biop);
    AddMissingSourceInfo (biop);
  }

  if (BioseqGetSeqDescr (bsp, Seq_descr_molinfo, NULL) == NULL) {
    mip = MolInfoNew ();
    if (mip != NULL) {
      mip->biomol = MOLECULE_TYPE_GENOMIC;
      SeqDescrAddPointer (&(bsp->descr), Seq_descr_molinfo, (Pointer) mip);
    }
  }

  sep = GetBestTopParentForData (entityID, bsp);
  genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);

  if (findorf) {
    AnnotateBestOrf (bsp, genCode, altstart);
  }

  fp = OpenOneFile (directory, base, ".tbl");
  if (fp != NULL) {

    while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
      if (datatype == OBJ_SEQANNOT) {

        sap = (SeqAnnotPtr) dataptr;
        bsp = AttachSeqAnnotEntity (entityID, sap);
        if (bsp != NULL) {

          /* if existing accession, coerce all SeqIds */

          if (! StringHasNoText (accn)) {
            sip = SeqIdFromAccession (accn, 0, NULL);
            if (sip != NULL) {
              bsp->id = SeqIdSetFree (bsp->id);
              bsp->id = sip;
              SeqMgrReplaceInBioseqIndex (bsp);
              VisitFeaturesOnBsp (bsp, (Pointer) bsp->id, CorrectFeatureSeqIds);
            }
          }

          /* for parsed in features or best ORF, promote CDS products to protein bioseq */

          for (sap = bsp->annot; sap != NULL; sap = sap->next) {
            if (sap->type == 1) {
              SetEmptyGeneticCodes (sap, genCode);
              sfp = (SeqFeatPtr) sap->data;
              PromoteXrefs (sfp, bsp, entityID);
            }
          }
        }

      } else {
        ObjMgrFree (datatype, dataptr);
      }
    }
    FileClose (fp);
  }

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    WriteOneFile (results, base, ".sqn", sep, sbp);
    if (validate || flatfile) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    if (validate) {
      Message (MSG_POST, "Validating %s\n", base);
      ValidateOneFile (results, base, ".val", sep);
    }
    if (flatfile) {
      Message (MSG_POST, "Flatfile %s\n", base);
      sep = FindNucSeqEntry (sep);
      FlatfileOneFile (results, base, ".gbf", sep);
    }
  }

  ObjMgrFreeByEntityID (entityID);
}

/* command-line argument list */

#define p_argInputPath  0
#define r_argOutputPath 1
#define f_argSingleFile 2
#define x_argSuffix     3
#define t_argTemplate   4
#define a_argAccession  5
#define n_argOrgName    6
#define c_argFindOrf    7
#define m_argAltStart   8
#define v_argValidate   9
#define b_argGenBank   10

Args myargs [] = {
  {"Path to files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Only this file", NULL, NULL, NULL,
    TRUE, 'f', ARG_FILE_IN, 0.0, 0, NULL},
  {"Suffix", ".fsa", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Template file", NULL, NULL, NULL,
    TRUE, 't', ARG_FILE_IN, 0.0, 0, NULL},
  {"Accession", NULL, NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Organism name", NULL, NULL, NULL,
    TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
  {"Annotate longest ORF", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Allow alternative starts", "F", NULL, NULL,
    TRUE, 'm', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Validate", "F", NULL, NULL,
    TRUE, 'v', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Generate GenBank file", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL}
};

Int2 Main (void)

{
  AsnIoPtr        aip;
  CharPtr         base, directory, results, suffix, accn, organism, ptr, tmplate;
  Boolean         altstart, findorf, flatfile, validate;
  ValNodePtr      head, vnp;
  SubmitBlockPtr  sbp = NULL;
  SeqEntryPtr     sep;
  BioSourcePtr    src = NULL;
  SeqSubmitPtr    ssp = NULL;

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
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
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }

  if (! GetArgs ("tbl2asn", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = directory;
  }
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  base = (CharPtr) myargs [f_argSingleFile].strvalue;
  tmplate = (CharPtr) myargs [t_argTemplate].strvalue;
  accn = (CharPtr) myargs [a_argAccession].strvalue;
  organism = (CharPtr) myargs [n_argOrgName].strvalue;
  findorf = (Boolean) myargs [c_argFindOrf].intvalue;
  altstart = (Boolean) myargs [m_argAltStart].intvalue;
  validate = (Boolean) myargs [v_argValidate].intvalue;
  flatfile = (Boolean) myargs [b_argGenBank].intvalue;

  if (StringHasNoText (base) && (! StringHasNoText (accn))) {
    Message (MSG_FATAL, "Accession can be entered only for a single record");
    return 1;
  }

  if (! StringHasNoText (tmplate)) {
    aip = AsnIoOpen (tmplate, "r");
    if (aip != NULL) {
      ssp = SeqSubmitAsnRead (aip, NULL);
      AsnIoClose (aip);
    }

    if (ssp == NULL) {
      Message (MSG_FATAL, "Unable to read required template file");
      return 1;
    }

    sbp = ssp->sub;
    if (sbp == NULL) {
      Message (MSG_FATAL, "Unable to read submit block within required template file");
      ssp = SeqSubmitFree (ssp);
      return 1;
    }

    if (sbp != NULL) {
      sbp->tool = MemFree (sbp->tool);
      sbp->tool = StringSave ("tbl2asn");
      sbp->hup = FALSE;
      sbp->reldate = DateFree (sbp->reldate);
    }
    if (ssp->datatype == 1) {
      sep = (SeqEntryPtr) ssp->data;
      if (sep != NULL) {
        SeqEntryToBioSource (sep, NULL, NULL, 0, &src);
      }
    }
  }

  if (! StringHasNoText (base)) {
    ptr = StringStr (base, suffix);
    if (ptr != NULL) {
      *ptr = '\0';
      ProcessOneRecord (sbp, src, directory, results, base, suffix, accn, organism, findorf, altstart, validate, flatfile);
    }
  } else {
    head = DirCatalog (directory);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == 0) {
        base = (CharPtr) vnp->data.ptrvalue;
        if (! StringHasNoText (base)) {
          ptr = StringStr (base, suffix);
          if (ptr != NULL) {
            *ptr = '\0';
            Message (MSG_POST, "Processing %s\n", base);
            ProcessOneRecord (sbp, src, directory, results, base, suffix, NULL, organism, findorf, altstart, validate, flatfile);
          }
        }
      }
    }
    ValNodeFreeData (head);
  }

  ssp = SeqSubmitFree (ssp);

  return 0;
}

