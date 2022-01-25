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

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
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

static Pointer ReadOneFile (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix,
  Uint2Ptr datatype
)

{
  Pointer  dataptr;
  Char     file [FILENAME_MAX], path [PATH_MAX];
  FILE*    fp;

  *datatype = 0;

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s.%s", base, suffix);
  FileBuildPath (path, NULL, file);

  fp = FileOpen (path, "r");
  if (fp == NULL) return NULL;

  /* ReadAsnFastaOrFlatFile can read ASN.1, FASTA, and feature table formats */

  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, NULL, FALSE, FALSE, TRUE, FALSE);

  FileClose (fp);

  return dataptr;
}

static void WriteOneFile (
  CharPtr directory,
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

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s.%s", base, suffix);
  FileBuildPath (path, NULL, file);

  aip = AsnIoOpen (path, "w");
  if (aip == NULL) return;

  SeqSubmitAsnWrite (&ssb, aip, NULL);

  AsnIoFlush (aip);
  AsnIoClose (aip);
}

static void ValidateOneFile (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix,
  SeqEntryPtr sep
)

{
  Char            file [FILENAME_MAX], path [PATH_MAX];
  ErrSev          oldErrSev;
  ValidStructPtr  vsp;

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s.%s", base, suffix);
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
  CharPtr directory,
  CharPtr base,
  CharPtr suffix,
  SeqEntryPtr sep
)

{
  Char    file [FILENAME_MAX], path [PATH_MAX];
  FILE    *fp;
  ErrSev  oldErrSev;

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s.%s", base, suffix);
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

static void ProcessOneRecord (
  SubmitBlockPtr sbp,
  BioSourcePtr src,
  CharPtr directory,
  CharPtr base,
  CharPtr accn,
  CharPtr organism,
  Boolean findorf,
  Boolean altstart,
  Boolean validate,
  Boolean flatfile
)

{
  BioSourcePtr  biop = NULL;
  BioseqPtr     bsp;
  Uint2         datatype, entityID;
  Int2          genCode;
  MolInfoPtr    mip;
  SeqAnnotPtr   sap;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SqnTagPtr     stp;
  CharPtr       ttl;
  ValNodePtr    vnp;

  Message (MSG_POST, "Processing %s\n", base);

  bsp = (BioseqPtr) ReadOneFile (directory, base, "fsa", &datatype);
  if (bsp == NULL || datatype != OBJ_BIOSEQ) {
    ObjMgrFree (datatype, (Pointer) bsp);
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

  mip = MolInfoNew ();
  if (mip != NULL) {
    mip->biomol = MOLECULE_TYPE_GENOMIC;
    SeqDescrAddPointer (&(bsp->descr), Seq_descr_molinfo, (Pointer) mip);
  }

  entityID = ObjMgrRegister (datatype, (Pointer) bsp);

  sep = GetBestTopParentForData (entityID, bsp);
  genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);

  if (findorf) {
    AnnotateBestOrf (bsp, genCode, altstart);
  }

  sap = (SeqAnnotPtr) ReadOneFile (directory, base, "tbl", &datatype);
  if (sap != NULL && datatype == OBJ_SEQANNOT && sap->type == 1) {
    sap->next = bsp->annot;
    bsp->annot = sap;
  } else {
    ObjMgrFree (datatype, (Pointer) sap);
  }

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

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    WriteOneFile (directory, base, "sqn", sep, sbp);
    if (validate || flatfile) {
      SeqMgrIndexFeatures (entityID, 0);
    }
    if (validate) {
      Message (MSG_POST, "Validating %s\n", base);
      ValidateOneFile (directory, base, "val", sep);
    }
    if (flatfile) {
      Message (MSG_POST, "Flatfile %s\n", base);
      sep = FindNucSeqEntry (sep);
      FlatfileOneFile (directory, base, "gbf", sep);
    }
  }

  ObjMgrFreeByEntityID (entityID);
}

/* command-line argument list */

#define i_argInputPath  0
#define i_argSingleFile 1
#define i_argTemplate   2
#define i_argAccession  3
#define i_argOrgName    4
#define i_argFindOrf    5
#define i_argAltStart   6
#define i_argValidate   7
#define i_argGenBank    8

Args myargs [] = {
  {"Path to files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Only this file", NULL, NULL, NULL,
    TRUE, 'f', ARG_FILE_IN, 0.0, 0, NULL},
  {"Template file", NULL, NULL, NULL,
    FALSE, 't', ARG_FILE_IN, 0.0, 0, NULL},
  {"Accession", NULL, NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Organism name", NULL, NULL, NULL,
    TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
  {"Annotate longest ORF", "F", NULL, NULL,
    TRUE, 'r', ARG_BOOLEAN, 0.0, 0, NULL},
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
  CharPtr         base, directory, accn, organism, ptr, tmplate;
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

  if (! GetArgs ("tbl2asn", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  directory = (CharPtr) myargs [i_argInputPath].strvalue;
  base = (CharPtr) myargs [i_argSingleFile].strvalue;
  tmplate = (CharPtr) myargs [i_argTemplate].strvalue;
  accn = (CharPtr) myargs [i_argAccession].strvalue;
  organism = (CharPtr) myargs [i_argOrgName].strvalue;
  findorf = (Boolean) myargs [i_argFindOrf].intvalue;
  altstart = (Boolean) myargs [i_argAltStart].intvalue;
  validate = (Boolean) myargs [i_argValidate].intvalue;
  flatfile = (Boolean) myargs [i_argGenBank].intvalue;

  if (StringHasNoText (base) && (! StringHasNoText (accn))) {
    Message (MSG_FATAL, "Accession can be entered only for a single record");
    return 1;
  }

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

  if (! StringHasNoText (base)) {
    ptr = StringStr (base, ".fsa");
    if (ptr != NULL) {
      *ptr = '\0';
      ProcessOneRecord (sbp, src, directory, base, accn, organism, findorf, altstart, validate, flatfile);
    }
  } else {
    head = DirCatalog (directory);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == 0) {
        base = (CharPtr) vnp->data.ptrvalue;
        if (! StringHasNoText (base)) {
          ptr = StringStr (base, ".fsa");
          if (ptr != NULL) {
            *ptr = '\0';
            ProcessOneRecord (sbp, src, directory, base, NULL, organism, findorf, altstart, validate, flatfile);
          }
        }
      }
    }
    ValNodeFreeData (head);
  }

  ssp = SeqSubmitFree (ssp);

  return 0;
}

