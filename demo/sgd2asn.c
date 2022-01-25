/*   sgd2asn.c
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
* File Name:  sgd2asn.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   8/20/99
*
* $Revision: 6.8 $
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

static Pointer ReadOneFile (CharPtr directory, CharPtr base, CharPtr suffix, Uint2Ptr datatype)

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

  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, NULL, FALSE, FALSE, TRUE, FALSE);

  FileClose (fp);

  return dataptr;
}

static void WriteOneFile (CharPtr directory, CharPtr base, CharPtr suffix, SeqEntryPtr sep)

{
  AsnIoPtr  aip;
  Char      file [FILENAME_MAX], path [PATH_MAX];

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s.%s", base, suffix);
  FileBuildPath (path, NULL, file);

  aip = AsnIoOpen (path, "w");
  if (aip == NULL) return;

  SeqEntryAsnWrite (sep, aip, NULL);

  AsnIoFlush (aip);
  AsnIoClose (aip);
}

static void ValidateOneFile (CharPtr directory, CharPtr base, CharPtr suffix, SeqEntryPtr sep)

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

static void FlatfileOneFile (CharPtr directory, CharPtr base, CharPtr suffix, SeqEntryPtr sep)

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

static BioSourcePtr ParseTitleIntoBioSource (SqnTagPtr stp, BioseqPtr bsp)

{
  BioSourcePtr  biop = NULL;
  ValNodePtr    db;
  DbtagPtr      dbt;
  ObjectIdPtr   oip;
  OrgModPtr     omp;
  OrgNamePtr    onp;
  OrgRefPtr     orp;
  SubSourcePtr  ssp;
  CharPtr       str;
  int           val;

  if (stp == NULL) return NULL;

  str = SqnTagFind (stp, "top");
  if (str != NULL && bsp != NULL) {
    if (StringICmp (str, "linear") == 0) {
      bsp->topology = TOPOLOGY_LINEAR;
    } else if (StringICmp (str, "circular") == 0) {
      bsp->topology = TOPOLOGY_CIRCULAR;
    }
  }
  str = SqnTagFind (stp, "org");
  if (str == NULL) return NULL;
   biop = BioSourceNew ();
  if (biop == NULL) return NULL;
  orp = OrgRefNew ();
  if (orp == NULL) return NULL;
  biop->org = orp;
  onp = OrgNameNew ();
  if (onp == NULL) return NULL;
  orp->orgname = onp;
  orp->taxname = StringSave (str);
  if (StringICmp (orp->taxname, "Saccharomyces cerevisiae") == 0) {
    onp->gcode = 1; /* standard */
    onp->mgcode = 3; /* yeast mitochondrial */
    onp->div = StringSave ("PLN");
    onp->lineage = StringSave ("Eukaryota; Fungi; Ascomycota; Hemiascomycetes; Saccharomycetales; Saccharomycetaceae; Saccharomyces");
    db = ValNodeNew (NULL);
    if (db != NULL) {
      dbt = DbtagNew ();
      if (dbt != NULL) {
        oip = ObjectIdNew ();
        if (oip != NULL) {
          oip->id = 4932;
          dbt->db = StringSave ("taxon");
          dbt->tag = oip;
          db->data.ptrvalue = (Pointer) dbt;
          orp->db = db;
        }
      }
    }
  } else if (StringICmp (orp->taxname, "Escherichia coli") == 0) {
    onp->gcode = 11; /* bacterial */
    onp->div = StringSave ("BCT");
    onp->lineage = StringSave ("Bacteria; Proteobacteria; gamma subdivision; Enterobacteriaceae; Escherichia");
    db = ValNodeNew (NULL);
    if (db != NULL) {
      dbt = DbtagNew ();
      if (dbt != NULL) {
        oip = ObjectIdNew ();
        if (oip != NULL) {
          oip->id = 562;
          dbt->db = StringSave ("taxon");
          dbt->tag = oip;
          db->data.ptrvalue = (Pointer) dbt;
          orp->db = db;
        }
      }
    }
  } else if (StringICmp (orp->taxname, "Helicobacter pylori") == 0) {
    onp->gcode = 11; /* bacterial */
    onp->div = StringSave ("BCT");
    onp->lineage = StringSave ("Bacteria; Proteobacteria; epsilon subdivision; Helicobacter group; Helicobacter");
    db = ValNodeNew (NULL);
    if (db != NULL) {
      dbt = DbtagNew ();
      if (dbt != NULL) {
        oip = ObjectIdNew ();
        if (oip != NULL) {
          oip->id = 210;
          dbt->db = StringSave ("taxon");
          dbt->tag = oip;
          db->data.ptrvalue = (Pointer) dbt;
          orp->db = db;
        }
      }
    }
  }

  str = SqnTagFind (stp, "gcode");
  if (str != NULL && sscanf (str, "%d", &val) == 1) {
    onp->gcode = (Uint1) val; /* cytoplasmic */
  }

  str = SqnTagFind (stp, "mgcode");
  if (str != NULL && sscanf (str, "%d", &val) == 1) {
    onp->mgcode = (Uint1) val; /* mitochondrial */
  }

  str = SqnTagFind (stp, "location");
  if (str != NULL && StringICmp (str, "mitochondrion") == 0) {
    biop->genome = GENOME_mitochondrion;
  } else {
    biop->genome = GENOME_genomic;
  }

  str = SqnTagFind (stp, "strain");
  if (str != NULL) {
    omp = OrgModNew ();
    if (omp != NULL) {
      omp->subtype = ORGMOD_strain;
      omp->subname = StringSave (str);
      onp->mod = omp;
    }
  }

  str = SqnTagFind (stp, "chromosome");
  if (str != NULL) {
    ssp = SubSourceNew ();
    if (ssp != NULL) {
      ssp->subtype = SUBSRC_chromosome;
      ssp->name = StringSave (str);
      biop->subtype = ssp;
    }
  }

  return biop;
}

static void ProcessOneRecord (CharPtr directory, CharPtr base,
                              Boolean validate, Boolean flatfile)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  Uint2         datatype, entityID;
  Int2          genCode;
  MolInfoPtr    mip;
  SeqAnnotPtr   sap;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp;
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

  vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
  if (vnp != NULL) {
    ttl = (CharPtr) vnp->data.ptrvalue;
    if (ttl != NULL) {
      stp = SqnTagParse (ttl);
      if (stp != NULL) {
        biop = ParseTitleIntoBioSource (stp, bsp);
      }
      SqnTagFree (stp);
      if (biop != NULL) {
        SeqDescrAddPointer (&(bsp->descr), Seq_descr_source, (Pointer) biop);
      }
    }
    ValNodeFreeData (vnp);
  }

  mip = MolInfoNew ();
  if (mip != NULL) {
    mip->biomol = MOLECULE_TYPE_GENOMIC;
    SeqDescrAddPointer (&(bsp->descr), Seq_descr_molinfo, (Pointer) mip);
  }

  entityID = ObjMgrRegister (datatype, (Pointer) bsp);

  sap = (SeqAnnotPtr) ReadOneFile (directory, base, "tbl", &datatype);
  if (sap != NULL && datatype == OBJ_SEQANNOT && sap->type == 1) {
    sfp = (SeqFeatPtr) sap->data;
    if (sfp != NULL) {
      sap->next = bsp->annot;
      bsp->annot = sap;
      sep = GetBestTopParentForData (entityID, bsp);
      genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
      SetEmptyGeneticCodes (sap, genCode);
      PromoteXrefs (sfp, bsp, entityID);
    }
    sep = GetTopSeqEntryForEntityID (entityID);
    if (sep != NULL) {
      SeriousSeqEntryCleanup (sep, NULL, NULL);
      WriteOneFile (directory, base, "sqn", sep);
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
  } else {
    ObjMgrFree (datatype, (Pointer) sap);
  }

  ObjMgrFreeByEntityID (entityID);
}

Args myargs [] = {
  {"Path to files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Only this file", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"Validate", "F", NULL, NULL,
    TRUE, 'v', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Generate GenBank file", "F", NULL, NULL,
    TRUE, 'g', ARG_BOOLEAN, 0.0, 0, NULL}
};

Int2 Main (void)

{
  CharPtr     base, directory, ptr;
  Boolean     flatfile, validate;
  ValNodePtr  head, vnp;

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

  if (! GetArgs ("sgd2asn", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  directory = (CharPtr) myargs [0].strvalue;
  base = (CharPtr) myargs [1].strvalue;
  validate = (Boolean) myargs [2].intvalue;
  flatfile = (Boolean) myargs [3].intvalue;

  if (! StringHasNoText (base)) {
    ptr = StringStr (base, ".fsa");
    if (ptr != NULL) {
      *ptr = '\0';
      ProcessOneRecord (directory, base, validate, flatfile);
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
            ProcessOneRecord (directory, base, validate, flatfile);
          }
        }
      }
    }
    ValNodeFreeData (head);
  }

  return 0;
}

