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
* $Revision: 6.3 $
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

static Pointer ReadYeastFile (CharPtr directory, CharPtr base, CharPtr suffix, Uint2Ptr datatype)

{
  Pointer  dataptr;
  Char     file [FILENAME_MAX];
  FILE*    fp;
  Char     path [PATH_MAX];

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

static void WriteYeastFile (CharPtr directory, CharPtr base, CharPtr suffix, SeqEntryPtr sep)

{
  AsnIoPtr  aip;
  Char      file [FILENAME_MAX];
  Char      path [PATH_MAX];

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s.%s", base, suffix);
  FileBuildPath (path, NULL, file);

  aip = AsnIoOpen (path, "w");
  if (aip == NULL) return;

  SeqEntryAsnWrite (sep, aip, NULL);

  AsnIoFlush (aip);
  AsnIoClose (aip);
}

static BioSourcePtr ParseTitleIntoBioSource (SqnTagPtr stp)

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

  if (stp == NULL) return NULL;

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
          oip->id = 7227;
          dbt->db = StringSave ("taxon");
          dbt->tag = oip;
          db->data.ptrvalue = (Pointer) dbt;
          orp->db = db;
        }
      }
    }
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

static void ProcessYeastRecord (CharPtr directory, CharPtr base)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  Uint2         datatype;
  Uint2         entityID;
  Int2          genCode;
  MolInfoPtr    mip;
  SeqAnnotPtr   sap;
  SeqEntryPtr   sep;
  SeqFeatPtr    sfp;
  SqnTagPtr     stp;
  CharPtr       ttl;
  ValNodePtr    vnp;

  Message (MSG_POST, "%s\n", base);

  bsp = (BioseqPtr) ReadYeastFile (directory, base, "fsa", &datatype);
  if (bsp == NULL || datatype != OBJ_BIOSEQ) {
    ObjMgrFree (datatype, (Pointer) bsp);
    return;
  }

  vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
  if (vnp != NULL) {
    ttl = (CharPtr) vnp->data.ptrvalue;
    if (ttl != NULL) {
      stp = SqnTagParse (ttl);
      biop = ParseTitleIntoBioSource (stp);
      SqnTagFree (stp);
      if (biop != NULL) {
        ValNodeAddPointer (&(bsp->descr), Seq_descr_source, (Pointer) biop);
      }
      mip = MolInfoNew ();
      if (mip != NULL) {
        mip->biomol = 1;
        ValNodeAddPointer (&(bsp->descr), Seq_descr_molinfo, (Pointer) mip);
      }
    }
    ValNodeFreeData (vnp);
  }

  entityID = ObjMgrRegister (datatype, (Pointer) bsp);

  sap = (SeqAnnotPtr) ReadYeastFile (directory, base, "tbl", &datatype);
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
      WriteYeastFile (directory, base, "sqn", sep);
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
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL}
};

Int2 Main (void)

{
  Int2     i;
  CharPtr  yeastFiles [18] = {
    "chr01", "chr02", "chr03", "chr04",
    "chr05", "chr06", "chr07", "chr08",
    "chr09", "chr10", "chr11", "chr12",
    "chr13", "chr14", "chr15", "chr16",
    "chrmt", NULL
  };

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

  if (! GetArgs ("yeastfeat", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  if (! StringHasNoText (myargs [1].strvalue)) {
    ProcessYeastRecord (myargs [0].strvalue, myargs [1].strvalue);
  } else {
    for (i = 0; yeastFiles [i] != NULL; i++) {
      ProcessYeastRecord (myargs [0].strvalue, yeastFiles [i]);
    }
  }

  return 0;
}

