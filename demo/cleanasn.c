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
* $Revision: 6.10 $
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
#include <sqnutils.h>
#include <explore.h>
#include <toasn3.h>
#include <pmfapi.h>
#include <tax3api.h>
#ifdef INTERNAL_NCBI_CLEANASN
#include <accpubseq.h>
#endif

#define CLEANASN_APP_VER "1.4"

CharPtr CLEANASN_APPLICATION = CLEANASN_APP_VER;

typedef struct cleanflags {
  Boolean       batch;
  Boolean       binary;
  Boolean       compressed;
  Int2          type;
  CharPtr       results;
  CharPtr       outfile;
  CharPtr       clean;
  CharPtr       link;
  CharPtr       feat;
  Boolean       taxon;
  AsnModulePtr  amp;
  AsnTypePtr    atp_bss;
  AsnTypePtr    atp_bsss;
  AsnTypePtr    atp_se;
  AsnTypePtr    atp_bsc;
  AsnTypePtr    bssp_atp;
  BioseqSet     bss;
  FILE          *logfp;
} CleanFlagData, PNTR CleanFlagPtr;

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

typedef struct dummysmfedata {
  Int4  max;
  Int4  num_at_max;
} DummySmfeData, PNTR DummySmfePtr;

static Boolean LIBCALLBACK CADummySMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)


{
  DummySmfePtr  dsp;
  Int4          len;

  if (sfp == NULL || context == NULL) return TRUE;
  dsp = context->userdata;
  if (dsp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < dsp->max) {
    dsp->max = len;
    dsp->num_at_max = 1;
  } else if (len == dsp->max) {
    (dsp->num_at_max)++;
  }

  return TRUE;
}

static void RemoveUnnecGeneXref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Int2                 count;
  SeqFeatXrefPtr       curr, next;
  DummySmfeData        dsd;
  SeqMgrFeatContext    fcontext;
  SeqFeatXrefPtr PNTR  last;
  GeneRefPtr           grp, grpx;
  SeqFeatPtr           sfpx;
  CharPtr              syn1, syn2;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
  sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE) return;
  grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  if (grpx == NULL) return;

  if ((StringDoesHaveText (grp->locus)) &&
       (StringDoesHaveText (grpx->locus))) {
    if ((StringICmp (grp->locus, grpx->locus) != 0)) return;
  } else if (StringDoesHaveText (grp->locus_tag) &&
             StringDoesHaveText (grpx->locus_tag)) {
    if ((StringICmp (grp->locus_tag, grpx->locus_tag) != 0)) return;
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if ((StringDoesHaveText (syn1)) && (StringDoesHaveText (syn2))) {
      if ((StringICmp (syn1, syn2) != 0)) return;
    }
  }

  MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
  dsd.max = INT4_MAX;
  dsd.num_at_max = 0;
  count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE,
                                           NULL, 0, LOCATION_SUBSET,
                                           (Pointer) &dsd, CADummySMFEProc);

  if (dsd.num_at_max < 2) {
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

static void DoCleeanup (
  SeqEntryPtr sep,
  Uint2 entityID,
  CleanFlagPtr cfp
)

{
  if (sep == NULL || cfp == NULL) return;

  if (StringChr (cfp->clean, 'b') != NULL) {
	BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->clean, 's') != NULL) {
	SeriousSeqEntryCleanup (sep, NULL, NULL);
  }

  if (cfp->taxon) {
	Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
  }

  if (StringChr (cfp->link, 'o') != NULL) {
	SeqMgrIndexFeatures (entityID, 0);
	LinkCDSmRNAbyOverlap (sep);
  }
  if (StringChr (cfp->link, 'p') != NULL) {
	SeqMgrIndexFeatures (entityID, 0);
	LinkCDSmRNAbyProduct (sep);
  }
  if (StringChr (cfp->link, 'r') != NULL) {
	SeqMgrIndexFeatures (entityID, 0);
	ReassignFeatureIDs (sep);
  }
  if (StringChr (cfp->link, 'c') != NULL) {
	ClearFeatureIDs (sep);
  }

  if (StringChr (cfp->feat, 'u') != NULL) {
	VisitFeaturesInSep (sep, NULL, RemoveFeatUser);
  }
  if (StringChr (cfp->feat, 'd') != NULL) {
	VisitFeaturesInSep (sep, NULL, RemoveFeatDbxref);
  }
  if (StringChr (cfp->feat, 'r') != NULL) {
	SeqMgrIndexFeatures (entityID, 0);
	VisitFeaturesInSep (sep, NULL, RemoveUnnecGeneXref);
  }
}

static void CleanupSingleRecord (
  CharPtr filename,
  CleanFlagPtr cfp
)

{
  AsnIoPtr      aip, aop;
  AsnTypePtr    atp = NULL;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype, entityID = 0;
  FILE          *fp;
  Char          path [PATH_MAX];
  CharPtr       ptr;
  SeqEntryPtr   sep;

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

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
        datatype = OBJ_SEQSUB;
        break;
      default :
        break;
    }

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
      if (sep != NULL && StringDoesHaveText (path)) {

        DoCleeanup (sep, entityID, cfp);

        aop = AsnIoOpen (path, "w");
        if (aop != NULL) {
          if (datatype == OBJ_SEQSUB) {
            SeqSubmitAsnWrite ((SeqSubmitPtr) dataptr, aop, NULL);
          } else {
            SeqEntryAsnWrite (sep, aop, NULL);
          }
          AsnIoFlush (aop);
          AsnIoClose (aop);
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
  AsnIoPtr     aip, aop;
  AsnTypePtr   atp;
  DataVal      av;
  BioseqPtr    bsp;
  Char         buf [41];
  Char         cmmd [256];
  Uint2        entityID;
  FILE         *fp;
  SeqEntryPtr  fsep;
  size_t       len;
  Char         longest [41];
  Int4         numrecords;
  Char         path [PATH_MAX];
  CharPtr      ptr;
  SeqEntryPtr  sep;
  time_t       starttime, stoptime, worsttime;
#ifdef OS_UNIX
  CharPtr      gzcatprog;
  int          ret;
  Boolean      usedPopen = FALSE;
#endif

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

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
  if (StringHasNoText (path)) return;

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
    fprintf (cfp->logfp, "%s\n\n", filename);
    fflush (cfp->logfp);
  }

  longest [0] = '\0';
  worsttime = 0;
  numrecords = 0;

  aop = AsnIoOpen (path, cfp->binary? "wb" : "w");
  if (aop != NULL) {

    AsnOpenStruct (aop, cfp->bssp_atp, (Pointer) &(cfp->bss));
    av.intvalue = 7;
    AsnWrite (aop, cfp->atp_bsc, &av);
    AsnOpenStruct (aop, cfp->atp_bsss, (Pointer) &(cfp->bss.seq_set));

    atp = cfp->atp_bss;

    while ((atp = AsnReadId (aip, cfp->amp, atp)) != NULL) {
      if (atp == cfp->atp_se) {

        sep = SeqEntryAsnRead (aip, atp);
        if (sep != NULL) {

          entityID = ObjMgrGetEntityIDForChoice (sep);

          fsep = FindNthBioseq (sep, 1);
          if (fsep != NULL && fsep->choice == 1) {
            bsp = (BioseqPtr) fsep->data.ptrvalue;
            if (bsp != NULL) {
              SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
              if (cfp->logfp != NULL) {
                fprintf (cfp->logfp, "%s\n", buf);
                fflush (cfp->logfp);
              }
            }
          }

          starttime = GetSecs ();
          DoCleeanup (sep, entityID, cfp);
          stoptime = GetSecs ();

          if (stoptime - starttime > worsttime) {
            worsttime = stoptime - starttime;
            StringCpy (longest, buf);
          }
          numrecords++;

          SeqEntryAsnWrite (sep, aop, cfp->atp_se);

          ObjMgrFreeByEntityID (entityID);
        }

      } else {

        AsnReadVal (aip, atp, NULL);
      }
    }

    AsnCloseStruct (aop, cfp->atp_bsss, (Pointer) &(cfp->bss.seq_set));
    AsnCloseStruct (aop, cfp->bssp_atp, (Pointer) &(cfp->bss));
  }

  AsnIoClose (aop);
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
  if (cfp->logfp != NULL && (! StringHasNoText (longest))) {
    fprintf (cfp->logfp, "Longest processing time %ld seconds on %s\n",
             (long) worsttime, longest);
    fprintf (cfp->logfp, "Total number of records %ld\n", (long) numrecords);
    fflush (cfp->logfp);
  }
}

static void CleanupOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;

  if (StringHasNoText (filename)) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  if (cfp->batch) {
    CleanupMultipleRecord (filename, cfp);
  } else {
    CleanupSingleRecord (filename, cfp);
  }
}

/* Args structure contains command-line arguments */

#define p_argInputPath     0
#define r_argOutputPath    1
#define i_argInputFile     2
#define o_argOutputFile    3
#define f_argFilter        4
#define x_argSuffix        5
#define a_argType          6
#define b_argBinary        7
#define c_argCompressed    8
#define l_argLogFile       9
#define R_argRemote       10
#define K_argClean        11
#define N_argLink         12
#define F_argFeat         13
#define T_argTaxonLookup  14

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
  {"ASN.1 Type\n"
   "      a Any\n"
   "      e Seq-entry\n"
   "      b Bioseq\n"
   "      s Bioseq-set\n"
   "      m Seq-submit\n"
   "      t Batch Processing", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log fFile", NULL, NULL, NULL,
    TRUE, 'l', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Cleanup\n"
   "      b BasicSeqEntryCleanup\n"
   "      s SeriousSeqEntryCleanup", NULL, NULL, NULL,
    TRUE, 'K', ARG_STRING, 0.0, 0, NULL},
  {"Link\n"
   "      o LinkCDSmRNAbyOverlap\n"
   "      p LinkCDSmRNAbyProduct\n"
   "      r ReassignFeatureIDs\n"
   "      c ClearFeatureIDs", NULL, NULL, NULL,
    TRUE, 'N', ARG_STRING, 0.0, 0, NULL},
  {"Feature\n"
   "      u Remove User Object\n"
   "      d Remove db_xref\n"
   "      r Remove Redundant Gene xref", NULL, NULL, NULL,
    TRUE, 'F', ARG_STRING, 0.0, 0, NULL},
  {"Taxonomy Lookup", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char           app [64], type;
  CleanFlagData  cfd;
  CharPtr        directory, filter, infile, logfile, outfile, results, str, suffix;
  Boolean        remote;
  time_t         runtime, starttime, stoptime;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
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
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  cfd.batch = FALSE;
  cfd.binary = (Boolean) myargs [b_argBinary].intvalue;
  cfd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  cfd.type = 1;

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
      cfd.type = 1;
      cfd.batch = TRUE;
      break;
    default :
      cfd.type = 1;
      break;
  }

  remote = (Boolean) myargs [R_argRemote].intvalue;

  cfd.clean = myargs [K_argClean].strvalue;
  cfd.link = myargs [N_argLink].strvalue;
  cfd.feat = myargs [F_argFeat].strvalue;
  cfd.taxon = (Boolean) myargs [T_argTaxonLookup].intvalue;

  cfd.amp = AsnAllModPtr ();
  cfd.atp_bss = AsnFind ("Bioseq-set");
  cfd.atp_bsss = AsnFind ("Bioseq-set.seq-set");
  cfd.atp_se = AsnFind ("Bioseq-set.seq-set.E");
  cfd.atp_bsc = AsnFind ("Bioseq-set.class");
  cfd.bssp_atp = AsnLinkType (NULL, cfd.atp_bss);

  logfile = (CharPtr) myargs [l_argLogFile].strvalue;
  if (StringDoesHaveText (logfile)) {
    cfd.logfp = FileOpen (logfile, "w");
  }

  if (remote) {
#ifdef INTERNAL_NCBI_CLEANASN
    if (! PUBSEQBioseqFetchEnable ("cleanasn", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
#else
    PubSeqFetchEnable ();
#endif
  }

  starttime = GetSecs ();

  if (StringDoesHaveText (directory)) {

    cfd.results = results;

    DirExplore (directory, NULL, suffix, FALSE, CleanupOneRecord, (Pointer) &cfd);

  } else if (StringDoesHaveText (infile) && StringDoesHaveText (outfile)) {

    cfd.outfile = outfile;

    CleanupOneRecord (infile, (Pointer) &cfd);
  }

  stoptime = GetSecs ();
  runtime = stoptime - starttime;
  if (cfd.logfp != NULL) {
    fprintf (cfd.logfp, "Finished in %ld seconds\n", (long) runtime);
    FileClose (cfd.logfp);
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

