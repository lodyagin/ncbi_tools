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

#define CLEANASN_APP_VER "1.1"

CharPtr CLEANASN_APPLICATION = CLEANASN_APP_VER;

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
             StringDoesHaveText (grp->locus_tag)) {
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

static void CleanupOneRecord (
  CharPtr directory,
  CharPtr results,
  CharPtr filename,
  CharPtr clean,
  CharPtr link,
  CharPtr feat,
  Boolean taxon
)

{
  AsnIoPtr     aip;
  Pointer      dataptr;
  Uint2        datatype;
  Uint2        entityID;
  FILE*        fp;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;

  Message (MSG_POST, "%s\n", filename);

  StringNCpy_0 (path, directory, sizeof (path));
  FileBuildPath (path, NULL, filename);

  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE,
                                    FALSE, TRUE, FALSE);

  FileClose (fp);

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep != NULL) {

    if (StringChr (clean, 'b') != NULL) {
      BasicSeqEntryCleanup (sep);
    }
    if (StringChr (clean, 's') != NULL) {
      SeriousSeqEntryCleanup (sep, NULL, NULL);
    }

    if (taxon) {
      Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
    }

    if (StringChr (link, 'o') != NULL) {
      SeqMgrIndexFeatures (entityID, 0);
      LinkCDSmRNAbyOverlap (sep);
    }
    if (StringChr (link, 'p') != NULL) {
      SeqMgrIndexFeatures (entityID, 0);
      LinkCDSmRNAbyProduct (sep);
    }
    if (StringChr (link, 'r') != NULL) {
      SeqMgrIndexFeatures (entityID, 0);
      ReassignFeatureIDs (sep);
    }

    if (StringChr (feat, 'u') != NULL) {
      VisitFeaturesInSep (sep, NULL, RemoveFeatUser);
    }
    if (StringChr (feat, 'd') != NULL) {
      VisitFeaturesInSep (sep, NULL, RemoveFeatDbxref);
    }
    if (StringChr (feat, 'r') != NULL) {
      SeqMgrIndexFeatures (entityID, 0);
      VisitFeaturesInSep (sep, NULL, RemoveUnnecGeneXref);
    }

    StringNCpy_0 (path, results, sizeof (path));
    FileBuildPath (path, NULL, filename);

    aip = AsnIoOpen (path, "w");
    if (aip != NULL) {
      if (datatype == OBJ_SEQSUB) {
        SeqSubmitAsnWrite ((SeqSubmitPtr) dataptr, aip, NULL);
      } else {
        SeqEntryAsnWrite (sep, aip, NULL);
      }
      AsnIoFlush (aip);
      AsnIoClose (aip);
    }
  }

  ObjMgrFreeByEntityID (entityID);
}

/* Args structure contains command-line arguments */

#define p_argInputPath      0
#define r_argOutputPath     1
#define c_argClean          2
#define l_argLink           3
#define f_argFeat           4
#define t_argTaxonLookup    5
#define R_argRemote         6

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Cleanup (b BasicSeqEntryCleanup, s SeriousSeqEntryCleanup)", NULL, NULL, NULL,
    TRUE, 'c', ARG_STRING, 0.0, 0, NULL},
  {"Link (o LinkCDSmRNAbyOverlap, p LinkCDSmRNAbyProduct, r ReassignFeatureIDs)", NULL, NULL, NULL,
    TRUE, 'l', ARG_STRING, 0.0, 0, NULL},
  {"Feature (u Remove User Object, d Remove db_xref, r Remove Redundant Gene xref)", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"Taxonomy Lookup", "F", NULL, NULL,
    TRUE, 't', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char        app [64];
  CharPtr     clean, feat, link;
  CharPtr     directory, results;
  ValNodePtr  head, vnp;
  Boolean     remote, taxon;

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

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  if (StringHasNoText (directory)) {
    Message (MSG_FATAL, "You must supply an input directory (-p).\nUse -p . to specify the current directory.\n\n");
    return 1;
  }
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = directory;
  }

  clean = myargs [c_argClean].strvalue;
  link = myargs [l_argLink].strvalue;
  feat = myargs [f_argFeat].strvalue;

  taxon = (Boolean) myargs [t_argTaxonLookup].intvalue;
  remote = (Boolean) myargs [R_argRemote].intvalue;

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

  head = DirCatalog (directory);
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      CleanupOneRecord (directory, results,
                        (CharPtr) vnp->data.ptrvalue,
                        clean, link, feat, taxon);
    }
  }
  ValNodeFreeData (head);

  if (remote) {
#ifdef INTERNAL_NCBI_CLEANASN
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
  }

  return 0;
}

