/*   sqn2agp.c
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
* File Name:  sqn2agp.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   3/14/11
*
* $Revision: 1.5 $
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
#include <gather.h>
#include <seqport.h>
#include <tofasta.h>

#define SQN2AGP_APP_VER "1.2"

CharPtr SQN2AGP_APPLICATION = SQN2AGP_APP_VER;

typedef struct s2aflags {
  CharPtr  results;
  Int2     known;
  Boolean  unknown;
  FILE     *afp;
  FILE     *ffp;
} S2AFlagData, PNTR S2AFlagPtr;

static void DoOneBioseq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Int4         cumulative = 0;
  DeltaSeqPtr  dsp;
  Int4         gap_sizes [2];
  Char         id [64];
  Int4         len;
  SeqLitPtr    lit;
  Int4         part = 0;
  Int4         seg = 0;
  S2AFlagPtr   sfp;

  if (bsp == NULL) return;
  sfp = (S2AFlagPtr) userdata;
  if (sfp == NULL) return;

  if (sfp->afp == NULL || sfp->ffp == NULL) return;

  if (! ISA_na (bsp->mol)) return;

  if (sfp->unknown) {
    gap_sizes [0] = 100;
  } else {
    gap_sizes [0] = 0;
  }
  gap_sizes [1] = -(sfp->known);

  ConvertNsToGaps (bsp, (Pointer) gap_sizes);

  SeqIdWrite (bsp->id, id, PRINTID_REPORT, sizeof (id) - 1);

  if (bsp->repr == Seq_repr_raw) {
    seg++;
    fprintf (sfp->ffp, ">%s_%ld\n", id, (long) seg);
    BioseqFastaStream (bsp, sfp->ffp, 0, 60, 0, 0, FALSE);

    fprintf (sfp->afp, "%s\t1\t%ld\t1\tW\t%s_%ld\t1\t%ld\t+\n", id,
             (long) bsp->length,
             id, (long) seg, (long) bsp->length);

    return;
  }

  if (bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return;

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next, cumulative += len) {

    len = 0;

    if (dsp->choice != 2) continue;
    lit = (SeqLitPtr) dsp->data.ptrvalue;
    if (lit == NULL) continue;

    if (lit->length < 1) {
      /* unknown length */
      continue;
    }

    len = lit->length;

    if (sfp->unknown && len == 100) {
      /* designated unknown length */
      part++;
      fprintf (sfp->afp, "%s\t%ld\t%ld\t%ld\tU\t%ld\tfragment\tyes\t\n", id,
               (long) cumulative + 1, (long) cumulative + len,
               (long) part, (long) len);
      continue;
    }

    if (lit->seq_data == NULL || lit->seq_data_type == Seq_code_gap) {
      /* known length */
      part++;
      fprintf (sfp->afp, "%s\t%ld\t%ld\t%ld\tN\t%ld\tfragment\tyes\t\n", id,
               (long) cumulative + 1, (long) cumulative + len,
               (long) part, (long) len);
      continue;
    }

    seg++;
    fprintf (sfp->ffp, ">%s_%ld\n", id, (long) seg);
    SeqLitFastaStream (lit, sfp->ffp, 0, 60, 0, 0);

    part++;
    fprintf (sfp->afp, "%s\t%ld\t%ld\t%ld\tW\t%s_%ld\t1\t%ld\t+\n", id,
             (long) cumulative + 1, (long) cumulative + len,
             (long) part, id, (long) seg, (long) len);
  }
}

static void DoWriteAgpAndFsa (
  SeqEntryPtr sep,
  Uint2 entityID,
  S2AFlagPtr sfp,
  CharPtr agppath,
  CharPtr fsapath
)

{
  if (sep == NULL || sfp == NULL) return;

  sfp->afp = FileOpen (agppath, "w");
  sfp->ffp = FileOpen (fsapath, "w");

  VisitBioseqsInSep (sep, (Pointer) sfp, DoOneBioseq);

  FileClose (sfp->afp);
  FileClose (sfp->ffp);
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  Char          agppath [PATH_MAX], fsapath [PATH_MAX];
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype, entityID = 0;
  FILE          *fp;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  S2AFlagPtr    sfp;

  if (StringHasNoText (filename)) return;
  sfp = (S2AFlagPtr) userdata;
  if (sfp == NULL) return;

  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    Message (MSG_POSTERR, "Failed to open '%s'", filename);
    return;
  }

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

  FileClose (fp);

  entityID = ObjMgrRegister (datatype, dataptr);

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

      ptr = StringRChr (filename, '.');
      if (ptr != NULL) {
        *ptr = '\0';
      }

      agppath [0] = '\0';
      fsapath [0] = '\0';
      if (StringDoesHaveText (sfp->results)) {

        ptr = StringRChr (filename, DIRDELIMCHR);
        if (ptr != NULL) {
          StringNCpy_0 (agppath, sfp->results, sizeof (agppath));
          StringNCpy_0 (fsapath, sfp->results, sizeof (fsapath));
          ptr++;
          filename = ptr;
        }
      }
      FileBuildPath (agppath, NULL, filename);
      FileBuildPath (fsapath, NULL, filename);
      StringCat (agppath, ".agp");
      StringCat (fsapath, ".contigs.fsa");

      sep = GetTopSeqEntryForEntityID (entityID);
      if (sep != NULL) {

        AssignIDsInEntity (entityID, 0, NULL);
        DoWriteAgpAndFsa (sep, entityID, sfp, agppath, fsapath);

      }

      ObjMgrFreeByEntityID (entityID);
    }

  } else {

    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }
}

/* Args structure contains command-line arguments */

typedef enum {
  p_argInputPath = 0,
  r_argOutputPath,
  i_argInputFile,
  f_argFilter,
  x_argSuffix,
  k_argKnown,
  u_argUnknown,
} Arguments;

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".sqn", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Known gap length",  "20", NULL, NULL,
    TRUE, 'k', ARG_INT, 0.0, 0, NULL},
  {"Unknown gap length 100",  "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char         app [64];
  CharPtr      directory, filter, infile, results, suffix;
  S2AFlagData  sfd;

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

  sprintf (app, "sqn2agp %s", SQN2AGP_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &sfd, 0, sizeof (S2AFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = directory;
  }
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  sfd.known = (Int2) myargs [k_argKnown].intvalue;
  sfd.unknown = (Boolean) myargs [u_argUnknown].intvalue;

  if (StringDoesHaveText (directory)) {

    sfd.results = results;
    DirExplore (directory, filter, suffix, FALSE, ProcessOneRecord, (Pointer) &sfd);

  } else if (StringDoesHaveText (infile)) {

    sfd.results = results;
    ProcessOneRecord (infile, (Pointer) &sfd);
  }

  return 0;
}

