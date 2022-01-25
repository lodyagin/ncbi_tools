/*   gi2accn.c
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
* File Name:  gi2accn.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/15/02
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
#include <objfdef.h>
#include <sequtil.h>
#include <gather.h>
#include <sqnutils.h>
#include <subutil.h>
#include <explore.h>
#include <pmfapi.h>

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
  SeqIdStripLocus (newsip);
  sip->choice = newsip->choice;
  sip->data.ptrvalue = newsip->data.ptrvalue;
  newsip->choice = SEQID_NOT_SET;
  newsip->data.ptrvalue = NULL;
  SeqIdFree (newsip);
}

static void UpdateAligns (
  SeqAlignPtr sap,
  Pointer userdata
)

{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  SeqIdPtr      sip;
  SeqLocPtr     slp = NULL;
  StdSegPtr     ssp;

  if (sap == NULL) return;

  if (sap->bounds != NULL) {
    sip = SeqLocId (sap->bounds);
    ConvertGiToAccn (sip);
  }
  if (sap->segs == NULL) return;

  switch (sap->segtype) {
    case SAS_DENDIAG :
      ddp = (DenseDiagPtr) sap->segs;
      if (ddp != NULL) {
        for (sip = ddp->id; sip != NULL; sip = sip->next) {
          ConvertGiToAccn (sip);
        }
      }
      break;
    case SAS_DENSEG :
      dsp = (DenseSegPtr) sap->segs;
      if (dsp != NULL) {
        for (sip = dsp->ids; sip != NULL; sip = sip->next) {
          ConvertGiToAccn (sip);
        }
      }
      break;
    case SAS_STD :
      ssp = (StdSegPtr) sap->segs;
      for (slp = ssp->loc; slp != NULL; slp = slp->next) {
        sip = SeqLocId (slp);
        ConvertGiToAccn (sip);
      }
      break;
    case SAS_DISC :
      /* recursive */
      for (sap = (SeqAlignPtr) sap->segs; sap != NULL; sap = sap->next) {
        UpdateAligns (sap, userdata);
      }
      break;
    default :
      break;
  }
}

/* Args structure contains command-line arguments */

#define i_argInputFile  0
#define o_argOutputFile 1

Args myargs [] = {
  {"Input File", "stdin", NULL, NULL,
    FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File", "stdout", NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
};

Int2 Main (void)

{
  AsnIoPtr     aip;
  CharPtr      infile, outfile;
  SeqEntryPtr  sep;

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

  if (! GetArgs ("gi2accn", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;

  aip = AsnIoOpen (infile, "r");
  if (aip == NULL) {
    Message (MSG_FATAL, "AsnIoOpen failed");
    return 1;
  }

  sep = SeqEntryAsnRead (aip, NULL);
  AsnIoClose (aip);
  if (sep == NULL) {
    Message (MSG_FATAL, "SeqEntryAsnRead failed");
    return 1;
  }

  PubSeqFetchEnable ();

  LookupFarSeqIDs (sep, FALSE, FALSE, FALSE, TRUE, FALSE);

  VisitAlignmentsInSep (sep, NULL, UpdateAligns);

  PubSeqFetchDisable ();

  BasicSeqEntryCleanup (sep);

  aip = AsnIoOpen (outfile, "w");
  if (aip == NULL) {
    Message (MSG_FATAL, "AsnIoOpen failed");
    return 1;
  }

  SeqEntryAsnWrite (sep, aip, NULL);
  AsnIoClose (aip);

  return 0;
}

