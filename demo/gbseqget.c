/*   gbseqget.c
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
* File Name:  gbseqget.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   11/4/02
*
* $Revision: 6.3 $
*
* File Description:  Demo to fetch by accession, write GBSet XML
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <objgbseq.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <explore.h>
#include <pmfapi.h>
#include <asn2gnbp.h>

static CharPtr ReadALine (CharPtr str, size_t size, FILE *fp)

{
  Char     ch;
  CharPtr  ptr;
  CharPtr  rsult;

  if (str == NULL || size < 1 || fp == NULL) return NULL;
  *str = '\0';
  rsult = FileGets (str, size, fp);
  if (rsult != NULL) {
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch != '\n' && ch != '\r') {
      ptr++;
      ch = *ptr;
    }
    *ptr = '\0';
  }
  return rsult;
}

static void ProcessAccession (CharPtr accn, XtraPtr extra, Boolean only_new)

{
  Char         ch;
  Int4         gi = 0;
  Char         id [41];
  Boolean      is_numeric = TRUE;
  Int4         newgi = 0;
  CharPtr      ptr;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  Char         tmp [41];
  long         val;

  ptr = accn;
  ch = *ptr;
  while (ch != '\0' && is_numeric) {
    if (! IS_DIGIT (ch)) {
      is_numeric = FALSE;
    }
    ptr++;
    ch = *ptr;
  }

  if (is_numeric) {
    if (sscanf (accn, "%ld", &val) == 1) {
      gi = (Int4) val;
      if (gi < 1) return;
      if (only_new) {
        sip = GetSeqIdForGI (gi);
        if (sip != NULL) {
          SeqIdWrite (sip, tmp, PRINTID_TEXTID_ACC_VER, sizeof (tmp));
          SeqIdFree (sip);
          ptr = StringChr (tmp, '.');
          if (ptr != NULL) {
            *ptr = '\0';
            sip = SeqIdFromAccessionDotVersion (tmp);
            newgi = GetGIForSeqId (sip);
            SeqIdFree (sip);
            if (newgi == gi) return;
          }
        }
      }
    }
  } else {
    sip = SeqIdFromAccessionDotVersion (accn);
    gi = GetGIForSeqId (sip);
    SeqIdFree (sip);
    if (only_new) {
      sip = GetSeqIdForGI (gi);
      if (sip != NULL) {
        SeqIdWrite (sip, id, PRINTID_TEXTID_ACC_VER, sizeof (id));
        SeqIdFree (sip);
        if (StringICmp (accn, id) == 0) return;
      }
    }
  }
  if (gi < 1) return;

  sep = PubSeqSynchronousQuery (gi, 0, 0);
  if (sep == NULL) return;

  SeqEntryToGnbk (sep, NULL, GENBANK_FMT, ENTREZ_MODE, SEGMENT_STYLE,
                  SHOW_FAR_TRANSLATION, LOCK_FAR_COMPONENTS, extra, NULL);

  SeqEntryFree (sep);
}

#define i_argInputFile  0
#define o_argOutputFile 1
#define n_argNewRecords 2

Args myargs [] = {
  {"Input File Name", "stdin", NULL, NULL,
    FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File Name", "stdout", NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"New Records Only", "F", NULL, NULL,
    TRUE, 'n', ARG_BOOLEAN, 0.0, 0, NULL},
};

NLM_EXTERN void AsnPrintNewLine PROTO((AsnIoPtr aip));

Int2 Main (void)

{
  AsnIoPtr    aip;
  AsnTypePtr  atp;
  XtraPtr     extra;
  FILE        *fp;
  GBSeq       gbsq;
  GBSet       gbst;
  Char        line [256];
  Boolean     only_new;
  CharPtr     str;
  Char        xmlbuf [128];
  XtraBlock   xtra;

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
  if (! objgbseqAsnLoad ()) {
    Message (MSG_POSTERR, "objgbseqAsnLoad failed");
    return 1;
  }

  if (! GetArgs ("fetchgbseq", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  fp = FileOpen (myargs [i_argInputFile].strvalue, "r");
  if (fp == NULL) {
    return 1;
  }

  if (GetAppParam ("NCBI", "SETTINGS", "XMLPREFIX", NULL, xmlbuf, sizeof (xmlbuf))) {
    AsnSetXMLmodulePrefix (StringSave (xmlbuf));
  }

  MemSet ((Pointer) &xtra, 0, sizeof (XtraBlock));
  MemSet ((Pointer) &gbsq, 0, sizeof (GBSeq));
  xtra.gbseq = &gbsq;
  aip = AsnIoOpen (myargs [o_argOutputFile].strvalue, "wx");

  if (aip == NULL) {
    Message (MSG_POSTERR, "AsnIoOpen failed");
    FileClose (fp);
    return 1;
  }

  only_new = (Boolean) myargs [n_argNewRecords].intvalue;

  PubSeqFetchEnable ();

  xtra.aip = aip;
  atp = AsnLinkType (NULL, AsnFind ("GBSet"));
  xtra.atp = AsnLinkType (NULL, AsnFind ("GBSet.E"));
  if (atp == NULL || xtra.atp == NULL) {
    Message (MSG_POSTERR, "AsnLinkType or AsnFind failed");
    return 1;
  }
  extra = &xtra;
  MemSet ((Pointer) &gbst, 0, sizeof (GBSet));
  AsnOpenStruct (aip, atp, (Pointer) &gbst);

  str = ReadALine (line, sizeof (line), fp);
  while (str != NULL) {
    if (! StringHasNoText (str)) {
      ProcessAccession (str, extra, only_new);
    }
    str = ReadALine (line, sizeof (line), fp);
  }

  AsnCloseStruct (aip, atp, NULL);
  AsnPrintNewLine (aip);
  AsnIoClose (aip);

  FileClose (fp);

  PubSeqFetchDisable ();

  return 0;
}

