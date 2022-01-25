/*   fetchent.c
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
*  any work or product based on this material.
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
* File Name:  fetchent.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/10/98
*
* $Revision: 6.1 $
*
* File Description: 
*
*   Sample program to demonstrate fetching MEDLINE or Sequence records from
*   Entrez, using the string query evaluation functions of <accutils.h>.  In
*   this format, terms have the term name in double quotes followed by the
*   field name in square brackets.  For example:
*
*     "Perutz MF" [AUTH]
*
*   Field names from all Entrez databases, including nucleotide, protein,
*   genome, and structure are:
*
*     [ACCN], [AFFL], [ALL],  [AUTH], [ECNO], [EDAT], [FKEY], [GENE], [ISS],
*     [JOUR], [KYWD], [LANG], [MAJR], [MDAT], [MESH], [ORGN], [PACC], [PAGE],
*     [PDAT], [PROP], [PROT], [PTYP], [SLEN], [SQID], [SuBH], [SUBS], [TITL],
*     [VOL],  [WORD].
*
*     [*] or [ALL] will search all fields. 
*
*   Operators are:
*
*     & (and), | (or), - (butnot), and : (range).
*
*   A more complicated example is shown below:
*
*     (("glucagon" [WORD] | "insulin" [MESH]) & ("1995" : "1996" [PDAT]))
*
*   At some point in the future, a new Entrez network access API will use
*   strings, not hard-coded numbers, to refer to the database.  For now,
*   the database is passed in as a string (ML, AA, or NT), which map to
*   TYP_ML, TYP_AA, and TYP_NT.
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*/

#include <ncbi.h>
#include <accentr.h>
#include <accutils.h>
#include <objmedli.h>
#include <objsset.h>
#include <objacces.h>
#include <tomedlin.h>
#include <asn2ff.h>
#include <sqnutils.h>

static void ReadPubMedRecords (LinkSetPtr lsp, FILE *fp)

{
  Int4                  count;
  Int2                  num;
  MedlineEntryPtr PNTR  list;  /* see <objmedli.h> */
  MedlineEntryPtr       mep;

  if (lsp == NULL || lsp->num == 0 || lsp->uids == NULL) return;
  list = (MedlineEntryPtr PNTR) MemNew (lsp->num * sizeof (MedlineEntryPtr));
  if (list != NULL) {

    /* EntrezMedlineEntryListGet get a maximum of 32767 records at once */
    num = EntrezMedlineEntryListGet (list, lsp->num, lsp->uids, FALSE);

    for (count = 0; count < num; count++) {
      mep = list [count];
      if (mep != NULL) {
        /* the following call saves the record in traditional MEDLINE format */
        if (MedlineEntryToDataFile (mep, fp)) {
          fprintf (fp, "\n\n");
        }
      }
    }

    for (count = 0; count < lsp->num; count++) {
      list [count] = MedlineEntryFree (list [count]);
    }
    MemFree (list);
  }
}

static void ReadPubSeqRecords (LinkSetPtr lsp, Int2 db, FILE *fp)

{
  Int4              count;
  Uint1             format = TYP_NT;
  Int2              num;
  SeqEntryPtr PNTR  list;  /* see <objsset.h> */
  SeqEntryPtr       sep;

  if (lsp == NULL || lsp->num == 0 || lsp->uids == NULL) return;
  list = (SeqEntryPtr PNTR) MemNew (lsp->num * sizeof (SeqEntryPtr));
  if (list != NULL) {

    /* EntrezSeqEntryListGet get a maximum of 32767 records at once */
    num = EntrezSeqEntryListGet (list, lsp->num, lsp->uids, 0, FALSE);

    if (db == TYP_AA) {
      format = GENPEPT_FMT;
    } else if (db == TYP_NT) {
      format = GENBANK_FMT;
    }

    for (count = 0; count < num; count++) {
      sep = list [count];
      if (sep != NULL) {
        /* the following call saves the record in GenBank or GenPept format */
        if (SeqEntryToFlat (sep, fp, format, RELEASE_MODE)) {
          fprintf (fp, "\n\n");
        }
      }
    }

    for (count = 0; count < lsp->num; count++) {
      list [count] = SeqEntryFree (list [count]);
    }
    MemFree (list);
  }
}

static Int2 ProcessQuery (Int2 db, CharPtr query, FILE *fp)

{
  Int4        count;
  LinkSetPtr  lsp;    /* see <objacces.h> */

  if (query == NULL || fp == NULL) return 1;

  /* check query for proper syntax */
  if (! EntrezTLParseString (query, db, -1, NULL, NULL)) {
    Message (MSG_FATAL, "Query string is not well formed");
    return 1;
  }

  /* calculate number of documents that satisfy the query */
  count = EntrezTLEvalCountString (query, db, -1, NULL, NULL);
  if (count > 32000) {
    Message (MSG_FATAL, "Too many documents");
    return 1;
  }

  /* EntrezTLEvalXString returns a ByteStore that can have > 32767 uids */
  lsp = EntrezTLEvalString (query, db, -1, NULL, NULL);

  if (db == TYP_ML) {
    ReadPubMedRecords (lsp, fp);
  } else if (db == TYP_AA || db == TYP_NT) {
    ReadPubSeqRecords (lsp, db, fp);
  }

  LinkSetFree (lsp);
  return 0;
}

#ifdef NUMARG
#undef NUMARG
#endif
#define NUMARG 3

Args myargs [NUMARG] = {
  {"Database (ML/AA/NT)", "ML", NULL, NULL,
    FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"Entrez Query String", "\"Perutz MF\" [AUTH]", NULL, NULL,
    FALSE, 'q', ARG_STRING, 0.0, 0, NULL},
  {"Output File Name", "stdout", NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
};

static CharPtr databases [] = {
  "ML", "AA", "NT", NULL
};

Int2 Main (void)

{
  Int2     db = -1;
  Int2     i;
  Char     path [PATH_MAX];
  CharPtr  progname;
  FILE     *fp;
  Int2     rsult;

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }

  ProgramPath (path, sizeof (path));
  progname = StringRChr (path, DIRDELIMCHR);
  if (progname != NULL) {
    progname++;
  } else {
    progname = "fetchent";
  }

  /* GetArgs is a portable way of obtaining arguments */
  if (! GetArgs (progname, NUMARG, myargs)) {
    Message (MSG_FATAL, "GetArgs failed");
    return 1;
  }

  /* Map database argument to TYP_XX value */
  for (i = 0; databases [i] != NULL; i++) {
    if (StringICmp (myargs [0].strvalue, databases [i]) == 0) {
      db = i;
    }
  }
  if (db < 0 || db > 2) {
    Message (MSG_FATAL, "Database must be ML, AA, or NT");
    return 1;
  }

  if (! EntrezInit (progname, FALSE, NULL)) {
    Message (MSG_FATAL, "EntrezInit failed");
    return 1;
  }

  fp = FileOpen (myargs [2].strvalue, "w");
  if (fp == NULL) {
    Message (MSG_FATAL, "FileOpen failed");
    return 1;
  }

  rsult = ProcessQuery (db, myargs [1].strvalue, fp);

  FileClose (fp);
  EntrezFini ();
  return rsult;
}

